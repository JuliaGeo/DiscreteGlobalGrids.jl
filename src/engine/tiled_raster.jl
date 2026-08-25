# A grid whose cells are the pixels of raster tiles gets a tree in two layers: a
# packed tree over the tiles' caps, and a bisection quadtree inside each tile.
# Leaf indices are the grid's own indices at every node.

"""
    LEAF_CELLS

Cells at or below which a raster node stops splitting and yields its cells
directly. It is the capacity of [`LeafCells`](@ref).
"""
const LEAF_CELLS = 9

"""
    LeafCells(entries, len)

A leaf's cells as `(grid index, cap)` pairs, in a fixed inline buffer.

The return value of [`STI.child_indices_extents`](@ref
GeometryOps.SpatialTreeInterface.child_indices_extents) for a raster tree: a
read-only `AbstractVector` that indexes, iterates and `collect`s like the
`Vector` it replaces, but is `isbits`, so it lives in the caller's frame and
never reaches the heap.

[`STI.isleaf`](@ref GeometryOps.SpatialTreeInterface.isleaf) splits any node
holding more than [`LEAF_CELLS`](@ref) cells, so the entries fit an
`NTuple{LEAF_CELLS}`. Every call must return its own value rather than share a
buffer: callers keep the result past the call, and a dual-tree self-join loops
over two leaves' entries at once.

The unused tail repeats entry one; an empty leaf repeats a placeholder no
`size` reaches. Entry order is row-major: the column index varies fastest.
"""
struct LeafCells <: AbstractVector{Tuple{Int,Cap}}
    entries::NTuple{LEAF_CELLS,Tuple{Int,Cap}}
    len::Int
end

Base.size(e::LeafCells) = (e.len,)
Base.IndexStyle(::Type{LeafCells}) = IndexLinear()

Base.@propagate_inbounds function Base.getindex(e::LeafCells, k::Int)
    @boundscheck checkbounds(e, k)
    return @inbounds e.entries[k]
end

# The placeholder an empty leaf's unused slots hold; no `size` reaches it.
const _EMPTY_CAP = SphericalCap(USPoint(0.0, 0.0, 1.0), 0.0)

# Build a leaf from a callable that answers entry `k` of the row-major order.
@inline function leaf_cells(entry, n::Int)
    0 <= n <= LEAF_CELLS || throw(ArgumentError(
        "a leaf holds 0 to $LEAF_CELLS cells, not $n"))
    n == 0 && return LeafCells(ntuple(_ -> (0, _EMPTY_CAP), Val(LEAF_CELLS)), 0)
    head = entry(1)
    entries = ntuple(Val(LEAF_CELLS)) do k
        k == 1 ? head : (k <= n ? entry(k) : head)
    end
    return LeafCells(entries, n)
end

# Part `t` of a non-empty, near-equal partition of `lo:lo+len-1`.
@inline function rect_part(lo::Int, len::Int, parts::Int, t::Int)
    base, rem = divrem(len, parts)
    off = (t - 1) * base + min(t - 1, rem)
    return (lo + off, lo + off + base + (t <= rem ? 1 : 0) - 1)
end

# How many parts each axis of a `(nrows, ncols)` rectangle splits into: the
# longer axis in two, and neither when the rectangle is one cell wide and tall.
@inline bisect_parts(nr::Int, nc::Int) = nr >= nc ? (min(2, nr), 1) : (1, min(2, nc))

# ===========================================================================
# The packed tree over the tiles
# ===========================================================================

"""
    RASTER_TILE_ARITY

Children of a packed tile node. The tiles are sorted by the Morton key of their
cap centres and split into this many near-equal blocks per level, down to one
tile per node, which is where the per-tile quadtree begins. Four is
[`IndexTree`](@ref)'s arity, the other packed tree here.

The alternative is a flat root — one node with every tile as a child — and it
measures the same on the tile counts a holding of tens has, for a point query
and for a dual-tree join alike, because a shallow packed tree tests about as
many caps as there are tiles. It diverges as the holding grows: a flat root
tests every tile on every query, so a 1024-tile holding's point query costs
3.4 times a packed root's, and the ratio grows with the tile count. A global
Copernicus holding is 26 450 tiles.
"""
const RASTER_TILE_ARITY = 4

"""
    RasterTileTree(grid, inds; arity = RASTER_TILE_ARITY)

The tile layer of the tree [`treeify`](@ref) builds for a grid that answers
[`raster_tiles`](@ref): the tiles covering grid indices `inds`, sorted by the
Morton key of their cap centres, and a stored cap and pixel count per node.

Node `1` is the root. A node holding one tile is where the packed tree stops and
the tile's own quadtree begins, so its stored cap is that tile's.

Every cap is a merge of the caps below it, so the packed layer's extents nest.
"""
struct RasterTileTree{G<:AbstractGrid,T}
    grid::G
    tiles::Vector{T}
    tilecaps::Vector{Cap}
    shapes::Vector{Tuple{Int,Int}}
    node_lo::Vector{Int}            # node -> first tile slot
    node_hi::Vector{Int}            # node -> last tile slot
    node_children::Vector{Vector{Int}}
    node_cap::Vector{Cap}
    node_weight::Vector{Int}        # pixels beneath the node
end

function RasterTileTree(grid::AbstractGrid, inds::AbstractUnitRange{<:Integer};
        arity::Integer = RASTER_TILE_ARITY)
    Int(arity) >= 2 || throw(ArgumentError("arity must be at least 2, got $arity"))
    listed = raster_tiles(grid, inds)
    listed === nothing && throw(ArgumentError(
        "$(typeof(grid)) does not answer `raster_tiles`, so it has no tiled raster tree"))
    tiles = collect(listed)
    n = length(tiles)
    caps = Vector{Cap}(undef, n)
    shapes = Vector{Tuple{Int,Int}}(undef, n)
    keys = Vector{UInt64}(undef, n)
    for k in 1:n
        nrows, ncols = raster_shape(grid, tiles[k])
        (nrows > 0 && ncols > 0) || throw(ArgumentError(
            "tile $(tiles[k]) reports a $(nrows)x$(ncols) raster; both must be positive"))
        shapes[k] = (Int(nrows), Int(ncols))
        cap = raster_cap(grid, tiles[k], 0, Int(nrows) - 1, 0, Int(ncols) - 1)
        caps[k] = cap
        keys[k] = _morton_key(cap.point)
    end
    order = sortperm(keys)
    tree = RasterTileTree{typeof(grid),eltype(tiles)}(grid, tiles[order], caps[order],
        shapes[order], Int[], Int[], Vector{Int}[], Cap[], Int[])
    _build_tile_node!(tree, 1, n, Int(arity))
    return tree
end

# Depth-first recursion records children by node index, as `IndexTree` does.
function _build_tile_node!(tree::RasterTileTree, lo::Int, hi::Int, arity::Int)
    index = length(tree.node_lo) + 1
    push!(tree.node_lo, lo)
    push!(tree.node_hi, hi)
    push!(tree.node_children, Int[])
    push!(tree.node_cap, _EMPTY_CAP)
    push!(tree.node_weight, 0)
    count = hi - lo + 1
    if count <= 0
        tree.node_cap[index] = _EMPTY_CAP
        return index
    elseif count == 1
        tree.node_cap[index] = tree.tilecaps[lo]
        tree.node_weight[index] = prod(tree.shapes[lo])
        return index
    end
    per = cld(count, arity)
    start = lo
    cap = nothing
    weight = 0
    while start <= hi
        stop = min(hi, start + per - 1)
        child = _build_tile_node!(tree, start, stop, arity)
        push!(tree.node_children[index], child)
        cap = cap === nothing ? tree.node_cap[child] :
              Extents.union(cap, tree.node_cap[child])
        weight += tree.node_weight[child]
        start = stop + 1
    end
    tree.node_cap[index] = cap === nothing ? _EMPTY_CAP : cap
    tree.node_weight[index] = weight
    return index
end

Base.show(io::IO, tree::RasterTileTree) =
    print(io, "RasterTileTree(", typeof(tree.grid).name.name, ", ntiles=",
        length(tree.tiles), ", nodes=", length(tree.node_lo), ")")

# ===========================================================================
# The cursor
# ===========================================================================

"""
    TiledRasterCursor(tree)

`GeometryOps.SpatialTreeInterface` cursor over a grid whose cells are the pixels
of raster tiles. Prefer [`treeify(grid)`](@ref treeify) to direct construction.

A node is either

  - a **tile node**: a node of the packed [`RasterTileTree`](@ref), covering
    every pixel of a block of tiles, or
  - a **raster node**: a rectangle of rows and columns inside one tile.

Descent runs through the packed tree to a single tile and then bisects that
tile's rectangle — the longer axis in two — until a node holds
[`LEAF_CELLS`](@ref) pixels or fewer, which is a leaf. A leaf names its pixels
by [`raster_localindex`](@ref), so leaf indices are the grid's own indices and a
tile's pixels are one contiguous block of them.

Tile node extents are stored and nest. Raster node extents are derived from
[`raster_cap`](@ref) and memoized per task; like every
[`node_extent`](@ref) here they bound cell *geometry*, so a node's cap covers
the boundaries of the pixels beneath it but not necessarily their caps.
"""
struct TiledRasterCursor{G<:AbstractGrid,T}
    tree::RasterTileTree{G,T}
    node::Int       # packed node index; 0 inside a tile
    slot::Int       # tile slot; 0 at a multi-tile node
    j0::Int
    j1::Int
    i0::Int
    i1::Int         # the rectangle inside tile `slot`, inclusive
    inraster::Bool
end

TiledRasterCursor(tree::RasterTileTree) =
    TiledRasterCursor(tree, 1, 0, 0, 0, 0, 0, false)

function Base.show(io::IO, c::TiledRasterCursor)
    if c.inraster
        print(io, "TiledRasterCursor(rows ", c.j0, ":", c.j1, " x cols ", c.i0, ":",
            c.i1, " of tile ", c.slot, "/", length(c.tree.tiles), ")")
    else
        lo, hi = c.tree.node_lo[c.node], c.tree.node_hi[c.node]
        print(io, "TiledRasterCursor(tiles ", lo, ":", hi, " of ",
            length(c.tree.tiles), ")")
    end
end

Base.show(io::IO, ::MIME"text/plain", c::TiledRasterCursor) = show(io, c)

# What a single-tile node's children partition: the parts each axis splits into,
# the tile, and the rectangle. A node holding one tile stands for that tile's
# whole raster, so the tree spends no level on the tile itself. Zero parts mean
# the node is a leaf. Every accessor below reads this once.
@inline function _raster_split(c::TiledRasterCursor)
    if c.inraster
        slot, j0, j1, i0, i1 = c.slot, c.j0, c.j1, c.i0, c.i1
    else
        slot = @inbounds c.tree.node_lo[c.node]
        nrows, ncols = @inbounds c.tree.shapes[slot]
        j0, j1, i0, i1 = 0, nrows - 1, 0, ncols - 1
    end
    (j1 - j0 + 1) * (i1 - i0 + 1) <= LEAF_CELLS &&
        return (0, 0, slot, j0, j1, i0, i1)
    pr, pc = bisect_parts(j1 - j0 + 1, i1 - i0 + 1)
    # A stalled split must terminate.
    pr * pc <= 1 && return (0, 0, slot, j0, j1, i0, i1)
    return (pr, pc, slot, j0, j1, i0, i1)
end

@inline _issingletile(c::TiledRasterCursor) =
    c.inraster || @inbounds(c.tree.node_lo[c.node]) == @inbounds(c.tree.node_hi[c.node])

@inline _isempty(c::TiledRasterCursor) =
    !c.inraster && @inbounds(c.tree.node_hi[c.node]) < @inbounds(c.tree.node_lo[c.node])

STI.isspatialtree(::Type{<:TiledRasterCursor}) = true

# Tile extents are stored and raster extents memoized, so a hit is a load and
# the search should not pay a vector per visited node to avoid one.
STI.node_extent_is_expensive(::Type{<:TiledRasterCursor}) = false

function STI.isleaf(c::TiledRasterCursor)
    _isempty(c) && return true
    _issingletile(c) || return false
    return first(_raster_split(c)) == 0
end

function STI.nchild(c::TiledRasterCursor)
    _isempty(c) && return 0
    _issingletile(c) || return length(@inbounds c.tree.node_children[c.node])
    pr, pc, _, _, _, _, _ = _raster_split(c)
    return pr * pc
end

function STI.getchild(c::TiledRasterCursor, k::Int)
    if !_issingletile(c)
        kids = @inbounds c.tree.node_children[c.node]
        1 <= k <= length(kids) || throw(BoundsError(c, k))
        return TiledRasterCursor(c.tree, @inbounds(kids[k]), 0, 0, 0, 0, 0, false)
    end
    pr, pc, slot, j0, j1, i0, i1 = _raster_split(c)
    1 <= k <= pr * pc || throw(BoundsError(c, k))
    tr, tc = (k - 1) ÷ pc + 1, (k - 1) % pc + 1
    a0, a1 = rect_part(j0, j1 - j0 + 1, pr, tr)
    b0, b1 = rect_part(i0, i1 - i0 + 1, pc, tc)
    return TiledRasterCursor(c.tree, 0, slot, a0, a1, b0, b1, true)
end

STI.getchild(c::TiledRasterCursor) = (STI.getchild(c, k) for k in 1:STI.nchild(c))

STI.node_extent(c::TiledRasterCursor) =
    c.inraster ? _memo_extent(c) : @inbounds c.tree.node_cap[c.node]

"""
    STI.child_indices_extents(cursor::TiledRasterCursor) -> LeafCells

A leaf's pixels as `(grid index, cap)` pairs. See [`LeafCells`](@ref).
"""
function STI.child_indices_extents(c::TiledRasterCursor)
    STI.isleaf(c) ||
        throw(ArgumentError("child_indices_extents is only valid for leaf nodes"))
    _isempty(c) && return leaf_cells(k -> (0, _EMPTY_CAP), 0)
    tree = c.tree
    _, _, slot, j0, j1, i0, i1 = _raster_split(c)
    tile = @inbounds tree.tiles[slot]
    nfast = i1 - i0 + 1
    return leaf_cells(nfast * (j1 - j0 + 1)) do k
        slow, fast = divrem(k - 1, nfast)
        j = j0 + slow
        i = i0 + fast
        return (raster_localindex(tree.grid, tile, j, i),
            raster_cap(tree.grid, tile, j, j, i, i))
    end
end

# ===========================================================================
# ConservativeRegridding.Trees
# ===========================================================================

GOCore.best_manifold(c::TiledRasterCursor) = GOCore.best_manifold(c.tree.grid)

# Leaf indices are grid indices at every node, so these answer about the whole
# grid regardless of which node holds them.
Trees.ncells(c::TiledRasterCursor) = ncells(c.tree.grid)
Trees.getcell(c::TiledRasterCursor, i::Int) = getcell(c.tree.grid, i)
Trees.getcell(c::TiledRasterCursor) = getcell(c.tree.grid)

# `Trees.ncells` answers for the whole grid, so the frontier's default estimate
# would be wrong here; a node knows its own pixel count.
function Trees.split_weight(c::TiledRasterCursor)
    c.inraster && return (c.j1 - c.j0 + 1) * (c.i1 - c.i0 + 1)
    return @inbounds c.tree.node_weight[c.node]
end

# ===========================================================================
# The per-task extent memo
# ===========================================================================

"""
    RasterExtentMemo()

One task's direct-mapped cache of derived raster node extents, so a repeat ask
is a key compare and a load rather than a fresh [`raster_cap`](@ref).

A node's rectangle hashes to exactly one slot, and a miss derives the extent and
overwrites whatever sat there. Hence bounded memory per task whatever the tile
size, no eviction policy, no lock — tasks sharing a tree have separate tables —
and a collision costs a re-derive, never a wrong extent. The hit rate comes from
revisits: the dual-tree join asks each node's extent once per opposing node.

The slots describe one tree, and are cleared whenever the task turns to another.
Tile nodes are not cached: their extents are stored in the tree itself.
"""
mutable struct RasterExtentMemo
    tree::Any
    const keys::Vector{NTuple{5,Int}}
    const vals::Vector{Cap}
end

# A power of two, so the slot is a mask rather than a modulo.
const _RASTER_MEMO_SLOTS = 1024

# No node has a negative tile slot, so this can never equal a real key.
const _NO_RASTER_NODE = (-1, -1, -1, -1, -1)

RasterExtentMemo() = RasterExtentMemo(nothing,
    fill(_NO_RASTER_NODE, _RASTER_MEMO_SLOTS),
    Vector{Cap}(undef, _RASTER_MEMO_SLOTS))

function _taskrastermemo(tree)
    memo = get!(RasterExtentMemo, task_local_storage(),
        :_dgg_raster_extent_memo)::RasterExtentMemo
    if memo.tree !== tree
        fill!(memo.keys, _NO_RASTER_NODE)
        memo.tree = tree
    end
    return memo
end

# The slot holds the whole key and compares it, so a collision is a miss that
# overwrites, never a wrong extent.
function _memo_extent(c::TiledRasterCursor)
    key = (c.slot, c.j0, c.j1, c.i0, c.i1)
    memo = _taskrastermemo(c.tree)
    s = Int(hash(key) & UInt(_RASTER_MEMO_SLOTS - 1)) + 1
    @inbounds memo.keys[s] == key && return @inbounds memo.vals[s]
    tile = @inbounds c.tree.tiles[c.slot]
    extent = raster_cap(c.tree.grid, tile, c.j0, c.j1, c.i0, c.i1)
    @inbounds memo.keys[s] = key
    @inbounds memo.vals[s] = extent
    return extent
end

# ===========================================================================
# treeify and subcursor
# ===========================================================================

"""
    tiled_raster_tree(grid, inds = 1:ncells(grid); arity = RASTER_TILE_ARITY)

The [`TiledRasterCursor`](@ref) root over the cells at grid indices `inds`, or
`nothing` when `grid` answers no [`raster_tiles`](@ref) for them.

Leaf indices are `grid`'s own, so a window's tree names its cells exactly as the
whole grid's does, which is what [`subcursor`](@ref) needs.
"""
function tiled_raster_tree(grid::AbstractGrid,
        inds::AbstractUnitRange{<:Integer} = 1:ncells(grid);
        arity::Integer = RASTER_TILE_ARITY)
    raster_tiles(grid, inds) === nothing && return nothing
    return TiledRasterCursor(RasterTileTree(grid, inds; arity))
end

treeify(::GOCore.Manifold, c::TiledRasterCursor) = c
treeify(c::TiledRasterCursor) = c
