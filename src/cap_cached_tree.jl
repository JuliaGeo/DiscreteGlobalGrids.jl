# A chunk's tree with per-cell work done once, not once per visit.

"""
    CapCachedTree(node, caps)

Wrap a grid cursor so a leaf's extent is `caps[p]` for grid index `p`, rather
than re-derived from the cell boundary on every visit. `caps` is indexed by raw
grid index; a chunk's vector covers only its own window and is wrapped in
`_ShiftedCaps`.
"""
struct CapCachedTree{C<:HierarchicalGridCursor,V<:AbstractVector}
    node::C
    caps::V
end

"""
    _ShiftedCaps(data, offset)

`data` addressed by grid index: element `i` is `data[i - offset]`, with axes
`offset .+ axes(data, 1)`.
"""
struct _ShiftedCaps{E} <: AbstractVector{E}
    data::Vector{E}
    offset::Int
end

Base.size(v::_ShiftedCaps) = size(v.data)
Base.axes(v::_ShiftedCaps) = (v.offset .+ axes(v.data, 1),)
Base.IndexStyle(::Type{<:_ShiftedCaps}) = IndexLinear()
Base.parent(v::_ShiftedCaps) = v.data
Base.@propagate_inbounds function Base.getindex(v::_ShiftedCaps, i::Int)
    @boundscheck checkbounds(v, i)
    return @inbounds v.data[i - v.offset]
end

Base.show(io::IO, t::CapCachedTree) =
    print(io, "CapCachedTree(", t.node, ", ", length(t.caps), " caps)")

STI.isspatialtree(::Type{<:CapCachedTree}) = true
STI.node_extent_is_expensive(::Type{<:CapCachedTree{C}}) where {C} =
    STI.node_extent_is_expensive(C)
STI.isleaf(t::CapCachedTree) = STI.isleaf(t.node)
STI.nchild(t::CapCachedTree) = STI.nchild(t.node)
STI.getchild(t::CapCachedTree) =
    Iterators.map(n -> CapCachedTree(n, t.caps), STI.getchild(t.node))
STI.getchild(t::CapCachedTree, i::Int) =
    CapCachedTree(STI.getchild(t.node, i), t.caps)

# A window node at the leaf level is exactly one stored cell; its index locates
# it in the cache. Everything else keeps the cursor's own extent logic.
function STI.node_extent(t::CapCachedTree)
    c = t.node
    if !Engine._issynthetic(c) && c.level >= c.leaf_level
        return @inbounds t.caps[c.first_index]
    end
    return STI.node_extent(c)
end

"""
    _CapCachedLeafEntries(cursor, caps)

A read-only view of one [`CapCachedTree`](@ref) leaf's `(index, cap)` entries.

The cached tree deliberately groups up to [`_CACHED_BUCKET_SIZE`](@ref) cells in
one leaf. A dual-tree search may pair that leaf with many source leaves, so
materialising its entries on every pairing turns the cache into hundreds of
megabytes of short-lived vectors. This view keeps the entries lazy: their
indices come from the cursor and their already-computed caps come from the
cache, with no per-visit buffer.
"""
struct _CapCachedLeafEntries{C,V<:AbstractVector,E} <: AbstractVector{Tuple{Int,E}}
    cursor::C
    caps::V
end

_CapCachedLeafEntries(cursor::C, caps::V) where {C,E,V<:AbstractVector{E}} =
    _CapCachedLeafEntries{C,V,E}(cursor, caps)

Base.size(entries::_CapCachedLeafEntries) = (Engine._stored_count(entries.cursor),)
Base.IndexStyle(::Type{<:_CapCachedLeafEntries}) = IndexLinear()

Base.@propagate_inbounds function Base.getindex(entries::_CapCachedLeafEntries, k::Int)
    @boundscheck checkbounds(entries, k)
    index = Engine._stored_index(entries.cursor, k)
    return (index, @inbounds entries.caps[index])
end

function STI.child_indices_extents(t::CapCachedTree)
    c = t.node
    STI.isleaf(c) ||
        throw(ArgumentError("child_indices_extents is only valid for leaf nodes"))
    return _CapCachedLeafEntries(c, t.caps)
end

getcell(t::CapCachedTree, i::Int) = getcell(t.node, i)
getcell(t::CapCachedTree) = getcell(t.node)
split_weight(t::CapCachedTree) = split_weight(t.node)
subtree_indices(t::CapCachedTree) = subtree_indices(t.node)
GOCore.best_manifold(t::CapCachedTree) = GOCore.best_manifold(t.node)

"""
    _CACHED_BUCKET_SIZE

Leaf size for a cursor whose cell caps are precomputed: stop descent at `7^2`
stored cells, two IGeo7 refinement levels above the grid's own resolution. It
is a cell budget, not a level count, so a system of another aperture stops
wherever `49` cells falls.

A cursor's default leaf is one cell, so the dual search descends every level and
spends most of its visits on the bottom one — for an IGeo7 level-12 column
rooted at level 5, 823543 of the tree's 960799 nodes. Stopping two levels early
deletes that layer; the leaf then hands back its `49` cells through
`child_indices_extents`, which loads their cached caps.

Measured on the production CopDEM GLO-90 -> IGeo7 L12 column regrid, core-seconds
per column, single-threaded:

| column         | leaf 1 | leaf 8 | **leaf 49** | leaf 343 |
|:---------------|-------:|-------:|------------:|---------:|
| 728            |  20.39 |  13.16 |   **12.69** |    18.67 |
| 98241 (all NaN)|  18.27 |  11.66 |   **10.97** |    16.60 |
| 115426 (polar) |  20.36 |  13.30 |   **12.34** |    18.89 |

!!! warning "This is not the analytical leaf size"
    A bigger leaf revisits each cell cap more often. Systems with a cheap
    analytical cap still evaluate it per visit and share one measured optimum,
    [`_ANALYTICAL_BUCKET_SIZE`](@ref), calibrated on IGeo7. A generic `cell_boundary` cap remains at
    the cursor's own leaf size.
"""
const _CACHED_BUCKET_SIZE = 49

"""
    _ANALYTICAL_BUCKET_SIZE

Leaf size for analytical caps. The cursor constructs one inline leaf when the
traversal first reaches it. The traversal carries that value while the opposing
tree descends, so this can use the same broad-search optimum as a precomputed
cap load without retaining a cache.
"""
const _ANALYTICAL_BUCKET_SIZE = Engine.ANALYTICAL_LEAF_CAPACITY

# Give a cursor the optimized seam's leaf size, leaving an explicit caller
# choice alone. `0` is the grid default ("descend to single cells"), not a
# request.
function _bucketed(c::HierarchicalGridCursor, bucket_size::Int = _CACHED_BUCKET_SIZE)
    c.bucket_size == 0 || return c
    return typeof(c)(c.grid, c.system, c.top_level, c.leaf_level,
        bucket_size, c.level, c.id, c.first_index, c.last_index,
        c.complete_subtree, c.selection)
end

# The whole-grid sweep tree, `treeify(grid; leaf_capacity = ForSweeps())`, is
# the engine's: the hierarchy cursor with `leaf_capacity(grid, ForSweeps())`
# cells per leaf, O(1) to build, its caps derived per visit. No cap vector is
# kept for a whole grid: that is an O(ncells) fill before the first query, and
# `ForSweeps` must stay the tree of its integer capacity. Caps are cached per
# chunk only (`_cachedchunktree`), where the chunk bounds the vector.

# Above this a chunk's cap vector costs more to fill than the revisits it saves
# (2 MiB of caps at the limit).
const _CHUNK_CAP_CACHE_MAX = 2^16

# A chunk's cursor with analytical uncached caps where available, otherwise its own
# cap vector when small enough. Analytical caps take the same bucketed leaf
# shape without retaining a vector proportional to the chunk.
function _cachedchunktree(cursor::HierarchicalGridCursor,
        inds::AbstractUnitRange{<:Integer})
    Fallbacks.cell_cap_is_cheap(cursor.grid) isa Val{true} &&
        return _bucketed(cursor, _ANALYTICAL_BUCKET_SIZE)
    # Past the limit the cursor goes back bare, so it keeps its own leaf size.
    length(inds) > _CHUNK_CAP_CACHE_MAX && return cursor
    caps = _ShiftedCaps(_leafcaps(cursor.grid, inds), Int(first(inds)) - 1)
    return CapCachedTree(_bucketed(cursor), caps)
end

# One tight cap per cell of `inds`; entry `k` is index `first(inds) + k - 1`.
function _leafcaps(grid::AbstractGrid, inds::AbstractUnitRange{<:Integer})
    n = length(inds)
    n == 0 && return [_cellcap(grid, Int(i)) for i in inds]
    off = Int(first(inds)) - 1
    caps = Vector{typeof(_cellcap(grid, off + 1))}(undef, n)
    nt = min(Threads.nthreads(), max(1, n >> 14))
    if nt > 1
        tasks = Vector{Task}(undef, nt)
        for t in 1:nt
            lo = (n * (t - 1)) ÷ nt + 1
            hi = (n * t) ÷ nt
            tasks[t] = Threads.@spawn _fillcaps!(caps, grid, off, lo, hi)
        end
        foreach(wait, tasks)
    else
        _fillcaps!(caps, grid, off, 1, n)
    end
    return caps
end

_cellcap(grid::AbstractGrid, i::Int) = Fallbacks.cell_cap(grid, cellindex(grid, i))

function _fillcaps!(caps::Vector, grid::AbstractGrid, off::Int, lo::Int, hi::Int)
    for k in lo:hi
        @inbounds caps[k] = _cellcap(grid, off + k)
    end
    return nothing
end

# --- GridInterface's `subtree` on a rooted chunk of a level grid ---------------
# `chunks(grid)` hands out rooted subsets whose positions are one ancestor's
# descendant range. Their tree is this package's hierarchy cursor started at
# that ancestor (O(1)), with the leaf shape and cap treatment the regridding
# face measured (`_cachedchunktree`); the leaves emit the level grid's positions.
# Any other range gets GridInterface's default.
function _rangecursor(grid::HierarchicalLevelGrid, inds::AbstractUnitRange{<:Integer})
    sys = system(grid)
    (isempty(inds) || !has_sorted_subtrees(sys)) && return nothing
    root = treeify(grid)
    (root isa HierarchicalGridCursor && root.selection === nothing) || return nothing
    c = cellindex(grid, Int(first(inds)))
    for l in levels(sys)
        l > root.leaf_level && break
        a = ancestor(sys, c, l)
        r = descendant_range(sys, a, root.leaf_level)
        (Int(first(r)) == Int(first(inds)) && Int(last(r)) == Int(last(inds))) || continue
        return typeof(root)(grid, sys, root.top_level, root.leaf_level, root.bucket_size,
            l, a, Int(first(inds)), Int(last(inds)), true, nothing)
    end
    return nothing
end

function GridInterface.subtree(grid::HierarchicalLevelGrid, inds::AbstractUnitRange{<:Integer})
    window = subcursor(grid, inds)      # dispatch-only; CopernicusDEM answers it
    window === nothing || return window
    cursor = _rangecursor(grid, inds)
    cursor === nothing || return _cachedchunktree(cursor, inds)
    return invoke(GridInterface.subtree, Tuple{Any,AbstractVector{<:Integer}}, grid, inds)
end
function GridInterface.subtree(grid::HierarchicalLevelGrid, s::GridInterface.Subset)
    pos = GridInterface.parentpositions(s)
    (parent(s) === grid && pos isa AbstractUnitRange) && return GridInterface.subtree(grid, pos)
    return invoke(GridInterface.subtree, Tuple{Any,GridInterface.Subset}, grid, s)
end
# A holding windows its own tree where it can (a Copernicus holding's tile
# chunks keep the tiled raster cursor); GridInterface's chunk route reaches
# this through the chunk's position range.
function GridInterface.subtree(grid::PartialGrid, inds::AbstractUnitRange{<:Integer})
    window = subcursor(grid, inds)
    window === nothing || return window
    return invoke(GridInterface.subtree, Tuple{Any,AbstractVector{<:Integer}}, grid, inds)
end
# A rooted chunk of a holding is the holding's hierarchy cursor started at the
# chunk's root, O(1), with the same cap treatment as a level grid's chunk.
# (GridInterface roots the cursor only for a grid without a parent.)
function GridInterface.subtree(grid::PartialGrid, s::GridInterface.Subset)
    pos = GridInterface.parentpositions(s)
    if parent(s) === grid && s.root !== nothing && pos isa AbstractUnitRange
        window = subcursor(grid, pos)
        window === nothing || return window
        cursor = _rootedcursor(grid, s.root, pos)
        cursor === nothing || return _cachedchunktree(cursor, pos)
    end
    return invoke(GridInterface.subtree, Tuple{Any,GridInterface.Subset}, grid, s)
end

function _rootedcursor(grid::AbstractGrid, root_id, inds::AbstractUnitRange{<:Integer})
    isempty(inds) && return nothing
    root = treeify(grid)
    (root isa HierarchicalGridCursor && root.selection === nothing) || return nothing
    complete_subtree = length(inds) ==
        length(descendant_range(root.system, root_id, root.leaf_level))
    return typeof(root)(grid, root.system, root.top_level, root.leaf_level,
        root.bucket_size, level(root_id), root_id, Int(first(inds)), Int(last(inds)),
        complete_subtree, nothing)
end
