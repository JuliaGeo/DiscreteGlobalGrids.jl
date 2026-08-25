# A lazy spatial tree over tile and pixel index rectangles, for the complete
# lattice; a holding of tiles gets the engine's tiled raster tree instead.

# The leaf size and the leaf container are the engine's, so both trees over this
# lattice hand a leaf back the same way.
const LEAF_CELLS = Engine.LEAF_CELLS
const LeafCells = Engine.LeafCells

"""
    BlockStrategy

How a [`BlockCursor`](@ref) partitions rectangles: [`Bisected`](@ref) by default,
or [`Blocked`](@ref). Both produce the same intersections;
`scripts/bench_copdem_cursor.jl` times them.
"""
abstract type BlockStrategy end

"""
    Blocked{K}()

Split both axes into `K` near-equal parts, giving up to `K²` children. Edge blocks
may be one cell narrower. Bisection benchmarked faster.
"""
struct Blocked{K} <: BlockStrategy end

"""
    Bisected()

Split the longer axis in two, the rule both trees over this lattice use. This
is what [`treeify`](@ref) builds.
"""
struct Bisected <: BlockStrategy end

"The strategy `treeify` builds."
const DEFAULT_STRATEGY = Bisected()

"""
    BlockCursor(grid::HierarchicalLevelGrid{<:CopernicusDEMSystem};
                strategy = DEFAULT_STRATEGY)

A spatial-tree cursor over a rectangle of the complete Copernicus DEM lattice.
Prefer [`treeify`](@ref), which wraps it in a [`MemoBlockCursor`](@ref).

  - A node is either a tile rectangle (`inpixels == false`) or a raster
    rectangle within tile `(r0, q0)` (`inpixels == true`).
  - Its addressing law is `index = id - origin`, so it holds only where the
    cells it covers are one contiguous id run forming one rectangle — the whole
    level, and the windows [`subcursor`](@ref) cuts out of it.
  - A holding of tiles satisfies neither in general and gets the engine's
    [`TiledRasterCursor`](@ref), whose law is a tile offset plus a row-major
    pixel index.
"""
struct BlockCursor{G<:DGG.AbstractGrid,S<:CopernicusDEMSystem,F<:BlockStrategy}
    grid::G
    sys::S
    strategy::F
    level::Int          # the GRID's level: 0 = tiles are cells, 1 = pixels are
    origin::Int64       # grid index of lattice id `x` is `x - origin`
    r0::Int
    r1::Int
    q0::Int
    q1::Int             # tile rectangle, inclusive
    j0::Int
    j1::Int
    i0::Int
    i1::Int             # raster rectangle inside tile `(r0, q0)`, inclusive
    inpixels::Bool
end

function Base.show(io::IO, c::BlockCursor)
    print(io, "BlockCursor(", c.inpixels ? "pixels " : "tiles ",
        c.inpixels ? "$(c.j0):$(c.j1) x $(c.i0):$(c.i1) of tile ($(c.r0),$(c.q0))" :
        "$(c.r0):$(c.r1) x $(c.q0):$(c.q1)", ", level=", c.level, ")")
end

# ===========================================================================
# The generic seam: four functions that know about this lattice
# ===========================================================================

"""
    _node_box(cursor) -> (west, east, south, north)

The node's longitude/latitude box in degrees. It contains every descendant box;
tile rectangles use the westernmost band offset.
"""
_node_box(c::BlockCursor) = c.inpixels ?
    _pixel_box(c.sys, c.r0, c.q0, c.j0, c.j1, c.i0, c.i1) :
    _tile_box(c.sys, c.r0, c.r1, c.q0, c.q1)

"""
    _pixel_box(sys, r, q, j0, j1, i0, i1) -> (west, east, south, north)

The longitude/latitude box in degrees of raster rows `j0:j1` and columns
`i0:i1` of tile `(r, q)`. Pole rows extend to the pole.
"""
function _pixel_box(sys::CopernicusDEMSystem, r::Int, q::Int, j0::Int, j1::Int,
        i0::Int, i1::Int)
    N = Int(lat_intervals(sys))
    half_dlat = (1 / N) / 2
    lat_s = _lat_s(r)
    lon_w = _lon_w(q)
    nc = ncols(sys, r)
    half_dlon = (1 / nc) / 2
    west = (lon_w + i0 / nc) - half_dlon
    east = (lon_w + (i1 + 1) / nc) - half_dlon
    north = (lat_s + 1 - j0 / N) + half_dlat
    south = (lat_s + 1 - (j1 + 1) / N) + half_dlat
    j0 == 0 && lat_s == 89 && (north = 90.0)
    j1 == N - 1 && lat_s == -90 && (south = -90.0)
    return (west, east, south, north)
end

"""
    _tile_box(sys, r0, r1, q0, q1) -> (west, east, south, north)

The longitude/latitude box in degrees of the tile rectangle `r0:r1` x `q0:q1`,
taken over the westernmost half-pixel offset of the bands it spans.
"""
function _tile_box(sys::CopernicusDEMSystem, r0::Int, r1::Int, q0::Int, q1::Int)
    N = Int(lat_intervals(sys))
    half_dlat = (1 / N) / 2
    lat_n = _lat_s(r0) + 1
    lat_s = _lat_s(r1)
    north = _lat_s(r0) == 89 ? 90.0 : lat_n + half_dlat
    south = lat_s == -90 ? -90.0 : lat_s + half_dlat
    widest = 0.0
    for r in r0:r1
        widest = max(widest, (1 / ncols(sys, r)) / 2)
    end
    return (_lon_w(q0) - widest, Float64(_lon_w(q1) + 1), south, north)
end

"""
    _leaf_pad(cursor) -> Float64

Radians of cap headroom at a leaf: `Δλ²/16`, the poleward bow of an undensified
geodesic quad edge, for the coarsest `Δλ` any leaf beneath the node can have.
"""
function _leaf_pad(c::BlockCursor)
    c.level == 0 && return deg2rad(1.0)^2 / 16
    if c.inpixels
        return deg2rad(1 / ncols(c.sys, c.r0))^2 / 16
    end
    coarsest = 0.0
    for r in c.r0:c.r1
        coarsest = max(coarsest, 1 / ncols(c.sys, r))
    end
    return deg2rad(coarsest)^2 / 16
end

"""
    _index(cursor, r, q, j, i) -> Int

The grid index of the lattice cell at tile `(r, q)` and — at level 1 —
raster `(j, i)`. Closed form via `tilebase`, offset by the grid's first id.
"""
@inline function _index(c::BlockCursor, r::Int, q::Int, j::Int, i::Int)
    c.level == 0 && return Int(Int64(tileordinal(r, q)) - c.origin)
    return Int(tilebase(c.sys, r, q) + Int64(j) * ncols(c.sys, r) + Int64(i) - c.origin)
end

"""
    _childspace(cursor) -> (row_lo, nrows, col_lo, ncols, inpixels)

The rectangle a node's children partition, and at which scale. A node covering
exactly one tile of a level-1 grid returns that tile's whole raster, so the
tree spends no level on a one-tile block.
"""
function _childspace(c::BlockCursor)
    c.inpixels && return (c.j0, c.j1 - c.j0 + 1, c.i0, c.i1 - c.i0 + 1, true)
    if c.level == 1 && c.r0 == c.r1 && c.q0 == c.q1
        return (0, Int(lat_intervals(c.sys)), 0, Int(ncols(c.sys, c.r0)), true)
    end
    return (c.r0, c.r1 - c.r0 + 1, c.q0, c.q1 - c.q0 + 1, false)
end

# ===========================================================================
# Rectangles: splitting, counting, and the cap
# ===========================================================================

# Part `t` of a non-empty, near-equal partition of `lo:lo+len-1`.
const _part = Engine.rect_part

# How many parts each axis of a `(nrows, ncols)` rectangle splits into.
@inline _parts(::Blocked{K}, nr::Int, nc::Int) where {K} = (min(K, nr), min(K, nc))
@inline _parts(::Bisected, nr::Int, nc::Int) = Engine.bisect_parts(nr, nc)

"Cells in this node — pixels at level 1, tiles at level 0."
function _node_cells(c::BlockCursor)
    c.inpixels && return Int64(c.j1 - c.j0 + 1) * Int64(c.i1 - c.i0 + 1)
    width = Int64(c.q1 - c.q0 + 1)
    c.level == 0 && return Int64(c.r1 - c.r0 + 1) * width
    N = Int64(lat_intervals(c.sys))
    total = Int64(0)
    for r in c.r0:c.r1
        total += width * ncols(c.sys, r) * N
    end
    return total
end

const NODE_NORTH_POLE = GO.UnitSphericalPoint(0.0, 0.0, 1.0)
const NODE_SOUTH_POLE = GO.UnitSphericalPoint(0.0, 0.0, -1.0)

"""
    _box_cap(west, east, south, north, pad) -> SphericalCap

A cap containing the lon/lat box plus `pad` radians. Boxes wider than 180° use
a polar cap.
"""
function _box_cap(west::Float64, east::Float64, south::Float64, north::Float64,
        pad::Float64)
    if east - west > 180.0
        rn = deg2rad(90.0 - south)
        rs = deg2rad(north + 90.0)
        return rn <= rs ? SphericalCap(NODE_NORTH_POLE, min(Float64(π), rn + pad)) :
               SphericalCap(NODE_SOUTH_POLE, min(Float64(π), rs + pad))
    end
    centre = TO_SPHERE(((west + east) / 2, (south + north) / 2))
    rmax = 0.0
    for (lon, lat) in ((west, south), (east, south), (east, north), (west, north))
        rmax = max(rmax, US.spherical_distance(centre, TO_SPHERE((lon, lat))))
    end
    return SphericalCap(centre, nextfloat(min(Float64(π), rmax + pad)))
end

# ===========================================================================
# SpatialTreeInterface
# ===========================================================================

STI.isspatialtree(::Type{<:BlockCursor}) = true

# Derive caps on demand to keep nodes compact.
STI.node_extent_is_expensive(::Type{<:BlockCursor}) = true

function STI.isleaf(c::BlockCursor)
    # A tile block of a level-1 grid is never a leaf: its cells are pixels.
    c.level == 1 && !c.inpixels && return false
    _node_cells(c) <= LEAF_CELLS && return true
    # A stalled split must terminate.
    _, nr, _, nc, _ = _childspace(c)
    pr, pc = _parts(c.strategy, nr, nc)
    return pr * pc <= 1
end

function STI.nchild(c::BlockCursor)
    STI.isleaf(c) && return 0
    _, nr, _, nc, _ = _childspace(c)
    pr, pc = _parts(c.strategy, nr, nc)
    return pr * pc
end

function STI.getchild(c::BlockCursor, k::Int)
    n = STI.nchild(c)
    1 <= k <= n || throw(BoundsError(c, k))
    rlo, nr, clo, nc, pixels = _childspace(c)
    pr, pc = _parts(c.strategy, nr, nc)
    tr, tc = (k - 1) ÷ pc + 1, (k - 1) % pc + 1
    a0, a1 = _part(rlo, nr, pr, tr)
    b0, b1 = _part(clo, nc, pc, tc)
    pixels && return BlockCursor(c.grid, c.sys, c.strategy, c.level, c.origin,
        c.r0, c.r0, c.q0, c.q0, a0, a1, b0, b1, true)
    return BlockCursor(c.grid, c.sys, c.strategy, c.level, c.origin,
        a0, a1, b0, b1, 0, 0, 0, 0, false)
end

function STI.node_extent(c::BlockCursor)
    west, east, south, north = _node_box(c)
    return _box_cap(west, east, south, north, _leaf_pad(c))
end

# (cells along the fast axis, cells in the leaf).
@inline function _leafshape(c::BlockCursor)
    c.inpixels && return (c.i1 - c.i0 + 1, (c.j1 - c.j0 + 1) * (c.i1 - c.i0 + 1))
    return (c.q1 - c.q0 + 1, (c.r1 - c.r0 + 1) * (c.q1 - c.q0 + 1))
end

# Entry `k` of the leaf's row-major cell order.
@inline function _leafcell(c::BlockCursor, nfast::Int, k::Int)
    slow, fast = divrem(k - 1, nfast)
    if c.inpixels
        j = c.j0 + slow
        i = c.i0 + fast
        leaf = BlockCursor(c.grid, c.sys, c.strategy, c.level, c.origin,
            c.r0, c.r0, c.q0, c.q0, j, j, i, i, true)
        return (_index(c, c.r0, c.q0, j, i), STI.node_extent(leaf))
    end
    r = c.r0 + slow
    q = c.q0 + fast
    leaf = BlockCursor(c.grid, c.sys, c.strategy, c.level, c.origin,
        r, r, q, q, 0, 0, 0, 0, false)
    return (_index(c, r, q, 0, 0), STI.node_extent(leaf))
end

"""
    LeafCells(cursor::BlockCursor)

A [`BlockCursor`](@ref) leaf's cells as `(grid index, cap)` pairs. See
[`LeafCells`](@ref DiscreteGlobalGrids.Engine.LeafCells) for the container.
"""

function LeafCells(c::BlockCursor)
    nfast, n = _leafshape(c)
    1 <= n <= LEAF_CELLS || throw(ArgumentError(
        "a leaf holds 1 to $LEAF_CELLS cells; $c reports $n"))
    return Engine.leaf_cells(k -> _leafcell(c, nfast, k), n)
end

"""
    STI.child_indices_extents(cursor) -> LeafCells

A leaf's cells as `(grid index, cap)` pairs. See [`LeafCells`](@ref).
"""
function STI.child_indices_extents(c::BlockCursor)
    STI.isleaf(c) ||
        throw(ArgumentError("child_indices_extents is only valid for leaf nodes"))
    return LeafCells(c)
end

# ===========================================================================
# ConservativeRegridding.Trees
# ===========================================================================

GOCore.best_manifold(c::BlockCursor) = GOCore.best_manifold(c.grid)

# Leaf indices are grid indices, so these answer about the whole grid
# regardless of which node holds them.
Trees.ncells(c::BlockCursor) = DGG.ncells(c.grid)
Trees.getcell(c::BlockCursor, i::Int) = DGG.getcell(c.grid, i)
Trees.getcell(c::BlockCursor) = DGG.getcell(c.grid)

# `Trees.ncells` answers for the whole grid, so the frontier's default estimate
# would be wrong here; a block knows its own cell count.
Trees.split_weight(c::BlockCursor) = Int(_node_cells(c))

# ===========================================================================
# The lattice as raster tiles
# ===========================================================================

"""
    RasterTile(r, q, j0, i0, nrows, ncols, offset, pixels)

One rectangle of a Copernicus DEM holding, as [`raster_tiles`](@ref) names it:
rows `j0:j0+nrows-1` and columns `i0:i0+ncols-1` of tile `(r, q)`, whose pixels
occupy grid indices `offset+1 : offset+nrows*ncols` in row-major order.

  - `pixels` is `false` on a level-0 grid, where the cell is the tile itself and
    the rectangle is the tile's own `1x1`.
  - The engine treats this as an opaque handle, so it carries everything the
    three geometry hooks need without decoding an id again.
"""
struct RasterTile
    r::Int
    q::Int
    j0::Int
    i0::Int
    nrows::Int
    ncols::Int
    offset::Int
    pixels::Bool
end

Base.show(io::IO, t::RasterTile) = print(io, "RasterTile(",
    t.pixels ? "rows $(t.j0):$(t.j0 + t.nrows - 1) x cols $(t.i0):$(t.i0 + t.ncols - 1) of " : "",
    "tile ($(t.r),$(t.q)), offset=", t.offset, ")")

"""
    raster_tiles(grid::PartialGrid{<:CopernicusDEMSystem}, inds)

The rectangles of the lattice that hold grid indices `inds`, or `nothing` when
the grid names an id this lattice does not have.

  - A whole tile is one rectangle, the shape a holding of tiles is made of.
  - A run that starts or ends inside a raster row splits into at most three: a
    part row, the whole rows below it, and a part row.
  - Every rectangle is therefore one contiguous block of grid indices, and the
    row-major index law holds on it.
"""
function DGG.raster_tiles(grid::DGG.PartialGrid{<:CopernicusDEMSystem},
        inds::AbstractUnitRange{<:Integer})
    isempty(inds) && return RasterTile[]
    (1 <= first(inds) && last(inds) <= DGG.ncells(grid)) || return nothing
    sys = DGG.system(grid)
    l = DGG.level(grid)
    N = Int(lat_intervals(sys))
    tiles = RasterTile[]
    p, stop = Int(first(inds)), Int(last(inds))
    while p <= stop
        c = DGG.cellindex(grid, p)
        # `PartialGrid` does not range-check its ids; an id off this lattice has
        # no rectangle, and the grid keeps the generic cursor.
        DGG.globalindex(sys, c) === nothing && return nothing
        r, q, j, i = decode(sys, c)
        if l == 0
            push!(tiles, RasterTile(r, q, 0, 0, 1, 1, p - 1, false))
            p += 1
            continue
        end
        nc = Int(ncols(sys, r))
        len = _idrun(grid, p, min(stop - p + 1, nc * N - (j * nc + i)))
        jj, rest = j, len
        if i > 0
            head = min(len, nc - i)
            push!(tiles, RasterTile(r, q, j, i, 1, head, p - 1, true))
            p += head
            jj, rest = j + 1, len - head
        end
        full = rest ÷ nc
        if full > 0
            push!(tiles, RasterTile(r, q, jj, 0, full, nc, p - 1, true))
            p += full * nc
            jj += full
            rest -= full * nc
        end
        if rest > 0
            push!(tiles, RasterTile(r, q, jj, 0, 1, rest, p - 1, true))
            p += rest
        end
    end
    return tiles
end

# The longest run of consecutive ids starting at grid index `p`, at most
# `maxlen`. Ids ascend strictly, so `id(p+k) - id(p) - k` never decreases and a
# binary search finds the last `k` where it is still zero.
function _idrun(grid, p::Int, maxlen::Int)
    base = DGG.cellindex(grid, p).index
    lo, hi = 0, maxlen - 1
    while lo < hi
        mid = (lo + hi + 1) >> 1
        if DGG.cellindex(grid, p + mid).index - base == mid
            lo = mid
        else
            hi = mid - 1
        end
    end
    return lo + 1
end

"The rectangle's extent in pixel rows and columns."
DGG.raster_shape(::DGG.PartialGrid{<:CopernicusDEMSystem}, t::RasterTile) =
    (t.nrows, t.ncols)

"""
    raster_localindex(grid, tile, j, i)

The rectangle's offset in the grid plus the pixel's row-major index within it.
"""
DGG.raster_localindex(::DGG.PartialGrid{<:CopernicusDEMSystem}, t::RasterTile,
    j::Int, i::Int) = t.offset + j * t.ncols + i + 1

"""
    raster_cap(grid, tile, j0, j1, i0, i1)

The cap of the sub-rectangle's [`_pixel_box`](@ref), padded for edge bow by
`Δλ²/16` — the same cap [`BlockCursor`](@ref) derives for the same pixels.
"""
function DGG.raster_cap(grid::DGG.PartialGrid{<:CopernicusDEMSystem}, t::RasterTile,
        j0::Int, j1::Int, i0::Int, i1::Int)
    sys = DGG.system(grid)
    if !t.pixels
        west, east, south, north = _tile_box(sys, t.r, t.r, t.q, t.q)
        return _box_cap(west, east, south, north, deg2rad(1.0)^2 / 16)
    end
    west, east, south, north = _pixel_box(sys, t.r, t.q, t.j0 + j0, t.j0 + j1,
        t.i0 + i0, t.i0 + i1)
    return _box_cap(west, east, south, north, deg2rad(1 / ncols(sys, t.r))^2 / 16)
end

# ===========================================================================
# treeify
# ===========================================================================

"""
    treeify(grid::HierarchicalLevelGrid{<:CopernicusDEMSystem})
    treeify(grid::PartialGrid{<:CopernicusDEMSystem})

  - The complete lattice is one rectangle, and gets a [`BlockCursor`](@ref)
    whose node boxes are closed-form in the lattice coordinates: no table of
    tile caps, and `O(1)` to build.
  - A holding is a collection of tiles in no particular arrangement, and gets
    the engine's [`TiledRasterCursor`](@ref): a packed tree over the tiles' caps
    with a bisection quadtree inside each tile. Any tile set qualifies, one
    crossing a latitude row included.
  - Only an id this lattice does not name keeps the generic
    [`HierarchicalGridCursor`](@ref), such an id having no rectangle to place.
  - The block cursor comes back wrapped in a [`MemoBlockCursor`](@ref), which
    memoizes derived node extents per task; `BlockCursor(grid)` gives the bare
    cursor. The tiled raster cursor memoizes its own.
"""
DGG.treeify(::GOCore.Manifold, grid::LevelGrid) = _memoized(BlockCursor(grid))
DGG.treeify(::GOCore.Manifold, c::BlockCursor) = c
DGG.treeify(c::BlockCursor) = c

function DGG.treeify(::GOCore.Manifold,
        grid::DGG.PartialGrid{<:CopernicusDEMSystem})
    tree = Engine.tiled_raster_tree(grid)
    return tree === nothing ? DGG.HierarchicalGridCursor(grid) : tree
end

BlockCursor(grid::LevelGrid; strategy::BlockStrategy=DEFAULT_STRATEGY) =
    _level_cursor(grid, strategy)

function _level_cursor(grid::LevelGrid, strategy::BlockStrategy)
    sys = DGG.system(grid)
    l = DGG.level(grid)
    return BlockCursor(grid, sys, strategy, l, Int64(-1),
        0, NROWS - 1, 0, NCOLS_TILES - 1, 0, 0, 0, 0, false)
end

"""
    subcursor(grid, inds) -> tree or `nothing`

The node covering grid indices `inds`, with leaf indices still `grid`'s own.

  - On the complete lattice: a [`MemoBlockCursor`](@ref) when `inds` is one
    contiguous id run forming one lattice rectangle — a segment of a tile row or
    whole tile rows at level 0, one tile's whole raster rows or a run of whole
    tiles at level 1 — and `nothing` otherwise. The run test is over the window,
    so a tile-sized chunk qualifies.
  - On a holding: the [`TiledRasterCursor`](@ref) over the rectangles that
    window holds, which every window has.
"""
DGG.subcursor(grid::LevelGrid, inds::AbstractUnitRange{<:Integer}) =
    _memoized(_window_cursor(grid, DEFAULT_STRATEGY, inds))

DGG.subcursor(grid::DGG.PartialGrid{<:CopernicusDEMSystem},
    inds::AbstractUnitRange{<:Integer}) = Engine.tiled_raster_tree(grid, inds)

# `nothing` unless `inds` is one contiguous id run forming one rectangle.
function _window_cursor(grid::LevelGrid, strategy::BlockStrategy,
        inds::AbstractUnitRange{<:Integer})
    isempty(inds) && return nothing
    lo_p, hi_p = Int(first(inds)), Int(last(inds))
    (1 <= lo_p && hi_p <= DGG.ncells(grid)) || return nothing
    lo = DGG.cellindex(grid, lo_p)
    hi = DGG.cellindex(grid, hi_p)
    # Strictly ascending ids form a run iff their span equals their count.
    hi.index - lo.index + 1 == hi_p - lo_p + 1 || return nothing
    return _run_cursor(grid, strategy, Int64(lo.index) - lo_p, lo, hi)
end

# The node covering the id run `lo:hi`, whose grid index is `id - origin`.
function _run_cursor(grid::LevelGrid, strategy::BlockStrategy, origin::Int64,
        lo::DGG.LevelIndex, hi::DGG.LevelIndex)
    sys = DGG.system(grid)
    l = DGG.level(grid)
    ra, qa, ja, ia = decode(sys, lo)
    rb, qb, jb, ib = decode(sys, hi)
    # Rectangular tile runs are one row segment or whole rows.
    tilerect = ra == rb || (qa == 0 && qb == NCOLS_TILES - 1)
    if l == 0
        tilerect || return nothing
        return BlockCursor(grid, sys, strategy, 0, origin, ra, rb, qa, qb,
            0, 0, 0, 0, false)
    end
    if ra == rb && qa == qb
        # One tile, whole raster rows.
        nc = Int(ncols(sys, ra))
        (ia == 0 && ib == nc - 1) || return nothing
        return BlockCursor(grid, sys, strategy, 1, origin, ra, ra, qa, qa,
            ja, jb, 0, nc - 1, true)
    end
    # A run of whole tiles, represented as a level-1 tile node.
    N = Int(lat_intervals(sys))
    (ja == 0 && ia == 0 && jb == N - 1 && ib == Int(ncols(sys, rb)) - 1) ||
        return nothing
    tilerect || return nothing
    return BlockCursor(grid, sys, strategy, 1, origin, ra, rb, qa, qb,
        0, 0, 0, 0, false)
end
