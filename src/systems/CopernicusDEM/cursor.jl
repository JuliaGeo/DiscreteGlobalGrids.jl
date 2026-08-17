# A lazy spatial tree over tile and pixel index rectangles.

# Cells at or below which a node stops splitting and yields its cells directly.
const LEAF_CELLS = 9

"""
    BlockStrategy

How a [`BlockCursor`](@ref) partitions rectangles: [`Bisected`](@ref) by default,
or [`Blocked`](@ref). Both produce the same intersections; see the benchmark script.
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

Split the longer axis in two. This is what [`treeify`](@ref) builds.
"""
struct Bisected <: BlockStrategy end

"The strategy `treeify` builds."
const DEFAULT_STRATEGY = Bisected()

"""
    BlockCursor(grid; strategy = DEFAULT_STRATEGY)

A spatial-tree cursor over a Copernicus DEM rectangle. Prefer [`treeify`](@ref),
which falls back for grids this cursor cannot represent.

A node is either:

  - a tile rectangle (`inpixels == false`), or
  - a raster rectangle within tile `(r0, q0)` (`inpixels == true`).

`position = id - origin`; therefore partial grids must be one contiguous id run
forming one rectangle.
"""
struct BlockCursor{G<:DGG.AbstractGrid,S<:CopernicusDEMSystem,F<:BlockStrategy}
    grid::G
    sys::S
    strategy::F
    level::Int          # the GRID's level: 0 = tiles are cells, 1 = pixels are
    origin::Int64       # grid position of lattice id `x` is `x - origin`
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
function _node_box(c::BlockCursor{G,S}) where {G,S}
    N = lat_intervals(c.sys)
    half_dlat = (1 / N) / 2
    if c.inpixels
        lat_s = _lat_s(c.r0)
        lon_w = _lon_w(c.q0)
        nc = ncols(c.sys, c.r0)
        half_dlon = (1 / nc) / 2
        west = (lon_w + c.i0 / nc) - half_dlon
        east = (lon_w + (c.i1 + 1) / nc) - half_dlon
        north = (lat_s + 1 - c.j0 / N) + half_dlat
        south = (lat_s + 1 - (c.j1 + 1) / N) + half_dlat
        c.j0 == 0 && lat_s == 89 && (north = 90.0)
        c.j1 == N - 1 && lat_s == -90 && (south = -90.0)
        return (west, east, south, north)
    end
    lat_n = _lat_s(c.r0) + 1
    lat_s = _lat_s(c.r1)
    north = _lat_s(c.r0) == 89 ? 90.0 : lat_n + half_dlat
    south = lat_s == -90 ? -90.0 : lat_s + half_dlat
    # Westernmost half-pixel offset over the block's bands.
    widest = 0.0
    for r in c.r0:c.r1
        widest = max(widest, (1 / ncols(c.sys, r)) / 2)
    end
    return (_lon_w(c.q0) - widest, Float64(_lon_w(c.q1) + 1), south, north)
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
    _position(cursor, r, q, j, i) -> Int

The grid position of the lattice cell at tile `(r, q)` and — at level 1 —
raster `(j, i)`. Closed form via `tilebase`, offset by the grid's first id.
"""
@inline function _position(c::BlockCursor, r::Int, q::Int, j::Int, i::Int)
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
@inline function _part(lo::Int, len::Int, parts::Int, t::Int)
    base, rem = divrem(len, parts)
    off = (t - 1) * base + min(t - 1, rem)
    return (lo + off, lo + off + base + (t <= rem ? 1 : 0) - 1)
end

# How many parts each axis of a `(nrows, ncols)` rectangle splits into.
@inline _parts(::Blocked{K}, nr::Int, nc::Int) where {K} = (min(K, nr), min(K, nc))
@inline _parts(::Bisected, nr::Int, nc::Int) =
    nr >= nc ? (min(2, nr), 1) : (1, min(2, nc))

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

"""
    STI.child_indices_extents(cursor) -> Vector{Tuple{Int,SphericalCap{Float64}}}

A leaf's cells as `(grid position, cap)` pairs.
"""
function STI.child_indices_extents(c::BlockCursor)
    STI.isleaf(c) ||
        throw(ArgumentError("child_indices_extents is only valid for leaf nodes"))
    entries = Tuple{Int,Cap}[]
    if c.inpixels
        for j in c.j0:c.j1, i in c.i0:c.i1
            leaf = BlockCursor(c.grid, c.sys, c.strategy, c.level, c.origin,
                c.r0, c.r0, c.q0, c.q0, j, j, i, i, true)
            push!(entries, (_position(c, c.r0, c.q0, j, i), STI.node_extent(leaf)))
        end
    else
        for r in c.r0:c.r1, q in c.q0:c.q1
            leaf = BlockCursor(c.grid, c.sys, c.strategy, c.level, c.origin,
                r, r, q, q, 0, 0, 0, 0, false)
            push!(entries, (_position(c, r, q, 0, 0), STI.node_extent(leaf)))
        end
    end
    return entries
end

# ===========================================================================
# ConservativeRegridding.Trees
# ===========================================================================

GOCore.best_manifold(c::BlockCursor) = GOCore.best_manifold(c.grid)

# Leaf indices are grid positions, so these answer about the whole grid
# regardless of which node holds them.
Trees.ncells(c::BlockCursor) = DGG.ncells(c.grid)
Trees.getcell(c::BlockCursor, i::Int) = DGG.getcell(c.grid, i)
Trees.getcell(c::BlockCursor) = DGG.getcell(c.grid)

function Trees.should_parallelize(c::BlockCursor, ::US.SphericalCap)
    threshold = max(Int64(1), Int64(DGG.ncells(c.grid)) ÷
                              (Int64(Threads.nthreads()) * DGG.PARALLELIZE_CHUNKS_PER_THREAD))
    return _node_cells(c) <= threshold
end

# ===========================================================================
# treeify
# ===========================================================================

"""
    treeify(grid::HierarchicalLevelGrid{<:CopernicusDEMSystem})
    treeify(grid::PartialGrid{<:CopernicusDEMSystem})

Use a [`BlockCursor`](@ref) when the grid has a rectangular contiguous id run.

A `PartialGrid` qualifies only for:

  * level 0 — one segment of a tile row, or whole tile rows;
  * level 1 — one tile's whole raster rows, or a run of whole tiles that is a
    tile rectangle by the level-0 rule.

Anything else falls back to [`HierarchicalGridCursor`](@ref).
"""
DGG.treeify(::GOCore.Manifold, grid::LevelGrid) = BlockCursor(grid)
DGG.treeify(::GOCore.Manifold, c::BlockCursor) = c
DGG.treeify(c::BlockCursor) = c

function DGG.treeify(::GOCore.Manifold,
        grid::DGG.PartialGrid{<:CopernicusDEMSystem})
    cursor = _block_cursor(grid, DEFAULT_STRATEGY)
    return cursor === nothing ? DGG.HierarchicalGridCursor(grid) : cursor
end

BlockCursor(grid::LevelGrid; strategy::BlockStrategy=DEFAULT_STRATEGY) =
    _level_cursor(grid, strategy)

function BlockCursor(grid::DGG.PartialGrid{<:CopernicusDEMSystem};
        strategy::BlockStrategy=DEFAULT_STRATEGY)
    cursor = _block_cursor(grid, strategy)
    cursor === nothing && throw(ArgumentError(
        "this partial grid's cells are not one rectangle of the Copernicus DEM " *
        "lattice held as one contiguous id run, so `position = id - origin` does " *
        "not hold; `treeify` falls back to HierarchicalGridCursor for it"))
    return cursor
end

function _level_cursor(grid::LevelGrid, strategy::BlockStrategy)
    sys = DGG.system(grid)
    l = DGG.level(grid)
    return BlockCursor(grid, sys, strategy, l, Int64(-1),
        0, NROWS - 1, 0, NCOLS_TILES - 1, 0, 0, 0, 0, false)
end

# `nothing` when the grid is not one rectangle held as one contiguous id run.
function _block_cursor(grid::DGG.PartialGrid{<:CopernicusDEMSystem},
        strategy::BlockStrategy)
    sys = DGG.system(grid)
    l = DGG.level(grid)
    n = DGG.ncells(grid)
    n == 0 && return nothing
    lo = DGG.cellindex(grid, 1)
    hi = DGG.cellindex(grid, n)
    # Strictly ascending ids form a run iff their span equals their count.
    hi.index - lo.index + 1 == n || return nothing
    # `PartialGrid` does not range-check ids.
    (DGG.cellposition(sys, lo) === nothing || DGG.cellposition(sys, hi) === nothing) &&
        return nothing
    origin = Int64(lo.index) - 1
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
