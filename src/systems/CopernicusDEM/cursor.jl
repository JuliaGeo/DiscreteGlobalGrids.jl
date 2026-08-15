# An interior tree for a lattice that has none.
#
# `levels(sys) == 0:1` gives the generic `HierarchicalGridCursor` nothing to
# descend: a level-0 node's children are its 960 000 (GLO-90) or 12 960 000
# (GLO-30) pixels, and the whole-globe root's children are 64 800 tiles. A dual
# tree search against such a tree costs O(n_src x dst-depth) instead of
# O(log n_src x log n_dst), because one side never narrows.
#
# Nothing about the lattice forces that. A Copernicus DEM cell family is a
# RECTANGLE in two integer axes at both scales — 180 x 360 tiles, then N x ncols
# pixels inside a tile — and any axis-aligned sub-rectangle has a lon/lat box and
# a spherical cap in closed form from the band tables. So this file gives the
# system a real tree by RECURSIVELY SPLITTING those rectangles, down to single
# cells, with nothing materialised. A node is eight integers, and no node is
# ever built that the search does not visit.
#
# Building the whole 960 000-pixel `N50_00_E006_00` tile's matched-resolution
# regridder goes from 105.7 s to 11.5 s onto IGEO7 and from 83.9 s to 9.2 s onto
# HEALPix, with the intersection matrix identical to the last bit. How a node
# splits changes only the speed, and [`BlockStrategy`](@ref) has that table and
# the script that reproduces it.
#
# # The generic seam
#
# Only four functions know that the rectangle is a Copernicus DEM one:
#
#   * [`_node_box`](@ref)     — the node's lon/lat box, closed form;
#   * [`_leaf_pad`](@ref)     — how far a leaf's published ring may bow outside it;
#   * [`_position`](@ref)     — (tile, raster row, column) -> grid position;
#   * [`_childspace`](@ref)   — which rectangle a node's children partition.
#
# Everything else — the split arithmetic, the cap, and the thirteen
# `SpatialTreeInterface`/`Trees` methods — is about rectangles and would move to
# `ConservativeRegridding.Trees` unchanged as a cursor over an "abstract
# curvilinear grid" whose extent function is closed form. That is where it
# belongs long-term; it lives here because CR is pinned by `[sources]` in this
# repo's `Project.toml` and must not be edited from here.

# Cells at or below which a node stops splitting and yields its cells directly.
# Nine — one 3x3 block's worth: small enough that a leaf's cells are genuinely
# adjacent, large enough that the traversal is not all node overhead. It is the
# same trade `POSITION_TREE_LEAF_SIZE` makes for the fallback tree, one notch
# tighter because these cells are pixels.
const LEAF_CELLS = 9

"""
    BlockStrategy

How a [`BlockCursor`](@ref) node partitions its rectangle: [`Bisected`](@ref),
which `treeify` uses, or [`Blocked`](@ref), which it was measured against.

# Which one, and why it was measured rather than argued

Every strategy here produces the SAME intersection matrix, bit for bit — they
differ only in how fast the dual tree search reaches it. Building the
matched-resolution regridder for the whole `N50_00_E006_00` GLO-90 tile
(960 000 pixels), by `scripts/bench_copdem_cursor.jl`, which is the reproducer
for every number in this table:

| source tree                        | onto IGEO7 12 | onto HEALPix 16 |
|:-----------------------------------|--------------:|----------------:|
| generic `HierarchicalGridCursor`   |     105.7 s   |        83.9 s   |
| `Blocked{3}` — 9 children a node   |      52.3 s   |        51.9 s   |
| `Blocked{2}` — 4 children a node   |      24.1 s   |        25.1 s   |
| `Bisected`   — 2 children a node   |      11.5 s   |         9.2 s   |

MACHINE-LOCAL, and measured once: an M-series laptop on 8 threads, in the `docs`
project, which is the environment `examples/copernicus_dem.jl` runs in and quotes
its own two builds from. The ratios are the claim; the seconds are this machine's.

Monotone in the fanout, and the mechanism is the dual search's shape rather than
this tree's: when neither node is a leaf it tests `nchild(src) x nchild(dst)`
extent pairs before it can prune any of them
(`GeometryOps/src/utils/SpatialTreeInterface/dual_depth_first_search.jl:71-86`).
Nine source children against IGEO7's seven is 63 tests where two is 14, and a
coarse split also over-commits: a block nine times too big for the opposing node
cannot be pruned against it at all, so the two trees stop narrowing together.
Halving descends in the smallest step that still separates, so every level
prunes. The 20-level depth that costs is free — a node is eight integers with no
allocation.

The multithreaded walk that dominates an 8-thread build makes that argument
stronger than the single-threaded one does. `ConservativeRegridding`'s
`multithreaded_dual_depth_first_search`
(`src/utils/MultithreadedDualDepthFirstSearch.jl:27-31`) derives the opposing
node's extent INSIDE the loop over this node's children, so it calls
`node_extent` `nchild(src) x nchild(dst)` times per pair of interior nodes,
unconditionally: it has no expense cache, and `node_extent_is_expensive` — which
this cursor sets, because its cap is four inverse projections and four distances
derived on demand — buys nothing there. Every child this strategy does not
create is a row of that product not paid for.
"""
abstract type BlockStrategy end

"""
    Blocked{K}()

Split both axes into `K` near-equal parts, giving up to `K²` children. `K = 3`
reaches a single pixel of a 1 200 x 800 tile in seven levels.

The parts are `ceil`-divided, so an edge block can be one cell narrower than its
siblings. Nothing here needs them uniform: a node carries a box and a cap, not a
shape.

Not the default — see [`BlockStrategy`](@ref) for the measurement that decided
against it. Kept because it is the natural quadtree-style shape, because it is
the baseline the default is quoted against, and because a different destination
system with a much wider fanout could turn the table back.
"""
struct Blocked{K} <: BlockStrategy end

"""
    Bisected()

Split the LONGER axis in two, so the axis alternates by construction and blocks
stay near-square — in INDEX space, which is the only space this arithmetic sees:
a block of equal row and column counts spans `N/ncols` times as many degrees of
longitude as of latitude, ten times as many in the 10x reduction band above
latitude 85, of which the cosine of the latitude then takes most back on the
ground. What [`treeify`](@ref) builds; [`BlockStrategy`](@ref) carries the
measurement that chose it.
"""
struct Bisected <: BlockStrategy end

"The strategy `treeify` builds. See [`BlockStrategy`](@ref) for why this one."
const DEFAULT_STRATEGY = Bisected()

"""
    BlockCursor(grid; strategy = DEFAULT_STRATEGY)

A `GeometryOps.SpatialTreeInterface` cursor over a Copernicus DEM grid, built by
recursive splitting of the lattice rectangle. Prefer [`treeify`](@ref) to direct
construction; it falls back to the generic cursor for grids this cannot
represent. [`BlockStrategy`](@ref) is how a node splits, and carries the
measurement that picked the default.

A node is a rectangle at one of two scales, and `inpixels` says which:

  - `false`: tile rows `r0:r1` by tile columns `q0:q1`. At a level-0 grid its
    cells are those tiles; at a level-1 grid it is an interior node, and a node
    covering ONE tile descends straight into that tile's pixel rectangle rather
    than spending a level on itself.
  - `true`: raster rows `j0:j1` by raster columns `i0:i1` of the single tile
    `(r0, q0)`.

`origin` turns a lattice id into a grid position: `position = id - origin`. It is
`-1` for a complete level grid and `first_id - 1` for a partial one, which is why
this cursor is only built for grids whose ids are one contiguous run over one
rectangle — see [`treeify`](@ref).
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

The node's lon/lat box in degrees, containing every cell box beneath it.

A pixel rectangle is exactly [`cell_box`](@ref)'s expression over a column and a
row interval, so a 1x1 node reproduces `cell_box` bit for bit. A tile rectangle
is the union of its tiles' boxes, and those are NOT flush across a band
boundary: each tile's west edge is offset by half of ITS OWN `Δlon`, which steps
at latitude 50/60/70/80/85. The union takes the westernmost of them and drops
the east edge's offset entirely — the result over-covers by at most half a pixel
and a node box only has to COVER.
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
    # The westernmost half-pixel offset over the block's bands; the east edge
    # keeps the whole degree, which is the same over-covering by omission.
    widest = 0.0
    for r in c.r0:c.r1
        widest = max(widest, (1 / ncols(c.sys, r)) / 2)
    end
    return (_lon_w(c.q0) - widest, Float64(_lon_w(c.q1) + 1), south, north)
end

"""
    _leaf_pad(cursor) -> Float64

Radians of cap headroom at a LEAF, sized by the bow of a leaf's published ring
outside its own box. [`cell_boundary`](@ref) emits undensified geodesic quads,
so a cell's north and south edges bow poleward by about
`(Δλ²/8)·sin φ·cos φ ≤ Δλ²/16`, and this returns that bound for the leaf's own
`Δλ` — at level 1 a pixel of the block's COARSEST band, the largest `Δλ` any
leaf beneath the node can have, and at level 0 a whole tile.

# What the pad is actually doing, which is not bounding that bow

The cap's radius is set by the box's CORNERS ([`_box_cap`](@ref)), and a corner
is farther from the box centre than any point on a bowed edge — over leaves
sampled in every band and both pole rows, the corner clears the farthest
interior ring point by a measured median of 1.7 of these pads on a 1° tile and
about 1.8e3 of them on a pixel. So with `pad = 0` a leaf's ring is still inside
its own cap: no point in the INTERIOR of a ring edge — which is where the bow
lives — leaves the unpadded cap at all (worst `-2.4e-15` rad), and the only
overshoot is at the ring VERTICES, which are the corners the radius was built
from, at `+3.2e-17` rad. That is float rounding in `spherical_distance`, not
geometry, and it is asserted rather than merely measured in
`test/systems/CopernicusDEM/runtests.jl`, "the block cursor is a tree over the
lattice".

So the pad is belt-and-braces against rounding at the corners, and the `Δλ²/16`
expression is an ESTIMATE rather than a bound, in exactly the sense
`cell_boundary`'s own bow figure is (`system.jl:302-305`) — small angle, and it
errs low. Nothing here relies on it being tight; the covering law it serves is
asserted by sampling in that same testset, not derived from this number.
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

The grid position of a lattice cell named by tile `(r, q)` and — at level 1 —
raster `(j, i)`. Closed form through the same `tilebase` prefix sum the id codec
uses, offset by the grid's own first id.
"""
@inline function _position(c::BlockCursor, r::Int, q::Int, j::Int, i::Int)
    c.level == 0 && return Int(Int64(tileordinal(r, q)) - c.origin)
    return Int(tilebase(c.sys, r, q) + Int64(j) * ncols(c.sys, r) + Int64(i) - c.origin)
end

"""
    _childspace(cursor) -> (row_lo, nrows, col_lo, ncols, inpixels)

The rectangle a node's children partition, and at which scale.

A node covering exactly one tile of a LEVEL-1 grid returns that tile's whole
raster rather than itself: the tile scale and the pixel scale meet at one node,
so the tree spends no level on a one-tile block.
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

# Part `t` of `parts` near-equal pieces of `lo:lo+len-1`. The first `len % parts`
# pieces are one longer, so the pieces differ by at most one cell and none is
# empty.
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

A cap containing the lon/lat box, plus `pad` radians.

Two regimes, because the corner argument has a hypothesis. Write `h` for the
box's HALF-width in longitude. When `h ≤ 90°` the farthest point of the box from
its midpoint is a CORNER — fix `φ` and
`cos d = sin φ_c sin φ + cos φ_c cos φ cos(λ − λ_c)` falls as `|λ − λ_c|` grows,
so the extreme is at `λ ∈ {W, E}`; fix that `λ` and the expression is
`A sin φ + B cos φ` with `B = cos φ_c cos h ≥ 0` exactly when `h ≤ 90°`, a
sinusoid whose trough lies outside `[−90°, 90°]`, hence unimodal in `φ` with its
extreme at `φ ∈ {S, N}`. That is the same argument [`node_extent`](@ref) makes
for one cell, and it is why this cursor never samples.

A box wider than 180° breaks that hypothesis, and it is covered by a POLAR cap
instead: every point of a box reaching latitude `S` is within `90° − S` of the
north pole. Loose, and it costs nothing — bisection halves the span at every
level, so only the whole-globe root of a tile tree ever takes this branch.
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

# The cap is derived on demand — four inverse projections and four distances —
# rather than stored, so a node stays eight integers and the dual search should
# cache a node's child extents instead of re-deriving them per opposing child.
STI.node_extent_is_expensive(::Type{<:BlockCursor}) = true

function STI.isleaf(c::BlockCursor)
    # A tile block of a LEVEL-1 grid is never a leaf, however few tiles it holds:
    # its cells are pixels, and there are `ncols * N` of them per tile.
    c.level == 1 && !c.inpixels && return false
    _node_cells(c) <= LEAF_CELLS && return true
    # A node that would be its own only child is a leaf, whatever its cell count.
    # Unreachable for the two shipped strategies — a rectangle above `LEAF_CELLS`
    # always has an axis longer than one — and cheap insurance against a future
    # strategy that stalls, which would otherwise recurse forever.
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

A leaf's cells as (grid position, cap) pairs, each cap the 1x1 case of
[`_box_cap`](@ref) and so the same cap [`node_extent`](@ref) gives that cell.
Materialized because the dual search binds it once and iterates it against every
cell of the opposing leaf.
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

# Leaf indices are GRID positions at every node — the same index space
# `child_indices_extents` yields — so these two answer about the whole grid
# regardless of which node holds them, exactly as the generic position tree does.
Trees.ncells(c::BlockCursor) = DGG.ncells(c.grid)
Trees.getcell(c::BlockCursor, i::Int) = DGG.getcell(c.grid, i)
Trees.getcell(c::BlockCursor) = DGG.getcell(c.grid)

function Trees.should_parallelize(c::BlockCursor, ::US.SphericalCap)
    threshold = max(Int64(1), Int64(DGG.ncells(c.grid)) ÷
                              (Int64(Threads.nthreads()) * DGG.Fallbacks.PARALLELIZE_CHUNKS_PER_THREAD))
    return _node_cells(c) <= threshold
end

# ===========================================================================
# treeify
# ===========================================================================

"""
    treeify(grid::HierarchicalLevelGrid{<:CopernicusDEMSystem})
    treeify(grid::PartialGrid{<:CopernicusDEMSystem})

Give a Copernicus DEM grid the [`BlockCursor`](@ref) tree instead of the generic
[`HierarchicalGridCursor`](@ref), which on a two-level system is a root with
64 800 children over tile nodes with up to 12 960 000 each.

A `PartialGrid` gets it only when its cells are exactly one axis-aligned
rectangle of the lattice, held as one contiguous id run, which is what makes
`position = id - origin` a closed form. Three shapes qualify, and the id order
is why there are only three:

  * level 0 — one segment of a tile row, or whole tile rows;
  * level 1 — one tile's whole raster rows, which is `PartialGrid(sys, tile, 1)`
    and any band of it;
  * level 1 — a run of WHOLE tiles that is a tile rectangle by the level-0 rule,
    which descends as a tile node until it reaches one tile.

Anything else (a scattered id list, a window that starts mid-row, a multi-tile
window with a partial tile at either end) falls back to the generic cursor,
which is correct for every grid and merely slower here.
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
    # CONTIGUITY, and two endpoints are the whole of it: `PartialGrid` rejects an
    # id vector that is not STRICTLY ascending at construction
    # (`src/fallbacks/partial_grid.jl:80-81`), so `n` ids spanning `hi - lo + 1 == n`
    # values are exactly the run `lo:hi`, with none missing and none repeated. That
    # invariant is what makes the shape test below a test on `lo` and `hi` alone,
    # and what makes `position = id - origin` a closed form.
    hi.index - lo.index + 1 == n || return nothing
    # `decode` throws on an id no cell has and a `PartialGrid` does not range-check
    # its ids, while `treeify` owes an answer for every grid: a miss falls back to
    # the generic cursor here rather than escaping as an `ArgumentError`.
    # `cellposition` is the decoder that answers `nothing` instead of throwing.
    (DGG.cellposition(sys, lo) === nothing || DGG.cellposition(sys, hi) === nothing) &&
        return nothing
    origin = Int64(lo.index) - 1
    ra, qa, ja, ia = decode(sys, lo)
    rb, qb, jb, ib = decode(sys, hi)
    # One tile row segment, or whole tile rows: those are the contiguous tile runs
    # that are rectangles of the tile lattice, and there are no others.
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
    # Several tiles, each of them whole: the run starts at the first pixel of its
    # first tile and ends at the last pixel of its last, so it is every pixel of
    # that tile rectangle. The node is then a level-1 TILE node, which is the
    # interior shape `_childspace` and `_node_box` already build for a level-1
    # level grid — one tile down, it descends straight into the raster.
    N = Int(lat_intervals(sys))
    (ja == 0 && ia == 0 && jb == N - 1 && ib == Int(ncols(sys, rb)) - 1) ||
        return nothing
    tilerect || return nothing
    return BlockCursor(grid, sys, strategy, 1, origin, ra, rb, qa, qb,
        0, 0, 0, 0, false)
end
