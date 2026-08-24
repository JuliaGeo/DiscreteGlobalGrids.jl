# ISEA4R uses 0-based Morton identifiers over ten diamonds. Parent/child and
# subtree operations are radix-4 arithmetic and grid positions are identifier
# + 1, so they are the quad-face family's; this file writes the ten-diamond
# chart's share.

# ===========================================================================
# Types
# ===========================================================================

"""
    ISEA4RSystem() <: AbstractQuadFaceGridSystem

Ten equal-area ISEA rhombus charts refined by aperture-4 subdivision. Level
`l in 0:29` contains `10*4^l` cells of solid angle `4π/(10*4^l)`.

Canonical [`LevelIndex`](@ref) values store
`diamond*4^level + morton(ix,iy)`. Morton order gives contiguous subtrees.
Vertex connectivity has at most 9 neighbors and edge connectivity at most 4.
Chart edges are curved and [`cell_boundary`](@ref) densifies them.
"""
struct ISEA4RSystem <: DGG.AbstractQuadFaceGridSystem end

# Grid descriptor for all `10 * 4^l` cells in Morton order.
const LevelGrid = DGG.HierarchicalLevelGrid{ISEA4RSystem}

"The deepest level at which `10 * 4^level` still fits a signed 64-bit integer."
const MAX_LEVEL = 29

"The number of diamonds — the ten rhombi of `diamonds.jl`, and the root count."
const NDIAMONDS = 10

# ===========================================================================
# System interface
# ===========================================================================

DGG.nbasefaces(::ISEA4RSystem) = NDIAMONDS
DGG.systemname(::ISEA4RSystem) = "ISEA4R"
DGG.idname(::ISEA4RSystem) = "ISEA4R id"

DGG.levels(::ISEA4RSystem) = 0:MAX_LEVEL

"""
    maxneighbors(ISEA4RSystem(), connectivity) -> Int

`9` under `Vertex()`, `4` under `Edge()`.

Interior cells have eight vertex neighbours. At vertices 0 and 11, five
diamond corners meet and a cell can have nine. Other valence-3 icosahedron
vertices give seven. At level zero every diamond has six vertex neighbours.
"""
DGG.maxneighbors(::ISEA4RSystem, ::DGG.Vertex) = 9
DGG.maxneighbors(::ISEA4RSystem, ::DGG.Edge) = 4

# `CustomOrder`, not `CounterClockwise`, and no `maxring`: both omissions have
# the same cause. The diamond lattice's flat laws would be `8k` per ring and a
# faithfully rotational shell, but the 5-valent icosahedral vertices break both
# — rings reach 25 and 35 against the flat 24 and 32 at k = 3 and 4, and the
# distorted shells are not rotations of their one-rings, so azimuth is not
# monotone around them. The one-rings themselves do measure counter-clockwise;
# what cannot be carried outward is the order, which is what `CustomOrder` says.
DGG.winding(::ISEA4RSystem, ::DGG.Connectivity) = DGG.CustomOrder()

# Row-major codecs remain internal; only the canonical Morton index is exposed.
DGG.cellindextypes(::ISEA4RSystem) = (DGG.LevelIndex,)

# ===========================================================================
# Geometry
# ===========================================================================

# Eight great-circle segments approximate each curved chart edge. Shared points
# are bit-identical within a diamond; cross-diamond incidence requires tolerance.
#
# Hashed coordinate incidence is valid only within one diamond. Cross-diamond
# boundaries require a tolerance, and boundary coordinates do not contain
# negative zero.
const BOUNDARY_SEGMENTS = 8

"""
    _perimeter_points(ix, iy, diamond, nside, nseg) -> Vector{UnitSphericalPoint}

The cell's boundary walked in `nseg` equal chart steps per edge, starting at the
`(x+, y+)` corner and running in the [`cell_corners`](@ref) order
`(x+,y+) → (x-,y+) → (x-,y-) → (x+,y-)`. `nseg == 1` reproduces
[`cell_corners`](@ref) exactly.
"""
_perimeter_points(ix::Integer, iy::Integer, diamond::Integer,
        nside::Integer, nseg::Integer) =
    DGG.chart_perimeter(xyd_to_point, ix, iy, diamond, nside, nseg)

"""
    cell_boundary(grid, c) -> Vector{UnitSphericalPoint}

Return the implicitly closed boundary counterclockwise from `(x+,y+)`, seen from
outside the sphere. Each curved chart edge uses `BOUNDARY_SEGMENTS`
great-circle segments. Use
[`cell_area`](@ref) for the exact equal-area value.
"""
function DGG.cell_boundary(sys::ISEA4RSystem, c::DGG.LevelIndex)
    nside = DGG.nside(DGG.level(c))
    ix, iy, d = morton_to_xyd(DGG.checked_id(sys, c), nside)
    return _perimeter_points(ix, iy, d, nside, BOUNDARY_SEGMENTS)
end

"""
    cell_area(grid, c) -> Float64

Exact cell area in steradians: `4π/(10*4^level)`. This `O(1)` value is
independent of the approximate boundary polygon.
"""
DGG.cell_area(g::LevelGrid, c::DGG.LevelIndex) =
    (DGG.checked_id(g, c); 4 * Float64(π) / DGG.ncells(g))

"""
    cell_centroid(grid, c) -> UnitSphericalPoint

The cell centre: the chart evaluated at the lattice cell's midpoint
`((ix + 0.5)/nside, (iy + 0.5)/nside)`.

This is the centre *by definition* — the chart is exactly equal-area, so the
midpoint of the chart rectangle is the canonical centre — and it is strictly
interior, as [`cell_centroid`](@ref) requires. It is not the spherical centroid
of the published 4-gon; no equal-area DGGS claims that of its cell centres.
"""
function DGG.cell_centroid(sys::ISEA4RSystem, c::DGG.LevelIndex)
    nside = DGG.nside(DGG.level(c))
    ix, iy, d = morton_to_xyd(DGG.checked_id(sys, c), nside)
    return cell_center(ix, iy, d, nside)
end

# ===========================================================================
# node_extent — the subtree cap
# ===========================================================================

# How many chart samples per edge the subtree cap is built from. See
# `_subtree_cap` for why this number, and not the boundary's, sets the bound.
const CAP_EDGE_SEGMENTS = 8

"""
    _subtree_cap(ix, iy, diamond, nside) -> SphericalCap

Bounding cap for a cell and its subtree. Children tile the parent chart
rectangle exactly, so a cap covering that rectangle covers all descendants;
[`DGG.sampled_cap`](@ref) turns the corner-inclusive perimeter samples into the
radius.
"""
_subtree_cap(ix::Integer, iy::Integer, diamond::Integer, nside::Integer) =
    DGG.sampled_cap(cell_center(ix, iy, diamond, nside),
        _perimeter_points(ix, iy, diamond, nside, CAP_EDGE_SEGMENTS))

"""
    node_extent(ISEA4RSystem(), c) -> SphericalCap

Return the uninflated subtree cap. Exact chart nesting makes the cell rectangle
a bound for every descendant's geometry; [`_subtree_cap`](@ref) supplies a
conservative sampled radius. It bounds no descendant's *cap*: a child cap is
recentred and may reach outside its parent's extent without violating the
covering law. All returned caps are geodesically convex: the widest in the
system is a level-0 diamond's, 62.3°, against the 90° bound.
"""
function DGG.node_extent(::ISEA4RSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    nside = DGG.nside(l)
    ix, iy, d = morton_to_xyd(c.index, nside)
    return _subtree_cap(ix, iy, d, nside)
end

# ===========================================================================
# Location
# ===========================================================================

"""
    cellat(grid, p::UnitSphericalPoint) -> LevelIndex

Return the cell containing `p` via the analytic chart inverse. Complete grids
never return `nothing`. Boundary ties use Snyder's face choice and the
higher-side lattice cell, deterministically per floating-point platform.
"""
DGG.cellat(g::LevelGrid, p::GO.UnitSphericalPoint) =
    DGG.LevelIndex(g.level, point_to_morton(p, DGG.nside(g.level)))

# ===========================================================================
# Topology
# ===========================================================================

"""
    one_ring(grid, c, connectivity) -> SmallVector{9,LevelIndex}

The immediate neighbours of `c` in **counter-clockwise rotational order seen
from outside the sphere, starting at the `(+1, 0)` chart direction**, from
[`lattice_neighbors`](@ref).

The subject's id is validated once, by `DGG.checked_id`; the neighbours are then
encoded through [`xyd_to_morton_unchecked`](@ref), because
[`lattice_neighbors`](@ref) derives them from that validated cell by table
lookup and they cannot be off the lattice. Re-deriving `ispow2(nside)` and two
range tests per neighbour is the whole cost of this function otherwise, and this
is the inner loop of every breadth-first shell walk.
"""
function DGG.one_ring(g::LevelGrid, c::DGG.LevelIndex, connectivity::DGG.Connectivity)
    nside = DGG.nside(g.level)
    ix, iy, d = morton_to_xyd(DGG.checked_id(g, c), nside)
    out = SmallVector{9,DGG.LevelIndex}()
    for (jx, jy, jd) in lattice_neighbors(ix, iy, d, nside, connectivity)
        out = SmallCollections.push(out,
            DGG.LevelIndex(g.level, xyd_to_morton_unchecked(jx, jy, jd, nside)))
    end
    return out
end

"""
    neighbors(grid, c, k = 1; connectivity = Vertex())

Cells within `k` lattice steps, excluding `c`, with rings concatenated outward.
Each ring is counterclockwise seen from outside the sphere, starting at chart
direction `(+1,0)` — a chart direction, not a compass one, so the same slot
points a different way on each diamond. Vertex connectivity has degree 7–9 at
corner cells and 8 in the interior; edge connectivity uses
the four axis offsets. Missing corner slots are omitted and multi-cell vertex
slots follow fan order. Outer rings use the same starting azimuth.

At level zero each diamond has six vertex-neighbors; the 7–9 counts apply to
corner cells at finer levels.

`k == 0` returns an empty container; `k == 1` returns a
`SmallCollections.SmallVector` sized by [`maxneighbors`](@ref).
"""
Base.@constprop :aggressive function DGG.neighbors(g::LevelGrid, c::DGG.LevelIndex, k::Integer = 1;
        connectivity::DGG.Connectivity = DGG.Vertex())
    steps = DGG.checked_steps(k)
    steps == 0 && return SmallVector{9,DGG.LevelIndex}()
    steps == 1 && return DGG.one_ring(g, c, connectivity)
    return DGG.shell_disc(g, c, steps, connectivity)
end

"""
    ring(grid, c, k; connectivity = Vertex())

Cells at lattice distance exactly `k`, counterclockwise from the first ring-1
direction. `k == 0` returns `[c]`; outer-ring azimuth ties use canonical order.
"""
Base.@constprop :aggressive function DGG.ring(g::LevelGrid, c::DGG.LevelIndex, k::Integer;
        connectivity::DGG.Connectivity = DGG.Vertex())
    steps = DGG.checked_steps(k)
    steps == 0 && return DGG.LevelIndex[c]
    steps == 1 && return DGG.one_ring(g, c, connectivity)
    return DGG.shell_ring(g, c, steps, connectivity)
end
