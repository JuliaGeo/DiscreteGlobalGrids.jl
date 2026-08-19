# Nested HEALPix uses dense 0-based Morton ids: level `l` is
# `0:12*4^l-1`, children are `4p:4p+3`, and the parent is `p ÷ 4` — the quad-face
# family's arithmetic, which this file declares itself a member of rather than
# rewriting. Grid positions and alternate RING indices are 1-based. Location,
# topology, and subtree extents have direct chart or lattice implementations.

# ===========================================================================
# Types
# ===========================================================================

"""
    HEALPixSystem() <: AbstractQuadFaceGridSystem

HEALPix on the unit sphere in nested Morton order. Its twelve equal-area base
pixels refine by aperture 4, giving `12 * 4^l` cells of area
`4π / (12 * 4^l)` at level `l`. Chart edges are not great circles, so
[`cell_boundary`](@ref) densifies them.

Canonical [`LevelIndex`](@ref) ids are 0-based; grid positions and alternate
[`HEALPixRingIndex`](@ref) ids are 1-based. Supported levels are `0:29`, the
largest range whose cell count fits `Int64`. Nested ordering makes subtrees
contiguous. Maximum neighbour counts are 8 for `Vertex()` and 4 for `Edge()`.
[`node_extent`](@ref) returns a direct subtree cap without generic inflation.
"""
struct HEALPixSystem <: DGG.AbstractQuadFaceGridSystem end

# Grid descriptor for all `12 * 4^level` pixels in nested order.
const LevelGrid = DGG.HierarchicalLevelGrid{HEALPixSystem}

"""
    HEALPixRingIndex(level, index) <: AbstractCellIndex

A 1-based HEALPix RING index, ordered north-to-south by iso-latitude ring and
west-to-east within each ring. This alternate scheme is compatible with
ring-ordered data vectors but does not preserve subtree order. Convert with
[`reindex`](@ref).
"""
struct HEALPixRingIndex <: DGG.AbstractCellIndex
    level::Int32
    index::Int64
end

DGG.level(c::HEALPixRingIndex) = Int(c.level)
DGG.rawid(c::HEALPixRingIndex) = c.index
Base.isless(a::HEALPixRingIndex, b::HEALPixRingIndex) =
    isless((a.level, a.index), (b.level, b.index))
Base.show(io::IO, c::HEALPixRingIndex) =
    print(io, "HEALPixRingIndex(", c.level, ", ", c.index, ")")

const MAX_LEVEL = 29

# ===========================================================================
# System interface
# ===========================================================================

# The twelve base pixels, which `rootcells` names `LevelIndex(0, 0:11)`.
DGG.nbasefaces(::HEALPixSystem) = 12
DGG.systemname(::HEALPixSystem) = "HEALPix"
DGG.idname(::HEALPixSystem) = "nested id"

DGG.cellindextypes(::HEALPixSystem) = (DGG.LevelIndex, HEALPixRingIndex)
DGG.levels(::HEALPixSystem) = 0:MAX_LEVEL

DGG.maxneighbors(::HEALPixSystem, ::DGG.Vertex) = 8
DGG.maxneighbors(::HEALPixSystem, ::DGG.Edge) = 4

# ===========================================================================
# Alternate index scheme
# ===========================================================================

function DGG.reindex(::Type{HEALPixRingIndex}, ::HEALPixSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    return HEALPixRingIndex(l, nested_to_ring(c.index, DGG.nside(l)))
end

function DGG.reindex(::Type{DGG.LevelIndex}, ::HEALPixSystem, c::HEALPixRingIndex)
    l = DGG.level(c)
    return DGG.LevelIndex(l, ring_to_nested(c.index, DGG.nside(l)))
end

# ===========================================================================
# Geometry
# ===========================================================================

# HEALPix chart edges are not great circles. Eight great-circle segments per
# edge give 32 boundary vertices and about 0.18% level-independent relative
# area error. A power-of-two count also makes shared-edge sample arguments, and
# therefore vertices, bit-identical across adjacent cells. The relative error
# falls as `nseg^-2`.
const BOUNDARY_SEGMENTS = 8

"""
    _perimeter_points(ix, iy, face, nside, nseg) -> Vector{UnitSphericalPoint}

Return `nseg` chart samples per edge in north → west → south → east order.
The implicitly closed ring omits duplicated edge endpoints; `nseg == 1`
reproduces [`pixel_corners`](@ref).
"""
_perimeter_points(ix::Integer, iy::Integer, face::Integer,
        nside::Integer, nseg::Integer) =
    DGG.chart_perimeter(xyf_to_point, ix, iy, face, nside, nseg)

"""
    cell_boundary(grid, c) -> Vector{UnitSphericalPoint}

Return the implicitly closed boundary counter-clockwise from the north corner,
in north → west → south → east order. Each chart edge is represented by
`BOUNDARY_SEGMENTS` great-circle segments; the four corners occur at vertices
1, 9, 17, and 25. Use [`cell_area`](@ref) for the exact equal-area value.
"""
function DGG.cell_boundary(sys::HEALPixSystem, c::DGG.LevelIndex)
    nside = DGG.nside(DGG.level(c))
    ix, iy, face = nested_to_xyf(DGG.checked_id(sys, c), nside)
    return _perimeter_points(ix, iy, face, nside, BOUNDARY_SEGMENTS)
end

"""
    cell_area(grid, c) -> Float64

Return the exact equal-area solid angle `4π / (12 * 4^level)` in O(1). This is
independent of boundary densification; the 32-vertex polygon underestimates it
by about 0.18%.
"""
DGG.cell_area(g::LevelGrid, c::DGG.LevelIndex) =
    (DGG.checked_id(g, c); 4 * Float64(π) / DGG.ncells(g))

"""
    cell_centroid(grid, c) -> UnitSphericalPoint

Return the canonical, strictly interior pixel centre: the chart evaluated at
the lattice-cell midpoint.
"""
function DGG.cell_centroid(sys::HEALPixSystem, c::DGG.LevelIndex)
    nside = DGG.nside(DGG.level(c))
    ix, iy, face = nested_to_xyf(DGG.checked_id(sys, c), nside)
    return pixel_center(ix, iy, face, nside)
end

# ===========================================================================
# node_extent — the subtree cap
# ===========================================================================

# Chart samples per edge used to bound a subtree.
const CAP_EDGE_SEGMENTS = 8

"""
    _subtree_cap(ix, iy, face, nside) -> SphericalCap

Return a cap for the pixel and its complete subtree. Nested refinement exactly
subdivides the parent's chart square, so bounding that square bounds every
descendant in O(1); [`DGG.sampled_cap`](@ref) turns the corner-inclusive
perimeter samples into the radius.
"""
_subtree_cap(ix::Integer, iy::Integer, face::Integer, nside::Integer) =
    DGG.sampled_cap(pixel_center(ix, iy, face, nside),
        _perimeter_points(ix, iy, face, nside, CAP_EDGE_SEGMENTS))

"""
    node_extent(HEALPixSystem(), c) -> SphericalCap

Return the pixel's subtree cap. Nested children exactly partition the parent,
so no generic inflation is required. `_subtree_cap` derives its radius from
the corner-inclusive perimeter samples plus conservative slack.
"""
function DGG.node_extent(::HEALPixSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    nside = DGG.nside(l)
    ix, iy, face = nested_to_xyf(c.index, nside)
    return _subtree_cap(ix, iy, face, nside)
end

# ===========================================================================
# Location
# ===========================================================================

"""
    cellat(grid, p::UnitSphericalPoint) -> LevelIndex

Return the pixel containing `p` via the chart's closed-form inverse. A complete
level grid covers the sphere, so this never returns `nothing`. Boundary ties
use `point_to_xyf`'s deterministic higher-side `floor` convention. Other
HEALPix implementations may choose a different valid cell at shared borders.
"""
DGG.cellat(g::LevelGrid, p::GO.UnitSphericalPoint) =
    DGG.LevelIndex(g.level, point_to_nested(p, g.level))

# ===========================================================================
# Topology
# ===========================================================================

"""
    one_ring(grid, c, connectivity) -> SmallVector{8,LevelIndex}

Return immediate neighbours counter-clockwise from `SW`, as seen from outside
the sphere. Missing entries are omitted and level-0 duplicates keep their
first occurrence to preserve the cycle.
"""
function DGG.one_ring(g::LevelGrid, c::DGG.LevelIndex, connectivity::DGG.Connectivity)
    DGG.checked_id(g, c)
    raw = nested_neighbors(c.index, g.level)
    out = SmallVector{8,DGG.LevelIndex}()
    for m in _neighbor_cycle(connectivity)
        id = raw[m]
        id < 0 && continue
        nb = DGG.LevelIndex(g.level, id)
        nb in out && continue
        out = SmallCollections.push(out, nb)
    end
    return out
end

"""
    neighbors(grid, c, k = 1; connectivity = Vertex())

Return rings `1:k` concatenated outward, excluding `c`. Each ring is ordered
counter-clockwise as seen from outside the sphere. `Vertex()` uses all eight
lattice directions, except the missing neighbour at degree-3 vertices;
`Edge()` uses the edge-sharing `SW, NW, NE, SE` directions. Every ring starts on
the `SW` spoke; outer rings use azimuth about the cell centre.

`k == 0` returns an empty container. `k == 1` returns a fixed-capacity
`SmallVector` without allocation.
"""
function DGG.neighbors(g::LevelGrid, c::DGG.LevelIndex, k::Integer = 1;
        connectivity::DGG.Connectivity = DGG.Vertex())
    steps = DGG.checked_steps(k)
    steps == 0 && return SmallVector{8,DGG.LevelIndex}()
    steps == 1 && return DGG.one_ring(g, c, connectivity)
    shells = DGG.adjacency_shells(g, c, steps, connectivity)
    isempty(shells) && return DGG.LevelIndex[]
    return reduce(vcat, shells)
end

"""
    ring(grid, c, k; connectivity = Vertex())

Return cells at lattice distance exactly `k`, counter-clockwise as seen from
outside the sphere. `k == 0` returns `[c]`; `k == 1` uses the lattice cycle.
Outer rings are sorted by azimuth about the cell centre from the spoke through
the `SW` neighbour, with canonical ids breaking ties.
"""
function DGG.ring(g::LevelGrid, c::DGG.LevelIndex, k::Integer;
        connectivity::DGG.Connectivity = DGG.Vertex())
    steps = DGG.checked_steps(k)
    steps == 0 && return DGG.LevelIndex[c]
    steps == 1 && return DGG.one_ring(g, c, connectivity)
    shells = DGG.adjacency_shells(g, c, steps, connectivity)
    steps <= length(shells) || return DGG.LevelIndex[]
    return shells[steps]
end
