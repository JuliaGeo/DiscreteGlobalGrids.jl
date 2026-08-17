# Nested HEALPix uses dense 0-based Morton ids: level `l` is
# `0:12*4^l-1`, children are `4p:4p+3`, and the parent is `p ÷ 4`.
# Grid positions and alternate RING indices are 1-based. Hierarchy operations
# use integer arithmetic; location, topology, and subtree extents have direct
# chart or lattice implementations.

# ===========================================================================
# Types
# ===========================================================================

"""
    HEALPixSystem() <: AbstractHierarchicalGridSystem

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
struct HEALPixSystem <: DGG.AbstractHierarchicalGridSystem end

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

# Convert a validated level to its face side and pixel count.
@inline _nside(level::Integer) = Int64(1) << Int(level)
@inline _npix(level::Integer) = 12 * (Int64(1) << (2 * Int(level)))

const MAX_LEVEL = 29

# ===========================================================================
# System interface
# ===========================================================================

DGG.cellindextype(::HEALPixSystem) = DGG.LevelIndex
DGG.cellindextypes(::HEALPixSystem) = (DGG.LevelIndex, HEALPixRingIndex)
DGG.levels(::HEALPixSystem) = 0:MAX_LEVEL
DGG.has_sorted_subtrees(::HEALPixSystem) = true

DGG.max_neighbors(::HEALPixSystem, ::DGG.Vertex) = 8
DGG.max_neighbors(::HEALPixSystem, ::DGG.Edge) = 4

DGG.rootcells(::HEALPixSystem) = [DGG.LevelIndex(0, i) for i in 0:11]

"""
    parent(HEALPixSystem(), c) -> LevelIndex

The nested parent: `index ÷ 4`, one level up. Throws an `ArgumentError` on a
level-0 cell, which has no parent.
"""
function Base.parent(::HEALPixSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    l > 0 || throw(ArgumentError(
        "level-0 HEALPix cell $c is a root and has no parent"))
    return DGG.LevelIndex(l - 1, c.index >> 2)
end

"""
    children(HEALPixSystem(), c)

Return the four ascending nested children `4*index .+ (0:3)`. Throws an
`ArgumentError` at `max_level`.
"""
function DGG.children(sys::HEALPixSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    l < DGG.max_level(sys) || throw(ArgumentError(
        "HEALPix cell $c is at max_level $(DGG.max_level(sys)) and has no children"))
    base = c.index << 2
    return [DGG.LevelIndex(l + 1, base + k) for k in 0:3]
end

"""
    ancestor(HEALPixSystem(), c, l) -> LevelIndex

Return the level-`l` ancestor by dropping `2Δ` low Morton bits:
`index >> 2Δ`.
"""
function DGG.ancestor(sys::HEALPixSystem, c::DGG.LevelIndex, l::Integer)
    target = Int(l)
    lc = DGG.level(c)
    target <= lc || throw(ArgumentError(
        "ancestor level $target is deeper than the cell's own level $lc"))
    target >= 0 || throw(ArgumentError(
        "ancestor level $target is above the root level 0"))
    return DGG.LevelIndex(target, c.index >> (2 * (lc - target)))
end

"""
    descendant_range(HEALPixSystem(), c, l) -> UnitRange{Int}

Return the contiguous 1-based positions of `c`'s level-`l` descendants. For
`Δ = l - level(c)`, their 0-based ids are
`index * 4^Δ : (index + 1) * 4^Δ - 1`. The range is exact and hole-free.
"""
function DGG.descendant_range(sys::HEALPixSystem, c::DGG.LevelIndex, l::Integer)
    target = Int(l)
    lc = DGG.level(c)
    target >= lc || throw(ArgumentError(
        "descendant level $target is above the cell's own level $lc"))
    target <= DGG.max_level(sys) || throw(ArgumentError(
        "descendant level $target is past max_level $(DGG.max_level(sys))"))
    shift = 2 * (target - lc)
    lo = c.index << shift
    hi = ((c.index + 1) << shift) - 1
    return Int(lo + 1):Int(hi + 1)
end

"""
    descendants(HEALPixSystem(), c, l)

Return every level-`l` descendant of `c` in ascending nested order. Dense,
subtree-contiguous ids permit direct conversion from [`descendant_range`](@ref).
"""
function DGG.descendants(sys::HEALPixSystem, c::DGG.LevelIndex, l::Integer)
    r = DGG.descendant_range(sys, c, l)      # validates `l` both ways
    target = Int(l)
    return [DGG.LevelIndex(target, i - 1) for i in r]
end

# ===========================================================================
# The level grid: size, and positions <-> ids
# ===========================================================================

DGG.ncells(::HEALPixSystem, l::Integer) = Int(_npix(l))

# The grid bounds-checks `i`, so this is the bijection and nothing else.
DGG.cellindex(::HEALPixSystem, l::Integer, i::Int) = DGG.LevelIndex(l, i - 1)

"""
    cellposition(HEALPixSystem(), c) -> Union{Int,Nothing}

Return `index + 1` for an in-range nested id, or `nothing` otherwise. The grid
must reject cells from another level and convert [`HEALPixRingIndex`](@ref)
values before calling this method.
"""
function DGG.cellposition(::HEALPixSystem, c::DGG.LevelIndex)
    0 <= c.index < _npix(DGG.level(c)) || return nothing
    return Int(c.index + 1)
end

# ===========================================================================
# Alternate index scheme
# ===========================================================================

function DGG.reindex(::Type{HEALPixRingIndex}, ::HEALPixSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    return HEALPixRingIndex(l, nested_to_ring(c.index, _nside(l)))
end

function DGG.reindex(::Type{DGG.LevelIndex}, ::HEALPixSystem, c::HEALPixRingIndex)
    l = DGG.level(c)
    return DGG.LevelIndex(l, ring_to_nested(c.index, _nside(l)))
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
function _perimeter_points(ix::Integer, iy::Integer, face::Integer,
        nside::Integer, nseg::Integer)
    n = nside
    x0 = Int64(ix)
    y0 = Int64(iy)
    pts = Vector{GO.UnitSphericalPoint{Float64}}(undef, 4 * nseg)
    k = 0
    for i in 0:(nseg - 1)          # north -> west, along y = (iy+1)/n
        t = i / nseg
        pts[k += 1] = xyf_to_point((x0 + 1 - t) / n, (y0 + 1) / n, face)
    end
    for i in 0:(nseg - 1)          # west -> south, along x = ix/n
        t = i / nseg
        pts[k += 1] = xyf_to_point(x0 / n, (y0 + 1 - t) / n, face)
    end
    for i in 0:(nseg - 1)          # south -> east, along y = iy/n
        t = i / nseg
        pts[k += 1] = xyf_to_point((x0 + t) / n, y0 / n, face)
    end
    for i in 0:(nseg - 1)          # east -> north, along x = (ix+1)/n
        t = i / nseg
        pts[k += 1] = xyf_to_point((x0 + 1) / n, (y0 + t) / n, face)
    end
    return pts
end

"""
    cell_boundary(grid, c) -> Vector{UnitSphericalPoint}

Return the implicitly closed boundary counter-clockwise from the north corner,
in north → west → south → east order. Each chart edge is represented by
`BOUNDARY_SEGMENTS` great-circle segments; the four corners occur at vertices
1, 9, 17, and 25. Use [`cell_area`](@ref) for the exact equal-area value.
"""
function DGG.cell_boundary(::HEALPixSystem, c::DGG.LevelIndex)
    nside = _nside(DGG.level(c))
    ix, iy, face = nested_to_xyf(_checked_index(c), nside)
    return _perimeter_points(ix, iy, face, nside, BOUNDARY_SEGMENTS)
end

"""
    cell_area(grid, c) -> Float64

Return the exact equal-area solid angle `4π / (12 * 4^level)` in O(1). This is
independent of boundary densification; the 32-vertex polygon underestimates it
by about 0.18%.
"""
DGG.cell_area(g::LevelGrid, c::DGG.LevelIndex) =
    (_checked_index(g, c); 4 * Float64(π) / _npix(g.level))

"""
    cell_centroid(grid, c) -> UnitSphericalPoint

Return the canonical, strictly interior pixel centre: the chart evaluated at
the lattice-cell midpoint.
"""
function DGG.cell_centroid(::HEALPixSystem, c::DGG.LevelIndex)
    nside = _nside(DGG.level(c))
    ix, iy, face = nested_to_xyf(_checked_index(c), nside)
    return pixel_center(ix, iy, face, nside)
end

# Validate a nested id before decoding it into face coordinates.
@inline function _checked_index(c::DGG.LevelIndex)
    l = DGG.level(c)
    0 <= c.index < _npix(l) || throw(ArgumentError(
        "nested id $(c.index) is out of range 0:$(_npix(l) - 1) at level $l"))
    return c.index
end

# Also require the cell and grid to have the same level.
@inline function _checked_index(g::LevelGrid, c::DGG.LevelIndex)
    DGG.level(c) == g.level || throw(ArgumentError(
        "cell $c is at level $(DGG.level(c)), not the grid's level $(g.level)"))
    return _checked_index(c)
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
descendant in O(1).

The cap is centred on the pixel centre. The maximum distance over the square
occurs at a corner, and all corners are sampled, so `rmax` supplies the bound.
`gap/2` adds conservative measured slack. It is not itself a formal Lipschitz
bound because `gap` is a geodesic chord rather than chart-edge arc length.
"""
function _subtree_cap(ix::Integer, iy::Integer, face::Integer, nside::Integer)
    center = pixel_center(ix, iy, face, nside)
    pts = _perimeter_points(ix, iy, face, nside, CAP_EDGE_SEGMENTS)
    rmax = 0.0
    gap = 0.0
    prev = pts[end]
    for p in pts
        rmax = max(rmax, US.spherical_distance(center, p))
        gap = max(gap, US.spherical_distance(prev, p))
        prev = p
    end
    radius = min(Float64(π), rmax + gap / 2)
    return SphericalCap(center, nextfloat(radius))
end

"""
    node_extent(HEALPixSystem(), c) -> SphericalCap

Return the pixel's subtree cap. Nested children exactly partition the parent,
so no generic inflation is required. `_subtree_cap` derives its radius from
the corner-inclusive perimeter samples plus conservative slack.
"""
function DGG.node_extent(::HEALPixSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    nside = _nside(l)
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
    _checked_index(g, c)
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
