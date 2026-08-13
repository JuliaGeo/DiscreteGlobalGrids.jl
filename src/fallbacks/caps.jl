# ---------------------------------------------------------------------------
# Spherical caps: the one extent vocabulary
#
# Every node of every tree in this package reports a `SphericalCap` (radians,
# on the unit sphere), which is what lets one predicate serve all of them —
# `ConservativeRegridding` dispatches its dual-tree predicate on exactly that
# type for a `Spherical` manifold.
#
# Two properties are load-bearing and are why the guards below exist:
#
#   * a cap of angular radius <= pi/2 is geodesically convex, so a cap that
#     contains a ring's *vertices* contains the great-circle arcs between them.
#     Past pi/2 that argument fails, and a "bounding" cap built from vertices
#     alone can miss the geometry it claims to bound — so we return the full
#     sphere there instead, which is loose but never wrong.
#   * the radius is nudged outward by a relative slack before it is stored, so
#     that a point on the boundary of the cap survives `_contains`' closed-cap
#     arithmetic under rounding.
# ---------------------------------------------------------------------------

# Relative outward nudge on every derived cap radius. 1e-4 rad of a cell radius
# is far below any pruning consequence and six orders of magnitude above the
# ~1e-14 noise of the dot products the radius comes from.
const CAP_SLACK = 1.0001

# How many cells a union cap will look at before giving up and answering the
# full sphere. A batch this big is a coarse node, and a coarse node's own
# `node_extent` is the O(1) answer that exists for exactly this reason.
const UNION_CAP_BATCH_LIMIT = 2048

"""
    full_sphere_cap() -> SphericalCap

The cap that contains everything: the sound answer wherever a tighter bound
cannot be derived. Loose extents cost time; wrong ones cost correctness.
"""
full_sphere_cap() = SphericalCap(USPoint(0.0, 0.0, 1.0), nextfloat(Float64(pi)))

"""
    intersects_cap(a::SphericalCap, b::SphericalCap) -> Bool
    intersects_cap(a::SphericalCap) -> predicate

Whether two caps meet, and — curried — the one-argument predicate the
`SpatialTreeInterface` traversals take. Exists so call sites stop reaching into
`GeometryOps.UnitSpherical._intersects`, which is private.
"""
intersects_cap(a::US.SphericalCap, b::US.SphericalCap) = US._intersects(a, b)
intersects_cap(a::US.SphericalCap) = Base.Fix1(intersects_cap, a)

"""
    cap_contains(cap::SphericalCap, p::UnitSphericalPoint) -> Bool

Whether a point lies in the closed cap.
"""
cap_contains(cap::US.SphericalCap, p) = US._contains(cap, p)

"""
    inflate_cap(cap::SphericalCap, factor::Real) -> SphericalCap

`cap` with its angular radius multiplied by `factor` — the default
[`node_extent`](@ref)'s covering headroom, clamped at the full sphere.
"""
function inflate_cap(cap::US.SphericalCap, factor::Real)
    p = cap.point
    return SphericalCap(USPoint(p[1], p[2], p[3]),
        nextfloat(min(Float64(pi), Float64(cap.radius) * Float64(factor) + 1e-9)))
end

"""
    points_cap(points) -> SphericalCap

Bounding cap of a set of unit-sphere points: centred on the normalised mean
direction, with the radius the farthest point needs.

Because the result is convex only up to `pi/2`, a point set that spreads wider
than that gets the full sphere — the vertex-based radius would no longer bound
the great-circle arcs between the points, and this function's callers all use
it to bound *rings*, not just their vertices.
"""
function points_cap(points)
    isempty(points) && return full_sphere_cap()
    sx = sy = sz = 0.0
    n = 0
    for p in points
        sx += p[1]; sy += p[2]; sz += p[3]
        n += 1
    end
    norm = sqrt(sx * sx + sy * sy + sz * sz)
    norm <= eps(Float64) && return full_sphere_cap()
    center = USPoint(sx / norm, sy / norm, sz / norm)
    radius = 0.0
    for p in points
        radius = max(radius, US.spherical_distance(center, USPoint(p[1], p[2], p[3])))
    end
    radius > Float64(pi) / 2 && return full_sphere_cap()
    return SphericalCap(center, nextfloat(min(Float64(pi), radius * CAP_SLACK + 1e-12)))
end

"""
    cell_cap(grid, c) -> SphericalCap

The **tight** bounding cap of one cell — its boundary ring's own cap, with no
covering headroom. This is a cell's extent as a *leaf*: it bounds the cell and
nothing else, which is the tightest sound thing a tree can prune a leaf with.

It is deliberately not [`node_extent`](@ref), which bounds the cell's whole
subtree and is therefore inflated in systems where children overhang their
parent.
"""
cell_cap(grid::AbstractGrid, c::AbstractCellIndex) = points_cap(cell_boundary(grid, c))

"""
    cells_cap(grid, ids) -> SphericalCap

Bounding cap of a batch of cells: the exact union cap of their boundary
vertices, or the full sphere once the batch passes `UNION_CAP_BATCH_LIMIT`.

This is what a *sparse* tree node uses — one that stores far fewer leaves than
its cell's subtree has, where a cap around what it really owns is much tighter
than the cell's own, and the boundary calls are bounded by the caller.
"""
function cells_cap(grid::AbstractGrid, ids)
    points = USPoint[]
    count = 0
    for c in ids
        count += 1
        count > UNION_CAP_BATCH_LIMIT && return full_sphere_cap()
        for p in cell_boundary(grid, c)
            push!(points, USPoint(p[1], p[2], p[3]))
        end
    end
    return points_cap(points)
end

"""
    merge_caps(a::SphericalCap, b::SphericalCap) -> SphericalCap

A cap containing both. Used to build the fallback tree's internal extents
bottom-up, where the leaf caps are in hand but their vertices are not.
"""
merge_caps(a::US.SphericalCap, b::US.SphericalCap) = US._merge(a, b)

"""
    unit_point(lon, lat) -> UnitSphericalPoint

Longitude/latitude **in degrees** onto the unit sphere. The only place in this
package where degrees turn into geometry; everything above it is `xyz`.
"""
function unit_point(lon::Real, lat::Real)
    # `Float64` unconditionally: a Float32-coordinate geometry (GeoJSON reads
    # Natural Earth that way) would otherwise carry its element type into every
    # downstream container and miss the concrete struct fields below.
    lon64, lat64 = Float64(lon), Float64(lat)
    coslat = cosd(lat64)
    return USPoint(coslat * cosd(lon64), coslat * sind(lon64), sind(lat64))
end

"""
    lonlat(p) -> (lon, lat)

A unit-sphere point back to longitude/latitude **in degrees**.
"""
lonlat(p) = (atand(p[2], p[1]), asind(clamp(p[3], -1.0, 1.0)))

# A GeoInterface point as a unit-sphere point: 3-D coordinates are xyz as-is,
# 2-D coordinates are lon/lat degrees. The same convention the spherical
# RelateNG kernel ingests vertices with.
query_point(p) = GI.is3d(p) ?
                 USPoint(Float64(GI.x(p)), Float64(GI.y(p)), Float64(GI.z(p))) :
                 unit_point(GI.x(p), GI.y(p))
