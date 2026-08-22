# Tree extents are unit-sphere `SphericalCap`s. Vertex-derived ring caps are
# valid only through radius `π/2`; larger candidates become full-sphere caps.
# Derived radii are nudged outward to preserve closed-cap containment.

# Relative outward nudge above dot-product rounding noise.
const CAP_SLACK = 1.0001

# How many cells a union cap will look at before giving up and answering the
# full sphere. A batch this big is a coarse node, and a coarse node's own
# `node_extent` is the O(1) answer that exists for exactly this reason.
const UNION_CAP_BATCH_LIMIT = 2048

"""
    full_sphere_cap() -> SphericalCap

The cap containing the full sphere.
"""
full_sphere_cap() = SphericalCap(USPoint(0.0, 0.0, 1.0), nextfloat(Float64(pi)))

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

Bound unit-sphere points by a cap centred on their normalized mean. Return the
full sphere if the vertex radius exceeds `π/2`, where it would not necessarily
contain the great-circle arcs between vertices.
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

The tight cap of one cell boundary, without subtree headroom. Unlike
[`node_extent`](@ref), it does not cover descendants that overhang the cell.

It bounds `c`'s own geometry and nothing else. In particular it is not itself
bounded by the `node_extent` of an ancestor of `c` — the covering law is a
statement about descendant geometry, not about descendant caps — so comparing
this cap against an ancestor's extent tests nothing.
"""
cell_cap(grid::AbstractGrid, c::AbstractCellIndex) = points_cap(cell_boundary(grid, c))

"""
    cells_cap(grid, ids) -> SphericalCap

Bound a batch of cell boundaries. Return the full sphere after
`UNION_CAP_BATCH_LIMIT` cells.

The result bounds the geometry of the listed cells only. It covers neither their
descendants nor any cap derived from them, so it is not a [`node_extent`](@ref)
for any cell whose subtree reaches below `grid`'s level.
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
    unit_point(lon, lat) -> UnitSphericalPoint

Longitude/latitude **in degrees** onto the unit sphere. The only place in this
package where degrees turn into geometry; everything above it is `xyz`.
"""
function unit_point(lon::Real, lat::Real)
    # Normalize input coordinate types for concrete downstream containers.
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
