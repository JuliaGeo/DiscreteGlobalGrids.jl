# Geometry uses A5's geographic lon/lat frame. Do not mix it with the native
# internal Cartesian frame, which omits longitude and latitude conversions:
# `cell_boundary_cartesian` stops in that frame, while `cell_boundary`,
# `cell_to_lonlat` and `lonlat_to_cell` are geographic. Mixing the two puts ring
# and centroid 93 degrees apart in longitude.

# The one degrees -> unit-sphere conversion this module makes.
function _unit_point(lon::Real, lat::Real)
    λ = deg2rad(Float64(lon))
    φ = deg2rad(Float64(lat))
    cosφ = cos(φ)
    return USPoint(cosφ * cos(λ), cosφ * sin(λ), sin(φ))
end

"""
    cell_boundary(::A5System, c::A5Cell) -> Vector{UnitSphericalPoint}

The A5 boundary on the unit sphere, implicitly closed and counter-clockwise seen
from outside.

Edges are straight in the face plane, not great circles. A5's `:auto`
subdivision uses `2^(6-level)` segments per edge through level 6 and one segment
afterward.

So a ring has `corners × segments` vertices: 5 corners at level 0, 3 at level 1
(a quintant is a triangular slice of a face) and 5 below, times 64, 32, 16, …
down to 1 from level 6.
"""
function cell_boundary(::A5System, c::A5Cell)
    ring = A5Native.cell_boundary(c.id; closed_ring=false, segments=:auto)
    out = Vector{USPoint}(undef, length(ring))
    for (i, p) in enumerate(ring)
        @inbounds out[i] = _unit_point(p[1], p[2])
    end
    return out
end

"""
    cell_centroid(::A5System, c::A5Cell) -> UnitSphericalPoint

The interior face-plane polygon centre, inverse-projected with A5's
`cell_to_lonlat`, strictly interior to the cell. Two level-0 centres are
geographic poles, where a lon/lat east/north tangent frame does not exist.
"""
function cell_centroid(::A5System, c::A5Cell)
    lon, lat = A5Native.cell_to_lonlat(c.id)
    return _unit_point(lon, lat)
end

# ===========================================================================
# Location — the marquee fast path
# ===========================================================================

"""
    cellat(grid::LevelGrid, p::UnitSphericalPoint) -> A5Cell
    cellat(grid::LevelGrid, lon::Real, lat::Real) -> A5Cell

The cell containing a point, computed by A5's O(1) `lonlat_to_cell` inverse.

Never `nothing`: A5 tessellates the whole sphere, so every point is in a cell.

Shared-boundary ties follow A5's deterministic inverse-projection rule.

The `(lon, lat)` overload takes degrees.
"""
function cellat(grid::LevelGrid, lon::Real, lat::Real)
    return A5Cell(A5Native.lonlat_to_cell(lon, lat, grid.level))
end

function cellat(grid::LevelGrid, p::GO.UnitSphericalPoint)
    lon = atand(p[2], p[1])
    lat = asind(clamp(p[3], -1.0, 1.0))
    return cellat(grid, lon, lat)
end
