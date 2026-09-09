# `cell_boundary_cartesian` uses A5's internal Cartesian frame. The public
# geometry methods use geographic longitude and geodetic latitude.

# The projection itself is the authalic one: `A5Native._from_lonlat` converts
# geodetic latitude to authalic before projecting, and `_to_lonlat` converts
# back. That conversion is already inside every method here, so the geometry
# below is geodetic and `AuthalicSystem` must refuse it rather than convert a
# second time.
DGG.Fallbacks.publishes_geodetic_geometry(::A5System) = true

# Convert geographic longitude and latitude in degrees to a unit vector.
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

The returned vertex count is `corners × segments`: five corners at level 0,
three at level 1, and five at deeper levels.
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
# Location
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
