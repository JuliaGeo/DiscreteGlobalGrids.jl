# ---------------------------------------------------------------------------
# Geometry and location
#
# Geometry methods preserve libh3's cell definitions while converting results
# to unit-sphere points.
# ---------------------------------------------------------------------------

"""
    cell_boundary(::H3System, c::H3Cell) -> Helpers.SmallList{10,UnitSphericalPoint}

The libh3 `cellToBoundary` ring on the unit sphere, implicitly closed and
counter-clockwise seen from outside.

Cells crossing an icosahedron edge include distortion vertices, so the ring has
5 to 10 vertices. These vertices are part of the exact boundary and are retained.
"""
function cell_boundary(::H3System, c::H3Cell)
    verts, n = H3Native.boundary_verts(c.id)
    out = DGG.Helpers.empty_small_list(Val(10), USPoint(1.0, 0.0, 0.0))
    for i in 1:n
        v = @inbounds verts[i]
        out = DGG.Helpers.small_push(out, USPoint(v[1], v[2], v[3]))
    end
    return out
end

"""
    cell_centroid(::H3System, c::H3Cell) -> UnitSphericalPoint

The cell centre from libh3's `cellToLatLng`, converted to a unit-sphere point.
It is not the mean of the boundary vertices.
"""
function cell_centroid(::H3System, c::H3Cell)
    x, y, z = H3Native.cell_center_cartesian(c.id)
    return USPoint(x, y, z)
end

# ===========================================================================
# Location
# ===========================================================================

"""
    cellat(grid::LevelGrid, p::UnitSphericalPoint) -> H3Cell
    cellat(grid::LevelGrid, lon::Real, lat::Real) -> H3Cell

The cell containing a point, computed by libh3's O(1) `latLngToCell` inverse.

Never `nothing`: H3 tessellates the whole sphere, so every point is in a cell.

Shared-edge ties follow libh3's deterministic rule.

The `(lon, lat)` overload takes degrees.
"""
function cellat(grid::LevelGrid, p::GO.UnitSphericalPoint)
    lon = atand(p[2], p[1])
    lat = asind(clamp(p[3], -1.0, 1.0))
    return H3Cell(H3Native.lonlat_to_cell(lon, lat, grid.level))
end

cellat(grid::LevelGrid, lon::Real, lat::Real) =
    H3Cell(H3Native.lonlat_to_cell(lon, lat, grid.level))
