# ---------------------------------------------------------------------------
# Geometry and location
#
# All three primitives are libh3's own answers, converted to the unit sphere and
# otherwise untouched. That is deliberate: H3's geometry is defined by libh3, so
# anything this file "improved" would be a different grid wearing the same name.
# ---------------------------------------------------------------------------

"""
    cell_boundary(::H3System, c::H3Cell) -> SmallVector{10,UnitSphericalPoint}

The exact boundary ring of `c`: libh3's `cellToBoundary`, on the unit sphere.

Implicitly closed and counter-clockwise seen from outside, as the base
interface requires — both are libh3's own conventions, verified in this
system's test suite rather than assumed.

# Why up to ten vertices

A hexagon has six vertices and a pentagon five, *until* the cell crosses an
edge of the underlying icosahedron. There the gnomonic chart changes face and
libh3 emits extra **distortion vertices** at the crossing, up to ten in all.

They are kept. They are where the cell's boundary actually bends, so dropping
them — as this port's predecessor did, to keep every ring convex for a
spherical clipper that no longer participates — moves the boundary and costs
real area: the old cleanup admitted errors up to 12% of a cell, while the
untouched ring reproduces `cellAreaRads2` to about 1e-15 sr. The container is
sized for the ten libh3 can produce and carries however many it did.
"""
function cell_boundary(::H3System, c::H3Cell)
    verts, n = H3Native.boundary_verts(c.id)
    out = SmallVector{10,USPoint}()
    for i in 1:n
        v = @inbounds verts[i]
        out = SmallCollections.push(out, USPoint(v[1], v[2], v[3]))
    end
    return out
end

"""
    cell_centroid(::H3System, c::H3Cell) -> UnitSphericalPoint

The cell centre from libh3's `cellToLatLng`, straight onto the unit sphere.

This is the centre H3's hierarchy is built around — the point child cells are
arranged about — rather than the mean of the boundary vertices, and it is
strictly interior to the cell for every valid index.
"""
function cell_centroid(::H3System, c::H3Cell)
    x, y, z = H3Native.cell_center_cartesian(c.id)
    return USPoint(x, y, z)
end

# ===========================================================================
# Location — the marquee fast path
# ===========================================================================

"""
    cellat(grid::LevelGrid, p::UnitSphericalPoint) -> H3Cell
    cellat(grid::LevelGrid, lon::Real, lat::Real) -> H3Cell

The cell of `grid` containing a point, by libh3's `latLngToCell`.

This is a closed-form inverse projection — pick the icosahedron face, project,
round to the nearest cell centre in that face's hex lattice — so it replaces
the generic tree-descent-plus-point-in-polygon fallback with an O(1) call that
does not touch a single boundary.

Never `nothing`: H3 tessellates the whole sphere, so every point is in a cell.

**Ties.** A point exactly on a shared edge resolves the way `latLngToCell`
resolves it, which is deterministic and is the same cell every other H3
implementation names for that point. This system deliberately does not impose
its own tie rule on top: agreeing with libh3 is worth more than agreeing with
a convention no other H3 binding follows.

The `(lon, lat)` method takes **degrees** and is the primitive here — it hands
the coordinates to libh3 directly rather than routing them through `xyz` and
back, which is both faster and one rounding step shorter.
"""
function cellat(grid::LevelGrid, p::GO.UnitSphericalPoint)
    lon = atand(p[2], p[1])
    lat = asind(clamp(p[3], -1.0, 1.0))
    return H3Cell(H3Native.lonlat_to_cell(lon, lat, grid.level))
end

cellat(grid::LevelGrid, lon::Real, lat::Real) =
    H3Cell(H3Native.lonlat_to_cell(lon, lat, grid.level))
