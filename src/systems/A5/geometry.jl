# ---------------------------------------------------------------------------
# Geometry and location
#
# All three primitives are a5's own answers, converted to the unit sphere and
# otherwise untouched. That is deliberate: A5's geometry is defined by upstream
# a5, so anything this file "improved" would be a different grid wearing the
# same name.
#
# THE FRAME, which is the one thing this file could get wrong and did once
# before. `A5Native` has two spherical frames: the internal one `_to_cartesian`
# produces, which applies neither the `LONGITUDE_OFFSET` un-rotation nor the
# authalic -> geodetic latitude conversion, and the geographic one `_to_lonlat`
# produces, which applies both. `cell_boundary_cartesian` stops in the first;
# `cell_boundary` and `cell_to_lonlat` land in the second. Wiring the ring to
# the former while the centroid used the latter put the two 93 degrees apart in
# longitude — leaf caps at res 3 came out at 1.168 rad instead of 0.144, subtree
# containment became meaningless, and every cross-system regrid was silently
# misaligned, while every test that looked at the boundary *alone* still passed.
# Both ends of this file therefore go through the lon/lat surface, which is also
# the frame `lonlat_to_cell` reads, so `cellat` and the geometry agree.
# ---------------------------------------------------------------------------

# The one degrees -> unit-sphere conversion this module makes.
function _unit_point(lon::Real, lat::Real)
    λ = deg2rad(Float64(lon))
    φ = deg2rad(Float64(lat))
    cosφ = cos(φ)
    return USPoint(cosφ * cos(λ), cosφ * sin(λ), sin(φ))
end

"""
    cell_boundary(::A5System, c::A5Cell) -> Vector{UnitSphericalPoint}

The exact boundary ring of `c`: a5's own `cell_boundary`, on the unit sphere.

Implicitly closed and counter-clockwise seen from outside, as the base interface
requires — the counter-clockwise winding is a5's own (its ring is built
clockwise in the face plane and reversed on the way out) and is verified in this
system's test suite rather than assumed.

# Why the ring is densified

An A5 cell's edges are straight **in the dodecahedral face plane**, and the
equal-area projection does not map those straight lines to great circles. A
five-chord ring would therefore describe a different region from the cell —
1.6% wider in area spread at level 2 than the subdivided ring. a5's `:auto`
subdivision splits each edge into `2^(6-resolution)` pieces down to resolution 6
and one above it, and *that* polyline is the cell: it is what the area, the
caps, the regridding polygon and every containment test see.

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

The cell centre from a5's `cell_to_lonlat`, straight onto the unit sphere.

This is the centroid of the cell's polygon **in the face plane**, pulled back
through the equal-area projection — the point the lattice is built around,
rather than the mean of the boundary vertices — and it is strictly interior to
the cell for every valid index.

One consequence worth naming: the twelve res-0 faces include two whose centre is
exactly a geographic pole, so nothing downstream may build a tangent frame from
a lon/lat east/north pair. The neighbour ordering in `neighbors.jl` measures its
frame from a neighbour's direction instead, which is well-conditioned there.
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

The cell of `grid` containing a point, by a5's own `lonlat_to_cell`.

This replaces the generic tree-descent-plus-point-in-polygon fallback, which for
A5 would be worse than slow: with no [`has_sorted_subtrees`](@ref) the fallback's
tree materialises `1:ncells(grid)` at its root, so the generic `cellat` is
impossible past about level 12 while this one is O(1) at every level.

Never `nothing`: A5 tessellates the whole sphere, so every point is in a cell.

**Ties, and how the search works.** The point is projected onto the nearest
dodecahedron face and the lattice cell is *estimated* in closed form; if that
cell's face-plane polygon contains the point the estimate is the answer. Where
it does not — near a cell corner, or across a face seam — a5 widens over a fixed
24-sample spiral scaled to the cell size, then over the neighbours of the three
closest misses, and finally returns the nearest miss. Every step of that is
deterministic, so a point exactly on a shared boundary resolves the same way on
every call and in every other a5 implementation; this system deliberately does
not impose a tie rule of its own on top, because agreeing with a5 is worth more
than agreeing with a convention no other a5 binding follows.

The `(lon, lat)` method takes **degrees** and is the primitive here — it hands
the coordinates to the a5 arithmetic directly rather than routing them through
`xyz` and back, which is one rounding step shorter.
"""
function cellat(grid::LevelGrid, lon::Real, lat::Real)
    return A5Cell(A5Native.lonlat_to_cell(lon, lat, grid.level))
end

function cellat(grid::LevelGrid, p::GO.UnitSphericalPoint)
    lon = atand(p[2], p[1])
    lat = asind(clamp(p[3], -1.0, 1.0))
    return cellat(grid, lon, lat)
end
