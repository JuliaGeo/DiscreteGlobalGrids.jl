# Geometry fallbacks derived from `cell_boundary`. Boundaries are implicitly
# closed, counter-clockwise unit-sphere rings with great-circle edges. A system
# may override any of these for speed, none of them for a different answer.

"""
    closed_ring(points) -> AbstractVector{UnitSphericalPoint{Float64}}

Return a 1-based, explicitly closed copy of a boundary. An already repeated
closing vertex is not duplicated. A `Helpers.SmallList` boundary closes into a
`SmallList` one longer, so inline storage survives; anything else copies into a
`Vector`.
"""
function closed_ring(points)
    n = length(points)
    n == 0 && return USPoint[]
    out = Vector{USPoint}(undef, n + 1)
    # `enumerate` rather than `points[i]`: the boundary contract allows any
    # `AbstractVector`, including one that is not 1-based.
    for (i, p) in enumerate(points)
        @inbounds out[i] = USPoint(p[1], p[2], p[3])
    end
    @inbounds out[n+1] = out[1]
    # Already closed: drop the duplicate we just introduced.
    n > 1 && out[n] == out[1] && return out[1:n]
    return out
end

# Inline storage in, inline storage out: a system that answers its boundary in a
# `SmallList` must not have it spilled onto the heap just to close it.
function closed_ring(points::Helpers.SmallList{N}) where {N}
    out = Helpers.empty_small_list(Val(N + 1), USPoint(1.0, 0.0, 0.0))
    for p in points
        out = Helpers.small_push(out, USPoint(p[1], p[2], p[3]))
    end
    n = length(out)
    n == 0 && return out
    # Already closed: nothing to append.
    n > 1 && @inbounds(out[n] == out[1]) && return out
    return Helpers.small_push(out, @inbounds out[1])
end

"""
    open_ring(points) -> (ring, n)

Return a 1-based boundary and its open length `n`, with distinct vertices in
`ring[1:n]` and the closing edge implied. A repeated closing vertex is excluded
from `n` to avoid a degenerate edge in parity predicates.
"""
function open_ring(points)
    ring = firstindex(points) == 1 ? points : collect(points)
    n = length(ring)
    n > 1 && ring[n] == ring[1] && (n -= 1)
    return ring, n
end

"""
    point_in_cell(boundary, p) -> Union{Bool,Nothing}

Test whether unit-sphere point `p` is inside or on `boundary`. It tries a
conditioned exterior-anchor parity test, a local-wedge parity test, then winding
number. Return `nothing` for fewer than three distinct vertices or an undefined
antipodal bearing.
"""
function point_in_cell(boundary, p)
    ring, n = open_ring(boundary)
    n >= 3 || return nothing
    # Computed here rather than left to the `encloses` default so that the
    # conditioning test and the algorithm see the same anchor, and so the
    # O(n) vertex-mass sum is walked once rather than twice.
    anchor = US.spherical_exterior_anchor(ring, n)
    if anchor_arc_is_conditioned(ring, n, p, anchor)
        verdict = US.spherical_ring_encloses(ring, n, p; anchor)
        verdict === nothing || return verdict
    end
    verdict = US.spherical_ring_contains(ring, n, p)
    verdict === nothing || return verdict
    return ring_winding_verdict(ring, n, p)
end

"""
    anchor_arc_is_conditioned(ring, n, q, anchor) -> Bool

Return whether the exterior-anchor parity arc is sufficiently separated from
antipodality. For ring radius `R` about `q` and complementary-arc length `δ`,
the sufficient condition `δ ≥ R` is

    -(q · anchor)  ≤  min over vertices v of (q · v)

Near-antipodal arcs fail the underlying between-ness test even before exact
antipodality, so this geometric bound is required in addition to its tolerance.
A `q` inside the cell — its centroid above all — fails the bound and falls
through to the wedge test; without it such points are reported outside the cell.
"""
function anchor_arc_is_conditioned(ring, n, q, anchor)
    anchor === nothing && return false
    cos_slack = -(q[1] * anchor[1] + q[2] * anchor[2] + q[3] * anchor[3])
    cos_radius = 1.0
    for i in 1:n
        v = ring[i]
        d = q[1] * v[1] + q[2] * v[2] + q[3] * v[3]
        d < cos_radius && (cos_radius = d)
    end
    return cos_slack <= cos_radius
end

"""
    ring_winding_verdict(ring, n, p) -> Union{Bool,Nothing}

Test containment by the absolute winding number around `p`. Tangent-plane
bearing steps folded into `(-π, π]` sum to `±2π` for a point the ring encircles
and `0` for one it does not. The last resort of [`point_in_cell`](@ref), and
the only one of its algorithms that picks no bootstrap point to be degenerate.

Encircling is not containment: winding about `p` is the negation of the winding
about `-p`, so the antipode of an interior point also sums to `±2π`. A
sub-hemispheric guard therefore runs first — every vertex must lie strictly
within a quarter turn of `p`, else the verdict abstains. Without it correctness
would rest on this running last in the cascade, an ordering invariant nothing
enforces.

Magnitude, not sign, so a ring wound either way reports the region it bounds.

Return `nothing` when any vertex is a quarter turn or more from `p`, or when a
vertex is antipodal to `p` and its bearing is undefined.
"""
function ring_winding_verdict(ring, n, p)
    # The sub-hemispheric guard. See the docstring: without it this function
    # answers `true` for the antipode of an interior point just as readily as
    # for the point itself.
    for i in 1:n
        v = ring[i]
        p[1] * v[1] + p[2] * v[2] + p[3] * v[3] > 0.0 || return nothing
    end

    # A right-handed tangent frame at `p` seen from outside the sphere; which
    # reference direction seeds it does not matter, only the handedness.
    a = abs(p[3]) < 0.9 ? (0.0, 0.0, 1.0) : (1.0, 0.0, 0.0)
    radial = a[1] * p[1] + a[2] * p[2] + a[3] * p[3]
    t = (a[1] - radial * p[1], a[2] - radial * p[2], a[3] - radial * p[3])
    tn = sqrt(t[1]^2 + t[2]^2 + t[3]^2)
    tn < 1e-12 && return nothing
    e1 = (t[1] / tn, t[2] / tn, t[3] / tn)
    e2 = (p[2] * e1[3] - p[3] * e1[2],
          p[3] * e1[1] - p[1] * e1[3],
          p[1] * e1[2] - p[2] * e1[1])

    first_bearing = _ring_bearing(p, e1, e2, ring[1])
    first_bearing === nothing && return nothing
    total = 0.0
    previous = first_bearing
    for i in 2:n
        b = _ring_bearing(p, e1, e2, ring[i])
        b === nothing && return nothing
        total += _fold_turn(b - previous)
        previous = b
    end
    # The closing edge `ring[n] -> ring[1]`, which the loop never took.
    total += _fold_turn(first_bearing - previous)
    return abs(total) > π
end

# The bearing of `v` seen from `p` in the tangent frame `(e1, e2)`, or
# `nothing` where `v` is (anti)podal to `p` and no bearing exists.
function _ring_bearing(p, e1, e2, v)
    d = (v[1] - p[1], v[2] - p[2], v[3] - p[3])
    x = d[1] * e1[1] + d[2] * e1[2] + d[3] * e1[3]
    y = d[1] * e2[1] + d[2] * e2[2] + d[3] * e2[3]
    return (x == 0.0 && y == 0.0) ? nothing : atan(y, x)
end

# One bearing step folded into (-pi, pi].
function _fold_turn(d)
    while d <= -π
        d += 2π
    end
    while d > π
        d -= 2π
    end
    return d
end

"""
    cell_polygon(grid, c) -> GI.Polygon

See the interface docstring. One `LinearRing`, explicitly closed, unit-sphere
`(x, y, z)` coordinates. The ring carries whatever [`closed_ring`](@ref) gives
it, so a system with an inline boundary gets an `isbits` polygon.
"""
cell_polygon(grid::AbstractGrid, c::AbstractCellIndex) =
    GI.Polygon(_ring_list(GI.LinearRing(closed_ring(cell_boundary(grid, c)))))

# `[ring]` would wrap a ring that is itself inline in a one-element heap vector.
_ring_list(ring::T) where {T} = Helpers.SmallList{1,T}(1, (ring,))

"""
    getcell(grid, i) -> GI.Polygon

`ConservativeRegridding.Trees.getcell`: position -> unit-sphere polygon.
Implemented once, here, as `cell_polygon(grid, cellindex(grid, i))`; grid
authors never write it.
"""
getcell(grid::AbstractGrid, i::Int) = cell_polygon(grid, cellindex(grid, i))

# `Trees` addresses a grid as a tree source too, and its lazy all-cells form is
# part of the surface every grid gets.
getcell(grid::AbstractGrid) = (getcell(grid, i) for i in 1:ncells(grid))

"""
    cell_area(grid, c) -> Float64

Return the cell's spherical area in steradians, computed on the exact ring with
a unit-radius spherical manifold.
"""
cell_area(grid::AbstractGrid, c::AbstractCellIndex) =
    Float64(GO.area(GO.Spherical(; radius=1.0), cell_polygon(grid, c)))

# ===========================================================================
# Longitude/latitude extent
# ===========================================================================

"""
    cell_extent(grid, c) -> Extents.Extent{(:X, :Y)}

Return a conservative lon/lat bounding box in degrees. It includes
great-circle latitude extrema and expands longitude to `(-180, 180)` for
antimeridian crossings or enclosed poles.
"""
function cell_extent(grid::AbstractGrid, c::AbstractCellIndex)
    points, n = open_ring(cell_boundary(grid, c))
    n == 0 && throw(ArgumentError("cell $c has an empty boundary"))

    latmin = latmax = 0.0
    lonmin = 180.0
    lonmax = -180.0
    winding = 0.0
    crosses = false
    first_lon = 0.0
    previous_lon = 0.0
    for i in 1:n
        p = points[i]
        lon, lat = lonlat(p)
        if i == 1
            latmin = latmax = lat
            first_lon = lon
        else
            latmin = min(latmin, lat)
            latmax = max(latmax, lat)
            delta = lon - previous_lon
            if delta > 180.0
                delta -= 360.0
                crosses = true
            elseif delta < -180.0
                delta += 360.0
                crosses = true
            end
            winding += delta
        end
        previous_lon = lon
        lonmin = min(lonmin, lon)
        lonmax = max(lonmax, lon)
        # Arc bulge: the closing edge is included by the wrap-around index.
        q = points[i == n ? 1 : i + 1]
        elat = _arc_lat_extremes(p, q)
        if elat !== nothing
            latmin = min(latmin, elat[1])
            latmax = max(latmax, elat[2])
        end
    end
    # The closing edge's longitude step, which the loop above never took.
    delta = first_lon - previous_lon
    if delta > 180.0
        delta -= 360.0
        crosses = true
    elseif delta < -180.0
        delta += 360.0
        crosses = true
    end
    winding += delta

    if abs(winding) > 180.0
        # The ring winds around a pole: no longitude bound exists, and the
        # enclosed pole is at the latitude limit. An undecidable ring gives up
        # on both poles rather than guessing one.
        north = point_in_cell(points, USPoint(0.0, 0.0, 1.0))
        south = point_in_cell(points, USPoint(0.0, 0.0, -1.0))
        ymin = south === false ? latmin : -90.0
        ymax = north === false ? latmax : 90.0
        return Extents.Extent(X=(-180.0, 180.0), Y=(ymin, ymax))
    end
    crosses && return Extents.Extent(X=(-180.0, 180.0), Y=(latmin, latmax))
    return Extents.Extent(X=(lonmin, lonmax), Y=(latmin, latmax))
end

# Great-circle latitude extrema on the minor arc `a -> b`, if present. For unit
# normal `n`, candidates are `±normalize(z - (z ⋅ n)n)`.
function _arc_lat_extremes(a, b)
    nx = a[2] * b[3] - a[3] * b[2]
    ny = a[3] * b[1] - a[1] * b[3]
    nz = a[1] * b[2] - a[2] * b[1]
    nn = nx * nx + ny * ny + nz * nz
    nn <= 1e-24 && return nothing            # degenerate or antipodal edge
    inv = 1.0 / sqrt(nn)
    ux, uy, uz = nx * inv, ny * inv, nz * inv
    # z - (z . u) u, the tangent direction of steepest latitude gain.
    vx, vy, vz = -uz * ux, -uz * uy, 1.0 - uz * uz
    vn = sqrt(vx * vx + vy * vy + vz * vz)
    # `u` parallel to the pole means the arc's great circle IS the equator,
    # where latitude is constant and there is no extreme to find.
    vn <= 1e-12 && return nothing
    px, py, pz = vx / vn, vy / vn, vz / vn
    lo = nothing
    hi = nothing
    for s in (1.0, -1.0)
        qx, qy, qz = s * px, s * py, s * pz
        # `a -> q -> b` in the orientation of `n`.
        ((a[2] * qz - a[3] * qy) * nx + (a[3] * qx - a[1] * qz) * ny +
         (a[1] * qy - a[2] * qx) * nz) > 0 || continue
        ((qy * b[3] - qz * b[2]) * nx + (qz * b[1] - qx * b[3]) * ny +
         (qx * b[2] - qy * b[1]) * nz) > 0 || continue
        lat = asind(clamp(qz, -1.0, 1.0))
        lo = lo === nothing ? lat : min(lo, lat)
        hi = hi === nothing ? lat : max(hi, lat)
    end
    lo === nothing && return nothing
    return (lo, hi)
end

# ===========================================================================
# The default node extent
# ===========================================================================

"""
    node_extent(sys, c) -> SphericalCap

Return the cell cap inflated by [`cap_inflation(sys)`](@ref cap_inflation).
This `O(1)` fallback is sound when descendant overhang stays within that factor,
where overhang is measured on descendant *boundary geometry* against the cell's
cap. Descendant caps are separate bounds and are not covered by the result — see
the covering law in [`node_extent`](@ref).
"""
node_extent(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex) =
    inflate_cap(cell_cap(levelgrid(sys, level(c)), c), cap_inflation(sys))

# ===========================================================================
# Manifold
# ===========================================================================

"""
    GeometryOpsCore.best_manifold(grid::AbstractGrid) -> GO.Spherical

Return the unit-sphere compute manifold. Areas are therefore steradians; callers
apply physical scaling once with `R²`.
"""
GOCore.best_manifold(grid::AbstractGrid) = GO.Spherical(; radius=1.0)
