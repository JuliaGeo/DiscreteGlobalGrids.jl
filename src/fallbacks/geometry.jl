# ---------------------------------------------------------------------------
# Geometry generics
#
# Everything a grid can say about a cell's shape, derived from the one required
# primitive `cell_boundary`. A system overrides any of these for speed; none of
# them may be overridden for a different answer.
#
# The boundary contract these are all written against: implicitly closed (the
# first vertex is NOT repeated), counter-clockwise seen from outside the
# sphere, unit-norm points, great-circle arcs between consecutive vertices.
# ---------------------------------------------------------------------------

"""
    closed_ring(points) -> Vector{UnitSphericalPoint{Float64}}

The implicitly closed [`cell_boundary`](@ref) ring made explicitly closed: the
first vertex appended at the end, which is what a GeoInterface `LinearRing`
wants. A ring that already repeats its first vertex is passed through unchanged
— boundaries are implicitly closed by contract, but a system that closes its
own is not thereby made to describe a doubled edge.
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

"""
    open_ring(points) -> (ring, n)

A boundary as a **1-based** vector together with its **open** length: `ring[1:n]`
are the distinct vertices, with the closing edge `ring[n] -> ring[1]` implied.

Two conversions in one, both needed by the spherical ring predicates, which
index `pts[1:n]` positionally and count the closing edge themselves:

  - a boundary that is not 1-based is collected (the contract allows any
    `AbstractVector`);
  - a boundary that repeats its first vertex at the end — which the contract
    says it should not, but which [`closed_ring`](@ref) also tolerates — has
    that vertex dropped from the count. Left in, it becomes a retraced
    zero-length edge, and a retraced edge silently breaks the crossing parity
    those predicates decide containment by.
"""
function open_ring(points)
    ring = firstindex(points) == 1 ? points : collect(points)
    n = length(ring)
    n > 1 && ring[n] == ring[1] && (n -= 1)
    return ring, n
end

"""
    point_in_cell(boundary, p) -> Union{Bool,Nothing}

Whether the unit-sphere point `p` lies in the cell bounded by `boundary`
(boundary points count as inside), or `nothing` when the question is
genuinely undecidable.

Three algorithms, tried in order, because each has a degenerate case the
others do not.

 1. `spherical_ring_encloses` — even-odd parity along the arc from `p` to a
    definitionally exterior anchor, the antipode of the ring's vertex mass.
    Robust for a `p` well away from the cell, and the only one of the three
    that survives a ring which self-intersects *on the sphere*. Consulted
    only when [`anchor_arc_is_conditioned`](@ref) says its test arc is, which
    is the whole of the fix described there.
 2. `spherical_ring_contains` — parity from a wedge at one ring edge, valid
    here because the boundary contract fixes the winding (counter-clockwise
    seen from outside). Well conditioned exactly where (1) is not: its
    bootstrap is a *local* edge midpoint rather than a near-antipode.
 3. the winding number, [`ring_winding_verdict`](@ref) — no bootstrap point at
    all, so no bootstrap can be degenerate. This is what decides the exactly
    symmetric rings on which both (1) and (2) hit exact vertex incidences and
    decline.

`nothing` now means only that `p` is (anti)podal to a boundary vertex, or that
the ring has fewer than three distinct vertices.
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

Whether `spherical_ring_encloses`' test arc `q -> anchor` is long enough for
its own between-ness test to mean anything. **This guard is load-bearing:
without it the generic point-in-cell test reports a cell's own centroid as
outside the cell.**

`spherical_ring_encloses` counts how often the arc from `q` to `anchor`
crosses the ring, and `anchor` is the antipode of the ring's vertex mass. So
for any `q` *inside* a small cell — the centroid above all — `q` is within a
whisker of the vertex mass and the test arc is a **near-half-turn**:
`q · anchor ≈ -1`.

Every "does this point lie between the arc's endpoints" decision in that
parity walk is `point_on_spherical_arc`, whose between-ness test is

    (q · v ≥ q · anchor - tol)  ∧  (anchor · v ≥ q · anchor - tol)

With `q · anchor ≈ -1`, and every dot product on the sphere `≥ -1`, **both
inequalities hold for every point on the sphere**: the arc stops behaving like
an arc and behaves like the entire great circle. That circle meets a cell's
boundary twice, so the parity comes out even and the verdict is a confident
`false` for a point that is not merely inside the cell but is its centroid.

Writing `δ = π - length(arc)`, a point on the *complementary* arc is wrongly
admitted exactly when it lies more than `2δ` behind `q`. So the walk is sound
only while the ring stays nearer to `q` than that, and requiring `δ ≥ R`, for
`R` the ring's angular radius about `q`, is the cheap sufficient condition —
`cos δ ≤ cos R`, i.e.

    -(q · anchor)  ≤  min over vertices v of (q · v)

No tuned constant: the bound is the failure mode written down. A `q` well away
from the cell passes it comfortably and still gets algorithm (1); a `q` inside
the cell fails it and goes straight to the wedge bootstrap, which is exactly
the division of labour [`point_in_cell`](@ref) always claimed to have.

The upstream predicate does carry a guard of its own — it returns `nothing`
when `q · anchor < -1 + 1e-9` — but that is scaled to catch *exact*
antipodality (`δ ≲ 4.5e-5` rad) rather than the near-antipodality that breaks
the between-ness test. That is why a perfectly symmetric cell whose centroid
coincides with its vertex mass to machine precision (IGeo7's twelve pentagons)
used to return `nothing`, while every slightly asymmetric cell — 282 of 3072
HEALPix level-4 diamonds, and cells of H3, S2 and ISEA4R besides — returned a
wrong `false`: one bug, two symptoms either side of one badly-scaled
threshold.
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

Containment by **winding number**: `p` is inside the region the ring bounds
when the ring winds around it a nonzero number of times.

The last resort of [`point_in_cell`](@ref), and the only one of its three
algorithms that picks no bootstrap point — which is the point. Both parity
algorithms decide containment relative to a reference (an exterior anchor, an
edge-midpoint wedge), and an exactly symmetric ring can put that reference in
an exactly degenerate position: a regular pentagon's centre-to-edge-midpoint
arc runs precisely through the opposite vertex, so every wedge
`spherical_ring_contains` tries grazes a vertex and it declines on all of
them. The winding number has no such reference to be unlucky with.

Measured as the net turning of the ring's bearing seen from `p`: each vertex
is projected into the tangent plane at `p`, and the bearing steps — each
folded into `(-π, π]`, which is exact as long as no single edge subtends half
a turn about `p` — sum to `±2π` for a point the ring encircles and `0` for one
it does not.

**Encircling is not containment, because it cannot tell `p` from `-p`.**
Azimuthal winding about a point is the negation of the winding about its
antipode, so `±2π` is reached both inside the cell and at the antipode of the
inside, and a bare magnitude test calls both contained. Hence the guard, which
is the first thing this function does: **every vertex must lie strictly within
a quarter turn of `p`**, which no sub-hemispheric ring can satisfy from the far
side of the sphere, and the verdict abstains otherwise.

Without that guard the correctness of this function would rest on its position
*last* in [`point_in_cell`](@ref)'s cascade — an antipodal probe is in practice
answered by one of the two parity algorithms and never arrives here — which is
an ordering invariant nothing states and nothing enforces. The guard closes the
mode outright instead.

Magnitude, not sign: the boundary contract fixes the winding
counter-clockwise from outside, but a ring that arrives wound the other way
still has an inside, and this reports the region it bounds either way.

Returns `nothing` when any vertex is a quarter turn or more from `p`, and when
a vertex is (anti)podal to `p`, where the bearing is undefined.
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
`(x, y, z)` coordinates.
"""
cell_polygon(grid::AbstractGrid, c::AbstractCellIndex) =
    GI.Polygon([GI.LinearRing(closed_ring(cell_boundary(grid, c)))])

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

The spherical area of the cell in **steradians**, from the exact ring.

Delegates to `GeometryOps.area` on a unit-radius `Spherical` manifold, which
fan-triangulates the ring and sums signed spherical excess — right for cells
that span a large solid angle, where a planar formula is wrong in sign rather
than merely inaccurate.
"""
cell_area(grid::AbstractGrid, c::AbstractCellIndex) =
    Float64(GO.area(GO.Spherical(; radius=1.0), cell_polygon(grid, c)))

# ===========================================================================
# Longitude/latitude extent
# ===========================================================================

"""
    cell_extent(grid, c) -> Extents.Extent{(:X, :Y)}

The lon/lat bounding box of a cell, in degrees, conservative where a rectangle
cannot be exact (see the interface docstring). This is an interoperability
convenience — tree pruning uses [`node_extent`](@ref), a `SphericalCap`.

Three things a naive vertex bounding box gets wrong, and what is done instead:

  - **Arc bulge.** A great-circle edge between two vertices reaches a higher
    latitude than either endpoint. Each edge's extreme-latitude point is
    computed and included when it actually falls on the arc.
  - **The antimeridian.** An edge whose endpoints straddle ±180° makes the
    vertex longitude range meaningless, so `X` becomes `(-180, 180)`.
  - **Poles.** A cell enclosing a pole has no longitude bound at all, and its
    latitude range runs to ±90. Enclosure is decided by the ring's longitude
    winding, and which pole by a spherical containment test.
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

# The latitude extremes of the great-circle arc `a -> b`, or `nothing` when
# neither extreme falls on the arc (the common case, where the endpoints
# already bound it).
#
# The great circle with unit normal `n` attains its extreme latitudes at
# `+-normalize(z - (z . n) n)`; that point is on the minor arc `a -> b` exactly
# when it lies strictly between them, which is the same pair of triple products
# the boundary-arc distance uses.
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

The generic default: the cell's own bounding cap inflated by
[`cap_inflation(sys)`](@ref cap_inflation), which is O(1) at every depth and
sound for every system whose measured child overhang stays under that factor.

See the covering law in the interface docstring — this default is the reason
`cap_inflation` exists, and a system that overrides `node_extent` outright
(HEALPix's exact subtree cap) ignores it.
"""
node_extent(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex) =
    inflate_cap(cell_cap(levelgrid(sys, level(c)), c), cap_inflation(sys))

# ===========================================================================
# Manifold
# ===========================================================================

"""
    GeometryOpsCore.best_manifold(grid::AbstractGrid) -> GO.Spherical

The compute manifold of every grid in this package: the **unit** sphere, which
is the frame all its geometry is expressed in.

`radius = 1` is not a placeholder. Cell coordinates are unit-sphere `(x, y, z)`,
so an area computed against this manifold comes out in steradians — see
[`cell_area`](@ref) — and scaling to a physical sphere is the caller's one
multiplication by `R^2`. Handing `ConservativeRegridding` a manifold whose
radius disagreed with the coordinates would put that factor in twice.
"""
GOCore.best_manifold(grid::AbstractGrid) = GO.Spherical(; radius=1.0)
