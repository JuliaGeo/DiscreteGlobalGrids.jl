# Spherical polygon area and vector first moment, and the pair operator that
# measures both for one overlap; `conservative_second_order.jl` reads them.

"""
    PolygonMoments(area, moment)

A spherical polygon's area and vector first moment, the two integrals a
second-order conservative weight is built from:

    area   = ∫_P dA
    moment = ∫_P x dA

`x` is the polygon's position on the *unit* sphere, so `moment ./ area` is a
non-unit mean position vector — the direction of the polygon's centroid, with a
length below one that measures how spread out the polygon is. Both integrals
carry the manifold's `radius^2`, exactly as `GO.intersection_area` does, so the
scale cancels in that ratio and `area` alone answers in the manifold's units.

Addition sums both fields, which is what makes this type an `eltype` a sparse
matrix can be assembled with: `SparseArrays.sparse(rows, cols, vals, m, n)`
combines duplicate coordinates with `+`, so a destination-source pair split
across several overlap polygons accumulates into one entry whose `area` and
`moment` are the union's.
"""
struct PolygonMoments
    area::Float64
    moment::NTuple{3,Float64}
end

Base.zero(::Type{PolygonMoments}) = PolygonMoments(0.0, (0.0, 0.0, 0.0))
Base.zero(m::PolygonMoments) = zero(typeof(m))

Base.:+(a::PolygonMoments, b::PolygonMoments) =
    PolygonMoments(a.area + b.area,
        (a.moment[1] + b.moment[1], a.moment[2] + b.moment[2],
            a.moment[3] + b.moment[3]))

# Both fields are integrals over the same region, so scaling the region's
# measure scales both. `sum` over a sparse column of moments needs this.
Base.:*(m::PolygonMoments, c::Real) =
    PolygonMoments(m.area * c,
        (m.moment[1] * c, m.moment[2] * c, m.moment[3] * c))
Base.:*(c::Real, m::PolygonMoments) = m * c

# An empty overlap: what a sparse assembly may drop, and what `zero` returns.
Base.iszero(m::PolygonMoments) = iszero(m.area) && all(iszero, m.moment)

Base.show(io::IO, m::PolygonMoments) =
    print(io, "PolygonMoments(area = ", m.area, ", moment = ", m.moment, ")")

# Ring moments

"""
    polygonmoments(m::GO.Spherical, pts; closed::Bool) -> PolygonMoments

Measure the ring `pts` — a vector of points on the unit sphere, converted the
way `GO._ring_area` converts them — returning its area and vector first moment.

`closed` says whether the ring repeats its first vertex last, and is handled
exactly as `GO._ring_area` handles it: a closed ring's repeated last vertex is
dropped, an open one's is kept and the closing edge supplied. Passing `true`
for an open ring would swallow the last vertex of any sliver whose ends fall
within `isapprox`'s default tolerance of each other.

The area comes from `GO._ring_area`, absolute value and radius scaling applied
in the same order as in `GO.intersection_area`, so it is bit for bit the area
[`IntersectionAreaOperator`](@ref) and `Conservative()` measure for the same
ring. That equality is the point: the second-order weights must not disagree
with the coverage areas they correct.

The moment uses the closed form for a great-circle ring. The vector area of a
surface is `½∮ x × dx`, and on the unit sphere the outward normal is `x`
itself, so `∫_P x dA = ½ Σ_edges θ_e (a_e × b_e) / |a_e × b_e|` with `θ_e` the
arc length of the edge from `a_e` to `b_e`. A degenerate edge contributes
nothing and is skipped.

That sum is evaluated shifted to the first vertex `r`, as

    ½ Σ (a_e − r) × (b_e − r) + ½ Σ (θ_e / sin θ_e − 1) (a_e × b_e)

the same number in exact arithmetic — around a closed ring the terms in `r`
cancel — and a different one in floating point:

  - **Unshifted**, each `θ_e n̂_e` is of order `θ` while their sum is of order
    `θ²`, and the cross product of two nearly parallel unit vectors carries an
    absolute error near machine epsilon whatever its length.
  - **Shifted to `r`**, every term is of order `θ²` with relative error near
    epsilon, and the spherical correction is `O(θ²)` smaller again, so cells of
    any size measure to full precision.

For a 30 m Copernicus pixel, `θ ≈ 5e-6`: the error of the plain sum is the size
of the cell, and the gradient fitted from such positions is noise.

The formula is signed by the ring's orientation. The area is not, so a
clockwise ring's moment is negated to match: both fields then describe the same
enclosed region whichever way its vertices were wound.
"""
function polygonmoments(m::GO.Spherical, pts::AbstractVector; closed::Bool)
    signedarea = GO._ring_area(m, pts, Float64; closed)
    scale = GO._area_scale(m)
    n = _ringvertexcount(pts, closed)
    n < 3 && return zero(PolygonMoments)
    mx = my = mz = 0.0
    firstpoint = US.UnitSphericalPoint(GI.PointTrait(), pts[1])
    r = firstpoint
    a = firstpoint
    for i in 1:n
        # The last edge closes the ring, whether or not `pts` spelled it out.
        b = i == n ? firstpoint : US.UnitSphericalPoint(GI.PointTrait(), pts[i + 1])
        cx = a[2] * b[3] - a[3] * b[2]
        cy = a[3] * b[1] - a[1] * b[3]
        cz = a[1] * b[2] - a[2] * b[1]
        s = sqrt(cx * cx + cy * cy + cz * cz)
        # `s == 0` is a repeated or antipodal vertex: no arc, no contribution.
        if s > 0
            # The chord term, from the vertices shifted to `r`; the first edge
            # starts at `r` and contributes nothing here.
            if i > 1
                ax, ay, az = a[1] - r[1], a[2] - r[2], a[3] - r[3]
                bx, by, bz = b[1] - r[1], b[2] - r[2], b[3] - r[3]
                mx += ay * bz - az * by
                my += az * bx - ax * bz
                mz += ax * by - ay * bx
            end
            # The spherical correction: `atan(s, a⋅b)` is the arc length, and
            # `θ / sin θ - 1` is the excess of the arc over its chord, `O(θ²)`.
            g = atan(s, a[1] * b[1] + a[2] * b[2] + a[3] * b[3]) / s - 1
            mx += cx * g
            my += cy * g
            mz += cz * g
        end
        a = b
    end
    orientation = signbit(signedarea) ? -0.5 : 0.5
    w = orientation * scale
    return PolygonMoments(abs(signedarea) * scale, (mx * w, my * w, mz * w))
end

# The vertex count `GO._ring_area` fans over: a ring the caller called closed
# drops its repeated first vertex, by the same `isapprox` test. Kept in step
# with that function so the moment and the area sum over the same edges.
function _ringvertexcount(pts::AbstractVector, closed::Bool)
    n = length(pts)
    n < 3 && return 0
    if closed && US.UnitSphericalPoint(GI.PointTrait(), pts[n]) ≈
                 US.UnitSphericalPoint(GI.PointTrait(), pts[1])
        n -= 1
    end
    return n < 3 ? 0 : n
end

"""
    cellmoments(space::RegridSpace, i::Int) -> PolygonMoments

Measure the space's own cell `i`: its area and vector first moment, over the
closed exterior ring `getcell` returns, on the space's [`manifold`](@ref).

The area is bit for bit `GO.area(manifold(space), getcell(space, i))`, which is
what makes a cell's moment comparable with the overlap moments that tile it.
"""
cellmoments(space::RegridSpace, i::Int) = polygonmoments(
    manifold(space), _exteriorpoints(getcell(space, i)); closed = true)

_exteriorpoints(geom) = collect(GI.getpoint(GI.getexterior(geom)))

# Pair operator

# The one call into a GeometryOps internal. `GO.intersection_area` is itself a
# five-line wrapper over this same seam — clip into the cache's buffer, then
# `GO._ring_area` that buffer — so measuring the buffer here is exactly what
# keeps the two operators' areas identical. A public GO clip-to-buffer API
# would replace this one line.
_clipring!(cache, subject, clip) =
    GO._sh_clip_spherical!(cache, subject, clip, Float64)

"""
    IntersectionMomentOperator(manifold)

Pair-moment operator: the moment counterpart of
[`IntersectionAreaOperator`](@ref), returning a [`PolygonMoments`](@ref) for
the overlap of a source cell and a destination cell rather than an area alone.
The clipped ring is integrated where it is clipped, without materializing the
result polygon.

Both inputs must be convex and wound counter-clockwise seen from outside the
sphere, which is what `GO.ConvexConvexSutherlandHodgman` requires; the returned
buffer is open, so the ring is measured with `closed = false`.

`.area` is bit for bit what [`IntersectionAreaOperator`](@ref) returns for the
same pair, so second-order weights and conservative coverage are measured
against the same overlap.

Carries the clipping cache that keeps the call allocation-free.
[`task_local_operator`](@ref ConservativeRegridding.task_local_operator) gives
each assembly task an operator with a cache of its own.
"""
struct IntersectionMomentOperator{M<:GO.Spherical,C}
    manifold::M
    cache::C
end

IntersectionMomentOperator(m::GO.Spherical) =
    IntersectionMomentOperator(m, GO.SutherlandHodgmanCache(m))

ConservativeRegridding.task_local_operator(op::IntersectionMomentOperator) =
    IntersectionMomentOperator(op.manifold)

# The two assembly hooks that stop defaulting once the operator stops returning
# a number: `output_eltype` defaults to `Float64`, and `should_store_result`
# errors outright on anything but a number rather than guess. An empty overlap
# is dropped on its area, exactly as the area-only path drops a zero area.
ConservativeRegridding.output_eltype(::IntersectionMomentOperator) = PolygonMoments

ConservativeRegridding.should_store_result(
    ::IntersectionMomentOperator, pm::PolygonMoments) = pm.area > 0

(op::IntersectionMomentOperator)(p1, p2) = polygonmoments(
    op.manifold, _clipring!(op.cache, p1, p2); closed = false)
