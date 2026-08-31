# Spherical polygon area and first moment, and the pair operator measuring both
# for one overlap. Read by `conservative_second_order.jl`.

"""
    PolygonMoments(area, moment)

A spherical polygon's area `∫_P dA` and first moment `∫_P x dA`, with `x` on the
*unit* sphere. `moment ./ area` gives the centroid's direction, its length below
one measuring the polygon's spread; both integrals carry the manifold's
`radius^2`, so that ratio is scale-free and `area` answers in manifold units.

Addition sums both fields, so `sparse` folds a destination-source pair split
across several overlaps into one entry.
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

# Both fields integrate the same region, so its measure scales both.
Base.:*(m::PolygonMoments, c::Real) =
    PolygonMoments(m.area * c,
        (m.moment[1] * c, m.moment[2] * c, m.moment[3] * c))
Base.:*(c::Real, m::PolygonMoments) = m * c

Base.iszero(m::PolygonMoments) = iszero(m.area) && all(iszero, m.moment)

Base.show(io::IO, m::PolygonMoments) =
    print(io, "PolygonMoments(area = ", m.area, ", moment = ", m.moment, ")")

# Ring moments

"""
    polygonmoments(m::GO.Spherical, pts; closed::Bool) -> PolygonMoments

Area and first moment of the ring `pts`, converted and closed exactly as
`GO._ring_area` does. Passing `closed = true` for an open ring swallows the last
vertex of any sliver whose ends fall within `isapprox`'s tolerance.

The area is `GO._ring_area`'s, scaled in `GO.intersection_area`'s order, so it
matches [`IntersectionAreaOperator`](@ref) bit for bit: second-order weights
must agree with the coverage areas they correct.

The moment is the great-circle form `∫_P x dA = ½ Σ θ_e n̂_e`, summed shifted to
the first vertex `r`:

    ½ Σ (a_e − r) × (b_e − r) + ½ Σ (θ_e / sin θ_e − 1) (a_e × b_e)

Unshifted, terms of order `θ` cancel to `θ²` while their cross products carry
absolute error near epsilon. Shifted, every term is `O(θ²)` with relative error
near epsilon, so a 30 m pixel (`θ ≈ 5e-6`) keeps full precision where the plain
sum loses the whole cell.

Orientation signs the moment but not the area, so a clockwise ring's moment is
negated to match.
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

# The vertex count `GO._ring_area` fans over, by the same `isapprox` test, so
# the moment and the area sum over the same edges.
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

Area and first moment of cell `i`, over the closed exterior ring `getcell`
returns. The area is bit for bit `GO.area(manifold(space), getcell(space, i))`,
which makes a cell's moment comparable with the overlap moments tiling it.
"""
cellmoments(space::RegridSpace, i::Int) = polygonmoments(
    manifold(space), _exteriorpoints(getcell(space, i)); closed = true)

_exteriorpoints(geom) = collect(GI.getpoint(GI.getexterior(geom)))

# Pair operator

# The one call into a GeometryOps internal. `GO.intersection_area` wraps this
# same seam, so measuring its buffer here keeps both operators' areas identical.
_clipring!(cache, subject, clip) =
    GO._sh_clip_spherical!(cache, subject, clip, Float64)

"""
    IntersectionMomentOperator(manifold)

The moment counterpart of [`IntersectionAreaOperator`](@ref), returning a
[`PolygonMoments`](@ref) for one source-destination overlap. The clipped ring is
integrated in place, and `.area` is bit for bit the area operator's, so weights
and coverage measure the same overlap.

`GO.ConvexConvexSutherlandHodgman` requires both inputs convex and wound
counter-clockwise from outside the sphere, and returns an open buffer.

The clipping cache keeps the call allocation-free;
[`task_local_operator`](@ref ConservativeRegridding.task_local_operator) gives
each assembly task its own.
"""
struct IntersectionMomentOperator{M<:GO.Spherical,C}
    manifold::M
    cache::C
end

IntersectionMomentOperator(m::GO.Spherical) =
    IntersectionMomentOperator(m, GO.SutherlandHodgmanCache(m))

ConservativeRegridding.task_local_operator(op::IntersectionMomentOperator) =
    IntersectionMomentOperator(op.manifold)

# Both hooks default to numbers only: `output_eltype` to `Float64`, and
# `should_store_result` errors on anything else. Empty overlaps drop on area.
ConservativeRegridding.output_eltype(::IntersectionMomentOperator) = PolygonMoments

ConservativeRegridding.should_store_result(
    ::IntersectionMomentOperator, pm::PolygonMoments) = pm.area > 0

(op::IntersectionMomentOperator)(p1, p2) = polygonmoments(
    op.manifold, _clipring!(op.cache, p1, p2); closed = false)
