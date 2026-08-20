# Area-only spherical intersection

"""
    IntersectionAreaOperator(manifold)

Pair-area operator over `GO.intersection_area`: the clipped ring is integrated
where it is clipped, without materializing the result polygon or the area
pass's point collection. Equal bit for bit to
`GO.area(manifold, GO.intersection(alg, subject, clip))`.

Carries the clipping cache that keeps the call allocation-free.
[`task_local_operator`](@ref ConservativeRegridding.task_local_operator) gives
each assembly task an operator with a cache of its own.
"""
struct IntersectionAreaOperator{M<:GO.Spherical,C}
    manifold::M
    cache::C
end

IntersectionAreaOperator(m::GO.Spherical) =
    IntersectionAreaOperator(m, GO.SutherlandHodgmanCache(m))

ConservativeRegridding.task_local_operator(op::IntersectionAreaOperator) =
    IntersectionAreaOperator(op.manifold)

# Spherical clips take the area-only kernel; other manifolds keep the default.
_intersectionoperator(m::GO.Spherical) = IntersectionAreaOperator(m)
_intersectionoperator(m::GOCore.Manifold) =
    ConservativeRegridding.DefaultIntersectionOperator(m)

(op::IntersectionAreaOperator)(p1, p2) = GO.intersection_area(
    GO.ConvexConvexSutherlandHodgman(op.manifold), p1, p2; cache = op.cache)
