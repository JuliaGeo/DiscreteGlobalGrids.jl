# Area-only spherical intersection

"""
    IntersectionAreaOperator(manifold, cache)

Pair-area operator over `GO.intersection_area`: the clipped ring is integrated
where it is clipped, without materializing the result polygon or the area
pass's point collection. Equal bit for bit to
`GO.area(manifold, GO.intersection(alg, subject, clip))`.
[`task_local_operator`](@ref ConservativeRegridding.task_local_operator) hands
each assembly task a private cache, which keeps the call allocation-free.
"""
struct IntersectionAreaOperator{M<:GO.Spherical,C}
    manifold::M
    cache::C
end

IntersectionAreaOperator(m::GO.Spherical) = IntersectionAreaOperator(m, nothing)

ConservativeRegridding.task_local_operator(op::IntersectionAreaOperator) =
    IntersectionAreaOperator(op.manifold, GO.SutherlandHodgmanCache(op.manifold))

# Spherical clips take the area-only kernel; other manifolds keep the default.
_intersectionoperator(m::GO.Spherical) = IntersectionAreaOperator(m)
_intersectionoperator(m::GOCore.Manifold) =
    ConservativeRegridding.DefaultIntersectionOperator(m)

(op::IntersectionAreaOperator)(p1, p2) = _clip_area(op.manifold, op.cache, p1, p2)

_clip_area(m, cache, p1, p2) =
    GO.intersection_area(GO.ConvexConvexSutherlandHodgman(m), p1, p2; cache)

# No task-local cache yet: let GeometryOps allocate its own.
_clip_area(m, ::Nothing, p1, p2) =
    GO.intersection_area(GO.ConvexConvexSutherlandHodgman(m), p1, p2)
