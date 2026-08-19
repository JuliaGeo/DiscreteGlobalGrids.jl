# Area-only spherical intersection

"""
    SphericalClipAreaOperator(manifold, cache)

Compute `GO.area(manifold, GO.intersection(SutherlandHodgman, subject, clip))`
with the clipped ring integrated inside the task-local cache: the same
clipping and the same triangle fan, without materializing the result polygon
or the area pass's point collection. [`task_local_operator`](@ref
ConservativeRegridding.task_local_operator) hands each assembly task a private
cache.
"""
struct SphericalClipAreaOperator{M<:GO.Spherical,C}
    manifold::M
    cache::C
end

SphericalClipAreaOperator(m::GO.Spherical) = SphericalClipAreaOperator(m, nothing)

ConservativeRegridding.task_local_operator(op::SphericalClipAreaOperator) =
    SphericalClipAreaOperator(op.manifold, GO.SutherlandHodgmanCache(op.manifold))

# Spherical clips take the streaming kernel; other manifolds keep the default.
_intersectionoperator(m::GO.Spherical) = SphericalClipAreaOperator(m)
_intersectionoperator(m::GOCore.Manifold) =
    ConservativeRegridding.DefaultIntersectionOperator(m)

(op::SphericalClipAreaOperator)(p1, p2) =
    _clip_area(op.manifold, something(op.cache, GO.SutherlandHodgmanCache(op.manifold)), p1, p2)

# Mirrors GeometryOps' spherical `_intersection_sutherland_hodgman` bit for bit,
# ending in the fan integral instead of the copied-out result ring.
function _clip_area(m::GO.Spherical, cache::GO.SutherlandHodgmanCache, poly_a, poly_b)
    T = Float64
    ring_a = GI.getexterior(poly_a)
    ring_b = GI.getexterior(poly_b)

    first_pt = GI.getpoint(ring_a, 1)
    first_pt isa US.UnitSphericalPoint || throw(ArgumentError(
        "spherical clipping requires UnitSphericalPoint coordinates, got $(typeof(first_pt))"))

    # Clip ring vertices, deduplicating the closing point.
    clip_points = cache.clip
    empty!(clip_points)
    for point in GI.getpoint(ring_b)
        if !isempty(clip_points) && point ≈ clip_points[1]
            continue
        end
        push!(clip_points, US.UnitSphericalPoint{T}(point))
    end

    # Subject ring, ping-ponging between the cache buffers as it is clipped.
    buf_in, buf_out = cache.input, cache.output
    empty!(buf_in)
    for point in GI.getpoint(ring_a)
        if !isempty(buf_in) && point == buf_in[1]
            continue
        end
        push!(buf_in, US.UnitSphericalPoint{T}(point))
    end

    original_subject = cache.subject
    empty!(original_subject)
    append!(original_subject, buf_in)

    n_clip = length(clip_points)
    for i in 1:n_clip
        isempty(buf_in) && break
        edge_start = clip_points[i]
        edge_end = clip_points[mod1(i + 1, n_clip)]
        edge_start == edge_end && continue
        GO._sh_clip_to_edge_spherical!(buf_out, buf_in, edge_start, edge_end, T)
        buf_in, buf_out = buf_out, buf_in
    end

    if isempty(buf_in)
        # The clip ring may lie entirely inside the subject; its own fan is the area.
        if !isempty(clip_points) &&
           all(p -> GO._point_in_convex_spherical_polygon(p, original_subject), clip_points)
            return _openring_area(m, clip_points)
        end
        return 0.0
    end
    # 1-2 surviving points cannot form a ring; the degenerate polygon's area is zero.
    length(buf_in) < 3 && return 0.0
    return _openring_area(m, buf_in)
end

# `GO.area` of the closed result ring, read off the open buffer: the same
# Eriksson triangle fan in the same order, `abs`, times the radius squared.
function _openring_area(m::GO.Spherical, buf::Vector{US.UnitSphericalPoint{Float64}})
    n = length(buf)
    n < 3 && return 0.0
    p1 = @inbounds buf[1]
    area = 0.0
    for i in 2:(n - 1)
        area += GO._spherical_triangle_area(GO.Eriksson(), p1,
            @inbounds(buf[i]), @inbounds(buf[i + 1]))
    end
    return abs(area) * m.radius^2
end
