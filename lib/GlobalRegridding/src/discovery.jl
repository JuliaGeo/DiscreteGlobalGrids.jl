# Discover connected source and destination chunks: the source grid's chunk
# extents (`grid_extent` of each chunk, read lazily) packed in an R-tree and
# tested against each destination chunk's extent dilated by the method's
# support radius.

"""
    DilatedIntersects(radius)

Test whether two caps are within `radius` radians. Dilation applies at every
tree level, so pruning remains valid for methods with nonzero support.

# Example

```julia
cap_a = SphericalCap(USPoint(1.0, 0.0, 0.0), 0.1)
cap_b = SphericalCap(USPoint(cos(0.2), sin(0.2), 0.0), 0.1)
intersects = DilatedIntersects(0.05)
intersects(cap_a, cap_b)
```
"""
struct DilatedIntersects
    radius::Float64
end

@inline (p::DilatedIntersects)(a::Cap, b::Cap) =
    Extents.intersects(_dilatedcap(a, p.radius), b)

# --------------------------------------------------------------------------
# The chunk index
# --------------------------------------------------------------------------

struct EmptyChunkIndex end

function _packedchunkindex(caps::AbstractVector{<:SphericalCap})
    isempty(caps) && return EmptyChunkIndex()
    data = collect(Cap, caps)
    boxes = map(cap -> convert(Extents.Extent, cap), data)
    return RTree(HPR(), data; extents = boxes)
end

"""
    chunkindex(grid) -> index

The source-chunk query object [`candidatechunks!`](@ref) searches: the extents
of `chunks(grid)` (`grid_extent` of each chunk) packed in a GeometryOps
`FlexibleRTree`. Built per plan, never stored on the grid.
"""
chunkindex(grid) = _packedchunkindex(_chunkextents(grid))

@inline function _dilatedcap(cap::Cap, radius::Float64)
    radius == 0.0 && return cap
    return SphericalCap(cap.point, min(Float64(pi), cap.radius + radius))
end

"""
    candidatechunks!(out::Vector{Int}, index, dstcap; radius = 0.0) -> out

Replace `out` with ascending unique chunk numbers from `index` that may lie
within `radius` radians of `dstcap`: the chunks whose extent meets the
dilated cap exactly, so the relation is the cap join. False positives are
allowed; false negatives are not.
"""
function candidatechunks!(out::Vector{Int}, ::EmptyChunkIndex, ::Cap;
        radius::Real = 0.0)
    empty!(out)
    return out
end

function candidatechunks!(out::Vector{Int}, index::RTree, dstcap::Cap;
        radius::Real = 0.0)
    r = Float64(radius)
    empty!(out)
    querycap = _dilatedcap(dstcap, r)
    broad = Base.Fix1(Extents.intersects, querycap)
    exact = DilatedIntersects(r)
    STI.depth_first_search(broad, index) do chunk
        exact(dstcap, index.data[chunk]) && push!(out, chunk)
    end
    sort!(out)
    unique!(out)
    return out
end

# Any other SpatialTreeInterface tree over chunk numbers with cap extents.
function candidatechunks!(out::Vector{Int}, index, dstcap::Cap; radius::Real = 0.0)
    empty!(out)
    predicate = Base.Fix1(DilatedIntersects(Float64(radius)), dstcap)
    STI.depth_first_search(predicate, index) do chunk
        push!(out, chunk)
    end
    sort!(out)
    unique!(out)
    return out
end

# Where the chunk-to-chunk relation is spelled:
#
#   - one destination chunk's sources: `sourcesof(dependencies(plan), d)`, or
#     `sourcesof(chunk_dependency_graph(dst, src; radius), d)` without a plan;
#   - several destinations' union: `_unionrows!` over those rows, which is what
#     a derived lazy tile takes;
#   - one prebuilt index against one cap: `candidatechunks!` itself, which is
#     the seam all of the above are built from and the only query implementation
#     that defines a graph edge.
#
# `_chunkextents` answers none of them: it hands out the caps as *values*, for
# `spacestamp`, `_builddependencies`, `subspace_dependencies` and the
# `chunkindex` that packs them.
