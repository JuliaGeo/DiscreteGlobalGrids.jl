# Discover connected source and destination chunks through each space's native
# spatial index. Unlike cell trees, chunk indexes need not expose one common
# node-extent representation.

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
# Chunk-candidate index implementations
# --------------------------------------------------------------------------

struct EmptyChunkIndex end

function _packedchunkindex(caps::AbstractVector{<:SphericalCap})
    isempty(caps) && return EmptyChunkIndex()
    data = collect(Cap, caps)
    boxes = map(cap -> convert(Extents.Extent, cap), data)
    return FlexibleRTrees.RTree(FlexibleRTrees.HPR(), data; extents = boxes)
end

chunkindex(space::RegridSpace) = _packedchunkindex(chunkextents(space))
chunkindex(space::RasterGrid) = _rasterchunkcursor(space)

@inline function _dilatedcap(cap::Cap, radius::Float64)
    radius == 0.0 && return cap
    return SphericalCap(cap.point, min(Float64(pi), cap.radius + radius))
end

function candidatechunks!(out::Vector{Int}, ::EmptyChunkIndex, ::Cap;
        radius::Real = 0.0)
    empty!(out)
    return out
end

function candidatechunks!(out::Vector{Int}, index::FlexibleRTrees.RTree, dstcap::Cap;
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

# A RasterGridView is indexed in array storage order. Convert its cursor ranges
# back to X/Y before looking up the owning DiskArrays chunks.
@inline function _rastercursorbox(node::Trees.TopDownQuadtreeCursor{<:RasterGridView})
    a, b = node.leafranges
    return node.grid.space.xfast ? (a, b) : (b, a)
end

@inline function _rasterchunkspan(node::Trees.TopDownQuadtreeCursor{<:RasterGridView})
    space = node.grid.space
    xr, yr = _rastercursorbox(node)
    return (_chunkofindex(space.xchunks, first(xr)),
        _chunkofindex(space.xchunks, last(xr)),
        _chunkofindex(space.ychunks, first(yr)),
        _chunkofindex(space.ychunks, last(yr)))
end

function _rastercandidates!(out::Vector{Int},
        node::Trees.TopDownQuadtreeCursor{<:RasterGridView}, intersects)
    intersects(STI.node_extent(node)) || return out
    cx0, cx1, cy0, cy1 = _rasterchunkspan(node)
    space = node.grid.space
    if cx0 == cx1 && cy0 == cy1
        push!(out, chunkposition(space, cx0, cy0))
    elseif STI.isleaf(node)
        for cy in cy0:cy1, cx in cx0:cx1
            push!(out, chunkposition(space, cx, cy))
        end
    else
        for child in STI.getchild(node)
            _rastercandidates!(out, child, intersects)
        end
    end
    return out
end

function candidatechunks!(out::Vector{Int},
        index::Trees.TopDownQuadtreeCursor{<:RasterGridView}, dstcap::Cap;
        radius::Real = 0.0)
    empty!(out)
    # `index` may be retained by a LazyRegridArray. Traverse a private cursor
    # copy so a task-owned chart wrapper never escapes into shared index state.
    localindex = _task_prepared_raster_tree(index)
    _rastercandidates!(out, localindex,
        Base.Fix1(DilatedIntersects(Float64(radius)), dstcap))
    sort!(out)
    unique!(out)
    return out
end

# Compatibility fallback for existing extension trees. New structured indexes
# specialize `candidatechunks!` directly and need not impersonate cap trees.
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

chunkextent(space::RegridSpace, chunk::Integer) = chunkextents(space)[Int(chunk)]

# Three spellings of one relation lived below this line until Phase 4.
#
# `connectedchunkpairs(f, dst_space, src_space; radius)` went in Task G4:
# since PR #69 it was line for line the loop `_chunkgraph` runs to fill its
# destination-major rows.
#
# `connectedchunks(dst_space, dstchunk, src_space; radius)` and its
# `connectedchunks!` in-place forms went in Task E2, together with the
# `chunktree`-collecting `chunkextents` fallback that gave a space a *second*
# way to answer a chunk query. Since Task E1 the lazy executor takes a tile's
# sources from `sourcesof(dependencies(plan), d)`, so the last caller of a
# one-off chunk query was a test. What they did is now spelled:
#
#   - one destination chunk's sources: `sourcesof(dependencies(plan), d)`, or
#     `sourcesof(chunk_dependency_graph(dst, src; radius), d)` without a plan;
#   - several destinations' union: `_unionrows!` over those rows, which is what
#     a derived lazy tile takes;
#   - one prebuilt index against one cap: `candidatechunks!` itself, which is
#     the seam all of the above are built from and the only query implementation
#     that defines a graph edge.
#
# `chunkextents` remains, because it has real consumers that are not queries:
# `spacestamp`, `_builddependencies` and `subspace_dependencies` need the caps
# as *values*, and the generic `chunkindex` packs them. It is a required hook
# now rather than a fallback over a compatibility tree.
