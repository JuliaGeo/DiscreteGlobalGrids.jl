# Discover connected source and destination chunks through each space's native
# spatial index. Chunk indexes are private query objects; unlike cell trees,
# they need not expose one common node-extent representation.

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

@inline (p::DilatedIntersects)(a, b) =
    US.spherical_distance(a.point, b.point) <= a.radius + b.radius + p.radius

# --------------------------------------------------------------------------
# Private chunk-candidate index
# --------------------------------------------------------------------------

struct EmptyChunkIndex end

"""
    cap_xyz_extent(cap) -> Extents.Extent{(:X, :Y, :Z)}

Return an outward-rounded Cartesian bounding box for every unit-sphere point in
`cap`. This is the leaf extent used by the generic packed R-tree.
"""
function cap_xyz_extent(cap::Cap)
    r = cap.radius
    if !(r < Float64(pi))
        return Extents.Extent(X = (-1.0, 1.0), Y = (-1.0, 1.0), Z = (-1.0, 1.0))
    end
    cr, sr = cos(r), sin(r)
    guard = 8 * eps(Float64)
    bounds = ntuple(Val(3)) do i
        x = clamp(Float64(cap.point[i]), -1.0, 1.0)
        transverse = sqrt(max(0.0, (1.0 - x) * (1.0 + x)))
        lo = x <= -cr ? -1.0 : x * cr - transverse * sr
        hi = x >= cr ? 1.0 : x * cr + transverse * sr
        (max(-1.0, prevfloat(lo - guard)), min(1.0, nextfloat(hi + guard)))
    end
    return Extents.Extent(X = bounds[1], Y = bounds[2], Z = bounds[3])
end

function _packedchunkindex(caps::AbstractVector{<:SphericalCap})
    isempty(caps) && return EmptyChunkIndex()
    data = collect(Cap, caps)
    boxes = map(cap_xyz_extent, data)
    return FlexibleRTrees.RTree(FlexibleRTrees.HPR(), data; extents = boxes)
end

"""
    chunkindex(space)

Build the private source-chunk query object. Structured spaces specialize this
to expose their implicit hierarchy; the fallback packs their covering caps in a
GeometryOps `FlexibleRTree`.
"""
chunkindex(space::RegridSpace) = _packedchunkindex(chunkextents(space))
chunkindex(space::RasterGrid) = _rasterchunkcursor(space)

@inline function _dilatedcap(cap::Cap, radius::Float64)
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
    _rastercandidates!(out, index, Base.Fix1(DilatedIntersects(Float64(radius)), dstcap))
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

"""
    chunkextents(space::RegridSpace) -> Vector{SphericalCap}

Collect each chunk's spherical-cap extent from [`chunktree`](@ref). Structured
spaces should specialize this when their query index is not a cap tree.
"""
function chunkextents(space::RegridSpace)
    caps = Vector{Cap}(undef, Int(nchunks(space)))
    filled = falses(length(caps))
    _collectextents!(caps, filled, chunktree(space))
    all(filled) || throw(ArgumentError(
        "the chunk tree of $(typeof(space)) does not reach every chunk in " *
        "1:$(length(caps)); chunk tree leaf indices must be chunk numbers"))
    return caps
end

function _collectextents!(caps::Vector{Cap}, filled::BitVector, node)
    if STI.isleaf(node)
        for (i, extent) in STI.child_indices_extents(node)
            1 <= i <= length(caps) || throw(ArgumentError(
                "chunk tree leaf index $i is outside 1:$(length(caps))"))
            caps[i] = extent
            filled[i] = true
        end
    else
        for child in STI.getchild(node)
            _collectextents!(caps, filled, child)
        end
    end
    return caps
end

"""
    chunkextent(space::RegridSpace, chunk::Integer) -> SphericalCap

Return one chunk extent. Spaces may specialize this to avoid walking the tree.
"""
chunkextent(space::RegridSpace, chunk::Integer) = chunkextents(space)[Int(chunk)]

"""
    connectedchunks(dst_space, dstchunk, src_space; radius = 0.0) -> Vector{Int}

Return ascending source chunks within `radius` radians of `dstchunk`. The result
may include false positives, but must include every contributing chunk.
"""
connectedchunks(dst_space::RegridSpace, dstchunk::Integer, src_space::RegridSpace;
    radius::Real = 0.0) =
    connectedchunks!(Int[], chunkextent(dst_space, dstchunk), src_space; radius)

"""
    connectedchunks!(out, dstcap::SphericalCap, src_space; radius = 0.0) -> out
    connectedchunks!(out, dstcap::SphericalCap, srcindex; radius = 0.0) -> out

Write connected source chunks for `dstcap` into `out`. Passing a prebuilt source
index avoids rebuilding it for repeated queries.
"""
connectedchunks!(out::Vector{Int}, dstcap::Cap, src_space::RegridSpace;
    radius::Real = 0.0) =
    candidatechunks!(out, chunkindex(src_space), dstcap; radius)

connectedchunks!(out::Vector{Int}, dstcap::Cap, srcindex; radius::Real = 0.0) =
    candidatechunks!(out, srcindex, dstcap; radius)

"""
    connectedchunks!(out, dstcaps::AbstractVector{<:SphericalCap}, srcindex; radius = 0.0)

Write the union of connected chunks for several destination extents into `out`.
Querying caps separately avoids the loose bound from merging distant caps.
"""
function connectedchunks!(out::Vector{Int}, dstcaps::AbstractVector{<:SphericalCap},
    srcindex; radius::Real = 0.0)
    empty!(out)
    buffer = Int[]
    for cap in dstcaps
        candidatechunks!(buffer, srcindex, cap; radius)
        append!(out, buffer)
    end
    sort!(out)
    unique!(out)
    return out
end

"""
    connectedchunkpairs(f, dst_space, src_space; radius = 0.0)

Call `f(dstchunk, srcchunk)` for every potentially contributing pair through
the source space's native chunk index.
"""
function connectedchunkpairs(f::F, dst_space::RegridSpace, src_space::RegridSpace;
    radius::Real = 0.0) where {F}
    index = chunkindex(src_space)
    out = Int[]
    for (d, cap) in pairs(chunkextents(dst_space))
        candidatechunks!(out, index, cap; radius)
        for s in out
            f(d, s)
        end
    end
    return nothing
end
