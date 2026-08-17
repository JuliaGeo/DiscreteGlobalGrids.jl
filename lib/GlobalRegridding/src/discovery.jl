# Discover connected source and destination chunks from their cap trees.

"""
    CapQuery(cap)

A one-leaf spatial tree for querying one spherical cap. Its leaf index is `1`.

# Example

```julia
cap = SphericalCap(USPoint(0.0, 0.0, 1.0), 0.1)
query = CapQuery(cap)
STI.node_extent(query) == cap
```
"""
struct CapQuery
    cap::Cap
end

STI.isspatialtree(::Type{CapQuery}) = true
STI.node_extent_is_expensive(::Type{CapQuery}) = false
STI.isleaf(::CapQuery) = true
STI.nchild(::CapQuery) = 0
STI.getchild(::CapQuery) = ()
STI.node_extent(q::CapQuery) = q.cap
STI.child_indices_extents(q::CapQuery) = ((1, q.cap),)

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

"""
    chunkextents(space::RegridSpace) -> Vector{SphericalCap}

Collect each chunk's spherical-cap extent from [`chunktree`](@ref).
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
    connectedchunks!(out, dstcap::SphericalCap, srctree; radius = 0.0) -> out

Write connected source chunks for `dstcap` into `out`. Passing a prebuilt source
tree avoids rebuilding it for repeated queries.
"""
connectedchunks!(out::Vector{Int}, dstcap::Cap, src_space::RegridSpace;
    radius::Real = 0.0) =
    connectedchunks!(out, dstcap, chunktree(src_space); radius)

function connectedchunks!(out::Vector{Int}, dstcap::Cap, srctree; radius::Real = 0.0)
    empty!(out)
    _descend!(out, dstcap, srctree, Float64(radius))
    sort!(out)
    unique!(out)
    return out
end

"""
    connectedchunks!(out, dstcaps::AbstractVector{<:SphericalCap}, srctree; radius = 0.0)

Write the union of connected chunks for several destination extents into `out`.
Querying caps separately avoids the loose bound from merging distant caps.
"""
function connectedchunks!(out::Vector{Int}, dstcaps::AbstractVector{<:SphericalCap},
    srctree; radius::Real = 0.0)
    empty!(out)
    r = Float64(radius)
    for cap in dstcaps
        _descend!(out, cap, srctree, r)
    end
    sort!(out)
    unique!(out)
    return out
end

function _descend!(out::Vector{Int}, dstcap::Cap, srctree, radius::Float64)
    STI.dual_depth_first_search(DilatedIntersects(radius),
        CapQuery(dstcap), srctree) do _, s
        push!(out, s)
    end
    return out
end

"""
    connectedchunkpairs(f, dst_space, src_space; radius = 0.0)

Call `f(dstchunk, srcchunk)` for every potentially contributing pair using one
dual-tree descent.
"""
function connectedchunkpairs(f::F, dst_space::RegridSpace, src_space::RegridSpace;
    radius::Real = 0.0) where {F}
    STI.dual_depth_first_search(f, DilatedIntersects(Float64(radius)),
        chunktree(dst_space), chunktree(src_space))
    return nothing
end
