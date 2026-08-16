# Connected chunk-pair discovery over the two chunk trees. Owned by task T7.
#
# The lazy path needs one question answered: which source chunks can contribute
# to this destination chunk? It is answered by a dual depth-first descent over
# the two chunk trees, whose node extents are spherical caps, with the
# destination side dilated by the method's support radius. Dilation is what
# keeps an interpolation stencil that reaches past a source chunk's boundary
# discoverable; without it the pair is dropped and the missing policy silently
# renormalizes a truncated stencil.

"""
    CapQuery(cap)

A one-leaf spatial tree holding a single [`SphericalCap`] extent, so that a
query about one chunk is the same dual descent as a query about a whole tree.

Its leaf index is always `1`; callers care only about the opposing tree's
indices.
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

The descent predicate: whether two spherical caps come within `radius` angular
radians of each other, i.e. whether the first cap **dilated by `radius`** meets
the second.

Dilating rather than testing bare intersection is the whole of the support-radius
contract, and it is applied at every level of the descent rather than only at the
leaves: a node's cap dilated by `radius` contains every child's cap dilated by
`radius`, so a pair pruned here has no descendant pair that would have survived.
Over-reporting is free — an extra chunk pair yields an all-zero block — while
under-reporting is silent data loss.
"""
struct DilatedIntersects
    radius::Float64
end

@inline (p::DilatedIntersects)(a, b) =
    US.spherical_distance(a.point, b.point) <= a.radius + b.radius + p.radius

"""
    chunkextents(space::RegridSpace) -> Vector{SphericalCap}

Every chunk's extent, indexed by chunk number, collected from
[`chunktree`](@ref)'s leaves.

Geometry only: no cell polygon is built and no data is read.
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

One chunk's extent. The generic implementation walks [`chunktree`](@ref); a
space that can answer in O(1) should say so by defining this.
"""
chunkextent(space::RegridSpace, chunk::Integer) = chunkextents(space)[Int(chunk)]

"""
    connectedchunks(dst_space, dstchunk, src_space; radius = 0.0) -> Vector{Int}

The source chunks that can contribute to `dstchunk`: those whose extent comes
within `radius` angular radians of `dstchunk`'s, ascending.

`radius` is `support_radius(method, src_space)` — zero for a method whose
weights are pure overlap, positive for one whose stencil reaches beyond the cell
it is centred on.

The answer is a superset, never a subset: extents bound their geometry, so a
pair that is dropped cannot contribute, while a pair reported in error costs one
empty block. Every method is linear and every block accumulates, so a superset
gives the same numbers as the exact set.
"""
connectedchunks(dst_space::RegridSpace, dstchunk::Integer, src_space::RegridSpace;
    radius::Real = 0.0) =
    connectedchunks!(Int[], chunkextent(dst_space, dstchunk), src_space; radius)

"""
    connectedchunks!(out, dstcap::SphericalCap, src_space; radius = 0.0) -> out
    connectedchunks!(out, dstcap::SphericalCap, srctree; radius = 0.0) -> out

[`connectedchunks`](@ref) against a destination extent already in hand, into a
reused vector.

The descent is `GeometryOps.SpatialTreeInterface`'s dual depth-first search with
the destination as a one-node tree, so a hierarchical source chunk tree is
pruned level by level and a flat one degrades to the pairwise cap scan that is
all a flat tree can offer.

The second form takes the source's chunk tree directly, for a caller that queries
one source repeatedly: [`chunktree`](@ref) may cost O(nchunks) to build, and that
cost belongs once per plan and not once per destination chunk.
"""
connectedchunks!(out::Vector{Int}, dstcap::Cap, src_space::RegridSpace;
    radius::Real = 0.0) =
    connectedchunks!(out, dstcap, chunktree(src_space); radius)

function connectedchunks!(out::Vector{Int}, dstcap::Cap, srctree; radius::Real = 0.0)
    empty!(out)
    STI.dual_depth_first_search(DilatedIntersects(Float64(radius)),
        CapQuery(dstcap), srctree) do _, s
        push!(out, s)
    end
    sort!(out)
    unique!(out)
    return out
end

"""
    connectedchunkpairs(f, dst_space, src_space; radius = 0.0)

Call `f(dstchunk, srcchunk)` for every chunk pair that can contribute, in one
dual descent over both chunk trees.

The batched form of [`connectedchunks`](@ref), for callers that want every
destination chunk's connections at once: two hierarchical chunk trees prune each
other here, which a per-destination query cannot do.
"""
function connectedchunkpairs(f::F, dst_space::RegridSpace, src_space::RegridSpace;
    radius::Real = 0.0) where {F}
    STI.dual_depth_first_search(f, DilatedIntersects(Float64(radius)),
        chunktree(dst_space), chunktree(src_space))
    return nothing
end
