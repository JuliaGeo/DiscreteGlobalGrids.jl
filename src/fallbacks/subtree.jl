# The eager subtree verbs. Both are `collect` of the lazy iterators in
# `subtree_iterators.jl`, which is also where the per-system fast paths hang —
# so a system writes its walker once and both faces of it get faster.

# `collect`, plus the one check it cannot do for us. Given `HasLength`, `collect`
# sizes its vector from `length` and fills by iterating; a walk that emits fewer
# cells than its closed-form count claims therefore returns a tail of `undef` —
# arbitrary integers handed back as cell ids, silently. HEALPix guarded exactly
# this for its own rim; every automaton with a counted rim needs it, so it lives
# here. One comparison against Θ(rim) work.
#
# Takes any counted walk, not just the two iterators — `border_descendants` runs
# an engine directly — so the message names the walk by `show` rather than by
# fields only the iterators have.
#
# `hint` is an APPROXIMATE size for a walk that declares none, and it is a
# second positional argument rather than anything `IteratorSize` can see. That
# separation is the point: `counted` still comes from `IteratorSize` alone, so
# the guard below arms on exactly the walks that promise a count, and a number
# that may be wrong reaches `sizehint!` and nothing else. A hint that is high
# costs a roomier `Vector`; a hint that is low costs one reallocation. Neither
# can put an `undef` slot in a caller's hands, which is what a wrong `length`
# does. See `halo_sizehint`.
function collect_subtree(it, hint::Union{Int,Nothing} = nothing)
    counted = Base.IteratorSize(typeof(it)) isa Base.HasLength
    out = eltype(it)[]
    if counted
        sizehint!(out, length(it))
    elseif hint !== nothing
        sizehint!(out, hint)
    end
    for c in it
        push!(out, c)
    end
    counted && length(out) != length(it) && error(
        "$(nameof(typeof(it))) walked $(length(out)) cells but counts \
         $(length(it)): $it")
    return out
end

"""
    subtree_border(sys, c, l; connectivity = Vertex()) -> Vector

Return level-`l` descendants of `c` with a neighbor outside the subtree, in
ascending canonical order. `collect` of [`EdgeCellIterator`](@ref); the generic
walk costs `O(subtree · degree)`, a system's automaton `O(rim)`.
"""
subtree_border(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
    l::Integer; connectivity::Connectivity=Vertex()) =
    collect_subtree(EdgeCellIterator(sys, c, l; connectivity))

"""
    subtree_interior(sys, c, l; connectivity = Vertex()) -> Vector

Return level-`l` descendants of `c` excluding [`subtree_border`](@ref), in
ascending canonical order. `collect` of [`InnerCellIterator`](@ref), which
generates the interior from the rim walk's pruned branches rather than by
subtracting a border set.
"""
subtree_interior(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
    l::Integer; connectivity::Connectivity=Vertex()) =
    collect_subtree(InnerCellIterator(sys, c, l; connectivity))
