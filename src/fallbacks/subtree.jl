# Eager subtree operations collect the lazy, system-specializable iterators.

# Collect an iterator while validating any declared length. An incorrect
# `HasLength` count would otherwise leave uninitialized output slots. `hint` is
# used only by `sizehint!` for iterators whose size is unknown.
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
walk costs `O(subtree · degree)`, a system's automaton `O(border)`.
"""
subtree_border(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
    l::Integer; connectivity::Connectivity=Vertex()) =
    collect_subtree(EdgeCellIterator(sys, c, l; connectivity))

"""
    subtree_interior(sys, c, l; connectivity = Vertex()) -> Vector

Return level-`l` descendants of `c` excluding [`subtree_border`](@ref), in
ascending canonical order. `collect` of [`InnerCellIterator`](@ref), which
generates the interior from the border walk's pruned branches rather than by
subtracting a border set.
"""
subtree_interior(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
    l::Integer; connectivity::Connectivity=Vertex()) =
    collect_subtree(InnerCellIterator(sys, c, l; connectivity))
