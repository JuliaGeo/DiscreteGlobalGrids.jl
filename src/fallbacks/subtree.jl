# Collecting the lazy, system-specializable walks.

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
