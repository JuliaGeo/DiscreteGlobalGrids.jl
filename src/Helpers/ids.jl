"""
    sorted_index(values, value) -> Int

Return the index of `value` in sorted `values`, or zero when it is absent.
"""
@inline function sorted_index(values, value)
    i = searchsortedfirst(values, value)
    return i <= lastindex(values) && values[i] == value ? i : 0
end

"""
    strictly_increasing(values) -> Bool

Check sortedness and uniqueness in one allocation-free pass.
"""
function strictly_increasing(values)
    first = iterate(values)
    first === nothing && return true
    previous, state = first
    while true
        next = iterate(values, state)
        next === nothing && return true
        value, state = next
        isless(previous, value) || return false
        previous = value
    end
end

to_uint64_id(id::UInt64) = id
to_uint64_id(id::Unsigned) = UInt64(id)
to_uint64_id(id::Integer) = UInt64(id)

function to_uint64_id(id::AbstractString)
    prefixed = startswith(id, "0x") || startswith(id, "0X")
    clean = prefixed ? SubString(id, 3) : id
    return parse(UInt64, clean; base=16)
end
