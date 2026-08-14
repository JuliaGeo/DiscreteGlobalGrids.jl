"""
    SmallList{N,T}

An immutable `AbstractVector{T}` with inline storage for at most `N` elements.
It is intended for small, compile-time capacities and isbits element types on
allocation-sensitive paths.
"""
struct SmallList{N,T} <: AbstractVector{T}
    len::Int
    data::NTuple{N,T}
end

Base.IndexStyle(::Type{<:SmallList}) = IndexLinear()
Base.size(list::SmallList) = (list.len,)
Base.length(list::SmallList) = list.len
Base.eltype(::Type{SmallList{N,T}}) where {N,T} = T
Base.similar(::SmallList, ::Type{T}, dims::Dims) where {T} =
    Array{T}(undef, dims)

Base.@propagate_inbounds function Base.getindex(list::SmallList, i::Int)
    @boundscheck checkbounds(list, i)
    @inbounds return list.data[i]
end

@inline Base.iterate(list::SmallList, state::Int=1) =
    state > list.len ? nothing : (@inbounds list.data[state], state + 1)

"""
    empty_small_list(Val(capacity), filler)

Construct an empty `SmallList`. `filler` initializes unused inline storage.
`capacity` should be a compile-time `Val` for the result to remain fully
inferred.
"""
empty_small_list(::Val{N}, filler::T) where {N,T} =
    SmallList{N,T}(0, ntuple(_ -> filler, Val(N)))

"""
    small_push(list, item)

Return a new list with `item` appended. Throws `BoundsError` when full.
"""
@inline function small_push(list::SmallList{N,T}, item::T) where {N,T}
    list.len < N || throw(BoundsError(list, list.len + 1))
    next_index = list.len + 1
    data = ntuple(i -> i == next_index ? item : list.data[i], Val(N))
    return SmallList{N,T}(next_index, data)
end

small_push(list::SmallList{N,T}, item) where {N,T} =
    small_push(list, convert(T, item))

@inline tuple_set(data::NTuple{N,T}, item::T, index::Int) where {N,T} =
    ntuple(i -> i == index ? item : data[i], Val(N))

"""
    small_pop(list)

Return `list` without its last element. The dropped slot keeps its old contents,
which `len` puts out of reach. Throws `BoundsError` when empty.
"""
@inline function small_pop(list::SmallList{N,T}) where {N,T}
    list.len > 0 || throw(BoundsError(list, 0))
    return SmallList{N,T}(list.len - 1, list.data)
end

"""
    small_setlast(list, item)

Return `list` with its last element replaced — the in-place update an immutable
stack needs to advance the frame on top of it. Throws `BoundsError` when empty.
"""
@inline function small_setlast(list::SmallList{N,T}, item::T) where {N,T}
    list.len > 0 || throw(BoundsError(list, 0))
    return SmallList{N,T}(list.len, tuple_set(list.data, item, list.len))
end

small_setlast(list::SmallList{N,T}, item) where {N,T} =
    small_setlast(list, convert(T, item))

"""
    small_sort(list)

Return an ascending copy of a small list without heap allocation.
"""
function small_sort(list::SmallList{N,T}) where {N,T}
    out = SmallList{N,T}(0, list.data)
    used = ntuple(_ -> false, Val(N))
    for _ in 1:list.len
        best_i = 0
        best = list.data[1]
        for i in 1:list.len
            if !used[i] && (best_i == 0 || isless(list.data[i], best))
                best_i = i
                best = list.data[i]
            end
        end
        out = small_push(out, best)
        used = tuple_set(used, true, best_i)
    end
    return out
end
