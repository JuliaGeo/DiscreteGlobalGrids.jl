# ---------------------------------------------------------------------------
# Cell fields: a vector over a collection that knows some of its entries and
# computes the rest.
#
# A field is an ordinary `AbstractVector` laid out on the collection's cell
# axis, so anything that accepts such a vector accepts one — including the
# `Value` field request. Reading it is pure and the object is never mutated,
# which is what lets one field be shared by every task of a threaded sweep;
# the caching is a per-task reader built around it, below.
# ---------------------------------------------------------------------------

# What the caller already knows, when it covers only part of the collection:
# the values, and the collection they are laid out on. Membership is that
# collection's own `localindex`, which answers `nothing` off it.
struct _SubsetKnown{S,V<:AbstractVector}
    on::S
    values::V
end

@inline function _knownat(s::_SubsetKnown, c)
    j = localindex(s.on, c)
    j === nothing && return nothing
    return @inbounds s.values[j]
end

struct CellField{T,F,CV<:CellVector,K} <: AbstractVector{T}
    f::F
    cv::CV
    known::K
end

"""
    cellfield(f, cv::CellVector; known = nothing) -> CellField

A vector over `cv`'s cell axis whose `k`th entry is `f(cv.grid, cv[k])` —
`cellfield(cell_centroid, cv)` is the collection's centroid field — with
`known` naming entries the caller has already computed, which are read instead
of called for.

`known` accepts

  - `nothing`, the default: every entry is computed on read;
  - an `AbstractVector` with axis `1:length(cv)`: the field is complete and
    nothing is ever computed;
  - a one-dimensional cube over a subset of `cv` — a `DimensionalData` array
    on a `Cells` dimension at the same system and level — whose cells
    are read from it and whose absent cells are computed.

The element type is `f`'s return type for one cell of `cv`, inferred once when
the field is built; an `f` whose return type cannot be inferred is an
`ArgumentError`, as is a `known` whose elements are not of that type or whose
shape or extent is none of the three above.

Reading is pure: the field never mutates and never remembers, so one field is
safely shared by concurrent readers and reading the same index twice costs
what reading it once did. Remembering belongs to the reader — a sweep handed
`Value(field)` gives each of its tasks a bounded window over the field and
computes an entry on its first read inside that window, which is why the field
must be over the collection being swept and is rejected when it is not.
"""
function cellfield(f::F, cv::CellVector; known = nothing) where {F}
    T = Base.promote_op(f, typeof(cv.grid), eltype(cv))
    isconcretetype(T) || _field_eltype(f, cv, T)
    kn = _cellknown(known, cv, T)
    return CellField{T,F,typeof(cv),typeof(kn)}(f, cv, kn)
end

Base.size(a::CellField) = (length(a.cv),)
Base.IndexStyle(::Type{<:CellField}) = Base.IndexLinear()

# The three knowledge shapes, one method each rather than one branch, so a
# reader of a complete field carries no probe and a reader of an unknown field
# carries no lookup. Each takes the cell as well as its index: whoever asks
# usually holds it already — a sweep has just decoded it to find its
# neighbours — and decoding it again is a lattice walk this layer does not
# need to repeat. `_fieldat(a, k)` is the form for a caller that holds only
# the index, and it is the one place the cell is decoded.
@inline _fieldat(a::CellField{T,F,CV,Nothing}, k::Int, c) where {T,F,CV} =
    a.f(a.cv.grid, c)::T

@inline _fieldat(a::CellField{T,F,CV,K}, k::Int, c) where {T,F,CV<:CellVector,
    K<:AbstractVector} = @inbounds a.known[k]

@inline function _fieldat(a::CellField{T,F,CV,K}, k::Int,
        c) where {T,F,CV,K<:_SubsetKnown}
    v = _knownat(a.known, c)
    v === nothing || return v::T
    return a.f(a.cv.grid, c)::T
end

@inline _fieldat(a::CellField, k::Int) = _fieldat(a, k, @inbounds a.cv[k])

Base.@propagate_inbounds function Base.getindex(a::CellField, k::Int)
    @boundscheck checkbounds(a, k)
    return _fieldat(a, k)
end

# --- what `known` may be ---------------------------------------------------

_cellknown(::Nothing, ::CellVector, ::Type) = nothing

function _cellknown(v::AbstractVector, cv::CellVector, ::Type{T}) where {T}
    axes(v) == (Base.OneTo(length(cv)),) || _field_known_axes(v, cv)
    eltype(v) <: T || _field_known_eltype(eltype(v), T)
    return v
end

# The cube shape is named by the DimensionalData layer, which sits above this
# one; it adds its own method (`src/dimensionaldata.jl`).
_cellknown(x, ::CellVector, ::Type) = _field_known_kind(x)

@noinline _field_eltype(f, cv::CellVector, ::Type{T}) where {T} =
    throw(ArgumentError(
        "cellfield needs one element type for the whole field, but $f " *
        "applied to $(typeof(cv.grid)) and $(eltype(cv)) is inferred as $T; " *
        "annotate the function's return type or pass a complete `known`"))

@noinline _field_known_axes(v, cv::CellVector) = throw(ArgumentError(
    "a vector `known` must be laid out on the whole collection — index k of " *
    "the vector is index k of the collection — so it needs axis " *
    "1:$(length(cv)), got $(axes(v, 1)); a partial one is a cube on a Cells " *
    "dimension"))

@noinline _field_known_eltype(::Type{K}, ::Type{T}) where {K,T} =
    throw(ArgumentError(
        "`known` holds $K where the field's element type is $T"))

@noinline _field_known_kind(x) = throw(ArgumentError(
    "`known` is nothing, a vector on the collection's cell axis, or a cube " *
    "on a Cells dimension over a subset of it, got a $(typeof(x))"))

# --- the per-task reader ---------------------------------------------------

# A cell's field entry is read `1 + degree` times over a sweep — once as the
# visited cell, once from each neighbour that names it — and those reads are
# close together in local index wherever the traversal has locality. So one
# task keeps a direct-mapped cache of `W` recent entries keyed by local index:
# slot `k & (W - 1) + 1`, a hit iff that slot's tag is `k`. No replacement
# policy and no recency bookkeeping — a fully associative LRU of the same size
# was measured against it and buys at most 0.008 computations per cell.
#
# The key is the local index, so the hit rate is the visit order's property,
# not the window's: storage order and any locality-preserving `order` keep it,
# while a random permutation destroys it. Local indices start at 1, so a zero
# tag is an empty slot and no separate valid bit is needed.
#
# The window is per task and never shared: two tasks writing one buffer would
# race, and their ranges have no locality with each other anyway. It is the
# only mutable thing in this file, and the field it reads from stays pure.
struct _FieldWindow{T,A<:CellField} <: AbstractVector{T}
    field::A
    tags::Vector{Int}
    vals::Vector{T}
    mask::Int
end

# `W` is a power of two so the slot is a mask, no larger than the task itself
# so a short sweep carries a short window, and capped at 16384 slots because
# every registered system is already down to 1.008-1.043 computations per cell
# there and a bigger window is memory spent on nothing.
function _FieldWindow(a::CellField{T}, n::Int) where {T}
    W = n >= 16384 ? 16384 : nextpow(2, max(n, 1))
    return _FieldWindow{T,typeof(a)}(a, zeros(Int, W), Vector{T}(undef, W),
        W - 1)
end

Base.size(w::_FieldWindow) = size(w.field)
Base.IndexStyle(::Type{<:_FieldWindow}) = Base.IndexLinear()

@inline function Base.getindex(w::_FieldWindow{T}, k::Int) where {T}
    slot = (k & w.mask) + 1
    if (@inbounds w.tags[slot]) == k
        return @inbounds w.vals[slot]
    end
    return _window_decode(w, k, slot)
end

# The miss is deliberately out of line. A hit is a tag compare and a load, and
# it is what six of every seven reads are; inlining the field's own call — a
# handful of trig calls, for a centroid — beside it bloats the ring loop the
# sweep spends its time in for the sake of the seventh. What `known` already
# answered is inserted too, so a cell the caller precomputed is probed once
# rather than once per touch.
@noinline function _window_miss(w::_FieldWindow{T}, k::Int, slot::Int,
        c) where {T}
    v = _fieldat(w.field, k, c)::T
    @inbounds w.tags[slot] = k
    @inbounds w.vals[slot] = v
    return v
end

# The miss taken by a reader that was handed only an index: decoding the cell
# is part of the miss and belongs out of line with the rest of it.
@noinline _window_decode(w::_FieldWindow, k::Int, slot::Int) =
    _window_miss(w, k, slot, @inbounds w.field.cv[k])

# --- reading with the cell in hand -----------------------------------------

# What a sweep reads a need through. It holds both names for the cell — the
# local index the reader is keyed by, and the cell itself, decoded once to
# find its neighbours — so nothing below has to decode it a second time. A
# stored vector ignores the cell; a window uses it only on a miss.
@inline _read(a::AbstractVector, k::Int, c) = @inbounds a[k]

@inline function _read(w::_FieldWindow{T}, k::Int, c) where {T}
    slot = (k & w.mask) + 1
    if (@inbounds w.tags[slot]) == k
        return @inbounds w.vals[slot]
    end
    return _window_miss(w, k, slot, c)
end

# One reader per task, over the range that task sweeps. An ordinary vector is
# its own reader; a complete field is the vector it was built from; anything
# that computes gets the window.
_taskreader(a::AbstractVector, ::UnitRange{Int}) = a
_taskreader(a::CellField{T,F,CV,K}, ::UnitRange{Int}) where {T,F,CV<:CellVector,
    K<:AbstractVector} = a.known
_taskreader(a::CellField, r::UnitRange{Int}) = _FieldWindow(a, length(r))
