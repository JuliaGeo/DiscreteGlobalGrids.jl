# ---------------------------------------------------------------------------
# Field requests: the per-neighbour quantities a sweep streams.
#
# A request is a tuple of `AbstractNeed`s, resolved once per cell against the
# single membership clip `_indexed` already produced. Every function here
# either dispatches on one need or walks the tuple with `map`/`Base.tail`, so
# a heterogeneous request stays inferable and the rings stay stack-allocated
# wherever the plain sweep's do.
# ---------------------------------------------------------------------------

"""
    AbstractNeed

One per-neighbour quantity a neighbourhood sweep streams: [`Cell`](@ref),
[`Index`](@ref DiscreteGlobalGrids.Index), [`Value`](@ref) or
[`Centroid`](@ref).

A tuple of these is the `needs` keyword of [`mapneighbors`](@ref) and
[`foreachneighbors`](@ref). The callback is then `f(center, rings)`, where
`center` holds one entry per need for the visited cell and `rings` holds one
ring per need — field-major, so `rings[j]` is need `j`'s value for every
clipped neighbour and slot `i` of every ring names the same neighbour.
"""
abstract type AbstractNeed end

"""
    Cell()

Request each neighbour's cell identity, in the system's canonical id scheme.
The ring's element type is the collection's own (`eltype(cv)`), and the
center entry is the visited cell. Use [`Index`](@ref DiscreteGlobalGrids.Index)
to ask for the same cell in another scheme.
"""
struct Cell <: AbstractNeed end

"""
    Local()

The index space of the collection the sweep was called on: `1:length(cv)`,
the same numbers `localindex` answers with. An argument of
[`Index`](@ref DiscreteGlobalGrids.Index).
"""
struct Local end

"""
    Global()

The index space of the complete grid at the collection's level: the numbers
`globalindex` answers with. An argument of
[`Index`](@ref DiscreteGlobalGrids.Index).
"""
struct Global end

"""
    Index(Local())
    Index(Global())
    Index(T::Type{<:AbstractCellIndex})

Request each neighbour's index in one named space: [`Local`](@ref) for the
collection's own `1:length(cv)`, [`Global`](@ref) for the complete grid at
that level, or an id type listed in `cellindextypes(system(cv))` for the same
cell re-encoded, as `reindex` answers.

`Index(Local())` names the collection the caller passed — the `CellVector`,
or the cube's cell axis — not any chunk a sweep splits it into. The two
integer spaces coincide on a complete grid and differ on every subset.
"""
struct Index{S} <: AbstractNeed
    space::S
    # An inner constructor suppresses the permissive default outer one, so the
    # three spaces named below are the only arguments that build an `Index`
    # and anything else is an `ArgumentError` at the call site rather than a
    # `MethodError` from inside the sweep.
    Index{S}(space::S) where {S} = new{S}(space)
end

Index(space::Local) = Index{Local}(space)
Index(space::Global) = Index{Global}(space)
Index(::Type{T}) where {T<:AbstractCellIndex} = Index{Type{T}}(T)
Index(x) = _need_space(x)

"""
    Value(data::AbstractVector)

Request each neighbour's entry in `data`, which must be laid out against the
collection: index `k` of the vector is index `k` of the collection, so
`axes(data) == (Base.OneTo(length(cv)),)`. The ring's element type is
`eltype(data)`.

Any number of `Value`s may appear in one request, each with its own element
type; the sweep reads them all from the one membership clip it already made
for the cell.
"""
struct Value{A<:AbstractVector} <: AbstractNeed
    data::A
end

Value(x) = _need_vector(x)

"""
    Centroid()

Request each neighbour's `cell_centroid` on the unit sphere, as a
`GO.UnitSphericalPoint{Float64}`. Distances and bearings need a radius and
stay downstream of the sweep.

`Centroid()` is `Value(`[`cellfield`](@ref)`(cell_centroid, cv))` asked for by
name, and the two answer identically: the sweep gives each task a bounded
window over the field, keyed by local index, and computes an entry on its
first read inside that window. The values are exactly what `cell_centroid`
answers; what the window changes is how often it is called — about once per
cell wherever the visit order keeps neighbours close in local index, which
storage order and a locality-preserving `order` do and a random permutation
does not.

Spell the field out to hand the sweep what you have already computed:
`Value(cellfield(cell_centroid, cv; known = table))` for the whole
collection, read straight through with no window at all, or a cube over a
subset for the part of it you have.
"""
struct Centroid <: AbstractNeed end

# --- resolution ------------------------------------------------------------

# A request is resolved once, at the sweep's entry: every need that names a
# quantity rather than a place to read it becomes the `Value` of a field, and
# nothing below this line knows that centroids are special. The field is built
# here, once per sweep, and shared by every task — it is pure, so that is safe.
_resolve(::Centroid, cv::CellVector) = Value(cellfield(cell_centroid, cv))
_resolve(n::AbstractNeed, ::CellVector) = n

_resolveneeds(needs::Tuple, cv::CellVector) = map(n -> _resolve(n, cv), needs)

# Each task reads through its own reader: a plain vector is its own, and a
# field gets the bounded window built inside that task from its range.
_taskneed(n::Value, r::UnitRange{Int}) = Value(_taskreader(n.data, r))
_taskneed(n::AbstractNeed, ::UnitRange{Int}) = n

_taskneeds(needs::Tuple, r::UnitRange{Int}) = map(n -> _taskneed(n, r), needs)

# --- element types ---------------------------------------------------------

# The element type one need's ring carries, and the type of its center entry.
_needtype(::Cell, cv::CellVector) = eltype(cv)
_needtype(::Index{Local}, ::CellVector) = Int
_needtype(::Index{Global}, ::CellVector) = Int
_needtype(::Index{Type{T}}, ::CellVector) where {T<:AbstractCellIndex} = T
_needtype(n::Value, ::CellVector) = eltype(n.data)

_centertype(needs::Tuple, cv::CellVector) =
    Tuple{map(n -> _needtype(n, cv), needs)...}

_ringstype(needs::Tuple, cv::CellVector, cap) =
    Tuple{map(n -> _ringtype(cap, _needtype(n, cv)), needs)...}

# --- per-need answers ------------------------------------------------------

# The visited cell, named by its local index `k` and its handle `c`.
@inline _center(::Cell, ::CellVector, k::Int, c) = cellid(c)
@inline _center(::Index{Local}, ::CellVector, k::Int, c) = k
@inline _center(::Index{Global}, cv::CellVector, k::Int, c) =
    globalindex(cv.grid, c)::Int
@inline _center(::Index{Type{T}}, cv::CellVector, k::Int, c) where {T} =
    reindex(T, system(cv), c)
# Both names for the visited cell reach the reader: `k` is what it is keyed
# by, and `cellid(c)` is the cell a field would otherwise decode again.
@inline _center(n::Value, ::CellVector, k::Int, c) = _read(n.data, k, cellid(c))

# One ring slot, named by the clipped handle the sweep produced. The handle
# carries both of the neighbour's names — its local index, which is the key a
# field's reader is built on, so a neighbour and the visit that later lands on
# it share a slot, and the cell itself, which the clip already resolved.
@inline _slot(::Cell, ::CellVector, h::SubsetIndexedCell) = h.cell
@inline _slot(::Index{Local}, ::CellVector, h::SubsetIndexedCell) = h.index
@inline _slot(::Index{Global}, cv::CellVector, h::SubsetIndexedCell) =
    globalindex(cv.grid, h.cell)::Int
@inline _slot(::Index{Type{T}}, cv::CellVector,
    h::SubsetIndexedCell) where {T} = reindex(T, system(cv), h.cell)
@inline _slot(n::Value, ::CellVector, h::SubsetIndexedCell) =
    _read(n.data, h.index, h.cell)

# --- rings -----------------------------------------------------------------

# One need's answers for the whole clipped ring, in the container the capacity
# witness names — the same witness `_ringstype` derives the declared type from,
# so the ring the callback receives and the type the output was sized for
# cannot disagree. A declared bound keeps the row stack-allocated.
@inline function _ring(n::AbstractNeed, cv::CellVector, nbrs, cap::Val)
    out = _newring(cap, _needtype(n, cv))
    for h in nbrs
        out = _pushring(out, _slot(n, cv, h))
    end
    return out
end

@inline function _ring(n::AbstractNeed, cv::CellVector, nbrs, ::Nothing)
    out = Vector{_needtype(n, cv)}(undef, length(nbrs))
    for (i, h) in enumerate(nbrs)
        @inbounds out[i] = _slot(n, cv, h)
    end
    return out
end

# The two tuples the callback receives. `map` over the request, never a loop:
# the entries have different types and the tuple must stay inferable.
@inline _centers(needs::Tuple, cv::CellVector, k::Int, c) =
    map(n -> _center(n, cv, k, c), needs)

@inline _rings(needs::Tuple, cv::CellVector, nbrs, cap) =
    map(n -> _ring(n, cv, nbrs, cap), needs)

# --- validation ------------------------------------------------------------

@noinline _need_empty() = throw(ArgumentError(
    "needs must name at least one field, got ()"))

@noinline _need_kind(n) = throw(ArgumentError(
    "each entry of needs must be an AbstractNeed — Cell(), Index(Local()), " *
    "Index(Global()), Index(T), Value(data) or Centroid() — got $(typeof(n))"))

# A lone request behind a missing comma is not a tuple. Name the value, not
# just its type, and spell the form that would have worked.
_tuple_hint(n::AbstractNeed) = " — did you mean `($(repr(n)),)`?"
_tuple_hint(_) = ""

@noinline _need_tuple(needs) = throw(ArgumentError(
    "needs must be a tuple of field requests, got $(repr(needs))" *
    _tuple_hint(needs)))

# `Index(Local)` names the type where the space was meant — the typo the
# documented `Index(T)` form invites — so the message points at the instance.
_space_hint(::Type{Local}) = "; `Local` is the type, `Local()` the space"
_space_hint(::Type{Global}) = "; `Global` is the type, `Global()` the space"
_space_hint(_) = ""

@noinline _need_space(x) = throw(ArgumentError(
    "Index takes one index space: Index(Local()), Index(Global()), or " *
    "Index(T) for a cell id type T in cellindextypes(system(cv)) — got " *
    "$(repr(x))$(_space_hint(x))"))

@noinline _need_vector(x) = throw(ArgumentError(
    "Value takes one vector laid out on the collection's cell axis — index " *
    "k of the vector is index k of the collection — got a $(typeof(x))"))

@noinline _need_idtype(::Type{T}, cv::CellVector) where {T} = throw(ArgumentError(
    "Index(T) takes an id type the collection's system declares: " *
    "$(system(cv)) accepts $(join(cellindextypes(system(cv)), ", ")), " *
    "got $T"))

_checkneed(n::Value, cv::CellVector) =
    (_check_data(n.data, length(cv)); _checkover(n.data, cv))

# A stored vector says which cell an entry is for by its axis index alone, so
# the layout check above is the whole contract. A field also carries the
# collection it computes against, and a sweep reads it by local index: over
# any other collection those indices name other cells and every answer would
# be silently wrong. Identity settles it in one comparison for the field a
# sweep built itself; otherwise `CellVector`'s own equality does, which is
# system, level and window bounds — O(number of windows), and O(n) over leaf
# indices only when the two window representations differ.
_checkover(::AbstractVector, ::CellVector) = nothing
_checkover(a::CellField, cv::CellVector) =
    (a.cv === cv || a.cv == cv) ? nothing : _field_collection(a.cv, cv)

# Two collections of the same length are the case the layout check cannot
# see, so the message names the cells rather than the count.
_cvbrief(cv::CellVector) = isempty(cv) ? "no cells of $(system(cv)) level $(level(cv))" :
                           "$(length(cv)) cells of $(system(cv)) level " *
                           "$(level(cv)), $(cv[1]) to $(cv[end])"

@noinline _field_collection(from::CellVector, cv::CellVector) =
    throw(ArgumentError(
        "a cell field is read by local index, so it must be over the " *
        "collection being swept, and this one is not: the field is over " *
        "$(_cvbrief(from)); the sweep is over $(_cvbrief(cv))"))
# The scheme is checked here, before the sweep starts, rather than left to the
# first `reindex` call to fail on.
_checkneed(::Index{Type{T}}, cv::CellVector) where {T<:AbstractCellIndex} =
    T in cellindextypes(system(cv)) ? nothing : _need_idtype(T, cv)
_checkneed(::AbstractNeed, ::CellVector) = nothing
_checkneed(n, ::CellVector) = _need_kind(n)

# Recursion, not iteration: the entries have different types.
_checkeach(::Tuple{}, ::CellVector) = nothing
_checkeach(needs::Tuple, cv::CellVector) =
    (_checkneed(first(needs), cv); _checkeach(Base.tail(needs), cv))

_checkneeds(needs::Tuple, cv::CellVector) = _checkeach(needs, cv)
_checkneeds(::Tuple{}, ::CellVector) = _need_empty()
_checkneeds(needs, ::CellVector) = _need_tuple(needs)

# --- entry points ----------------------------------------------------------

# The `needs` half of `mapneighbors`/`foreachneighbors`; the `nothing` half is
# in neighborhood.jl, and `needs = nothing` reaches it by dispatch.
function _mapneighbors(f::F, cv::CellVector, needs, order, threaded,
        connectivity::Connectivity) where {F}
    # Check what the caller wrote, then resolve it: an error names the request
    # the caller made, not the field it would have become.
    _checkneeds(needs, cv)
    rn = _resolveneeds(needs, cv)
    cap = _capacity(system(cv), connectivity)
    T = Base.promote_op(f, _centertype(rn, cv), _ringstype(rn, cv, cap))
    outs = _outputs(T, length(cv))
    return _mapstore!(f, outs, cv, rn, connectivity, order,
        GOCore.booltype(threaded), cap)
end

# Built after the output container type is known, like the other `_mapstore!`s.
# The callback is built once per task rather than once per sweep, so the
# readers it closes over belong to that task alone.
function _mapstore!(f::F, outs::O, cv::CellVector, needs::Tuple,
        conn::Connectivity, order, thr, cap::CAP) where {F,O,CAP}
    _runeach!(cv, conn, order, thr, cap) do r
        tn = _taskneeds(needs, r)
        (k, c, nbrs) -> _store!(outs, k,
            f(_centers(tn, cv, k, c), _rings(tn, cv, nbrs, cap)))
    end
    return outs
end

# Which visited cells the callback is called for. A sweep over one chunk of a
# larger axis computes rings for the chunk's halo too — those cells are read to
# complete the owned cells' rings, and their own rings are incomplete — and a
# callback that runs for its effects has no result to discard afterwards, so
# the chunk route hands its owned range down here. Every other caller visits
# every cell, through a singleton the branch folds away against, so the plain
# sweep keeps its zero-allocation shape.
struct _EveryCell end

@inline (::_EveryCell)(::Int) = true

function _foreachneighbors(f::F, cv::CellVector, needs, order, threaded,
        connectivity::Connectivity, keep::K=_EveryCell()) where {F,K}
    _checkneeds(needs, cv)
    rn = _resolveneeds(needs, cv)
    # One capacity witness for the whole sweep: the clip, the rings it feeds
    # and the ring type all come from this value.
    cap = _capacity(system(cv), connectivity)
    _runeach!(cv, connectivity, order, GOCore.booltype(threaded), cap) do r
        tn = _taskneeds(rn, r)
        (k, c, nbrs) -> (keep(k) && f(_centers(tn, cv, k, c),
            _rings(tn, cv, nbrs, cap)); nothing)
    end
    return nothing
end
