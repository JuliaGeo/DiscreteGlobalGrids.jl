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

A cell's centroid is touched once as the visited cell and once from each
neighbour that names it, so the sweep keeps a bounded working set of recent
centroids per task and computes each one on its first touch inside that
window. The values are exactly what `cell_centroid` answers; what the window
changes is how often it is called.

The window is keyed by the local index, so it hits wherever the visit order
keeps neighbours close in that index: storage order does, and so does a
locality-preserving `order` permutation; a random one leaves the window
missing on nearly every touch. What survives the window is about one
computation per cell, and most of that is what a full precomputed table
handed in as `Value(table)` removes — faster still, at O(ncells) memory plus
its own build, which is why the window is the default and the table the
caller's opt-in.
"""
struct Centroid <: AbstractNeed end

# --- the centroid working set ----------------------------------------------

# A cell's centroid is touched `1 + degree` times over a sweep — once as the
# visited cell, once from each neighbour that names it — and those touches are
# close together in local index wherever the traversal has locality. So one
# task keeps a direct-mapped cache of `W` recent centroids keyed by local
# index: slot `k & (W - 1) + 1`, a hit iff that slot's tag is `k`. No
# replacement policy and no recency bookkeeping — a fully associative LRU of
# the same size was measured and buys at most 0.008 computations per cell.
#
# The key is the local index, so the hit rate is the visit order's property,
# not the window's: storage order and any locality-preserving `order` keep it,
# while a random permutation destroys it (see the permuted rows of
# `benchmark/needs_centroid.jl`). What is left after the window is about one
# computation per cell — 1.026-1.043 on the registered systems — and that
# residual is most of what a full precomputed table removes on top, so
# `Value(table)` stays worth its O(ncells) memory for a caller who can pay it.
#
# Local indices start at 1, so a zero tag is an empty slot and no separate
# valid bit is needed. The window is per task and never shared: two tasks that
# wrote one buffer would race, and their ranges have no locality with each
# other anyway.
struct _CentroidWindow
    tags::Vector{Int}
    vals::Vector{USPoint}
    mask::Int
end

# `W` is a power of two so the slot is a mask, no larger than the task itself
# so a short sweep carries a short window, and capped at 16384 slots (512 KiB
# of tags and points) because every registered system is already down to
# 1.008-1.043 computations per cell there and a bigger window is memory spent
# on nothing.
function _CentroidWindow(n::Int)
    W = n >= 16384 ? 16384 : nextpow(2, max(n, 1))
    return _CentroidWindow(zeros(Int, W), Vector{USPoint}(undef, W), W - 1)
end

@inline function _centroid(ws::_CentroidWindow, grid, k::Int, c)
    slot = (k & ws.mask) + 1
    if (@inbounds ws.tags[slot]) == k
        return @inbounds ws.vals[slot]
    end
    return _centroid_miss(ws, grid, c, slot, k)
end

# The miss is deliberately out of line, and it is worth an eighth of what the
# centroids cost the sweep. A hit is a tag compare and a load, and it is what
# six of every seven touches are; inlining `cell_centroid` — a lattice decode
# and a handful of trig calls — beside it bloats the ring loop the sweep
# spends its time in for the sake of the seventh.
@noinline function _centroid_miss(ws::_CentroidWindow, grid, c, slot::Int,
        k::Int)
    p = cell_centroid(grid, c)::USPoint
    @inbounds ws.tags[slot] = k
    @inbounds ws.vals[slot] = p
    return p
end

# No window: every touch computes. This is the path a request that never asks
# for a centroid takes, and it costs it nothing.
@inline _centroid(::Nothing, grid, k::Int, c) = cell_centroid(grid, c)::USPoint

# One window per task, built only when the request names a centroid — several
# `Centroid()`s in one request share it. Recursion over the tuple rather than
# a loop, so the answer is a type the compiler already knows.
_taskcache(::Tuple{}, ::UnitRange{Int}) = nothing
_taskcache(::Tuple{Centroid,Vararg{Any}}, r::UnitRange{Int}) =
    _CentroidWindow(length(r))
_taskcache(needs::Tuple, r::UnitRange{Int}) = _taskcache(Base.tail(needs), r)

# --- element types ---------------------------------------------------------

# The element type one need's ring carries, and the type of its center entry.
_needtype(::Cell, cv::CellVector) = eltype(cv)
_needtype(::Index{Local}, ::CellVector) = Int
_needtype(::Index{Global}, ::CellVector) = Int
_needtype(::Index{Type{T}}, ::CellVector) where {T<:AbstractCellIndex} = T
_needtype(n::Value, ::CellVector) = eltype(n.data)
_needtype(::Centroid, ::CellVector) = USPoint

_centertype(needs::Tuple, cv::CellVector) =
    Tuple{map(n -> _needtype(n, cv), needs)...}

_ringstype(needs::Tuple, cv::CellVector, cap) =
    Tuple{map(n -> _ringtype(cap, _needtype(n, cv)), needs)...}

# --- per-need answers ------------------------------------------------------

# The visited cell, named by its local index `k` and its handle `c`. `ws` is
# the task's centroid working set, or `nothing` where no centroid was asked
# for; every need is handed it so the request can stay a plain `map`.
@inline _center(::Cell, ::CellVector, k::Int, c, ws) = cellid(c)
@inline _center(::Index{Local}, ::CellVector, k::Int, c, ws) = k
@inline _center(::Index{Global}, cv::CellVector, k::Int, c, ws) =
    globalindex(cv.grid, c)::Int
@inline _center(::Index{Type{T}}, cv::CellVector, k::Int, c, ws) where {T} =
    reindex(T, system(cv), c)
@inline _center(n::Value, ::CellVector, k::Int, c, ws) = @inbounds n.data[k]
@inline _center(::Centroid, cv::CellVector, k::Int, c, ws) =
    _centroid(ws, cv.grid, k, c)

# One ring slot, named by the clipped handle the sweep produced. The handle
# carries the neighbour's own local index, which is the key the working set is
# built on, so a neighbour and the visit that later lands on it share a slot.
@inline _slot(::Cell, ::CellVector, h::SubsetIndexedCell, ws) = h.cell
@inline _slot(::Index{Local}, ::CellVector, h::SubsetIndexedCell, ws) = h.index
@inline _slot(::Index{Global}, cv::CellVector, h::SubsetIndexedCell, ws) =
    globalindex(cv.grid, h.cell)::Int
@inline _slot(::Index{Type{T}}, cv::CellVector, h::SubsetIndexedCell,
    ws) where {T} = reindex(T, system(cv), h.cell)
@inline _slot(n::Value, ::CellVector, h::SubsetIndexedCell, ws) =
    @inbounds n.data[h.index]
@inline _slot(::Centroid, cv::CellVector, h::SubsetIndexedCell, ws) =
    _centroid(ws, cv.grid, h.index, h.cell)

# --- rings -----------------------------------------------------------------

# One need's answers for the whole clipped ring, in the container the capacity
# witness names — the same witness `_ringstype` derives the declared type from,
# so the ring the callback receives and the type the output was sized for
# cannot disagree. A declared bound keeps the row stack-allocated.
@inline function _ring(n::AbstractNeed, cv::CellVector, nbrs, cap::Val, ws)
    out = _newring(cap, _needtype(n, cv))
    for h in nbrs
        out = _pushring(out, _slot(n, cv, h, ws))
    end
    return out
end

@inline function _ring(n::AbstractNeed, cv::CellVector, nbrs, ::Nothing, ws)
    out = Vector{_needtype(n, cv)}(undef, length(nbrs))
    for (i, h) in enumerate(nbrs)
        @inbounds out[i] = _slot(n, cv, h, ws)
    end
    return out
end

# The two tuples the callback receives. `map` over the request, never a loop:
# the entries have different types and the tuple must stay inferable.
@inline _centers(needs::Tuple, cv::CellVector, k::Int, c, ws) =
    map(n -> _center(n, cv, k, c, ws), needs)

@inline _rings(needs::Tuple, cv::CellVector, nbrs, cap, ws) =
    map(n -> _ring(n, cv, nbrs, cap, ws), needs)

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

_checkneed(n::Value, cv::CellVector) = _check_data(n.data, length(cv))
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
    _checkneeds(needs, cv)
    cap = _capacity(system(cv), connectivity)
    T = Base.promote_op(f, _centertype(needs, cv), _ringstype(needs, cv, cap))
    outs = _outputs(T, length(cv))
    return _mapstore!(f, outs, cv, needs, connectivity, order,
        GOCore.booltype(threaded), cap)
end

# Built after the output container type is known, like the other `_mapstore!`s.
# The callback is built once per task rather than once per sweep, so the
# working set it closes over belongs to that task alone.
function _mapstore!(f::F, outs::O, cv::CellVector, needs::Tuple,
        conn::Connectivity, order, thr, cap::CAP) where {F,O,CAP}
    _runeach!(cv, conn, order, thr, cap) do r
        ws = _taskcache(needs, r)
        (k, c, nbrs) -> _store!(outs, k,
            f(_centers(needs, cv, k, c, ws), _rings(needs, cv, nbrs, cap, ws)))
    end
    return outs
end

function _foreachneighbors(f::F, cv::CellVector, needs, order, threaded,
        connectivity::Connectivity) where {F}
    _checkneeds(needs, cv)
    # One capacity witness for the whole sweep: the clip, the rings it feeds
    # and the ring type all come from this value.
    cap = _capacity(system(cv), connectivity)
    _runeach!(cv, connectivity, order, GOCore.booltype(threaded), cap) do r
        ws = _taskcache(needs, r)
        (k, c, nbrs) -> (f(_centers(needs, cv, k, c, ws),
            _rings(needs, cv, nbrs, cap, ws)); nothing)
    end
    return nothing
end
