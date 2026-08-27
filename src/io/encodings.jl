# The two halves of a stored cell axis.
#
#   * A GRID owns id ARITHMETIC: rank, select, and the count in an id interval,
#     over the complete level in canonical order. Three functions per grid.
#   * An ENCODING owns the LAYOUT the store keeps that axis in: dense ids, an
#     inclusive-range array, or nothing at all. One set of verbs per encoding.
#
# Neither knows the other, so a new grid gains every encoding by implementing
# the arithmetic, and a new encoding works on every grid that has it. The verbs
# that need the axis and manifest TYPES are declared here and implemented in
# `chunked_lookup.jl`, which defines them.

"""
    Encodings

Grid id arithmetic and the pluggable cell-axis encodings a DGGS store can be
written in.

[`idrank`](@ref), [`idselect`](@ref) and [`idcount_between`](@ref) are the grid
side: closed-form rank/select over one complete level, defined on every
`UInt64`/`Int64` rather than only on ids that name cells, so a reader can count
what a stored interval holds without enumerating it.

[`CellEncoding`](@ref) and its four singletons are the layout side, resolved
from a store's vocabulary through [`ENCODING_REGISTRY`](@ref).
"""
module Encodings

import ..DiscreteGlobalGrids: AbstractGrid, AbstractCellIndex,
    AbstractHierarchicalGridSystem, HierarchicalLevelGrid,
    IGeo7System, HEALPixSystem, Z7Cell, LevelIndex,
    ncells, cellindex, level, levels, levelgrid, rawid
# The compacted layout's axis IS the mixed-level container; nothing else here
# touches a collection type.
import ..DiscreteGlobalGrids.Engine: MultiOrderVector
# The Z7 digit tables, validity predicate and inverse-rank walk are the system's
# own; this file adds only the rank that is total on ids naming no cell.
import ..DiscreteGlobalGrids: IGeo7
# From `errors.jl`, which the including module reads before this file.
import ..DGGSFormatError

# WIRING: rejections here throw `DGGSFormatError` with no store context, because
# this layer sees ids and lengths and never a store. The extension that opened
# the store adds the identifier and the conventions that fired on the way out,
# with `with_store_context`.

export CellEncoding, DenseEncoding, RangesEncoding, ImplicitEncoding,
    CompactedEncoding
export ENCODING_REGISTRY, encodingname, register_encoding!
export idrank, idselect, idcount_between, idvalid, idcell, idtype, idlevel
export idranges, write_eligible, validate_ranges, cellaxis, storedid

# ===========================================================================
# The grid side: id arithmetic over one complete level
# ===========================================================================

"""
    idrank(grid::AbstractGrid, id::Integer) -> Int

The number of cells of `grid` whose raw id is strictly less than `id` — a
COUNT, so it is a zero-based rank and `idrank(grid, rawid(c)) + 1` is
[`globalindex`](@ref)`(grid, c)`.

**Total on the integer type.** `id` need not name a cell: an id above every
cell of the level answers `ncells(grid)`, one below every cell answers `0`, and
one that is well formed but names nothing — a Z7 phantom on a pentagon's
deleted branch — answers where it would sit. That totality is the whole point:
a stored `[start, stop]` interval is counted by subtracting two ranks, and an
interval's endpoints are not required to be cells.

`grid` is a complete level grid; rank is meaningless against a subset, which
has [`localindex`](@ref) instead.

**Required** of a grid that is to be read from a store, together with
[`idselect`](@ref) and [`idcount_between`](@ref).
"""
function idrank end

"""
    idselect(grid::AbstractGrid, r::Integer) -> Integer

The raw id of the cell at zero-based rank `r` in `grid`'s canonical order — the
inverse of [`idrank`](@ref):

    idrank(grid, idselect(grid, r)) == r          for r in 0:ncells(grid)-1
    idselect(grid, idrank(grid, x)) == x          for every cell id x

Equivalently `rawid(cellindex(grid, r + 1))`, without constructing the typed id.
`r` outside `0:ncells(grid)-1` throws a `BoundsError`.
"""
function idselect end

"""
    idcount_between(grid::AbstractGrid, lo::Integer, hi::Integer) -> Int

The number of cells of `grid` whose raw id lies in the INCLUSIVE interval
`[lo, hi]`; zero when `hi < lo`. Neither endpoint need name a cell.

This is what makes a `ranges` store readable without touching a data array: the
axis length is the sum of this over the stored intervals, and the same sum
prefixed gives the index of every interval's first cell.

The generic implementation is the rank difference, and is correct for any grid
that implements [`idrank`](@ref).
"""
function idcount_between end

"""
    idvalid(grid::AbstractGrid, id::Integer) -> Bool

Whether `id` names a cell **of this grid's level**. Never throws.

The ingest-time check for a dense axis: a Z7 id can be well formed and still
name nothing (the twelve pentagons delete one child digit each), and DGGRID
itself does not reject those, so a reader must.
"""
function idvalid end

"""
    idlevel(grid::AbstractGrid, id::Integer) -> Union{Int,Nothing}

The level `id` names a cell of, read out of the id ITSELF, or `nothing` where
the id scheme carries no level to read.

**Optional**, and informational only: a stored axis is read at the level the
store declares, and this never overrides it. What it is for is the error a
level-lying store produces — "names no cell of level 3" is more useful with
"and is a cell of level 4" after it, and that second half is the id's own claim,
not a second opinion about the store's.

Z7 pads its digit slots to a fixed width, so the level is the count of leading
non-pad digits and every well-formed id answers. HEALPix nested ids carry
nothing: `0:12*4^L-1` is a prefix of `0:12*4^(L+1)-1`, so a level-2 id is a
structurally valid level-3 id and the honest answer is `nothing` rather than a
guess. A grid that does not implement this gets the default, which is the same
answer.
"""
idlevel(::AbstractGrid, ::Integer) = nothing

"""
    idcell(grid::AbstractGrid, id::Integer) -> AbstractCellIndex

The typed cell id for a raw stored integer at this grid's level. The inverse of
`rawid` given the level, which the raw integer does not always carry.
"""
function idcell end

"""
    idtype(grid::AbstractGrid) -> Type{<:Integer}

The integer type this grid's raw ids are stored in. The generic implementation
asks the first cell, which is exact but not free; grids with a fixed width
should say so directly.
"""
idtype(grid::AbstractGrid) = typeof(rawid(cellindex(grid, 1)))

"""
    storedid(::Type{I}, x) -> I

A value out of a store's cell coordinate as the grid's own id type `I`.

**A stored id is REINTERPRETED, not converted.** Integer width and signedness
are a writer's choice and no part of a cell's identity: an IGEO7 id in base cell
8 or above sets the top bit, so a store whose coordinate is `Int64` — which is
what a producer gets by letting the ids through any signed dtype — holds those
cells as negative integers whose bits are exactly the same cells. A signed value
of `I`'s own width is therefore read as those bits.

Nothing downstream takes that on trust: a bit pattern naming no cell is what
[`idvalid`](@ref) rejects in the scan, and this function decides only whether
the bits fit. Anything that fits neither the reinterpretation nor a plain
conversion — a fractional float, a value of another width — is a malformed
coordinate and raises `DGGSFormatError(check = :coordinate_width)` naming the
value, rather than the `InexactError` a bare `convert` raises with no store in
it.

A persisted chunk manifest is deliberately stricter: one written at another
width is DECLINED rather than reinterpreted, because declining it costs a scan
and nothing else, where declining the coordinate would make the store
unreadable.
"""
function storedid(::Type{I}, x::Integer) where {I<:Integer}
    typemin(I) <= x <= typemax(I) && return x % I
    # The same width read with the other signedness: the bits are the id.
    (x isa Signed && sizeof(x) == sizeof(I)) && return x % I
    return _coordinate_width(I, x)
end

# A coordinate stored as something other than an integer. Zarr keeps whatever
# dtype the writer chose, and a float array is an id array only where each of
# its values is exactly one id.
function storedid(::Type{I}, x) where {I<:Integer}
    (isinteger(x) && typemin(I) <= x <= typemax(I)) && return convert(I, x)
    return _coordinate_width(I, x)
end

@noinline _coordinate_width(::Type{I}, x) where {I<:Integer} =
    throw(DGGSFormatError(check=:coordinate_width, declared=I, observed=x,
        detail="the stored value $x is no cell id of type $I. A cell coordinate " *
               "holds the grid's own ids: at $I, or at the same width read with " *
               "the other signedness, which is read as the same bits."))

@noinline function _no_arithmetic(grid)
    name = string(nameof(typeof(grid)))
    throw(DGGSFormatError(check=:no_id_arithmetic, observed=name,
        detail="$name has no id arithmetic: reading a DGGS store needs " *
               "`idrank`, `idselect` and `idcount_between` on its complete " *
               "level grid."))
end

idrank(grid::AbstractGrid, ::Integer) = _no_arithmetic(grid)
idselect(grid::AbstractGrid, ::Integer) = _no_arithmetic(grid)
idcell(grid::AbstractGrid, ::Integer) = _no_arithmetic(grid)

function idcount_between(grid::AbstractGrid, lo::Integer, hi::Integer)
    I = idtype(grid)
    a = convert(I, lo)
    b = convert(I, hi)
    b < a && return 0
    # `b + 1` would wrap at the top of the integer type; every cell is <= b there.
    upper = b == typemax(I) ? ncells(grid) : idrank(grid, b + one(I))
    return upper - idrank(grid, a)
end

# --- IGEO7 / Z7: the existence model ---------------------------------------
#
# A level-L Z7 id is a base cell and L base-7 digits, tail-padded with the 7
# sentinel, so ascending integers are ascending mixed-radix numerals and rank is
# a digit walk. The correction is the pentagon chain: at every all-zero-prefix
# node one child digit is deleted (2 for bases 0:5, 5 for bases 6:11), which
# leaves `10*7^L + 2` cells rather than `12*7^L` well-formed strings. Both walks
# below carry `azero`, the "prefix so far is all zeros" flag that is exactly the
# condition under which the deletion applies.

const _Z7Grid = HierarchicalLevelGrid{IGeo7System}

idtype(::_Z7Grid) = UInt64
idcell(::_Z7Grid, id::Integer) = Z7Cell(UInt64(id))

function idrank(grid::_Z7Grid, id::Integer)
    L = grid.level
    x = UInt64(id)
    b = Int(x >> 60)
    # Above every cell: the four high bits name a base cell that does not exist.
    b >= IGeo7.Z7_NUM_BASES && return Int(IGeo7.num_cells(L))
    K = @inbounds IGeo7.Z7_DELETED_DIGIT[b+1]
    cnt = b * (@inbounds IGeo7.PENT_COUNT[L+1])
    azero = true
    @inbounds for j in 1:L
        e = IGeo7._z7_digit(x, j)
        w = IGeo7.POW7[L-j+1]
        if azero
            if e > 0
                # digit 0 keeps the pentagon chain and precedes every other
                # digit; the hexagon siblings below `e` are 1:(e-1) less `K`.
                cnt += IGeo7.PENT_COUNT[L-j+1] + (e - 1 - (K < e ? 1 : 0)) * w
            end
        else
            cnt += e * w
        end
        # A pad digit here, or the deleted digit on the chain, means `x` is not
        # a cell and every cell sharing this prefix has already been counted.
        (e <= 6 && !(azero && e == K)) || return Int(cnt)
        azero &= e == 0
    end
    # Exact prefix match: the one cell with this prefix has an all-7 tail, which
    # is maximal, so it is >= x and contributes nothing.
    return Int(cnt)
end

# The package's own inverse-rank walk is this level's select, one-based.
idselect(grid::_Z7Grid, r::Integer) = IGeo7.index_to_cell(Int(r) + 1, grid.level)

idvalid(grid::_Z7Grid, id::Integer) =
    id >= 0 && IGeo7.is_valid_z7(UInt64(id)) &&
    IGeo7._z7_leading_resolution(UInt64(id)) == grid.level

# The leading-resolution walk is the level: the padding tail says where the
# digits stop. A structurally invalid id has none to report — a broken padding
# tail or a deleted-digit prefix names no cell at ANY level, so there is nothing
# to tell the reader about it.
function idlevel(::_Z7Grid, id::Integer)
    id >= 0 && IGeo7.is_valid_z7(UInt64(id)) || return nothing
    return Int(IGeo7._z7_leading_resolution(UInt64(id)))
end

# --- HEALPix nested: the contiguous model ----------------------------------

const _HPXGrid = HierarchicalLevelGrid{HEALPixSystem}

idtype(::_HPXGrid) = Int64
idcell(grid::_HPXGrid, id::Integer) = LevelIndex(grid.level, Int64(id))

# Nested ids are exactly `0:12*4^L-1` in canonical order, so rank is the id
# itself, clamped to the level at both ends.
function idrank(grid::_HPXGrid, id::Integer)
    n = ncells(grid)
    id <= 0 && return 0
    return id >= n ? n : Int(id)
end

function idselect(grid::_HPXGrid, r::Integer)
    n = ncells(grid)
    # The rank as given: `idselect` is zero-based, and reporting `r + 1` would
    # name an index the caller never passed.
    0 <= r < n || throw(BoundsError(grid, Int(r)))
    return Int64(r)
end

idvalid(grid::_HPXGrid, id::Integer) = 0 <= id < ncells(grid)

# No `idlevel`: the default `nothing` is the correct answer here and not a gap.
# Every level-L nested id is also a well-formed level-(L+1) id, so an id names
# no level of its own and a reader cannot be told one.

# ===========================================================================
# The encoding side
# ===========================================================================

"""
    CellEncoding

How a store lays its cell axis out. The four shipped layouts are
[`DenseEncoding`](@ref), [`RangesEncoding`](@ref), [`ImplicitEncoding`](@ref)
and [`CompactedEncoding`](@ref); a downstream package adds its own by subtyping
this and registering an instance in [`ENCODING_REGISTRY`](@ref).

An encoding implements [`cellaxis`](@ref) (build the axis, and with it the
chunk manifest), its own validation, and [`write_eligible`](@ref), which
`encoding = :auto` consults. It never touches id arithmetic: everything it
needs about the ids themselves is [`idrank`](@ref) / [`idselect`](@ref) /
[`idcount_between`](@ref) on the grid.
"""
abstract type CellEncoding end

"""
    DenseEncoding()

One stored id per cell. Universal and the interop escape hatch; the axis costs
a chunked pass over the id array to verify and to build the manifest from.
"""
struct DenseEncoding <: CellEncoding end

"""
    RangesEncoding()

An `(n, 2)` array of INCLUSIVE `[start, stop]` raw-id intervals at one level.
The axis is computed from the intervals by rank/select, so it needs no id
storage at all and no data IO: the length, every chunk's first and last id, and
every selector are closed-form.
"""
struct RangesEncoding <: CellEncoding end

"""
    ImplicitEncoding()

No stored axis: index `k` is the cell at rank `k - 1` of the level. The
whole-level case, as written by the DKRZ-style conventions.
"""
struct ImplicitEncoding <: CellEncoding end

"""
    CompactedEncoding()

One stored cell per position at MIXED refinement levels: two aligned columns,
a level and a raw id at that level, in [`MultiOrderVector`](@ref) container
order — ascending by subtree-interval start. Such a store declares
`refinement_level: null` and `compression: "compacted"`, names its level
column in `refinement_levels`, and its axis reads back as a
`MultiOrderVector` rather than a single-level vector.

This is an EXTENSION, not a conforming layout: v1 of `zarr-conventions/dggs`
takes the vocabulary word `compacted` but requires `compression: "none"`
wherever `refinement_level` is null, and for healpix it requires a `*uniq`
indexing scheme there — a self-describing id with no level column at all. A
conforming reader will reject these stores; only this package reads them.
"""
struct CompactedEncoding <: CellEncoding end

"""
    ENCODING_REGISTRY

Store vocabulary (`"none"`, `"ranges"`, `"implicit"`, `"compacted"`) to
[`CellEncoding`](@ref) instance. Conventions resolve an attribute's string
through this table, and keyword symbols are sugar over the same entries, so a
downstream encoding becomes usable by adding one pair —
[`register_encoding!`](@ref).
"""
const ENCODING_REGISTRY = Dict{String,CellEncoding}(
    "none" => DenseEncoding(),
    "ranges" => RangesEncoding(),
    "implicit" => ImplicitEncoding(),
    "compacted" => CompactedEncoding(),
)

"""
    register_encoding!(name::AbstractString, enc::CellEncoding) -> ENCODING_REGISTRY

Register `enc` under the vocabulary string `name`, the way a store spells it.

A store's `compression` attribute, and `dggwrite`'s `encoding` keyword, are both
resolved through [`ENCODING_REGISTRY`](@ref), so this is what makes a downstream
encoding reachable by name. Reading it also needs [`cellaxis`](@ref) and writing
it the write path's own verbs; an encoding that registers without them is
refused by name rather than by `MethodError`.
"""
register_encoding!(name::AbstractString, enc::CellEncoding) =
    setindex!(ENCODING_REGISTRY, enc, String(name))

"""
    encodingname(enc::CellEncoding) -> String

The [`ENCODING_REGISTRY`](@ref) key `enc` is registered under — what a writer
stamps into the store's attributes.
"""
encodingname(::DenseEncoding) = "none"
encodingname(::RangesEncoding) = "ranges"
encodingname(::ImplicitEncoding) = "implicit"
encodingname(::CompactedEncoding) = "compacted"

"""
    cellaxis(enc::CellEncoding, grid_or_system, source...; kw...) -> AbstractVector

The stored axis as an `AbstractVector` of typed cell ids. `source` is whatever
the encoding reads:

| encoding | `source` | cost |
|---|---|---|
| [`RangesEncoding`](@ref) | the `(n, 2)` inclusive-range array | arithmetic, no IO |
| [`ImplicitEncoding`](@ref) | the axis length | arithmetic, no IO |
| [`DenseEncoding`](@ref) | the id vector, lazy or not | one chunked pass |
| [`CompactedEncoding`](@ref) | the level and id columns | one validating pass |

The single-level encodings take a level grid and answer a `ChunkedCellVector`
(built in `chunked_lookup.jl`, where that type lives); the compacted one takes
the SYSTEM and answers a [`MultiOrderVector`](@ref).
"""
function cellaxis end

# --- write-side eligibility and range construction -------------------------

"""
    write_eligible(enc::CellEncoding, grid::AbstractGrid, ids::AbstractVector) -> Bool

Whether `ids` can be written in `enc`. This is what `encoding = :auto` asks:
[`RangesEncoding`](@ref) is eligible exactly when the ids are sorted, unique,
and all cells of `grid`'s single level; [`ImplicitEncoding`](@ref) additionally
requires them to be the whole level; [`DenseEncoding`](@ref) always is, which
is what makes `:auto` total.

`ids` may be raw integers or typed cell ids. A downstream encoding that states
no restriction is eligible: the fallback is `true`, and the write path refuses
it by name if it implements nothing else.
"""
write_eligible(::CellEncoding, ::AbstractGrid, ::AbstractVector) = true

write_eligible(::DenseEncoding, ::AbstractGrid, ::AbstractVector) = true

write_eligible(::RangesEncoding, grid::AbstractGrid, ids::AbstractVector) =
    _sorted_unique_cells(grid, ids)

# A single-level axis is written dense or as ranges; the compacted layout is
# reached from a `MultiOrderLookup`, never from one level's ids.
write_eligible(::CompactedEncoding, ::AbstractGrid, ::AbstractVector) = false

function write_eligible(::ImplicitEncoding, grid::AbstractGrid, ids::AbstractVector)
    length(ids) == ncells(grid) || return false
    _sorted_unique_cells(grid, ids) || return false
    return _rawid(grid, first(ids)) == idselect(grid, 0)
end

function _sorted_unique_cells(grid::AbstractGrid, ids::AbstractVector)
    previous = nothing
    for c in ids
        x = _rawid(grid, c)
        (x === nothing || !idvalid(grid, x)) && return false
        previous === nothing || x > previous || return false
        previous = x
    end
    return true
end

# A typed id from another level is not a cell of this grid, and says so by
# answering `nothing` rather than by handing over a raw id that would validate.
_rawid(grid::AbstractGrid, x::Integer) = convert(idtype(grid), x)
_rawid(grid::AbstractGrid, c::AbstractCellIndex) =
    level(c) == level(grid) ? convert(idtype(grid), rawid(c)) : nothing

"""
    idranges(grid::AbstractGrid, ids::AbstractVector; merge = :rank) -> Matrix

The `(n, 2)` inclusive-range array for a sorted, unique, single-level `ids`, one
row per run. `merge` chooses what counts as a run:

  - `:rank` — maximal runs of CONSECUTIVE CELLS. An interval is read back by
    counting the cells inside it, not by stepping through the integers, so a run
    may span integers that name nothing: at a Z7 digit rollover it always does,
    since a 3-bit slot holds eight values and a base-7 digit uses seven, and
    across a pentagon's deleted branch it may.
  - `:step` — runs of ids ADJACENT AS INTEGERS, at the level's unit increment.
    No interval can then span an id that names nothing, so a reader that counts
    well-formed digit strings rather than cells gets the same answer. It costs
    rows: the whole res-3 earth is 504 runs of at most 7 cells this way, and one
    single run the other.

Both encode the same axis, and [`cellaxis`](@ref) reads either back to it. The
difference is interop: `:step` is what a structural-count reader needs, and what
the published IGEO7 range stores hold.

Throws a `DGGSFormatError` when `ids` is not [`write_eligible`](@ref) for
[`RangesEncoding`](@ref).

!!! note "Runs across a pentagon's deleted branch"
    Under `:rank` a structural-count reader overcounts an interval by the
    phantoms inside it — 4 116 against 3 432 on the whole res-3 level. The three
    published IGEO7 range stores contain no such interval, because they lie off
    the pentagon chains entirely.
"""
function idranges(grid::AbstractGrid, ids::AbstractVector; merge::Symbol=:rank)
    merge in (:rank, :step) || throw(ArgumentError(
        "the ranges merge rule is `:rank` or `:step`, not $(repr(merge))"))
    write_eligible(RangesEncoding(), grid, ids) || throw(DGGSFormatError(
        check=:not_range_eligible, declared=level(grid),
        detail="ranges encoding needs sorted, unique cell ids from level $(level(grid))."))
    I = idtype(grid)
    isempty(ids) && return Matrix{I}(undef, 0, 2)
    unit = merge === :step ? _idunit(grid) : nothing
    starts = I[]
    stops = I[]
    previous = _rawid(grid, first(ids))
    rank = idrank(grid, previous)
    push!(starts, previous)
    for k in 2:length(ids)
        x = _rawid(grid, ids[k])
        r = idrank(grid, x)
        run = r == rank + 1 && (unit === nothing || x - previous == unit)
        run || (push!(stops, previous); push!(starts, x))
        previous = x
        rank = r
    end
    push!(stops, previous)
    return [starts stops]
end

# The id distance between two cells that are adjacent in canonical order with
# nothing rolling over between them — the level's unit increment, `2^(3(20-L))`
# for Z7 and 1 for HEALPix nested. Read off the first two cells of the level,
# which every grid with a hierarchical id orders as siblings.
function _idunit(grid::AbstractGrid)
    ncells(grid) >= 2 || return one(idtype(grid))
    return idselect(grid, 1) - idselect(grid, 0)
end

"""
    rangeindex(grid::AbstractGrid, ranges::AbstractMatrix)
        -> (starts, stops, rstart, offsets)

Validate an `(n, 2)` inclusive-range array and reduce it to the rank/select
dictionary the axis is read through: `rstart[i]` is the global rank of interval
`i`'s first cell, and `offsets[i]` the number of cells before it, with
`offsets[end]` the total.

Each row must be a non-empty interval and the rows must ascend and be disjoint;
anything else is a malformed store and throws a `DGGSFormatError` naming the row.
The whole thing is `O(n * level)` digit arithmetic over an array that is
kilobytes even for a store of tens of millions of cells, and it reads no data.
"""
function rangeindex(grid::AbstractGrid, ranges::AbstractMatrix)
    size(ranges, 2) == 2 || throw(DGGSFormatError(check=:invalid_ranges_shape,
        observed=size(ranges),
        detail="a ranges array is (n, 2) inclusive [start, stop] pairs, got " *
               "$(size(ranges, 1))x$(size(ranges, 2))."))
    I = idtype(grid)
    n = size(ranges, 1)
    starts = I[storedid(I, ranges[i, 1]) for i in 1:n]
    stops = I[storedid(I, ranges[i, 2]) for i in 1:n]
    rstart = Vector{Int}(undef, n)
    offsets = Vector{Int}(undef, n + 1)
    offsets[1] = 0
    for i in 1:n
        stops[i] >= starts[i] || throw(DGGSFormatError(check=:empty_range_row,
            declared=starts[i], observed=stops[i],
            detail="ranges row $i is empty: stop $(stops[i]) is below start $(starts[i])."))
        i == 1 || starts[i] > stops[i-1] || throw(DGGSFormatError(
            check=:overlapping_ranges, declared=stops[i-1], observed=starts[i],
            detail="ranges row $i starts at $(starts[i]), which does not follow " *
                   "the previous row's stop $(stops[i-1]); rows must ascend and " *
                   "be disjoint."))
        rstart[i] = idrank(grid, starts[i])
        offsets[i+1] = offsets[i] + idcount_between(grid, starts[i], stops[i])
    end
    return starts, stops, rstart, offsets
end

"""
    validate_ranges(grid::AbstractGrid, ranges::AbstractMatrix, declared::Integer)

Check an `(n, 2)` inclusive-range array against the length the store declares
for its data arrays, and throw a `DGGSFormatError` naming the failure otherwise.

[`rangeindex`](@ref)'s structural checks, and then the normative one: the total
cell count is `declared`. That last comparison is what catches an
expansion-semantics disagreement — an exclusive `stop`, a wrong step, a
structural rather than existence count — before a single data byte is read.
"""
function validate_ranges(grid::AbstractGrid, ranges::AbstractMatrix, declared::Integer)
    offsets = last(rangeindex(grid, ranges))
    checkcount(grid, offsets[end], declared)
    return nothing
end

"""
    checkcount(grid::AbstractGrid, total::Integer, declared::Integer)

Throw unless a closed-form cell count matches the length the store declares.
"""
function checkcount(grid::AbstractGrid, total::Integer, declared::Integer)
    total == declared || throw(DGGSFormatError(check=:count_mismatch,
        declared=Int(declared), observed=Int(total),
        detail="the cell axis holds $total cells of level $(level(grid)) but " *
               "the store declares $declared; the axis and the data disagree."))
    return nothing
end

# ===========================================================================
# The compacted axis
# ===========================================================================

"""
    cellaxis(CompactedEncoding(), sys::AbstractHierarchicalGridSystem,
             cell_levels::AbstractVector{<:Integer}, cell_ids::AbstractVector;
             declared_length = nothing) -> MultiOrderVector

Build the mixed-level axis of a `compression: "compacted"` store from its two
aligned columns: position `k` holds the cell whose raw id is `cell_ids[k]` at
level `cell_levels[k]`.

Every pair is checked — the level must be one of the system's, the id must name
a cell of that level — and the whole set must already BE a
[`MultiOrderVector`](@ref): pairwise disjoint as subtrees and ascending by
subtree-interval start. Container order is required rather than restored,
because the store's data arrays are laid out against these positions and a
sort here would silently misalign every value.
"""
function cellaxis(::CompactedEncoding, sys::AbstractHierarchicalGridSystem,
    cell_levels::AbstractVector{<:Integer}, cell_ids::AbstractVector;
    declared_length::Union{Integer,Nothing}=nothing)
    n = length(cell_ids)
    length(cell_levels) == n || throw(DGGSFormatError(
        check=:compacted_column_mismatch, declared=n, observed=length(cell_levels),
        detail="a compacted axis holds one level per id: $(length(cell_levels)) " *
               "levels against $n ids."))
    declared_length === nothing || declared_length == n || throw(DGGSFormatError(
        check=:count_mismatch, declared=Int(declared_length), observed=n,
        detail="the cell axis holds $n cells but the store declares " *
               "$declared_length; the axis and the data disagree."))
    grids = Dict{Int,Any}()
    cells = map(1:n) do k
        # Membership before narrowing: `Int` of a wild unsigned is an
        # InexactError, not the format error a malformed store deserves.
        raw = cell_levels[k]
        raw in levels(sys) || throw(DGGSFormatError(check=:invalid_stored_level,
            declared=Int(first(levels(sys))):Int(last(levels(sys))), observed=raw,
            detail="position $k of the cell axis declares level $raw, which " *
                   "$(nameof(typeof(sys))) does not have."))
        l = Int(raw)
        grid = get!(() -> levelgrid(sys, l), grids, l)
        id = storedid(idtype(grid), cell_ids[k])
        idvalid(grid, id) || throw(DGGSFormatError(check=:id_names_no_cell,
            declared=l, observed=id,
            detail="the id $id at position $k of the cell axis names no cell " *
                   "of level $l."))
        idcell(grid, id)
    end
    mov = try
        MultiOrderVector(sys, cells)
    catch err
        err isa ArgumentError || rethrow()
        throw(DGGSFormatError(check=:invalid_compacted_axis,
            observed=nameof(typeof(sys)),
            detail="the stored cells do not form a multi-order axis: " * err.msg))
    end
    all(k -> mov[k] == cells[k], 1:n) || throw(DGGSFormatError(
        check=:compacted_axis_order,
        detail="the compacted cell axis is not in container order — ascending " *
               "by subtree-interval start — so the store's data arrays do not " *
               "align with its cells. Rewrite the store with `dggwrite`."))
    return mov
end

end # module Encodings
