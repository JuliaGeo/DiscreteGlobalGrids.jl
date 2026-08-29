# Separate grid-arithmetic and storage-layout interfaces make every combination reusable.

"""
    Encodings

Grid id arithmetic and pluggable DGGS cell-axis encodings.

[`idrank`](@ref), [`idselect`](@ref) and [`idcount_between`](@ref) provide
closed-form arithmetic over one complete level. Their integer-domain definitions
let readers count stored intervals without enumerating their cells.

[`CellEncoding`](@ref) and its four shipped layouts form the storage side,
resolved from a store's vocabulary through [`ENCODING_REGISTRY`](@ref).
"""
module Encodings

import ..DiscreteGlobalGrids: AbstractGrid, AbstractCellIndex,
    AbstractHierarchicalGridSystem, HierarchicalLevelGrid,
    IGeo7System, HEALPixSystem, Z7Cell, LevelIndex,
    ncells, cellindex, level, levels, levelgrid, rawid
import ..DiscreteGlobalGrids.Engine: MultiOrderVector
# IGeo7 owns the Z7 tables and validity walk; storage adds total rank arithmetic.
import ..DiscreteGlobalGrids: IGeo7
# From `errors.jl`, which the including module reads before this file.
import ..DGGSFormatError

# The store boundary adds identifiers and convention names to these layer-neutral errors.

export CellEncoding, DenseEncoding, RangesEncoding, ImplicitEncoding,
    CompactedEncoding
export ENCODING_REGISTRY, encodingname, register_encoding!
export idrank, idselect, idcount_between, idvalid, idcell, idtype, idlevel
export idranges, write_eligible, validate_ranges, cellaxis, storedid

"""
    idrank(grid::AbstractGrid, id::Integer) -> Int

Return the number of `grid` cells whose raw id is less than `id`. This zero-based
rank satisfies `idrank(grid, rawid(c)) + 1 == globalindex(grid, c)`.

The function accepts every value of the grid's integer type:

  - values below the level return `0`;
  - values above the level return `ncells(grid)`; and
  - non-cell values return their insertion rank.

This total definition counts stored `[start, stop]` intervals by rank difference,
including intervals whose endpoints are not cells.

`grid` must be a complete level grid. Subsets use [`localindex`](@ref).

**Required** of a grid that is to be read from a store, together with
[`idselect`](@ref) and [`idcount_between`](@ref).
"""
function idrank end

"""
    idselect(grid::AbstractGrid, r::Integer) -> Integer

Return the raw id at zero-based rank `r` in `grid`'s canonical order. This is the
inverse of [`idrank`](@ref):

    idrank(grid, idselect(grid, r)) == r          for r in 0:ncells(grid)-1
    idselect(grid, idrank(grid, x)) == x          for every cell id x

The result equals `rawid(cellindex(grid, r + 1))` while avoiding a typed-cell
allocation. Ranks outside `0:ncells(grid)-1` throw `BoundsError`.
"""
function idselect end

"""
    idcount_between(grid::AbstractGrid, lo::Integer, hi::Integer) -> Int

Return the number of `grid` cells in the inclusive raw-id interval `[lo, hi]`.
The result is zero when `hi < lo`, and either endpoint may be a non-cell value.

Ranges stores use these counts to derive the axis length and each interval's
starting index without reading a data array.

The generic implementation is the rank difference, and is correct for any grid
that implements [`idrank`](@ref).
"""
function idcount_between end

"""
    idvalid(grid::AbstractGrid, id::Integer) -> Bool

Return whether `id` names a cell at `grid`'s level. This predicate never throws.

Dense-axis ingestion uses this check to reject Z7 phantom ids from the deleted
pentagon branches.
"""
function idvalid end

"""
    idlevel(grid::AbstractGrid, id::Integer) -> Union{Int,Nothing}

Return the level encoded by `id`, or `nothing` when the id scheme carries no
level. This optional information enriches validation errors; the store's declared
level remains authoritative.

Z7 derives the level from the leading non-padding digits. HEALPix nested ids
return `nothing` because the id ranges of successive levels overlap. Grids
without a specialized method use the same `nothing` default.
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

Interpret a stored coordinate value as the grid's id type `I`.

Integer signedness is a storage choice. A same-width signed value therefore
preserves its bit pattern when read as `I`; this supports IGEO7 ids whose top bit
appears negative in `Int64`. [`idvalid`](@ref) separately verifies that the
resulting bits name a cell.

Exact numeric values within `I`'s range convert normally. Other widths and
fractional values raise `DGGSFormatError(check = :coordinate_width)` with the
offending value.

Persisted manifests use stricter width matching because the reader can safely
discard an incompatible manifest and scan the coordinate.
"""
function storedid(::Type{I}, x::Integer) where {I<:Integer}
    typemin(I) <= x <= typemax(I) && return x % I
    # The same width read with the other signedness: the bits are the id.
    (x isa Signed && sizeof(x) == sizeof(I)) && return x % I
    return _coordinate_width(I, x)
end

# Non-integer coordinate dtypes are valid only for exact in-range integer values.
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

# A Z7 id contains a base cell, `L` base-7 digits and a 7-filled tail. Integer
# order therefore follows the mixed-radix digits. Pentagon chains delete digit 2
# for bases 0:5 and digit 5 for bases 6:11 whenever the prefix remains all zero;
# `azero` tracks that condition during the rank walk.

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

# The padding boundary gives a valid Z7 id's level; invalid digit paths return nothing.
function idlevel(::_Z7Grid, id::Integer)
    id >= 0 && IGeo7.is_valid_z7(UInt64(id)) || return nothing
    return Int(IGeo7._z7_leading_resolution(UInt64(id)))
end

# --- HEALPix nested: the contiguous model ----------------------------------

const _HPXGrid = HierarchicalLevelGrid{HEALPixSystem}

idtype(::_HPXGrid) = Int64
idcell(grid::_HPXGrid, id::Integer) = LevelIndex(grid.level, Int64(id))

# HEALPix nested rank is the id clamped to its level's `0:12*4^L-1` range.
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

# Overlapping HEALPix level ranges make the default `idlevel == nothing` exact.

"""
    CellEncoding

Describe how a store lays out its cell axis. The shipped layouts are
[`DenseEncoding`](@ref), [`RangesEncoding`](@ref), [`ImplicitEncoding`](@ref)
and [`CompactedEncoding`](@ref). A downstream package can subtype this type and
register an instance in [`ENCODING_REGISTRY`](@ref).

An encoding implements [`cellaxis`](@ref), its validation, and
[`write_eligible`](@ref), which `encoding = :auto` consults. Single-level
encodings build a chunk-aware axis; the compacted encoding builds a
[`MultiOrderVector`](@ref). Encodings delegate id arithmetic to
[`idrank`](@ref), [`idselect`](@ref) and [`idcount_between`](@ref) on the grid.
"""
abstract type CellEncoding end

"""
    DenseEncoding()

Store one id per cell. Reading verifies the id array in a chunked pass and
builds its manifest. This universal layout provides broad interoperability.
"""
struct DenseEncoding <: CellEncoding end

"""
    RangesEncoding()

Store an `(n, 2)` array of inclusive `[start, stop]` raw-id intervals at one
level. Rank/select arithmetic derives the axis length, chunk endpoints and
selectors without reading an id array.
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

Store one raw id and one level per mixed-level cell. The aligned columns follow
[`MultiOrderVector`](@ref) order, ascending by subtree-interval start.

The store declares `refinement_level: null` and `compression: "compacted"`;
its `refinement_levels` attribute names the level column. Reading restores a
`MultiOrderVector` axis.

This package defines the layout as an extension to v1 of
`zarr-conventions/dggs`. Version 1 specifies `compression: "none"` when
`refinement_level` is null and a `*uniq` scheme without a HEALPix level column.
This package supplies the reader for the extended layout.
"""
struct CompactedEncoding <: CellEncoding end

"""
    ENCODING_REGISTRY

Map store vocabulary (`"none"`, `"ranges"`, `"implicit"`, `"compacted"`) to
[`CellEncoding`](@ref) instances. Conventions and keyword symbols resolve
through the same table. [`register_encoding!`](@ref) adds a downstream
encoding.
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

Both a store's `compression` attribute and `dggwrite`'s `encoding` keyword
resolve through [`ENCODING_REGISTRY`](@ref). Reading also requires a
[`cellaxis`](@ref) implementation; writing requires the write-path verbs.
Incomplete registrations raise a named unsupported-encoding error.
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

Build the stored axis as an `AbstractVector` of typed cell ids. `source` depends
on the encoding:

| encoding | `source` | cost |
|---|---|---|
| [`RangesEncoding`](@ref) | the `(n, 2)` inclusive-range array | arithmetic, no IO |
| [`ImplicitEncoding`](@ref) | the axis length | arithmetic, no IO |
| [`DenseEncoding`](@ref) | the id vector, lazy or not | one chunked pass |
| [`CompactedEncoding`](@ref) | the level and id columns | one validating pass |

Single-level encodings take a level grid and return a `ChunkedCellVector`.
Compacted encoding takes a system and returns a [`MultiOrderVector`](@ref).
"""
function cellaxis end

# --- write-side eligibility and range construction -------------------------

"""
    write_eligible(enc::CellEncoding, grid::AbstractGrid, ids::AbstractVector) -> Bool

Return whether `ids` can be written with `enc`. `encoding = :auto` uses these
rules:

  - [`RangesEncoding`](@ref) requires sorted, unique cells from `grid`'s level.
  - [`ImplicitEncoding`](@ref) additionally requires the complete level.
  - [`DenseEncoding`](@ref) accepts every axis, making `:auto` total.

`ids` may contain raw integers or typed cells. The downstream fallback returns
`true`; the write path separately checks for the required implementation verbs.
"""
write_eligible(::CellEncoding, ::AbstractGrid, ::AbstractVector) = true

write_eligible(::DenseEncoding, ::AbstractGrid, ::AbstractVector) = true

write_eligible(::RangesEncoding, grid::AbstractGrid, ids::AbstractVector) =
    _sorted_unique_cells(grid, ids)

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

# Returning `nothing` prevents a typed cell from validating against another level.
_rawid(grid::AbstractGrid, x::Integer) = convert(idtype(grid), x)
_rawid(grid::AbstractGrid, c::AbstractCellIndex) =
    level(c) == level(grid) ? convert(idtype(grid), rawid(c)) : nothing

"""
    idranges(grid::AbstractGrid, ids::AbstractVector; merge = :rank) -> Matrix

Build one inclusive `[start, stop]` row per run in sorted, unique, single-level
`ids`. `merge` defines a run:

  - `:rank` forms maximal runs of consecutive cells. These intervals may span
    unused integers at Z7 digit rollovers and deleted pentagon branches.
  - `:step` forms runs of ids separated by the level's integer unit. Every
    integer in such an interval names a cell, so structural and cell-aware
    readers produce the same count.

Both modes encode the same axis, and [`cellaxis`](@ref) reads either mode.
Published IGEO7 range stores use interoperable `:step` runs.

Throws a `DGGSFormatError` when `ids` is not [`write_eligible`](@ref) for
[`RangesEncoding`](@ref).

!!! note "Runs across a pentagon's deleted branch"
    At resolution 3, a structural reader counts 4,116 digit strings where a
    `:rank` interval contains 3,432 cells. The published IGEO7 stores use
    `:step` intervals that exclude those phantom ids.
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

# The first sibling pair exposes the level's raw-id unit without rollover effects.
function _idunit(grid::AbstractGrid)
    ncells(grid) >= 2 || return one(idtype(grid))
    return idselect(grid, 1) - idselect(grid, 0)
end

"""
    rangeindex(grid::AbstractGrid, ranges::AbstractMatrix)
        -> (starts, stops, rstart, offsets)

Validate an `(n, 2)` inclusive-range array and build its rank/select index.
`rstart[i]` gives interval `i`'s global starting rank; `offsets[i]` gives the
number of preceding cells; `offsets[end]` gives the total.

Rows must contain non-empty, ascending, disjoint intervals. A malformed row
raises `DGGSFormatError` with its index. The method performs `O(n * level)` digit
arithmetic and reads no data arrays.
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

Validate an `(n, 2)` inclusive-range array against the store's declared data
length. [`rangeindex`](@ref) first checks the rows, then the method requires a
total cell count of `declared`. This catches exclusive stops, incorrect steps
and structural-count mismatches before reading data.
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

"""
    cellaxis(CompactedEncoding(), sys::AbstractHierarchicalGridSystem,
             cell_levels::AbstractVector{<:Integer}, cell_ids::AbstractVector;
             declared_length = nothing) -> MultiOrderVector

Build the mixed-level axis of a `compression: "compacted"` store. Index `k`
combines raw id `cell_ids[k]` with level `cell_levels[k]`.

The method validates that:

  - both columns have the declared length;
  - every level belongs to the system;
  - every id names a cell at its paired level; and
  - the cells form disjoint subtrees in ascending subtree-interval order.

The method preserves the stored order because sorting the cells would misalign
them with the store's data values.
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
        # Check membership before narrowing so an out-of-range UInt gets a format error.
        raw = cell_levels[k]
        raw in levels(sys) || throw(DGGSFormatError(check=:invalid_stored_level,
            declared=Int(first(levels(sys))):Int(last(levels(sys))), observed=raw,
            detail="cell-axis index $k declares level $raw outside the " *
                   "$(nameof(typeof(sys))) level range."))
        l = Int(raw)
        grid = get!(() -> levelgrid(sys, l), grids, l)
        id = storedid(idtype(grid), cell_ids[k])
        idvalid(grid, id) || throw(DGGSFormatError(check=:id_names_no_cell,
            declared=l, observed=id,
            detail="the id $id at index $k of the cell axis names no cell " *
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
        detail="the compacted cell axis must ascend by subtree-interval start; " *
               "its current order misaligns the store's cells and data arrays. " *
               "Rewrite the store with `dggwrite`."))
    return mov
end

end # module Encodings
