# z7_ranges.jl — the Z7 monotonic number line, and the lazy id vectors on it.
#
# Pure integer code, like `z7.jl` on which it depends and nothing else: no
# geometry, no `ISEA`, no I/O. It exists because a Z7 archive that stores its
# cell ids as *ranges* (the `compression: "ranges"` form of the Zarr DGGS
# convention) needs a compact, order-preserving projection of the id onto a
# dense line before "contiguous" means anything — and because that projection
# is **not** the one `grid.jl` already provides.
#
# The two number lines, and why the difference matters
# ----------------------------------------------------
# `cell_to_index` (grid.jl) is the *pentagon-aware dense rank*: the 1-based
# position of a cell among the `10·7^res + 2` cells that actually exist at its
# resolution. It skips the deleted digit, so consecutive ranks are consecutive
# cells.
#
# `z7_to_monotonic` (here) is the *base-7 positional value* of the digit
# string, `base·7^res + Σ dₖ·7^(res−k)`, over `[0, 12·7^res)`. It reserves a
# slot for every digit combination, including the deleted-digit paths that name
# no cell.
#
# Both are strictly increasing in the packed `UInt64` — they agree on *order*,
# which is why either would look correct on a spot check. They disagree on
# *gaps*, and the range table is a statement about gaps: `[start, end]` means
# "every cell whose monotonic value lies between these". Expanding a stored
# range with `cell_to_index` therefore yields the wrong cell count, silently.
# The archives are written against the base-7 line, so that is what this file
# implements; `test/DGGSZarr/test_monotonic.jl` pins the two apart.
#
# Provenance: this is the ordinary base-7 place-value reading of a Z7 digit
# string — first-principles arithmetic, matching the `z7_to_monotonic_int` used
# by the archive writer. Nothing here derives from the AGPL reference
# implementation; see AGENTS.md "Provenance constraint".

"""
    z7_to_monotonic(z7) -> UInt64
    z7_to_monotonic(z7, res) -> UInt64

Position of `z7` on the **base-7 monotonic number line**:
`base·7^res + Σₖ dₖ·7^(res−k)`, the place value of the digit string read in
base 7 with the base cell as the leading digit.

This is the projection the `compression: "ranges"` archive convention is
written against — all `7^d` descendants of a prefix occupy one contiguous
block, so a spatially compact cell set collapses to few ranges.

It is **not** [`cell_to_index`](@ref). That one ranks a cell among the cells
that exist (`10·7^res + 2` of them, deleted digit skipped); this one indexes
all `12·7^res` digit combinations, deleted-digit paths included. The two induce
the same order and different spacing — see the note at the top of
`src/IGeo7/z7_ranges.jl` — so they are never interchangeable, however similar
they look on a contiguous sample.

The one-argument form reads the resolution off the id. The two-argument form
requires `res` to equal it, which is what makes a mixed-resolution id vector
fail loudly rather than project onto a line it does not belong to.

Throws [`InvalidZ7Error`](@ref) for a malformed id (`:bad_padding`, via
[`z7_resolution`](@ref)), for `res ∉ 0:$Z7_MAX_RESOLUTION`
(`:resolution_range`) and for `res` disagreeing with the id's own resolution
(`:monotonic_res`).
"""
@inline function z7_to_monotonic(z7::UInt64, res::Integer)
    r = Int(res)
    0 <= r <= Z7_MAX_RESOLUTION || throw(InvalidZ7Error(
        :resolution_range, z7, _z7_int(res), Z7_MAX_RESOLUTION))
    own = z7_resolution(z7)
    own == r || throw(InvalidZ7Error(:monotonic_res, z7, _z7_int(own), r))
    value = UInt64(z7_base_cell(z7))
    for k in 1:r
        value = value * UInt64(7) + UInt64(_z7_digit(z7, k))
    end
    return value
end

@inline z7_to_monotonic(z7::UInt64) = z7_to_monotonic(z7, z7_resolution(z7))
@inline z7_to_monotonic(z7::Unsigned, res::Integer) = z7_to_monotonic(UInt64(z7), res)
@inline z7_to_monotonic(z7::Unsigned) = z7_to_monotonic(UInt64(z7))

"`MONOTONIC_LIMIT[res+1] = 12·7^res` — one past the last valid monotonic value."
const MONOTONIC_LIMIT = ntuple(i -> UInt64(12) * UInt64(7)^(i - 1), Z7_MAX_RESOLUTION + 1)

"""
    z7_from_monotonic(value, res) -> UInt64

Inverse of [`z7_to_monotonic`](@ref): the packed id at position `value` of the
base-7 number line at resolution `res`, padding slots past `res` with the
sentinel `7`.

The result is a structurally well-formed id but not necessarily a *cell*:
`value` may name a deleted-digit path, which [`is_valid_z7`](@ref) rejects.
That is the honest shape of this line — it has more positions than the grid has
cells — and expanding a range table is exactly where it shows up. Callers that
need cells (rather than digit strings) validate.

Throws [`InvalidZ7Error`](@ref) for `res ∉ 0:$Z7_MAX_RESOLUTION`
(`:resolution_range`) and for `value ≥ 12·7^res` (`:monotonic_range`).
"""
function z7_from_monotonic(value::UInt64, res::Integer)
    r = Int(res)
    0 <= r <= Z7_MAX_RESOLUTION || throw(InvalidZ7Error(
        :resolution_range, zero(UInt64), _z7_int(res), Z7_MAX_RESOLUTION))
    limit = @inbounds MONOTONIC_LIMIT[r+1]
    value < limit ||
        throw(InvalidZ7Error(:monotonic_range, zero(UInt64), Int(value), Int(limit - UInt64(1))))
    # Peel the digits low-to-high, writing each into its slot as we go; the
    # remainder left over after `r` divisions is the base cell.
    z = Z7_PAD_MASK
    remainder = value
    for k in r:-1:1
        d = remainder % UInt64(7)
        remainder ÷= UInt64(7)
        z = _z7_set_digit(z, k, d)
    end
    return (remainder << Z7_BASE_SHIFT) | z
end

z7_from_monotonic(value::Unsigned, res::Integer) = z7_from_monotonic(UInt64(value), res)

# ---------------------------------------------------------------------------
# Lazy id vectors
#
# Both types below are `AbstractVector{UInt64}` that an `IGeo7Lookup` can wrap
# without the ids ever existing as a vector — the same trick `DGGSGlobeIds`
# (core/globe_ids.jl) plays for a globe-complete dimension, for the same reason:
# a `DimensionalData` dimension otherwise costs one word per cell, and the
# whole point of a ranges archive is that its index cost is O(R), not O(N).
#
# They share a supertype only so the `IGeo7Lookup` constructor fast path and the
# `strictly_increasing` short circuit can be stated once each. Neither is a new
# lookup family; `IGeo7Lookup{<:Z7RangeIds}` is an ordinary `IGeo7Lookup`, with
# the same selectors and the same `rebuild`.
# ---------------------------------------------------------------------------

"""
    Z7LazyIds <: AbstractVector{UInt64}

Supertype of the Z7 id vectors that compute or defer their contents rather than
storing them: [`Z7RangeIds`](@ref) (arithmetic, O(R) state) and
[`Z7CachedIds`](@ref) (deferred read, O(1) until pulled).

Every subtype guarantees strictly ascending ids, which is what lets
`Helpers.strictly_increasing` short-circuit instead of walking the vector — see
the note on that method.
"""
abstract type Z7LazyIds <: AbstractVector{UInt64} end

# The O(N) pass the `<X>Lookup` constructors run to check sortedness is, for
# these two, a check of something already established: `Z7RangeIds` validated
# its range table at construction (ascending, non-overlapping, non-empty) and
# derives every element by adding a non-negative offset to it, while
# `Z7CachedIds` runs the real check once, at materialization, where the ids
# actually arrive. Walking the vector here would instead *cause* the
# materialization this type exists to defer.
Helpers.strictly_increasing(::Z7LazyIds) = true

Base.IndexStyle(::Type{<:Z7LazyIds}) = Base.IndexLinear()

"""
    z7_level(ids::Z7LazyIds) -> Int

Resolution every id in `ids` is at.
"""
z7_level(g::Z7LazyIds) = g.level

# --- Z7RangeIds ------------------------------------------------------------

"""
    Z7RangeIds(range_table::AbstractMatrix{UInt64}, level)
    Z7RangeIds(starts::AbstractVector{UInt64}, ends::AbstractVector{UInt64}, level)

Every Z7 cell covered by a table of `R` inclusive `[start, end]` ranges on the
base-7 monotonic line, as an `AbstractVector{UInt64}` of length `N` — the
decoded form of a `compression: "ranges"` DGGS archive coordinate.

`range_table` is `R × 2`, column 1 the range starts and column 2 the inclusive
ends, both as **packed Z7 ids** at `level` (that is the on-disk layout; note
that Zarr stores it C-ordered, so a reader must transpose). The ranges must be
ascending, non-empty and non-overlapping; that is checked here, once, in O(R).

State is `O(R)` and independent of `N`: two vectors of `R` words. At Estonia
res-12 (`R ≈ 39,697`, `N ≈ 11.9M`) that is ~1 MB in place of ~95 MB, which is
the property the whole ranges form exists to buy.

```julia
ids = Z7RangeIds(table, 10)          # 136 ranges -> 3101 cells
length(ids)                          # 3101
ids[1]                               # first cell, O(log R), no allocation
```

Element access and the reverse lookup (`DiscreteGlobalGrids.cell_position`) are
both `O(log R)` binary searches, so an `IGeo7Lookup` over one of these answers
`At(id)` without ever materializing the dimension.
"""
struct Z7RangeIds <: Z7LazyIds
    # Monotonic value of each range's first cell; ascending, non-overlapping.
    starts::Vector{UInt64}
    # Cumulative cell counts, length R+1, `offsets[1] == 0`. Range `r` owns
    # positions `offsets[r]+1 : offsets[r+1]`, so the position -> range step is
    # one `searchsortedlast` and the range -> position step is a subtraction.
    offsets::Vector{Int}
    level::Int

    function Z7RangeIds(starts::Vector{UInt64}, offsets::Vector{Int}, level::Int)
        return new(starts, offsets, level)
    end
end

function Z7RangeIds(starts_z7::AbstractVector, ends_z7::AbstractVector, level::Integer)
    r = Int(level)
    0 <= r <= Z7_MAX_RESOLUTION || throw(InvalidZ7Error(
        :resolution_range, zero(UInt64), _z7_int(level), Z7_MAX_RESOLUTION))
    n = length(starts_z7)
    n == length(ends_z7) || throw(ArgumentError(
        "Z7 range starts and ends must have equal length, got $n and $(length(ends_z7))"))

    starts = Vector{UInt64}(undef, n)
    offsets = Vector{Int}(undef, n + 1)
    @inbounds offsets[1] = 0
    previous_end = zero(UInt64)
    @inbounds for i in 1:n
        s = z7_to_monotonic(UInt64(starts_z7[i]), r)
        e = z7_to_monotonic(UInt64(ends_z7[i]), r)
        s <= e || throw(ArgumentError(
            "Z7 range $i ends before it starts (monotonic $s > $e)"))
        # `>` and not `>=`: two ranges that merely touch would still be a
        # writer bug (they should have been one range), but an *overlap* is the
        # one that breaks the position arithmetic, by giving a cell two
        # positions. Only that is rejected.
        i == 1 || s > previous_end || throw(ArgumentError(
            "Z7 range $i starts at monotonic $s, overlapping the previous range's end $previous_end"))
        starts[i] = s
        previous_end = e
        offsets[i+1] = offsets[i] + Int(e - s) + 1
    end
    return Z7RangeIds(starts, offsets, r)
end

function Z7RangeIds(range_table::AbstractMatrix, level::Integer)
    size(range_table, 2) == 2 || throw(ArgumentError(
        "a Z7 range table must be R × 2 (start, end), got $(size(range_table))"))
    return Z7RangeIds(view(range_table, :, 1), view(range_table, :, 2), level)
end

Base.size(g::Z7RangeIds) = (@inbounds(g.offsets[end]),)

"Number of ranges backing `g` — the `R` whose independence from `N` is the point."
z7_nranges(g::Z7RangeIds) = length(g.starts)

@inline function Base.getindex(g::Z7RangeIds, i::Int)
    @boundscheck checkbounds(g, i)
    # `offsets` is strictly increasing (every range holds at least one cell),
    # so the last offset `<= i-1` is the range that owns position `i`.
    r = searchsortedlast(g.offsets, i - 1)
    @inbounds return z7_from_monotonic(g.starts[r] + UInt64(i - 1 - g.offsets[r]), g.level)
end

"""
    z7_range_position(g::Z7RangeIds, id) -> Union{Nothing,Int}

Position of `id` in `g`, or `nothing` when `g` does not cover it — the inverse
of `getindex`, in `O(log R)` and without materializing anything.

The answer is exact: an id at the wrong resolution, a malformed id, and an id
falling in a gap between two ranges are all `nothing` rather than a wrong
position or a throw. This is what `DiscreteGlobalGrids.cell_position` dispatches
to, and hence what makes `At`/`Contains` on a ranges-backed lookup range-aware.
"""
function z7_range_position(g::Z7RangeIds, id::UInt64)
    # A membership test, so every way `id` can fail to be a cell of this level
    # answers "not here" rather than propagating — the same contract
    # `cell_position` states for the globe branch.
    value = try
        z7_to_monotonic(id, g.level)
    catch
        return nothing
    end
    r = searchsortedlast(g.starts, value)
    r == 0 && return nothing
    @inbounds offset = g.offsets[r]
    @inbounds span = g.offsets[r+1] - offset
    delta = value - @inbounds(g.starts[r])
    delta < UInt64(span) || return nothing      # past this range's end: in a gap
    return offset + Int(delta) + 1
end

z7_range_position(g::Z7RangeIds, id::Unsigned) = z7_range_position(g, UInt64(id))

# `AbstractArray`'s `show` prints elements, which on a large archive is the
# materialization these types exist to avoid — `DimensionalData`'s non-compact
# `show(::Lookup)` takes exactly that path. The three numbers below are the
# whole content anyway: every id follows from them.
Base.show(io::IO, g::Z7RangeIds) = print(io, "Z7RangeIds(level=", g.level,
    ", ranges=", z7_nranges(g), ", cells=", length(g), ")")
Base.show(io::IO, ::MIME"text/plain", g::Z7RangeIds) = show(io, g)

# --- Z7CachedIds -----------------------------------------------------------

"""
    Z7CachedIds(source, level)
    Z7CachedIds(source, level, length)

A dense Z7 id array that is **not read until something needs it** — the
`compression: "none"` counterpart of [`Z7RangeIds`](@ref), where the archive
really does store one id per cell and there is no arithmetic to replace them
with.

`source` is any indexable, sized store of `UInt64` — in practice the lazy
`Zarr.ZArray` of the archive's `cell_ids` coordinate. Opening a dataset,
printing its structure, and reading its resolution all stay O(1) in `N`; the
first operation that needs the ids as a whole (a cell-id selection, an
iteration) pulls them once and caches the result.

The deferred read is also where the sortedness contract is *checked*: the
convention requires the coordinate be ascending and unique, and
[`z7_materialize!`](@ref) verifies that at the moment the ids arrive, throwing
if the archive lied. Until then the guarantee is asserted, not tested — which
is the honest cost of not reading `N` words at open time.

Materialization is idempotent and races benignly: two threads that pull at once
each build an identical vector and one wins, so no lock is taken.
"""
mutable struct Z7CachedIds{S} <: Z7LazyIds
    const source::S
    const level::Int
    const length::Int
    cache::Union{Nothing,Vector{UInt64}}
end

function Z7CachedIds(source, level::Integer, length::Integer)
    r = Int(level)
    0 <= r <= Z7_MAX_RESOLUTION || throw(InvalidZ7Error(
        :resolution_range, zero(UInt64), _z7_int(level), Z7_MAX_RESOLUTION))
    n = Int(length)
    n >= 0 || throw(ArgumentError("Z7CachedIds length must be non-negative, got $n"))
    return Z7CachedIds{typeof(source)}(source, r, n, nothing)
end

Z7CachedIds(source, level::Integer) = Z7CachedIds(source, level, length(source))

Base.size(g::Z7CachedIds) = (g.length,)

"""
    z7_materialize!(g::Z7CachedIds) -> Vector{UInt64}

Read `g`'s ids from its source (once) and return the cached vector, checking on
the way that they are strictly ascending and all at `g.level` — the validation
[`Z7CachedIds`](@ref) defers from construction to first use.
"""
function z7_materialize!(g::Z7CachedIds)
    cached = g.cache
    cached === nothing || return cached
    ids = Vector{UInt64}(undef, g.length)
    copyto!(ids, g.source[1:g.length])
    Helpers.strictly_increasing(ids) || throw(ArgumentError(
        "the archive's dense cell_ids coordinate is not strictly ascending"))
    for id in ids
        z7_resolution(id) == g.level || throw(InvalidZ7Error(
            :monotonic_res, id, _z7_int(z7_resolution(id)), g.level))
    end
    g.cache = ids
    return ids
end

"`true` once `g` has been read from its source."
z7_is_materialized(g::Z7CachedIds) = g.cache !== nothing

# Element access materializes — except at the two endpoints.
#
# Reading `source[i]` for a single element looks cheaper than a bulk read and
# is catastrophically not: the source is a *chunked, compressed* store, so one
# element costs one whole decompressed chunk, and walking N of them costs N
# chunks. Measured on the res-12 reference archive (N = 158,430, chunk =
# 39,608), elementwise reads spent 33 s and allocated 49 GiB to deliver 1.3 MB
# of ids. So any index that suggests someone is walking the vector triggers the
# single bulk read instead.
#
# The first and last elements are the exception because of who asks for them.
# Both are probes by callers that only want to know what kind of ids these are
# and never touch the rest: the `IGeo7Lookup` constructor checks the endpoints'
# resolution at open time, and `DimensionalData`'s compact lookup display
# prints `v[begin]`, `…`, `v[end]` — note that it indexes the *parent* vector,
# so forwarding `first`/`last` on the lookup would not reach it. Serving those
# two from the source is one chunk each and keeps "open an archive and look at
# its structure" free of reading the coordinate.
@inline function Base.getindex(g::Z7CachedIds, i::Int)
    @boundscheck checkbounds(g, i)
    cached = g.cache
    cached === nothing || return @inbounds cached[i]
    (i == 1 || i == g.length) && return UInt64(g.source[i])
    return @inbounds z7_materialize!(g)[i]
end

"""
    z7_cached_position(g::Z7CachedIds, id) -> Union{Nothing,Int}

Position of `id` in `g`, or `nothing`. Binary search over the dense ids, so this
is the call that pulls them: `O(log N)` once the read has happened, and the read
itself the first time.
"""
function z7_cached_position(g::Z7CachedIds, id::UInt64)
    ids = z7_materialize!(g)
    i = Helpers.sorted_index(ids, id)
    return iszero(i) ? nothing : i
end

z7_cached_position(g::Z7CachedIds, id::Unsigned) = z7_cached_position(g, UInt64(id))

Base.show(io::IO, g::Z7CachedIds) = print(io, "Z7CachedIds(level=", g.level,
    ", cells=", g.length, z7_is_materialized(g) ? ", read" : ", unread", ")")
Base.show(io::IO, ::MIME"text/plain", g::Z7CachedIds) = show(io, g)
