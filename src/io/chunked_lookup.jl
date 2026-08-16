# The stored cell axis and the chunk manifest a read is planned against.
#
# One shape, three backings. The axis is always an `AbstractVector` of typed
# cell ids; what differs is where element `k` comes from:
#
#   ranges   -> rank/select arithmetic over the stored intervals   (no IO)
#   implicit -> rank/select arithmetic over the whole level        (no IO)
#   dense    -> one cached chunk of the stored id array            (one read)
#
# The manifest is the single source of truth about the chunk grid: per chunk the
# first id, the last id, and the length. Selection is two-level — binary search
# the manifest to name the one chunk an id could be in, then search inside it —
# so resolving a selector costs at most one chunk, never a scan.

"""
    ChunkedLookups

The lazy cube axis of a DGGS store: [`ChunkManifest`](@ref),
[`ChunkedCellVector`](@ref), and the `DimensionalData` lookup over them,
[`ChunkedCellLookup`](@ref).

Where `CellLookup` compresses a cell set into position windows it computes
ids from, this lookup describes a set someone else has already written down,
chunk by chunk, and resolves selectors without reading more of it than the one
chunk an answer can be in.
"""
module ChunkedLookups

import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractCellIndex,
    AbstractHierarchicalGridSystem, ncells, cellindex, cellposition, cellat,
    level, system, cellindextype, rawid, query, descendants,
    has_sorted_subtrees, level_ranges, MultiOrderCoverage, CellVector,
    CellLookup, Covering, covering_positions

import ..Encodings
using ..Encodings: CellEncoding, DenseEncoding, RangesEncoding, ImplicitEncoding,
    cellaxis, idrank, idselect, idcount_between, idvalid, idcell, idtype, idlevel,
    rangeindex, checkcount
# From `errors.jl`, which the including module reads before this file.
import ..DGGSFormatError

# WIRING: the scan's rejections throw `DGGSFormatError` with no store context —
# a scan sees ids, not stores. The extension that opened the store enriches them
# with the identifier and the conventions that fired, with `with_store_context`.

import DimensionalData as DD
using DimensionalData: Lookups, Dimensions

export ChunkManifest, nchunks, chunkof, chunkbounds
export ChunkedCellVector, axisposition, chunkmanifest
export ChunkedCellLookup

# ===========================================================================
# The manifest
# ===========================================================================

"""
    ChunkManifest(axis::ChunkedCellVector, chunklength::Integer)

The chunk grid of a stored array, described in cells: for chunk `c`, the first
and last cell id it holds, how many cells that is, and how many precede it.

    manifest.firstids[c]   manifest.lastids[c]
    manifest.lengths[c]    manifest.offsets[c]     # cells before chunk c

This is what the lookup owns and what every IO decision is made from — which
chunks a selection touches, where a halo's neighbours live, which chunk an id
resolves in. Chunking is a property of an ARRAY, not of the axis, so a store
whose data variables are chunked differently from its coordinate has one
manifest per chunk length; each is built by the same call.

The constructor reads only the axis's chunk boundaries, so on a computed axis
([`RangesEncoding`](@ref), [`ImplicitEncoding`](@ref)) it is closed-form
arithmetic over the stored intervals and touches no store at all, however many
cells the axis holds. On a dense axis it costs the two boundary reads per chunk
that the cache usually already holds.

The final chunk is short whenever the length is not a multiple of
`chunklength`; Zarr pads it on disk and the manifest does not. Every other
chunk holds exactly `chunklength` cells, and the constructor refuses a manifest
that says otherwise: [`chunkof`](@ref) divides by `chunklength` rather than
searching `offsets`, so a manifest read back from a sidecar under a different
chunk length would resolve every position into the wrong chunk.
"""
struct ChunkManifest{I<:Integer}
    firstids::Vector{I}
    lastids::Vector{I}
    lengths::Vector{Int}
    offsets::Vector{Int}
    chunklength::Int

    function ChunkManifest{I}(firstids::Vector{I}, lastids::Vector{I},
        lengths::Vector{Int}, offsets::Vector{Int}, chunklength::Integer) where {I}
        cl = Int(chunklength)
        _checkmanifest(firstids, lastids, lengths, offsets, cl)
        return new{I}(firstids, lastids, lengths, offsets, cl)
    end
end

ChunkManifest(firstids::Vector{I}, lastids::Vector{I}, lengths::Vector{Int},
    offsets::Vector{Int}, chunklength::Integer) where {I<:Integer} =
    ChunkManifest{I}(firstids, lastids, lengths, offsets, chunklength)

function _checkmanifest(firstids, lastids, lengths, offsets, cl::Int)
    nc = length(lengths)
    all(==(nc), (length(firstids), length(lastids), length(offsets))) ||
        throw(DGGSFormatError(check=:inconsistent_manifest,
            observed=(length(firstids), length(lastids), nc, length(offsets)),
            detail="a chunk manifest holds one firstid, lastid, length and " *
                   "offset per chunk."))
    cl > 0 || throw(DGGSFormatError(check=:inconsistent_manifest, observed=cl,
        detail="a chunk manifest's chunk length is positive."))
    for c in 1:nc
        want = c == nc ? (1 <= lengths[c] <= cl) : lengths[c] == cl
        (want && offsets[c] == (c - 1) * cl) || throw(DGGSFormatError(
            check=:inconsistent_manifest, declared=cl,
            observed=(chunk=c, length=lengths[c], offset=offsets[c]),
            detail="chunk $c holds $(lengths[c]) cells at offset $(offsets[c]), " *
                   "which is not what uniform chunks of $cl look like; only the " *
                   "last chunk may be short."))
    end
    return nothing
end

"""
    nchunks(m::ChunkManifest) -> Int

The number of chunks the axis spans.
"""
nchunks(m::ChunkManifest) = length(m.lengths)

Base.length(m::ChunkManifest) = isempty(m.lengths) ? 0 : m.offsets[end] + m.lengths[end]

Base.:(==)(a::ChunkManifest, b::ChunkManifest) =
    a.firstids == b.firstids && a.lastids == b.lastids &&
    a.lengths == b.lengths && a.offsets == b.offsets &&
    a.chunklength == b.chunklength

Base.show(io::IO, m::ChunkManifest) = print(io, "ChunkManifest(", nchunks(m),
    " chunks of ", m.chunklength, ", ", length(m), " cells)")

"""
    chunkof(m::ChunkManifest, k::Integer) -> Int

The chunk holding axis position `k`.
"""
function chunkof(m::ChunkManifest, k::Integer)
    1 <= k <= length(m) || throw(BoundsError(m, k))
    return (Int(k) - 1) ÷ m.chunklength + 1
end

"""
    chunkbounds(m::ChunkManifest, c::Integer) -> UnitRange{Int}

The axis positions chunk `c` holds.
"""
chunkbounds(m::ChunkManifest, c::Integer) =
    (m.offsets[c]+1):(m.offsets[c]+m.lengths[c])

# The prune half of the two-level search: the only chunk whose id interval can
# contain `x`, or `nothing`. `lastids` ascends because the axis is verified
# sorted, which is what makes one binary search enough.
function _prunechunk(m::ChunkManifest, x)
    c = searchsortedfirst(m.lastids, x)
    c > nchunks(m) && return nothing
    return x >= m.firstids[c] ? c : nothing
end

# ===========================================================================
# The axis
# ===========================================================================

"""
    ChunkedCellVector(grid, source, length)

The stored cell axis as a lazy `AbstractVector` of typed cell ids at one level.

Semantically it **is** the id vector the store wrote: `length` is the number of
cells, `axis[k]` is the `k`th of them, `collect(axis)` is the vector itself.
What backs it is the encoding's business — closed-form rank/select over stored
intervals, over the whole level, or a cached chunk of a stored id array — and
[`axisposition`](@ref) is the inverse in every case.

Build one with [`cellaxis`](@ref), never directly.
"""
struct ChunkedCellVector{ID,G<:AbstractGrid,S} <: AbstractVector{ID}
    grid::G
    source::S
    length::Int
end

function ChunkedCellVector(grid::AbstractGrid, source, n::Integer)
    ID = cellindextype(system(grid))
    return ChunkedCellVector{ID,typeof(grid),typeof(source)}(grid, source, Int(n))
end

Base.size(axis::ChunkedCellVector) = (axis.length,)
Base.IndexStyle(::Type{<:ChunkedCellVector}) = Base.IndexLinear()

Base.@propagate_inbounds function Base.getindex(axis::ChunkedCellVector, k::Int)
    @boundscheck checkbounds(axis, k)
    return idcell(axis.grid, rawcell(axis, k))
end

"""
    rawcell(axis::ChunkedCellVector, k::Integer) -> Integer

The RAW stored id at axis position `k` — `axis[k]` without the typed wrapper,
which is what a comparison against stored bytes wants.
"""
rawcell(axis::ChunkedCellVector, k::Integer) = _rawcell(axis.source, axis, Int(k))

"""
    axisposition(axis::ChunkedCellVector, id::Integer) -> Union{Int,Nothing}

The position of raw id `id` in the axis, or `nothing` when the axis does not
hold it. The inverse of [`rawcell`](@ref), and the half of the bijection every
selector ends at.

Resolution is two-level wherever the ids are stored rather than computed: the
manifest names the one chunk `id` could be in, and only that chunk is read.
"""
axisposition(axis::ChunkedCellVector, id::Integer) =
    _axisposition(axis.source, axis, id)

DGG.system(axis::ChunkedCellVector) = system(axis.grid)
DGG.level(axis::ChunkedCellVector) = level(axis.grid)

"""
    cellposition(axis::ChunkedCellVector, c::AbstractCellIndex) -> Union{Int,Nothing}

The position of a typed cell id in the axis, or `nothing` when the axis does not
hold it — including when `c` is at another level.
"""
function DGG.cellposition(axis::ChunkedCellVector, c::AbstractCellIndex)
    level(c) == level(axis.grid) || return nothing
    return axisposition(axis, convert(idtype(axis.grid), rawid(c)))
end

Base.:(==)(a::ChunkedCellVector, b::ChunkedCellVector) =
    system(a) == system(b) && level(a) == level(b) && length(a) == length(b) &&
    all(k -> rawcell(a, k) == rawcell(b, k), eachindex(a))

Base.show(io::IO, axis::ChunkedCellVector) =
    print(io, "ChunkedCellVector(", typeof(system(axis)).name.name, ", level=",
        level(axis), ", ncells=", length(axis), ", ",
        Encodings.encodingname(encoding(axis)), ")")

Base.show(io::IO, ::MIME"text/plain", axis::ChunkedCellVector) = show(io, axis)

# --- ranges: arithmetic over the stored intervals --------------------------

# `rstart[i]` is the global rank of interval `i`'s first cell and `offsets[i]`
# the number of cells before it, so both directions are one binary search over
# a prefix-sum array plus one O(level) digit walk. Nothing else is stored.
struct RangesSource{I<:Integer}
    starts::Vector{I}
    stops::Vector{I}
    rstart::Vector{Int}
    offsets::Vector{Int}
end

function _rawcell(s::RangesSource, axis::ChunkedCellVector, k::Int)
    i = searchsortedlast(s.offsets, k - 1)
    return idselect(axis.grid, s.rstart[i] + (k - 1 - s.offsets[i]))
end

function _axisposition(s::RangesSource, axis::ChunkedCellVector, id::Integer)
    x = convert(idtype(axis.grid), id)
    i = searchsortedlast(s.starts, x)
    (i == 0 || x > s.stops[i]) && return nothing
    # An interval may span integers that name no cell, so membership of the
    # interval is not membership of the axis.
    idvalid(axis.grid, x) || return nothing
    return s.offsets[i] + (idrank(axis.grid, x) - s.rstart[i]) + 1
end

encoding(::ChunkedCellVector{<:Any,<:AbstractGrid,<:RangesSource}) = RangesEncoding()

"""
    cellaxis(RangesEncoding(), grid, ranges::AbstractMatrix; declared_length=nothing)

Build the axis of a `compression: "ranges"` store from its `(n, 2)` array of
inclusive `[start, stop]` raw ids alone. Zero data IO: the length, every id, and
every position are closed-form rank/select arithmetic.

The rows are checked for shape and disjointness as they are indexed, and
`declared_length`, when given, against the closed-form count — the normative
check that the axis and the data agree, and the one that catches an
expansion-semantics disagreement before a byte is read.
"""
function Encodings.cellaxis(::RangesEncoding, grid::AbstractGrid,
    ranges::AbstractMatrix; declared_length::Union{Integer,Nothing}=nothing)
    starts, stops, rstart, offsets = rangeindex(grid, ranges)
    total = offsets[end]
    declared_length === nothing || checkcount(grid, total, declared_length)
    return ChunkedCellVector(grid, RangesSource(starts, stops, rstart, offsets), total)
end

# --- implicit: position is rank --------------------------------------------

struct ImplicitSource end

_rawcell(::ImplicitSource, axis::ChunkedCellVector, k::Int) = idselect(axis.grid, k - 1)

function _axisposition(::ImplicitSource, axis::ChunkedCellVector, id::Integer)
    idvalid(axis.grid, id) || return nothing
    r = idrank(axis.grid, convert(idtype(axis.grid), id))
    return r < length(axis) ? r + 1 : nothing
end

encoding(::ChunkedCellVector{<:Any,<:AbstractGrid,ImplicitSource}) = ImplicitEncoding()

"""
    cellaxis(ImplicitEncoding(), grid, n::Integer)

Build the axis of a store that writes no cell coordinate at all: position `k` is
the cell at rank `k - 1` of the level. `n` is the declared length, which must be
a prefix of the level — the whole of it, for a global store.
"""
function Encodings.cellaxis(::ImplicitEncoding, grid::AbstractGrid, n::Integer)
    0 <= n <= ncells(grid) || throw(DGGSFormatError(check=:invalid_axis_length,
        declared=Int(n), observed=ncells(grid),
        detail="an implicit axis of $n cells does not fit level $(level(grid)), " *
               "which has $(ncells(grid))."))
    return ChunkedCellVector(grid, ImplicitSource(), n)
end

# --- dense: one cached chunk of the stored ids -----------------------------

# The cache holds one decoded chunk. A read walks positions in order and a
# selector resolves inside one chunk, so a single slot serves both without a
# policy.
#
# `verify` is set when the manifest was READ rather than built: the ids behind a
# persisted sidecar may have been rewritten since, so every chunk that is
# decoded is checked against the row that described it. `store` is the label
# that error carries, because a lazy axis is touched long after the boundary
# that would otherwise have named the store had returned, and `sidecar` the name
# of the array the manifest came out of — a manifest is found by its marker, so
# the array it was found in is the only spelling the message can name.
struct DenseSource{V,I<:Integer}
    ids::V
    manifest::ChunkManifest{I}
    cached::Base.RefValue{Int}
    block::Base.RefValue{Vector{I}}
    verify::Bool
    store::Union{String,Nothing}
    sidecar::Union{String,Nothing}
end

DenseSource(ids, manifest::ChunkManifest{I}, cached, block) where {I} =
    DenseSource(ids, manifest, cached, block, false, nothing, nothing)

function _chunkblock(s::DenseSource{V,I}, c::Int) where {V,I}
    s.cached[] == c && return s.block[]
    r = chunkbounds(s.manifest, c)
    block = I[convert(I, x) for x in s.ids[first(r):last(r)]]
    # The spot check that bounds a trusted manifest: the two ids the manifest
    # published for this chunk are the two the chunk begins and ends with. Two
    # comparisons, on data that was being read anyway.
    if s.verify
        m = s.manifest
        (length(block) == m.lengths[c] && first(block) == m.firstids[c] &&
         last(block) == m.lastids[c]) || _stalemanifest(s, c, block)
    end
    s.cached[] = c
    s.block[] = block
    return block
end

@noinline function _stalemanifest(s::DenseSource, c::Int, block)
    m = s.manifest
    throw(DGGSFormatError(check=:stale_manifest, store=s.store,
        declared=(m.firstids[c], m.lastids[c], m.lengths[c]),
        observed=isempty(block) ? (nothing, nothing, 0) :
                 (first(block), last(block), length(block)),
        detail="chunk $c of the cell axis does not hold the ids the store's " *
               "persisted chunk manifest says it does: they have been rewritten " *
               "since the manifest was written. Rewrite the store with " *
               "`dggwrite`, delete its " *
               (s.sidecar === nothing ? "chunk-manifest array" : "`$(s.sidecar)` array") *
               " to make the reader scan the ids instead, or read it once with " *
               "`validate = :scan`."))
end

function _rawcell(s::DenseSource, axis::ChunkedCellVector, k::Int)
    c = chunkof(s.manifest, k)
    return _chunkblock(s, c)[k-s.manifest.offsets[c]]
end

function _axisposition(s::DenseSource, axis::ChunkedCellVector, id::Integer)
    x = convert(idtype(axis.grid), id)
    c = _prunechunk(s.manifest, x)
    c === nothing && return nothing
    block = _chunkblock(s, c)
    j = searchsortedfirst(block, x)
    (j <= length(block) && block[j] == x) || return nothing
    return s.manifest.offsets[c] + j
end

encoding(::ChunkedCellVector{<:Any,<:AbstractGrid,<:DenseSource}) = DenseEncoding()

"""
    cellaxis(DenseEncoding(), grid, ids::AbstractVector;
             chunklength, declared_length = length(ids), samples = nothing)

Build the axis of a store that writes one id per cell, from a single chunked
pass over the id array — not over the data. The pass does four things at once,
because it is the only pass anyone wants to pay for:

  - verifies the ids are strictly ascending, and REJECTS duplicates, naming
    both how many there are and where the first one is;
  - verifies that every id names a cell of this level, which is what catches an
    id set written at another level, or a Z7 id on a pentagon's deleted branch;
  - truncates to `declared_length`, because Zarr pads the final chunk with the
    fill value and a naive concatenation would read that padding as ids;
  - builds the [`ChunkManifest`](@ref) the lookup then resolves against.

The validation policy is a pair. By default — `validate = :strict` — every id is
checked, which is what the pass is for: the ids are all in hand already and the
digit walk is `O(level)` each. `samples = n` is the explicit opt-out, checking
`n` ids per chunk (both ends, which the manifest publishes, plus an even
interior spread) and admitting a phantom anywhere between them; that is what
`validate = :lazy` selects when a caller would rather have the axis now.

`ids` may be lazy: it is read one `chunklength` block at a time, in order.
"""
function Encodings.cellaxis(::DenseEncoding, grid::AbstractGrid,
    ids::AbstractVector; chunklength::Integer=length(ids),
    declared_length::Integer=length(ids),
    samples::Union{Integer,Nothing}=nothing)
    I = idtype(grid)
    n = Int(declared_length)
    cl = max(Int(chunklength), 1)
    0 <= n <= length(ids) || throw(DGGSFormatError(check=:declared_length_mismatch,
        declared=n, observed=length(ids),
        detail="the store declares $n cells but its id array holds $(length(ids))."))
    nc = cld(n, cl)
    firstids = Vector{I}(undef, nc)
    lastids = Vector{I}(undef, nc)
    lengths = Vector{Int}(undef, nc)
    offsets = Vector{Int}(undef, nc)
    previous = nothing
    duplicates = 0
    firstduplicate = 0
    for c in 1:nc
        lo = (c - 1) * cl + 1
        hi = min(c * cl, n)
        block = I[convert(I, x) for x in ids[lo:hi]]
        length(block) == hi - lo + 1 || throw(DGGSFormatError(
            check=:short_chunk_read, declared=hi - lo + 1, observed=length(block),
            detail="reading chunk $c of the cell axis returned $(length(block)) " *
                   "ids, not the $(hi - lo + 1) its chunk grid says it holds."))
        for (j, x) in pairs(block)
            if previous !== nothing && x <= previous
                if x == previous
                    duplicates += 1
                    firstduplicate == 0 && (firstduplicate = lo + j - 1)
                else
                    throw(DGGSFormatError(check=:unsorted_cell_axis,
                        declared=previous, observed=x,
                        detail="the cell axis is not sorted: id $x at position " *
                               "$(lo + j - 1) is below $previous at the position " *
                               "before it."))
                end
            end
            previous = x
        end
        for j in _checkslots(length(block), samples)
            idvalid(grid, block[j]) || throw(DGGSFormatError(
                check=:id_names_no_cell, declared=level(grid), observed=block[j],
                detail="the id $(block[j]) at position $(lo + j - 1) of the cell " *
                       "axis names no cell of level $(level(grid))." *
                       _ownlevel(grid, block[j])))
        end
        firstids[c] = first(block)
        lastids[c] = last(block)
        lengths[c] = length(block)
        offsets[c] = lo - 1
    end
    duplicates == 0 || throw(DGGSFormatError(check=:duplicate_ids,
        observed=duplicates, declared=firstduplicate,
        detail="the cell axis holds $duplicates duplicate ids, the first at " *
               "position $firstduplicate; a DGGS store names each cell once."))
    manifest = ChunkManifest(firstids, lastids, lengths, offsets, cl)
    source = DenseSource(ids, manifest, Ref(0), Ref(I[]))
    return ChunkedCellVector(grid, source, n)
end

"""
    cellaxis(DenseEncoding(), grid, ids::AbstractVector, manifest::ChunkManifest;
             store = nothing, sidecar = nothing)

Build the axis of a dense store whose chunk grid was PERSISTED with it, from
that manifest alone. No id is read: this is the constructor that lets a store of
tens of millions of cells open in constant time, and the reason a writer
persists a manifest at all.

Everything the scanning method PROVES — strictly ascending, no duplicates, every
id a cell of this level — this one takes on the writer's word, which is what the
marker's `validated = "strict"` attests to. **Nothing here checks any of it**,
not in a chunk that is decoded and not in one that is not: what the spot check
below compares is a chunk's first id, its last id and its length, so ids
duplicated or unsorted or naming no cell strictly INSIDE a chunk pass through it
untouched. Declining the sidecar — `dggread(store; validate = :scan)` — is what
puts those checks back.

The word is not taken forever. Every chunk this axis later decodes — a lazy id
access, the fine search inside a selector — has those three numbers checked
against the row that described it, and a store whose ids were rewritten behind
the sidecar throws `DGGSFormatError(check = :stale_manifest)` on the first chunk
it is wrong about rather than answering out of a stale map. `store` labels that
error, because a lazy axis is touched long after the reader that opened it has
returned and no boundary is left to name the store; `sidecar` names the array
the manifest was read from, which the message tells the caller to delete.
"""
function Encodings.cellaxis(::DenseEncoding, grid::AbstractGrid,
    ids::AbstractVector, manifest::ChunkManifest;
    store::Union{AbstractString,Nothing}=nothing,
    sidecar::Union{AbstractString,Nothing}=nothing)
    I = idtype(grid)
    n = length(manifest)
    n <= length(ids) || throw(DGGSFormatError(check=:declared_length_mismatch,
        declared=n, observed=length(ids),
        detail="the persisted chunk manifest describes $n cells but the id " *
               "array holds $(length(ids))."))
    source = DenseSource(ids, manifest, Ref(0), Ref(I[]), true,
        store === nothing ? nothing : String(store),
        sidecar === nothing ? nothing : String(sidecar))
    return ChunkedCellVector(grid, source, n)
end

# What the offending id says about its OWN level, where the id scheme carries
# one ([`idlevel`](@ref)). It is information and not a ruling: the store names
# the level and this reader reads at that one, so the sentence reports a
# disagreement rather than resolving it. Empty where the scheme says nothing,
# which is most of them.
function _ownlevel(grid, x)
    own = idlevel(grid, x)
    (own === nothing || own == level(grid)) && return ""
    return " That id is itself well formed at level $own, so the store's " *
           "attributes and its ids name different levels. Ids kept at their " *
           "own level against a coarser declared one are the reference level " *
           "pattern, which a single-level axis does not express; a store " *
           "meaning level $own says so in its attributes, or is read through " *
           "an asserted `description`."
end

# Which slots of a block have their validity checked: all of them, or the
# sampled subset the opt-out asks for.
_checkslots(n::Int, ::Nothing) = Base.OneTo(n)
_checkslots(n::Int, samples::Integer) = _sampleslots(n, Int(samples))

# Spot-sample positions: both ends of the block, which are the ids the manifest
# itself publishes, and an even interior spread between them.
function _sampleslots(n::Int, samples::Int)
    n == 0 && return Int[]
    s = clamp(samples, 1, n)
    s == 1 && return [1]
    return unique!([1 + ((j - 1) * (n - 1)) ÷ (s - 1) for j in 1:s])
end

# --- the manifest, from any axis -------------------------------------------

ChunkManifest(axis::ChunkedCellVector, chunklength::Integer) =
    _boundary_manifest(axis, chunklength)

# A dense axis already scanned its own chunk grid; asking for it again is free.
function ChunkManifest(axis::ChunkedCellVector{<:Any,<:AbstractGrid,<:DenseSource},
    chunklength::Integer)
    stored = axis.source.manifest
    return Int(chunklength) == stored.chunklength ? stored :
           _boundary_manifest(axis, chunklength)
end

function _boundary_manifest(axis::ChunkedCellVector, chunklength::Integer)
    I = idtype(axis.grid)
    n = length(axis)
    cl = max(Int(chunklength), 1)
    nc = cld(n, cl)
    firstids = Vector{I}(undef, nc)
    lastids = Vector{I}(undef, nc)
    lengths = Vector{Int}(undef, nc)
    offsets = Vector{Int}(undef, nc)
    for c in 1:nc
        lo = (c - 1) * cl + 1
        hi = min(c * cl, n)
        firstids[c] = convert(I, rawcell(axis, lo))
        lastids[c] = convert(I, rawcell(axis, hi))
        lengths[c] = hi - lo + 1
        offsets[c] = lo - 1
    end
    return ChunkManifest(firstids, lastids, lengths, offsets, cl)
end

# ===========================================================================
# Region selection
# ===========================================================================

"""
    covering_positions(axis::ChunkedCellVector, target) -> Vector{Int}

The positions in `axis` of the cells a [`MultiOrderCoverage`](@ref) of `target`
names, ascending — the position-space form of the [`Covering`](@ref) selector.

The coverage is walked in ascending position order, so a chunk-backed axis
touches each chunk the region meets once.

This is the same verb [`CellVector`](@ref) answers, on a stored axis instead of
a computed one.
"""
function covering_positions(axis::ChunkedCellVector, target)
    sys = system(axis.grid)
    l = level(axis.grid)
    set = query(sys, MultiOrderCoverage(target); level=l)
    out = Int[]
    _each_leaf(sys, set, axis.grid, l) do p
        k = axisposition(axis, convert(idtype(axis.grid), rawid(cellindex(axis.grid, p))))
        k === nothing || push!(out, k)
    end
    return issorted(out) ? out : sort!(out)
end

# Ranges where the system has sorted subtrees, explicit descendants otherwise —
# `CellVector` makes the same distinction for the same reason.
function _each_leaf(f, sys::AbstractHierarchicalGridSystem, set, grid, l::Int)
    if has_sorted_subtrees(sys)
        for r in level_ranges(set, l), p in r
            f(p)
        end
    else
        positions = Int[]
        for c in set, d in descendants(sys, c, l)
            p = cellposition(grid, d)
            p === nothing || push!(positions, p)
        end
        for p in sort!(positions)
            f(p)
        end
    end
    return nothing
end

# ===========================================================================
# The lookup
# ===========================================================================

"""
    ChunkedCellLookup(axis::ChunkedCellVector)

The `DimensionalData` lookup over a STORED cell axis. Pair it with `Cells` to
make a cube axis, exactly as [`CellLookup`](@ref) is paired:

```julia
axis = cellaxis(RangesEncoding(), grid, ranges)
A    = DimensionalData.DimArray(data, Cells(ChunkedCellLookup(axis)))
A[Cells(DimensionalData.At(cell))]
A[Cells(Covering(basin))]
```

It answers the same selectors as [`CellLookup`](@ref) — `At` and `Contains` on
a cell id, `Contains` on a lon/lat point, [`Covering`](@ref) on a region — and
differs in where the answer comes from: the axis is what a store wrote, so a
selector is resolved by the manifest first and by at most one chunk of ids
after, never by a scan.

`ForwardOrdered` is claimed on every axis, and how far it is EARNED depends on
where the axis came from. A scanned store proves it: sorted, unique and
single-level are verified id by id at open. A store opened on a persisted
manifest ([`cellaxis`](@ref)) does not — there the three are the writer's
attestation, and all that is verified per chunk, as chunks are decoded, is the
first id, the last id and the length. `dggread(store; validate = :scan)` declines
the attestation and scans.

A SUBSET is no longer a stored axis — indexing or selecting materialises the
cells it names. A sorted, unique subset becomes the package's own compressed
[`CellLookup`](@ref), which every later operation then treats as any other cell
axis; one that is neither becomes an `Unordered` `Categorical` lookup over the
same cells, since a cell axis is sorted by definition and this one is not.
`Base.reverse` is the everyday way to reach the second case.
"""
struct ChunkedCellLookup{ID,A<:ChunkedCellVector} <: Lookups.Lookup{ID,1}
    axis::A
end

ChunkedCellLookup(axis::ChunkedCellVector{ID}) where {ID} =
    ChunkedCellLookup{ID,typeof(axis)}(axis)

ChunkedCellLookup(lk::ChunkedCellLookup) = lk

Base.parent(lk::ChunkedCellLookup) = lk.axis
Base.IndexStyle(::Type{<:ChunkedCellLookup}) = Base.IndexLinear()

Base.@propagate_inbounds Base.getindex(lk::ChunkedCellLookup, k::Int) = parent(lk)[k]
Base.@propagate_inbounds Base.getindex(lk::ChunkedCellLookup, k::CartesianIndex{1}) =
    parent(lk)[k[1]]

for f in (:getindex, :view, :dotview)
    @eval Base.$f(lk::ChunkedCellLookup, ::Colon) = lk
    @eval Base.$f(lk::ChunkedCellLookup, i::AbstractVector{<:Integer}) = _subset(lk, i)
end

Base.reverse(lk::ChunkedCellLookup) = lk[lastindex(lk):-1:firstindex(lk)]

function _subset(lk::ChunkedCellLookup, mask::AbstractArray{Bool})
    axes(mask) == axes(lk) || throw(BoundsError(lk, (mask,)))
    return _subset(lk, findall(mask))
end

# A subset leaves the store behind: its cells are named explicitly, which is
# what `CellVector` compresses.
function _subset(lk::ChunkedCellLookup, idx)
    axis = parent(lk)
    ids = [axis[Int(k)] for k in idx]
    issorted(ids) && allunique(ids) || return Lookups.Categorical(ids;
        order=Lookups.Unordered())
    return CellLookup(CellVector(system(axis), level(axis), ids))
end

"""
    chunkmanifest(lk::ChunkedCellLookup, chunklength) -> ChunkManifest
    chunkmanifest(axis::ChunkedCellVector, chunklength) -> ChunkManifest

The [`ChunkManifest`](@ref) of an array chunked in blocks of `chunklength`
cells along this axis. The lookup owns the chunk pattern, so this is where an
IO plan starts.
"""
chunkmanifest(axis::ChunkedCellVector, chunklength::Integer) =
    ChunkManifest(axis, chunklength)
chunkmanifest(lk::ChunkedCellLookup, chunklength::Integer) =
    ChunkManifest(parent(lk), chunklength)

DGG.system(lk::ChunkedCellLookup) = system(parent(lk))
DGG.level(lk::ChunkedCellLookup) = level(parent(lk))
DGG.cellposition(lk::ChunkedCellLookup, c::AbstractCellIndex) =
    cellposition(parent(lk), c)

encoding(lk::ChunkedCellLookup) = encoding(parent(lk))

# --- DimensionalData plumbing ----------------------------------------------

Lookups.order(::ChunkedCellLookup) = Lookups.ForwardOrdered()
Lookups.metadata(::ChunkedCellLookup) = Lookups.NoMetadata()
Lookups.bounds(lk::ChunkedCellLookup) =
    isempty(lk) ? (nothing, nothing) : (first(lk), last(lk))

Lookups.reducelookup(::ChunkedCellLookup) = Lookups.NoLookup(Base.OneTo(1))

Dimensions.format(lk::ChunkedCellLookup, ::Type, values, axis::AbstractRange) = lk

function Lookups.rebuild(lk::ChunkedCellLookup; data=nothing, kw...)
    (data === nothing || data === lk || data === parent(lk)) && return lk
    return _rebuild(lk, data)
end

_rebuild(lk::ChunkedCellLookup, ids::AbstractVector{<:AbstractCellIndex}) =
    issorted(ids) && allunique(ids) ?
    CellLookup(CellVector(system(lk), level(lk), collect(ids))) :
    Lookups.Categorical(collect(ids); order=Lookups.Unordered())

@noinline _rebuild(lk::ChunkedCellLookup, data) = throw(ArgumentError(
    "a ChunkedCellLookup holds the cell ids a store wrote at one level; it " *
    "cannot be rebuilt around $(typeof(data)). Replace an axis wholesale with " *
    "`set(A, Cells => NoLookup())`."))

Base.:(==)(a::ChunkedCellLookup, b::ChunkedCellLookup) = parent(a) == parent(b)

function Base.show(io::IO, lk::ChunkedCellLookup)
    axis = parent(lk)
    print(io, "ChunkedCellLookup(", typeof(system(axis)).name.name, ", level=",
        level(axis), ", ncells=", length(axis), ", ",
        Encodings.encodingname(encoding(axis)), ")")
end

Base.show(io::IO, ::MIME"text/plain", lk::ChunkedCellLookup) = show(io, lk)

# --- selectors -------------------------------------------------------------

Lookups.hasselection(lk::ChunkedCellLookup, sel::Lookups.At{<:AbstractCellIndex}) =
    cellposition(parent(lk), Lookups.val(sel)) !== nothing

Lookups.hasselection(lk::ChunkedCellLookup, sel::Lookups.Contains{<:AbstractCellIndex}) =
    cellposition(parent(lk), Lookups.val(sel)) !== nothing

Lookups.hasselection(lk::ChunkedCellLookup, sel::Lookups.Contains{<:Tuple{Real,Real}}) =
    _pointposition(parent(lk), Lookups.val(sel)...) !== nothing

Lookups.selectindices(lk::ChunkedCellLookup, sel::Lookups.At{<:AbstractCellIndex}; kw...) =
    _found(lk, cellposition(parent(lk), Lookups.val(sel)), sel)

Lookups.selectindices(lk::ChunkedCellLookup,
    sel::Lookups.Contains{<:AbstractCellIndex}; kw...) =
    _found(lk, cellposition(parent(lk), Lookups.val(sel)), sel)

Lookups.selectindices(lk::ChunkedCellLookup,
    sel::Lookups.Contains{<:Tuple{Real,Real}}; kw...) =
    _found(lk, _pointposition(parent(lk), Lookups.val(sel)...), sel)

Lookups.selectindices(lk::ChunkedCellLookup, sel::Lookups.At{<:Tuple{Real,Real}}; kw...) =
    _found(lk, _pointposition(parent(lk), Lookups.val(sel)...), sel)

Lookups.selectindices(lk::ChunkedCellLookup, sel::Covering; kw...) =
    covering_positions(parent(lk), Lookups.val(sel))

Lookups.selectindices(lk::ChunkedCellLookup, sel::Covering{<:AbstractVector}; kw...) =
    covering_positions(parent(lk), Lookups.val(sel))

# Resolve a point through the complete level and then look the cell up here, so
# a point inside the level but outside the store answers `nothing`.
function _pointposition(axis::ChunkedCellVector, lon::Real, lat::Real)
    c = cellat(axis.grid, lon, lat)
    c === nothing && return nothing
    return cellposition(axis, c)
end

_found(::ChunkedCellLookup, k::Int, sel) = k
_found(lk::ChunkedCellLookup, ::Nothing, sel) = throw(Lookups.SelectorError(lk, sel))

end # module ChunkedLookups
