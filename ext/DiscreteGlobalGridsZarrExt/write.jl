# The write half: a DimensionalData cube over a cell axis to a Zarr v2 store.
#
# Everything format-semantic is borrowed, not restated. The encoding decides the
# layout of the cell coordinate (`idranges`, `cellaxis`), the conventions decide
# the attributes (`encode!` on a `StoreSnapshot`), and the lookup layer decides
# the manifest (`chunkmanifest`). What lives here is the order those are called
# in, the chunk plan, and the Zarr calls themselves.
#
# The pipeline, once:
#
#   cube -> cell axis -> encoding -> chunk plan -> snapshot -> encode! -> bytes
#
# The snapshot is built BEFORE anything is created, so the group and its arrays
# are written with their final attributes in one pass and no store is ever left
# half-stamped.
#
# A submodule so that this file's helper names cannot collide with the read
# half's, which shares the extension's namespace.

module DGGSZarrWrite

import DiscreteGlobalGrids as DGG
using ..DiscreteGlobalGridsZarrExt: storeidentifier,
    MANIFEST_MARKER, MANIFEST_WRITER, MANIFEST_FORMAT, MANIFEST_VALIDATED,
    COMPACTED_LEVELS_ARRAY
using DiscreteGlobalGrids: AbstractCellIndex, ArrayEntry,
    CellEncoding, CellLookup, ChunkManifest, ChunkedCellLookup, DGGSFormatError,
    DEFAULT_WRITE_CONVENTIONS, DGGSConvention, DenseEncoding, ENCODING_REGISTRY,
    CompactedEncoding, ImplicitEncoding, MultiOrderLookup, MultiOrderVector,
    RangesEncoding, StoreDescription, StoreSnapshot, XdggsConvention,
    ancestor, cellaxis, chunkmanifest, encodingname, has_sorted_subtrees,
    idcell, idranges, idtype, level, levels, levelgrid, rawid, system,
    write_eligible
import DimensionalData as DD
import Zarr

# --- the on-disk contract ---------------------------------------------------
#
# Names follow the published IGEO7 stores (dggs-storage-landscape.md 4.1-4.3):
# `spatial_dimension` names the dimension the DATA variables share and
# `coordinate` names the array that encodes it, and in the ranges case those are
# deliberately different — the coordinate is `(n, 2)` and carries neither.

const CELL_IDS_ARRAY = "cell_ids"
const CELL_RANGES_ARRAY = "cell_id_ranges"
const SPATIAL_DIMENSION = "cell_ids"
const RANGES_DIMS = ["ranges", "bounds"]
const MANIFEST_ARRAY = "cell_chunk_manifest"
const MANIFEST_DIMS = ["chunks", "bounds"]
const ARRAY_DIMENSIONS = "_ARRAY_DIMENSIONS"

"""
    DEFAULT_CHUNK_TARGET

The number of ELEMENTS `chunks = :auto` aims a chunk at — cells times the
extents of every other dimension, not cells alone. Roughly a million, which is
the size the published stores' aperture-7 chunking lands on for a
one-value-per-cell layer; `chunk_target` overrides it, and a test axis of a few
hundred cells collapses to one chunk.
"""
const DEFAULT_CHUNK_TARGET = 1_000_000

# Keyword sugar over `ENCODING_REGISTRY`: the same table a convention resolves a
# store's vocabulary through, so a downstream encoding needs no keyword of its own.
const ENCODING_KEYWORDS = Dict(:dense => "none", :ranges => "ranges",
    :implicit => "implicit", :compacted => "compacted")

const REMOTE_SCHEMES = ("gs://", "s3://", "http://", "https://", "az://", "abfs://")

# ===========================================================================
# Entry points
# ===========================================================================

const Cube = Union{DD.AbstractDimArray,DD.AbstractDimStack}

"""
    dggwrite(dest, stack_or_array; encoding = :auto,
             conventions = DEFAULT_WRITE_CONVENTIONS, chunks = :auto,
             merge = :step, chunk_target = DEFAULT_CHUNK_TARGET) -> dest

Write a `DimStack` or `DimArray` over a cell axis to a **Zarr v2 directory
store**, consolidated metadata included. `dest` is a local path or an open
writeable `Zarr.ZGroup`; a `gs://`/`s3://`/`https://` URL is refused rather than
half-written.

The cell dimension carries a `CellLookup` or `ChunkedCellLookup` (a
single-level axis, sorted and unique — `reverse` and friends degrade it to a
`Categorical`, which is refused) or a `MultiOrderLookup` (a mixed-level axis,
written as `compacted`).

  - `encoding = :auto` writes a mixed-level axis as compacted, ranges where a
    single-level axis is eligible, and dense otherwise. `:dense` is the interop
    escape for readers that cannot expand ranges, `:ranges` forces the compact
    form, `:implicit` writes no cell coordinate at all — position is the cell —
    which needs a whole level, and `:compacted` names the mixed-level layout:
    two aligned columns, `cell_ids` and `cell_levels`, under
    `refinement_level: null`. A single-level encoding requested for a
    mixed-level axis is refused; `expand(A, level)` is the bridge.
  - `merge = :step` merges only ids adjacent as integers, so no interval can
    enclose an id that names no cell — what a structural-count reader needs, and
    what the published IGEO7 range stores hold. `merge = :rank` merges runs of
    consecutive CELLS instead, giving the fewest rows, and is read back correctly
    only by a rank-aware reader such as this package; see [`idranges`](@ref).
  - `chunks = :auto` groups whole coarse-ancestor subtree runs into chunks of
    about `chunk_target` elements; an integer is a fixed chunk length in CELLS.
    See `ChunkPlan` for what that guarantees and what it only aims at.
    `chunk_target` counts the elements of a chunk — cells times the extents of
    the non-cell dimensions, which are one chunk each — so a layer with a
    40-step time axis gets a fortieth of the cells per chunk.
  - `conventions` stamps the store, `zarr-conventions/dggs` plus xdggs by
    default, so both a convention-aware reader and xdggs can open it.

The chunk grid is persisted as a `(n_chunks, 2)` sidecar array of per-chunk
first and last id, so a reader need not scan the axis to rebuild it.

**Attributes.** Each layer's `metadata` is written as its array attributes and
the stack's `metadata["attrs"]` as the group's — the two places `dggread` puts
them, so a store read and rewritten keeps its `units`, its `long_name` and its
group vocabulary. Convention-generated keys are stamped OVER the producer's: a
`_ARRAY_DIMENSIONS` or a `dggs` object carried in from another layout would
describe this store wrongly. Other stack metadata is not written; the cell
coordinate is regenerated by the encoding and carries no producer attributes.

**Two documented normalizations of a round trip.** Layers are written in
alphabetical order, so that is the order they come back in whatever order went
in; and a layer's attributes include the `_ARRAY_DIMENSIONS` this writer stamps,
so a stack read back carries it in each layer's `metadata`.

Layers are never overwritten: a `ZGroup` destination that already holds an array
this write would create raises before anything is stamped.
"""
function DGG.dggwrite(dest::AbstractString, src::Cube; kw...)
    path = String(dest)
    _reject_remote(path)
    # `zgroup` refuses a store that is not empty, so a path needs no name guard.
    _write(path, (attrs, names) -> Zarr.zgroup(path; attrs=attrs), src; kw...)
    return dest
end

function DGG.dggwrite(dest::Zarr.ZGroup, src::Cube; kw...)
    dest.writeable || throw(ArgumentError(
        "dggwrite needs a writeable group; this one was opened read-only."))
    _write(storeidentifier(dest), (attrs, names) -> _stamp(dest, attrs, names),
        src; kw...)
    return dest
end

@noinline function _reject_remote(path)
    for scheme in REMOTE_SCHEMES
        startswith(path, scheme) || continue
        throw(ArgumentError("dggwrite writes local directory stores only; " *
                            "$(repr(path)) names a $(rstrip(scheme, ['/', ':'])) store. " *
                            "Write locally and upload, or open the remote group yourself " *
                            "and pass the ZGroup."))
    end
    return nothing
end

# Stamping is the first irreversible step, and `zcreate` throws on the first
# name it cannot take: a group checked afterwards would be left carrying
# attributes for an encoding its arrays do not have, which no reader can open.
function _stamp(g::Zarr.ZGroup, attrs, names)
    taken = sort!(String[n for n in names if haskey(g.arrays, n)])
    isempty(taken) || throw(DGGSFormatError(check=:destination_not_empty,
        declared=taken, observed=sort!(collect(keys(g.arrays))),
        detail="this group already holds " * join(taken, ", ") *
               ", and dggwrite does not overwrite an array. Write to a new " *
               "group, or delete these first."))
    merge!(g.attrs, attrs)
    Zarr.writeattrs(g.zarr_format, g.storage, g.path, g.attrs)
    return g
end

# ===========================================================================
# The pipeline
# ===========================================================================

function _write(identifier, opengroup, src; encoding=:auto,
    conventions=DEFAULT_WRITE_CONVENTIONS, chunks=:auto, merge::Symbol=:step,
    chunk_target::Integer=DEFAULT_CHUNK_TARGET)

    # Both keywords are checked whatever the encoding: `merge` only reaches the
    # ranges coordinate and `chunk_target` only the automatic plan, but a
    # misspelling is a misspelling and silently doing something else is how a
    # store of tens of millions of cells ends up with one chunk per cell.
    merge in (:rank, :step) || throw(ArgumentError(
        "the ranges merge rule is `:rank` or `:step`, not $(repr(merge))"))
    chunk_target >= 1 || throw(ArgumentError(
        "chunk_target counts the elements of a chunk and is at least one, " *
        "not $chunk_target"))

    names = String[DGG.conventionname(c) for c in conventions]
    return DGG.with_store_context(identifier; conventions=names) do
        mixed = _mixedaxis(src)
        mixed === nothing || return _writemixed(opengroup, identifier, src,
            mixed..., encoding, conventions, chunks, Int(chunk_target))
        celldim, grid, cells = _cellaxis(src)
        isempty(cells) && throw(ArgumentError(
            "dggwrite has nothing to write: the cell axis is empty."))
        enc = _encoding(encoding, grid, cells)
        layers = _layers(src, celldim)
        target = _celltarget(Int(chunk_target), layers, celldim)
        plan = _chunkplan(chunks, grid, cells, target)

        # The coordinate is computed once and used twice: it is the array that
        # goes on disk, and it is what the axis is rebuilt from. `cellaxis` is
        # the READER's constructor, so the length, the ids and the manifest a
        # reader will derive are all checked by the reader's own code before a
        # byte is committed.
        coord = _coordinate(enc, grid, cells, merge)
        axis = _axis(enc, grid, coord, plan, length(cells))
        manifest = chunkmanifest(axis, plan.chunklength)

        # The description first: the manifest marker records the level and the
        # grid name the axis was validated at, and those are the description's
        # to say, not the plan's.
        desc = _description(system(grid), level(grid), enc, layers)
        arrays = _arrayplan(enc, coord, layers, celldim, plan, manifest, desc)
        return _commit(opengroup, identifier, src, desc, conventions, arrays)
    end
end

# The mixed-level pipeline: the same shape as the single-level one, with the
# level grid replaced by the container and the coordinate by its two columns.
function _writemixed(opengroup, identifier, src, celldim, mov,
    encoding, conventions, chunks, chunk_target::Int)
    isempty(mov) && throw(ArgumentError(
        "dggwrite has nothing to write: the cell axis is empty."))
    enc = _mixedencoding(encoding)
    layers = _layers(src, celldim)
    target = _celltarget(chunk_target, layers, celldim)
    plan = _movchunkplan(chunks, mov, target)
    lv = Int8[level(c) for c in mov]
    I = idtype(levelgrid(system(mov), Int(maximum(lv))))
    ids = I[convert(I, rawid(c)) for c in mov]
    # The reader's constructor validates the columns — pairwise disjoint,
    # container order — before a byte is committed, as `_axis` does for the
    # single-level encodings.
    cellaxis(enc, system(mov), lv, ids; declared_length=length(mov))
    manifest = _movmanifest(mov, plan.chunklength)
    desc = _description(system(mov), nothing, enc, layers)
    arrays = _arrayplan(enc, (lv, ids), layers, celldim, plan, manifest, desc)
    return _commit(opengroup, identifier, src, desc, conventions, arrays;
        reference_level=DGG.Engine.reference_level(mov))
end

# Stamp, create and fill: everything after the arrays are planned, shared by
# both pipelines. `reference_level` rides into the manifest marker for a
# compacted store, whose description carries no single level.
function _commit(opengroup, identifier, src, desc, conventions, arrays;
    reference_level::Union{Int,Nothing}=nothing)
    reference_level === nothing || _markerreference!(arrays, reference_level)
    snapshot = StoreSnapshot(identifier=identifier, attrs=_groupattrs(src),
        arrays=[a.entry for a in arrays])
    for c in conventions
        _stampable(c, desc) || continue
        DGG.encode!(c, snapshot, desc)
    end
    reference_level === nothing || _declarelevels!(snapshot.attrs)

    group = opengroup(snapshot.attrs, String[a.entry.name for a in arrays])
    for a in arrays
        z = Zarr.zcreate(a.entry.eltype, group, a.entry.name,
            reverse(a.entry.shape)...; chunks=a.chunks, attrs=a.entry.attrs)
        _fill!(z, a.source, first(a.chunks))
    end
    Zarr.consolidate_metadata(group)
    return group
end

# xdggs attributes describe a dense single-level coordinate; stamped on a
# mixed-level store they would make its readers decode one, so they are skipped
# there.
_stampable(::DGGSConvention, desc) = true
_stampable(::XdggsConvention, desc) = desc.level !== nothing

# `coordinate` names the id array; nothing in v1 of zarr-conventions/dggs names
# the level array beside it, so the store says it itself — the one key an
# independent reader needs to align `cell_levels[k]` with `cell_ids[k]`.
function _declarelevels!(attrs)
    dggs = get(attrs, "dggs", nothing)
    dggs isa AbstractDict || return attrs
    dggs["refinement_levels"] = COMPACTED_LEVELS_ARRAY
    return attrs
end

function _markerreference!(arrays, reference_level::Int)
    for a in arrays
        marker = get(a.entry.attrs, MANIFEST_MARKER, nothing)
        marker === nothing || (marker["reference_level"] = reference_level)
    end
    return arrays
end

# `chunk_target` counts the ELEMENTS of a chunk, so the cell-count target is
# what is left after the non-cell dimensions have taken their share — a
# (1M, 40) Float32 chunk is 160 MB and nobody asked for that. The widest layer
# sets it for all of them: one chunk length is shared by the axis and every
# array on it.
function _celltarget(target::Int, layers, celldim)
    trailing = 1
    for (_, A) in layers
        width = 1
        for d in DD.otherdims(A, celldim)
            width *= length(d)
        end
        trailing = max(trailing, width)
    end
    return max(1, target ÷ trailing)
end

# The group attributes a `dggread` stack carries verbatim under `"attrs"`, which
# is where a round trip finds the producer's own vocabulary. The conventions
# stamp on top of these.
_groupattrs(src) = _attrs(get(_attrs(DD.metadata(src)), "attrs", nothing))

_attrs(x) = Dict{String,Any}()
_attrs(md::AbstractDict) = Dict{String,Any}(String(k) => deepcopy(v) for (k, v) in md)
_attrs(md::NamedTuple) = Dict{String,Any}(String(k) => deepcopy(v) for (k, v) in pairs(md))
_attrs(md::DD.Metadata) = _attrs(DD.val(md))

# ===========================================================================
# The cell axis
# ===========================================================================

"""
    _cellaxis(src) -> (dim, grid, ids)

The cube's cell dimension, the complete level grid its cells live at, and the
RAW ids on it.

The ids are read out of the lookup exactly ONCE, and as the integers a store
holds rather than as typed cells, because that vector is the one the dense
coordinate writes: everything downstream — eligibility, the chunk plan, the
coordinate itself — works from this array, and a write of tens of millions of
cells never holds the axis twice. Where a typed cell is wanted, `idcell` puts
the wrapper back on for the one id in hand.

Only this package's own single-level cell lookups are accepted, and that is
the canonicity check rather than a restriction: an ascending, unique subset of
a cell axis is a `CellLookup` again, and one that is neither is exactly what
DimensionalData degrades to a `Categorical` — so a cell dimension that is not
a cell lookup is a cell dimension that is no longer sorted and unique. A
`MultiOrderLookup` never reaches this: `_mixedaxis` claims it first.
"""
function _cellaxis(src)
    for d in DD.dims(src)
        lk = DD.val(d)
        lk isa Union{CellLookup,ChunkedCellLookup} || continue
        grid = levelgrid(system(lk), level(lk))
        return d, grid, _rawids(grid, lk)
    end
    return _nocellaxis(src)
end

function _rawids(grid, cells)
    I = idtype(grid)
    return I[convert(I, rawid(c)) for c in cells]
end

# The mixed-level cell dimension and its container, or `nothing` when the cube
# has none and the single-level path applies.
function _mixedaxis(src)
    for d in DD.dims(src)
        lk = DD.val(d)
        lk isa MultiOrderLookup && return d, parent(lk)
    end
    return nothing
end

@noinline function _nocellaxis(src)
    for d in DD.dims(src)
        lk = DD.val(d)
        eltype(lk) <: AbstractCellIndex || continue
        throw(DGGSFormatError(check=:noncanonical_cell_axis,
            declared=string(DD.name(d)), observed=nameof(typeof(DD.val(d))),
            detail="the $(DD.name(d)) dimension holds cell ids in a lookup this " *
                   "package uses for a cell axis that has stopped being sorted and " *
                   "unique — `reverse` is the everyday way to get one. A DGGS store " *
                   "names each cell once in canonical order, and one written out of " *
                   "order could not be opened again; sort the cube along its cell " *
                   "dimension first."))
    end
    throw(ArgumentError("dggwrite needs a cube with a cell dimension: none of " *
                        join(map(d -> string(DD.name(d)), DD.dims(src)), ", ") *
                        " carries a CellLookup or a ChunkedCellLookup."))
end

# ===========================================================================
# Encoding choice
# ===========================================================================

"""
    _encoding(spec, grid, cells) -> CellEncoding

`:auto` is [`RangesEncoding`](@ref) where the axis is eligible for it — sorted,
unique, one level — and [`DenseEncoding`](@ref) otherwise, which is what makes
it total. `:dense`, `:ranges` and `:implicit` are sugar over
`ENCODING_REGISTRY`, and an encoding instance takes the same eligibility check
its keyword would.
"""
function _encoding(spec::Symbol, grid, cells)
    spec === :auto && return write_eligible(RangesEncoding(), grid, cells) ?
                             RangesEncoding() : DenseEncoding()
    vocab = get(ENCODING_KEYWORDS, spec, nothing)
    vocab === nothing && throw(ArgumentError(
        "unknown encoding $(repr(spec)); it is :auto, or one of " *
        join(sort!([repr(k) for k in keys(ENCODING_KEYWORDS)]), ", ") *
        ", or a CellEncoding instance."))
    return _encoding(ENCODING_REGISTRY[vocab], grid, cells)
end

# An instance is held to the same eligibility as its keyword: `write_eligible`
# is part of the encoding interface, not of the keyword sugar.
function _encoding(enc::CellEncoding, grid, cells)
    write_eligible(enc, grid, cells) || throw(DGGSFormatError(
        check=:not_write_eligible, declared=encodingname(enc), observed=length(cells),
        detail="this axis cannot be written as $(encodingname(enc)): " *
               _ineligible(enc, grid)))
    return enc
end

_ineligible(::ImplicitEncoding, grid) =
    "an implicit axis is the whole of one level, in order, and this one is not."
_ineligible(::RangesEncoding, grid) =
    "ranges need sorted, unique ids from level $(level(grid))."
_ineligible(::CompactedEncoding, grid) =
    "compacted stores cells at MIXED levels, reached from a `MultiOrderLookup` " *
    "axis; a level-$(level(grid)) axis is written dense or as ranges."
_ineligible(::CellEncoding, grid) = "its `write_eligible` answered false."

"""
    _mixedencoding(spec) -> CompactedEncoding

The encoding for a `MultiOrderLookup` axis: `:auto` and `:compacted` both mean
[`CompactedEncoding`](@ref), the one layout that stores mixed levels. Any
single-level encoding, requested by name, is refused with the `expand` bridge.
"""
function _mixedencoding(spec)
    (spec === :auto || spec === :compacted || spec isa CompactedEncoding) &&
        return CompactedEncoding()
    known = spec isa CellEncoding ||
            (spec isa Symbol && haskey(ENCODING_KEYWORDS, spec))
    known || throw(ArgumentError(
        "unknown encoding $(repr(spec)); it is :auto, or one of " *
        join(sort!([repr(k) for k in keys(ENCODING_KEYWORDS)]), ", ") *
        ", or a CellEncoding instance."))
    label = spec isa Symbol ? ENCODING_KEYWORDS[spec] : _encodinglabel(spec)
    throw(DGGSFormatError(check=:mixed_level_axis, declared=label,
        observed=:MultiOrderLookup,
        detail="the cell axis is a `MultiOrderLookup` — cells at mixed " *
               "refinement levels — which only the `compacted` encoding " *
               "writes. Write it with `encoding = :auto`, or present the cube " *
               "at one level with `expand(A, level)` and write that as " *
               "`$label`."))
end

"""
    _coordinate(enc, grid, cells, merge)

What the encoding stores for the axis: the `(n, 2)` inclusive-range array, the
raw ids, or — for an implicit axis, which stores nothing — the length alone.
Computed once, and both written and read back through `_axis`.

The dense coordinate is the id vector `_cellaxis` already materialized,
handed on rather than copied: one vector of ids exists between the lookup and
the bytes on disk.
"""
_coordinate(::RangesEncoding, grid, cells, merge) = idranges(grid, cells; merge=merge)
_coordinate(::DenseEncoding, grid, cells, merge) = cells
_coordinate(::ImplicitEncoding, grid, cells, merge) = length(cells)
_coordinate(enc::CellEncoding, grid, cells, merge) = _nowritepath(enc)

"""
    _axis(enc, grid, coord, plan, n) -> ChunkedCellVector

The axis a reader would rebuild from `coord`, checked against the `n` cells that
went in. This is where a writer's mistake surfaces: an expansion-semantics
disagreement fails the normative count here rather than in someone else's reader.
"""
_axis(::RangesEncoding, grid, coord, plan, n) =
    cellaxis(RangesEncoding(), grid, coord; declared_length=n)

_axis(::DenseEncoding, grid, coord, plan, n) =
    cellaxis(DenseEncoding(), grid, coord; chunklength=plan.chunklength)

_axis(::ImplicitEncoding, grid, coord, plan, n) =
    cellaxis(ImplicitEncoding(), grid, coord)

_axis(enc::CellEncoding, grid, coord, plan, n) = _nowritepath(enc)

"""
    _nowritepath(enc)

Refuse an encoding this writer has no verbs for, by name.

The write pipeline is four private verbs on the encoding — `_coordinate`,
`_axis`, `_coordinate!` and `_coordinatename` — and a downstream encoding
that registers itself without them would otherwise fall off the end of dispatch
into a `MethodError` about a private function. This is the write-side twin of
the read half's `storedaxis` fallback: same check, same shape of message, from
the other direction.
"""
@noinline function _nowritepath(enc::CellEncoding)
    registered = sort!(collect(keys(ENCODING_REGISTRY)))
    throw(DGGSFormatError(check=:unsupported_encoding, observed=_encodinglabel(enc),
        detail="dggwrite writes the dense (`none`), `ranges`, `implicit` and " *
               "`compacted` layouts; `$(_encodinglabel(enc))` names an encoding it has no " *
               "write path for. A downstream encoding is written by implementing " *
               "this extension's `_coordinate`, `_axis`, `_coordinate!` and " *
               "`_coordinatename` for it. Registered encodings: " *
               join(registered, ", ") * "."))
end

# An encoding's own name where it has one, and its registry key where it does
# not: `encodingname` is the encoding's to define, and a type that has skipped
# the write verbs may well have skipped that too.
function _encodinglabel(enc::CellEncoding)
    for (name, registered) in ENCODING_REGISTRY
        registered === enc && return name
    end
    return string(nameof(typeof(enc)))
end

# ===========================================================================
# The chunk plan
# ===========================================================================

"""
    ChunkPlan(chunklength, ancestor_level, aligned)

What `chunks = :auto` decided, and how much of the coarse-ancestor property it
could keep.

**Zarr v2 chunks are uniform by format**: every chunk but the last holds exactly
`chunklength` cells. "Chunk on ancestor boundaries" is therefore a property of
one integer, not of individual boundaries, and `:auto` chooses that integer as
the largest whole number of level-`ancestor_level` subtree runs that fits under
the target.

  - **Guaranteed**: `chunklength` is a whole number of complete
    level-`ancestor_level` runs of the axis as written, so the FIRST chunk
    boundary is always a subtree boundary, and the persisted
    [`ChunkManifest`](@ref) describes the chunk grid exactly.
  - **`aligned`**: whether EVERY chunk boundary lands on a run boundary. It does
    whenever the runs are equal — full coverage of the coarse level, the case
    aggregation cares about, and the whole-level case always — and it can fail
    where they are not, because no uniform length lands on unequal boundaries.
    `:auto` does not trade the target away to chase it.

`ancestor_level` is `nothing` when no coarse level helped: an integer `chunks`,
a system whose subtrees are not contiguous in canonical order, or an axis whose
level-`(L-1)` runs already exceed the target. The chunk length is then the
target itself, clamped to the axis.

The target here is a CELL count: `dggwrite`'s `chunk_target` counts elements,
and the non-cell extents have already been divided out of it.
"""
struct ChunkPlan
    chunklength::Int
    ancestor_level::Union{Int,Nothing}
    aligned::Bool
end

function _chunkplan(chunks::Integer, grid, cells, target)
    chunks >= 1 || throw(ArgumentError(
        "a chunk length is at least one cell, not $chunks"))
    # Longer than the axis is not wrong, but it is not what was written either,
    # and the manifest must describe the chunks that exist.
    return ChunkPlan(min(Int(chunks), length(cells)), nothing, false)
end

function _chunkplan(chunks::Symbol, grid, cells, target)
    chunks === :auto || throw(ArgumentError(
        "chunks is :auto or a positive integer, not $(repr(chunks))"))
    n = length(cells)
    sys = system(grid)
    L = level(grid)
    plain = ChunkPlan(clamp(target, 1, n), nothing, false)
    (L < 1 || !has_sorted_subtrees(sys)) && return plain

    # Runs at level L-1 cost one `ancestor` per cell; every coarser level is then
    # a merge over the runs already found, so the whole descent is one pass.
    ends = _runends(grid, cells, L - 1, 1:n)
    best, bestlevel = nothing, nothing
    A = L - 1
    while _maxrun(ends) <= target
        best, bestlevel = ends, A
        A == 0 && break
        A -= 1
        ends = _runends(grid, cells, A, _runstarts(ends))
    end
    best === nothing && return plain

    # The largest whole number of runs under the target. `best` was only taken
    # while its longest run fitted, and `best[1]` is the first run's length, so
    # at least one always does.
    j = searchsortedlast(best, target)
    @assert j >= 1
    cl = best[j]
    return ChunkPlan(cl, bestlevel, _allaligned(best, cl, n))
end

# End positions of the maximal runs of cells sharing a level-`A` ancestor.
# `starts` names the positions to look at: every cell for the first pass, one
# representative per known run for each coarsening after it. The axis is raw
# ids, so the typed cell `ancestor` wants is put back together one at a time —
# a wrapper around an integer already in hand, and nothing is allocated.
function _runends(grid, cells, A::Int, starts)
    sys = system(grid)
    ends = Int[]
    previous = nothing
    for k in starts
        a = ancestor(sys, idcell(grid, cells[k]), A)
        previous === nothing || a == previous || push!(ends, k - 1)
        previous = a
    end
    isempty(cells) || push!(ends, length(cells))
    return ends
end

_runstarts(ends) = [i == 1 ? 1 : ends[i-1] + 1 for i in eachindex(ends)]

_maxrun(ends) = isempty(ends) ? 0 :
                maximum(i -> i == 1 ? ends[1] : ends[i] - ends[i-1], eachindex(ends))

# Whether every interior chunk boundary is also a run boundary. The last one is
# the axis length, which always is.
function _allaligned(ends, cl::Int, n::Int)
    boundaries = Set(ends)
    for b in cl:cl:(n-1)
        b in boundaries || return false
    end
    return true
end

# --- the mixed-level plan ---------------------------------------------------

_movchunkplan(chunks::Integer, mov, target) = _chunkplan(chunks, nothing, mov, target)

function _movchunkplan(chunks::Symbol, mov, target)
    chunks === :auto || throw(ArgumentError(
        "chunks is :auto or a positive integer, not $(repr(chunks))"))
    n = length(mov)
    plain = ChunkPlan(clamp(target, 1, n), nothing, false)
    floor = Int(first(levels(system(mov))))
    A = maximum(level, mov) - 1
    A < floor && return plain
    ends = _movrunends(mov, A, 1:n)
    best, bestlevel = nothing, nothing
    while _maxrun(ends) <= target
        best, bestlevel = ends, A
        A == floor && break
        A -= 1
        ends = _movrunends(mov, A, _runstarts(ends))
    end
    best === nothing && return plain
    j = searchsortedlast(best, target)
    @assert j >= 1
    cl = best[j]
    return ChunkPlan(cl, bestlevel, _allaligned(best, cl, n))
end

# Runs of cells sharing a level-`A` ancestor, on the container itself. A cell
# at level `A` or shallower is a complete subtree already and keys its own run.
function _movrunends(mov, A::Int, starts)
    sys = system(mov)
    ends = Int[]
    previous = nothing
    for k in starts
        c = mov[k]
        a = level(c) <= A ? c : ancestor(sys, c, A)
        previous === nothing || a == previous || push!(ends, k - 1)
        previous = a
    end
    isempty(mov) || push!(ends, length(mov))
    return ends
end

# The compacted manifest is keyed at the reference level: per chunk, the
# interval start of its first cell and the interval stop of its last. Both
# ascend because the container is sorted by interval start and disjoint.
function _movmanifest(mov::MultiOrderVector, chunklength::Int)
    n = length(mov)
    cl = max(chunklength, 1)
    nc = cld(n, cl)
    firstids = Vector{Int}(undef, nc)
    lastids = Vector{Int}(undef, nc)
    lengths = Vector{Int}(undef, nc)
    offsets = Vector{Int}(undef, nc)
    for c in 1:nc
        lo = (c - 1) * cl + 1
        hi = min(c * cl, n)
        firstids[c] = mov.starts[lo]
        lastids[c] = mov.stops[hi]
        lengths[c] = hi - lo + 1
        offsets[c] = lo - 1
    end
    return ChunkManifest(firstids, lastids, lengths, offsets, cl)
end

# ===========================================================================
# The arrays
# ===========================================================================

# One `ArrayEntry` (what the conventions stamp) plus what it takes to create the
# array and fill it. `source` is either a materialized array — the coordinate,
# the dimension values and the manifest, all kilobytes — or a `CellStream`,
# which is the layer as it was handed in.
struct ArrayWrite{S}
    entry::ArrayEntry
    source::S
    chunks::Tuple{Vararg{Int}}
end

# A layer still in whatever array it arrived in: a lazy `ZArray` straight out of
# `dggread` as readily as an `Array`. `celldim` is the position of the cell
# dimension in it and `perm` the permutation that puts cells first, `nothing`
# where they already are. Nothing bigger than one chunk is ever taken from it.
struct CellStream{A,P}
    data::A
    celldim::Int
    perm::P
end

# The coordinate, the dimension values and the manifest go over in one piece.
_fill!(z, values::AbstractArray, chunklength::Int) =
    (z[ntuple(_ -> Colon(), ndims(values))...] = values; z)

# A layer goes over a chunk at a time, along the cell axis and on the chunk
# boundaries the array was created with, so a lazy source reads chunk-sized
# pieces and the whole axis is never in memory at once. The premise is tens of
# millions of cells: materializing a layer to write it costs the axis twice.
function _fill!(z, s::CellStream, chunklength::Int)
    n = size(s.data, s.celldim)
    rest = ntuple(_ -> Colon(), ndims(s.data) - 1)
    for lo in 1:chunklength:n
        r = lo:min(lo + chunklength - 1, n)
        z[r, rest...] = _block(s, r)
    end
    return z
end

function _block(s::CellStream, r)
    block = s.data[ntuple(i -> i == s.celldim ? r : Colon(), ndims(s.data))...]
    return s.perm === nothing ? block : permutedims(block, s.perm)
end

# The names this writer owns, whatever the encoding chosen: the dense
# coordinate, the range array, and the manifest sidecar. They are reserved
# together rather than one encoding at a time, so that renaming a layer is not
# something `encoding = :dense` asks for and `:auto` does not.
const RESERVED_ARRAYS = (CELL_IDS_ARRAY, CELL_RANGES_ARRAY, MANIFEST_ARRAY,
    COMPACTED_LEVELS_ARRAY)

function _arrayplan(enc, coord, layers, celldim, plan, manifest, desc)
    _checklayernames(layers)
    out = ArrayWrite[]
    _coordinate!(out, enc, coord, plan)
    seen = Set{String}()
    for (name, A) in layers
        push!(out, _layerwrite(name, A, celldim, plan))
        for d in DD.otherdims(A, celldim)
            key = string(DD.name(d))
            key in seen && continue
            push!(seen, key)
            w = _dimcoordinate(d, key)
            w === nothing || push!(out, w)
        end
    end
    # Every store gets a manifest, including the ranges and implicit ones whose
    # axis is arithmetic and whose reader never opens the sidecar. That is
    # deliberate: it is a kilobyte-scale table that records the chunk grid and
    # the level and grid the axis was validated at, which is what a later
    # encoding, an aggregation reading chunk boundaries, or a reader of a
    # partially rewritten store would otherwise have to recompute or assume.
    push!(out, _manifestwrite(manifest, plan, desc))
    sort!(out; by=a -> a.entry.name)
    _checkunique(out)
    return out
end

# One Zarr array per name, checked before the destination is touched. `zcreate`
# would find the collision too, but only on the second array of that name — by
# which time the group has been created and stamped and some of its arrays
# written, and what is on disk is neither the old store nor the new one.
@noinline function _checklayernames(layers)
    taken = sort!(String[String(n) for (n, _) in layers if String(n) in RESERVED_ARRAYS])
    isempty(taken) && return nothing
    throw(DGGSFormatError(check=:reserved_array_name, declared=taken,
        observed=collect(RESERVED_ARRAYS),
        detail="the layer " * join(taken, ", ") * " would take an array name this " *
               "writer owns: a DGGS store keeps its cell coordinate and its chunk " *
               "manifest in " * join(RESERVED_ARRAYS, ", ") * ", whichever encoding " *
               "it is written in. Rename the layer."))
end

@noinline function _checkunique(out)
    for i in 2:length(out)
        name = out[i].entry.name
        name == out[i-1].entry.name || continue
        throw(DGGSFormatError(check=:duplicate_array_name, declared=name,
            observed=String[a.entry.name for a in out],
            detail="two arrays of this store would be called `$name`. A layer and " *
                   "a dimension of one name are one Zarr array and not two; rename " *
                   "the layer, or the dimension."))
    end
    return nothing
end

# The cell coordinate. Grid attributes are stamped onto it later by the flat
# conventions, so it carries only its dimension names here. Zarr.jl reverses
# shape between the JSON and Julia, so the `(n, 2)` range array is handed over
# transposed to land as `(n, 2)` in the store, as the published stores have it.
function _coordinate!(out, ::RangesEncoding, R::AbstractMatrix, plan)
    return _push!(out, CELL_RANGES_ARRAY, permutedims(R), (size(R, 1), 2),
        RANGES_DIMS, (2, size(R, 1)))
end

_coordinate!(out, ::DenseEncoding, ids::AbstractVector, plan) =
    _push!(out, CELL_IDS_ARRAY, ids, (length(ids),), [SPATIAL_DIMENSION],
        (plan.chunklength,))

# An implicit axis stores nothing: position IS the cell.
_coordinate!(out, ::ImplicitEncoding, n::Integer, plan) = out

# The compacted coordinate is the id column on the spatial dimension plus the
# level column on a dimension of its own, so the levels never read as a data
# variable.
function _coordinate!(out, ::CompactedEncoding, coord::Tuple, plan)
    lv, ids = coord
    _push!(out, CELL_IDS_ARRAY, ids, (length(ids),), [SPATIAL_DIMENSION],
        (plan.chunklength,))
    return _push!(out, COMPACTED_LEVELS_ARRAY, lv, (length(lv),),
        [COMPACTED_LEVELS_ARRAY], (plan.chunklength,))
end

_coordinate!(out, enc::CellEncoding, coord, plan) = _nowritepath(enc)

function _push!(out, name, values, shape, dims, chunks)
    entry = ArrayEntry(name=name,
        attrs=Dict{String,Any}(ARRAY_DIMENSIONS => copy(dims)),
        shape=shape, eltype=eltype(values), dims=copy(dims))
    push!(out, ArrayWrite(entry, values, chunks))
    return out
end

# Cells become the LAST Zarr dimension and so the fastest-varying one, which is
# what makes a chunk of cells contiguous. Zarr.jl reverses shape and chunks
# between the JSON and Julia, so cells go first here — as a permutation applied
# to each block rather than to the cube, which would materialize it.
function _layerwrite(name, A, celldim, plan)
    ds = DD.dims(A)
    cd = findfirst(d -> DD.name(d) === DD.name(celldim), ds)
    perm = (cd, filter(!=(cd), ntuple(identity, length(ds)))...)
    juliadims = [_dimname(ds[i], celldim) for i in perm]
    shape = map(i -> size(A, i), perm)
    return ArrayWrite(
        ArrayEntry(name=String(name), attrs=_layerattrs(A, reverse(juliadims)),
            shape=reverse(shape), eltype=eltype(A), dims=reverse(juliadims)),
        CellStream(DD.data(A), cd, cd == 1 ? nothing : perm),
        (plan.chunklength, Base.tail(shape)...))
end

# The producer's attributes ride through to disk, with the ones this writer
# generates stamped OVER them: a `_ARRAY_DIMENSIONS` carried in from another
# layout would describe this array wrongly.
function _layerattrs(A, dims)
    attrs = _attrs(DD.metadata(A))
    attrs[ARRAY_DIMENSIONS] = dims
    return attrs
end

_dimname(d, celldim) = DD.name(d) === DD.name(celldim) ? SPATIAL_DIMENSION :
                       string(DD.name(d))

# A non-cell dimension's own values, where they are something Zarr holds
# directly. Anything else (a DateTime axis, a lookup of structs) becomes a bare
# Zarr dimension with no coordinate array; encoding those is not in v1.
function _dimcoordinate(d, key)
    lk = DD.val(d)
    lk isa DD.Lookups.NoLookup && return nothing
    values = collect(lk)
    eltype(values) <: Union{Real,AbstractString} || return nothing
    entry = ArrayEntry(name=key,
        attrs=Dict{String,Any}(ARRAY_DIMENSIONS => [key]),
        shape=(length(values),), eltype=eltype(values), dims=[key])
    return ArrayWrite(entry, values, (max(length(values), 1),))
end

# The manifest as the design's `(n_chunks, 2)` sidecar: first and last id per
# chunk, with the chunk length and the axis length that make it interpretable,
# and — when the plan came from a coarse level — the ancestor level it grouped
# by and whether every boundary really landed on one of its subtrees. Its
# dimensions are its own, so it is invisible to every convention and to
# `arrays_on` — an extra variable an xarray reader ignores.
#
# `writer`/`format`/`validated` are what a consumer decides how far to trust a
# stale sidecar by: `validated` is always `"strict"`, because the axis is
# rebuilt through the reader's own `cellaxis` before a byte is committed.
# `grid`/`level` say what it was validated AGAINST — a claim of strictness is
# empty without them, since the same ids are a clean axis at the level they were
# written at and name no cell at all at any other.
function _manifestwrite(manifest, plan, desc)
    rows = permutedims(hcat(manifest.firstids, manifest.lastids))
    n = size(rows, 2)
    marker = Dict{String,Any}(
        "writer" => MANIFEST_WRITER,
        "format" => MANIFEST_FORMAT,
        "validated" => MANIFEST_VALIDATED,
        "grid" => desc.gridname,
        "level" => desc.level,
        "spatial_dimension" => SPATIAL_DIMENSION,
        "chunk_length" => manifest.chunklength,
        "length" => length(manifest))
    if plan.ancestor_level !== nothing
        marker["ancestor_level"] = plan.ancestor_level
        marker["ancestor_aligned"] = plan.aligned
    end
    entry = ArrayEntry(name=MANIFEST_ARRAY,
        attrs=Dict{String,Any}(ARRAY_DIMENSIONS => copy(MANIFEST_DIMS),
            MANIFEST_MARKER => marker),
        shape=(n, 2), eltype=eltype(rows), dims=copy(MANIFEST_DIMS))
    return ArrayWrite(entry, rows, (2, n))
end

# ===========================================================================
# Layers and the description
# ===========================================================================

_layers(A::DD.AbstractDimArray, celldim) = _checklayers([(_layername(A), A)], celldim)
_layers(s::DD.AbstractDimStack, celldim) =
    _checklayers([(String(k), s[k]) for k in keys(DD.layers(s))], celldim)

function _checklayers(layers, celldim)
    isempty(layers) && throw(ArgumentError("dggwrite has no layers to write."))
    for (name, A) in layers
        DD.hasdim(A, celldim) || throw(ArgumentError(
            "layer $(repr(name)) has no cell dimension; every layer of a DGGS " *
            "store shares one cell axis."))
    end
    return sort!(layers; by=first)
end

function _layername(A)
    n = DD.name(A)
    n isa Symbol && n !== Symbol("") && return String(n)
    return "layer"
end

function _description(sys, lev, enc, layers)
    name, ref = _gridname(sys)
    return StoreDescription(gridname=name, system=ref.system, idscheme=ref.idscheme,
        level=lev, encoding=enc, coordinate=_coordinatename(enc),
        spatial_dimension=SPATIAL_DIMENSION,
        variables=String[n for (n, _) in layers])
end

_coordinatename(::RangesEncoding) = CELL_RANGES_ARRAY
_coordinatename(::DenseEncoding) = CELL_IDS_ARRAY
_coordinatename(::ImplicitEncoding) = nothing
_coordinatename(::CompactedEncoding) = CELL_IDS_ARRAY
_coordinatename(enc::CellEncoding) = _nowritepath(enc)

# The reference table read backwards: a grid name pins the id packing, so a
# system with no registered name has no store spelling either.
function _gridname(sys)
    for (name, ref) in DGG.GRID_REFERENCE
        ref.system == sys && return name, ref
    end
    throw(DGGSFormatError(check=:unknown_grid_name, observed=nameof(typeof(sys)),
        detail="$(nameof(typeof(sys))) has no canonical store name; registered " *
               "names are " * join(sort!(collect(keys(DGG.GRID_REFERENCE))), ", ") *
               ". Add one with `register_grid!(name, GridReference(...))` before " *
               "writing."))
end

end # module DGGSZarrWrite
