# `dggwrite` follows cube → cell axis → encoding → chunk plan → snapshot →
# bytes. Core IO code supplies encoding, convention and manifest semantics; this
# extension plans Zarr arrays and writes them. Planning the complete snapshot
# before creation keeps store attributes aligned with the arrays. A submodule
# isolates these helpers from the read half of the extension.

module DGGSZarrWrite

import DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGridsZarrExt
using ..DiscreteGlobalGridsZarrExt: storeidentifier,
    MANIFEST_MARKER, MANIFEST_WRITER, MANIFEST_FORMAT, MANIFEST_VALIDATED,
    COMPACTED_LEVELS_ARRAY
using DiscreteGlobalGrids: AbstractCellIndex, ArrayEntry,
    AbstractCellLookup, CellEncoding, CellLookup, ChunkManifest, CompactedEncoding,
    DGGSFormatError, DEFAULT_WRITE_CONVENTIONS, DGGSConvention, DenseEncoding,
    ENCODING_REGISTRY, ImplicitEncoding, MultiOrderLookup, MultiOrderVector,
    RangesEncoding, StoreDescription, StoreSnapshot, XdggsConvention,
    ancestor, cellaxis, chunkmanifest, encodingname, has_sorted_subtrees,
    idcell, idranges, idtype, level, levels, levelgrid, rawid, system,
    write_eligible
import DimensionalData as DD
import Zarr

# Data variables share `spatial_dimension`; `coordinate` names the array that
# encodes it. A ranges coordinate has its own `(n, 2)` dimensions.

const CELL_IDS_ARRAY = "cell_ids"
const CELL_RANGES_ARRAY = "cell_id_ranges"
const SPATIAL_DIMENSION = "cell_ids"
const RANGES_DIMS = ["ranges", "bounds"]
const MANIFEST_ARRAY = "cell_chunk_manifest"
const MANIFEST_DIMS = ["chunks", "bounds"]
const ARRAY_DIMENSIONS = "_ARRAY_DIMENSIONS"

"""
    DEFAULT_CHUNK_TARGET

The element-count target for `chunks = :auto`, including cells and every other
dimension. `chunk_target` overrides this default.
"""
const DEFAULT_CHUNK_TARGET = 1_000_000

# Keyword aliases resolve through the same encoding registry as store metadata.
const ENCODING_KEYWORDS = Dict(:dense => "none", :ranges => "ranges",
    :implicit => "implicit", :compacted => "compacted")

const REMOTE_SCHEMES = ("gs://", "s3://", "http://", "https://", "az://", "abfs://")

const Cube = Union{DD.AbstractDimArray,DD.AbstractDimStack}

"""
    dggwrite(dest, stack_or_array; encoding = :auto,
             conventions = DEFAULT_WRITE_CONVENTIONS, chunks = :auto,
             merge = :step, chunk_target = DEFAULT_CHUNK_TARGET) -> dest

Write a `DimStack` or `DimArray` over a cell axis to a **Zarr v2 directory
store**, consolidated metadata included. `dest` is a local path or an open,
writable `Zarr.ZGroup`. String URL destinations are rejected; write locally and
upload, or pass an already-open writable remote group.

The cell dimension carries one of two lookup families:

  - `AbstractCellLookup` represents a sorted, unique, single-level axis.
    Operations such as `reverse` degrade it to a `Categorical`, which the writer
    rejects as noncanonical.
  - `MultiOrderLookup` represents a mixed-level axis and writes as compacted.

`encoding` selects the cell-axis layout:

  - `:auto` selects compacted for a mixed-level axis, ranges for an eligible
    single-level axis, and dense otherwise.
  - `:dense` writes every id for broad reader compatibility.
  - `:ranges` writes the compact single-level range representation.
  - `:implicit` writes a complete level with no cell coordinate.
  - `:compacted` writes `cell_ids` and `cell_levels` as aligned columns under
    `refinement_level: null`.

Single-level encodings require `expand(A, level)` to present mixed-level data
at one level.

Other options are:

  - `merge = :step` merges integer-adjacent ids for structural-reader
    compatibility. `merge = :rank` merges consecutive cells for fewer rows and
    requires a rank-aware reader; see [`idranges`](@ref).
  - `chunks = :auto` groups complete coarse-ancestor subtree runs near
    `chunk_target`. An integer fixes the chunk length in cells.
  - `chunk_target` counts every element in a chunk, including non-cell
    dimensions.
  - `conventions` stamps a single-level store with `zarr-conventions/dggs` and
    xdggs by default. Compacted stores carry only the compatible DGGS
    convention metadata because xdggs describes a single-level coordinate.

The writer persists each chunk's first and last id in an `(n_chunks, 2)`
sidecar, allowing readers to rebuild the chunk grid without an axis scan.

**Attributes.** Layer metadata becomes array attributes, and
`metadata["attrs"]` becomes group attributes. Convention-generated keys take
precedence over producer values. The encoding regenerates the cell coordinate,
so it carries fresh layout attributes.

**Round-trip normalizations.** Layers return in alphabetical order, and each
layer's metadata gains the writer's `_ARRAY_DIMENSIONS` attribute.

The writer raises before stamping a `ZGroup` that already contains a planned
array name.
"""
function DGG.dggwrite(dest::AbstractString, src::Cube; layout::Symbol=:cells, kw...)
    path = String(dest)
    _reject_remote(path)
    layout === :cells || return (_otherlayout(layout, path, src; kw...); dest)
    # `zgroup` refuses a store that is not empty, so a path needs no name guard.
    _write(path, (attrs, names) -> Zarr.zgroup(path; attrs=attrs), src; kw...)
    return dest
end

function DGG.dggwrite(dest::Zarr.ZGroup, src::Cube; layout::Symbol=:cells, kw...)
    dest.writeable || throw(ArgumentError(
        "dggwrite needs a writeable group; this one was opened read-only."))
    layout === :cells || return (_otherlayout(layout, dest, src; kw...); dest)
    _write(storeidentifier(dest), (attrs, names) -> _stamp(dest, attrs, names),
        src; kw...)
    return dest
end

# `layout` selects the store shape; `encoding` selects its cell-coordinate shape.
# The subzone module is referenced by name because it loads after this file.
@noinline function _otherlayout(layout::Symbol, dest, src; kw...)
    layout === :subzones && return DiscreteGlobalGridsZarrExt.DGGSZarrSubzones.write_subzones(
        dest, src; kw...)
    throw(ArgumentError(
        "dggwrite writes the `:cells` layout — one cell dimension, the default — " *
        "and the `:subzones` layout, which is the two-dimensional " *
        "ancestor-subzone store and takes an `ancestor_level`. " *
        "$(repr(layout)) is neither."))
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

# Validate every destination name before the first irreversible attribute write.
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
    # Reject invalid alignment before creating a partial store.
    cellaxis(enc, system(mov), lv, ids; declared_length=length(mov))
    manifest = _movmanifest(mov, plan.chunklength)
    desc = _description(system(mov), nothing, enc, layers)
    arrays = _arrayplan(enc, (lv, ids), layers, celldim, plan, manifest, desc)
    return _commit(opengroup, identifier, src, desc, conventions, arrays;
        reference_level=DGG.reference_level(mov))
end

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

# xdggs readers decode coordinates at one declared level.
_stampable(::DGGSConvention, desc) = true
_stampable(::XdggsConvention, desc) = desc.level !== nothing

# v1 of zarr-conventions/dggs lacks a level-column key, so publish the extension.
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

# Divide the element target by the widest non-cell extent because every layer
# and the cell axis share one chunk length.
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

# Restore producer group attributes before conventions stamp authoritative keys.
_groupattrs(src) = _attrs(get(_attrs(DD.metadata(src)), "attrs", nothing))

_attrs(x) = Dict{String,Any}()
_attrs(md::AbstractDict) = Dict{String,Any}(String(k) => deepcopy(v) for (k, v) in md)
_attrs(md::NamedTuple) = Dict{String,Any}(String(k) => deepcopy(v) for (k, v) in pairs(md))
_attrs(md::DD.Metadata) = _attrs(DD.val(md))

"""
    _cellaxis(src) -> (dim, grid, ids)

Return the cube's cell dimension, its complete level grid and its raw ids.

Materializing raw ids once supplies eligibility checks, chunk planning and the
dense coordinate without a second axis copy. `idcell` reconstructs an individual
typed cell when needed.

[`AbstractCellLookup`](@ref) guarantees a sorted, unique single-level axis.
DimensionalData represents noncanonical cell axes as `Categorical`, which this
writer rejects. `_mixedaxis` handles `MultiOrderLookup` before this method runs.
"""
function _cellaxis(src)
    for d in DD.dims(src)
        lk = DD.val(d)
        lk isa AbstractCellLookup || continue
        grid = levelgrid(system(lk), level(lk))
        return d, grid, _rawids(grid, lk)
    end
    return _nocellaxis(src)
end

function _rawids(grid, cells)
    I = idtype(grid)
    return I[convert(I, rawid(c)) for c in cells]
end

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
                        " carries a cell lookup."))
end

"""
    _encoding(spec, grid, cells) -> CellEncoding

Resolve the encoding for a single-level axis. `:auto` selects
[`RangesEncoding`](@ref) for an eligible sorted, unique axis and
[`DenseEncoding`](@ref) otherwise. Keyword aliases resolve through
`ENCODING_REGISTRY`; direct instances use the same eligibility check.
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

# Apply the encoding's eligibility rule equally to instances and keyword names.
function _encoding(enc::CellEncoding, grid, cells)
    write_eligible(enc, grid, cells) || throw(DGGSFormatError(
        check=:not_write_eligible, declared=encodingname(enc), observed=length(cells),
        detail="this axis cannot be written as $(encodingname(enc)): " *
               _ineligible(enc, grid)))
    return enc
end

_ineligible(::ImplicitEncoding, grid) =
    "an implicit axis requires every cell of one level in canonical order."
_ineligible(::RangesEncoding, grid) =
    "ranges require sorted, unique ids from level $(level(grid))."
_ineligible(::CompactedEncoding, grid) =
    "compacted requires a mixed-level `MultiOrderLookup`; write this " *
    "level-$(level(grid)) axis as dense or ranges."
_ineligible(::CellEncoding, grid) = "its `write_eligible` method returned false."

"""
    _mixedencoding(spec) -> CompactedEncoding

Resolve the encoding for a `MultiOrderLookup` axis. `:auto` and `:compacted`
select [`CompactedEncoding`](@ref); single-level encodings report how to use
`expand` first.
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
        detail="a `MultiOrderLookup` axis requires the `compacted` encoding. " *
               "Use `encoding = :auto`, or present the cube at one level with " *
               "`expand(A, level)` before writing it as `$label`."))
end

"""
    _coordinate(enc, grid, cells, merge)

Compute the axis representation once:

  - ranges use an `(n, 2)` inclusive-range array;
  - dense uses the materialized raw-id vector directly; and
  - implicit uses the axis length.

`_axis` validates the same representation before writing.
"""
_coordinate(::RangesEncoding, grid, cells, merge) = idranges(grid, cells; merge=merge)
_coordinate(::DenseEncoding, grid, cells, merge) = cells
_coordinate(::ImplicitEncoding, grid, cells, merge) = length(cells)
_coordinate(enc::CellEncoding, grid, cells, merge) = _nowritepath(enc)

"""
    _axis(enc, grid, coord, plan, n) -> ChunkedCellVector

Rebuild and validate the reader-visible axis from `coord`. The declared count
`n` catches expansion or counting disagreements before the store is committed.
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

Raise a named unsupported-encoding error when an encoding lacks any required
write verb: `_coordinate`, `_axis`, `_coordinate!` or `_coordinatename`.
"""
@noinline function _nowritepath(enc::CellEncoding)
    registered = sort!(collect(keys(ENCODING_REGISTRY)))
    throw(DGGSFormatError(check=:unsupported_encoding, observed=_encodinglabel(enc),
        detail="dggwrite writes the dense (`none`), `ranges`, `implicit` and " *
               "`compacted` layouts; `$(_encodinglabel(enc))` names an " *
               "encoding it has no " *
               "write path for. A downstream encoding is written by implementing " *
               "this extension's `_coordinate`, `_axis`, `_coordinate!` and " *
               "`_coordinatename` for it. Registered encodings: " *
               join(registered, ", ") * "."))
end

# Prefer a registry key so incomplete downstream encodings still have a useful label.
function _encodinglabel(enc::CellEncoding)
    for (name, registered) in ENCODING_REGISTRY
        registered === enc && return name
    end
    return string(nameof(typeof(enc)))
end

"""
    WriteChunkPlan(chunklength, ancestor_level, aligned)

Record the chunk length selected by `chunks = :auto` and its ancestor alignment.

Zarr v2 uses a uniform `chunklength` except for the final chunk. Automatic
planning chooses the largest whole number of level-`ancestor_level` subtree runs
within the cell target.

  - `chunklength` always contains complete runs from that ancestor level, so the
    first chunk boundary aligns and the [`ChunkManifest`](@ref) describes the
    exact grid.
  - `aligned` reports whether every interior chunk boundary aligns. Equal run
    lengths guarantee this property; unequal runs may prevent it.
  - `ancestor_level = nothing` records a fixed integer request, a system without
    contiguous subtrees, or ancestor runs larger than the target. In these cases
    the clamped target becomes the chunk length.

The target here counts cells because `_celltarget` has already divided out the
non-cell extents.
"""
struct WriteChunkPlan
    chunklength::Int
    ancestor_level::Union{Int,Nothing}
    aligned::Bool
end

function _chunkplan(chunks::Integer, grid, cells, target)
    chunks >= 1 || throw(ArgumentError(
        "a chunk length is at least one cell, not $chunks"))
    # Longer than the axis is not wrong, but it is not what was written either,
    # and the manifest must describe the chunks that exist.
    return WriteChunkPlan(min(Int(chunks), length(cells)), nothing, false)
end

function _chunkplan(chunks::Symbol, grid, cells, target)
    chunks === :auto || throw(ArgumentError(
        "chunks is :auto or a positive integer, not $(repr(chunks))"))
    n = length(cells)
    sys = system(grid)
    L = level(grid)
    plain = WriteChunkPlan(clamp(target, 1, n), nothing, false)
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
    return WriteChunkPlan(cl, bestlevel, _allaligned(best, cl, n))
end

# `starts` samples every cell on the first pass and one representative per run on
# later coarsening passes; typed wrappers are reconstructed one at a time.
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

# Check whether every interior uniform-chunk boundary ends an ancestor run.
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
    plain = WriteChunkPlan(clamp(target, 1, n), nothing, false)
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
    return WriteChunkPlan(cl, bestlevel, _allaligned(best, cl, n))
end

# Cells at or above level `A` key their own complete-subtree run.
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

# Reference-level interval endpoints give each compacted chunk searchable bounds.
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

struct ArrayWrite{S}
    entry::ArrayEntry
    source::S
    chunks::Tuple{Vararg{Int}}
end

# Retain the source so writing materializes only one cell chunk at a time.
struct CellStream{A,P}
    data::A
    celldim::Int
    perm::P
end

# The coordinate, the dimension values and the manifest go over in one piece.
_fill!(z, values::AbstractArray, chunklength::Int) =
    (z[ntuple(_ -> Colon(), ndims(values))...] = values; z)

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

# Reserve every writer-owned name independently of the selected encoding.
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

# Detect array-name collisions before creating or stamping the destination.
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

# Transpose the ranges coordinate to preserve its `(n, 2)` on-disk Zarr shape.
function _coordinate!(out, ::RangesEncoding, R::AbstractMatrix, plan)
    return _push!(out, CELL_RANGES_ARRAY, permutedims(R), (size(R, 1), 2),
        RANGES_DIMS, (2, size(R, 1)))
end

_coordinate!(out, ::DenseEncoding, ids::AbstractVector, plan) =
    _push!(out, CELL_IDS_ARRAY, ids, (length(ids),), [SPATIAL_DIMENSION],
        (plan.chunklength,))

_coordinate!(out, ::ImplicitEncoding, n::Integer, plan) = out

# A separate dimension keeps the level column out of data-variable discovery.
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

# Put cells first in Julia so reversed Zarr metadata makes them the contiguous
# final on-disk dimension; permute each block to preserve lazy streaming.
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

# Stamp authoritative dimensions over producer attributes from an earlier layout.
function _layerattrs(A, dims)
    attrs = _attrs(DD.metadata(A))
    attrs[ARRAY_DIMENSIONS] = dims
    return attrs
end

_dimname(d, celldim) = DD.name(d) === DD.name(celldim) ? SPATIAL_DIMENSION :
                       string(DD.name(d))

# Write directly representable dimension values; use a bare dimension for other types.
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

# The independent `(n_chunks, 2)` sidecar records chunk bounds, axis geometry,
# validation provenance and optional ancestor alignment. Rebuilding the axis
# through `cellaxis` before commit justifies its `validated = "strict"` marker.
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

# Reverse the grid registry because a canonical store name also fixes id packing.
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
