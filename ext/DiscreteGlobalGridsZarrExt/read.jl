# `dggread` follows open → snapshot → describe → cell axis → `DimStack`.
# Core IO code owns format semantics; this extension supplies Zarr access and
# adds store and convention context to errors at the boundary.

import DimensionalData as DD
using DimensionalData: Lookups

using DiscreteGlobalGrids: StoreDescription, DGGSFormatError, Cells,
    ChunkManifest, ChunkedCellLookup, MultiOrderLookup, with_store_context,
    arrays_on, getarray
using DiscreteGlobalGrids.Encodings: CellEncoding, DenseEncoding, RangesEncoding,
    ImplicitEncoding, CompactedEncoding, cellaxis, encodingname, idtype

"""
    LAZY_SAMPLES

How many ids per chunk `validate = :lazy` checks for existence. Both ends of
the chunk — the ids the manifest publishes — plus an even interior spread. A
phantom between two samples survives, which is what the opt-out buys.
"""
const LAZY_SAMPLES = 8

"What `dggread` will open: an already-open group, a store, or a path or URL."
const StoreLike = Union{Zarr.ZGroup,Zarr.AbstractStore,AbstractString}

# Shared marker constants keep the reader and writer on one manifest contract.

"""
    MANIFEST_MARKER

The identifying attribute on a chunk-manifest sidecar. Its fields record:

  - `writer`, `format` and `validated` for provenance and validation strength;
  - `grid` and `level` for the axis used during validation; and
  - `spatial_dimension`, `chunk_length` and `length` for the described chunk grid.

`persistedmanifest` trusts a sidecar only when all eight fields match the open
store.
"""
const MANIFEST_MARKER = "dggs_chunk_manifest"

"The only `writer` this reader builds an axis from without scanning it."
const MANIFEST_WRITER = "DiscreteGlobalGrids.jl"

"The sidecar layout version this reader understands."
const MANIFEST_FORMAT = 1

"The `validated` level that attests every id was checked when it was written."
const MANIFEST_VALIDATED = "strict"

"""
    COMPACTED_LEVELS_ARRAY

The second column of a `compacted` cell axis: one integer level per stored
cell, aligned with the id coordinate. It sits on a dimension of its own so
that no convention reads it as a data variable.
"""
const COMPACTED_LEVELS_ARRAY = "cell_levels"

"""
    dggread(store; vars = All(), lazy = true, validate = :strict,
            conventions = CONVENTION_REGISTRY, description = nothing) -> DimStack
    dggread(store, var::Symbol; kwargs...) -> DimArray

Read a DGGS store into plain DimensionalData with one `Cells` dimension shared
by every layer. The dimension carries one of two lookups:

  - [`ChunkedCellLookup`](@ref) resolves cells on a stored single-level axis
    without scanning it.
  - [`MultiOrderLookup`](@ref) represents the mixed-level cells of a
    `compacted` store.

Every other store dimension, such as time or bands, becomes an ordinary `Dim`
and uses the values of its like-named coordinate array when present.

`store` accepts a `Zarr.ZGroup`, `Zarr.AbstractStore`, local path or URL.
`gs://BUCKET/PATH` maps to the public
`https://storage.googleapis.com/BUCKET/PATH` endpoint. `s3://` support comes from
Zarr's AWSS3 extension, loaded with `using AWSS3`.

  - `vars`: `All()` reads every data variable; symbols select specific variables.
    An unknown name raises an error that lists the available variables.
  - `lazy`: `true` keeps the store's chunked arrays; `false` materializes them.
  - `validate`: controls cell-id validation.
    - `:strict` checks every stored id unless a trusted package-written manifest
      supplies the axis.
    - `:lazy` samples `LAZY_SAMPLES` ids per chunk.
    - `:scan` ignores a manifest and checks every id.

    Every mode checks sortedness, uniqueness and declared lengths when it scans.
    A trusted manifest instead verifies each chunk's first id, last id and
    length when that chunk is read; see `persistedmanifest`.
  - `conventions`: the conventions to try, in order.
  - `ancestors`: restricts an ancestor-subzone store to selected ancestor cells
    or column indices. The default returns the complete level, including unwritten
    columns as fill values.
  - `description`: supplies a [`StoreDescription`](@ref) and skips detection.
    The id and closed-form count checks still run, enabling attribute-free stores.

The stack's `metadata` carries the provenance a value-identical rewrite needs:

| key | |
|---|---|
| `"source"` | the store URL or path |
| `"conventions"` | the conventions that fired, in order |
| `"encoding"` | the cell-axis layout, as the store spells it |
| `"attrs"` | the group attributes verbatim |
| `"description"` | the [`StoreDescription`](@ref) everything was read through |

Each layer keeps its own array attributes as its metadata.

**Ancestor-subzone stores** use their own attributes and read path. Their layers
retain a `Cells` dimension backed by a lazy [`SubzoneCellArray`](@ref), with one
subtree per chunk. Their stack metadata records a [`SubzoneLayout`](@ref) under
`"layout"`.
"""
function DiscreteGlobalGrids.dggread(store::StoreLike; vars=DD.All(), lazy::Bool=true,
    validate::Symbol=:strict, conventions=DiscreteGlobalGrids.CONVENTION_REGISTRY,
    description::Union{StoreDescription,Nothing}=nothing, ancestors=nothing)
    samples = validation_samples(validate)
    group, identifier = opengroup(store)
    return with_store_context(identifier) do
        snap = snapshot(group; identifier)
        # The two-dimensional layout is recognized before the conventions are
        # asked, and by its own attribute rather than by one of theirs: its data
        # arrays share no cell dimension, so there is nothing for a convention to
        # describe and no `StoreDescription` that could describe it.
        if DiscreteGlobalGrids.issubzonestore(snap.attrs)
            description === nothing || throw(ArgumentError(
                "`description` asserts a one-dimensional cell axis; this store is " *
                "the ancestor-subzone layout, whose shape is in its own attributes."))
            return DGGSZarrSubzones.assemble(group, snap, identifier, vars, lazy, ancestors)
        end
        ancestors === nothing || throw(ArgumentError(
            "`ancestors` selects the columns of an ancestor-subzone store, and " *
            "this store is not one."))
        desc, names = describe(snap, conventions, description)
        with_store_context(identifier; conventions=names) do
            assemble(group, snap, desc, names, vars, lazy, validate, samples)
        end
    end
end

function DiscreteGlobalGrids.dggread(store::StoreLike, var::Symbol; kwargs...)
    haskey(kwargs, :vars) && throw(ArgumentError(
        "`dggread(store, var)` reads the one variable named; `vars` is the " *
        "DimStack method's keyword."))
    return DiscreteGlobalGrids.dggread(store; vars=(var,), kwargs...)[var]
end

"""
    validation_samples(validate) -> Union{Int,Nothing}

Return the per-chunk validation sample count. `:lazy` returns `LAZY_SAMPLES`;
`:strict` and `:scan` return `nothing` to request every id. `:scan` additionally
disables the trusted-manifest path in `storedaxis`.
"""
function validation_samples(validate::Symbol)
    (validate === :strict || validate === :scan) && return nothing
    validate === :lazy && return LAZY_SAMPLES
    throw(ArgumentError(
        "`validate` is `:strict` (every id checked, and a chunk manifest this " *
        "package wrote taken on its word), `:scan` (the same checks, on every " *
        "store, sidecar or not) or `:lazy` (sampled), not $(repr(validate))"))
end

"""
    normalize_store_url(url) -> String

Normalize `gs://BUCKET/PATH` to the anonymous Google Cloud Storage endpoint
`https://storage.googleapis.com/BUCKET/PATH`. Return every other string
unchanged. Private buckets require a credentialed store supplied by the caller.
"""
function normalize_store_url(url::AbstractString)
    startswith(url, "gs://") || return String(url)
    return string(GCS_BASE, lstrip(SubString(url, 6), '/'))
end

"""
    opengroup(store) -> (ZGroup, identifier)

The group `store` names, and the URL or path to report it by.
"""
opengroup(g::Zarr.ZGroup) = (g, storeidentifier(g))

function opengroup(s::Zarr.AbstractStore)
    g = Zarr.zopen(s, "r")
    return requiregroup(g, storeurl(s))
end

function opengroup(url::AbstractString)
    u = normalize_store_url(url)
    startswith(u, "s3://") && !s3_supported() && throw(DGGSFormatError(
        check=:unsupported_store_scheme, store=u, declared="s3://",
        detail="Zarr.jl reads s3:// through its own AWSS3 extension: run " *
               "`using AWSS3` first, or open the store by its https:// URL."))
    # Zarr's own error says what is wrong with a path or a host; a store that
    # does not open is not a DGGS format failure and is not reported as one.
    return requiregroup(Zarr.zopen(u, "r"), u)
end

requiregroup(g::Zarr.ZGroup, identifier) = (g, String(identifier))

requiregroup(x, identifier) = throw(DGGSFormatError(check=:not_a_group,
    store=String(identifier), observed=string(nameof(typeof(x))),
    detail="a DGGS store is a GROUP of arrays sharing one cell dimension; " *
           "this path names a single array."))

s3_supported() = Base.get_extension(Zarr, :ZarrAWSS3Ext) !== nothing

# A separate detection pass preserves convention names for provenance and errors.
function describe(snap, conventions, ::Nothing)
    names = String[DiscreteGlobalGrids.conventionname(det.convention)
                   for det in DiscreteGlobalGrids.detections(snap; conventions)]
    return DiscreteGlobalGrids.describe_store(snap; conventions), names
end

# A supplied description asserts complete layout metadata and skips detection.
function describe(snap, conventions, description::StoreDescription)
    desc = DiscreteGlobalGrids.applydefaults(description)
    return DiscreteGlobalGrids.requirecomplete(desc; store=snap.identifier), String[]
end

function assemble(group, snap, desc, names, vars, lazy, validate, samples)
    dim = desc.spatial_dimension
    n = axislength(snap, dim)
    lookup = storedlookup(desc.encoding, group, snap, desc, n, validate, samples)

    available = desc.variables === nothing ?
                arrays_on(snap, dim; exclude=(desc.coordinate,)) : desc.variables
    selected = selectvars(available, vars)

    otherdims = Dict{String,Any}()
    layers = map(selected) do name
        entry = getarray(snap, name)
        z = group[name]
        dims = layerdims(entry, dim, lookup, group, snap, otherdims)
        DD.DimArray(lazy ? z : Array(z), dims;
            name=Symbol(name), metadata=deepcopy(entry.attrs))
    end

    metadata = Dict{String,Any}(
        "source" => snap.identifier,
        "conventions" => names,
        "encoding" => encodingname(desc.encoding),
        "attrs" => deepcopy(snap.attrs),
        "description" => desc)
    return DD.DimStack(NamedTuple{Tuple(Symbol.(selected))}(Tuple(layers)); metadata)
end

# --- the grid ---------------------------------------------------------------

function describedsystem(desc::StoreDescription)
    ref = get(DiscreteGlobalGrids.GRID_REFERENCE, desc.gridname, nothing)
    # A description that names the system and the scheme itself needs no table
    # entry: that is what the escape hatch is for. Anything less is looked up,
    # and an unregistered name raises there with the registry listed.
    ref === nothing && (desc.system === nothing || desc.idscheme === nothing) &&
        (ref = DiscreteGlobalGrids.gridreference(desc.gridname))
    system = desc.system === nothing ? ref.system : desc.system
    scheme = desc.idscheme === nothing ? ref.idscheme : desc.idscheme
    # Ids of another scheme are the same integers in a different order, so
    # reading them as the canonical one would misplace every cell silently.
    ref === nothing || scheme === ref.idscheme || throw(DGGSFormatError(
        check=:unsupported_indexing_scheme, declared=scheme, observed=ref.idscheme,
        detail="the id arithmetic of `$(desc.gridname)` is implemented for the " *
               "`$(ref.idscheme)` scheme; ids written in `$scheme` have to be " *
               "reindexed before they can be read as a cell axis."))
    return system
end

function describedgrid(desc::StoreDescription)
    system = describedsystem(desc)
    desc.level === nothing && throw(DGGSFormatError(check=:missing_level,
        declared=desc.gridname,
        detail="the store declares no refinement level, and a stored cell axis " *
               "is read at one level; pass `description` to assert it."))
    return DiscreteGlobalGrids.levelgrid(system, desc.level)
end

# --- the axis length --------------------------------------------------------

# Compare per-array Zarr shapes to establish one spatial-axis length.
function axislength(snap, dim)
    onaxis = [a for a in snap.arrays if dim in a.dims]
    isempty(onaxis) && throw(DGGSFormatError(check=:unknown_axis_length,
        declared=dim,
        detail="no array in this group has a dimension named `$dim`, so the " *
               "cell axis has no length to read."))
    lengths = unique(Int[a.shape[findfirst(==(dim), a.dims)] for a in onaxis])
    length(lengths) == 1 || throw(DGGSFormatError(
        check=:dimension_length_disagreement, declared=dim,
        observed=Dict(a.name => a.shape[findfirst(==(dim), a.dims)] for a in onaxis),
        detail="the arrays on dimension `$dim` declare $(join(lengths, " and ")) " *
               "cells; one dimension has one length."))
    return only(lengths)
end

# --- the axis ---------------------------------------------------------------

"""
    storedlookup(encoding, group, snapshot, desc, n, validate, samples) -> Lookup

The cube axis of the store: a [`ChunkedCellLookup`](@ref) over the
single-level axis `storedaxis` builds, or — for a `compacted` store — a
`MultiOrderLookup` over the validated `MultiOrderVector` its two columns name.
"""
storedlookup(enc::CellEncoding, group, snap, desc, n, validate, samples) =
    ChunkedCellLookup(storedaxis(enc, describedgrid(desc), group, snap, desc, n,
        validate, samples))

function storedlookup(enc::CompactedEncoding, group, snap, desc, n, validate, samples)
    desc.level === nothing || throw(DGGSFormatError(check=:level_with_compacted,
        declared=desc.level,
        detail="compacted cells declare one level per coordinate and require " *
               "`refinement_level: null`; found refinement level " *
               "$(desc.level)."))
    sys = describedsystem(desc)
    z = coordinatearray(group, snap, desc, 1)
    lv = levelscolumn(group, snap, n)
    # Validate every compacted pair regardless of the single-level scan policy.
    mov = cellaxis(enc, sys, lv, Array(z); declared_length=n)
    return MultiOrderLookup(mov)
end

function levelscolumn(group, snap, n)
    entry = getarray(snap, COMPACTED_LEVELS_ARRAY)
    (entry === nothing || !haskey(group, COMPACTED_LEVELS_ARRAY)) &&
        throw(DGGSFormatError(check=:missing_levels_array,
            declared=COMPACTED_LEVELS_ARRAY,
            observed=sort!(String[a.name for a in snap.arrays]),
            detail="a compacted store keeps each cell's level in " *
                   "`$COMPACTED_LEVELS_ARRAY`, and this one has no such array."))
    (entry.shape == (n,) && entry.eltype <: Integer) || throw(DGGSFormatError(
        check=:invalid_levels_array, declared=(n,), observed=entry.shape,
        detail="`$COMPACTED_LEVELS_ARRAY` holds one integer level per cell " *
               "of the axis; expected $n of them, found shape $(entry.shape) " *
               "of $(entry.eltype)."))
    return Int.(Array(group[COMPACTED_LEVELS_ARRAY]))
end

"""
    storedaxis(encoding, grid, group, snapshot, desc, n, validate, samples)
        -> ChunkedCellVector

Build a single-level stored axis according to its encoding:

  - Dense encoding scans the id array in chunks.
  - Ranges encoding reads the small range array once.
  - Implicit encoding uses only the axis length.
  - Dense encoding with a trusted persisted manifest defers id reads until
    selection.

`validate = :scan` disables the trusted-manifest path. Arithmetic encodings
require no id scan under any validation mode.

This function is the read-side extension point for single-level layouts. A
downstream encoding registers an instance, implements `cellaxis`, and adds a
`storedaxis` method that supplies its source data. Encodings with another axis
type add `storedlookup`, as compacted encoding does. Missing implementations
raise a named unsupported-encoding error.
"""
function storedaxis(enc::DenseEncoding, grid, group, snap, desc, n, validate, samples)
    z = coordinatearray(group, snap, desc, 1)
    cl = chunklength(z)
    # A manifest this package persisted is the chunk grid, and skipping the scan
    # for it is what makes a store of tens of millions of cells open at all —
    # which is also exactly what `:scan` is for declining.
    if validate !== :scan
        found = persistedmanifest(group, snap, desc, grid, Int(n), Int(cl))
        found === nothing || return cellaxis(enc, grid, z, found.manifest;
            store=snap.identifier, sidecar=found.sidecar)
    end
    # One block per stored chunk, so the scan reads each chunk exactly once and
    # the lookup's later per-chunk reads land on the same boundaries.
    return cellaxis(enc, grid, z; chunklength=cl, declared_length=n, samples)
end

function storedaxis(enc::RangesEncoding, grid, group, snap, desc, n, validate, samples)
    z = coordinatearray(group, snap, desc, 2)
    # Kilobytes even for a store of tens of millions of cells, and the axis is
    # arithmetic over it from there: `(n, 2)` as Zarr declares it, which is
    # `(2, n)` in Julia's own order.
    return cellaxis(enc, grid, permutedims(Array(z)); declared_length=n)
end

storedaxis(enc::ImplicitEncoding, grid, group, snap, desc, n, validate, samples) =
    cellaxis(enc, grid, n)

function storedaxis(enc::CellEncoding, grid, group, snap, desc, n, validate, samples)
    throw(DGGSFormatError(check=:unsupported_encoding, observed=encodingname(enc),
        detail="the Zarr reader implements the dense, ranges and implicit " *
               "layouts; `$(encodingname(enc))` names an encoding it has no " *
               "read path for."))
end

function coordinatearray(group, snap, desc, n::Int)
    name = desc.coordinate
    name === nothing && throw(DGGSFormatError(check=:missing_coordinate,
        declared=encodingname(desc.encoding),
        detail="the $(encodingname(desc.encoding)) layout stores its cell ids " *
               "in an array, and the description names none."))
    haskey(group, name) || throw(DGGSFormatError(check=:missing_coordinate,
        declared=name, observed=sort!(collect(keys(group.arrays))),
        detail="the store's cell coordinate `$name` is not one of its arrays."))
    z = group[name]
    ndims(z) == n || throw(DGGSFormatError(check=:invalid_coordinate_shape,
        declared=n, observed=reverse(size(z)),
        detail="the $(encodingname(desc.encoding)) layout's coordinate `$name` " *
               "is $n-dimensional; this one has $(ndims(z)) dimensions."))
    return z
end

chunklength(z::Zarr.ZArray) = first(z.metadata.chunks)

# --- the persisted manifest -------------------------------------------------

"""
    persistedmanifest(group, snapshot, desc, grid, n, chunklength)
        -> Union{NamedTuple, Nothing}

Return a trusted persisted [`ChunkManifest`](@ref) and its sidecar name. Return
`nothing` when any trust check fails; the caller then scans the coordinate.

The marker must match the reader's writer and format, strict validation level,
grid, refinement level, spatial dimension, chunk length and axis length. The
sidecar rows must also form ascending, disjoint chunk intervals.

Trust covers the writer's earlier proof that coordinate ids are sorted, unique
and valid at the declared level. Later chunk reads verify each chunk's first id,
last id and length. They do not rescan interior ids; `validate = :scan` requests
that stronger check.
"""
function persistedmanifest(group, snap, desc, grid, n::Int, cl::Int)
    (n >= 1 && cl >= 1) || return nothing
    entry = manifestentry(snap, desc, n, cl)
    entry === nothing && return nothing
    nc = cld(n, cl)
    # Zarr's `(n_chunks, 2)` as the snapshot reports it, which is the shape the
    # array really has: `arrayentry` reads it off the same `ZArray` the rows are
    # then read from.
    (entry.shape == (nc, 2) && entry.eltype <: Integer) || return nothing
    haskey(group, entry.name) || return nothing
    rows = group[entry.name][:, :]      # `(2, n_chunks)` in Julia's own order
    I = idtype(grid)
    firstids = storedids(I, view(rows, 1, :))
    lastids = storedids(I, view(rows, 2, :))
    (firstids === nothing || lastids === nothing) && return nothing
    disjointrows(firstids, lastids) || return nothing
    # Zarr chunks are uniform by format, so the shape of the grid follows from
    # the two numbers the marker was believed on.
    lengths = Int[c == nc ? n - (nc - 1) * cl : cl for c in 1:nc]
    offsets = Int[(c - 1) * cl for c in 1:nc]
    return (manifest=ChunkManifest(firstids, lastids, lengths, offsets, cl),
        sidecar=entry.name)
end

# Search every marked array so name ordering cannot let a foreign sidecar shadow ours.
function manifestentry(snap, desc, n::Int, cl::Int)
    for a in snap.arrays
        marker = get(a.attrs, MANIFEST_MARKER, nothing)
        marker === nothing && continue
        trustedmarker(marker, desc, n, cl) && return a
    end
    return nothing
end

function trustedmarker(marker, desc, n::Int, cl::Int)
    marker isa AbstractDict || return false
    for (key, want) in (("writer", MANIFEST_WRITER), ("format", MANIFEST_FORMAT),
        ("validated", MANIFEST_VALIDATED),
        ("spatial_dimension", desc.spatial_dimension),
        ("grid", desc.gridname), ("level", desc.level),
        ("chunk_length", cl), ("length", n))
        get(marker, key, nothing) == want || return false
    end
    return true
end

# Reject a manifest whose ids do not fit the grid's exact integer width.
function storedids(::Type{I}, values) where {I<:Integer}
    out = Vector{I}(undef, length(values))
    for (i, x) in pairs(values)
        typemin(I) <= x <= typemax(I) || return nothing
        out[i] = x % I
    end
    return out
end

# Binary search requires ascending, disjoint chunk intervals.
function disjointrows(firstids, lastids)
    for c in eachindex(firstids)
        firstids[c] <= lastids[c] || return false
        c == 1 || lastids[c-1] < firstids[c] || return false
    end
    return true
end

# --- the layers -------------------------------------------------------------

function selectvars(available, vars)
    selected = varnames(available, vars)
    isempty(selected) && throw(DGGSFormatError(check=:no_data_variables,
        observed=available,
        detail="no data variable was selected; the store's cell axis carries " *
               (isempty(available) ? "no arrays." : join(available, ", ") * ".")))
    return selected
end

varnames(available, ::DD.All) = available
varnames(available, ::Colon) = available
varnames(available, ::Nothing) = available
varnames(available, v::Union{Symbol,AbstractString}) = varnames(available, (v,))

function varnames(available, vars)
    out = String[]
    for v in vars
        name = String(v)
        name in available || throw(DGGSFormatError(check=:unknown_variable,
            declared=name, observed=available,
            detail="the store has no data variable `$name`; it has " *
                   join(available, ", ") * "."))
        push!(out, name)
    end
    return out
end

# Reversing Zarr dimension order keeps cells contiguous in Julia and enforces one
# shared length for each named dimension.
function layerdims(entry, spatial, lookup, group, snap, cache)
    return map(reverse(Tuple(entry.dims)), reverse(entry.shape)) do name, len
        name == spatial && return Cells(lookup)
        d = get!(cache, name) do
            otherdim(name, len, group, snap)
        end
        length(d) == len || throw(DGGSFormatError(
            check=:dimension_length_disagreement, declared=name,
            observed=Dict(entry.name => len),
            detail="array `$(entry.name)` gives dimension `$name` $len steps, " *
                   "which the arrays read before it made $(length(d)) long."))
        return d
    end
end

# A matching one-dimensional coordinate supplies values; other dimensions use indices.
function otherdim(name, len, group, snap)
    entry = getarray(snap, name)
    D = DD.Dim{Symbol(name)}
    if entry !== nothing && entry.dims == [name] && entry.shape == (len,)
        return D(Array(group[name]))
    end
    return D(Lookups.NoLookup(Base.OneTo(len)))
end
