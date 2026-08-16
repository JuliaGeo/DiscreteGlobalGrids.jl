# `dggread`: open a store, describe it, build its cell axis, and hand back a
# DimensionalData cube over the arrays that share it.
#
# The pipeline is
#
#   open -> snapshot -> describe_store -> cellaxis -> DimStack
#
# and every step but the first two lives in `src/io/`: this file supplies the
# store access and nothing format-semantic. Errors from those layers cross the
# API boundary through `with_store_context`, which fills in the store
# identifier and the conventions that fired without changing the error.

import DimensionalData as DD
using DimensionalData: Lookups

using DiscreteGlobalGrids: StoreDescription, DGGSFormatError, Cells,
    ChunkedCellLookup, with_store_context, arrays_on, getarray
using DiscreteGlobalGrids.Encodings: CellEncoding, DenseEncoding, RangesEncoding,
    ImplicitEncoding, cellaxis, encodingname

"""
    LAZY_SAMPLES

How many ids per chunk `validate = :lazy` checks for existence. Both ends of
the chunk — the ids the manifest publishes — plus an even interior spread. A
phantom between two samples survives, which is what the opt-out buys.
"""
const LAZY_SAMPLES = 8

"What `dggread` will open: an already-open group, a store, or a path or URL."
const StoreLike = Union{Zarr.ZGroup,Zarr.AbstractStore,AbstractString}

# ===========================================================================
# The API surface
# ===========================================================================

"""
    dggread(store; vars = All(), lazy = true, validate = :strict,
            conventions = CONVENTION_REGISTRY, description = nothing) -> DimStack
    dggread(store, var::Symbol; kwargs...) -> DimArray

Read a DGGS store into plain DimensionalData: one `Cells` dimension shared by
every layer, carrying a [`ChunkedCellLookup`](@ref) that resolves a cell to a
position without scanning the axis. Any other dimension of the store — time,
bands — becomes an ordinary `Dim`, with the values of the like-named coordinate
array where the store has one.

`store` is a `Zarr.ZGroup`, a `Zarr.AbstractStore`, a local path, or a URL. A
`gs://BUCKET/PATH` URL is read as `https://storage.googleapis.com/BUCKET/PATH`,
which works for public-read buckets; an `s3://` URL needs Zarr's own AWSS3
extension (`using AWSS3`) and says so otherwise.

  - `vars`: `All()` for every data variable, or the `Symbol`s to read. An
    unknown name raises, listing what the store holds.
  - `lazy`: leave the data as the store's own chunked arrays (the default), or
    materialize them.
  - `validate`: `:strict` checks that every stored id names a cell of the
    declared level; `:lazy` checks [`LAZY_SAMPLES`](@ref) per chunk instead.
    Sortedness, uniqueness and the length checks are not optional and run
    either way.
  - `conventions`: the conventions to try, in order.
  - `description`: a [`StoreDescription`](@ref) that bypasses detection. The
    caller then asserts grid, level, encoding and array names, and only the
    mechanical checks — the id scan, the closed-form counts — still run. This
    is how an attribute-less store is read.

The stack's `metadata` carries the provenance a value-identical rewrite needs:

| key | |
|---|---|
| `"source"` | the store URL or path |
| `"conventions"` | the conventions that fired, in order |
| `"encoding"` | the cell-axis layout, as the store spells it |
| `"attrs"` | the group attributes verbatim |
| `"description"` | the [`StoreDescription`](@ref) everything was read through |

Each layer keeps its own array attributes as its metadata.
"""
function DiscreteGlobalGrids.dggread(store::StoreLike; vars=DD.All(), lazy::Bool=true,
    validate::Symbol=:strict, conventions=DiscreteGlobalGrids.CONVENTION_REGISTRY,
    description::Union{StoreDescription,Nothing}=nothing)
    samples = validation_samples(validate)
    group, identifier = opengroup(store)
    return with_store_context(identifier) do
        snap = snapshot(group; identifier)
        desc, names = describe(snap, conventions, description)
        with_store_context(identifier; conventions=names) do
            assemble(group, snap, desc, names, vars, lazy, samples)
        end
    end
end

function DiscreteGlobalGrids.dggread(store::StoreLike, var::Symbol; kwargs...)
    haskey(kwargs, :vars) && throw(ArgumentError(
        "`dggread(store, var)` reads the one variable named; `vars` is the " *
        "DimStack method's keyword."))
    return DiscreteGlobalGrids.dggread(store; vars=(var,), kwargs...)[var]
end

function validation_samples(validate::Symbol)
    validate === :strict && return nothing
    validate === :lazy && return LAZY_SAMPLES
    throw(ArgumentError(
        "`validate` is `:strict` (every id checked) or `:lazy` (sampled), not $(repr(validate))"))
end

# ===========================================================================
# Opening
# ===========================================================================

"""
    normalize_store_url(url) -> String

`gs://BUCKET/PATH` as `https://storage.googleapis.com/BUCKET/PATH`, and every
other string unchanged.

Both forms reach the same anonymous Google Cloud Storage endpoint, so this is
a spelling, not a capability: a bucket that is not public-read is unreadable
either way, and needs credentials that this package does not carry.
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
    store=String(identifier), observed=typeof(x),
    detail="a DGGS store is a GROUP of arrays sharing one cell dimension; " *
           "this path names a single array."))

s3_supported() = Base.get_extension(Zarr, :ZarrAWSS3Ext) !== nothing

# ===========================================================================
# Description
# ===========================================================================

# Detection runs twice: once for the names the errors and the provenance
# report, once inside `describe_store`, which reconciles but does not say who
# spoke. Both passes only read attribute dictionaries.
function describe(snap, conventions, ::Nothing)
    names = String[DiscreteGlobalGrids.conventionname(det.convention)
                   for det in DiscreteGlobalGrids.detections(snap; conventions)]
    return DiscreteGlobalGrids.describe_store(snap; conventions), names
end

# A supplied description is an assertion, so nothing is detected and nothing is
# merged; it still has to name a layout and a data dimension.
function describe(snap, conventions, description::StoreDescription)
    desc = DiscreteGlobalGrids.applydefaults(description)
    return DiscreteGlobalGrids.requirecomplete(desc; store=snap.identifier), String[]
end

# ===========================================================================
# Assembly
# ===========================================================================

function assemble(group, snap, desc, names, vars, lazy, samples)
    grid = describedgrid(desc)
    dim = desc.spatial_dimension
    n = axislength(snap, dim)
    axis = storedaxis(desc.encoding, grid, group, snap, desc, n, samples)
    lookup = ChunkedCellLookup(axis)

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

function describedgrid(desc::StoreDescription)
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
    desc.level === nothing && throw(DGGSFormatError(check=:missing_level,
        declared=desc.gridname,
        detail="the store declares no refinement level, and a stored cell axis " *
               "is read at one level; pass `description` to assert it."))
    return DiscreteGlobalGrids.levelgrid(system, desc.level)
end

# --- the axis length --------------------------------------------------------

# The spatial dimension's extent, which every array carrying it must agree on:
# Zarr stores a shape per array and nothing anywhere says they match.
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
    storedaxis(encoding, grid, group, snapshot, desc, n, samples) -> ChunkedCellVector

The store's cell axis, read the way its encoding keeps it: a chunked pass over
the id array, the small range array read whole, or no read at all.
"""
function storedaxis(enc::DenseEncoding, grid, group, snap, desc, n, samples)
    z = coordinatearray(group, snap, desc, 1)
    # One block per stored chunk, so the scan reads each chunk exactly once and
    # the lookup's later per-chunk reads land on the same boundaries.
    return cellaxis(enc, grid, z; chunklength=chunklength(z),
        declared_length=n, samples)
end

function storedaxis(enc::RangesEncoding, grid, group, snap, desc, n, samples)
    z = coordinatearray(group, snap, desc, 2)
    # Kilobytes even for a store of tens of millions of cells, and the axis is
    # arithmetic over it from there: `(n, 2)` as Zarr declares it, which is
    # `(2, n)` in Julia's own order.
    return cellaxis(enc, grid, permutedims(Array(z)); declared_length=n)
end

storedaxis(enc::ImplicitEncoding, grid, group, snap, desc, n, samples) =
    cellaxis(enc, grid, n)

function storedaxis(enc::CellEncoding, grid, group, snap, desc, n, samples)
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

# Julia's dimension order is the reverse of the store's, so the cell dimension
# of a `(time, cell)` array is the FIRST here and its data stays contiguous.
#
# One dimension has one length, whatever it is named: `axislength` says so for
# the cell axis, and this is the same check for every other one, since two
# layers disagreeing about `time` is the same lie about the same store.
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

# A dimension with a like-named one-dimensional array is that array's values —
# the CF coordinate-variable convention, and small enough to read. Anything
# else is a bare axis: DimensionalData indexes it by position.
function otherdim(name, len, group, snap)
    entry = getarray(snap, name)
    D = DD.Dim{Symbol(name)}
    if entry !== nothing && entry.dims == [name] && entry.shape == (len,)
        return D(Array(group[name]))
    end
    return D(Lookups.NoLookup(Base.OneTo(len)))
end
