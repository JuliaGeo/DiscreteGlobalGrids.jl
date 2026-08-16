# The user-facing verbs. Bare `regrid` builds an ephemeral plan and drops it;
# `plan_regrid` hands one back. Applying a plan takes no keyword arguments
# because the plan already answers every question they would.

# Whether a source is chunked storage rather than a plain in-memory array.
# `DiskArrays.haschunks` is total, so this is total.
isdiskbacked(data) = DiskArrays.haschunks(data) isa DiskArrays.Chunked

# The metadata keys a nodata sentinel is written under, in the order they are
# consulted. `_FillValue` and `missing_value` are the CF spellings a NetCDF or
# GRIB reader carries through; `missingval` is what a reader that has already
# normalized them says.
const MISSINGVAL_KEYS = ("missingval", "_FillValue", "missing_value")

"""
    sourcemissingval(data) -> value or nothing

The nodata sentinel `data` declares of itself, or `nothing` when it declares
none. This is the default of [`regrid`](@ref)'s `missingval`, so a source that
carries its sentinel is handled without the caller repeating it; an explicit
`missingval =` always wins.

A `DimensionalData` array is asked for its metadata and answers with the first
of `$(MISSINGVAL_KEYS)` it carries, under a `String` or a `Symbol` key.
Anything else answers `nothing`.

The extension point for a reader whose sentinel lives somewhere else — a struct
field rather than metadata, as `Rasters.missingval` does — is a method of this
function on that reader's array type.
"""
sourcemissingval(::Any) = nothing

function sourcemissingval(data::DD.AbstractDimArray)
    md = DD.metadata(data)
    md isa DD.NoMetadata && return nothing
    for k in MISSINGVAL_KEYS
        haskey(md, k) && return md[k]
        s = Symbol(k)
        haskey(md, s) && return md[s]
    end
    return nothing
end

"""
    regrid(data; to, from = nothing, method = Conservative(),
           missingpolicy = Weighted(0.5), missingval = sourcemissingval(data),
           lazy = isdiskbacked(data), chunks = nothing, budget = 2^30,
           storage = nothing)
    regrid(data, plan::AbstractRegriddingPlan)

Regrid `data` onto the cells named by `to`.

Spatial dimensions are replaced by the destination's; every other dimension
passes through unchanged, and the output dimensions are
`(destination spatial dim, other dims...)`. One spatial plan is built and
reused across the non-spatial slices, so an N-D call costs one weight
construction, not one per field.

`data` flattens over its spatial dimensions ([`spatialdims`](@ref)) in memory
order, which must be the order the source space numbers its cells in. The
spatial dimensions must come first; `permutedims` a source whose do not. The
result's cell axis is a plain `1:ncells` axis, which a package supplying spaces
replaces with its own.

Values are floating point, and destinations blanked by [`Weighted`](@ref) are
`missing` when the source element type carries `missing` and NaN otherwise.

# Keyword arguments

  - `to`: the destination, a [`RegridSpace`](@ref). Packages that supply spaces
    accept their own target spellings and resolve them into one.
  - `from`: the source space. `nothing` builds a [`RasterGrid`](@ref) from
    `data`.
  - `method`: an [`AbstractRegriddingMethod`](@ref). [`Conservative`](@ref) by
    default, the only exactly conservative choice.
  - `missingpolicy`: [`Weighted`](@ref) — coverage-normalized mean, blanking
    destinations below the threshold — or [`Extensive`](@ref) — raw sums.
  - `missingval`: the source's nodata sentinel, invalid on top of `missing` and
    NaN, which are always invalid. `nothing` means there is none; the default
    is whatever the source declares of itself ([`sourcemissingval`](@ref)), so a
    reader that carries its sentinel needs no keyword here, and passing one
    overrides what it says. Apply-time only: weights stay geometry-blind, and a
    sentinel field gives exactly what the same field with its sentinels replaced
    by NaN gives.
  - `lazy`: return a [`LazyRegridArray`](@ref) instead of a materialized array.
    Defaults to whether `data` is chunked storage. Constructing one reads no
    source data at all.
  - `chunks`: the destination's chunking, as a `DiskArrays.GridChunks` or a
    tuple of chunk sizes. `nothing` derives it from the destination space, or
    from `budget` when the destination space has no chunking of its own. It is
    also the destination tiling the lazy path computes one tile at a time — a
    lazy-path knob, ignored by the eager path, which produces one whole array.
  - `budget`: the bytes of transient memory a lazy read aims to stay within,
    shared between the weights it holds and the source data it loads. Purely a
    performance knob — connected source chunks are held together when they fit
    and streamed one at a time when they do not, to the same answer either way —
    and likewise ignored by the eager path.
  - `storage`: where a lazy plan keeps built weights — [`PerChunk`](@ref) in
    memory, [`Spilled`](@ref) on disk. `nothing` takes a `PerChunk` bounded by
    the weight share of `budget`.

The one-argument-plus-plan form takes no keyword arguments: a plan already
carries the method, both spaces, the missing policy, and the storage and budget
policy. Bare `regrid` is exactly `plan_regrid`, apply, drop — a one-shot call
retains nothing.
"""
function regrid end

function regrid(data; to, from = nothing,
    method::AbstractRegriddingMethod = Conservative(),
    missingpolicy::AbstractMissingPolicy = Weighted(0.5),
    missingval = sourcemissingval(data),
    lazy::Bool = isdiskbacked(data), chunks = nothing, budget::Integer = 2^30,
    storage::Union{Nothing,AbstractBlockStorage} = nothing)
    plan = plan_regrid(data; to, from, method, missingpolicy, missingval, lazy,
        chunks, budget, storage)
    return regrid(data, plan)
end

function regrid(data, plan::DirectPlan)
    sd, othersizes, src = _flatten(data, plan)
    ndst = size(plan.block, 1)
    out = Array{outputeltype(eltype(data))}(undef, ndst, othersizes...)
    applyplan!(reshape(out, ndst, prod(othersizes)), plan, src)
    return wrapoutput(out, data, sd)
end

regrid(data, plan::AbstractRegriddingPlan) =
    error("$(typeof(plan).name.name) defines no `regrid` application")

"""
    regrid!(dest, data; to, from = nothing, method = Conservative(),
            missingpolicy = Weighted(0.5), missingval = sourcemissingval(data),
            lazy = isdiskbacked(data), chunks = nothing, budget = 2^30,
            storage = nothing)
    regrid!(dest, data, plan::AbstractRegriddingPlan)

Regrid `data` into the preallocated `dest` and return `dest`.

`dest`'s first dimension is the destination's cells and the rest are `data`'s
non-spatial dimensions, in order.

Keyword arguments are [`regrid`](@ref)'s. The plan form takes none.
"""
function regrid! end

function regrid!(dest, data; to, from = nothing,
    method::AbstractRegriddingMethod = Conservative(),
    missingpolicy::AbstractMissingPolicy = Weighted(0.5),
    missingval = sourcemissingval(data),
    lazy::Bool = isdiskbacked(data), chunks = nothing, budget::Integer = 2^30,
    storage::Union{Nothing,AbstractBlockStorage} = nothing)
    plan = plan_regrid(data; to, from, method, missingpolicy, missingval, lazy,
        chunks, budget, storage)
    return regrid!(dest, data, plan)
end

function regrid!(dest, data, plan::DirectPlan)
    _, othersizes, src = _flatten(data, plan)
    ndst = size(plan.block, 1)
    size(dest) == (ndst, othersizes...) || throw(DimensionMismatch(
        "destination of size $(size(dest)) cannot hold a regrid of size " *
        "$((ndst, othersizes...))"))
    raw = dest isa DD.AbstractDimArray ? parent(dest) : dest
    applyplan!(reshape(raw, ndst, prod(othersizes)), plan, src)
    return dest
end

regrid!(dest, data, plan::AbstractRegriddingPlan) =
    error("$(typeof(plan).name.name) defines no `regrid!` application")

"""
    plan_regrid(data; to, from = nothing, method = Conservative(),
                missingpolicy = Weighted(0.5), missingval = sourcemissingval(data),
                lazy = isdiskbacked(data), chunks = nothing, budget = 2^30,
                storage = nothing) -> AbstractRegriddingPlan

Build the regridding operator from `data`'s space onto `to` without applying it,
for reuse across fields, slices, and repeated reads.

Keyword arguments are [`regrid`](@ref)'s, and the plan retains all of them:
`regrid(data, plan)` and `regrid!(dest, data, plan)` accept none.

Building a plan reads **no source data**. Weights are geometry-only, and a
chunked plan's blocks are built on first touch rather than up front, so planning
a regrid of a terabyte source is as cheap as planning a megabyte one. Reusing
the plan is what makes the second and later applications free of weight
construction — the reason to reach for `plan_regrid` over bare
[`regrid`](@ref).

An in-memory source plans to a [`DirectPlan`](@ref): one whole-domain block over
both spaces entire.

A lazy plan's default `storage` is a [`PerChunk`](@ref) bounded by the weight
share of `budget` ([`weightbudget`](@ref)) — so a plan's resident weights are
bounded by what the caller asked for rather than by how much of the destination
has been touched. Pass `storage = PerChunk()` for the unbounded cache, or
`storage = Spilled(dir)` to keep the weights on disk instead.
"""
function plan_regrid(data; to, from = nothing,
    method::AbstractRegriddingMethod = Conservative(),
    missingpolicy::AbstractMissingPolicy = Weighted(0.5),
    missingval = sourcemissingval(data),
    lazy::Bool = isdiskbacked(data), chunks = nothing, budget::Integer = 2^30,
    storage::Union{Nothing,AbstractBlockStorage} = nothing)
    dst_space = _asspace(to, "to")
    src_space = from === nothing ? _sourcespace(data) : _asspace(from, "from")
    _checkchunks(chunks)
    budget > 0 || throw(ArgumentError("budget must be positive, got $budget"))
    manifold(dst_space) == manifold(src_space) || throw(ArgumentError(
        "the two sides of a regrid must live on one manifold, but the source " *
        "is on $(manifold(src_space)) and the destination on $(manifold(dst_space))"))
    lazy && return ChunkedPlan(method, missingpolicy, dst_space, src_space;
        storage, budget, chunks, missingval)
    storage === nothing || throw(ArgumentError(
        "`storage` is a lazy-path knob: an eager plan holds one whole-domain " *
        "block and has nothing to store it in. Pass `lazy = true` to use it."))
    return DirectPlan(method, missingpolicy, dst_space, src_space,
        wholeblock(method, dst_space, src_space), missingval)
end

"""
    wholeblock(method, dst_space, src_space) -> WeightBlock

The [`WeightBlock`](@ref) of one unrestricted pair: every cell of `dst_space`
against every cell of `src_space`, in one [`build_weights!`](@ref) call. Its
chunk-local indices are therefore cell positions.
"""
function wholeblock(method::AbstractRegriddingMethod, dst_space::RegridSpace,
    src_space::RegridSpace)
    ndst = Int(ncells(dst_space))
    nsrc = Int(ncells(src_space))
    coo = WeightCOO(ndst)
    build_weights!(coo, method, dst_space, 1:ndst, src_space, 1:nsrc)
    return WeightBlock(coo, ndst, nsrc)
end

# The source space when `from` was not given. Only a dimensional raster carries
# enough geometry to derive one; a plain array is numbers with no coordinates on
# them, and guessing an extent for it would be worse than saying so.
_sourcespace(data::DD.AbstractDimArray) = RasterGrid(data)
_sourcespace(data) = throw(ArgumentError(
    "`from` was not given, and a source space cannot be derived from a " *
    "$(typeof(data)): it carries no coordinates. Pass `from = ` a RegridSpace, " *
    "or a DimensionalData array whose dimensions describe the raster."))

function _asspace(space, name)
    space isa RegridSpace || throw(ArgumentError(
        "`$name` must be a RegridSpace, got $(typeof(space)). A package that " *
        "supplies spaces resolves its own target spellings into one."))
    return space
end

_checkchunks(::Nothing) = nothing
_checkchunks(chunks::Tuple{Vararg{Integer}}) =
    all(>(0), chunks) ? nothing :
    throw(ArgumentError("chunk sizes must be positive, got $chunks"))
_checkchunks(::DiskArrays.GridChunks) = nothing
_checkchunks(chunks) = throw(ArgumentError(
    "`chunks` must be a tuple of chunk sizes, a DiskArrays.GridChunks, or " *
    "nothing, got $(typeof(chunks))"))

# The source as an `ncells × nslices` matrix, plus the shape the pass-through
# dimensions keep.
function _flatten(data, plan::AbstractRegriddingPlan)
    nsrc = Int(ncells(plan.src_space))
    sd = resolvespatialdims(data, nsrc)
    othersizes = _otherdimsizes(data, sd)
    return sd, othersizes, flatsource(data, nsrc, prod(othersizes))
end
