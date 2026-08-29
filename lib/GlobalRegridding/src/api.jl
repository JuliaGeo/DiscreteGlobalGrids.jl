# Public regridding API.

declareschunks(data) = DiskArrays.haschunks(data) isa DiskArrays.Chunked

# Nodata metadata keys, in precedence order.
const MISSINGVAL_KEYS = ("missingval", "_FillValue", "missing_value")

"""
    sourcemissingval(data) -> value or nothing

Return the nodata sentinel declared by `data`, or `nothing`. Dimensional arrays
use the first matching `$(MISSINGVAL_KEYS)` metadata key. Readers may extend
this function when their sentinel is stored elsewhere.
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
           lazy = declareschunks(data), chunks = nothing, budget = nothing,
           storage = nothing, sampling = nothing)
    regrid(data, plan::AbstractRegriddingPlan)

Regrid `data` onto `to`. Spatial dimensions lead in source-cell order;
non-spatial dimensions retain their order. One plan serves every non-spatial
slice.

Dimensional results use the destination axes followed by the unchanged
non-spatial dimensions. A destination without axes uses one flat `Cell` axis.
Lazy results preserve the same labels and shape over a disk-backed array.

Results are floating point. [`Weighted`](@ref) writes `missing` when supported
by the source element type, and `NaN` otherwise.

# Keyword arguments

  - `to`: destination [`RegridSpace`](@ref) or a package-specific target.
  - `from`: source space. `nothing` infers a self-describing cell axis or a
    [`RasterGrid`](@ref).
  - `method`: weight-building method; defaults to [`Conservative`](@ref).
  - `missingpolicy`: [`Weighted`](@ref) means or [`Extensive`](@ref) sums.
  - `missingval`: additional nodata sentinel. `missing` and `NaN` remain invalid.
  - `lazy`: compute on demand with [`LazyRegridArray`](@ref); defaults to chunked
    sources.
  - `chunks`: lazy destination tiling. `nothing` derives it automatically.
  - `budget`: target bytes for lazy reads and weights, default `2^31`.
  - `storage`: lazy weight storage, [`PerChunk`](@ref) or [`Spilled`](@ref).
  - `sampling`: destination lookup sampling. `nothing` follows the method —
    area-based methods give `Intervals`, point samples give `Points`
    ([`outputsampling`](@ref)).

Keyword applicability:

  - Lazy regrids accept `chunks`, `budget`, and `storage`.
  - Eager regrids accept `sampling`.
  - The plan form accepts no keywords because its plan contains all settings.

The keyword form delegates to [`plan_regrid`](@ref). Dependency-relation
keywords belong to reusable plans and are accepted directly by `plan_regrid`.
"""
function regrid end

function regrid(data; kwargs...)
    _rejectplankeywords(kwargs, "regrid")
    return regrid(data, plan_regrid(data; kwargs...))
end

function regrid(data, plan::DirectPlan)
    sd, othersizes, src = _flatten(data, plan)
    ndst = size(plan.block, 1)
    out = Array{outputeltype(eltype(data))}(undef, ndst, othersizes...)
    applyplan!(reshape(out, ndst, prod(othersizes)), plan, src)
    return wrapoutput(out, data, sd, destinationdims(plan))
end

regrid(data, plan::AbstractRegriddingPlan) =
    error("$(typeof(plan).name.name) defines no `regrid` application")

"""
    destinationdims(plan::DirectPlan) -> Tuple or nothing
    destinationdims(plan::ChunkedPlan) -> Tuple or nothing

Return dimensions for results produced by `plan`. A direct plan uses its
sampling override when present, then [`outputsampling`](@ref). A chunked plan
uses `outputsampling`.
"""
destinationdims(plan::DirectPlan) = destinationdims(plan.dst_space,
    something(plan.sampling, outputsampling(plan.method)))

destinationdims(plan::ChunkedPlan) =
    destinationdims(plan.dst_space, outputsampling(plan.method))

"""
    regrid!(dest, data; to, from = nothing, method = Conservative(),
            missingpolicy = Weighted(0.5), missingval = sourcemissingval(data),
            lazy = declareschunks(data), chunks = nothing, budget = nothing,
            storage = nothing, sampling = nothing)
    regrid!(dest, data, plan::AbstractRegriddingPlan)

Regrid `data` into the preallocated `dest` and return `dest`.

`dest` starts with the destination axes or one flat cell dimension, followed by
the non-spatial dimensions of `data`. Keywords match [`regrid`](@ref) and
delegate to [`plan_regrid`](@ref). The plan form accepts no keywords.
"""
function regrid! end

function regrid!(dest, data; kwargs...)
    _rejectplankeywords(kwargs, "regrid!")
    return regrid!(dest, data, plan_regrid(data; kwargs...))
end

function regrid!(dest, data, plan::DirectPlan)
    _, othersizes, src = _flatten(data, plan)
    ndst = size(plan.block, 1)
    dstdims = destinationdims(plan)
    shaped = dstdims === nothing ? (ndst, othersizes...) :
             (map(length, dstdims)..., othersizes...)
    size(dest) == shaped || size(dest) == (ndst, othersizes...) ||
        throw(DimensionMismatch(
            "destination of size $(size(dest)) cannot hold a regrid of size $shaped"))
    raw = dest isa DD.AbstractDimArray ? parent(dest) : dest
    applyplan!(reshape(raw, ndst, prod(othersizes)), plan, src)
    return dest
end

regrid!(dest, data, plan::AbstractRegriddingPlan) =
    error("$(typeof(plan).name.name) defines no `regrid!` application")

"""
    plan_regrid(data; to, from = nothing, method = Conservative(),
                missingpolicy = Weighted(0.5), missingval = sourcemissingval(data),
                lazy = declareschunks(data), chunks = nothing, budget = nothing,
                storage = nothing, sampling = nothing, dependencies = nothing,
                refine = nothing, narrow = nothing) -> AbstractRegriddingPlan

Build a reusable regridding plan without reading source values.

  - Eager plans use one whole-domain [`DirectPlan`](@ref).
  - Lazy plans build blocks on demand and default to a budget-limited
    [`PerChunk`](@ref) cache.

Use `PerChunk()` for an unlimited memory cache or `Spilled(dir)` for disk
storage. Keywords match [`regrid`](@ref). Lazy plans accept `chunks`, `budget`,
`storage`, `dependencies`, `refine`, and `narrow`; eager plans accept `sampling`.

# The chunk dependency relation

A lazy plan owns one chunk dependency relation. `dependencies` selects its
origin:

  - `nothing` or `true` builds a relation;
  - a [`ChunkDependencyGraph`](@ref) adopts and validates that relation;
  - `false` omits the relation.

Every [`LazyRegridArray`](@ref) requires a relation for source selection, tile
order, wave costing, reference counts, and prefetching. `refine(dstchunk,
srcchunk) -> Bool` supplies a conservative narrow phase; `narrow` names that
phase in the relation identity. `refine` must reject only pairs proven
disconnected because a false rejection corrupts results.

[`dependencies`](@ref)`(plan)` returns the relation. The relation remains fixed
for the plan's lifetime; build another plan to use a different one.
"""
function plan_regrid(data; to, from = nothing,
    method::AbstractRegriddingMethod = Conservative(),
    missingpolicy::AbstractMissingPolicy = Weighted(0.5),
    missingval = sourcemissingval(data),
    lazy::Bool = declareschunks(data), chunks = nothing,
    budget::Union{Nothing,Integer} = nothing,
    storage::Union{Nothing,AbstractBlockStorage} = nothing,
    sampling::Union{Nothing,DD.Lookups.Sampling} = nothing,
    dependencies = nothing, refine = nothing,
    narrow::Union{Nothing,Symbol} = nothing)
    src_space = from === nothing ? _sourcespace(data, method) :
                sourcespacefor(from, method)
    from === nothing || checksource(from, data, src_space)
    dst_space = _asspace(to, "to", src_space)
    manifold(dst_space) == manifold(src_space) || throw(ArgumentError(
        "the two sides of a regrid must live on one manifold, but the source " *
        "is on $(manifold(src_space)) and the destination on $(manifold(dst_space))"))
    if !lazy
        _rejectlazykeywords(chunks, budget, storage, dependencies, refine, narrow)
        return eagerplan(method, missingpolicy, dst_space, src_space,
            missingval, sampling)
    end
    sampling === nothing || throw(ArgumentError(
        "a lazy regrid returns an unlabelled disk array, so there is no lookup " *
        "for `sampling` to describe; pass `lazy = false` to label the destination."))
    _checkchunks(chunks)
    budget === nothing || budget > 0 ||
        throw(ArgumentError("budget must be positive, got $budget"))
    return ChunkedPlan(method, missingpolicy, dst_space, src_space;
        storage, budget = something(budget, DEFAULT_BUDGET), chunks, missingval,
        dependencies, refine, narrow)
end

function _rejectlazykeywords(chunks, budget, storage, dependencies, refine, narrow)
    named = String[]
    chunks === nothing || push!(named, "`chunks`")
    budget === nothing || push!(named, "`budget`")
    storage === nothing || push!(named, "`storage`")
    dependencies === nothing || push!(named, "`dependencies`")
    refine === nothing || push!(named, "`refine`")
    narrow === nothing || push!(named, "`narrow`")
    isempty(named) && return nothing
    throw(ArgumentError(
        "an eager plan holds one whole-domain block and takes no " *
        "$(join(named, ", ", " or ")); pass `lazy = true` for the chunked path."))
end

function _rejectplankeywords(kwargs, name::AbstractString)
    named = String[]
    for k in (:dependencies, :refine, :narrow)
        k in keys(kwargs) && push!(named, "`$k`")
    end
    isempty(named) && return nothing
    throw(ArgumentError(
        "`$name` builds a plan, applies it and drops it, so it takes no " *
        "$(join(named, ", ", " or ")): a chunk dependency relation is settled " *
        "when the plan is built and is worth supplying only to a plan that is " *
        "reused. Build it with `plan_regrid` and pass the plan to `$name`."))
end

"""
    wholeblock(method, dst_space, src_space) -> WeightBlock

Build one [`WeightBlock`](@ref) over the full source and destination domains.
The shared [`weightblock`](@ref) path keeps eager and chunk-pair construction
equivalent. A task-local memo is cheaper here because one block offers no later
build with which to share prepared destination geometry.
"""
wholeblock(method::AbstractRegriddingMethod, dst_space::RegridSpace,
    src_space::RegridSpace) =
    weightblock(method, dst_space, 1:Int(ncells(dst_space)),
        src_space, 1:Int(ncells(src_space)))

function _sourcespace(data::DD.AbstractDimArray, method)
    for d in DD.dims(data)
        lookup = DD.lookup(d)
        view = sourceview(lookup, data, method)
        view === nothing || return _presentedspace(view, method)
        named = dimsource(lookup)
        named === nothing && continue
        return sourcespacefor(named, method)
    end
    return RasterGrid(data)
end

function _presentedspace(view::DD.AbstractDimArray, method)
    for d in DD.dims(view)
        named = dimsource(DD.lookup(d))
        named === nothing || return sourcespacefor(named, method)
    end
    throw(ArgumentError(
        "a presented source must name the cells it is written against, but " *
        "$(DD.dims(view)) names none"))
end

_presentedspace(view) = throw(ArgumentError(
    "a presented source must be a dimensional array naming its own cells, " *
    "got a $(typeof(view))"))

_sourcespace(data, method) = throw(ArgumentError(
    "a $(typeof(data)) carries no coordinates, so no source space can be " *
    "derived from it; pass `from = ` a RegridSpace."))

function _asspace(space, name)
    space isa RegridSpace || throw(ArgumentError(
        "`$name` must be a RegridSpace, got $(typeof(space)). A package that " *
        "supplies spaces resolves its own target spellings into one."))
    return space
end

_asspace(space, name, src_space) = _asspace(space, name)

_checkchunks(::Nothing) = nothing
_checkchunks(chunks::Tuple{Vararg{Integer}}) =
    all(>(0), chunks) ? nothing :
    throw(ArgumentError("chunk sizes must be positive, got $chunks"))
_checkchunks(::DiskArrays.GridChunks) = nothing
_checkchunks(chunks) = throw(ArgumentError(
    "`chunks` must be a tuple of chunk sizes, a DiskArrays.GridChunks, or " *
    "nothing, got $(typeof(chunks))"))

# Apply the method-specific source presentation before spatial flattening.
function _flatten(data, plan::AbstractRegriddingPlan)
    data = sourceview(data, plan.method)
    nsrc = Int(ncells(plan.src_space))
    sd = resolvespatialdims(data, nsrc)
    othersizes = _otherdimsizes(data, sd)
    return sd, othersizes, flatsource(data, nsrc, prod(othersizes))
end
