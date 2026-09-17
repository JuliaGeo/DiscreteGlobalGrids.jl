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
           storage = nothing, sampling = nothing, locator = TreeLocator())
    regrid(data, plan::AbstractRegriddingPlan; missingval = outputmissingval(data))

Regrid `data` onto `to`. Spatial dimensions lead in source-cell order;
non-spatial dimensions retain their order. One plan serves every non-spatial
slice.

Dimensional results use the destination axes — a `RasterGrid` echoes the
dimension order it was constructed with — followed by the unchanged non-spatial
dimensions. A destination without axes uses one flat `Cell` axis. Lazy results
preserve the same labels and shape over a disk-backed array.

Results are floating point. [`Weighted`](@ref) blanks uncovered destination
cells with the source's own nodata sentinel (`outputmissingval`): a
`Rasters.AbstractRaster` comes back as a raster declaring the `missingval` it
was handed, and every other array takes `missing` when its element type holds
it and `NaN` otherwise.

# Keyword arguments

  - `to`: destination `RegridSpace`, a dimensional raster or tuple of
    dimensions naming a `RasterGrid`, or a package-specific target.
  - `from`: source space, spelled any of those ways. `nothing` infers a
    self-describing cell axis or a `RasterGrid`.
  - `method`: weight-building method; defaults to [`Conservative`](@ref).
  - `locator`: how candidate cells are found while weights are built, a
    `CandidateLocator`; defaults to [`TreeLocator`](@ref). The weights
    do not depend on it.
  - `missingpolicy`: [`Weighted`](@ref) means or [`Extensive`](@ref) sums.
  - `missingval`: the nodata sentinel of the regrid — invalid in the source, and
    written into blanked destination cells. Left out, the source's comes from
    `sourcemissingval(data)` and the destination's from
    `outputmissingval(data)`, a raster's own `missingval`. `missing`
    and `NaN` are always invalid whatever it is. Give it a value the destination
    element type holds and the result stays concrete: `missingval = NaN` regrids
    a `Union{Missing,Float64}` raster into a `Float64` one.
  - `lazy`: compute on demand (`LazyRegridArray`); defaults to chunked sources.
  - `chunks`: lazy destination tiling. `nothing` derives it automatically.
  - `budget`: target bytes for lazy reads and weights, default `2^31`.
  - `storage`: lazy weight storage, [`PerChunk`](@ref) or [`Spilled`](@ref).
  - `sampling`: destination lookup sampling. `nothing` follows the method —
    area-based methods give `Intervals`, point samples give `Points`
    (`outputsampling`).

Keyword applicability:

  - Lazy regrids accept `chunks`, `budget`, and `storage`.
  - Eager regrids accept `sampling`.
  - The plan form takes `missingval` alone: a plan settles how weights are
    built, and the sentinel is what the caller does with them.

The keyword form delegates to [`plan_regrid`](@ref). Dependency-relation
keywords belong to reusable plans and are accepted directly by `plan_regrid`.
"""
function regrid end

function regrid(data; kwargs...)
    _rejectplankeywords(kwargs, "regrid")
    haskey(kwargs, :missingval) || return regrid(data, plan_regrid(data; kwargs...))
    mv = kwargs[:missingval]
    plan = plan_regrid(data; kwargs..., missingval = _sourcesentinel(mv))
    return regrid(data, plan; missingval = mv)
end

function regrid(data, plan::DirectPlan; missingval = outputmissingval(data))
    sd, othersizes, src = _flatten(data, plan)
    ndst = size(plan.block, 1)
    out = Array{outputeltype(eltype(data), missingval)}(undef, ndst, othersizes...)
    applyplan!(reshape(out, ndst, prod(othersizes)), plan, src, missingval)
    return wrapoutput(out, data, sd, destinationdims(plan), missingval)
end

regrid(data, plan::AbstractRegriddingPlan; missingval = outputmissingval(data)) =
    error("$(typeof(plan).name.name) defines no `regrid` application")

# The source half of a `missingval`. `missing` is invalid wherever it appears, so
# declaring it as a sentinel would only cost `anyinvalid` a scan it can skip.
_sourcesentinel(missingval) = missingval
_sourcesentinel(::Missing) = nothing

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
            missingpolicy = Weighted(0.5), missingval = destinationmissingval(dest),
            lazy = declareschunks(data), chunks = nothing, budget = nothing,
            storage = nothing, sampling = nothing)
    regrid!(dest, data, plan::AbstractRegriddingPlan;
            missingval = destinationmissingval(dest))

Regrid `data` into the preallocated `dest` and return `dest`.

`dest` starts with the destination axes or one flat cell dimension, followed by
the non-spatial dimensions of `data`; either leading shape is accepted. Keywords
match [`regrid`](@ref) and delegate to [`plan_regrid`](@ref).

`dest` declares the destination's nodata convention here, so `missingval`
defaults to `destinationmissingval(dest)` — its own `missingval` for a
`Rasters.AbstractRaster`, and `missing` or NaN for a plain array. Passing one
names the sentinel on both sides, as it does for [`regrid`](@ref), and it must
be a value `eltype(dest)` holds.
"""
function regrid! end

function regrid!(dest, data; kwargs...)
    _rejectplankeywords(kwargs, "regrid!")
    haskey(kwargs, :missingval) ||
        return regrid!(dest, data, plan_regrid(data; kwargs...))
    mv = kwargs[:missingval]
    plan = plan_regrid(data; kwargs..., missingval = _sourcesentinel(mv))
    return regrid!(dest, data, plan; missingval = mv)
end

function regrid!(dest, data, plan::DirectPlan;
    missingval = destinationmissingval(dest))
    _, othersizes, src = _flatten(data, plan)
    ndst = size(plan.block, 1)
    dstdims = destinationdims(plan)
    shaped = dstdims === nothing ? (ndst, othersizes...) :
             (map(length, dstdims)..., othersizes...)
    size(dest) == shaped || size(dest) == (ndst, othersizes...) ||
        throw(DimensionMismatch(
            "destination of size $(size(dest)) cannot hold a regrid of size $shaped"))
    raw = dest isa DD.AbstractDimArray ? parent(dest) : dest
    applyplan!(reshape(raw, ndst, prod(othersizes)), plan, src, missingval)
    return dest
end

regrid!(dest, data, plan::AbstractRegriddingPlan;
    missingval = destinationmissingval(dest)) =
    error("$(typeof(plan).name.name) defines no `regrid!` application")

"""
    plan_regrid(data; to, from = nothing, method = Conservative(),
                missingpolicy = Weighted(0.5), missingval = sourcemissingval(data),
                lazy = declareschunks(data), chunks = nothing, budget = nothing,
                storage = nothing, sampling = nothing, dependencies = nothing,
                refine = nothing, narrow = nothing) -> AbstractRegriddingPlan

Build a reusable regridding plan without reading source values. `missingval` is
the source sentinel alone here — a plan reads data and never writes it, so the
destination's sentinel belongs to [`regrid`](@ref).

  - Eager plans use one whole-domain `DirectPlan`.
  - Lazy plans build blocks on demand and default to a budget-limited
    [`PerChunk`](@ref) cache.

Use `PerChunk()` for an unlimited memory cache or `Spilled(dir)` for disk
storage. Keywords match [`regrid`](@ref). Lazy plans accept `chunks`, `budget`,
`storage`, `dependencies`, `refine`, and `narrow`; eager plans accept `sampling`.

# The chunk dependency relation

A lazy plan owns one chunk dependency relation. `dependencies` selects its
origin:

  - `nothing` or `true` builds a relation;
  - a `ChunkDependencyGraph` adopts and validates that relation;
  - `false` omits the relation.

Every `LazyRegridArray` requires a relation for source selection, tile
order, wave costing, reference counts, and prefetching. `refine(dstchunk,
srcchunk) -> Bool` supplies a conservative narrow phase; `narrow` names that
phase in the relation identity. `refine` must reject only pairs proven
disconnected because a false rejection corrupts results.

`dependencies(plan)` returns the relation. The relation remains fixed
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
    narrow::Union{Nothing,Symbol} = nothing,
    locator::CandidateLocator = TreeLocator())
    if from === nothing
        src_space = _sourcespace(data, method)
    else
        src_space = sourcespacefor(from, method)
        checksource(from, data, src_space)
    end
    dst_space = _asspace(to, "to", src_space)
    manifold(dst_space) == manifold(src_space) || throw(ArgumentError(
        "the two sides of a regrid must live on one manifold, but the source " *
        "is on $(manifold(src_space)) and the destination on $(manifold(dst_space))"))
    if !lazy
        _rejectlazykeywords(chunks, budget, storage, dependencies, refine, narrow)
        return eagerplan(method, missingpolicy, dst_space, src_space,
            missingval, sampling, locator)
    end
    sampling === nothing || throw(ArgumentError(
        "a lazy regrid returns an unlabelled disk array, so there is no lookup " *
        "for `sampling` to describe; pass `lazy = false` to label the destination."))
    _checkchunks(chunks)
    budget === nothing || budget > 0 ||
        throw(ArgumentError("budget must be positive, got $budget"))
    return ChunkedPlan(method, missingpolicy, dst_space, src_space;
        storage, budget = something(budget, DEFAULT_BUDGET), chunks, missingval,
        dependencies, refine, narrow, locator)
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
    wholeblock(method, dst_space, src_space, locator = TreeLocator()) -> WeightBlock

Build one [`WeightBlock`](@ref) over the full source and destination domains.
The shared [`weightblock`](@ref) path keeps eager and chunk-pair construction
equivalent. A task-local memo is cheaper here because one block offers no later
build with which to share prepared destination geometry.
"""
wholeblock(method::AbstractRegriddingMethod, dst_space::RegridSpace,
    src_space::RegridSpace, locator::CandidateLocator = TreeLocator()) =
    weightblock(method, dst_space, 1:Int(ncells(dst_space)),
        src_space, 1:Int(ncells(src_space)), locator)

# Only dimensional arrays carry enough geometry to infer a source space. A
# presented view or a cell-naming dimension wins over raster inference.
function _sourcespace(data::DD.AbstractDimArray, method)
    for d in DD.dims(data)
        lookup = DD.lookup(d)
        view = sourceview(lookup, data, method)
        view === nothing || return _presentedspace(view, method)
        named = dimsource(lookup)
        named === nothing || return sourcespacefor(named, method)
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

_presentedspace(view, method) = throw(ArgumentError(
    "a presented source must be a dimensional array naming its own cells, " *
    "got a $(typeof(view))"))

_sourcespace(data, method) = throw(ArgumentError(
    "a $(typeof(data)) carries no coordinates, so no source space can be " *
    "derived from it; pass `from = ` a RegridSpace."))

# A dimension that already names cells ([`dimsource`](@ref)) is not a raster
# axis, so name the source it holds rather than ask for the raster axis it does
# not have. Every route that reads a lattice off dimensions comes through here:
# `context` says which route it was, and `name` is the keyword whose spelling
# the caller has to change.
function _checkrasterdims(ds, name::AbstractString, context::AbstractString)
    for d in ds
        named = dimsource(DD.lookup(d))
        named === nothing && continue
        throw(ArgumentError("""
        $context, but its $(DD.name(d)) dimension names cells rather than a \
        raster lattice. Pass `$name = $(named)`.
        """))
    end
    return nothing
end

# A dimensional raster, or the bare dimensions of one, stands for the lattice it
# carries. Either side of a regrid may be spelled that way, so putting a result
# back on the axes it came from needs no `RasterGrid` written out by hand.
function _asspace(A::DD.AbstractDimArray, name::AbstractString)
    _checkrasterdims(DD.dims(A), name, "`$name` was given a dimensional raster")
    return RasterGrid(A)
end

function _asspace(ds::Tuple{Vararg{DD.Dimension}}, name::AbstractString)
    _checkrasterdims(ds, name, "`$name` was given a dimension tuple")
    return RasterGrid(ds)
end

function _asspace(space, name)
    space isa RegridSpace || throw(ArgumentError(
        "`$name` must be a RegridSpace, a dimensional raster or a tuple of " *
        "dimensions, got $(typeof(space)). A package that supplies spaces " *
        "resolves its own target spellings into one."))
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
