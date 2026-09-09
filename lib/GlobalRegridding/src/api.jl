# Public regridding API.

# Whether `data` declares a chunking of its own: what makes a regrid lazy by
# default, what `SourceChunking` reads, and what `flatsource` materializes before
# reshaping. It says nothing about residence; `_isdisksource` tests that.
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
    regrid(data, plan::AbstractRegriddingPlan; missingval = outputmissingval(data))

Regrid `data` onto `to`. Spatial dimensions must come first and flatten in the
source space's cell order. Non-spatial dimensions retain their order. One plan
is reused for all non-spatial slices.

A dimensional source comes back labelled with the destination's own axes — a
[`RasterGrid`](@ref) echoes the dimension order it was constructed with —
followed by its unchanged non-spatial dimensions. Destinations without axes of
their own keep a flat `Cell` axis. Lazy results carry the same labels and
shape over a disk-backed array.

Results are floating point. [`Weighted`](@ref) blanks uncovered destination
cells with the source's own nodata sentinel ([`outputmissingval`](@ref)): a
`Rasters.AbstractRaster` comes back as a raster declaring the `missingval` it
was handed, and every other array takes `missing` when its element type holds
it and `NaN` otherwise.

# Keyword arguments

  - `to`: destination [`RegridSpace`](@ref), a dimensional raster or tuple of
    dimensions naming a [`RasterGrid`](@ref), or a package-specific target.
  - `from`: source space, spelled any of those ways; `nothing` derives a
    [`RasterGrid`](@ref) from `data`.
  - `method`: weight-building method; defaults to [`Conservative`](@ref).
  - `missingpolicy`: [`Weighted`](@ref) means or [`Extensive`](@ref) sums.
  - `missingval`: the nodata sentinel of the regrid — invalid in the source, and
    written into blanked destination cells. Left out, the source's comes from
    [`sourcemissingval`](@ref)`(data)` and the destination's from
    [`outputmissingval`](@ref)`(data)`, a raster's own `missingval`. `missing`
    and `NaN` are always invalid whatever it is. Give it a value the destination
    element type holds and the result stays concrete: `missingval = NaN` regrids
    a `Union{Missing,Float64}` raster into a `Float64` one.
  - `lazy`: compute on demand ([`LazyRegridArray`](@ref)); defaults to chunked sources.
  - `chunks`: lazy destination tiling. `nothing` derives it automatically.
  - `budget`: target bytes for lazy reads and weights, default `2^31`.
  - `storage`: lazy weight storage, [`PerChunk`](@ref) or [`Spilled`](@ref).
  - `sampling`: destination lookup sampling. `nothing` follows the method —
    area-based methods give `Intervals`, point samples give `Points`
    ([`outputsampling`](@ref)).

`chunks`, `budget` and `storage` apply only to `lazy = true`, and `sampling`
only to `lazy = false`. The plan form takes `missingval` alone: a plan settles
how weights are built, and the sentinel is what the caller does with them.

Every other keyword above is [`plan_regrid`](@ref)'s and is forwarded to it, so
each default and each check is stated there once. The relation keywords
`dependencies`, `refine` and `narrow` describe a plan that is kept and are
refused here.
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

Return the dimensions labelling this plan's results, under the plan's own
`sampling` when it declares one and the method's otherwise. Chunked plans
declare no sampling, so the method's always applies.
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

`dest` starts with the destination's own axes, or one flat cell dimension,
followed by `data`'s non-spatial dimensions; either leading shape is accepted.
Keywords match [`regrid`](@ref) and are forwarded to [`plan_regrid`](@ref).

`dest` declares the destination's nodata convention here, so `missingval`
defaults to [`destinationmissingval`](@ref)`(dest)` — its own `missingval` for a
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
destination's sentinel belongs to [`regrid`](@ref). In-memory data uses one
whole-domain [`DirectPlan`](@ref). Lazy plans build blocks on demand
and default to a budget-limited [`PerChunk`](@ref) cache. Use `PerChunk()` for
an unlimited cache or `Spilled(dir)` for disk storage. Keywords match
[`regrid`](@ref); `chunks`, `budget`, `storage`, `dependencies`, `refine` and
`narrow` apply only to `lazy = true`, and `sampling` only to `lazy = false`.

# The chunk dependency relation

A lazy plan is the sole owner of its chunk dependency relation, and this is the
only place a narrow phase may be supplied. `dependencies` chooses whether the
plan builds one (`nothing`, the default, or `true`), adopts and validates one
somebody else built (a [`ChunkDependencyGraph`](@ref)), or holds none (`false`).
Every lazy read needs one — for tile order, wave costing, refcounts and
prefetch, and on the chunk-pair route for the source chunks themselves — so a
plan that holds none cannot back a [`LazyRegridArray`](@ref). `refine` is the
conservative narrow phase to apply while building, `refine(dstchunk, srcchunk)
-> Bool`, and `narrow` the `Symbol` that names it in the relation's identity. A
`refine` must only ever reject pairs it can *prove* disconnected; a wrong one
silently corrupts results. [`dependencies`](@ref) documents each branch.

[`dependencies`](@ref)`(plan)` reads the relation back and builds nothing. It is
deliberately impossible to narrow, replace or rebuild a plan's relation once the
plan exists: [`regrid`](@ref) and [`regrid!`](@ref) forward every other keyword
here but refuse `dependencies`, `refine` and `narrow`, and
[`chunk_dependency_graph`](@ref) has no `plan` method. A caller that wants a
different relation makes a different plan.
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
    src_space = from === nothing ? _sourcespace(data) : _asspace(from, "from")
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

# The three keywords that describe a plan somebody keeps: a relation to adopt,
# the narrow phase to build it with, and the name that phase goes by. A one-shot
# regrid builds its plan and drops it, so there is nothing for them to describe.
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

Build one [`WeightBlock`](@ref) over all source and destination cells. The build
path is [`weightblock`](@ref)'s, so the eager domain and a chunk pair are built
the same way.

The whole domain is one block, so it prepares no destination geometry
([`preparedestination`](@ref)): with no second block to share it, a task-local
memo is cheaper than a slot per destination cell.
"""
wholeblock(method::AbstractRegriddingMethod, dst_space::RegridSpace,
    src_space::RegridSpace) =
    weightblock(method, dst_space, 1:Int(ncells(dst_space)),
        src_space, 1:Int(ncells(src_space)))

# Only dimensional arrays carry enough geometry to infer a source space.
function _sourcespace(data::DD.AbstractDimArray)
    _checkrasterdims(DD.dims(data), "from",
        "no `from` was given, so the source space was derived from the data")
    return RasterGrid(data)
end

_sourcespace(data) = throw(ArgumentError(
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

# Flatten the source to `ncells × nslices` and retain pass-through sizes.
function _flatten(data, plan::AbstractRegriddingPlan)
    nsrc = Int(ncells(plan.src_space))
    sd = resolvespatialdims(data, nsrc)
    othersizes = _otherdimsizes(data, sd)
    return sd, othersizes, flatsource(data, nsrc, prod(othersizes))
end
