# Public regridding API.

# Whether the source exposes chunked storage.
isdiskbacked(data) = DiskArrays.haschunks(data) isa DiskArrays.Chunked

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
           lazy = isdiskbacked(data), chunks = nothing, budget = nothing,
           storage = nothing, sampling = nothing)
    regrid(data, plan::AbstractRegriddingPlan)

Regrid `data` onto `to`. Spatial dimensions must come first and flatten in the
source space's cell order. Non-spatial dimensions retain their order. One plan
is reused for all non-spatial slices.

A dimensional source comes back labelled with the destination's own axes — a
[`RasterGrid`](@ref) echoes the dimension order it was constructed with —
followed by its unchanged non-spatial dimensions. Destinations without axes of
their own keep a flat `Cell` axis. Lazy results carry the same labels and
shape over a disk-backed array.

Results are floating point. [`Weighted`](@ref) writes `missing` when supported
by the source element type, and `NaN` otherwise.

# Keyword arguments

  - `to`: destination [`RegridSpace`](@ref) or a package-specific target.
  - `from`: source space; `nothing` derives a [`RasterGrid`](@ref) from `data`.
  - `method`: weight-building method; defaults to [`Conservative`](@ref).
  - `missingpolicy`: [`Weighted`](@ref) means or [`Extensive`](@ref) sums.
  - `missingval`: additional nodata sentinel. `missing` and `NaN` are always invalid.
  - `lazy`: compute on demand ([`LazyRegridArray`](@ref)); defaults to chunked sources.
  - `chunks`: lazy destination tiling. `nothing` derives it automatically.
  - `budget`: target bytes for lazy reads and weights, default `2^30`.
  - `storage`: lazy weight storage, [`PerChunk`](@ref) or [`Spilled`](@ref).
  - `sampling`: destination lookup sampling. `nothing` follows the method —
    area-based methods give `Intervals`, point samples give `Points`
    ([`outputsampling`](@ref)).

`chunks`, `budget` and `storage` apply only to `lazy = true`, and `sampling`
only to `lazy = false`. The plan form accepts no keywords because the plan
contains all settings.
"""
function regrid end

function regrid(data; to, from = nothing,
    method::AbstractRegriddingMethod = Conservative(),
    missingpolicy::AbstractMissingPolicy = Weighted(0.5),
    missingval = sourcemissingval(data),
    lazy::Bool = isdiskbacked(data), chunks = nothing,
    budget::Union{Nothing,Integer} = nothing,
    storage::Union{Nothing,AbstractBlockStorage} = nothing,
    sampling::Union{Nothing,DD.Lookups.Sampling} = nothing)
    plan = plan_regrid(data; to, from, method, missingpolicy, missingval, lazy,
        chunks, budget, storage, sampling)
    return regrid(data, plan)
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
            missingpolicy = Weighted(0.5), missingval = sourcemissingval(data),
            lazy = isdiskbacked(data), chunks = nothing, budget = nothing,
            storage = nothing, sampling = nothing)
    regrid!(dest, data, plan::AbstractRegriddingPlan)

Regrid `data` into the preallocated `dest` and return `dest`.

`dest` starts with the destination's own axes, or one flat cell dimension,
followed by `data`'s non-spatial dimensions; either leading shape is accepted.
Keywords match [`regrid`](@ref); the plan form takes none.
"""
function regrid! end

function regrid!(dest, data; to, from = nothing,
    method::AbstractRegriddingMethod = Conservative(),
    missingpolicy::AbstractMissingPolicy = Weighted(0.5),
    missingval = sourcemissingval(data),
    lazy::Bool = isdiskbacked(data), chunks = nothing,
    budget::Union{Nothing,Integer} = nothing,
    storage::Union{Nothing,AbstractBlockStorage} = nothing,
    sampling::Union{Nothing,DD.Lookups.Sampling} = nothing)
    plan = plan_regrid(data; to, from, method, missingpolicy, missingval, lazy,
        chunks, budget, storage, sampling)
    return regrid!(dest, data, plan)
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
                lazy = isdiskbacked(data), chunks = nothing, budget = nothing,
                storage = nothing, sampling = nothing) -> AbstractRegriddingPlan

Build a reusable regridding plan without reading source values. In-memory data
uses one whole-domain [`DirectPlan`](@ref). Lazy plans build blocks on demand
and default to a budget-limited [`PerChunk`](@ref) cache. Use `PerChunk()` for
an unlimited cache or `Spilled(dir)` for disk storage. Keywords match
[`regrid`](@ref); `chunks`, `budget` and `storage` apply only to `lazy = true`,
and `sampling` only to `lazy = false`.
"""
function plan_regrid(data; to, from = nothing,
    method::AbstractRegriddingMethod = Conservative(),
    missingpolicy::AbstractMissingPolicy = Weighted(0.5),
    missingval = sourcemissingval(data),
    lazy::Bool = isdiskbacked(data), chunks = nothing,
    budget::Union{Nothing,Integer} = nothing,
    storage::Union{Nothing,AbstractBlockStorage} = nothing,
    sampling::Union{Nothing,DD.Lookups.Sampling} = nothing)
    src_space = from === nothing ? _sourcespace(data) : _asspace(from, "from")
    dst_space = _asspace(to, "to", src_space)
    manifold(dst_space) == manifold(src_space) || throw(ArgumentError(
        "the two sides of a regrid must live on one manifold, but the source " *
        "is on $(manifold(src_space)) and the destination on $(manifold(dst_space))"))
    if !lazy
        _rejectlazykeywords(chunks, budget, storage)
        return DirectPlan(method, missingpolicy, dst_space, src_space,
            wholeblock(method, dst_space, src_space), missingval, sampling)
    end
    sampling === nothing || throw(ArgumentError(
        "a lazy regrid returns an unlabelled disk array, so there is no lookup " *
        "for `sampling` to describe; pass `lazy = false` to label the destination."))
    _checkchunks(chunks)
    budget === nothing || budget > 0 ||
        throw(ArgumentError("budget must be positive, got $budget"))
    return ChunkedPlan(method, missingpolicy, dst_space, src_space;
        storage, budget = something(budget, 2^30), chunks, missingval)
end

function _rejectlazykeywords(chunks, budget, storage)
    named = String[]
    chunks === nothing || push!(named, "`chunks`")
    budget === nothing || push!(named, "`budget`")
    storage === nothing || push!(named, "`storage`")
    isempty(named) && return nothing
    throw(ArgumentError(
        "an eager plan holds one whole-domain block and takes no " *
        "$(join(named, ", ", " or ")); pass `lazy = true` for the chunked path."))
end

"""
    wholeblock(method, dst_space, src_space) -> WeightBlock

Build one [`WeightBlock`](@ref) over all source and destination cells.
"""
function wholeblock(method::AbstractRegriddingMethod, dst_space::RegridSpace,
    src_space::RegridSpace)
    ndst = Int(ncells(dst_space))
    nsrc = Int(ncells(src_space))
    coo = WeightCOO(ndst)
    build_weights!(coo, method, dst_space, 1:ndst, src_space, 1:nsrc)
    return WeightBlock(coo, ndst, nsrc)
end

# Only dimensional arrays carry enough geometry to infer a source space.
# A dimension that already names cells ([`dimsource`](@ref)) is not a raster
# axis, so name it rather than ask for the raster axis it does not have.
function _sourcespace(data::DD.AbstractDimArray)
    for d in DD.dims(data)
        named = dimsource(DD.lookup(d))
        named === nothing && continue
        throw(ArgumentError("""
        no `from` was given, so the source space was derived from the data, but \
        its $(DD.name(d)) dimension names cells rather than a raster lattice. \
        Pass `from = $(named)`.
        """))
    end
    return RasterGrid(data)
end

_sourcespace(data) = throw(ArgumentError(
    "a $(typeof(data)) carries no coordinates, so no source space can be " *
    "derived from it; pass `from = ` a RegridSpace."))

"""
    _asspace(space, name) -> RegridSpace
    _asspace(space, name, src_space) -> RegridSpace

Resolve a `to` or `from` argument into a [`RegridSpace`](@ref). Packages that
supply spaces extend the two-argument form for their own target spellings, and
the three-argument form when the destination depends on the resolved source
space. `name` names the keyword in error messages.
"""
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

# Flatten the source to `ncells × nslices` and retain pass-through sizes.
function _flatten(data, plan::AbstractRegriddingPlan)
    nsrc = Int(ncells(plan.src_space))
    sd = resolvespatialdims(data, nsrc)
    othersizes = _otherdimsizes(data, sd)
    return sd, othersizes, flatsource(data, nsrc, prod(othersizes))
end
