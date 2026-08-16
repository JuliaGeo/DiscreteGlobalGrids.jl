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
           lazy = isdiskbacked(data), chunks = nothing, budget = 2^30,
           storage = nothing)
    regrid(data, plan::AbstractRegriddingPlan)

Regrid `data` onto `to`. Spatial dimensions must come first and flatten in the
source space's cell order. Non-spatial dimensions retain their order. One plan
is reused for all non-spatial slices.

Results are floating point. [`Weighted`](@ref) writes `missing` when supported
by the source element type, and `NaN` otherwise.

# Keyword arguments

  - `to`: destination [`RegridSpace`](@ref) or a package-specific target.
  - `from`: source space; `nothing` derives a [`RasterGrid`](@ref) from `data`.
  - `method`: weight-building method; defaults to [`Conservative`](@ref).
  - `missingpolicy`: [`Weighted`](@ref) means or [`Extensive`](@ref) sums.
  - `missingval`: additional nodata sentinel. `missing` and `NaN` are always invalid.
  - `lazy`: return a [`LazyRegridArray`](@ref); defaults to chunked sources.
  - `chunks`: lazy destination tiling. `nothing` derives it automatically.
  - `budget`: target bytes for transient lazy-read data and weights.
  - `storage`: lazy weight storage, [`PerChunk`](@ref) or [`Spilled`](@ref).

The plan form accepts no keywords because the plan contains all settings.
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

`dest` starts with the destination cell dimension, followed by `data`'s
non-spatial dimensions. Keywords match [`regrid`](@ref); the plan form takes none.
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

Build a reusable regridding plan without reading source values. In-memory data
uses one whole-domain [`DirectPlan`](@ref). Lazy plans build blocks on demand
and default to a budget-limited [`PerChunk`](@ref) cache. Use `PerChunk()` for
an unlimited cache or `Spilled(dir)` for disk storage. Keywords match
[`regrid`](@ref).
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

# Flatten the source to `ncells × nslices` and retain pass-through sizes.
function _flatten(data, plan::AbstractRegriddingPlan)
    nsrc = Int(ncells(plan.src_space))
    sd = resolvespatialdims(data, nsrc)
    othersizes = _otherdimsizes(data, sd)
    return sd, othersizes, flatsource(data, nsrc, prod(othersizes))
end
