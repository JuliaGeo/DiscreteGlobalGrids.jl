# Shared accumulation, finalization, and array-shaping logic.

# Validity

# `missing` and NaN are invalid by default.
@inline _isvalid(::Missing) = false
@inline _isvalid(x::AbstractFloat) = !isnan(x)
@inline _isvalid(::Any) = true

# Invalid entries contribute zero.
@inline _value(::Missing) = 0.0
@inline _value(x::AbstractFloat) = isnan(x) ? 0.0 : Float64(x)
@inline _value(x::Real) = Float64(x)

# Valid entries contribute one to coverage; booleans act as explicit masks.
@inline _validity(::Missing) = 0.0
@inline _validity(x::Bool) = x ? 1.0 : 0.0
@inline _validity(x::AbstractFloat) = isnan(x) ? 0.0 : 1.0
@inline _validity(::Real) = 1.0

# `isequal` keeps sentinel comparisons total for `missing` values.
@inline _isvalid(x, ::Nothing) = _isvalid(x)
@inline _isvalid(x, missingval) = _isvalid(x) && !isequal(x, missingval)
@inline _value(x, ::Nothing) = _value(x)
@inline _value(x, missingval) = isequal(x, missingval) ? 0.0 : _value(x)
@inline _validity(x, ::Nothing) = _validity(x)
@inline _validity(x, missingval) = isequal(x, missingval) ? 0.0 : _validity(x)

"""
    knownempty(data, ndchunk::Tuple{Vararg{AbstractUnitRange}}) -> Bool

Return whether N-D storage chunk `ndchunk` contains no valid values, without
reading it. Defaults to `false`. A `true` result permits skipping reads and
matvecs; zeros are valid data and must not be reported empty. Backend adapters
may specialize this hook.
"""
knownempty(data, ndchunk::Tuple{Vararg{AbstractUnitRange}}) = false

# Skip validity scans for element types that cannot hold invalid values.
_canbeinvalid(::Type{T}) where {T} = !(T <: Real) || T <: AbstractFloat || Missing <: T

"""
    anyinvalid(src, missingval = nothing) -> Bool

Return whether `src` contains an invalid value. Avoid scanning when its element
type cannot contain one and no sentinel is declared.
"""
anyinvalid(src::AbstractArray) =
    _canbeinvalid(eltype(src)) ? any(!_isvalid, src) : false

anyinvalid(src::AbstractArray, ::Nothing) = anyinvalid(src)
anyinvalid(src::AbstractArray, missingval) = any(x -> !_isvalid(x, missingval), src)

# Use `missing` when supported and NaN otherwise.
@inline function _maskedvalue(::Type{T}) where {T}
    if Missing <: T
        return missing
    elseif T <: AbstractFloat
        return T(NaN)
    else
        throw(ArgumentError(
            "`Weighted` blanks destinations below its coverage threshold and " *
            "needs a sentinel to do it, but $T holds neither `missing` nor NaN"))
    end
end

# Accumulation

"""
    blockreference!(ref, block::WeightBlock) -> ref

Write each destination's reference weight to `ref`, using stored denominators or
row sums. The result is data-independent and reusable.
"""
function blockreference!(ref::AbstractVector{Float64}, block::WeightBlock)
    fill!(ref, 0.0)
    return addreference!(ref, block)
end

"""
    addreference!(ref, block::WeightBlock) -> ref

Add `block`'s per-destination reference weight to `ref`.
"""
function addreference!(ref::AbstractVector{Float64}, block::WeightBlock)
    d = block.denom
    if d === nothing
        _addrowsums!(ref, block.weights)
    else
        @inbounds for j in eachindex(ref, d)
            ref[j] += d[j]
        end
    end
    return ref
end

function _addrowsums!(ref::AbstractVector{Float64}, W::SparseMatrixCSC)
    rows = SparseArrays.rowvals(W)
    vals = SparseArrays.nonzeros(W)
    @inbounds for p in eachindex(rows, vals)
        ref[rows[p]] += vals[p]
    end
    return ref
end

function _addrowsums!(ref::AbstractVector{Float64}, W::AbstractMatrix)
    @inbounds for k in axes(W, 2), j in axes(W, 1)
        ref[j] += W[j, k]
    end
    return ref
end

"""
    applyblock!(num, cover, block::WeightBlock, src, valid = nothing, ref = nothing,
                missingval = nothing)

Accumulate one source chunk into weighted sums `num` and valid coverage `cover`.
`src` follows source-chunk cell order. `valid = nothing` asserts all values are
valid and uses `ref` when supplied; otherwise `valid` is a same-length mask or
the source itself. `missingval` adds a nodata sentinel. This function allocates
nothing.
"""
function applyblock!(num::AbstractVector{Float64}, cover::AbstractVector{Float64},
    block::WeightBlock, src::AbstractVector,
    valid::Union{Nothing,AbstractVector} = nothing,
    ref::Union{Nothing,AbstractVector{Float64}} = nothing,
    missingval = nothing)
    W = block.weights
    Base.require_one_based_indexing(num, cover, src)
    size(W, 1) == length(num) == length(cover) || throw(DimensionMismatch(
        "block has $(size(W, 1)) destinations, accumulators have " *
        "$(length(num)) and $(length(cover))"))
    size(W, 2) == length(src) || throw(DimensionMismatch(
        "block expects $(size(W, 2)) source cells, got $(length(src))"))
    if valid === nothing
        _accumulate!(num, W, src, missingval)
        if ref === nothing
            addreference!(cover, block)
        else
            @inbounds for j in eachindex(cover, ref)
                cover[j] += ref[j]
            end
        end
    else
        Base.require_one_based_indexing(valid)
        length(valid) == length(src) || throw(DimensionMismatch(
            "validity mask has $(length(valid)) entries, source has $(length(src))"))
        _accumulate!(num, cover, W, src, valid, missingval)
    end
    return num
end

# Very sparse blocks locate columns from nonzeros; denser blocks scan columns.
@inline _walknonzeros(W::SparseMatrixCSC) =
    SparseArrays.nnz(W) * 256 < size(W, 2)

# Locate the column owning nonzero position `p`.
@inline _columnof(cols::Vector{Int}, p::Int) = searchsortedlast(cols, p)

# Accumulate values when all entries are valid.
function _accumulate!(num::AbstractVector{Float64}, W::SparseMatrixCSC,
    src::AbstractVector, missingval = nothing)
    rows = SparseArrays.rowvals(W)
    vals = SparseArrays.nonzeros(W)
    cols = SparseArrays.getcolptr(W)
    if _walknonzeros(W)
        p = 1
        nz = SparseArrays.nnz(W)
        @inbounds while p <= nz
            k = _columnof(cols, p)
            p2 = cols[k+1] - 1
            x = _value(src[k], missingval)
            if !iszero(x)
                for q in p:p2
                    num[rows[q]] += vals[q] * x
                end
            end
            p = p2 + 1
        end
        return num
    end
    @inbounds for k in axes(W, 2)
        p1, p2 = cols[k], cols[k+1] - 1
        p1 > p2 && continue
        x = _value(src[k], missingval)
        iszero(x) && continue
        for p in p1:p2
            num[rows[p]] += vals[p] * x
        end
    end
    return num
end

function _accumulate!(num::AbstractVector{Float64}, W::AbstractMatrix,
    src::AbstractVector, missingval = nothing)
    @inbounds for k in axes(W, 2)
        x = _value(src[k], missingval)
        iszero(x) && continue
        for j in axes(W, 1)
            num[j] += W[j, k] * x
        end
    end
    return num
end

# Accumulate values and coverage in one pass.
function _accumulate!(num::AbstractVector{Float64}, cover::AbstractVector{Float64},
    W::SparseMatrixCSC, src::AbstractVector, valid::AbstractVector, missingval = nothing)
    rows = SparseArrays.rowvals(W)
    vals = SparseArrays.nonzeros(W)
    cols = SparseArrays.getcolptr(W)
    if _walknonzeros(W)
        p = 1
        nz = SparseArrays.nnz(W)
        @inbounds while p <= nz
            k = _columnof(cols, p)
            p2 = cols[k+1] - 1
            x = _value(src[k], missingval)
            m = _validity(valid[k], missingval)
            if !(iszero(x) && iszero(m))
                for q in p:p2
                    j = rows[q]
                    v = vals[q]
                    num[j] += v * x
                    cover[j] += v * m
                end
            end
            p = p2 + 1
        end
        return num
    end
    @inbounds for k in axes(W, 2)
        p1, p2 = cols[k], cols[k+1] - 1
        p1 > p2 && continue
        x = _value(src[k], missingval)
        m = _validity(valid[k], missingval)
        (iszero(x) && iszero(m)) && continue
        for p in p1:p2
            j = rows[p]
            v = vals[p]
            num[j] += v * x
            cover[j] += v * m
        end
    end
    return num
end

function _accumulate!(num::AbstractVector{Float64}, cover::AbstractVector{Float64},
    W::AbstractMatrix, src::AbstractVector, valid::AbstractVector, missingval = nothing)
    @inbounds for k in axes(W, 2)
        x = _value(src[k], missingval)
        m = _validity(valid[k], missingval)
        (iszero(x) && iszero(m)) && continue
        for j in axes(W, 1)
            v = W[j, k]
            num[j] += v * x
            cover[j] += v * m
        end
    end
    return num
end

# Finalization

"""
    finalize!(out, num, cover, total, policy) -> out

Finalize one destination chunk after all source blocks have accumulated.

`num` is the weighted sum with invalid sources contributing zero, `cover` the
weight of the sources that were valid, and `total` the accumulated reference
weight from [`blockreference!`](@ref) — the denominator the method reported, or
the row sums when it reported none.

  - [`Extensive`](@ref): `num`.
  - [`Weighted`](@ref)`(t)`: `num / cover`, blank where `cover ≤ 0` or
    `cover < t · total`.

`Weighted` always normalizes by valid coverage and applies its threshold against
`total`. Blanked values use `missing` when supported and NaN otherwise.
`Extensive` returns raw sums and never blanks.
"""
function finalize! end

"""
    usesreference(policy::AbstractMissingPolicy) -> Bool

Return whether `policy` reads the accumulated reference weight. When false,
known-empty chunks can be skipped before their weights are built. When true the
skip is unsafe: a dropped pair's reference weight would be missing from the
coverage [`Weighted`](@ref) blanks against, changing answers.
"""
usesreference(::Extensive) = false
usesreference(policy::Weighted) = policy.threshold > 0

function finalize!(out::AbstractVector, num::AbstractVector{Float64},
    cover::AbstractVector{Float64}, total::AbstractVector{Float64},
    ::Extensive)
    @inbounds for j in eachindex(out, num)
        out[j] = num[j]
    end
    return out
end

function finalize!(out::AbstractVector, num::AbstractVector{Float64},
    cover::AbstractVector{Float64}, total::AbstractVector{Float64},
    policy::Weighted)
    blank = _maskedvalue(eltype(out))
    t = policy.threshold
    @inbounds for j in eachindex(out, num, cover, total)
        c = cover[j]
        if c <= 0 || c < t * total[j]
            out[j] = blank
        else
            out[j] = num[j] / c
        end
    end
    return out
end

# Cell areas

"""
    cellarea(space::RegridSpace, i::Integer) -> Float64

Return cell `i`'s area on the space manifold. The fallback builds the polygon;
spaces with a closed-form area should specialize it.
"""
cellarea(space::RegridSpace, i::Integer) =
    Float64(GO.area(manifold(space), getcell(space, Int(i))))

# N-D shaping

"""
    spatialdims(A) -> Tuple{Vararg{Int}}

Return the dimension numbers replaced by regridding. Dimensional arrays use
their X/Y dimensions; other arrays use the first two dimensions, or one for a
vector. The source cell count validates the result.
"""
spatialdims(A::AbstractArray) = ndims(A) <= 1 ? (1,) : (1, 2)

function spatialdims(A::DD.AbstractDimArray)
    ds = DD.dims(A)
    pos = Int[]
    for (i, d) in enumerate(ds)
        (d isa DD.Dimensions.XDim || d isa DD.Dimensions.YDim) && push!(pos, i)
    end
    isempty(pos) && return ndims(A) <= 1 ? (1,) : (1, 2)
    length(pos) == 1 && return (pos[1],)
    length(pos) == 2 && return (pos[1], pos[2])
    throw(ArgumentError("$(length(pos)) spatial dimensions found in $(DD.dims(A)); " *
                        "a regrid replaces at most two"))
end

# Resolve spatial dimensions against the source space's cell count.
function resolvespatialdims(data::AbstractArray, nsrc::Integer)
    sd = spatialdims(data)
    _spatialsize(data, sd) == nsrc && return _checkleading(data, sd)
    alt = length(sd) == 1 ? (1, 2) : (1,)
    ndims(data) >= length(alt) && _spatialsize(data, alt) == nsrc &&
        return _checkleading(data, alt)
    throw(DimensionMismatch(
        "source data of size $(size(data)) does not flatten to the source " *
        "space's $nsrc cells over dimensions $sd"))
end

_spatialsize(data, sd) = prod(ntuple(i -> size(data, sd[i]), length(sd)))

function _checkleading(data, sd)
    sd == ntuple(identity, length(sd)) || throw(ArgumentError(
        "spatial dimensions $sd of the source must come first; " *
        "`permutedims` the source so they do"))
    return sd
end

# Sizes of pass-through dimensions.
_otherdimsizes(data, sd) =
    ntuple(i -> size(data, length(sd) + i), ndims(data) - length(sd))

# Reshape the source to `ncells × nslices` in memory order.
function flatsource(data::AbstractArray, nsrc::Integer, nslices::Integer)
    raw = data isa DD.AbstractDimArray ? parent(data) : data
    isdiskbacked(raw) && (raw = Array(raw))
    return reshape(raw, Int(nsrc), Int(nslices))
end

# Output is floating point and preserves the source's ability to hold `missing`.
function outputeltype(::Type{Tin}) where {Tin}
    T = nonmissingtype(Tin)
    T <: Real || throw(ArgumentError("cannot regrid data of element type $Tin"))
    F = float(T)
    return Missing <: Tin ? Union{Missing,F} : F
end

"""
    wrapoutput(out, data, sd, dstdims) -> array

Put the destination's own axes `dstdims`, or one flat `Cell` axis when it has
none ([`destinationdims`](@ref)), before the source's unchanged non-spatial
dimensions. Sources that are not dimensional are returned unlabelled.
"""
function wrapoutput(out::AbstractArray, data, sd, dstdims)
    data isa DD.AbstractDimArray || return out
    ds = DD.dims(data)
    others = Tuple(ds[i] for i in eachindex(ds) if !(i in sd))
    dstdims === nothing &&
        return DD.DimArray(out, (DD.Dim{:Cell}(1:size(out, 1)), others...))
    shaped = reshape(out, map(length, dstdims)..., Base.tail(size(out))...)
    return DD.DimArray(shaped, (dstdims..., others...))
end

# Whole-domain apply

"""
    applyplan!(out, plan::DirectPlan, src) -> out

Apply the whole-domain plan to each source column. Accumulators and the reference
vector are reused across slices.
"""
function applyplan!(out::AbstractMatrix, plan::DirectPlan, src::AbstractMatrix)
    block = plan.block
    ndst, nsrc = size(block)
    size(src, 1) == nsrc || throw(DimensionMismatch(
        "plan expects $nsrc source cells, got $(size(src, 1))"))
    size(out, 1) == ndst || throw(DimensionMismatch(
        "plan produces $ndst destination cells, output holds $(size(out, 1))"))
    size(out, 2) == size(src, 2) || throw(DimensionMismatch(
        "$(size(src, 2)) source slices into $(size(out, 2)) output slices"))
    num = zeros(Float64, ndst)
    cover = zeros(Float64, ndst)
    ref = blockreference!(Vector{Float64}(undef, ndst), block)
    policy = plan.missingpolicy
    mv = plan.missingval
    for s in axes(src, 2)
        fill!(num, 0.0)
        fill!(cover, 0.0)
        x = view(src, :, s)
        applyblock!(num, cover, block, x, anyinvalid(x, mv) ? x : nothing, ref, mv)
        finalize!(view(out, :, s), num, cover, ref, policy)
    end
    return out
end
