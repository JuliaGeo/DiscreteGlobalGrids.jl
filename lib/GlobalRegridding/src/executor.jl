# The apply core: sparse accumulation and missing-policy finalize, plus the N-D
# shaping around them. A block, a loaded source vector, three accumulators and a
# policy — nothing here knows about IO, chunk discovery, or which plan it serves,
# which is why the eager path and the lazy one call the same functions.

# ===========================================================================
# Validity
# ===========================================================================

# A value is invalid when it is `missing` or a NaN; anything else counts, so an
# integer or boolean source can never be invalid.
@inline _isvalid(::Missing) = false
@inline _isvalid(x::AbstractFloat) = !isnan(x)
@inline _isvalid(::Any) = true

# The contribution an entry makes to a weighted sum: itself, or zero when
# invalid.
@inline _value(::Missing) = 0.0
@inline _value(x::AbstractFloat) = isnan(x) ? 0.0 : Float64(x)
@inline _value(x::Real) = Float64(x)

# The contribution an entry makes to coverage: one when valid, zero when not.
# A `Bool` is read as the flag it is, so an explicit validity mask works as well
# as the data itself.
@inline _validity(::Missing) = 0.0
@inline _validity(x::Bool) = x ? 1.0 : 0.0
@inline _validity(x::AbstractFloat) = isnan(x) ? 0.0 : 1.0
@inline _validity(::Real) = 1.0

# Whether an element type can hold an invalid value at all. `false` skips the
# scan and the coverage matvec outright.
_canbeinvalid(::Type{T}) where {T} = !(T <: Real) || T <: AbstractFloat || Missing <: T

"""
    anyinvalid(src) -> Bool

Whether `src` holds any `missing` or NaN. `false` for an element type that
cannot hold one, without reading `src`.
"""
anyinvalid(src::AbstractArray) =
    _canbeinvalid(eltype(src)) ? any(!_isvalid, src) : false

# The sentinel a masked destination is written as: `missing` when the output
# element type admits it, NaN otherwise.
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

# ===========================================================================
# Accumulation
# ===========================================================================

"""
    blockreference!(ref, block::WeightBlock) -> ref

Overwrite `ref` with `block`'s per-destination reference weight: its stored
denominator, or its row sums when it carries none.

This is what a destination's accumulated valid coverage is compared against, and
it is data-independent — computed once per block and reused across every slice
and every read the block serves.
"""
function blockreference!(ref::AbstractVector{Float64}, block::WeightBlock)
    fill!(ref, 0.0)
    return addreference!(ref, block)
end

"""
    addreference!(ref, block::WeightBlock) -> ref

Add `block`'s per-destination reference weight to `ref`. Blocks over one
destination chunk sum, so accumulating over the connected source chunks gives
the destination's whole reference.
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
    applyblock!(num, cover, block::WeightBlock, src, valid = nothing, ref = nothing)

Accumulate one source chunk's contribution into the destination accumulators.

`num[j] += Σₖ W[j,k] · src[k]`, with invalid `src` entries contributing zero, and
`cover[j] += Σₖ W[j,k] · valid[k]`, the weight of the source cells that were not
missing. Both accumulate: applying every block that meets a destination chunk
sums to the whole, which is what lets any linear method stream.

`src` is one source chunk's cells, flattened in `cellindices` order; its length
must be the block's column count.

`valid` is the validity mask. `nothing` asserts that every entry of `src` is
valid, and takes the shortcut of crediting the block's full reference weight
instead of a second matvec. Otherwise it is any vector the same length as `src`
read entrywise as valid-or-not — `missing` and NaN are invalid, a `Bool` is
itself, anything else is valid — so passing `src` itself derives the mask from
the data without materializing one.

`ref` is `block`'s reference weight from [`blockreference!`](@ref), used only on
the `valid === nothing` path; passing it avoids recomputing row sums per slice.

Allocates nothing.
"""
function applyblock!(num::AbstractVector{Float64}, cover::AbstractVector{Float64},
    block::WeightBlock, src::AbstractVector,
    valid::Union{Nothing,AbstractVector} = nothing,
    ref::Union{Nothing,AbstractVector{Float64}} = nothing)
    W = block.weights
    Base.require_one_based_indexing(num, cover, src)
    size(W, 1) == length(num) == length(cover) || throw(DimensionMismatch(
        "block has $(size(W, 1)) destinations, accumulators have " *
        "$(length(num)) and $(length(cover))"))
    size(W, 2) == length(src) || throw(DimensionMismatch(
        "block expects $(size(W, 2)) source cells, got $(length(src))"))
    if valid === nothing
        _accumulate!(num, W, src)
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
        _accumulate!(num, cover, W, src, valid)
    end
    return num
end

# Values only — every source entry is valid.
function _accumulate!(num::AbstractVector{Float64}, W::SparseMatrixCSC, src::AbstractVector)
    rows = SparseArrays.rowvals(W)
    vals = SparseArrays.nonzeros(W)
    @inbounds for k in axes(W, 2)
        x = _value(src[k])
        iszero(x) && continue
        for p in SparseArrays.nzrange(W, k)
            num[rows[p]] += vals[p] * x
        end
    end
    return num
end

function _accumulate!(num::AbstractVector{Float64}, W::AbstractMatrix, src::AbstractVector)
    @inbounds for k in axes(W, 2)
        x = _value(src[k])
        iszero(x) && continue
        for j in axes(W, 1)
            num[j] += W[j, k] * x
        end
    end
    return num
end

# Values and coverage in one pass over the nonzeros.
function _accumulate!(num::AbstractVector{Float64}, cover::AbstractVector{Float64},
    W::SparseMatrixCSC, src::AbstractVector, valid::AbstractVector)
    rows = SparseArrays.rowvals(W)
    vals = SparseArrays.nonzeros(W)
    @inbounds for k in axes(W, 2)
        x = _value(src[k])
        m = _validity(valid[k])
        (iszero(x) && iszero(m)) && continue
        for p in SparseArrays.nzrange(W, k)
            j = rows[p]
            v = vals[p]
            num[j] += v * x
            cover[j] += v * m
        end
    end
    return num
end

function _accumulate!(num::AbstractVector{Float64}, cover::AbstractVector{Float64},
    W::AbstractMatrix, src::AbstractVector, valid::AbstractVector)
    @inbounds for k in axes(W, 2)
        x = _value(src[k])
        m = _validity(valid[k])
        (iszero(x) && iszero(m)) && continue
        for j in axes(W, 1)
            v = W[j, k]
            num[j] += v * x
            cover[j] += v * m
        end
    end
    return num
end

# ===========================================================================
# Finalize — the missing-data semantics
# ===========================================================================

"""
    finalize!(out, num, cover, total, policy, hasdenom::Bool) -> out
    finalize!(out, num, cover, total, block::WeightBlock, policy) -> out

Turn the accumulators of one destination chunk into values, applying the missing
policy. Called once per destination, after every connected source chunk has been
accumulated.

`num` is the weighted sum with invalid sources contributing zero, `cover` the
weight of the sources that were valid, and `total` the accumulated reference
weight from [`blockreference!`](@ref) — the denominator the method reported, or
the row sums when it reported none. `hasdenom` says which of those two `total`
is, and is the whole of the distinction below; the second form reads it off a
block.

| block | [`Extensive`](@ref) | [`Weighted`](@ref)`(t)` |
|---|---|---|
| with a denominator | `num` | `num / cover`, blank where `cover ≤ 0` or `cover < t · total` |
| without one | `num` | `num`, blank where `cover ≤ 0` or `cover < t · total` |

A method's denominator is its declaration of what full coverage weighs, and only
a block that declares one is normalized: weights that already sum to one are
returned as they are rather than silently renormalized over the part of a
stencil that survived. A point method that wants a truncated stencil rescaled
rather than blanked reports a denominator ([`adddenom!`](@ref)) and gets the
first row.

Blanked destinations are written as `missing` when `eltype(out)` admits it and
NaN otherwise. `Extensive` blanks nothing: a destination no valid source reached
is a raw sum of zero, which is the lost mass rather than an absence of one.
"""
function finalize! end

finalize!(out::AbstractVector, num::AbstractVector{Float64},
    cover::AbstractVector{Float64}, total::AbstractVector{Float64},
    block::WeightBlock, policy::AbstractMissingPolicy) =
    finalize!(out, num, cover, total, policy, block.denom !== nothing)

function finalize!(out::AbstractVector, num::AbstractVector{Float64},
    cover::AbstractVector{Float64}, total::AbstractVector{Float64},
    ::Extensive, ::Bool)
    @inbounds for j in eachindex(out, num)
        out[j] = num[j]
    end
    return out
end

function finalize!(out::AbstractVector, num::AbstractVector{Float64},
    cover::AbstractVector{Float64}, total::AbstractVector{Float64},
    policy::Weighted, hasdenom::Bool)
    blank = _maskedvalue(eltype(out))
    t = policy.threshold
    @inbounds for j in eachindex(out, num, cover, total)
        c = cover[j]
        if c <= 0 || c < t * total[j]
            out[j] = blank
        else
            out[j] = hasdenom ? num[j] / c : num[j]
        end
    end
    return out
end

# ===========================================================================
# Cell areas
# ===========================================================================

"""
    cellarea(space::RegridSpace, i::Integer) -> Float64

The area of cell `i` of `space` on its manifold.

Generic and correct for any space, at the cost of building the cell polygon; a
space with a closed-form area should say so by defining this.

Nothing in the executor calls it: a destination's coverage threshold is measured
against the denominator its method reported, not against the destination cell's
own area. It exists for methods and callers that need an area the contract does
not otherwise expose.
"""
cellarea(space::RegridSpace, i::Integer) =
    Float64(GO.area(manifold(space), getcell(space, Int(i))))

# ===========================================================================
# N-D shaping
# ===========================================================================

"""
    spatialdims(A) -> Tuple{Vararg{Int}}

The dimensions of `A` that the regrid replaces, as dimension numbers.

A `DimensionalData` array answers with the positions of its `XDim`/`YDim`
dimensions; anything else answers `(1, 2)`, or `(1,)` when it is a vector. Every
other dimension passes through a regrid untouched.

The answer is a shape claim, not a certainty — [`regrid`](@ref) checks it
against the source space's cell count and falls back to the other reading before
erroring.
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

# The spatial dimensions the source data is actually flattened over: the trait's
# answer when it matches the space, the other reading when that one does.
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

# The sizes of the dimensions a regrid passes through, in order.
_otherdimsizes(data, sd) =
    ntuple(i -> size(data, length(sd) + i), ndims(data) - length(sd))

# The source as an `ncells × nslices` matrix: `vec` over the spatial dimensions
# in memory order, which is the order the space numbers its cells in.
function flatsource(data::AbstractArray, nsrc::Integer, nslices::Integer)
    raw = data isa DD.AbstractDimArray ? parent(data) : data
    isdiskbacked(raw) && (raw = Array(raw))
    return reshape(raw, Int(nsrc), Int(nslices))
end

# The element type a regrid of `Tin` produces: floating point, since weights are,
# and carrying `missing` exactly when the source does.
function outputeltype(::Type{Tin}) where {Tin}
    T = nonmissingtype(Tin)
    T <: Real || throw(ArgumentError("cannot regrid data of element type $Tin"))
    F = float(T)
    return Missing <: Tin ? Union{Missing,F} : F
end

# Destination dimensions: the destination's own cell axis, then the source's
# non-spatial dimensions unchanged. A plain array in, a plain array out.
function wrapoutput(out::AbstractArray, data, sd)
    data isa DD.AbstractDimArray || return out
    ds = DD.dims(data)
    others = Tuple(ds[i] for i in eachindex(ds) if !(i in sd))
    return DD.DimArray(out, (DD.Dim{:Cell}(1:size(out, 1)), others...))
end

# ===========================================================================
# Whole-domain apply
# ===========================================================================

"""
    applyplan!(out, plan::DirectPlan, src) -> out

Apply `plan`'s whole-domain block to every column of the `ncells × nslices`
source matrix `src`, writing the `ndst × nslices` result into `out`.

One plan, one set of accumulators, one reference vector, reused across the
slices: an N-D regrid costs one weight construction and one reference pass, not
one per field.
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
    hasdenom = block.denom !== nothing
    policy = plan.missingpolicy
    for s in axes(src, 2)
        fill!(num, 0.0)
        fill!(cover, 0.0)
        x = view(src, :, s)
        applyblock!(num, cover, block, x, anyinvalid(x) ? x : nothing, ref)
        finalize!(view(out, :, s), num, cover, ref, policy, hasdenom)
    end
    return out
end
