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

# The blank an element type carries on its own, or `nothing` where it carries
# none: `Extensive` writes every cell, so a destination that could never be
# blanked is still a destination.
@inline _elementblank(::Type{T}) where {T} =
    Missing <: T ? missing : T <: AbstractFloat ? T(NaN) : nothing

# Use `missing` when supported and NaN otherwise.
@inline function _maskedvalue(::Type{T}) where {T}
    blank = _elementblank(T)
    blank === nothing && throw(ArgumentError(
        "`Weighted` blanks destinations below its coverage threshold and " *
        "needs a sentinel to do it, but $T holds neither `missing` nor NaN"))
    return blank
end

# The sentinel a destination of element type `T` is blanked with. A declared
# `missingval` must be a value `T` can hold, so a mismatch is the caller's to fix.
@inline _blankvalue(::Type{T}, ::Nothing) where {T} = _maskedvalue(T)
@inline _blankvalue(::Type{T}, ::Missing) where {T} =
    Missing <: T ? missing : throw(ArgumentError(
        "`missingval = missing` needs a destination that holds `missing`, " *
        "but this one holds $T"))
@inline _blankvalue(::Type{T}, missingval) where {T} = convert(T, missingval)

# Accumulation

"""
    addreference!(ref, block::WeightBlock) -> ref

Add `block`'s per-destination reference weight to `ref`, from the vector the
block stores. Blocks over one destination accumulate into one `ref`, which is
the reference the destination's values are normalized against.
"""
function addreference!(ref::AbstractVector{Float64}, block::WeightBlock)
    r = block.reference
    @inbounds for j in eachindex(ref, r)
        ref[j] += r[j]
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

Under a mask, coverage comes from the block's own coverage weights where it
declares them and from its value weights otherwise ([`WeightBlock`](@ref)): a
block whose values are signed cannot measure coverage with them.

`fallback` and `taint` are the accumulators [`degradetainted!`](@ref) reads, and
are filled only for a block that declares coverage, under a mask. Pass both or
neither.
"""
function applyblock!(num::AbstractVector{Float64}, cover::AbstractVector{Float64},
    block::WeightBlock, src::AbstractVector,
    valid::Union{Nothing,AbstractVector} = nothing,
    ref::Union{Nothing,AbstractVector{Float64}} = nothing,
    missingval = nothing;
    fallback::Union{Nothing,AbstractVector{Float64}} = nothing,
    taint::Union{Nothing,AbstractVector{Float64}} = nothing)
    W = block.weights
    Base.require_one_based_indexing(num, cover, src)
    size(W, 1) == length(num) == length(cover) || throw(DimensionMismatch(
        "block has $(size(W, 1)) destinations, accumulators have " *
        "$(length(num)) and $(length(cover))"))
    size(W, 2) == length(src) || throw(DimensionMismatch(
        "block expects $(size(W, 2)) source cells, got $(length(src))"))
    if valid === nothing
        _accumulate!(num, W, src, missingval)
        # A clean chunk still owes the first-order sum: another chunk's hole can
        # taint a destination this one also feeds, and the degrade then reads
        # `fallback` for the whole destination, not just the tainted chunk.
        if fallback !== nothing && block.coverage !== nothing
            _accumulate!(fallback, block.coverage, src, missingval)
        end
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
        C = block.coverage
        if C === nothing
            _accumulate!(num, cover, W, src, valid, missingval)
        else
            _accumulate!(num, W, src, missingval)
            _accumulatecoverage!(cover, C, valid, missingval)
            if fallback !== nothing && taint !== nothing
                Base.require_one_based_indexing(fallback, taint)
                length(fallback) == length(taint) == length(num) ||
                    throw(DimensionMismatch(
                        "degrade accumulators hold $(length(fallback)) and " *
                        "$(length(taint)) entries, destinations $(length(num))"))
                _accumulate!(fallback, C, src, missingval)
                _accumulatetainted!(taint, W, valid, missingval)
            end
        end
    end
    return num
end

# Very sparse blocks locate columns from nonzeros; denser blocks scan columns.
@inline _walknonzeros(W::SparseMatrixCSC) =
    SparseArrays.nnz(W) * 256 < size(W, 2)

# Locate the column owning nonzero index `p`.
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

# Accumulate coverage alone, from the block's separate coverage operator.
function _accumulatecoverage!(cover::AbstractVector{Float64}, C::SparseMatrixCSC,
    valid::AbstractVector, missingval = nothing)
    rows = SparseArrays.rowvals(C)
    vals = SparseArrays.nonzeros(C)
    cols = SparseArrays.getcolptr(C)
    if _walknonzeros(C)
        p = 1
        nz = SparseArrays.nnz(C)
        @inbounds while p <= nz
            k = _columnof(cols, p)
            p2 = cols[k+1] - 1
            m = _validity(valid[k], missingval)
            if !iszero(m)
                for q in p:p2
                    cover[rows[q]] += vals[q] * m
                end
            end
            p = p2 + 1
        end
        return cover
    end
    @inbounds for k in axes(C, 2)
        p1, p2 = cols[k], cols[k+1] - 1
        p1 > p2 && continue
        m = _validity(valid[k], missingval)
        iszero(m) && continue
        for p in p1:p2
            cover[rows[p]] += vals[p] * m
        end
    end
    return cover
end

function _accumulatecoverage!(cover::AbstractVector{Float64}, C::AbstractMatrix,
    valid::AbstractVector, missingval = nothing)
    @inbounds for k in axes(C, 2)
        m = _validity(valid[k], missingval)
        iszero(m) && continue
        for j in axes(C, 1)
            cover[j] += C[j, k] * m
        end
    end
    return cover
end

# Accumulate, per destination, the value weight that reached an invalid source.
# Nonzero exactly where a hole entered a signed sum — including through a
# gradient stencil, whose terms sit on a neighbour's column and so leave the
# coverage accounting untouched.
function _accumulatetainted!(taint::AbstractVector{Float64}, W::SparseMatrixCSC,
    valid::AbstractVector, missingval = nothing)
    rows = SparseArrays.rowvals(W)
    vals = SparseArrays.nonzeros(W)
    cols = SparseArrays.getcolptr(W)
    if _walknonzeros(W)
        p = 1
        nz = SparseArrays.nnz(W)
        @inbounds while p <= nz
            k = _columnof(cols, p)
            p2 = cols[k+1] - 1
            if iszero(_validity(valid[k], missingval))
                for q in p:p2
                    taint[rows[q]] += abs(vals[q])
                end
            end
            p = p2 + 1
        end
        return taint
    end
    @inbounds for k in axes(W, 2)
        p1, p2 = cols[k], cols[k+1] - 1
        p1 > p2 && continue
        iszero(_validity(valid[k], missingval)) || continue
        for p in p1:p2
            taint[rows[p]] += abs(vals[p])
        end
    end
    return taint
end

function _accumulatetainted!(taint::AbstractVector{Float64}, W::AbstractMatrix,
    valid::AbstractVector, missingval = nothing)
    @inbounds for k in axes(W, 2)
        iszero(_validity(valid[k], missingval)) || continue
        for j in axes(W, 1)
            taint[j] += abs(W[j, k])
        end
    end
    return taint
end

# Finalization

"""
    degradetainted!(num, fallback, taint) -> num

Take the first-order sum wherever a signed weight reached an invalid source.

`num` is the sum a signed block accumulated, `fallback` the same sum over the
block's coverage weights alone, and `taint` the weight that reached a hole
([`applyblock!`](@ref)). A destination with `taint > 0` takes `fallback`, which
is what [`Conservative`](@ref) would have given it.

The correction a signed method adds is fitted from a source cell's neighbours,
so one hole biases every neighbour it has by a share of the *field's* value, not
of its gradient — and, because a stencil term sits on a neighbour's column, that
bias never reaches `cover` for [`Weighted`](@ref) to threshold against. Dropping
the correction where it was fitted from a hole is the conservative reading: a
destination away from every hole keeps the better answer.

Accuracy near a hole, not conservation, is what this trades. The correction sums
to zero over the destinations covering one source cell, so replacing it for some
of them and not others leaves an [`Extensive`](@ref) total short by the
difference. Build the plan against a known mask to keep both.
"""
function degradetainted!(num::AbstractVector{Float64},
    fallback::AbstractVector{Float64}, taint::AbstractVector{Float64})
    @inbounds for j in eachindex(num, fallback, taint)
        taint[j] > 0 && (num[j] = fallback[j])
    end
    return num
end

"""
    finalize!(out, num, cover, total, policy, missingval = _maskedvalue(eltype(out))) -> out

Finalize one destination chunk after all source blocks have accumulated.

`num` is the weighted sum with invalid sources contributing zero, `cover` the
weight of the sources that were valid — the applied blocks' coverage weights, or
their value weights where they report none — and `total` the reference weight
the applied blocks accumulated ([`addreference!`](@ref)) — the denominator the
method reported, or the row sums when it reported none.

  - [`Extensive`](@ref): `num`.
  - [`Weighted`](@ref)`(t)`: `num / cover`, blank where `cover ≤ 0` or
    `cover < t · total`.

`Weighted` always normalizes by valid coverage and applies its threshold against
`total`. `missingval` is the sentinel blanked cells receive, and defaults to
`missing` when `out` holds it and NaN otherwise ([`outputmissingval`](@ref)).
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

finalize!(out::AbstractVector, num::AbstractVector{Float64},
    cover::AbstractVector{Float64}, total::AbstractVector{Float64},
    policy::AbstractMissingPolicy) =
    finalize!(out, num, cover, total, policy, _maskedvalue(eltype(out)))

function finalize!(out::AbstractVector, num::AbstractVector{Float64},
    cover::AbstractVector{Float64}, total::AbstractVector{Float64},
    ::Extensive, missingval)
    @inbounds for j in eachindex(out, num)
        out[j] = num[j]
    end
    return out
end

function finalize!(out::AbstractVector, num::AbstractVector{Float64},
    cover::AbstractVector{Float64}, total::AbstractVector{Float64},
    policy::Weighted, missingval)
    blank = _blankvalue(eltype(out), missingval)
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
    declareschunks(raw) && (raw = Array(raw))
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
    outputeltype(Tin, missingval) -> Type

The floating point element type a regrid of a `Tin` source writes, wide enough
to hold `missingval`. `missing` unions itself in; every other sentinel is a
number the float type already holds, so the result stays concrete.
"""
outputeltype(::Type{Tin}, missingval) where {Tin} =
    nonmissingtype(outputeltype(Tin))

outputeltype(::Type{Tin}, ::Missing) where {Tin} =
    Union{Missing,nonmissingtype(outputeltype(Tin))}

outputeltype(::Type{Tin}, ::Nothing) where {Tin} = outputeltype(Tin)

"""
    outputmissingval(data) -> value

The sentinel a regrid of `data` blanks uncovered destination cells with, and
the `missingval` a labelled result declares. Arrays take `missing` when their
element type holds it and NaN otherwise; a `Rasters.AbstractRaster` inherits
its own `missingval`, so a regrid keeps the nodata convention it was handed.

Pass `missingval` to [`regrid`](@ref) to override it.
"""
outputmissingval(data) = _maskedvalue(outputeltype(eltype(data)))

"""
    destinationmissingval(dest) -> value

The sentinel [`regrid!`](@ref) blanks uncovered cells of `dest` with. A
preallocated destination declares its own nodata convention, so this reads
`dest` where [`outputmissingval`](@ref) reads the source.
"""
destinationmissingval(dest) = _elementblank(eltype(dest))

"""
    wrapoutput(out, data, sd, dstdims, missingval = _maskedvalue(eltype(out))) -> array

Put the destination's own axes `dstdims`, or one flat `Cell` axis when it has
none ([`destinationdims`](@ref)), before the source's unchanged non-spatial
dimensions. Sources that are not dimensional are returned unlabelled.

A destination naming one axis is already the shape the cells were written in,
so the result is labelled over the written array itself; two or more axes split
the leading cell axis by reshaping it. [`rebuildoutput`](@ref) decides what the
labelled result is, and `missingval` is the sentinel it declares.
"""
function wrapoutput(out::AbstractArray, data, sd, dstdims,
    missingval = _maskedvalue(eltype(out)))
    data isa DD.AbstractDimArray || return out
    shaped, newdims = _outputlabels(out, data, sd, dstdims)
    # Declare the value actually written, so a sentinel given in the source's
    # own type reaches the result in the destination's.
    return rebuildoutput(data, shaped, newdims, _blankvalue(eltype(out), missingval))
end

# The destination's axes and the source's pass-through dimensions, over an array
# reshaped to match them.
function _outputlabels(out::AbstractArray, data, sd, dstdims)
    ds = DD.dims(data)
    others = Tuple(ds[i] for i in eachindex(ds) if !(i in sd))
    dstdims === nothing &&
        return out, (DD.Dim{:Cell}(1:size(out, 1)), others...)
    length(dstdims) == 1 && return out, (dstdims..., others...)
    shaped = reshape(out, map(length, dstdims)..., Base.tail(size(out))...)
    return shaped, (dstdims..., others...)
end

"""
    rebuildoutput(data, out, dims, missingval) -> array

Label `out` with `dims` as the result of regridding `data`. A dimensional
source that carries nodata of its own — a `Rasters.AbstractRaster` — rebuilds
into its own array type declaring `missingval`; everything else becomes a plain
`DimArray`, which has no nodata to declare.

Both the eager and the lazy path label through here, so a reader and a
`LazyRegridArray` come back wrapped the same way.
"""
rebuildoutput(data, out::AbstractArray, dims::Tuple, missingval) =
    DD.DimArray(out, dims)

# Whole-domain apply

"""
    applyplan!(out, plan::DirectPlan, src, missingval = _maskedvalue(eltype(out))) -> out

Apply the whole-domain plan to each source column, blanking uncovered cells with
`missingval`. Accumulators are reused across slices, and the reference weights
are the block's own, so a second application of one plan allocates nothing beyond
its accumulators.
"""
function applyplan!(out::AbstractMatrix, plan::DirectPlan, src::AbstractMatrix,
    missingval = _maskedvalue(eltype(out)))
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
    # A block that declares coverage has signed weights, so a hole in the source
    # can reach a destination without reaching its coverage.
    signed = block.coverage !== nothing
    fallback = signed ? zeros(Float64, ndst) : nothing
    taint = signed ? zeros(Float64, ndst) : nothing
    # The block's own reference, not a copy: applying one plan again allocates
    # no second reference vector.
    ref = block.reference
    policy = plan.missingpolicy
    mv = plan.missingval
    for s in axes(src, 2)
        fill!(num, 0.0)
        fill!(cover, 0.0)
        signed && (fill!(fallback, 0.0); fill!(taint, 0.0))
        x = view(src, :, s)
        applyblock!(num, cover, block, x, anyinvalid(x, mv) ? x : nothing, ref, mv;
            fallback, taint)
        signed && degradetainted!(num, fallback, taint)
        finalize!(view(out, :, s), num, cover, ref, policy, missingval)
    end
    return out
end
