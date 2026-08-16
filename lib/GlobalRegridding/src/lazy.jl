# `LazyRegridArray`, the lazy destination array. Owned by task T7.
#
# The whole lazy path is one function — `readblock!` — and it is generic over
# methods, spaces and storage: requested range → covering destination chunks →
# connected source chunks (dilated discovery) → per pair, a block from the plan's
# storage → load that one source chunk → accumulate → finalize once per
# destination chunk. Nothing here knows what a method is; the blocks column-
# partition the operator by source chunk, and that is what makes every linear
# method stream without a line of per-method code.

# ===========================================================================
# Source addressing
# ===========================================================================

"""
    chunkranges(space::RegridSpace, chunk, spatialsize::NTuple{NS,Int})
        -> NTuple{NS,UnitRange{Int}}

The index ranges of `chunk` in the **array's** spatial dimensions, in array
dimension order, for an array whose spatial dimensions have size `spatialsize`.

This is what turns a chunk number into one contiguous read. Its contract is that
`vec(A[chunkranges(space, chunk, size)...])` is the chunk's cells in
`cellindices` order — the same flattening the executor uses, column-major over
the array's spatial dimensions.

The generic implementation recovers the ranges from [`cellindices`](@ref) and
verifies that the chunk really is a rectangle of the lattice; a space that knows
its own rectangle should say so by defining this, and one whose chunks are not
lattice rectangles must.
"""
function chunkranges end

function chunkranges(space::RegridSpace, chunk::Integer, ::NTuple{1,Int})
    inds = cellindices(space, Int(chunk))
    n = length(inds)
    n == 0 && return (1:0,)
    lo, hi = Int(first(inds)), Int(last(inds))
    hi - lo + 1 == n || _notrectangular(space, chunk)
    return (lo:hi,)
end

function chunkranges(space::RegridSpace, chunk::Integer, ssz::NTuple{2,Int})
    inds = cellindices(space, Int(chunk))
    n = length(inds)
    n == 0 && return (1:0, 1:0)
    n1 = ssz[1]
    c1, r1 = fldmod1(Int(first(inds)), n1)
    c2, r2 = fldmod1(Int(last(inds)), n1)
    rows, cols = r1:r2, c1:c2
    length(rows) * length(cols) == n || _notrectangular(space, chunk)
    k = 0
    for c in cols, r in rows
        k += 1
        @inbounds Int(inds[k]) == r + (c - 1) * n1 || _notrectangular(space, chunk)
    end
    return (rows, cols)
end

chunkranges(space::RegridSpace, chunk::Integer, ssz::Tuple) = throw(ArgumentError(
    "a regrid flattens over one or two spatial dimensions, got $(length(ssz))"))

_notrectangular(space, chunk) = throw(ArgumentError(
    "chunk $chunk of $(typeof(space)) is not a rectangle of the array lattice, " *
    "so it cannot be read in one block; define " *
    "`GlobalRegridding.chunkranges(::$(typeof(space)), chunk, spatialsize)`"))

# `RasterGrid` holds its chunk rectangle outright; `chunkbox` is (X, Y) and the
# array's spatial dimensions are in `xfast` order.
chunkranges(space::RasterGrid, chunk::Integer, ::NTuple{2,Int}) =
    space.xfast ? chunkbox(space, Int(chunk)) : reverse(chunkbox(space, Int(chunk)))

# ===========================================================================
# The array
# ===========================================================================

"""
    LazyRegridArray(data, plan::ChunkedPlan)

The regrid of `data` under `plan`, as a chunked `DiskArrays.AbstractDiskArray`
that computes each block on demand.

Dimensions are `(destination cells, data's non-spatial dims...)`, and chunks are
the destination space's chunks against the source's own chunking of the
pass-through dimensions — so one chunk of this array is one destination chunk of
one slice group, which is exactly the unit `readblock!` computes.

**Constructing one reads no source data.** Weights are built on first touch and
kept in the plan's storage, so the second read of a destination chunk, and every
further non-spatial slice, builds nothing.

Reading a block loads only the source chunks discovered to reach it, each once,
and holds only one of them at a time: accumulators are sized to the destination
**chunk**, never to the destination space.
"""
struct LazyRegridArray{T,N,NS,A,P<:ChunkedPlan,K,C} <: DiskArrays.AbstractDiskArray{T,N}
    source::A
    plan::P
    srcsize::NTuple{NS,Int}
    size::NTuple{N,Int}
    srctree::K
    dstcaps::Vector{Cap}
    dstspans::Vector{UnitRange{Int}}
    radius::Float64
    chunks::C
end

function LazyRegridArray(data, plan::ChunkedPlan)
    src_space, dst_space = plan.src_space, plan.dst_space
    nsrc = Int(ncells(src_space))
    ndst = Int(ncells(dst_space))
    sd = resolvespatialdims(data, nsrc)
    othersizes = _otherdimsizes(data, sd)
    source = data isa DD.AbstractDimArray ? parent(data) : data
    nspatial = length(sd)
    srcsize = ntuple(i -> size(source, i), nspatial)
    spans, contiguous = _chunkspans(dst_space)
    caps = chunkextents(dst_space)
    radius = Float64(support_radius(plan.method, src_space))
    radius >= 0 || throw(ArgumentError(
        "support_radius must be a non-negative angular radius, got $radius"))
    chunks = _outputchunks(plan, source, ndst, spans, contiguous, nspatial, othersizes)
    # The source's chunk tree is descended once per destination chunk per read;
    # building it is O(nchunks) for a flat tree, so it is built once here.
    srctree = chunktree(src_space)
    T = outputeltype(eltype(data))
    return LazyRegridArray{T,length(othersizes) + 1,nspatial,typeof(source),
        typeof(plan),typeof(srctree),typeof(chunks)}(
        source, plan, srcsize, (ndst, othersizes...), srctree, caps, spans, radius, chunks)
end

Base.size(A::LazyRegridArray) = A.size

Base.show(io::IO, ::MIME"text/plain", A::LazyRegridArray) = show(io, A)
Base.show(io::IO, A::LazyRegridArray{T}) where {T} =
    print(io, "LazyRegridArray{", T, "}(", join(size(A), "×"), ", ",
        typeof(A.plan.method).name.name, ")")

DiskArrays.haschunks(::LazyRegridArray) = DiskArrays.Chunked()
DiskArrays.eachchunk(A::LazyRegridArray) = A.chunks

# The span of cell positions each destination chunk owns, plus whether every
# chunk is contiguous — which decides whether the cell axis can report the
# destination's chunks as its own chunking.
function _chunkspans(space::RegridSpace)
    nc = Int(nchunks(space))
    spans = Vector{UnitRange{Int}}(undef, nc)
    contiguous = true
    for c in 1:nc
        inds = cellindices(space, c)
        if isempty(inds)
            spans[c] = 1:0
            contiguous = false
        else
            spans[c] = Int(first(inds)):Int(last(inds))
            contiguous &= length(spans[c]) == length(inds)
        end
    end
    return spans, contiguous
end

function _outputchunks(plan::ChunkedPlan, source, ndst::Int,
    spans::Vector{UnitRange{Int}}, contiguous::Bool, nspatial::Int, othersizes::Tuple)
    declared = plan.chunks
    nd = length(othersizes) + 1
    if declared isa DiskArrays.GridChunks
        length(declared.chunks) == nd || throw(ArgumentError(
            "the plan's chunking has $(length(declared.chunks)) dimensions, but " *
            "the regrid of this source has $nd"))
        return declared
    elseif declared isa Tuple{Vararg{Integer}}
        length(declared) == nd || throw(ArgumentError(
            "the plan's chunk sizes $declared do not cover the regrid's $nd dimensions"))
        return DiskArrays.GridChunks((ndst, othersizes...), declared)
    end
    return DiskArrays.GridChunks(_cellchunks(spans, contiguous, ndst),
        _passthroughchunks(source, nspatial, othersizes)...)
end

# The cell axis reports the destination's own chunks when they partition it into
# ascending contiguous runs, and one whole chunk otherwise — a chunking is a
# read-granularity hint, and `readblock!` covers any request either way.
function _cellchunks(spans::Vector{UnitRange{Int}}, contiguous::Bool, ndst::Int)
    whole = DiskArrays.RegularChunks(max(ndst, 1), 0, ndst)
    contiguous || return whole
    sizes = Vector{Int}(undef, length(spans))
    next = 1
    for (c, sp) in enumerate(spans)
        first(sp) == next || return whole
        sizes[c] = length(sp)
        next = last(sp) + 1
    end
    next == ndst + 1 || return whole
    return DiskArrays.IrregularChunks(; chunksizes = sizes)
end

# Pass-through dimensions keep the source's chunking, so a slice group is read
# whole rather than in pieces.
function _passthroughchunks(source, nspatial::Int, othersizes::Tuple)
    if DiskArrays.haschunks(source) isa DiskArrays.Chunked
        ec = DiskArrays.eachchunk(source)
        if ec isa DiskArrays.GridChunks && length(ec.chunks) == nspatial + length(othersizes)
            return ntuple(i -> ec.chunks[nspatial+i], length(othersizes))
        end
    end
    return ntuple(i -> DiskArrays.RegularChunks(max(othersizes[i], 1), 0, othersizes[i]),
        length(othersizes))
end

# ===========================================================================
# readblock!
# ===========================================================================

function DiskArrays.readblock!(A::LazyRegridArray, aout::AbstractArray,
    r::AbstractUnitRange...)
    length(r) == ndims(A) || throw(DimensionMismatch(
        "$(length(r)) ranges requested from a $(ndims(A))-dimensional regrid"))
    cellr = r[1]
    others = Base.tail(r)
    nslices = prod(map(length, others); init = 1)
    return _readdestination!(reshape(aout, length(cellr), nslices), A, cellr, others, nslices)
end

"""
    _readdestination!(out, A, cellr, others, nslices)

The lazy path's one loop, into an `ncells × nslices` view of the output block.

Residency at any moment: one source chunk over the requested slices, one
`WeightBlock` (plus whatever the storage's bound keeps), and the destination
**chunk's** accumulators — never the destination space's, and never the source
array's.
"""
function _readdestination!(out::AbstractMatrix, A::LazyRegridArray{T,N,NS},
    cellr::AbstractUnitRange, others::Tuple, nslices::Int) where {T,N,NS}
    plan = A.plan
    dst_space, src_space = plan.dst_space, plan.src_space
    policy = plan.missingpolicy
    olens = map(length, others)
    srcbuf = _emptybuffer(A)
    srcchunks = Int[]
    num = Matrix{Float64}(undef, 0, nslices)
    cover = Matrix{Float64}(undef, 0, nslices)
    total = Float64[]
    vals = T[]
    for c in _coveringchunks(A, cellr)
        dinds = cellindices(dst_space, c)
        nd = length(dinds)
        num = _fitmatrix(num, nd, nslices)
        cover = _fitmatrix(cover, nd, nslices)
        length(total) == nd || (total = Vector{Float64}(undef, nd))
        length(vals) == nd || (vals = Vector{T}(undef, nd))
        fill!(num, 0.0)
        fill!(cover, 0.0)
        fill!(total, 0.0)
        connectedchunks!(srcchunks, A.dstcaps[c], A.srctree; radius = A.radius)
        denominated = false
        for s in srcchunks
            entry = blockfor(plan, c, s)
            addreference!(total, entry.block)
            denominated |= hasdenom(entry.block)
            sr = chunkranges(src_space, s, A.srcsize)
            slens = map(length, sr)
            srcbuf = _fitbuffer(srcbuf, slens, olens)
            _readsource!(srcbuf, A.source, sr, others)
            _applyslices!(num, cover, entry.block,
                reshape(srcbuf, prod(slens), nslices), entry.ref, nslices)
        end
        _writechunk!(out, vals, num, cover, total, policy, denominated, dinds, cellr)
    end
    return out
end

# The destination chunks whose cells can fall in `cellr`. Spans bound a chunk's
# positions, so this over-reports only for a space whose chunks interleave, and
# never misses one.
function _coveringchunks(A::LazyRegridArray, cellr::AbstractUnitRange)
    out = Int[]
    lo, hi = first(cellr), last(cellr)
    @inbounds for c in eachindex(A.dstspans)
        sp = A.dstspans[c]
        (isempty(sp) || last(sp) < lo || first(sp) > hi) && continue
        push!(out, c)
    end
    return out
end

# One block against every requested slice of one loaded source chunk. This is
# the function barrier: the block's concrete type is recovered here, once per
# chunk pair, so the accumulation kernels below it are statically dispatched.
function _applyslices!(num::Matrix{Float64}, cover::Matrix{Float64},
    block::WeightBlock, src::AbstractMatrix, ref::Vector{Float64}, nslices::Int)
    dirty = anyinvalid(src)
    for t in 1:nslices
        col = view(src, :, t)
        applyblock!(view(num, :, t), view(cover, :, t), block, col,
            dirty ? col : nothing, ref)
    end
    return num
end

# Finalize once per destination chunk per slice, then scatter the cells the
# request actually asked for.
function _writechunk!(out::AbstractMatrix, vals::Vector, num::Matrix{Float64},
    cover::Matrix{Float64}, total::Vector{Float64}, policy::AbstractMissingPolicy,
    denominated::Bool, dinds, cellr::AbstractUnitRange)
    lo, hi = first(cellr), last(cellr)
    off = lo - 1
    for t in axes(num, 2)
        finalize!(vals, view(num, :, t), view(cover, :, t), total, policy, denominated)
        @inbounds for (j, p) in enumerate(dinds)
            lo <= p <= hi || continue
            out[p-off, t] = vals[j]
        end
    end
    return out
end

# ===========================================================================
# Source loading
# ===========================================================================

_emptybuffer(::LazyRegridArray{T,N,NS,A}) where {T,N,NS,A} =
    Array{eltype(A),NS + N - 1}(undef, ntuple(_ -> 0, NS + N - 1))

_fitmatrix(M::Matrix{Float64}, n::Int, m::Int) =
    size(M) == (n, m) ? M : Matrix{Float64}(undef, n, m)

function _fitbuffer(buf::Array{S,M}, slens::Tuple, olens::Tuple) where {S,M}
    shape = (slens..., olens...)::NTuple{M,Int}
    return size(buf) == shape ? buf : Array{S,M}(undef, shape)
end

# One chunk-aligned read into the reused buffer. A disk-backed source is asked
# through `readblock!`, which is one read per source chunk by construction.
function _readsource!(buf::Array, source, sr::Tuple, others::Tuple)
    if _isdisksource(source)
        DiskArrays.readblock!(source, buf, sr..., others...)
    else
        copyto!(buf, view(source, sr..., others...))
    end
    return buf
end

_isdisksource(x) = x isa DiskArrays.AbstractDiskArray || DiskArrays.isdisk(x)

# ===========================================================================
# API
# ===========================================================================

"""
    regrid(data, plan::ChunkedPlan) -> LazyRegridArray

The lazy application of a chunked plan: a disk array that computes destination
chunks on demand and reads no source data until one is asked for.
"""
regrid(data, plan::ChunkedPlan) = LazyRegridArray(data, plan)

"""
    regrid!(dest, data, plan::ChunkedPlan) -> dest

Materialize a chunked plan into `dest`, one destination chunk at a time.

Each chunk's source chunks are loaded, applied and dropped before the next
chunk's, so peak residency is a chunk's worth of source and not the source
array.
"""
function regrid!(dest, data, plan::ChunkedPlan)
    A = LazyRegridArray(data, plan)
    size(dest) == size(A) || throw(DimensionMismatch(
        "destination of size $(size(dest)) cannot hold a regrid of size $(size(A))"))
    raw = dest isa DD.AbstractDimArray ? parent(dest) : dest
    for ranges in DiskArrays.eachchunk(A)
        raw[ranges...] = A[ranges...]
    end
    return dest
end
