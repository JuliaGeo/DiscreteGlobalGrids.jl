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
# Destination tiling
# ===========================================================================

# The destination side of a read costs, per cell and per slice, four Float64
# accumulators and one output value; an eighth of the budget is what a derived
# tile is allowed to spend on them.
const DEST_BUDGET_SHARE = 8
const DEST_BYTES_PER_CELL = 40

"""
    DestTiling(runs, capsof, spacetiled)

How the destination is cut up for reading: tile `t` owns the cells `runs[t]`
names, and is bounded by the extents of the destination chunks `capsof[t]`.

`spacetiled` says tile `t` **is** destination chunk `t`, in which case `runs[t]`
is only the span of its cell positions and its cells come from
[`cellindices`](@ref) — the case a destination space with real chunking gets,
and the one whose blocks are keyed the way [`ChunkedPlan`](@ref) describes.
Otherwise a tile is exactly the contiguous run of cell positions `runs[t]`, cut
from the requested chunking or derived from the budget.
"""
struct DestTiling
    runs::Vector{UnitRange{Int}}
    capsof::Vector{Vector{Int}}
    spacetiled::Bool
end

ntiles(tiling::DestTiling) = length(tiling.runs)

# A destination space can tile the read only when its chunks partition the cell
# axis into ascending contiguous runs — and only when it has more than one, since
# a single chunk is a declaration that the space has no chunking to speak of and
# is exactly the case that must not collapse a whole destination into one read.
function _spacetileable(spans::Vector{UnitRange{Int}}, contiguous::Bool, ndst::Int)
    (contiguous && length(spans) > 1) || return false
    next = 1
    for sp in spans
        (!isempty(sp) && first(sp) == next) || return false
        next = last(sp) + 1
    end
    return next == ndst + 1
end

function _desttiling(dst_space::RegridSpace, ndst::Int, spans::Vector{UnitRange{Int}},
    contiguous::Bool, cellsizes::Union{Nothing,Vector{Int}}, budget::Int)
    if cellsizes === nothing && _spacetileable(spans, contiguous, ndst)
        return DestTiling(copy(spans), [[c] for c in eachindex(spans)], true)
    end
    sizes = cellsizes === nothing ?
            _defaulttilesizes(ndst, Int(nchunks(dst_space)), budget) : cellsizes
    runs = _runs(sizes, ndst)
    return DestTiling(runs, [_overlappingchunks(spans, r) for r in runs], false)
end

# A derived tile is budget-sized, but never coarser than the destination's own
# chunks: a space whose chunks interleave in position space cannot tile the read
# itself, and giving it one tile per chunk's worth of cells keeps the read
# granularity it meant to offer.
function _defaulttilesizes(ndst::Int, nchunk::Int, budget::Int)
    frombudget = max(1, fld(budget, DEST_BUDGET_SHARE * DEST_BYTES_PER_CELL))
    fromchunks = cld(ndst, max(nchunk, 1))
    return [clamp(min(frombudget, fromchunks), 1, ndst)]
end

# Cell-position runs from chunk sizes: one size repeats to cover the axis, a
# full list is taken as given.
function _runs(sizes::Vector{Int}, ndst::Int)
    all(>(0), sizes) || throw(ArgumentError(
        "destination chunk sizes must be positive, got $sizes"))
    out = UnitRange{Int}[]
    next = 1
    i = 0
    while next <= ndst
        i += 1
        n = sizes[min(i, length(sizes))]
        last_i = min(ndst, next + n - 1)
        push!(out, next:last_i)
        next = last_i + 1
    end
    isempty(out) && push!(out, 1:ndst)
    return out
end

# The destination chunks a run of cell positions can draw cells from. Spans
# bound a chunk's positions, so this is a superset for a space whose chunks
# interleave and exact for one whose chunks are contiguous.
function _overlappingchunks(spans::Vector{UnitRange{Int}}, run::UnitRange{Int})
    out = Int[]
    for (c, sp) in enumerate(spans)
        (isempty(sp) || last(sp) < first(run) || first(sp) > last(run)) && continue
        push!(out, c)
    end
    isempty(out) && append!(out, eachindex(spans))
    return out
end

# ===========================================================================
# Residency accounting
# ===========================================================================

"""
    LazyStats()

What one [`LazyRegridArray`](@ref) has done: source chunk loads issued, loads
answered from a held chunk instead, `(spatial × non-spatial)` chunk combinations
skipped as known-empty, chunk pairs dropped before their weights were built, and
the source bytes resident now and at their peak.

Accounting only — nothing here changes an answer, and nothing is updated per
element. It is how the bounded-residency law is checked without reading
`Sys.maxrss()`.
"""
mutable struct LazyStats
    loads::Int
    hits::Int
    skipped::Int
    dropped::Int
    bytes::Int
    peakbytes::Int
end

LazyStats() = LazyStats(0, 0, 0, 0, 0, 0)

Base.show(io::IO, s::LazyStats) =
    print(io, "LazyStats(loads=", s.loads, ", hits=", s.hits, ", skipped=", s.skipped,
        ", dropped=", s.dropped, ", peak=", s.peakbytes, " bytes)")

# The loaded source chunks one `readblock!` call may keep. Holding is what
# `budget` buys: a source chunk reached by several destination tiles of one
# request is read once rather than once per tile. Nothing survives the call —
# weights are cached across reads, data never is.
mutable struct SourceHold{A<:AbstractArray}
    held::Dict{Tuple{Int,Int},A}
    used::Dict{Tuple{Int,Int},Int}
    bytes::Int
    limit::Int
    clock::Int
    scratch::A
    stats::LazyStats
end

SourceHold(scratch::A, limit::Integer, stats::LazyStats) where {A<:AbstractArray} =
    SourceHold{A}(Dict{Tuple{Int,Int},A}(), Dict{Tuple{Int,Int},Int}(), 0,
        Int(limit), 0, scratch, stats)

_arraybytes(a::AbstractArray) = sizeof(eltype(a)) * length(a)

function _notepeak!(hold::SourceHold)
    s = hold.stats
    s.bytes = hold.bytes + _arraybytes(hold.scratch)
    s.peakbytes = max(s.peakbytes, s.bytes)
    return nothing
end

function _holdtake!(hold::SourceHold, key::Tuple{Int,Int})
    a = get(hold.held, key, nothing)
    a === nothing && return nothing
    hold.used[key] = (hold.clock += 1)
    hold.stats.hits += 1
    return a
end

function _holdstore!(hold::SourceHold, key::Tuple{Int,Int}, a::AbstractArray)
    n = _arraybytes(a)
    n <= hold.limit || return a
    while hold.bytes + n > hold.limit && !isempty(hold.held)
        _evictoldest!(hold)
    end
    hold.held[key] = a
    hold.used[key] = (hold.clock += 1)
    hold.bytes += n
    _notepeak!(hold)
    return a
end

function _evictoldest!(hold::SourceHold)
    victim = first(keys(hold.held))
    oldest = typemax(Int)
    for (k, u) in hold.used
        if u < oldest && haskey(hold.held, k)
            oldest = u
            victim = k
        end
    end
    hold.bytes -= _arraybytes(hold.held[victim])
    delete!(hold.held, victim)
    delete!(hold.used, victim)
    return hold
end

# ===========================================================================
# The array
# ===========================================================================

"""
    LazyRegridArray(data, plan::ChunkedPlan)

The regrid of `data` under `plan`, as a chunked `DiskArrays.AbstractDiskArray`
that computes each block on demand.

Dimensions are `(destination cells, data's non-spatial dims...)`, and chunks are
the destination tiling against the source's own chunking of the pass-through
dimensions — so one chunk of this array is one destination tile of one slice
group, which is exactly the unit `readblock!` computes.

The **destination tiling** is the destination space's own chunks when they
partition the cell axis into more than one ascending contiguous run, the plan's
`chunks` when it declares one, and a budget-sized partition of the cell axis
otherwise. A destination that reports a single whole-domain chunk therefore
still streams, which it cannot do if the space's chunking is taken as the last
word.

**Constructing one reads no source data.** Weights are built on first touch and
kept in the plan's storage, so the second read of a destination tile, and every
further non-spatial slice, builds nothing.

Reading a block loads only the source chunks discovered to reach it, and holds
several of them at once only while the plan's [`databudget`](@ref) allows;
otherwise it streams load-apply-drop, one source chunk resident at a time.
Accumulators are sized to the destination **tile**, never to the destination
space, and the source array is never materialized.
"""
struct LazyRegridArray{T,N,NS,NO,A,P<:ChunkedPlan,K,C} <: DiskArrays.AbstractDiskArray{T,N}
    source::A
    plan::P
    srcsize::NTuple{NS,Int}
    size::NTuple{N,Int}
    srctree::K
    dstcaps::Vector{Cap}
    tiling::DestTiling
    radius::Float64
    chunks::C
    otherchunks::NTuple{NO,Vector{UnitRange{Int}}}
    othergroups::Vector{NTuple{NO,UnitRange{Int}}}
    emptymemo::Vector{Int8}
    dropempty::Bool
    stats::LazyStats
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
    chunks, tiling = _outputgrid(plan, source, ndst, spans, contiguous, nspatial, othersizes)
    otherchunks = _sourceotherchunks(source, nspatial, othersizes)
    othergroups = _groupgrid(otherchunks)
    # The source's chunk tree is descended once per destination tile per read;
    # building it is O(nchunks) for a flat tree, so it is built once here.
    srctree = chunktree(src_space)
    T = outputeltype(eltype(data))
    return LazyRegridArray{T,length(othersizes) + 1,nspatial,length(othersizes),
        typeof(source),typeof(plan),typeof(srctree),typeof(chunks)}(
        source, plan, srcsize, (ndst, othersizes...), srctree, caps, tiling, radius,
        chunks, otherchunks, othergroups, zeros(Int8, Int(nchunks(src_space))),
        !usesreference(plan.missingpolicy), LazyStats())
end

Base.size(A::LazyRegridArray) = A.size

Base.show(io::IO, ::MIME"text/plain", A::LazyRegridArray) = show(io, A)
Base.show(io::IO, A::LazyRegridArray{T}) where {T} =
    print(io, "LazyRegridArray{", T, "}(", join(size(A), "×"), ", ",
        typeof(A.plan.method).name.name, ")")

DiskArrays.haschunks(::LazyRegridArray) = DiskArrays.Chunked()
DiskArrays.eachchunk(A::LazyRegridArray) = A.chunks

"""
    residency(A::LazyRegridArray) -> LazyStats

What `A` has read, held, skipped and dropped so far. Accounting only.
"""
residency(A::LazyRegridArray) = A.stats

# The span of cell positions each destination chunk owns, plus whether every
# chunk is contiguous — which decides whether the destination space can tile the
# read itself.
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

# The chunk grid the lazy array reports, and the destination tiling that has to
# agree with it: one chunk of the output is one tile of one slice group.
function _outputgrid(plan::ChunkedPlan, source, ndst::Int, spans::Vector{UnitRange{Int}},
    contiguous::Bool, nspatial::Int, othersizes::Tuple)
    nd = length(othersizes) + 1
    declared = plan.chunks
    cellsizes = nothing
    others = nothing
    if declared isa DiskArrays.GridChunks
        length(declared.chunks) == nd || throw(ArgumentError(
            "the plan's chunking has $(length(declared.chunks)) dimensions, but " *
            "the regrid of this source has $nd"))
        cellsizes = [length(r) for r in declared.chunks[1]]
        others = ntuple(i -> declared.chunks[i+1], nd - 1)
    elseif declared isa Tuple{Vararg{Integer}}
        length(declared) == nd || throw(ArgumentError(
            "the plan's chunk sizes $declared do not cover the regrid's $nd dimensions"))
        all(>(0), declared) || throw(ArgumentError(
            "chunk sizes must be positive, got $declared"))
        cellsizes = [Int(declared[1])]
        others = ntuple(i -> DiskArrays.RegularChunks(Int(declared[i+1]), 0, othersizes[i]),
            nd - 1)
    elseif declared !== nothing
        throw(ArgumentError(
            "`chunks` must be a tuple of chunk sizes, a DiskArrays.GridChunks, or " *
            "nothing, got $(typeof(declared))"))
    end
    tiling = _desttiling(plan.dst_space, ndst, spans, contiguous, cellsizes, plan.budget)
    cells = DiskArrays.IrregularChunks(; chunksizes = [length(r) for r in tiling.runs])
    passthrough = others === nothing ?
                  _passthroughchunks(source, nspatial, othersizes) : others
    return DiskArrays.GridChunks(cells, passthrough...), tiling
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

# The source's own chunking of its non-spatial dimensions, as index ranges. This
# is the grid `knownempty` is keyed on, and the grid a read is split along so a
# known-empty combination can be skipped without splitting anything else.
function _sourceotherchunks(source, nspatial::Int, othersizes::NTuple{NO,Int}) where {NO}
    if DiskArrays.haschunks(source) isa DiskArrays.Chunked
        ec = DiskArrays.eachchunk(source)
        if ec isa DiskArrays.GridChunks && length(ec.chunks) == nspatial + NO
            return ntuple(i -> UnitRange{Int}[UnitRange{Int}(r) for r in ec.chunks[nspatial+i]],
                NO)
        end
    end
    return ntuple(i -> UnitRange{Int}[1:othersizes[i]], NO)
end

# Every combination of one chunk per non-spatial dimension, in column-major
# order. A purely spatial source has exactly one, the empty tuple.
function _groupgrid(splits::NTuple{NO,Vector{UnitRange{Int}}}) where {NO}
    counts = map(length, splits)
    out = Vector{NTuple{NO,UnitRange{Int}}}(undef, prod(counts; init = 1))
    k = 0
    for I in CartesianIndices(counts)
        out[k+=1] = ntuple(d -> splits[d][I[d]], NO)
    end
    return out
end

# ===========================================================================
# readblock!
# ===========================================================================

function DiskArrays.readblock!(A::LazyRegridArray, aout::AbstractArray,
    r::AbstractUnitRange...)
    length(r) == ndims(A) || throw(DimensionMismatch(
        "$(length(r)) ranges requested from a $(ndims(A))-dimensional regrid"))
    cellr = UnitRange{Int}(r[1])
    others = map(UnitRange{Int}, Base.tail(r))
    nslices = prod(map(length, others); init = 1)
    return _readdestination!(reshape(aout, length(cellr), nslices), A, cellr, others, nslices)
end

"""
    _readdestination!(out, A, cellr, others, nslices)

The lazy path's one loop, into an `ncells × nslices` view of the output block.

Residency at any moment: the source chunks the budget lets the call hold, one
streamed source chunk beyond them, one `WeightBlock` (plus whatever the
storage's bound keeps), and the destination **tile's** accumulators — never the
destination space's, and never the source array's.
"""
function _readdestination!(out::AbstractMatrix, A::LazyRegridArray{T,N,NS,NO},
    cellr::UnitRange{Int}, others::NTuple{NO,UnitRange{Int}}, nslices::Int) where {T,N,NS,NO}
    plan = A.plan
    policy = plan.missingpolicy
    mv = plan.missingval
    groups = _slicegroups(A, others)
    strides = _slicestrides(others)
    hold = SourceHold(_emptybuffer(A), databudget(plan.budget), A.stats)
    srcchunks = Int[]
    srcranges = NTuple{NS,UnitRange{Int}}[]
    num = Matrix{Float64}(undef, 0, nslices)
    cover = Matrix{Float64}(undef, 0, nslices)
    total = Float64[]
    vals = T[]
    for t in _coveringtiles(A, cellr)
        dinds = _tileindices(A, t)
        nd = length(dinds)
        num = _fitmatrix(num, nd, nslices)
        cover = _fitmatrix(cover, nd, nslices)
        length(total) == nd || (total = Vector{Float64}(undef, nd))
        length(vals) == nd || (vals = Vector{T}(undef, nd))
        fill!(num, 0.0)
        fill!(cover, 0.0)
        fill!(total, 0.0)
        _connectedsource!(srcchunks, A, t)
        _sourceranges!(srcranges, A, srcchunks)
        keep = _canhold(hold, srcranges, groups, sizeof(eltype(hold.scratch)))
        denominated = false
        for (i, s) in enumerate(srcchunks)
            entry = blockfor(plan, (t, s), dinds)
            addreference!(total, entry.block)
            denominated |= hasdenom(entry.block)
            sr = srcranges[i]
            ncell = prod(map(length, sr))
            for (gi, pos) in enumerate(groups)
                gr = _grouprange(others, pos)
                if knownempty(A.source, (sr..., gr...))
                    A.stats.skipped += 1
                    continue
                end
                buf = _sourcefor!(hold, A, (s, gi), sr, gr, pos, keep)
                _applygroup!(num, cover, entry.block, entry.ref, buf, ncell, pos, strides, mv)
            end
        end
        _writechunk!(out, vals, num, cover, total, policy, denominated, dinds, cellr)
    end
    return out
end

# The destination tiles whose cells can fall in `cellr`. Runs bound a tile's
# positions, so this over-reports only for a destination space whose chunks
# interleave, and never misses one.
function _coveringtiles(A::LazyRegridArray, cellr::UnitRange{Int})
    out = Int[]
    lo, hi = first(cellr), last(cellr)
    runs = A.tiling.runs
    @inbounds for t in eachindex(runs)
        sp = runs[t]
        (isempty(sp) || last(sp) < lo || first(sp) > hi) && continue
        push!(out, t)
    end
    return out
end

# The cells a tile owns: the destination chunk's own, or the run itself.
_tileindices(A::LazyRegridArray, t::Int) =
    A.tiling.spacetiled ? cellindices(A.plan.dst_space, t) : A.tiling.runs[t]

# The source chunks that can reach a tile: the union of the descents against
# every destination chunk bounding it, minus the ones the source asserts hold no
# valid data at any non-spatial chunk.
function _connectedsource!(out::Vector{Int}, A::LazyRegridArray, t::Int)
    caps = A.tiling.capsof[t]
    if length(caps) == 1
        connectedchunks!(out, A.dstcaps[caps[1]], A.srctree; radius = A.radius)
    else
        connectedchunks!(out, view(A.dstcaps, caps), A.srctree; radius = A.radius)
    end
    if A.dropempty
        before = length(out)
        filter!(s -> !_allempty(A, s), out)
        A.stats.dropped += before - length(out)
    end
    return out
end

# Whether the source asserts that a spatial chunk holds no valid data at **any**
# non-spatial chunk — the only case in which a pair may be dropped before its
# weights are built. Memoized per source chunk: it is a property of the data's
# storage, not of the request.
function _allempty(A::LazyRegridArray, s::Int)
    memo = A.emptymemo
    m = memo[s]
    m == 0 || return m == 1
    sr = chunkranges(A.plan.src_space, s, A.srcsize)
    empty = true
    for g in A.othergroups
        knownempty(A.source, (sr..., g...)) || (empty = false; break)
    end
    memo[s] = empty ? Int8(1) : Int8(2)
    return empty
end

function _sourceranges!(out::Vector{NTuple{NS,UnitRange{Int}}}, A::LazyRegridArray{T,N,NS},
    srcchunks::Vector{Int}) where {T,N,NS}
    resize!(out, length(srcchunks))
    @inbounds for i in eachindex(srcchunks)
        out[i] = chunkranges(A.plan.src_space, srcchunks[i], A.srcsize)
    end
    return out
end

# Hold or stream. Holding a tile's whole connected set is what lets the next
# tile of the same request reuse it; a set that does not fit is streamed
# load-apply-drop instead, whose floor is one source chunk and one block.
function _canhold(hold::SourceHold, srcranges::Vector, groups::Vector, elbytes::Int)
    total = 0
    for sr in srcranges
        cells = prod(map(length, sr))
        for pos in groups
            total += elbytes * cells * prod(map(length, pos); init = 1)
            total > hold.limit && return false
        end
    end
    return true
end

# One source chunk over one non-spatial chunk, from the hold if it is there and
# from storage otherwise.
function _sourcefor!(hold::SourceHold, A::LazyRegridArray{T,N,NS,NO}, key::Tuple{Int,Int},
    sr::NTuple{NS,UnitRange{Int}}, gr::NTuple{NO,UnitRange{Int}},
    pos::NTuple{NO,UnitRange{Int}}, keep::Bool) where {T,N,NS,NO}
    got = _holdtake!(hold, key)
    got === nothing || return got
    shape = (map(length, sr)..., map(length, pos)...)
    A.stats.loads += 1
    if keep
        buf = similar(hold.scratch, shape)
        _readsource!(buf, A.source, sr, gr)
        _holdstore!(hold, key, buf)
        return buf
    end
    hold.scratch = _fitbuffer(hold.scratch, shape)
    _readsource!(hold.scratch, A.source, sr, gr)
    _notepeak!(hold)
    return hold.scratch
end

# One block against the slices of one loaded chunk combination. This is the
# function barrier: the block's concrete type is recovered here, once per
# combination, so the accumulation kernels below it are statically dispatched.
function _applygroup!(num::Matrix{Float64}, cover::Matrix{Float64}, block::WeightBlock,
    ref::Vector{Float64}, buf::AbstractArray, ncell::Int,
    pos::NTuple{NO,UnitRange{Int}}, strides::NTuple{NO,Int}, missingval) where {NO}
    src = reshape(buf, ncell, :)
    dirty = anyinvalid(buf, missingval)
    k = 0
    for I in CartesianIndices(pos)
        k += 1
        col = 1
        @inbounds for d in 1:NO
            col += (I[d] - 1) * strides[d]
        end
        x = view(src, :, k)
        applyblock!(view(num, :, col), view(cover, :, col), block, x,
            dirty ? x : nothing, ref, missingval)
    end
    return num
end

# Finalize once per destination tile per slice, then scatter the cells the
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
# Slice groups
# ===========================================================================

# The requested slices, split along the source's own non-spatial chunking, as
# positions within the requested ranges. Splitting there and nowhere else is
# what makes a known-empty chunk skippable without changing how anything else is
# read: a request that lies inside one source chunk is one group, as before.
_slicegroups(A::LazyRegridArray{T,N,NS,NO}, others::NTuple{NO,UnitRange{Int}}) where {T,N,NS,NO} =
    _groupgrid(ntuple(d -> _splitpositions(others[d], A.otherchunks[d]), NO))

function _splitpositions(r::UnitRange{Int}, chunks::Vector{UnitRange{Int}})
    out = UnitRange{Int}[]
    isempty(r) && return push!(out, 1:0)
    lo, hi = first(r), last(r)
    covered = 0
    for c in chunks
        a, b = max(first(c), lo), min(last(c), hi)
        a <= b || continue
        push!(out, (a-lo+1):(b-lo+1))
        covered += b - a + 1
    end
    # A grid that does not cover the request is not a partition of it; read the
    # request whole rather than silently dropping the uncovered part.
    covered == length(r) || return [1:length(r)]
    return out
end

# The linear stride of each non-spatial dimension in the flattened slice index.
_slicestrides(others::NTuple{NO,UnitRange{Int}}) where {NO} =
    ntuple(d -> prod(ntuple(i -> length(others[i]), d - 1); init = 1), NO)

# A group's positions back to index ranges of the source array.
_grouprange(others::NTuple{NO,UnitRange{Int}}, pos::NTuple{NO,UnitRange{Int}}) where {NO} =
    ntuple(d -> (first(others[d])+first(pos[d])-1):(first(others[d])+last(pos[d])-1), NO)

# ===========================================================================
# Source loading
# ===========================================================================

_emptybuffer(::LazyRegridArray{T,N,NS,NO,A}) where {T,N,NS,NO,A} =
    Array{eltype(A),NS + NO}(undef, ntuple(_ -> 0, NS + NO))

_fitmatrix(M::Matrix{Float64}, n::Int, m::Int) =
    size(M) == (n, m) ? M : Matrix{Float64}(undef, n, m)

function _fitbuffer(buf::Array{S,M}, shape::NTuple{M,Int}) where {S,M}
    return size(buf) == shape ? buf : Array{S,M}(undef, shape)
end

# One chunk-aligned read into the given buffer. A disk-backed source is asked
# through `readblock!`, which is one read per source chunk combination by
# construction.
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
tiles on demand and reads no source data until one is asked for.
"""
regrid(data, plan::ChunkedPlan) = LazyRegridArray(data, plan)

"""
    regrid!(dest, data, plan::ChunkedPlan) -> dest

Materialize a chunked plan into `dest`, one destination tile at a time.

Each tile's source chunks are loaded, applied and dropped before the next
tile's, so peak residency is a tile's worth of source and not the source
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
