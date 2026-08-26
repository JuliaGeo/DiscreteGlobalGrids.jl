# Lazy destination array and chunked execution.

# Source addressing

function chunkranges(space::RegridSpace, chunk::Integer, ::NTuple{1,Int})
    inds = ownedindices(space, Int(chunk))
    n = length(inds)
    n == 0 && return (1:0,)
    lo, hi = Int(first(inds)), Int(last(inds))
    hi - lo + 1 == n || _notrectangular(space, chunk)
    return (lo:hi,)
end

function chunkranges(space::RegridSpace, chunk::Integer, ssz::NTuple{2,Int})
    inds = ownedindices(space, Int(chunk))
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

# `chunkbox` is `(X, Y)`; array dimensions follow `xfast` order.
chunkranges(space::RasterGrid, chunk::Integer, ::NTuple{2,Int}) =
    space.xfast ? chunkbox(space, Int(chunk)) : reverse(chunkbox(space, Int(chunk)))

# Source chunk metadata

"""
    SourceChunking(source, ::Val{NS}, othersizes::NTuple{NO,Int})

What `source` declares about the chunking of its `NO` pass-through dimensions,
read once at construction.

  - Usable when `declareschunks(source)` holds and `DiskArrays.eachchunk`
    answers a `GridChunks` over the `NS` flattened dimensions followed by the
    `NO` pass-through ones; anything else makes every pass-through dimension
    one whole chunk.
  - Spatial chunking is not carried; a lazy read addresses it through
    [`chunkranges`](@ref).
  - `NS` arrives as a `Val` so each pass-through dimension keeps its chunk type.

# Fields

  - `passthrough`: the chunks as declared, which the output reports.
  - `splits`: the same chunks as index ranges, which a read splits a request along.
  - `groups`: every combination of those chunks, column-major.
"""
struct SourceChunking{NO,P<:Tuple}
    passthrough::P
    splits::NTuple{NO,Vector{UnitRange{Int}}}
    groups::Vector{NTuple{NO,UnitRange{Int}}}
end

function SourceChunking(source, ::Val{NS}, othersizes::NTuple{NO,Int}) where {NS,NO}
    if declareschunks(source)
        declared = DiskArrays.eachchunk(source)
        if declared isa DiskArrays.GridChunks && length(declared.chunks) == NS + NO
            passthrough = ntuple(i -> declared.chunks[NS+i], Val(NO))
            splits = ntuple(i -> UnitRange{Int}[UnitRange{Int}(r) for r in passthrough[i]],
                Val(NO))
            return SourceChunking{NO,typeof(passthrough)}(passthrough, splits,
                _groupgrid(splits))
        end
    end
    # Nothing usable: one whole chunk per pass-through dimension, said in both
    # of the spellings above.
    whole = ntuple(i -> DiskArrays.RegularChunks(max(othersizes[i], 1), 0, othersizes[i]),
        Val(NO))
    splits = ntuple(i -> UnitRange{Int}[1:othersizes[i]], Val(NO))
    return SourceChunking{NO,typeof(whole)}(whole, splits, _groupgrid(splits))
end

# Cartesian product of non-spatial chunks in column-major order.
function _groupgrid(splits::NTuple{NO,Vector{UnitRange{Int}}}) where {NO}
    counts = map(length, splits)
    out = Vector{NTuple{NO,UnitRange{Int}}}(undef, prod(counts; init = 1))
    k = 0
    for I in CartesianIndices(counts)
        out[k+=1] = ntuple(d -> splits[d][I[d]], NO)
    end
    return out
end

# Destination tiling

# Reserve one eighth of the budget for the per-cell state of one destination
# tile: the executor's accumulators, and whatever a build prepares per cell.
const DEST_BUDGET_SHARE = 8
const DEST_BYTES_PER_CELL = 40

"""
    destcellbudget(budget::Integer) -> Int

The bytes one destination tile's per-cell state may hold.
"""
destcellbudget(budget::Integer) = max(1, Int(budget) ÷ DEST_BUDGET_SHARE)

"""
    destcellbytes(method, dst_space, ndst::Int) -> Int

The bytes one destination cell holds while its tile is live: the executor's
per-cell accumulators, plus one prepared cell polygon and its flag where
`method` keeps them ([`preparedestination`](@ref)).

One cell is probed, because a space's cells are the same shape.
"""
function destcellbytes(method::AbstractRegriddingMethod, dst_space::RegridSpace,
    ndst::Int)
    (ndst > 0 && preparesdestination(method, dst_space)) ||
        return DEST_BYTES_PER_CELL
    return DEST_BYTES_PER_CELL + Int(Base.summarysize(getcell(dst_space, 1))) + 1
end

"""
    destcellsfit(method, dst_space, ndst::Integer, budget::Integer) -> Bool

Whether a tile of `ndst` destination cells holds its per-cell state, prepared
geometry included, inside [`destcellbudget`](@ref).

  - Tiling sizes tiles by the executor's own per-cell cost alone; preparing
    takes what is left of the share.
  - Prepared geometry that does not fit is refused here rather than charged to a
    budget that does not hold it.
"""
destcellsfit(method::AbstractRegriddingMethod, dst_space::RegridSpace,
    ndst::Integer, budget::Integer) =
    Int(ndst) * destcellbytes(method, dst_space, Int(ndst)) <= destcellbudget(budget)

"""
    DestTiling(runs, chunksof, spacetiled)

Describe destination tiles by cell-index runs and bounding space chunks.
When `spacetiled` is true, tile `t` is destination chunk `t`; otherwise the tile
is the corresponding contiguous run.
"""
struct DestTiling
    runs::Vector{UnitRange{Int}}
    chunksof::Vector{Vector{Int}}
    spacetiled::Bool
end

# Use space chunks as tiles only when multiple chunks partition the cell axis.
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

# Derived tiles honor both the budget and the destination's chunk granularity.
function _defaulttilesizes(ndst::Int, nchunk::Int, budget::Int)
    frombudget = max(1, fld(destcellbudget(budget), DEST_BYTES_PER_CELL))
    fromchunks = cld(ndst, max(nchunk, 1))
    return [clamp(min(frombudget, fromchunks), 1, ndst)]
end

# Convert repeated or explicit chunk sizes to cell-index runs.
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

# Return space chunks whose index spans overlap `run`.
function _overlappingchunks(spans::Vector{UnitRange{Int}}, run::UnitRange{Int})
    out = Int[]
    for (c, sp) in enumerate(spans)
        (isempty(sp) || last(sp) < first(run) || first(sp) > last(run)) && continue
        push!(out, c)
    end
    isempty(out) && append!(out, eachindex(spans))
    return out
end

# Residency accounting

"""
    LazyStats()

Track source loads, cache hits, skipped empty chunks, dropped pairs, and current
and peak source bytes. Statistics do not affect results.
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

# Source chunks held during one `readblock!` call. Data is not cached across reads.
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

# Prefetching destination tiles

"""
    TilePrefetch(ntiles, limit)

Builds the weights of upcoming destination tiles in background tasks while the
current tile's sources are read, so a sweep over a lazy array overlaps the two.

  - Starts after two consecutive tiles have been served; a single read never
    spawns anything. A read that jumps waits for the queued tasks and empties
    the queue.
  - Queue depth is `min(nthreads() - 1, limit ÷ largest - 1)`: `limit` is the
    storage's [`weightlimit`](@ref), `largest` the biggest tile built so far, so
    queued tiles plus the one in hand never exceed what the cache may hold.
    Depth is 0 on one thread, under [`OUTER_PARALLEL`](@ref), and before the
    first tile is built.
  - Tasks write only into the plan's storage under its lock; the queue is
    guarded by `lock`, so concurrent readers of one array are safe.
  - Source reads, their order, every `LazyStats` counter and the values are
    unchanged.
"""
mutable struct TilePrefetch
    ntiles::Int
    queue::Vector{StableTask{CachedTile}}   # `queue[k]` is tile `served + k`
    served::Int            # the tile last served, `0` before the first
    largest::Int           # bytes of the largest tile built so far
    limit::Int
    lock::ReentrantLock    # guards every field above
end

TilePrefetch(ntiles::Integer, limit::Integer) = TilePrefetch(
    Int(ntiles), StableTask{CachedTile}[], 0, 0, Int(limit), ReentrantLock())

Base.show(io::IO, p::TilePrefetch) =
    print(io, "TilePrefetch(", length(p.queue), " tiles queued after ", p.served, ")")

# The queue depth, in one place.
function _prefetchdepth(p::TilePrefetch)
    (p.largest > 0 && !OUTER_PARALLEL[]) || return 0
    return max(0, min(max(Threads.nthreads() - 1, 0), p.limit ÷ p.largest - 1))
end

# Lazy array

"""
    LazyRegridArray(data, plan::ChunkedPlan)

Return a chunked disk array that computes destination tiles on demand.

  - Dimensions are destination cells followed by the source's non-spatial
    dimensions. Tiling uses compatible destination chunks, declared `chunks`, or
    a budget-derived fallback.
  - Pass-through dimensions report the plan's declared chunking, and the
    source's own where there is none.
  - Construction reads no source data: one [`SourceChunking`](@ref) declaration
    and nothing else. Reads keep what they load within [`databudget`](@ref).
  - A tile's source chunks come from the plan's dependency relation alone — a
    destination chunk's row, [`sourcesof`](@ref)`(dependencies(plan), d)`, or the
    ascending union of a derived tile's rows. A read performs no geometric
    dependency discovery and issues no [`candidatechunks!`](@ref) query.
  - Dropping the chunks [`knownempty`](@ref) reports empty comes after that
    selection, never instead of it.
  - A plan built with `dependencies = false` cannot back a lazy array.
"""
struct LazyRegridArray{T,N,NS,NO,A,P<:ChunkedPlan,G<:ChunkDependencyGraph,C,
    S<:SourceChunking{NO}} <: DiskArrays.AbstractDiskArray{T,N}
    source::A
    plan::P
    srcsize::NTuple{NS,Int}
    size::NTuple{N,Int}
    # The plan's relation by reference (`graph === dependencies(plan)`): tile
    # adjacency and the per-chunk caps wave costing weighs, neither copied.
    graph::G
    tiling::DestTiling
    chunks::C
    # What the source declares about its own chunking, read once and read here.
    chunking::S
    emptymemo::Vector{Int8}
    dropempty::Bool
    stats::LazyStats
    # Carries prefetch state between reads: a sweep asks for one tile per read.
    prefetch::TilePrefetch
end

function LazyRegridArray(data, plan::ChunkedPlan)
    src_space, dst_space = plan.src_space, plan.dst_space
    graph = _lazygraph(plan)
    nsrc = Int(ncells(src_space))
    ndst = Int(ncells(dst_space))
    sd = resolvespatialdims(data, nsrc)
    othersizes = _otherdimsizes(data, sd)
    source = data isa DD.AbstractDimArray ? parent(data) : data
    nspatial = length(sd)
    srcsize = ntuple(i -> size(source, i), nspatial)
    spans, contiguous = _chunkspans(dst_space)
    chunking = SourceChunking(source, Val(nspatial), othersizes)
    chunks, tiling = _outputgrid(plan, chunking, ndst, spans, contiguous, othersizes)
    T = outputeltype(eltype(data))
    return LazyRegridArray{T,length(othersizes) + 1,nspatial,length(othersizes),
        typeof(source),typeof(plan),typeof(graph),typeof(chunks),typeof(chunking)}(
        source, plan, srcsize, (ndst, othersizes...), graph, tiling,
        chunks, chunking, zeros(Int8, Int(nchunks(src_space))),
        !usesreference(plan.missingpolicy), LazyStats(),
        TilePrefetch(length(tiling.runs), weightlimit(plan.storage)))
end

# Refuse a plan that owns no relation at construction, naming the keyword that
# caused it: every lazy read needs one, and on the pair route it is the read.
function _lazygraph(plan::ChunkedPlan)
    g = dependencies(plan)
    g === nothing && throw(ArgumentError(
        "a lazy regrid orders, costs and bounds its source reads with the " *
        "plan's dependency relation, but this plan owns none. It was built with " *
        "`dependencies = false`; drop that keyword to have the plan build one, " *
        "or pass a `ChunkDependencyGraph` for it to adopt."))
    hasextents(g) || throw(ArgumentError(
        "$g carries no chunk extents, so a lazy read cannot weigh a wave's " *
        "blocks. Build the relation from the spaces — `dependencies = true`, or " *
        "`chunk_dependency_graph(dst, src; radius)` — rather than assembling it " *
        "from bare CSR arrays."))
    # Restate the plan's validation in the numbering the tiling indexes.
    ndestinationchunks(g) == Int(nchunks(plan.dst_space)) || throw(ArgumentError(
        "$g holds $(ndestinationchunks(g)) destination rows for a " *
        "$(nchunks(plan.dst_space))-chunk destination space"))
    nsourcechunks(g) == Int(nchunks(plan.src_space)) || throw(ArgumentError(
        "$g holds $(nsourcechunks(g)) source chunks for a " *
        "$(nchunks(plan.src_space))-chunk source space"))
    return g
end

"""
    dependencies(A::LazyRegridArray) -> ChunkDependencyGraph

Return the relation `A` reads against, the identical object its plan holds:
`dependencies(A) === dependencies(A.plan)`.

  - It orders tiles, costs waves, holds the caps, and carries the refcounts and
    prefetch.
  - On the chunk-pair route it is also the source selection itself.
  - A tile built as one [`TileWeights`](@ref) selects from its own manifest
    instead; the relation is then a superset that decides no read of it.
"""
dependencies(A::LazyRegridArray) = A.graph

Base.size(A::LazyRegridArray) = A.size

Base.show(io::IO, ::MIME"text/plain", A::LazyRegridArray) = show(io, A)
Base.show(io::IO, A::LazyRegridArray{T}) where {T} =
    print(io, "LazyRegridArray{", T, "}(", join(size(A), "×"), ", ",
        typeof(A.plan.method).name.name, ")")

DiskArrays.haschunks(::LazyRegridArray) = DiskArrays.Chunked()
DiskArrays.eachchunk(A::LazyRegridArray) = A.chunks

"""
    residency(A::LazyRegridArray) -> LazyStats

Return `A`'s load, cache, skip, drop, and residency statistics.
"""
residency(A::LazyRegridArray) = A.stats

# Return chunk spans and whether every chunk is contiguous.
function _chunkspans(space::RegridSpace)
    nc = Int(nchunks(space))
    spans = Vector{UnitRange{Int}}(undef, nc)
    contiguous = true
    for c in 1:nc
        inds = ownedindices(space, c)
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

# Build the reported chunk grid and matching destination tiling. Pass-through
# dimensions report the plan's declared chunking, or the source's where none.
function _outputgrid(plan::ChunkedPlan, chunking::SourceChunking, ndst::Int,
    spans::Vector{UnitRange{Int}}, contiguous::Bool, othersizes::Tuple)
    nd = length(othersizes) + 1
    declared = plan.chunks
    # Plans built without `plan_regrid` have not passed the API-boundary check.
    _checkchunks(declared)
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
        cellsizes = [Int(declared[1])]
        others = ntuple(i -> DiskArrays.RegularChunks(Int(declared[i+1]), 0, othersizes[i]),
            nd - 1)
    end
    tiling = _desttiling(plan.dst_space, ndst, spans, contiguous, cellsizes, plan.budget)
    cells = DiskArrays.IrregularChunks(; chunksizes = [length(r) for r in tiling.runs])
    passthrough = others === nothing ? chunking.passthrough : others
    return DiskArrays.GridChunks(cells, passthrough...), tiling
end

# Block reads

function DiskArrays.readblock!(A::LazyRegridArray, aout::AbstractArray,
    r::AbstractUnitRange...)
    length(r) == ndims(A) || throw(DimensionMismatch(
        "$(length(r)) ranges requested from a $(ndims(A))-dimensional regrid"))
    cellr = UnitRange{Int}(r[1])
    others = map(UnitRange{Int}, Base.tail(r))
    nslices = prod(map(length, others); init = 1)
    return _readdestination!(reshape(aout, length(cellr), nslices), A, cellr, others, nslices)
end

# Serving and queueing tiles. The reading loop below is their only caller; they
# name the array because a queued tile takes its cells off it.

# Serve tile `t` from the queue where it holds it and by building here where it
# does not, then queue the tiles behind it. Both guarded regions leave `queue[k]`
# naming tile `served + k`, so a take pops the tile it asked for however takes
# interleave; the wait for the tile is between them, holding no lock.
function _taketile!(p::TilePrefetch, A::LazyRegridArray, plan::ChunkedPlan, t::Int,
    dinds, smp)
    prev = 0
    task = nothing
    @lock p.lock begin
        prev = p.served
        if t == prev + 1 && !isempty(p.queue)
            task = popfirst!(p.queue)
        else
            _drain!(p)
        end
        p.served = t
    end
    tile = task === nothing ? tilefor(plan, t, dinds, smp) : _fetchtile(task)
    @lock p.lock begin
        p.largest = max(p.largest, tile.bytes)
        # Two tiles in succession are the evidence that a sweep is running.
        (prev >= 1 && t == prev + 1) && _prefetch!(p, A, plan, smp)
    end
    return tile
end

# Raise what the build raised, not the wrapper a task failure arrives in.
function _fetchtile(task::StableTask{CachedTile})
    try
        return fetch(task)
    catch err
        err isa TaskFailedException && throw(err.task.result)
        rethrow()
    end
end

# Queue the tiles behind the one served, up to the depth. Called under `p.lock`.
function _prefetch!(p::TilePrefetch, A::LazyRegridArray, plan::ChunkedPlan, smp)
    depth = _prefetchdepth(p)
    while length(p.queue) < depth
        next = p.served + length(p.queue) + 1
        next <= p.ntiles || break
        # The reading loop reads the destination space too, so a tile's cells
        # are resolved here rather than in the task.
        dinds = _tileindices(A, next)
        # A build inside the task sees the scope and keeps its own loops serial.
        # The typeassert is what the task's return type is inferred from, so the
        # queue holds a `StableTask{CachedTile}` whatever `tilefor` infers to.
        push!(p.queue, @with OUTER_PARALLEL => true StableTasks.@spawn tilefor(
            plan, next, dinds, smp)::CachedTile)
    end
    return p
end

# Wait for every queued task and forget them. A build that failed raises where
# its tile is served, not here.
function _drain!(p::TilePrefetch)
    @lock p.lock begin
        for task in p.queue
            try
                wait(task)
            catch
            end
        end
        empty!(p.queue)
    end
    return p
end

"""
    _readdestination!(out, A, cellr, others, nslices)

Read destination cells and slices into `out`.

  - Source residency is limited to held chunks plus one streamed chunk.
  - Blocks are built in waves but applied in chunk order, so results are
    independent of thread count.
  - A point method with a [`sampler`](@ref) builds one [`TileWeights`](@ref) per
    tile before any source read, because its manifest is what selects the
    chunks; its blocks are applied in the same ascending chunk order.
  - Tiles not yet asked for are built concurrently, through the array's
    [`TilePrefetch`](@ref).
"""
function _readdestination!(out::AbstractMatrix, A::LazyRegridArray{T,N,NS,NO},
    cellr::UnitRange{Int}, others::NTuple{NO,UnitRange{Int}}, nslices::Int) where {T,N,NS,NO}
    plan = A.plan
    policy = plan.missingpolicy
    mv = plan.missingval
    # One sampler serves every tile of this read; a plan on the pair route has
    # none and builds a block at a time.
    smp = tilesampler(plan)
    groups = _slicegroups(A, others)
    strides = _slicestrides(others)
    hold = SourceHold(_emptybuffer(A), databudget(plan.budget), A.stats)
    srcchunks = Int[]
    srcranges = NTuple{NS,UnitRange{Int}}[]
    wave = CachedBlock[]
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
        # On the tile route the weights come before the selection, because the
        # chunks that carry one are exactly the chunks read.
        tile = smp === nothing ? nothing : _taketile!(A.prefetch, A, plan, t, dinds, smp)
        _connectedsource!(srcchunks, A, t, tile)
        _sourceranges!(srcranges, A, srcchunks)
        # Prepare the tile's destination geometry once where more than one block
        # reads it; a single block uses its own task-local memo instead. The
        # choice crosses an inference barrier so that only the assembly this run
        # takes is compiled; `_fillwave!` specializes on what it is handed.
        destination = Base.inferencebarrier(
            (tile === nothing && length(srcchunks) > 1) ?
            preparedestination(plan.method, plan.dst_space, dinds, plan.budget) :
            dinds)
        keep = _canhold(hold, srcranges, groups, sizeof(eltype(hold.scratch)))
        w = tile === nothing ?
            _wavesize(plan, nd, srcchunks, srcranges, A.graph, A.tiling.chunksof[t]) : 1
        i = 1
        while i <= length(srcchunks)
            j = min(i + w - 1, length(srcchunks))
            if tile === nothing
                _fillwave!(wave, plan, t, srcchunks, i, j, dinds, destination)
            else
                _tilewave!(wave, tile, srcchunks, i, j)
            end
            for k in i:j
                entry = wave[k-i+1]
                s = srcchunks[k]
                addreference!(total, entry.block)
                sr = srcranges[k]
                ncell = prod(map(length, sr))
                for (gi, pos) in enumerate(groups)
                    gr = _grouprange(others, pos)
                    if knownempty(A.source, (sr..., gr...))
                        A.stats.skipped += 1
                        continue
                    end
                    buf = _sourcefor!(hold, A, (s, gi), sr, gr, pos, keep)
                    _applygroup!(num, cover, entry.block, buf, ncell, pos, strides, mv)
                end
            end
            # Retain at most one local wave; storage manages its own cache.
            empty!(wave)
            i = j + 1
        end
        _writechunk!(out, vals, num, cover, total, policy, dinds, cellr)
    end
    return out
end

"""
    _wavesize(plan, nd, srcchunks, srcranges, graph, rows) -> Int

Return the number of chunk-pair blocks to build concurrently; one on a
single-threaded session.

  - A wave spawns one task per block and declares [`OUTER_PARALLEL`](@ref),
    which turns off threading *inside* each build, so the two are alternatives
    rather than additions.
  - Under an outer loop that already declared `OUTER_PARALLEL` the width is
    bounded only by thread count, block count and the weight budget.
  - At top level the wave must beat inner threading — [`_waveideal`](@ref) of
    [`_blockcosts!`](@ref) against `innerspeedup * nthreads()` — or this returns
    one and the build gets the threads.
  - `innerspeedup` is the per-thread speedup one weight build reaches on its own.
  - The caps the estimate needs are read off `graph`'s
    [`chunkextents`](@ref) and never copied, so a wave's tasks share them;
    `rows` are the tile's destination rows.
"""
function _wavesize(plan::ChunkedPlan, nd::Int, srcchunks::Vector{Int}, srcranges::Vector,
    graph::ChunkDependencyGraph, rows::Vector{Int}; innerspeedup::Float64 = 0.73)
    n = length(srcchunks)
    nt = Threads.nthreads()
    (nt > 1 && n > 1) || return 1
    ncols = 0
    for sr in srcranges
        ncols = max(ncols, prod(map(length, sr)))
    end
    floorbytes = 8 * (ncols + 1) + 16 * nd + 64
    fits = max(1, weightbudget(plan.budget) ÷ floorbytes)
    w = Int(min(nt, n, fits))
    # An outer loop already owns the cores; the wave is all that is left here.
    OUTER_PARALLEL[] && return w
    costs = _blockcosts!(Vector{Float64}(undef, n), srcchunks, srcranges, graph, rows)
    return _waveideal(costs, w) >= innerspeedup * nt ? w : 1
end

"""
    _blockcosts!(costs, srcchunks, srcranges, graph, rows) -> costs

Estimate each chunk pair's build cost before building it, as the source cells
the destination tile can reach: the chunk's cell count scaled by the fraction of
its extent a tile extent overlaps.

  - Cell counts alone do not track build cost, because chunks of one source are
    near enough the same size while their overlap with a tile is not.
  - Extents are bounding caps, so the estimate is an upper bound on how evenly
    the work spreads: a chunk touching only the tile's cap, and so contributing
    nothing, still scores above zero.
  - That bias favours keeping the wave.
"""
function _blockcosts!(costs::Vector{Float64}, srcchunks::Vector{Int}, srcranges::Vector,
    graph::ChunkDependencyGraph, rows::Vector{Int})
    srccaps = sourceextents(graph)
    dstcaps = destinationextents(graph)
    @inbounds for k in eachindex(srcchunks)
        scap = srccaps[srcchunks[k]]
        area = _caparea(Float64(scap.radius))
        ncell = prod(map(length, srcranges[k]))
        shared = 0.0
        for d in rows
            dcap = dstcaps[destinationchunk(graph, d)]
            # Tile extents may overlap each other, so take the widest reach
            # rather than summing and double counting the shared part.
            shared = max(shared, _capoverlap(Float64(dcap.radius),
                Float64(scap.radius), _capdistance(dcap, scap)))
        end
        costs[k] = area > 0 ? ncell * min(shared / area, 1.0) : 0.0
    end
    return costs
end

# The best speedup a wave of `w` can reach: waves run back to back and each one
# takes as long as its slowest block.
function _waveideal(costs::Vector{Float64}, w::Int)
    total = sum(costs)
    total > 0 || return 1.0
    serial = 0.0
    i = 1
    @inbounds while i <= length(costs)
        j = min(i + w - 1, length(costs))
        serial += maximum(view(costs, i:j))
        i = j + 1
    end
    return serial > 0 ? total / serial : 1.0
end

# The tile route's wave is a lookup: the tile's block for each source chunk, in
# the ascending order the loop walks them. `srcchunks` came from the manifest,
# so every one of them has a block.
function _tilewave!(wave::Vector{CachedBlock}, tile::CachedTile, srcchunks::Vector{Int},
    i::Int, j::Int)
    empty!(wave)
    for k in i:j
        push!(wave, tileblock(tile, srcchunks[k])::CachedBlock)
    end
    return wave
end

# Spawn one block build. A task inherits the scope, so a nested build sees the
# wave and stays serial; the typeassert is what its return type is inferred from,
# so a wave's handles are `StableTask{CachedBlock}`s whatever `blockfor` infers
# to.
_spawnblock(plan::ChunkedPlan, key::Tuple{Int,Int}, dinds, destination) =
    @with OUTER_PARALLEL => true StableTasks.@spawn blockfor(
        plan, key, dinds, destination)::CachedBlock

# Build one wave concurrently and preserve chunk order in `wave`.
function _fillwave!(wave::Vector{CachedBlock}, plan::ChunkedPlan, t::Int,
    srcchunks::Vector{Int}, i::Int, j::Int, dinds, destination)
    empty!(wave)
    if i == j
        # No spawn and no declaration: a wave of one leaves the threads to the
        # build itself, which is what `_wavesize` returning one is asking for.
        push!(wave, blockfor(plan, (t, srcchunks[i]), dinds, destination))
        return wave
    end
    tasks = map(k -> _spawnblock(plan, (t, srcchunks[k]), dinds, destination), i:j)
    # Every spawned task is waited for, whatever happens, so no build outlives
    # the read that spawned it. The first exception is the one raised.
    err = nothing
    for task in tasks
        try
            block = fetch(task)
            err === nothing && push!(wave, block)
        catch e
            err === nothing && (err = e)
        end
    end
    err === nothing || throw(err)
    return wave
end

# Return destination tiles whose index spans overlap `cellr`.
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

# Return cells owned by a tile.
_tileindices(A::LazyRegridArray, t::Int) =
    A.tiling.spacetiled ? ownedindices(A.plan.dst_space, t) : A.tiling.runs[t]

"""
    _connectedsource!(out, A, t[, tile]) -> out

Write tile `t`'s source chunks into `out`, ascending. This is the whole of the
executor's source selection, and the tile's build unit decides it.

  - A tile with no [`TileWeights`](@ref) — every area method, and a point method
    supplying no [`sampler`](@ref) — reads the plan's relation, with no index,
    no query and no cap test: row `d` for destination chunk `d`, or the
    ascending union of the rows a derived tile spans, which depends on the
    tiling and not on the order the rows are visited.
  - A tile with `TileWeights` takes its manifest exactly — the sorted union of
    the source chunks owning its nonzero stencil entries. The relation's row is
    a superset of it and decides nothing here.
  - A manifest naming a chunk no row of the tile holds is a declared
    [`supportradius`](@ref) too small to bound the method's stencils, and is
    refused rather than silently dropped.
  - `knownempty` filtering comes **after**, on both routes: it may drop a chunk
    the selection holds, and may never add one the selection does not.
"""
_connectedsource!(out::Vector{Int}, A::LazyRegridArray, t::Int) =
    _connectedsource!(out, A, t, nothing)

function _connectedsource!(out::Vector{Int}, A::LazyRegridArray, t::Int, tile)
    _selectsource!(out, A, t, tile)
    if A.dropempty
        before = length(out)
        filter!(s -> !_allempty(A, s), out)
        A.stats.dropped += before - length(out)
    end
    return out
end

# The relation's rows: one row, or a derived tile's ascending union of rows.
function _selectsource!(out::Vector{Int}, A::LazyRegridArray, t::Int, ::Nothing)
    rows = A.tiling.chunksof[t]
    if length(rows) == 1
        _copyrow!(out, sourcesof(A.graph, @inbounds rows[1]))
    else
        _unionrows!(out, A.graph, rows)
    end
    return out
end

# The tile's own manifest, checked against the rows it has to sit inside: a
# chunk the relation never named would be read outside every refcount.
function _selectsource!(out::Vector{Int}, A::LazyRegridArray, t::Int, tile::CachedTile)
    _copyrow!(out, tile.sourcechunks)
    s = _outsiderows(A, t, out)
    s == 0 || throw(ArgumentError(_reachmessage(A, t, s)))
    return out
end

# The first chunk of `manifest` that no row of tile `t` holds, or `0` when every
# one of them is inside. Rows are ascending, so each test is a binary search.
function _outsiderows(A::LazyRegridArray, t::Int, manifest::Vector{Int})
    rows = A.tiling.chunksof[t]
    @inbounds for s in manifest
        inside = false
        for d in rows
            _inrow(sourcesof(A.graph, d), s) && (inside = true; break)
        end
        inside || return s
    end
    return 0
end

function _inrow(row, s::Int)
    k = searchsortedfirst(row, s)
    return k <= length(row) && Int(@inbounds row[k]) == s
end

function _reachmessage(A::LazyRegridArray, t::Int, s::Int)
    method, src_space = A.plan.method, A.plan.src_space
    return "destination tile $t takes weights from source chunk $s, which the " *
           "plan's dependency relation does not name for it. A method that " *
           "supplies a `sampler` must declare a `supportradius` bounding every " *
           "stencil it emits, so that the relation stays a superset of what is " *
           "read; `supportradius` of $(typeof(method)) on $(typeof(src_space)) " *
           "answers $(supportradius(method, src_space)) radians, which does not " *
           "reach chunk $s."
end

# One ascending list of chunk numbers, widened to the `Int` every consumer
# downstream uses: a relation row, or a tile's manifest.
function _copyrow!(out::Vector{Int}, row)
    resize!(out, length(row))
    @inbounds for i in eachindex(row)
        out[i] = Int(row[i])
    end
    return out
end

# The ascending union of several rows: concatenate, then sort in place, so the
# answer is the same for any row order.
function _unionrows!(out::Vector{Int}, g::ChunkDependencyGraph, rows)
    empty!(out)
    @inbounds for d in rows
        row = sourcesof(g, d)
        n = length(out)
        resize!(out, n + length(row))
        for i in eachindex(row)
            out[n+i] = Int(row[i])
        end
    end
    sort!(out)
    unique!(out)
    return out
end

# Memoize whether a spatial chunk is empty across all non-spatial chunks.
function _allempty(A::LazyRegridArray, s::Int)
    memo = A.emptymemo
    m = memo[s]
    m == 0 || return m == 1
    sr = chunkranges(A.plan.src_space, s, A.srcsize)
    empty = true
    for g in A.chunking.groups
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

# Return whether all needed source groups fit in the hold budget.
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

# Load one spatial/non-spatial source chunk combination, using the hold when possible.
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

# Apply one block to every slice in a loaded chunk combination. Coverage comes
# off the block's own reference weights, which every slice shares.
function _applygroup!(num::Matrix{Float64}, cover::Matrix{Float64}, block::WeightBlock,
    buf::AbstractArray, ncell::Int,
    pos::NTuple{NO,UnitRange{Int}}, strides::NTuple{NO,Int}, missingval) where {NO}
    ref = block.reference
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

# Finalize each tile slice and scatter requested cells.
function _writechunk!(out::AbstractMatrix, vals::Vector, num::Matrix{Float64},
    cover::Matrix{Float64}, total::Vector{Float64}, policy::AbstractMissingPolicy,
    dinds, cellr::AbstractUnitRange)
    lo, hi = first(cellr), last(cellr)
    off = lo - 1
    for t in axes(num, 2)
        finalize!(vals, view(num, :, t), view(cover, :, t), total, policy)
        @inbounds for (j, p) in enumerate(dinds)
            lo <= p <= hi || continue
            out[p-off, t] = vals[j]
        end
    end
    return out
end

# Slice groups

# Split requested slices along source non-spatial chunk boundaries.
_slicegroups(A::LazyRegridArray{T,N,NS,NO}, others::NTuple{NO,UnitRange{Int}}) where {T,N,NS,NO} =
    _groupgrid(ntuple(d -> _splitindices(others[d], A.chunking.splits[d]), NO))

function _splitindices(r::UnitRange{Int}, chunks::Vector{UnitRange{Int}})
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
    # Fall back to one group when the chunk grid does not cover the request.
    covered == length(r) || return [1:length(r)]
    return out
end

# Linear strides for flattened non-spatial slices.
_slicestrides(others::NTuple{NO,UnitRange{Int}}) where {NO} =
    ntuple(d -> prod(ntuple(i -> length(others[i]), d - 1); init = 1), NO)

# Convert group-relative indices to source ranges.
_grouprange(others::NTuple{NO,UnitRange{Int}}, pos::NTuple{NO,UnitRange{Int}}) where {NO} =
    ntuple(d -> (first(others[d])+first(pos[d])-1):(first(others[d])+last(pos[d])-1), NO)

# Source loading

_emptybuffer(::LazyRegridArray{T,N,NS,NO,A}) where {T,N,NS,NO,A} =
    Array{eltype(A),NS + NO}(undef, ntuple(_ -> 0, NS + NO))

_fitmatrix(M::Matrix{Float64}, n::Int, m::Int) =
    size(M) == (n, m) ? M : Matrix{Float64}(undef, n, m)

function _fitbuffer(buf::Array{S,M}, shape::NTuple{M,Int}) where {S,M}
    return size(buf) == shape ? buf : Array{S,M}(undef, shape)
end

# Read one aligned source chunk combination into `buf`.
function _readsource!(buf::Array, source, sr::Tuple, others::Tuple)
    if _isdisksource(source)
        DiskArrays.readblock!(source, buf, sr..., others...)
    else
        copyto!(buf, view(source, sr..., others...))
    end
    return buf
end

# Whether the source's values are reached through a DiskArrays read. This is not
# `declareschunks`: an in-memory array may declare a chunking and is copied from.
_isdisksource(x) = x isa DiskArrays.AbstractDiskArray || DiskArrays.isdisk(x)

# Shaped view

"""
    ShapedRegridArray(parent::LazyRegridArray, shape)

Lazily split `parent`'s leading cell axis into `shape`, column-major, so cell
indices keep their order. Reads forward to `parent` without materializing it.
"""
struct ShapedRegridArray{T,N,ND,P<:LazyRegridArray} <: DiskArrays.AbstractDiskArray{T,N}
    parent::P
    shape::NTuple{ND,Int}
    size::NTuple{N,Int}
end

function ShapedRegridArray(A::LazyRegridArray{T}, shape::Tuple{Vararg{Int}}) where {T}
    prod(shape) == size(A, 1) || throw(DimensionMismatch(
        "$(size(A, 1)) cells do not reshape to $shape"))
    sz = (shape..., Base.tail(size(A))...)
    return ShapedRegridArray{T,length(sz),length(shape),typeof(A)}(A, shape, sz)
end

Base.size(A::ShapedRegridArray) = A.size
Base.parent(A::ShapedRegridArray) = A.parent
residency(A::ShapedRegridArray) = residency(A.parent)

Base.show(io::IO, ::MIME"text/plain", A::ShapedRegridArray) = show(io, A)
Base.show(io::IO, A::ShapedRegridArray{T}) where {T} =
    print(io, "ShapedRegridArray{", T, "}(", join(size(A), "×"), ")")

DiskArrays.haschunks(::ShapedRegridArray) = DiskArrays.Chunked()

function DiskArrays.eachchunk(A::ShapedRegridArray{T,N,ND}) where {T,N,ND}
    pc = DiskArrays.eachchunk(A.parent).chunks
    return DiskArrays.GridChunks(_shapedchunks(pc[1], A.shape)..., Base.tail(pc)...)
end

# Full leading axes; the parent's cell tiling projects onto the last spatial
# axis where its run boundaries land on whole columns.
function _shapedchunks(cells, shape::NTuple{ND,Int}) where {ND}
    lead = ntuple(d -> DiskArrays.RegularChunks(shape[d], 0, shape[d]), ND - 1)
    ncol = prod(Base.front(shape); init = 1)
    offs = Int[0]
    for r in cells
        q, rem = divrem(last(r), ncol)
        rem == 0 && last(offs) < q < shape[ND] && push!(offs, q)
    end
    push!(offs, shape[ND])
    return (lead..., DiskArrays.IrregularChunks(; chunksizes = diff(offs)))
end

# A request is one contiguous cell run when partial axes trail only singletons.
function _iscontiguous(shape::NTuple{ND,Int}, rs::NTuple{ND,UnitRange{Int}}) where {ND}
    partial = false
    for d in 1:ND
        partial && length(rs[d]) > 1 && return false
        partial |= length(rs[d]) != shape[d]
    end
    return true
end

function DiskArrays.readblock!(A::ShapedRegridArray{T,N,ND}, aout::AbstractArray,
    r::AbstractUnitRange...) where {T,N,ND}
    length(r) == N || throw(DimensionMismatch(
        "$(length(r)) ranges requested from a $N-dimensional regrid"))
    rs = ntuple(d -> UnitRange{Int}(r[d]), ND)
    others = ntuple(d -> UnitRange{Int}(r[ND+d]), N - ND)
    strides = ntuple(d -> prod(ntuple(i -> A.shape[i], d - 1); init = 1), ND)
    lo = 1 + sum(ntuple(d -> (first(rs[d]) - 1) * strides[d], ND); init = 0)
    hi = 1 + sum(ntuple(d -> (last(rs[d]) - 1) * strides[d], ND); init = 0)
    if _iscontiguous(A.shape, rs)
        DiskArrays.readblock!(A.parent,
            reshape(aout, hi - lo + 1, map(length, others)...), lo:hi, others...)
        return aout
    end
    # Rectangles that span partial columns read the covering cell run once and
    # keep the requested lattice rows.
    tmp = Array{T}(undef, hi - lo + 1, map(length, others)...)
    DiskArrays.readblock!(A.parent, tmp, lo:hi, others...)
    tmp2 = reshape(tmp, hi - lo + 1, :)
    out2 = reshape(aout, prod(map(length, rs)), :)
    row = 0
    for I in CartesianIndices(rs)
        row += 1
        flat = 1 + sum(ntuple(d -> (I[d] - 1) * strides[d], ND); init = 0)
        out2[row, :] = view(tmp2, flat - lo + 1, :)
    end
    return aout
end

# Mirror `wrapoutput` without materializing: label a dimensional source's regrid
# with the destination's own axes, or one flat `Cell` axis. A destination naming
# one axis is labelled directly, with no shaped view between it and its reads.
function wraplazy(A::LazyRegridArray{T,N,NS}, data, dstdims) where {T,N,NS}
    data isa DD.AbstractDimArray || return A
    ds = DD.dims(data)
    others = ntuple(i -> ds[NS+i], ndims(data) - NS)
    dstdims === nothing &&
        return DD.DimArray(A, (DD.Dim{:Cell}(1:size(A, 1)), others...))
    length(dstdims) == 1 && return DD.DimArray(A, (dstdims..., others...))
    shaped = ShapedRegridArray(A, map(length, dstdims))
    return DD.DimArray(shaped, (dstdims..., others...))
end

# API

"""
    regrid(data, plan::ChunkedPlan) -> lazy array

Return a disk-backed array that computes destination tiles on demand.

  - Labels and shape match the eager result ([`wrapoutput`](@ref)): a
    dimensional source comes back as a `DimArray` carrying the destination's own
    axes in their construction order, or a flat `Cell` axis.
  - Several destination axes are a [`ShapedRegridArray`](@ref) view of the cell
    axis; one axis, or none, is the [`LazyRegridArray`](@ref) itself.
  - Other sources stay a flat [`LazyRegridArray`](@ref).
"""
function regrid(data, plan::ChunkedPlan)
    A = LazyRegridArray(data, plan)
    return wraplazy(A, data, destinationdims(plan))
end

"""
    regrid!(dest, data, plan::ChunkedPlan) -> dest

Materialize a chunked plan into `dest`, one destination tile at a time. As for
[`DirectPlan`](@ref), `dest` may lead with the destination's own axes or one
flat cell dimension.
"""
function regrid!(dest, data, plan::ChunkedPlan)
    A = LazyRegridArray(data, plan)
    dstdims = destinationdims(plan)
    shaped = dstdims === nothing ? size(A) :
             (map(length, dstdims)..., Base.tail(size(A))...)
    size(dest) == shaped || size(dest) == size(A) || throw(DimensionMismatch(
        "destination of size $(size(dest)) cannot hold a regrid of size $shaped"))
    raw = dest isa DD.AbstractDimArray ? parent(dest) : dest
    flat = reshape(raw, size(A))
    for ranges in DiskArrays.eachchunk(A)
        flat[ranges...] = A[ranges...]
    end
    return dest
end
