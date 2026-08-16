# Plans. A plan carries everything execution needs — method, spaces, missing
# policy, weight storage — which is why applying one takes no keyword arguments.
# `DirectPlan` is the whole-domain case; `ChunkedPlan` and the storage policies
# below are declared here and implemented with the lazy path.

"""
    AbstractRegriddingPlan

A complete, reusable regridding operator: method, both spaces, missing policy,
and the weights themselves.

A plan is **self-contained** — applying one takes no keyword arguments, because
there is nothing left to decide. It is also immutable in everything that
affects the answer; only weight *storage* changes as blocks are built on first
touch, so one plan may be applied concurrently to different slices.

[`DirectPlan`](@ref) holds one whole-domain block; [`ChunkedPlan`](@ref) holds
blocks per chunk pair.
"""
abstract type AbstractRegriddingPlan end

"""
    WeightBlock(weights, denom)
    WeightBlock(coo::WeightCOO, ndst::Integer, nsrc::Integer)

The weights of one destination-chunk/source-chunk pair, immutable once built.

`weights` is an `ndst × nsrc` sparse matrix in **chunk-local** indices:
`weights[j, k]` is the weight of the `k`-th cell of the source chunk in the
`j`-th cell of the destination chunk. `denom` is either a length-`ndst` vector
of per-destination denominators — covered area, for [`Conservative`](@ref) —
or `nothing` when the method reported none, in which case the block finalizes
as a raw sum.

Building from a [`WeightCOO`](@ref) sums duplicate entries.

Blocks column-partition the operator by source chunk. That is what lets every
method stream: applying a block to one loaded source chunk yields that chunk's
partial contribution, and summing the partials over the connected source chunks
reconstructs the whole, with no per-method streaming code and no mosaic
container anywhere in the contract.
"""
struct WeightBlock{M<:AbstractMatrix{Float64},D<:Union{Nothing,Vector{Float64}}}
    weights::M
    denom::D
end

function WeightBlock(coo::WeightCOO, ndst::Integer, nsrc::Integer)
    weights = sparse(coo.rows, coo.cols, coo.vals, Int(ndst), Int(nsrc))
    return WeightBlock(weights, coo.hasdenom ? copy(coo.denom) : nothing)
end

Base.size(block::WeightBlock) = size(block.weights)
Base.size(block::WeightBlock, d::Integer) = size(block.weights, d)

"""
    hasdenom(block::WeightBlock) -> Bool

Whether `block` carries a per-destination denominator, which is what decides
whether [`Weighted`](@ref) normalizes its destinations or only blanks the
under-covered ones. See `finalize!`.
"""
hasdenom(block::WeightBlock) = block.denom !== nothing

Base.show(io::IO, block::WeightBlock) =
    print(io, "WeightBlock(", size(block, 1), "×", size(block, 2),
        block.denom === nothing ? "" : ", denom", ")")

"""
    DirectPlan(method, missingpolicy, dst_space, src_space, block)

A plan holding one [`WeightBlock`](@ref) over the whole of both spaces.

The in-memory case: no chunk machinery, no discovery, no storage policy — the
source is read once and the block applied once. This is what bare
[`regrid`](@ref) builds for an in-memory source, and what
[`plan_regrid`](@ref) returns when the source is not chunked.

`block` spans both spaces entire: its size is
`(ncells(dst_space), ncells(src_space))` and its indices are cell positions,
the chunk-local addressing of a whole-domain pair. `chunks` and `budget` are
absent because a whole-domain plan has nothing to decide with them.
"""
struct DirectPlan{M<:AbstractRegriddingMethod,P<:AbstractMissingPolicy,
                  D<:RegridSpace,S<:RegridSpace,B<:WeightBlock} <: AbstractRegriddingPlan
    method::M
    missingpolicy::P
    dst_space::D
    src_space::S
    block::B
end

Base.show(io::IO, plan::DirectPlan) =
    print(io, "DirectPlan(", typeof(plan.method).name.name, ", ",
        size(plan.block, 1), " ← ", size(plan.block, 2), " cells)")

"""
    AbstractBlockStorage

Where a [`ChunkedPlan`](@ref) keeps the blocks it has built:
[`PerChunk`](@ref) in memory, [`Spilled`](@ref) on disk.
"""
abstract type AbstractBlockStorage end

"""
    CachedBlock(block, ref, bytes)

A built [`WeightBlock`](@ref) together with the data-independent quantities worth
keeping beside it: its reference vector (`blockreference!`) and its approximate
size in bytes. `used` is the storage's recency stamp.
"""
mutable struct CachedBlock
    block::WeightBlock
    ref::Vector{Float64}
    bytes::Int
    used::Int
end

# Approximate resident bytes of a block and its reference vector. Sparse storage
# is counted from its own arrays rather than by `summarysize`, which walks.
function _blockbytes(block::WeightBlock, ref::Vector{Float64})
    W = block.weights
    w = W isa SparseMatrixCSC ?
        16 * SparseArrays.nnz(W) + 8 * (size(W, 2) + 1) :
        8 * length(W)
    d = block.denom === nothing ? 0 : 8 * length(block.denom)
    return w + d + 8 * length(ref) + 64
end

"""
    PerChunk(capacity)
    PerChunk(; capacity = typemax(Int), maxbytes = typemax(Int))

Keep built [`WeightBlock`](@ref)s in memory, keyed by
`(destination chunk, source chunk)`, evicting least-recently-used blocks once
either bound is exceeded.

`capacity` bounds the number of resident blocks and `maxbytes` their approximate
total size; the most recently touched block is never evicted, so a bound of one
still applies. Eviction costs correctness nothing — an evicted block is rebuilt
from geometry on its next touch — so the bounds are a residency knob and never
change the answer.

Concurrent readers of **different** keys are safe: the dictionary is guarded by a
lock that is never held across a weight build. Two tasks racing on one key may
both build it, and the loser's block is discarded.
"""
mutable struct PerChunk <: AbstractBlockStorage
    capacity::Int
    maxbytes::Int
    blocks::Dict{Tuple{Int,Int},CachedBlock}
    bytes::Int
    clock::Int
    builds::Int
    lock::ReentrantLock
end

function PerChunk(; capacity::Integer = typemax(Int), maxbytes::Integer = typemax(Int))
    capacity >= 1 || throw(ArgumentError("PerChunk capacity must be at least one block, got $capacity"))
    maxbytes >= 1 || throw(ArgumentError("PerChunk maxbytes must be positive, got $maxbytes"))
    return PerChunk(Int(capacity), Int(maxbytes), Dict{Tuple{Int,Int},CachedBlock}(),
        0, 0, 0, ReentrantLock())
end

PerChunk(capacity::Integer) = PerChunk(; capacity)

Base.show(io::IO, s::PerChunk) =
    print(io, "PerChunk(", length(s.blocks), " blocks, ", s.bytes, " bytes)")

"""
    nblocks(storage) -> Int

How many [`WeightBlock`](@ref)s `storage` currently holds.
"""
nblocks(s::PerChunk) = length(s.blocks)

"""
    storagebytes(storage) -> Int

The approximate resident size in bytes of the blocks `storage` holds.
"""
storagebytes(s::PerChunk) = s.bytes

"""
    getblock!(storage, key, build) -> CachedBlock

The cached block for `key`, building it with `build()` on first touch.

`build` runs **outside** the storage lock, so a slow weight construction never
blocks a reader of another key; a duplicate build racing on one key is resolved
by keeping the block already inserted.
"""
function getblock!(storage::PerChunk, key::Tuple{Int,Int}, build::F) where {F}
    hit = _touch!(storage, key)
    hit === nothing || return hit
    block = build()
    ref = blockreference!(Vector{Float64}(undef, size(block, 1)), block)
    entry = CachedBlock(block, ref, _blockbytes(block, ref), 0)
    return _insert!(storage, key, entry)
end

function _touch!(storage::PerChunk, key::Tuple{Int,Int})
    @lock storage.lock begin
        entry = get(storage.blocks, key, nothing)
        entry === nothing && return nothing
        entry.used = (storage.clock += 1)
        return entry
    end
end

function _insert!(storage::PerChunk, key::Tuple{Int,Int}, entry::CachedBlock)
    @lock storage.lock begin
        existing = get(storage.blocks, key, nothing)
        if existing !== nothing
            existing.used = (storage.clock += 1)
            return existing
        end
        entry.used = (storage.clock += 1)
        storage.blocks[key] = entry
        storage.bytes += entry.bytes
        storage.builds += 1
        _evict!(storage, key)
        return entry
    end
end

# Drop least-recently-used entries until both bounds hold, never the one just
# touched. Called under the lock.
function _evict!(storage::PerChunk, keep::Tuple{Int,Int})
    while length(storage.blocks) > 1 &&
          (length(storage.blocks) > storage.capacity || storage.bytes > storage.maxbytes)
        victim = keep
        oldest = typemax(Int)
        for (k, e) in storage.blocks
            k == keep && continue
            if e.used < oldest
                oldest = e.used
                victim = k
            end
        end
        victim == keep && break
        storage.bytes -= storage.blocks[victim].bytes
        delete!(storage.blocks, victim)
    end
    return storage
end

"""
    Spilled(dir)

Serialize built [`WeightBlock`](@ref)s to scratch directory `dir` and reload
them on demand, so a plan whose weights exceed memory still applies.

The on-disk format is private to this package and is not a weight-exchange
format.
"""
struct Spilled <: AbstractBlockStorage
    dir::String
    Spilled(::AbstractString) = error("Spilled is not yet implemented")
end

"""
    ChunkedPlan(method, missingpolicy, dst_space, src_space, storage, budget, chunks)
    ChunkedPlan(method, missingpolicy, dst_space, src_space;
                storage = PerChunk(), budget = 2^30, chunks = nothing)

A plan whose weights are [`WeightBlock`](@ref)s keyed by
`(destination chunk, source chunk)` and built on first touch.

The lazy and streaming case: reading a destination chunk discovers its
connected source chunks, takes each pair's block from `storage`, and
accumulates. `budget` is a performance knob only — it decides whether the
connected source chunks are held together or streamed one at a time — and never
changes the answer. `chunks` is the chunking the destination array reports, as a
`DiskArrays.GridChunks`, a tuple of chunk sizes, or `nothing` to derive it from
the destination space's own chunks and the source's non-spatial chunking.

Constructing a plan builds no weights and reads no data.
"""
struct ChunkedPlan{M<:AbstractRegriddingMethod,P<:AbstractMissingPolicy,
                   D<:RegridSpace,S<:RegridSpace,T<:AbstractBlockStorage,C} <: AbstractRegriddingPlan
    method::M
    missingpolicy::P
    dst_space::D
    src_space::S
    storage::T
    budget::Int
    chunks::C
end

ChunkedPlan(method::AbstractRegriddingMethod, missingpolicy::AbstractMissingPolicy,
    dst_space::RegridSpace, src_space::RegridSpace;
    storage::AbstractBlockStorage = PerChunk(), budget::Integer = 2^30, chunks = nothing) =
    ChunkedPlan(method, missingpolicy, dst_space, src_space, storage, Int(budget), chunks)

Base.show(io::IO, plan::ChunkedPlan) =
    print(io, "ChunkedPlan(", typeof(plan.method).name.name, ", ",
        ncells(plan.dst_space), " cells / ", nchunks(plan.dst_space), " chunks ← ",
        ncells(plan.src_space), " cells / ", nchunks(plan.src_space), " chunks)")

"""
    blockfor(plan::ChunkedPlan, dstchunk::Integer, srcchunk::Integer) -> CachedBlock

The [`WeightBlock`](@ref) of one chunk pair, from the plan's storage, built on
first touch by one [`build_weights!`](@ref) call over the two chunks' cell
indices.

Building is geometry-only and reads no data, so a block is as cheap to rebuild
after eviction as it was to build.
"""
function blockfor(plan::ChunkedPlan, dstchunk::Integer, srcchunk::Integer)
    key = (Int(dstchunk), Int(srcchunk))
    return getblock!(plan.storage, key, () -> buildblock(plan, key[1], key[2]))
end

"""
    buildblock(plan::ChunkedPlan, dstchunk::Integer, srcchunk::Integer) -> WeightBlock

One chunk pair's weights, built unconditionally — the storage decides whether
this is called.
"""
function buildblock(plan::ChunkedPlan, dstchunk::Integer, srcchunk::Integer)
    dinds = cellindices(plan.dst_space, Int(dstchunk))
    sinds = cellindices(plan.src_space, Int(srcchunk))
    coo = WeightCOO(length(dinds))
    build_weights!(coo, plan.method, plan.dst_space, dinds, plan.src_space, sinds)
    return WeightBlock(coo, length(dinds), length(sinds))
end
