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
    DirectPlan(method, missingpolicy, dst_space, src_space, block, missingval = nothing)

A plan holding one [`WeightBlock`](@ref) over the whole of both spaces.

The in-memory case: no chunk machinery, no discovery, no storage policy — the
source is read once and the block applied once. This is what bare
[`regrid`](@ref) builds for an in-memory source, and what
[`plan_regrid`](@ref) returns when the source is not chunked.

`block` spans both spaces entire: its size is
`(ncells(dst_space), ncells(src_space))` and its indices are cell positions,
the chunk-local addressing of a whole-domain pair. `chunks` and `budget` are
absent because a whole-domain plan has nothing to decide with them.

`missingval` is the source's nodata sentinel, or `nothing` when it has none;
see [`isvalidvalue`](@ref).
"""
struct DirectPlan{M<:AbstractRegriddingMethod,P<:AbstractMissingPolicy,
                  D<:RegridSpace,S<:RegridSpace,B<:WeightBlock,V} <: AbstractRegriddingPlan
    method::M
    missingpolicy::P
    dst_space::D
    src_space::S
    block::B
    missingval::V
end

DirectPlan(method::AbstractRegriddingMethod, missingpolicy::AbstractMissingPolicy,
    dst_space::RegridSpace, src_space::RegridSpace, block::WeightBlock) =
    DirectPlan(method, missingpolicy, dst_space, src_space, block, nothing)

Base.show(io::IO, plan::DirectPlan) =
    print(io, "DirectPlan(", typeof(plan.method).name.name, ", ",
        size(plan.block, 1), " ← ", size(plan.block, 2), " cells)")

"""
    AbstractBlockStorage

Where a [`ChunkedPlan`](@ref) keeps the blocks it has built:
[`PerChunk`](@ref) in memory, [`Spilled`](@ref) on disk.
"""
abstract type AbstractBlockStorage end

# ===========================================================================
# The budget split
# ===========================================================================

# A plan's `budget` is its whole transient-memory target, and two consumers
# share it: the weights built while a destination chunk is produced, and the
# source data loaded to apply them. The split is fixed rather than adaptive —
# an adaptive one would make a plan's residency depend on the order chunks were
# touched in, which is exactly the property a budget exists to remove.
const WEIGHT_BUDGET_SHARE = 0.25

"""
    weightbudget(budget::Integer) -> Int

The bytes of a plan's `budget` given to built [`WeightBlock`](@ref)s — the
default `maxbytes` of the [`PerChunk`](@ref) a plan builds for itself.

A quarter of the budget: weights are rebuilt from geometry when evicted, while
source data must be re-read from storage, so the cheaper-to-lose consumer gets
the smaller share.
"""
weightbudget(budget::Integer) = max(1, floor(Int, WEIGHT_BUDGET_SHARE * budget))

"""
    databudget(budget::Integer) -> Int

The bytes of a plan's `budget` available for loaded source chunks — the rest of
it, after [`weightbudget`](@ref).

This is what decides hold-versus-stream while a destination chunk is produced:
a connected set that fits here may be held and reused, one that does not is
streamed load-apply-drop, one chunk resident at a time.
"""
databudget(budget::Integer) = max(1, Int(budget) - weightbudget(budget))

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
    Spilled(dir; capacity = typemax(Int), maxbytes = typemax(Int))

Serialize built [`WeightBlock`](@ref)s to scratch directory `dir` and reload
them on demand, so a plan whose weights exceed memory still applies.

A [`PerChunk`](@ref) cache sits in front, bounded by `capacity` and `maxbytes`,
so a block in memory is answered without touching the filesystem; a block only
on disk is deserialized and enters the cache; a block on neither is built,
written, and cached. Correctness never depends on which of the three happens.

Files are **not** shared between plans. Each `Spilled` stamps its own random tag
into every file name at construction, so two plans pointed at one directory
never read each other's weights — a plan's blocks are meaningless to a plan with
different spaces, a different method, or different chunking, and there is no
key that could tell those apart. Files are left behind when the plan is dropped;
`dir` is a scratch directory the caller owns.

The on-disk format is private, minimal, and **not durable**: raw CSC arrays with
a version byte, valid for this package version on this machine and nothing more.
It is not a weight-exchange format.
"""
struct Spilled <: AbstractBlockStorage
    dir::String
    tag::String
    memory::PerChunk
end

function Spilled(dir::AbstractString; capacity::Integer = typemax(Int),
    maxbytes::Integer = typemax(Int))
    d = String(dir)
    isdir(d) || mkpath(d)
    # Plan identity: a fresh tag per storage, so a directory reused by a second
    # plan is a cold cache rather than a wrong one.
    tag = string(rand(UInt64); base = 16, pad = 16)
    return Spilled(d, tag, PerChunk(; capacity, maxbytes))
end

Base.show(io::IO, s::Spilled) =
    print(io, "Spilled(", repr(s.dir), ", ", length(s.memory.blocks), " cached)")

"""
    blockpath(storage::Spilled, key) -> String

Where one chunk pair's block is spilled. The plan's tag is part of the name, so
this can never name another plan's file.
"""
blockpath(storage::Spilled, key::Tuple{Int,Int}) =
    joinpath(storage.dir, string("gr-", storage.tag, "-", key[1], "-", key[2], ".blk"))

"""
    spilledfiles(storage::Spilled) -> Vector{String}

The files this storage has written, by name. For inspection and tests; nothing
in the executor calls it.
"""
function spilledfiles(storage::Spilled)
    prefix = string("gr-", storage.tag, "-")
    isdir(storage.dir) || return String[]
    return sort!([f for f in readdir(storage.dir)
                  if startswith(f, prefix) && endswith(f, ".blk")])
end

nblocks(s::Spilled) = nblocks(s.memory)
storagebytes(s::Spilled) = storagebytes(s.memory)

function getblock!(storage::Spilled, key::Tuple{Int,Int}, build::F) where {F}
    return getblock!(storage.memory, key, function ()
        path = blockpath(storage, key)
        isfile(path) && return readblockfile(path)
        block = build()
        writeblockfile(path, block)
        return block
    end)
end

# The spill format. Four bytes of magic, a version byte, then the CSC arrays and
# the optional denominator, all little-endian native. `nothing` denom is a
# distinct flag rather than a zero-length vector, because a block without a
# denominator finalizes differently from one whose denominator is empty.
const SPILL_MAGIC = 0x42575247  # "GRWB"
const SPILL_VERSION = 0x01

"""
    writeblockfile(path, block::WeightBlock) -> path

Write `block` to `path` in the private spill format. Written to a temporary
name and moved into place, so a reader never sees a half-written block.
"""
function writeblockfile(path::AbstractString, block::WeightBlock)
    W = block.weights
    W isa SparseMatrixCSC || throw(ArgumentError(
        "Spilled serializes sparse weight blocks; this block's weights are a " *
        "$(typeof(W)). Use PerChunk storage for a method that builds dense blocks."))
    nz = SparseArrays.nnz(W)
    # Unique per writer, not merely per process: two tasks racing one key both
    # build and both write, and a shared temporary name would interleave them.
    tmp = string(path, ".", getpid(), ".", rand(UInt32), ".tmp")
    open(tmp, "w") do io
        write(io, SPILL_MAGIC)
        write(io, SPILL_VERSION)
        write(io, Int64(size(W, 1)), Int64(size(W, 2)), Int64(nz))
        write(io, Int64.(SparseArrays.getcolptr(W)))
        write(io, Int64.(view(SparseArrays.rowvals(W), 1:nz)))
        write(io, Float64.(view(SparseArrays.nonzeros(W), 1:nz)))
        d = block.denom
        write(io, UInt8(d === nothing ? 0 : 1))
        d === nothing || write(io, d)
    end
    mv(tmp, path; force = true)
    return path
end

"""
    readblockfile(path) -> WeightBlock

Read back a block [`writeblockfile`](@ref) wrote. The reference vector is not
stored — it is derived from the block on reload, which costs one pass over the
nonzeros and cannot disagree with the weights.
"""
function readblockfile(path::AbstractString)
    return open(path, "r") do io
        read(io, UInt32) == SPILL_MAGIC || throw(ArgumentError(
            "$path is not a GlobalRegridding weight spill"))
        read(io, UInt8) == SPILL_VERSION || throw(ArgumentError(
            "$path was written by a different version of the spill format; " *
            "the format is private and not durable, so delete the scratch directory"))
        m = Int(read(io, Int64))
        n = Int(read(io, Int64))
        nz = Int(read(io, Int64))
        colptr = Vector{Int}(read!(io, Vector{Int64}(undef, n + 1)))
        rowval = Vector{Int}(read!(io, Vector{Int64}(undef, nz)))
        nzval = read!(io, Vector{Float64}(undef, nz))
        W = SparseMatrixCSC{Float64,Int}(m, n, colptr, rowval, nzval)
        if read(io, UInt8) == 0x00
            return WeightBlock(W, nothing)
        end
        return WeightBlock(W, read!(io, Vector{Float64}(undef, m)))
    end
end

"""
    ChunkedPlan(method, missingpolicy, dst_space, src_space, storage, budget, chunks,
                missingval = nothing)
    ChunkedPlan(method, missingpolicy, dst_space, src_space;
                storage = nothing, budget = 2^30, chunks = nothing, missingval = nothing)

A plan whose weights are [`WeightBlock`](@ref)s keyed by
`(destination tile, source chunk)` and built on first touch.

The lazy and streaming case: reading a destination tile discovers its connected
source chunks, takes each pair's block from `storage`, and accumulates. A
destination tile is the destination space's own chunk whenever those chunks
partition the cell axis into ascending contiguous runs, and a run of cell
positions from `chunks` otherwise — see [`LazyRegridArray`](@ref).

`budget` is the plan's whole transient-memory target in bytes, and a performance
knob only: it bounds the resident weights ([`weightbudget`](@ref)) and decides
whether a destination tile's connected source chunks are held together or
streamed one at a time ([`databudget`](@ref)). Neither changes the answer.
`storage` defaults to a `PerChunk` bounded by the weight share of the budget;
pass one explicitly to override, including an unbounded `PerChunk()`.

`chunks` is the chunking the destination array reports — a
`DiskArrays.GridChunks`, a tuple of chunk sizes, or `nothing` to derive it — and
its cell axis is also the destination tiling when the destination space has no
usable one. `missingval` is the source's nodata sentinel; see
[`isvalidvalue`](@ref).

Constructing a plan builds no weights and reads no data.
"""
struct ChunkedPlan{M<:AbstractRegriddingMethod,P<:AbstractMissingPolicy,
                   D<:RegridSpace,S<:RegridSpace,T<:AbstractBlockStorage,C,V} <: AbstractRegriddingPlan
    method::M
    missingpolicy::P
    dst_space::D
    src_space::S
    storage::T
    budget::Int
    chunks::C
    missingval::V
end

ChunkedPlan(method::AbstractRegriddingMethod, missingpolicy::AbstractMissingPolicy,
    dst_space::RegridSpace, src_space::RegridSpace, storage::AbstractBlockStorage,
    budget::Integer, chunks) =
    ChunkedPlan(method, missingpolicy, dst_space, src_space, storage, Int(budget),
        chunks, nothing)

ChunkedPlan(method::AbstractRegriddingMethod, missingpolicy::AbstractMissingPolicy,
    dst_space::RegridSpace, src_space::RegridSpace;
    storage::Union{Nothing,AbstractBlockStorage} = nothing, budget::Integer = 2^30,
    chunks = nothing, missingval = nothing) =
    ChunkedPlan(method, missingpolicy, dst_space, src_space,
        storage === nothing ? PerChunk(; maxbytes = weightbudget(budget)) : storage,
        Int(budget), chunks, missingval)

Base.show(io::IO, plan::ChunkedPlan) =
    print(io, "ChunkedPlan(", typeof(plan.method).name.name, ", ",
        ncells(plan.dst_space), " cells / ", nchunks(plan.dst_space), " chunks ← ",
        ncells(plan.src_space), " cells / ", nchunks(plan.src_space), " chunks)")

"""
    blockfor(plan::ChunkedPlan, key::Tuple{Int,Int}, dinds) -> CachedBlock
    blockfor(plan::ChunkedPlan, dstchunk::Integer, srcchunk::Integer) -> CachedBlock

The [`WeightBlock`](@ref) of one pair, from the plan's storage, built on first
touch by one [`build_weights!`](@ref) call over the two sides' cell indices.

`key` is `(destination tile, source chunk)` and `dinds` the destination cells
that tile owns; the second form is the special case where the tile is the
destination space's chunk of the same number. The key is what the storage
addresses, so a plan's tiling must not change over its lifetime.

Building is geometry-only and reads no data, so a block is as cheap to rebuild
after eviction as it was to build.
"""
function blockfor(plan::ChunkedPlan, key::Tuple{Int,Int}, dinds)
    return getblock!(plan.storage, key,
        () -> buildblock(plan, dinds, cellindices(plan.src_space, key[2])))
end

blockfor(plan::ChunkedPlan, dstchunk::Integer, srcchunk::Integer) =
    blockfor(plan, (Int(dstchunk), Int(srcchunk)),
        cellindices(plan.dst_space, Int(dstchunk)))

"""
    buildblock(plan::ChunkedPlan, dinds, sinds) -> WeightBlock
    buildblock(plan::ChunkedPlan, dstchunk::Integer, srcchunk::Integer) -> WeightBlock

One pair's weights, built unconditionally — the storage decides whether this is
called.
"""
function buildblock(plan::ChunkedPlan, dinds, sinds)
    coo = WeightCOO(length(dinds))
    build_weights!(coo, plan.method, plan.dst_space, dinds, plan.src_space, sinds)
    return WeightBlock(coo, length(dinds), length(sinds))
end

buildblock(plan::ChunkedPlan, dstchunk::Integer, srcchunk::Integer) =
    buildblock(plan, cellindices(plan.dst_space, Int(dstchunk)),
        cellindices(plan.src_space, Int(srcchunk)))
