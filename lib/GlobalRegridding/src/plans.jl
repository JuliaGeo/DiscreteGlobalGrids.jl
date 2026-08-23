# Reusable plans and weight-block storage.

"""
    AbstractRegriddingPlan

A self-contained regridding operator with method, spaces, missing policy, and
weights. [`DirectPlan`](@ref) stores one whole-domain block;
[`ChunkedPlan`](@ref) builds blocks by chunk pair. Only weight storage mutates.
"""
abstract type AbstractRegriddingPlan end

"""
    WeightBlock(weights, denom)
    WeightBlock(coo::WeightCOO, ndst::Integer, nsrc::Integer)

Store one chunk pair's immutable weights. `weights` uses chunk-local indices;
`denom` contains optional per-destination denominators. Construction from
[`WeightCOO`](@ref) sums duplicate entries. Summing blocks across source chunks
reconstructs the full operator.
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

Return whether the block carries per-destination denominators.
"""
hasdenom(block::WeightBlock) = block.denom !== nothing

Base.show(io::IO, block::WeightBlock) =
    print(io, "WeightBlock(", size(block, 1), "×", size(block, 2),
        block.denom === nothing ? "" : ", denom", ")")

"""
    DirectPlan(method, missingpolicy, dst_space, src_space, block, missingval = nothing,
               sampling = nothing)

Store one [`WeightBlock`](@ref) over both complete spaces. The block size is
`(ncells(dst_space), ncells(src_space))`; its local indices are cell positions.
`missingval` is an optional source nodata sentinel. `sampling` overrides the
destination lookup sampling the method would imply ([`outputsampling`](@ref)).
"""
struct DirectPlan{M<:AbstractRegriddingMethod,P<:AbstractMissingPolicy,
                  D<:RegridSpace,S<:RegridSpace,B<:WeightBlock,V,
                  L<:Union{Nothing,DD.Lookups.Sampling}} <: AbstractRegriddingPlan
    method::M
    missingpolicy::P
    dst_space::D
    src_space::S
    block::B
    missingval::V
    sampling::L
end

DirectPlan(method::AbstractRegriddingMethod, missingpolicy::AbstractMissingPolicy,
    dst_space::RegridSpace, src_space::RegridSpace, block::WeightBlock,
    missingval = nothing) =
    DirectPlan(method, missingpolicy, dst_space, src_space, block, missingval, nothing)

Base.show(io::IO, plan::DirectPlan) =
    print(io, "DirectPlan(", typeof(plan.method).name.name, ", ",
        size(plan.block, 1), " ← ", size(plan.block, 2), " cells)")

"""
    AbstractBlockStorage

Where a [`ChunkedPlan`](@ref) keeps the blocks it has built:
[`PerChunk`](@ref) in memory, [`Spilled`](@ref) on disk.
"""
abstract type AbstractBlockStorage end

# Memory-budget split

# Fixed split between cached weights and loaded source data.
const WEIGHT_BUDGET_SHARE = 0.25

"""
    weightbudget(budget::Integer) -> Int

Return the part of `budget` assigned to resident weight blocks.
"""
weightbudget(budget::Integer) = max(1, floor(Int, WEIGHT_BUDGET_SHARE * budget))

"""
    databudget(budget::Integer) -> Int

Return the part of `budget` available for loaded source chunks.
"""
databudget(budget::Integer) = max(1, Int(budget) - weightbudget(budget))

"""
    CachedBlock(block, ref, bytes)

A weight block with its reference vector, approximate size, and recency stamp.
"""
mutable struct CachedBlock
    block::WeightBlock
    ref::Vector{Float64}
    bytes::Int
    used::Int
end

# Estimate resident bytes from the block's backing arrays.
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

Cache blocks by `(destination chunk, source chunk)` and evict least-recently-used
entries when `capacity` or `maxbytes` is exceeded. The newest block is retained.
Builds run outside the lock; duplicate concurrent builds keep the first result.
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

# The three observers below take the lock like every other reader of `blocks`.
# A wave builds concurrently, and `length(::Dict)` read while `_insert!` rehashes
# is the one way a caller merely *looking* at the cache can see a torn table.
Base.show(io::IO, s::PerChunk) =
    @lock s.lock print(io, "PerChunk(", length(s.blocks), " blocks, ", s.bytes, " bytes)")

"""
    nblocks(storage) -> Int

How many [`WeightBlock`](@ref)s `storage` currently holds.
"""
nblocks(s::PerChunk) = @lock s.lock length(s.blocks)

"""
    storagebytes(storage) -> Int

The approximate resident size in bytes of the blocks `storage` holds.
"""
storagebytes(s::PerChunk) = @lock s.lock s.bytes

"""
    getblock!(storage, key, build) -> CachedBlock

Return the cached block for `key`, calling `build()` on a miss. Builds run
outside the storage lock.
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

# Evict under the lock until both bounds hold; retain the latest entry.
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

Store blocks in scratch directory `dir` behind a [`PerChunk`](@ref) cache. The
caller owns `dir` and its lifetime. Filenames carry a per-instance tag, so only
this storage can read them: once the plan is dropped they are unreadable
garbage, and deleting the directory is the caller's job.
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
    # Unique storage tag prevents cross-plan cache hits.
    tag = string(rand(UInt64); base = 16, pad = 16)
    return Spilled(d, tag, PerChunk(; capacity, maxbytes))
end

Base.show(io::IO, s::Spilled) =
    print(io, "Spilled(", repr(s.dir), ", ", length(s.memory.blocks), " cached)")

"""
    blockpath(storage::Spilled, key) -> String

Return the spill path for one chunk pair.
"""
blockpath(storage::Spilled, key::Tuple{Int,Int}) =
    joinpath(storage.dir, string("gr-", storage.tag, "-", key[1], "-", key[2], ".blk"))

"""
    spilledfiles(storage::Spilled) -> Vector{String}

Return the names of files written by this storage.
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

# Private format: magic, version, CSC arrays, then an optional denominator.
const SPILL_MAGIC = 0x42575247  # "GRWB"
const SPILL_VERSION = 0x01

"""
    writeblockfile(path, block::WeightBlock) -> path

Atomically write `block` to `path` in the private spill format.
"""
function writeblockfile(path::AbstractString, block::WeightBlock)
    W = block.weights
    W isa SparseMatrixCSC || throw(ArgumentError(
        "Spilled serializes sparse weight blocks; this block's weights are a " *
        "$(typeof(W)). Use PerChunk storage for a method that builds dense blocks."))
    nz = SparseArrays.nnz(W)
    # Concurrent writers need distinct temporary files.
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

Read a block written by [`writeblockfile`](@ref). Reference vectors are rebuilt
from the weights.
"""
function readblockfile(path::AbstractString)
    return open(path, "r") do io
        read(io, UInt32) == SPILL_MAGIC || throw(ArgumentError(
            "$path is not a GlobalRegridding weight spill"))
        read(io, UInt8) == SPILL_VERSION || throw(ArgumentError(
            "$path was written by another version of the private spill format; " *
            "delete the scratch directory"))
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
                missingval = nothing, dependencies = nothing)
    ChunkedPlan(method, missingpolicy, dst_space, src_space;
                storage = nothing, budget = 2^30, chunks = nothing, missingval = nothing,
                dependencies = nothing, refine = nothing, narrow = nothing)

Build and store [`WeightBlock`](@ref)s by `(destination tile, source chunk)` on
first use. `budget` limits transient weight and source-data residency without
changing results. `chunks` sets lazy destination tiling; `nothing` derives it.
`missingval` is an optional source nodata sentinel. Construction reads no data
and builds no weights.

The plan is also the sole owner of its chunk dependency relation. `dependencies`
selects which one it holds, once, at construction; [`dependencies`](@ref) reads
it back and builds nothing. See that accessor for the whole contract, and
[`plan_regrid`](@ref) for the keywords as an API user meets them.
"""
struct ChunkedPlan{M<:AbstractRegriddingMethod,P<:AbstractMissingPolicy,
                   D<:RegridSpace,S<:RegridSpace,T<:AbstractBlockStorage,C,V,
                   G<:Union{Nothing,ChunkDependencyGraph}} <: AbstractRegriddingPlan
    method::M
    missingpolicy::P
    dst_space::D
    src_space::S
    storage::T
    budget::Int
    chunks::C
    missingval::V
    dependencies::G
end

ChunkedPlan(method::AbstractRegriddingMethod, missingpolicy::AbstractMissingPolicy,
    dst_space::RegridSpace, src_space::RegridSpace, storage::AbstractBlockStorage,
    budget::Integer, chunks) =
    ChunkedPlan(method, missingpolicy, dst_space, src_space, storage, Int(budget),
        chunks, nothing, nothing)

ChunkedPlan(method::AbstractRegriddingMethod, missingpolicy::AbstractMissingPolicy,
    dst_space::RegridSpace, src_space::RegridSpace, storage::AbstractBlockStorage,
    budget::Integer, chunks, missingval) =
    ChunkedPlan(method, missingpolicy, dst_space, src_space, storage, Int(budget),
        chunks, missingval, nothing)

ChunkedPlan(method::AbstractRegriddingMethod, missingpolicy::AbstractMissingPolicy,
    dst_space::RegridSpace, src_space::RegridSpace;
    storage::Union{Nothing,AbstractBlockStorage} = nothing, budget::Integer = 2^30,
    chunks = nothing, missingval = nothing,
    dependencies = nothing, refine = nothing, narrow = nothing) =
    ChunkedPlan(method, missingpolicy, dst_space, src_space,
        storage === nothing ? PerChunk(; maxbytes = weightbudget(budget)) : storage,
        Int(budget), chunks, missingval,
        _plandependencies(dependencies, refine, narrow, method, dst_space, src_space))

function Base.show(io::IO, plan::ChunkedPlan)
    print(io, "ChunkedPlan(", typeof(plan.method).name.name, ", ",
        ncells(plan.dst_space), " cells / ", nchunks(plan.dst_space), " chunks ← ",
        ncells(plan.src_space), " cells / ", nchunks(plan.src_space), " chunks")
    g = plan.dependencies
    g === nothing || print(io, ", ", Graphs.ne(g), " dependency edges")
    print(io, ")")
end

# --------------------------------------------------------------------------
# The one relation a plan owns
# --------------------------------------------------------------------------

"""
    dependencies(plan) -> Union{Nothing,ChunkDependencyGraph}

Return the chunk dependency relation `plan` owns, or `nothing` when it owns
none. **This accessor builds nothing**: it is a field read, it queries neither
space, and it returns the identical object every time it is called.

A `ChunkedPlan` holds *at most one* relation, fixed when the plan was
constructed. There is no way to obtain a second one from a plan, and no way to
apply a narrow phase to a plan that already exists — `refine` is a keyword of
[`plan_regrid`](@ref) (and of the [`ChunkedPlan`](@ref) constructor it forwards
to) and of nothing else, and [`chunk_dependency_graph`](@ref) has no `plan`
method. That is what makes "the plan's relation" a well-defined phrase: whoever
schedules, evicts, prefetches or validates against a plan is looking at one
object.

An eager [`DirectPlan`](@ref) never has one: it holds a single whole-domain
block, so there is no chunk pairing to describe.

# What the plan holds, and how it is chosen

The `dependencies` keyword on lazy `plan_regrid` / `ChunkedPlan`:

- `nothing` (the default) — **the plan owns no relation**, and construction does
  no graph work at all. `dependencies(plan)` returns `nothing`. This is the
  right choice for a plan whose destination is one chunk, or a few: the lazy
  executor discovers a tile's sources by querying the source index directly, and
  a whole graph object for a handful of rows costs more than the queries it
  would replace, because [`restrict`](@ref) and a rebuild both pay
  `O(nsourcechunks)` for the source-major direction that a one-row view barely
  uses ([measured]: `regrid-notes/2026-08-23-g4-plan-owns-graph.md`).
- `true` — build it now, once, from the plan's own spaces at the plan's own
  radius, [`support_radius`](@ref)`(method, src_space)`. Passing `refine` or
  `narrow` implies this.
- a [`ChunkDependencyGraph`](@ref) — adopt a relation somebody else built, after
  [`validate_dependencies`](@ref) certifies it against *these* spaces, *this*
  radius and the `narrow` phase the caller claims it carries. An invalid reuse
  therefore fails at plan construction rather than as a wrong answer later.
  `refine` cannot be combined with a supplied graph: a narrow phase applies as a
  relation is built, and reapplying it afterwards would renumber nothing and
  prove nothing. Name the phase the graph already carries with `narrow` instead.
- `false` — hold none, explicitly. Only differs from `nothing` in that it
  rejects `refine`/`narrow` rather than acting on them.

Whichever branch runs, it runs **once**, at construction, and reads no source
data, builds no weights and issues no network metadata request.

# Example

```julia
global_plan = plan_regrid(data; to = dst, from = src, lazy = true,
                          dependencies = true)
g = dependencies(global_plan)
sourcesof(g, 12)                  # what destination chunk 12 may read
consumerdegree(g, 7)              # initial refcount for source chunk 7

# A second plan over the same pair reuses the relation instead of rebuilding it.
other = plan_regrid(data; to = dst, from = src, lazy = true, dependencies = g)
dependencies(other) === g
```
"""
dependencies(plan::ChunkedPlan) = plan.dependencies
dependencies(::DirectPlan) = nothing

# The default: no relation asked for, so no relation work. This method exists
# so that constructing a per-destination-chunk plan — which production does
# 66 175 times — pays nothing at all, not even `support_radius`.
_plandependencies(::Nothing, ::Nothing, ::Nothing, ::AbstractRegriddingMethod,
    ::RegridSpace, ::RegridSpace) = nothing

function _plandependencies(dependencies, refine, narrow,
        method::AbstractRegriddingMethod, dst_space::RegridSpace,
        src_space::RegridSpace)
    if dependencies isa ChunkDependencyGraph
        refine === nothing || throw(ArgumentError(
            "`refine` narrows a relation while it is being built, so it cannot " *
            "be applied to a `dependencies` graph the plan did not build. Pass " *
            "`narrow` to name the narrow phase the supplied graph already " *
            "carries, or drop `dependencies` to build a narrowed one here."))
        return validate_dependencies(dependencies, dst_space, src_space;
            radius = Float64(support_radius(method, src_space)),
            narrow = narrow === nothing ? :none : narrow)
    elseif dependencies === true ||
           (dependencies === nothing && (refine !== nothing || narrow !== nothing))
        return _builddependencies(dst_space, src_space,
            support_radius(method, src_space), refine, narrow)
    elseif dependencies === false
        (refine === nothing && narrow === nothing) || throw(ArgumentError(
            "`dependencies = false` asks the plan to hold no relation, but " *
            "`refine`/`narrow` describe one it would have to build; pass " *
            "`dependencies = true` or drop them."))
        return nothing
    elseif dependencies === nothing
        return nothing
    end
    throw(ArgumentError(
        "`dependencies` must be `nothing`, `true`, `false`, or a " *
        "ChunkDependencyGraph to adopt, got $(typeof(dependencies))"))
end

"""
    blockfor(plan::ChunkedPlan, key::Tuple{Int,Int}, dinds[, dst_space]) -> CachedBlock
    blockfor(plan::ChunkedPlan, dstchunk::Integer, srcchunk::Integer) -> CachedBlock

Return one cached chunk-pair block, building it on first use. `key` identifies
the destination tile and source chunk; `dinds` lists the tile's destination
cells. A shared `dst_space` can reuse tile geometry across several pairs.
"""
blockfor(plan::ChunkedPlan, key::Tuple{Int,Int}, dinds) =
    blockfor(plan, key, dinds, TileCells(plan.dst_space, dinds))

function blockfor(plan::ChunkedPlan, key::Tuple{Int,Int}, dinds, dst_space::RegridSpace)
    return getblock!(plan.storage, key,
        () -> buildblock(plan, dinds, cellindices(plan.src_space, key[2]), dst_space))
end

blockfor(plan::ChunkedPlan, dstchunk::Integer, srcchunk::Integer) =
    blockfor(plan, (Int(dstchunk), Int(srcchunk)),
        cellindices(plan.dst_space, Int(dstchunk)))

"""
    buildblock(plan::ChunkedPlan, dinds, sinds[, dst_space]) -> WeightBlock
    buildblock(plan::ChunkedPlan, dstchunk::Integer, srcchunk::Integer) -> WeightBlock

Build one chunk pair's weights without consulting storage.
"""
buildblock(plan::ChunkedPlan, dinds, sinds) =
    buildblock(plan, dinds, sinds, TileCells(plan.dst_space, dinds))

function buildblock(plan::ChunkedPlan, dinds, sinds, dst_space::RegridSpace)
    coo = WeightCOO(length(dinds))
    build_weights!(coo, plan.method, dst_space, dinds, plan.src_space, sinds)
    return WeightBlock(coo, length(dinds), length(sinds))
end

buildblock(plan::ChunkedPlan, dstchunk::Integer, srcchunk::Integer) =
    buildblock(plan, cellindices(plan.dst_space, Int(dstchunk)),
        cellindices(plan.src_space, Int(srcchunk)))
