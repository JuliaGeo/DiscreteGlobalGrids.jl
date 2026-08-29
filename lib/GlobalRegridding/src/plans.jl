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
    WeightBlock(weights, denom, coverage)
    WeightBlock(coo::WeightCOO, ndst::Integer, nsrc::Integer)

Store one chunk pair's immutable weights.

  - `weights` uses chunk-local indices; `denom` contains optional
    per-destination denominators. Summing blocks across source chunks
    reconstructs the full operator.
  - `weights` may be signed: a method whose stencil corrects as well as averages
    subtracts from some sources.
  - `coverage` is the non-negative per-(destination, source) weight a valid
    source contributes to its destination's coverage, over the same indices as
    `weights`, and `nothing` when `weights` serves as its own coverage — which
    it does exactly when every entry is non-negative. Execution reads it only
    where a validity mask says some sources are missing.
  - Construction from [`WeightCOO`](@ref) sums duplicate entries, and assembles
    the coverage matrix when the COO carries a coverage list of its own.
  - `reference` is the per-destination weight the values are normalized against,
    held once and reusable by every application: `denom` itself on a denominated
    block (`reference === denom`), the row sums of `weights` otherwise, computed
    when the block is built.
  - Nothing outside the block recomputes or copies the reference.
"""
struct WeightBlock{M<:AbstractMatrix{Float64},D<:Union{Nothing,Vector{Float64}},
                   C<:Union{Nothing,AbstractMatrix{Float64}}}
    weights::M
    denom::D
    coverage::C
    reference::Vector{Float64}
end

WeightBlock(weights::AbstractMatrix{Float64}, denom::Vector{Float64},
    coverage::Union{Nothing,AbstractMatrix{Float64}} = nothing) =
    WeightBlock{typeof(weights),Vector{Float64},typeof(coverage)}(
        weights, denom, _checkcoverage(weights, coverage), denom)

WeightBlock(weights::AbstractMatrix{Float64}, ::Nothing,
    coverage::Union{Nothing,AbstractMatrix{Float64}} = nothing) =
    WeightBlock{typeof(weights),Nothing,typeof(coverage)}(
        weights, nothing, _checkcoverage(weights, coverage), _rowsums(weights))

function WeightBlock(coo::WeightCOO, ndst::Integer, nsrc::Integer)
    weights = sparse(coo.rows, coo.cols, coo.vals, Int(ndst), Int(nsrc))
    d = coo.denom
    c = coo.coverage
    coverage = c === nothing ? nothing :
               sparse(c.rows, c.cols, c.vals, Int(ndst), Int(nsrc))
    return WeightBlock(weights, d === nothing ? nothing : copy(d), coverage)
end

# Coverage is indexed exactly like the values it stands in for, so a mismatch is
# caught where the block is built rather than in an `@inbounds` accumulation.
_checkcoverage(::AbstractMatrix{Float64}, ::Nothing) = nothing

function _checkcoverage(weights::AbstractMatrix{Float64},
    coverage::AbstractMatrix{Float64})
    size(coverage) == size(weights) || throw(DimensionMismatch(
        "coverage is $(size(coverage, 1))×$(size(coverage, 2)) over " *
        "$(size(weights, 1))×$(size(weights, 2)) weights"))
    return coverage
end

Base.size(block::WeightBlock) = size(block.weights)
Base.size(block::WeightBlock, d::Integer) = size(block.weights, d)

Base.show(io::IO, block::WeightBlock) =
    print(io, "WeightBlock(", size(block, 1), "×", size(block, 2),
        block.denom === nothing ? "" : ", denom",
        block.coverage === nothing ? "" : ", coverage", ")")

# A block with no denominator references its row sums, which is what the
# accumulated weight of every source cell a destination draws on comes to.
_rowsums(W::AbstractMatrix) = _addrowsums!(zeros(Float64, size(W, 1)), W)

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
    DirectPlan(method, missingpolicy, dst_space, src_space, block, missingval = nothing,
               sampling = nothing)

Store one [`WeightBlock`](@ref) over both complete spaces. The block size is
`(ncells(dst_space), ncells(src_space))`, so its chunk-local indices are the
spaces' own local indices.
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

# Bytes a chunked plan may hold in transient weights and source data when the
# caller names no budget; `plan_regrid` resolves its API keyword against this.
const DEFAULT_BUDGET = 2^31

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
    CachedBlock(block, bytes, used)

A weight block with its approximate size and recency stamp. The block carries
its own reference vector, so caching one copies no numerical state.
"""
mutable struct CachedBlock
    block::WeightBlock
    bytes::Int
    used::Int
end

CachedBlock(block::WeightBlock) = CachedBlock(block, _blockbytes(block), 0)

# Estimate resident bytes from the block's backing arrays. A denominated block's
# reference is its denominator, so either way one vector is counted; a block
# that separates coverage from values holds a second matrix.
function _blockbytes(block::WeightBlock)
    w = _matrixbytes(block.weights)
    c = block.coverage === nothing ? 0 : _matrixbytes(block.coverage)
    d = block.denom === nothing ? 0 : 8 * length(block.denom)
    r = block.reference === block.denom ? 0 : 8 * length(block.reference)
    return w + c + d + r + 64
end

_matrixbytes(W::SparseMatrixCSC) =
    16 * SparseArrays.nnz(W) + 8 * (size(W, 2) + 1)
_matrixbytes(W::AbstractMatrix) = 8 * length(W)

"""
    TileWeights(sourcechunks::Vector{Int}, blocks::Vector{WeightBlock})

One destination tile's weights, split by the source chunk each entry belongs to.

  - `sourcechunks` is the exact ascending list of the source chunks the tile's
    stencils name; `blocks[k]` holds the weights `sourcechunks[k]` supplies.
  - A block's rows are chunk-local destination indices within the tile and its
    columns chunk-local source indices within that chunk, the convention
    [`WeightCOO`](@ref) documents.
  - A chunk no stencil names has no block and does not appear, so the manifest
    bounds the tile's reads exactly.
"""
struct TileWeights
    sourcechunks::Vector{Int}
    blocks::Vector{WeightBlock}
end

Base.show(io::IO, tw::TileWeights) =
    print(io, "TileWeights(", length(tw.sourcechunks), " source chunks, ",
        sum(SparseArrays.nnz, (b.weights for b in tw.blocks); init = 0), " entries)")

"""
    CachedTile(weights::TileWeights)

One destination tile's cached weights: its source-chunk manifest and one
[`CachedBlock`](@ref) per chunk in the same order. `bytes` counts the whole tile,
manifest included, so a budget bounds a tile as it bounds a chunk pair.
"""
mutable struct CachedTile
    sourcechunks::Vector{Int}
    entries::Vector{CachedBlock}
    bytes::Int
    used::Int
end

function CachedTile(weights::TileWeights)
    entries = Vector{CachedBlock}(undef, length(weights.blocks))
    bytes = 8 * length(weights.sourcechunks) + 64
    for k in eachindex(weights.blocks)
        entry = CachedBlock(weights.blocks[k])
        entries[k] = entry
        bytes += entry.bytes
    end
    return CachedTile(copy(weights.sourcechunks), entries, bytes, 0)
end

Base.show(io::IO, tile::CachedTile) =
    print(io, "CachedTile(", length(tile.sourcechunks), " source chunks, ",
        tile.bytes, " bytes)")

"""
    tileblock(tile::CachedTile, chunk::Integer) -> Union{Nothing,CachedBlock}

The tile's cached block for source `chunk`, or `nothing` when no stencil in the
tile names it.
"""
function tileblock(tile::CachedTile, chunk::Integer)
    s = Int(chunk)
    k = searchsortedfirst(tile.sourcechunks, s)
    (k <= length(tile.sourcechunks) && @inbounds(tile.sourcechunks[k]) == s) ||
        return nothing
    return @inbounds tile.entries[k]
end

"""
    PerChunk(capacity)
    PerChunk(; capacity = typemax(Int), maxbytes = typemax(Int))

Cache weight blocks, evicting least-recently-used entries when `capacity` or
`maxbytes` is exceeded.

  - A chunk pair is keyed by `(destination chunk, source chunk)`; a point method
    whose build unit is a destination tile is keyed by tile number instead,
    through [`gettile!`](@ref).
  - The newest entry is retained. Tiles and chunk pairs share the recency clock,
    the entry count and the byte budget, and evict each other by recency alone.
  - Builds run outside the lock. Duplicate concurrent builds of a pair keep the
    first result; one tile is built once, a request meeting a build already in
    flight waiting for it.
"""
mutable struct PerChunk <: AbstractBlockStorage
    capacity::Int
    maxbytes::Int
    blocks::Dict{Tuple{Int,Int},CachedBlock}
    tiles::Dict{Int,CachedTile}
    building::Set{Int}
    bytes::Int
    clock::Int
    builds::Int
    lock::ReentrantLock
    ready::Threads.Condition
end

function PerChunk(; capacity::Integer = typemax(Int), maxbytes::Integer = typemax(Int))
    capacity >= 1 || throw(ArgumentError("PerChunk capacity must be at least one block, got $capacity"))
    maxbytes >= 1 || throw(ArgumentError("PerChunk maxbytes must be positive, got $maxbytes"))
    guard = ReentrantLock()
    return PerChunk(Int(capacity), Int(maxbytes), Dict{Tuple{Int,Int},CachedBlock}(),
        Dict{Int,CachedTile}(), Set{Int}(), 0, 0, 0, guard, Threads.Condition(guard))
end

PerChunk(capacity::Integer) = PerChunk(; capacity)

# The three observers below take the lock like every other reader of `blocks`: a
# wave builds concurrently, and an unlocked `length(::Dict)` can read a table
# `_insert!` is rehashing.
Base.show(io::IO, s::PerChunk) =
    @lock s.lock print(io, "PerChunk(", _entrycount(s), " blocks, ", s.bytes, " bytes)")

# Chunk pairs and destination tiles are both cache entries, and both bounds
# count them together.
_entrycount(s::PerChunk) = length(s.blocks) + length(s.tiles)

"""
    nblocks(storage) -> Int

How many entries `storage` currently holds: one per chunk-pair
[`WeightBlock`](@ref), and one per destination tile whose weights it keeps.
"""
nblocks(s::PerChunk) = @lock s.lock _entrycount(s)

"""
    storagebytes(storage) -> Int

The approximate resident size in bytes of the blocks `storage` holds.
"""
storagebytes(s::PerChunk) = @lock s.lock s.bytes

"""
    weightlimit(storage) -> Int

The bytes of weights `storage` holds before it evicts, or zero where it states
no bound. This is what bounds a lazy read's [`TilePrefetch`](@ref) queue.
"""
weightlimit(::AbstractBlockStorage) = 0

weightlimit(s::PerChunk) = s.maxbytes

"""
    getblock!(storage, key, build) -> CachedBlock

Return the cached block for `key`, calling `build()` on a miss. Builds run
outside the storage lock.
"""
function getblock!(storage::PerChunk, key::Tuple{Int,Int}, build::F) where {F}
    hit = _touch!(storage, key)
    hit === nothing || return hit
    return _insert!(storage, key, CachedBlock(build()))
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

# Evict under the lock until both bounds hold, retaining `keep`, the key just
# inserted: a chunk pair or a tile number, the two competing by recency alone.
function _evict!(storage::PerChunk, keep)
    while _entrycount(storage) > 1 &&
          (_entrycount(storage) > storage.capacity || storage.bytes > storage.maxbytes)
        victim = nothing
        oldest = typemax(Int)
        for (k, e) in storage.blocks
            k == keep && continue
            if e.used < oldest
                oldest = e.used
                victim = k
            end
        end
        for (k, e) in storage.tiles
            k == keep && continue
            if e.used < oldest
                oldest = e.used
                victim = k
            end
        end
        victim === nothing && break
        _drop!(storage, victim)
    end
    return storage
end

function _drop!(storage::PerChunk, key::Tuple{Int,Int})
    storage.bytes -= storage.blocks[key].bytes
    delete!(storage.blocks, key)
    return storage
end

function _drop!(storage::PerChunk, tile::Int)
    storage.bytes -= storage.tiles[tile].bytes
    delete!(storage.tiles, tile)
    return storage
end

"""
    gettile!(storage, tile::Int, build) -> CachedTile

Return the cached weights of destination `tile`, calling `build()` on a miss.

  - One tile is built once: a request meeting a build already in flight waits
    for it instead of starting a second.
  - The build itself runs outside the storage lock.
"""
function gettile!(storage::PerChunk, tile::Int, build::F) where {F}
    @lock storage.lock begin
        while true
            entry = get(storage.tiles, tile, nothing)
            if entry !== nothing
                entry.used = (storage.clock += 1)
                return entry
            end
            in(tile, storage.building) || break
            wait(storage.ready)
        end
        push!(storage.building, tile)
    end
    local built::CachedTile
    try
        built = CachedTile(build())
    catch
        _releasetile!(storage, tile)
        rethrow()
    end
    return _inserttile!(storage, tile, built)
end

# Give up a build claim nothing will finish, and wake whoever waited on it.
function _releasetile!(storage::PerChunk, tile::Int)
    @lock storage.lock begin
        delete!(storage.building, tile)
        notify(storage.ready)
    end
    return storage
end

function _inserttile!(storage::PerChunk, tile::Int, entry::CachedTile)
    @lock storage.lock begin
        delete!(storage.building, tile)
        notify(storage.ready)
        existing = get(storage.tiles, tile, nothing)
        if existing !== nothing
            existing.used = (storage.clock += 1)
            return existing
        end
        entry.used = (storage.clock += 1)
        storage.tiles[tile] = entry
        storage.bytes += entry.bytes
        storage.builds += 1
        _evict!(storage, tile)
        return entry
    end
end

"""
    Spilled(dir; capacity = typemax(Int), maxbytes = typemax(Int))

Store blocks in scratch directory `dir` behind a [`PerChunk`](@ref) cache.

  - The caller owns `dir` and its lifetime, deleting it included.
  - Filenames carry a per-instance tag, so only this storage can read them: once
    the plan is dropped they are unreadable garbage.
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
    tilepath(storage::Spilled, tile::Integer) -> String

Return the spill path for one destination tile's weights. It is a file of its
own, beside the chunk-pair files and in the same private format.
"""
tilepath(storage::Spilled, tile::Integer) =
    joinpath(storage.dir, string("gr-", storage.tag, "-t", Int(tile), ".tile"))

"""
    spilledfiles(storage::Spilled) -> Vector{String}

Return the names of files written by this storage.
"""
function spilledfiles(storage::Spilled)
    prefix = string("gr-", storage.tag, "-")
    isdir(storage.dir) || return String[]
    return sort!([f for f in readdir(storage.dir)
                  if startswith(f, prefix) &&
                     (endswith(f, ".blk") || endswith(f, ".tile"))])
end

nblocks(s::Spilled) = nblocks(s.memory)
storagebytes(s::Spilled) = storagebytes(s.memory)
weightlimit(s::Spilled) = weightlimit(s.memory)

function getblock!(storage::Spilled, key::Tuple{Int,Int}, build::F) where {F}
    return getblock!(storage.memory, key, function ()
        path = blockpath(storage, key)
        isfile(path) && return readblockfile(path)
        block = build()
        writeblockfile(path, block)
        return block
    end)
end

# A spilled tile keeps its manifest beside its blocks, so an evicted dependency
# list comes back off disk rather than from another destination pass.
function gettile!(storage::Spilled, tile::Int, build::F) where {F}
    return gettile!(storage.memory, tile, function ()
        path = tilepath(storage, tile)
        isfile(path) && return readtilefile(path)
        weights = build()
        writetilefile(path, weights)
        return weights
    end)
end

# Private format: magic, version, CSC arrays, then an optional denominator and
# an optional coverage matrix in the same CSC layout as the weights.
const SPILL_MAGIC = 0x42575247  # "GRWB"
const TILE_MAGIC = 0x54575247  # "GRWT"
const SPILL_VERSION = 0x02

"""
    writeblockfile(path, block::WeightBlock) -> path

Atomically write `block` to `path` in the private spill format.
"""
function writeblockfile(path::AbstractString, block::WeightBlock)
    return _atomicwrite(path) do io
        write(io, SPILL_MAGIC)
        write(io, SPILL_VERSION)
        _writeblock(io, block)
    end
end

"""
    readblockfile(path) -> WeightBlock

Read a block written by [`writeblockfile`](@ref). The file holds the weights,
the optional denominator, and the optional coverage matrix; the block's
reference vector is not stored but reconstructed here, as the denominator itself
or as the weights' row sums.
"""
function readblockfile(path::AbstractString)
    return open(path, "r") do io
        _readheader(io, SPILL_MAGIC, path)
        return _readblock(io)
    end
end

"""
    writetilefile(path, weights::TileWeights) -> path

Atomically write one destination tile's weights to `path`: the source-chunk
manifest, then one block per chunk in the same order.
"""
function writetilefile(path::AbstractString, weights::TileWeights)
    return _atomicwrite(path) do io
        write(io, TILE_MAGIC)
        write(io, SPILL_VERSION)
        write(io, Int64(length(weights.sourcechunks)))
        write(io, Int64.(weights.sourcechunks))
        for block in weights.blocks
            _writeblock(io, block)
        end
    end
end

"""
    readtilefile(path) -> TileWeights

Read a tile written by [`writetilefile`](@ref), manifest and all, without
rebuilding a stencil. Each block's reference vector is reconstructed as it is
read, the same way [`readblockfile`](@ref) does.
"""
function readtilefile(path::AbstractString)
    return open(path, "r") do io
        _readheader(io, TILE_MAGIC, path)
        n = Int(read(io, Int64))
        chunks = Vector{Int}(read!(io, Vector{Int64}(undef, n)))
        blocks = [_readblock(io) for _ in 1:n]
        return TileWeights(chunks, blocks)
    end
end

# Write through a uniquely named temporary: concurrent writers must not share
# one, and the reader must never see a partial file.
function _atomicwrite(body::F, path::AbstractString) where {F}
    tmp = string(path, ".", getpid(), ".", rand(UInt32), ".tmp")
    open(tmp, "w") do io
        body(io)
    end
    mv(tmp, path; force = true)
    return path
end

function _readheader(io::IO, magic::UInt32, path::AbstractString)
    read(io, UInt32) == magic || throw(ArgumentError(
        "$path is not a GlobalRegridding weight spill"))
    read(io, UInt8) == SPILL_VERSION || throw(ArgumentError(
        "$path was written by another version of the private spill format; " *
        "delete the scratch directory"))
    return io
end

function _writeblock(io::IO, block::WeightBlock)
    _writesparse(io, block.weights)
    d = block.denom
    write(io, UInt8(d === nothing ? 0 : 1))
    d === nothing || write(io, d)
    c = block.coverage
    write(io, UInt8(c === nothing ? 0 : 1))
    c === nothing || _writesparse(io, c)
    return io
end

function _readblock(io::IO)
    W = _readsparse(io)
    d = read(io, UInt8) == 0x00 ? nothing :
        read!(io, Vector{Float64}(undef, size(W, 1)))
    c = read(io, UInt8) == 0x00 ? nothing : _readsparse(io)
    return WeightBlock(W, d, c)
end

function _writesparse(io::IO, W::AbstractMatrix)
    throw(ArgumentError(
        "Spilled serializes sparse weight blocks; this block holds a " *
        "$(typeof(W)). Use PerChunk storage for a method that builds dense blocks."))
end

function _writesparse(io::IO, W::SparseMatrixCSC)
    nz = SparseArrays.nnz(W)
    write(io, Int64(size(W, 1)), Int64(size(W, 2)), Int64(nz))
    write(io, Int64.(SparseArrays.getcolptr(W)))
    write(io, Int64.(view(SparseArrays.rowvals(W), 1:nz)))
    write(io, Float64.(view(SparseArrays.nonzeros(W), 1:nz)))
    return io
end

function _readsparse(io::IO)
    m = Int(read(io, Int64))
    n = Int(read(io, Int64))
    nz = Int(read(io, Int64))
    colptr = Vector{Int}(read!(io, Vector{Int64}(undef, n + 1)))
    rowval = Vector{Int}(read!(io, Vector{Int64}(undef, nz)))
    nzval = read!(io, Vector{Float64}(undef, nz))
    return SparseMatrixCSC{Float64,Int}(m, n, colptr, rowval, nzval)
end

"""
    ChunkedPlan(method, missingpolicy, dst_space, src_space, storage, budget, chunks,
                missingval = nothing, dependencies = nothing)
    ChunkedPlan(method, missingpolicy, dst_space, src_space;
                storage = nothing, budget = 2^31, chunks = nothing, missingval = nothing,
                dependencies = nothing, refine = nothing, narrow = nothing)

Build and store [`WeightBlock`](@ref)s by `(destination tile, source chunk)` on
first use. Construction reads no data and builds no weights.

  - `budget` limits transient weight and source-data residency without changing
    results.
  - `chunks` sets lazy destination tiling; `nothing` derives it. `missingval` is
    an optional source nodata sentinel.
  - The plan is the sole owner of its chunk dependency relation, which
    `dependencies` selects once, at construction. See [`dependencies`](@ref) for
    the whole contract and [`plan_regrid`](@ref) for the keywords as an API user
    meets them.
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

# The positional forms build the relation too; the nine-argument field
# constructor is the one way to assemble a plan around an existing relation, or
# around none.
ChunkedPlan(method::AbstractRegriddingMethod, missingpolicy::AbstractMissingPolicy,
    dst_space::RegridSpace, src_space::RegridSpace, storage::AbstractBlockStorage,
    budget::Integer, chunks) =
    ChunkedPlan(method, missingpolicy, dst_space, src_space, storage, Int(budget),
        chunks, nothing,
        _plandependencies(nothing, nothing, nothing, method, dst_space, src_space))

ChunkedPlan(method::AbstractRegriddingMethod, missingpolicy::AbstractMissingPolicy,
    dst_space::RegridSpace, src_space::RegridSpace, storage::AbstractBlockStorage,
    budget::Integer, chunks, missingval) =
    ChunkedPlan(method, missingpolicy, dst_space, src_space, storage, Int(budget),
        chunks, missingval,
        _plandependencies(nothing, nothing, nothing, method, dst_space, src_space))

ChunkedPlan(method::AbstractRegriddingMethod, missingpolicy::AbstractMissingPolicy,
    dst_space::RegridSpace, src_space::RegridSpace;
    storage::Union{Nothing,AbstractBlockStorage} = nothing,
    budget::Integer = DEFAULT_BUDGET, chunks = nothing, missingval = nothing,
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

A `ChunkedPlan` holds *at most one* relation, fixed at construction and never
replaced or narrowed afterwards: `refine` is a keyword of [`plan_regrid`](@ref)
and the [`ChunkedPlan`](@ref) constructor it forwards to and of nothing else,
and [`chunk_dependency_graph`](@ref) has no `plan` method. An eager
[`DirectPlan`](@ref), holding one whole-domain block, never has one. For a
chunk-pair build unit a row *is* the read; for a point method with a
[`sampler`](@ref) the rows only have to stay a superset of the tile's
[`TileWeights`](@ref) manifest, which is what [`supportradius`](@ref) is for.

The `dependencies` keyword on lazy `plan_regrid` / `ChunkedPlan` selects which
one is held. Whichever branch runs, it runs **once**, at construction, and reads
no source data, builds no weights and issues no network metadata request.

  - `nothing` (the default), `true`, or any use of `refine`/`narrow`: build it
    from the plan's own spaces at the plan's own radius,
    [`supportradius`](@ref)`(method, src_space)`. This is the default because
    every lazy read needs a relation.
  - a [`ChunkDependencyGraph`](@ref): adopt one somebody else built, once
    [`validate_dependencies`](@ref) certifies it against *these* spaces, *this*
    radius and the `narrow` phase the caller claims it carries, so an invalid
    reuse fails at construction rather than as a wrong answer later. `refine` is
    refused here — name the phase the graph already carries with `narrow`
    instead. A plan over part of a bigger destination adopts the bigger relation
    through [`subspace_dependencies`](@ref).
  - `false`: hold none, explicitly, and reject `refine`/`narrow` rather than
    acting on them. Such a plan cannot back a [`LazyRegridArray`](@ref), which
    needs the rows to read; it is for a caller that wants nothing but
    [`blockfor`](@ref) and the plan's spaces.

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
            radius = Float64(supportradius(method, src_space)),
            narrow = narrow === nothing ? :none : narrow)
    elseif dependencies === true || dependencies === nothing
        return _builddependencies(dst_space, src_space,
            supportradius(method, src_space), refine, narrow)
    elseif dependencies === false
        (refine === nothing && narrow === nothing) || throw(ArgumentError(
            "`dependencies = false` asks the plan to hold no relation, but " *
            "`refine`/`narrow` describe one it would have to build; pass " *
            "`dependencies = true` or drop them."))
        return nothing
    end
    throw(ArgumentError(
        "`dependencies` must be `nothing`, `true`, `false`, or a " *
        "ChunkDependencyGraph to adopt, got $(typeof(dependencies))"))
end

"""
    blockfor(plan::ChunkedPlan, key::Tuple{Int,Int}, dinds[, destination]) -> CachedBlock
    blockfor(plan::ChunkedPlan, dstchunk::Integer, srcchunk::Integer) -> CachedBlock

Return the cached block destination tile `key[1]` takes from source chunk
`key[2]`, building what it comes from on first use.

  - `dinds` lists the tile's destination cells; a shared `destination`, what
    [`preparedestination`](@ref) answered for the tile, lets several pairs reuse
    one preparation of its geometry.
  - What is built, cached and locked is the method's build unit. An area method,
    and a point method with no [`sampler`](@ref), builds and caches the one pair.
  - A point method that supplies a sampler builds the whole destination tile at
    once — one [`TileWeights`](@ref) keyed by tile number, built once however
    many pairs ask for it — and answers `key[2]` from it, or with an empty block
    when no stencil in the tile names that chunk.
"""
blockfor(plan::ChunkedPlan, key::Tuple{Int,Int}, dinds) =
    blockfor(plan, key, dinds, dinds)

blockfor(plan::ChunkedPlan, key::Tuple{Int,Int}, dinds, destination) =
    _blockfor(tilesampler(plan), plan, key, dinds, destination)

_blockfor(::Nothing, plan::ChunkedPlan, key::Tuple{Int,Int}, dinds, destination) =
    getblock!(plan.storage, key,
        () -> buildblock(plan, dinds, ownedindices(plan.src_space, key[2]), destination))

# A tile route reads destination sample sites, not cell geometry, so it takes
# the tile's cells and never the preparation of their polygons.
function _blockfor(smp, plan::ChunkedPlan, key::Tuple{Int,Int}, dinds, destination)
    entry = tileblock(tilefor(plan, key[1], dinds, smp), key[2])
    entry === nothing || return entry
    return _emptyblock(length(dinds), length(ownedindices(plan.src_space, key[2])))
end

blockfor(plan::ChunkedPlan, dstchunk::Integer, srcchunk::Integer) =
    blockfor(plan, (Int(dstchunk), Int(srcchunk)),
        ownedindices(plan.dst_space, Int(dstchunk)))

"""
    tilefor(plan::ChunkedPlan, tile::Integer, dinds, smp) -> CachedTile

Return destination tile `tile`'s cached [`TileWeights`](@ref), building them on
first use from the sampler `smp`.

  - The tile is the build, cache and locking unit: concurrent requests for one
    tile wait on one build rather than starting a second.
  - The whole tile, manifest included, counts against the plan's
    [`weightbudget`](@ref).
"""
tilefor(plan::ChunkedPlan, tile::Integer, dinds, smp) =
    gettile!(plan.storage, Int(tile),
        () -> tileweights(plan.method, plan.dst_space, dinds, plan.src_space, smp))

# A source chunk no stencil in the tile names contributes nothing. This block is
# not cached: it holds no weights to keep, and nothing reads it twice.
_emptyblock(nd::Int, ns::Int) =
    CachedBlock(WeightBlock(sparse(Int[], Int[], Float64[], nd, ns), nothing))

"""
    buildblock(plan::ChunkedPlan, dinds, sinds[, destination]) -> WeightBlock
    buildblock(plan::ChunkedPlan, dstchunk::Integer, srcchunk::Integer) -> WeightBlock

Build one chunk pair's weights without consulting storage.

  - This is always the single pair, through [`weightblock`](@ref).
  - The whole-tile build a point method with a [`sampler`](@ref) uses belongs to
    [`blockfor`](@ref), where the tile is kept.
"""
buildblock(plan::ChunkedPlan, dinds, sinds) =
    buildblock(plan, dinds, sinds, dinds)

buildblock(plan::ChunkedPlan, dinds, sinds, destination) =
    weightblock(plan.method, plan.dst_space, destination, plan.src_space, sinds)

buildblock(plan::ChunkedPlan, dstchunk::Integer, srcchunk::Integer) =
    buildblock(plan, ownedindices(plan.dst_space, Int(dstchunk)),
        ownedindices(plan.src_space, Int(srcchunk)))

"""
    weightblock(method, dst_space, dst_inds, src_space, src_inds) -> WeightBlock

Build the weights destination cells `dst_inds` take from source cells
`src_inds`.

  - Every builder, eager or chunked, goes through here and dispatches on
    [`outputsampling`](@ref), so a sampling may specialise the assembly and no
    concrete method type takes part.
  - `dst_inds` names the destination cells either as their index set or as the
    geometry [`preparedestination`](@ref) prepared for them; nothing here reads
    it either way.
  - Every sampling assembles one pair through [`pairblock`](@ref), whose generic
    route fills one [`WeightCOO`](@ref) through [`buildweights!`](@ref) — all a
    point method that builds no weights of its own needs. The eager whole domain
    and an un-cached [`buildblock`](@ref) go through it.
  - The whole-tile route a point method with a [`sampler`](@ref) takes is not
    chosen here: a chunked plan selects it once, by [`tilesampler`](@ref), in
    [`blockfor`](@ref).
"""
weightblock(method::AbstractRegriddingMethod, dst_space::RegridSpace, dst_inds,
    src_space::RegridSpace, src_inds) =
    weightblock(outputsampling(method), method, dst_space, dst_inds, src_space, src_inds)

weightblock(::DD.Lookups.Sampling, method::AbstractRegriddingMethod,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds) =
    pairblock(method, dst_space, dst_inds, src_space, src_inds)

"""
    pairblock(method, dst_space, dst_inds, src_space, src_inds) -> WeightBlock

Assemble one `(destination cells, source chunk)` pair. Rows are chunk-local
within `dst_inds`, columns chunk-local within `src_inds`.

  - The generic route fills a single [`WeightCOO`](@ref) through
    [`buildweights!`](@ref), which is all a method has to supply.
  - A method that assembles a block of its own specializes here instead, as
    [`Conservative`](@ref) does.
  - A method that wraps another takes the inner method's build by forwarding
    `pairblock` — and [`sampler`](@ref) too, for a point method. Forwarding
    [`buildweights!`](@ref) alone reaches the generic route whatever the inner
    method assembles for itself.
"""
function pairblock(method::AbstractRegriddingMethod, dst_space::RegridSpace, dst_inds,
    src_space::RegridSpace, src_inds)
    coo = WeightCOO(length(dst_inds))
    buildweights!(coo, method, dst_space, dst_inds, src_space, src_inds)
    return WeightBlock(coo, length(dst_inds), length(src_inds))
end

# --------------------------------------------------------------------------
# Point weights, one destination tile at a time
# --------------------------------------------------------------------------

"""
    sampler(method::AbstractRegriddingMethod, space::RegridSpace) -> Nothing

A method with no sampler builds its weights through [`buildweights!`](@ref), one
chunk pair at a time.
"""
sampler(::AbstractRegriddingMethod, ::RegridSpace) = nothing

"""
    tilesampler(plan::ChunkedPlan) -> Union{Nothing,Sampler}

The [`sampler`](@ref) this plan's method prepares for its source space, or
`nothing` when the plan builds one chunk pair at a time.

  - A method whose destination sampling is `Points()` and which supplies a
    sampler builds by destination tile.
  - Every other method, and every point method that supplies no sampler, keeps
    the chunk-pair build.
  - No concrete method type takes part in the choice.
"""
tilesampler(plan::ChunkedPlan) = _tilesampler(outputsampling(plan.method), plan)

_tilesampler(::DD.Lookups.Sampling, ::ChunkedPlan) = nothing
_tilesampler(::DD.Lookups.Points, plan::ChunkedPlan) =
    sampler(plan.method, plan.src_space)

# One source chunk's share of a destination tile: the chunk-local column map and
# the entries that actually landed there. Nothing here is sized by the tile.
struct ChunkAccumulator
    chunk::Int
    ncols::Int
    map::Union{OffsetIndexMap,LookupIndexMap}
    rows::Vector{Int}
    cols::Vector{Int}
    vals::Vector{Float64}
end

# `indexmap` is `i - first(inds) + 1` on a contiguous chunk and a lookup table
# on any other, which is the chunk-local conversion both kinds of space need.
ChunkAccumulator(chunk::Int, inds) =
    ChunkAccumulator(chunk, length(inds), indexmap(inds), Int[], Int[], Float64[])

@inline function _fileentry!(a::ChunkAccumulator, dst_local::Int, src_index::Int,
    w::Float64)
    col = localindex(a.map, src_index)
    col == 0 && throw(ArgumentError(
        "`chunkat` places source cell $src_index in chunk $(a.chunk), whose " *
        "`ownedindices` do not list it; chunks must partition 1:ncells(space)"))
    push!(a.rows, dst_local)
    push!(a.cols, col)
    push!(a.vals, w)
    return a
end

"""
    tileweights(method, dst_space, dinds, src_space, smp) -> TileWeights

Build destination tile `dinds`' weights in one pass over its sample sites.

  - Each destination cell gets one [`weightsat!`](@ref) query whatever the
    source chunking is, so no chunk boundary changes a stencil.
  - Every nonzero entry is filed under the chunk that owns it,
    [`chunkat`](@ref)`(src_space, i)`, by that chunk's chunk-local index;
    entries naming one source cell twice are summed.
  - Blocks are finalized in ascending source-chunk order, and a chunk no stencil
    names produces none, so the manifest is exact.
  - Temporary storage is one [`WeightRow`](@ref), one slot per source chunk
    number, and one accumulator per contributing chunk — never destination cells
    times candidate chunks.
  - Nothing here reads source data or depends on field values or execution order.
"""
function tileweights(method::AbstractRegriddingMethod, dst_space::RegridSpace, dinds,
    src_space::RegridSpace, smp)
    sites = samplesites(dst_space)
    row = WeightRow()
    slots = zeros(Int, Int(nchunks(src_space)))
    accums = ChunkAccumulator[]
    for (j, i) in enumerate(dinds)
        ismapped(weightsat!(row, smp, sites[Int(i)])) || continue
        indices, weights = row.indices, row.weights
        for k in eachindex(indices, weights)
            src_index = indices[k]
            chunk = Int(chunkat(src_space, src_index))
            slot = slots[chunk]
            if slot == 0
                push!(accums, ChunkAccumulator(chunk, ownedindices(src_space, chunk)))
                slot = slots[chunk] = length(accums)
            end
            _fileentry!(accums[slot], j, src_index, weights[k])
        end
    end
    sort!(accums; by = a -> a.chunk)
    nd = length(dinds)
    chunks = Vector{Int}(undef, length(accums))
    blocks = Vector{WeightBlock}(undef, length(accums))
    for k in eachindex(accums)
        a = accums[k]
        chunks[k] = a.chunk
        # `sparse` sums the duplicate entries a stencil naming one source cell
        # twice leaves, so a destination keeps one entry of the summed weight.
        blocks[k] = WeightBlock(sparse(a.rows, a.cols, a.vals, nd, a.ncols), nothing)
    end
    return TileWeights(chunks, blocks)
end
