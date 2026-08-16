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
    PerChunk(capacity)

Keep built [`WeightBlock`](@ref)s in memory, evicting least-recently-used
blocks past `capacity`.
"""
struct PerChunk <: AbstractBlockStorage
    capacity::Int
    PerChunk(::Integer) = error("PerChunk is not yet implemented")
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
    ChunkedPlan(method, missingpolicy, dst_space, src_space, storage, budget)

A plan whose weights are [`WeightBlock`](@ref)s keyed by
`(destination chunk, source chunk)` and built on first touch.

The lazy and streaming case: reading a destination chunk discovers its
connected source chunks, takes each pair's block from `storage`, and
accumulates. `budget` is a performance knob only — it decides whether the
connected source chunks are held together or streamed one at a time — and never
changes the answer.
"""
struct ChunkedPlan{M<:AbstractRegriddingMethod,P<:AbstractMissingPolicy,
                   D<:RegridSpace,S<:RegridSpace,T<:AbstractBlockStorage} <: AbstractRegriddingPlan
    method::M
    missingpolicy::P
    dst_space::D
    src_space::S
    storage::T
    budget::Int
    ChunkedPlan(args...) = error("ChunkedPlan is not yet implemented")
end
