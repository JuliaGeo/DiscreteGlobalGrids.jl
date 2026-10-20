# Conservative weights from spherical intersection areas.

# Candidate discovery descends the grids' trees: the whole grid's sweep tree,
# or the grid's `subtree` over a chunk or any other index set (`_subtree`,
# spaces.jl): its own tree restricted, or a typed chunk's own tree. Leaves
# emit the grids' positions; the block operator maps them to its rows and
# columns through the index maps.

# Prepared destination geometry

"""
    DestinationCache(space, inds) -> DestinationCache or `nothing`

Prepared destination geometry for one tile: its index set, its chunk-local
index map, its restricted tree, and one polygon slot per tile row.

  - Every block the tile takes from a source chunk shares one cache: the
    restricted tree is built once, and each destination polygon synthesized at
    most once, for the whole tile.
  - Slots fill from the candidate pairs a block is about to measure, so a
    destination cell no source overlaps is never synthesized. A fill holds the
    lock; the clipping loop only reads the slot named by its own block row.
  - Not a grid: the destination grid still answers every question about the
    destination.
  - `nothing` when one probe cannot name a concrete polygon type, dynamic
    dispatch in the clipping loop costing more than prepared geometry saves.
  - The polygon type is a parameter and the tree is not, so a cache has one type
    per destination grid whatever tree `GridInterface.subtree` returned.
"""
struct DestinationCache{I,M,P}
    inds::I
    map::M
    # Untyped on purpose: a compile barrier. `subtree` returns one of several
    # tree types for a grid, and the whole-tile assemblies that read this field
    # would otherwise compile once per type. The clipping loop never reads it.
    tree::Any
    # The grid the polygons come from, behind the same barrier.
    space::Any
    polygons::Vector{P}
    # A flag per row rather than `isassigned`, which cannot answer for an
    # isbits polygon: those are stored inline and leave no empty slot.
    filled::Vector{Bool}
    lock::ReentrantLock
end

DestinationCache(space, inds) =
    isempty(inds) ? nothing :
    DestinationCache(space, inds, indexmap(inds),
        _subtree(space, inds; leaf_capacity = GridInterface.ForSweeps()))

# The map and the tree are the tile's, built once by `preparedestination` and
# handed to whichever preparation it settles on.
function DestinationCache(space, inds, imap, tree)
    isempty(inds) && return nothing
    P = typeof(GridInterface.cell_geometry(space, Int(first(inds))))
    isconcretetype(P) || return nothing
    n = length(inds)
    return DestinationCache{typeof(inds),typeof(imap),P}(inds, imap,
        tree, space, Vector{P}(undef, n), fill(false, n), ReentrantLock())
end

Base.show(io::IO, c::DestinationCache) = print(io, "DestinationCache(",
    length(c.inds), " cells, ", count(c.filled), " prepared)")

"""
    DestinationTree(inds, map, tree)

One tile's restricted tree and index map, with no polygon slots: the
preparation a tile gets when its cells' geometry does not fit
[`destcellbudget`](@ref), or when the method reads no cell polygon at all.

Building the tree is the part of preparing that every tile repays. A tile's
tree is a function of its index set alone, so one build serves every block the
tile takes from a source chunk, and on a whole-space DGGS destination that
build is a cap per cell — 823,543 of them for one production chunk, which is
most of the block's cost when each block repeats it. The polygons are the part
the budget may refuse; refusing them is not a reason to rebuild the tree.

Blocks addressed this way synthesize destination geometry through their own
task-local [`CellMemo`](@ref), exactly as an index set does.
"""
struct DestinationTree{I,M}
    inds::I
    map::M
    # Untyped for the same reason `DestinationCache.tree` is.
    tree::Any
end

Base.show(io::IO, d::DestinationTree) =
    print(io, "DestinationTree(", length(d.inds), " cells)")

"""
    preparesdestination(method, dst_grid) -> Bool

Whether a build of `method` over `dst_grid` keeps the tile's destination cell
polygons, given room for them.

An area method reads a destination cell's polygon once per overlapping source
leaf, so it keeps them where synthesizing one is expensive
(`GridInterface.cell_geometry_is_expensive`). A point method evaluates at a
destination's sample site and reads no cell polygon, so it keeps none.
"""
preparesdestination(::AbstractRegriddingMethod, ::Any) = false
preparesdestination(::Conservative, dst_space) =
    GridInterface.cell_geometry_is_expensive(dst_space)

"""
    preparedestination(method, dst_space, dst_inds, budget) -> DestinationCache, DestinationTree or dst_inds

The destination a tile's builds address: the tile's tree with prepared polygons
where `method` keeps them ([`preparesdestination`](@ref)) and `budget` holds
them ([`destcellsfit`](@ref)), the tile's tree alone otherwise, and `dst_inds`
itself for an empty tile.

  - All three answers reach the same weights, value for value and entry for entry.
  - Only a tile with more than one block repays a preparation, so a single-block
    build — the eager whole domain included — takes the index set and its
    task-local [`CellMemo`](@ref) instead.
  - The tree is hoisted whether or not the polygons are. It is built from the
    index set alone and every block of the tile clips against the same one, so
    the budget's refusal applies to the polygon slots and to nothing else. A
    [`DestinationTree`](@ref) carries the tree without them.
"""
function preparedestination(method::AbstractRegriddingMethod,
    dst_space, dst_inds, budget::Integer)
    isempty(dst_inds) && return dst_inds
    imap = indexmap(dst_inds)
    tree = _subtree(dst_space, dst_inds; leaf_capacity = GridInterface.ForSweeps())
    if preparesdestination(method, dst_space) &&
       destcellsfit(method, dst_space, length(dst_inds), budget)
        cache = DestinationCache(dst_space, dst_inds, imap, tree)
        cache === nothing || return cache
    end
    return DestinationTree(dst_inds, imap, tree)
end

# Fill every row the pairs name and no other, once per row for the whole tile.
# Concurrent blocks of one tile fill under the lock, and a block's own reads
# happen after it releases the lock, so a slot another block filled is visible.
function _preparepolygons!(cache::DestinationCache, pairs)
    lock(cache.lock)
    try
        _fillpolygons!(cache.polygons, cache.filled, cache.map, cache.space, pairs)
    finally
        unlock(cache.lock)
    end
    return cache
end

function _fillpolygons!(polygons::Vector{P}, filled::Vector{Bool}, map, space,
    pairs) where {P}
    @inbounds for (_, i2) in pairs
        row = localindex(map, i2)
        (row == 0 || filled[row]) && continue
        polygons[row] = GridInterface.cell_geometry(space, i2)
        filled[row] = true
    end
    return polygons
end

"""
    _destinationtree(dst_space, dst_inds)
    _destinationtree(prepared::Union{DestinationCache,DestinationTree})

The tree a block clips its destination against, held opaque to inference.

`GridInterface.subtree` answers one of several tree types for a grid. Crossing an
inference barrier compiles a block's assembly for the one type the run reaches
rather than for all of them, behind a single dynamic dispatch per block; the
clipping loop runs inside that specialization and pays nothing.
"""
_destinationtree(dst_space, dst_inds) =
    Base.inferencebarrier(_subtree(dst_space, dst_inds; leaf_capacity = GridInterface.ForSweeps()))
# A preparation holds its tree opaque already.
_destinationtree(cache::DestinationCache) = cache.tree
_destinationtree(prepared::DestinationTree) = prepared.tree
_destinationtree(::Any, prepared::Union{DestinationCache,DestinationTree}) =
    _destinationtree(prepared)

# The index set a block's destination stands for, however it was named.
_destinationinds(dst_inds) = dst_inds
_destinationinds(prepared::Union{DestinationCache,DestinationTree}) = prepared.inds

# Intersection operator

"""
    CellMemo(::Type{P}; slots)

A task-local direct-mapped memo of cell polygons by tree index. Candidate
pairs leave a leaf pairing in index runs — the destination repeats within a
pairing and the source's few leaf cells cycle across consecutive pairings — so
a handful of slots serves most `cell_geometry` calls with the value already
built. The memo reads cells from the grid the operator is handed
(`GridInterface.cell_geometry(grid, i)`).
"""
struct CellMemo{P}
    keys::Vector{Int}
    vals::Vector{P}
end

CellMemo(::Type{P}; slots::Int = 64) where {P} =
    CellMemo{P}(zeros(Int, slots), Vector{P}(undef, slots))

# A fresh, empty memo of the same shape, for one assembly task.
_fresh(memo::CellMemo{P}) where {P} = CellMemo(P; slots = length(memo.keys))
_fresh(::Nothing) = nothing

@inline function _memocell(memo::CellMemo, space, i::Int)
    slot = (i & (length(memo.keys) - 1)) + 1
    @inbounds memo.keys[slot] == i && return @inbounds memo.vals[slot]
    cell = GridInterface.cell_geometry(space, i)
    @inbounds memo.keys[slot] = i
    @inbounds memo.vals[slot] = cell
    return cell
end

@inline _memocell(::Nothing, space, i::Int) = GridInterface.cell_geometry(space, i)

# A concretely typed memo where one probe names the polygon type; otherwise none.
function _cellmemo(space, inds)
    P = typeof(GridInterface.cell_geometry(space, Int(first(inds))))
    return isconcretetype(P) ? CellMemo(P) : nothing
end

"""
    BlockAreaOperator(inner, dstmap, srcmap, srcmemo, dstcells)

Measure intersections with `inner` and map the grids' positions to chunk-local
block rows and columns (the engine hands the operator the two grids; their
trees only drive the candidate search). Each assembly task receives its own
mutable clipping cache and its own source cell memo.

`dstcells` is where the destination's geometry comes from: a
[`DestinationCache`](@ref) the whole tile shares, a task-local
[`CellMemo`](@ref), or `nothing` for a memo-free build. A cache answers by the
block row the operator has already computed; the other two synthesize.
"""
struct BlockAreaOperator{O,DM,SM,MS,MD}
    inner::O
    dstmap::DM
    srcmap::SM
    srcmemo::MS
    dstcells::MD
end

# Memo-free form: every `cell_geometry` synthesizes. The production path passes memos.
BlockAreaOperator(inner, dstmap, srcmap) =
    BlockAreaOperator(inner, dstmap, srcmap, nothing, nothing)

ConservativeRegridding.IntersectionReturnStyle(::BlockAreaOperator) =
    ConservativeRegridding.InPlace()

ConservativeRegridding.output_matrix_size(op::BlockAreaOperator, src, dst) =
    (length(op.dstmap), length(op.srcmap))

ConservativeRegridding.task_local_operator(op::BlockAreaOperator) =
    BlockAreaOperator(ConservativeRegridding.task_local_operator(op.inner),
        op.dstmap, op.srcmap, _fresh(op.srcmemo), _fresh(op.dstcells))

# One tile's cache is read-only for the length of an assembly and is shared by
# every task, rather than copied per task as a memo is.
_fresh(cache::DestinationCache) = cache

# The candidate pairs are complete and ordered before any assembly task starts,
# so this is where the block learns exactly which destination rows it measures.
ConservativeRegridding.work_items(
    op::BlockAreaOperator{O,DM,SM,MS,<:DestinationCache}, pairs) where {O,DM,SM,MS} =
    (_preparepolygons!(op.dstcells, pairs); pairs)

# The source is the subject and the destination is the clip ring. A block emits
# weights only for the pairs both of its chunks contain.
@inline function (op::BlockAreaOperator)(rows, cols, vals, item, src_space, dst_space)
    i1, i2 = item
    row = localindex(op.dstmap, i2)
    col = localindex(op.srcmap, i1)
    (row == 0 || col == 0) && return nothing
    area = op.inner(_memocell(op.srcmemo, src_space, i1),
        _destinationcell(op.dstcells, dst_space, i2, row))
    area > 0 || return nothing
    push!(rows, row)
    push!(cols, col)
    push!(vals, area)
    return nothing
end

@inline _destinationcell(dstcells, dst_space, i::Int, ::Int) =
    _memocell(dstcells, dst_space, i)
@inline _destinationcell(cache::DestinationCache, dst_space, ::Int, row::Int) =
    @inbounds cache.polygons[row]

# Conservative method

"""
    buildweights!(coo, ::Conservative, dst_space, dst_inds, src_space, src_inds,
                  locator = TreeLocator())

Append spherical intersection areas for the two chunks. Each destination
denominator accumulates its covered area. A conservative block always carries a
denominator, including when coverage is zero. Source and destination manifolds
must match. Intersection discovery and clipping may run in parallel. `locator`
chooses how candidate pairs are discovered ([`CandidateLocator`](@ref)); the
weights do not depend on it.

This is the generic [`WeightCOO`](@ref) route, which any method may build
through. [`pairblock`](@ref)`(::Conservative, …)` is the route a conservative
plan takes, and reaches the same weights without a coordinate list.
"""
function buildweights!(coo::WeightCOO, ::Conservative,
    dst_space, dst_inds, src_space, src_inds,
    locator::CandidateLocator = TreeLocator())
    isempty(dst_inds) && return coo
    markdenominated!(coo)
    isempty(src_inds) && return coo

    m = _sharedmanifold(dst_space, src_space)
    op = BlockAreaOperator(_intersectionoperator(m),
        indexmap(dst_inds), indexmap(src_inds),
        _cellmemo(src_space, src_inds), _cellmemo(dst_space, dst_inds))
    block = _intersectionareas(m, locator, dst_space, dst_inds, src_space, src_inds, op)

    return _fillcoo!(coo, block)
end

"""
    pairblock(::Conservative, dst_space, dst_inds, src_space, src_inds[, locator]) -> WeightBlock
    pairblock(::Conservative, dst_space, dst_cache::DestinationCache, src_space, src_inds[, locator])
    pairblock(::Conservative, dst_space, prepared::DestinationTree, src_space, src_inds[, locator])

Adopt the assembled sparse matrix of intersection areas as the block's weights,
reading each destination's denominator off it once. Every conservative block —
the eager whole domain and a chunk pair alike — is built here. `locator`
([`CandidateLocator`](@ref), default [`TreeLocator`](@ref)) chooses how the
candidate pairs are discovered and leaves the weights unchanged.

  - The values, their CSC layout, and the denominator accumulation order are the
    generic [`WeightCOO`](@ref) route's, bit for bit.
  - Naming the destination by index set builds its restricted tree here and
    synthesizes its geometry through a task-local [`CellMemo`](@ref); naming it
    by a [`DestinationTree`](@ref) reuses the tile's tree and keeps the memo;
    naming it by a [`DestinationCache`](@ref) reuses the tile's tree and its
    prepared polygons. The weights are the same, entry for entry.
"""
function pairblock(::Conservative, dst_space, dst_inds,
    src_space, src_inds, locator::CandidateLocator = TreeLocator())
    ndst = length(dst_inds)
    # A degenerate side keeps the generic route's exact semantics, including
    # which of the two sides reports a denominator.
    (ndst == 0 || isempty(src_inds)) && return invoke(pairblock,
        Tuple{AbstractRegriddingMethod,Any,Any,Any,Any},
        Conservative(), dst_space, dst_inds, src_space, src_inds)

    m = _sharedmanifold(dst_space, src_space)
    op = BlockAreaOperator(_intersectionoperator(m),
        indexmap(dst_inds), indexmap(src_inds),
        _cellmemo(src_space, src_inds), _cellmemo(dst_space, dst_inds))
    block = _intersectionareas(m, locator, dst_space, dst_inds, src_space, src_inds, op)

    return WeightBlock(block, _blockdenom(block, ndst))
end

function pairblock(::Conservative, dst_space,
    dst_cache::DestinationCache, src_space, src_inds,
    locator::CandidateLocator = TreeLocator())
    ndst = length(dst_cache.inds)
    (ndst == 0 || isempty(src_inds)) && return invoke(pairblock,
        Tuple{AbstractRegriddingMethod,Any,Any,Any,Any},
        Conservative(), dst_space, dst_cache.inds, src_space, src_inds)

    m = _sharedmanifold(dst_space, src_space)
    op = BlockAreaOperator(_intersectionoperator(m),
        dst_cache.map, indexmap(src_inds),
        _cellmemo(src_space, src_inds), dst_cache)
    block = _intersectionareas(m, locator, dst_space, dst_cache, src_space, src_inds, op)

    return WeightBlock(block, _blockdenom(block, ndst))
end

# The tree-only preparation: the tile's tree, and a task-local memo for the
# geometry it holds no slots for. Same assembly as the index-set spelling, one
# tree build per tile rather than one per block.
function pairblock(::Conservative, dst_space,
    prepared::DestinationTree, src_space, src_inds,
    locator::CandidateLocator = TreeLocator())
    ndst = length(prepared.inds)
    (ndst == 0 || isempty(src_inds)) && return invoke(pairblock,
        Tuple{AbstractRegriddingMethod,Any,Any,Any,Any},
        Conservative(), dst_space, prepared.inds, src_space, src_inds)

    m = _sharedmanifold(dst_space, src_space)
    op = BlockAreaOperator(_intersectionoperator(m),
        prepared.map, indexmap(src_inds),
        _cellmemo(src_space, src_inds), _cellmemo(dst_space, prepared.inds))
    block = _intersectionareas(m, locator, dst_space, prepared, src_space, src_inds, op)

    return WeightBlock(block, _blockdenom(block, ndst))
end

# A method that reads no prepared geometry still takes the destination a tile
# prepared, and names its cells.
pairblock(method::AbstractRegriddingMethod, dst_space,
    dst_cache::DestinationCache, src_space, src_inds) =
    pairblock(method, dst_space, dst_cache.inds, src_space, src_inds)

pairblock(method::AbstractRegriddingMethod, dst_space,
    prepared::DestinationTree, src_space, src_inds) =
    pairblock(method, dst_space, prepared.inds, src_space, src_inds)

pairblock(method::AbstractRegriddingMethod, dst_space,
    dst_cache::DestinationCache, src_space, src_inds,
    locator::CandidateLocator) =
    pairblock(method, dst_space, dst_cache.inds, src_space, src_inds, locator)

pairblock(method::AbstractRegriddingMethod, dst_space,
    prepared::DestinationTree, src_space, src_inds,
    locator::CandidateLocator) =
    pairblock(method, dst_space, prepared.inds, src_space, src_inds, locator)

function _sharedmanifold(dst_space, src_space)
    m = GOCore.best_manifold(dst_space)
    m == GOCore.best_manifold(src_space) || throw(ArgumentError(
        "conservative weights need one manifold on both sides: destination is " *
        "$(m), source is $(GOCore.best_manifold(src_space))"))
    return m
end

# `_fillcoo!`'s denominator pass, over an adopted matrix. Assembly types the
# block only as `SparseMatrixCSC`, so the loop needs its own dispatch to specialise.
function _blockdenom(block::SparseArrays.AbstractSparseMatrixCSC, ndst::Int)
    denom = zeros(Float64, ndst)
    rows = SparseArrays.rowvals(block)
    vals = SparseArrays.nonzeros(block)
    @inbounds for col in axes(block, 2), t in SparseArrays.nzrange(block, col)
        w = vals[t]
        w > 0 || continue
        denom[rows[t]] += w
    end
    return denom
end

"""
    _intersectionareas(manifold, locator, dst_grid, dst, src_grid, src_inds, op)
    _intersectionareas(manifold, dst_grid, dst_tree, src_grid, src_tree, op) -> SparseMatrixCSC

Compute intersection areas, threaded when more than one thread is available
and no outer loop is already parallel ([`OUTER_PARALLEL`](@ref)). `dst` is the
destination's index set or its [`preparedestination`](@ref) result.

Both routes run `ConservativeRegridding.intersection_areas` over the two grids
with the block operator. A [`TreeLocator`](@ref) hands the engine the two
grids' restricted trees, whose leaves emit the grids' positions, and the
engine's dual-tree search names the candidate pairs. An
[`AnalyticLocator`](@ref) discovers the pairs itself with
[`overlappairs`](@ref) and hands them to the engine's positional-pairs form.
Either way the operator reads cells through `GridInterface.cell_geometry(grid,
i)` and maps positions to block rows and columns. The two routes assemble the
same matrix, entry for entry.
"""
function _intersectionareas(m::GOCore.Manifold, ::TreeLocator, dst_space,
    dst, src_space, src_inds, op)
    return _intersectionareas(m, dst_space, _destinationtree(dst_space, dst),
        src_space, _subtree(src_space, src_inds; leaf_capacity = GridInterface.ForSweeps()), op)
end

function _intersectionareas(m::GOCore.Manifold, locator::AnalyticLocator,
    dst_space, dst, src_space, src_inds, op)
    threaded = _innerthreaded()
    pairs = overlappairs(locator, dst_space, _destinationinds(dst), src_space, src_inds;
        threaded = threaded isa GOCore.True)
    ValType = ConservativeRegridding.output_eltype(op, src_space, dst_space)
    return _with_sparse_scratch(ValType) do cache
        ConservativeRegridding.intersection_areas(
            m, threaded, dst_space, src_space, pairs;
            intersection_operator = op, npartitions = Threads.nthreads() * 4, cache)
    end
end

function _intersectionareas(m::GOCore.Manifold, dst_space, dst_tree,
    src_space, src_tree, op)
    threaded = _innerthreaded()
    ValType = ConservativeRegridding.output_eltype(op, src_space, dst_space)
    return _with_sparse_scratch(ValType) do cache
        ConservativeRegridding.intersection_areas(
            m, threaded, dst_space, src_space;
            intersection_operator = op, dst_tree, src_tree, cache)
    end
end

# Assembly scratch is module-owned because the block-build call chain has no
# operation context.  The freelist grows on an empty acquire rather than
# waiting for another build, and its retained size is therefore bounded by the
# peak number of overlapping builds for each value type.
struct _SparseScratchKey{T} end

mutable struct _SparseScratchPool{T}
    lock::ReentrantLock
    free::Vector{ConservativeRegridding.SparseMatrixAssemblyCache{T}}
end

_SparseScratchPool(::Type{T}) where {T} =
    _SparseScratchPool{T}(ReentrantLock(),
        ConservativeRegridding.SparseMatrixAssemblyCache{T}[])

const _SPARSE_SCRATCH_POOL_LOCK = ReentrantLock()
const _SPARSE_SCRATCH_POOLS = Dict{DataType,Any}()

function _sparse_scratch_pool(::Type{T}) where {T}
    key = _SparseScratchKey{T}
    lock(_SPARSE_SCRATCH_POOL_LOCK)
    try
        return get!(_SPARSE_SCRATCH_POOLS, key) do
            _SparseScratchPool(T)
        end::_SparseScratchPool{T}
    finally
        unlock(_SPARSE_SCRATCH_POOL_LOCK)
    end
end

function _acquire_sparse_scratch(pool::_SparseScratchPool{T}) where {T}
    lock(pool.lock)
    try
        return isempty(pool.free) ?
            ConservativeRegridding.SparseMatrixAssemblyCache(T) : pop!(pool.free)
    finally
        unlock(pool.lock)
    end
end

function _release_sparse_scratch!(pool::_SparseScratchPool{T},
    cache::ConservativeRegridding.SparseMatrixAssemblyCache{T}) where {T}
    lock(pool.lock)
    try
        push!(pool.free, cache)
    finally
        unlock(pool.lock)
    end
    return cache
end

function _with_sparse_scratch(f, ::Type{T}) where {T}
    pool = _sparse_scratch_pool(T)
    cache = _acquire_sparse_scratch(pool)
    try
        return f(cache)
    finally
        _release_sparse_scratch!(pool, cache)
    end
end

# Copy stored entries directly to avoid `findnz` allocations.
function _fillcoo!(coo::WeightCOO, block::SparseArrays.AbstractSparseMatrixCSC)
    rows = SparseArrays.rowvals(block)
    vals = SparseArrays.nonzeros(block)
    @inbounds for col in axes(block, 2), t in SparseArrays.nzrange(block, col)
        w = vals[t]
        w > 0 || continue
        addweight!(coo, rows[t], col, w)
        adddenom!(coo, rows[t], w)
    end
    return coo
end

function _fillcoo!(coo::WeightCOO, block::AbstractMatrix)
    @inbounds for col in axes(block, 2), row in axes(block, 1)
        w = block[row, col]
        w > 0 || continue
        addweight!(coo, row, col, w)
        adddenom!(coo, row, w)
    end
    return coo
end
