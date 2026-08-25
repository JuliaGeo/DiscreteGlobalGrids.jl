# Conservative weights from spherical intersection areas.

# Chunk-restricted cell trees

function subtree(space::RegridSpace, inds)
    _iswholespace(space, inds) && return celltree(space)
    return CellSpaceRTree(space, inds)
end

_iswholespace(space::RegridSpace, inds::AbstractUnitRange{<:Integer}) =
    first(inds) == 1 && length(inds) == ncells(space)
_iswholespace(::RegridSpace, inds) = false

"""
    CellSpaceRTree(space, inds; nodecapacity = 16)

The unstructured cell fallback: a thin cell-space adapter over GeometryOps'
packed `FlexibleRTrees.RTree`. The R-tree owns outward-rounded Cartesian
extents, while the adapter retains `space` for geometry access and maps every
packed leaf back to its requested local index in `space`.

Structured spaces should return a native restricted cursor before reaching
this fallback. `nodecapacity` is exposed for packing-invariance checks.
"""
struct CellSpaceRTree{S,N,V<:AbstractVector{Int},C<:AbstractVector{<:Cap}}
    space::S
    node::N
    indices::V
    caps::C
    leafcount::Int
    height::Int
    nodecapacity::Int
end

function CellSpaceRTree(space::S, inds; nodecapacity::Int = 16) where {S<:RegridSpace}
    indices = collect(Int, inds)
    isempty(indices) && throw(ArgumentError(
        "cannot build a cell tree from an empty index set"))
    caps = [_packedcellcap(space, i) for i in indices]
    boxes = map(cap -> convert(Extents.Extent, cap), caps)
    packed = FlexibleRTrees.RTree(FlexibleRTrees.HPR(), indices;
        nodecapacity, extents = boxes)
    return CellSpaceRTree(space, packed, indices, caps, length(indices),
        length(packed.levels), nodecapacity)
end

# A cap of at most one hemisphere is geodesically convex, so it covers the cell
# edges as well as the vertices. Wider vertex caps conservatively become the
# whole sphere. Construction and merging use GeometryOps' public cap protocol;
# the small fractional growth preserves the previous fallback's headroom.
function _packedcellcap(space::RegridSpace, i::Int)
    points = GI.getpoint(getcell(space, i))
    isempty(points) && return _WHOLE_SPHERE
    firstpoint = first(points)
    cap = SphericalCap(USPoint(
        GI.x(firstpoint), GI.y(firstpoint), GI.z(firstpoint)), 0.0)
    for point in Iterators.drop(points, 1)
        cap = Extents.union(cap, SphericalCap(
            USPoint(GI.x(point), GI.y(point), GI.z(point)), 0.0))
        cap.radius > Float64(pi) / 2 && return _WHOLE_SPHERE
    end
    return Extents.grow(cap, 5e-5)
end

Base.show(io::IO, tree::CellSpaceRTree) =
    print(io, "CellSpaceRTree(", tree.leafcount, " cells, ", tree.node, ")")

STI.isspatialtree(::Type{<:CellSpaceRTree}) = true
STI.node_extent_is_expensive(::Type{<:CellSpaceRTree}) = false
STI.isleaf(tree::CellSpaceRTree) = STI.isleaf(tree.node)
STI.nchild(tree::CellSpaceRTree) = STI.nchild(tree.node)
STI.node_extent(tree::CellSpaceRTree) = STI.node_extent(tree.node)

function STI.getchild(tree::CellSpaceRTree, i::Int)
    childheight = tree.height - 1
    childcapacity = tree.nodecapacity^childheight
    childcount = clamp(tree.leafcount - (i - 1) * childcapacity, 0, childcapacity)
    return CellSpaceRTree(tree.space, STI.getchild(tree.node, i), tree.indices,
        tree.caps, childcount, childheight, tree.nodecapacity)
end

STI.getchild(tree::CellSpaceRTree) =
    (STI.getchild(tree, i) for i in 1:STI.nchild(tree))

# FlexibleRTrees reports slots in its payload vector. Translate those back to
# the space's stable local indices, and expose their exact cap rather than the
# cap's broad-phase XYZ box.
STI.child_indices_extents(tree::CellSpaceRTree) =
    ((@inbounds(tree.indices[i]), @inbounds(tree.caps[i]))
     for (i, _) in STI.child_indices_extents(tree.node))

# `ncells` and `split_weight` describe this node's restricted population;
# matrix assembly separately needs the owning space's complete cell count.
ncells(tree::CellSpaceRTree) = tree.leafcount
Trees.cell_index_count(tree::CellSpaceRTree) = ncells(tree.space)
Trees.split_weight(tree::CellSpaceRTree) = tree.leafcount
getcell(tree::CellSpaceRTree, i::Int) = getcell(tree.space, i)

function _packedcellindices!(out::Vector{Int}, tree::CellSpaceRTree)
    if STI.isleaf(tree)
        append!(out, first(entry) for entry in STI.child_indices_extents(tree))
    else
        for child in STI.getchild(tree)
            _packedcellindices!(out, child)
        end
    end
    return out
end

getcell(tree::CellSpaceRTree) =
    (getcell(tree.space, i) for i in _packedcellindices!(Int[], tree))
GOCore.best_manifold(tree::CellSpaceRTree) = manifold(tree.space)

# Prepared destination geometry

"""
    DestinationCache(space, inds) -> DestinationCache or `nothing`

Prepared destination geometry for one tile: its index set, its chunk-local
index map, its restricted tree, and one polygon slot per tile row.

Every block the tile takes from a source chunk shares one cache, so the
restricted tree is built once for the tile and each destination polygon is
synthesized at most once across all of the tile's blocks. Slots fill from the
candidate pairs a block is about to measure, so a destination cell no source
overlaps is never synthesized. A fill holds the lock; the clipping loop only
reads the slot named by the block row it has already computed.

This is not a `RegridSpace` and answers no question about the destination:
the destination space still does. Construction returns `nothing` when one
probe cannot name a concrete polygon type, because dynamic dispatch in the
clipping loop costs more than prepared geometry saves.

The polygon type is a parameter and the tree is not, so a cache has one type
per destination space whatever [`subtree`](@ref) returned.
"""
struct DestinationCache{I,M,P}
    inds::I
    map::M
    # Untyped on purpose: a compile barrier. `subtree` returns one of several
    # tree types for a space, and both of the calls that read this field are
    # whole assemblies over the tile. Inferring the field would compile each
    # of them for every tree type the space can produce; leaving it opaque
    # compiles each for the type a run actually reaches, behind one dynamic
    # dispatch per block. The clipping loop never reads it.
    tree::Any
    polygons::Vector{P}
    # A flag per row rather than `isassigned`, which cannot answer for an
    # isbits polygon: those are stored inline and leave no empty slot.
    filled::Vector{Bool}
    lock::ReentrantLock
end

function DestinationCache(space::RegridSpace, inds)
    isempty(inds) && return nothing
    P = typeof(getcell(space, Int(first(inds))))
    isconcretetype(P) || return nothing
    n = length(inds)
    imap = indexmap(inds)
    return DestinationCache{typeof(inds),typeof(imap),P}(inds, imap,
        subtree(space, inds), Vector{P}(undef, n), fill(false, n), ReentrantLock())
end

Base.show(io::IO, c::DestinationCache) = print(io, "DestinationCache(",
    length(c.inds), " cells, ", count(c.filled), " prepared)")

"""
    preparesdestination(method, dst_space::RegridSpace) -> Bool

Whether a build of `method` over `dst_space` keeps the tile's destination cell
polygons, given room for them.

An area method reads a destination cell's polygon once per overlapping source
leaf, so it keeps them where synthesizing one is expensive
([`expensivecellgeometry`](@ref)). A point method evaluates at a destination's
sample site and reads no cell polygon, so it keeps none.
"""
preparesdestination(::AbstractRegriddingMethod, ::RegridSpace) = false
preparesdestination(::Conservative, dst_space::RegridSpace) =
    expensivecellgeometry(dst_space)

"""
    preparedestination(method, dst_space, dst_inds, budget) -> dst_inds or DestinationCache

The destination a tile's builds address: prepared geometry where `method`
keeps it ([`preparesdestination`](@ref)) and `budget` holds it
([`destcellsfit`](@ref)), and `dst_inds` itself otherwise.

Both answers reach the same weights, value for value and entry for entry.
Preparing is worth its memory only where more than one block reads it — one
tile prepares once and every block of that tile shares the result — so a build
of a single block, the eager whole domain included, takes the index set and
its task-local [`CellMemo`](@ref) instead.
"""
function preparedestination(method::AbstractRegriddingMethod,
    dst_space::RegridSpace, dst_inds, budget::Integer)
    preparesdestination(method, dst_space) || return dst_inds
    destcellsfit(method, dst_space, length(dst_inds), budget) || return dst_inds
    cache = DestinationCache(dst_space, dst_inds)
    return cache === nothing ? dst_inds : cache
end

# Fill every row the pairs name and no other, once per row for the whole tile.
# Concurrent blocks of one tile fill under the lock, and a block's own reads
# happen after it releases the lock, so a slot another block filled is visible.
function _preparepolygons!(cache::DestinationCache, pairs)
    lock(cache.lock)
    try
        _fillpolygons!(cache.polygons, cache.filled, cache.map, cache.tree, pairs)
    finally
        unlock(cache.lock)
    end
    return cache
end

function _fillpolygons!(polygons::Vector{P}, filled::Vector{Bool}, map, tree,
    pairs) where {P}
    @inbounds for (_, i2) in pairs
        row = localindex(map, i2)
        (row == 0 || filled[row]) && continue
        polygons[row] = Trees.getcell(tree, i2)
        filled[row] = true
    end
    return polygons
end

"""
    _destinationtree(dst_space, dst_inds)
    _destinationtree(cache::DestinationCache)

The tree a block clips its destination against, held opaque to inference.

[`subtree`](@ref) answers one of several tree types for a space, and the
assembly a block runs is compiled for the tree it is handed. Inferring the
answer compiles that assembly for every tree type the space can produce;
crossing a barrier compiles it for the one the run reaches, behind a single
dynamic dispatch per block. The clipping loop runs inside that specialization
and pays nothing.
"""
_destinationtree(dst_space::RegridSpace, dst_inds) =
    Base.inferencebarrier(subtree(dst_space, dst_inds))
# The cache holds its tree opaque already.
_destinationtree(cache::DestinationCache) = cache.tree

# Intersection operator

"""
    CellMemo(::Type{P}; slots)

A task-local direct-mapped memo of cell polygons by tree index. Candidate
pairs leave a leaf pairing in index runs — the destination repeats within a
pairing and the source's few leaf cells cycle across consecutive pairings — so
a handful of slots serves most `getcell` calls with the value already built.
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

@inline function _memocell(memo::CellMemo, tree, i::Int)
    slot = (i & (length(memo.keys) - 1)) + 1
    @inbounds memo.keys[slot] == i && return @inbounds memo.vals[slot]
    cell = Trees.getcell(tree, i)
    @inbounds memo.keys[slot] = i
    @inbounds memo.vals[slot] = cell
    return cell
end

@inline _memocell(::Nothing, tree, i::Int) = Trees.getcell(tree, i)

# A concretely typed memo where one probe names the polygon type; otherwise none.
function _cellmemo(space::RegridSpace, inds)
    P = typeof(getcell(space, Int(first(inds))))
    return isconcretetype(P) ? CellMemo(P) : nothing
end

"""
    BlockAreaOperator(inner, dstmap, srcmap, srcmemo, dstcells)

Measure intersections with `inner` and map the trees' local indices to
chunk-local block rows and columns. Each assembly task receives its own mutable
clipping cache and its own source cell memo.

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

# Memo-free form: every `getcell` synthesizes. The production path passes memos.
BlockAreaOperator(inner, dstmap, srcmap) =
    BlockAreaOperator(inner, dstmap, srcmap, nothing, nothing)

ConservativeRegridding.IntersectionReturnStyle(::BlockAreaOperator) =
    ConservativeRegridding.InPlace()

ConservativeRegridding.output_matrix_size(op::BlockAreaOperator, src_tree, dst_tree) =
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
@inline function (op::BlockAreaOperator)(rows, cols, vals, item, src_tree, dst_tree)
    i1, i2 = item
    row = localindex(op.dstmap, i2)
    col = localindex(op.srcmap, i1)
    (row == 0 || col == 0) && return nothing
    area = op.inner(_memocell(op.srcmemo, src_tree, i1),
        _destinationcell(op.dstcells, dst_tree, i2, row))
    area > 0 || return nothing
    push!(rows, row)
    push!(cols, col)
    push!(vals, area)
    return nothing
end

@inline _destinationcell(dstcells, dst_tree, i::Int, ::Int) =
    _memocell(dstcells, dst_tree, i)
@inline _destinationcell(cache::DestinationCache, dst_tree, ::Int, row::Int) =
    @inbounds cache.polygons[row]

# Conservative method

"""
    buildweights!(coo, ::Conservative, dst_space, dst_inds, src_space, src_inds)

Append spherical intersection areas for the two chunks. Each destination
denominator accumulates its covered area. A conservative block always carries a
denominator, including when coverage is zero. Source and destination manifolds
must match. Intersection discovery and clipping may run in parallel.

This is the generic [`WeightCOO`](@ref) route, which any method may build
through. [`pairblock`](@ref)`(::Conservative, …)` is the route a conservative
plan takes, and reaches the same weights without a coordinate list.
"""
function buildweights!(coo::WeightCOO, ::Conservative,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    isempty(dst_inds) && return coo
    markdenominated!(coo)
    isempty(src_inds) && return coo

    m = _sharedmanifold(dst_space, src_space)
    op = BlockAreaOperator(_intersectionoperator(m),
        indexmap(dst_inds), indexmap(src_inds),
        _cellmemo(src_space, src_inds), _cellmemo(dst_space, dst_inds))
    block = _intersectionareas(m, _destinationtree(dst_space, dst_inds),
        subtree(src_space, src_inds), op)

    return _fillcoo!(coo, block)
end

"""
    pairblock(::Conservative, dst_space, dst_inds, src_space, src_inds) -> WeightBlock
    pairblock(::Conservative, dst_space, dst_cache::DestinationCache, src_space, src_inds)

Adopt the assembled sparse matrix of intersection areas as the block's weights,
reading each destination's denominator off it once. Every conservative block —
the eager whole domain and a chunk pair alike — is built here.

The generic route copies every entry into a [`WeightCOO`](@ref) and assembles a
second, identical CSC from it. The values, their CSC layout, and the denominator
accumulation order here are that route's, bit for bit.

Naming the destination by index set builds its restricted tree here and
synthesizes its geometry through a task-local [`CellMemo`](@ref). Naming it by
a [`DestinationCache`](@ref) reuses the tile's tree and its prepared polygons
instead; the weights are the same, entry for entry.
"""
function pairblock(::Conservative, dst_space::RegridSpace, dst_inds,
    src_space::RegridSpace, src_inds)
    ndst = length(dst_inds)
    # A degenerate side keeps the generic route's exact semantics, including
    # which of the two sides reports a denominator.
    (ndst == 0 || isempty(src_inds)) && return invoke(pairblock,
        Tuple{AbstractRegriddingMethod,RegridSpace,Any,RegridSpace,Any},
        Conservative(), dst_space, dst_inds, src_space, src_inds)

    m = _sharedmanifold(dst_space, src_space)
    op = BlockAreaOperator(_intersectionoperator(m),
        indexmap(dst_inds), indexmap(src_inds),
        _cellmemo(src_space, src_inds), _cellmemo(dst_space, dst_inds))
    block = _intersectionareas(m, _destinationtree(dst_space, dst_inds),
        subtree(src_space, src_inds), op)

    return WeightBlock(block, _blockdenom(block, ndst))
end

function pairblock(::Conservative, dst_space::RegridSpace,
    dst_cache::DestinationCache, src_space::RegridSpace, src_inds)
    ndst = length(dst_cache.inds)
    (ndst == 0 || isempty(src_inds)) && return invoke(pairblock,
        Tuple{AbstractRegriddingMethod,RegridSpace,Any,RegridSpace,Any},
        Conservative(), dst_space, dst_cache.inds, src_space, src_inds)

    m = _sharedmanifold(dst_space, src_space)
    op = BlockAreaOperator(_intersectionoperator(m),
        dst_cache.map, indexmap(src_inds),
        _cellmemo(src_space, src_inds), dst_cache)
    block = _intersectionareas(m, _destinationtree(dst_cache),
        subtree(src_space, src_inds), op)

    return WeightBlock(block, _blockdenom(block, ndst))
end

# A method that reads no prepared geometry still takes the destination a tile
# prepared, and names its cells.
pairblock(method::AbstractRegriddingMethod, dst_space::RegridSpace,
    dst_cache::DestinationCache, src_space::RegridSpace, src_inds) =
    pairblock(method, dst_space, dst_cache.inds, src_space, src_inds)

function _sharedmanifold(dst_space::RegridSpace, src_space::RegridSpace)
    m = manifold(dst_space)
    m == manifold(src_space) || throw(ArgumentError(
        "conservative weights need one manifold on both sides: destination is " *
        "$(m), source is $(manifold(src_space))"))
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
    _intersectionareas(manifold, dst_tree, src_tree, op) -> SparseMatrixCSC

Compute intersection areas, threaded when more than one thread is available
and no outer loop is already parallel ([`OUTER_PARALLEL`](@ref)).
"""
function _intersectionareas(m::GOCore.Manifold, dst_tree, src_tree, op)
    threaded = _innerthreaded()
    dst_tree, src_tree = _task_prepared_intersection_trees(
        threaded, dst_tree, src_tree)
    ValType = ConservativeRegridding.output_eltype(op, src_tree, dst_tree)
    return _with_sparse_assembly_cache(ValType) do cache
        ConservativeRegridding.intersection_areas(
            m, threaded, dst_tree, src_tree; intersection_operator = op, cache)
    end
end

# Assembly scratch is module-owned because the block-build call chain has no
# operation context.  The freelist grows on an empty acquire rather than
# waiting for another build, and its retained size is therefore bounded by the
# peak number of overlapping builds for each value type.
struct _SparseAssemblyCacheKey{T} end

mutable struct _SparseAssemblyCachePool{T}
    lock::ReentrantLock
    free::Vector{ConservativeRegridding.SparseMatrixAssemblyCache{T}}
end

_SparseAssemblyCachePool(::Type{T}) where {T} =
    _SparseAssemblyCachePool{T}(ReentrantLock(),
        ConservativeRegridding.SparseMatrixAssemblyCache{T}[])

const _SPARSE_ASSEMBLY_POOL_LOCK = ReentrantLock()
const _SPARSE_ASSEMBLY_POOLS = Dict{DataType,Any}()

function _sparse_assembly_cache_pool(::Type{T}) where {T}
    key = _SparseAssemblyCacheKey{T}
    lock(_SPARSE_ASSEMBLY_POOL_LOCK)
    try
        return get!(_SPARSE_ASSEMBLY_POOLS, key) do
            _SparseAssemblyCachePool(T)
        end::_SparseAssemblyCachePool{T}
    finally
        unlock(_SPARSE_ASSEMBLY_POOL_LOCK)
    end
end

function _acquire_sparse_assembly_cache(pool::_SparseAssemblyCachePool{T}) where {T}
    lock(pool.lock)
    try
        return isempty(pool.free) ?
            ConservativeRegridding.SparseMatrixAssemblyCache(T) : pop!(pool.free)
    finally
        unlock(pool.lock)
    end
end

function _release_sparse_assembly_cache!(pool::_SparseAssemblyCachePool{T},
    cache::ConservativeRegridding.SparseMatrixAssemblyCache{T}) where {T}
    lock(pool.lock)
    try
        push!(pool.free, cache)
    finally
        unlock(pool.lock)
    end
    return cache
end

function _with_sparse_assembly_cache(f, ::Type{T}) where {T}
    pool = _sparse_assembly_cache_pool(T)
    cache = _acquire_sparse_assembly_cache(pool)
    try
        return f(cache)
    finally
        _release_sparse_assembly_cache!(pool, cache)
    end
end

# CR's threaded frontier hands cursor nodes to spawned tasks, so its roots must
# retain safe chart callables. A serial search and assembly stay in this task
# and may hoist each raster chart once on private cursor copies.
_task_prepared_intersection_trees(::GOCore.True, dst_tree, src_tree) =
    (dst_tree, src_tree)
_task_prepared_intersection_trees(::GOCore.False, dst_tree, src_tree) =
    (_task_prepared_raster_tree(dst_tree), _task_prepared_raster_tree(src_tree))

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
