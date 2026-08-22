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
packed leaf back to its requested global cell position.

Structured spaces should return a native restricted cursor before reaching
this fallback. `nodecapacity` is exposed for packing-invariance checks.
"""
struct CellSpaceRTree{S,N,V<:AbstractVector{Int},C<:AbstractVector{<:Cap}}
    space::S
    node::N
    positions::V
    caps::C
    leafcount::Int
    height::Int
    nodecapacity::Int
end

function CellSpaceRTree(space::S, inds; nodecapacity::Int = 16) where {S<:RegridSpace}
    positions = collect(Int, inds)
    isempty(positions) && throw(ArgumentError(
        "cannot build a cell tree from an empty index set"))
    caps = [_packedcellcap(space, i) for i in positions]
    boxes = map(cap -> convert(Extents.Extent, cap), caps)
    packed = FlexibleRTrees.RTree(FlexibleRTrees.HPR(), positions;
        nodecapacity, extents = boxes)
    return CellSpaceRTree(space, packed, positions, caps, length(positions),
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
    return CellSpaceRTree(tree.space, STI.getchild(tree.node, i), tree.positions,
        tree.caps, childcount, childheight, tree.nodecapacity)
end

STI.getchild(tree::CellSpaceRTree) =
    (STI.getchild(tree, i) for i in 1:STI.nchild(tree))

# FlexibleRTrees reports positions in its payload vector. Translate those local
# positions back to stable global cell positions, and expose their exact cap
# rather than the cap's broad-phase XYZ box.
STI.child_indices_extents(tree::CellSpaceRTree) =
    ((@inbounds(tree.positions[i]), @inbounds(tree.caps[i]))
     for (i, _) in STI.child_indices_extents(tree.node))

# `ncells` and `split_weight` describe this node's restricted population;
# matrix assembly separately needs the owning space's complete global domain.
ncells(tree::CellSpaceRTree) = tree.leafcount
Trees.cell_index_count(tree::CellSpaceRTree) = ncells(tree.space)
Trees.split_weight(tree::CellSpaceRTree) = tree.leafcount
getcell(tree::CellSpaceRTree, i::Int) = getcell(tree.space, i)

function _packedcellpositions!(out::Vector{Int}, tree::CellSpaceRTree)
    if STI.isleaf(tree)
        append!(out, first(entry) for entry in STI.child_indices_extents(tree))
    else
        for child in STI.getchild(tree)
            _packedcellpositions!(out, child)
        end
    end
    return out
end

getcell(tree::CellSpaceRTree) =
    (getcell(tree.space, i) for i in _packedcellpositions!(Int[], tree))
GOCore.best_manifold(tree::CellSpaceRTree) = manifold(tree.space)

# Destination-cell geometry cache

# Avoid caching tiles large enough to create excessive temporary geometry.
const _TILE_CELL_CACHE_MAX = 1 << 16

"""
    TileCells(space, inds)

Wrap `space` and cache geometry for `inds` on first subtree access. Concurrent
block builds share the cache. Tiles above $(_TILE_CELL_CACHE_MAX) cells keep
on-demand geometry.
"""
mutable struct TileCells{S<:RegridSpace,I,M} <: RegridSpace
    space::S
    inds::I
    map::M
    initialized::Bool
    "the tile's cached geometry, or `nothing` when caching is off"
    cells::Union{Nothing,Vector}
    "the tile's restricted tree, built once and shared by every block build"
    tree::Any
    lock::ReentrantLock
end

TileCells(space::RegridSpace, inds) =
    TileCells(space, inds, indexmap(inds), false, nothing, nothing, ReentrantLock())

Base.show(io::IO, tc::TileCells) =
    print(io, "TileCells(", tc.space, ", ", length(tc.inds), " cells)")

# Forward the space interface; only the restricted tree uses cached geometry.
ncells(tc::TileCells) = ncells(tc.space)
getcell(tc::TileCells, i::Int) = getcell(tc.space, i)
manifold(tc::TileCells) = manifold(tc.space)
nchunks(tc::TileCells) = nchunks(tc.space)
cellindices(tc::TileCells, chunk::Int) = cellindices(tc.space, chunk)
chunkat(tc::TileCells, i::Integer) = chunkat(tc.space, i)
chunkat(tc::TileCells, p::US.UnitSphericalPoint) = chunkat(tc.space, p)
cellat(tc::TileCells, p::US.UnitSphericalPoint) = cellat(tc.space, p)
cellcentroid(tc::TileCells, i::Int) = cellcentroid(tc.space, i)
hascellchart(tc::TileCells) = hascellchart(tc.space)
celltree(tc::TileCells) = celltree(tc.space)
chunktree(tc::TileCells) = chunktree(tc.space)

"""
    subtree(tc::TileCells, inds)

Return the memoized restricted tree for the tile's exact index set, wrapped
with cached cell geometry where available. Initialization runs once under a
lock; the immutable tree is then shared by concurrent block builds.
"""
function subtree(tc::TileCells, inds)
    inds == tc.inds || return subtree(tc.space, inds)
    tree, cells = _tiletree!(tc)
    cells === nothing && return tree
    return _cachedtree(tree, tc.map, cells)
end

function _tiletree!(tc::TileCells)
    lock(tc.lock)
    try
        if !tc.initialized
            tc.tree = subtree(tc.space, tc.inds)
            tc.cells = length(tc.inds) > _TILE_CELL_CACHE_MAX ? nothing :
                       _synthesizecells(tc.space, tc.inds)
            tc.initialized = true
        end
        return tc.tree, tc.cells
    finally
        unlock(tc.lock)
    end
end

function _synthesizecells(space::RegridSpace, inds)
    cells = [getcell(space, Int(i)) for i in inds]
    # Dynamic dispatch in the clipping loop costs more than the cache saves.
    return isconcretetype(eltype(cells)) ? cells : nothing
end

# Function barrier for the cache field's abstract element type.
_cachedtree(tree::T, map::M, cells::Vector{P}) where {T,M,P} =
    CachedCellTree{T,M,P}(tree, map, cells)

"""
    CachedCellTree(tree, map, cells)

Wrap `tree` so `getcell` uses `cells`. Spatial-tree methods and extents are
unchanged. The immutable wrapper is safe to share across assembly tasks.
"""
struct CachedCellTree{T,M,P}
    tree::T
    map::M
    cells::Vector{P}
end

Base.show(io::IO, t::CachedCellTree) =
    print(io, "CachedCellTree(", t.tree, ", ", length(t.cells), " cells)")

STI.isspatialtree(::Type{<:CachedCellTree{T}}) where {T} = STI.isspatialtree(T)
STI.node_extent_is_expensive(::Type{<:CachedCellTree{T}}) where {T} =
    STI.node_extent_is_expensive(T)
STI.isleaf(t::CachedCellTree) = STI.isleaf(t.tree)
STI.nchild(t::CachedCellTree) = STI.nchild(t.tree)
STI.getchild(t::CachedCellTree) = STI.getchild(t.tree)
STI.getchild(t::CachedCellTree, i::Int) = STI.getchild(t.tree, i)
STI.node_extent(t::CachedCellTree) = STI.node_extent(t.tree)
STI.child_indices_extents(t::CachedCellTree) = STI.child_indices_extents(t.tree)

GOCore.best_manifold(t::CachedCellTree) = GOCore.best_manifold(t.tree)
ncells(t::CachedCellTree) = ncells(t.tree)
Trees.cell_index_count(t::CachedCellTree) = Trees.cell_index_count(t.tree)
getcell(t::CachedCellTree) = getcell(t.tree)

# A CachedCellTree returned from a shared TileCells owns only a safe retained
# cursor. Serial traversal may use a private wrapper copy with a prepared inner
# cursor; the stored tile tree and the original wrapper remain unchanged.
function _task_prepared_raster_tree(t::CachedCellTree)
    tree = _task_prepared_raster_tree(t.tree)
    tree === t.tree && return t
    return CachedCellTree(tree, t.map, t.cells)
end

@inline function getcell(t::CachedCellTree, i::Int)
    k = localindex(t.map, i)
    k == 0 && return getcell(t.tree, i)
    return @inbounds t.cells[k]
end

# Intersection operator

"""
    CellMemo(::Type{P}; slots)

A task-local direct-mapped memo of cell polygons by tree position. Candidate
pairs leave a leaf pairing in position runs — the destination repeats within a
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
    BlockAreaOperator(inner, dstmap, srcmap, srcmemo, dstmemo)

Measure intersections with `inner` and map global tree positions to block-local
indices. Each assembly task receives its own mutable clipping cache and its own
cell memos.
"""
struct BlockAreaOperator{O,DM,SM,MS,MD}
    inner::O
    dstmap::DM
    srcmap::SM
    srcmemo::MS
    dstmemo::MD
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
        op.dstmap, op.srcmap, _fresh(op.srcmemo), _fresh(op.dstmemo))

# The source is the subject and the destination is the clip ring. A block emits
# weights only for the pairs both of its chunks contain.
@inline function (op::BlockAreaOperator)(rows, cols, vals, item, src_tree, dst_tree)
    i1, i2 = item
    row = localindex(op.dstmap, i2)
    col = localindex(op.srcmap, i1)
    (row == 0 || col == 0) && return nothing
    area = op.inner(_memocell(op.srcmemo, src_tree, i1),
        _memocell(op.dstmemo, dst_tree, i2))
    area > 0 || return nothing
    push!(rows, row)
    push!(cols, col)
    push!(vals, area)
    return nothing
end

# Conservative method

"""
    build_weights!(coo, ::Conservative, dst_space, dst_inds, src_space, src_inds)

Append spherical intersection areas for the two chunks. Each destination
denominator accumulates its covered area. A conservative block always carries a
denominator, including when coverage is zero. Source and destination manifolds
must match. Intersection discovery and clipping may run in parallel.
"""
function build_weights!(coo::WeightCOO, ::Conservative,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    isempty(dst_inds) && return coo
    markdenominated!(coo)
    isempty(src_inds) && return coo

    m = manifold(dst_space)
    m == manifold(src_space) || throw(ArgumentError(
        "conservative weights need one manifold on both sides: destination is " *
        "$(m), source is $(manifold(src_space))"))

    op = BlockAreaOperator(_intersectionoperator(m),
        indexmap(dst_inds), indexmap(src_inds),
        _cellmemo(src_space, src_inds), _cellmemo(dst_space, dst_inds))
    block = _intersectionareas(m, subtree(dst_space, dst_inds),
        subtree(src_space, src_inds), op)

    return _fillcoo!(coo, block)
end

"""
    wholeblock(::Conservative, dst_space, src_space) -> WeightBlock

The eager whole-domain block, adopting the assembled sparse matrix directly.
The generic path copies every entry into a [`WeightCOO`](@ref) and builds a
second, identical CSC from it; here only the denominators are read off. The
values, their CSC layout, and the denominator accumulation order are the
generic path's, bit for bit.
"""
function wholeblock(::Conservative, dst_space::RegridSpace, src_space::RegridSpace)
    ndst = Int(ncells(dst_space))
    nsrc = Int(ncells(src_space))
    # Degenerate sides keep the generic path's exact semantics.
    (ndst == 0 || nsrc == 0) &&
        return invoke(wholeblock, Tuple{AbstractRegriddingMethod,RegridSpace,RegridSpace},
            Conservative(), dst_space, src_space)

    m = manifold(dst_space)
    m == manifold(src_space) || throw(ArgumentError(
        "conservative weights need one manifold on both sides: destination is " *
        "$(m), source is $(manifold(src_space))"))

    op = BlockAreaOperator(_intersectionoperator(m),
        indexmap(1:ndst), indexmap(1:nsrc),
        _cellmemo(src_space, 1:nsrc), _cellmemo(dst_space, 1:ndst))
    block = _intersectionareas(m, subtree(dst_space, 1:ndst),
        subtree(src_space, 1:nsrc), op)

    return WeightBlock(block, _blockdenom(block, ndst))
end

# `_fillcoo!`'s denominator pass, without the COO round trip. Assembly types the
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
    return ConservativeRegridding.intersection_areas(
        m, threaded, dst_tree, src_tree; intersection_operator = op)
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
