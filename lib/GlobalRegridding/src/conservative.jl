# Conservative weights from spherical intersection areas.

# Chunk-local index maps

# Contiguous chunks use arithmetic; other chunks use a lookup table.

struct OffsetIndexMap
    offset::Int
    n::Int
end

struct LookupIndexMap
    lookup::Dict{Int,Int}
    n::Int
end

indexmap(inds::AbstractUnitRange{<:Integer}) =
    OffsetIndexMap(Int(first(inds)) - 1, length(inds))
indexmap(inds) =
    LookupIndexMap(Dict{Int,Int}(Int(p) => k for (k, p) in enumerate(inds)), length(inds))

@inline localindex(m::OffsetIndexMap, i::Int) = i - m.offset
@inline localindex(m::LookupIndexMap, i::Int) = m.lookup[i]

Base.length(m::OffsetIndexMap) = m.n
Base.length(m::LookupIndexMap) = m.n

# Chunk-restricted cell trees

const _FULL_SPHERE = SphericalCap(USPoint(0.0, 0.0, 1.0), nextfloat(Float64(pi)))

# Maximum cells per fallback-tree leaf.
const _SUBTREE_LEAFSIZE = 16

"""
    subtree(space::RegridSpace, inds) -> tree

Return a spatial tree over `inds`, with leaves addressed by global cell
position. The fallback builds a bounding-cap hierarchy and one polygon per
cell. Spaces with a cheaper restricted tree should specialize this function.
"""
function subtree(space::RegridSpace, inds)
    _iswholespace(space, inds) && return celltree(space)
    return CellCapTree(space, inds)
end

_iswholespace(space::RegridSpace, inds::AbstractUnitRange{<:Integer}) =
    first(inds) == 1 && length(inds) == ncells(space)
_iswholespace(::RegridSpace, inds) = false

"""
    CellCapTree(space, inds)

Build the fallback balanced cap tree for `inds`. Nodes store their extents and
leaves contain at most $(_SUBTREE_LEAFSIZE) cells.
"""
struct CellCapTree{S}
    space::S
    inds::Vector{Int}
    caps::Vector{Cap}
    lo::Int
    hi::Int
    extent::Cap
    children::Vector{CellCapTree{S}}
end

function CellCapTree(space::S, inds) where {S<:RegridSpace}
    ix = collect(Int, inds)
    caps = [_cellcap(space, i) for i in ix]
    return _cellcapnode(space, ix, caps, 1, length(ix))
end

function _cellcapnode(space::S, ix::Vector{Int}, caps::Vector{Cap},
    lo::Int, hi::Int) where {S}
    children = CellCapTree{S}[]
    if hi - lo + 1 > _SUBTREE_LEAFSIZE
        mid = (lo + hi) >> 1
        push!(children, _cellcapnode(space, ix, caps, lo, mid))
        push!(children, _cellcapnode(space, ix, caps, mid + 1, hi))
        extent = _mergecaps([child.extent for child in children])
    else
        extent = _mergecaps(view(caps, lo:hi))
    end
    return CellCapTree{S}(space, ix, caps, lo, hi, extent, children)
end

Base.show(io::IO, tree::CellCapTree) =
    print(io, "CellCapTree(", tree.hi - tree.lo + 1, " cells)")

# Merge caps around their mean centre. Use the full sphere beyond the convex range.
function _mergecaps(caps)
    sx = sy = sz = 0.0
    n = 0
    for c in caps
        sx += c.point[1]
        sy += c.point[2]
        sz += c.point[3]
        n += 1
    end
    n == 0 && return _FULL_SPHERE
    norm = sqrt(sx^2 + sy^2 + sz^2)
    norm <= eps(Float64) && return _FULL_SPHERE
    centre = USPoint(sx / norm, sy / norm, sz / norm)
    radius = 0.0
    for c in caps
        radius = max(radius, US.spherical_distance(centre, c.point) + c.radius)
    end
    radius > Float64(pi) / 2 && return _FULL_SPHERE
    # Nudged outward past dot-product rounding noise, so containment stays closed.
    return SphericalCap(centre, nextfloat(radius * 1.0001 + 1e-12))
end

_cellcap(space::RegridSpace, i::Int) = _mergecaps(
    [SphericalCap(USPoint(GI.x(p), GI.y(p), GI.z(p)), 0.0)
     for p in GI.getpoint(getcell(space, i))])

STI.isspatialtree(::Type{<:CellCapTree}) = true
STI.node_extent_is_expensive(::Type{<:CellCapTree}) = false
STI.isleaf(tree::CellCapTree) = isempty(tree.children)
STI.nchild(tree::CellCapTree) = length(tree.children)
STI.getchild(tree::CellCapTree) = tree.children
STI.getchild(tree::CellCapTree, i::Int) = tree.children[i]
STI.node_extent(tree::CellCapTree) = tree.extent
STI.child_indices_extents(tree::CellCapTree) =
    ((tree.inds[k], tree.caps[k]) for k in tree.lo:tree.hi)

# Tree leaves and cell access both use global space positions.
ncells(tree::CellCapTree) = ncells(tree.space)
getcell(tree::CellCapTree, i::Int) = getcell(tree.space, i)
getcell(tree::CellCapTree) =
    (getcell(tree.space, i) for i in view(tree.inds, tree.lo:tree.hi))
GOCore.best_manifold(tree::CellCapTree) = manifold(tree.space)

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
    # `nothing` before initialization, `false` when caching is disabled.
    cells::Any
    lock::ReentrantLock
end

TileCells(space::RegridSpace, inds) =
    TileCells(space, inds, indexmap(inds), nothing, ReentrantLock())

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

Return the wrapped subtree, using cached cell geometry for the tile's exact
index set. Initialization runs once under a lock.
"""
function subtree(tc::TileCells, inds)
    inds == tc.inds || return subtree(tc.space, inds)
    tree = subtree(tc.space, inds)
    cells = _tilecells!(tc)
    cells === false && return tree
    return _cachedtree(tree, tc.map, cells)
end

function _tilecells!(tc::TileCells)
    lock(tc.lock)
    try
        if tc.cells === nothing
            tc.cells = length(tc.inds) > _TILE_CELL_CACHE_MAX ? false :
                       _synthesizecells(tc.space, tc.inds)
        end
        return tc.cells
    finally
        unlock(tc.lock)
    end
end

function _synthesizecells(space::RegridSpace, inds)
    cells = [getcell(space, Int(i)) for i in inds]
    # Avoid dynamic dispatch in the clipping loop.
    return isconcretetype(eltype(cells)) ? cells : false
end

# Function barrier for the untyped cache field.
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
getcell(t::CachedCellTree) = getcell(t.tree)

@inline function getcell(t::CachedCellTree, i::Int)
    k = _cachedposition(t.map, i)
    k == 0 && return getcell(t.tree, i)
    return @inbounds t.cells[k]
end

@inline _cachedposition(m::OffsetIndexMap, i::Int) =
    (k = i - m.offset; 1 <= k <= m.n ? k : 0)
@inline _cachedposition(m::LookupIndexMap, i::Int) = get(m.lookup, i, 0)

# Intersection operator

"""
    BlockAreaOperator(inner, dstmap, srcmap)

Measure intersections with `inner` and map global tree positions to block-local
indices. Each assembly task receives its own mutable clipping cache.
"""
struct BlockAreaOperator{O,DM,SM}
    inner::O
    dstmap::DM
    srcmap::SM
end

ConservativeRegridding.IntersectionReturnStyle(::BlockAreaOperator) =
    ConservativeRegridding.InPlace()

ConservativeRegridding.output_matrix_size(op::BlockAreaOperator, src_tree, dst_tree) =
    (length(op.dstmap), length(op.srcmap))

ConservativeRegridding.task_local_operator(op::BlockAreaOperator) =
    BlockAreaOperator(ConservativeRegridding.task_local_operator(op.inner),
        op.dstmap, op.srcmap)

# The source is the subject and the destination is the clip ring.
@inline function (op::BlockAreaOperator)(rows, cols, vals, item, src_tree, dst_tree)
    i1, i2 = item
    area = op.inner(Trees.getcell(src_tree, i1), Trees.getcell(dst_tree, i2))
    if area > 0
        push!(rows, localindex(op.dstmap, i2))
        push!(cols, localindex(op.srcmap, i1))
        push!(vals, area)
    end
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
    # Zero coverage still requires a denominator.
    adddenom!(coo, 1, 0.0)
    isempty(src_inds) && return coo

    m = manifold(dst_space)
    m == manifold(src_space) || throw(ArgumentError(
        "conservative weights need one manifold on both sides: destination is " *
        "$(m), source is $(manifold(src_space))"))

    op = BlockAreaOperator(ConservativeRegridding.DefaultIntersectionOperator(m),
        indexmap(dst_inds), indexmap(src_inds))
    block = _intersectionareas(m, subtree(dst_space, dst_inds),
        subtree(src_space, src_inds), op)

    return _fillcoo!(coo, block)
end

# Work around ConservativeRegridding.jl#132 until the pinned version includes it.
# Other failures recur in the serial retry and are rethrown.
_isemptyreduction(err) = err isa ArgumentError &&
                         occursin("reducing over an empty collection", err.msg)

"""
    _intersectionareas(manifold, dst_tree, src_tree, op) -> SparseMatrixCSC

Compute intersection areas using all available threads. When the threaded API
reports an empty task reduction for disjoint trees, retry serially to obtain the
empty block. Other errors are rethrown.
"""
function _intersectionareas(m::GOCore.Manifold, dst_tree, src_tree, op)
    Threads.nthreads() > 1 || return ConservativeRegridding.intersection_areas(
        m, GOCore.False(), dst_tree, src_tree; intersection_operator = op)
    try
        return ConservativeRegridding.intersection_areas(
            m, GOCore.True(), dst_tree, src_tree; intersection_operator = op)
    catch err
        _isemptyreduction(err) || rethrow()
        return ConservativeRegridding.intersection_areas(
            m, GOCore.False(), dst_tree, src_tree; intersection_operator = op)
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
