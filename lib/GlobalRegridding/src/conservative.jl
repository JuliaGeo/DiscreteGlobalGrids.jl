# `build_weights!` for `Conservative()`. Owned by task T3.
#
# The weights are spherical intersection areas, and they come from
# `ConservativeRegridding`: a dual-tree descent over the two chunk-restricted
# cell trees finds the candidate pairs, and its `DefaultIntersectionOperator`
# clips and measures each one. What this file adds is the chunk-local
# addressing the `WeightCOO` seam wants — CR's trees speak cell positions, a
# block speaks positions *within* `dst_inds` and `src_inds` — and the
# restricted trees themselves.

# ===========================================================================
# Chunk-local index maps
# ===========================================================================

# A chunk's cell positions are a `UnitRange` when the space's cell order makes
# the chunk contiguous and an arbitrary ascending vector otherwise, so the
# global→local map has two shapes. Both are concrete and the operator is
# parameterized on them, so the lookup in the assembly loop is not a dynamic
# dispatch.

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

# ===========================================================================
# Chunk-restricted cell trees
# ===========================================================================

const _FULL_SPHERE = SphericalCap(USPoint(0.0, 0.0, 1.0), nextfloat(Float64(pi)))

# Cells per leaf of the fallback tree. Small enough that a leaf-against-leaf
# comparison is cheap, large enough that the tree stays shallow.
const _SUBTREE_LEAFSIZE = 16

"""
    subtree(space::RegridSpace, inds) -> tree

A `SpatialTreeInterface` tree over just the cells `inds` of `space`, with
`SphericalCap` node extents.

Leaf indices and `getcell` indices are **cell positions of `space`**, not
positions within `inds` — the same convention [`celltree`](@ref) uses, which is
why the whole-space case can return `celltree(space)` unchanged. Localization to
block-local indices happens in the weight builder, not here.

The fallback builds a bounding-cap hierarchy over `inds` in the order given,
one cap per cell from [`getcell`](@ref). That costs one polygon per cell and
assumes the space's cell order carries spatial locality, so a space that can
restrict its own tree — a raster chunk, a DGGS subtree — should define a method
of its own rather than inherit this one.
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

A balanced bounding-cap hierarchy over the cells `inds` of `space`, in the
order given — the generic [`subtree`](@ref).

Every node stores its extent, so `node_extent` is free and the dual descent
does not cache. Leaves hold up to $(_SUBTREE_LEAFSIZE) cells and report their
extents as `(cell position, cap)` pairs.
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

# A cap centred on the normalized mean of the given caps' centres, wide enough
# to contain all of them. Past π/2 a cap is no longer convex and no longer
# guaranteed to contain the great-circle arcs between the points it covers, so
# the whole sphere is reported instead — pessimistic, never wrong.
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

# `ncells`/`getcell` here are `ConservativeRegridding.Trees`' own bindings, so
# the tree is addressable as a cell source without a wrapper. Indices are cell
# positions of the space, matching what the leaves emit.
ncells(tree::CellCapTree) = ncells(tree.space)
getcell(tree::CellCapTree, i::Int) = getcell(tree.space, i)
getcell(tree::CellCapTree) =
    (getcell(tree.space, i) for i in view(tree.inds, tree.lo:tree.hi))
GOCore.best_manifold(tree::CellCapTree) = manifold(tree.space)

# ===========================================================================
# The intersection operator
# ===========================================================================

"""
    BlockAreaOperator(inner, dstmap, srcmap)

A `ConservativeRegridding` intersection operator that measures with `inner` and
stores in block-local indices.

CR's default operator stores its result at the tree indices the descent found —
cell positions, which would size the assembled matrix by the whole space. This
one runs the same measurement and relabels through `dstmap`/`srcmap`, so the
matrix it assembles is exactly one chunk pair's block.

The relabelling maps are immutable; the only mutable state is the clipping
cache inside `inner`, which `task_local_operator` hands out per assembly task.
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

# Argument order is CR's: the source cell is the subject, the destination cell
# the clip ring. Keeping it means the spherical clip behaves here exactly as it
# does in `Regridder`, including the convexity precondition on the destination.
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

# ===========================================================================
# The hook
# ===========================================================================

"""
    build_weights!(coo, ::Conservative, dst_space, dst_inds, src_space, src_inds)

Append the spherical area of every nonempty overlap between a destination cell
of `dst_inds` and a source cell of `src_inds`, and return `coo`.

Entry `(j, k)` is the area, on the two spaces' shared [`manifold`](@ref), of
the intersection of destination cell `dst_inds[j]` with source cell
`src_inds[k]`;
`denom[j]` accumulates the same areas, so it is the part of destination cell
`j` this source chunk covers. Summed over the source chunks of one destination
chunk, that is the destination cell's covered area — which is why a partly
covered destination is recoverable and [`Weighted`](@ref) can normalize by it.

Only pairs whose source cell lies in `src_inds` are emitted; the descent runs
over trees restricted to the two chunks ([`subtree`](@ref)), so the partition
invariant holds by construction rather than by filtering. A conservative block
always carries a denominator, including when the two chunks turn out not to
meet at all: zero coverage is an answer.

Both spaces must report the same [`manifold`](@ref); a mismatch is an error
rather than a silent reprojection.

Geometry only — no data, no IO, no missing-value logic — and no state outside
the call, so concurrent builds of different chunk pairs are independent.

One block is built on every thread of the session: the dual-tree descent and
the clip both run in parallel, so a whole-domain plan — which is one block — is
parallel too. The weights do not depend on the thread count; see
`_intersectionareas`.
"""
function build_weights!(coo::WeightCOO, ::Conservative,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    isempty(dst_inds) && return coo
    # Flip the block into carrying a denominator before anything can return
    # early: a conservative block with no coverage still means zero coverage,
    # not "finalize as a raw sum".
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

# The empty-reduction `ArgumentError` `Base.reduce` raises with no `init`, which
# is how CR's threaded dual query reports a chunk pair whose cells do not meet.
# Matched on the message because that is what distinguishes it: the type alone
# is `ArgumentError`, which a genuine argument fault also is.
_isemptyreduction(err) = err isa ArgumentError &&
                         occursin("reducing over an empty collection", err.msg)

"""
    _intersectionareas(manifold, dst_tree, src_tree, op) -> SparseMatrixCSC

`ConservativeRegridding.intersection_areas` over the two trees, on every thread
of the session.

Threading is inside one block build — the dual-tree descent that finds the
candidate pairs and the clip that measures them — so a whole-domain plan, which
is one block, is threaded too. The candidate pairs come back in descent order
either way and no pair is emitted twice, so the assembled matrix is the serial
matrix bit for bit.

**The empty pair.** CR's threaded dual query ends in `reduce(vcat, map(fetch,
tasks))` with no `init`, and `tasks` is empty exactly when no pair of nodes
survives the descent — two chunks that do not meet, which is an answer here and
not a fault. That reduction's `ArgumentError` is caught and the pair rebuilt
serially. The classification is a routing hint and never a correctness
argument: the serial rebuild computes the same block from the same trees, so a
misread exception costs one cheap recomputation, and anything the threaded run
failed on for another reason is raised again by the serial one. It is cheap
because an empty reduction means no descent work was done.
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

# CR's assembled block, read straight into the `WeightCOO`.
#
# Column-major over the stored entries, which is the order `findnz` would have
# produced, so the weights and the accumulated denominators are the same
# `Float64`s to the bit. Reading the matrix in place is what keeps them the same
# *and* keeps the triple `findnz` allocates — three vectors as long as the
# block has nonzeros, live at once with both matrices and the COO — out of the
# build's peak.
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
