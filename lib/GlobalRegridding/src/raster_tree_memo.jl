# Per-task memos for the lazy raster tree's per-visit geometry.

"""
    MemoRasterTree(tree::RasterCellTree)

[`RasterCellTree`](@ref) wrapped so node extents and leaf `(position, cap)`
entries consult a per-task direct-mapped memo. The dual-tree join re-derives a
node's extent once per opposing node and a leaf's cell caps once per opposing
leaf; each search task walks one spatially coherent subtree pair, so a few
hundred slots serve most revisits with the value already built. Memory is
`O(slots)` per task regardless of raster size, so the wrap holds for chunked
larger-than-memory sources too. Same [`_rectcap`](@ref) and
[`_rastercellcap`](@ref) values, bit for bit.
"""
struct MemoRasterTree{T<:RasterCellTree}
    tree::T
end

Base.show(io::IO, t::MemoRasterTree) = print(io, "Memo", t.tree)

# Slot counts are powers of two; extents are memoized at every level, leaf
# entry vectors only at the bottom, where almost all revisits land.
const _MEMO_EXTENT_SLOTS = 1024
const _MEMO_ENTRY_SLOTS = 256

const _NO_RECT = (0, 0, 0, 0)

mutable struct RasterTreeMemo
    space::Any
    const extkeys::Vector{NTuple{4,Int}}
    const extvals::Vector{Cap}
    const entkeys::Vector{NTuple{4,Int}}
    const entvals::Vector{Vector{Tuple{Int,Cap}}}
end

RasterTreeMemo() = RasterTreeMemo(nothing,
    fill(_NO_RECT, _MEMO_EXTENT_SLOTS), Vector{Cap}(undef, _MEMO_EXTENT_SLOTS),
    fill(_NO_RECT, _MEMO_ENTRY_SLOTS), Vector{Vector{Tuple{Int,Cap}}}(undef, _MEMO_ENTRY_SLOTS))

# The calling task's memo, reset when the task turns to a different raster.
function _taskmemo(space::RasterGrid)
    memo = get!(RasterTreeMemo, task_local_storage(), :_gr_raster_tree_memo)::RasterTreeMemo
    if memo.space !== space
        fill!(memo.extkeys, _NO_RECT)
        fill!(memo.entkeys, _NO_RECT)
        memo.space = space
    end
    return memo
end

@inline _rectslot(key::NTuple{4,Int}, n::Int) = Int(hash(key) & UInt(n - 1)) + 1

@inline _rect(r::RasterCellTree) = (r.ix0, r.ix1, r.iy0, r.iy1)

STI.isspatialtree(::Type{<:MemoRasterTree}) = true
STI.node_extent_is_expensive(::Type{<:MemoRasterTree}) = false
STI.isleaf(t::MemoRasterTree) = STI.isleaf(t.tree)
STI.nchild(t::MemoRasterTree) = STI.nchild(t.tree)
STI.getchild(t::MemoRasterTree) = map(MemoRasterTree, STI.getchild(t.tree))
STI.getchild(t::MemoRasterTree, i::Int) = MemoRasterTree(STI.getchild(t.tree, i))

function STI.node_extent(t::MemoRasterTree)
    r = t.tree
    key = _rect(r)
    memo = _taskmemo(r.space)
    s = _rectslot(key, _MEMO_EXTENT_SLOTS)
    @inbounds memo.extkeys[s] == key && return @inbounds memo.extvals[s]
    extent = STI.node_extent(r)
    @inbounds memo.extkeys[s] = key
    @inbounds memo.extvals[s] = extent
    return extent
end

# A miss collects the lazy generator in its own order; the vector is never
# mutated after the fill, so handing the slot's reference out stays safe.
function STI.child_indices_extents(t::MemoRasterTree)
    r = t.tree
    key = _rect(r)
    memo = _taskmemo(r.space)
    s = _rectslot(key, _MEMO_ENTRY_SLOTS)
    @inbounds memo.entkeys[s] == key && return @inbounds memo.entvals[s]
    entries = vec(collect(Tuple{Int,Cap}, STI.child_indices_extents(r)))
    @inbounds memo.entkeys[s] = key
    @inbounds memo.entvals[s] = entries
    return entries
end

GOCore.best_manifold(t::MemoRasterTree) = GOCore.best_manifold(t.tree)
Trees.ncells(t::MemoRasterTree) = Trees.ncells(t.tree)
Trees.split_weight(t::MemoRasterTree) = Trees.split_weight(t.tree)
Trees.getcell(t::MemoRasterTree, i::Int) = Trees.getcell(t.tree, i)
Trees.getcell(t::MemoRasterTree) = Trees.getcell(t.tree)
