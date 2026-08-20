# A per-task memo for the block cursor's derived node extents.

"""
    MemoBlockCursor(node::BlockCursor)

[`BlockCursor`](@ref) wrapped so [`node_extent`](@ref
GeometryOps.SpatialTreeInterface.node_extent) consults a per-task direct-mapped
memo instead of re-deriving `_node_box` → [`_box_cap`](@ref) — five `sincosd`
pairs and four spherical distances — every time it is asked.

The dual-tree join derives a node's extent once per opposing node, and a tile's
interior tree is walked again for every source-block build, every destination
column and every worker, for a lattice that never changes. Each search task
walks one spatially coherent subtree pair, so a thousand slots hold the working
set of the nodes it revisits.

Memory is `O(slots)` per task, not `O(tiles)`, so the wrap holds for the whole
26 475-tile lattice with no eviction policy and no ceiling. The memo lives in
`task_local_storage`, so concurrent readers of one shared grid never touch each
other's slots and no lock is taken. Same [`_box_cap`](@ref) values, bit for bit.

Only interior geometry is memoized here: a leaf's per-cell
`child_indices_extents` entries are rebuilt on every call, as for the bare
cursor.
"""
struct MemoBlockCursor{C<:BlockCursor}
    node::C
end

Base.show(io::IO, t::MemoBlockCursor) = print(io, "Memo", t.node)

# A power of two, so the slot is a mask rather than a modulo.
const _MEMO_EXTENT_SLOTS = 1024

# No node has a negative row, so this can never equal a real key.
const _NO_NODE = (-1, -1, -1, -1, -1, -1, -1, -1, -1)

"""
    BlockCursorMemo()

One task's direct-mapped extent slots, plus the lattice they describe. A node's
extent is fixed by its system, its grid level and its rectangle, so the slots
are cleared whenever the task turns to a grid with a different one.
"""
mutable struct BlockCursorMemo
    grid::Any
    sys::Any
    level::Int
    const keys::Vector{NTuple{9,Int}}
    const vals::Vector{Cap}
end

BlockCursorMemo() = BlockCursorMemo(nothing, nothing, -1,
    fill(_NO_NODE, _MEMO_EXTENT_SLOTS), Vector{Cap}(undef, _MEMO_EXTENT_SLOTS))

# The calling task's memo, reset when the task turns to a different lattice.
function _taskmemo(c::BlockCursor)
    memo = get!(BlockCursorMemo, task_local_storage(),
        :_dgg_copdem_cursor_memo)::BlockCursorMemo
    if memo.grid !== c.grid || memo.sys !== c.sys || memo.level != c.level
        fill!(memo.keys, _NO_NODE)
        memo.grid = c.grid
        memo.sys = c.sys
        memo.level = c.level
    end
    return memo
end

# Everything `_node_box` and `_leaf_pad` read beyond the lattice identity the
# memo is keyed on. `origin` and `strategy` change no extent, so they are out.
@inline _nodekey(c::BlockCursor) =
    (c.r0, c.r1, c.q0, c.q1, c.j0, c.j1, c.i0, c.i1, c.inpixels ? 1 : 0)

@inline _nodeslot(key::NTuple{9,Int}, n::Int) = Int(hash(key) & UInt(n - 1)) + 1

STI.isspatialtree(::Type{<:MemoBlockCursor}) = true

# A hit is a compare and a load, so the search should not pay a vector per
# visited node to avoid it.
STI.node_extent_is_expensive(::Type{<:MemoBlockCursor}) = false

STI.isleaf(t::MemoBlockCursor) = STI.isleaf(t.node)
STI.nchild(t::MemoBlockCursor) = STI.nchild(t.node)
STI.getchild(t::MemoBlockCursor) =
    Iterators.map(MemoBlockCursor, STI.getchild(t.node))
STI.getchild(t::MemoBlockCursor, k::Int) = MemoBlockCursor(STI.getchild(t.node, k))

# The slot holds the whole key and compares it, so a collision is a miss that
# overwrites, never a wrong extent.
function STI.node_extent(t::MemoBlockCursor)
    c = t.node
    key = _nodekey(c)
    memo = _taskmemo(c)
    s = _nodeslot(key, _MEMO_EXTENT_SLOTS)
    @inbounds memo.keys[s] == key && return @inbounds memo.vals[s]
    extent = STI.node_extent(c)
    @inbounds memo.keys[s] = key
    @inbounds memo.vals[s] = extent
    return extent
end

STI.child_indices_extents(t::MemoBlockCursor) = STI.child_indices_extents(t.node)

GOCore.best_manifold(t::MemoBlockCursor) = GOCore.best_manifold(t.node)
Trees.ncells(t::MemoBlockCursor) = Trees.ncells(t.node)
Trees.getcell(t::MemoBlockCursor, i::Int) = Trees.getcell(t.node, i)
Trees.getcell(t::MemoBlockCursor) = Trees.getcell(t.node)
Trees.split_weight(t::MemoBlockCursor) = Trees.split_weight(t.node)

DGG.treeify(::GOCore.Manifold, t::MemoBlockCursor) = t
DGG.treeify(t::MemoBlockCursor) = t

# `treeify` and `subcursor` hand the tree to the regridder, so they wrap; the
# `BlockCursor` constructors stay bare for callers that want the raw cursor.
_memoized(::Nothing) = nothing
_memoized(c::BlockCursor) = MemoBlockCursor(c)
