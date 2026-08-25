# A per-task memo for the block cursor's derived node extents.

"""
    MemoBlockCursor(node::BlockCursor)

[`BlockCursor`](@ref) wrapped so [`node_extent`](@ref
GeometryOps.SpatialTreeInterface.node_extent) reads a cached cap instead of
re-deriving `_node_box` → [`_box_cap`](@ref) on every ask.

  - The values are unchanged, bit for bit.
  - The cache is [`ExtentTable`](@ref DiscreteGlobalGrids.Engine.ExtentTable):
    one direct-mapped table per lattice, a few tables per task, so alternating
    between lattices hits in both.
  - Interior nodes only: a leaf's `child_indices_extents` entries come back as
    the bare cursor's [`LeafCells`](@ref DiscreteGlobalGrids.Engine.LeafCells),
    which never reaches the heap.
"""
struct MemoBlockCursor{C<:BlockCursor}
    node::C
end

Base.show(io::IO, t::MemoBlockCursor) = print(io, "Memo", t.node)

# What fixes a node's extent beyond its rectangle: the pad and the box both read
# the system, and the level decides whether a rectangle names tiles or pixels.
@inline _nodeid(c::BlockCursor) = (c.grid, c.sys, c.level)

# The calling task's table for the lattice this cursor addresses.
@inline _taskmemo(c::BlockCursor) = Engine.extent_table(NTuple{9,Int}, _nodeid(c))

# Everything `_node_box` and `_leaf_pad` read beyond the lattice identity the
# memo is keyed on. `origin` and `strategy` change no extent, so they are out.
@inline _nodekey(c::BlockCursor) =
    (c.r0, c.r1, c.q0, c.q1, c.j0, c.j1, c.i0, c.i1, c.inpixels ? 1 : 0)

STI.isspatialtree(::Type{<:MemoBlockCursor}) = true

# A memoized extent is a compare and a load.
STI.node_extent_is_expensive(::Type{<:MemoBlockCursor}) = false

STI.isleaf(t::MemoBlockCursor) = STI.isleaf(t.node)
STI.nchild(t::MemoBlockCursor) = STI.nchild(t.node)
STI.getchild(t::MemoBlockCursor) =
    Iterators.map(MemoBlockCursor, STI.getchild(t.node))
STI.getchild(t::MemoBlockCursor, k::Int) = MemoBlockCursor(STI.getchild(t.node, k))

function STI.node_extent(t::MemoBlockCursor)
    c = t.node
    return Engine.memoized_extent(_taskmemo(c), _nodekey(c)) do
        STI.node_extent(c)
    end
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
