# Whole-space destination tree with per-cell work done once, not once per visit.

"""
    CapCachedTree(node, caps)

Wrap a grid cursor so leaf extents read a position-indexed cap vector computed
once per build. The dual-tree search otherwise re-derives each leaf's cap — an
inverse projection over its boundary — once per opposing leaf, and then a
second time in `child_indices_extents` for the same visit.
"""
struct CapCachedTree{C<:HierarchicalGridCursor}
    node::C
    caps::Vector{_Cap}
end

Base.show(io::IO, t::CapCachedTree) =
    print(io, "CapCachedTree(", t.node, ", ", length(t.caps), " caps)")

STI.isspatialtree(::Type{<:CapCachedTree}) = true
STI.node_extent_is_expensive(::Type{CapCachedTree{C}}) where {C} =
    STI.node_extent_is_expensive(C)
STI.isleaf(t::CapCachedTree) = STI.isleaf(t.node)
STI.nchild(t::CapCachedTree) = STI.nchild(t.node)
STI.getchild(t::CapCachedTree) =
    Iterators.map(Base.Fix2(CapCachedTree, t.caps), STI.getchild(t.node))
STI.getchild(t::CapCachedTree, i::Int) = CapCachedTree(STI.getchild(t.node, i), t.caps)

# A window node at the leaf level is exactly one stored cell; its position
# indexes the cache. Everything else keeps the cursor's own extent logic.
function STI.node_extent(t::CapCachedTree)
    c = t.node
    if !Engine._issynthetic(c) && c.level >= c.leaf_level
        return @inbounds t.caps[c.first_index]
    end
    return STI.node_extent(c)
end

function STI.child_indices_extents(t::CapCachedTree)
    c = t.node
    STI.isleaf(c) ||
        throw(ArgumentError("child_indices_extents is only valid for leaf nodes"))
    count = Engine._stored_count(c)
    entries = Vector{Tuple{Int,_Cap}}(undef, count)
    for k in 1:count
        index = Engine._stored_index(c, k)
        entries[k] = (index, @inbounds t.caps[index])
    end
    return entries
end

Trees.ncells(t::CapCachedTree) = Trees.ncells(t.node)
Trees.getcell(t::CapCachedTree, i::Int) = Trees.getcell(t.node, i)
Trees.getcell(t::CapCachedTree) = Trees.getcell(t.node)
GOCore.best_manifold(t::CapCachedTree) = GOCore.best_manifold(t.node)

# The whole-space tree, or the plain cursor where the wrap does not apply
# (selection cursors index leaves by selection slot, not grid position).
function _cachedcelltree(space::DGGSpace)
    root = treeify(_decodedgrid(space.grid))
    (root isa HierarchicalGridCursor && root.selection === nothing) ||
        return GR.celltree(space)
    return CapCachedTree(root, _leafcaps(root.grid))
end

# Decode a compressed id vector once so descent and geometry read O(1) ids.
_decodedgrid(grid::PartialGrid{<:AbstractHierarchicalGridSystem,<:CellVector}) =
    PartialGrid(system(grid), level(grid), collect(grid.ids);
        bucket_size = grid.bucket_size,
        root = Engine._is_rooted(grid) ? grid.root_id : nothing)
_decodedgrid(grid::AbstractGrid) = grid

# One tight cap per cell, in parallel at top level ("outer parallelism wins").
function _leafcaps(grid::AbstractGrid)
    n = ncells(grid)
    caps = Vector{_Cap}(undef, n)
    nt = GR._innerthreaded() isa GOCore.True ?
        min(Threads.nthreads(), max(1, n >> 14)) : 1
    if nt > 1
        tasks = Vector{Task}(undef, nt)
        for t in 1:nt
            lo = (n * (t - 1)) ÷ nt + 1
            hi = (n * t) ÷ nt
            tasks[t] = Threads.@spawn _fillcaps!(caps, grid, lo, hi)
        end
        foreach(wait, tasks)
    else
        _fillcaps!(caps, grid, 1, n)
    end
    return caps
end

function _fillcaps!(caps::Vector{_Cap}, grid::AbstractGrid, lo::Int, hi::Int)
    for i in lo:hi
        @inbounds caps[i] = Fallbacks.cell_cap(grid, cellindex(grid, i))
    end
    return nothing
end
