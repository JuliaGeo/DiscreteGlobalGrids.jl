# Whole-space destination tree with per-cell work done once, not once per visit.

"""
    CapCachedTree(node, caps, offset = 0)

Wrap a grid cursor so leaf extents read a cap vector computed once per build:
grid position `p` is `caps[p - offset]`. The dual-tree search otherwise
re-derives each leaf's cap — an inverse projection over its boundary — once per
opposing leaf, and then a second time in `child_indices_extents` for the same
visit.

`offset` is `0` for a whole-space tree and `first(inds) - 1` for a chunk's,
whose cap vector covers only that chunk's positions.
"""
struct CapCachedTree{C<:HierarchicalGridCursor,E}
    node::C
    caps::Vector{E}
    offset::Int
end

CapCachedTree(node::HierarchicalGridCursor, caps::Vector) =
    CapCachedTree(node, caps, 0)

Base.show(io::IO, t::CapCachedTree) =
    print(io, "CapCachedTree(", t.node, ", ", length(t.caps), " caps)")

STI.isspatialtree(::Type{<:CapCachedTree}) = true
STI.node_extent_is_expensive(::Type{<:CapCachedTree{C}}) where {C} =
    STI.node_extent_is_expensive(C)
STI.isleaf(t::CapCachedTree) = STI.isleaf(t.node)
STI.nchild(t::CapCachedTree) = STI.nchild(t.node)
STI.getchild(t::CapCachedTree) =
    Iterators.map(n -> CapCachedTree(n, t.caps, t.offset), STI.getchild(t.node))
STI.getchild(t::CapCachedTree, i::Int) =
    CapCachedTree(STI.getchild(t.node, i), t.caps, t.offset)

# A window node at the leaf level is exactly one stored cell; its position
# indexes the cache. Everything else keeps the cursor's own extent logic.
function STI.node_extent(t::CapCachedTree)
    c = t.node
    if !Engine._issynthetic(c) && c.level >= c.leaf_level
        return @inbounds t.caps[c.first_index - t.offset]
    end
    return STI.node_extent(c)
end

function STI.child_indices_extents(t::CapCachedTree{<:Any,E}) where {E}
    c = t.node
    STI.isleaf(c) ||
        throw(ArgumentError("child_indices_extents is only valid for leaf nodes"))
    count = Engine._stored_count(c)
    entries = Vector{Tuple{Int,E}}(undef, count)
    for k in 1:count
        index = Engine._stored_index(c, k)
        entries[k] = (index, @inbounds t.caps[index - t.offset])
    end
    return entries
end

Trees.ncells(t::CapCachedTree) = Trees.ncells(t.node)
Trees.getcell(t::CapCachedTree, i::Int) = Trees.getcell(t.node, i)
Trees.getcell(t::CapCachedTree) = Trees.getcell(t.node)
GOCore.best_manifold(t::CapCachedTree) = GOCore.best_manifold(t.node)

"""
    _CACHED_BUCKET_SIZE

Leaf size for a cursor whose caps this file precomputes: stop descent at `7^2`
stored cells, two IGeo7 refinement levels above the grid's own resolution. It is
a cell budget, not a level count, so a system of another aperture stops wherever
`49` cells falls — the swept optimum below is flat enough to carry that.

A cursor's default leaf is one cell, so the dual search descends every level and
spends most of its visits on the bottom one — for an IGeo7 level-12 column
rooted at level 5, 823543 of the tree's 960799 nodes. Stopping two levels early
deletes that layer; the leaf then hands back its `49` cells through
`child_indices_extents`, which reads them straight out of `caps`.

Measured on the production CopDEM GLO-90 -> IGeo7 L12 column regrid, core-seconds
per column, single-threaded:

| column         | leaf 1 | leaf 8 | **leaf 49** | leaf 343 |
|:---------------|-------:|-------:|------------:|---------:|
| 728            |  20.39 |  13.16 |   **12.69** |    18.67 |
| 98241 (all NaN)|  18.27 |  11.66 |   **10.97** |    16.60 |
| 115426 (polar) |  20.36 |  13.30 |   **12.34** |    18.89 |

!!! warning "This constant belongs to the cached seam, not to the cursor"
    A big leaf is only cheap because `caps` already holds its cells' extents. A
    bare [`HierarchicalGridCursor`](@ref) re-derives them on every visit, and the
    same sweep against one costs `+25%` at leaf 50 and `+441%` at leaf 350 — the
    change flips sign. So it is applied at the two sites that return a
    `CapCachedTree`, and never on a path that hands back a plain cursor:
    `GR.celltree`, a [`subcursor`](@ref) window, the oversized-chunk return
    below, or the selection-cursor fallback in `_cachedcelltree`.
"""
const _CACHED_BUCKET_SIZE = 49

# Give a cursor the cached seam's leaf size, leaving an explicit caller choice
# alone. `0` is the grid default ("descend to single cells"), not a request.
function _bucketed(c::HierarchicalGridCursor)
    c.bucket_size == 0 || return c
    return typeof(c)(c.grid, c.system, c.top_level, c.leaf_level,
        _CACHED_BUCKET_SIZE, c.level, c.id, c.first_index, c.last_index,
        c.selection)
end

# The whole-space tree, or the plain cursor where the wrap does not apply
# (selection cursors index leaves by selection slot, not grid position).
function _cachedcelltree(space::DGGSpace)
    root = treeify(_decodedgrid(space.grid))
    (root isa HierarchicalGridCursor && root.selection === nothing) ||
        return GR.celltree(space)
    caps = _leafcaps(root.grid, 1:ncells(root.grid))
    return CapCachedTree(_bucketed(root), caps)
end

# Above this a chunk's cap vector costs more to fill than the revisits it saves.
# A chunk at the limit holds 2 MiB of caps.
const _CHUNK_CAP_CACHE_MAX = 1 << 16

# A chunk's cursor with its own cap vector. `subtree` calls this once per block
# build, and the destination side reaches it through `GR.TileCells`, which keeps
# one tree per tile — so the caps are filled once and read by every source chunk
# paired with that tile.
function _cachedchunktree(cursor::HierarchicalGridCursor,
        inds::AbstractUnitRange{<:Integer})
    # Past the limit the cursor goes back bare, so it keeps its own leaf size.
    length(inds) > _CHUNK_CAP_CACHE_MAX && return cursor
    caps = _leafcaps(cursor.grid, inds)
    return CapCachedTree(_bucketed(cursor), caps, Int(first(inds)) - 1)
end

# Decode a compressed id vector once so descent and geometry read O(1) ids.
_decodedgrid(grid::PartialGrid{<:AbstractHierarchicalGridSystem,<:CellVector}) =
    PartialGrid(system(grid), level(grid), collect(grid.ids);
        bucket_size = grid.bucket_size,
        root = Engine._is_rooted(grid) ? grid.root_id : nothing)
_decodedgrid(grid::AbstractGrid) = grid

# One tight cap per cell of `inds`, in parallel at top level ("outer
# parallelism wins"). Entry `k` holds the cap of position `first(inds) + k - 1`.
function _leafcaps(grid::AbstractGrid, inds::AbstractUnitRange{<:Integer})
    n = length(inds)
    n == 0 && return [_cellcap(grid, Int(i)) for i in inds]
    off = Int(first(inds)) - 1
    caps = Vector{typeof(_cellcap(grid, off + 1))}(undef, n)
    nt = GR._innerthreaded() isa GOCore.True ?
        min(Threads.nthreads(), max(1, n >> 14)) : 1
    if nt > 1
        tasks = Vector{Task}(undef, nt)
        for t in 1:nt
            lo = (n * (t - 1)) ÷ nt + 1
            hi = (n * t) ÷ nt
            tasks[t] = Threads.@spawn _fillcaps!(caps, grid, off, lo, hi)
        end
        foreach(wait, tasks)
    else
        _fillcaps!(caps, grid, off, 1, n)
    end
    return caps
end

_cellcap(grid::AbstractGrid, i::Int) = Fallbacks.cell_cap(grid, cellindex(grid, i))

function _fillcaps!(caps::Vector, grid::AbstractGrid, off::Int, lo::Int, hi::Int)
    for k in lo:hi
        @inbounds caps[k] = _cellcap(grid, off + k)
    end
    return nothing
end
