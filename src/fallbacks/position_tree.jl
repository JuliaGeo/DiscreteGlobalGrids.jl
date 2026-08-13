# ---------------------------------------------------------------------------
# The fallback tree: position space, built from cell caps
#
# What a grid with no hierarchy gets. There is no parent/child arithmetic to
# lean on, so the tree is built once, eagerly, over the grid's positions:
#
#   1. every cell's tight cap (one `cell_boundary` per cell — the O(n) part),
#   2. positions sorted by a Morton key on the cap centre's lon/lat, so that
#      contiguous blocks of the sorted order are spatially compact whatever the
#      grid's own position order happens to be,
#   3. a 4-ary tree over that order, with each node's extent the merged cap of
#      its children.
#
# Node extents are therefore *stored*, not derived: `node_extent_is_expensive`
# is `false` here, unlike the hierarchical cursor.
#
# Leaf indices are grid positions throughout — `child_indices_extents` yields
# them and `Trees.getcell` takes them — so a `Regridder` built on this tree
# lines up with data laid out in the grid's own order.
# ---------------------------------------------------------------------------

# Cells per leaf block. Small enough that a leaf's cells are genuinely close
# together, large enough that the per-node overhead of the traversal is
# amortised — the same trade `QUERY_BUCKET_SIZE` makes for the cursor.
const POSITION_TREE_LEAF_SIZE = 16
const POSITION_TREE_ARITY = 4

"""
    PositionTree(grid)

A spatial tree over the positions of any [`AbstractGrid`](@ref), built from the
cells' bounding caps. This is what [`treeify`](@ref) returns for a grid with no
hierarchical system; grids that have one get a [`HierarchicalGridCursor`](@ref)
instead, which needs no build step at all.

Construction is O(ncells) in time and memory (one cap per cell). The tree is
immutable and its nodes are addressed by [`PositionTreeNode`](@ref).
"""
struct PositionTree{G<:AbstractGrid}
    grid::G
    order::Vector{Int}              # grid positions, in tree order
    caps::Vector{Cap}               # cap of `order[k]`
    node_first::Vector{Int}         # node -> first index into `order`
    node_last::Vector{Int}
    node_children::Vector{Vector{Int}}
    node_cap::Vector{Cap}
end

"""
    PositionTreeNode(tree, index)

One node of a [`PositionTree`](@ref) — the `SpatialTreeInterface` cursor over
it. `PositionTreeNode(tree, 1)` is the root, which is what `treeify` returns.
"""
struct PositionTreeNode{G<:AbstractGrid}
    tree::PositionTree{G}
    index::Int
end

function PositionTree(grid::AbstractGrid)
    n = ncells(grid)
    caps = Vector{Cap}(undef, n)
    keys = Vector{UInt64}(undef, n)
    for i in 1:n
        cap = cell_cap(grid, cellindex(grid, i))
        caps[i] = cap
        keys[i] = _morton_key(cap.point)
    end
    order = sortperm(keys)
    tree = PositionTree{typeof(grid)}(grid, order, caps[order],
        Int[], Int[], Vector{Int}[], Cap[])
    _build_node!(tree, 1, n)
    return tree
end

# Depth-first build; children are recorded by index because the recursion
# interleaves them (a breadth-first build would make them contiguous and buy
# nothing — the child list is read once per visited node either way).
function _build_node!(tree::PositionTree, lo::Int, hi::Int)
    index = length(tree.node_first) + 1
    push!(tree.node_first, lo)
    push!(tree.node_last, hi)
    push!(tree.node_children, Int[])
    push!(tree.node_cap, full_sphere_cap())
    count = hi - lo + 1
    if count <= POSITION_TREE_LEAF_SIZE
        tree.node_cap[index] = _merge_range(tree.caps, lo, hi)
        return index
    end
    per = cld(count, POSITION_TREE_ARITY)
    start = lo
    cap = nothing
    while start <= hi
        stop = min(hi, start + per - 1)
        child = _build_node!(tree, start, stop)
        push!(tree.node_children[index], child)
        cap = cap === nothing ? tree.node_cap[child] : merge_caps(cap, tree.node_cap[child])
        start = stop + 1
    end
    tree.node_cap[index] = cap === nothing ? full_sphere_cap() : cap
    return index
end

function _merge_range(caps::Vector{Cap}, lo::Int, hi::Int)
    hi >= lo || return full_sphere_cap()
    cap = caps[lo]
    for k in (lo+1):hi
        cap = merge_caps(cap, caps[k])
    end
    return cap
end

# A Morton key on the cap centre's lon/lat, 16 bits per axis. Only the ordering
# matters: contiguous runs of the sorted keys are spatially local, which is what
# makes a block-split tree prune at all on a grid whose own position order is
# arbitrary.
function _morton_key(p)
    lon, lat = lonlat(p)
    x = UInt64(clamp(floor(Int, (lon + 180.0) / 360.0 * 65535), 0, 65535))
    y = UInt64(clamp(floor(Int, (lat + 90.0) / 180.0 * 65535), 0, 65535))
    key = UInt64(0)
    for b in 0:15
        key |= ((x >> b) & 0x1) << (2b)
        key |= ((y >> b) & 0x1) << (2b + 1)
    end
    return key
end

Base.show(io::IO, tree::PositionTree) =
    print(io, "PositionTree(", typeof(tree.grid).name.name, ", ncells=",
        length(tree.order), ", nodes=", length(tree.node_first), ")")

Base.show(io::IO, node::PositionTreeNode) =
    print(io, "PositionTreeNode(node=", node.index, ", ncells=",
        node.tree.node_last[node.index] - node.tree.node_first[node.index] + 1, ")")

# --------------------------------------------------------------------------
# SpatialTreeInterface
# --------------------------------------------------------------------------

STI.isspatialtree(::Type{<:PositionTreeNode}) = true

# Extents are stored, so there is nothing for the dual search to cache.
STI.node_extent_is_expensive(::Type{<:PositionTreeNode}) = false

STI.isleaf(node::PositionTreeNode) = isempty(node.tree.node_children[node.index])
STI.nchild(node::PositionTreeNode) = length(node.tree.node_children[node.index])
STI.getchild(node::PositionTreeNode) =
    (PositionTreeNode(node.tree, i) for i in node.tree.node_children[node.index])
STI.getchild(node::PositionTreeNode, i::Int) =
    PositionTreeNode(node.tree, node.tree.node_children[node.index][i])
STI.node_extent(node::PositionTreeNode) = node.tree.node_cap[node.index]

function STI.child_indices_extents(node::PositionTreeNode)
    STI.isleaf(node) ||
        throw(ArgumentError("child_indices_extents is only valid for leaf nodes"))
    tree = node.tree
    lo, hi = tree.node_first[node.index], tree.node_last[node.index]
    return Tuple{Int,Cap}[(tree.order[k], tree.caps[k]) for k in lo:hi]
end

# --------------------------------------------------------------------------
# ConservativeRegridding.Trees
# --------------------------------------------------------------------------

GOCore.best_manifold(node::PositionTreeNode) = GOCore.best_manifold(node.tree.grid)
GOCore.best_manifold(tree::PositionTree) = GOCore.best_manifold(tree.grid)

# Leaf indices are grid positions here, at every node — see the file header.
Trees.ncells(node::PositionTreeNode) = ncells(node.tree.grid)
Trees.getcell(node::PositionTreeNode, i::Int) = getcell(node.tree.grid, i)
Trees.getcell(node::PositionTreeNode) = getcell(node.tree.grid)

function Trees.should_parallelize(node::PositionTreeNode, ::US.SphericalCap)
    tree = node.tree
    threshold = max(1, length(tree.order) ÷
                       (Threads.nthreads() * PARALLELIZE_CHUNKS_PER_THREAD))
    return tree.node_last[node.index] - tree.node_first[node.index] + 1 <= threshold
end

# --------------------------------------------------------------------------
# treeify
# --------------------------------------------------------------------------

"""
    treeify(grid::AbstractGrid)
    treeify(manifold, grid::AbstractGrid)

A spatial tree over the grid — total on [`AbstractGrid`](@ref).

A grid from a hierarchical system gets a [`HierarchicalGridCursor`](@ref): the
hierarchy already is the tree, so there is nothing to build and construction is
O(1). Any other grid gets a [`PositionTree`](@ref) built from its cell caps in
O(ncells).

The manifold argument carries no information — this package's geometry is on
the unit sphere by construction — and exists so that
`ConservativeRegridding.Regridder`, which always passes one, resolves here.
"""
treeify(grid::AbstractGrid) = treeify(GOCore.best_manifold(grid), grid)
treeify(::GOCore.Manifold, grid::AbstractGrid) = _grid_tree(grid)

_grid_tree(grid::AbstractGrid) = system(grid) === nothing ?
                                 PositionTreeNode(PositionTree(grid), 1) :
                                 HierarchicalGridCursor(grid)

# Idempotent on the trees themselves, as `Trees.treeify` is for its own cursors.
treeify(::GOCore.Manifold, cursor::HierarchicalGridCursor) = cursor
treeify(::GOCore.Manifold, node::PositionTreeNode) = node
treeify(cursor::HierarchicalGridCursor) = cursor
treeify(node::PositionTreeNode) = node
