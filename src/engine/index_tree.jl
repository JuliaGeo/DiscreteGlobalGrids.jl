# Eager fallback tree for grids without a hierarchy. It sorts tight cell caps by
# Morton key and builds a 4-ary tree with stored extents. Leaf indices remain
# grid indices.

# Cells per leaf block. Small enough that a leaf's cells are genuinely close
# together, large enough that the per-node overhead of the traversal is
# amortised — the same trade `QUERY_BUCKET_SIZE` makes for the cursor.
const INDEX_TREE_LEAF_SIZE = 16
const INDEX_TREE_ARITY = 4

"""
    IndexTree(grid)

A spatial tree over grid indices, built from cell caps in `O(ncells)` time and
memory. [`treeify`](@ref) uses it only for grids without a hierarchy.

Its extents nest — every node's cap is a merge of its children's — because the
tree bottoms out at the grid's own cells and has nothing below them. That is a
property of this tree, not of extents generally; a system's
[`node_extent`](@ref) hierarchy covers descendant geometry down to `maxlevel`
and its caps do not nest.
"""
struct IndexTree{G<:AbstractGrid}
    grid::G
    order::Vector{Int}              # grid indices, in tree order
    caps::Vector{Cap}               # cap of `order[k]`
    node_first::Vector{Int}         # node -> first index into `order`
    node_last::Vector{Int}
    node_children::Vector{Vector{Int}}
    node_cap::Vector{Cap}
end

"""
    IndexTreeNode(tree, index)

One node of a [`IndexTree`](@ref) — the `SpatialTreeInterface` cursor over
it. `IndexTreeNode(tree, 1)` is the root, which is what `treeify` returns.
"""
struct IndexTreeNode{G<:AbstractGrid}
    tree::IndexTree{G}
    index::Int
end

function IndexTree(grid::AbstractGrid)
    n = ncells(grid)
    caps = Vector{Cap}(undef, n)
    keys = Vector{UInt64}(undef, n)
    for i in 1:n
        cap = cell_cap(grid, cellindex(grid, i))
        caps[i] = cap
        keys[i] = _morton_key(cap.point)
    end
    order = sortperm(keys)
    tree = IndexTree{typeof(grid)}(grid, order, caps[order],
        Int[], Int[], Vector{Int}[], Cap[])
    _build_node!(tree, 1, n)
    return tree
end

# Depth-first recursion records children by node index.
function _build_node!(tree::IndexTree, lo::Int, hi::Int)
    index = length(tree.node_first) + 1
    push!(tree.node_first, lo)
    push!(tree.node_last, hi)
    push!(tree.node_children, Int[])
    push!(tree.node_cap, full_sphere_cap())
    count = hi - lo + 1
    if count <= INDEX_TREE_LEAF_SIZE
        tree.node_cap[index] = _merge_range(tree.caps, lo, hi)
        return index
    end
    per = cld(count, INDEX_TREE_ARITY)
    start = lo
    cap = nothing
    while start <= hi
        stop = min(hi, start + per - 1)
        child = _build_node!(tree, start, stop)
        push!(tree.node_children[index], child)
        cap = cap === nothing ? tree.node_cap[child] :
              Extents.union(cap, tree.node_cap[child])
        start = stop + 1
    end
    tree.node_cap[index] = cap === nothing ? full_sphere_cap() : cap
    return index
end

function _merge_range(caps::Vector{Cap}, lo::Int, hi::Int)
    hi >= lo || return full_sphere_cap()
    cap = caps[lo]
    for k in (lo+1):hi
        cap = Extents.union(cap, caps[k])
    end
    return cap
end

# A Morton key on the cap centre's lon/lat, 16 bits per axis. Only the ordering
# matters: contiguous runs of the sorted keys are spatially local, which is what
# makes a block-split tree prune at all on a grid whose own index order is
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

Base.show(io::IO, tree::IndexTree) =
    print(io, "IndexTree(", typeof(tree.grid).name.name, ", ncells=",
        length(tree.order), ", nodes=", length(tree.node_first), ")")

Base.show(io::IO, node::IndexTreeNode) =
    print(io, "IndexTreeNode(node=", node.index, ", ncells=",
        node.tree.node_last[node.index] - node.tree.node_first[node.index] + 1, ")")

# --------------------------------------------------------------------------
# SpatialTreeInterface
# --------------------------------------------------------------------------

STI.isspatialtree(::Type{<:IndexTreeNode}) = true

# Extents are stored, so there is nothing for the dual search to cache.
STI.node_extent_is_expensive(::Type{<:IndexTreeNode}) = false

STI.isleaf(node::IndexTreeNode) = isempty(node.tree.node_children[node.index])
STI.nchild(node::IndexTreeNode) = length(node.tree.node_children[node.index])
STI.getchild(node::IndexTreeNode) =
    (IndexTreeNode(node.tree, i) for i in node.tree.node_children[node.index])
STI.getchild(node::IndexTreeNode, i::Int) =
    IndexTreeNode(node.tree, node.tree.node_children[node.index][i])
STI.node_extent(node::IndexTreeNode) = node.tree.node_cap[node.index]

function STI.child_indices_extents(node::IndexTreeNode)
    STI.isleaf(node) ||
        throw(ArgumentError("child_indices_extents is only valid for leaf nodes"))
    tree = node.tree
    lo, hi = tree.node_first[node.index], tree.node_last[node.index]
    return Tuple{Int,Cap}[(tree.order[k], tree.caps[k]) for k in lo:hi]
end

# --------------------------------------------------------------------------
# ConservativeRegridding.Trees
# --------------------------------------------------------------------------

GOCore.best_manifold(node::IndexTreeNode) = GOCore.best_manifold(node.tree.grid)
GOCore.best_manifold(tree::IndexTree) = GOCore.best_manifold(tree.grid)

# Leaf indices are grid indices here, at every node — see the file header.
Trees.ncells(node::IndexTreeNode) = ncells(node.tree.grid)
Trees.getcell(node::IndexTreeNode, i::Int) = getcell(node.tree.grid, i)
Trees.getcell(node::IndexTreeNode) = getcell(node.tree.grid)

# `Trees.ncells` answers for the whole grid, so the frontier's default estimate
# would be wrong here; the node's stored leaf window is exact.
Trees.split_weight(node::IndexTreeNode) =
    node.tree.node_last[node.index] - node.tree.node_first[node.index] + 1

# --------------------------------------------------------------------------
# treeify
# --------------------------------------------------------------------------

"""
    treeify(grid::AbstractGrid)
    treeify(manifold, grid::AbstractGrid)

Return a spatial tree for any [`AbstractGrid`](@ref). Hierarchical grids receive
an `O(1)` [`HierarchicalGridCursor`](@ref); other grids receive an `O(ncells)`
[`IndexTree`](@ref). The manifold overload supports callers that pass the
grid's unit-sphere manifold.
"""
treeify(grid::AbstractGrid) = treeify(GOCore.best_manifold(grid), grid)
treeify(::GOCore.Manifold, grid::AbstractGrid) = _grid_tree(grid)

_grid_tree(grid::AbstractGrid) = system(grid) === nothing ?
                                 IndexTreeNode(IndexTree(grid), 1) :
                                 HierarchicalGridCursor(grid)

# Idempotent on the trees themselves, as `Trees.treeify` is for its own cursors.
treeify(::GOCore.Manifold, cursor::HierarchicalGridCursor) = cursor
treeify(::GOCore.Manifold, node::IndexTreeNode) = node
treeify(cursor::HierarchicalGridCursor) = cursor
treeify(node::IndexTreeNode) = node
