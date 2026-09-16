# Eager fallback tree for grids without a hierarchy: a packed STR R-tree over
# the boxes of the cells' caps. Leaf indices remain grid indices.

const XYZExtent = Extents.Extent{(:X, :Y, :Z),NTuple{3,Tuple{Float64,Float64}}}

# Cells per leaf, the same trade `QUERY_BUCKET_SIZE` makes for the cursor.
const INDEX_TREE_LEAF_SIZE = 16

const CellTree = GO.FlexibleRTrees.RTree{GO.FlexibleRTrees.STR,XYZExtent,Base.OneTo{Int},Vector{Int}}
const CellNode = GO.FlexibleRTrees.RTreeNode{CellTree,XYZExtent}

"""
    IndexTree(grid)

A spatial tree over grid indices, built from cell caps in `O(ncells)` time and
memory. [`treeify`](@ref ConservativeRegridding.Trees.treeify) uses it only for grids without a hierarchy.

- `caps[i]` is the cap of grid index `i`; `rtree` is a
  `GeometryOps.FlexibleRTrees.RTree` over their boxes in STR order, and leaf
  slot `k` holds grid index `rtree.indices[k]`.
- [`node_extent`](@ref) is the cap bounding a node's box and a leaf lists its
  cells with their stored caps, so every consumer sees caps.
"""
struct IndexTree{G<:AbstractGrid}
    grid::G
    caps::Vector{Cap}
    rtree::CellTree
end

"""
    IndexTreeNode(tree)
    IndexTreeNode(tree, node)

One node of an [`IndexTree`](@ref) — the `SpatialTreeInterface` cursor over
it. `IndexTreeNode(tree)` is the root, which is what `treeify` returns.
"""
struct IndexTreeNode{G<:AbstractGrid}
    tree::IndexTree{G}
    node::CellNode
end

IndexTreeNode(tree::IndexTree) = IndexTreeNode(tree, _root_node(tree.rtree))

function IndexTree(grid::AbstractGrid)
    n = ncells(grid)
    caps = Cap[cell_cap(grid, cellindex(grid, i)) for i in 1:n]
    rtree = n == 0 ? _empty_cell_tree() :
            GO.FlexibleRTrees.RTree(GO.FlexibleRTrees.STR(), Base.OneTo(n);
                nodecapacity=INDEX_TREE_LEAF_SIZE,
                extents=XYZExtent[convert(Extents.Extent, cap) for cap in caps])
    return IndexTree{typeof(grid)}(grid, caps, rtree)
end

# One empty leaf: `RTree` refuses an empty collection, and a search over an
# empty grid must find nothing rather than throw.
_empty_cell_tree() = CellTree(GO.FlexibleRTrees.STR(), INDEX_TREE_LEAF_SIZE,
    Extents.Extent(X=(0.0, 0.0), Y=(0.0, 0.0), Z=(1.0, 1.0)), [XYZExtent[]], Int[], Base.OneTo(0))

_root_node(tree::GO.FlexibleRTrees.RTree) =
    GO.FlexibleRTrees.RTreeNode(tree, 0, 1, Extents.extent(tree))

# The leaf slots under an R-tree node: packing unions consecutive runs of the
# capacity, so every subtree is one contiguous run of slots.
function _leaf_range(node::GO.FlexibleRTrees.RTreeNode)
    tree = node.tree
    span = tree.nodecapacity^(length(tree.levels) - node.level)
    return ((node.index - 1) * span + 1):min(node.index * span, length(tree.indices))
end

@inline function _box_half_diagonal(ext::XYZExtent)
    dx = ext.X[2] - ext.X[1]
    dy = ext.Y[2] - ext.Y[1]
    dz = ext.Z[2] - ext.Z[1]
    return 0.5 * sqrt(dx * dx + dy * dy + dz * dz)
end

# A box's part of the sphere lies within its half diagonal of the box centre,
# plus that centre's distance from the sphere along the centre's direction.
function _box_cap(ext::XYZExtent)
    cx = (ext.X[1] + ext.X[2]) / 2
    cy = (ext.Y[1] + ext.Y[2]) / 2
    cz = (ext.Z[1] + ext.Z[2]) / 2
    len = sqrt(cx * cx + cy * cy + cz * cz)
    chord = _box_half_diagonal(ext) + abs(1 - len)
    (len > 0 && chord < 2) || return full_sphere_cap()
    return SphericalCap(USPoint(cx / len, cy / len, cz / len), 2 * asin(chord / 2))
end

Base.show(io::IO, tree::IndexTree) =
    print(io, "IndexTree(", typeof(tree.grid).name.name, ", ncells=",
        length(tree.rtree.indices), ", levels=", length(tree.rtree.levels), ")")

Base.show(io::IO, node::IndexTreeNode) =
    print(io, "IndexTreeNode(level=", node.node.level, ", index=", node.node.index,
        ", ncells=", length(_leaf_range(node.node)), ")")

# --------------------------------------------------------------------------
# SpatialTreeInterface
# --------------------------------------------------------------------------

STI.isspatialtree(::Type{<:IndexTreeNode}) = true

# A node cap is a few flops from its stored box; nothing for the dual search to cache.
STI.node_extent_is_expensive(::Type{<:IndexTreeNode}) = false

STI.isleaf(node::IndexTreeNode) = STI.isleaf(node.node)
STI.nchild(node::IndexTreeNode) = STI.nchild(node.node)
STI.getchild(node::IndexTreeNode) =
    (IndexTreeNode(node.tree, child) for child in STI.getchild(node.node))
STI.getchild(node::IndexTreeNode, i::Int) = IndexTreeNode(node.tree, STI.getchild(node.node, i))
STI.node_extent(node::IndexTreeNode) = _box_cap(node.node.extent)
function STI.child_indices_extents(node::IndexTreeNode)
    STI.isleaf(node) ||
        throw(ArgumentError("child_indices_extents is only valid for leaf nodes"))
    caps = node.tree.caps
    return ((index, caps[index]) for (index, _) in STI.child_indices_extents(node.node))
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
# would be wrong here; the node's leaf run is exact.
Trees.split_weight(node::IndexTreeNode) = length(_leaf_range(node.node))

# --------------------------------------------------------------------------
# treeify
# --------------------------------------------------------------------------

"""
    treeify(grid::AbstractGrid)
    treeify(manifold, grid::AbstractGrid)

Return a spatial tree for any [`AbstractGrid`](@ref). Hierarchical grids receive
an `O(1)` [`HierarchicalGridCursor`](@ref); other grids receive an
[`IndexTree`](@ref) built over every cell. The manifold overload supports
callers that pass the grid's unit-sphere manifold.
"""
treeify(grid::AbstractGrid) = treeify(GOCore.best_manifold(grid), grid)
treeify(::GOCore.Manifold, grid::AbstractGrid) = _grid_tree(grid)

_grid_tree(grid::AbstractGrid) = system(grid) === nothing ?
                                 IndexTreeNode(IndexTree(grid)) :
                                 HierarchicalGridCursor(grid)

# Idempotent on the trees themselves, as `Trees.treeify` is for its own cursors.
treeify(::GOCore.Manifold, cursor::HierarchicalGridCursor) = cursor
treeify(::GOCore.Manifold, node::IndexTreeNode) = node
treeify(cursor::HierarchicalGridCursor) = cursor
treeify(node::IndexTreeNode) = node
