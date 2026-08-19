# A hierarchical cursor uses grid positions as leaf indices. Sorted-subtree
# systems descend through contiguous windows; others materialize selections.
# Sparse partial-grid nodes may tighten their system-provided extent.

"""
    HierarchicalGridCursor(grid; bucket_size = nothing)

`GeometryOps.SpatialTreeInterface` cursor over a hierarchical grid. Prefer
[`treeify(grid)`](@ref treeify) to direct construction.

The tree is the system's own hierarchy: a node is one cell, its children are
that cell's [`children`](@ref), and its leaves are the grid's cells beneath it.
A synthetic root above [`rootcells`](@ref) covers the sphere. A rooted
[`PartialGrid`](@ref) starts at its stored root.

`bucket_size` stops descent when a node has at most that many stored cells;
`nothing` uses the grid default.

Node extents are system [`node_extent`](@ref) caps and are marked expensive.

!!! note "What a system has to get right for this to work"
    A non-`PartialGrid` grid with a system must be its complete level grid in
    canonical position order. `descendant_range(sys, c, level(c))` must return
    the cell's one-element position range.
"""
struct HierarchicalGridCursor{G<:AbstractGrid,S<:AbstractHierarchicalGridSystem,ID,X}
    grid::G
    system::S
    top_level::Int       # `first(levels(sys))`; `level < top_level` is the synthetic root
    leaf_level::Int
    bucket_size::Int
    level::Int
    id::ID
    first_index::Int
    last_index::Int
    selection::X
end

function HierarchicalGridCursor(grid::AbstractGrid; bucket_size::Union{Nothing,Integer}=nothing)
    sys = system(grid)
    sys === nothing && throw(ArgumentError(
        "$(typeof(grid)) has no hierarchical system; treeify builds a position-space tree for it"))
    leaf = level(grid)
    leaf === nothing && throw(ArgumentError(
        "$(typeof(grid)) reports a system but no level, so it names no leaf level"))
    top = first(levels(sys))
    bs = bucket_size === nothing ? _grid_bucket_size(grid) : Int(bucket_size)
    bs >= 0 || throw(ArgumentError("bucket_size must be non-negative"))
    n = ncells(grid)
    lvl, id = _tree_root(grid, sys, top)
    # Selection mode materializes the root position space once.
    selection = has_sorted_subtrees(sys) ? nothing : collect(1:n)
    return HierarchicalGridCursor{typeof(grid),typeof(sys),typeof(id),typeof(selection)}(
        grid, sys, top, Int(leaf), bs, lvl, id, 1, n, selection)
end

_grid_bucket_size(grid::PartialGrid) = grid.bucket_size
_grid_bucket_size(::AbstractGrid) = 0

# Where descent starts: the chunk's own root where the grid names one, the
# synthetic whole-sphere root otherwise.
function _tree_root(grid::PartialGrid, sys, top)
    _is_rooted(grid) && return (grid.root_level, grid.root_id)
    return (top - 1, _placeholder_root(sys))
end
_tree_root(::AbstractGrid, sys, top) = (top - 1, _placeholder_root(sys))

const WindowCursor{G,S,ID} = HierarchicalGridCursor{G,S,ID,Nothing}
const SelectionCursor{G,S,ID} = HierarchicalGridCursor{G,S,ID,Vector{Int}}

# The synthetic root stands for the whole sphere and its `id` is a placeholder:
# check the level, never the id.
_issynthetic(cursor::HierarchicalGridCursor) = cursor.level < cursor.top_level

_rebuild(cursor::HierarchicalGridCursor{G,S,ID,X}, level::Int, id, lo::Int, hi::Int,
    selection::X) where {G,S,ID,X} =
    HierarchicalGridCursor{G,S,ID,X}(cursor.grid, cursor.system, cursor.top_level,
        cursor.leaf_level, cursor.bucket_size, level, id, lo, hi, selection)

function Base.show(io::IO, cursor::HierarchicalGridCursor)
    print(io, "HierarchicalGridCursor(", typeof(cursor.system).name.name, ", level=")
    _issynthetic(cursor) ? print(io, "root") : print(io, cursor.level, ", id=", cursor.id)
    print(io, ", ncells=", _stored_count(cursor), ")")
end

Base.show(io::IO, ::MIME"text/plain", cursor::HierarchicalGridCursor) = show(io, cursor)

# --------------------------------------------------------------------------
# Node accessors — every one O(1); the dual search calls them per visited node
# --------------------------------------------------------------------------

_stored_count(cursor::WindowCursor) = max(0, cursor.last_index - cursor.first_index + 1)
_stored_count(cursor::SelectionCursor) = length(cursor.selection)

_stored_index(cursor::WindowCursor, i::Int) = cursor.first_index + i - 1
_stored_index(cursor::SelectionCursor, i::Int) = cursor.selection[i]

_stored_id(cursor::HierarchicalGridCursor, i::Int) =
    cellindex(cursor.grid, _stored_index(cursor, i))

"""
    node_indices(cursor) -> AbstractVector{Int}

Return the ascending grid positions owned by this node. A window cursor returns
an `O(1)` range; a selection cursor returns its materialized vector.
"""
node_indices(cursor::WindowCursor) = cursor.first_index:cursor.last_index
node_indices(cursor::SelectionCursor) = cursor.selection

"""
    node_cell(cursor) -> Union{AbstractCellIndex,Nothing}

The cell this node stands for, or `nothing` at the synthetic whole-sphere root.
"""
node_cell(cursor::HierarchicalGridCursor) = _issynthetic(cursor) ? nothing : cursor.id

# --------------------------------------------------------------------------
# Child enumeration
# --------------------------------------------------------------------------

# Treat bucketed and maximum-level nodes uniformly as leaves.
function _child_ids(cursor::HierarchicalGridCursor)
    STI.isleaf(cursor) && return cellindextype(cursor.system)[]
    _issynthetic(cursor) && return rootcells(cursor.system)
    return children(cursor.system, cursor.id)
end

# Position window of a child on a COMPLETE level grid: `descendant_range` is
# already in this grid's position space, by its own contract.
_child_window(cursor::HierarchicalGridCursor, child_id) =
    _range_bounds(descendant_range(cursor.system, child_id, cursor.leaf_level))

_range_bounds(r::AbstractUnitRange) = (Int(first(r)), Int(last(r)))

# A partial-grid child window is the exact intersection with its descendant-id
# bounds, found by two binary searches.
function _child_window(cursor::HierarchicalGridCursor{<:PartialGrid}, child_id)
    cursor.last_index >= cursor.first_index || return (cursor.first_index, cursor.first_index - 1)
    range = descendant_range(cursor.system, child_id, cursor.leaf_level)
    complete = cursor.grid.complete
    lo_id = cellindex(complete, first(range))
    hi_id = cellindex(complete, last(range))
    ids = cursor.grid.ids
    lo = searchsortedfirst(ids, lo_id, cursor.first_index, cursor.last_index, Base.Order.Forward)
    hi = searchsortedlast(ids, hi_id, cursor.first_index, cursor.last_index, Base.Order.Forward)
    return (Int(lo), Int(hi))
end

function _window_child(cursor::WindowCursor, child_id)
    lo, hi = _child_window(cursor, child_id)
    return _rebuild(cursor, cursor.level + 1, child_id, lo, hi, nothing)
end

_nonempty(cursor::HierarchicalGridCursor) = _stored_count(cursor) > 0

"""
    _selection_children(cursor) -> Vector{<:HierarchicalGridCursor}

Return nonempty children of a selection-mode node. Bucketing by child ancestor
costs `O(selection * log nchild)`.
"""
function _selection_children(cursor::SelectionCursor)
    child_level = cursor.level + 1
    child_ids = collect(_child_ids(cursor))
    buckets = [Int[] for _ in eachindex(child_ids)]
    for index in cursor.selection
        c = cellindex(cursor.grid, index)
        anc = child_level == cursor.leaf_level ? c : ancestor(cursor.system, c, child_level)
        slot = searchsortedfirst(child_ids, anc)
        (slot <= length(child_ids) && child_ids[slot] == anc) || continue
        push!(buckets[slot], index)
    end
    out = typeof(cursor)[]
    for slot in eachindex(buckets)
        isempty(buckets[slot]) && continue
        push!(out, _rebuild(cursor, child_level, child_ids[slot],
            1, length(buckets[slot]), buckets[slot]))
    end
    return out
end

# --------------------------------------------------------------------------
# SpatialTreeInterface
# --------------------------------------------------------------------------

STI.isspatialtree(::Type{<:HierarchicalGridCursor}) = true

# Every extent below is derived from `cell_boundary` — an inverse projection
# per cell — rather than read off the node, so the dual depth-first search
# should cache a node's child extents instead of re-deriving them per opposing
# child.
STI.node_extent_is_expensive(::Type{<:HierarchicalGridCursor}) = true

function STI.isleaf(cursor::HierarchicalGridCursor)
    count = _stored_count(cursor)
    count == 0 && return true
    cursor.level >= cursor.leaf_level && return true
    return cursor.bucket_size > 0 && count <= cursor.bucket_size
end

function STI.nchild(cursor::WindowCursor)
    count = 0
    for child_id in _child_ids(cursor)
        lo, hi = _child_window(cursor, child_id)
        count += hi >= lo
    end
    return count
end

STI.nchild(cursor::SelectionCursor) = length(_selection_children(cursor))

STI.getchild(cursor::WindowCursor) = Iterators.filter(_nonempty,
    (_window_child(cursor, child_id) for child_id in _child_ids(cursor)))

STI.getchild(cursor::SelectionCursor) = _selection_children(cursor)

# The indexed form is the interface's convenience accessor, not the traversal
# path (both depth-first searches iterate `getchild(node)`), so walking the
# iterator keeps one definition of "which children exist".
function STI.getchild(cursor::HierarchicalGridCursor, i::Int)
    # A leaf has no children to index, and asking `children(sys, id)` for them
    # at `maxlevel` would raise the system's own error instead of a BoundsError.
    (i >= 1 && !STI.isleaf(cursor)) || throw(BoundsError(cursor, i))
    seen = 0
    for child in STI.getchild(cursor)
        seen += 1
        seen == i && return child
    end
    throw(BoundsError(cursor, i))
end

# Maximum per-node boundary calls used to tighten a sparse node extent. 64
# clears one child row at every aperture in scope (4, 7, 9) and two rows of an
# aperture-7 grid (49).
const STORED_UNION_CAP_LIMIT = 64

function STI.node_extent(cursor::HierarchicalGridCursor)
    _issynthetic(cursor) && return full_sphere_cap()
    # A node at the leaf level IS one cell, and its tight cap is both sound and
    # strictly better than the subtree-covering `node_extent`.
    cursor.level >= cursor.leaf_level && return cell_cap(cursor.grid, cursor.id)
    count = _stored_count(cursor)
    if 0 < count <= STORED_UNION_CAP_LIMIT && count < _subtree_count(cursor)
        # Tighten only proper sparse subsets; complete subtrees use the system cap.
        return cells_cap(cursor.grid, (_stored_id(cursor, i) for i in 1:count))
    end
    return node_extent(cursor.system, cursor.id)
end

# Complete grids report their stored count; partial grids use the descendant
# range when available and otherwise conservatively qualify as sparse.
_subtree_count(cursor::HierarchicalGridCursor) = _stored_count(cursor)

function _subtree_count(cursor::HierarchicalGridCursor{<:PartialGrid})
    has_sorted_subtrees(cursor.system) || return typemax(Int)
    return length(descendant_range(cursor.system, cursor.id, cursor.leaf_level))
end

"""
    STI.child_indices_extents(cursor) -> Vector{Tuple{Int,SphericalCap{Float64}}}

Return leaf grid positions and tight cell caps. The result is materialized to
avoid recomputing caps during repeated dual-tree passes.
"""
function STI.child_indices_extents(cursor::HierarchicalGridCursor)
    STI.isleaf(cursor) ||
        throw(ArgumentError("child_indices_extents is only valid for leaf nodes"))
    count = _stored_count(cursor)
    entries = Vector{Tuple{Int,Cap}}(undef, count)
    for i in 1:count
        index = _stored_index(cursor, i)
        entries[i] = (index, cell_cap(cursor.grid, cellindex(cursor.grid, index)))
    end
    return entries
end

# --------------------------------------------------------------------------
# ConservativeRegridding.Trees
# --------------------------------------------------------------------------

GOCore.best_manifold(cursor::HierarchicalGridCursor) = GOCore.best_manifold(cursor.grid)

# The root cursor's window is the whole grid, so one definition serves both the
# tree-level (`1:ncells(tree)`) and node-level index spaces.
Trees.ncells(cursor::HierarchicalGridCursor) = _stored_count(cursor)

# `i` is a GRID position, at every node — the same index space
# `child_indices_extents` yields, so that the index pairs a dual tree walk
# collects are valid arguments here. (At the root, where every consumer calls
# it, the node-local and grid position spaces coincide anyway.)
function Trees.getcell(cursor::HierarchicalGridCursor, i::Int)
    1 <= i <= ncells(cursor.grid) || throw(BoundsError(cursor.grid, i))
    return cell_polygon(cursor.grid, cellindex(cursor.grid, i))
end

Trees.getcell(cursor::HierarchicalGridCursor) =
    (Trees.getcell(cursor, i) for i in node_indices(cursor))

"""
    PARALLELIZE_CHUNKS_PER_THREAD

Work chunks per thread that `Trees.should_parallelize` aims for: a dual-tree
walk keeps splitting while a node holds more than
`ncells(grid) ÷ (nthreads * PARALLELIZE_CHUNKS_PER_THREAD)` stored leaves.

An internal extension point, re-exported as
`DiscreteGlobalGrids.PARALLELIZE_CHUNKS_PER_THREAD`, so a system writing its own
cursor splits on the same threshold.
"""
const PARALLELIZE_CHUNKS_PER_THREAD = 32

# Mirrors `Trees.AbstractQuadtreeCursor`'s policy. Without this method
# ConservativeRegridding's `::Any` fallback decides on cap area, not work.
function Trees.should_parallelize(cursor::HierarchicalGridCursor, ::US.SphericalCap)
    threshold = max(1, ncells(cursor.grid) ÷
                       (Threads.nthreads() * PARALLELIZE_CHUNKS_PER_THREAD))
    return _stored_count(cursor) <= threshold
end
