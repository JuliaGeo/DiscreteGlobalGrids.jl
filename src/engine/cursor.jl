# A hierarchical cursor uses grid indices as leaf indices. Sorted-subtree
# systems descend through contiguous windows; others materialize selections.

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

Node extents are marked expensive. They are system [`node_extent`](@ref) caps,
except at the leaf level, where a node is one cell and reports its own
[`cell_cap`](@ref), which covers nothing below it.

!!! note "What a system has to get right for this to work"
    A non-`PartialGrid` grid with a system must be its complete level grid in
    canonical index order. `descendant_range(sys, c, level(c))` must return
    the cell's one-element index range.
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
    complete_subtree::Bool # every leaf below `id` is stored
    selection::X
end

function HierarchicalGridCursor(grid::AbstractGrid; bucket_size::Union{Nothing,Integer}=nothing)
    sys = system(grid)
    sys === nothing && throw(ArgumentError(
        "$(typeof(grid)) has no hierarchical system; treeify builds an index-space tree for it"))
    leaf = level(grid)
    leaf === nothing && throw(ArgumentError(
        "$(typeof(grid)) reports a system but no level, so it names no leaf level"))
    top = first(levels(sys))
    bs = bucket_size === nothing ? _grid_bucket_size(grid) : Int(bucket_size)
    bs >= 0 || throw(ArgumentError("bucket_size must be non-negative"))
    n = ncells(grid)
    lvl, id = _tree_root(grid, sys, top)
    complete_subtree = lvl >= top && has_sorted_subtrees(sys) &&
        n == length(descendant_range(sys, id, leaf))
    # Selection mode materializes the root index space once.
    selection = has_sorted_subtrees(sys) ? nothing : collect(1:n)
    return HierarchicalGridCursor{typeof(grid),typeof(sys),typeof(id),typeof(selection)}(
        grid, sys, top, Int(leaf), bs, lvl, id, 1, n, complete_subtree, selection)
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
    complete_subtree::Bool, selection::X) where {G,S,ID,X} =
    HierarchicalGridCursor{G,S,ID,X}(cursor.grid, cursor.system, cursor.top_level,
        cursor.leaf_level, cursor.bucket_size, level, id, lo, hi,
        complete_subtree, selection)

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

Return the ascending grid indices owned by this node. A window cursor returns
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

# Index window of a child on a COMPLETE level grid: `descendant_range` is
# already in this grid's index space, by its own contract. Every such child is
# complete as well.
function _child_window_state(cursor::HierarchicalGridCursor, child_id)
    lo, hi = _range_bounds(descendant_range(cursor.system, child_id, cursor.leaf_level))
    return (lo, hi, true)
end

_range_bounds(r::AbstractUnitRange) = (Int(first(r)), Int(last(r)))

# A partial-grid child window is the exact intersection with its descendant-id
# bounds. Once a node is known complete, its child windows are offsets into the
# parent's contiguous descendant range; only boundary nodes need two searches.
function _child_window_state(cursor::HierarchicalGridCursor{<:PartialGrid}, child_id)
    cursor.last_index >= cursor.first_index ||
        return (cursor.first_index, cursor.first_index - 1, false)
    range = descendant_range(cursor.system, child_id, cursor.leaf_level)
    if cursor.complete_subtree
        parent = descendant_range(cursor.system, cursor.id, cursor.leaf_level)
        lo = cursor.first_index + Int(first(range) - first(parent))
        return (lo, lo + length(range) - 1, true)
    end
    complete = cursor.grid.complete
    lo, hi = _partial_child_window(cursor.grid.ids, complete, range,
        cursor.first_index, cursor.last_index)
    full = hi >= lo && hi - lo + 1 == length(range)
    return (Int(lo), Int(hi), full)
end

@inline function _child_window(cursor::HierarchicalGridCursor, child_id)
    lo, hi, _ = _child_window_state(cursor, child_id)
    return (lo, hi)
end

@inline function _partial_child_window(ids, complete, range, first_index::Int, last_index::Int)
    lo_id = cellindex(complete, first(range))
    hi_id = cellindex(complete, last(range))
    lo = searchsortedfirst(ids, lo_id, first_index, last_index, Base.Order.Forward)
    hi = searchsortedlast(ids, hi_id, first_index, last_index, Base.Order.Forward)
    return (lo, hi)
end

function _window_child(cursor::WindowCursor, child_id)
    lo, hi, complete_subtree = _child_window_state(cursor, child_id)
    return _rebuild(cursor, cursor.level + 1, child_id, lo, hi,
        complete_subtree, nothing)
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
            1, length(buckets[slot]), false, buckets[slot]))
    end
    return out
end

# --------------------------------------------------------------------------
# SpatialTreeInterface
# --------------------------------------------------------------------------

STI.isspatialtree(::Type{<:HierarchicalGridCursor}) = true

# Every extent below is computed rather than read from the cursor. Depending on
# the system and node it is a subtree cap, an analytical cell-cap enclosure, or
# a boundary-derived enclosure, so the dual search should carry child extents
# instead of recomputing them for each opposing child.
STI.node_extent_is_expensive(::Type{<:HierarchicalGridCursor}) = true

function STI.isleaf(cursor::HierarchicalGridCursor)
    count = _stored_count(cursor)
    count == 0 && return true
    cursor.level >= cursor.leaf_level && return true
    return cursor.bucket_size > 0 && count <= cursor.bucket_size
end

function STI.nchild(cursor::WindowCursor)
    cursor.complete_subtree && return length(_child_ids(cursor))
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

"""
    STI.node_extent(cursor::HierarchicalGridCursor) -> SphericalCap

The cap this node is pruned by. A node at or below the leaf level is one cell
and bounds that cell's own geometry (`cell_cap`); every other node bounds the
system's whole subtree (`node_extent(system, id)`, the covering law). Only the
latter reaches past the cursor's leaf level, so code that needs a bound over
what lies *below* a leaf — a chunk index, where aperture-7 children overhang
their parent — must call `node_extent(system, id)` itself rather than reuse
the cursor's cap.
"""
function STI.node_extent(cursor::HierarchicalGridCursor)
    _issynthetic(cursor) && return full_sphere_cap()
    # A node at the leaf level is one cell, so its own `cell_cap` is sound: a
    # boundary-derived cap is tighter than the subtree cap, and an analytical
    # one may coincide with it (HEALPix).
    cursor.level >= cursor.leaf_level && return cell_cap(cursor.grid, cursor.id)
    return node_extent(cursor.system, cursor.id)
end

"""
    STI.child_indices_extents(cursor) -> Vector{Tuple{Int,SphericalCap{Float64}}}

The leaf's grid indices and tight cell caps, built once per visit: the pair
search reads a leaf's entries once per opposing cell, so deriving them lazily
would recompute every cap that many times (measured 2.6x slower on HEALPix and
IGeo7, 22x on H3).

Each cap is `cell_cap` of that one leaf cell: a bound over that cell's own
geometry, with no subtree headroom and no coverage of anything below the leaf
level. A consumer that needs the subtree bound must call
`node_extent(system, id)`, whose covering law runs to `maxlevel`.
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
# The tree as a cell source, and the SpatialTreeInterface work verbs
# --------------------------------------------------------------------------

GOCore.best_manifold(cursor::HierarchicalGridCursor) = GOCore.best_manifold(cursor.grid)

# `i` is a GRID index, at every node — the same index space
# `child_indices_extents` yields, so that the index pairs a dual tree walk
# collects are valid arguments here. (At the root, where every consumer calls
# it, the node-local and grid index spaces coincide anyway.)
function getcell(cursor::HierarchicalGridCursor, i::Int)
    1 <= i <= ncells(cursor.grid) || throw(BoundsError(cursor.grid, i))
    return cell_polygon(cursor.grid, cellindex(cursor.grid, i))
end

getcell(cursor::HierarchicalGridCursor) =
    (getcell(cursor, i) for i in node_indices(cursor))

# The node's own work estimate is the cells it stores, and the positions its
# leaves emit are its window (a range) or its selection (a vector): both lazy,
# neither copied.
split_weight(cursor::HierarchicalGridCursor) = _stored_count(cursor)
subtree_indices(cursor::HierarchicalGridCursor) = node_indices(cursor)
# A window's span is the window; a selection stays ascending (children are
# bucketed in the parent's order), so its span is its ends.
GridInterface.position_span(cursor::WindowCursor) = cursor.first_index:cursor.last_index
GridInterface.position_span(cursor::SelectionCursor) =
    isempty(cursor.selection) ? (1:0) : (first(cursor.selection):last(cursor.selection))
