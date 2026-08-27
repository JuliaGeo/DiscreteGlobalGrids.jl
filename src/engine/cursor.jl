# A hierarchical cursor uses grid indices as leaf indices. Sorted-subtree
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

Node extents are marked expensive. They are system [`node_extent`](@ref) caps,
except at the leaf level and on sparse nodes, which report tight caps over the
cells they actually hold and therefore cover nothing below the leaf level.

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

# Maximum per-node boundary calls used to tighten a sparse node extent. 64
# clears one child row at every aperture in scope (4, 7, 9) and two rows of an
# aperture-7 grid (49).
const STORED_UNION_CAP_LIMIT = 64

"""
    STI.node_extent(cursor::HierarchicalGridCursor) -> SphericalCap

The cap this node is pruned by. It bounds cell *geometry* in every case, but not
the same geometry in every case: a node at or below the leaf level bounds its one
cell (`cell_cap`), a sparse node encloses the stored cells beneath it from either
their boundaries or their analytical caps, and every other node bounds the
system's whole subtree (`node_extent(system, id)`). Only the last reaches past
the cursor's leaf level, so code that needs a bound over what lies *below* a leaf
— a chunk index, where aperture-7 children overhang their parent — must call
`node_extent(system, id)` itself rather than reuse the cursor's cap.
"""
function STI.node_extent(cursor::HierarchicalGridCursor)
    _issynthetic(cursor) && return full_sphere_cap()
    # A node at the leaf level IS one cell, and its tight cap is both sound and
    # strictly better than the subtree-covering `node_extent`.
    cursor.level >= cursor.leaf_level && return cell_cap(cursor.grid, cursor.id)
    count = _stored_count(cursor)
    if 0 < count <= STORED_UNION_CAP_LIMIT && count < _subtree_count(cursor)
        # Tighten only proper sparse subsets; complete subtrees use the system cap.
        return _stored_cells_cap(cursor, count, cell_cap_is_cheap(cursor.grid))
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

# Boundary-derived caps remain the tight generic fallback.  A grid with a cheap
# analytical cell cap can instead enclose those caps directly: no vertex vector,
# no boundary construction, and the cursor still has one concrete extent type.
_stored_cells_cap(cursor, count, ::Val{false}) =
    cells_cap(cursor.grid, (_stored_id(cursor, i) for i in 1:count))

function _stored_cells_cap(cursor, count, ::Val{true})
    return Fallbacks._caps_cap(count) do i
        cell_cap(cursor.grid, _stored_id(cursor, i))
    end
end

"""
    LazyCellCapEntries(cursor)

A read-only, allocation-free leaf view for a grid with a closed-form
`cell_cap`. Each index and cap is derived only when the dual-tree walk reads
that entry; no per-cell cap vector is retained.
"""
struct LazyCellCapEntries{C} <: AbstractVector{Tuple{Int,Cap}}
    cursor::C
end

Base.size(entries::LazyCellCapEntries) = (_stored_count(entries.cursor),)
Base.IndexStyle(::Type{<:LazyCellCapEntries}) = IndexLinear()

Base.@propagate_inbounds function Base.getindex(entries::LazyCellCapEntries, i::Int)
    @boundscheck checkbounds(entries, i)
    cursor = entries.cursor
    index = _stored_index(cursor, i)
    return (index, cell_cap(cursor.grid, cellindex(cursor.grid, index)))
end

"""
    ANALYTICAL_LEAF_CAPACITY

Capacity of an inline analytical leaf. The value is one aperture-7 grandchild
block, matching the broad-search leaf-size optimum when the dual-tree traversal
derives a fixed leaf once and carries it while the opposing tree descends.
"""
const ANALYTICAL_LEAF_CAPACITY = 49

"""
    AnalyticalLeafEntries

One immutable analytical leaf's entries in inline storage. The dual-tree
traversal owns this value only for the fixed leaf's recursive descent, so cap
memory is `O(ANALYTICAL_LEAF_CAPACITY)` rather than `O(ncells(grid))`.
"""
struct AnalyticalLeafEntries <: AbstractVector{Tuple{Int,Cap}}
    entries::NTuple{ANALYTICAL_LEAF_CAPACITY,Tuple{Int,Cap}}
    len::Int
end

Base.size(entries::AnalyticalLeafEntries) = (entries.len,)
Base.IndexStyle(::Type{AnalyticalLeafEntries}) = IndexLinear()

Base.@propagate_inbounds function Base.getindex(entries::AnalyticalLeafEntries, i::Int)
    @boundscheck checkbounds(entries, i)
    return @inbounds entries.entries[i]
end

const _EMPTY_ANALYTICAL_CAP = SphericalCap(USPoint(0.0, 0.0, 1.0), 0.0)
const _EMPTY_ANALYTICAL_ENTRIES = AnalyticalLeafEntries(
    ntuple(_ -> (0, _EMPTY_ANALYTICAL_CAP), Val(ANALYTICAL_LEAF_CAPACITY)), 0)

@inline function _analytical_leaf_entries_indirect(cursor::HierarchicalGridCursor)
    lazy = LazyCellCapEntries(cursor)
    n = length(lazy)
    return if n == 0
        _EMPTY_ANALYTICAL_ENTRIES
    else
        head = lazy[1]
        AnalyticalLeafEntries(ntuple(Val(ANALYTICAL_LEAF_CAPACITY)) do i
            i == 1 ? head : (i <= n ? lazy[i] : head)
        end, n)
    end
end

# A complete partial-grid node occupies one contiguous interval of the complete
# leaf grid. Resolve that interval once, then decode its cells directly. This
# avoids the `CellVector` run search that a logical subset index otherwise pays
# for every entry.
@inline function _analytical_leaf_entries(cursor::HierarchicalGridCursor{<:PartialGrid})
    cursor.complete_subtree || return _analytical_leaf_entries_indirect(cursor)
    range = descendant_range(cursor.system, cursor.id, cursor.leaf_level)
    first_complete_index = Int(first(range))
    n = _stored_count(cursor)
    n == 0 && return _EMPTY_ANALYTICAL_ENTRIES
    first_entry = _direct_analytical_entry(cursor, first_complete_index, 1)
    return AnalyticalLeafEntries(ntuple(Val(ANALYTICAL_LEAF_CAPACITY)) do i
        i == 1 ? first_entry :
        (i <= n ? _direct_analytical_entry(cursor, first_complete_index, i) : first_entry)
    end, n)
end

_analytical_leaf_entries(cursor::HierarchicalGridCursor) =
    _analytical_leaf_entries_indirect(cursor)

@inline function _direct_analytical_entry(cursor, first_complete_index::Int, i::Int)
    index = _stored_index(cursor, i)
    id = cellindex(cursor.grid.complete, first_complete_index + i - 1)
    return (index, cell_cap(cursor.grid, id))
end

function _leaf_indices_extents(cursor::HierarchicalGridCursor, ::Val{false})
    count = _stored_count(cursor)
    entries = Vector{Tuple{Int,Cap}}(undef, count)
    for i in 1:count
        index = _stored_index(cursor, i)
        entries[i] = (index, cell_cap(cursor.grid, cellindex(cursor.grid, index)))
    end
    return entries
end

function _leaf_indices_extents(cursor::HierarchicalGridCursor, ::Val{true})
    _stored_count(cursor) <= ANALYTICAL_LEAF_CAPACITY ||
        return LazyCellCapEntries(cursor)
    return _analytical_leaf_entries(cursor)
end

"""
    STI.child_indices_extents(cursor) -> AbstractVector{Tuple{Int,SphericalCap{Float64}}}

Return leaf grid indices and tight cell caps. Systems whose cap is boundary-
derived materialize the entries once per visit. A system with a closed-form cap
returns an inline [`AnalyticalLeafEntries`](@ref) value. The dual-tree traversal
carries that value while the opposing tree descends, avoiding both repeated cap
construction and retained leaf or grid-sized caches.

Each cap is `cell_cap` of that one leaf cell: a bound over that cell's own
geometry, with no subtree headroom and no coverage of anything below the leaf
level. A consumer that needs the subtree bound must call
`node_extent(system, id)`, whose covering law runs to `maxlevel`.
"""
function STI.child_indices_extents(cursor::HierarchicalGridCursor)
    STI.isleaf(cursor) ||
        throw(ArgumentError("child_indices_extents is only valid for leaf nodes"))
    return _leaf_indices_extents(cursor, cell_cap_is_cheap(cursor.grid))
end

# --------------------------------------------------------------------------
# ConservativeRegridding.Trees
# --------------------------------------------------------------------------

GOCore.best_manifold(cursor::HierarchicalGridCursor) = GOCore.best_manifold(cursor.grid)

# The root cursor's window is the whole grid, so one definition serves both the
# tree-level (`1:ncells(tree)`) and node-level index spaces.
Trees.ncells(cursor::HierarchicalGridCursor) = _stored_count(cursor)

# `i` is a GRID index, at every node — the same index space
# `child_indices_extents` yields, so that the index pairs a dual tree walk
# collects are valid arguments here. (At the root, where every consumer calls
# it, the node-local and grid index spaces coincide anyway.)
function Trees.getcell(cursor::HierarchicalGridCursor, i::Int)
    1 <= i <= ncells(cursor.grid) || throw(BoundsError(cursor.grid, i))
    return cell_polygon(cursor.grid, cellindex(cursor.grid, i))
end

Trees.getcell(cursor::HierarchicalGridCursor) =
    (Trees.getcell(cursor, i) for i in node_indices(cursor))

# No `Trees.split_weight` method: the frontier's default reads `Trees.ncells`,
# which is `_stored_count` here — already the node's own work estimate.
