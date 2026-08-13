# ---------------------------------------------------------------------------
# `HierarchicalGridCursor` — the hierarchy IS the tree
#
# One cursor type for every hierarchical system; systems never define cursors.
# A cursor is simultaneously the tree and a node in it: `treeify` returns the
# root, whose leaf index space is `1:ncells(grid)` — the space
# `ConservativeRegridding` addresses through `Trees.getcell`, and, for a
# `PartialGrid`, the position space of the caller's own id vector.
#
# Two descent modes, distinguished by the `selection` type parameter so that no
# node pays a runtime branch:
#
#   * window mode (`has_sorted_subtrees(sys)`, `selection === nothing`): a node
#     owns a contiguous *position* window. On a complete level grid that window
#     is `descendant_range` itself — positions, not ids, so there is nothing to
#     convert. On a `PartialGrid` it is two binary searches of the stored ids
#     against the child's descendant id bounds: O(log n) per node, no per-node
#     vector.
#   * selection mode (`selection::Vector{Int}`): a node materialises the
#     positions it owns, split from its parent's by one `ancestor` pass. The
#     fallback for a system that does not declare sorted subtrees; the root
#     selection is `1:ncells(grid)` materialised, which is the price of not
#     having the trait.
#
# Node extents delegate to the system's `node_extent` — the covering law is the
# system's to state, and the cursor has nothing to add to it except where a
# node stores a *proper* subset of its subtree. There it may tighten to the
# union cap of what it really owns, for a bounded number of boundary calls.
# ---------------------------------------------------------------------------

"""
    HierarchicalGridCursor(grid; bucket_size = nothing)

`GeometryOps.SpatialTreeInterface` cursor over a grid of a hierarchical
system — build one with [`treeify(grid)`](@ref treeify) rather than by calling
this constructor.

The tree is the system's own hierarchy: a node is one cell, its children are
that cell's [`children`](@ref), and its leaves are the grid's cells beneath it.
A synthetic root one level above [`rootcells`](@ref) stands for the whole
sphere; a [`PartialGrid`](@ref) built over a subtree starts at its own root
cell instead, so descent stays windowed over the chunk.

`bucket_size` stops descent once a node covers that few stored cells (they are
then scanned, with one cap prune first); `nothing` takes the grid's own — a
`PartialGrid` field, otherwise `0`, which descends to single cells.

Node extents are `SphericalCap`s from the system's [`node_extent`](@ref), and
`node_extent_is_expensive` is `true`, so GeometryOps' dual depth-first search
caches a node's child extents rather than re-deriving them per opposing child.

!!! note "What a system has to get right for this to work"
    A grid that reports a system and is *not* a [`PartialGrid`](@ref) is taken
    to be the complete `levelgrid(sys, level(grid))`, so that a child's window
    is [`descendant_range`](@ref) read straight off — positions, no conversion.
    A system whose `levelgrid` returned some other subset with a different
    position order would descend into the wrong cells; that is precisely what
    the position contract of `descendant_range` forbids.

    Window descent also asks for `descendant_range(sys, c, level(c))` at the
    level above the leaves, so the `l == level(c)` case must answer the cell's
    own one-element position range rather than throwing.
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
    # Selection mode materialises the whole position space once, at the root.
    # That is what a system without `has_sorted_subtrees` costs; every system
    # in scope declares it and takes the window path.
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

The grid **positions** this node owns, ascending — the same index space
`Trees.getcell(tree, i)` addresses, and, for a [`PartialGrid`](@ref), position
`i` of the id vector the grid was built from.

This is how a traversal accepts a whole subtree without walking it: for a
window cursor it is the O(1) range, for a selection cursor the materialised
vector.
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

_child_ids(cursor::HierarchicalGridCursor) =
    _issynthetic(cursor) ? rootcells(cursor.system) : children(cursor.system, cursor.id)

# Position window of a child on a COMPLETE level grid: `descendant_range` is
# already in this grid's position space, by its own contract.
_child_window(cursor::HierarchicalGridCursor, child_id) =
    _range_bounds(descendant_range(cursor.system, child_id, cursor.leaf_level))

_range_bounds(r::AbstractUnitRange) = (Int(first(r)), Int(last(r)))

# ... and on a PartialGrid: two binary searches inside the parent's window,
# against the child's descendant id bounds. Ids inside those bounds that are
# not descendants cannot exist (that is the two-sided range contract), so the
# window is exact rather than a superset.
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

Children of a selection-mode node, non-empty ones only. One pass over the
parent's selection buckets every stored position under its ancestor at the
child level, so the whole child row costs `O(selection * log nchild)` rather
than `O(selection * nchild)`.
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
    STI.isleaf(cursor) && return 0
    count = 0
    for child_id in _child_ids(cursor)
        lo, hi = _child_window(cursor, child_id)
        count += hi >= lo
    end
    return count
end

STI.nchild(cursor::SelectionCursor) =
    STI.isleaf(cursor) ? 0 : length(_selection_children(cursor))

STI.getchild(cursor::WindowCursor) = Iterators.filter(_nonempty,
    (_window_child(cursor, child_id) for child_id in _child_ids(cursor)))

STI.getchild(cursor::SelectionCursor) = _selection_children(cursor)

# The indexed form is the interface's convenience accessor, not the traversal
# path (both depth-first searches iterate `getchild(node)`), so walking the
# iterator keeps one definition of "which children exist".
function STI.getchild(cursor::HierarchicalGridCursor, i::Int)
    # A leaf has no children to index, and asking `children(sys, id)` for them
    # at `max_level` would raise the system's own error instead of a BoundsError.
    (i >= 1 && !STI.isleaf(cursor)) || throw(BoundsError(cursor, i))
    seen = 0
    for child in STI.getchild(cursor)
        seen += 1
        seen == i && return child
    end
    throw(BoundsError(cursor, i))
end

# How many stored leaves an internal node will spend `cell_boundary` calls on
# to tighten its extent past the system's own `node_extent`. A constant,
# deliberately, not a fraction of the subtree: it is the per-node work bound,
# and the defect it replaces was a limit that scaled with the subtree instead.
# 64 clears one child row at every aperture in scope (4, 7, 9) and two rows of
# an aperture-7 grid (49); past that a node is not "a handful of stored leaves"
# any more and the system's own bound is the honest one.
const STORED_UNION_CAP_LIMIT = 64

function STI.node_extent(cursor::HierarchicalGridCursor)
    _issynthetic(cursor) && return full_sphere_cap()
    # A node at the leaf level IS one cell, and its tight cap is both sound and
    # strictly better than the subtree-covering `node_extent`.
    cursor.level >= cursor.leaf_level && return cell_cap(cursor.grid, cursor.id)
    count = _stored_count(cursor)
    if 0 < count <= STORED_UNION_CAP_LIMIT && count < _subtree_count(cursor)
        # A node storing a handful of a big cell's leaves — the sparse chunk a
        # partial grid exists for — is bounded far more tightly by a cap around
        # those leaves, and that tightness is pruning power. A node that stores
        # its WHOLE subtree gains at most the covering headroom from the union
        # and pays one `cell_boundary` per leaf for it, which is pure loss and
        # falls on the nodes nearest the leaves — the many.
        return cells_cap(cursor.grid, (_stored_id(cursor, i) for i in 1:count))
    end
    return node_extent(cursor.system, cursor.id)
end

# How many leaves the node's cell has in total, against which `_stored_count`
# decides whether the node is sparse. A complete level grid stores every one of
# them by construction, so the union cap can never win there and the count is
# reported as the stored count itself. Where the subtree size is not O(1) to
# know (no sorted subtrees) the node is assumed sparse: it is a partial grid on
# the slow path, which is exactly the case the union cap exists for.
_subtree_count(cursor::HierarchicalGridCursor) = _stored_count(cursor)

function _subtree_count(cursor::HierarchicalGridCursor{<:PartialGrid})
    has_sorted_subtrees(cursor.system) || return typemax(Int)
    return length(descendant_range(cursor.system, cursor.id, cursor.leaf_level))
end

"""
    STI.child_indices_extents(cursor) -> Vector{Tuple{Int,SphericalCap{Float64}}}

Positions and tight cell caps of a leaf node's cells. Materialised
deliberately: the dual depth-first search re-iterates this once per opposing
leaf, so a lazy generator would recompute every cap — and its boundary — on
each pass.
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

function Trees.getcell(cursor::HierarchicalGridCursor, i::Int)
    count = _stored_count(cursor)
    1 <= i <= count || throw(BoundsError(1:count, i))
    return cell_polygon(cursor.grid, _stored_id(cursor, i))
end

Trees.getcell(cursor::HierarchicalGridCursor) =
    (Trees.getcell(cursor, i) for i in 1:Trees.ncells(cursor))

# Mirrors `Trees.AbstractQuadtreeCursor`'s policy: spawn once a subtree's leaf
# count drops below `total / (nthreads * 32)`, which lands a few hundred tasks
# for the scheduler to balance. Without this method ConservativeRegridding's
# `::Any` fallback applies — a cap-area heuristic, not a work estimate.
const PARALLELIZE_CHUNKS_PER_THREAD = 32

function Trees.should_parallelize(cursor::HierarchicalGridCursor, ::US.SphericalCap)
    threshold = max(1, ncells(cursor.grid) ÷
                       (Threads.nthreads() * PARALLELIZE_CHUNKS_PER_THREAD))
    return _stored_count(cursor) <= threshold
end
