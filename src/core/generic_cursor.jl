# ---------------------------------------------------------------------------
# Generic DGGS spatial-tree cursor
#
# One cursor type over the grid types in `grid_types.jl`, replacing the four
# per-system tree/node clones this package used to carry. The tree is the
# DGGS hierarchy itself: a synthetic level −1 root stands for the whole sphere,
# a node at `level` for one cell, and its children for that cell's immediate
# children. Leaves are the grid's cells — all of them for a `DGGSGrid`, the
# stored `ids` for a `DGGSPartialGrid`.
#
# A cursor is simultaneously the tree and a node in it: `treeify` returns the
# root cursor, whose leaf index space (`1:ncells`) is the one
# `ConservativeRegridding` addresses through `Trees.getcell`, and every
# descendant addresses the same leaves through a local `1:ncells(node)` window.
# For a partial grid those indices are positions in `grid.ids`, which is what
# makes a `Regridder` line up with the `DimensionalData` lookup the ids came
# from.
#
# Three descent modes, distinguished by type so no node pays a runtime branch:
#
# - dense (`grid::DGGSGrid`): a node owns a contiguous *ordinal* interval.
#   `has_descendant_ranges` makes it O(1) — `descendant_range`'s endpoints are
#   themselves valid descendants, so the interval is
#   `cell_to_ordinal(lo):cell_to_ordinal(hi)` (ordinal monotonicity plus the
#   two-sided range contract). Without the trait, non-root internal nodes carry
#   no interval and fall back to `subtree_leaf_count`/`cell_descendants`.
# - partial + `has_descendant_ranges` (`selection === nothing`): a node stores
#   `(first_index, last_index)` into `grid.ids`, found by two binary searches
#   against `descendant_range` — O(log n) per node, no per-node vector. HEALPix,
#   H3 and IGeo7 all take this path.
# - partial without the trait (`selection::Vector{Int}`): a node materializes
#   the indices its cell owns, split from the parent's selection by one
#   `cell_parent` pass. A5-style fallback.
# ---------------------------------------------------------------------------

"""
    DGGSCursor(grid, level, id, first_index, last_index, selection)

`GeometryOps.SpatialTreeInterface` cursor over a [`DGGSGrid`](@ref) or
[`DGGSPartialGrid`](@ref); build one with [`treeify(grid)`](@ref) rather than
calling this constructor.

`level == -1` is the synthetic whole-sphere root; otherwise the cursor is the
cell `(level, id)`. `first_index:last_index` is the leaf index interval the
node owns (ordinals for a dense grid, positions in `grid.ids` for a partial
one) and `selection` is the materialized index vector of the non-range partial
path — exactly one of the two carries the node's leaf set.

What a traversal has any business reading is the cell the node stands for —
[`node_level`](@ref) and [`node_id`](@ref) — plus the leaf-index set it owns,
through [`node_indices`](@ref). The raw fields behind them are index
bookkeeping whose representation is a dispatch detail (see the three descent
modes at the top of this file).

Construction is guarded: `level` must be in `-1:grid.level` (with `-1` only in
the whole-tree root form), exactly one of the index window and `selection` may
carry the node's leaf set, and the window must lie inside the grid's leaf index
space. See the inner constructor for the full list and the reasons.
"""
struct DGGSCursor{G,ID,X}
    grid::G
    level::Int
    id::ID
    first_index::Int
    last_index::Int
    selection::X
    # This constructor is pure guard: `treeify` and the three descent modes
    # build only states that satisfy it, so nothing on a valid path changes.
    # It exists because `DGGSCursor` is *exported* and its raw six-argument
    # form was open, and every field below is read by the traversal without a
    # second thought:
    #
    #   * `level` is what `node_level` reports and what every kernel call a
    #     traversal makes is keyed on, so a forged level is a wrong answer, not
    #     an error — `-1 <= level <= grid.level`, and `-1` (the synthetic
    #     whole-sphere root) only where the node really does own the whole leaf
    #     set, which is the only shape `treeify` produces.
    #   * the index window is what `Trees.getcell`'s own bounds check is
    #     computed FROM (`_stored_count`), so a window past the leaf count
    #     cannot be caught downstream: it *is* the downstream authority. Hence
    #     `1 <= first_index <= last_index <= total` here, once.
    #   * `selection` is the other representation of the same leaf set, and
    #     `_stored_count` / `_stored_ids` pick between them by *type*: a dense
    #     cursor carrying a selection would silently ignore it, and a selection
    #     whose length disagrees with the window would report one size and
    #     iterate another. Exactly one of the two carries the set.
    #
    # Both empty forms the file documents stay legal: the dense "no interval"
    # marker `(0, -1)` of `_dense_child`, and a partial node's empty window,
    # which the binary searches always return as `last_index == first_index - 1`
    # (`_stored_bounds`' own `(1, 0)` included).
    function DGGSCursor(grid::Union{DGGSGrid,DGGSPartialGrid}, level::Integer, id,
            first_index::Integer, last_index::Integer,
            selection::Union{Nothing,Vector{Int}})
        lvl = Int(level)
        lo = Int(first_index)
        hi = Int(last_index)
        -1 <= lvl <= grid.level || throw(ArgumentError(
            "cursor level must be in -1:$(grid.level) (-1 is the synthetic root), got $lvl"))
        # `length(grid.ids)` for a partial grid, `num_cells` for a dense one —
        # O(1) either way, and the dense count is the number `treeify` already
        # computes for the root window.
        total = _total_cells(grid)
        if selection === nothing
            if hi == lo - 1                     # the two empty/sentinel forms
                dense = grid isa DGGSGrid
                (dense ? lo == 0 : 1 <= lo <= total + 1) || throw(ArgumentError(
                    "empty cursor window $lo:$hi is not one of the documented empty forms"))
            else
                1 <= lo <= hi <= total || throw(ArgumentError(
                    "cursor window $lo:$hi is not inside the grid's leaf indices 1:$total"))
            end
        else
            grid isa DGGSPartialGrid || throw(ArgumentError(
                "only a partial grid's cursor carries a materialized selection"))
            (lo == 1 && hi == length(selection)) || throw(ArgumentError(
                "a selection cursor's window must be 1:$(length(selection)), got $lo:$hi"))
            # One comparison per stored index, against the `O(selection)`
            # bucketing `_child_cursors` does per node anyway. Out-of-range
            # entries would index `grid.ids` out of bounds later.
            all(index -> 1 <= index <= total, selection) || throw(ArgumentError(
                "cursor selection indices must be inside the grid's leaf indices 1:$total"))
        end
        lvl == -1 && hi != total && throw(ArgumentError(
            "the synthetic root (level -1) must own the whole grid, 1:$total"))
        cid = cell_id_type(grid.system)(id)
        return new{typeof(grid),typeof(cid),typeof(selection)}(grid, lvl, cid, lo, hi, selection)
    end
end

"""
    node_level(cursor::DGGSCursor) -> Int

Refinement level of the cell this node stands for, or `-1` at the synthetic
whole-sphere root. Paired with [`node_id`](@ref) it is the `(level, id)` every
kernel operation takes, so a traversal can ask the kernel about the node it is
sitting on:

```julia
STI.isleaf(node) || subtree_cap(system, node_level(node), node_id(node), leaf)
```
"""
node_level(cursor::DGGSCursor) = cursor.level

"""
    node_id(cursor::DGGSCursor) -> id

Canonical cell id of the node, of the system's [`cell_id_type`](@ref). Only
meaningful once `node_level(cursor) >= 0`: the synthetic root stands for the
whole sphere and reports the system's zero id, which is also a real level-0
cell — check the level, not the id.
"""
node_id(cursor::DGGSCursor) = cursor.id

"""
    node_indices(cursor::DGGSCursor) -> AbstractVector{Int}

The leaf indices this node owns, in the *tree's* leaf index space — positions
in `grid.ids` for a partial grid, dense ordinals for a [`DGGSGrid`](@ref) —
so they line up with `Trees.getcell(tree, i)` and, for a lookup-backed grid,
with positions in the `DimensionalData` lookup.

This is the public accessor a traversal uses to accept a whole subtree
without walking it (the alternative, an `STI.depth_first_search` collect,
re-runs the per-node binary searches over every node underneath): for the
range-backed cursors it is the O(1) window `first_index:last_index` —
ascending, and empty ranges for the documented empty node forms — and for a
selection cursor the node's materialized index vector. Only the dense
no-interval fallback (a system without [`has_descendant_ranges`](@ref))
enumerates its leaves.
"""
node_indices(cursor::DGGSCursor{<:DGGSPartialGrid,<:Any,Nothing}) =
    cursor.first_index:cursor.last_index
node_indices(cursor::DGGSCursor{<:DGGSPartialGrid,<:Any,Vector{Int}}) = cursor.selection
function node_indices(cursor::DGGSCursor{<:DGGSGrid})
    _has_interval(cursor) && return cursor.first_index:cursor.last_index
    return Int[first(_leaf_entry(cursor, i)) for i in 1:_stored_count(cursor)]
end

function Base.show(io::IO, cursor::DGGSCursor)
    print(io, "DGGSCursor(", system_name(_system(cursor)), ", level=", cursor.level)
    cursor.level >= 0 && print(io, ", id=", cursor.id)
    print(io, ", ncells=", Trees.ncells(cursor), ")")
end

# --------------------------------------------------------------------------
# Node accessors
#
# Every one of these is O(1): the dual depth-first search calls them once per
# visited node, and `should_parallelize` needs a leaf count per node.
# --------------------------------------------------------------------------

_system(cursor::DGGSCursor) = cursor.grid.system
_leaf_level(cursor::DGGSCursor) = cursor.grid.level
# Only partial grids carry a bucket size — a dense node always splits to a cell.
_bucket_size(cursor::DGGSCursor{<:DGGSPartialGrid}) = cursor.grid.bucket_size

# Dense internal nodes of a system without `descendant_range` carry no ordinal
# interval; `first_index == 0` is the marker (leaf indices are 1-based).
_has_interval(cursor::DGGSCursor) = cursor.first_index >= 1

# Total leaves of the whole tree — `should_parallelize`'s denominator. For a
# subtree-rooted partial grid this is the chunk's cell count, not the globe's.
_total_cells(grid::DGGSGrid) = Int(num_cells(grid.system, grid.level))
_total_cells(grid::DGGSPartialGrid) = length(grid.ids)

function _stored_count(cursor::DGGSCursor{<:DGGSGrid})
    _has_interval(cursor) && return cursor.last_index - cursor.first_index + 1
    return Int(subtree_leaf_count(_system(cursor), cursor.level, cursor.id, _leaf_level(cursor)))
end
_stored_count(cursor::DGGSCursor{<:DGGSPartialGrid,<:Any,Nothing}) =
    max(0, cursor.last_index - cursor.first_index + 1)
_stored_count(cursor::DGGSCursor{<:DGGSPartialGrid,<:Any,Vector{Int}}) =
    length(cursor.selection)

# The node's stored ids, as a view — `cells_cap` only iterates them.
_stored_ids(cursor::DGGSCursor{<:DGGSPartialGrid,<:Any,Nothing}) =
    view(cursor.grid.ids, cursor.first_index:cursor.last_index)
_stored_ids(cursor::DGGSCursor{<:DGGSPartialGrid,<:Any,Vector{Int}}) =
    view(cursor.grid.ids, cursor.selection)

"""
    _leaf_entry(cursor, i) -> (index, id)

The `i`-th leaf of a node (1-based, node-local): its index in the *tree's* leaf
index space and its canonical cell id. At the root the two index spaces
coincide, which is what keeps `Trees.getcell(tree, i)` and node-level indexing
consistent.
"""
function _leaf_entry(cursor::DGGSCursor{<:DGGSGrid}, i::Int)
    system = _system(cursor)
    leaf_level = _leaf_level(cursor)
    if _has_interval(cursor)
        index = cursor.first_index + i - 1
        return (index, ordinal_to_cell(system, leaf_level, index))
    end
    # Fallback for systems without `descendant_range`: only reachable on a
    # dense *internal* node, whose leaves nothing but `Trees.getcell` asks for.
    id = cell_descendants(system, cursor.level, cursor.id, leaf_level)[i]
    return (Int(cell_to_ordinal(system, leaf_level, id)), id)
end

function _leaf_entry(cursor::DGGSCursor{<:DGGSPartialGrid,<:Any,Nothing}, i::Int)
    index = cursor.first_index + i - 1
    return (index, cursor.grid.ids[index])
end

function _leaf_entry(cursor::DGGSCursor{<:DGGSPartialGrid,<:Any,Vector{Int}}, i::Int)
    index = cursor.selection[i]
    return (index, cursor.grid.ids[index])
end

# --------------------------------------------------------------------------
# Child enumeration
# --------------------------------------------------------------------------

"""
    _child_ids(cursor) -> AbstractVector

Canonical ids of the node's immediate children, ascending, empty for a leaf.
Ordinal-id systems get the radix block as a `UnitRange`: the child list of an
ordinal node *is* `id * radix .+ (0:radix-1)` by the definition of
[`has_ordinal_ids`](@ref) — exactly what `cell_children`'s default collects —
so descent allocates nothing per node, as the old per-system HEALPix tree
managed. Structural-id systems go through the wired enumeration.
"""
@inline function _child_ids(cursor::DGGSCursor)
    system = _system(cursor)
    idtype = cell_id_type(system)
    if has_ordinal_ids(system)
        STI.isleaf(cursor) && return one(idtype):zero(idtype)
        cursor.level < 0 && return zero(idtype):(idtype(root_count(system)) - one(idtype))
        base = idtype(cursor.id) * idtype(radix(system))
        return base:(base + idtype(radix(system)) - one(idtype))
    end
    STI.isleaf(cursor) && return Vector{idtype}()
    cursor.level < 0 && return convert(Vector{idtype}, root_ids(system))
    return convert(Vector{idtype}, cell_children(system, cursor.level, cursor.id))
end

function _dense_child(cursor::DGGSCursor{<:DGGSGrid}, child_id)
    system = _system(cursor)
    leaf_level = _leaf_level(cursor)
    child_level = cursor.level + 1
    if has_descendant_ranges(system)
        lo, hi = descendant_range(system, child_level, child_id, leaf_level)
        return DGGSCursor(cursor.grid, child_level, child_id,
            Int(cell_to_ordinal(system, leaf_level, lo)),
            Int(cell_to_ordinal(system, leaf_level, hi)), nothing)
    elseif child_level == leaf_level
        ordinal = Int(cell_to_ordinal(system, leaf_level, child_id))
        return DGGSCursor(cursor.grid, child_level, child_id, ordinal, ordinal, nothing)
    end
    return DGGSCursor(cursor.grid, child_level, child_id, 0, -1, nothing)
end

# Index bounds of a child's stored ids: two binary searches inside the parent's
# window, against the child's two-sided `descendant_range`. Invalid ids inside
# the range can never be stored (that is the range contract), so the bounds are
# exact, not a superset.
function _stored_bounds(cursor::DGGSCursor{<:DGGSPartialGrid}, child_level::Int, child_id)
    cursor.last_index >= cursor.first_index || return (1, 0)
    system = _system(cursor)
    ids = cursor.grid.ids
    lo, hi = descendant_range(system, child_level, child_id, _leaf_level(cursor))
    first_index = searchsortedfirst(ids, lo, cursor.first_index, cursor.last_index, Base.Order.Forward)
    last_index = searchsortedlast(ids, hi, cursor.first_index, cursor.last_index, Base.Order.Forward)
    return (Int(first_index), Int(last_index))
end

function _range_child(cursor::DGGSCursor{<:DGGSPartialGrid,<:Any,Nothing}, child_id)
    child_level = cursor.level + 1
    first_index, last_index = _stored_bounds(cursor, child_level, child_id)
    return DGGSCursor(cursor.grid, child_level, child_id, first_index, last_index, nothing)
end

_nonempty(cursor::DGGSCursor{<:DGGSPartialGrid,<:Any,Nothing}) =
    cursor.last_index >= cursor.first_index

"""
    _child_cursors(cursor) -> Vector{<:DGGSCursor}

Children of a partial node on the `cell_parent` fallback path, non-empty ones
only. One pass over the parent's selection buckets every stored index under its
ancestor at the child level, so the whole child row costs
`O(selection * log radix)` rather than `O(selection * radix)`.
"""
function _child_cursors(cursor::DGGSCursor{<:DGGSPartialGrid,<:Any,Vector{Int}})
    system = _system(cursor)
    leaf_level = _leaf_level(cursor)
    child_level = cursor.level + 1
    ids = cursor.grid.ids
    child_ids = _child_ids(cursor)
    buckets = [Int[] for _ in eachindex(child_ids)]
    for index in cursor.selection
        ancestor = child_level == leaf_level ? ids[index] :
                   cell_parent(system, leaf_level, ids[index], child_level)
        slot = searchsortedfirst(child_ids, ancestor)
        (slot <= length(child_ids) && child_ids[slot] == ancestor) || continue
        push!(buckets[slot], index)
    end
    children = DGGSCursor{typeof(cursor.grid),typeof(cursor.id),Vector{Int}}[]
    for slot in eachindex(buckets)
        isempty(buckets[slot]) && continue
        push!(children, DGGSCursor(cursor.grid, child_level, child_ids[slot],
            1, length(buckets[slot]), buckets[slot]))
    end
    return children
end

# --------------------------------------------------------------------------
# SpatialTreeInterface
# --------------------------------------------------------------------------

STI.isspatialtree(::Type{<:DGGSCursor}) = true

# A dense node splits until it *is* a cell; a partial node also stops when it
# covers nothing, or once `bucket_size` says the remaining leaves are cheaper
# scanned than descended.
STI.isleaf(cursor::DGGSCursor{<:DGGSGrid}) = cursor.level >= _leaf_level(cursor)

function STI.isleaf(cursor::DGGSCursor{<:DGGSPartialGrid})
    count = _stored_count(cursor)
    count == 0 && return true
    cursor.level >= _leaf_level(cursor) && return true
    bucket_size = _bucket_size(cursor)
    return bucket_size > 0 && count <= bucket_size
end

STI.nchild(cursor::DGGSCursor{<:DGGSGrid}) = length(_child_ids(cursor))

function STI.nchild(cursor::DGGSCursor{<:DGGSPartialGrid,<:Any,Nothing})
    count = 0
    child_level = cursor.level + 1
    for child_id in _child_ids(cursor)
        first_index, last_index = _stored_bounds(cursor, child_level, child_id)
        count += last_index >= first_index
    end
    return count
end

STI.nchild(cursor::DGGSCursor{<:DGGSPartialGrid,<:Any,Vector{Int}}) =
    length(_child_cursors(cursor))

STI.getchild(cursor::DGGSCursor{<:DGGSGrid}) =
    (_dense_child(cursor, child_id) for child_id in _child_ids(cursor))

STI.getchild(cursor::DGGSCursor{<:DGGSPartialGrid,<:Any,Nothing}) =
    Iterators.filter(_nonempty,
        (_range_child(cursor, child_id) for child_id in _child_ids(cursor)))

STI.getchild(cursor::DGGSCursor{<:DGGSPartialGrid,<:Any,Vector{Int}}) =
    _child_cursors(cursor)

# The indexed form is the interface's convenience accessor, not the traversal
# path (both depth-first searches iterate `getchild(node)`), so walking the
# iterator is fine and keeps one definition of "which children exist".
function STI.getchild(cursor::DGGSCursor, i::Int)
    i >= 1 || throw(BoundsError(cursor, i))
    seen = 0
    for child in STI.getchild(cursor)
        seen += 1
        seen == i && return child
    end
    throw(BoundsError(cursor, i))
end

# Every extent below is derived from `cell_boundary` — an inverse projection per
# cell — rather than read off the node, so the dual depth-first search should
# cache a node's child extents instead of re-deriving them per opposing child.
STI.node_extent_is_expensive(::Type{<:DGGSCursor}) = true

function STI.node_extent(cursor::DGGSCursor{<:DGGSGrid})
    cursor.level < 0 && return full_sphere_extent()
    system = _system(cursor)
    leaf_level = _leaf_level(cursor)
    cursor.level == leaf_level && return cell_cap(system, leaf_level, cursor.id)
    return subtree_cap(system, cursor.level, cursor.id, leaf_level)
end

# How many stored leaves an internal node will spend `cell_boundary` calls on
# to tighten its extent past its own cell's cap (`STI.node_extent` below). A
# constant, deliberately, not a fraction of the subtree: it is the per-node
# work bound, and the whole defect it replaces was a limit that scaled with the
# subtree instead. 64 clears one child row at every aperture this package wires
# (4, 7, 9) and two rows of an aperture-7 grid (49) — past that a node is not
# "a handful of stored leaves" any more and its own cell is the honest bound.
const STORED_UNION_CAP_LIMIT = 64

function STI.node_extent(cursor::DGGSCursor{<:DGGSPartialGrid})
    cursor.level < 0 && return full_sphere_extent()
    system = _system(cursor)
    leaf_level = _leaf_level(cursor)
    cursor.level == leaf_level && return cell_cap(system, leaf_level, cursor.id)
    # Where a parent geographically contains its descendants, its own O(1) cap
    # already bounds any stored subset as tightly as a union cap would, so
    # paying `cells_cap`'s O(stored) boundary calls per node buys nothing at
    # all. HEALPix is the case; that is what the trait means.
    has_exact_subtree_cap(system) &&
        return subtree_cap(system, cursor.level, cursor.id, leaf_level)
    # Otherwise the node *is* the cell `(level, id)`, and the hierarchy already
    # bounds that cell's whole subtree in O(1): `cell_cap_inflation` exists
    # precisely so a cell's cap covers the descendants that overhang it (1.048
    # measured for IGEO7 and 1.052 for H3 against the shared 1.2, 1.469 for A5
    # against the 1.75 it wires — `test/<System>/test_*_kernel.jl` measures all
    # of them). So the O(stored) union cap has to earn its keep, and it can
    # only do that where the node stores a *proper* subset of its subtree:
    #
    #   * a node that stores the whole subtree gains at most the inflation
    #     factor in radius from the union and pays one `cell_boundary` per leaf
    #     for it. That was the old rule's pure loss, and it fell on the nodes
    #     nearest the leaves — the many. Aperture 7 puts 7^3 = 343 and
    #     7^4 = 2401 either side of the old `SUBTREE_CAP_EXACT_LIMIT`, so every
    #     leaf's boundary was recomputed by each of its three nearest ancestors
    #     even under a traversal that visits each node once: 5.8M
    #     `cell_boundary` calls for 117,649 cells in the profile that prompted
    #     this, 49x redundancy, 80% of a `Regridder`'s construction time.
    #   * a node storing a handful of a big cell's leaves — the sparse chunk a
    #     partial grid exists for — really is bounded far more tightly by a cap
    #     around those leaves, and that tightness is pruning power. It keeps
    #     the union, for a bounded `STORED_UNION_CAP_LIMIT` boundary calls
    #     rather than for as many as the subtree happens to have.
    #
    # `_stored_count` is O(1) on every path, and `subtree_leaf_count` is asked
    # only once the count is already known small — it is O(1) in every wiring
    # (the generic trees lean on that elsewhere too, see its docstring). An
    # empty node takes the cell cap as well: `cells_cap` answers the full
    # sphere for an empty batch, which is sound and is also the worst extent
    # there is, and the documented empty node forms do reach here.
    stored = _stored_count(cursor)
    if 0 < stored <= STORED_UNION_CAP_LIMIT &&
       stored < subtree_leaf_count(system, cursor.level, cursor.id, leaf_level)
        return cells_cap(system, leaf_level, _stored_ids(cursor))
    end
    return cell_cap(system, cursor.level, cursor.id)
end

"""
    STI.child_indices_extents(cursor) -> Vector{Tuple{Int,SphericalCap{Float64}}}

Leaf indices and leaf caps of a leaf node. Materialized deliberately: the dual
depth-first search re-iterates this once per opposing leaf, so a lazy generator
would recompute every `cell_cap` (and its boundary) on each pass, and a `Tuple`
is type-unstable past a handful of entries.
"""
function STI.child_indices_extents(cursor::DGGSCursor)
    STI.isleaf(cursor) ||
        throw(ArgumentError("child_indices_extents is only valid for leaf nodes"))
    system = _system(cursor)
    leaf_level = _leaf_level(cursor)
    count = _stored_count(cursor)
    entries = Vector{Tuple{Int,GO.UnitSpherical.SphericalCap{Float64}}}(undef, count)
    for i in 1:count
        index, id = _leaf_entry(cursor, i)
        entries[i] = (index, cell_cap(system, leaf_level, id))
    end
    return entries
end

# --------------------------------------------------------------------------
# ConservativeRegridding.Trees
#
# `best_manifold` is what makes the one-argument `treeify(grid)` work: it is
# `Trees`' own generic `treeify(grid) = treeify(best_manifold(grid), grid)`
# resolving through these methods, not a second function. DGGS geometry is on
# the unit sphere by construction, so the manifold argument carries no
# information at a call site and naming it is pure ceremony — the package
# re-exports `treeify` (and `ncells`/`getcell`) so the short form is the one
# users reach for. The manifold-explicit methods stay available and are what
# `ConservativeRegridding.Regridder` itself calls.
# --------------------------------------------------------------------------

# Always a `Spherical`, never the declared manifold: the tessellation is on the
# authalic sphere, and `ConservativeRegridding` imports only `Planar`/
# `Spherical`, so handing it a `Geodesic` would be a `MethodError`. A grid that
# declares no ellipsoid resolves to itself, so this is a no-op on the default
# path. See `grid_manifold` for the declared frame.
GOCore.best_manifold(grid::DGGSGrid) = authalic_sphere(grid.manifold)
GOCore.best_manifold(grid::DGGSPartialGrid) = authalic_sphere(grid.manifold)
GOCore.best_manifold(cursor::DGGSCursor) = GOCore.best_manifold(cursor.grid)

function Trees.treeify(::GO.Spherical, grid::DGGSGrid)
    return DGGSCursor(grid, -1, zero(cell_id_type(grid.system)),
        1, Int(num_cells(grid.system, grid.level)), nothing)
end

function Trees.treeify(::GO.Spherical, grid::DGGSPartialGrid)
    # `grid.ids` is stored by reference, never collected or reordered: leaf
    # index i is position i of the lookup vector the grid was built from.
    count = length(grid.ids)
    has_descendant_ranges(grid.system) &&
        return DGGSCursor(grid, grid.root_level, grid.root_id, 1, count, nothing)
    return DGGSCursor(grid, grid.root_level, grid.root_id, 1, count,
        Int[i for i in eachindex(grid.ids)])
end

Trees.treeify(::GO.Spherical, cursor::DGGSCursor) = cursor

# The root cursor's window is the whole grid, so one definition serves both the
# tree-level (`1:ncells(tree)`) and node-level index spaces.
Trees.ncells(cursor::DGGSCursor) = _stored_count(cursor)

function Trees.getcell(cursor::DGGSCursor, i::Int)
    count = _stored_count(cursor)
    1 <= i <= count || throw(BoundsError(1:count, i))
    _, id = _leaf_entry(cursor, i)
    return cell_polygon_unitsphere(_system(cursor), _leaf_level(cursor), id)
end

Trees.getcell(cursor::DGGSCursor) = (Trees.getcell(cursor, i) for i in 1:Trees.ncells(cursor))

# --------------------------------------------------------------------------
# Parallelization policy
#
# Mirrors `Trees.AbstractQuadtreeCursor`'s: spawn once a subtree's leaf count
# drops below `total / (nthreads * 32)`, which lands a few hundred tasks for
# the scheduler to balance. Without this method `ConservativeRegridding`'s
# `::Any` fallback (spawn when the cap covers less than a quarter sphere)
# silently applies — a cap-area heuristic, not a work estimate.
# --------------------------------------------------------------------------

const PARALLELIZE_CHUNKS_PER_THREAD = 32

function Trees.should_parallelize(cursor::DGGSCursor, ::GO.UnitSpherical.SphericalCap)
    threshold = max(1, _total_cells(cursor.grid) ÷
                       (Threads.nthreads() * PARALLELIZE_CHUNKS_PER_THREAD))
    return _stored_count(cursor) <= threshold
end
