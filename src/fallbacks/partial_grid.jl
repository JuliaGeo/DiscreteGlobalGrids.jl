# ---------------------------------------------------------------------------
# `PartialGrid` — the one subset-grid type
#
# A subset of one level of one system, addressed by position in its sorted id
# vector. This is what a regional workflow, a chunk, or a `DimensionalData`
# lookup hands to the tree and the regridder: leaf index `i` of the resulting
# tree is position `i` of `grid.ids`, so a `Regridder` lines up with the data
# array without a permutation.
#
# It also absorbs the old lazy-subtree type. A partial grid built over one
# cell's subtree keeps `root_level`/`root_id`, so `treeify` starts descent at
# that node instead of at the whole sphere, and its ids need not be stored at
# all — `SubtreeIds` computes them from the level grid.
# ---------------------------------------------------------------------------

"""
    SubtreeIds(grid, first, n) <: AbstractVector

The ids of `n` consecutive positions of a complete level grid, computed on
demand rather than stored — the descendants of one cell at one level, when the
system has [`has_sorted_subtrees`](@ref) and the interval therefore exists.

A handful of words however many cells it names, which is what makes a
subtree-shaped [`PartialGrid`](@ref) free to build.
"""
struct SubtreeIds{G<:AbstractGrid,ID} <: AbstractVector{ID}
    grid::G
    first::Int
    n::Int
end

function SubtreeIds(grid::AbstractGrid, first::Integer, n::Integer)
    ID = typeof(cellindex(grid, Int(first)))
    return SubtreeIds{typeof(grid),ID}(grid, Int(first), Int(n))
end

Base.size(v::SubtreeIds) = (v.n,)
Base.IndexStyle(::Type{<:SubtreeIds}) = Base.IndexLinear()
Base.getindex(v::SubtreeIds, i::Int) = cellindex(v.grid, v.first + i - 1)

# Positions of a complete level grid ascend in canonical id order, so the O(n)
# verification `PartialGrid` runs on an arbitrary vector has nothing to find.
Helpers.strictly_increasing(::SubtreeIds) = true

"""
    PartialGrid(sys, level, ids; bucket_size = 0, root = nothing)
    PartialGrid(sys, c::AbstractCellIndex, level; bucket_size = 0)

A subset of `levelgrid(sys, level)`: the cells named by `ids`, which must be
**strictly ascending** canonical ids at `level`. Positions run `1:length(ids)`,
in that order — the vector is stored by reference, never copied or reordered,
so position `i` of a caller's id vector stays position `i` of the grid.

The second form is the subtree of one cell: every level-`level` descendant of
`c`, with `c` remembered as the tree root so descent starts there. Where the
system has [`has_sorted_subtrees`](@ref) the ids are computed rather than
materialised ([`SubtreeIds`](@ref)), so building one is O(1) whatever the
subtree's size.

# Keywords

  - `bucket_size` stops tree descent once a node covers that few stored cells
    and scans them instead. `0` (the default) descends to single cells.
  - `root` is a cell all the ids are descendants of. Passing it is what keeps
    cursor descent *windowed* over a chunk instead of restarting from the
    system's root cells.

# What is checked

`level` against [`levels`](@ref) (through [`levelgrid`](@ref)), `eltype(ids)`
against [`cellindextype`](@ref), strict ascent, that the **endpoints** live at
`level`, and — for a rooted grid — that the endpoints really are descendants of
`root` (complete and O(1) where `has_sorted_subtrees` holds, since the sorted
endpoints then bound every id; an [`ancestor`](@ref) check on the endpoints
otherwise).

Interior ids are trusted. Validating each one costs a hierarchy walk per cell,
which is the expense a chunk built from a validated lookup has already paid,
and an invalid id surfaces at the first geometry call that decodes it.
"""
struct PartialGrid{S<:AbstractHierarchicalGridSystem,V<:AbstractVector,ID,G<:AbstractGrid} <: AbstractGrid
    system::S
    level::Int
    ids::V
    complete::G          # `levelgrid(sys, level)`, for geometry and id bounds
    bucket_size::Int
    root_level::Int      # `< first(levels(sys))` means "not rooted"
    root_id::ID          # meaningless, and never read, when not rooted

    function PartialGrid(sys::AbstractHierarchicalGridSystem, lvl::Integer,
            ids::AbstractVector; bucket_size::Integer=0, root=nothing)
        l = Int(lvl)
        # `levelgrid` is the level validator, and its result is the grid every
        # geometry call below delegates to, so it is asked once and kept.
        complete = levelgrid(sys, l)
        ID = cellindextype(sys)
        eltype(ids) === ID || throw(ArgumentError(
            "ids must be a $ID vector for $(typeof(sys)), got eltype $(eltype(ids))"))
        Int(bucket_size) >= 0 || throw(ArgumentError("bucket_size must be non-negative"))
        Helpers.strictly_increasing(ids) ||
            throw(ArgumentError("ids must be strictly ascending"))
        if !isempty(ids)
            (level(first(ids)) == l && level(last(ids)) == l) || throw(ArgumentError(
                "ids must all be level-$l cells, got levels " *
                "$(level(first(ids))) and $(level(last(ids))) at the endpoints"))
        end
        root_level = first(levels(sys)) - 1
        root_id = _placeholder_root(sys)
        if root !== nothing
            # The root's type becomes the cursor's cell-id type parameter, and
            # every descent step is stored into that field — so a non-canonical
            # root fails on the first `children` call rather than here, where
            # the caller can still see what they passed.
            root isa ID || throw(ArgumentError(
                "root must be a $ID for $(typeof(sys)), got $(typeof(root))"))
            root_level = level(root)
            root_level <= l || throw(ArgumentError(
                "root cell is at level $root_level, deeper than the grid's level $l"))
            root_id = root
            isempty(ids) || _check_rooted(sys, complete, ids, root, l)
        end
        return new{typeof(sys),typeof(ids),typeof(root_id),typeof(complete)}(
            sys, l, ids, complete, Int(bucket_size), root_level, root_id)
    end
end

function PartialGrid(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex,
        lvl::Integer; bucket_size::Integer=0)
    l = Int(lvl)
    l >= level(c) || throw(ArgumentError(
        "subtree level $l is above the root cell's own level $(level(c))"))
    if has_sorted_subtrees(sys)
        range = descendant_range(sys, c, l)
        ids = SubtreeIds(levelgrid(sys, l), first(range), length(range))
        return PartialGrid(sys, l, ids; bucket_size, root=c)
    end
    return PartialGrid(sys, l, descendants(sys, c, l); bucket_size, root=c)
end

_placeholder_root(sys::AbstractHierarchicalGridSystem) = first(rootcells(sys))

function _check_rooted(sys, complete, ids, root, l)
    if has_sorted_subtrees(sys)
        # The sorted endpoints bound every id, so two comparisons decide the
        # whole vector — and they decide it completely, not as a guard.
        range = descendant_range(sys, root, l)
        lo = cellindex(complete, first(range))
        hi = cellindex(complete, last(range))
        (!isless(first(ids), lo) && !isless(hi, last(ids))) || throw(ArgumentError(
            "ids must all be descendants of the root cell $root"))
    else
        rl = level(root)
        (ancestor(sys, first(ids), rl) == root && ancestor(sys, last(ids), rl) == root) ||
            throw(ArgumentError("ids must all be descendants of the root cell $root"))
    end
    return nothing
end

# --- the base grid interface ----------------------------------------------

ncells(grid::PartialGrid) = length(grid.ids)
cellindex(grid::PartialGrid, i::Int) = grid.ids[i]
system(grid::PartialGrid) = grid.system
level(grid::PartialGrid) = grid.level
cell_boundary(grid::PartialGrid, c::AbstractCellIndex) = cell_boundary(grid.complete, c)
cell_centroid(grid::PartialGrid, c::AbstractCellIndex) = cell_centroid(grid.complete, c)

# Not derived, forwarded. `cell_area`'s generic is the ring's polygon area,
# which is the right answer only when the ring IS the cell: HEALPix and ISEA4R
# have curvilinear edges and override it on their level grid with the exact
# `4pi/ncells`. A `PartialGrid` is a different type, so the generic used to win
# here and a subset of a level grid reported areas its parent grid did not
# agree with. The subset changes which cells there are and nothing about their
# geometry, so the complete grid is the authority for every one of them.
cell_area(grid::PartialGrid, c::AbstractCellIndex) = cell_area(grid.complete, c)

# The ids are sorted, so the O(n) generic scan is two comparisons here.
function cellposition(grid::PartialGrid, c::AbstractCellIndex)
    target = _canonical(grid, c)
    target === nothing && return nothing
    i = searchsortedfirst(grid.ids, target)
    (i <= length(grid.ids) && grid.ids[i] == target) || return nothing
    return i
end

_is_rooted(grid::PartialGrid) = grid.root_level >= first(levels(grid.system))

function Base.show(io::IO, grid::PartialGrid)
    print(io, "PartialGrid(", typeof(grid.system).name.name, ", level=", grid.level,
        ", ncells=", ncells(grid))
    _is_rooted(grid) && print(io, ", root=", grid.root_id)
    grid.bucket_size > 0 && print(io, ", bucket_size=", grid.bucket_size)
    print(io, ")")
end

Base.show(io::IO, ::MIME"text/plain", grid::PartialGrid) = show(io, grid)
