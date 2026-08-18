# A `PartialGrid` is a sorted subset of one system level. Grid positions match
# positions in its id vector. Rooted subsets start tree descent at their root.

"""
    SubtreeIds(grid, first, n) <: AbstractVector

Lazy ids for `n` consecutive positions of a complete level grid. This provides
`O(1)` construction for rooted subsets in sorted-subtree systems.
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

# The check is written out rather than left to Base: a `getindex(::T, ::Int)`
# defined directly is the whole method, so without it a position past the end
# resolves an id from the NEXT subtree instead of throwing — and `PartialGrid`'s
# `cellindex`, which the position contract says bounds-checks, is this call.
Base.@propagate_inbounds function Base.getindex(v::SubtreeIds, i::Int)
    @boundscheck checkbounds(v, i)
    return cellindex(v.grid, v.first + i - 1)
end

# Positions of a complete level grid ascend in canonical id order, so the O(n)
# verification `PartialGrid` runs on an arbitrary vector has nothing to find.
Helpers.strictly_increasing(::SubtreeIds) = true

"""
    PartialGrid(sys, level, ids; bucket_size = 0, root = nothing)
    PartialGrid(sys, c::AbstractCellIndex, level; bucket_size = 0)

A subset of `levelgrid(sys, level)`. `ids` must be strictly ascending canonical
ids at `level` and is stored by reference without copying or reordering.

The second form contains all level-`level` descendants of `c` and stores `c` as
the tree root. Sorted-subtree systems use [`SubtreeIds`](@ref), making
construction `O(1)`.

# Keywords

  - `bucket_size` stops descent at that many stored cells; `0` reaches cells.
  - `root` declares a common ancestor and starts descent there.

# What is checked

Validation covers the level, id type, strict ascent, endpoint levels, and rooted
endpoint ancestry. For sorted subtrees, endpoint checks cover the entire vector.

Interior ids are not individually validated; a bad one surfaces at the first
geometry call that decodes it.
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
            # Validate before the root becomes the cursor's concrete id type.
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

# `AuthalicGrid`'s wrappability check, stated here because the subset type is.
Fallbacks._check_wrappable(::PartialGrid) = throw(ArgumentError(
    "wrap the SYSTEM, not the subset: `PartialGrid(AuthalicSystem(sys), level, ids)`. \
A subset is a property of the id set and the warp is a property of the system, and \
only that order keeps the tree cursor's position windows correct."))

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

# Forwarded, not derived. The generic `cell_area` is the ring's polygon area,
# right only where the ring IS the cell: HEALPix and ISEA4R have curvilinear
# edges and override it on their level grid with the exact `4pi/ncells`. A
# subset changes which cells exist, never their geometry, so the complete grid
# is the authority.
cell_area(grid::PartialGrid, c::AbstractCellIndex) = cell_area(grid.complete, c)

# Membership as a predicate, which is what the subset law in `neighbors`'
# contract is written with: `filter(in(sub), ring(complete, c, k))` has to RUN,
# and Base's fallback would need a grid to be iterable. `CellVector` carries the
# same method for the same reason.
Base.in(c::AbstractCellIndex, grid::PartialGrid) = cellposition(grid, c) !== nothing

# The ids are sorted, so the O(n) generic scan is two comparisons here.
function cellposition(grid::PartialGrid, c::AbstractCellIndex)
    target = _canonical(grid, c)
    target === nothing && return nothing
    i = searchsortedfirst(grid.ids, target)
    (i <= length(grid.ids) && grid.ids[i] == target) || return nothing
    return i
end

# [`subset_span`](@ref) over a sorted id vector. Ids ascend with positions on
# every system here — the same fact `_check_rooted` decides a whole vector's
# ancestry by — so the block `lo:hi` maps to the id interval `[idlo, idhi]` and
# one `searchsortedfirst` answers all three verdicts. The `ALL` case needs no
# scan: the ids are strictly ascending, so `hi - lo + 1` of them between the two
# endpoints inclusive is exactly the block with nothing missing.
function subset_span(grid::PartialGrid, lo::Int, hi::Int)
    ids = grid.ids
    isempty(ids) && return _SPAN_NONE
    idlo = cellindex(grid.complete, lo)
    i = searchsortedfirst(ids, idlo)
    i <= length(ids) || return _SPAN_NONE
    idhi = cellindex(grid.complete, hi)
    isless(idhi, @inbounds ids[i]) && return _SPAN_NONE
    if @inbounds(ids[i]) == idlo
        j = i + (hi - lo)
        (j <= length(ids) && @inbounds(ids[j]) == idhi) && return _SPAN_ALL
    end
    return _SPAN_SOME
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
