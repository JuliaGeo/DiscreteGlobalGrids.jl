# Identity and hierarchy fallbacks. A bare `Int` is a grid position; a typed
# `AbstractCellIndex` is a cell identifier.

"""
    cellindextypes(grid) -> Tuple

The index schemes a grid can name its cells in. Asks the grid's system; a
standalone grid reports the one scheme it demonstrably uses.
"""
function cellindextypes(grid::AbstractGrid)
    sys = system(grid)
    sys === nothing || return cellindextypes(sys)
    ncells(grid) == 0 && return ()
    return (typeof(cellindex(grid, 1)),)
end

"""
    cellindex(grid, i, T) -> T

The id of the cell at position `i` in the requested scheme. Generic:
`reindex(T, system(grid), cellindex(grid, i))`.
"""
function cellindex(grid::AbstractGrid, i::Int, ::Type{T}) where {T<:AbstractCellIndex}
    sys = system(grid)
    c = cellindex(grid, i)
    sys === nothing || return reindex(T, sys, c)
    c isa T && return c
    throw(ArgumentError("$(typeof(grid)) has no system, so it can only name cells " *
                        "as $(typeof(c)); requested $T"))
end

"""
    reindex(T, sys, c) -> T

Convert a cell id between the schemes of one system. The generic method knows
only the identity conversion; a system with more than one scheme in
[`cellindextypes`](@ref) wires the rest.
"""
function reindex(::Type{T}, sys::AbstractHierarchicalGridSystem,
        c::AbstractCellIndex) where {T<:AbstractCellIndex}
    # `cellindextypes` first, deliberately: the identity case must still be an
    # error for a scheme the system does not claim, and a system that has not
    # declared its canonical type has not answered the question at all.
    supported = cellindextypes(sys)
    T in supported || throw(ArgumentError(
        "$(typeof(sys)) cannot name cells as $T; it supports $(supported)"))
    c isa T && return c
    throw(ArgumentError(
        "no conversion from $(typeof(c)) to $T is wired for $(typeof(sys))"))
end

"""
    cellposition(grid, c) -> Union{Int,Nothing}

Return the position of `c`, or `nothing` if absent. The generic fallback scans
`1:ncells(grid)` linearly; searchable grids should override it.
"""
function cellposition(grid::AbstractGrid, c::AbstractCellIndex)
    target = _canonical(grid, c)
    target === nothing && return nothing
    for i in 1:ncells(grid)
        cellindex(grid, i) == target && return i
    end
    return nothing
end

# `BoundsError` on a grid names the valid position range instead of the grid's
# type parameters.
function Base.summary(io::IO, grid::AbstractGrid)
    sys = system(grid)
    sys === nothing || print(io, nameof(typeof(sys)), " ")
    l = level(grid)
    l === nothing || print(io, "level-", l, " ")
    print(io, nameof(typeof(grid)), " over positions 1:", ncells(grid))
    return nothing
end

# The registered systems that name their cells with `typeof(c)`. Empty for a
# cell type no shipped system claims. Only the error path below calls it.
_owning_systems(c::AbstractCellIndex) =
    filter(s -> typeof(c) in cellindextypes(s), collect(DGG.systems()))

"""
    _foreign_cell(f, sys, c)

Report `c` as belonging to another system, or re-raise the `MethodError` when it
does not: a scheme `sys` declares in [`cellindextypes`](@ref) but has no method
for is unimplemented, not foreign.
"""
@noinline function _foreign_cell(f, sys::AbstractHierarchicalGridSystem,
        c::AbstractCellIndex)
    typeof(c) in cellindextypes(sys) && throw(MethodError(f, (sys, c)))
    owners = _owning_systems(c)
    owner = isempty(owners) ? "a system not registered here" :
            join((string(nameof(typeof(s))) for s in owners), " or ")
    throw(ArgumentError(
        "$c is a cell of $owner, not of $(nameof(typeof(sys))), which names " *
        "cells as $(nameof(cellindextype(sys))). Cell ids do not convert " *
        "between systems: locate the cell by coordinate with `cellat`, or " *
        "move data across systems with `regrid`."))
end

# The system verbs a cell id reaches directly. Each shipped system defines a
# more specific method, so these are hit only by a cell the system cannot name.
for f in (:cell_boundary, :cell_centroid, :cellposition, :children)
    @eval $f(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex) =
        _foreign_cell($f, sys, c)
end

Base.parent(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex) =
    _foreign_cell(Base.parent, sys, c)

# Convert to the grid's canonical id, returning `nothing` for a wrong level,
# unsupported input scheme, or rejected value. Only `ArgumentError` from
# `reindex` denotes an absent value; other exceptions remain visible.
function _canonical(grid::AbstractGrid, c::AbstractCellIndex)
    l = level(grid)
    l === nothing || level(c) == l || return nothing
    sys = system(grid)
    sys === nothing && return c
    T = cellindextype(sys)
    c isa T && return c
    typeof(c) in cellindextypes(sys) || return nothing
    return try
        reindex(T, sys, c)
    catch err
        err isa ArgumentError || rethrow()
        nothing
    end
end

# ===========================================================================
# Derived hierarchy walks
# ===========================================================================

"""
    ancestor(sys, c, l) -> AbstractCellIndex

The ancestor of `c` at level `l`. Generic: repeated [`parent`](@ref). A system
that can drop `level(c) - l` digits in one operation should override.
"""
function ancestor(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex, l::Integer)
    target = Int(l)
    lc = level(c)
    target <= lc || throw(ArgumentError(
        "ancestor level $target is deeper than the cell's own level $lc"))
    target >= first(levels(sys)) || throw(ArgumentError(
        "ancestor level $target is above the root level $(first(levels(sys)))"))
    current = c
    while level(current) > target
        current = Base.parent(sys, current)
    end
    return current
end

"""
    descendants(sys, c, l)

Return all level-`l` descendants of `c`, ascending and materialized in
`O(subtree)`. Sorted-subtree systems resolve one [`descendant_range`](@ref);
others expand [`children`](@ref) level by level and sort.
"""
function descendants(sys::AbstractHierarchicalGridSystem, c::AbstractCellIndex, l::Integer)
    target = Int(l)
    lc = level(c)
    target >= lc || throw(ArgumentError(
        "descendant level $target is above the cell's own level $lc"))
    target <= max_level(sys) || throw(ArgumentError(
        "descendant level $target is past max_level $(max_level(sys))"))
    T = cellindextype(sys)
    target == lc && return T[c]
    if has_sorted_subtrees(sys)
        grid = levelgrid(sys, target)
        return T[cellindex(grid, i) for i in descendant_range(sys, c, target)]
    end
    current = T[c]
    for _ in lc:(target-1)
        next = T[]
        for cell in current
            append!(next, children(sys, cell))
        end
        current = next
    end
    # Sibling order is ascending by contract, but without sorted subtrees the
    # concatenation of two siblings' subtrees need not be.
    return sort!(current)
end
