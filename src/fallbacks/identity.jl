# ---------------------------------------------------------------------------
# Identity: index schemes, positions, and the derived hierarchy walks
#
# Position vs identity is the one convention everything here turns on: a bare
# `Int` is a position in `1:ncells(grid)`, a typed id is a name. These are the
# generics that convert between the two, and between one naming scheme and
# another.
# ---------------------------------------------------------------------------

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

The position of `c` in the grid's dense order, or `nothing` if it is not there.

The generic fallback is a **linear scan** over `1:ncells(grid)`: with only
[`cellindex`](@ref) in hand there is nothing else to do, since a grid's
canonical order is its own choice and need not be searchable. Every grid that
can answer faster should override this — `PartialGrid` binary-searches its
sorted ids, and a system's complete level grid answers in closed form.
"""
function cellposition(grid::AbstractGrid, c::AbstractCellIndex)
    target = _canonical(grid, c)
    target === nothing && return nothing
    for i in 1:ncells(grid)
        cellindex(grid, i) == target && return i
    end
    return nothing
end

# The id as the grid's own scheme, or `nothing` when it cannot be one of the
# grid's cells at all (wrong scheme, or a different level).
function _canonical(grid::AbstractGrid, c::AbstractCellIndex)
    l = level(grid)
    l === nothing || level(c) == l || return nothing
    sys = system(grid)
    sys === nothing && return c
    T = cellindextype(sys)
    c isa T && return c
    T in cellindextypes(sys) || return nothing
    return reindex(T, sys, c)
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

Every descendant of `c` at level `l`, ascending. O(subtree) and materialising —
see the interface docstring, and reach for [`descendant_range`](@ref) instead
wherever the system has it.

Two implementations, picked by [`has_sorted_subtrees`](@ref): with the trait
the answer is one `descendant_range` read off as positions in `levelgrid`;
without it, [`children`](@ref) expanded level by level and sorted.
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
