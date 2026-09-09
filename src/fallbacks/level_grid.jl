# The one package-provided complete-level grid, and the `levelgrid` default that
# builds it. A system that needs no state beyond `(system, level)` — which is
# every system shipped here — implements the base grid interface as five
# system-level methods and lets this type do the forwarding.

"""
    HierarchicalLevelGrid(sys, level) <: AbstractGrid

The complete grid of `sys` at `level`: every cell the system has there, in the
system's canonical dense order. This is what [`levelgrid`](@ref) returns unless
a system overrides it, and it is to a complete level what
[`PartialGrid`](@ref DiscreteGlobalGrids.Engine.PartialGrid) is to a subset
of one.

It stores the system and the level and nothing else, so constructing one is
O(1). The base grid interface is answered by forwarding to system-level
counterparts, which are the implementor surface:

| grid method | system method |
|---|---|
| `ncells(grid)` | `ncells(sys, level)` |
| `cellindex(grid, i)` | `cellindex(sys, level, i)` |
| `localindex(grid, c)` | `globalindex(sys, c)` |
| `cell_boundary(grid, c)` | `cell_boundary(sys, c)` |
| `cell_centroid(grid, c)` | `cell_centroid(sys, c)` |

The geometry pair takes no level: an [`AbstractCellIndex`](@ref) carries its own,
so the level would be redundant and could disagree with the id. The grid methods
supply what the level *is* needed for — the bounds check on an index, and the
rejection of an id from another level, which a complete-level grid must not
answer about.

Fast paths stay available: a system attaches them to
`HierarchicalLevelGrid{TheSystem}`, so `cellat`, `neighbors`, `ring`,
`cell_area` and the rest dispatch on the type parameter.

    julia> grid = levelgrid(HEALPixSystem(), 2);

    julia> grid isa HierarchicalLevelGrid{HEALPixSystem}
    true
"""
struct HierarchicalLevelGrid{S<:AbstractHierarchicalGridSystem} <: AbstractGrid
    system::S
    level::Int
end

HierarchicalLevelGrid(sys::AbstractHierarchicalGridSystem, l::Integer) =
    HierarchicalLevelGrid{typeof(sys)}(sys, Int(l))

"""
    levelgrid(sys, l) -> HierarchicalLevelGrid

The complete level grid of `sys` at `l`, after checking `l` against
[`levels`](@ref). A system overrides this only to return a grid type of its own.
"""
function levelgrid(sys::AbstractHierarchicalGridSystem, l::Integer)
    lvl = Int(l)
    lvl in levels(sys) || throw(ArgumentError(
        "level $lvl is outside levels($(nameof(typeof(sys)))()) = $(levels(sys))"))
    return HierarchicalLevelGrid(sys, lvl)
end

# --- the base grid interface ----------------------------------------------

system(grid::HierarchicalLevelGrid) = grid.system
level(grid::HierarchicalLevelGrid) = grid.level

ncells(grid::HierarchicalLevelGrid) = ncells(grid.system, grid.level)

function cellindex(grid::HierarchicalLevelGrid, i::Int)
    1 <= i <= ncells(grid) || throw(BoundsError(grid, i))
    return cellindex(grid.system, grid.level, i)
end

# `_canonical` is the whole "is this cell of this grid" question: wrong level,
# a scheme the system does not claim, and a scheme it claims but cannot convert
# all answer `nothing`, leaving the system method to see canonical ids only.
function localindex(grid::HierarchicalLevelGrid, c::AbstractCellIndex)
    target = _canonical(grid, c)
    target === nothing && return nothing
    return globalindex(grid.system, target)
end

cell_boundary(grid::HierarchicalLevelGrid, c::AbstractCellIndex) =
    cell_boundary(grid.system, _at_level(grid, c))

cell_centroid(grid::HierarchicalLevelGrid, c::AbstractCellIndex) =
    cell_centroid(grid.system, _at_level(grid, c))

# Geometry is a function of the id alone, so forwarding an id from another level
# would return the geometry of a cell this grid does not contain. That is an
# error, not an answer.
@inline function _at_level(grid::HierarchicalLevelGrid, c::AbstractCellIndex)
    level(c) == grid.level || throw(ArgumentError(
        "cell $c is at level $(level(c)), not the grid's level $(grid.level)"))
    return c
end

function Base.show(io::IO, grid::HierarchicalLevelGrid)
    print(io, "HierarchicalLevelGrid(", typeof(grid.system).name.name,
        ", level=", grid.level, ", ncells=", ncells(grid), ")")
end

Base.show(io::IO, ::MIME"text/plain", grid::HierarchicalLevelGrid) = show(io, grid)
