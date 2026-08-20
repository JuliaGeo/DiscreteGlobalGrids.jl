# # Naming a set of cells
#
# Everything the mesh builder needs is a list of cell ids plus something that
# can hand back a boundary ring for one of them.  `DiscreteGlobalGrids` spells
# the second one `cell_boundary(source, cell)` and implements it for systems,
# grids and multi-order sets alike, so the builder can stay indifferent to which
# of those the user passed.

"""
    CellSet(source, cells)

A set of DGGS cells to draw: the ids in `cells`, whose geometry is read as
`DiscreteGlobalGrids.cell_boundary(source, cell)`.

`source` may be a system, a grid or a `MultiOrderCellSet` — anything that
implements `cell_boundary`.  Cells need not share a level when `source` does not
require it, which is what lets a `MultiOrderCellSet` be drawn as it is stored.

Build one with [`cellset`](@ref) rather than by hand.
"""
struct CellSet{S, C <: AbstractVector}
    source::S
    cells::C
end

Base.length(cs::CellSet) = length(cs.cells)

"""
    GridCells(grid)

The cells of `grid` as a lazy vector, `cellindex(grid, i)` at position `i`.

Every grid can be enumerated this way, including ones outside a hierarchy that
have no `CellVector`.
"""
struct GridCells{G, ID} <: AbstractVector{ID}
    grid::G
end

function GridCells(grid)
    ID = DGG.cellindextype(DGG.system(grid))
    return GridCells{typeof(grid), ID}(grid)
end

Base.size(gc::GridCells) = (DGG.ncells(gc.grid),)
Base.@propagate_inbounds Base.getindex(gc::GridCells, i::Int) = DGG.cellindex(gc.grid, i)

"""
    cellset(x) -> CellSet
    cellset(system, ids) -> CellSet

Read `x` as a set of cells to draw.

Accepts a `CellSet`, an `AbstractGrid`, a `CellVector`, a `CellLookup`, a
`MultiOrderCellSet` or `MultiOrderCoverage`, one of the subtree iterators, or a
system together with a vector of ids.  Supporting a new container is a method
here.
"""
function cellset end

cellset(cs::CellSet) = cs
cellset(grid::DGG.AbstractGrid) = CellSet(grid, GridCells(grid))
cellset(set::DGG.MultiOrderCellSet) = CellSet(set, set.cells)
cellset(coverage::DGG.MultiOrderCoverage) = cellset(parent(coverage))

# A `CellVector` knows its system and its level but is not itself a geometry
# source, so pair it with the level grid it indexes into.
cellset(cv::DGG.CellVector) =
    CellSet(DGG.levelgrid(DGG.system(cv), DGG.level(cv)), cv)

cellset(lookup::DGG.CellLookup) = cellset(parent(lookup))

# The subtree iterators carry the system and the leaf level that a bare vector
# of ids lacks.
const SubtreeIterator = Union{DGG.EdgeCellIterator, DGG.InnerCellIterator}

cellset(iterator::SubtreeIterator) =
    CellSet(DGG.levelgrid(iterator.system, iterator.level), collect(iterator))

# `cell_boundary` on a system covers ids of any level, so a system paired with a
# bare vector of ids needs no level agreement.
cellset(source, ids::AbstractVector) = CellSet(source, ids)
