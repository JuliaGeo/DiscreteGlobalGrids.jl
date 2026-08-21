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

# # Naming a set of cells with adjacency
#
# A [`CellSet`](@ref) is enough to draw cells one boundary at a time.  A surface
# also needs which cells touch which, so its set must be one `adjacency` can
# answer about: a single level, positions doubling as data indices.

"""
    CellRegion(region, source, cells)

A set of DGGS cells whose adjacency can be read.

  * `region` — what `DiscreteGlobalGrids.adjacency` is asked about; its positions
    `1:length(region)` index a colour vector.
  * `source` — what `cell_centroid` is asked of.
  * `cells` — `cells[p]` is the cell at position `p`.

Build one with [`cellregion`](@ref) rather than by hand.
"""
struct CellRegion{R, S, C <: AbstractVector}
    region::R
    source::S
    cells::C
end

Base.length(cr::CellRegion) = length(cr.cells)

"""
    cellregion(x) -> CellRegion

Read `x` as a set of cells with adjacency, for [`triangulate`](@ref).

Accepts a `CellRegion`, an `AbstractGrid` — `PartialGrid` included — a
`CellVector`, or a `CellLookup`.

The rest of what [`cellset`](@ref) takes has no adjacency to read: a
`MultiOrderCellSet` spans several levels, so there is no one level to measure it
on, and a bare vector of ids names no set for a neighbour to be inside or
outside of.
"""
function cellregion end

cellregion(cr::CellRegion) = cr
cellregion(grid::DGG.AbstractGrid) = CellRegion(grid, grid, GridCells(grid))

# A `CellVector` is the region, but not a geometry source, so centroids come
# from the level grid it indexes into.
cellregion(cv::DGG.CellVector) =
    CellRegion(cv, DGG.levelgrid(DGG.system(cv), DGG.level(cv)), cv)

cellregion(lookup::DGG.CellLookup) = cellregion(parent(lookup))

cellregion(x) = throw(ArgumentError("$(typeof(x)) names cells but not their \
    adjacency, so it has no surface. `dggsurface` takes a grid, a `PartialGrid`, \
    a `CellVector` or a `CellLookup`; `dggpoly` draws anything."))
