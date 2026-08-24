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
    ==(a::CellSet, b::CellSet)

Do the two name the same cells, read out of the same source?

The default for a struct is `===`, which compares a field that is not `isbits`
by pointer — and a set rebuilt from the same cells is a different object with
the same contents, which is what a plot handed its arguments a second time gets.
Comparing the fields with `==` instead is what lets it see that nothing changed:
a window-backed `CellVector` answers in about a microsecond however many cells
it holds, and a [`GridCells`](@ref) in one comparison.

A bare vector of ids has no such shortcut and is compared element by element, so
`===` is tried first and the walk is only the price of an honest "no".
"""
Base.:(==)(a::CellSet, b::CellSet) =
    a === b || (a.source == b.source && a.cells == b.cells)

Base.hash(cs::CellSet, h::UInt) = hash(cs.cells, hash(cs.source, hash(:CellSet, h)))

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

# The grid is the whole of it, so two of them agree or not in one comparison
# rather than in `ncells` of them.  Without this the `AbstractVector` fallback
# walks every cell of both, which at a whole level is millions of `cellindex`
# calls to answer a question the grid already answers.
Base.:(==)(a::GridCells, b::GridCells) = a.grid == b.grid
Base.hash(gc::GridCells, h::UInt) = hash(gc.grid, hash(:GridCells, h))

"""
    cellset(x) -> CellSet
    cellset(system, ids) -> CellSet

Read `x` as a set of cells to draw.

Accepts a `CellSet`, an `AbstractGrid`, any `AbstractCellVector` or
`AbstractCellLookup` — so a cube axis read from a store as readily as one built
in memory — a `MultiOrderCellSet` or `MultiOrderCoverage`, one of the subtree
iterators, or a system together with a vector of ids.  Supporting a new
container is a method here.
"""
function cellset end

cellset(cs::CellSet) = cs
cellset(grid::DGG.AbstractGrid) = CellSet(grid, GridCells(grid))
cellset(set::DGG.MultiOrderCellSet) = CellSet(set, set.cells)
cellset(coverage::DGG.MultiOrderCoverage) = cellset(parent(coverage))

# A cell vector knows its system and its level but is not itself a geometry
# source, so pair it with the level grid it indexes into.
cellset(cv::DGG.AbstractCellVector) =
    CellSet(DGG.levelgrid(DGG.system(cv), DGG.level(cv)), cv)

cellset(lookup::DGG.AbstractCellLookup) = cellset(parent(lookup))

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
    ==(a::CellRegion, b::CellRegion)

Do the two name the same cells, with the same adjacency to read?

See `==(::CellSet, ::CellSet)` for why the default will not do.
"""
Base.:(==)(a::CellRegion, b::CellRegion) =
    a === b || (a.region == b.region && a.source == b.source && a.cells == b.cells)

Base.hash(cr::CellRegion, h::UInt) =
    hash(cr.cells, hash(cr.source, hash(cr.region, hash(:CellRegion, h))))

"""
    cellregion(x) -> CellRegion

Read `x` as a set of cells with adjacency, for [`triangulate`](@ref).

Accepts a `CellRegion`, an `AbstractGrid` — `PartialGrid` included — or any
`AbstractCellVector` or `AbstractCellLookup`.

The rest of what [`cellset`](@ref) takes has no adjacency to read: a
`MultiOrderCellSet` spans several levels, so there is no one level to measure it
on, and a bare vector of ids names no set for a neighbour to be inside or
outside of.
"""
function cellregion end

cellregion(cr::CellRegion) = cr
cellregion(grid::DGG.AbstractGrid) = CellRegion(grid, grid, GridCells(grid))

# A cell vector is the region, but not a geometry source, so centroids come
# from the level grid it indexes into.
cellregion(cv::DGG.AbstractCellVector) =
    CellRegion(cv, DGG.levelgrid(DGG.system(cv), DGG.level(cv)), cv)

cellregion(lookup::DGG.AbstractCellLookup) = cellregion(parent(lookup))

cellregion(x) = throw(ArgumentError("$(typeof(x)) names cells but not their \
    adjacency, so it has no surface. `dggsurface` takes a grid, a `PartialGrid`, \
    or a cell vector or lookup; `dggpoly` draws anything."))

# # A container, or a source and its ids
#
# A recipe with one argument never needs to tell those apart: it converts what
# it was given and dispatches on what came out.  One taking `(cells, zs)` does —
# `dggresample(cells, heights)` and `dggresample(system, ids)` are the same two
# slots — and the question is answered by the first argument's type.

"""
    CellContainer

Everything [`cellset`](@ref) reads on its own, as a type.

What it leaves out is the two-argument form, `cellset(system, ids)`, where the
cells are named by a pair rather than by one object.
"""
const CellContainer = Union{CellSet, CellRegion, DGG.AbstractGrid, DGG.CellVector,
    DGG.CellLookup, DGG.MultiOrderCellSet, DGG.MultiOrderCoverage, SubtreeIterator}
