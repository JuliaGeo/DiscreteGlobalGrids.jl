# ## Grids

# The abstract interface covers wrappers and downstream grid implementations.
function Makie.convert_arguments(P::Makie.PointBased, grid::DGG.AbstractGrid)
    Makie.convert_arguments(P,
        GO.transform(GO.GeographicFromUnitSphere(), DGG.getcell(grid)))
end

function Makie.convert_arguments(P::Type{<:Makie.Poly}, grid::DGG.AbstractGrid)
    Makie.convert_arguments(P,
        GO.transform(GO.GeographicFromUnitSphere(), DGG.getcell(grid)))
end

# ## Multi-order queries and results

# A coverage plots its query target; a resolved cell set plots selected cells.
function Makie.convert_arguments(P::Makie.PointBased, coverage::DGG.MultiOrderCoverage)
    Makie.convert_arguments(P, parent(coverage))
end

function Makie.convert_arguments(P::Type{<:Makie.Poly}, coverage::DGG.MultiOrderCoverage)
    Makie.convert_arguments(P, parent(coverage))
end

function Makie.convert_arguments(P::Makie.PointBased, set::DGG.MultiOrderCellSet)
    Makie.convert_arguments(P, GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygons(set)))
end

function Makie.convert_arguments(P::Type{<:Makie.Poly}, set::DGG.MultiOrderCellSet)
    Makie.convert_arguments(P, GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygons(set)))
end

# ## Cell collections

function Makie.convert_arguments(P::Makie.PointBased, vector::DGG.AbstractCellVector)
    Makie.convert_arguments(P, DGG.PartialGrid(vector))
end

function Makie.convert_arguments(P::Type{<:Makie.Poly}, vector::DGG.AbstractCellVector)
    Makie.convert_arguments(P, DGG.PartialGrid(vector))
end

function Makie.convert_arguments(P::Makie.PointBased, lookup::DGG.AbstractCellLookup)
    Makie.convert_arguments(P, parent(lookup))
end

function Makie.convert_arguments(P::Type{<:Makie.Poly}, lookup::DGG.AbstractCellLookup)
    Makie.convert_arguments(P, parent(lookup))
end

# Resolve mixed-level polygons through each cell's own level grid.
function Makie.convert_arguments(P::Makie.PointBased, vector::DGG.MultiOrderVector)
    Makie.convert_arguments(P,
        GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygons(vector)))
end

function Makie.convert_arguments(P::Type{<:Makie.Poly}, vector::DGG.MultiOrderVector)
    Makie.convert_arguments(P,
        GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygons(vector)))
end

function Makie.convert_arguments(P::Makie.PointBased, lookup::DGG.MultiOrderLookup)
    Makie.convert_arguments(P, parent(lookup))
end

function Makie.convert_arguments(P::Type{<:Makie.Poly}, lookup::DGG.MultiOrderLookup)
    Makie.convert_arguments(P, parent(lookup))
end

# Iterator metadata makes the partial-grid interpretation unambiguous.
const SubtreeIterator = Union{DGG.EdgeCellIterator,DGG.InnerCellIterator}

_partial_grid(iterator::SubtreeIterator) =
    DGG.PartialGrid(iterator.system, iterator.level, collect(iterator))

function Makie.convert_arguments(P::Makie.PointBased, iterator::SubtreeIterator)
    Makie.convert_arguments(P, _partial_grid(iterator))
end

function Makie.convert_arguments(P::Type{<:Makie.Poly}, iterator::SubtreeIterator)
    Makie.convert_arguments(P, _partial_grid(iterator))
end
