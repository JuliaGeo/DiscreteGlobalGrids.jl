# ## Grids

# Every grid publishes the same positional geometry interface.  Defining the
# conversion at that interface also covers wrappers such as `AuthalicGrid` and
# grid implementations supplied by downstream packages.
function Makie.convert_arguments(P::Makie.PointBased, grid::DGG.AbstractGrid)
    Makie.convert_arguments(P,
        GO.transform(GO.GeographicFromUnitSphere(), DGG.getcell(grid)))
end

function Makie.convert_arguments(P::Type{<:Makie.Poly}, grid::DGG.AbstractGrid)
    Makie.convert_arguments(P,
        GO.transform(GO.GeographicFromUnitSphere(), DGG.getcell(grid)))
end

# ## Multi-order queries and results

# A coverage is the query specification, so plotting it means plotting its
# target.  The selected DGGS cells belong to the `MultiOrderCellSet` returned
# after a system and a resolution or budget have been supplied.
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

# The subtree iterators carry the system and leaf level which a plain vector of
# cell ids lacks, so they can be read unambiguously as partial grids.
const SubtreeIterator = Union{DGG.EdgeCellIterator,DGG.InnerCellIterator}

_partial_grid(iterator::SubtreeIterator) =
    DGG.PartialGrid(iterator.system, iterator.level, collect(iterator))

function Makie.convert_arguments(P::Makie.PointBased, iterator::SubtreeIterator)
    Makie.convert_arguments(P, _partial_grid(iterator))
end

function Makie.convert_arguments(P::Type{<:Makie.Poly}, iterator::SubtreeIterator)
    Makie.convert_arguments(P, _partial_grid(iterator))
end
