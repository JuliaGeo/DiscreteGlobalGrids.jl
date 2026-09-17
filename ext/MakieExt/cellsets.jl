_lonlat(x) = GO.transform(GO.GeographicFromUnitSphere(), x)

# Iterator metadata makes the partial-grid interpretation unambiguous.
const SubtreeIterator = Union{DGG.EdgeCellIterator,DGG.InnerCellIterator}

_partial_grid(iterator::SubtreeIterator) =
    DGG.PartialGrid(iterator.system, iterator.level, collect(iterator))

# Point-based recipes and `poly` convert identically; two signatures keep each
# method more specific than Makie's own `AbstractVector` conversions.
for P in (:(Makie.PointBased), :(Type{<:Makie.Poly}))
    @eval begin
        # The abstract interface covers wrappers and downstream grid implementations.
        Makie.convert_arguments(P::$P, grid::DGG.AbstractGrid) =
            Makie.convert_arguments(P, _lonlat(DGG.getcell(grid)))
        # A coverage plots its query target; a resolved cell set plots selected cells.
        Makie.convert_arguments(P::$P, coverage::DGG.MultiOrderCoverage) =
            Makie.convert_arguments(P, parent(coverage))
        Makie.convert_arguments(P::$P, set::DGG.MultiOrderCellSet) =
            Makie.convert_arguments(P, _lonlat(DGG.cell_polygons(set)))
        Makie.convert_arguments(P::$P, vector::DGG.AbstractCellVector) =
            Makie.convert_arguments(P, DGG.PartialGrid(vector))
        # A cell axis, single-level or mixed, plots the container it wraps.
        Makie.convert_arguments(P::$P, lookup::DGG.CellLookups.AbstractCellAxis) =
            Makie.convert_arguments(P, parent(lookup))
        # Mixed-level polygons resolve through each cell's own level grid.
        Makie.convert_arguments(P::$P, vector::DGG.MultiOrderVector) =
            Makie.convert_arguments(P, _lonlat(DGG.cell_polygons(vector)))
        Makie.convert_arguments(P::$P, iterator::SubtreeIterator) =
            Makie.convert_arguments(P, _partial_grid(iterator))
    end
end
