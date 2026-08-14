Makie.convert_arguments(P::Makie.PointBased, set::DGG.MultiOrderCellSet) =
    Makie.convert_arguments(P,
        GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygons(set)))
