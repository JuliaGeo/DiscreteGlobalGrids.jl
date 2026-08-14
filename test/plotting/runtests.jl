module PlottingTests

using Test
import DiscreteGlobalGrids as DGG
import GeometryOps as GO
import Makie

@testset "multi-order Makie conversion" begin
    system = DGG.HEALPixSystem()
    cell = first(DGG.rootcells(system))
    set = DGG.MultiOrderCellSet(system, [cell], [1], trues(1), DGG.level(cell))
    polygons = GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygons(set))

    @test Makie.convert_arguments(Makie.PointBased(), set) ==
        Makie.convert_arguments(Makie.PointBased(), polygons)
end

end
