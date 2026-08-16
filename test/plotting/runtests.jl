module PlottingTests

using Test
import DiscreteGlobalGrids as DGG
import GeoInterface as GI
import GeometryOps as GO
import Makie

geographic_polygons(grid) =
    GO.transform(GO.GeographicFromUnitSphere(), DGG.getcell(grid))

function test_conversions(actual, expected)
    @test isequal(Makie.convert_arguments(Makie.PointBased(), actual),
        Makie.convert_arguments(Makie.PointBased(), expected))
    @test isequal(Makie.convert_arguments(Makie.Poly, actual),
        Makie.convert_arguments(Makie.Poly, expected))
end

@testset "multi-order Makie conversion" begin
    system = DGG.HEALPixSystem()
    cell = first(DGG.rootcells(system))
    set = DGG.MultiOrderCellSet(system, [cell], [1], trues(1), DGG.level(cell))
    polygons = GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygons(set))

    @test Makie.convert_arguments(Makie.PointBased(), set) ==
        Makie.convert_arguments(Makie.PointBased(), polygons)
    @test Makie.convert_arguments(Makie.Poly, set) ==
        Makie.convert_arguments(Makie.Poly, polygons)
end

@testset "grid Makie conversion" begin
    system = DGG.HEALPixSystem()
    grid = DGG.levelgrid(system, 1)
    root = first(DGG.rootcells(system))
    partial = DGG.PartialGrid(system, root, 1)
    authalic = DGG.levelgrid(DGG.AuthalicSystem(system), 1)

    for candidate in (grid, partial, authalic)
        test_conversions(candidate, geographic_polygons(candidate))
    end
end

@testset "cell collection Makie conversion" begin
    system = DGG.HEALPixSystem()
    root = first(DGG.rootcells(system))
    set = DGG.MultiOrderCellSet(system, [root], [1], trues(1), 2)
    vectors = (
        DGG.CellVector(DGG.levelgrid(system, 1)),
        DGG.CellVector(DGG.levelgrid(DGG.AuthalicSystem(system), 1)),
        DGG.CellVector(set),
    )

    for vector in vectors
        expected = geographic_polygons(DGG.PartialGrid(vector))
        test_conversions(vector, expected)
        test_conversions(DGG.CellLookup(vector), expected)
    end

    # Plot the cells represented by the expanded vector, not its mixed-level
    # provenance.
    expanded = last(vectors)
    @test length(expanded) == 16
    polygons = only(Makie.convert_arguments(Makie.Poly, expanded))
    @test length(polygons) == length(expanded)
end

@testset "multi-order container Makie conversion" begin
    system = DGG.HEALPixSystem()
    root = first(DGG.rootcells(system))
    # One child refined further than its siblings, so a conversion reading the
    # container through a single level's grid would name the wrong cells.
    kids = DGG.children(system, root)
    mixed = vcat(DGG.children(system, first(kids)), kids[2:end])
    vector = DGG.MultiOrderVector(system, mixed)
    @test length(unique(DGG.level, collect(vector))) > 1

    polygons = GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygons(vector))
    @test length(polygons) == length(vector)

    test_conversions(vector, polygons)
    test_conversions(DGG.MultiOrderLookup(vector), polygons)
end

@testset "multi-order coverage Makie conversion" begin
    target = GI.Polygon([GI.LinearRing([
        (-10.0, 40.0), (10.0, 40.0), (10.0, 50.0),
        (-10.0, 50.0), (-10.0, 40.0),
    ])])
    test_conversions(DGG.MultiOrderCoverage(target), target)
end

@testset "subtree iterator Makie conversion" begin
    system = DGG.HEALPixSystem()
    root = first(DGG.rootcells(system))

    for iterator in (DGG.EdgeCellIterator(system, root, 2),
                     DGG.InnerCellIterator(system, root, 2))
        grid = DGG.PartialGrid(iterator.system, iterator.level, collect(iterator))
        test_conversions(iterator, geographic_polygons(grid))
    end
end

end
