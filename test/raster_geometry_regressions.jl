module RasterGeometryRegressions

using Test
using Random
import DiscreteGlobalGrids as DGG
import GeoInterface as GI
import GeometryOps as GO

# The oracle is GeometryOps' spherical RelateNG asked of every cell, no tree.
prepare(geometry) = (; prepared = DGG.Engine._query_target(geometry).prepared,
    polygon = GI.trait(geometry) isa Union{GI.PolygonTrait,GI.MultiPolygonTrait})

function oracle_indices(prepared, grid, boundary)
    out = Int[]
    for i in 1:DGG.ncells(grid)
        predicate = prepared.polygon && boundary === :inside ?
            GO.pred_contains() : GO.pred_intersects()
        c = DGG.cellindex(grid, i)
        geometry = prepared.polygon && boundary === :center ?
            DGG.cell_centroid(grid, c) : DGG.cell_polygon(grid, c)
        GO.relate_predicate(prepared.prepared, predicate, geometry) && push!(out, i)
    end
    return out
end

function check_selection(grid, geometry, boundary; label="")
    prepared = prepare(geometry)
    prefix = "$label $(typeof(DGG.system(grid))) level $(DGG.level(grid)) $boundary"
    actual = try
        DGG._raster_indices(grid, geometry; boundary)
    catch exception
        error("$prefix accelerated selection threw $(sprint(showerror, exception))")
    end
    expected = try
        oracle_indices(prepared, grid, boundary)
    catch exception
        error("$prefix oracle threw $(sprint(showerror, exception))")
    end
    if actual != expected
        missing_indices = setdiff(expected, actual)
        extra_indices = setdiff(actual, expected)
        error("$prefix: " *
              "missing=$(missing_indices[1:min(end, 8)]) extra=$(extra_indices[1:min(end, 8)])")
    end
    return true
end

closed_ring(points) = GI.LinearRing([points; first(points)])
polygon(points) = GI.Polygon([closed_ring(points)])

function deterministic_geometries()
    geometries = Any[
        polygon([(-0.0001, -0.0001), (0.0001, -0.0001),
                 (0.0001, 0.0001), (-0.0001, 0.0001)]),
        polygon([(179.9999, -0.0001), (180.0001, -0.0001),
                 (180.0001, 0.0001), (179.9999, 0.0001)]),
        polygon([(-160.0, 82.0), (-40.0, 82.0), (80.0, 82.0),
                 (160.0, 82.0)]),
        polygon([(-170.0, -70.0), (170.0, -70.0), (170.0, 70.0),
                 (-170.0, 70.0)]),
        GI.Polygon([
            closed_ring([(-35.0, -30.0), (35.0, -30.0),
                         (35.0, 30.0), (-35.0, 30.0)]),
            closed_ring([(-0.001, -0.001), (-0.001, 0.001),
                         (0.001, 0.001), (0.001, -0.001)]),
        ]),
        GI.MultiPolygon([
            polygon([(-175.0, -5.0), (-165.0, -5.0),
                     (-165.0, 5.0), (-175.0, 5.0)]),
            polygon([(165.0, -5.0), (175.0, -5.0),
                     (175.0, 5.0), (165.0, 5.0)]),
        ]),
        GI.LineString([(-179.999, 0.0), (179.999, 0.0)]),
        GI.LineString([(0.0, 89.999), (120.0, 89.999), (-120.0, 89.999)]),
        GI.LineString([(12.0, 20.0), (12.0, 20.0)]),
        GI.LineString([(0.0, 0.0), (179.999999, 0.0)]),
    ]
    return geometries
end

function random_geometries(rng, count)
    geometries = Any[]
    for _ in 1:count
        lon = rand(rng) * 360 - 180
        lat = rand(rng) * 150 - 75
        dx = 10.0 ^ (rand(rng) * 4 - 3)
        dy = 10.0 ^ (rand(rng) * 3 - 3)
        push!(geometries, polygon([(lon-dx, lat-dy), (lon+dx, lat-dy),
                                   (lon+dx, lat+dy), (lon-dx, lat+dy)]))
        lon2 = lon + (rand(rng) * 100 - 50)
        lat2 = clamp(lat + (rand(rng) * 80 - 40), -89.99, 89.99)
        push!(geometries, GI.LineString([(lon, lat), (lon2, lat2)]))
    end
    return geometries
end

@testset "spherical crossing projection regressions" begin
    equator = GI.LineString([(0.0, 0.0), (10.0, 0.0)])
    meridian = GI.LineString([(5.0, -1.0), (5.0, 1.0)])
    equator_prepared = GO.prepare(GO.RelateNG(; manifold=GO.Spherical()), equator)
    @test GO.relate_predicate(equator_prepared, GO.pred_intersects(), meridian)

    # The minor arc runs east along the equator; Z7Cell("1006") (local index
    # 416) straddles it, spanning ±6.39° latitude at 163.98°-172.79° longitude.
    grid = DGG.levelgrid(DGG.IGeo7System(), 2)
    almost_antipodal = GI.LineString([(0.0, 0.0), (179.999999, 0.0)])
    almost_antipodal_prepared = prepare(almost_antipodal)
    @test GO.relate_predicate(almost_antipodal_prepared.prepared,
        GO.pred_intersects(), DGG.cell_polygon(grid, DGG.cellindex(grid, 416)))
    @test check_selection(grid, almost_antipodal, :intersects;
        label="almost antipodal equator line")

    # Cells 71 and 74 have (0, 0) as a vertex, a boundary point of the tiny hole,
    # so neither is inside the polygon. Their edge chords project singularly
    # into the XY plane, the proper-crossing case under test.
    grid = DGG.levelgrid(DGG.HEALPixSystem(), 2)
    with_equatorial_hole = deterministic_geometries()[5]
    with_hole_prepared = prepare(with_equatorial_hole)
    for i in (71, 74)
        @test !GO.relate_predicate(with_hole_prepared.prepared,
            GO.pred_contains(), DGG.cell_polygon(grid, DGG.cellindex(grid, i)))
    end
    @test check_selection(grid, with_equatorial_hole, :inside;
        label="equatorial hole")
end

@testset "raster geometry regressions" begin
    rng = MersenneTwister(0x59)
    systems = (DGG.H3System(), DGG.IGeo7System(), DGG.HEALPixSystem(),
               DGG.ISEA4RSystem(), DGG.A5System(), DGG.S2System())
    geometries = [deterministic_geometries(); random_geometries(rng, 12)]
    for system in systems
        grid = DGG.levelgrid(system, first(DGG.levels(system)) + 2)
        for (j, geometry) in pairs(geometries)
            for boundary in (:inside, :center, :intersects)
                @test check_selection(grid, geometry, boundary;
                    label="geometry $j")
            end
        end
    end
end

end
