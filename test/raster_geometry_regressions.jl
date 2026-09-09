module RasterGeometryRegressions

using Test
using Random
import DiscreteGlobalGrids as DGG
import GeoInterface as GI
import GeometryOps as GO

function oracle_indices(prepared, grid, boundary)
    out = Int[]
    for i in 1:DGG.ncells(grid)
        predicate = prepared.polygon && boundary === :inside ?
            GO.pred_contains() : GO.pred_intersects()
        c = DGG.cellindex(grid, i)
        geometry = prepared.polygon && boundary === :center ?
            DGG.cell_centroid(grid, c) : DGG.cell_polygon(grid, c)
        GO.relate_predicate(prepared.target.prepared, predicate, geometry) && push!(out, i)
    end
    return out
end

function check_selection(grid, geometry, boundary; label="")
    prepared = DGG._prepare_raster_geometry(geometry)
    prefix = "$label $(typeof(DGG.system(grid))) level $(DGG.level(grid)) $boundary"
    actual = try
        DGG._raster_indices(grid, prepared; boundary)
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
    # The unique minor arc runs east along the equator. Z7Cell("1006")
    # (local index 416) straddles it: its vertices range from -6.39° to
    # +6.39° latitude and 163.98° to 172.79° longitude.
    grid = DGG.levelgrid(DGG.IGeo7System(), 2)
    almost_antipodal = GI.LineString([(0.0, 0.0), (179.999999, 0.0)])
    @test DGG._raster_indices(grid, almost_antipodal; boundary=:intersects) ==
        [162, 163, 165, 166, 168, 170, 173, 174, 185, 187, 188, 191,
         234, 237, 238, 385, 389, 406, 411, 412, 416, 432, 435, 436, 450]

    # The cells at indices 71 and 74 have (0, 0) as a vertex. Each therefore
    # contains a boundary point of this tiny hole and cannot be wholly inside
    # the polygon. This exercises spherical proper crossings whose endpoint
    # chords have a singular projection into the XY plane.
    grid = DGG.levelgrid(DGG.HEALPixSystem(), 2)
    with_equatorial_hole = deterministic_geometries()[5]
    @test DGG._raster_indices(grid, with_equatorial_hole; boundary=:inside) ==
        [66, 67, 69, 72, 73, 76, 78, 79]
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
                # Covered by the focused tests above. On the affected
                # GeometryOps release these throw before an oracle exists.
                system isa DGG.IGeo7System && j == 10 && continue
                system isa DGG.HEALPixSystem && j == 5 && boundary === :inside && continue
                @test check_selection(grid, geometry, boundary;
                    label="geometry $j")
            end
        end
    end
end

end
