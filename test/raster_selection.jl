module RasterSelectionTests
using Test
import DiscreteGlobalGrids as DGG
import GeoInterface as GI
import GeometryOps as GO

struct EmptyRasterLine end
GI.geomtrait(::EmptyRasterLine) = GI.LineStringTrait()
GI.getpoint(::GI.LineStringTrait, ::EmptyRasterLine) = ()

# The oracle is GeometryOps' spherical RelateNG asked of every cell, no tree.
prepare(geometry) = (; prepared = DGG.Engine._query_target(geometry).prepared,
    polygon = GI.trait(geometry) isa Union{GI.PolygonTrait,GI.MultiPolygonTrait})

function oracle_match(prep, grid, i, boundary)
    c = DGG.cellindex(grid,i)
    geom = prep.polygon && boundary === :center ? DGG.cell_centroid(grid,c) : DGG.cell_polygon(grid,c)
    predicate = prep.polygon && boundary === :inside ? GO.pred_contains() : GO.pred_intersects()
    return GO.relate_predicate(prep.prepared, predicate, geom)
end

box(x, y, r) = GI.Polygon([GI.LinearRing([(x-r,y-r),(x+r,y-r),
    (x+r,y+r),(x-r,y+r),(x-r,y-r)])])

@testset "Spherical raster selection" begin
    for sys in (DGG.H3System(), DGG.IGeo7System(), DGG.HEALPixSystem(),
            DGG.ISEA4RSystem(), DGG.A5System())
        lev = first(DGG.levels(sys)) + 1
        grid = DGG.levelgrid(sys, lev)
        @testset "$(typeof(sys))" begin
            for geometry in (box(12.0, 20.0, 24.0), box(178.0, 8.0, 12.0),
                    GI.Polygon([GI.LinearRing([(20cos(t),20sin(t)) for t in range(0,2pi;length=129)])]),
                    GI.Polygon([GI.LinearRing([(-120.,75.),(0.,75.),(120.,75.),(-120.,75.)])]),
                    GI.Polygon([GI.LinearRing([(-40.,-40.),(40.,-40.),(40.,40.),(-40.,40.),(-40.,-40.)]),
                        GI.LinearRing([(-5.,-5.),(-5.,5.),(5.,5.),(5.,-5.),(-5.,-5.)])]))
                prep = prepare(geometry)
                selections = Dict{Symbol,Vector{Int}}()
                for boundary in (:inside, :center, :intersects)
                    actual = DGG._raster_indices(grid, geometry; boundary)
                    oracle = [i for i in 1:DGG.ncells(grid)
                        if oracle_match(prep, grid, i, boundary)]
                    @test actual == oracle
                    selections[boundary] = actual
                end
                @test issubset(selections[:inside], selections[:center])
                @test issubset(selections[:center], selections[:intersects])
                @test DGG._raster_indices(grid, geometry; boundary=:touches) == selections[:intersects]
            end
            partial = DGG.PartialGrid(sys, lev, [DGG.cellindex(grid,i) for i in 1:2:DGG.ncells(grid)])
            geom = box(12.,20.,24.)
            prep = prepare(geom)
            @test DGG._raster_indices(partial, geom) == [i for i in 1:DGG.ncells(partial)
                if oracle_match(prep, partial, i, :center)]
            point = GI.Point((12.,20.))
            c = DGG.cellat(grid, 12., 20.)
            @test DGG._raster_indices(grid, point) == [DGG.localindex(grid, c)]
            @test DGG._raster_indices(grid, GI.MultiPoint([point,point])) == [DGG.localindex(grid,c)]
            @test isempty(DGG._raster_indices(grid, missing))
            @test isempty(DGG._raster_indices(grid, EmptyRasterLine()))
            @test_throws ArgumentError DGG._raster_indices(grid, point; boundary=:invalid)
            line = GI.LineString([(-40.,0.),(40.,0.)])
            prep = prepare(line)
            @test DGG._raster_indices(grid, line) == [i for i in 1:DGG.ncells(grid)
                if oracle_match(prep, grid, i, :intersects)]
            # A coincident boundary is covered, and tiny enclosed polygons
            # intersect their host even when they miss the host centroid.
            cell = DGG.cell_polygon(grid,c)
            @test DGG.localindex(grid,c) in DGG._raster_indices(grid, cell; boundary=:inside)
            @test DGG.localindex(grid,c) in DGG._raster_indices(grid, box(12.,20.,0.001); boundary=:intersects)
        end
    end
    @testset "Regional CopDEM tiled cursor" begin
        CD = DGG.CopernicusDEM
        sys = CD.CopernicusDEMSystem{30}()
        ids = sort!(reduce(vcat, [collect(DGG.children(sys, CD.tilecell(sys, lat, lon)))
            for (lat,lon) in ((49,6),(50,6),(50,7))]))
        grid = DGG.PartialGrid(sys, 1, ids)
        @test DGG.treeify(grid) isa DGG.TiledRasterCursor
        geom = box(6.9, 50.0, 0.25)
        prep = prepare(geom)
        for boundary in (:center,:inside,:intersects)
            @test DGG._raster_indices(grid, geom; boundary) ==
                [i for i in 1:DGG.ncells(grid)
                 if oracle_match(prep, grid, i, boundary)]
        end
        @test isempty(DGG._raster_indices(grid, GI.Point((-120.,0.))))
    end

end
end
