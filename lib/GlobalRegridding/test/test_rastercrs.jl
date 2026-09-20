# CRS-derived RasterGrid charts through the Rasters and Rasters+Proj extensions.

import Rasters
import Rasters: EPSG, ProjString, X, Y

const GRRastersExt = Base.get_extension(GlobalRegridding, :GlobalRegriddingRastersExt)
const GRRastersProjExt = Base.get_extension(GlobalRegridding, :GlobalRegriddingRastersProjExt)

@testset "CRS-derived raster charts" begin
    @test GRRastersExt !== nothing
    @test GRRastersProjExt !== nothing

    values = reshape(collect(1.0:9.0), 3, 3)
    mercator = (X(-1_000_000.0:1_000_000.0:1_000_000.0),
        Y(4_000_000.0:1_000_000.0:6_000_000.0))
    geographic = (X(-1.0:1.0:1.0), Y(0.0:1.0:2.0))
    target = RasterGrid(DD.DimArray(zeros(36, 18), (X(-180.0:10.0:170.0), Y(-90.0:10.0:80.0))))

    @testset "projected CRS builds the Proj chart" begin
        raster = Rasters.Raster(values, mercator; crs = EPSG(3857))
        @test Rasters.crs(raster) == EPSG(3857)
        @test DD.lookup(raster, X) isa Rasters.Projected

        space = RasterGrid(raster)
        @test space.native_to_unit_sphere isa GRProjExt._NativeToUnitSphere
        @test space.unit_sphere_to_native isa GRProjExt._UnitSphereToNative
        @test hascellchart(space)
        @test all(cellat(space, cellcentroid(space, i)) == i for i in 1:ncells(space))

        template = Proj.Transformation("EPSG:3857", "EPSG:4326"; always_xy = true)
        explicit = RasterGrid(DD.DimArray(values, mercator);
            native_to_unit_sphere = template)
        derived = GR.wholeblock(Conservative(), target, space).weights
        reference = GR.wholeblock(Conservative(), target, explicit).weights
        @test derived.colptr == reference.colptr
        @test derived.rowval == reference.rowval
        @test all(derived.nzval .=== reference.nzval)

        # `regrid` reaches the same constructor with no keyword.
        out = regrid(raster; to = target, method = Conservative())
        @test out isa Rasters.AbstractRaster
        @test size(out) == (36, 18)

        # A dimension tuple carrying the CRS resolves the same way.
        @test RasterGrid(DD.dims(raster)).native_to_unit_sphere isa
              GRProjExt._NativeToUnitSphere

        # A bare projected PROJ string is classified without the `+type=crs`
        # tag and charted through Proj.
        merc = RasterGrid(Rasters.Raster(values, mercator;
            crs = ProjString("+proj=merc +datum=WGS84")))
        @test merc.native_to_unit_sphere isa GRProjExt._NativeToUnitSphere
        @test all(cellat(merc, cellcentroid(merc, i)) == i for i in 1:ncells(merc))

        # An explicit chart wins over CRS metadata.
        kept = RasterGrid(raster; native_to_unit_sphere = template)
        @test kept.native_to_unit_sphere.state.template === template
    end

    @testset "geographic CRS keeps the built-in chart" begin
        for raster in (Rasters.Raster(values, geographic; crs = EPSG(4326)),
                       Rasters.Raster(values, geographic))
            space = RasterGrid(raster)
            @test space.native_to_unit_sphere isa GO.UnitSpherical.UnitSphereFromGeographic
            @test space.xperiod == 360.0
            @test space.tables !== nothing
        end

        # Proj classifies a geographic CRS the Proj-free rules do not recognize.
        @test !GRRastersExt._geographic_without_proj(EPSG(4269))
        nad83 = RasterGrid(Rasters.Raster(values, geographic; crs = EPSG(4269)))
        @test nad83.native_to_unit_sphere isa GO.UnitSpherical.UnitSphereFromGeographic

        # `Mapped` lookup values live in `mappedcrs`.
        mapped = Rasters.Raster(values,
            (X(Rasters.Mapped(-1.0:1.0:1.0; crs = EPSG(3857), mappedcrs = EPSG(4326))),
             Y(Rasters.Mapped(0.0:1.0:2.0; crs = EPSG(3857), mappedcrs = EPSG(4326)))))
        @test RasterGrid(mapped).native_to_unit_sphere isa
              GO.UnitSpherical.UnitSphereFromGeographic

        # A `Mapped` lookup naming a projected `crs` but no `mappedcrs` leaves
        # its value CRS unknown.
        unmapped = Rasters.Raster(values,
            (X(Rasters.Mapped(-1.0:1.0:1.0; crs = EPSG(3857), mappedcrs = nothing)),
             Y(Rasters.Mapped(0.0:1.0:2.0; crs = EPSG(3857), mappedcrs = nothing))))
        @test_throws ArgumentError RasterGrid(unmapped)
        @test_throws "mappedcrs" RasterGrid(unmapped)

        mixed = Rasters.Raster(values,
            (X(Rasters.Projected(-1.0:1.0:1.0; crs = EPSG(3857))),
             Y(Rasters.Projected(0.0:1.0:2.0; crs = EPSG(4326)))))
        @test_throws ArgumentError RasterGrid(mixed)
    end

    @testset "prime meridians and angular units preserve latitude convention" begin
        paris = 2 + 20 / 60 + 14.025 / 3600
        coordinates = (X(-1.0:1.0:1.0), Y(49.0:1.0:51.0))
        for (crs, scale) in (
                (ProjString("+proj=longlat +datum=WGS84 +pm=paris"), 1.0),
                (ProjString("+proj=longlat +datum=WGS84 +pm=paris +type=crs"), 1.0),
                (ProjString("+proj=longlat +ellps=clrk80ign +towgs84=-168,-60,320 +pm=paris"), 1.0),
                (EPSG(4807), 0.9),
                (Rasters.GeoFormatTypes.WellKnownText(Proj.CRS("EPSG:4807")), 0.9),
                (Rasters.GeoFormatTypes.WellKnownText2(Proj.CRS("EPSG:4807")), 0.9))
            @test !GRRastersExt._geographic_without_proj(crs)
            raster = Rasters.Raster(values, coordinates; crs)
            space = RasterGrid(raster)
            @test space.native_to_unit_sphere isa GRRastersProjExt._GeographicChart
            for lon in (-1.0, 0.0, 1.0), lat in (49.0, 50.0, 51.0)
                point = space.native_to_unit_sphere((lon, lat))
                expected = GEOGRAPHIC_TO_UNIT_SPHERE((scale * lon + paris, scale * lat))
                @test all(isapprox.(Tuple(point), Tuple(expected); atol = 1e-9))
                @test all(isapprox.(space.unit_sphere_to_native(point), (lon, lat); atol = 1e-10))
            end
            @test all(cellat(space, cellcentroid(space, i)) == i for i in 1:ncells(space))
            @test space.xperiod ≈ 360 / scale
            @test all(isapprox.(GR._spherical_step_bounds_radians(space.native_to_unit_sphere, 1.0, 1.0),
                (deg2rad(scale), deg2rad(scale))))
            world = RasterGrid(Rasters.Raster(ones(36, 18),
                (X((-175.0:10.0:175.0) ./ scale), Y((-85.0:10.0:85.0) ./ scale)); crs))
            @test all(cellat(world, cellcentroid(world, i)) == i for i in 1:ncells(world))
            seam = RasterGrid(Rasters.Raster(ones(4, 3),
                (X((177.0:1.0:180.0) ./ scale), Y((49.0:1.0:51.0) ./ scale)); crs))
            @test all(cellat(seam, cellcentroid(seam, i)) == i for i in 1:ncells(seam))
            shifted = RasterGrid(DD.DimArray(values,
                (X(scale .* (-1.0:1.0:1.0) .+ paris), Y(scale .* (49.0:1.0:51.0)))))
            @test parent(regrid(raster; to = shifted, method = DirectNearest(), lazy = false)) ≈ values
        end
        for pm in ("greenwich", "0")
            crs = ProjString("+proj=longlat +datum=WGS84 +pm=$pm")
            @test RasterGrid(Rasters.Raster(values, geographic; crs)).native_to_unit_sphere isa
                GO.UnitSpherical.UnitSphereFromGeographic
        end
    end

    @testset "Proj-free classification and error" begin
        @test GRRastersExt._geographic_without_proj(EPSG(4326))
        @test !GRRastersExt._geographic_without_proj(EPSG(3857))
        @test GRRastersExt._geographic_without_proj(ProjString("+proj=longlat +datum=WGS84 +no_defs"))
        @test !GRRastersExt._geographic_without_proj(ProjString("+proj=merc +datum=WGS84"))
        @test !GRRastersExt._geographic_without_proj(
            Rasters.WellKnownText(Rasters.GeoFormatTypes.CRS(), "GEOGCS[\"WGS 84\",DATUM[]]"))
        @test !GRRastersExt._geographic_without_proj(
            Rasters.WellKnownText(Rasters.GeoFormatTypes.CRS(), "PROJCS[\"WGS 84 / Pseudo-Mercator\"]"))

        # With Proj loaded, the GeoFormat method of the hook shadows the
        # Proj-free throw, so the generic method is exercised directly.
        @test_throws ArgumentError GR._projected_crs_native_to_unit_sphere("EPSG:3857")
        err = try
            GR._projected_crs_native_to_unit_sphere(EPSG(3857).val)
        catch e
            e
        end
        @test occursin("Load Proj", err.msg)
        @test occursin("native_to_unit_sphere", err.msg)
    end
end
