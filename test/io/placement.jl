# Where a store's cells are: the IGEO7 placement (vertex 0 at 11.2°E) and the
# geodetic latitude a store can declare. Reading applies both, writing declares
# both, and a placement that is not IGEO7's is refused rather than read as if
# it were.

module DGGSIOPlacementTests

using Test
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: IGeo7System, HEALPixSystem, AuthalicSystem, AuthalicGrid,
    Z7Cell, CellVector, CellLookup, Cells, DGGSFormatError, GridOrientation,
    StoreDescription, levelgrid, descendants, dggread
import DimensionalData as DD
import GeometryOpsCore

const SYS = IGeo7System()
const LEVEL = 3
const ROOTS = [Z7Cell(DGG.idselect(levelgrid(SYS, 1), r)) for r in 1:3]
const CELLS = sort!(reduce(vcat, [collect(descendants(SYS, c, LEVEL)) for c in ROOTS]))

igeo7(; kw...) = StoreDescription(; gridname="igeo7", system=SYS, idscheme=:z7int,
    level=LEVEL, kw...)

@testset "placement and geodetic latitude" begin

    @testset "describedsystem: the declared placement is checked" begin
        @test DGG.describedsystem(igeo7(), SYS) === SYS
        # the published IGEO7 stores print the latitude to 8 decimals
        pori = GridOrientation(vert0_lon=11.2, vert0_lat=58.28252559, vert0_azimuth=0.0)
        @test DGG.describedsystem(igeo7(orientation=pori), SYS) === SYS
        @test DGG.describedsystem(igeo7(orientation=GridOrientation(vert0_lon=11.2)), SYS) === SYS
        for o in (GridOrientation(vert0_lon=11.25),
                  GridOrientation(vert0_lon=11.2, vert0_lat=58.3),
                  GridOrientation(vert0_lon=11.2, vert0_azimuth=10.0))
            err = try
                DGG.describedsystem(igeo7(orientation=o), SYS)
                nothing
            catch e
                e
            end
            @test err isa DGGSFormatError
            @test err.check === :orientation_mismatch
        end
        # other systems carry no placement here and pass through
        @test DGG.describedsystem(StoreDescription(gridname="healpix",
            orientation=GridOrientation(vert0_lon=0.0)), HEALPixSystem()) === HEALPixSystem()
    end

    @testset "describedsystem: declared geodetic latitude wraps the system" begin
        @test DGG.describedsystem(igeo7(geodetic_conversion=false), SYS) === SYS
        sys = DGG.describedsystem(igeo7(geodetic_conversion=true), SYS)
        @test sys isa AuthalicSystem
        @test parent(sys) === SYS
        @test levelgrid(sys, LEVEL) isa AuthalicGrid
    end

    @testset "storedsystem is the inverse on write" begin
        @test DGG.storedsystem(SYS) == (SYS, DGG.IGEO7_GRID_ORIENTATION, false)
        @test DGG.storedsystem(AuthalicSystem(SYS)) == (SYS, DGG.IGEO7_GRID_ORIENTATION, true)
        @test DGG.storedsystem(HEALPixSystem()) == (HEALPixSystem(), nothing, nothing)
        @test DGG.IGEO7_GRID_ORIENTATION.vert0_lon == DGG.IGeo7.IGEO7_VERT0_LON == 11.2
        # the flag says WGS84, so another ellipsoid has no spelling
        other = AuthalicSystem(SYS,
            GeometryOpsCore.Geodesic(; semimajor_axis=6378137.0, inv_flattening=300.0))
        @test_throws ArgumentError DGG.storedsystem(other)
        # describedsystem undoes storedsystem
        for s in (SYS, AuthalicSystem(SYS))
            base, o, g = DGG.storedsystem(s)
            back = DGG.describedsystem(igeo7(orientation=o, geodetic_conversion=g), base)
            @test back == s
        end
    end

    @testset "the authalic grid answers its parent's id arithmetic" begin
        g = levelgrid(SYS, LEVEL)
        a = levelgrid(AuthalicSystem(SYS), LEVEL)
        @test DGG.idtype(a) === DGG.idtype(g)
        for c in CELLS[1:7:end]
            x = DGG.rawid(c)
            @test DGG.idrank(a, x) == DGG.idrank(g, x)
            @test DGG.idselect(a, DGG.idrank(a, x)) == x
            @test DGG.idvalid(a, x)
            @test DGG.idcell(a, x) == c
            @test DGG.Encodings.idlevel(a, x) == LEVEL
        end
        @test DGG.idranges(a, CELLS; merge=:step) == DGG.idranges(g, CELLS; merge=:step)
    end

    ZARR_LOADED = try
        @eval using Zarr
        true
    catch
        false
    end

    if !ZARR_LOADED
        @info "Zarr.jl is unavailable: the placement round trips are skipped."
    else
        @testset "geodetic IGEO7 round-trips through a store" begin
            vals = Float32.(1:length(CELLS))
            for (sys, flag) in ((SYS, false), (AuthalicSystem(SYS), true))
                lookup = CellLookup(CellVector(sys, LEVEL, CELLS))
                stack = DD.DimStack((elevation=copy(vals),), (Cells(lookup),))
                path = joinpath(mktempdir(), "placement.zarr")
                DGG.dggwrite(path, stack)
                z = Zarr.zopen(path)
                @test z.attrs["dggs"]["dggs_vert0_lon"] == 11.2
                @test z.attrs["dggs"]["dggs_vert0_lat"] ≈ 58.282525588538995
                @test z["cell_id_ranges"].attrs["igeo7_wgs84_geodetic_conversion"] === flag
                back = dggread(path)
                grid = parent(DD.lookup(back, DD.dims(back)[1])).grid
                @test DGG.system(grid) == sys
                @test collect(DD.lookup(back, DD.dims(back)[1])) == CELLS
                @test Array(back[:elevation]) == vals
            end
        end

        @testset "a store at another placement is refused on read" begin
            lookup = CellLookup(CellVector(SYS, LEVEL, CELLS))
            stack = DD.DimStack((elevation=Float32.(1:length(CELLS)),), (Cells(lookup),))
            path = joinpath(mktempdir(), "isea7h.zarr")
            DGG.dggwrite(path, stack)
            # stamp the DGGRID ISEA7H default by hand, in both conventions, so
            # that they agree with each other and disagree with IGEO7
            for f in (joinpath(path, ".zattrs"), joinpath(path, "cell_id_ranges", ".zattrs"))
                s = read(f, String)
                t = replace(s, r"(\"(?:igeo7_)?dggs_vert0_lon\":\s*)11\.2\b" => s"\g<1>11.25")
                @test t != s
                write(f, t)
            end
            err = try
                dggread(path)
                nothing
            catch e
                e
            end
            @test err isa DGGSFormatError
            @test err.check === :orientation_mismatch
        end
    end
end

end # module DGGSIOPlacementTests
