using Test
using CopernicusUtils
import DiscreteGlobalGrids as DGG
import GeometryOps as GO
import DiskArrays

const CD = DGG.CopernicusDEM
const sys = DGG.CopernicusDEMSystem(90)

@testset "tile stems" begin
    for stem in ("Copernicus_DSM_COG_30_N00_00_E006_00_DEM",
                 "Copernicus_DSM_COG_30_S90_00_W180_00_DEM")
        @test tilestem(sys, stemtile(sys, stem)) == stem
    end
    @test_throws ArgumentError stemtile(sys, "not_a_stem")
    @test startswith(tilestem(DGG.CopernicusDEMSystem(30), CD.tilecell(sys, 0, 6)),
        "Copernicus_DSM_COG_10_")
end

@testset "tile list" begin
    mktempdir() do dir
        path = joinpath(dir, "list.txt")
        write(path, """
            Copernicus_DSM_COG_30_N46_00_E010_00_DEM

            Copernicus_DSM_COG_30_N00_00_E006_00_DEM
            Copernicus_DSM_COG_30_N00_00_E006_00_DEM
            """)
        all = landtiles(sys, path)
        @test length(all) == 2 && issorted(all)
        @test landtiles(sys, path, [(5.0, 7.0, -1.0, 1.0)]) ==
              [Int(CD.tilecell(sys, 0, 6).index)]
    end
end

@testset "real tiles: outside landtiles is ocean, missing is an error offline" begin
    tiles = [Int(CD.tilecell(sys, 0, 6).index)]
    ocean = Int(CD.tilecell(sys, 0, -30).index)
    mktempdir() do dir
        src = CopernicusTiles(sys, tiles; cachedir = dir, download = false)
        @test tilepath!(src, ocean) === nothing
        @test all(isnan, loadtile(src, ocean))
        @test_throws ErrorException loadtile(src, tiles[1])
        @test src.ndownloads[] == 0
        @test endswith(tileurl(src, tiles[1]),
            "Copernicus_DSM_COG_30_N00_00_E006_00_DEM.tif")
    end
end

@testset "tile locators" begin
    lat, lon = 0, 6
    ordinal = Int(CD.tilecell(sys, lat, lon).index)
    stem = tilestem(90, lat, lon)
    @test stem == tilestem(sys, CD.tilecell(sys, lat, lon))
    mktempdir() do dir
        locate = TileDirectory(dir, 90)
        flat = joinpath(dir, stem * ".tif")
        nested = joinpath(dir, stem, stem * ".tif")
        @test locate(lat, lon) == flat
        @test isempty(landtiles(sys, locate))
        mkpath(dirname(nested)); touch(nested)
        @test locate(lat, lon) == nested
        @test landtiles(sys, locate) == [ordinal]
        @test isempty(landtiles(sys, locate, [(20.0, 30.0, 20.0, 30.0)]))
        touch(flat)
        @test locate(lat, lon) == flat

        custom(lat, lon) = lat < 0 ? nothing : joinpath(dir, "n$(lat)e$(lon).tif")
        south = Int(CD.tilecell(sys, -1, lon).index)
        src = CopernicusTiles(sys, [ordinal, south]; locate = custom, download = false)
        @test tilecachepath(src, ordinal) == joinpath(dir, "n0e6.tif")
        @test tilepath!(src, south) === nothing
        @test all(isnan, loadtile(src, south))
        @test_throws ErrorException tilepath!(src, ordinal)
        touch(joinpath(dir, "n0e6.tif"))
        @test tilepath!(src, ordinal) == joinpath(dir, "n0e6.tif")
        @test_throws ArgumentError CopernicusTiles(sys, [ordinal])
        @test_throws ArgumentError CopernicusTiles(sys, [ordinal]; cachedir = dir, locate = custom)
    end
end

@testset "synthetic tiles" begin
    tile = CD.tilecell(sys, 46, 10)
    vals = synthetic_tile(sys, tile, NOMASK)
    nc = Int(CD.ncols_at(sys, 46))
    @test length(vals) == nc * 1200 && all(isfinite, vals)
    # Post (row j, column i), both zero-based, sits at (10 + i/nc, 47 - j/1200).
    @test vals[1] ≈ synthetic_elevation(10.0, 47.0)
    @test vals[3 * nc + 8] ≈ synthetic_elevation(10 + 7 / nc, 47 - 3 / 1200)

    # One half-degree lattice cell of land, so the tile is part ocean.
    bits = falses(720, 360)
    bits[381, 87] = true                       # 10..10.5 E, 46.5..47 N
    mask = LandMask(720, 360, 0.5, 0.5, bits)
    @test island(mask, 10.25, 46.75) && !island(mask, 10.75, 46.75)
    masked = loadtile(SyntheticTiles(sys, mask), Int(tile.index))
    @test any(isfinite, masked) && any(isnan, masked)
    @test landmask("unused.shp", 0) === NOMASK
end

@testset "SubtreeIds and TiledDEM" begin
    tiles = sort!([Int(CD.tilecell(sys, lat, 10).index) for lat in (46, 60, 85)])
    ids = TileIds(sys, tiles)
    widths = [length(DGG.descendant_range(sys, DGG.LevelIndex(0, t), 1)) for t in tiles]
    @test length(ids) == sum(widths)
    @test tileat(ids, 1) == (1, 1)
    @test tileat(ids, widths[1]) == (1, widths[1])
    @test tileat(ids, widths[1] + 1) == (2, 1)
    @test issorted([ids[1], ids[widths[1]], ids[widths[1] + 1], ids[end]])

    dem = TiledDEM(SyntheticTiles(sys), tiles; slots = 2, stripes = 1)
    @test size(dem) == (sum(widths),)
    @test [length(c[1]) for c in DiskArrays.eachchunk(dem)] == widths
    # A read spanning the first tile boundary stitches both tiles.
    r = (widths[1] - 2):(widths[1] + 3)
    a = loadtile(SyntheticTiles(sys), tiles[1])
    b = loadtile(SyntheticTiles(sys), tiles[2])
    @test dem[r] == vcat(a[(end - 2):end], b[1:3])

    grid = DGG.PartialGrid(sys, 1, dem.ids)
    @test DGG.ncells(grid) == length(dem)
end

@testset "StripedLRUCache" begin
    calls = Threads.Atomic{Int}(0)
    cache = StripedLRUCache{Vector{Float32}}(
        k -> (Threads.atomic_add!(calls, 1); fill(Float32(k), 4)); slots = 2, stripes = 1)
    @test cache(1) == fill(1.0f0, 4)
    cache(1); cache(2)
    @test calls[] == 2
    cache(3)                                   # evicts 1, the least recently used
    cache(1)
    @test calls[] == 4
    s = CopernicusUtils.cachestats(cache)
    @test s.loads == 4 && s.hits == 1 && s.live == 2
end

@testset "covering_chunks" begin
    sys7 = DGG.IGeo7System()
    tiles = [Int(CD.tilecell(sys, 46, 10).index)]
    chunks = covering_chunks(sys7, sys, tiles, 5)
    @test !isempty(chunks) && issorted(chunks)
    g5 = DGG.levelgrid(sys7, 5)
    centre = DGG.cellat(g5, GO.UnitSpherical.UnitSphereFromGeographic()((10.5, 46.5)))
    @test DGG.localindex(g5, centre) in chunks
    @test covering_chunks(sys7, sys, tiles, 5; nthreads = 2) == chunks
end
