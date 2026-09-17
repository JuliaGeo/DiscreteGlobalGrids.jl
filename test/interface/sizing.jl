# Cell size at a level, and the level a target resolution wants.
#
# Two systems carry the whole file. HEALPix is exactly equal-area, so its sizes
# are closed form and an area of interest cannot move them. CopernicusDEM's
# level-0 cells are 1° lon/lat boxes whose area falls off as cos(latitude), so
# it is the system where an area of interest has to move them.

module SizingTests

using Test
using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
import GeometryOpsCore as GOCore
import DimensionalData as DD
import Extents
import GlobalRegridding as GR

const HP = HEALPixSystem()
const COP = DGG.CopernicusDEMSystem(90)
const R = DGG.authalic_sphere(GOCore.Geodesic()).radius

_axis(D, centres, step) = D(DD.Sampled(collect(centres); span = DD.Regular(step),
    sampling = DD.Intervals(DD.Center()), order = DD.ForwardOrdered()))

function globalraster(step)
    lon = (-180 + step / 2):step:180
    lat = (-90 + step / 2):step:90
    return DD.DimArray([x + y for x in lon, y in lat],
        (_axis(DD.X, lon, step), _axis(DD.Y, lat, step)))
end

@testset "cellsize is the cell's own width, on the authalic sphere" begin
    # HEALPix level `l` is exactly 12·4^l cells of 4π/(12·4^l) steradians each.
    for l in (0, 3, 9)
        @test cellsize(HP, l) ≈ sqrt(4pi / (12 * 4^l)) * R
    end
    # A grid and its `(system, level)` spelling are the same grid.
    @test cellsize(levelgrid(HP, 4)) == cellsize(HP, 4)
    # `radius` scales the answer and nothing else.
    @test cellsize(HP, 4; radius = 1.0) * R ≈ cellsize(HP, 4)
end

@testset "levelfor takes the level bracketing the target" begin
    r = globalraster(1.0)
    # A 1° raster's median cell is the box at ±45°, ~93.5 km on a side. That
    # falls between HEALPix level 7 (~51 km) and level 6 (~102 km) and is the
    # nearer to level 6 in ratio; either neighbouring level is a factor of two
    # off.
    @test levelfor(HP, r) == 6
    @test cellsize(HP, 7) < cellsize(r) < cellsize(HP, 6)
    # A size in metres and a raster of that size choose the same level, and a
    # level's own size chooses that level back.
    @test levelfor(HP, cellsize(r)) == 6
    for l in (2, 6, 11)
        @test levelfor(HP, cellsize(HP, l)) == l
    end
end

const ARCTIC = Extents.Extent(X = (-10.0, 10.0), Y = (75.0, 85.0))
const EQUATOR = Extents.Extent(X = (-10.0, 10.0), Y = (-5.0, 5.0))

@testset "an area of interest moves the size where cell area varies" begin
    arctic, equator = ARCTIC, EQUATOR
    # 1° boxes shrink as cos(latitude): ~111 km at the equator and ~47 km at
    # 80° N, against a ~93.5 km global median (the box at ±45°).
    @test cellsize(COP, 0; over = arctic) < 0.6 * cellsize(COP, 0)
    @test cellsize(COP, 0; over = equator) > 1.1 * cellsize(COP, 0)
    # Equal-area cells cannot move, so the same calls on HEALPix return the
    # global number exactly: the area of interest selects cells, it does not
    # measure itself.
    @test cellsize(HP, 5; over = arctic) == cellsize(HP, 5)
    @test cellsize(HP, 5; over = equator) == cellsize(HP, 5)
end

# A 1° lon/lat box at latitude φ has area (π/180)^2·cos φ steradians, so the
# HEALPix level nearest a box centred on φ follows from cos φ alone.
boxlevel(φ) = argmin(l -> abs(log(4pi / (12 * 4^l)) - log((pi / 180)^2 * cosd(φ))), 0:29)

@testset "levelfor measures a spatial target within the area of interest" begin
    r = globalraster(1.0)
    grid = levelgrid(COP, 0)
    # The arctic boxes are ~47 km across and the equatorial ones ~111 km, one
    # HEALPix level apart, and each region takes the level its own boxes want.
    @test boxlevel(80) == 7 > boxlevel(0) == 6
    @test levelfor(HP, r; over = ARCTIC) == 7
    @test levelfor(HP, r; over = EQUATOR) == 6
    @test levelfor(HP, grid; over = ARCTIC) == 7
    @test levelfor(HP, grid; over = EQUATOR) == 6
    # A `RegridSpace` target — a raster space or a DGGSpace over the same
    # cells — agrees with the raster and grid spellings.
    for space in (GR.RasterGrid(r), DGG.DGGSpace(grid))
        @test levelfor(HP, space; over = ARCTIC) == 7
        @test levelfor(HP, space; over = EQUATOR) == 6
    end
    @test cellsize(r; over = ARCTIC) < cellsize(r; over = EQUATOR)
    @test cellsize(GR.RasterGrid(r); over = ARCTIC) == cellsize(r; over = ARCTIC)
    # A `DGGSpace` reads the grid it wraps, with an area of interest and without.
    @test cellsize(DGG.DGGSpace(grid); over = ARCTIC) == cellsize(grid; over = ARCTIC)
    @test cellsize(DGG.DGGSpace(grid)) == cellsize(grid)
end

@testset "a target in metres has one size everywhere" begin
    r = globalraster(1.0)
    for over in (ARCTIC, EQUATOR)
        @test levelfor(HP, 25_000; over) == levelfor(HP, 25_000) == 8
        @test levelfor(HP, cellsize(r); over) == levelfor(HP, cellsize(r)) == 6
    end
    # Without an area of interest the global median still decides.
    @test levelfor(HP, r) == levelfor(HP, GR.RasterGrid(r)) == 6
end

@testset "an area of interest the target misses is an error" begin
    step = 1.0
    lon = (-180 + step / 2):step:180
    lat = (step / 2):step:90
    north = DD.DimArray([x + y for x in lon, y in lat],
        (_axis(DD.X, lon, step), _axis(DD.Y, lat, step)))
    south = Extents.Extent(X = (-10.0, 10.0), Y = (-60.0, -50.0))
    @test_throws ArgumentError levelfor(HP, north; over = south)
    @test_throws "the target does not meet the area of interest" levelfor(
        HP, north; over = south)
    @test_throws "the target does not meet the area of interest" cellsize(
        north; over = south)
    @test levelfor(HP, north; over = ARCTIC) == 7
end

end # module SizingTests
