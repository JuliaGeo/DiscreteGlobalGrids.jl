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

@testset "an area of interest moves the size where cell area varies" begin
    arctic = Extents.Extent(X = (-10.0, 10.0), Y = (75.0, 85.0))
    equator = Extents.Extent(X = (-10.0, 10.0), Y = (-5.0, 5.0))
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

end # module SizingTests
