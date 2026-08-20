# Run from an environment containing the sibling Geomorphometry checkout:
#
# julia --project=<integration-env> test/integration/geomorphometry_igeo7.jl DEM.tif

using Test
using Statistics
import ArchGDAL
import ConservativeRegridding as CR
import DiscreteGlobalGrids as DGG
import Geomorphometry as GM
import GeometryOps as GO
using Rasters

isempty(ARGS) && error("pass the Copernicus DEM GeoTIFF path as the first argument")
dem_path = abspath(only(ARGS))
isfile(dem_path) || error("DEM does not exist: $dem_path")

source_dem = Raster(dem_path; lazy=true)
dem = aggregate(mean, source_dem, 60)

(west, east), (south, north) = bounds(dem, X), bounds(dem, Y)
to_sphere = GO.UnitSpherical.UnitSphereFromGeographic()
corners = [
    to_sphere((lon, lat))
    for lon in range(west, east; length=size(dem, X) + 1),
        lat in range(south, north; length=size(dem, Y) + 1)
]

sys = DGG.IGeo7System()
center = to_sphere(((west + east) / 2, (south + north) / 2))
root = DGG.cellat(DGG.levelgrid(sys, 5), center)
grid = DGG.subtree(sys, root, 8)

regridder = CR.Regridder(GO.Spherical(; radius=1.0), grid, corners)
source = vec(Float64.(reverse(parent(dem); dims=2)))
raw = zeros(DGG.ncells(grid))
cover = zeros(DGG.ncells(grid))
CR.regrid!(raw, regridder, source)
CR.regrid!(cover, regridder, ones(length(source)))

@test all(>(0), cover)
elevation = raw ./ cover
cells = DGG.CellVector(grid)
igeo7_dem = Raster(
    elevation,
    (DGG.Cells(DGG.CellLookup(cells)),);
    name=:height,
)

@testset "Geomorphometry on regridded IGeo7 DEM" begin
    indices = eachindex(igeo7_dem)
    @test indices === cells
    @test indices isa DGG.CellVector
    @test first(indices) in indices
    complete = DGG.levelgrid(sys, DGG.level(first(indices)))
    candidate = DGG.cellindex(complete, 1)
    candidate in indices &&
        (candidate = DGG.cellindex(complete, DGG.ncells(complete)))
    @test candidate ∉ indices

    tpi = GM.topographic_position_index(igeo7_dem)
    @test size(tpi) == size(igeo7_dem)
    @test all(isfinite, parent(tpi))

    cell_area = GM.cellarea(igeo7_dem, first(eachindex(igeo7_dem)))
    threshold = 5cell_area
    for method in (GM.D8(), GM.DInf(), GM.FD8())
        accumulation, directions = GM.flowaccumulation(igeo7_dem; method)
        hand = GM.height_above_nearest_drainage(
            igeo7_dem;
            method,
            threshold,
        )
        @test size(accumulation) == size(directions) == size(igeo7_dem)
        @test size(hand) == size(igeo7_dem)
        @test all(isfinite, parent(accumulation))
        @test all(isfinite, parent(hand))
        @test all(>=(0.99cell_area), parent(accumulation))
        @test all(>=(0), parent(hand))
    end
end

println((
    source_size=size(source_dem),
    aggregated_size=size(dem),
    igeo7_cells=length(igeo7_dem),
    elevation=extrema(elevation),
))
