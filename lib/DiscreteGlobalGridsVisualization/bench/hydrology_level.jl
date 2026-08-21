# Run the hydrology tutorial's pipeline at a chosen IGEO7 level, stage by stage,
# reporting wall time and peak resident memory after each one.  This is the
# question the recipe exists to answer: at which level does the tutorial stop
# fitting on a CI runner, and is plotting still what decides it?
#
#     xvfb-run -a julia --project=docs -t8 \
#         lib/DiscreteGlobalGridsVisualization/bench/hydrology_level.jl 13
#
# Pass `poly` as a second argument to draw with Makie's `poly` instead of
# `dggpoly`, for the before-and-after, or `resample` to draw through
# `dggresample`, which shows the level the screen can carry rather than the
# level the data is stored at.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH",
    joinpath(tempdir(), "rasterdatasources")))

import DiscreteGlobalGrids as DGG
import Geomorphometry as GM
import GeoInterface as GI
using DiscreteGlobalGridsVisualization
using Rasters, RasterDataSources
import ArchGDAL
using Statistics
using GLMakie, GeoMakie
using Printf

GLMakie.activate!()

const LEVEL = isempty(ARGS) ? 12 : parse(Int, ARGS[1])
const DRAW = any(==("poly"), ARGS) ? :poly :
    any(==("resample"), ARGS) ? :dggresample : :dggpoly
const OUT = mktempdir()

stage_start = time()

"Report what a stage cost and how much memory the process has ever held."
function stage(name)
    elapsed = time() - stage_start
    @printf("%-34s %8.2f s   peak RSS %6.2f GiB\n", name, elapsed, Sys.maxrss() / 2^30)
    global stage_start = time()
    return nothing
end

function draw!(axis, cells; kwargs...)
    DRAW === :poly && return poly!(axis, cells; kwargs...)
    DRAW === :dggresample && return dggresample!(axis, cells; kwargs...)
    return dggpoly!(axis, cells; kwargs...)
end

# `dggresample` draws no outlines at all, so it has no `strokewidth` to set.
const stroke = DRAW === :dggresample ? (;) : (; strokewidth = 0)

println("level $LEVEL, drawing with $DRAW, $(Threads.nthreads()) threads")

path = only(skipmissing(RasterDataSources.getraster(CopernicusDEM; extent = GI.extent((10.5, 46.5)))))
dem = Raster(path; lazy = false)
stage("read the DEM tile")

sys = DGG.IGeo7System()
region = DGG.query(sys, DGG.MultiOrderCoverage(Rasters.extent(dem)); level = LEVEL)
cells_in_region = length(DGG.CellVector(region))
stage("cover the tile ($cells_in_region cells)")

igeo7_dem = DGG.regrid(dem; to = region)
stage("regrid the DEM onto the cells")

covered = .!isnan.(igeo7_dem)
grid = DGG.PartialGrid(lookup(igeo7_dem, DGG.Cells))
shown = findall(covered)
cells = DGG.CellVector(grid)[shown]
stage("select the covered cells ($(length(shown)))")

figure = Figure(size = (900, 430))
axis = GeoAxis(figure[1, 1]; dest = "+proj=longlat +datum=WGS84", title = "elevation (m)")
plot = draw!(axis, cells; color = vec(igeo7_dem[shown]), colormap = :terrain, stroke...)
stage("build the elevation figure")

save(joinpath(OUT, "elevation.png"), figure)
stage("save the elevation figure")

nbrs = [filter(p -> covered[p], row) for row in DGG.adjacency(grid)]
stage("adjacency and the covered filter")

function downhill(i)
    isempty(nbrs[i]) && return 0
    j = nbrs[i][argmin(igeo7_dem[nbrs[i]])]
    return igeo7_dem[j] < igeo7_dem[i] ? j : 0
end

flow = [covered[i] ? downhill(i) : 0 for i in 1:DGG.ncells(grid)]
drop = [flow[i] == 0 ? NaN : igeo7_dem[i] - igeo7_dem[flow[i]] for i in eachindex(flow)]
stage("flow direction")

terrain = Raster(igeo7_dem; name = :height)
tpi = GM.topographic_position_index(terrain)
stage("topographic position index")

accumulation, directions = GM.flowaccumulation(terrain; method = GM.D8())
cell_area = GM.cellarea(terrain, first(eachindex(terrain)))
log_cells = log10.(accumulation ./ cell_area)
stage("D8 flow accumulation")

figure = Figure(size = (900, 430))
axis = GeoAxis(figure[1, 1]; dest = "+proj=longlat +datum=WGS84", title = "topographic position index (m)")
draw!(axis, cells; color = vec(tpi[shown]), colorrange = (-25, 25), colormap = :delta, stroke...)
axis = GeoAxis(figure[1, 2]; dest = "+proj=longlat +datum=WGS84", title = "D8 flow accumulation")
draw!(axis, cells; color = vec(log_cells[shown]), colormap = :devon, stroke...)
save(joinpath(OUT, "terrain.png"), figure)
stage("build and save the terrain figures")

@printf("TOTAL peak RSS %.2f GiB\n", Sys.maxrss() / 2^30)
