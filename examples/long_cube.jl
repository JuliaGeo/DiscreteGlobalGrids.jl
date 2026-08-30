#=
# Operations on long cubes

Here, we'll show operations on a long datacube (long in time) - in this case,
an ERA5 dataset in my downloads folder.
=#

using Rasters, RasterDataSources
import NCDatasets, ArchGDAL
using Dates

using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG

using GLMakie, GeoMakie, DiscreteGlobalGridsVisualization
GLMakie.activate!(inline=true)

# First, let's load our raster:
ERA5_FILE = joinpath(homedir(), "Downloads", "t2m.hh.an.era5.06.2000.nc")
if isfile(ERA5_FILE)
    ras = Raster(ERA5_FILE; lazy = true)
else
    ras = Rasters.combine(RasterSeries(WorldClim{Climate}, :tavg; month = 1:12))
end
# Then, we can get a global IGeo7 grid of this.  First, we establish the closest resolution:
sys = DGG.IGeo7System()
level = DGG.levelfor(sys, ras)

# Then, we construct a grid from the IGeo7 system and the level:
grid = DGG.levelgrid(sys, level)

# Now, we can create a regridding plan:
plan = @time DGG.plan_regrid(ras; to = grid, method = DGG.Conservative())
# You can also use a faster regridding method like `DGG.Bilinear()` or `DGG.NearestNeighbor()`, 
# but for this example, we'll use the conservative method.  It is a bit slower when planning,
# but the time to execute the plan is about the same between the three methods.

# Let's suppose we get a single slice of this raster and implement the plan:
ras_slice = Rasters.read(ras[:, :, 1])
#
igeo7_ras = @time Raster(DGG.regrid(ras_slice, plan))
# You can also regrid into an existing buffer:
@time DGG.regrid!(igeo7_ras, ras[:, :, 2], plan)

# Now, let's plot these side by side:
fig = Figure()
ras_ax, ras_plt = heatmap(fig[1, 1], ras_slice; axis = (; title = "ERA5 Raster Slice", aspect = DataAspect()))
igeo7_ax, igeo7_plt = DiscreteGlobalGridsVisualization.dggpoly(fig[2, 1], lookup(igeo7_ras, DGG.Cells); color = vec(Rasters.replace_missing(igeo7_ras; missingval = NaN)), axis = (; title = "IGeo7 Regridded Slice", aspect = DataAspect()))
GLMakie.tightlimits!(igeo7_ax)
fig
