# # Regridding: getting data onto a grid
#
# Use `regrid` to bring a longitude/latitude raster onto a global grid. The
# result keeps its time dimension and works with the package's spatial
# selectors, neighbourhood operations and plotting tools.
#
# Here we move twelve months of soil moisture onto HEALPix, choose how to
# handle missing ocean values, and map the result back to a raster. We use the
# default `Conservative()` method, which weights source values by their overlap
# with each destination cell. [Moving between DGGS](between_grids.md) compares
# area averaging with point interpolation.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
using Rasters, RasterDataSources
import NCDatasets
import Dates
using Statistics
using GLMakie, GeoMakie
using DiscreteGlobalGridsVisualization: dggpoly, dggpoly!
GLMakie.activate!(inline = true)

# ## The source: CPC soil moisture on a half-degree raster
#
# The data is a monthly climatology from the NOAA Climate Prediction Center:
# one soil moisture value per half-degree lon/lat cell per month.

soil = Rasters.set(Raster(RasterDataSources.getraster(CPCSoil; period = "1981-2010");
    name = :soilw), Ti => 1:12)
soil = DD.rebuild(soil; metadata = DD.NoMetadata())

# Every `missing` cell is ocean — the source measures land — and the ocean is
# most of the raster:

count(ismissing, soil) / length(soil)

# The January map shows the source on its `0 … 360` longitude axis. `regrid`
# handles that longitude convention directly; the ocean remains missing.

fig, ax, plt = heatmap(soil[Ti = 1]; colormap = :viridis,
    axis = (; title = "January", xlabel = "longitude", ylabel = "latitude"))
Colorbar(fig[1, 2], plt; label = "soil moisture (mm)")
fig

# ## Regrid the raster onto HEALPix
#
# HEALPix gives each cell the same area, which simplifies spatial averages.
# Choose a resolution close to the source with `levelfor`, then build that
# level with `levelgrid`. [Choosing a grid](choosing_a_grid.md) explains cell
# sizes and coordinate conventions.

sys = DGG.HEALPixSystem()
grid = DGG.levelgrid(sys, DGG.levelfor(sys, soil))

# Metres across a source cell, and across a destination cell:

DGG.cellsize(soil), DGG.cellsize(grid)

# One call regrids all twelve months:

onhealpix = DGG.regrid(soil; to = grid)

# The result is a `Raster` with dimensions `Cells` and `Ti`: one value per
# cell per month. You can still select January with `[Ti = 1]`.
#
# ## Draw the result on the sphere
#
# `replace_missing` turns the ocean cells into `NaN`, and `dggpoly!` leaves a
# `NaN` cell unpainted:

january = Rasters.replace_missing(onhealpix[Ti = 1], NaN)

fig = Figure(size = (820, 400))
ax = GeoAxis(fig[1, 1]; dest = "+proj=moll", title = "January")
plt = dggpoly!(ax, january; color = january, colormap = :viridis)
Colorbar(fig[1, 2], plt; label = "soil moisture (mm)")
fig

# Read January's value at a few locations with `DD.Contains((lon, lat))`.
# The mid-Pacific query also checks how missing ocean data appears:

[place => onhealpix[DGG.Cells(DD.Contains((lon, lat))), Ti = DD.At(1)]
 for (place, lon, lat) in (("Amazon", -60.0, -5.0), ("Sahara", 15.0, 24.0),
                           ("Kansas", -100.0, 38.0), ("mid-Pacific", -150.0, 0.0))]

# ## Choose what a coastal cell holds
#
# A coastal cell can overlap both valid land values and missing ocean values.
# The default `Weighted(0.5)` averages the valid contribution and requires at
# least half of the total weight to be valid. A cell below that threshold gets
# a missing value.
#
# Lowering the threshold keeps more coastal cells. Raising it asks for more
# complete coverage. Compare the number of valid January cells:

covered(t) = count(!ismissing, DGG.regrid(soil; to = grid,
    missingpolicy = DGG.Weighted(t))[:, 1])
[t => covered(t) for t in (0.01, 0.5, 1.0)]

# The difference between these counts shows how much of the result depends on
# partial coverage. This raster uses `missing` for those cells. The
# [regridding reference](../api/regridding-methods.md) explains other missing
# value conventions and the `Extensive()` policy.
#
# ## Reuse one plan across the twelve months
#
# Reuse a regridding plan when several variables share the same source and
# destination grids. The plan computes the spatial weights once:

plan = DGG.plan_regrid(soil; to = grid)

# Apply it to the full monthly cube:

monthly = DGG.regrid(soil, plan)

# Use `regrid!` to reuse an output array as well:

buffer = similar(monthly)
DGG.regrid!(buffer, soil, plan)
isequal(buffer, monthly)

# ## Animate the seasonal cycle
#
# A shared colour scale makes the months comparable. Animate the values on
# the same grid to see the seasonal cycle:

seasonal = Rasters.replace_missing(monthly, NaN)
crange = extrema(filter(!isnan, seasonal))
colors = Observable(seasonal[Ti = 1])
title = Observable(Dates.monthname(1))

fig = Figure(size = (820, 400))
ax = GeoAxis(fig[1, 1]; dest = "+proj=moll", title)
plt = dggpoly!(ax, seasonal[Ti = 1]; color = colors,
    colormap = :viridis, colorrange = crange)
Colorbar(fig[1, 2], plt; label = "soil moisture (mm)")

record(fig, "seasonal_soilw.mp4", 1:12; framerate = 2) do m
    colors[] = seasonal[Ti = m]
    title[] = Dates.monthname(m)
end
nothing #hide

#md # ```@raw html
#md # <video autoplay loop muted playsinline controls src="./seasonal_soilw.mp4" style="max-width: 100%;"/>
#md # ```

# ## Regrid back onto a lon/lat raster
#
# Use an existing raster as the destination template when a downstream tool
# expects longitude and latitude axes. Pass the source DGGS grid as `from`:

back = DGG.regrid(monthly; to = soil, from = grid)

# Compare the original and returned values where both are valid. These plain
# means are a useful check of the change in the field; they do not test area
# conservation, because longitude/latitude pixels have unequal areas. The
# round trip also retains the smoothing introduced by regridding.

both = .!ismissing.(soil) .&& .!ismissing.(back)
mean(soil[both]), mean(back[both])

# You can also choose new destination coordinates. A `-180 … 180` longitude
# axis fits the following world map:

shown = DGG.regrid(monthly[Ti = 1]; from = grid,
    to = Rasters.set(soil, X => -179.75:0.5:179.75))

fig = Figure(size = (820, 400))
ax = GeoAxis(fig[1, 1]; dest = "+proj=moll", title = "January, regridded twice")
plt = surface!(ax, Rasters.replace_missing(shown, NaN);
    colormap = :viridis, colorrange = crange, shading = NoShading)
Colorbar(fig[1, 2], plt; label = "soil moisture (mm)")
fig

# ## Regrid onto any other system
#
# Choose another system and match its cell size to the same source:

igeo7 = DGG.IGeo7System()
DGG.regrid(soil; to = DGG.levelgrid(igeo7, DGG.levelfor(igeo7, soil)))

# The result supports the same operations whichever system you choose. Try
# [Zonal statistics](zonal.md) to summarize it by region, or
# [Stencil operations](stencils.md) to compute from neighbouring cells.
#
# For a different interpretation of the source values, see
# [Choosing a regridding method](../api/regridding-methods.md).
