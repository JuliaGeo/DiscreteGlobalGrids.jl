# # Moving between DGGS
#
# Move data between grid systems when you need to combine datasets or use a
# different cell shape. `regrid` handles both a change of system and a change
# of resolution.
#
# This example moves monthly soil moisture from HEALPix to IGeo7. We compare
# area averaging and point interpolation, first on a much finer destination
# and then on grids of similar cell size. The setup repeats the data loading
# from [Regridding](regridding.md), so you can run this page on its own.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
using Rasters, RasterDataSources
import NCDatasets
import Extents
using Statistics
using GLMakie, GeoMakie
using DiscreteGlobalGridsVisualization: dggpoly, dggpoly!
GLMakie.activate!(inline = true)

soil = Rasters.set(Raster(RasterDataSources.getraster(CPCSoil; period = "1981-2010");
    name = :soilw), Ti => 1:12)
soil = DD.rebuild(soil; metadata = DD.NoMetadata())
#
healpix = DGG.levelgrid(DGG.HEALPixSystem(), 7)
soilonhealpix = DGG.regrid(soil; to = healpix)

# ## Move the cube from HEALPix onto IGeo7
#
# `to` names the destination grid and `from` names the source grid:

igeo7 = DGG.levelgrid(DGG.IGeo7System(), 5)
crossed = DGG.regrid(soilonhealpix; to = igeo7, from = healpix)

# The result has one value per IGeo7 cell for each of the twelve months.
# Plot January to see the soil moisture field on the new cells; missing ocean
# values remain blank.

january = Rasters.replace_missing(crossed[Ti = 1], NaN)

fig = Figure(size = (820, 400))
ax = GeoAxis(fig[1, 1]; dest = "+proj=natearth2", xticks = -120:60:120,
    yticks = -60:30:60, title = "January on IGeo7 level 5")
plt = dggpoly!(ax, january; color = january, colormap = :viridis)
Colorbar(fig[1, 2], plt; label = "soil moisture (mm)")
fig

# ## Check the cell size and the land mean
#
# Metres across a cell, on either grid:

DGG.cellsize(healpix), DGG.cellsize(igeo7)

# These levels have similar cell widths. Comparing their plain means gives
# a quick check of the effect on the data. For an area-weighted comparison,
# use `cell_area` weights on IGeo7 and account for changes in coastal coverage;
# its cells are only approximately equal in area.

mean(skipmissing(soilonhealpix)), mean(skipmissing(crossed))

# ## Refine a coarse field onto a finer grid
#
# A finer destination makes the choice of method visible. Choose according
# to what a source value represents: an average over a cell, or a measurement
# at a point.
#
# | method | what a destination cell gets |
# |---|---|
# | `Conservative()` (default) | the area-weighted mean of the source cells under it |
# | `BarycentricPoint()` | a sample at its centre, interpolated from surrounding source centres |
#
# HEALPix level 4 cells are 407 km across, against 55 km on IGeo7 level 5, so
# roughly fifty destination cells sit under each source cell:

coarse = DGG.levelgrid(DGG.HEALPixSystem(), 4)
oncoarse = DGG.regrid(soil[Ti = 1]; to = coarse)

# Both methods run across that pair:

blocky = DGG.regrid(oncoarse; to = igeo7, from = coarse)
smooth = DGG.regrid(oncoarse; to = igeo7, from = coarse,
    method = DGG.BarycentricPoint())

# Index by an extent to draw Europe alone:

europe = Extents.Extent(X = (-12.0, 42.0), Y = (35.0, 65.0))
crange = extrema(filter(!isnan, january))

fig = Figure(size = (900, 430))
plots = map(enumerate(("Conservative()" => blocky,
                       "BarycentricPoint()" => smooth))) do (k, (name, field))
    ax = GeoAxis(fig[1, k]; dest = "+proj=laea +lat_0=52 +lon_0=15",
        xticks = -10:10:40, yticks = 35:10:65, title = name)
    sub = Rasters.replace_missing(field, NaN)[DGG.Cells(DGG.Covering(europe))]
    dggpoly!(ax, sub; color = sub, colormap = :viridis, colorrange = crange)
end
Colorbar(fig[1, 3], first(plots); label = "soil moisture (mm)")
fig

# `Conservative()` repeats a source value in destination cells wholly inside
# that source cell, and combines values where cells overlap. The coarse
# HEALPix pattern remains visible. `BarycentricPoint()` treats the values as
# samples at the source centres and interpolates a smoother surface. Neither
# method adds measurements at the finer resolution.
#
# ## Compare the two methods at equal resolution
#
# HEALPix level 7 and IGeo7 level 5 hold cells of nearly the same size, so each
# destination cell draws on one source cell and its immediate neighbours. Run
# `BarycentricPoint()` across that pair:

pointwise = Rasters.replace_missing(
    DGG.regrid(soilonhealpix[Ti = 1]; to = igeo7, from = healpix,
        method = DGG.BarycentricPoint()), NaN)

# The 99th percentile of the absolute difference, beside the standard deviation
# of the field itself, in millimetres:

q99 = quantile(abs.(filter(!isnan, january .- pointwise)), 0.99)
(q99 = q99, sigma = std(filter(!isnan, january)))

# Compare `q99` with `sigma` to judge the method difference relative to the
# field's own variation. The third panel shows where the methods differ:

conservative = january[DGG.Cells(DGG.Covering(europe))]
barycentric = pointwise[DGG.Cells(DGG.Covering(europe))]

fig = Figure(size = (1240, 430))
plots = map(enumerate(("Conservative()" => conservative,
                       "BarycentricPoint()" => barycentric))) do (k, (name, field))
    ax = GeoAxis(fig[1, k]; dest = "+proj=laea +lat_0=52 +lon_0=15",
        xticks = -10:10:40, yticks = 35:10:65, title = name)
    dggpoly!(ax, field; color = field, colormap = :viridis, colorrange = crange)
end
Colorbar(fig[1, 3], first(plots); label = "soil moisture (mm)")
ax = GeoAxis(fig[1, 4]; dest = "+proj=laea +lat_0=52 +lon_0=15",
    xticks = -10:10:40, yticks = 35:10:65, title = "difference")
plt = dggpoly!(ax, conservative; color = conservative .- barycentric,
    colormap = :balance, colorrange = (-q99, q99))
Colorbar(fig[1, 5], plt; label = "Conservative() − BarycentricPoint() (mm)")
fig

# Coastal cells are particularly sensitive to the method: area averaging uses
# overlap weights, while point interpolation uses surrounding sample sites.
# Missing source values can therefore affect different destination cells.
# Compare how many cells each method leaves missing:

count(isnan, pointwise) - count(isnan, january)

# ## Reuse one plan across the twelve months
#
# `plan_regrid` builds the weights for the pair of grids once. Every field that
# crosses the same pair — the twelve months here — reuses them:

plan = DGG.plan_regrid(soilonhealpix; to = igeo7, from = healpix)

# Reuse the output buffer while computing a monthly mean. Each call replaces
# the previous month's values:

dest = DGG.regrid(soilonhealpix[Ti = 1], plan)
seasonal = map(1:12) do m
    DGG.regrid!(dest, soilonhealpix[Ti = m], plan)
    mean(skipmissing(dest))
end

# Plot the monthly cell means to see the seasonal cycle:

fig = Figure(size = (600, 340))
ax = Axis(fig[1, 1]; xticks = 1:12, xlabel = "month",
    ylabel = "mean soil moisture (mm)", title = "global land mean, IGeo7 level 5")
scatterlines!(ax, 1:12, seasonal)
fig

# ## Let the destination match the source resolution
#
# Pass a system as `to` to choose the destination level automatically. Here it
# selects the IGeo7 level closest in cell size to HEALPix level 7:

DGG.regrid(soilonhealpix[Ti = 1]; to = DGG.IGeo7System(), from = healpix)

# For a regional destination, pass a cell collection or coverage as `to`;
# [Multi-order coverage](multiorder.md) shows that workflow.
# [Choosing a regridding method](../api/regridding-methods.md) provides the
# method and missing-data reference.
