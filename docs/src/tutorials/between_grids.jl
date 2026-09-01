# # Moving between DGGS
#
# `regrid` also moves data between DGGS grids. This page starts from the CPC
# soil moisture cube that [Regridding](regridding.md) put on HEALPix, and moves
# it onto IGeo7.

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
soil = DD.rebuild(soil; metadata = DD.NoMetadata())   # the NetCDF header is noise here

healpix = DGG.levelgrid(DGG.HEALPixSystem(), 7)
onhealpix = DGG.regrid(soil; to = healpix)

# ## Move the cube from HEALPix onto IGeo7
#
# `to` names the destination grid and `from` names the source grid:

igeo7 = DGG.levelgrid(DGG.IGeo7System(), 5)
crossed = DGG.regrid(onhealpix; to = igeo7, from = healpix)

# The `Cells` axis now carries IGeo7 cells, and `Ti` came through untouched. A
# `Cells` lookup carries its own grid, and `regrid` still asks for it by name:
# leave `from` out and the error prints the grid to pass.
#
# `replace_missing` turns the ocean cells into `NaN`, which `dggpoly!` leaves
# unpainted:

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

# The two agree to 8 %, so the move changes the tiling and leaves the
# resolution where it was. Both grids are equal-area, so the unweighted land
# mean is the same on either:

mean(skipmissing(onhealpix)), mean(skipmissing(crossed))

# ## Refine a coarse field onto a finer grid
#
# A move onto a much finer grid is a refinement, and the method decides how one
# coarse cell is spread over the finer cells under it:
#
# | method | what a destination cell gets |
# |---|---|
# | `Conservative()` (default) | the area-weighted mean of the source cells under it |
# | `BarycentricPoint()` | a sample at its centre, interpolated between the three source centres around it |
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

# `Conservative()` holds each coarse value flat across the fine cells under it,
# so the blocks on the left are the HEALPix cells the numbers came from.
# `BarycentricPoint()` interpolates between the coarse centroids, and reads the
# values as samples of a continuous surface, which suits data measured at
# points.
#
# ## Compare the two methods at equal resolution
#
# HEALPix level 7 and IGeo7 level 5 hold cells of nearly the same size, so each
# destination cell draws on one source cell and its immediate neighbours. Run
# `BarycentricPoint()` across that pair:

pointwise = Rasters.replace_missing(
    DGG.regrid(onhealpix[Ti = 1]; to = igeo7, from = healpix,
        method = DGG.BarycentricPoint()), NaN)

# The 99th percentile of the absolute difference, beside the standard deviation
# of the field itself, in millimetres:

q99 = quantile(abs.(filter(!isnan, january .- pointwise)), 0.99)
(q99 = q99, sigma = std(filter(!isnan, january)))

# At matched resolution the two methods agree to a few millimetres, on a field
# whose own spread is 177 mm. The third panel maps what is left:

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

# What is left sits on the coastlines: an area mean weights the land fraction
# under a cell, and a point sample reads the source centres around the cell's
# centre.
#
# `BarycentricPoint()` also blanks a destination cell whose centre lies outside
# every triangle of source centres, which along a coastline costs a handful of
# cells:

count(isnan, pointwise) - count(isnan, january)

# ## Reuse one plan across the twelve months
#
# `plan_regrid` builds the weights for the pair of grids once. Every field that
# crosses the same pair — the twelve months here — reuses them:

plan = DGG.plan_regrid(onhealpix; to = igeo7, from = healpix)

# `regrid!` writes the result into a buffer you own, and a blanked cell takes
# the buffer's `missingval`:

dest = DGG.regrid(onhealpix[Ti = 1], plan)
seasonal = map(1:12) do m
    DGG.regrid!(dest, onhealpix[Ti = m], plan)
    mean(skipmissing(dest))
end

# Twelve applies of one plan:

fig = Figure(size = (600, 340))
ax = Axis(fig[1, 1]; xticks = 1:12, xlabel = "month",
    ylabel = "mean soil moisture (mm)", title = "global land mean, IGeo7 level 5")
scatterlines!(ax, 1:12, seasonal)
fig

# ## What `to` and `from` accept
#
# | spelling | names |
# |---|---|
# | a grid | itself |
# | a `Cells` lookup, a `CellVector`, a `MultiOrderCellSet` | the partial grid of their cells |
# | a bare system, as `to` | the level whose cells come closest to the source's in size |
# | a raster, or a tuple of dimensions | its lon/lat lattice |
#
# Any of those sits on either end, so a move between two subsets of two
# different systems is the same call. The bare-system spelling picks IGeo7
# level 5 here, to match HEALPix level 7:

DGG.regrid(onhealpix[Ti = 1]; to = DGG.IGeo7System(), from = healpix)

# [Choosing a regridding method](../api/regridding-methods.md) compares the
# methods side by side.
