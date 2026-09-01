# # Regridding: getting data onto a grid
#
# Regridding is how you get data that is not on a DGGS onto a DGGS.
#
# Three methods are available today:
#
# - `Conservative()` — the area-weighted mean of the source values over the
#   ground each destination cell covers. It conserves the integral, and it is
#   the default.
# - `BarycentricPoint()` — interpolation between the source sample sites around
#   each destination cell's centre: bilinear on a raster's quadrilaterals,
#   barycentric on a triangle.
# - `NearestCell()` — the value of the one source cell that contains each
#   destination cell's centre.
#
# More will come. [Choosing a regridding method](../api/regridding-methods.md)
# compares the three side by side.

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
soil = DD.rebuild(soil; metadata = DD.NoMetadata())   # the NetCDF header is noise here

# Every `missing` cell is ocean — the source measures land — and the ocean is
# most of the raster:

count(ismissing, soil) / length(soil)

# Two properties of this lattice the regridder reads for itself:
#
# - longitude runs `0 … 360`, and a longitude axis is periodic, so the cut
#   falls where the data puts it;
# - the lookups are `Points`, so cell edges are the midpoints between centres.
#
# Here is January. `heatmap` on a plain `Axis` draws the raster as it is
# stored, and leaves the ocean blank:

fig, ax, plt = heatmap(soil[Ti = 1]; colormap = :viridis,
    axis = (; title = "January", xlabel = "longitude", ylabel = "latitude"))
Colorbar(fig[1, 2], plt; label = "soil moisture (mm)")
fig

# ## Regrid the raster onto HEALPix
#
# HEALPix cells are exactly equal-area, so the unweighted mean of a field over
# them is its areal mean. `levelfor` picks the level whose cells come closest
# in size to the source's, and `levelgrid` builds the grid at that level:

sys = DGG.HEALPixSystem()
grid = DGG.levelgrid(sys, DGG.levelfor(sys, soil))

# Metres across a source cell, and across a destination cell:

DGG.cellsize(soil), DGG.cellsize(grid)

# `to` names the destination:

onhealpix = DGG.regrid(soil; to = grid)

# A raster in is a raster out: the two spatial dimensions have become one
# `Cells` dimension carrying the grid, and `Ti` passed through untouched.
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

# The values land where they should. A `Cells` lookup takes a lon/lat position
# through `DD.Contains`:

[place => onhealpix[DGG.Cells(DD.Contains((lon, lat))), Ti = DD.At(1)]
 for (place, lon, lat) in (("Amazon", -60.0, -5.0), ("Sahara", 15.0, 24.0),
                           ("Kansas", -100.0, 38.0), ("mid-Pacific", -150.0, 0.0))]

# ## Choose what a coastal cell holds
#
# A destination cell on a coastline is part land, part ocean. `missingpolicy`
# decides its value:
#
# | policy | value | blanked when |
# |---|---|---|
# | `Weighted(t)`, default `t = 0.5` | mean over the covered part | coverage below `t` |
# | `Extensive()` | sum over the covered part | never |
#
# `Weighted` normalizes by the covered area, so a coastal cell keeps the
# magnitude of the land it saw. A lower `t` admits more of the coastline:

covered(t) = count(!ismissing, DGG.regrid(soil; to = grid,
    missingpolicy = DGG.Weighted(t))[:, 1])
[t => covered(t) for t in (0.01, 0.5, 1.0)]

# The gap between the first count and the last is the coastline: cells that saw
# some land and some sea.
#
# A blanked cell holds the destination's own sentinel:
#
# - a `Raster` keeps its `missingval`, which is `missing` here because that is
#   what `soil` declares;
# - a plain array takes `missing` or `NaN`, by element type;
# - `missingval = NaN` names the sentinel explicitly;
# - `regrid!` uses the buffer's.
#
# ## Reuse one plan across the twelve months
#
# `plan_regrid` builds the weights from the two lattices alone, reading no
# array values:

plan = DGG.plan_regrid(soil; to = grid)

# The plan covers the spatial dimensions only, so `Ti` passes through it: one
# plan regrids all twelve months, and every other field on the same lattice.

monthly = DGG.regrid(soil, plan)

# `regrid!` applies the plan into a buffer you already hold, blanking with the
# buffer's own sentinel:

buffer = similar(monthly)
DGG.regrid!(buffer, soil, plan)
isequal(buffer, monthly)

# ## Animate the seasonal cycle
#
# Only the colours change from frame to frame, so they go in an `Observable`
# and the mesh is built once:

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
# The same call runs the other way:
#
# - `to` takes a raster, and its lookups become the destination lattice;
# - `from` names the source grid, because a `Cells` axis carries cell ids
#   rather than a lattice.

back = DGG.regrid(monthly; to = soil, from = grid)

# A round trip conserves the mean over the cells covered in both directions.
# Detail that the coarser grid averaged away stays averaged.

both = .!ismissing.(soil) .&& .!ismissing.(back)
mean(soil[both]), mean(back[both])

# Any lattice is a valid destination, and a `-180 … 180` one draws cleanly on a
# map cut at ±180:

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
# The same call takes any system:

igeo7 = DGG.IGeo7System()
DGG.regrid(soil; to = DGG.levelgrid(igeo7, DGG.levelfor(igeo7, soil)))

# Two things depend on the system:
#
# - an areal mean over cells that differ in area needs `cell_area` weights;
# - conservation into a destination holds for systems whose cell rings are
#   convex, IGeo7 among them. On HEALPix, `missingpolicy = DGG.Weighted(t)`
#   normalizes by covered area and carries the result.
#
# Every regrid above is conservative, an area operation. For data whose values
# sit at points rather than over footprints — a DEM's posts, a station network —
# see [Choosing a regridding method](../api/regridding-methods.md), and
# [Moving between DGGS](between_grids.md) for the DGGS-to-DGGS case.
