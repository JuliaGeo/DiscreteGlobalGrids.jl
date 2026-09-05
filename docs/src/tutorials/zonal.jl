# # Zonal statistics
#
# A zonal statistic reduces the cells in each region to one value. This example
# estimates each country's mean July temperature from a HEALPix grid.

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import NaturalEarth
import NCDatasets  # the netCDF backend Rasters reads CRU through
using Rasters, RasterDataSources
using Statistics
using DiscreteGlobalGridsVisualization: dggpoly, dggpoly!
using GLMakie, GeoMakie
GLMakie.activate!(inline = true)

# ## Load the July temperature raster
#
# CRU CL 2.0 is a station climatology of 1961–1990 over land, at 10 arcmin.
# Its `:tmp` variable holds mean temperature in °C with one layer per month,
# and July is layer 7.

tavg = read(Raster(RasterDataSources.getraster(CRUCL2); name = :tmp, lazy = true)[Ti = 7])

# Ocean cells hold `missing`, and the source record has limited Antarctic
# coverage.

fig, ax, plt = heatmap(tavg; colormap = :thermal,
    axis = (; aspect = DataAspect(), title = "CRU CL 2.0 mean July temperature"))
Colorbar(fig[1, 2], plt; label = "°C")
fig

# ## Regrid the raster onto HEALPix level 6
#
# [Regridding](regridding.md) returns a `Raster` whose single dimension is
# `Cells`: the grid's cells, in the grid's order. That axis is what every
# selector below indexes.

grid = DGG.levelgrid(DGG.HEALPixSystem(), 6)
field = DGG.regrid(tavg; to = grid)

# ## Average the field over every country
#
# `field[Cells(Covering(geom))]` selects a cell coverage for a polygon. Applying
# `mean` to that selection gives a grid-based estimate; any reduction can be
# used in its place.

countries = NaturalEarth.naturalearth("admin_0_countries", 50)

#

percountry = [mean(skipmissing(field[DGG.Cells(DGG.Covering(g))]))
              for g in countries.geometry]

# At this resolution, small islands may have no selected land cell, and the
# source has no observations for some Antarctic regions. Those means are
# `NaN`; the ranking places them at the end and the map draws them grey.

ranked = sort(countries.NAME .=> round.(percountry; digits = 1); by = last)
last(filter(!isnan ∘ last, ranked), 5)

# The five coldest:

first(ranked, 5)

#

crange = extrema(skipmissing(field))
fnan = Rasters.replace_missing(field, NaN)

fig = Figure(size = (900, 900))
ax1 = GeoAxis(fig[1, 1]; dest = "+proj=eqearth",
    title = "July mean temperature on HEALPix level 6",
    xgridcolor = (:black, 0.15), ygridcolor = (:black, 0.15))
dggpoly!(ax1, fnan; color = fnan, colormap = :thermal, colorrange = crange)
lines!(ax1, GeoMakie.coastlines(); color = (:black, 0.45), linewidth = 0.5)
Colorbar(fig[1, 2]; colormap = :thermal, colorrange = crange, label = "°C")
ax2 = GeoAxis(fig[2, 1]; dest = "+proj=eqearth",
    title = "mean July temperature per country",
    xgridcolor = (:black, 0.15), ygridcolor = (:black, 0.15))
poly!(ax2, countries.geometry; color = percountry, colormap = :thermal,
    colorrange = crange, strokecolor = (:black, 0.4), strokewidth = 0.4,
    nan_color = :lightgray)
Colorbar(fig[2, 2]; colormap = :thermal, colorrange = crange, label = "°C")
fig

# `Covering` returns a cell set that contains the polygon and may include a rim
# extending beyond the border. The country value is therefore an approximation,
# not an exact polygon statistic. [Count only cells wholly inside the outline](@ref)
# is a stricter alternative.
#
# ## Select the cells covering one region
#
# The same selector on one state gives a `Raster` over the cells covering it.

states = NaturalEarth.naturalearth("admin_1_states_provinces", 50)
texas = states.geometry[findfirst(==("Texas"), states.name)]

tx = field[DGG.Cells(DGG.Covering(texas))]

#

mean(skipmissing(tx))

# HEALPix cells have equal area, so the plain mean is area-weighted for this
# grid. On a system with unequal cells, weight by
# `DGG.cell_area.(grid, DD.lookup(tx, DGG.Cells))`.
#
# ## Draw the covering cells of Texas
#
# `replace_missing` turns the Gulf-coast cells into `NaN`, and `dggpoly!`
# leaves those undrawn.

txn = Rasters.replace_missing(tx, NaN)

fig = Figure(size = (760, 620))
ax = GeoAxis(fig[1, 1]; dest = "+proj=longlat +datum=WGS84",
    limits = ((-110.0, -91.0), (23.5, 39.0)), xticks = -108:2:-92,
    yticks = 24:2:38, xgridcolor = (:black, 0.15), ygridcolor = (:black, 0.15),
    title = "HEALPix level-6 cells covering Texas")
dggpoly!(ax, txn; color = txn, colormap = :thermal,
    strokecolor = (:white, 0.6), strokewidth = 0.5)
poly!(ax, texas; color = :transparent, strokecolor = :black, strokewidth = 2)
Colorbar(fig[1, 2]; colormap = :thermal, colorrange = extrema(skipmissing(tx)),
    label = "°C")
fig

# ## Count only cells wholly inside the outline
#
# These boundary rules select progressively narrower sets of cells. A raster
# zonal tool commonly uses the centre-in-zone rule for its pixels.
#
# | rule | cells kept | spelling |
# |---|---|---|
# | `Covering` | a cell set containing the outline, possibly with an outer rim | `field[Cells(Covering(geom))]` |
# | centre-in-zone | every cell whose centre is inside | — |
# | `Within` | every cell wholly inside the outline | `field[Cells(Within(geom))]` |

inside = field[DGG.Cells(DGG.Within(texas))]

#

mean(skipmissing(tx)), mean(skipmissing(inside))

# The two means differ because boundary cells change which ground contributes
# to the statistic.
#
# ## Read the cell holding a point
#
# `DD.Contains` takes a position in degrees and selects the cell holding it.
# Austin sits at 97.74° W, 30.27° N.

field[DGG.Cells(DD.Contains((-97.74, 30.27)))]

# ## Run the same statistic on another system
#
# The same selectors work with another grid system. Here `levelfor` chooses an
# IGEO7 level with cells approximately 100 km across:

igeo7 = DGG.IGeo7System()
level = DGG.levelfor(igeo7, 100_000)

#

field7 = DGG.regrid(tavg; to = DGG.levelgrid(igeo7, level))

#

mean(skipmissing(field7[DGG.Cells(DGG.Covering(texas))]))

# The width of the covering rim depends on the grid's refinement scheme;
# [Multi-order coverage](multiorder.md) shows where that width comes from.
#
# ## Without the DimArray
#
# `CellVector` exposes the grid as a vector of cell ids. The two functions below
# show the index operations used by `Covering` and `Contains` on the `Cells`
# axis:

cells = DGG.CellVector(grid)
DGG.covering_indices(cells, texas)

#

DGG.localindex(cells, -97.74, 30.27)
