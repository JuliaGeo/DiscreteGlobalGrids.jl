# # Hydrology: a DEM on an IGEO7 grid
#
# Slope, flow accumulation and catchment area are areal quantities. An IGEO7
# cell covers the same area in the Alps as on the equator, so a count of cells
# is an area; the same count on a lon/lat raster is a function of latitude.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))

import DiscreteGlobalGrids as DGG
import Geomorphometry as GM
import GeoInterface as GI
using Rasters, RasterDataSources
import ArchGDAL
import Extents
using Statistics
using GLMakie, GeoMakie
using DiscreteGlobalGridsVisualization: dggsurface, dggsurface!, dggpoly, dggpoly!
GLMakie.activate!(inline = true)

# ## Regrid a Copernicus DEM tile onto IGEO7 hexagons
#
# One 1°×1° tile of the [Copernicus DEM](https://dataspace.copernicus.eu/explore-data/data-collections/copernicus-contributing-missions/collections-description/COP-DEM)
# at 30 m over the Alps, fetched from its
# [AWS Open Data bucket](https://registry.opendata.aws/copernicus-dem/) by
# RasterDataSources:

centre = GI.extent((10.5, 46.5))
path = only(skipmissing(RasterDataSources.getraster(CopernicusDEM; extent = centre)))
Sys.isapple() && Rasters.checkmem!(false) # needed for Apple systems
dem = Raster(path; lazy = false)

# Four lines put it on a grid. `levelfor` picks a level from the raster,
# `MultiOrderCoverage` names the cells the tile touches, and `regrid` fills
# them:

sys = DGG.IGeo7System()
leaf_level = DGG.levelfor(sys, dem)
region = DGG.query(sys, DGG.MultiOrderCoverage(Rasters.extent(dem)); level = leaf_level)
elevation = DGG.regrid(dem; to = region)

# `regrid` returns a `Raster` over a `Cells` axis. Cells outside the tile hold
# its `missingval` — `NaN`, since Copernicus DEM declares none; `missingval =`
# chooses another.
#
# The level `levelfor` chose, and what its cells measure across:

leaf_level, DGG.cellsize(sys, leaf_level)

# Pixels on the left, cells on the right, the same terrain in both:

tile = Rasters.extent(dem)
lims = (tile.X, tile.Y)
shape = AxisAspect(cosd(46.5))            # a degree of longitude is cos(lat) as wide as a degree of latitude
crange = extrema(dem)

fig = Figure(size = (900, 470))
ax1 = Axis(fig[1, 1]; aspect = shape, limits = lims,
    xlabel = "longitude", ylabel = "latitude",
    title = "Copernicus DEM, $(size(dem, 1))×$(size(dem, 2)) pixels")
heatmap!(ax1, dem; colormap = :terrain, colorrange = crange)
ax2 = Axis(fig[1, 2]; aspect = shape, limits = lims,
    xlabel = "longitude", yticklabelsvisible = false,
    title = "IGEO7 level $leaf_level, $(round(length(elevation) / 1e6; digits = 1)) M cells")
p = dggsurface!(ax2, elevation; color = elevation,
    colormap = :terrain, colorrange = crange)
Colorbar(fig[1, 3], p; label = "elevation (m)")
fig

# ## Store sixteen million leaves as a multi-order coverage
#
# `region` is the coarsest set of cells that covers the tile: one big cell
# where the tile is solidly covered, smaller ones towards the border where a
# big cell would overhang. It stands for sixteen million leaves:

length(region), length(DGG.CellLookup(region))

# Coloured by level: a few tens of thousands of coverage cells in place of
# sixteen million ids.

levels = DGG.level.(collect(region))

fig = Figure(size = (620, 540))
ax = Axis(fig[1, 1]; aspect = shape, limits = lims,
    xlabel = "longitude", ylabel = "latitude",
    title = "$(length(region)) coverage cells over $(length(DGG.CellLookup(region))) leaves")
heatmap!(ax, dem; colormap = :grays)
p = dggpoly!(ax, region; color = levels, alpha = 0.75,
    colormap = cgrad(:managua, length(unique(levels)); categorical = true),
    colorrange = (minimum(levels) - 0.5, maximum(levels) + 0.5),
    strokewidth = 0.15, strokecolor = (:black, 0.35))
Colorbar(fig[1, 2], p; label = "cell level", ticks = minimum(levels):maximum(levels))
fig

# A coverage contains the extent, so its border cells overhang the tile by up
# to a kilometre or two. Those leaves have no pixels under them and hold `NaN`:

count(isnan, elevation)

# Explicit `limits` keep the axes on the tile; the overhanging cells would
# otherwise widen them.

# ## Accumulate flow with D8
#
# Geomorphometry's verbs take this raster as it is: they read neighbours and
# cell geometry off its `Cells` axis. One cell's area, in square metres:

area_per_cell = GM.cellarea(elevation, first(eachindex(elevation)))

# `flowaccumulation` returns the upstream area of every cell in square metres
# and the flow directions it used to get there. Its default method is D8, which
# sends all of a cell's water to one neighbour:

accumulation, directions = GM.flowaccumulation(elevation)

# A cell here is `area_per_cell` square metres, so dividing by it counts cells
# upstream:

log_cells = log10.(accumulation ./ area_per_cell)

# A window a dozen kilometres across, taken by indexing the cube with an
# extent. A channel is one cell wide, and over the full tile that is under a
# screen pixel:

window = Extents.Extent(X = (10.30, 10.45), Y = (46.55, 46.68))
here = log_cells[DGG.Cells(DGG.Covering(window))]

fig = Figure(size = (640, 560))
ax = Axis(fig[1, 1]; aspect = shape, limits = (window.X, window.Y),
    xlabel = "longitude", ylabel = "latitude",
    title = "D8 flow accumulation, $(length(here)) cells")
p = dggsurface!(ax, here; color = here,
    colormap = :devon, colorrange = (1, 3.5), highclip = :white)
Colorbar(fig[1, 2], p; label = "log₁₀ cells upstream")
fig

# ## Split the flow between neighbours with D∞
#
# D∞ ([Tarboton 1997](https://doi.org/10.1029/96WR03137)) splits each cell's
# outflow between the two neighbours that bracket its downslope direction, so
# hillslope flow fans out where D8 threads it through single cells. On a
# hexagon those two neighbours are ring slots, two of the six around each cell.
#
# D∞ over the whole tile takes a few minutes; this run covers a padded box
# around the same window, wide enough that the channels entering it arrive with
# their upstream area:

padded = Extents.Extent(X = (10.20, 10.55), Y = (46.45, 46.78))
around = elevation[DGG.Cells(DGG.Covering(padded))]
accumulation_dinf, _ = GM.flowaccumulation(around; method = GM.DInf())
here_dinf = log10.(accumulation_dinf ./ area_per_cell)[DGG.Cells(DGG.Covering(window))]

fig = Figure(size = (640, 560))
ax = Axis(fig[1, 1]; aspect = shape, limits = (window.X, window.Y),
    xlabel = "longitude", ylabel = "latitude",
    title = "D∞ flow accumulation, $(length(here_dinf)) cells")
p = dggsurface!(ax, here_dinf; color = here_dinf,
    colormap = :devon, colorrange = (1, 3.5), highclip = :white)
Colorbar(fig[1, 2], p; label = "log₁₀ cells upstream")
fig

# ## Find each cell's downhill neighbour by hand
#
# Each cell drains to its lowest neighbour. `mapneighbors` visits every cell
# with its ring of neighbours; `needs` names what the kernel receives for each
# of them — here the neighbour's elevation and its index in the cube:

function downhill((z, _), (zs, ids))
    isnan(z) && return (0, NaN32)
    dest, fall = 0, 0.0f0
    for k in eachindex(zs)
        isnan(zs[k]) && continue
        z - zs[k] > fall && ((dest, fall) = (ids[k], z - zs[k]))
    end
    return (dest, dest == 0 ? NaN32 : fall)
end

flow, drop = DGG.mapneighbors(downhill, elevation;
    needs = (DGG.Value(vec(elevation)), DGG.Index(DGG.Local())))
flow

# `flow` is the index each cell drains to and `drop` is the fall to it, one
# raster per component of the tuple the kernel returns. Cells that drain
# nowhere — a pit, or a cell the tile does not cover:

count(==(0), flow)

# Valley floors have little fall and come out dark; headwalls come out bright.

fig = Figure(size = (620, 540))
ax = Axis(fig[1, 1]; aspect = shape, limits = lims,
    xlabel = "longitude", ylabel = "latitude", title = "drop to the downhill neighbour")
p = dggsurface!(ax, drop; color = drop,
    colormap = :magma, colorrange = (0, 30), highclip = :white)
Colorbar(fig[1, 2], p; label = "metres")
fig

# ## Compare the two sets of directions
#
# `directions` from the D8 run holds one downstream neighbour per cell in
# Geomorphometry's local-drain-direction (LDD) code. Encoding `flow` the same
# way makes the two comparable:

cells = collect(lookup(elevation, DGG.Cells))
sends = vec(flow)
downstream(i) = sends[i] == 0 ? cells[i] : cells[sends[i]]
mine = [GM.FlowDirection{GM.LDD}(downstream(i) - cells[i]) for i in eachindex(cells)]

mean(mine .== directions)

# The disagreements are pits, where a one-ring rule stops. The library's D8 is
# a priority flood: it fills each depression and routes the water out over its
# rim. Where the kernel did find a lower neighbour, both pick the same one:

routed = findall(!=(0), sends)
mean(mine[routed] .== vec(directions)[routed])

# ## Which systems each method runs on
#
# `MultiOrderCoverage`, `regrid` and `mapneighbors` are interface methods and
# run on any system: set `sys = DGG.HEALPixSystem()` and nothing else above
# changes. Geomorphometry's verbs on a DGGS raster:
#
# | verb | systems |
# |---|---|
# | `flowaccumulation` with `D8()`, the default | any |
# | `flowaccumulation` with `DInf()` or `FD8()`, and `height_above_nearest_drainage` | IGEO7 — these need relative-cell arithmetic, which the IGEO7 backend provides |
