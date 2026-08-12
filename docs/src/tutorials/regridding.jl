# # Regridding a time series
#
# WorldClim ships monthly climatologies on a regular lon/lat grid, where a cell
# near the pole covers far less ground than one at the equator. Here we move all
# twelve months of mean temperature onto an equal-area HEALPix grid with
# first-order conservative regridding, animate the seasonal cycle, and pull out
# a monthly time series for Texas with `zonal`.
#
# ## Source rasters and destination grid
#
# RasterDataSources downloads WorldClim through the Rasters extension (about
# 35 MB on first run, cached afterwards); ArchGDAL is the backend that reads the
# GeoTIFFs. The series holds one global raster per month, 2160×1080 cells at 10
# arc-minutes, in °C, `missing` over the oceans.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))

using DiscreteGlobalGrids
using DiscreteGlobalGrids.HEALPix.HealpixLookups
import ConservativeRegridding as CR
using Rasters, RasterDataSources
import ArchGDAL
import NaturalEarth
import DimensionalData as DD
import GeometryOps as GO
using Statistics
import Dates
using CairoMakie, GeoMakie
CairoMakie.activate!()

ser = RasterSeries(WorldClim{Climate}, :tavg; month = 1:12, res = "10m")
nothing

# The destination is the whole globe at HEALPix level 6 — 49152 equal-area
# cells, a hair under 1° across. The source grid is described to
# ConservativeRegridding as the matrix of cell-corner points on the unit sphere,
# built from the raster's cell edges; `CR.Regridder(dst, src)` takes the
# destination first, and its `.intersections` matrix is destination cells ×
# source cells. Building it intersects every pair of overlapping cells and takes
# a few seconds.

l = HealpixLookup(DGGSGlobeIds(HEALPixDGGS(), 6))

r = first(ser)
(west, east), (south, north) = bounds(r, X), bounds(r, Y)
lon_edges = range(west, east; length = size(r, X) + 1)
lat_edges = range(south, north; length = size(r, Y) + 1)
to_sphere = GO.UnitSpherical.UnitSphereFromGeographic()
corners = [to_sphere((lon, lat)) for lon in lon_edges, lat in lat_edges]

regridder = CR.Regridder(DGGSPartialGrid(l), corners)
nothing

# ## Regridding all twelve months
#
# A conservative mean with gaps needs the standard coverage trick: regrid the
# field with `NaN`s zeroed *and* regrid a 0/1 data-coverage indicator, then
# divide — otherwise coastal cells get dragged down by the empty ocean fraction.
# WorldClim's rows run north to south, so each matrix is flipped to match the
# ascending `lat_edges` before flattening; cells with under 1% coverage stay
# `missing`. The result is one `DimArray` over `(Cells, month)`.

months = 1:12
tavg = Matrix{Union{Float64,Missing}}(missing, length(l), length(months))
field = zeros(length(l))
cover = zeros(length(l))
for m in months
    v = vec(reverse(parent(replace_missing(ser[m], NaN)); dims = 2))
    CR.regrid!(field, regridder, replace(v, NaN => 0.0))
    CR.regrid!(cover, regridder, Float64.(.!isnan.(v)))
    covered = cover .> 0.01
    tavg[covered, m] .= field[covered] ./ cover[covered]
end
A = DD.DimArray(tavg, (Cells(l), DD.Dim{:month}(months)); name = :tavg)
nothing

# ## The seasonal cycle, animated
#
# `cell_polygons` gives one lon/lat polygon per cell, ready for `poly!`;
# clamping their longitudes to ±180° keeps the handful of cells that straddle
# the antimeridian from smearing across the map. One figure, one `Observable`
# of colors on a fixed color range, and `record` walks it through the months.

polys = GO.transform(p -> (clamp(p[1], -180.0, 180.0), p[2]), cell_polygons(l))
crange = extrema(skipmissing(A))
month_colors(m) = collect(replace(parent(A[month = m]), missing => NaN))

colors = Observable(month_colors(1))
month_name = Observable(Dates.monthname(1))
fig = Figure(size = (800, 450))
ax = GeoAxis(fig[1, 1]; dest = "+proj=moll", title = month_name)
plt = poly!(ax, polys; color = colors, colormap = Reverse(:RdYlBu),
    colorrange = crange, strokewidth = 0)
Colorbar(fig[1, 2], plt; label = "mean temperature (°C)")

record(fig, "seasonal_tavg.mp4", months; framerate = 2) do m
    colors[] = month_colors(m)
    month_name[] = Dates.monthname(m)
end
nothing

#md # ```@raw html
#md # <video autoplay loop muted playsinline controls src="./seasonal_tavg.mp4" style="max-width: 100%;"/>
#md # ```

# ## A Texas time series
#
# `zonal` is qualified here because Rasters exports one too. HEALPix cells are
# equal-area, so the unweighted mean over the cells whose centers fall in Texas
# *is* the areal mean, and `missing` cells are skipped by default.

states = NaturalEarth.naturalearth("admin_1_states_provinces", 50)
texas = states.geometry[findfirst(==("Texas"), states.name)]
ts = [only(HealpixLookups.zonal(mean, A[month = m]; of = texas)) for m in months]
nothing

# About 60 level-6 cells land in Texas, and their mean traces the expected
# cycle: winter around 7 °C, high summer pushing 28 °C.

fig = Figure(size = (600, 350))
ax = Axis(fig[1, 1]; xticks = (months, Dates.monthabbr.(months)),
    ylabel = "mean temperature (°C)", title = "Texas monthly mean temperature")
scatterlines!(ax, months, Float64.(ts))
fig
