# # Regridding a time series
#
# WorldClim ships monthly mean temperature on a regular lon/lat grid, where a
# cell near the pole covers far less ground than one at the equator. This page
# moves all twelve months onto an equal-area HEALPix grid with first-order
# conservative regridding, animates the seasonal cycle, and ends with a monthly
# time series regridded onto just the cells that cover Texas.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))

import DiscreteGlobalGrids as DGG
import ConservativeRegridding as CR
using Rasters, RasterDataSources
import ArchGDAL
import NaturalEarth
import GeometryOps as GO
import GeoInterface as GI
using Statistics
import Dates
using CairoMakie, GeoMakie
CairoMakie.activate!()

# ## Source and destination
#
# RasterDataSources downloads WorldClim (about 35 MB once, cached afterwards):
# one global raster per month, 2160×1080 cells at 10 arc-minutes, in °C,
# `missing` over the oceans.

ser = RasterSeries(WorldClim{Climate}, :tavg; month = 1:12, res = "10m")
r = first(ser)
size(r)

# The destination is HEALPix level 6 — 49152 equal-area cells, a hair under 1°
# across. A grid from this package is a regridding endpoint as it stands; the
# lon/lat source is described as its matrix of cell-corner points on the unit
# sphere.

grid = DGG.levelgrid(DGG.HEALPixSystem(), 6)

(west, east), (south, north) = bounds(r, X), bounds(r, Y)
to_sphere = GO.UnitSpherical.UnitSphereFromGeographic()
corners = [to_sphere((lon, lat))
           for lon in range(west, east; length = size(r, X) + 1),
               lat in range(south, north; length = size(r, Y) + 1)]
size(corners)

# The manifold is passed explicitly: every grid here computes on the **unit**
# sphere, and a guessed manifold would be a WGS84 sphere — a factor of `R²` in
# every area. `Regridder` takes the destination first; building it clips every
# overlapping pair of cells and takes a few seconds.

manifold = GO.Spherical(; radius = 1.0)
regridder = CR.Regridder(manifold, grid, corners)
size(regridder.intersections)

# ## Twelve months, with gaps
#
# A field with gaps needs one trick: regrid the field with its `NaN`s zeroed
# *and* a 0/1 coverage indicator, then divide — otherwise coastal cells are
# dragged down by the empty ocean fraction. The ratio is a genuine weighted
# mean of the covered source values whatever each row of the intersection
# matrix sums to, and that robustness matters here: a DGGS as a regridding
# *source* is conservative on every system, but as a *destination* only where
# its cells' rings are convex — HEALPix's are not — pending an upstream clipper
# fix in GeometryOps. The README and
# `test/systems/crosssystem/regridding_conservation.jl` carry the full account.
#
# WorldClim's rows run north to south, so each month is flipped to the corner
# mesh's ascending latitudes before flattening; cells with under 1% coverage
# stay `missing`.

months = 1:12
month_field(m) = vec(reverse(parent(replace_missing(ser[m], NaN)); dims = 2))

tavg = Matrix{Union{Float64, Missing}}(missing, DGG.ncells(grid), length(months))
field = zeros(DGG.ncells(grid))
cover = zeros(DGG.ncells(grid))
for m in months
    v = month_field(m)
    CR.regrid!(field, regridder, replace(v, NaN => 0.0))
    CR.regrid!(cover, regridder, Float64.(.!isnan.(v)))
    ok = cover .> 0.01
    tavg[ok, m] .= field[ok] ./ cover[ok]
end
count(!ismissing, tavg[:, 1])

# ## The seasonal cycle
#
# `cell_polygon` is the unit-sphere polygon of one cell; one `GO.transform`
# takes the whole vector to lon/lat. A few hundred cells straddle the
# antimeridian and would smear into horizontal bands in planar lon/lat, so they
# are masked like the ocean. One `Observable` of colours, and `record` walks it
# through the months.

cells = [DGG.cellindex(grid, i) for i in 1:DGG.ncells(grid)]
polys = GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygon.(Ref(grid), cells))
lonspan(p) = (ring = GI.coordinates(p)[1]; maximum(first, ring) - minimum(first, ring))
seam = lonspan.(polys) .> 180

crange = extrema(skipmissing(tavg))
month_colors(m) = ifelse.(seam, NaN, replace(tavg[:, m], missing => NaN))

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
nothing #hide

#md # ```@raw html
#md # <video autoplay loop muted playsinline controls src="./seasonal_tavg.mp4" style="max-width: 100%;"/>
#md # ```

# ## A Texas time series
#
# For one region there is no need to regrid the globe. `covering` selects the
# cells of the grid that Texas' coverage names, as a lazy `CellVector`, and
# `PartialGrid` reads that selection back as a grid — in O(1) — so it is a
# regridding destination as it stands.

states = NaturalEarth.naturalearth("admin_1_states_provinces", 50)
texas = states.geometry[findfirst(==("Texas"), states.name)]
tx = DGG.covering(DGG.CellVector(grid), texas)
length(tx)

# The same divide as above, on the small regridder. HEALPix cells are
# equal-area, so the unweighted mean over them is the areal mean.

rgtx = CR.Regridder(manifold, DGG.PartialGrid(tx), corners)
f, c = zeros(length(tx)), zeros(length(tx))
ts = map(months) do m
    v = month_field(m)
    CR.regrid!(f, rgtx, replace(v, NaN => 0.0))
    CR.regrid!(c, rgtx, Float64.(.!isnan.(v)))
    mean(f ./ c)
end
round.(ts; digits = 1)

# Winter around 7 °C, high summer pushing 28 °C.

fig = Figure(size = (600, 350))
ax = Axis(fig[1, 1]; xticks = (months, Dates.monthabbr.(months)),
    ylabel = "mean temperature (°C)", title = "Texas monthly mean temperature")
scatterlines!(ax, months, ts)
fig
