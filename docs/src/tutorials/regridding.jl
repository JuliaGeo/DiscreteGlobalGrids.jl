# # Regridding a time series
#
# WorldClim ships monthly climatologies on a regular lon/lat grid, where a cell
# near the pole covers far less ground than one at the equator. Here we move all
# twelve months of mean temperature onto an equal-area HEALPix grid with
# first-order conservative regridding, animate the seasonal cycle, and pull out
# a monthly time series for Texas.
#
# The regridder needs no adapter. `treeify`, `ncells` and `getcell` are
# `ConservativeRegridding.Trees`' own bindings, extended for every
# `AbstractGrid` in this package, so a level grid is a regridding endpoint as it
# stands.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))

import DiscreteGlobalGrids as DGG
import ConservativeRegridding as CR
using Rasters, RasterDataSources
import ArchGDAL
import NaturalEarth
import GeometryOps as GO
using Statistics
import Dates
using CairoMakie, GeoMakie
CairoMakie.activate!()

# ## Source rasters and destination grid
#
# RasterDataSources downloads WorldClim through the Rasters extension (about
# 35 MB on first run, cached afterwards); ArchGDAL reads the GeoTIFFs. The
# series holds one global raster per month, 2160×1080 cells at 10 arc-minutes,
# in °C, `missing` over the oceans.

ser = RasterSeries(WorldClim{Climate}, :tavg; month = 1:12, res = "10m")
r = first(ser)
size(r)

# The destination is the whole globe at HEALPix level 6 — 49152 equal-area
# cells, a hair under 1° across. The source is described to
# ConservativeRegridding as the matrix of cell-corner points on the unit sphere.

grid = DGG.levelgrid(DGG.HEALPixSystem(), 6)

(west, east), (south, north) = bounds(r, X), bounds(r, Y)
to_sphere = GO.UnitSpherical.UnitSphereFromGeographic()
corners = [to_sphere((lon, lat))
           for lon in range(west, east; length = size(r, X) + 1),
               lat in range(south, north; length = size(r, Y) + 1)]

# Every grid here computes on the **unit** sphere — `cell_boundary` returns
# `UnitSphericalPoint`s and `cell_area` returns steradians — so the manifold is
# named once, explicitly. A bare corner matrix carries no manifold of its own,
# and mixing in a radius is a factor of `R²` in every area.
#
# `CR.Regridder(manifold, dst, src)` takes the destination first, so its
# `.intersections` matrix is destination cells × source cells. Building it
# intersects every overlapping pair and takes a few seconds.

manifold = GO.Spherical(; radius = 1.0)
regridder = CR.Regridder(manifold, grid, corners)
size(regridder.intersections)

# ## Regridding all twelve months
#
# A conservative mean with gaps needs the standard coverage trick: regrid the
# field with `NaN`s zeroed *and* regrid a 0/1 data-coverage indicator, then
# divide — otherwise coastal cells get dragged down by the empty ocean fraction.
# WorldClim's rows run north to south, so each matrix is flipped to match the
# ascending latitudes before flattening; cells with under 1% coverage stay
# `missing`.

months = 1:12
tavg = Matrix{Union{Float64, Missing}}(missing, DGG.ncells(grid), length(months))
field = zeros(DGG.ncells(grid))
cover = zeros(DGG.ncells(grid))
for m in months
    v = vec(reverse(parent(replace_missing(ser[m], NaN)); dims = 2))
    CR.regrid!(field, regridder, replace(v, NaN => 0.0))
    CR.regrid!(cover, regridder, Float64.(.!isnan.(v)))
    covered = cover .> 0.01
    tavg[covered, m] .= field[covered] ./ cover[covered]
end
count(!ismissing, tavg[:, 1])

# ## The seasonal cycle, animated
#
# `cell_polygon` is the unit-sphere polygon of one cell; one `GO.transform`
# takes the whole vector to lon/lat, and clamping the longitudes to ±180° keeps
# the handful of cells that straddle the antimeridian from smearing across the
# map. One figure, one `Observable` of colours on a fixed colour range, and
# `record` walks it through the months.

cells = [DGG.cellindex(grid, i) for i in 1:DGG.ncells(grid)]
polys = GO.transform(p -> (clamp(p[1], -180.0, 180.0), p[2]),
    GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygon.(Ref(grid), cells)))

crange = extrema(skipmissing(tavg))
month_colors(m) = replace(tavg[:, m], missing => NaN)

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
# Zonal statistics are a `query` away: `Intersects(texas)` returns the cells
# that meet the outline as sorted typed ids, and `cellposition` turns each into
# the row of `tavg` it names. HEALPix cells are equal-area, so the unweighted
# mean over them *is* the areal mean.

states = NaturalEarth.naturalearth("admin_1_states_provinces", 50)
texas = states.geometry[findfirst(==("Texas"), states.name)]
rows = [DGG.cellposition(grid, c) for c in DGG.query(grid, DGG.Intersects(texas))]
ts = [mean(skipmissing(tavg[rows, m])) for m in months]
length(rows)

# About 60 level-6 cells land in Texas, and their mean traces the expected
# cycle: winter around 7 °C, high summer pushing 28 °C.

fig = Figure(size = (600, 350))
ax = Axis(fig[1, 1]; xticks = (months, Dates.monthabbr.(months)),
    ylabel = "mean temperature (°C)", title = "Texas monthly mean temperature")
scatterlines!(ax, months, ts)
fig

# ## The same regridder on every system
#
# The destination grid is the only system-specific token in this page. Swapping
# it changes the cell count and the cell shape, but not the answer: the cell
# covering a given point gets the same January temperature on every system, to
# the resolution of its cells. `cellat` finds that cell and `cellposition` finds
# its row, so the comparison is three interface calls.
#
# Dividing by the coverage matters here for more than the coastline. Both
# numerator and denominator are the *same* row of the intersection matrix
# applied to two source fields, so their ratio is an honest weighted mean of the
# source values whatever the row sums to — which is what makes this the robust
# quantity to compare across grids.

january = vec(reverse(parent(replace_missing(ser[1], NaN)); dims = 2))
probes = ((2.35, 48.86, "Paris"), (-97.7, 30.3, "Austin"), (151.2, -33.9, "Sydney"))

println(rpad("system", 24), "level  cells   ", join(rpad.(last.(probes), 10)))
for sys in (DGG.systems()..., DGG.AuthalicSystem(DGG.HEALPixSystem()))
    base = sys isa DGG.AuthalicSystem ? parent(sys) : sys
    l = base isa Union{DGG.IGeo7System, DGG.H3System} ? 4 : 5
    g = DGG.levelgrid(sys, l)
    rg = CR.Regridder(manifold, g, corners)
    out, cov = zeros(DGG.ncells(g)), zeros(DGG.ncells(g))
    CR.regrid!(out, rg, replace(january, NaN => 0.0))
    CR.regrid!(cov, rg, Float64.(.!isnan.(january)))
    at(lon, lat) = (i = DGG.cellposition(g, DGG.cellat(g, lon, lat));
                    round(out[i] / cov[i]; digits = 2))
    name = sys isa DGG.AuthalicSystem ? "Authalic($(nameof(typeof(base))))" :
           string(nameof(typeof(sys)))
    println(rpad(name, 24), lpad(l, 5), lpad(DGG.ncells(g), 8), "   ",
            join(rpad.([at(p[1], p[2]) for p in probes], 10)))
end
