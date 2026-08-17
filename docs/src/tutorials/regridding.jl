# # Regridding a time series
#
# A temperature field on a regular lon/lat grid has cells that cover far less
# ground near the pole than at the equator. This page moves twelve deterministic
# monthly fields onto an equal-area HEALPix grid with first-order
# conservative regridding, animates the seasonal cycle, and ends with a monthly
# time series regridded onto just the cells that cover Texas.

import DiscreteGlobalGrids as DGG
using Rasters
import NaturalEarth
import GeometryOps as GO
import GeoInterface as GI
import Dates
using Statistics
using GLMakie, GeoMakie
GLMakie.activate!(inline = true)

# ## Source and destination
#
# A one-degree analytic raster keeps the example reproducible and offline. Its
# latitude trend and phase-shifted seasonal cycle are temperature-like; a few
# fixed gaps exercise the same missing-data path as an observational field. The
# twelve months are a dimension of one raster, not twelve rasters.

lon = -179.5:1.0:179.5
lat = 89.5:-1.0:-89.5
months = 1:12
isgap(x, y) = abs(y) > 80 || (20 < x < 80 && -25 < y < 20)
temperature(m, x, y) = isgap(x, y) ? NaN :
    28 - 0.35abs(y) + 2sind(x) + 10sign(y) * cospi((m - 7) / 6)
temps = Raster(
    [temperature(m, x, y) for x in lon, y in lat, m in months],
    (X(lon; sampling = Rasters.Intervals(Rasters.Center())),
     Y(lat; sampling = Rasters.Intervals(Rasters.Center())),
     Dim{:month}(months)),
)
size(temps)

# The destination is HEALPix level 6 — 49152 equal-area cells, a hair under 1°
# across. A grid from this package is a regridding destination as it stands.

grid = DGG.levelgrid(DGG.HEALPixSystem(), 6)

# ## Twelve months, one operator
#
# `plan_regrid` builds the operator and stops there. The weights are geometry —
# how much of each source cell falls in each destination cell — so building them
# is the expensive half and reads no data at all; the twelve months are one
# cheap apply, and a plan applies with no keyword arguments because it already
# carries every answer they would give.
#
# `missingpolicy` says what a partly covered destination cell holds.
# `Weighted(t)` is the coverage-normalized mean of the source values that were
# there — so a coastal cell is not dragged down by the empty gap fraction — and
# blanks a cell covered less than `t`. `Extensive()` is the other choice: raw
# conservative sums, no normalization.

plan = @time DGG.plan_regrid(temps; to = grid, missingpolicy = DGG.Weighted(0.01))
tavg = @time DGG.regrid(temps, plan)

# The spatial dimensions are replaced by the destination's cells and every other
# dimension passes through, so the result is a cube over `(cells, month)`, and
# the cells with under 1% coverage are `NaN`.

dims(tavg)
#
count(!isnan, view(tavg, :, 1))

# ## The seasonal cycle
#
# `cell_polygon` is the unit-sphere polygon of one cell; one `GO.transform`
# takes the whole vector to lon/lat. A few hundred cells straddle the
# antimeridian and would smear into horizontal bands in planar lon/lat, so they
# are masked like the ocean. One `Observable` of colours, and `record` walks it
# through the months.

cells = DGG.CellVector(grid)
polys = GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygon.(Ref(grid), cells))
lonspan(p) = (ring = GI.coordinates(p)[1]; maximum(first, ring) - minimum(first, ring))
seam = lonspan.(polys) .> 180

monthly = parent(tavg)
crange = extrema(filter(!isnan, monthly))
month_colors(m) = ifelse.(seam, NaN, view(monthly, :, m))

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
# cells of the grid that Texas' coverage names, as a lazy `CellVector`, and `to`
# takes that selection directly.

states = NaturalEarth.naturalearth("admin_1_states_provinces", 50)
texas = states.geometry[findfirst(==("Texas"), states.name)]
tx = DGG.covering(DGG.CellVector(grid), texas)
f, a, p = poly(tx; axis = (; aspect = DataAspect()))
lines!(a, texas)
f

# Bare `regrid` builds a plan, applies it and drops it — the one-shot form, for
# an operator with nothing later to reuse it. HEALPix cells are equal-area, so
# the mean over the selected cells is the mean over the ground they cover.

txavg = @time DGG.regrid(temps; to = tx)
ts = vec(mean(parent(txavg); dims = 1))
round.(ts; digits = 1)

# Winter around 5 °C, high summer around 25 °C.

fig = Figure(size = (600, 350))
ax = Axis(fig[1, 1]; xticks = (months, Dates.monthabbr.(months)),
    ylabel = "mean temperature (°C)", title = "Texas monthly mean temperature")
scatterlines!(ax, months, ts)
fig

# Nothing above is HEALPix-specific — any system slots into `levelgrid` — but
# the choice was deliberate: equal-area cells make the means above areal means
# without weights. Conservation does differ by system: a DGGS as a regridding
# *source* is conservative everywhere, but as a *destination* only where its
# cells' rings are convex. IGEO7 and S2 conserve; HEALPix does not, pending an
# upstream clipper fix in GeometryOps. The README and
# `test/systems/crosssystem/regridding_conservation.jl` carry the full account,
# and `missingpolicy = DGG.Weighted(t)` is the normalization that survives it.
