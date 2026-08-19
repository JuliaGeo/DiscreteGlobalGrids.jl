# # Hydrology: a DEM on an IGEO7 grid
#
# Terrain analysis wants an equal-area grid: slope, flow accumulation and
# catchment area are all areal quantities, and on a lon/lat raster every one of
# them is a function of latitude. This page moves a Copernicus 30 m DEM tile
# over the Alps onto IGEO7 — hexagons, equal-area by construction — and does
# the first step of a flow-routing model on it. The worked example averages
# the source to 240 m and works at a level whose cells are a little coarser
# again, so a tile's worth of them fits comfortably on a standard CI runner.
#
# Three calls carry the page: `MultiOrderCoverage` names the cells the tile
# touches, `regrid` fills them from the raster, and `adjacency` routes water
# out of every cell.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))

import DiscreteGlobalGrids as DGG
import Geomorphometry as GM
import GeoInterface as GI, GeometryOps as GO
using Rasters, RasterDataSources
import ArchGDAL
import Extents
using Statistics
using GLMakie, GeoMakie
GLMakie.activate!(inline = true)

# ## Acquiring data

# Let's first get a DEM tile and regrid it to IGeo7.
# We'll use a Copernicus DEM 30m-tile over the Alps for this example.

centre = GI.extent((10.5, 46.5))
path = only(skipmissing(RasterDataSources.getraster(CopernicusDEM; extent = centre)))
Sys.isapple() && Rasters.checkmem!(false)
dem = Raster(path; lazy = false)
## dem = aggregate(mean, dem, 8; progress = false)
# This is what the raster looks like:
plot(dem; axis = (; aspect = DataAspect()))
# Now, let's regrid it to IGeo7.  First
sys = DGG.IGeo7System()
# We can get the grid of IGeo7 cells that would cover the DEM tile,
# using a `MultiOrderCoverage` query:
leaf_level = 13
region = DGG.query(sys, DGG.MultiOrderCoverage(Rasters.extent(dem)); level = leaf_level)
# Here's what this looks like:
f, a, p = plot(dem; axis = (; aspect = DataAspect()))
poly!(a, region; color = :transparent, strokewidth = 1, strokecolor = (:black, 0.5))
f
# This is a nice way to compress the set of cells that would be covered in memory.
# Note that `region` says it has ~44,424 cells.  But when you look at the number of 
# cells at level 13,
DGG.CellLookup(region) |> length
# That's a lot of cells!  This optimization helps to decrease memory pressure,
# especially on datasets that don't fit in memory in the first place.
# 
# If you want, you can also choose the largest cell in the region to run this
# analysis on:
largest_contained_cell = DGG.coarsest_contained(region)

# `CellLookup` reads the set as a one-level cell axis, and `PartialGrid` reads
# that as an ordinary grid: positions run `1:ncells`, so a data vector indexes
# straight through it. This is where the leaf level is paid for: the axis is
# every level-`leaf` cell under the set, so its length grows sevenfold per
# level and it, not the set, sizes everything below.

grid = DGG.PartialGrid(DGG.CellLookup(region))

# ## Regridding the DEM
#
# The download asks for a point inside the tile rather than the tile itself:
# `getraster` fetches every 1° tile an extent touches, and a closed 1° box
# touches four. Averaging the native 30 m source to 240 m keeps this worked
# example compact and leaves the source a little finer than the destination,
# so every hexagon averages several pixels.

# [`regrid`](@ref) takes the grid as its destination and the raster as its source, and
# hands back a cube whose axis is the cells. The coverage overhangs the tile by
# a bit at the border.

igeo7_dem = @time DGG.regrid(dem; to = region)
# Let's now plot this too:
poly(lookup(igeo7_dem, DGG.Cells); color = vec(igeo7_dem); axis = (; aspect = DataAspect()))
# and we can compute the elevation of the cells:
covered = .!isnan.(igeo7_dem)
extrema(igeo7_dem[covered])
# Finally, let's extract the grid backing the lookup,
# which we can pass to the DGG neighbors API.
grid = DGG.PartialGrid(lookup(igeo7_dem, DGG.Cells))
# ## Flow direction
#
# Each cell sends its water to the lowest of its neighbours — the first step of
# every flow-routing model. `adjacency(grid)` is every cell's neighbours, as
# positions into `elevation`, in one call; keeping only the covered ones is
# this page's filter, not the grid's. A cell with no lower neighbour is a pit.

nbrs = [filter(p -> covered[p], row) for row in DGG.adjacency(grid)]

function downhill(i)
    isempty(nbrs[i]) && return 0
    j = nbrs[i][argmin(elevation[nbrs[i]])]
    return elevation[j] < elevation[i] ? j : 0
end

flow = [covered[i] ? downhill(i) : 0 for i in 1:DGG.ncells(grid)]
drop = [flow[i] == 0 ? NaN : igeo7_dem[i] - igeo7_dem[flow[i]] for i in eachindex(flow)]
(; n_pits = count(i -> covered[i] && flow[i] == 0, eachindex(flow)),
   max_drop = maximum(filter(!isnan, drop)))

# Elevation, and the drop to the downhill neighbour — the drop map picks out
# valley floors as the flat regions and headwalls as the steep ones.
# `CellVector(grid)[shown]` is the covered cells as a vector, which `poly!`
# draws directly, in the same order `elevation[shown]` follows.

shown = findall(covered)
cells = DGG.CellVector(grid)[shown]

fig = Figure(size = (900, 430))
ax1 = GeoAxis(fig[1, 1]; dest = "+proj=longlat +datum=WGS84", title = "elevation (m)")
p1 = poly!(ax1, cells; color = igeo7_dem[shown], colormap = :terrain, strokewidth = 0)
Colorbar(fig[2, 1], p1; vertical = false)
ax2 = GeoAxis(fig[1, 2]; dest = "+proj=longlat +datum=WGS84", title = "drop to downhill neighbour (m)")
p2 = poly!(ax2, cells; color = drop[shown], colormap = :magma,
    nan_color = :gray80, strokewidth = 0)
Colorbar(fig[2, 2], p2; vertical = false)
fig

# ## Terrain analysis with Geomorphometry
#
# The cube's axis is already a `CellLookup`: canonical IGEO7 cell identities
# over positional vector storage. Geomorphometry reads the DGGS neighbourhood
# and physical cell geometry off it directly. TPI is a single local call; flow
# accumulation uses D8 here, with each result expressed as upstream area in
# square metres.

terrain = Raster(igeo7_dem; name = :height)

tpi = @time GM.topographic_position_index(terrain)
accumulation, directions = @time GM.flowaccumulation(terrain; method = GM.D8())

cell_area = @time GM.cellarea(terrain, first(eachindex(terrain)))
log_cells = log10.(accumulation ./ cell_area)

fig = Figure(size = (900, 430))
ax1 = GeoAxis(fig[1, 1]; dest = "+proj=longlat +datum=WGS84",
    title = "topographic position index (m)")
p1 = poly!(ax1, cells; color = vec(tpi[shown]), colorrange = (-25, 25),
    colormap = :delta, strokewidth = 0)
Colorbar(fig[2, 1], p1; vertical = false)
ax2 = GeoAxis(fig[1, 2]; dest = "+proj=longlat +datum=WGS84",
    title = "D8 flow accumulation")
p2 = poly!(ax2, cells; color = vec(log_cells[shown]),
    colormap = :devon, strokewidth = 0)
Colorbar(fig[2, 2], p2; vertical = false,
    label = "log₁₀(upstream cell equivalents)")
save("geomorphometry_igeo7.png", fig)
fig

# `MultiOrderCoverage`, `subtree`, `regrid` and `adjacency` are interface
# methods, so the regridding and routing portions can use another system.
# `RelativeZ7Cell` and the lazy cell-index iterator provide the IGEO7-specific
# Geomorphometry integration.
