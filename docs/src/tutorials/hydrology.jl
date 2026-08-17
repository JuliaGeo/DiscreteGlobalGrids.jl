# # Hydrology: a DEM on an IGEO7 grid
#
# Terrain analysis wants an equal-area grid: slope, flow accumulation and
# catchment area are all areal quantities, and on a lon/lat raster every one of
# them is a function of latitude. This page moves a Copernicus 30 m DEM tile
# over the Alps onto IGEO7 — hexagons, equal-area by construction — and does
# the first step of a flow-routing model on it. The worked example averages
# the source to 120 m so it also fits comfortably on a standard CI runner.
#
# Three calls carry the page: `MultiOrderCoverage` names the cells the tile
# touches, `regrid` fills them from the raster, and `halo_table` routes water
# out of every cell.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))

import DiscreteGlobalGrids as DGG
import Geomorphometry as GM
using Rasters, RasterDataSources
import ArchGDAL
import Extents
using Statistics
using GLMakie, GeoMakie
GLMakie.activate!(inline = true)

sys = DGG.IGeo7System()

# ## The coarsest cell inside a DEM tile
#
# Copernicus DEM ships in 1°×1° tiles. Which single IGEO7 cell fits inside one?
# Cover the tile with `MultiOrderCoverage` — the coarsest cells that cover it,
# refined only where its outline cuts through — and `coarsest_contained` reads
# off the shallowest cell the traversal proved inside.

tile = Extents.Extent(X = (10.0, 11.0), Y = (46.0, 47.0))
root = DGG.coarsest_contained(DGG.query(sys, DGG.MultiOrderCoverage(tile); level = 10))
f, a, p = poly(Rect2f([tile.X[1], tile.Y[1]], [-(-)(tile.X...), -(-)(tile.Y...)]))
poly!(DGG.PartialGrid(sys, root, DGG.level(root)); strokewidth = 2)
f
# The cell, and its level:
root, DGG.level(root)

# ## The tile's coverage as a grid
#
# One contained cell gives up the tile's rim. The destination that keeps it is
# the tile's own coverage at the working level: every cell the tile touches.

leaf = 12                                          # ≈ 65 m cells
region = DGG.query(sys, DGG.MultiOrderCoverage(tile); level = leaf)
f, a, p = poly(Rect2f([tile.X[1], tile.Y[1]], [-(-)(tile.X...), -(-)(tile.Y...)]))
poly!(region; color = :transparent, strokewidth = 1)
f

# `CellLookup` reads the set as a one-level cell axis, and `PartialGrid` reads
# that as an ordinary grid: positions run `1:ncells`, so a data vector indexes
# straight through it.

grid = DGG.PartialGrid(DGG.CellLookup(region))
DGG.ncells(grid)

# ## Regridding the DEM
#
# The download asks for a point inside the tile rather than the tile itself:
# `getraster` fetches every 1° tile an extent touches, and a closed 1° box
# touches four. Averaging the native 30 m source to 120 m keeps this worked
# example compact.

centre = Extents.Extent(X = (10.5, 10.5), Y = (46.5, 46.5))
path = only(skipmissing(RasterDataSources.getraster(CopernicusDEM; extent = centre)))
dem = Raster(path; lazy = false)
dem = aggregate(mean, dem, 4; progress = false)

# `regrid` takes the grid as its destination and the raster as its source, and
# hands back a cube whose axis is the cells. The coverage overhangs the tile at
# the rim; a cell the raster covers less than half of comes back `NaN` rather
# than as a number standing for ground that was never seen.

igeo7_dem = @time DGG.regrid(dem; to = grid)
elevation = parent(igeo7_dem)
covered = .!isnan.(elevation)
extrema(elevation[covered])

# ## Flow direction
#
# Each cell sends its water to the lowest of its neighbours — the first step of
# every flow-routing model. `halo_table(grid)` is every cell's neighbours, as
# positions into `elevation`, in one call; keeping only the covered ones is
# this page's filter, not the grid's. A cell with no lower neighbour is a pit.

halo = [filter(p -> covered[p], row) for row in DGG.halo_table(grid)]

function downhill(i)
    isempty(halo[i]) && return 0
    j = halo[i][argmin(elevation[halo[i]])]
    return elevation[j] < elevation[i] ? j : 0
end

flow = [covered[i] ? downhill(i) : 0 for i in 1:DGG.ncells(grid)]
drop = [flow[i] == 0 ? NaN : elevation[i] - elevation[flow[i]] for i in eachindex(flow)]
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
p1 = poly!(ax1, cells; color = elevation[shown], colormap = :terrain, strokewidth = 0)
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

tpi = GM.topographic_position_index(terrain)
accumulation, directions = GM.flowaccumulation(terrain; method = GM.D8())

cell_area = GM.cellarea(terrain, first(eachindex(terrain)))
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

# `MultiOrderCoverage`, `PartialGrid`, `regrid` and `halo_table` are interface
# methods, so the regridding and routing portions can use another system.
# `RelativeZ7Cell` and the lazy cell-index iterator provide the IGEO7-specific
# Geomorphometry integration.
