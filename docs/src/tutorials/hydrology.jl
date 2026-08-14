# # Hydrology: a DEM on an IGEO7 grid
#
# Terrain analysis wants an equal-area grid: slope, flow accumulation and
# catchment area are all areal quantities, and on a lon/lat raster every one of
# them is a function of latitude. This page moves a Copernicus 30 m DEM tile
# over the Alps onto IGEO7 — hexagons, equal-area by construction — and does
# the first step of a flow-routing model on it.
#
# Three calls carry the page: `coarsest_contained` picks the cell to work in,
# `PartialGrid` names its subtree as a grid, and `halo_table` routes water out
# of every cell.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))

import DiscreteGlobalGrids as DGG
import ConservativeRegridding as CR
import Geomorphometry as GM
using Rasters, RasterDataSources
import ArchGDAL
import GeometryOps as GO, Extents
using Statistics
using CairoMakie, GeoMakie
CairoMakie.activate!()

sys = DGG.IGeo7System()

# ## The coarsest cell inside a DEM tile
#
# Copernicus DEM ships in 1°×1° tiles. Which single IGEO7 cell fits inside one?
# Cover the tile with `MultiOrderCoverage` — the coarsest cells that cover it,
# refined only where its outline cuts through — and `coarsest_contained` reads
# off the shallowest cell the traversal proved inside.

tile = Extents.Extent(X = (10.0, 11.0), Y = (46.0, 47.0))
root = DGG.coarsest_contained(DGG.query(sys, DGG.MultiOrderCoverage(tile); level = 10))
root, DGG.level(root)

# ## One subtree as a grid
#
# `PartialGrid(sys, cell, level)` is that cell's subtree at a working level, as
# an ordinary grid: positions run `1:ncells`, so a data vector indexes straight
# through it.

leaf = 10                                          # ≈ 430 m cells
grid = DGG.PartialGrid(sys, root, leaf)
DGG.ncells(grid)

# ## The DEM, conservatively regridded
#
# The download asks for a point inside the tile rather than the tile itself:
# `getraster` fetches every 1° tile an extent touches, and a closed 1° box
# touches four. At 30 m the tile is far finer than the destination cells, so
# `aggregate` averages it down before the intersection matrix is built.

Rasters.checkmem!(false)                           # the tile is bigger than free RAM
centre = Extents.Extent(X = (10.5, 10.5), Y = (46.5, 46.5))
path = only(skipmissing(RasterDataSources.getraster(CopernicusDEM; extent = centre)))
dem = aggregate(mean, Raster(path; lazy = true), 16)

# ConservativeRegridding wants the source's cell corners on the unit sphere.
# The destination needs no adapter: `treeify`, `ncells` and `getcell` are
# extended for every `AbstractGrid`.

(west, east), (south, north) = bounds(dem, X), bounds(dem, Y)
to_sphere = GO.UnitSpherical.UnitSphereFromGeographic()
corners = [to_sphere((lon, lat))
           for lon in range(west, east; length = size(dem, X) + 1),
               lat in range(south, north; length = size(dem, Y) + 1)]
regridder = CR.Regridder(GO.Spherical(; radius = 1.0), grid, corners)

# The tile need not cover every rim cell of the subtree, so regrid the field
# and a 0/1 indicator with the same matrix and divide. The ratio is a weighted
# mean of the source values — right for a partly covered cell, where the raw
# number is not.

source = vec(Float64.(reverse(parent(dem); dims = 2)))
raw = zeros(DGG.ncells(grid))
cover = zeros(DGG.ncells(grid))
CR.regrid!(raw, regridder, source)
CR.regrid!(cover, regridder, ones(length(source)))
covered = cover .> 0.5
elevation = fill(NaN, DGG.ncells(grid))
elevation[covered] .= raw[covered] ./ cover[covered]
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
# `getcell(grid, i)` is the polygon at position `i`; without an argument it is
# all of them.

lonlat(g) = GO.transform(GO.GeographicFromUnitSphere(), g)
polys = map(lonlat, DGG.getcell(grid))
shown = findall(covered)

fig = Figure(size = (900, 430))
ax1 = GeoAxis(fig[1, 1]; dest = "+proj=longlat +datum=WGS84", title = "elevation (m)")
p1 = poly!(ax1, polys[shown]; color = elevation[shown], colormap = :terrain, strokewidth = 0)
Colorbar(fig[2, 1], p1; vertical = false)
ax2 = GeoAxis(fig[1, 2]; dest = "+proj=longlat +datum=WGS84", title = "drop to downhill neighbour (m)")
p2 = poly!(ax2, polys[shown]; color = drop[shown], colormap = :magma,
    nan_color = :gray80, strokewidth = 0)
Colorbar(fig[2, 2], p2; vertical = false)
fig

# ## Terrain analysis with Geomorphometry
#
# A `CellLookup` gives the raster canonical IGEO7 cell identities while keeping
# its vector storage positional. Geomorphometry can then use the DGGS
# neighbourhood and physical cell geometry directly. TPI is a single local
# call; flow accumulation uses D8 here, with each result expressed as upstream
# area in square metres.

igeo7_dem = Raster(elevation,
    (DGG.Cells(DGG.CellLookup(DGG.CellVector(grid))),);
    name = :height)

tpi = GM.topographic_position_index(igeo7_dem)
accumulation, directions = GM.flowaccumulation(igeo7_dem; method = GM.D8())

cell_area = GM.cellarea(igeo7_dem, first(eachindex(igeo7_dem)))
log_cells = log10.(accumulation ./ cell_area)

fig = Figure(size = (900, 430))
ax1 = GeoAxis(fig[1, 1]; dest = "+proj=longlat +datum=WGS84",
    title = "topographic position index (m)")
p1 = poly!(ax1, polys[shown]; color = tpi[shown], colorrange = (-25, 25),
    colormap = :delta, strokewidth = 0)
Colorbar(fig[2, 1], p1; vertical = false)
ax2 = GeoAxis(fig[1, 2]; dest = "+proj=longlat +datum=WGS84",
    title = "D8 flow accumulation")
p2 = poly!(ax2, polys[shown]; color = log_cells[shown],
    colormap = :devon, strokewidth = 0)
Colorbar(fig[2, 2], p2; vertical = false,
    label = "log₁₀(upstream cell equivalents)")
save("geomorphometry_igeo7.png", fig)
fig

# `MultiOrderCoverage`, `PartialGrid`, the regridder and `halo_table` are
# interface methods, so the regridding and routing portions can use another
# system. The relative indexing methods added for IGEO7 provide the
# Geomorphometry integration.
