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
using Rasters, RasterDataSources
import ArchGDAL
import GeometryOps as GO, Extents
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
# This is the actual cell we have:
root, DGG.level(root)

# ## One subtree as a grid
#
# `PartialGrid(sys, cell, level)` is that cell's subtree at a working level, as
# an ordinary grid: positions run `1:ncells`, so a data vector indexes straight
# through it.

leaf = 13                                          # ≈ 430 m cells
grid = DGG.PartialGrid(sys, root, leaf)
#
DGG.ncells(grid)

# ## Regridding the DEM
#
# The download asks for a point inside the tile rather than the tile itself:
# `getraster` fetches every 1° tile an extent touches, and a closed 1° box
# touches four. The source stays at its native 30 m resolution; the source
# quadtree keeps intersection work localized while the matrix is built.

Rasters.checkmem!(false)                           # the tile is bigger than free RAM
centre = Extents.Extent(X = (10.5, 10.5), Y = (46.5, 46.5))
path = only(skipmissing(RasterDataSources.getraster(CopernicusDEM; extent = centre)))
dem = Raster(path; lazy = false)
dem = set(dem, X => Rasters.Intervals(Rasters.Start()), Y => Rasters.Intervals(Rasters.Start()))
# ConservativeRegridding wants the source's cell corners on the unit sphere.
# The destination needs no adapter: `treeify`, `ncells` and `getcell` are
# extended for every `AbstractGrid`.
#
# The source could skip the adapter too: `DGG.CopernicusDEMSystem(30)` is this
# exact lattice as a grid system, so `DGG.PartialGrid(sys, tile, 1)` is a whole
# tile as an ordinary grid — no corner matrix, no `CellBasedGrid`. This page
# keeps the hand-built route because it is the one that works for *any* raster;
# `examples/copernicus_dem.jl` is the tile-native one, end to end.

xbounds, ybounds = Rasters.intervalbounds(dem, (X, Y))
# `intervalbounds` follows array-index order while each pair remains `(low,
# high)`. The raster's X lookup runs west-to-east, so its corner sequence starts
# at the first low edge. Its Y lookup runs north-to-south, so that sequence
# starts at the first high edge. Cell `(i, j)` then describes `dem[i, j]`
# directly, without copying or reordering the DEM.
xrange = [first(first(xbounds)); last.(xbounds)]
yrange = [last(first(ybounds)); first.(ybounds)]
latlong_point_matrix = [GO.UnitSphericalPoint((x, y)) for x in xrange, y in yrange]
latlong_grid = CR.Trees.CellBasedGrid(GO.Spherical(), latlong_point_matrix) |> CR.Trees.TopDownQuadtreeCursor

regridder = @time CR.Regridder(GO.Spherical(; radius = 1.0), grid, latlong_grid)
# Just for fun, let's look at the sparsity pattern of the regridder matrix:
spy(regridder.intersections; axis = (; aspect = DataAspect(), title = "Sparsity pattern of lat-long -> IGEO7 regridder"))

# The tile need not cover every rim cell of the subtree, so regrid the field
# and a 0/1 indicator with the same matrix and divide. The ratio is a weighted
# mean of the source values — right for a partly covered cell, where the raw
# number is not.

source = vec(dem)
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

# Only the singleton on the first line was IGEO7: `MultiOrderCoverage`,
# `PartialGrid`, the regridder and `halo_table` are interface methods, so
# swapping the system reruns the page unchanged.
