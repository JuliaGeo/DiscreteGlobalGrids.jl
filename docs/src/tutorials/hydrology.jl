# # Hydrology: a DEM on an IGEO7 grid
#
# Terrain analysis wants an equal-area grid: slope, flow accumulation and
# catchment area are all areal quantities, and on a lon/lat raster every one of
# them is a function of latitude. This page moves a Copernicus 30 m DEM tile
# over the Alps onto IGEO7 — hexagons, equal-area by construction — and then
# does the first step of a flow-routing model on it.
#
# Three interface calls carry the whole page: `MultiOrderCoverage` to find the
# cells that cover a region, `PartialGrid` to name one subtree as a grid, and
# `halo_table` to route water out of each cell.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))

import DiscreteGlobalGrids as DGG
import ConservativeRegridding as CR
using Rasters, RasterDataSources
import ArchGDAL, NaturalEarth
import GeometryOps as GO, GeoInterface as GI, Extents
using Statistics
using CairoMakie, GeoMakie
CairoMakie.activate!()

# CairoMakie rather than GLMakie or WGLMakie: this page saves a PNG, and Cairo
# is the backend that renders `poly!` to a file without a display.

sys = DGG.IGeo7System()

# ## Covering a region, at mixed resolution
#
# `MultiOrderCoverage` returns the **coarsest** cells that cover a target,
# refining only where the outline cuts through: the interior arrives as a few
# big cells and the boundary as many small ones. That is the whole of "which
# cells does Switzerland touch", in one call.

countries = NaturalEarth.naturalearth("admin_0_countries", 10)
switzerland = countries.geometry[findfirst(==("Switzerland"), countries.NAME)]

coverage = DGG.query(sys, DGG.MultiOrderCoverage(switzerland); level = 8)
(; n_cells = length(coverage),
   levels = extrema(DGG.level, coverage),
   n_leaves = sum(length, DGG.level_ranges(coverage, 8)))

# Twenty-odd thousand level-8 cells, held as a few thousand mixed-level ones.
# Drawing it, coloured by level, is the picture of what the compression buys:

cellpoly(c) = GO.transform(GO.GeographicFromUnitSphere(),
    DGG.cell_polygon(DGG.levelgrid(sys, DGG.level(c)), c))

fig = Figure(size = (760, 520))
ax = GeoAxis(fig[1, 1]; dest = "+proj=longlat +datum=WGS84",
    limits = ((5.5, 10.8), (45.6, 48.0)), title = "IGEO7 multi-order coverage of Switzerland")
plt = poly!(ax, cellpoly.(collect(coverage)); color = DGG.level.(collect(coverage)),
    colormap = :viridis, strokecolor = (:white, 0.6), strokewidth = 0.3)
poly!(ax, switzerland; color = :transparent, strokecolor = :black, strokewidth = 2)
Colorbar(fig[1, 2], plt; label = "cell level")
fig

# ## The coarsest cell inside a DEM tile
#
# Copernicus DEM ships in 1°×1° tiles. Which single IGEO7 cell fits inside one?
# The coverage of the tile answers it: its shallowest member is the coarsest
# cell the traversal could emit whole, which is the coarsest cell contained by
# the tile.

tile = Extents.Extent(X = (10.0, 11.0), Y = (46.0, 47.0))
fitting = argmin(DGG.level, DGG.query(sys, DGG.MultiOrderCoverage(tile); level = 12))
fitting, DGG.level(fitting)

# It really is inside — `Within` is the predicate that says so, and asking the
# grid at that cell's own level gives the same answer the traversal did. (The
# shallowest member is a contained cell only when the traversal got to emit one
# whole; a target smaller than a max-depth cell yields only boundary cells, so
# check rather than assume.)

fitting in DGG.query(sys, DGG.Within(tile); level = DGG.level(fitting))

# ## One subtree as a grid
#
# `PartialGrid(sys, cell, level)` is that cell's subtree at a working level, as
# an ordinary `AbstractGrid`: positions run `1:ncells`, so a data vector indexes
# straight through it, and on a sorted-subtree system like IGEO7 the ids are a
# lazy window over the level grid, so building it is O(1).

leaf = 10                                          # ≈ 430 m cells
grid = DGG.PartialGrid(sys, fitting, leaf)
(; n = DGG.ncells(grid), first = DGG.cellindex(grid, 1),
   last = DGG.cellindex(grid, DGG.ncells(grid)))

# ## The DEM, conservatively regridded
#
# The tile is 3600×3600 at 30 m, which is far finer than the destination cells;
# `aggregate` averages it down to roughly one source pixel per DGGS cell before
# the intersection matrix is built.

# The download extent is a point inside the tile rather than the tile itself:
# `getraster` maps an extent to every 1° tile it touches, and a closed 1° box
# touches four of them.

Rasters.checkmem!(false)                           # the tile is bigger than free RAM

centre = Extents.Extent(X = (10.5, 10.5), Y = (46.5, 46.5))
path = only(skipmissing(RasterDataSources.getraster(CopernicusDEM; extent = centre)))
dem = aggregate(mean, Raster(path; lazy = true), 16)
size(dem)

# ConservativeRegridding wants the source's cell **corners** on the unit sphere.
# The destination is the DGGS grid, and it needs no adapter at all: `treeify`,
# `ncells` and `getcell` are `ConservativeRegridding.Trees`' own bindings,
# extended here for every `AbstractGrid`. The manifold is named explicitly
# because this package computes on the *unit* sphere.

(west, east), (south, north) = bounds(dem, X), bounds(dem, Y)
to_sphere = GO.UnitSpherical.UnitSphereFromGeographic()
corners = [to_sphere((lon, lat))
           for lon in range(west, east; length = size(dem, X) + 1),
               lat in range(south, north; length = size(dem, Y) + 1)]

manifold = GO.Spherical(; radius = 1.0)
regridder = CR.Regridder(manifold, grid, corners)

# The tile does not cover the whole subtree, so the standard coverage trick
# applies: regrid the field and a 0/1 indicator with the same matrix, then
# divide. The ratio is a weighted mean of the source values with the row of the
# intersection matrix as weights, which is what makes it right for a cell the
# tile only partly covers.

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
# The first step of every flow-routing model: each cell sends its water to the
# lowest of its neighbours. `halo_table` names them as positions in `elevation`
# in one call, and a cell with no lower neighbour is a pit. Adjacency on a
# subset is the complete level's clipped to membership, so a neighbour outside
# the subtree is simply absent from the row — and on a rooted subtree like this
# one the table takes the rim/interior split, which spares the interior any
# membership test at all. The DEM's own coverage is the second filter, and it is
# this page's rather than the grid's.

halo = [[p for p in row if covered[p]] for row in DGG.halo_table(grid)]

function downhill(i)
    isempty(halo[i]) && return 0
    j = halo[i][argmin(elevation[halo[i]])]
    return elevation[j] < elevation[i] ? j : 0
end

flow = [covered[i] ? downhill(i) : 0 for i in 1:DGG.ncells(grid)]
drop = [flow[i] == 0 ? NaN : elevation[i] - elevation[flow[i]] for i in eachindex(flow)]
(; n_covered = count(covered), n_pits = count(i -> covered[i] && flow[i] == 0, eachindex(flow)),
   max_drop = maximum(skipmissing(filter(!isnan, drop))))

# The two panels: elevation, and the drop to the downhill neighbour, which picks
# out the valley floors as the flat regions and the headwalls as the steep ones.

polys = [GO.transform(GO.GeographicFromUnitSphere(),
             DGG.cell_polygon(grid, DGG.cellindex(grid, i)))
         for i in 1:DGG.ncells(grid)]
shown = findall(covered)

fig = Figure(size = (900, 430))
ax1 = GeoAxis(fig[1, 1]; dest = "+proj=longlat +datum=WGS84", title = "elevation (m)")
p1 = poly!(ax1, polys[shown]; color = elevation[shown], colormap = :terrain, strokewidth = 0)
Colorbar(fig[2, 1], p1; vertical = false)
ax2 = GeoAxis(fig[1, 2]; dest = "+proj=longlat +datum=WGS84", title = "drop to downhill neighbour (m)")
p2 = poly!(ax2, polys[shown]; color = drop[shown], colormap = :magma,
    nan_color = :gray80, strokewidth = 0)
Colorbar(fig[2, 2], p2; vertical = false)
save("dem_igeo7.png", fig)
fig

# ## The same page on every system
#
# Only the singleton on the first line was IGEO7. `MultiOrderCoverage`,
# `PartialGrid` and `halo_table` are interface methods, so the coverage, the
# subtree grid and the halo table all come out the same way elsewhere — A5 alone
# has no `descendant_range`, so its coverage cannot be expanded to position
# ranges and its subtree ids are materialised instead of windowed.
#
# The `contained` column is the check from above: the shallowest emitted cell is
# a cell *inside* the tile only when the traversal was deep enough to emit one
# whole, which is why the maximum depth is chosen per system here.

for s in (DGG.systems()..., DGG.AuthalicSystem(DGG.IGeo7System()))
    base = s isa DGG.AuthalicSystem ? parent(s) : s
    l = base isa Union{DGG.IGeo7System, DGG.H3System} ? 8 : 10
    cov = DGG.query(s, DGG.MultiOrderCoverage(tile); level = l)
    top = argmin(DGG.level, cov)
    inside = top in DGG.query(s, DGG.Within(tile); level = DGG.level(top))
    g = DGG.PartialGrid(s, top, l)
    ranges = DGG.has_sorted_subtrees(s) ? length(DGG.level_ranges(cov, l)) : missing
    name = s isa DGG.AuthalicSystem ? "Authalic($(nameof(typeof(base))))" :
           string(nameof(typeof(s)))
    println(rpad(name, 24), "level $l: ", lpad(length(cov), 6), " coverage cells, ",
            "coarsest at level ", DGG.level(top), " (contained: ", inside, "), subtree ",
            lpad(DGG.ncells(g), 6), " cells, ranges ", ranges)
end
