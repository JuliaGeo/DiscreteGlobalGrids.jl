# # Stencil operations

import DiscreteGlobalGrids as DGG                                            # hide
import DimensionalData as DD                                                 # hide
import Extents                                                               # hide
using GLMakie, GeoMakie                                                      # hide
using DiscreteGlobalGridsVisualization: dggpoly!                             # hide
GLMakie.activate!(inline = true)                                             # hide
let                                                                          # hide
    grid = DGG.levelgrid(DGG.IGeo7System(), 4)                               # hide
    lookup = DGG.CellLookup(grid)                                            # hide
    role = DD.DimArray(zeros(length(lookup)), DGG.Cells(lookup))             # hide
    centre = DGG.cellat(grid, 10.0, 50.0)                                    # hide
    p = DGG.globalindex(grid, centre)                                        # hide
    for k in 1:3                                                             # hide
        parent(role)[DGG.ring(lookup, p, k)] .= 4 - k                        # hide
    end                                                                      # hide
    role[DGG.Cells(DD.At(centre))] = 4                                       # hide
    box = Extents.Extent(X = (-4.0, 24.0), Y = (41.0, 59.0))                 # hide
    sub = role[DGG.Cells(DGG.Covering(box))]                                 # hide
    shades = ["#dcf5d7", "#e5d6ec", "#c9a6db", "#8338b8", "#389826"]         # hide
    global fig = Figure(size = (900, 500))                                   # hide
    ax = GeoAxis(fig[1, 1]; dest = "+proj=laea +lat_0=50 +lon_0=10",         # hide
        limits = ((-4.0, 24.0), (41.0, 59.0)),                               # hide
        xgridvisible = false, ygridvisible = false,                          # hide
        xticklabelsvisible = false, yticklabelsvisible = false)              # hide
    dggpoly!(ax, sub; color = parent(sub), colorrange = (0, 4),              # hide
        colormap = shades, strokecolor = "#2c7a1e", strokewidth = 0.8)       # hide
    hi = role[sort(vcat(p, DGG.neighbors(lookup, p, 3)))]                    # hide
    dggpoly!(ax, hi; color = parent(hi), colorrange = (0, 4),                # hide
        colormap = shades, strokecolor = "#212529", strokewidth = 1.6)       # hide
    lines!(ax, GeoMakie.coastlines(); color = ("#212529", 0.7),              # hide
        linewidth = 0.7)                                                     # hide
    Legend(fig[1, 2],                                                        # hide
        [PolyElement(color = c, strokecolor = "#212529", strokewidth = 1)    # hide
         for c in reverse(shades)],                                          # hide
        ["the cell", "first-order neighbours", "second-order",               # hide
         "third-order", "the rest of the grid"]; framevisible = false)       # hide
end                                                                          # hide
fig                                                                          # hide

# A stencil recomputes every cell from its own value and its neighbours'. On a
# grid the neighbours are the cells that touch it — six on the IGeo7 hexagons
# above, eight on HEALPix — and a kernel that reads them once per cell is one
# `mapneighbors` call. Stacking passes reaches the second- and third-order
# neighbours in the picture.
#
# Every sweep on this page runs on a `DimArray` whose one dimension is `Cells`:
# a value per cell of a grid, indexed by the cells themselves.

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GeometryOps as GO
import Extents
using DataStructures: PriorityQueue, enqueue!, dequeue_pair!
using Statistics, Random
using DiscreteGlobalGridsVisualization: dggpoly, dggpoly!
using GLMakie, GeoMakie
GLMakie.activate!(inline = true)

# ## Building a field on the grid
#
# A HEALPix grid at level 5, and its cells as a lookup that a dimension can
# hold:

grid = DGG.levelgrid(DGG.HEALPixSystem(), 5)

#

lookup = DGG.CellLookup(grid)

# Each cell's centroid gives a longitude and a latitude to write the field in
# terms of. `cell_centroid` broadcasts over the lookup, and
# `GeographicFromUnitSphere` turns the unit-sphere points into lon/lat pairs.

lonlat = GO.UnitSpherical.GeographicFromUnitSphere()
centroids = lonlat.(DGG.cell_centroid.(grid, lookup))
lon = DD.DimArray(first.(centroids), DGG.Cells(lookup); name = :lon)
lat = DD.DimArray(last.(centroids), DGG.Cells(lookup); name = :lat)

# The field is a step at ±30°, a wave in longitude, and speckle on top. Any
# `DimArray` over the same `Cells` dimension replaces it — a [regridded
# raster](regridding.md), a column read out of [a store](store_io.md).

Random.seed!(42)
plateau = ifelse.(abs.(lat) .< 30, 8.0, 0.0)
wave = 3 .* sind.(2 .* lon) .* cosd.(lat)
field = DD.rebuild(plateau .+ wave .+ randn(length(lookup)); name = :value)

# `dggpoly` draws the cells of such an array as one mesh, coloured by its
# values.

crange = extrema(field)

#

fig = Figure(size = (820, 460))
ax = GeoAxis(fig[1, 1]; dest = "+proj=moll", title = "the field",
    xgridcolor = (:black, 0.15), ygridcolor = (:black, 0.15))
dggpoly!(ax, field; color = parent(field), colormap = :viridis,
    colorrange = crange)
lines!(ax, GeoMakie.coastlines(); color = (:black, 0.45), linewidth = 0.5)
Colorbar(fig[1, 2]; colormap = :viridis, colorrange = crange, label = "value")
fig

# ## Smoothing with mapneighbors
#
# `mapneighbors` calls a kernel once per cell, threaded. With
# `pass = Values()` the kernel is `f(cell, value, neighbours)`:
#
# - `value` — the cell's own entry;
# - `neighbours` — its neighbours' entries, counter-clockwise seen from outside
#   the sphere, in the order [`neighbors`](@ref) fixes.
#
# An oriented stencil (a gradient, an upwind scheme) relies on that order; an
# average ignores it. The result comes back as a `DimArray` over the same cells.

smooth(A) = DGG.mapneighbors((c, x, nbs) -> (x + sum(nbs)) / (1 + length(nbs)),
    A; pass = DGG.Values())

smoothed = smooth(field)

#

var(field), var(smoothed)

# Each pass mixes one hop further out, so iterating it is diffusion.

diffused = foldl((v, _) -> smooth(v), 1:10; init = field)
var(diffused)

# Close up, over Europe and North Africa. `Covering(box)` selects every cell
# over a lon/lat box, and the result keeps its cell lookup, so `dggpoly` draws
# it as readily as the whole grid.

box = Extents.Extent(X = (-30.0, 50.0), Y = (18.0, 68.0))
patch = field[DGG.Cells(DGG.Covering(box))]

#

zrange = extrema(patch)

#

fig = Figure(size = (900, 420))
for (k, (name, v)) in enumerate(("original" => field,
                                 "after 10 smoothing passes" => diffused))
    panel = GeoAxis(fig[1, k]; dest = "+proj=laea +lat_0=42 +lon_0=10",
        limits = ((-14.0, 34.0), (28.0, 58.0)), title = name,
        xgridvisible = false, ygridvisible = false,
        xticklabelsvisible = false, yticklabelsvisible = false)
    sub = v[DGG.Cells(DGG.Covering(box))]
    dggpoly!(panel, sub; color = parent(sub), colormap = :viridis,
        colorrange = zrange, strokecolor = (:white, 0.35), strokewidth = 0.3)
end
Colorbar(fig[1, 3]; colormap = :viridis, colorrange = zrange, label = "value")
fig

# ## Detecting edges with the Laplacian
#
# `mean(neighbours) - centre` is a discrete Laplacian, an edge detector.
# Smoothing first leaves the two step edges as the largest second differences
# in the field.

laplacian = DGG.mapneighbors((c, v, nbs) -> mean(nbs) - v, diffused;
    pass = DGG.Values())
extrema(laplacian)

# The colour range is the 95th percentile of |Laplacian|, so the field's
# texture stays visible; the two step edges exceed it and clip to `highclip`
# and `lowclip`.

q = quantile(abs.(laplacian), 0.95)

#

fig = Figure(size = (820, 460))
ax = GeoAxis(fig[1, 1]; dest = "+proj=moll",
    title = "Laplacian of the smoothed field",
    xgridcolor = (:black, 0.15), ygridcolor = (:black, 0.15))
dggpoly!(ax, laplacian; color = parent(laplacian), colormap = :balance,
    colorrange = (-q, q), highclip = :darkred, lowclip = :darkblue)
Colorbar(fig[1, 2]; colormap = :balance, colorrange = (-q, q),
    highclip = :darkred, lowclip = :darkblue, label = "mean(neighbours) - centre")
fig

# ## Reading more than one quantity per neighbour
#
# `needs` names what the kernel reads for each neighbour, and the sweep hands
# those quantities over directly. A finite-difference gradient reads two: the
# neighbour's value and its centroid, so each difference is taken over the real
# distance.

function steepest((value, centre), (values, centres))
    isempty(values) && return 0.0
    maximum(eachindex(values)) do k
        arc = GO.UnitSpherical.spherical_distance(centres[k], centre)
        abs(values[k] - value) / (6371.0 * arc)   # 6371 km: the Earth's radius
    end
end

gradient = DGG.mapneighbors(steepest, field;
    needs = (DGG.Value(field), DGG.Centroid()))
extrema(gradient)   # units per kilometre

# The callback is `f(center, rings)`; both are tuples with one entry per need,
# in `needs` order:
#
# | | `center` | `rings` |
# |---|---|---|
# | need 1, `Value(field)` | the cell's value | every neighbour's value |
# | need 2, `Centroid()` | the cell's centroid | every neighbour's centroid |
#
# Slot `k` of each ring is the same neighbour. Here is `rings[1]` for one cell,
# its eight neighbours' values:

valuerings = DGG.mapneighbors((centre, rings) -> first(rings), field;
    needs = (DGG.Value(field), DGG.Centroid()))
valuerings[1]

# A kernel that wants one record per neighbour writes `zip(rings...)` itself.
# [Requesting neighbour fields](../api/neighbor-fields.md) is the full account,
# including `Cell()` and `Index(Local())` and what each costs.
#
# ## Materialising the neighbour table with adjacency
#
# `adjacency` materialises every cell's neighbour list once, as one CSR table,
# for kernels that make many passes over the same grid. `Cells` names the
# dimension to walk; the sweeps above found it themselves.

table = DGG.adjacency(field, DGG.Cells)

# Row `i` is `neighbors(lookup, i)` in the same order the sweep uses, so the
# table gives the same answer the sweep does:

[mean(vcat(field[i], field[table[i]])) for i in eachindex(field)] ≈ smoothed

# The rows differ in length where the grid does — eight neighbours nearly
# everywhere, seven at the cells sitting on a degree-three vertex of the base
# tiling:

sort(unique(length.(table))), count(==(7), length.(table))

# `Vertex()`, the default connectivity, counts every cell that shares a corner;
# `Edge()` counts only the cells that share a side. Both are drawn around the
# cell over Zürich, which `Contains` resolves to one index along `Cells`:

p = only(DD.dims2indices(field, DGG.Cells(DD.Contains((8.5, 47.4)))))

#

disk = sort(vcat(p, DGG.neighbors(lookup, p, 3)))
shades = ["#dcf5d7", "#8338b8", "#389826"]

fig = Figure(size = (820, 420))
for (k, conn) in enumerate((DGG.Vertex(), DGG.Edge()))
    nbrs = DGG.neighbors(lookup, p; connectivity = conn)
    role = DD.DimArray(zeros(length(lookup)), DGG.Cells(lookup))
    parent(role)[nbrs] .= 1
    parent(role)[p] = 2
    panel = GeoAxis(fig[1, k]; dest = "+proj=laea +lat_0=47.4 +lon_0=8.5",
        title = "$(nameof(typeof(conn)))(): $(length(nbrs)) neighbours",
        xgridvisible = false, ygridvisible = false,
        xticklabelsvisible = false, yticklabelsvisible = false)
    around = role[disk]
    dggpoly!(panel, around; color = parent(around), colorrange = (0, 2),
        colormap = shades, strokecolor = "#2c7a1e", strokewidth = 1)
    hi = role[sort(vcat(p, nbrs))]
    dggpoly!(panel, hi; color = parent(hi), colorrange = (0, 2),
        colormap = shades, strokecolor = "#212529", strokewidth = 2.5)
end
fig

# ## Running a stencil on a subset of the grid
#
# A subset of the grid clips every neighbourhood to its membership: a cell on
# the subset's edge has fewer neighbours, and the same `smooth` call runs on
# it.

smoothed_patch = smooth(patch)
count(.!(parent(smoothed_patch) .≈ parent(smoothed[DGG.Cells(DGG.Covering(box))])))

# Those are the cells on the edge of `patch`, whose neighbourhoods lost
# members. [Out of core](out_of_core.md) loads the margin a stencil on a tile
# needs, so its edge cells agree with the whole-grid pass.
#
# ## Cost distance with Dijkstra
#
# Given a cost per cell, the cost distance to a cell is the smallest
# accumulated cost of any path from the seed. Dijkstra's algorithm computes it
# with a priority queue; `neighbors(lookup, i)` — a local index in, local
# indices out — is the only grid call it needs.

function costdistance(cost, seed)
    lookup = DD.lookup(cost, DGG.Cells)
    source = only(DD.dims2indices(cost, DGG.Cells(DD.Contains(seed))))
    dist = fill(Inf, length(cost))
    dist[source] = 0.0
    queue = PriorityQueue{Int,Float64}()
    enqueue!(queue, source => 0.0)
    while !isempty(queue)
        i, d = dequeue_pair!(queue)
        d > dist[i] && continue          # a stale entry, already improved
        for j in DGG.neighbors(lookup, i)
            alt = d + (cost[i] + cost[j]) / 2
            if alt < dist[j]
                dist[j] = alt
                queue[j] = alt           # insert, or lower the key
            end
        end
    end
    return DD.rebuild(cost; data = dist, name = :costdistance)
end

# Crossing a cell costs 1 everywhere except in a wall along the equator, which
# is impassable — cost `Inf` — apart from a gap between 0° and 25° E. The seed
# is Zürich. Every cell's cost distance is then what it takes to reach it
# around the wall, or through the gap.

wall = @. (abs(lat) < 5) & !(0 <= lon <= 25)
cost = DD.rebuild(lat; data = ifelse.(wall, Inf, 1.0), name = :cost)
travel = costdistance(cost, (8.5, 47.4))
reach = maximum(filter(isfinite, parent(travel)))

# The projection is azimuthal equidistant, centred on the seed: distance from
# the centre of the page is great-circle distance from the seed, so at uniform
# cost the bands are circles, and every departure from a circle is the wall's
# doing. The map covers a 172° cap around the seed, selected with
# `Covering(cap)`.

centre = DGG.cell_centroid(grid, DGG.cellat(grid, 8.5, 47.4))
cap = GO.UnitSpherical.SphericalCap(centre, deg2rad(172))
near = travel[DGG.Cells(DGG.Covering(cap))]

#

bands = cgrad(:viridis, 12; categorical = true)

fig = Figure(size = (960, 460))
ax = GeoAxis(fig[1, 1]; title = "cost of crossing a cell",
    dest = "+proj=aeqd +lat_0=47.4 +lon_0=8.5", xgridvisible = false,
    ygridvisible = false, xticklabelsvisible = false,
    yticklabelsvisible = false)
dggpoly!(ax, near; color = ifelse.(isinf.(parent(near)), 2.0, 1.0),
    colorrange = (1, 2), colormap = ["#eef3f7", "#b03030"])
lines!(ax, GeoMakie.coastlines(); color = (:black, 0.4), linewidth = 0.4)
scatter!(ax, [8.5], [47.4]; color = :dodgerblue, marker = :star5,
    markersize = 16, strokecolor = :black, strokewidth = 0.5)
ax = GeoAxis(fig[1, 2]; title = "cost distance from the seed",
    dest = "+proj=aeqd +lat_0=47.4 +lon_0=8.5", xgridvisible = false,
    ygridvisible = false, xticklabelsvisible = false,
    yticklabelsvisible = false)
dggpoly!(ax, near; color = parent(near), colormap = bands,
    colorrange = (0, reach), highclip = "#8c1b1b")
lines!(ax, GeoMakie.coastlines(); color = (:black, 0.4), linewidth = 0.4)
scatter!(ax, [8.5], [47.4]; color = :red, marker = :star5, markersize = 16,
    strokecolor = :black, strokewidth = 0.5)
Colorbar(fig[1, 3]; colormap = bands, colorrange = (0, reach),
    highclip = "#8c1b1b", label = "cost distance")
fig

# On the seed's side of the wall the bands are circles around the seed. Every
# cell beyond it is reached through the gap, so the bands there are centred on
# the gap. The wall itself is `Inf` and draws in dark red, above the top of the
# scale.
#
# ## Sweeping a plain vector of values
#
# The same sweeps take the cells and the values as two arguments, for code that
# holds them apart:

cells = DGG.CellVector(grid)
values = parent(field)
DGG.mapneighbors((c, x, nbs) -> (x + sum(nbs)) / (1 + length(nbs)), cells,
    values) ≈ parent(smoothed)

# ## Running the page on another grid system
#
# `levelgrid`, `mapneighbors`, `adjacency` and `neighbors` are interface
# methods, so the page runs unchanged on any system; only the neighbour count
# changes — six on IGeo7's hexagons, eight here. [Choosing a
# grid](choosing_a_grid.md) covers what that does to a kernel.
#
# The smoothing pass is also the graph convolution behind DeepSphere-style
# machine learning on spherical grids: stack `k` of them and every cell has a
# `k`-hop receptive field.
