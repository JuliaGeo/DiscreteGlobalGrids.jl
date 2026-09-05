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

# A stencil computes a new value from a cell and its neighbours. The figure
# shows successive rings of neighbours around one cell. Here we use stencils
# to smooth a field, detect its edges, request extra neighbour data, and
# compute a graph distance. Each example uses `mapneighbors` or the related
# `adjacency` table, so the same code can work with a different grid system.
#
# Every field here is a `DimArray` with one `Cells` dimension. That dimension
# holds the lookup that connects array positions to cells and their neighbours.

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

# We use cell centroids to define a reproducible test field. `cell_centroid`
# broadcasts over the lookup, and `GeographicFromUnitSphere` converts the
# returned unit-sphere points to longitude and latitude.

lonlat = GO.UnitSpherical.GeographicFromUnitSphere()
centroids = lonlat.(DGG.cell_centroid.(grid, lookup))
lon = DD.DimArray(first.(centroids), DGG.Cells(lookup); name = :lon)
lat = DD.DimArray(last.(centroids), DGG.Cells(lookup); name = :lat)

# The field combines a broad step, a longitude wave, and random noise. Any
# `DimArray` over the same `Cells` dimension can take its place, including a
# [regridded raster](regridding.md) or a column read from [a store](store_io.md).

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
# `mapneighbors` applies a kernel once per cell and can thread those calls.
# With `pass = Values()`, the kernel receives `f(cell, value, neighbours)`:
#
# - `value` — the cell's own entry;
# - `neighbours` — its neighbours' entries, counter-clockwise seen from outside
#   the sphere, in the order [`neighbors`](@ref) fixes.
#
# The neighbour order matters for oriented stencils such as gradients and
# upwind schemes. An average only needs the values. The result keeps the same
# `Cells` dimension.

smooth(A) = DGG.mapneighbors((c, x, nbs) -> (x + sum(nbs)) / (1 + length(nbs)),
    A; pass = DGG.Values())

smoothed = smooth(field)

#

var(field), var(smoothed)

# Repeating the pass increases the radius of the operation and produces a
# simple diffusion process.

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
# The discrete Laplacian compares a cell with the mean of its neighbours.
# Applying it after smoothing makes the step boundaries stand out as large
# second differences.

laplacian = DGG.mapneighbors((c, v, nbs) -> mean(nbs) - v, diffused;
    pass = DGG.Values())
extrema(laplacian)

# The colour range uses the 95th percentile of the absolute Laplacian. This
# keeps ordinary texture visible while larger edge values use the clip colours.

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
# `needs` declares the fields a kernel reads for each neighbour. This example
# requests values and centroids so it can measure a finite difference per unit
# great-circle distance.

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

# A kernel that wants one record per neighbour can combine the parallel rings
# with `zip(rings...)`. [Requesting neighbour fields](../api/neighbor-fields.md)
# documents the other available requests, including `Cell()` and
# `Index(Local())`.
#
# ## Materialising the neighbour table with adjacency
#
# `adjacency` materialises the neighbour lists once as a CSR table. Reuse it
# when many passes will traverse the same grid. `Cells` names the dimension to
# walk; `mapneighbors` inferred it in the earlier examples.

table = DGG.adjacency(field, DGG.Cells)

# Row `i` is `neighbors(lookup, i)` in the same order the sweep uses, so the
# table gives the same answer the sweep does:

[mean(vcat(field[i], field[table[i]])) for i in eachindex(field)] ≈ smoothed

# Rows can have different lengths because the grid has cells with different
# degrees. HEALPix has fewer neighbours at cells around a base-tiling vertex:

sort(unique(length.(table))), count(==(7), length.(table))

# `Vertex()` counts cells that share a corner, while `Edge()` keeps only cells
# that share a side. The next figure compares both neighbourhoods around the
# cell containing Zürich, selected with `Contains`:

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
# A subset is itself a valid input. Its boundary cells see only neighbours that
# belong to the subset, so the same `smooth` call runs with a clipped stencil.

smoothed_patch = smooth(patch)
count(.!(parent(smoothed_patch) .≈ parent(smoothed[DGG.Cells(DGG.Covering(box))])))

# The count identifies cells whose subset neighbourhood differs from the
# whole-grid result. [Out of core](out_of_core.md) shows how to load the margin
# a tile needs when boundary cells must agree with a full-grid pass.
#
# ## Cost distance with Dijkstra
#
# A cost field turns the cell graph into a travel network. The cost distance is
# the least accumulated cost from a seed, which Dijkstra's algorithm computes
# with a priority queue. `neighbors(lookup, i)` supplies the graph edges as
# local indices. This example measures cost per graph step; a travel cost per
# kilometre would also need the distance between cell centres.

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

# The example makes an equatorial wall impassable (`Inf`) except for a gap
# between 0° and 25° E. With Zürich as the seed, the resulting distance field
# shows how paths go around the wall or through the gap.

wall = @. (abs(lat) < 5) & !(0 <= lon <= 25)
cost = DD.rebuild(lat; data = ifelse.(wall, Inf, 1.0), name = :cost)
travel = costdistance(cost, (8.5, 47.4))
reach = maximum(filter(isfinite, parent(travel)))

# The azimuthal-equidistant projection centres the map on the seed, making
# great-circle distance easy to compare visually. The map uses a 172° cap
# selected with `Covering(cap)`.

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

# The bands bend around the wall and spread from the gap on its far side. The
# `Inf` wall is drawn in dark red above the distance scale.
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
# `levelgrid`, `mapneighbors`, `adjacency`, and `neighbors` share the same
# interface across grid systems. Switching the constructor changes the cell
# topology while the stencil code stays the same. [Choosing a
# grid](choosing_a_grid.md) compares those topology choices.
