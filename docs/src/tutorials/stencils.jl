# # Stencil operations
#
# A stencil operation recomputes every cell of a grid from its own value and the
# values of its immediate neighbours. Two calls set it up on any grid in this
# package: `neighbors` names the neighbours, and `cellposition` turns each name
# into an index into the data vector. That halo table is the stencil.
#
# This is also the primitive behind machine learning on spherical grids: a graph
# convolution in the style of DeepSphere is exactly such a pass, and stacking `k`
# of them gives each cell a `k`-hop receptive field.

import DiscreteGlobalGrids as DGG
import GeometryOps as GO
using Statistics, Random
using CairoMakie, GeoMakie
CairoMakie.activate!()

# ## A field on the whole globe
#
# HEALPix level 5 is `nside = 32`, 12288 cells. `levelgrid` is the complete
# level; it stores the system and the level and nothing else, so building it is
# free.

grid = DGG.levelgrid(DGG.HEALPixSystem(), 5)
cells = [DGG.cellindex(grid, i) for i in 1:DGG.ncells(grid)]

# Geometry is on the unit sphere throughout. `GeographicFromUnitSphere` is the
# one conversion this page needs, for sampling a lon/lat field at the cell
# centroids and, later, for plotting.

lonlat = GO.UnitSpherical.GeographicFromUnitSphere()
centers = [lonlat(DGG.cell_centroid(grid, c)) for c in cells]

# The test field is a hot zonal band with a wavy longitude signal on top, plus
# noise — sharp edges at ±30° latitude for the Laplacian to find, and speckle
# for smoothing to remove.

Random.seed!(42)
field(lon, lat) = (abs(lat) < 30 ? 8.0 : 0.0) + 3 * sind(2lon) * cosd(lat)
values = [field(c...) for c in centers] .+ 0.4 .* randn(length(cells))

# ## The halo table
#
# `neighbors(grid, c)` returns the ring-1 neighbours as typed cell ids, counter-
# clockwise seen from outside the sphere. `cellposition` maps each back to its
# position in the grid's dense order, which is the index into `values`. On a
# complete level every neighbour is present; on a `PartialGrid` the ones outside
# the coverage are simply absent, and the table is shorter there.

halo = [Int[DGG.cellposition(grid, nb) for nb in DGG.neighbors(grid, c)] for c in cells]
length.(halo) |> extrema

# HEALPix has 8 neighbours nearly everywhere; 24 cells per grid sit at a
# degree-3 vertex of the base tiling and have only 7.

count(==(7), length.(halo))

# ## Smoothing, edges, and diffusion
#
# With the table in hand a stencil is one comprehension: `f(centre, neighbours)`
# over every position. Averaging the neighbourhood smooths; the discrete
# Laplacian `mean(neighbours) - centre` is an edge detector.

stencil(f, v) = [f(v[i], v[halo[i]]) for i in eachindex(v)]

smoothed = stencil((c, nbs) -> mean(vcat(c, nbs)), values)
laplacian = stencil((c, nbs) -> mean(nbs) - c, values)
(var(values), var(smoothed))

# Each pass mixes a cell with its neighbours one hop out, so `k` repeated passes
# give a `k`-hop receptive field — iterated smoothing is diffusion.

diffused = foldl((v, _) -> stencil((c, nbs) -> mean(vcat(c, nbs)), v), 1:10; init = values)
var(diffused)

# ## Plotting
#
# `cell_polygon` returns a unit-sphere polygon; one `GO.transform` takes the
# whole vector to lon/lat. The `+over` flag keeps PROJ from wrapping the cells
# that straddle ±180°, a hairline stroke in the fill colour hides antialiasing
# seams, and the Laplacian panel gets a diverging colormap centred at 0 so the
# band edges at ±30° stand out against the noise.

polys = GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygon.(Ref(grid), cells))
crange = extrema(values)
lmax = maximum(abs, laplacian)

fig = Figure(size = (650, 1000))
for (row, (title, v, kw)) in enumerate((
        ("original", values, (; colorrange = crange)),
        ("after 10 smoothing passes", diffused, (; colorrange = crange)),
        ("Laplacian of the original", laplacian,
            (; colormap = :balance, colorrange = (-lmax, lmax)))))
    ax = GeoAxis(fig[row, 1]; dest = "+proj=moll +over", title)
    poly!(ax, polys; color = v, strokecolor = v, strokewidth = 0.7, kw...)
end
fig

# ## The same three lines on any system
#
# Nothing above named HEALPix except the singleton. `neighbors` and
# `cellposition` are interface methods, so the halo table is built the same way
# on every registered system, and on an `AuthalicSystem` wrap of one — only the
# degree changes.

for sys in (DGG.systems()..., DGG.AuthalicSystem(DGG.IGeo7System()))
    base = sys isa DGG.AuthalicSystem ? parent(sys) : sys
    name = sys isa DGG.AuthalicSystem ? "Authalic($(nameof(typeof(base))))" :
           string(nameof(typeof(sys)))
    l = base isa Union{DGG.IGeo7System, DGG.H3System} ? 3 : 4
    g = DGG.levelgrid(sys, l)
    degrees = [length(DGG.neighbors(g, DGG.cellindex(g, i))) for i in 1:DGG.ncells(g)]
    edges = [length(DGG.neighbors(g, DGG.cellindex(g, i); connectivity = DGG.Edge()))
             for i in 1:DGG.ncells(g)]
    println(rpad(name, 18), " level $l: ", DGG.ncells(g),
            " cells, vertex degree ", extrema(degrees), ", edge degree ", extrema(edges))
end
