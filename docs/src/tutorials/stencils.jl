# # Stencil operations
#
# A stencil operation recomputes every cell of a grid from its own value and the
# values of its immediate neighbours. One call sets it up on any grid in this
# package: `halo_table` gives, for every position, the positions of the cells
# within `k` steps of it. That table is the stencil.
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
# `halo_table(grid, k)` is the whole stencil: entry `p` is the positions of the
# cells within `k` steps of position `p`, ascending. It is exactly
# `[DGG.neighbors(grid, p, k) for p in 1:DGG.ncells(grid)]` — the same answer,
# built in one pass — and `neighbors(grid, c)` still names the same cells as
# typed ids, counter-clockwise seen from outside the sphere, when the identities
# rather than the indices are what is wanted.

halo = DGG.halo_table(grid)
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

# ## The same pass on a region
#
# A `PartialGrid` is a subset of one level, and adjacency on it is the complete
# level's **clipped to membership**: a neighbour outside the subset is omitted,
# not padded, and distance is still measured in the system rather than inside
# the subset. So the table is the same call, and the rows on the subset's edge
# are simply shorter.

sys = DGG.HEALPixSystem()
face = DGG.cellindex(DGG.levelgrid(sys, 0), 5)     # one of the twelve base cells
sub = DGG.PartialGrid(sys, face, 5)                # its level-5 subtree, 32 x 32
subcells = [DGG.cellindex(sub, i) for i in 1:DGG.ncells(sub)]
subhalo = DGG.halo_table(sub)
(; n = DGG.ncells(sub), interior = count(==(8), length.(subhalo)),
   edge = count(<(8), length.(subhalo)))

# Position `i` of the subset is position `i` of its data vector, so the stencil
# is the same comprehension as before.

subvalues = [field(lonlat(DGG.cell_centroid(sub, c))...) for c in subcells]
substencil(f, v) = [f(v[i], v[subhalo[i]]) for i in eachindex(v)]
subsmoothed = substencil((c, nbs) -> mean(vcat(c, nbs)), subvalues)
(var(subvalues), var(subsmoothed))

# The clipping is a law, not a convention. Ring 2 of a cell on the subset's edge
# is ring 2 of the complete level with the absent cells dropped — not a
# breadth-first walk of the subset, which would have to detour around what is
# missing and would name different cells:

edgecell = subcells[findfirst(<(8), length.(subhalo))]
DGG.ring(sub, edgecell, 2) ==
    [c for c in DGG.ring(grid, edgecell, 2) if DGG.cellposition(sub, c) !== nothing]

# A cell the subset does not hold has no neighbourhood here at all — asking is
# an `ArgumentError` rather than the complete level's answer about cells that
# are not in the data.

# ## The same three lines on any system
#
# Nothing above named HEALPix except the singleton. `neighbors`, `halo_table`
# and `cellposition` are interface methods, so the halo table is built the same
# way on every registered system, and on an `AuthalicSystem` wrap of one — only
# the degree changes.

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
