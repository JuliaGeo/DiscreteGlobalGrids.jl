# # Stencil operations
#
# A stencil operation recomputes every cell from its own value and its
# neighbours'. On any grid in this package that is three lines: `adjacency`
# gives, for every position, the positions of the cells it touches, and a
# stencil is a comprehension over its rows. The same pass is the graph
# convolution behind DeepSphere-style machine learning on spherical grids —
# stacking `k` passes gives each cell a `k`-hop receptive field.

import DiscreteGlobalGrids as DGG
import GeometryOps as GO
using Statistics, Random
using GLMakie, GeoMakie
GLMakie.activate!(inline = true)

# ## A field on the sphere
#
# HEALPix level 5 is 12288 cells. The test field is a hot band inside ±30°
# latitude with a wavy longitude signal on top, plus noise — sharp edges for a
# Laplacian to find, speckle for smoothing to remove. Geometry lives on the
# unit sphere, so sampling at the cell centroids takes one conversion to
# lon/lat; `CellVector(grid)` is the grid's cells as a lazy vector.

grid = DGG.levelgrid(DGG.HEALPixSystem(), 5)
cells = DGG.CellVector(grid)
lonlat = GO.UnitSpherical.GeographicFromUnitSphere()

Random.seed!(42)
field(lon, lat) = (abs(lat) < 30 ? 8.0 : 0.0) + 3 * sind(2lon) * cosd(lat)
values = [field(lonlat(DGG.cell_centroid(grid, c))...) for c in cells] .+
    0.4 .* randn(length(cells))

# ## Smoothing in three lines
#
# Row `p` of the table is `neighbors(grid, p)`: the in-set positions of `p`'s
# neighbours, counter-clockwise seen from outside the sphere, starting from the
# same neighbour every time the cell is asked. Every neighbour idiom in the
# package uses that one order, so on a complete grid like this one an oriented
# stencil — a gradient, an upwind scheme — can be written against the rows. On
# a subset, out-of-set members drop from a row without padding: the rotation
# stays the complete ring's, only the slots close up, so slot `j` alone no
# longer names a fixed direction. Averaging ignores the order and just smooths.

table = DGG.adjacency(grid)
stencil(f, v) = [f(v[i], v[table[i]]) for i in eachindex(v)]
smoothed = stencil((c, nbs) -> mean(vcat(c, nbs)), values)
(var(values), var(smoothed))

# The rows differ in length where the grid does: 8 neighbours nearly
# everywhere, 7 at the 24 cells that sit on a degree-3 vertex of the base
# tiling. Eight because adjacency is `connectivity = Vertex()` by default —
# touching at a corner counts — where `Edge()` would keep only the four cells
# that share a side.

extrema(length.(table)), count(==(7), length.(table))

# Two more stencils for free: the discrete Laplacian `mean(neighbours) - centre`
# is an edge detector, and each smoothing pass mixes one hop further out, so
# iterating it is diffusion.

laplacian = stencil((c, nbs) -> mean(nbs) - c, values)
diffused = foldl((v, _) -> stencil((c, nbs) -> mean(vcat(c, nbs)), v), 1:10; init = values)
var(diffused)

# The same pass without the intermediate table is `mapneighbors(f, cells)`,
# which calls `f` per cell and threads by default (`threaded = true`); its
# side-effect twin `foreachneighbors` defaults to `threaded = false`, because
# what it writes into is the caller's to make safe.
#
# ## Plotting
#
# `poly!` draws `cells` directly, as lon/lat polygons in position order.
# `+over` keeps PROJ from rewrapping the cells that straddle ±180°, a hairline
# stroke in the fill colour hides antialiasing seams, and the Laplacian gets a
# diverging colormap centred at 0.

crange = extrema(values)
lmax = maximum(abs, laplacian)

fig = Figure(size = (650, 1000))
for (row, (title, v, kw)) in enumerate((
        ("original", values, (; colorrange = crange)),
        ("after 10 smoothing passes", diffused, (; colorrange = crange)),
        ("Laplacian of the original", laplacian,
            (; colormap = :balance, colorrange = (-lmax, lmax)))))
    ax = GeoAxis(fig[row, 1]; dest = "+proj=moll +over", title)
    poly!(ax, cells; color = v, strokecolor = v, strokewidth = 0.7, kw...)
end
fig

# ## The same pass on a region
#
# A `PartialGrid` is a subset of one level, and adjacency on it is the complete
# level's clipped to membership: a neighbour outside the subset is omitted, not
# padded. The table is the same call; the rows on the subset's edge are simply
# shorter.

sys = DGG.system(grid)
face = DGG.cellindex(DGG.levelgrid(sys, 0), 5)     # one of the twelve base cells
sub = DGG.subtree(sys, face, 5)                # its level-5 subtree, 32 × 32
subtable = DGG.adjacency(sub)
(; n = DGG.ncells(sub), edge = count(<(8), length.(subtable)))

# Position `i` of the subset is position `i` of its data vector, so the stencil
# is the same comprehension.

subcells = DGG.CellVector(sub)
subvalues = [field(lonlat(DGG.cell_centroid(sub, c))...) for c in subcells]
substencil(f, v) = [f(v[i], v[subtable[i]]) for i in eachindex(v)]
subsmoothed = substencil((c, nbs) -> mean(vcat(c, nbs)), subvalues)
(var(subvalues), var(subsmoothed))

# The clipping is a law, not a convention. Ring 2 of a cell on the subset's
# edge is ring 2 of the complete level with the absent cells dropped — not a
# breadth-first walk of the subset, which would detour around what is missing
# and name different cells. (A cell the subset does not hold has no
# neighbourhood here at all: asking is an `ArgumentError`.)

edgecell = subcells[findfirst(<(8), length.(subtable))]
DGG.ring(sub, edgecell, 2) == filter(in(sub), DGG.ring(grid, edgecell, 2))

# What the clipping dropped has a name: `halo` is the cells just outside the
# subset that touch it — the extra fetch list a stencil on a tile needs, and
# the other half of `adjacency`'s answer. It is lazy, so ask for as much of it
# as you want.

DGG.ncells(sub), length(collect(DGG.halo(sub)))

# Nothing above named HEALPix except the singleton. `levelgrid`, `adjacency`
# and the position forms of `neighbors` and `ring` are interface methods, so
# the same three lines run unchanged on every registered system — only the
# degree changes.
