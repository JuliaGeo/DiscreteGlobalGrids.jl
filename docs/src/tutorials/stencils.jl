# # Stencil operations
#
# A stencil operation recomputes every cell of a grid from its own value and the
# values of its immediate neighbors. On HEALPix nearly every cell has the same
# 8-neighbor structure, so a whole-globe pass is one sweep over a precomputed
# neighbor table. This is also the primitive behind machine learning on
# spherical grids: a graph convolution in the style of DeepSphere is exactly
# such a pass, and stacking `k` of them gives each cell a `k`-hop receptive
# field.
#
# ## A field on the whole globe
#
# We work at level 5 (`nside = 32`, 12288 cells). `DGGSGlobeIds` names every
# cell at the level without materializing them, and the lookup built from it
# behaves like any other `HealpixLookup`.

using DiscreteGlobalGrids
using DiscreteGlobalGrids.HEALPix.HealpixLookups
import DimensionalData as DD
using Statistics
using Random
using CairoMakie, GeoMakie
CairoMakie.activate!()

l = HealpixLookup(DGGSGlobeIds(HEALPixDGGS(), 5))
nothing

# The test field is a hot zonal band with a wavy longitude signal on top, plus
# noise — sharp edges at ±30° latitude for the Laplacian to find, and speckle
# for smoothing to remove.

Random.seed!(42)
field(lon, lat) = (abs(lat) < 30 ? 8.0 : 0.0) + 3 * sind(2lon) * cosd(lat)
vals = [field(lon, lat) for (lon, lat) in cell_centers(l)] .+ 0.4 .* randn(length(l))
A = DD.DimArray(vals, Cells(l); name = :field)

# ## The neighbor table
#
# `neighbor_indices` builds the halo table once: for each cell, the positions
# (into the lookup) of its 8 HEALPix neighbors, in the order SW, W, NW, N, NE,
# E, SE, S. 24 pixels per grid sit at a degree-3 vertex of the base tiling and
# have only 7 — the missing slot holds 0.

nbi = HealpixLookups.neighbor_indices(l)
nbi[1]

# ## Smoothing, edges, and diffusion
#
# `stencil(f, A; nbidx)` applies `f(center, neighbors::Vector)` to every cell
# and returns a new array over the same lookup. Averaging the neighborhood
# smooths; the discrete Laplacian `mean(neighbors) - center` is an edge
# detector.

smoothed = stencil((c, nbs) -> mean(vcat(c, nbs)), A; nbidx = nbi)
lap = stencil((c, nbs) -> mean(nbs) - c, A; nbidx = nbi)
var(parent(A)), var(parent(smoothed))

# Each pass mixes a cell with its neighbors one hop out, so `k` repeated passes
# give a `k`-hop receptive field — iterated smoothing is diffusion.

smooth_pass(X) = stencil((c, nbs) -> mean(vcat(c, nbs)), X; nbidx = nbi)

function diffuse(X, passes)
    for _ in 1:passes
        X = smooth_pass(X)
    end
    return X
end

diffused = diffuse(A, 10)
nothing

# `cell_polygons` returns one lon/lat polygon per stored cell, ready to hand to
# `poly!`. The `+over` flag keeps PROJ from wrapping the cells that straddle
# ±180°, a hairline stroke in the fill color hides antialiasing seams between
# cells, and the Laplacian panel gets a diverging colormap centered at 0, so
# the band edges at ±30° stand out against the noise.

polys = cell_polygons(l)
crange = extrema(vals)
lmax = maximum(abs, parent(lap))

fig = Figure(size = (650, 1000))
ax1 = GeoAxis(fig[1, 1]; dest = "+proj=moll +over", title = "original")
poly!(ax1, polys; color = collect(parent(A)), colorrange = crange,
    strokecolor = collect(parent(A)), strokewidth = 0.7)
ax2 = GeoAxis(fig[2, 1]; dest = "+proj=moll +over", title = "after 10 smoothing passes")
poly!(ax2, polys; color = collect(parent(diffused)), colorrange = crange,
    strokecolor = collect(parent(diffused)), strokewidth = 0.7)
ax3 = GeoAxis(fig[3, 1]; dest = "+proj=moll +over", title = "Laplacian of the original")
poly!(ax3, polys; color = collect(parent(lap)), colormap = :balance,
    colorrange = (-lmax, lmax), strokecolor = collect(parent(lap)),
    strokewidth = 0.7)
fig
