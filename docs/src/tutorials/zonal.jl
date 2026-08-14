# # Zonal statistics
#
# Zonal statistics is rasterization's sibling: take a polygon, find the cells
# that cover it, aggregate data over them. Here the polygon is Texas, the grid
# is HEALPix at level 7, and the statistic is a mean.

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GeometryOps as GO
import NaturalEarth
using Statistics
using GLMakie, GeoMakie
GLMakie.activate!(inline = true)

# ## A polygon, a grid, and data
#
# Texas comes from Natural Earth. The grid is every HEALPix cell at level 7,
# read as a `CellVector` — the level's cell ids as a vector, so that position
# `k` in a data array means cell `cells[k]`. The data is a synthetic field
# sampled at each cell's centroid; only the sampling needs lon/lat.

fc = NaturalEarth.naturalearth("admin_1_states_provinces", 50)
texas = fc.geometry[findfirst(==("Texas"), fc.name)]

grid = DGG.levelgrid(DGG.HEALPixSystem(), 7)
cells = DGG.CellVector(grid)
#
f, a, p = poly(cells; axis = (; aspect = DataAspect()))
lines!(a, texas)
f
#
lonlat_tf = GO.UnitSpherical.GeographicFromUnitSphere()
field_f(lon, lat) = 20 - 0.5 * (lat - 25) + 2 * sind(3 * lon)
data = [field_f(lonlat_tf(DGG.cell_centroid(grid, c))...) for c in cells]

# ## The zonal mean
#
# `covering_positions` names the positions of every cell covering Texas,
# ready to index `data`. HEALPix cells are equal-area, so this unweighted
# mean *is* the areal mean — no latitude weights.

tx = DGG.covering_positions(cells, texas)
(; n = length(tx), tavg = mean(data[tx]))

# ## Touching or inside
#
# The covering keeps every cell that touches the outline. The opposite rule
# keeps only the cells wholly inside it — `query` with the `Within`
# predicate. Any boundary rule (a raster tool's centre-in-zone, say) lands
# between the two; at this resolution the bracket is a hundredth wide.

interior = DGG.query(grid, DGG.Within(texas))
(; touching = mean(data[tx]),
   inside = mean(data[DGG.cellposition(grid, c)] for c in interior))

# ## The picture
#
# The covering, coloured by the data — the polygon, rasterized onto the grid.

polys = GO.transform(lonlat, [DGG.cell_polygon(grid, cells[k]) for k in tx])

fig = Figure(size = (700, 620))
ax = GeoAxis(fig[1, 1]; dest = "+proj=longlat +datum=WGS84",
    limits = ((-107.5, -92.5), (25.0, 37.5)),
    title = "HEALPix level-7 cells covering Texas")
poly!(ax, polys; color = data[tx], colormap = :thermal,
    strokecolor = :white, strokewidth = 0.1)
poly!(ax, texas; color = :transparent, strokecolor = :black, strokewidth = 2)
fig

# ## The cube spelling
#
# The same selection inside a DimensionalData cube: `Cells` is the dimension,
# and `Covering` is `covering` wearing a selector hat.

A = DD.DimArray(data, DGG.Cells(DGG.CellLookup(cells)); name = :tavg)
mean(A[DGG.Cells(DGG.Covering(texas))])

# ## Any grid
#
# Nothing above is a HEALPix recipe: `covering_positions`, `query` and the
# selectors are interface methods, and any system in `DGG.systems()` runs this
# page unchanged. Two things do vary by system. Where cells are not equal-area,
# weight the mean by `cell_area`. And where a system's refinement is not
# congruent, a covering can name a few cells past the ones that touch —
# the multi-order coverage tutorial shows where that margin comes from.
