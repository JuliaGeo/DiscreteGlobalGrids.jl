# # Zonal statistics
#
# Which cells of a discrete global grid lie inside a polygon, and what is the
# mean of a field over them? `query` answers the first question directly, on any
# grid in this package, with a DE9IM predicate that says exactly which sense of
# "inside" is meant. This page answers it for Texas on a whole-globe HEALPix
# grid, then shows the multi-order form, which returns the same region as a few
# hundred mixed-level cells instead of a few thousand leaf ones.

import DiscreteGlobalGrids as DGG
import NaturalEarth
import GeometryOps as GO
using Statistics
using CairoMakie, GeoMakie
CairoMakie.activate!()

# ## Texas and a field on the globe
#
# Texas comes from Natural Earth's admin-1 dataset. The grid is every HEALPix
# cell at level 7 (`nside = 128`, 196608 cells), and the data is a smooth
# synthetic field sampled at the cell centroids.

fc = NaturalEarth.naturalearth("admin_1_states_provinces", 50)
texas = fc.geometry[findfirst(==("Texas"), fc.name)]

grid = DGG.levelgrid(DGG.HEALPixSystem(), 7)
lonlat = GO.UnitSpherical.GeographicFromUnitSphere()
field(lon, lat) = 20 - 0.5 * (lat - 25) + 2 * sind(3 * lon)

# ## The query
#
# `Within(texas)` keeps only the cells wholly inside the outline; `Intersects`
# keeps every cell that touches it. Both return **sorted typed cell ids**, and
# `cellposition` turns an id into the position that indexes a data vector — the
# two directions of the one bijection the interface is built on.
#
# The descent prunes whole subtrees against the conservative `node_extent` caps,
# so of 196608 cells only the handful straddling the outline is ever handed to
# an exact spherical predicate.

interior = DGG.query(grid, DGG.Within(texas))
touching = DGG.query(grid, DGG.Intersects(texas))
(; n_interior = length(interior), n_touching = length(touching))

# HEALPix cells are equal-area, so the unweighted mean over a cell set *is* the
# areal mean — no latitude weights. The two predicates bracket the answer, and
# at this resolution they agree to a hundredth of a degree.

zonal(cells) = mean(field(lonlat(DGG.cell_centroid(grid, c))...) for c in cells)
(; interior = zonal(interior), touching = zonal(touching))

# The rim is the difference of the two sets, and every predicate the engine
# implements is available the same way — `Covers`, `Touches`, `Overlaps`,
# `Disjoint`, and the rest.

rim = setdiff(touching, interior)
length(rim)

# ## The multi-order form
#
# `MultiOrderCoverage` asks the same question but keeps the answer compressed:
# the interior is emitted at the **coarsest level that is still wholly inside**,
# and only the boundary is refined down to the leaf. What comes back is a
# mixed-level `MultiOrderCellSet` in space-filling-curve order.

set = DGG.query(DGG.HEALPixSystem(), DGG.MultiOrderCoverage(texas); level = 7)
set

# `level_ranges` expands it to sorted, disjoint **position** ranges at any level
# no shallower than its coarsest cell. That is the handshake a lazy cell axis
# consumes: the ranges are what gets stored, and the leaf ids are computed from
# them on demand.

ranges = DGG.level_ranges(set, 7)
(; n_cells = length(set), n_ranges = length(ranges), n_leaves = sum(length, ranges))

# Coverage means *covering*, not *covered by*: the emitted set covers Texas, so
# it is a superset of `touching` rather than of `interior`.

leaves = [DGG.cellindex(grid, i) for r in ranges for i in r]
issubset(touching, leaves)

# ## Plotting
#
# Only the plot wants lon/lat; everything above stayed on the sphere.

topolys(cells) = GO.transform(GO.GeographicFromUnitSphere(),
    [DGG.cell_polygon(DGG.levelgrid(DGG.HEALPixSystem(), DGG.level(c)), c) for c in cells])

fig = Figure(size = (700, 620))
ax = GeoAxis(fig[1, 1]; dest = "+proj=longlat +datum=WGS84",
    limits = ((-107.5, -92.5), (25.0, 37.5)),
    title = "HEALPix multi-order coverage of Texas")
poly!(ax, topolys(collect(set)); color = DGG.level.(collect(set)),
    colormap = :viridis, strokecolor = :white, strokewidth = 0.5)
poly!(ax, texas; color = :transparent, strokecolor = :black, strokewidth = 2)
fig

# Cells are coloured by their level: the interior arrives as a few coarse cells,
# the outline as many fine ones.

# ## The DimensionalData layer, as it should read
#
# A cube over this region wants ONE dimension naming its cells at the leaf
# level, backed by the coverage rather than by the expanded id vector — memory
# `O(length(set))`, not `O(sum(length, ranges))`.
#
# !!! warning "Aspirational — T16"
#     Nothing implements the block below yet. It is the API this page wants,
#     written out as the acceptance test for the DimensionalData layer, and it
#     is a plain code fence rather than an `@example`, so the build never runs
#     it.
#
# ```julia
# lk = DGG.CellLookup(set)                     # MOC-backed, leaf level 7
# length(lk) == sum(length, ranges)            # the logical content
# lk[k]                                        # position -> typed leaf id, on demand
# DGG.cellposition(lk, c)                      # typed leaf id -> position, or nothing
#
# A = DD.DimArray(values, DGG.Cells(lk); name = :tavg)
# A[DGG.Cells(DD.At(c))]                       # a typed id
# A[DGG.Cells(DD.Contains(-97.7, 30.3))]       # a lon/lat point, via `cellat`
# A[DGG.Cells(DGG.Covering(travis_county))]    # a polygon: coverage ∩ backing
#
# mean(A[DGG.Cells(DGG.Covering(texas))])      # this whole page, in one line
# ```

# ## The same query on every system
#
# `query` is an interface method, so the zonal selection above is not a HEALPix
# recipe. Only the cell counts differ.

levels = Dict(DGG.IGeo7System => 4, DGG.H3System => 3)   # aperture 7; the rest are 4
for sys in (DGG.systems()..., DGG.AuthalicSystem(DGG.HEALPixSystem()))
    base = sys isa DGG.AuthalicSystem ? parent(sys) : sys
    l = get(levels, typeof(base), 6)
    g = DGG.levelgrid(sys, l)
    inner, outer = DGG.query(g, DGG.Within(texas)), DGG.query(g, DGG.Intersects(texas))
    name = sys isa DGG.AuthalicSystem ? "Authalic($(nameof(typeof(base))))" :
           string(nameof(typeof(sys)))
    println(rpad(name, 24), "level $l: ", lpad(length(inner), 6), " within, ",
            lpad(length(outer), 6), " touching")
end
