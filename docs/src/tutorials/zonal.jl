# # Zonal statistics
#
# Which cells of a discrete global grid lie inside a polygon, and what is the
# mean of a field over them? `query` answers the first question directly, on any
# grid in this package, with a DE9IM predicate that says exactly which sense of
# "inside" is meant. This page answers it for Texas on a whole-globe HEALPix
# grid, then shows the multi-order form, which returns the same region as a few
# hundred mixed-level cells instead of a few thousand leaf ones, and finally
# hands that compressed answer to a DimensionalData cell axis, which keeps it
# compressed.

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GeoInterface as GI
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

# ## The DimensionalData layer
#
# A cube over this region wants ONE dimension naming its cells at the leaf
# level, backed by the coverage rather than by the expanded id vector — memory
# `O(length(set))`, not `O(sum(length, ranges))`. `CellLookup` is that
# dimension's lookup: semantically the leaf id vector, structurally the ranges.

lk = DGG.CellLookup(set)
(; n_cells = length(lk), n_windows = length(ranges), bytes = Base.summarysize(lk))

# `lk[k]` is the `k`th leaf id, resolved on demand; `cellposition(lk, c)` is the
# inverse, and `nothing` for a cell the region does not hold. Those two are the
# whole of what a cube axis needs.

k = length(lk) ÷ 2
lk[k], DGG.cellposition(lk, lk[k])

# `Cells` is the dimension the lookup goes in.

A = DD.DimArray([field(lonlat(DGG.cell_centroid(grid, c))...) for c in lk],
    DGG.Cells(lk); name = :tavg)

# A typed id selects one cell, a lon/lat point selects the cell it falls in
# (through `cellat`), and a polygon selects through the same `MultiOrderCoverage`
# this page opened with, intersected with the axis.
#
# `At` and `Contains` are spelled `DD.At` and `DD.Contains` throughout: they are
# DimensionalData's selectors, and this package exports its own `Contains` — the
# DE9IM predicate about geometries — which would collide with the selector under
# a plain `using`.

A[DGG.Cells(DD.At(lk[k]))], A[DGG.Cells(DD.Contains(-97.7, 30.3))]

# Austin sits in Travis County; a box around it selects a few hundred of the
# region's cells, and the view keeps the compact backing.

travis = GI.Polygon([GI.LinearRing([(-98.2, 30.0), (-97.35, 30.0), (-97.35, 30.65),
    (-98.2, 30.65), (-98.2, 30.0)])])
austin = A[DGG.Cells(DGG.Covering(travis))]
(; n_cells = length(austin), lookup = DD.lookup(austin, DGG.Cells), tavg = mean(austin))

# And the zonal mean this page opened with is one line over the cube. It is the
# mean over the *coverage's* leaves, so it brackets with `touching` rather than
# with `interior` — the same distinction the two queries above drew.

mean(A[DGG.Cells(DGG.Covering(texas))])

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
