# # Multi-order coverage: one region, every resolution at once
#
# A discrete global grid has one cell size per level, and a region does not. Ask
# for California at a resolution fine enough to trace its coastline and you get
# tens of thousands of cells, almost all of them in the middle of the state
# where one cell a hundred times larger would have said the same thing.
#
# `MultiOrderCoverage` is the query that says it in mixed levels: descend from
# the roots, emit a cell whole as soon as the cell lies inside the target, and
# recurse only where the outline cuts through. The result is the picture below —
# the interior as a handful of giant hexagons, the coast and the Channel Islands
# at the leaf level — and it is also a data structure, because the cells come
# back in space-filling-curve order and expand to sorted, disjoint position
# ranges at any level.

import DiscreteGlobalGrids as DGG
import NaturalEarth
import GeometryOps as GO, GeoInterface as GI
using CairoMakie, GeoMakie
CairoMakie.activate!()

# CairoMakie rather than GLMakie or WGLMakie: this page draws thousands of
# polygons into a static figure, and Cairo renders `poly!` without a display.

sys = DGG.H3System()

fc = NaturalEarth.naturalearth("admin_1_states_provinces", 10)
california = fc.geometry[findfirst(==("California"), fc.name)]

# Natural Earth's California is a `MultiPolygon`: the mainland plus seven
# Channel Islands. Multipart targets need nothing special — the traversal sees
# one geometry.

GI.npolygon(california)

# ## The query

coverage = DGG.query(sys, DGG.MultiOrderCoverage(california); level = 7)

# Level 7 cells are about 2 km across. The set holds a few thousand cells
# spanning five levels, and stands for a leaf set twenty times larger:

(; n_cells = length(coverage),
   levels = extrema(DGG.level, coverage),
   n_leaves = sum(length, DGG.level_ranges(coverage, 7)),
   n_ranges = length(DGG.level_ranges(coverage, 7)))

# `level_ranges(coverage, 7)` is the compressed form: sorted, disjoint,
# adjacency-merged position ranges in `levelgrid(sys, 7)`. It is what a lookup
# layer slices arrays with, and it never materialises the leaf ids.

# ## The picture
#
# `cell_polygons` reads the geometry of a mixed-level set without the caller
# resolving a level grid per cell — the set knows which level each of its cells
# came from.

cells = collect(coverage)
polys = GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygons(coverage))

fig = Figure(size = (760, 820))
ax = GeoAxis(fig[1, 1]; dest = "+proj=longlat +datum=WGS84",
    limits = ((-125.5, -113.5), (32.0, 42.5)),
    title = "H3 multi-order coverage of California, to level 7")
plt = poly!(ax, polys; color = DGG.level.(cells), colormap = :viridis,
    strokecolor = (:black, 0.55), strokewidth = 0.25)
poly!(ax, california; color = :transparent, strokecolor = :black, strokewidth = 1.2)
Colorbar(fig[1, 2], plt; label = "cell level")
fig

# The interior is level 3 and 4 — cells a hundred kilometres across — and the
# coastline is level 7. Nothing chose those levels: the traversal emitted every
# cell at the coarsest level where the cell still fit inside the state.
#
# The white slivers between the big cells and their smaller neighbours are not a
# rendering artefact, and the next section is about them.

# ## The cells do not tile, on this system
#
# H3 has aperture 7: a cell's seven children are a rotated rosette that has the
# parent's area but not its footprint. Replacing a subtree by its root therefore
# replaces the subtree's shape by a different shape of the same size, and a
# mixed-level set has gaps in it. Aperture-4 systems whose four children tile
# their parent exactly — HEALPix, S2, ISEA4R — do not:

window = ((-121.6, -119.6), (35.6, 37.0))

pair = Figure(size = (900, 430))
for (k, (s, l, label)) in enumerate(((DGG.H3System(), 7, "H3 (aperture 7): slivers"),
                                     (DGG.S2System(), 9, "S2 (aperture 4): tiles")))
    cov = DGG.query(s, DGG.MultiOrderCoverage(california); level = l)
    ax = GeoAxis(pair[1, k]; dest = "+proj=longlat +datum=WGS84", limits = window,
        title = label, xticklabelsvisible = false, yticklabelsvisible = false)
    poly!(ax, GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygons(cov));
        color = DGG.level.(collect(cov)), colormap = :viridis,
        strokecolor = (:black, 0.5), strokewidth = 0.4)
end
pair

# What the set guarantees is a statement about the **leaf level**, not about the
# picture: every level-7 cell that meets California is a member of the set or a
# descendant of one, and no member is a descendant of another. That is the
# guarantee the compressed form rests on, and it holds on all seven systems.
# Expand the set before computing with it as a region, and read the figure as
# "which cells were chosen" rather than as the region itself.

# ## Which cells were proven to fit inside
#
# A coverage emits a cell for one of two reasons, and the difference matters.
# Most cells are emitted because the traversal asked whether they lie **inside**
# the target and the answer was yes. The rest are emitted at the deepest level,
# because the descent ran out of depth — they are in the set so that the set
# *covers* the target, and they are never asked the question at all.
#
# `is_contained` reports the first kind. It is a record of what was proven, not
# of what is true: a level-7 cell here may well sit entirely inside California
# and still read `false`, because nobody tested it. The containment predicate is
# the expensive one — some eighty times the allocation of an intersection test —
# and a deep coverage ends on thousands of those cells, so the traversal spends
# it only where the answer changes what it does next. `true` always means
# inside; `false` at the deepest level means untested.
#
# `argmin(level, coverage)` is therefore not "the coarsest cell inside
# California" — every emission can be unproven, and that idiom would silently
# return a boundary cell. `coarsest_contained` reads the flag instead, and
# answers `nothing` when nothing above the deepest level was proven to fit.

(; n_contained = count(i -> DGG.is_contained(coverage, i), eachindex(coverage)),
   coarsest = DGG.coarsest_contained(coverage))

# A target smaller than one cell is the clearest way to get `nothing`: every
# emission is at the deepest level, so nothing was proven and the accessor says
# so rather than handing back a boundary cell.

small = GI.Polygon([GI.LinearRing([(-122.42, 37.77), (-122.40, 37.77),
                                   (-122.40, 37.79), (-122.42, 37.79),
                                   (-122.42, 37.77)])])
tiny = DGG.query(sys, DGG.MultiOrderCoverage(small); level = 7)
(; n_cells = length(tiny), coarsest_contained = DGG.coarsest_contained(tiny))

# ## A budget instead of a depth
#
# Everything above fixes the *depth* and lets the cell count fall where it may.
# The other question is the one you ask when the cells are going on a slide, into
# a request payload, or into a coarse index: **give me ten cells that cover
# California, or a hundred.** That is the same query with the other keyword.

budgets = (10, 40, 100)
sets = [DGG.query(sys, DGG.MultiOrderCoverage(california); maxcells = n) for n in budgets]

for (n, set) in zip(budgets, sets)
    best = DGG.coarsest_contained(set)
    println("maxcells = ", lpad(n, 3), ": ", lpad(length(set), 4), " cells, levels ",
            extrema(DGG.level, set), ", coarsest proven inside: ",
            best === nothing ? "none" : "level $(DGG.level(best))")
end

# `maxcells` refines breadth first: it takes the coarsest cell the outline still
# crosses, replaces it by the children that meet California, and commits that
# swap only while the set it leaves behind still fits the budget. Contained cells
# are never touched — a cell already proven inside the state cannot be made a
# better statement about the state by splitting it.
#
# The budget is spent to the cell, and the depth is what it bought. At ten cells
# nothing has been refined far enough to sit inside California yet, so
# `coarsest_contained` says `nothing` rather than handing back a boundary cell.
# By forty, one has.

fig = Figure(size = (960, 400))
for (k, (n, set)) in enumerate(zip(budgets, sets))
    local panel = Axis(fig[1, k]; limits = ((-125.5, -113.5), (32.0, 42.5)),
        aspect = DataAspect(), title = "maxcells = $n",
        xticklabelsvisible = false, yticklabelsvisible = false)
    poly!(panel, GO.transform(GO.GeographicFromUnitSphere(), DGG.cell_polygons(set));
        color = DGG.level.(collect(set)), colormap = :viridis, colorrange = (1, 7),
        strokecolor = (:black, 0.6), strokewidth = 0.5)
    poly!(panel, california; color = :transparent, strokecolor = :black, strokewidth = 1.0)
end
fig

# Ten cells is a caricature of California, a hundred is recognisably the state,
# and the trade is visible in one direction only: more cells buy a tighter
# outline, never a different region.
#
# ### Does it actually cover?
#
# The covering statement is the same one the depth mode makes, read at the
# deepest level the budget happened to reach: every cell of that level which
# meets California is a member of the set or a descendant of one. Sample the
# state and ask it of each sample's own cell.

probes = [(lon, lat) for lon in -124.4:0.1:-114.2, lat in 32.6:0.1:42.0]

function uncovered(sys, set, points)
    leaf = maximum(DGG.level, set)
    grid = DGG.levelgrid(sys, leaf)
    members = Set(collect(set))
    count(points) do (lon, lat)
        c = DGG.cellat(grid, lon, lat)
        c === nothing && return true
        !(c in members || any(l -> DGG.ancestor(sys, c, l) in members,
                              (DGG.level(c)-1):-1:first(DGG.levels(sys))))
    end
end

interior = [p for p in vec(probes) if GO.contains(california, p)]
for (n, set) in zip(budgets, sets)
    println("maxcells = ", lpad(n, 3), ": deepest level ", maximum(DGG.level, set),
            ", ", uncovered(sys, set, interior), " of ", length(interior),
            " sampled points uncovered")
end

# Zero at every budget. That is a spot check and not the theorem: the statement
# is a theorem where four children tile their parent — HEALPix, S2, ISEA4R — and
# here, where they do not, refinement can drop a child that misses California
# while something beneath it does not. The package measures that rather than
# promising it, along the outline itself where it bites, in
# `test/systems/crosssystem/multiorder_budget.jl`.
#
# ## The second caveat: the compressed form is not universal
#
# The same non-congruence runs the other way too. A member's descendants are not
# confined to it, so expanding a set to the leaf level can name leaves the target
# does not touch — most visibly under a *hole* in the target, where a cell
# emitted whole outside the hole has descendants that land inside it. The set
# itself is exact; it is the expansion that over-covers, and only where the
# refinement over-covers.
#
# `level_ranges` exists because a
# cell's descendants occupy one contiguous interval of a deeper level's order.
# `has_sorted_subtrees(A5System())` is `false` — A5's descendants are scattered
# through their level — so `level_ranges` throws there rather than silently
# returning something that is not a range set. `descendants(sys, c, l)` is the
# always-available form, and generic code branches on the trait:

expand(sys, set, l) = DGG.has_sorted_subtrees(sys) ? DGG.cellindices(set, l) :
                      sort!(reduce(vcat, [DGG.descendants(sys, c, l) for c in set]))

# ## Every system, not just this one
#
# Multi-order coverage is generic. The same three lines run against all six
# registered systems and an `AuthalicSystem` wrap; levels are chosen per system
# so a leaf cell is roughly the same size on each, because the apertures differ.

sweep = [(DGG.IGeo7System(), 6), (DGG.H3System(), 5), (DGG.HEALPixSystem(), 9),
         (DGG.A5System(), 9), (DGG.S2System(), 9), (DGG.ISEA4RSystem(), 9),
         (DGG.AuthalicSystem(DGG.IGeo7System()), 6)]

for (s, l) in sweep
    cov = DGG.query(s, DGG.MultiOrderCoverage(california); level = l)
    name = s isa DGG.AuthalicSystem ? "Authalic($(nameof(typeof(parent(s)))))" :
           string(nameof(typeof(s)))
    ranges = DGG.has_sorted_subtrees(s) ? length(DGG.level_ranges(cov, l)) : missing
    best = DGG.coarsest_contained(cov)
    println(rpad(name, 24), "level $l: ", lpad(length(cov), 6), " cells, levels ",
            extrema(DGG.level, cov), ", ", lpad(length(expand(s, cov, l)), 7),
            " leaves, ", lpad(string(ranges), 7), " ranges, coarsest contained ",
            best === nothing ? "none" : string(DGG.level(best)))
end

# A5 is the one `missing` in the ranges column, and it is a stated exclusion
# rather than a failure: the trait says the compressed form does not exist
# there, and `expand` took the other branch to produce the same leaf set.

# ## Putting data on the cells
#
# A coverage names a region; it does not carry values. `CellLookup` reads the set
# as a one-level cell axis, `PartialGrid` reads that axis as a grid, and
# `ConservativeRegridding` fills it from anything with cell corners — here
# WorldClim's July mean temperature at 10 arc-minutes.
#
# The system is named explicitly for this section, and the reason is worth
# stating rather than hiding. Conservative regridding *onto* a DGGS is exact only
# where the destination cells' rings are convex: the regridder clips each source
# cell against the destination cell, and Sutherland-Hodgman is an intersection
# only when the clip window is convex. H3's rings are convex at even levels and
# not at odd ones, and a multi-order set is mixed-level by construction. IGEO7's
# are convex at every level, so the conservation check below comes out at machine
# precision. The regridding and hydrology pages take the other route — divide two
# fields regridded through the *same* matrix, which is an honest weighted mean
# whatever the rows sum to — and this page uses that trick too, for the oceans.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH",
    joinpath(tempdir(), "rasterdatasources")))

import ConservativeRegridding as CR
import DimensionalData as DD
import Extents
using Rasters, RasterDataSources
import ArchGDAL
using Statistics

igeo = DGG.IGeo7System()
region = DGG.query(igeo, DGG.MultiOrderCoverage(california); level = 6)
lk = DGG.CellLookup(region)
destination = DGG.PartialGrid(lk)

# 400-odd mixed-level entries stand for a thousand-odd level-6 cells, and the
# lookup stores the entries, not the cells.

(; entries = length(region), leaf_cells = length(lk), level = DGG.level(lk))

# The source is WorldClim cropped to a box comfortably larger than the coverage,
# described to the regridder as its matrix of cell-corner points on the unit
# sphere. `CR.Regridder(manifold, dst, src)` takes the destination first, so the
# intersection matrix is destination cells × source cells.

raster = Raster(WorldClim{Climate}, :tavg; month = 7, res = "10m")
box = raster[X(-128 .. -110), Y(30 .. 45)]

to_sphere = GO.UnitSpherical.UnitSphereFromGeographic()
(west, east), (south, north) = bounds(box, X), bounds(box, Y)
corners = [to_sphere((lon, lat))
           for lon in range(west, east; length = size(box, X) + 1),
               lat in range(south, north; length = size(box, Y) + 1)]

manifold = GO.Spherical(; radius = 1.0)
regridder = CR.Regridder(manifold, destination, corners)
size(regridder.intersections)

# Conservation, in the direction that is usually the broken one: regrid a field
# of ones and every destination cell must come back holding exactly one.

ones_out = zeros(DGG.ncells(destination))
CR.regrid!(ones_out, regridder, ones(size(regridder.intersections, 2)))
maximum(abs, ones_out .- 1)

# Machine precision, because these rings are convex. Now the real field. The
# ocean is `missing`, so the coverage trick from the regridding page applies:
# regrid the temperatures with the gaps zeroed *and* a 0/1 data indicator, then
# divide, so coastal cells are not dragged down by the empty half of themselves.

values = vec(reverse(parent(replace_missing(box, NaN)); dims = 2))
field = zeros(DGG.ncells(destination))
cover = zeros(DGG.ncells(destination))
CR.regrid!(field, regridder, replace(values, NaN => 0.0))
CR.regrid!(cover, regridder, Float64.(.!isnan.(values)))
tavg = field ./ cover
extrema(tavg)

# ## The set as a cube axis
#
# `Cells(lk)` is the DimensionalData dimension, and the array over it is an
# ordinary `DimArray` whose axis happens to be a compressed multi-order set —
# a thousand values whose labels are 115 stored windows. Its first few rows:

A = DD.DimArray(tavg, DGG.Cells(lk))
A[DGG.Cells(1:5)]

# Three selectors, three questions. `DD.At` takes a typed cell id, `DD.Contains`
# takes a lon/lat point and resolves it through `cellat`, and `Covering` takes a
# region and runs a coverage against the axis. `At` and `Contains` are
# DimensionalData's own and are reached through it: this package exports DE9IM's
# `Contains`, a predicate about two geometries, and the two names must not
# collide in a caller's namespace.

sanfrancisco = DGG.cellat(DGG.levelgrid(igeo, DGG.level(lk)), -122.42, 37.77)

(; at_id = A[DGG.Cells(DD.At(sanfrancisco))],
   at_point = A[DGG.Cells(DD.Contains((-118.24, 34.05)))])

# `Covering` answers with a view whose axis is a `CellLookup` again, so a
# subregion of the cube is still a cube over cells — and a mean over it is a
# zonal statistic in one line. IGEO7 cells are equal-area, so the unweighted mean
# over them is the areal mean.

for (name, ext) in (("Central Valley", Extents.Extent(X = (-121.5, -119.0), Y = (35.5, 38.0))),
                    ("north coast", Extents.Extent(X = (-123.5, -121.8), Y = (36.5, 39.0))))
    sub = A[DGG.Cells(DGG.Covering(ext))]
    println(rpad(name, 16), lpad(length(sub), 5), " cells, mean July temperature ",
            round(mean(sub); digits = 1), " °C")
end

# The valley runs some four degrees hotter than the coast at the same latitude,
# which is the answer the data has and the axis merely delivered.

fig = Figure(size = (620, 700))
ax = Axis(fig[1, 1]; limits = ((-125.0, -113.8), (32.2, 42.3)), aspect = DataAspect(),
    title = "July mean temperature on an IGEO7 level-$(DGG.level(lk)) coverage")
plt = poly!(ax, GO.transform(GO.GeographicFromUnitSphere(),
        [DGG.cell_polygon(DGG.levelgrid(igeo, DGG.level(lk)), c) for c in lk]);
    color = tavg, colormap = Reverse(:RdYlBu), strokewidth = 0)
poly!(ax, california; color = :transparent, strokecolor = :black, strokewidth = 1.0)
Colorbar(fig[1, 2], plt; label = "mean temperature (°C)")
fig

# The cells outside California are the coverage over-covering, exactly as the
# caveats above describe: the set is the coarsest cells that *meet* the state,
# and the leaf expansion of those is a superset of the state's own cells.
#
# A budget set backs the same axis. Nothing about `CellLookup` cares which mode
# produced the set, and a ten-cell coverage expanded to level 6 is a perfectly
# good — if generous — index into the same data.

budget_lk = DGG.CellLookup(DGG.query(igeo, DGG.MultiOrderCoverage(california);
        maxcells = 40, maxlevel = 5); level = DGG.level(lk))
(; budget_leaf_cells = length(budget_lk), exact_leaf_cells = length(lk))
