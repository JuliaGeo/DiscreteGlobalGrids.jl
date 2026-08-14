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

# ## Which cells fit inside, and which merely touch
#
# A coverage emits a cell for one of two reasons, and the difference matters.
# Most cells are emitted because they lie **inside** the target. The rest are
# emitted because the traversal ran out of depth with the outline still cutting
# through them — they are in the set so that the set *covers* the target.
#
# `argmin(level, coverage)` is therefore not "the coarsest cell inside
# California": when the target is smaller than one leaf cell, every emission is
# a crossing and that idiom silently returns one. `coarsest_contained` asks the
# question properly and answers `nothing` when there is no such cell.

(; n_contained = count(i -> DGG.is_contained(coverage, i), eachindex(coverage)),
   coarsest = DGG.coarsest_contained(coverage))

# A target smaller than a cell has no contained cell at all, and says so:

small = GI.Polygon([GI.LinearRing([(-122.42, 37.77), (-122.40, 37.77),
                                   (-122.40, 37.79), (-122.42, 37.79),
                                   (-122.42, 37.77)])])
tiny = DGG.query(sys, DGG.MultiOrderCoverage(small); level = 7)
(; n_cells = length(tiny), coarsest_contained = DGG.coarsest_contained(tiny))

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
