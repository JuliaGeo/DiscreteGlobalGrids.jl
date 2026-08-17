# # Multi-order coverage: one region, every resolution at once
#
# A discrete global grid has one cell size per level, and a region does not.
# Ask for California at a resolution fine enough to trace its coastline and you
# get tens of thousands of cells, almost all of them interior, where one cell a
# hundred times larger says the same thing.
#
# `MultiOrderCoverage` answers in mixed levels: descend from the roots, emit a
# cell whole as soon as it lies inside the target, recurse only where the
# outline cuts through. There are two ways to stop — fix the *depth* and let
# the cell count fall out (`level =`), or fix the *count* and let the depth
# fall out (`maxcells =`). This page shows both, then puts data on the result.

import DiscreteGlobalGrids as DGG
import NaturalEarth
import GeometryOps as GO, GeoInterface as GI
using GLMakie, GeoMakie
GLMakie.activate!(inline = true)

fc = NaturalEarth.naturalearth("admin_1_states_provinces", 10)
california = fc.geometry[findfirst(==("California"), fc.name)]

# California is a `MultiPolygon` — the mainland plus seven Channel Islands —
# and multipart targets need nothing special.

# ## Accuracy first: fix the level

sys = DGG.IGeo7System()
coverage = DGG.query(sys, DGG.MultiOrderCoverage(california); level = 7)

# Level 7 cells are about 6.8 km across (by √area, the convention on these
# pages). How many cells came back, and how many level-7 leaves they stand for:

(; n_cells = length(coverage),
   levels = extrema(DGG.level, coverage),
   n_leaves = sum(length, DGG.level_ranges(coverage, 7)))

# `level_ranges(coverage, 7)` is the compressed form — sorted, disjoint
# position ranges at level 7. It is what a lookup layer slices arrays with,
# and it never materialises the leaf ids. The query is generic — the same
# call runs on every registered system — but the compressed form needs
# sorted subtrees: on A5, `level_ranges` throws, and `descendants` is the
# always-available expansion.

cali_centroid = GO.centroid(california)

fig = Figure(size = (760, 820))
ax = GeoAxis(fig[1, 1]; dest = "+proj=ortho +datum=WGS84 +lon_0=$(GI.x(cali_centroid)) +lat_0=$(GI.y(cali_centroid))",
    limits = ((-125.5, -113.5), (32.0, 42.5)),
    title = "IGeo7 multi-order coverage of California, to level 7")
plt = poly!(ax, coverage; color = DGG.level.(coverage), colormap = :isoluminant_cgo_70_c39_n256,
    alpha = 0.7, strokecolor = (:black, 0.55), strokewidth = 0.25)
poly!(ax, california; color = :transparent, strokecolor = :black, strokewidth = 1.2)
Colorbar(fig[1, 2], plt; label = "cell level")
fig

# The interior is level 4 and 5 — cells up to some 130 km across — and the
# coastline is level 7. Nothing chose those levels: every cell was emitted at
# the coarsest level where it still fit inside the state.
#
# The white slivers between big cells and their smaller neighbours are real.
# IGEO7 has aperture 7: a cell's seven children are a rotated rosette with the
# parent's area but not its footprint, so a mixed-level set does not tile.
# What the set guarantees is a statement about the **leaf level** — every
# level-7 cell that meets California is a member or a descendant of one, and
# no member is a descendant of another — which is exactly the guarantee the
# compressed form rests on. Aperture-4 systems whose four children tile their
# parent (HEALPix, S2, ISEA4R) have no slivers. Expand the set before
# computing with it as a region; read the figure as "which cells were chosen".

# ## Budget first: fix the cell count
#
# The other question is the one you ask when the cells are going on a slide,
# into a request payload, or into a coarse index: **give me ten cells that
# cover California, or a hundred.** Same query, other keyword.

budgets = (10, 40, 100)
sets = [DGG.query(sys, DGG.MultiOrderCoverage(california); maxcells = n) for n in budgets]

fig = Figure(size = (960, 400))
for (k, (n, set)) in enumerate(zip(budgets, sets))
    local panel = Axis(fig[1, k]; limits = ((-125.5, -113.5), (32.0, 42.5)),
        aspect = DataAspect(), title = "maxcells = $n",
        xticklabelsvisible = false, yticklabelsvisible = false)
    poly!(panel, set;
        color = DGG.level.(set), colormap = :isoluminant_cgo_70_c39_n256, colorrange = (1, 7),
        strokecolor = (:black, 0.6), strokewidth = 0.5)
    poly!(panel, california; color = :transparent, strokecolor = :black, strokewidth = 1.0)
end
fig

# `maxcells` refines breadth first: take the coarsest cell the outline still
# crosses, replace it by the children that meet California, and commit the
# swap only while the set still fits the budget. Cells already proven inside
# the state are never split — splitting cannot make them a better statement
# about it. Ten cells is a caricature of California, a hundred is recognisably
# the state, and more cells only ever buy a tighter outline, never a different
# region. The covering statement holds at the deepest level the budget
# reached; it is a theorem on the aperture-4 systems and measured on the rest,
# in `test/systems/crosssystem/multiorder_budget.jl`.
#
# A budget set also knows which of its cells were *proven* to sit inside the
# state:

for (n, set) in zip(budgets, sets)
    best = DGG.coarsest_contained(set)
    println("maxcells = ", lpad(n, 3), ": ", lpad(length(set), 3), " cells, levels ",
            extrema(DGG.level, set), ", coarsest inside: ",
            best === nothing ? "none" : "level $(DGG.level(best))")
end

# At ten cells nothing has been refined far enough to fit inside California
# yet, so `coarsest_contained` answers `nothing` rather than handing back a
# boundary cell. By forty, one has. (On a `level =` set the flag is
# provenance, not geometry: `true` always means proven inside, but a cell
# emitted at the deepest level was never tested and reads `false` either way.)

# ## Putting data on the cells
#
# A coverage names a region; it does not carry values. `CellLookup` reads the
# set as a one-level cell axis, and `regrid` fills it from anything with cell
# corners — here a deterministic temperature-like field on a regular lon/lat
# raster.

import DimensionalData as DD
import Extents
using Rasters
using Statistics

region = DGG.query(sys, DGG.MultiOrderCoverage(california); level = 6)
lk = DGG.CellLookup(region)

# 400-odd mixed-level entries stand for a thousand-odd level-6 cells, and the
# lookup stores the entries, not the cells.

(; entries = length(region), leaf_cells = length(lk), level = DGG.level(lk))

# The source covers a box larger than the coverage.

lon = -127.95:0.1:-110.05
lat = 44.95:-0.1:30.05
july_temperature(x, y) = 30 - 0.45(y - 30) +
    4exp(-((x + 120) / 1.8)^2) + 1.5sind(3x)
box = Raster(
    [july_temperature(x, y) for x in lon, y in lat],
    (X(lon; sampling = Rasters.Intervals(Rasters.Center())),
     Y(lat; sampling = Rasters.Intervals(Rasters.Center()))),
)

# ## The set as a cube axis
#
# `to` takes the lookup as it stands. The result is an ordinary `DimArray` whose
# axis is a `Cells` dimension over the same cells, so nothing has to be looked
# up again to plot it or to slice it.

A = DGG.regrid(box; to = lk)
tavg = parent(A)
extrema(tavg)

# Slicing it is DimensionalData's business, and the selectors read cells.

A[DGG.Cells(DD.Contains((-118.24, 34.05)))]  # July mean at Los Angeles

# (`DD.Contains` is DimensionalData's point selector, reached through `DD`
# because this package exports DE9IM's `Contains`, a different thing.)
#
# `Covering` takes a region and runs a coverage against the axis, answering
# with a smaller `DimArray` whose axis is a `CellLookup` again — a subregion
# of the cube is still a cube over cells, and a mean over it is a zonal
# statistic in one line. IGEO7 cells are equal-area, so the unweighted mean is
# the areal mean.

for (name, ext) in (("Central Valley", Extents.Extent(X = (-121.5, -119.0), Y = (35.5, 38.0))),
                    ("north coast", Extents.Extent(X = (-123.5, -121.8), Y = (36.5, 39.0))))
    sub = A[DGG.Cells(DGG.Covering(ext))]
    println(rpad(name, 16), lpad(length(sub), 5), " cells, mean July temperature ",
            round(mean(sub); digits = 1), " °C")
end

# The analytic inland ridge runs hotter than the coast at the same latitude —
# the answer the data has, delivered by the axis.

fig = Figure(size = (620, 700))
ax = Axis(fig[1, 1]; limits = ((-125.0, -113.8), (32.2, 42.3)), aspect = DataAspect(),
    title = "July mean temperature on an IGEO7 level-$(DGG.level(lk)) coverage")
plt = poly!(ax, lk; color = tavg, colormap = Reverse(:RdYlBu), strokewidth = 0)
poly!(ax, california; color = :transparent, strokecolor = :black, strokewidth = 1.0)
Colorbar(fig[1, 2], plt; label = "mean temperature (°C)")
fig

# The cells outside California are the coverage over-covering, as the sliver
# section predicts from the other side: the set holds the coarsest cells that
# *meet* the state, and their leaf expansion is a superset of the state's own
# cells.
#
# A budget set backs the same axis. `CellLookup` does not care which mode
# produced the set — a forty-cell coverage expanded to level 6 is a perfectly
# good, if generous, index into the same data. (`maxlevel = 5` keeps this set
# coarser than the level it expands to; the forty-cell set above was
# unconstrained and reached level 6 itself.)

budget = DGG.query(sys, DGG.MultiOrderCoverage(california); maxcells = 40, maxlevel = 5)
budget_lk = DGG.CellLookup(budget; level = DGG.level(lk))
(; budget_leaf_cells = length(budget_lk), exact_leaf_cells = length(lk))
