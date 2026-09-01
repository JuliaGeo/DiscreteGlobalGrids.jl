# # Multi-order coverage

import DiscreteGlobalGrids as DGG                                                # hide
import NaturalEarth                                                              # hide
import GeometryOps as GO, GeoInterface as GI                                     # hide
using DiscreteGlobalGridsVisualization: dggpoly, dggpoly!                        # hide
using GLMakie, GeoMakie                                                          # hide
GLMakie.activate!(inline = true)                                                 # hide
fc = NaturalEarth.naturalearth("admin_1_states_provinces", 10)                   # hide
california = fc.geometry[findfirst(==("California"), fc.name)]                    # hide
sys = DGG.HEALPixSystem()                                                        # hide
moc = DGG.query(sys, DGG.MultiOrderCoverage(california); level = 9)              # hide
albers = "+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +datum=WGS84"     # hide
levels = DGG.level.(moc)                                                         # hide
lo, hi = extrema(levels)                                                         # hide
tints = cgrad(["#dcf5d7", "#9fd894", "#5cb84c", "#389826", "#2c7a1e"],           # hide
    hi - lo + 1; categorical = true)                                             # hide
fig = Figure(size = (720, 780))                                                  # hide
ax = GeoAxis(fig[1, 1]; dest = albers,                                           # hide
    limits = ((-125.5, -113.5), (32.0, 42.5)),                                   # hide
    title = "A multi-order coverage of California, HEALPix to level 9")          # hide
plt = dggpoly!(ax, moc; color = levels, colormap = tints,                        # hide
    colorrange = (lo - 0.5, hi + 0.5),                                           # hide
    strokecolor = "#2c7a1e", strokewidth = 0.4)                                  # hide
poly!(ax, california; color = :transparent,                                      # hide
    strokecolor = ("#212529", 0.9), strokewidth = 1.6)                           # hide
Colorbar(fig[1, 2], plt; ticks = lo:hi, label = "cell level")                    # hide
fig                                                                              # hide

# A multi-order coverage is a set of cells at mixed levels that covers a region:
# coarse cells wherever the region's interior fills them, finer cells along its
# boundary. Above, nineteen level-6 cells carry the interior of California and
# the coastline is traced at level 9.
#
# It is index compression. In a hierarchical system a parent's index stands for
# all of its children at the target level, so a subtree that lies wholly inside
# the region is stored as its root — and that root's parent, recursively upward,
# wherever the whole subtree fits. The 729 entries drawn above name 2721
# level-9 cells, and `level_ranges` expands them to sorted index ranges when the
# data is read.

# ## Query California as a coverage in HEALPix
#
# [`MultiOrderCoverage`](@ref) is a query target. Hand it a region and the
# finest level you will accept, and the query answers with the mixed-level set.

import DiscreteGlobalGrids as DGG
import NaturalEarth
import GeometryOps as GO, GeoInterface as GI
using DiscreteGlobalGridsVisualization: dggpoly, dggpoly!
using GLMakie, GeoMakie
GLMakie.activate!(inline = true)

fc = NaturalEarth.naturalearth("admin_1_states_provinces", 10)
california = fc.geometry[findfirst(==("California"), fc.name)]

#

sys = DGG.HEALPixSystem()
moc = DGG.query(sys, DGG.MultiOrderCoverage(california); level = 9)

# `level = 9` fixes the finest cell the set may hold. Its width in metres, by
# √area — the convention on these pages:

DGG.cellsize(sys, 9)

# Every entry carries the level it was emitted at. Counting them by level shows
# where the coverage spent its cells:

levels = DGG.level.(moc)
[l => count(==(l), levels) for l in minimum(levels):maximum(levels)]

# `level_ranges` reads the set as sorted, disjoint ranges of level-9 indices —
# the form a lookup layer slices arrays with, without ever materialising the ids
# one at a time:

ranges = DGG.level_ranges(moc, 9)

# Their total length is the number of level-9 cells the set stands for:

sum(length, ranges)

# ## Put an elevation raster onto the coverage
#
# A coverage names a region and `regrid` fills it with values. The source here
# is [CRU CL 2.0](https://zenodo.org/records/20754689), a 10-arcminute
# climatology of the world's land, whose `:elv` layer holds the mean elevation
# of each grid cell in metres.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))
using Rasters, RasterDataSources
import NCDatasets
import DimensionalData as DD
import Extents
using Statistics

world = Raster(RasterDataSources.getraster(CRUCL2); name = :elv, lazy = true)
elevation = read(world[X(-125.5 .. -113.5), Y(32.0 .. 42.5)])
elevation = DD.rebuild(elevation; metadata = DD.NoMetadata())

# The crop is a box, so it carries a margin of Nevada and Arizona along with
# California. `surface!` draws a raster on the `GeoAxis` these pages use for
# cells; pushing it back in `z` keeps the state outline on top:

albers = "+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +datum=WGS84"

fig = Figure(size = (700, 760))
ax = GeoAxis(fig[1, 1]; dest = albers, limits = ((-125.5, -113.5), (32.0, 42.5)),
    title = "CRU CL 2.0 elevation, 10 arcmin")
plt = surface!(ax, replace_missing(elevation, NaN);
    colormap = :batlow, shading = NoShading)
translate!(plt, 0, 0, -100)
poly!(ax, california; color = :transparent,
    strokecolor = ("#212529", 0.9), strokewidth = 1.5)
Colorbar(fig[1, 2], plt; label = "elevation (m)")
fig

# `regrid` takes the coverage as its destination:

A = DGG.regrid(elevation; to = moc)

# The result is a `Raster` over a single `Cells` dimension, one value per
# level-9 leaf cell. The cells offshore hold `missing`, because the source
# covers land:

count(ismissing, A)

# DimensionalData's selectors read that axis. `DD.Contains` picks the cell
# holding a point — spelled through `DD` because this package exports DE9IM's
# `Contains`, a geometry predicate:

A[DGG.Cells(DD.Contains((-118.24, 34.05)))]  # Los Angeles

# `Covering` takes a region and answers with a smaller `Raster` whose axis is a
# cell axis again, so a statistic over a subregion is one more line:

sierra = A[DGG.Cells(DGG.Covering(Extents.Extent(X = (-119.5, -118.0), Y = (36.5, 38.0))))]

#

mean(skipmissing(sierra))

# `replace_missing` puts `NaN` in the missing cells, and `dggpoly!` leaves those
# cells undrawn:

A = replace_missing(A, NaN)

fig = Figure(size = (700, 760))
ax = GeoAxis(fig[1, 1]; dest = albers, limits = ((-125.5, -113.5), (32.0, 42.5)),
    title = "Elevation on the coverage, HEALPix to level 9")
plt = dggpoly!(ax, A; color = A, colormap = :batlow)
poly!(ax, california; color = :transparent,
    strokecolor = ("#212529", 0.9), strokewidth = 1.5)
Colorbar(fig[1, 2], plt; label = "elevation (m)")
fig

# The Sierra Nevada runs down the eastern edge, Death Valley sits below sea
# level to its south-east, and the Central Valley is the low band between the
# Sierra and the Coast Ranges. Cells reach past the state line because the
# coverage holds the coarsest cells that *meet* California, so its leaf
# expansion is a superset of the state's own cells.

# ## Cap the cell count with `maxcells`
#
# `maxcells` asks the other question: cover California in ten cells, or forty,
# or a hundred. The query refines coarsest first.
#
# 1. Take the coarsest cell the outline crosses.
# 2. Replace it by the children that meet California.
# 3. Keep it whole when the replacement would exceed the budget.
#
# Cells already proven inside the state stay as they are.

budgets = (10, 40, 100)
sets = [DGG.query(sys, DGG.MultiOrderCoverage(california); maxcells = n) for n in budgets]

# Each set stops where its budget runs out, so each spans its own levels. The
# budget is the only limit here, so the hundred-cell set refines to level 11
# along the coast. The three panels share one colour scale, taken from all of
# them:

levelrange = extrema(DGG.level(c) for set in sets for c in set)

#

budgettints = cgrad(["#dcf5d7", "#9fd894", "#5cb84c", "#389826", "#2c7a1e"],
    levelrange[2] - levelrange[1] + 1; categorical = true)

fig = Figure(size = (1000, 430))
for (k, (n, set)) in enumerate(zip(budgets, sets))
    local panel = GeoAxis(fig[1, k]; dest = albers,
        limits = ((-127.0, -112.0), (30.0, 44.0)), title = "maxcells = $n",
        xticklabelsvisible = false, yticklabelsvisible = false,
        xgridvisible = false, ygridvisible = false)
    dggpoly!(panel, set; color = DGG.level.(set), colormap = budgettints,
        colorrange = (levelrange[1] - 0.5, levelrange[2] + 0.5),
        strokecolor = "#2c7a1e", strokewidth = 0.5)
    poly!(panel, california; color = :transparent,
        strokecolor = ("#212529", 0.9), strokewidth = 1.0)
end
Colorbar(fig[1, 4]; colormap = budgettints,
    colorrange = (levelrange[1] - 0.5, levelrange[2] + 0.5),
    ticks = levelrange[1]:levelrange[2], label = "cell level")
fig

# ## The same coverage on IGeo7
#
# IGeo7 level 6 cells come nearest the 10-arcminute source in width:

igeo7 = DGG.IGeo7System()
DGG.cellsize(igeo7, 6)

#

moc7 = DGG.query(igeo7, DGG.MultiOrderCoverage(california); level = 6)

#

levels7 = DGG.level.(moc7)
lo7, hi7 = extrema(levels7)
tints7 = cgrad(["#dcf5d7", "#9fd894", "#5cb84c", "#389826", "#2c7a1e"],
    hi7 - lo7 + 1; categorical = true)

fig = Figure(size = (720, 780))
ax = GeoAxis(fig[1, 1]; dest = albers, limits = ((-125.5, -113.5), (32.0, 42.5)),
    title = "A multi-order coverage of California, IGeo7 to level 6")
plt = dggpoly!(ax, moc7; color = levels7, colormap = tints7,
    colorrange = (lo7 - 0.5, hi7 + 0.5),
    strokecolor = "#2c7a1e", strokewidth = 0.4)
poly!(ax, california; color = :transparent,
    strokecolor = ("#212529", 0.9), strokewidth = 1.6)
Colorbar(fig[1, 2], plt; ticks = lo7:hi7, label = "cell level")
fig

# HEALPix children cover exactly their parent, so a HEALPix set at mixed levels
# tiles the region. IGeo7's seven children carry their parent's area with a
# different footprint, so where the level changes here the drawn cells overlap
# their coarser neighbour and leave slivers around it. The guarantee holds at
# the leaf level: every level-6 cell that meets California belongs to the set or
# descends from a member, and no member descends from another.
# [`MultiOrderCoverage`](@ref) tabulates how far each system departs from a
# tiling.
#
# The raster regrids onto it by the same call:

A7 = replace_missing(DGG.regrid(elevation; to = moc7), NaN)

#

fig = Figure(size = (700, 760))
ax = GeoAxis(fig[1, 1]; dest = albers, limits = ((-125.5, -113.5), (32.0, 42.5)),
    title = "Elevation on the coverage, IGeo7 to level 6")
plt = dggpoly!(ax, A7; color = A7, colormap = :batlow)
poly!(ax, california; color = :transparent,
    strokecolor = ("#212529", 0.9), strokewidth = 1.5)
Colorbar(fig[1, 2], plt; label = "elevation (m)")
fig

# ## What the query is doing
#
# The coverage walk is a depth-first search over the system's hierarchy.
# `DGG.treeify` exposes a grid as a `GeometryOps.SpatialTreeInterface` tree, and
# the search reads four things from each node: its bounding cap, its cell,
# whether it is a leaf, and its children. At each node:
#
# 1. drop it when its cap misses the target;
# 2. drop it when its cell polygon misses the target;
# 3. emit it when it is a leaf, or when its polygon lies inside the target;
# 4. otherwise descend into its children.

import GeometryOps.SpatialTreeInterface as STI
import GeometryOps.UnitSpherical as US

# The predicates run on the sphere, so the region moves to unit-sphere
# coordinates once. `Float64` because `RelateNG` works in double precision and
# Natural Earth ships `Float32`:

tosphere = US.UnitSphereFromGeographic()
onsphere(geom) = GO.transform(p -> tosphere((Float64(GI.x(p)), Float64(GI.y(p)))), geom)

# A node's extent is a `SphericalCap`, so the target needs one too: centred on
# its centroid, with the radius reaching its farthest vertex.

function boundingcap(geom)
    centre = tosphere(GO.centroid(geom))
    radius = maximum(US.spherical_distance(centre, tosphere(p))
                     for p in GO.flatten(GI.PointTrait, geom))
    return US.SphericalCap(centre, radius)
end

# `cell_boundary` returns a cell's corners on the unit sphere for a cell id at
# any level, so the walk builds its polygons straight from the system, with no
# per-level grid in sight:

cellpolygon(sys, c) = (ring = DGG.cell_boundary(sys, c); GI.Polygon([GI.LinearRing([ring; ring[1:1]])]))

#

function coverage(sys, region, maxlevel)
    target = onsphere(region)
    cap = boundingcap(region)
    prepared = GO.prepare(GO.RelateNG(; manifold = GO.Spherical()), target)
    meets(poly) = GO.relate_predicate(prepared, GO.pred_intersects(), poly)
    inside(poly) = GO.relate_predicate(prepared, GO.pred_contains(), poly)
    out = DGG.cellindextype(sys)[]
    function visit(node)
        Extents.intersects(cap, STI.node_extent(node)) || return nothing
        c = DGG.node_cell(node)   # `nothing` at the synthetic root above the base cells
        if c !== nothing
            poly = cellpolygon(sys, c)
            meets(poly) || return nothing
            if STI.isleaf(node) || inside(poly)
                push!(out, c)
                return nothing
            end
        end
        for child in STI.getchild(node)
            visit(child)
        end
        return nothing
    end
    visit(DGG.treeify(DGG.levelgrid(sys, maxlevel)))
    return out
end

# Compare it with `query` at level 7:

walked7 = coverage(sys, california, 7)
queried7 = DGG.query(sys, DGG.MultiOrderCoverage(california); level = 7)
length(walked7), length(queried7), Set(walked7) == Set(queried7)

# And at level 9, the coverage this page has been using:

walked9 = coverage(sys, california, 9)
length(walked9), length(moc), Set(walked9) == Set(moc)

# The same sets, both times. `query` adds two accelerations to the same walk: a
# proof that a cap lying on the far side of the target's boundary is wholly
# outside, which discards whole subtrees without building a polygon, and a
# cheaper accept for cells whose centroid falls inside.

# ## The cells behind the axis
#
# `regrid` and the selectors read the coverage through the calls below. Reach
# for them when you want the cell ids themselves.
#
# `CellLookup` reads a set as a one-level cell axis. Its window count is how
# many contiguous runs of leaf ids the set compresses to:

lk = DGG.CellLookup(moc)

# `expand` flattens the set to one level as a `CellVector`, the cells you would
# compute on:

flat = DGG.expand(moc, 9)

# `compact` merges every complete sibling group into its parent, as far up as it
# can go. It reads membership alone, so along the coast it merges groups the
# coverage kept apart and returns fewer entries than the coverage started with:

DGG.compact(flat)

# `level_ranges` needs contiguous subtrees (`has_sorted_subtrees`); `expand`
# works on every system.
#
# A budget set records which of its cells were proven to lie inside the target,
# and `coarsest_contained` returns the shallowest of those. At ten cells every
# cell still straddles the state line, so the answer is `nothing`:

DGG.coarsest_contained(sets[1])

# By forty cells, one is inside:

DGG.coarsest_contained(sets[2])

# Cell areas over the coverage's leaves, in steradians. HEALPix gives the same
# number twice, so an unweighted mean over these cells is the areal mean:

g9 = DGG.levelgrid(sys, 9)
extrema(DGG.cell_area(g9, c) for c in flat)

# A budget set backs the same axis: `CellLookup` expands its forty cells to
# level-9 leaves.

DGG.CellLookup(sets[2]; level = 9)
