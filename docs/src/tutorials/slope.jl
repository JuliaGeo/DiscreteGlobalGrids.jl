# # Terrain slope on the interior of a DGGS tile
#
# This continues `hydrology.jl`: same box over the Swiss Alps, same Copernicus
# DEM regridded onto the same IGEO7 subtree. What it adds is the reason to care
# which cells of a tile you are allowed to answer for.
#
# Slope at a cell is computed from its *neighbours*, so a cell whose
# neighbourhood is cut off by the tile edge does not have a slope — it has a
# slope estimated from a subset, biased toward whatever side of the cell
# survived. On a steepest-descent measure the bias has a direction: drop the
# downhill neighbour and the cell reads flat, which is indistinguishable from a
# genuine pit. Flow routing then terminates there, and the basin above it
# drains into a hole that is an artefact of the tiling.
#
# `edge_cells` and `interior_cells` name those two sets, so the sweep runs on
# the interior and leaves the rim undefined rather than quietly wrong.
#
# Not part of the docs build: it reads a DEM from a local path.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))

using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
import ConservativeRegridding as CR
using Rasters
import ArchGDAL, NaturalEarth
import DimensionalData as DD, GeometryOps as GO, GeoInterface as GI, Extents
import GeometryOps: SpatialTreeInterface as STI
using Statistics
using CairoMakie, GeoMakie
CairoMakie.activate!()

# The DEM tile that covers the box, and how fine to go. Level 11 puts roughly
# 25,800 m² in a cell (~160 m across) against the DEM's ~30 m posting, so each
# cell averages tens of DEM pixels — fine enough for alpine relief, coarse
# enough to stay a quick run. Level 12 is 3,700 m² and takes about 7x longer.
const DEM_PATH = expanduser("~/Downloads/Copernicus_DSM_10_N46_00_E010_00_DEM.tif")
const DESTINATION_LEVEL = 11
const R_EARTH = 6_371_008.8

# ## The tile
#
# Same descent as `hydrology.jl`: find the coarsest single IGEO7 cell that
# covers the box, and take its descendants at `DESTINATION_LEVEL`.

ext = Extents.Extent(X = (10, 11), Y = (46, 47))
box = GI.Polygon([GI.LinearRing([(ext.X[1], ext.Y[1]), (ext.X[2], ext.Y[1]),
                                 (ext.X[2], ext.Y[2]), (ext.X[1], ext.Y[2]),
                                 (ext.X[1], ext.Y[1])])])
prep = GO.prepare(GO.RelateNG(; manifold = GO.Spherical()), box)
center = GO.UnitSphereFromGeographic()((sum(ext.X) / 2, sum(ext.Y) / 2))
boxcap = GO.UnitSpherical.SphericalCap(center,
    maximum(GO.UnitSpherical.spherical_distance(center, GO.UnitSphereFromGeographic()(p))
            for p in GI.getpoint(box)))

function first_fitting_cell(sys, box_prep, box_cap; leaf = 12)
    nodes = [treeify(DGGSGrid(sys, leaf))]
    while !isempty(nodes)
        for node in nodes
            node_level(node) >= 0 || continue  # the root is synthetic: no cell
            GO.relate_predicate(box_prep, GO.pred_covers(),
                cell_polygon_unitsphere(sys, node_level(node), node_id(node))) && return node
        end
        nodes = [child for node in nodes for child in STI.getchild(node)
                    if intersects_cap(box_cap, STI.node_extent(child))]
    end
    return nothing
end

sys = IGEO7DGGS()
cell = first_fitting_cell(sys, prep, boxcap)
level, id = node_level(cell), node_id(cell)

# `DGGSSubtreeIds` *is* that subtree, lazily — the cell ids in ascending order,
# without materializing them.

tile = DGGSSubtreeIds(sys, level, id, DESTINATION_LEVEL)
cell_ids = collect(tile)
length(tile)

# ## The DEM, conservatively regridded onto it

Rasters.checkmem!(false)

dem_longlat_ras_points = Raster(DEM_PATH)
dem_longlat_ras = set(dem_longlat_ras_points, X => DD.Intervals(DD.Start()),
                                              Y => DD.Intervals(DD.Start()))

xintervalbounds = DD.intervalbounds(dims(dem_longlat_ras, X))
xgridcorners = vcat(first.(xintervalbounds), [last(last(xintervalbounds))])
yintervalbounds = DD.intervalbounds(dims(dem_longlat_ras, Y))
ygridcorners = vcat(first.(yintervalbounds), [last(last(yintervalbounds))])

longlat_cr_grid = CR.Trees.TopDownQuadtreeCursor(
    CR.Trees.CellBasedGrid(GO.Spherical(),
        GO.UnitSphereFromGeographic().(GI.Point.(xgridcorners, ygridcorners'))))

lookup = DGG.IGeo7.IGeo7Lookups.IGeo7Lookup(cell_ids, DESTINATION_LEVEL, Dict{String,Any}())
igeo7_cr_grid = DGG.DGGSPartialGrid(lookup; root_level = level, root_id = id)

dem = Raster(zeros(length(tile)), (Dim{:cells}(lookup),); name = :height)
regridder = @time CR.Regridder(igeo7_cr_grid, longlat_cr_grid; threaded = true)
@time CR.regrid!(dem, regridder, dem_longlat_ras |> vec)

heights = parent(dem)
extrema(heights)

# ## Slope, on the interior only
#
# `interior_cells(tile)` is the complement of the rim, allocating nothing: at
# this level it is 115,465 of the tile's 117,649 cells, and materializing it
# would cost ~0.9 MiB to say "almost all of them".
#
# The neighbour stepper hands back tile *positions*, which index `heights`
# directly. Since IGEO7's `cell_neighbors` is GBT digit arithmetic, this whole
# sweep touches no geometry beyond the cell centres it needs for distances.

edge = edge_cells(tile)
interior = interior_cells(tile, edge)
centers = [DGG.cell_center(sys, DESTINATION_LEVEL, cid) for cid in tile]
stepper = neighbor_stepper(tile)

# The tile is a rosette and the DEM is a 1°x1° square, so the tile pokes out
# past the data on every side — and the cells that straddle the edge are the
# dangerous ones. A conservative regrid divides by the *whole* cell's area, so a
# cell 5% covered by 2000 m terrain reports 100 m. That is not a missing value
# to test for; it is a plausible-looking elevation, and it manufactures an 85°
# cliff along the entire DEM boundary.
#
# The regridder already knows the answer: push a field of ones through it and
# each cell reports the fraction of itself the source covers. Anything short of
# full coverage is dropped, along with its neighbours — for the same reason the
# rim is, an estimate from a neighbourhood you do not fully have is not a slope.

coverage = Raster(zeros(length(tile)), (Dim{:cells}(lookup),); name = :coverage)
CR.regrid!(coverage, regridder, ones(length(dem_longlat_ras)))
covered = parent(coverage) .> 1 - 1e-6
count(!, covered)

# The neighbour stepper hands back tile *positions*, which index `heights`
# directly. Since IGEO7's `cell_neighbors` is GBT digit arithmetic, this whole
# sweep touches no geometry beyond the cell centres it needs for distances.

"Steepest-descent slope in degrees; `NaN` wherever the neighbourhood is incomplete."
function steepest_descent_slope(tile, heights, centers, stepper, interior, covered)
    slope = fill(NaN, length(tile))
    for i in interior
        covered[i] || continue
        hi = heights[i]
        best = 0.0
        complete = true
        for j in step_neighbors(stepper, i)
            if !covered[j]
                complete = false
                break
            end
            drop = hi - heights[j]
            drop > 0 || continue
            d = R_EARTH * GO.UnitSpherical.spherical_distance(centers[i], centers[j])
            g = drop / d
            g > best && (best = g)
        end
        complete && (slope[i] = atand(best))
    end
    return slope
end

slope = @time steepest_descent_slope(tile, heights, centers, stepper, interior, covered)
valid = [i for i in interior if !isnan(slope[i])]

# A cell with no downhill neighbour reads as slope 0 — a pit. On the interior
# these are real (lakes, valley floors, DEM noise). Had we swept the rim too,
# every edge cell whose downhill side lies outside the tile would have joined
# them, and a flow-routing pass could not tell the two apart.

slopes = slope[valid]
(; cells = length(tile), rim = length(edge), interior = length(interior),
   partial_coverage = length(interior) - length(valid), scored = length(valid),
   pits = count(iszero, slopes), median_slope = median(slopes),
   max_slope = maximum(slopes))

# ## The map

polys = [GO.transform(GO.GeographicFromUnitSphere(),
                      cell_polygon_unitsphere(sys, DESTINATION_LEVEL, tile[i]))
         for i in valid]

fig = Figure(size = (900, 800))
ax = GeoAxis(fig[1, 1];
    title = "Steepest-descent slope, IGEO7 level $DESTINATION_LEVEL",
    subtitle = "$(length(valid)) cells scored; $(length(edge)) rim + " *
               "$(length(interior) - length(valid)) partially covered omitted",
    dest = "+proj=ortho +lon_0=10.5 +lat_0=46.5")
p = poly!(ax, polys; color = slopes, colormap = :batlow, strokewidth = 0)
Colorbar(fig[1, 2], p; label = "slope (°)")
save(joinpath(@__DIR__, "slope_switzerland.png"), fig)
fig
