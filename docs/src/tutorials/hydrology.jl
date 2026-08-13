
ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))

using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
import ConservativeRegridding as CR
using Rasters, RasterDataSources
import ArchGDAL, GADM, NaturalEarth
import DimensionalData as DD, GeometryOps as GO, GeoInterface as GI, Extents
import GeometryOps: SpatialTreeInterface as STI
using Statistics, Dates
# using WGLMakie, GeoMakie, Makie
using GLMakie, GeoMakie, Makie
using ImageTransformations

sys = IGEO7DGGS()
DESTINATION_LEVEL = 13



# First understand which 1x1 tiles switzerland covers
ne_countries = NaturalEarth.naturalearth("admin_0_countries", 10)
switzerland_planar = ne_countries[findfirst(==("Switzerland"), ne_countries.NAME)].geometry
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
    # Descend into the nodes that intersect the box
    # Stop once a node is found that is covered by the box
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

cell = @time first_fitting_cell(sys, prep, boxcap)
level, id = node_level(cell), node_id(cell)
cp = cell_polygon_unitsphere(sys, level, id)
cp_longlat = GO.transform(GO.GeographicFromUnitSphere(), cp)

poly(cp_longlat)
poly!(box; alpha = 0.5)

# Just to check, that the children are in fact inside the box (they aren't, but let's ignore that :sweat_smile:)
DESTINATION_LEVEL = 13
lo, hi = descendant_range(sys, 5, id, DESTINATION_LEVEL)
ords = DGG.cell_to_ordinal(sys, DESTINATION_LEVEL, lo):DGG.cell_to_ordinal(sys, DESTINATION_LEVEL, hi)   # 9912026292:9912849834
cell_ids = DGG.ordinal_to_cell.((sys,), (DESTINATION_LEVEL,), ords)   # id of each
cell_polys = [cell_polygon_unitsphere(sys, level, id) for id in cell_ids] .|> x -> GO.transform(GO.GeographicFromUnitSphere(), x)

poly(cell_polys)
poly!(box; alpha = 0.5)

Rasters.checkmem!(false)

lookup = DGG.IGeo7.IGeo7Lookups.IGeo7Lookup(cell_ids; resolution = DESTINATION_LEVEL)
dem_igeo7_ras = Raster(zeros(length(lookup)), (Dim{:cells}(lookup),); name = :height)
dem_longlat_ras_points = Raster(joinpath(homedir(), "Downloads/Copernicus_DSM_10_N46_00_E010_00_DEM.tif"))
dem_longlat_ras = set(dem_longlat_ras_points, X => DD.Intervals(DD.Start()), Y => DD.Intervals(DD.Start()))

xintervalbounds = DD.intervalbounds(dims(dem_longlat_ras, X))
xgridcorners = vcat(first.(xintervalbounds), [last(last(xintervalbounds))])

yintervalbounds = DD.intervalbounds(dims(dem_longlat_ras, Y))
ygridcorners = vcat(first.(yintervalbounds), [last(last(yintervalbounds))])

longlat_cr_grid = CR.Trees.TopDownQuadtreeCursor(CR.Trees.CellBasedGrid(GO.Spherical(), GO.UnitSphereFromGeographic().(GI.Point.(xgridcorners, ygridcorners'))))

lookup = DGG.IGeo7.IGeo7Lookups.IGeo7Lookup(cell_ids, DESTINATION_LEVEL, Dict{String,Any}())
igeo7_cr_grid = DGG.DGGSPartialGrid(lookup; root_level = 5, root_id = id)

regridder = @time CR.Regridder(igeo7_cr_grid, longlat_cr_grid; threaded = true)
CR.regrid!(dem_igeo7_ras, regridder, dem_longlat_ras |> vec)

fap = poly(
    [GO.transform(GO.GeographicFromUnitSphere(), cell_polygon_unitsphere(sys, DESTINATION_LEVEL, id)) for id in DD.lookup(dem_igeo7_ras, Dim{:cells}())]; 
    color = vec(dem_igeo7_ras)
);
save("dem_igeo7_ras.png", fap)

using JLD2
# jldsave("dem_igeo7_ras.bin"; dem_igeo7_ras)
rl = load("dem_igeo7_ras.bin", "dem_igeo7_ras")

# `eachindex` yields canonical IGEO7 cell identities for a one-dimensional
# raster whose only dimension has an IGeo7Lookup.
x = eachindex(rl)
random_cell = rand(eachindex(rl))
random_height = rl[random_cell]

# The raster form omits neighbors outside the saved coverage. Use
# `DGG.neighbors(random_cell)` when all global neighbors are wanted.
stored_neighbors = DGG.neighbors(rl, random_cell)
neighbor_heights = [(cell = neighbor,
                     distance_m = DGG.celldistance(rl, random_cell, neighbor),
                     height = rl[neighbor])
                    for neighbor in stored_neighbors]

@info "Random IGEO7 cell and stored neighbors" random_cell random_height neighbor_heights

# Cells with at least one neighbor outside the saved raster coverage.
boundary_cells = DGG.edges(rl)



using Geomorphometry

# Extension
Geomorphometry.cellarea(r::Raster, i::DGG.IGEO7Index; cellsize) = DGG.cellarea(r, i)
Geomorphometry.neighbors(r::Raster, cell) = DGG.neighbors(r, cell)
Geomorphometry.celldistance(r::Raster, from, to; kwargs...) = DGG.celldistance(r, from, to)
Geomorphometry.outlets(r::Raster) = DGG.edges(r)

Geomorphometry.FlowDirection{C}(ci::RelativeIGEO7Index) where {C <: Geomorphometry.FlowDirectionConvention} =
    FlowDirection{C}(DGG.directioncode(ci))

tpi = topographic_position_index(rl)

slop = slope(rl)

hand = height_above_nearest_drainage(rl; threshold = 500*1000)

acc = flowaccumulation(rl, method=FD8())

polys = [GO.transform(GO.GeographicFromUnitSphere(), cell_polygon_unitsphere(sys, DESTINATION_LEVEL, id)) for id in DD.lookup(rl, Dim{:cells}())]; 
fap = poly(
    polys,
    color = vec(tpi),
    colorrange = (-5, 5),
    colormap = :delta,
)
save("tpi.png", fap, dpi=600);


scale = 1
width, height = 1500, 1500

fig = Figure(size = (scale*width, scale*height))
ax = Axis(fig[1, 1], title = "Flow accumulation on Ortler Massif | IGeo7 level 13", xlabel = "Longitude", ylabel = "Latitude")
fap = poly!(
    ax,
    polys,
    color = log10.(vec(acc)),
    colormap = :rain,
    strokewidth = 0,
);
Colorbar(
    fig[1, 2],
    fap,
    label = "Flow accumulation",
)
save("flow_accumulation_ortler.png", fig, px_per_unit = 2)

# img = GLMakie.colorbuffer(fig.scene)
# img_final = imresize(img, (height, width))
# save("slope_ortler_resized.png", img_final)
