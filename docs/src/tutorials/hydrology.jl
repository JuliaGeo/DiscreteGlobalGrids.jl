
ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))

using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
import ConservativeRegridding as CR
using Rasters, RasterDataSources
import ArchGDAL, GADM, NaturalEarth
import DimensionalData as DD, GeometryOps as GO, GeoInterface as GI, Extents
import GeometryOps: SpatialTreeInterface as STI
using Statistics, Dates
using WGLMakie, GeoMakie, Makie

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
dem_longlat_ras_points = Raster("/Users/anshul/Downloads/Copernicus_DSM_10_N46_00_E010_00_DEM.tif")
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

# fap = poly(
#     [GO.transform(GO.GeographicFromUnitSphere(), cell_polygon_unitsphere(sys, DESTINATION_LEVEL, id)) for id in DD.lookup(dem_igeo7_ras, Dim{:cells}())]; 
#     color = vec(dem_igeo7_ras)
# );
# save("dem_igeo7_ras.png", fap)