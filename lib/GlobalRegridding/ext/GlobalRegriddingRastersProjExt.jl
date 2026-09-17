"""
    GlobalRegriddingRastersProjExt

CRS-derived projected charts for `GlobalRegridding.RasterGrid`.

With Rasters and Proj both loaded, a raster whose X/Y lookups carry a CRS needs
no `native_to_unit_sphere` keyword: Proj classifies the CRS, a geographic one
keeps the built-in longitude/latitude chart, and a projected one becomes
`Proj.Transformation(crs, "EPSG:4326"; always_xy = true)`, which the Proj
extension turns into task-local forward and inverse charts.
"""
module GlobalRegriddingRastersProjExt

import GlobalRegridding
import Proj
import Rasters
import Rasters.GeoFormatTypes as GFT

"""
    GlobalRegridding._projected_crs_native_to_unit_sphere(crs::GeoFormat)

Return `_UNSET_RASTER_KEYWORD` when Proj classifies `crs` as geographic, and
otherwise a `Proj.Transformation` from `crs` to EPSG:4326 in X/Y order.
"""
function GlobalRegridding._projected_crs_native_to_unit_sphere(crs::GFT.GeoFormat)
    definition = convert(String, crs)
    Proj.is_geographic(Proj.CRS(definition)) &&
        return GlobalRegridding._UNSET_RASTER_KEYWORD
    return Proj.Transformation(definition, "EPSG:4326"; always_xy = true)
end

end # module GlobalRegriddingRastersProjExt
