"""
    GlobalRegriddingRastersProjExt

CRS-derived charts for `GlobalRegridding.RasterGrid`.

With Rasters and Proj both loaded, a raster whose X/Y lookups carry a CRS needs
no `native_to_unit_sphere` keyword: geographic coordinates are normalized to
Greenwich degrees without changing the datum or latitude convention. A
projected CRS becomes
`Proj.Transformation(crs, "EPSG:4326"; always_xy = true)`, which the Proj
extension turns into task-local forward and inverse charts.
"""
module GlobalRegriddingRastersProjExt

import GlobalRegridding
import GlobalRegridding as GR
import Proj
import Rasters
import Rasters.GeoFormatTypes as GFT

"""
    GlobalRegridding._projected_crs_native_to_unit_sphere(crs::GeoFormat)

Normalize geographic coordinates to Greenwich degrees, preserving their datum
and latitude convention. Projected CRSs use a transformation to EPSG:4326.
"""
function GlobalRegridding._projected_crs_native_to_unit_sphere(crs::GFT.GeoFormat)
    definition = crs isa GFT.ProjString ?
        GFT.ProjString(GFT.val(crs) * " +type=crs") : crs
    source = Proj.CRS(definition)
    while Proj.is_bound(source)
        native = Proj.proj_get_source_crs(source)
        native == C_NULL && throw(ArgumentError("cannot read the bound CRS source"))
        source = Proj.CRS(native)
    end
    if Proj.is_geographic(source)
        return _geographic_chart(source)
    end
    return Proj.Transformation(convert(String, crs), "EPSG:4326"; always_xy = true)
end

struct _GeographicChart
    xscale::Float64
    yscale::Float64
    offset::Float64
end

struct _GeographicInverse
    forward::_GeographicChart
end

@inline (t::_GeographicChart)(xy) = GR.US.UnitSphereFromGeographic()(
    (t.xscale * xy[1] + t.offset, t.yscale * xy[2]))

@inline function (t::_GeographicInverse)(point)
    lon, lat = GR.US.GeographicFromUnitSphere()(point)
    return ((lon - t.forward.offset) / t.forward.xscale, lat / t.forward.yscale)
end

GR._default_unit_sphere_to_native(t::_GeographicChart) = _GeographicInverse(t)
GR._native_coordinate_limits(t::_GeographicChart) =
    (360 / abs(t.xscale), (-90 / abs(t.yscale), 90 / abs(t.yscale)))
GR._spherical_step_bounds_radians(t::_GeographicChart, dx, dy) =
    (deg2rad(abs(t.xscale) * Float64(dx)), deg2rad(abs(t.yscale) * Float64(dy)))

function _geographic_chart(source)
    geographic = Proj.proj_crs_get_geodetic_crs(source)
    geographic == C_NULL && throw(ArgumentError("cannot read the geographic CRS"))
    meridian = coordinate_system = C_NULL
    try
        meridian = Proj.proj_get_prime_meridian(geographic)
        coordinate_system = Proj.proj_crs_get_coordinate_system(geographic)
        (meridian == C_NULL || coordinate_system == C_NULL) &&
            throw(ArgumentError("cannot read the geographic CRS axes and prime meridian"))
        longitude, unit = Ref{Cdouble}(), Ref{Cdouble}()
        Proj.proj_prime_meridian_get_parameters(meridian, longitude, unit, C_NULL) == 1 ||
            throw(ArgumentError("cannot read the geographic CRS prime meridian"))
        offset = rad2deg(longitude[] * unit[])
        xscale = yscale = NaN
        for axis in 0:(Proj.proj_cs_get_axis_count(coordinate_system) - 1)
            direction = Ref{Cstring}()
            Proj.proj_cs_get_axis_info(coordinate_system, axis, C_NULL, C_NULL,
                direction, unit, C_NULL, C_NULL, C_NULL) == 1 ||
                throw(ArgumentError("cannot read the geographic CRS angular units"))
            orientation = unsafe_string(direction[])
            scale = rad2deg(unit[])
            orientation == "east" && (xscale = scale)
            orientation == "west" && (xscale = -scale)
            orientation == "north" && (yscale = scale)
            orientation == "south" && (yscale = -scale)
        end
        (isfinite(xscale) && isfinite(yscale)) || throw(ArgumentError(
            "geographic raster axes must describe longitude and latitude"))
        if iszero(offset) && xscale == 1 && yscale == 1
            return GlobalRegridding._UNSET_RASTER_KEYWORD
        end
        return _GeographicChart(xscale, yscale, offset)
    finally
        coordinate_system == C_NULL || Proj.proj_destroy(coordinate_system)
        meridian == C_NULL || Proj.proj_destroy(meridian)
        Proj.proj_destroy(geographic)
    end
end

end # module GlobalRegriddingRastersProjExt
