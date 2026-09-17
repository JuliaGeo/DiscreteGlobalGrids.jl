"""
    GlobalRegriddingRastersExt

Teaches regridding the `Rasters.AbstractRaster` nodata convention on both
sides: [`GlobalRegridding.sourcemissingval`](@ref) reads a raster's sentinel,
[`GlobalRegridding.outputmissingval`](@ref) and
[`GlobalRegridding.destinationmissingval`](@ref) carry it to the destination,
and [`GlobalRegridding.rebuildoutput`](@ref) returns a raster declaring it.

It also reads CRS metadata off `Projected` and `Mapped` X/Y lookups so
[`GlobalRegridding.RasterGrid`](@ref) charts a projected raster through its
CRS. Proj classifies and transforms any CRS this extension cannot recognize as
geographic on its own; see `GlobalRegriddingRastersProjExt`.
"""
module GlobalRegriddingRastersExt

import GlobalRegridding
import DimensionalData as DD
import Rasters
import Rasters.GeoFormatTypes as GFT

"""
    GlobalRegridding.sourcemissingval(A::Rasters.AbstractRaster)

Return `Rasters.missingval(A)`, mapping `nothing` and `missing` to `nothing`.
Other sentinels, including NaN and integer values, pass through unchanged.
"""
GlobalRegridding.sourcemissingval(A::Rasters.AbstractRaster) =
    _sentinel(Rasters.missingval(A))

_sentinel(mv) = mv
_sentinel(::Nothing) = nothing
_sentinel(::Missing) = nothing

"""
    GlobalRegridding.outputmissingval(A::Rasters.AbstractRaster)

Return `Rasters.missingval(A)`: a regrid keeps the nodata convention its source
declares, so a raster of `Union{Missing,Float64}` comes back holding `missing`
and one with a `-9999.0` sentinel comes back holding `-9999.0`. A raster
declaring no sentinel falls back to the element type's own blank.

Pass `missingval` to [`GlobalRegridding.regrid`](@ref) to choose another one —
`missingval = NaN` is what makes a `Union{Missing,Float64}` source regrid into
a concrete `Float64` raster.
"""
GlobalRegridding.outputmissingval(A::Rasters.AbstractRaster) =
    _declared(Rasters.missingval(A),
        GlobalRegridding._maskedvalue(GlobalRegridding.outputeltype(eltype(A))))

"""
    GlobalRegridding.destinationmissingval(A::Rasters.AbstractRaster)

Return `Rasters.missingval(A)`, the sentinel [`GlobalRegridding.regrid!`](@ref)
blanks a preallocated raster with, falling back to the element type's own blank
when the raster declares none.
"""
GlobalRegridding.destinationmissingval(A::Rasters.AbstractRaster) =
    _declared(Rasters.missingval(A), GlobalRegridding._elementblank(eltype(A)))

_declared(mv, fallback) = mv
_declared(::Nothing, fallback) = fallback

"""
    GlobalRegridding.rebuildoutput(data::Rasters.AbstractRaster, out, dims, missingval)

Return a raster over `out` labelled with `dims`, declaring `missingval` and
keeping `data`'s name, metadata and reference dimensions.
"""
GlobalRegridding.rebuildoutput(data::Rasters.AbstractRaster, out::AbstractArray,
    dims::Tuple, missingval) =
    Rasters.rebuild(data; data = out, dims, missingval)

# CRS-derived charts

"""
    GlobalRegridding._crs_native_to_unit_sphere(A)

Return the forward chart implied by the CRS of `A`'s X/Y lookups, for a raster,
dimensional array or dimension tuple. Lookups without CRS metadata and CRS
values recognized as geographic return `_UNSET_RASTER_KEYWORD`; any other CRS
goes to [`GlobalRegridding._projected_crs_native_to_unit_sphere`](@ref).
"""
GlobalRegridding._crs_native_to_unit_sphere(A::DD.AbstractDimArray) =
    _lookup_crs_native_to_unit_sphere(DD.dims(A))
GlobalRegridding._crs_native_to_unit_sphere(ds::Tuple{Vararg{DD.Dimension}}) =
    _lookup_crs_native_to_unit_sphere(ds)

function _lookup_crs_native_to_unit_sphere(ds)
    crs = _valuecrs(ds)
    crs === nothing && return GlobalRegridding._UNSET_RASTER_KEYWORD
    _geographic_without_proj(crs) && return GlobalRegridding._UNSET_RASTER_KEYWORD
    return GlobalRegridding._projected_crs_native_to_unit_sphere(crs)
end

"""
    _valuecrs(dims) -> GeoFormat or nothing

Return the CRS the X and Y lookup values are expressed in. `Projected` values
live in `crs`; `Mapped` values live in `mappedcrs`, with `crs` describing the
projection the mapped values came from. Both lookups must agree when both
carry one.
"""
function _valuecrs(ds::Tuple)
    xcrs = _valuecrs(DD.dims(ds, DD.XDim))
    ycrs = _valuecrs(DD.dims(ds, DD.YDim))
    xcrs === nothing && return ycrs
    ycrs === nothing && return xcrs
    xcrs == ycrs || throw(ArgumentError(
        "RasterGrid: the X lookup is in $(repr(xcrs)) and the Y lookup in " *
        "$(repr(ycrs)); both spatial lookups must share one CRS"))
    return xcrs
end

# An unformatted dimension holds a bare vector, which carries no CRS.
_valuecrs(::Any) = nothing
_valuecrs(d::DD.Dimension) = _valuecrs(DD.lookup(d))
_valuecrs(l::Rasters.Projected) = Rasters.crs(l)
function _valuecrs(l::Rasters.Mapped)
    valuecrs = Rasters.mappedcrs(l)
    valuecrs === nothing && Rasters.crs(l) !== nothing && throw(ArgumentError(
        "RasterGrid: the Mapped lookup names the projected CRS " *
        "$(repr(Rasters.crs(l))) but no `mappedcrs`, so the CRS of its values is " *
        "unknown. Set it with `Rasters.setmappedcrs`, or pass " *
        "`native_to_unit_sphere` explicitly."))
    return valuecrs
end

"""
    _geographic_without_proj(crs::GeoFormat) -> Bool

Return whether `crs` is recognizably geographic from its GeoFormatTypes value
alone: `EPSG(4326)`, a PROJ string whose `+proj` is a longitude/latitude
family, or WKT whose root node is a geographic or geodetic CRS. Geographic
here means longitude/latitude in degrees from Greenwich. WKT2 `GEODCRS` also
covers geocentric CRSs, which this heuristic accepts. Anything else is left
for Proj to classify.
"""
_geographic_without_proj(crs::GFT.EPSG) = crs.val == (4326,)
_geographic_without_proj(crs::GFT.ProjString) =
    occursin(r"\+proj=(longlat|latlong|latlon|lonlat)\b", GFT.val(crs))
_geographic_without_proj(crs::GFT.AbstractWellKnownText) =
    startswith(lstrip(GFT.val(crs)), r"(GEOGCS|GEOGCRS|GEODCRS|GEODETICCRS|GEOGRAPHICCRS)\[")
_geographic_without_proj(::GFT.GeoFormat) = false

end # module GlobalRegriddingRastersExt
