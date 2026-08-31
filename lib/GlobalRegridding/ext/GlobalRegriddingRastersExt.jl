"""
    GlobalRegriddingRastersExt

Teaches regridding the `Rasters.AbstractRaster` nodata convention on both
sides: [`GlobalRegridding.sourcemissingval`](@ref) reads a raster's sentinel,
[`GlobalRegridding.outputmissingval`](@ref) and
[`GlobalRegridding.destinationmissingval`](@ref) carry it to the destination,
and [`GlobalRegridding.rebuildoutput`](@ref) returns a raster declaring it.
"""
module GlobalRegriddingRastersExt

import GlobalRegridding
import Rasters

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

end # module GlobalRegriddingRastersExt
