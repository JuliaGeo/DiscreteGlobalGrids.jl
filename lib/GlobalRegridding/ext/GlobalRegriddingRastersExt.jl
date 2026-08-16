"""
    GlobalRegriddingRastersExt

`Rasters.AbstractRaster`'s declared nodata sentinel, for
[`GlobalRegridding.sourcemissingval`](@ref).

A raster carries its sentinel in a **field**, not in metadata: a GeoTIFF's
nodata has already become `missingval` by the time the array exists, and
`DimensionalData.metadata` of a GDAL-backed raster carries only `"filepath"`.
The package's metadata route therefore answers `nothing` for every raster, and
this is what makes `regrid(raster; to)` see a sentinel the file declared without
the caller repeating it.
"""
module GlobalRegriddingRastersExt

import GlobalRegridding
import Rasters

"""
    GlobalRegridding.sourcemissingval(A::Rasters.AbstractRaster)

`Rasters.missingval(A)`, with both of its "nothing to compare against" spellings
mapped to `nothing`.

Rasters says a raster has no sentinel in two ways: `nothing`, and `missing` for
an array whose element type already carries `Missing`. Neither is a value to
test data against — `GlobalRegridding.isvalidvalue` rejects `missing` on its own
— and `nothing` is the spelling that compiles the comparison away, so an
array of `Union{Missing,Float32}` costs nothing per element for declaring one.
Every other value is passed through unchanged, including a `NaN` sentinel and an
integer one such as `-9999`.
"""
GlobalRegridding.sourcemissingval(A::Rasters.AbstractRaster) =
    _sentinel(Rasters.missingval(A))

_sentinel(mv) = mv
_sentinel(::Nothing) = nothing
_sentinel(::Missing) = nothing

end # module GlobalRegriddingRastersExt
