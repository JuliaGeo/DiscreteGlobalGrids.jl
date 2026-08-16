"""
    GlobalRegriddingRastersExt

Adds `Rasters.AbstractRaster` nodata support to
[`GlobalRegridding.sourcemissingval`](@ref).
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

end # module GlobalRegriddingRastersExt
