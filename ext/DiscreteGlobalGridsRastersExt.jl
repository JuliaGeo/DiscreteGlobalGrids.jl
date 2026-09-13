module DiscreteGlobalGridsRastersExt
import DiscreteGlobalGrids as DGG
import Rasters as RA
import DimensionalData as DD

DGG._raster_missingval(A::RA.AbstractRaster) = RA.missingval(A)
DGG._raster_rebuild(template::RA.AbstractRaster, data, dims; name=nothing, metadata=nothing, kw...) =
    RA.rebuild(template; data, dims, name=something(name, DD.name(template)), metadata=something(metadata, DD.metadata(template)), kw...)
DGG._raster_stack(template::RA.AbstractRaster, layers::NamedTuple; metadata=nothing) =
    RA.RasterStack(layers; metadata=something(metadata, DD.metadata(template)), refdims=DD.refdims(template))
end
