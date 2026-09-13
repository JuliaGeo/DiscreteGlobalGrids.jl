module DiscreteGlobalGridsRastersExt
import DiscreteGlobalGrids as DGG
import Rasters as RA
import DimensionalData as DD

DGG._raster_missingval(A::RA.AbstractRaster) = RA.missingval(A)
function DGG._raster_rebuild(template::RA.AbstractRaster, data, dims;
        name=nothing, metadata=nothing, missingval=DGG._RASTER_UNSET)
    return RA.rebuild(template; data, dims, name=something(name, DD.name(template)),
        metadata=something(metadata, DD.metadata(template)),
        missingval=missingval isa DGG._RasterUnset ? RA.missingval(template) : missingval)
end
function DGG._raster_stack(template::Union{RA.AbstractRasterStack,RA.AbstractRaster}, layers::NamedTuple; metadata=nothing)
    return RA.RasterStack(layers; metadata=something(metadata, DD.metadata(template)), refdims=DD.refdims(template))
end
end
