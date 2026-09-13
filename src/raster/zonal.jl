"""
    DiscreteGlobalGrids.zonal(f, A; of, spatialslices=true, skipmissing=true, kw...)

Reduce the selected cells of each zone; selection runs once per zone.

- `spatialslices=true` reduces the cell dimension for every other slice;
  `false` reduces the whole selected cube; a tuple names the dimensions to
  reduce and must include the cell dimension.
- Several zones append `Dim{:Zone}`; a single geometry returns its result.
- Zones outside the holding return `missing`. `emptyval` replaces empty
  slices; without it `f` receives an empty iterator.
- Means are unweighted cell means.
"""
function zonal(f, A::Union{DD.AbstractDimArray,DD.AbstractDimStack}; of,
        spatialslices=true, skipmissing=true, emptyval=_RASTER_UNSET,
        geometrycolumn=nothing, boundary=:center, shape=nothing,
        threaded=true, progress=true, kw...)
    _raster_check_keywords(kw)
    boundary = _raster_boundary(boundary)
    shape = _raster_shape(shape)
    geoms = _raster_geometries(of; geometrycolumn)
    single = _raster_issingle(of)
    if A isa DD.AbstractDimStack
        layers = _raster_eachlayer(A, geoms; boundary, shape) do _, layer, selections
            _raster_zonal(f, layer, geoms, selections, single; spatialslices, skipmissing, emptyval, shape, threaded)
        end
        return all(x -> x isa DD.AbstractDimArray, layers) ? _raster_stack(A, layers) : layers
    end
    selections = _raster_selections(_raster_grid(A), geoms; boundary, shape)
    return _raster_zonal(f, A, geoms, selections, single; spatialslices, skipmissing, emptyval, shape, threaded)
end

# Like VectorDataCubes' SpatialSliceify, applied to slices selected once per zone.
struct _RasterSpatialSliceify{F,E,M}
    f::F
    emptyval::E
    missingval::M
    skipmissing::Bool
end
function (s::_RasterSpatialSliceify)(x)
    vals = s.skipmissing ? Iterators.filter(v -> !_raster_ismissing(v, s.missingval), x) : x
    s.emptyval isa _RasterUnset || !isempty(vals) || return _raster_copy(s.emptyval)
    return s.f(vals)
end

function _raster_zonal(f, A, geoms, selections, single; spatialslices, skipmissing, emptyval, shape, threaded)
    d = _raster_dimnum(A)
    grid = _raster_grid(A)
    celldim = DD.dims(A)[d]
    reduced = spatialslices === true ? (celldim,) : spatialslices === false ? DD.dims(A) :
        spatialslices isa Tuple ? DD.dims(A, spatialslices) :
        throw(ArgumentError("spatialslices must be true, false, or a tuple of dimensions"))
    DD.hasdim(reduced, celldim) || throw(ArgumentError("spatialslices must include the cell dimension"))
    kept = DD.otherdims(A, reduced)
    reducer = _RasterSpatialSliceify(f, emptyval, _raster_missingval(A), skipmissing)
    zones = _raster_map(eachindex(geoms), threaded) do j
        indices = selections[j]
        # Empty membership is a sub-cell zone, unless the zone lies outside this holding.
        if isempty(indices) && !_raster_empty_geometry(geoms[j]) &&
                isempty(_raster_indices(grid, geoms[j]; boundary=:intersects, shape))
            return isempty(kept) ? missing : fill(missing, map(length, kept))
        end
        zone = view(A, ntuple(i -> i == d ? indices : Colon(), ndims(A))...)
        isempty(kept) && return reducer(zone)
        return map(reducer, eachslice(zone; dims=kept))
    end
    if single
        z = only(zones)
        return isempty(kept) ? z : _raster_zonal_rebuild(A, z isa DD.AbstractDimArray ? parent(z) : z, kept)
    end
    zonedim = DD.Dim{:Zone}(Base.OneTo(length(geoms)))
    if isempty(kept)
        return spatialslices === false ? zones : _raster_zonal_rebuild(A, zones, (zonedim,))
    end
    arrays = map(z -> z isa DD.AbstractDimArray ? parent(z) : z, zones)
    data = isempty(arrays) ? Array{Missing}(undef, map(length, kept)..., 0) : Base.stack(arrays)
    return _raster_zonal_rebuild(A, data, (kept..., zonedim))
end

_raster_zonal_rebuild(A, data, dims) =
    _raster_rebuild(A, data, dims; missingval=Missing <: eltype(data) ? missing : nothing)
