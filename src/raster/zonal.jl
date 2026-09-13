"""
    DiscreteGlobalGrids.zonal(f, A; of, spatialslices=true, skipmissing=true, kw...)

Reduce the selected cells of each zone; selection runs once per zone.

- `spatialslices`: `true` reduces the cell dimension per other slice, `false`
  the whole selected cube, a tuple those dimensions (cell dimension required).
- `spatialslices=false` returns a `Vector` (one entry per zone) for an array
  and a `NamedTuple` of such vectors for a stack.
- Several zones append `Dim{:Zone}`; a single geometry returns its result.
- Zones outside the holding give `missing`; `emptyval` replaces empty slices,
  else `f` receives an empty iterator.
- Means are unweighted cell means.
"""
function zonal(f, A::Union{DD.AbstractDimArray,DD.AbstractDimStack}; of,
        spatialslices=true, skipmissing=true, emptyval=_RASTER_UNSET, threaded=true, kw...)
    geoms, opts = _raster_inputs(of; kw...)
    single = _raster_issingle(of)
    if A isa DD.AbstractDimStack
        layers = _raster_eachlayer(A, geoms; opts.boundary, opts.shape) do _, layer, selections
            _raster_zonal(f, layer, geoms, selections, single; spatialslices, skipmissing, emptyval, opts.shape, threaded)
        end
        return all(x -> x isa DD.AbstractDimArray, layers) ? _raster_stack(A, layers) : layers
    end
    selections = _raster_selections(_raster_grid(A), geoms; opts.boundary, opts.shape)
    return _raster_zonal(f, A, geoms, selections, single; spatialslices, skipmissing, emptyval, opts.shape, threaded)
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
    mv = _raster_missingval(A)
    # One slice: drop missing values, then reduce or hand back `emptyval`.
    function reducer(x)
        vals = skipmissing ? Iterators.filter(v -> !_raster_ismissing(v, mv), x) : x
        return emptyval isa _RasterUnset || !isempty(vals) ? f(vals) : _raster_copy(emptyval)
    end
    # Only a partial holding can place a zone outside the data, so only then is
    # the second `:intersects` selection worth its cost.
    partial = _raster_ispartial(grid)
    zones = _raster_map(eachindex(geoms), threaded) do j
        indices = selections[j]
        # Empty membership is a sub-cell zone, unless the zone lies outside this holding.
        if isempty(indices) && partial && !_raster_empty_geometry(geoms[j]) &&
                isempty(_raster_indices(grid, geoms[j]; boundary=:intersects, shape))
            return isempty(kept) ? missing : fill(missing, map(length, kept))
        end
        zone = view(A, DD.rebuild(celldim, indices))
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
