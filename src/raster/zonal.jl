"""
    DiscreteGlobalGrids.zonal(f, A; of, spatialslices=true, skipmissing=true, kw...)

Reduce selected cells of each zone. By default, reduce the cell dimension and
preserve every other labelled dimension; multiple zones append `Dim{:Zone}`.
`spatialslices=false` reduces the whole selected cube; a tuple names dimensions
to reduce and must include the cell dimension. Select geometry once per zone.
`emptyval` is applied independently to empty slices. Without it, `f` receives
an empty iterator. Means are unweighted cell means.
"""
function zonal(f, A::Union{DD.AbstractDimArray,DD.AbstractDimStack}; of,
        spatialslices=true, skipmissing=true, emptyval=_RASTER_UNSET,
        geometrycolumn=nothing, boundary=:center, shape=nothing,
        threaded=true, progress=true, kw...)
    _raster_check_keywords(kw)
    boundary = _raster_boundary(boundary)
    shape = _raster_shape(shape)
    feats = _raster_features(of;geometrycolumn)
    single = _raster_issingle(of)
    if A isa DD.AbstractDimStack
        selections = Dict{Any,Vector{Vector{Int}}}()
        layers = map(keys(A)) do k
            layer = A[k]
            lookup = DD.lookup(layer,_raster_dimnum(layer))
            hits = get!(selections,lookup) do
                grid, _, _ = _raster_target(layer)
                [_raster_indices(grid,feat.geom;boundary,shape) for feat in feats]
            end
            _raster_zonal(f,layer,feats,single;spatialslices,skipmissing,emptyval,boundary,shape,threaded,selections=hits)
        end
        nt = NamedTuple{keys(A)}(layers)
        return all(x -> x isa DD.AbstractDimArray,layers) ? _raster_stack(A,nt) : nt
    end
    return _raster_zonal(f,A,feats,single;spatialslices,skipmissing,emptyval,boundary,shape,threaded)
end

# Like VectorDataCubes' SpatialSliceify, but selection is outside the slice loop.
struct _RasterSpatialSliceify{F,E,M}
    f::F
    emptyval::E
    missingval::M
    skipmissing::Bool
end
function (s::_RasterSpatialSliceify)(x)
    vals = s.skipmissing ? Iterators.filter(v -> !_raster_ismissing(v,s.missingval),x) : x
    if !(s.emptyval isa _RasterUnset) && isempty(vals)
        return deepcopy(s.emptyval)
    end
    return s.f(vals)
end
function _raster_zonal(f,A,feats,single;spatialslices,skipmissing,emptyval,boundary,shape,threaded,selections=nothing)
    d = _raster_dimnum(A)
    grid, dims, _ = _raster_target(A)
    reduceaxes = if spatialslices === false
        Tuple(1:ndims(A))
    elseif spatialslices === true
        (d,)
    elseif spatialslices isa Tuple
        reduceddims = DD.dims(A,spatialslices)
        Tuple(findfirst(==(dim),dims) for dim in reduceddims)
    else
        throw(ArgumentError("spatialslices must be true, false, or a tuple of dimensions"))
    end
    d in reduceaxes || throw(ArgumentError("spatialslices must include the cell dimension"))
    kept = Tuple(i for i in 1:ndims(A) if !(i in reduceaxes))
    otherdims = Tuple(dims[i] for i in kept)
    outshape = Tuple(length(dim) for dim in otherdims)
    zones = Vector{Any}(undef,length(feats))
    reducer = _RasterSpatialSliceify(f,emptyval,_raster_missingval(A),skipmissing)
    function onezone(j)
        indices = selections === nothing ? _raster_indices(grid,feats[j].geom;boundary,shape) : selections[j]
        # Empty membership may mean a sub-cell zone, or a zone outside this
        # regional holding altogether. Sub-cell and genuinely empty geometries
        # call f on empty data; outside-domain geometries return missing.
        if isempty(indices) && !_raster_empty_geometry(feats[j].geom) &&
                isempty(_raster_indices(grid,feats[j].geom;boundary=:intersects,shape))
            zones[j] = isempty(kept) ? missing : fill(missing,outshape)
            return
        end
        # Index in one operation so labels of custom reducers remain available.
        zoneview = view(A,ntuple(i -> i==d ? indices : Colon(),ndims(A))...)
        if isempty(kept)
            zones[j] = reducer(zoneview)
        else
            zoneresult = Array{Any}(undef,outshape)
            for I in CartesianIndices(zoneresult)
                inds = ntuple(ndims(A)) do i
                    k = findfirst(==(i),kept)
                    k === nothing ? Colon() : I[k]
                end
                zoneresult[I] = reducer(view(zoneview,inds...))
            end
            zones[j] = _raster_narrow(zoneresult)
        end
    end
    if threaded && Threads.nthreads()>1
        Threads.@threads for j in eachindex(feats)
            onezone(j)
        end
    else
        foreach(onezone,eachindex(feats))
    end
    if single
        isempty(kept) && return only(zones)
        return _raster_zonal_rebuild(A,only(zones),otherdims)
    end
    if spatialslices === false
        return _raster_narrow(zones)
    end
    result = Array{Any}(undef,(outshape...,length(feats)))
    for j in eachindex(zones)
        if isempty(kept)
            result[j] = zones[j]
        else
            copyto!(selectdim(result,ndims(result),j),zones[j])
        end
    end
    return _raster_zonal_rebuild(A,_raster_narrow(result),(otherdims...,DD.Dim{:Zone}(Base.OneTo(length(feats)))))
end
function _raster_narrow(data::AbstractArray)
    isempty(data) && return data
    T = mapreduce(typeof,promote_type,data)
    return T === Any ? data : map(x -> convert(T,x),data)
end

_raster_zonal_rebuild(A,data,dims) = _raster_rebuild(A,data,dims;missingval=Missing <: eltype(data) ? missing : nothing)
