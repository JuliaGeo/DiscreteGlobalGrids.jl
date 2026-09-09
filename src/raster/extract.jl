"""
    DiscreteGlobalGrids.extract(A, geometries; boundary=:center, kw...)

Extract cell values as named-tuple rows. Points use containing-cell lookup;
lines intersect cells; polygons use `:center`, `:intersects` (`:touches`), or
`:inside`. `index=true` reports the local cell-axis index. Unsampled dimensions
remain labelled slices. `flatten=false` groups non-point results by feature.
"""
function extract(A::Union{DD.AbstractDimArray,DD.AbstractDimStack}, data;
        geometrycolumn=nothing, names=nothing, name=names, skipmissing=false,
        flatten=true, id=false, geometry=true, index=false, boundary=:center,
        shape=nothing, atol=nothing, threaded=false, progress=true, kw...)
    _raster_check_keywords(kw)
    boundary = _raster_boundary(boundary)
    shape = _raster_shape(shape)
    # Cell membership is an interval lookup: atol cannot change containment.
    feats = _raster_features(data; geometrycolumn)
    layers = _raster_extract_layers(A, name)
    firstlayer = first(values(layers))
    grid, _, _ = _raster_target(firstlayer)
    for layer in values(layers)
        lk = DD.lookup(layer, _raster_dimnum(layer))
        lk == DD.lookup(firstlayer,_raster_dimnum(firstlayer)) || throw(ArgumentError("Extraction stack layers must share the cell axis"))
    end
    reservednames = ((id ? (:id,) : ())..., (geometry ? (:geometry,) : ())..., (index ? (:index,) : ())...)
    any(k -> k in reservednames, keys(layers)) && throw(ArgumentError("Layer names conflict with requested extraction fields"))
    groups = Vector{Vector{NamedTuple}}(undef,length(feats))
    function onefeature(j)
        feat = feats[j]
        point = GI.trait(feat.geom) isa GI.AbstractPointTrait
        absent = feat.geom === nothing || ismissing(feat.geom)
        indices = _raster_indices(grid, feat.geom; boundary,shape)
        rows = NamedTuple[]
        visits = isempty(indices) && (point || absent) ? (nothing,) : indices
        for i in visits
            vals = map(values(layers)) do layer
                i === nothing ? missing : _raster_extract_value(layer,i)
            end
            if skipmissing && any(t -> _raster_extract_missing(t[1],_raster_missingval(t[2])),zip(vals,values(layers)))
                continue
            end
            coord = absent ? missing : point ? _raster_point_tuple(feat.geom) : _raster_center(grid,i)
            row = merge(id ? (;id=j) : (;), geometry ? (;geometry=coord) : (;),
                index ? (;index=i === nothing ? missing : i) : (;), NamedTuple{keys(layers)}(vals))
            push!(rows,row)
        end
        groups[j] = rows
    end
    if threaded && Threads.nthreads()>1
        Threads.@threads for j in eachindex(feats)
            onefeature(j)
        end
    else
        foreach(onefeature,eachindex(feats))
    end
    allpoints = all(f -> f.geom === nothing || ismissing(f.geom) || GI.trait(f.geom) isa GI.AbstractPointTrait, feats)
    return flatten || allpoints ? reduce(vcat,groups; init=NamedTuple[]) : groups
end

function _raster_extract_layers(A::DD.AbstractDimArray, name)
    n = name === nothing ? DD.name(A) : name
    n = n isa Tuple ? only(n) : n
    n = n isa Symbol && !isempty(string(n)) ? n : :layer
    return NamedTuple{(n,)}((A,))
end
function _raster_extract_layers(A::DD.AbstractDimStack, name)
    names = name === nothing ? keys(A) : name isa Symbol ? (name,) : Tuple(name)
    isempty(names) && throw(ArgumentError("Select at least one stack layer"))
    return NamedTuple{names}(Tuple(A[k] for k in names))
end
function _raster_extract_value(A, i)
    d = _raster_dimnum(A)
    inds = ntuple(k -> k==d ? i : Colon(), ndims(A))
    return ndims(A)==1 ? A[i] : view(A,inds...)
end
_raster_extract_missing(v::DD.AbstractDimArray,mv) = any(x -> _raster_ismissing(x,mv),v)
_raster_extract_missing(v,mv) = _raster_ismissing(v,mv)
