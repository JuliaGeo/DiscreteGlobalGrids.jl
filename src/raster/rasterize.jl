import DimensionalData as RasterDD
import Statistics

function _raster_fillvalues(features, fill; single=false)
    fill isa Function && return fill
    if fill isa Symbol
        return [f.properties isa AbstractDict ? f.properties[fill] : Tables.getcolumn(f.properties, fill) for f in features]
    elseif single
        return [fill]
    elseif !(fill isa Union{Number,AbstractString,Missing,Nothing}) && applicable(iterate, fill)
        values = collect(fill)
        length(values) == 1 && return [only(values) for _ in features]
        length(values) == length(features) || throw(DimensionMismatch("fill must have one value per feature"))
        return values
    end
    return [fill for _ in features]
end

function _raster_validate_options(; kw...)
    isempty(kw) || throw(ArgumentError("unsupported DGGS rasterization keywords: $(keys(kw)); file output, res and size are not supported"))
end

function _raster_incidence(grid, features; boundary, shape, selections=nothing)
    hits = [Int[] for _ in 1:ncells(grid)]
    for (i, feature) in enumerate(features)
        for j in (isnothing(selections) ? _raster_indices(grid, feature.geom; boundary, shape) : selections[i])
            push!(hits[j], i)
        end
    end
    return hits
end

function _raster_reduce_values(reducer, op, values, ids, init, old, mv)
    if values isa Function
        state = _raster_ismissing(old, mv) ? (init isa _RasterUnset ? zero(Base.nonmissingtype(typeof(old))) : deepcopy(init)) : deepcopy(old)
        for _ in ids
            state = values(state)
        end
        return state
    elseif reducer === count
        return (_raster_ismissing(old, mv) ? 0 : old) + length(ids)
    elseif !isnothing(op)
        state = _raster_ismissing(old, mv) ? (init isa _RasterUnset ? deepcopy(values[first(ids)]) : deepcopy(init)) : deepcopy(old)
        start = _raster_ismissing(old, mv) && init isa _RasterUnset ? 2 : 1
        for k in start:length(ids)
            state = op(state, deepcopy(values[ids[k]]))
        end
        return state
    end
    xs = [deepcopy(values[i]) for i in ids]
    init isa _RasterUnset || pushfirst!(xs, deepcopy(init))
    return isnothing(reducer) ? only(xs) : reducer(xs)
end

"""
    rasterize([reducer,] data; to, fill, boundary=:center, ...)

Burn spherical geometries into a Cells array. `fill` is a scalar, one value per
feature, a property Symbol, a tuple of property Symbols, a NamedTuple of layers,
or a function updating each selected cell. Multiple features require a reducer,
binary `op`, or function fill. Values are combined in source order. Mutable
initial and missing values are copied independently for every cell.

`to` may be a grid, cell axis, dimensional array, or stack. Nonspatial dimensions
are retained and each spatial result is broadcast over their slices. `threaded` parallelizes independent destination writes. Arbitrary in-place user
operations use a serial fallback unless `threadsafe=true` is supplied.
"""
rasterize(reducer::Function, data; kw...) = rasterize(data; reducer, kw...)
function rasterize(data; to, fill=_RASTER_UNSET, reducer=nothing, op=nothing,
    init=_RASTER_UNSET, eltype=nothing, missingval=_RASTER_UNSET,
    boundary=:center, shape=nothing, geometrycolumn=nothing, level=nothing,
    threaded=false, threadsafe=false, progress=true, verbose=true,
    name=nothing, metadata=nothing, crs=nothing, mappedcrs=nothing,
    _features=nothing, _selections=nothing, kw...)
    _raster_validate_options(; kw...)
    boundary = _raster_boundary(boundary)
    shape = _raster_shape(shape)
    crs === nothing && mappedcrs === nothing || throw(ArgumentError("DGGS cell axes have intrinsic spherical coordinates; transform input geometries to longitude/latitude instead of setting crs/mappedcrs"))
    features = isnothing(_features) ? _raster_features(data; geometrycolumn) : _features
    fill isa Tuple && all(x -> x isa Symbol, fill) && (fill = NamedTuple{fill}(fill))
    if fill isa NamedTuple
        shared = if !(to isa RasterDD.AbstractDimStack)
            grid, _, _ = _raster_target(to; level)
            isnothing(_selections) ? [_raster_indices(grid, f.geom; boundary, shape) for f in features] : _selections
        else
            nothing
        end
        cache = Dict{Any,Vector{Vector{Int}}}()
        layers = map(keys(fill)) do key
            template = to isa RasterDD.AbstractDimStack ? to[key] : to
            selected = to isa RasterDD.AbstractDimStack ? _raster_cached_selections!(cache,template,features;boundary,shape) : shared
            rasterize(data; to=template, fill=getproperty(fill,key), reducer, op,
                init=init isa NamedTuple ? getproperty(init,key) : init,
                eltype=eltype isa NamedTuple ? getproperty(eltype,key) : eltype,
                missingval=missingval isa NamedTuple ? getproperty(missingval,key) : missingval,
                boundary, shape, geometrycolumn, level, threaded, threadsafe, progress, verbose,
                name=key, metadata, crs, mappedcrs, _features=features, _selections=selected)
        end
        return _raster_stack(to, NamedTuple{keys(fill)}(layers);metadata)
    elseif to isa RasterDD.AbstractDimStack
        cache = Dict{Any,Vector{Vector{Int}}}()
        layers = map(key -> rasterize(data; to=to[key], fill, reducer, op,
            init=init isa NamedTuple ? getproperty(init,key) : init,
            eltype=eltype isa NamedTuple ? getproperty(eltype,key) : eltype,
            missingval=missingval isa NamedTuple ? getproperty(missingval,key) : missingval,
            boundary, shape, geometrycolumn, level, threaded, threadsafe,
            progress, verbose, name=key, metadata, crs, mappedcrs, _features=features,
            _selections=_raster_cached_selections!(cache,to[key],features;boundary,shape)), keys(to))
        return _raster_stack(to, NamedTuple{keys(to)}(layers);metadata)
    end
    reducer === count && (fill = 1)
    fill isa _RasterUnset && throw(ArgumentError("fill is required except for count"))
    length(features) > 1 && isnothing(reducer) && isnothing(op) && !(fill isa Function) &&
        throw(ArgumentError("multiple features require a reducer, op, or function fill"))
    grid, dims, template = _raster_target(to; level)
    mv = missingval isa _RasterUnset ? (reducer === count ? 0 : isnothing(template) ? missing : _raster_missingval(template)) : missingval
    mv === nothing && (mv = reducer === count ? 0 : missing)
    values = _raster_fillvalues(features, fill;single=_raster_issingle(data))
    spatial = _raster_spatial(grid, features, values, reducer, op, init, mv; boundary, shape, selections=_selections)
    T = isnothing(eltype) ? (isempty(spatial) ? typeof(mv) : mapreduce(typeof, promote_type, spatial)) : promote_type(eltype, typeof(mv))
    dataout = Array{T}(undef, map(length, dims))
    axis = only(findall(d -> RasterDD.lookup(d) isa Union{AbstractCellLookup,_RasterGridLookup}, dims))
    _raster_foreach(CartesianIndices(dataout), threaded) do I
        dataout[I] = deepcopy(spatial[I[axis]])
    end
    return _raster_rebuild(template, dataout, dims; name, metadata, crs, mappedcrs, missingval=mv)
end

"""`rasterize!([reducer,] destination, data; fill, ...)` updates selected cells,
retaining untouched cells. Binary operations, associative built-in reducers, and
function fills use existing nonmissing values as their initial state. Other
iterable reducers replace selected values (including mean, whose previous count
is unavailable)."""
rasterize!(reducer::Function, A, data; kw...) = rasterize!(A, data; reducer, kw...)
function rasterize!(A::RasterDD.AbstractDimArray, data; fill=_RASTER_UNSET,
    reducer=nothing, op=nothing, init=_RASTER_UNSET, missingval=_raster_missingval(A),
    boundary=:center, shape=nothing, geometrycolumn=nothing,
    threaded=false, threadsafe=false, progress=true, verbose=true,
    _features=nothing, _selections=nothing, kw...)
    _raster_validate_options(; kw...)
    boundary = _raster_boundary(boundary)
    shape = _raster_shape(shape)
    features = isnothing(_features) ? _raster_features(data; geometrycolumn) : _features
    reducer === count && (fill = 1)
    fill isa _RasterUnset && throw(ArgumentError("fill is required except for count"))
    length(features) > 1 && isnothing(reducer) && isnothing(op) && !(fill isa Function) &&
        throw(ArgumentError("multiple features require a reducer, op, or function fill"))
    grid, _, _ = _raster_target(A)
    hits = _raster_incidence(grid, features; boundary, shape, selections=_selections)
    values = _raster_fillvalues(features, fill;single=_raster_issingle(data))
    axis = _raster_dimnum(A)
    if isnothing(op)
        op = reducer === sum ? Base.add_sum : reducer === prod ? Base.mul_prod :
            reducer === minimum ? min : reducer === maximum ? max :
            reducer === first ? ((a,b) -> a) : reducer === last ? ((a,b) -> b) : nothing
    end
    values isa Function && init isa _RasterUnset && (init = zero(Base.nonmissingtype(Base.eltype(A))))
    safe = threadsafe || reducer in (sum, prod, minimum, maximum, extrema, first, last, count, Statistics.mean) || op in (+, *, min, max)
    _raster_foreach(CartesianIndices(A), threaded && safe && !(parent(A) isa BitArray)) do I
        ids = hits[I[axis]]
        isempty(ids) && return
        A[I] = _raster_reduce_values(reducer, op, values, ids, init, A[I], missingval)
    end
    return A
end
function rasterize!(A::RasterDD.AbstractDimStack, data; fill=_RASTER_UNSET,
    init=_RASTER_UNSET, missingval=_RASTER_UNSET, geometrycolumn=nothing,
    boundary=:center,shape=nothing,kw...)
    features = _raster_features(data; geometrycolumn)
    cache = Dict{Any,Vector{Vector{Int}}}()
    fill isa Tuple && all(x -> x isa Symbol, fill) && (fill = NamedTuple{fill}(fill))
    for key in keys(A)
        rasterize!(A[key], data; fill=fill isa NamedTuple ? getproperty(fill,key) : fill,
            init=init isa NamedTuple ? getproperty(init,key) : init,
            missingval=missingval isa _RasterUnset ? _raster_missingval(A[key]) :
                missingval isa NamedTuple ? getproperty(missingval,key) : missingval,
            _features=features, _selections=_raster_cached_selections!(cache,A[key],features;boundary,shape), boundary,shape,kw...)
    end
    return A
end

"""`boolmask(data; to, invert=false, ...)` marks cells covered by any geometry."""
function boolmask(data; to, invert=false, kw...)
    A = rasterize(count, data; to, kw...)
    return _raster_rebuild(A, map(x -> invert ? iszero(x) : !iszero(x), parent(A)), RasterDD.dims(A); name=:boolmask, metadata=nothing, missingval=nothing)
end
function boolmask!(A::RasterDD.AbstractDimArray, data; invert=false, kw...)
    B = boolmask(data; to=A, invert, kw...)
    copyto!(parent(A), parent(B))
    return A
end
"""`missingmask(data; to, ...)` returns true on covered cells and missing elsewhere."""
function missingmask(data; to, invert=false, kw...)
    B = boolmask(data; to, invert, kw...)
    return _raster_rebuild(B, map(x -> x ? true : missing, parent(B)), RasterDD.dims(B); name=:missingmask, metadata=nothing, missingval=missing)
end
function missingmask!(A::RasterDD.AbstractDimArray, data; kw...)
    B = missingmask(data; to=A, kw...)
    copyto!(parent(A), parent(B))
    return A
end
"""`mask(A; with, missingval=missing, invert=false, ...)` replaces uncovered cells."""
function mask(A::RasterDD.AbstractDimArray; with, missingval=_raster_missingval(A), invert=false, kw...)
    T = promote_type(Base.eltype(A), typeof(missingval))
    B = _raster_rebuild(A, Array{T}(parent(A)), RasterDD.dims(A); name=nothing, metadata=nothing, missingval)
    return mask!(B; with, missingval, invert, kw...)
end
function mask!(A::RasterDD.AbstractDimArray; with, missingval=_raster_missingval(A), invert=false, kw...)
    B = boolmask(with; to=A, invert, kw...)
    for I in eachindex(A, B)
        B[I] || (A[I] = deepcopy(missingval))
    end
    return A
end

# Common reductions stream one accumulator per target cell. User iterable
# reducers use an ordered reverse incidence, never a dense cell-by-feature mask.
function _raster_spatial(grid, features, values, reducer, op, init, mv; boundary, shape, selections=nothing)
    combine = !isnothing(op) ? op : reducer === sum ? Base.add_sum : reducer === prod ? Base.mul_prod :
        reducer === minimum ? min : reducer === maximum ? max :
        reducer === extrema ? _raster_extrema_add : reducer === nothing ? ((a,b) -> b) :
        reducer === first ? ((a,b) -> a) : reducer === last ? ((a,b) -> b) :
        reducer === Statistics.mean ? Base.add_sum : op
    fast = !isnothing(combine) || reducer === count || values isa Function
    if !fast
        hits = _raster_incidence(grid, features; boundary, shape, selections)
        return map(hits) do ids
            isempty(ids) ? deepcopy(mv) : _raster_reduce_values(reducer, op, values, ids, init, mv, mv)
        end
    end
    firstvalue = x -> isnothing(op) && reducer in (sum, Statistics.mean) ? sum((x,)) :
        isnothing(op) && reducer === prod ? prod((x,)) :
        isnothing(op) && reducer === extrema ? (x,x) : x
    seedtype = values isa Function ? (init isa _RasterUnset ? Int : typeof(init)) :
        isempty(values) ? Base.eltype(values) : mapreduce(typeof, promote_type, values)
    valuetype = values isa Function ? Base.promote_op(values, seedtype) :
        reducer === count ? Int : Base.promote_op(firstvalue, seedtype)
    statetype = promote_type(valuetype, init isa _RasterUnset ? Union{} : typeof(init))
    if !isnothing(combine)
        statetype = promote_type(statetype, Base.promote_op(combine, statetype, seedtype))
    end
    reducer === Statistics.mean && isnothing(op) && (statetype = promote_type(statetype, Base.promote_op(/, statetype, Int)))
    T = promote_type(typeof(mv), statetype)
    state = T[deepcopy(mv) for _ in 1:ncells(grid)]
    counts = zeros(Int, ncells(grid))
    for (i, feature) in enumerate(features)
        for j in (isnothing(selections) ? _raster_indices(grid, feature.geom; boundary, shape) : selections[i])
            if reducer === count
                state[j] = counts[j] + 1
            elseif values isa Function
                old = counts[j] == 0 ? (init isa _RasterUnset ? 0 : deepcopy(init)) : state[j]
                state[j] = values(old)
            elseif counts[j] == 0
                value = deepcopy(values[i])
                state[j] = init isa _RasterUnset ? firstvalue(value) : combine(deepcopy(init), value)
            else
                state[j] = combine(state[j], deepcopy(values[i]))
            end
            counts[j] += 1
        end
    end
    if reducer === Statistics.mean && isnothing(op)
        for j in eachindex(state)
            counts[j] == 0 && continue
            state[j] /= counts[j] + (init isa _RasterUnset ? 0 : 1)
        end
    end
    return state
end

function mask(A::RasterDD.AbstractDimStack; with, missingval=_RASTER_UNSET, kw...)
    layers = map(keys(A)) do key
        mv = missingval isa _RasterUnset ? _raster_missingval(A[key]) :
            missingval isa NamedTuple ? getproperty(missingval, key) : missingval
        mask(A[key]; with, missingval=mv, kw...)
    end
    return _raster_stack(A, NamedTuple{keys(A)}(layers))
end
function mask!(A::RasterDD.AbstractDimStack; with, missingval=_RASTER_UNSET, kw...)
    for key in keys(A)
        mv = missingval isa _RasterUnset ? _raster_missingval(A[key]) :
            missingval isa NamedTuple ? getproperty(missingval, key) : missingval
        mask!(A[key]; with, missingval=mv, kw...)
    end
    return A
end

function _raster_foreach(f, indices, threaded)
    if threaded && Threads.nthreads() > 1
        Threads.@threads :dynamic for k in eachindex(indices)
            f(indices[k])
        end
    else
        foreach(f, indices)
    end
    return nothing
end

_raster_extrema_add(a::Tuple,b) = (min(a[1],b),max(a[2],b))
_raster_extrema_add(a,b) = (min(a,b),max(a,b))
