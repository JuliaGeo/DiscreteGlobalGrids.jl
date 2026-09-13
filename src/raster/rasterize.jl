# --- reducers ---------------------------------------------------------------

# Reducers with a binary form stream one accumulator per cell; the rest gather
# each cell's values in input order and reduce them at once.
_take_first(a, b) = a
_take_last(a, b) = b
_count_fill(x) = x + 1
_reduce_op(::typeof(sum)) = Base.add_sum
_reduce_op(::typeof(prod)) = Base.mul_prod
_reduce_op(::typeof(minimum)) = min
_reduce_op(::typeof(maximum)) = max
_reduce_op(::typeof(first)) = _take_first
_reduce_op(::typeof(last)) = _take_last
_reduce_op(::typeof(any)) = |
_reduce_op(::typeof(all)) = &
_reduce_op(::Nothing) = _take_last
_reduce_op(f) = nothing
const _RASTER_THREADSAFE_OPS = (Base.add_sum, Base.mul_prod, +, *, min, max, |, &, _take_first, _take_last)
const _RASTER_THREADSAFE_REDUCERS = (sum, prod, minimum, maximum, extrema, first, last, Statistics.mean)

struct _FunctionFill{F}
    f::F
end
struct _OpFill{O,V}
    op::O
    values::V
end
struct _GatherFill{R,V}
    reducer::R
    values::V
end

# Each cell folds its incident features into its current value; a value equal
# to `mv` counts as empty and is seeded from `init` or the first feature.
function (b::_FunctionFill)(ids, old, init, mv)
    state = _raster_copy(_raster_ismissing(old, mv) ? init : old)
    for _ in ids
        state = b.f(state)
    end
    return state
end
function (b::_OpFill)(ids, old, init, mv)
    seeded = !_raster_ismissing(old, mv)
    state = _raster_copy(seeded ? old : init isa _RasterUnset ? b.values[first(ids)] : init)
    for k in (seeded || !(init isa _RasterUnset) ? 1 : 2):length(ids)
        state = b.op(state, _raster_copy(b.values[ids[k]]))
    end
    return state
end
function (b::_GatherFill)(ids, old, init, mv)
    xs = [_raster_copy(b.values[i]) for i in ids]
    return b.reducer(init isa _RasterUnset ? xs : vcat([_raster_copy(init)], xs))
end

_raster_function_init(::Type{T}) where {T} = (S = Base.nonmissingtype(T); S === Union{} ? 0 : zero(S))
function _raster_burn(data, geoms, fill, reducer, op, init, ::Type{T}) where {T}
    fill isa Function && return _FunctionFill(fill), (init isa _RasterUnset ? _raster_function_init(T) : init)
    values = _raster_fills(data, fill, length(geoms))
    op === nothing && (op = _reduce_op(reducer))
    op === nothing || return _OpFill(op, values), init
    return _GatherFill(reducer, values), init
end
_raster_threadsafe(b::_OpFill, threadsafe) = threadsafe || b.op in _RASTER_THREADSAFE_OPS
_raster_threadsafe(b::_GatherFill, threadsafe) = threadsafe || b.reducer in _RASTER_THREADSAFE_REDUCERS
_raster_threadsafe(b::_FunctionFill, threadsafe) = threadsafe || b.f === _count_fill

# Output element types follow the fold, promoted with the missing value.
function _raster_eltype(b::_FunctionFill, init, mv)
    S = Base.promote_op(b.f, typeof(init))
    return promote_type(typeof(mv), S, Base.promote_op(b.f, S))
end
function _raster_eltype(b::_OpFill, init, mv)
    V = eltype(b.values)
    S = init isa _RasterUnset ? V : Base.promote_op(b.op, typeof(init), V)
    return promote_type(typeof(mv), S, Base.promote_op(b.op, S, V))
end
function _raster_eltype(b::_GatherFill, init, mv)
    V = init isa _RasterUnset ? eltype(b.values) : promote_type(eltype(b.values), typeof(init))
    return promote_type(typeof(mv), Base.promote_op(b.reducer, Vector{V}))
end
# Uninferrable reducers pay a first pass over the touched cells for an exact type.
function _raster_alloc_eltype(b, inc, init, mv)
    T = _raster_eltype(b, init, mv)
    T === Any || return T
    touched = filter(j -> _raster_touched(inc, j), 1:length(inc.offsets) - 1)
    return mapreduce(j -> typeof(b(_raster_hits(inc, j), mv, init, mv)), promote_type, touched; init=typeof(mv))
end

# The one write kernel: every destination element folds its cell's features.
function _raster_burn!(A, axis, inc, burn, init, mv, threaded)
    _raster_foreach(CartesianIndices(A), threaded && !(parent(A) isa BitArray)) do I
        ids = _raster_hits(inc, I[axis])
        isempty(ids) || (A[I] = burn(ids, A[I], init, mv))
    end
    return A
end

function _raster_check_fill(fill, reducer, op, n)
    fill isa _RasterUnset && throw(ArgumentError("fill is required except for count"))
    n > 1 && reducer === nothing && op === nothing && !(fill isa Function) &&
        throw(ArgumentError("multiple features require a reducer, op, or function fill"))
end

# --- rasterize --------------------------------------------------------------

"""
    rasterize([reducer,] data; to, fill, boundary=:center, kw...)

Burn spherical geometries into a Cells array.

- `fill`: a scalar, one value per feature, a property Symbol, a tuple of
  Symbols or NamedTuple for several layers, or a function updating each cell.
- Multiple features need a `reducer`, binary `op`, or function fill; values
  combine in input order. `count` needs no fill; `mean` is `sum ./ count`,
  with `init` added to the sum.
- `to`: grid, cell axis, dimensional array, or stack; the spatial result is
  broadcast over the other dimensions.
- `init`, `eltype`, `missingval` set the cell state; mutable values are copied per cell.
- `threaded` parallelises cell writes; a custom `op`, reducer, or fill function
  also needs `threadsafe=true`.
"""
rasterize(reducer::Function, data; kw...) = rasterize(data; reducer, kw...)
function rasterize(data; to, fill=_RASTER_UNSET, reducer=nothing, op=nothing,
        init=_RASTER_UNSET, eltype=nothing, missingval=_RASTER_UNSET, level=nothing,
        threadsafe=false, name=nothing, metadata=nothing, kw...)
    (; boundary, shape, geometrycolumn, threaded) = _raster_options(; kw...)
    geoms = _raster_geometries(data; geometrycolumn)
    fill isa Tuple && all(x -> x isa Symbol, fill) && (fill = NamedTuple{fill}(fill))
    layer(key, target, selections) = _rasterize_layer(data, geoms, selections, target...;
        fill=_raster_layerkw(fill, key), reducer, op, init=_raster_layerkw(init, key),
        outtype=_raster_layerkw(eltype, key), missingval=_raster_layerkw(missingval, key),
        threaded, threadsafe, name=key === nothing ? name : key, metadata)
    if to isa DD.AbstractDimStack
        level === nothing || throw(ArgumentError("level is only valid with a system destination"))
        layers = _raster_eachlayer((key, A, sel) -> layer(key, _raster_target(A), sel), to, geoms;
            boundary, shape, keys=fill isa NamedTuple ? keys(fill) : keys(to))
        return _raster_stack(to, layers; metadata)
    end
    target = _raster_target(to; level)
    selections = _raster_selections(target[1], geoms; boundary, shape)
    fill isa NamedTuple || return layer(nothing, target, selections)
    return _raster_stack(to, NamedTuple{keys(fill)}(map(key -> layer(key, target, selections), keys(fill))); metadata)
end

function _rasterize_layer(data, geoms, selections, grid, dims, template; fill, reducer, op, init,
        outtype, missingval, threaded, threadsafe, name, metadata)
    if reducer === Statistics.mean && op === nothing
        return _rasterize_mean(data, geoms, selections, grid, dims, template; fill, init, outtype, missingval, threaded, threadsafe, name, metadata)
    elseif reducer === count
        missingval isa Union{_RasterUnset,Nothing} && (missingval = 0)
        fill, init, reducer = _count_fill, 0, nothing
    end
    _raster_check_fill(fill, reducer, op, length(geoms))
    mv = _raster_missing(missingval, template)
    burn, init = _raster_burn(data, geoms, fill, reducer, op, init, typeof(mv))
    inc = _RasterIncidence(ncells(grid), selections)
    T = outtype === nothing ? _raster_alloc_eltype(burn, inc, init, mv) : promote_type(outtype, typeof(mv))
    dataout = Array{T}(undef, map(length, dims))
    _raster_foreach(eachindex(dataout), threaded) do k
        dataout[k] = _raster_copy(mv)
    end
    _raster_burn!(dataout, _raster_dimnum(dims), inc, burn, init, mv, threaded && _raster_threadsafe(burn, threadsafe))
    return _raster_rebuild(template, dataout, dims; name, metadata, missingval=mv)
end

function _rasterize_mean(data, geoms, selections, grid, dims, template; fill, init, outtype, missingval, threaded, threadsafe, name, metadata)
    sums = _rasterize_layer(data, geoms, selections, grid, dims, nothing; fill, reducer=sum, op=nothing, init, outtype=nothing, missingval=missing, threaded, threadsafe, name, metadata)
    counts = _rasterize_layer(data, geoms, selections, grid, dims, nothing; fill=_RASTER_UNSET, reducer=count, op=nothing, init=_RASTER_UNSET, outtype=nothing, missingval=0, threaded, threadsafe, name, metadata)
    mv = _raster_missing(missingval, template)
    S = Base.promote_op(/, Base.nonmissingtype(Base.eltype(sums)), Int)
    means = Array{promote_type(typeof(mv), something(outtype, S))}(undef, size(sums))
    map!((s, c) -> c == 0 ? mv : s / c, means, parent(sums), parent(counts))
    return _raster_rebuild(template, means, dims; name, metadata, missingval=mv)
end

"""
    rasterize!([reducer,] destination, data; fill, kw...)

Update selected cells in place; untouched cells keep their values.

- Binary operations, streaming reducers, and function fills fold onto an
  existing nonmissing value.
- Gathering reducers such as `mean` replace the selected values.
"""
rasterize!(reducer::Function, A, data; kw...) = rasterize!(A, data; reducer, kw...)
function rasterize!(A::DD.AbstractDimArray, data; fill=_RASTER_UNSET, reducer=nothing, op=nothing,
        init=_RASTER_UNSET, missingval=_raster_missingval(A), threadsafe=false, kw...)
    (; boundary, shape, geometrycolumn, threaded) = _raster_options(; kw...)
    geoms = _raster_geometries(data; geometrycolumn)
    selections = _raster_selections(_raster_grid(A), geoms; boundary, shape)
    return _rasterize_into!(A, data, geoms, selections; fill, reducer, op, init, missingval, threaded, threadsafe)
end
function rasterize!(st::DD.AbstractDimStack, data; fill=_RASTER_UNSET, reducer=nothing, op=nothing,
        init=_RASTER_UNSET, missingval=_RASTER_UNSET, threadsafe=false, kw...)
    (; boundary, shape, geometrycolumn, threaded) = _raster_options(; kw...)
    geoms = _raster_geometries(data; geometrycolumn)
    fill isa Tuple && all(x -> x isa Symbol, fill) && (fill = NamedTuple{fill}(fill))
    _raster_eachlayer(st, geoms; boundary, shape) do key, A, selections
        _rasterize_into!(A, data, geoms, selections; fill=_raster_layerkw(fill, key), reducer, op,
            init=_raster_layerkw(init, key), missingval=_raster_layermissing(A, missingval, key), threaded, threadsafe)
    end
    return st
end
function _rasterize_into!(A, data, geoms, selections; fill, reducer, op, init, missingval, threaded, threadsafe)
    reducer === count && ((fill, init, reducer) = (_count_fill, 0, nothing))
    _raster_check_fill(fill, reducer, op, length(geoms))
    burn, init = _raster_burn(data, geoms, fill, reducer, op, init, Base.eltype(A))
    axis = _raster_dimnum(A)
    inc = _RasterIncidence(size(A, axis), selections)
    return _raster_burn!(A, axis, inc, burn, init, missingval, threaded && _raster_threadsafe(burn, threadsafe))
end

# --- masks ------------------------------------------------------------------

# Coverage of `A`'s cell axis by `data`, shaped to broadcast over the whole array.
function _raster_covered(A, data; invert=false, kw...)
    (; boundary, shape, geometrycolumn) = _raster_options(; kw...)
    selections = _raster_selections(_raster_grid(A), _raster_geometries(data; geometrycolumn); boundary, shape)
    return _raster_covered(A, selections, invert)
end
function _raster_mask(to, level, ::Type{T}, name, missingval) where {T}
    _, dims, template = _raster_target(to; level)
    return _raster_rebuild(template, Array{T}(undef, map(length, dims)), dims; name, missingval)
end

"""`boolmask(data; to, invert=false, kw...)` marks cells covered by any geometry."""
boolmask(data; to, level=nothing, kw...) = boolmask!(_raster_mask(to, level, Bool, :boolmask, nothing), data; kw...)
boolmask!(A::DD.AbstractDimArray, data; kw...) = (parent(A) .= _raster_covered(A, data; kw...); A)

"""`missingmask(data; to, kw...)` is `true` on covered cells and `missing` elsewhere."""
missingmask(data; to, level=nothing, kw...) = missingmask!(_raster_mask(to, level, Union{Missing,Bool}, :missingmask, missing), data; kw...)
missingmask!(A::DD.AbstractDimArray, data; kw...) = (parent(A) .= ifelse.(_raster_covered(A, data; kw...), true, missing); A)

"""`mask(A; with, missingval=missing, invert=false, kw...)` replaces uncovered cells."""
function mask(A::DD.AbstractDimArray; with, missingval=_raster_missingval(A), kw...)
    T = promote_type(Base.eltype(A), typeof(missingval))
    return mask!(_raster_rebuild(A, Array{T}(parent(A)), DD.dims(A); missingval); with, missingval, kw...)
end
mask!(A::DD.AbstractDimArray; with, missingval=_raster_missingval(A), kw...) =
    _raster_mask!(A, _raster_covered(A, with; kw...), missingval)
_raster_mask!(A, keep, missingval) = (parent(A) .= _raster_keep.(keep, parent(A), Ref(missingval)); A)
_raster_keep(keep, x, missingval) = keep ? x : _raster_copy(missingval)
function mask(st::DD.AbstractDimStack; with, missingval=_RASTER_UNSET, kw...)
    copied = map(keys(st)) do key
        A = st[key]
        mv = _raster_layermissing(A, missingval, key)
        _raster_rebuild(A, Array{promote_type(Base.eltype(A), typeof(mv))}(parent(A)), DD.dims(A); missingval=mv)
    end
    return mask!(_raster_stack(st, NamedTuple{keys(st)}(copied)); with, missingval, kw...)
end
function mask!(st::DD.AbstractDimStack; with, missingval=_RASTER_UNSET, invert=false, kw...)
    (; boundary, shape, geometrycolumn) = _raster_options(; kw...)
    _raster_eachlayer(st, _raster_geometries(with; geometrycolumn); boundary, shape) do key, A, selections
        _raster_mask!(A, _raster_covered(A, selections, invert), _raster_layermissing(A, missingval, key))
    end
    return st
end
