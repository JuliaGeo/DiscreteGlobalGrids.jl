# Point-sampling methods and their source-chart interface.

# Chunk-local index lookup

# Use arithmetic for ranges and a dictionary for other index collections.
struct _RangeLocalIndex{R<:AbstractUnitRange{<:Integer}}
    inds::R
end

struct _DictLocalIndex
    map::Dict{Int,Int}
end

_localindexer(inds::AbstractUnitRange{<:Integer}) = _RangeLocalIndex(inds)
_localindexer(inds) = _DictLocalIndex(Dict{Int,Int}(Int(p) => k for (k, p) in enumerate(inds)))

"""
    _localindex(indexer, position::Integer) -> Int

Return the chunk-local index of `position`, or `0` when absent.
"""
_localindex(l::_RangeLocalIndex, position::Integer) =
    position in l.inds ? Int(position - first(l.inds)) + 1 : 0
_localindex(l::_DictLocalIndex, position::Integer) = get(l.map, Int(position), 0)

# Nearest-cell weights

"""
    build_weights!(coo, ::NearestCell, dst_space, dst_inds, src_space, src_inds)

Add weight 1 for the source cell containing each destination centroid. Emit no
entry when the point is outside coverage or the source belongs to another
chunk. Point samples have no coverage denominator.
"""
function build_weights!(coo::WeightCOO, method::NearestCell,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    _require_centroids(method, dst_space)
    _require_pointlocation(method, src_space)
    indexer = _localindexer(src_inds)
    for (j, i) in enumerate(dst_inds)
        source = cellat(src_space, cellcentroid(dst_space, Int(i)))
        source === nothing && continue
        k = _localindex(indexer, source)
        k == 0 && continue
        addweight!(coo, j, k, 1.0)
    end
    return coo
end

# Cell-chart interface

"""
    chartaxes(space::RegridSpace) -> (xs, ys)

Return strictly monotonic cell-centre coordinates for each separable lattice
axis. Required when [`hascellchart`](@ref) is `true`.
"""
function chartaxes end

"""
    chartcoords(space::RegridSpace, p) -> Union{Tuple{Real,Real},Nothing}

Convert `p` to native chart coordinates, or return `nothing` outside the chart.
Coordinates must use the same branch as [`chartaxes`](@ref), except on periodic
axes.
"""
function chartcoords end

"""
    chartposition(space::RegridSpace, ix::Int, iy::Int) -> Int

Return the cell position at lattice index `(ix, iy)`. Required when
[`hascellchart`](@ref) is `true`.
"""
function chartposition end

"""
    chartperiod(space::RegridSpace) -> (px, py)

Return each axis period in native coordinates, or `nothing` for no wrap.
Defaults to `(nothing, nothing)`.
"""
chartperiod(::RegridSpace) = (nothing, nothing)

"""
    chartspacing(space::RegridSpace) -> (Δx, Δy)

Return upper bounds, in radians, on adjacent-centre distance along each axis.
Required when [`hascellchart`](@ref) is `true`.
"""
function chartspacing end

function _chart_required(f::Symbol, method, space::RegridSpace)
    throw(ArgumentError(
        "$(nameof(typeof(method))) needs a cell chart, but " *
        "$(f)(::$(typeof(space))) is not defined. A space that answers " *
        "hascellchart = true must implement chartaxes, chartcoords, " *
        "chartposition and chartspacing (chartperiod defaults to no wrap)."))
end

chartaxes(space::RegridSpace) = _chart_required(:chartaxes, BilinearPoint(), space)
chartcoords(space::RegridSpace, _) = _chart_required(:chartcoords, BilinearPoint(), space)
chartposition(space::RegridSpace, ::Int, ::Int) =
    _chart_required(:chartposition, BilinearPoint(), space)
chartspacing(space::RegridSpace) = _chart_required(:chartspacing, BilinearPoint(), space)

function _require_chart(method, src_space::RegridSpace)
    hascellchart(src_space) || throw(ArgumentError(
        "$(nameof(typeof(method))) interpolates on the source chart, but " *
        "hascellchart(::$(typeof(src_space))) is false. Use Conservative() or " *
        "NearestCell() for a source with no chart."))
    return nothing
end

function _require_centroids(method, dst_space::RegridSpace)
    hasmethod(cellcentroid, Tuple{typeof(dst_space),Int}) || throw(ArgumentError(
        "$(nameof(typeof(method))) samples at destination centroids, but " *
        "cellcentroid(::$(typeof(dst_space)), ::Int) is not defined."))
    return nothing
end

function _require_pointlocation(method, src_space::RegridSpace)
    hasmethod(cellat, Tuple{typeof(src_space),USPoint}) || throw(ArgumentError(
        "$(nameof(typeof(method))) locates points in the source, but " *
        "cellat(::$(typeof(src_space)), ::UnitSphericalPoint) is not defined."))
    return nothing
end

# Coordinate location

"""
    _ChartAxis(values, period)

Prepare a monotonic chart axis for repeated location queries. Input may be
ascending or descending; `period` is `nothing` for a non-periodic axis.
"""
struct _ChartAxis
    values::Vector{Float64}
    n::Int
    reversed::Bool
    period::Union{Nothing,Float64}
end

function _ChartAxis(values, period)
    n = length(values)
    n >= 1 || throw(ArgumentError("a chart axis needs at least one cell centre"))
    v = Vector{Float64}(undef, n)
    @inbounds for (k, x) in enumerate(values)
        v[k] = Float64(x)
    end
    reversed = n > 1 && v[n] < v[1]
    reversed && reverse!(v)
    @inbounds for k in 2:n
        v[k] > v[k-1] || throw(ArgumentError(
            "chart axes must be strictly monotonic; got $(v[k-1]) then $(v[k])"))
    end
    p = period === nothing ? nothing : Float64(period)
    if p !== nothing
        p > v[n] - v[1] || throw(ArgumentError(
            "a chart axis period ($p) must exceed the span of its cell centres " *
            "($(v[n] - v[1]))"))
    end
    return _ChartAxis(v, n, reversed, p)
end

_latticeindex(ax::_ChartAxis, k::Int) = ax.reversed ? ax.n + 1 - k : k

"""
    _locate(ax::_ChartAxis, x) -> (i0, i1, w0, w1)

Return the two indices bracketing `x` and linear weights summing to one.
Non-periodic axes clamp to the nearest centre; periodic axes wrap across the seam.
"""
function _locate(ax::_ChartAxis, x::Float64)
    v = ax.values
    n = ax.n
    n == 1 && return (_latticeindex(ax, 1), _latticeindex(ax, 1), 1.0, 0.0)
    if ax.period === nothing
        x <= v[1] && return (_latticeindex(ax, 1), _latticeindex(ax, 1), 1.0, 0.0)
        x >= v[n] && return (_latticeindex(ax, n), _latticeindex(ax, n), 1.0, 0.0)
        k = searchsortedlast(v, x)
        t = (x - v[k]) / (v[k+1] - v[k])
        return (_latticeindex(ax, k), _latticeindex(ax, k + 1), 1.0 - t, t)
    end
    p = ax.period::Float64
    xw = v[1] + mod(x - v[1], p)
    k = searchsortedlast(v, xw)
    k = max(k, 1)
    if k < n
        t = (xw - v[k]) / (v[k+1] - v[k])
        return (_latticeindex(ax, k), _latticeindex(ax, k + 1), 1.0 - t, t)
    end
    seam = v[1] + p - v[n]
    t = (xw - v[n]) / seam
    return (_latticeindex(ax, n), _latticeindex(ax, 1), 1.0 - t, t)
end

# Bilinear weights

"""
    build_weights!(coo, ::BilinearPoint, dst_space, dst_inds, src_space, src_inds)

Build bilinear weights at destination centroids. Edges clamp instead of
extrapolating, while periodic axes wrap. Emit only stencil points in `src_inds`;
other chunks emit their own shares. Point samples have no denominator.
"""
function build_weights!(coo::WeightCOO, method::BilinearPoint,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    _require_chart(method, src_space)
    _require_centroids(method, dst_space)
    xs, ys = chartaxes(src_space)
    px, py = chartperiod(src_space)
    xax = _ChartAxis(xs, px)
    yax = _ChartAxis(ys, py)
    indexer = _localindexer(src_inds)
    for (j, i) in enumerate(dst_inds)
        coords = chartcoords(src_space, cellcentroid(dst_space, Int(i)))
        coords === nothing && continue
        x, y = Float64(coords[1]), Float64(coords[2])
        (isfinite(x) && isfinite(y)) || continue
        ix0, ix1, wx0, wx1 = _locate(xax, x)
        iy0, iy1, wy0, wy1 = _locate(yax, y)
        for (ix, wx) in ((ix0, wx0), (ix1, wx1)),
            (iy, wy) in ((iy0, wy0), (iy1, wy1))

            w = wx * wy
            iszero(w) && continue
            k = _localindex(indexer, chartposition(src_space, ix, iy))
            k == 0 && continue
            addweight!(coo, j, k, w)
        end
    end
    return coo
end

"""
    support_radius(::BilinearPoint, src_space) -> Float64

Return the larger chart-axis spacing, in radians, as a safe stencil bound for
chunk discovery.
"""
function support_radius(method::BilinearPoint, src_space::RegridSpace)
    _require_chart(method, src_space)
    dx, dy = chartspacing(src_space)
    return Float64(max(dx, dy))
end
