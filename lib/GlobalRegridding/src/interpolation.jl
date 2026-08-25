# Point-sampling methods, their samplers, and the source-chart interface.

# Nearest-cell weights

"""
    buildweights!(coo, ::NearestCell, dst_space, dst_inds, src_space, src_inds)

Add weight 1 for the source cell containing each destination centroid. Emit no
entry when the point is outside coverage or the source belongs to another
chunk. Point samples have no coverage denominator.
"""
function buildweights!(coo::WeightCOO, ::NearestCell,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    indexer = indexmap(src_inds)
    for (j, i) in enumerate(dst_inds)
        source = cellat(src_space, cellcentroid(dst_space, Int(i)))
        source === nothing && continue
        k = localindex(indexer, source)
        k == 0 && continue
        addweight!(coo, j, k, 1.0)
    end
    return coo
end

"""
    supportradius(::NearestCell, src_space) -> Float64

Zero. The stencil is the source cell the destination point already lies in, so
it reaches no further than that point.

A destination sample site lies inside its own cell, and so inside the covering
cap of the chunk that owns it; the source cell [`cellat`](@ref) names contains
that same point, so the source chunk owning that cell covers the point too. Two
covers sharing a point overlap, so cap overlap at radius zero already relates
every source chunk a destination tile reads. That is the containment every
method's chunk discovery already rests on, not a second assumption.
"""
supportradius(::NearestCell, ::RegridSpace) = 0.0

"""
    NearestSampler(space)

A source space prepared to answer [`NearestCell`](@ref)'s point queries.

It holds the space and nothing else: the stencil at a point is the one cell
containing it, which [`cellat`](@ref) answers by itself, so there is no state
to prepare and nothing for concurrent queries to share.
"""
struct NearestSampler{S<:RegridSpace}
    space::S
end

Base.show(io::IO, s::NearestSampler) = print(io, "NearestSampler(", s.space, ")")

"""
    sampler(::NearestCell, space::RegridSpace) -> NearestSampler

Prepare `space` to answer `NearestCell`'s point queries. This runs once per plan
or source space, never once per destination.
"""
sampler(::NearestCell, space::RegridSpace) = NearestSampler(space)

"""
    weightsat!(row, s::NearestSampler, p) -> WeightStatus

Give weight one to the source cell containing `p`.

The row is cleared on entry and left holding that single entry, named by the
source space's local index. A point the source covers nowhere leaves the row
empty and answers `WeightsOutside`; the missing policy then decides what the
destination becomes. Nothing here reads source values and nothing here is told
about chunks.
"""
function weightsat!(row::WeightRow, s::NearestSampler, p)
    empty!(row)
    i = cellat(s.space, p)
    i === nothing && return WeightsOutside
    _addentry!(row, i, 1.0)
    return WeightsMapped
end

# Cell-chart fallbacks

chartperiod(::RegridSpace) = (nothing, nothing)

function _chart_required(f::Symbol, space::RegridSpace)
    throw(ArgumentError(
        "$(typeof(space)) claims a cell chart but supplies no $(f), so " *
        "BilinearPoint cannot write a stencil on it."))
end

chartaxes(space::RegridSpace) = _chart_required(:chartaxes, space)
chartcoords(space::RegridSpace, _) = _chart_required(:chartcoords, space)
chartlocalindex(space::RegridSpace, ::Int, ::Int) = _chart_required(:chartlocalindex, space)
chartspacing(space::RegridSpace) = _chart_required(:chartspacing, space)

function _require_chart(method, src_space::RegridSpace)
    hascellchart(src_space) || throw(ArgumentError(
        "$(nameof(typeof(method))) interpolates on the source chart, but " *
        "hascellchart(::$(typeof(src_space))) is false; use Conservative() or " *
        "NearestCell() on a source with no chart."))
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
    buildweights!(coo, ::BilinearPoint, dst_space, dst_inds, src_space, src_inds)

Build bilinear weights at destination centroids. Edges clamp instead of
extrapolating, while periodic axes wrap. Emit only stencil points in `src_inds`;
other chunks emit their own shares. Point samples have no denominator.
"""
function buildweights!(coo::WeightCOO, method::BilinearPoint,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    _require_chart(method, src_space)
    xs, ys = chartaxes(src_space)
    px, py = chartperiod(src_space)
    xax = _ChartAxis(xs, px)
    yax = _ChartAxis(ys, py)
    indexer = indexmap(src_inds)
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
            k = localindex(indexer, chartlocalindex(src_space, ix, iy))
            k == 0 && continue
            addweight!(coo, j, k, w)
        end
    end
    return coo
end

"""
    supportradius(::BilinearPoint, src_space) -> Float64

Return the larger chart-axis spacing, in radians, as a safe stencil bound for
chunk discovery.
"""
function supportradius(method::BilinearPoint, src_space::RegridSpace)
    _require_chart(method, src_space)
    dx, dy = chartspacing(src_space)
    return Float64(max(dx, dy))
end
