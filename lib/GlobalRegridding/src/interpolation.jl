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
        "no chart stencil can be written on it."))
end

chartaxes(space::RegridSpace) = _chart_required(:chartaxes, space)
chartcoords(space::RegridSpace, _) = _chart_required(:chartcoords, space)
chartlocalindex(space::RegridSpace, ::Int, ::Int) = _chart_required(:chartlocalindex, space)
chartspacing(space::RegridSpace) = _chart_required(:chartspacing, space)

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

# Chart support radius

"""
    chartradius(src_space) -> Float64

Return the larger chart-axis spacing, in radians, as a safe stencil bound for
chunk discovery.

The source must answer `true` to [`hascellchart`](@ref); a source with no chart
has no chart spacing to bound.
"""
function chartradius(src_space::RegridSpace)
    hascellchart(src_space) || throw(ArgumentError(
        "a chart stencil interpolates on the source chart, but " *
        "hascellchart(::$(typeof(src_space))) is false; use Conservative() or " *
        "NearestCell() on a source with no chart."))
    dx, dy = chartspacing(src_space)
    return Float64(max(dx, dy))
end
