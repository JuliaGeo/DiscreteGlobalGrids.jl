# Point-sampling methods, their samplers, and the source-chart interface.

# Nearest-cell weights

"""
    buildweights!(coo, ::NearestCell, dst_space, dst_inds, src_space, src_inds)

Add weight 1 for the source cell containing each destination centroid. Emit no
entry when the point is outside coverage or the source belongs to another
chunk. Point samples have no coverage denominator.
"""
function buildweights!(coo::WeightCOO, ::NearestCell,
    dst_space, dst_inds, src_space, src_inds)
    indexer = indexmap(src_inds)
    for (j, i) in enumerate(dst_inds)
        source = GridInterface.cellat(src_space, GridInterface.cell_centroid(dst_space, Int(i)))
        source === nothing && continue
        k = localindex(indexer, Int(source))
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
cap of the chunk that owns it; the source cell `GridInterface.cellat` names contains
that same point, so the source chunk owning that cell covers the point too. Two
covers sharing a point overlap, so cap overlap at radius zero already relates
every source chunk a destination tile reads. That is the containment every
method's chunk discovery already rests on, not a second assumption.
"""
supportradius(::NearestCell, ::Any) = 0.0

"""
    NearestSampler(grid)

A source grid prepared to answer [`NearestCell`](@ref)'s point queries.

It holds the grid and nothing else: the stencil at a point is the one cell
containing it, which `GridInterface.cellat` answers by itself, so there is no
state to prepare and nothing for concurrent queries to share.
"""
struct NearestSampler{S}
    space::S
end

Base.show(io::IO, s::NearestSampler) = print(io, "NearestSampler(", s.space, ")")

"""
    sampler(::NearestCell, grid) -> NearestSampler

Prepare `grid` to answer `NearestCell`'s point queries. This runs once per plan
or source grid, never once per destination.
"""
sampler(::NearestCell, space) = NearestSampler(space)

"""
    weightsat!(row, s::NearestSampler, p) -> WeightStatus

Give weight one to the source cell containing `p`.

The row is cleared on entry and left holding that single entry, named by the
source grid's position. A point the source covers nowhere leaves the row
empty and answers `WeightsOutside`; the missing policy then decides what the
destination becomes. Nothing here reads source values and nothing here is told
about chunks.
"""
function weightsat!(row::WeightRow, s::NearestSampler, p)
    empty!(row)
    i = GridInterface.cellat(s.space, p)
    i === nothing && return WeightsOutside
    addentry!(row, Int(i), 1.0)
    return WeightsMapped
end

# Cell charts: the face axes of a one-face structured grid
#
# A grid whose cells sit on a separable two-dimensional lattice — a one-face
# `Structured` grid with face axes — carries a cell chart: its sample sites are
# the centres between consecutive axis coordinates, a point's chart coordinates
# are its face coordinates, and the cell at a lattice index is the face
# position. A transparent wrapper (`AbstractWrappedGrid`) has the chart of the
# grid it wraps, with the same positions.

_chartgrid(grid) = grid
_chartgrid(g::GridInterface.AbstractWrappedGrid) = _chartgrid(parent(g))

_chartface(grid) = (g = _chartgrid(grid);
    isstructured(g) && nfaces(g) == 1 && hasfaceaxes(g, 1))

# Whether a grid has a cell chart a point method interpolates on.
_hascellchart(grid) = _chartface(grid)

# Face axes hold vertex coordinates; the chart samples at cell centres.
_axiscentres(xs) = [(Float64(xs[k]) + Float64(xs[k+1])) / 2 for k in 1:(length(xs)-1)]
_axismaxstep(xs) = maximum(abs(Float64(xs[k+1]) - Float64(xs[k])) for k in 1:(length(xs)-1))

function _chartaxes(grid)
    xs, ys = faceaxes(_chartgrid(grid), 1)
    return (_axiscentres(xs), _axiscentres(ys))
end

function _chartcoords(grid, p)
    xy = facecoords(_chartgrid(grid), 1, p)
    xy === nothing && return nothing
    return (Float64(xy[1]), Float64(xy[2]))
end

_chartlocalindex(grid, ix::Int, iy::Int) = Int(faceposition(_chartgrid(grid), 1, ix, iy))

# A face period is the chart's (360 for longitude whatever the span); the chart
# wraps only where the axis spans it, as a regional raster does not.
function _chartperiod(grid)
    g = _chartgrid(grid)
    xs, ys = faceaxes(g, 1)
    px, py = faceperiod(g, 1)
    return (_spanperiod(xs, px), _spanperiod(ys, py))
end

_spanperiod(xs, ::Nothing) = nothing
_spanperiod(xs, p) =
    isapprox(abs(Float64(last(xs)) - Float64(first(xs))), Float64(p); rtol = 1e-9) ?
    Float64(p) : nothing

_chartspacing(grid) = _facespacing(_chartgrid(grid))

"""
    _facespacing(grid) -> (Δx, Δy)

Upper bounds, in radians, on the distance between adjacent cell centres along
each face axis of a one-face structured grid. In a geographic chart a step of
Δλ degrees bounds the distance along a parallel and Δφ the distance along a
meridian; any other chart names its own bound by extending this function for
its grid type.
"""
function _facespacing(grid)
    px, py = faceperiod(grid, 1)
    (px !== nothing && Float64(px) == 360.0 && py === nothing) || throw(ArgumentError(
        "the chart spacing of a $(typeof(grid)) needs a bound in radians on the " *
        "distance between adjacent cell centres along each face axis; only a " *
        "geographic chart (faceperiod (360, nothing)) derives one, so define " *
        "`GlobalRegridding._facespacing(::$(typeof(grid)))`"))
    xs, ys = faceaxes(grid, 1)
    return (deg2rad(_axismaxstep(xs)), deg2rad(_axismaxstep(ys)))
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

# Chart support radius

"""
    chartradius(src_grid) -> Float64

Return the larger chart-axis spacing, in radians, as a safe stencil bound for
chunk discovery.

The source must have a cell chart (a one-face `Structured` grid with face
axes); a source with no chart has no chart spacing to bound.
"""
function chartradius(src_space)
    _hascellchart(src_space) || throw(ArgumentError(
        "a chart stencil interpolates on the source chart, but a " *
        "$(typeof(src_space)) has no cell chart (a one-face Structured grid " *
        "with face axes); use Conservative() or NearestCell() on a source " *
        "with no chart."))
    dx, dy = _chartspacing(src_space)
    return Float64(max(dx, dy))
end
