# Cell size at a level, and the level that matches a target resolution. Both
# read `cell_area` through the grid interface and choose levels through
# `levels`/`levelgrid`, so every system answers them without a method of its
# own. Included after `regridding.jl` because a target resolution may be spelled
# as a raster or a `GlobalRegridding.RegridSpace`.

# WGS84 authalic radius: the sphere on which unit-sphere areas are ellipsoidal
# areas.
const _EARTH_RADIUS = authalic_sphere(GOCore.Geodesic()).radius

"""
    cellsize(grid::AbstractGrid; over, radius, samples) -> Float64
    cellsize(sys::AbstractHierarchicalGridSystem, l::Integer; over, radius, samples) -> Float64
    cellsize(x; radius, samples) -> Float64

The typical width of a cell, in **metres**: the side of a square whose area is
the median cell area of the grid, laid on a sphere of `radius`.

`x` may be a raster (a `DimensionalData.AbstractDimArray`, measured through its
X/Y lookups) or a `GlobalRegridding.RegridSpace`, which is how a source
resolution is put on the same footing as a system's levels.

The number is a median over a sample, not a bound: cell area varies within a
level on every system that is not equal-area, and with latitude on a system
whose cells are lon/lat boxes. Use [`cell_area`](@ref) for one exact cell.

# Keywords

  - `over = nothing`: measure only the cells that meet an area of interest — an
    `Extents.Extent` in lon/lat degrees, a GeoInterface geometry, or a
    `SphericalCap`, the targets [`query`](@ref) accepts. Where cell area varies
    with position this can differ from the global number by a factor of two.
  - `radius = $(_EARTH_RADIUS)`: the sphere the area is laid on, in metres. The
    default is the WGS84 authalic radius, so an equal-area system's answer is an
    ellipsoidal one.
  - `samples = 256`: how many cells the median is taken over.

See also [`levelfor`](@ref), which asks the question the other way round.
"""
function cellsize end

cellsize(grid::AbstractGrid; over=nothing, radius::Real=_EARTH_RADIUS,
    samples::Integer=256) =
    sqrt(_mediancellarea(grid, over, Int(samples))) * radius

cellsize(sys::AbstractHierarchicalGridSystem, l::Integer; kwargs...) =
    cellsize(levelgrid(sys, l); kwargs...)

cellsize(space::GR.RegridSpace; radius::Real=_EARTH_RADIUS, samples::Integer=256) =
    sqrt(_mediancellarea(space, Int(samples))) * radius

cellsize(A::DD.AbstractDimArray; kwargs...) = cellsize(GR.RasterGrid(A); kwargs...)

"""
    levelfor(sys::AbstractHierarchicalGridSystem, target; over, radius, samples) -> Int

The level of `sys` whose cells come closest to `target`, so that
`levelgrid(sys, levelfor(sys, target))` is the grid of `sys` nearest `target` in
resolution.

`target` is either a cell size in metres (a `Real`) or anything
[`cellsize`](@ref) measures: a raster, an [`AbstractGrid`](@ref), or a
`GlobalRegridding.RegridSpace`.

Levels are compared in ratio rather than in difference, so a target falling
between two levels takes the geometrically nearer of the two. A target coarser
than every level of `sys`, or finer than every level, takes that end of
[`levels`](@ref).

The keywords are [`cellsize`](@ref)'s and mean the same thing; `over` restricts
both sides of the comparison to an area of interest.
"""
function levelfor(sys::AbstractHierarchicalGridSystem, target;
        over=nothing, radius::Real=_EARTH_RADIUS, samples::Integer=256)
    k = Int(samples)
    goal = _targetarea(target, radius, k)
    goal > 0 || throw(ArgumentError("a target cell size must be positive, got $target"))
    best, bestscore = first(levels(sys)), Inf
    for l in levels(sys)
        area = _mediancellarea(levelgrid(sys, l), over, k)
        score = abs(log(area) - log(goal))
        score < bestscore && ((best, bestscore) = (l, score))
        # Areas shrink with depth; the first level at or below the target
        # brackets it, and nothing deeper can score better.
        area <= goal && break
    end
    return best
end

_targetarea(width::Real, radius::Real, ::Int) = (width / radius)^2
_targetarea(grid::AbstractGrid, ::Real, samples::Int) =
    _mediancellarea(grid, nothing, samples)
_targetarea(space::GR.RegridSpace, ::Real, samples::Int) =
    _mediancellarea(space, samples)
_targetarea(A::DD.AbstractDimArray, radius::Real, samples::Int) =
    _targetarea(GR.RasterGrid(A), radius, samples)

# Sample with an irrational stride to avoid aliasing regular raster columns.
function _sampleindices(n::Int, samples::Int)
    k = clamp(samples, 1, n)
    k == n && return collect(1:n)
    return [mod1(round(Int, j * n * 0.6180339887498949), n) for j in 1:k]
end

function _median!(v::Vector{Float64})
    isempty(v) && throw(ArgumentError("no cell to measure"))
    sort!(v)
    k = length(v)
    return isodd(k) ? v[(k + 1) ÷ 2] : (v[k ÷ 2] + v[k ÷ 2 + 1]) / 2
end

function _mediancellarea(space::GR.RegridSpace, samples::Int)
    n = Int(ncells(space))
    n > 0 || throw(ArgumentError("cannot measure cell areas against an empty space"))
    return _median!([GR.cellarea(space, i) for i in _sampleindices(n, samples)])
end

function _mediancellarea(grid::AbstractGrid, over, samples::Int)
    n = ncells(grid)
    n > 0 || throw(ArgumentError("an empty grid has no cell size"))
    over === nothing && return _median!(
        [cell_area(grid, cellindex(grid, i)) for i in _sampleindices(n, samples)])
    return _median!([cell_area(grid, c) for c in _aoicells(grid, over, samples)])
end

# The cells an area of interest meets, sampled rather than enumerated: query the
# coarsest level whose sample is worth taking, descend while the area of
# interest holds fewer cells than that, then take the grid's own cell under each
# sampled centroid. Cost is bounded by `samples` at every level, not by the
# number of cells the area of interest covers.
function _aoicells(grid::AbstractGrid, over, samples::Int)
    sys = system(grid)
    l = level(grid)
    (sys === nothing || l === nothing) && throw(ArgumentError(
        "`over` needs a grid that belongs to a grid system; got $(typeof(grid))"))
    k = max(samples, 1)
    lq = first(levels(sys))
    while lq < l && ncells(levelgrid(sys, lq)) < 4k
        lq += 1
    end
    probe = lq == l ? grid : levelgrid(sys, lq)
    ids = query(probe, Intersects(over))
    while lq < l && length(ids) < k
        lq += 1
        probe = lq == l ? grid : levelgrid(sys, lq)
        ids = query(probe, Intersects(over))
    end
    isempty(ids) && throw(ArgumentError(
        "no cell of the grid meets the area of interest"))
    sampled = length(ids) <= k ? ids :
              [ids[i] for i in _sampleindices(length(ids), k)]
    probe === grid && return sampled
    out = cellindextype(sys)[]
    for c in sampled
        d = cellat(grid, cell_centroid(probe, c))
        d === nothing || push!(out, d)
    end
    isempty(out) && throw(ArgumentError(
        "no cell of the grid meets the area of interest"))
    return out
end
