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
    cellsize(x; over, radius, samples) -> Float64

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
    A grid is sampled through its own system; a raster or regridding space is
    sampled at points spread evenly by area across the area of interest, each
    located by `GlobalRegridding.cellat`, so a projected raster needs the
    inverse chart that call reads. An area of interest that meets no cell of
    `x` throws an `ArgumentError`.
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

cellsize(space::GR.RegridSpace; over=nothing, radius::Real=_EARTH_RADIUS,
    samples::Integer=256) =
    sqrt(_mediancellarea(space, over, Int(samples))) * radius

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

`radius` and `samples` follow [`cellsize`](@ref). `over` restricts sampling to
an area of interest on both sides of the comparison: the candidate levels of
`sys` and a raster, grid, or regridding-space `target` are each measured within
`over`, so a global lon/lat raster asks for a finer level over the Arctic than
over the equator. A `target` in metres has one size everywhere and ignores
`over`. An `over` that meets no cell of the target throws an `ArgumentError`.
"""
function levelfor(sys::AbstractHierarchicalGridSystem, target;
        over=nothing, radius::Real=_EARTH_RADIUS, samples::Integer=256)
    k = Int(samples)
    goal = _targetarea(target, over, radius, k)
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

_targetarea(width::Real, over, radius::Real, ::Int) = (width / radius)^2
_targetarea(grid::AbstractGrid, over, ::Real, samples::Int) =
    _mediancellarea(grid, over, samples)
_targetarea(space::GR.RegridSpace, over, ::Real, samples::Int) =
    _mediancellarea(space, over, samples)
_targetarea(A::DD.AbstractDimArray, over, radius::Real, samples::Int) =
    _targetarea(GR.RasterGrid(A), over, radius, samples)

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

# A grid with a system samples through its own hierarchy; any other grid takes
# the probe route every space has.
function _mediancellarea(space::DGGSpace, over, samples::Int)
    g = space.grid
    (over === nothing || (system(g) !== nothing && level(g) !== nothing)) &&
        return _mediancellarea(g, over, samples)
    return invoke(_mediancellarea, Tuple{GR.RegridSpace,Any,Int}, space, over, samples)
end

function _mediancellarea(space::GR.RegridSpace, over, samples::Int)
    n = Int(ncells(space))
    n > 0 || throw(ArgumentError("cannot measure cell areas against an empty space"))
    over === nothing && return _median!(
        [GR.cellarea(space, i) for i in _sampleindices(n, samples)])
    return _median!([GR.cellarea(space, i) for i in _aoicells(space, over, samples)])
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
    lq, probe, sampled = _aoiprobe(sys, l, over, samples, grid)
    lq == l && return sampled
    return _locate(c -> cellat(grid, c), probe, sampled,
        "no cell of the grid meets the area of interest")
end

# A space has no hierarchy of its own to descend, so an equal-area system
# stands in as the probe and `GR.cellat` finds the space's cell under each point.
const _PROBE = HEALPixSystem()

function _aoicells(space::GR.RegridSpace, over, samples::Int)
    _, probe, sampled = _aoiprobe(_PROBE, last(levels(_PROBE)), over, samples, nothing)
    return unique!(_locate(p -> GR.cellat(space, p), probe, sampled,
        "the target does not meet the area of interest"))
end

# `deepest` is the grid standing at level `l`, queried in place of the level
# grid there so a regional grid samples only its own cells.
function _aoiprobe(sys::AbstractHierarchicalGridSystem, l::Int, over, samples::Int,
        deepest)
    k = max(samples, 1)
    lq = first(levels(sys))
    while lq < l && ncells(levelgrid(sys, lq)) < 4k
        lq += 1
    end
    probeat(lv) = lv == l && deepest !== nothing ? deepest : levelgrid(sys, lv)
    probe = probeat(lq)
    ids = query(probe, Intersects(over))
    while lq < l && length(ids) < k
        lq += 1
        probe = probeat(lq)
        ids = query(probe, Intersects(over))
    end
    isempty(ids) && throw(ArgumentError(
        "no cell of the grid meets the area of interest"))
    sampled = length(ids) <= k ? ids :
              [ids[i] for i in _sampleindices(length(ids), k)]
    return lq, probe, sampled
end

function _locate(at, probe::AbstractGrid, sampled, message::String)
    out = [d for d in (at(cell_centroid(probe, c)) for c in sampled) if d !== nothing]
    isempty(out) && throw(ArgumentError(message))
    return out
end
