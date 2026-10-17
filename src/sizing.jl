# Cell size at a level, and the level that matches a target resolution. Both
# read `cell_area` through the grid interface and choose levels through
# `levels`/`levelgrid`, so every system answers them without a method of its
# own. Included after `regrid_targets.jl` because a target resolution may be
# spelled as a raster or any GridInterface grid a regrid accepts.

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
X/Y lookups) or any GridInterface grid (a `GlobalRegridding.RasterGrid`, a
chunked wrapper of a grid here), which is how a source resolution is put on the
same footing as a system's levels.

The number is a median over a sample, not a bound: cell area varies within a
level on every system that is not equal-area, and with latitude on a system
whose cells are lon/lat boxes. Use [`cell_area`](@ref) for one exact cell.

# Keywords

  - `over = nothing`: measure only the cells that meet an area of interest — an
    `Extents.Extent` in lon/lat degrees, a GeoInterface geometry, or a
    `SphericalCap`, the targets [`query`](@ref) accepts. Where cell area varies
    with position this can differ from the global number by a factor of two.
    A grid is sampled through its own system; a raster or any other grid is
    sampled at points spread evenly by area across the area of interest, each
    located by `GridInterface.cellat`, so a projected raster needs the
    inverse chart that call reads. An area of interest that meets no cell of
    `x` throws an `ArgumentError`.
  - `radius = $(_EARTH_RADIUS)`: the sphere the area is laid on, in metres. The
    default is the WGS84 authalic radius, so an equal-area system's answer is an
    ellipsoidal one.
  - `samples = 256`: how many cells the median is taken over.

See also [`levelfor`](@ref), which asks the question the other way round.
"""
cellsize

cellsize(grid::AbstractGrid; over=nothing, radius::Real=_EARTH_RADIUS,
    samples::Integer=256) =
    sqrt(_mediancellarea(grid, over, Int(samples))) * radius

cellsize(sys::AbstractHierarchicalGridSystem, l::Integer; kwargs...) =
    cellsize(levelgrid(sys, l); kwargs...)

# `GridInterface.cellsize(grid)` is a one-argument radian answer for any grid;
# the metre form with its keywords is this package's, answered for its own
# grids above and for `GlobalRegridding`'s raster grid. A wrapper (a chunked
# grid) forwards to the grid it wraps, keywords and all, by GridInterface's
# own forwarding.
cellsize(grid::GR.RasterGrid; over=nothing, radius::Real=_EARTH_RADIUS,
    samples::Integer=256) =
    sqrt(_mediancellarea(grid, over, Int(samples))) * radius

cellsize(A::DD.AbstractDimArray; kwargs...) = cellsize(GR.RasterGrid(A); kwargs...)

"""
    levelfor(sys::AbstractHierarchicalGridSystem, target; over, radius, samples) -> Int

The level of `sys` whose cells come closest to `target`, so that
`levelgrid(sys, levelfor(sys, target))` is the grid of `sys` nearest `target` in
resolution.

`target` is either a cell size in metres (a `Real`) or anything
[`cellsize`](@ref) measures: a raster, an [`AbstractGrid`](@ref), or any other
GridInterface grid.

Levels are compared in ratio rather than in difference, so a target falling
between two levels takes the geometrically nearer of the two. A target coarser
than every level of `sys`, or finer than every level, takes that end of
[`levels`](@ref).

`radius` and `samples` follow [`cellsize`](@ref). `over` restricts the
comparison to an area of interest: the candidate levels of `sys` and a raster
or grid `target` are each measured within it.

  - A global lon/lat raster asks for a finer level over the Arctic than over
    the equator.
  - A `target` in metres has one size everywhere and ignores `over`.
  - An `over` meeting no cell of the target throws an `ArgumentError`.
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
function _targetarea(grid, over, ::Real, samples::Int)
    GridInterface.isgrid(grid) || throw(ArgumentError(
        "a target is a cell size in metres, a raster, a grid, or a system and " *
        "level; got $(typeof(grid))"))
    return _mediancellarea(grid, over, samples)
end
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

# A grid of this package samples through its own hierarchy; a wrapper of one
# (a chunked grid) is measured as the grid it wraps; any other grid takes the
# probe route every grid has.
_mediancellarea(g::GridInterface.AbstractWrappedGrid, over, samples::Int) =
    _mediancellarea(parent(g), over, samples)

function _mediancellarea(grid, over, samples::Int)
    n = Int(ncells(grid))
    n > 0 || throw(ArgumentError("cannot measure cell areas against an empty grid"))
    area(i) = Float64(GridInterface.cell_area(grid, Int(i)))
    over === nothing && return _median!([area(i) for i in _sampleindices(n, samples)])
    return _median!([area(i) for i in _aoicells(grid, over, samples)])
end

function _mediancellarea(grid::AbstractGrid, over, samples::Int)
    # Without a system there is no hierarchy to sample an area of interest
    # through, so such a grid takes the probe route every grid has.
    over === nothing || (system(grid) !== nothing && level(grid) !== nothing) ||
        return invoke(_mediancellarea, Tuple{Any,Any,Int}, grid, over, samples)  # positions
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
    # No system, no hierarchy: the probe route, which answers positions.
    (sys === nothing || l === nothing) &&
        return invoke(_aoicells, Tuple{Any,Any,Int}, grid, over, samples)
    lq, probe, sampled = _aoiprobe(sys, l, over, samples, grid)
    lq == l && return sampled
    return _locate(c -> cellat(grid, c), probe, sampled,
        "no cell of the grid meets the area of interest")
end

# A subset's coarse ancestors can contain mostly absent cells. Sample its own
# storage, then fill any remaining slots with exact spatial-query matches.
function _aoicells(grid::PartialGrid, over, samples::Int)
    target = Engine._query_target(over)
    k = min(max(samples, 1), ncells(grid))
    out = cellindextype(system(grid))[]
    pred = Intersects(nothing)
    indices = unique!(_sampleindices(ncells(grid), k))
    for i in indices
        c = cellindex(grid, i)
        Engine._matches(pred, target, grid, c) && push!(out, c)
    end
    if length(out) < k && length(indices) < ncells(grid)
        seen = Set(out)
        STI.depth_first_search(Base.Fix1(Extents.intersects, target.cap), treeify(grid)) do i
            c = cellindex(grid, i)
            c in seen && return nothing
            if Engine._matches(pred, target, grid, c)
                push!(out, c)
                push!(seen, c)
                length(out) == k && return GO.LoopStateMachine.Action(:full_return, nothing)
            end
            return nothing
        end
    end
    isempty(out) && throw(ArgumentError("no cell of the grid meets the area of interest"))
    return out
end

# A grid of another package has no hierarchy of its own to descend, so an
# equal-area system stands in as the probe and `GridInterface.cellat` finds the
# grid's cell under each point.
const _PROBE = HEALPixSystem()

function _aoicells(grid, over, samples::Int)
    _, probe, sampled = _aoiprobe(_PROBE, last(levels(_PROBE)), over, samples, nothing)
    return unique!(_locate(p -> GridInterface.cellat(grid, p), probe, sampled,
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
