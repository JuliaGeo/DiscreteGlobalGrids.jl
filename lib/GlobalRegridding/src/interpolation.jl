# `build_weights!` for the point-sampling methods — `NearestCell()` and
# `BilinearPoint()` — and the chart accessor a bilinear stencil is written
# against. Owned by task T4.
#
# Both methods sample the source *at a point* (the destination centroid) rather
# than integrating over the destination cell, so neither reports a denominator:
# there is no covered area to divide by, and the block finalizes as the raw
# weighted value.

# ===========================================================================
# Chunk-local index lookup
# ===========================================================================

# Membership in `src_inds` is asked once per destination cell (four times, for a
# bilinear stencil), so it must not be a linear scan. A contiguous chunk answers
# by arithmetic; anything else pays for one dictionary per block.
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

The chunk-local index of cell `position`, or `0` when the cell is not in this
chunk — the partition invariant's test, in O(1).
"""
_localindex(l::_RangeLocalIndex, position::Integer) =
    position in l.inds ? Int(position - first(l.inds)) + 1 : 0
_localindex(l::_DictLocalIndex, position::Integer) = get(l.map, Int(position), 0)

# ===========================================================================
# NearestCell
# ===========================================================================

"""
    build_weights!(coo, ::NearestCell, dst_space, dst_inds, src_space, src_inds)

One weight-1 entry per destination cell, at the source cell containing that
destination cell's centroid.

Requires [`cellcentroid`](@ref) of `dst_space` and [`cellat`](@ref) of
`src_space`. A destination cell emits nothing when its centroid falls outside
the source's coverage, and nothing *here* when the containing source cell
belongs to another chunk — that chunk's own block emits it. So a destination
row is empty in every block exactly when the centroid is a genuine miss, which
is what the missing policy then decides.

No denominator: a point sample is a value, not a coverage.
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

# ===========================================================================
# The cell chart
# ===========================================================================

"""
    chartaxes(space::RegridSpace) -> (xs, ys)

The cell-centre coordinates of `space`'s chart along each lattice axis, in the
space's own native coordinates.

Required of any space answering `true` to [`hascellchart`](@ref).

`xs` and `ys` are strictly monotonic vectors — either direction — and their
lengths are the lattice size, so cell `(ix, iy)` has its centre at
`(xs[ix], ys[iy])`. Separable axes are the whole assumption: the chart may be
irregularly spaced along either axis, but not curvilinear.

Called once per weight block, not once per cell, so materializing the vectors
is acceptable.
"""
function chartaxes end

"""
    chartcoords(space::RegridSpace, p) -> Union{Tuple{Real,Real},Nothing}

The unit-sphere point `p` in `space`'s native chart coordinates, or `nothing`
when `p` does not lie on the chart at all.

Required of any space answering `true` to [`hascellchart`](@ref).

The result is in the same coordinates as [`chartaxes`](@ref) and on the same
branch: a chart periodic in `x` may return any representative, since the
periodic reduction is done for it, but a non-periodic chart must return the
branch its axes are written on.
"""
function chartcoords end

"""
    chartposition(space::RegridSpace, ix::Int, iy::Int) -> Int

The cell position of lattice cell `(ix, iy)`, inverting the subscripting
[`chartaxes`](@ref) implies.

Required of any space answering `true` to [`hascellchart`](@ref). Indices are
always within `(1:length(xs), 1:length(ys))`; behaviour outside is the space's
own business.
"""
function chartposition end

"""
    chartperiod(space::RegridSpace) -> (px, py)

The period of each chart axis, or `nothing` for an axis that does not wrap.
Defaults to `(nothing, nothing)`.

A period is expressed in native chart coordinates — `360.0` for a global
longitude axis in degrees — and means the lattice closes on itself: the last
cell centre and the first are neighbours, separated by `xs[1] + px - xs[end]`.
A stencil then spans the seam instead of degrading against it, and a point
anywhere on the axis is reduced into range before it is located.

Reporting a period for an axis that does not in fact close wraps the stencil
onto the far edge of the domain, which is silent and wrong; reporting none for
one that does leaves a seam of one-sided stencils.
"""
chartperiod(::RegridSpace) = (nothing, nothing)

"""
    chartspacing(space::RegridSpace) -> (Δx, Δy)

An upper bound, in **radians of arc**, on the angular distance between the
centres of two lattice-adjacent cells, along each chart axis.

Required of any space answering `true` to [`hascellchart`](@ref). This is the
only place the chart's native units meet the sphere, and its only consumer is
[`support_radius`](@ref) — an upper bound is always safe, a lower bound loses
chunk pairs.
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

# ---------------------------------------------------------------------------
# Locating a coordinate between cell centres
# ---------------------------------------------------------------------------

"""
    _ChartAxis(values, period)

One chart axis prepared for repeated location queries: cell-centre coordinates
ascending, with the flip back to lattice indices remembered.

`values` may arrive in either direction — a raster whose latitudes run north to
south is the ordinary case — and `period` is the axis period or `nothing`.
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

The two lattice indices bracketing `x` and their linear weights, `w0 + w1 == 1`.

Inside the span of cell centres this is linear interpolation. Outside it the
stencil degrades to the nearest centre (`i0 == i1`, `w1 == 0`) rather than
extrapolating — so a stencil straddling one edge of the lattice is linear along
the other axis, and one at a corner is nearest-neighbour. On a periodic axis
there is no outside: `x` is reduced into range and the seam between the last
centre and the first is an interval like any other.
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

# ===========================================================================
# BilinearPoint
# ===========================================================================

"""
    build_weights!(coo, ::BilinearPoint, dst_space, dst_inds, src_space, src_inds)

Bilinear interpolation on the source cell-centre lattice, evaluated at each
destination cell's centroid.

The source must answer `true` to [`hascellchart`](@ref) and implement the chart
accessor ([`chartaxes`](@ref), [`chartcoords`](@ref), [`chartposition`](@ref),
[`chartperiod`](@ref)); the destination must provide [`cellcentroid`](@ref).

Each destination centroid is placed fractionally between the four surrounding
cell centres and given the product weights, which sum to `1`. At a lattice edge
the stencil degrades rather than extrapolates — linear along the axis that
still brackets the point, nearest-neighbour at a corner — and across a periodic
axis it spans the seam instead. A centroid the chart cannot place at all emits
nothing.

A stencil may straddle source chunks. Only the points falling in `src_inds` are
emitted here; the neighbouring chunk's block emits the rest, and accumulating
over source chunks sums the stencil back to `1`.

No denominator: a point sample is a value, not a coverage.
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

One source cell-centre spacing, in radians: the largest angular distance a
stencil point can lie from the source cell the destination centroid falls in.

A destination centroid inside source chunk A can need cell centres from chunk B
whenever it lies in the outermost half-cell of A, so chunk discovery must dilate
source extents by this much or that pair is never built and the stencil is
silently truncated. Taken as the larger of [`chartspacing`](@ref)'s two axes,
which over-covers the half-cell the stencil actually reaches.
"""
function support_radius(method::BilinearPoint, src_space::RegridSpace)
    _require_chart(method, src_space)
    dx, dy = chartspacing(src_space)
    return Float64(max(dx, dy))
end
