# Generic point interpolation between source sample sites.
#
# One destination sample site at a time: find the dual cell of source sample
# sites that contains it, then weight that cell's nodes with the coordinates its
# basis names. Nothing here knows about chunks. A row carries the source space's
# local indices and whoever asked for it partitions them.

# --- the reusable row ------------------------------------------------------

"""
    WeightRow()

The source cells one destination point takes its value from, and their weights.

`indices` and `weights` are parallel: entry `k` gives weight `weights[k]` to the
source cell the source space calls `indices[k]`. [`weightsat!`](@ref) clears the
row on entry, so one row serves a whole sweep and grows only to the widest
stencil it has met. A row belongs to one task.
"""
struct WeightRow
    indices::Vector{Int}
    weights::Vector{Float64}
end

WeightRow() = WeightRow(Int[], Float64[])

Base.length(row::WeightRow) = length(row.indices)
Base.isempty(row::WeightRow) = isempty(row.indices)

function Base.empty!(row::WeightRow)
    empty!(row.indices)
    empty!(row.weights)
    return row
end

Base.show(io::IO, row::WeightRow) = print(io, "WeightRow(", length(row), " entries)")

# A zero weight is dropped rather than stored, so a point on an edge or at a
# node emits only the nodes that carry it.
@inline function _addentry!(row::WeightRow, i::Integer, w::Float64)
    w > 0 || return row
    push!(row.indices, Int(i))
    push!(row.weights, w)
    return row
end

"""
    WeightStatus

What [`weightsat!`](@ref) made of one destination point.

`WeightsMapped` says the row holds a complete stencil. Every other value says
the row is empty and the destination takes no weights at all; they differ only
in the reason, which is diagnostic. Execution reads
[`ismapped`](@ref) and nothing else.

  - `WeightsOutside` — the point lies outside the dual complex, or outside the
    cell that was offered for it;
  - `WeightsRim` — a required sample site is not in the source collection, so no
    dual cell exists there. A space's own `weightsat!` answers this;
  - `WeightsDegenerate` — the cell is unusable: repeated nodes, no area, a fold,
    a reflex corner, or an inverse map that did not converge.
"""
@enum WeightStatus::UInt8 begin
    WeightsMapped
    WeightsOutside
    WeightsRim
    WeightsDegenerate
end

"""
    ismapped(status::WeightStatus) -> Bool

Whether `status` says the row carries a stencil.
"""
@inline ismapped(status::WeightStatus) = status === WeightsMapped

"""
    BasisKind

Which coordinates weight a dual cell's nodes: `Linear` for triangle barycentric
coordinates, `Bilinear` for inverse isoparametric Q1 on a quadrilateral, and
`MeanValue` for mean-value coordinates on a convex polygon.

The kind belongs to the cell, not to its node count: a four-node cell is
`Bilinear` only where the space that built it says so.
"""
@enum BasisKind::UInt8 begin
    Linear
    Bilinear
    MeanValue
end

# --- tolerances ------------------------------------------------------------

# A coordinate below this is outside its cell. This is a rounding margin, not a
# policy: a point further out than this is unmapped, and no point is ever
# clamped into a cell it does not lie in.
const COORD_TOL = 1e-12

# A cell is degenerate when its area, or the turn at a corner, is this small
# beside the square of its own size. Repeated nodes, collinear triangles,
# zero-area or folded quadrilaterals and reflex corners all arrive here.
const SHAPE_TOL = 1e-12

# The inverse bilinear map's Newton iteration: a step moving both coordinates
# less than this has converged, a cell needing more than these steps or leaving
# a residual above `NEWTON_RESIDUAL` times its size is treated as degenerate.
const NEWTON_TOL = 1e-12
const NEWTON_MAXITER = 20
const NEWTON_RESIDUAL = 1e-8

# --- plane arithmetic ------------------------------------------------------

@inline _sub(a::NTuple{2,Float64}, b::NTuple{2,Float64}) = (a[1] - b[1], a[2] - b[2])
@inline _cross(a::NTuple{2,Float64}, b::NTuple{2,Float64}) = a[1] * b[2] - a[2] * b[1]
@inline _dot(a::NTuple{2,Float64}, b::NTuple{2,Float64}) = a[1] * b[1] + a[2] * b[2]
@inline _norm2(a::NTuple{2,Float64}) = _dot(a, a)

@inline _next(i::Int, n::Int) = i == n ? 1 : i + 1

# A point at a node reproduces that node exactly, whatever the surrounding cell
# would otherwise make of it.
@inline function _nodehit(nodes, p::NTuple{2,Float64})
    for k in eachindex(nodes)
        n = nodes[k]
        (n[1] == p[1] && n[2] == p[2]) && return Int(k)
    end
    return 0
end

# The turning direction of a closed node ring: `1` counter-clockwise, `-1`
# clockwise, and `0` when some corner neither turns nor exists — a repeated
# node, a straight corner, or a reflex one. A ring answering `0` is degenerate
# for every basis here, since all three need a simple convex cell.
function _orientation(nodes)
    n = length(nodes)
    n >= 3 || return 0
    turn = 0
    for i in 1:n
        a, b, c = nodes[i], nodes[_next(i, n)], nodes[_next(_next(i, n), n)]
        e1, e2 = _sub(b, a), _sub(c, b)
        cr = _cross(e1, e2)
        abs(cr) > SHAPE_TOL * max(_norm2(e1), _norm2(e2)) || return 0
        t = cr > 0 ? 1 : -1
        turn == 0 ? (turn = t) : (turn == t || return 0)
    end
    return turn
end

# Coordinates already sum to one; dividing by the total they actually reached
# keeps a reported row inside the tolerance the acceptance law states.
function _normalize!(row::WeightRow)
    total = 0.0
    for w in row.weights
        total += w
    end
    (isfinite(total) && total > 0) || (empty!(row); return WeightsDegenerate)
    if total != 1.0
        for k in eachindex(row.weights)
            row.weights[k] /= total
        end
    end
    return WeightsMapped
end

# --- the three kernels -----------------------------------------------------

"""
    linearweights!(row, indices, nodes, p) -> WeightStatus

Weight a triangle's three nodes by the barycentric coordinates of `p`.

`nodes` holds the three node coordinates in one two-dimensional chart, `p` is a
point in that chart, and `indices` names the nodes' source cells in the same
order. The row is cleared on entry and left holding the nonzero weights, so a
point on an edge emits two entries and a point at a node one.

Weights reproduce every affine field on the triangle and sum to one. A point
outside the triangle is `WeightsOutside` and is never clamped in; a triangle
with repeated or collinear nodes is `WeightsDegenerate`.
"""
function linearweights!(row::WeightRow, indices, nodes, p::NTuple{2,Float64})
    empty!(row)
    length(nodes) == 3 || return WeightsDegenerate
    a, b, c = nodes[1], nodes[2], nodes[3]
    ab, ac, bc = _sub(b, a), _sub(c, a), _sub(c, b)
    twicearea = _cross(ab, ac)
    abs(twicearea) > SHAPE_TOL * max(_norm2(ab), _norm2(ac), _norm2(bc)) ||
        return WeightsDegenerate

    hit = _nodehit(nodes, p)
    if hit != 0
        _addentry!(row, indices[hit], 1.0)
        return WeightsMapped
    end

    wa = _cross(_sub(b, p), _sub(c, p)) / twicearea
    wb = _cross(_sub(c, p), _sub(a, p)) / twicearea
    wc = _cross(_sub(a, p), _sub(b, p)) / twicearea
    (wa < -COORD_TOL || wb < -COORD_TOL || wc < -COORD_TOL) && return WeightsOutside
    _addentry!(row, indices[1], wa)
    _addentry!(row, indices[2], wb)
    _addentry!(row, indices[3], wc)
    return _normalize!(row)
end

# The bilinear map of the unit square onto a quadrilateral, and its derivatives.
@inline _q1image(n1, n2, n3, n4, u::Float64, v::Float64) =
    ((1 - u) * (1 - v) * n1[1] + u * (1 - v) * n2[1] + u * v * n3[1] +
     (1 - u) * v * n4[1],
        (1 - u) * (1 - v) * n1[2] + u * (1 - v) * n2[2] + u * v * n3[2] +
        (1 - u) * v * n4[2])

# Roundoff at the edge of the parameter square is the square's edge: a point at
# a node or on an edge emits exactly the nodes that carry it.
@inline function _snapunit(t::Float64)
    abs(t) <= COORD_TOL && return 0.0
    abs(t - 1.0) <= COORD_TOL && return 1.0
    return t
end

"""
    bilinearweights!(row, indices, nodes, p) -> WeightStatus

Weight a quadrilateral's four nodes by tensor Q1 coordinates of `p`.

`nodes` holds the four node coordinates in cyclic order in one two-dimensional
chart, `p` is a point in that chart, and `indices` names the nodes' source cells
in the same order. The isoparametric map of the unit square onto the
quadrilateral is inverted by Newton iteration from its centre, to a step below
`$NEWTON_TOL` within `$NEWTON_MAXITER` iterations; the recovered coordinates
then give the tensor weights `(1-u)(1-v)`, `u(1-v)`, `uv`, `(1-u)v`.

Weights reproduce `a + bx + cy + dxy` on an axis-aligned rectangle, reproduce
every affine field on any convex quadrilateral, and sum to one. The row is
cleared on entry and left holding the nonzero weights. A point outside the
quadrilateral is `WeightsOutside`; a folded, zero-area, reflex or repeated-node
quadrilateral, and an inverse map that fails to converge, are
`WeightsDegenerate`.
"""
function bilinearweights!(row::WeightRow, indices, nodes, p::NTuple{2,Float64})
    empty!(row)
    length(nodes) == 4 || return WeightsDegenerate
    n1, n2, n3, n4 = nodes[1], nodes[2], nodes[3], nodes[4]
    _orientation(nodes) == 0 && return WeightsDegenerate

    hit = _nodehit(nodes, p)
    if hit != 0
        _addentry!(row, indices[hit], 1.0)
        return WeightsMapped
    end

    # x(u, v) = n1 + u B + v C + u v D.
    B = _sub(n2, n1)
    C = _sub(n4, n1)
    D = (n1[1] - n2[1] + n3[1] - n4[1], n1[2] - n2[2] + n3[2] - n4[2])
    size2 = max(_norm2(B), _norm2(C), _norm2(_sub(n3, n1)))

    u = v = 0.5
    converged = false
    for _ in 1:NEWTON_MAXITER
        image = _q1image(n1, n2, n3, n4, u, v)
        r = _sub(image, p)
        du = (B[1] + v * D[1], B[2] + v * D[2])
        dv = (C[1] + u * D[1], C[2] + u * D[2])
        det = _cross(du, dv)
        abs(det) > SHAPE_TOL * size2 || return WeightsDegenerate
        su = _cross(r, dv) / det
        sv = _cross(du, r) / det
        u -= su
        v -= sv
        if abs(su) <= NEWTON_TOL && abs(sv) <= NEWTON_TOL
            converged = true
            break
        end
    end
    converged || return WeightsDegenerate
    r = _sub(_q1image(n1, n2, n3, n4, u, v), p)
    _norm2(r) <= NEWTON_RESIDUAL^2 * size2 || return WeightsDegenerate

    u, v = _snapunit(u), _snapunit(v)
    (0.0 <= u <= 1.0 && 0.0 <= v <= 1.0) || return WeightsOutside
    _addentry!(row, indices[1], (1 - u) * (1 - v))
    _addentry!(row, indices[2], u * (1 - v))
    _addentry!(row, indices[3], u * v)
    _addentry!(row, indices[4], (1 - u) * v)
    return _normalize!(row)
end

# The tangent of half the angle `a` and `b` subtend, taken positive whichever
# way the ring turns.
@inline function _tanhalfangle(turn::Int, a::NTuple{2,Float64}, b::NTuple{2,Float64})
    return tan(atan(turn * _cross(a, b), _dot(a, b)) / 2)
end

"""
    meanvalueweights!(row, indices, nodes, p) -> WeightStatus

Weight a convex polygon's nodes by the mean-value coordinates of `p`.

`nodes` holds three or more node coordinates in cyclic order in one
two-dimensional chart, `p` is a point in that chart, and `indices` names the
nodes' source cells in the same order. Node `i` takes
`(tan(α₋/2) + tan(α₊/2)) / rᵢ`, normalized, where `rᵢ` is its distance from `p`
and `α∓` are the angles it subtends at `p` with its two neighbours.

Weights are positive, reproduce every affine field on the polygon, and sum to
one; on a triangle they are the barycentric coordinates. The row is cleared on
entry and left holding the nonzero weights, so a point on an edge emits that
edge's two nodes and a point at a node one. A point outside the polygon is
`WeightsOutside`; a polygon with repeated nodes, a straight corner or a reflex
corner is `WeightsDegenerate`.
"""
function meanvalueweights!(row::WeightRow, indices, nodes, p::NTuple{2,Float64})
    empty!(row)
    n = length(nodes)
    n >= 3 || return WeightsDegenerate
    turn = _orientation(nodes)
    turn == 0 && return WeightsDegenerate

    hit = _nodehit(nodes, p)
    if hit != 0
        _addentry!(row, indices[hit], 1.0)
        return WeightsMapped
    end

    # Containment. A convex polygon meets an edge's line exactly along that
    # edge, so a point on the line is on the edge or outside the polygon.
    for i in 1:n
        j = _next(i, n)
        a, b = nodes[i], nodes[j]
        e, d = _sub(b, a), _sub(p, a)
        cr = turn * _cross(e, d)
        scale = sqrt(_norm2(e) * _norm2(d))
        cr < -COORD_TOL * scale && return WeightsOutside
        cr <= COORD_TOL * scale || continue
        t = _dot(d, e) / _norm2(e)
        (-COORD_TOL <= t <= 1 + COORD_TOL) || return WeightsOutside
        t = clamp(t, 0.0, 1.0)
        _addentry!(row, indices[i], 1.0 - t)
        _addentry!(row, indices[j], t)
        return _normalize!(row)
    end

    # Mean-value coordinates. Each node's two angles are visited once, in one
    # pass that carries the previous node's tangent.
    d1 = _sub(nodes[1], p)
    d = d1
    tprev = _tanhalfangle(turn, _sub(nodes[n], p), d1)
    isfinite(tprev) || return WeightsDegenerate
    for i in 1:n
        dnext = i == n ? d1 : _sub(nodes[i+1], p)
        t = _tanhalfangle(turn, d, dnext)
        isfinite(t) || (empty!(row); return WeightsDegenerate)
        w = (tprev + t) / sqrt(_norm2(d))
        (isfinite(w) && w > 0) || (empty!(row); return WeightsDegenerate)
        _addentry!(row, indices[i], w)
        tprev = t
        d = dnext
    end
    return _normalize!(row)
end

# --- dual cells ------------------------------------------------------------

"""
    DualCell(indices, nodes, kind::BasisKind)

The polygon of source sample sites a destination point falls in.

`indices` names its nodes' cells by the source space's local index and `nodes`
gives their coordinates in one two-dimensional chart, both in cyclic order.
`kind` names the coordinates that weight them, which is the cell's own
statement and not a reading of its node count.

The cell is read, never written, so one prepared cell is shared by concurrent
queries. It promises no fixed node count.
"""
struct DualCell{I<:AbstractVector{Int},N<:AbstractVector{NTuple{2,Float64}}}
    indices::I
    nodes::N
    kind::BasisKind

    function DualCell(indices::I, nodes::N, kind::BasisKind) where {I<:AbstractVector{Int},
        N<:AbstractVector{NTuple{2,Float64}}}
        length(indices) == length(nodes) || throw(ArgumentError(
            "a dual cell needs one source index per node, got $(length(indices)) " *
            "for $(length(nodes)) nodes"))
        return new{I,N}(indices, nodes, kind)
    end
end

"""
    nodecount(cell::DualCell) -> Int

The number of sample sites the cell interpolates between.
"""
nodecount(cell::DualCell) = length(cell.indices)

Base.show(io::IO, cell::DualCell) =
    print(io, "DualCell(", nodecount(cell), " nodes, ", cell.kind, ")")

# The answer where no dual cell exists. Returning one prepared object costs a
# query nothing.
const NO_DUALCELL = DualCell(Int[], NTuple{2,Float64}[], Linear)

"""
    dualweights!(row, cell::DualCell, p) -> WeightStatus

Weight `cell`'s nodes at chart point `p` with the coordinates its kind names.
"""
function dualweights!(row::WeightRow, cell::DualCell, p::NTuple{2,Float64})
    kind = cell.kind
    kind === Linear && return linearweights!(row, cell.indices, cell.nodes, p)
    kind === Bilinear && return bilinearweights!(row, cell.indices, cell.nodes, p)
    return meanvalueweights!(row, cell.indices, cell.nodes, p)
end

# --- samplers --------------------------------------------------------------

"""
    samplesites(space::RegridSpace) -> AbstractVector

The point each source cell's value is taken to sit at, by local index:
`samplesites(space)[i] == cellcentroid(space, i)` for every space.

The fallback reads [`cellcentroid`](@ref) on demand and holds nothing. A space
whose sites are already computed, or which knows a site other than the
centroid, returns its own vector over `1:ncells(space)`.
"""
samplesites(space::RegridSpace) = CentroidSites(space)

struct CentroidSites{S<:RegridSpace} <: AbstractVector{USPoint}
    space::S
end

Base.size(sites::CentroidSites) = (Int(ncells(sites.space)),)
Base.IndexStyle(::Type{<:CentroidSites}) = Base.IndexLinear()
Base.@propagate_inbounds function Base.getindex(sites::CentroidSites, i::Int)
    @boundscheck checkbounds(sites, i)
    return cellcentroid(sites.space, i)
end

"""
    samplerstate(space::RegridSpace) -> state

Immutable query state a [`Sampler`](@ref) over `space` carries: axis locators,
topology tables, prepared dual cells. Whatever a space returns is read
concurrently and must not be written during a sweep.

The fallback prepares a [`ChartState`](@ref) where the space has a cell chart,
which is what puts it on the fused chart path in [`weightsat!`](@ref), and is
`nothing` otherwise. A space with a state of its own returns it here and takes
the path its type names.
"""
samplerstate(space::RegridSpace) =
    hascellchart(space) ? ChartState(space) : nothing

"""
    Sampler(space, sites, state)

A source space prepared to be asked for weights at points.

It holds the space, its [`samplesites`](@ref), and whatever
[`samplerstate`](@ref) the space prepared once. It is immutable and read
concurrently; per-point scratch belongs to the calling task's
[`WeightRow`](@ref).
"""
struct Sampler{S<:RegridSpace,V<:AbstractVector,T}
    space::S
    sites::V
    state::T
end

Base.show(io::IO, s::Sampler) = print(io, "Sampler(", s.space, ")")

"""
    sampler(method, space::RegridSpace) -> Sampler

Prepare `space` to answer `method`'s point queries. This runs once per plan or
source space, never once per destination.
"""
sampler(::BarycentricPoint, space::RegridSpace) =
    Sampler(space, samplesites(space), samplerstate(space))

"""
    chartat(sampler, p) -> Union{NTuple{2,Float64},Nothing}

The coordinates of `p` in the chart the sampler's dual cells are written in, or
`nothing` where the chart does not reach. The fallback is the source space's
cell chart, so a space with no chart answers `nothing`.
"""
function chartat(s::Sampler, p)
    hascellchart(s.space) || return nothing
    c = chartcoords(s.space, p)
    c === nothing && return nothing
    x, y = Float64(c[1]), Float64(c[2])
    (isfinite(x) && isfinite(y)) || return nothing
    return (x, y)
end

"""
    dualcellat(sampler, p) -> DualCell

The dual cell of source sample sites containing `p`, or a cell with no nodes
where none does.

The fallback has no construction and answers no nodes, so a source space that
does not implement it maps nothing. A space either answers here, in the chart
its [`chartat`](@ref) uses, or answers [`weightsat!`](@ref) whole.
"""
dualcellat(::Sampler, p) = NO_DUALCELL

"""
    weightsat!(row, sampler, p) -> WeightStatus

Fill `row` with the source cells and weights the point `p` takes its value
from, and return what became of it.

The row is cleared on entry, and is left empty for every status but
`WeightsMapped`. The default locates [`dualcellat`](@ref) and weights its nodes
with the coordinates the cell's kind names; a source space with a fused
algorithm of its own replaces this method. Nothing here reads source values,
and nothing here is told about chunks: `row.indices` are the source space's
local indices, `1:ncells(space)`.
"""
function weightsat!(row::WeightRow, s::Sampler, p)
    empty!(row)
    cell = dualcellat(s, p)
    nodecount(cell) == 0 && return WeightsOutside
    q = chartat(s, p)
    q === nothing && return WeightsOutside
    return dualweights!(row, cell, q)
end

# --- weights for one chunk pair --------------------------------------------

"""
    buildweights!(coo, ::BarycentricPoint, dst_space, dst_inds, src_space, src_inds)

Interpolate between source sample sites at each destination sample site.

Each destination gets one [`weightsat!`](@ref) query, whose stencil is the whole
source space's; only the entries `src_inds` owns are emitted, so the other
source chunks emit their own shares and no chunk boundary changes a stencil. A
destination the source cannot map emits no entry at all, and the missing policy
decides what it becomes. Point samples have no coverage denominator.
"""
function buildweights!(coo::WeightCOO, method::BarycentricPoint,
    dst_space::RegridSpace, dst_inds, src_space::RegridSpace, src_inds)
    smp = sampler(method, src_space)
    sites = samplesites(dst_space)
    indexer = indexmap(src_inds)
    row = WeightRow()
    for (j, i) in enumerate(dst_inds)
        ismapped(weightsat!(row, smp, sites[Int(i)])) || continue
        for k in 1:length(row)
            c = localindex(indexer, row.indices[k])
            c == 0 && continue
            addweight!(coo, j, c, row.weights[k])
        end
    end
    return coo
end

# --- chart spaces ----------------------------------------------------------
#
# A source whose cells sit on a two-dimensional lattice needs no dual cell
# object and no inverse map: its dual cells are the axis-aligned rectangles
# between neighbouring sample sites, and the two bracketing fractions are the
# rectangle's own Q1 coordinates. Preparation is the two axes; a query is one
# chart conversion, two bracket searches, and at most four products.

"""
    ChartState(space::RegridSpace)

The two prepared chart axes of a space with a cell chart.

Each axis holds its sample coordinates in ascending order, the lattice
orientation they arrived in, and the period where the source wraps, so a query
searches a sorted vector and never rebuilds an axis. A [`Sampler`](@ref)
carrying one takes the fused chart path in [`weightsat!`](@ref).
"""
struct ChartState
    x::_ChartAxis
    y::_ChartAxis
end

function ChartState(space::RegridSpace)
    xs, ys = chartaxes(space)
    px, py = chartperiod(space)
    return ChartState(_ChartAxis(xs, px), _ChartAxis(ys, py))
end

Base.show(io::IO, st::ChartState) =
    print(io, "ChartState(", st.x.n, "×", st.y.n, " sample sites)")

# The unmapped answer, whose indices and weights are never read.
@inline _nobracket(status::WeightStatus) = (0, 0, 0.0, 0.0, status)

# Both axes speak: a point outside the source is outside whatever the other
# axis made of it, and a rim only where nothing was outside.
@inline _bracketstatus(sx::WeightStatus, sy::WeightStatus) =
    (sx === WeightsOutside || sy === WeightsOutside) ? WeightsOutside : WeightsRim

"""
    _bracket(ax::_ChartAxis, x) -> (i0, i1, w0, w1, status)

Bracket `x` between two of an axis' sample coordinates, without clamping.

`i0` and `i1` are lattice indices as the space numbers the axis, so a descending
axis answers the indices its own lookup order gives; `w0` and `w1` are the two
linear fractions. A fraction within `$COORD_TOL` of an end of its interval is
snapped to it, so a sample site reached through roundoff takes that site alone
and never a second entry of weight `1e-17`.

A periodic axis wraps exactly and is never unmapped: past the last site it
brackets that site and the first across the seam. A one-cell axis has no
interval, so it brackets its single site with itself at weights `1` and `0` and
the repeat leaves with the zero; it cannot exclude a point, and the other axis
decides.

On a nonperiodic axis the Q1 domain ends at the outermost sample sites. Beyond
them `x` is not bracketed: within half of the end interval — the rim between the
outermost sites and the source's boundary — it is `WeightsRim`, and further out
`WeightsOutside`.
"""
@inline function _bracket(ax::_ChartAxis, x::Float64)
    v = ax.values
    n = ax.n
    n == 1 &&
        return (_latticeindex(ax, 1), _latticeindex(ax, 1), 1.0, 0.0, WeightsMapped)
    period = ax.period
    if period === nothing
        if x < v[1]
            t = (x - v[1]) / (v[2] - v[1])
            abs(t) <= COORD_TOL && return (_latticeindex(ax, 1), _latticeindex(ax, 2),
                1.0, 0.0, WeightsMapped)
            return _nobracket(t >= -0.5 ? WeightsRim : WeightsOutside)
        elseif x > v[n]
            t = (x - v[n]) / (v[n] - v[n-1])
            abs(t) <= COORD_TOL && return (_latticeindex(ax, n - 1), _latticeindex(ax, n),
                0.0, 1.0, WeightsMapped)
            return _nobracket(t <= 0.5 ? WeightsRim : WeightsOutside)
        end
        k = min(searchsortedlast(v, x), n - 1)
        t = _snapunit((x - v[k]) / (v[k+1] - v[k]))
        return (_latticeindex(ax, k), _latticeindex(ax, k + 1), 1.0 - t, t, WeightsMapped)
    end
    p = period::Float64
    xw = v[1] + mod(x - v[1], p)
    k = max(searchsortedlast(v, xw), 1)
    if k < n
        t = _snapunit((xw - v[k]) / (v[k+1] - v[k]))
        return (_latticeindex(ax, k), _latticeindex(ax, k + 1), 1.0 - t, t, WeightsMapped)
    end
    t = _snapunit((xw - v[n]) / (v[1] + p - v[n]))
    return (_latticeindex(ax, n), _latticeindex(ax, 1), 1.0 - t, t, WeightsMapped)
end

"""
    weightsat!(row, s::Sampler{<:RegridSpace,<:AbstractVector,ChartState}, p)

Weight the at most four source sample sites bracketing `p` on the source's cell
chart.

The point crosses into the chart once and each prepared axis brackets it once;
the two pairs of fractions multiply into the tensor Q1 weights of the
axis-aligned dual rectangle they name. That rectangle's Q1 coordinates are the
fractions themselves, so the basis is `Bilinear` by construction and no cell is
built and no inverse map solved. Zero products are dropped, so a point on a
lattice line takes two sites and a point at a sample site takes one with weight
one; the four fractions are a partition of unity, so the row sums to one to
roundoff.

Periodic axes wrap exactly. On a nonperiodic axis the domain ends at the
outermost sample sites: a point in the rim between those sites and the source's
boundary is `WeightsRim`, a point beyond it `WeightsOutside`, and neither takes
weights.
"""
function weightsat!(row::WeightRow,
    s::Sampler{<:RegridSpace,<:AbstractVector,ChartState}, p)
    empty!(row)
    q = chartat(s, p)
    q === nothing && return WeightsOutside
    st = s.state
    ix0, ix1, wx0, wx1, sx = _bracket(st.x, q[1])
    iy0, iy1, wy0, wy1, sy = _bracket(st.y, q[2])
    (ismapped(sx) && ismapped(sy)) || return _bracketstatus(sx, sy)
    space = s.space
    for (ix, wx) in ((ix0, wx0), (ix1, wx1)),
        (iy, wy) in ((iy0, wy0), (iy1, wy1))

        w = wx * wy
        w > 0 || continue
        _addentry!(row, chartlocalindex(space, ix, iy), w)
    end
    return WeightsMapped
end

"""
    supportradius(::BarycentricPoint, src_space) -> Float64

The larger chart-axis spacing in radians on a source the method queries through
its cell chart, and `0.0` on every other source.

A chart stencil reaches out to the sample sites bracketing a destination, up to
one source cell beyond the destination's own extent, so chunk discovery must
search that far or lose weights. Which sources those are is the state they
prepare and not their type, the same reading [`weightsat!`](@ref) makes: a
source answering dual cells of its own bounds them itself, and until it does its
stencils are found by cap overlap alone.
"""
supportradius(::BarycentricPoint, src_space::RegridSpace) =
    _chartradius(samplerstate(src_space), src_space)

_chartradius(_, ::RegridSpace) = 0.0
_chartradius(::ChartState, space::RegridSpace) = supportradius(BilinearPoint(), space)
