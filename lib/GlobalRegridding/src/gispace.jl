# A GridInterface grid as a `RegridSpace`: every hook is one verb call.
#
# The table in this file is the adapt step of the GridInterface migration. The
# methods, plans, executor and chunk graph run unchanged against it, and the
# parity suite (`test/test_gispace_parity.jl`) compares the wrapped path with
# the spaces it stands in for. Nothing here is precomputed per cell: the
# constructor reads `chunks(grid)` and holds it, and every other answer is
# derived when it is asked for.

"""
    GISpace(grid)

Wrap a GridInterface grid (`GridInterface.isgrid(grid)`) as a
[`RegridSpace`](@ref).

Every space hook is the GridInterface verb of the same meaning: cells are
`cell_geometry`, the cell tree is `treeify(grid; leaf_capacity = ForSweeps())`
and a chunk's restricted tree `subtree`, chunks are `chunks(grid)` (one chunk
over every cell when the grid declares none), point lookup is `cellat`, charts
are the `Structured` face axes of a one-face grid with separable coordinates,
and dual cells are read off `neighbors` and `cell_centroid` on a conforming
grid. Construction holds the grid and its chunk collection and computes nothing
per cell.
"""
struct GISpace{G,C} <: RegridSpace
    grid::G
    chunks::C
end

function GISpace(grid)
    GridInterface.isgrid(grid) || throw(ArgumentError(
        "GISpace wraps a GridInterface grid (`GridInterface.isgrid`), got " *
        "$(typeof(grid))"))
    return GISpace(grid, GridInterface.chunks(grid))
end

GISpace(space::GISpace) = space

Base.show(io::IO, s::GISpace) =
    print(io, "GISpace(", s.grid, ", ", nchunks(s), s.chunks === nothing ?
        " chunk)" : " chunks)")

# --------------------------------------------------------------------------
# Cell geometry
# --------------------------------------------------------------------------

ncells(s::GISpace) = Int(GridInterface.ncells(s.grid))
getcell(s::GISpace, i::Int) = GridInterface.cell_geometry(s.grid, i)
expensivecellgeometry(s::GISpace) = GridInterface.cell_geometry_is_expensive(s.grid)
manifold(s::GISpace) = GOCore.best_manifold(s.grid)

celltree(s::GISpace) =
    GridInterface.treeify(s.grid; leaf_capacity = GridInterface.ForSweeps())

# The whole space is its tree; a chunk's own positions take the chunk, so a
# rooted hierarchy chunk starts its cursor at the root in O(1); any other index
# vector is the grid's `subtree` (a tile's lazy index vector dispatches to its
# quadtree block there). Only an index set that is not a vector reaches the
# packed cell-space fallback.
function subtree(s::GISpace, inds::AbstractUnitRange{<:Integer})
    _iswholespace(s, inds) && return celltree(s)
    chunk = _ownchunk(s, inds)
    chunk === nothing || return GridInterface.subtree(s.grid, chunk)
    return GridInterface.subtree(s.grid, inds)
end
subtree(s::GISpace, inds::AbstractVector{<:Integer}) = GridInterface.subtree(s.grid, inds)
subtree(s::GISpace, inds) = CellSpaceRTree(s, inds)

# The chunk whose positions are exactly `inds`, or `nothing`.
function _ownchunk(s::GISpace, inds::AbstractUnitRange{<:Integer})
    ch = s.chunks
    (ch === nothing || isempty(inds)) && return nothing
    c = GridInterface.chunkat(ch, Int(first(inds)))
    chunk = ch[c]
    owned = GridInterface.parentpositions(chunk)
    owned isa AbstractUnitRange || return nothing
    return (Int(first(owned)) == Int(first(inds)) &&
            Int(last(owned)) == Int(last(inds))) ? chunk : nothing
end

# --------------------------------------------------------------------------
# Chunks
# --------------------------------------------------------------------------

nchunks(s::GISpace) = s.chunks === nothing ? 1 : length(s.chunks)

function ownedindices(s::GISpace, chunk::Int)
    s.chunks === nothing || return GridInterface.parentpositions(s.chunks[chunk])
    chunk == 1 || throw(BoundsError(s, chunk))
    return 1:ncells(s)
end

# One extent per chunk, read when asked: `grid_extent` of a chunk is cheap by
# contract (stored or O(1)), so nothing is held.
chunkextents(s::GISpace) = s.chunks === nothing ?
    Cap[GridInterface.grid_extent(s.grid)] :
    Cap[GridInterface.grid_extent(c) for c in s.chunks]

function chunkextent(s::GISpace, chunk::Integer)
    s.chunks === nothing ||
        return convert(Cap, GridInterface.grid_extent(s.chunks[Int(chunk)]))
    Int(chunk) == 1 || throw(BoundsError(s, Int(chunk)))
    return convert(Cap, GridInterface.grid_extent(s.grid))
end

function chunkat(s::GISpace, i::Integer)
    1 <= i <= ncells(s) || throw(BoundsError(s, Int(i)))
    s.chunks === nothing && return 1
    return Int(GridInterface.chunkat(s.chunks, Int(i)))
end

# `chunkindex` and `candidatechunks!` stay the generic ones: the chunk extents
# packed in an R-tree and tested with `DilatedIntersects`, so the relation is
# the cap join exactly.

# --------------------------------------------------------------------------
# Point lookup and neighbourhoods
# --------------------------------------------------------------------------

cellat(s::GISpace, p) = GridInterface.cellat(s.grid, p)
cellcentroid(s::GISpace, i::Int) = GridInterface.cell_centroid(s.grid, i)
hasanalyticlocation(s::GISpace) = GridInterface.location_is_analytic(s.grid)
cellneighbors(s::GISpace, i::Int) = GridInterface.neighbors(s.grid, i)
cellcap(s::GISpace, i::Int) = GridInterface.cell_extent(s.grid, i)
cellcorners(s::GISpace, i::Int) = GridInterface.cell_corners(s.grid, i)

# --------------------------------------------------------------------------
# Cell charts: the face axes of a one-face structured grid
# --------------------------------------------------------------------------

_chartface(grid) = isstructured(grid) && nfaces(grid) == 1 && hasfaceaxes(grid, 1)

hascellchart(s::GISpace) = _chartface(s.grid)

# Face axes hold vertex coordinates; the chart samples at cell centres.
_axiscentres(xs) = [(Float64(xs[k]) + Float64(xs[k+1])) / 2 for k in 1:(length(xs)-1)]
_axismaxstep(xs) = maximum(abs(Float64(xs[k+1]) - Float64(xs[k])) for k in 1:(length(xs)-1))

function chartaxes(s::GISpace)
    xs, ys = faceaxes(s.grid, 1)
    return (_axiscentres(xs), _axiscentres(ys))
end

function chartcoords(s::GISpace, p)
    xy = facecoords(s.grid, 1, p)
    xy === nothing && return nothing
    return (Float64(xy[1]), Float64(xy[2]))
end

chartlocalindex(s::GISpace, ix::Int, iy::Int) = Int(faceposition(s.grid, 1, ix, iy))

# A face period is the chart's (360 for longitude whatever the span); the chart
# wraps only where the axis spans it, as a regional raster does not.
function chartperiod(s::GISpace)
    xs, ys = faceaxes(s.grid, 1)
    px, py = faceperiod(s.grid, 1)
    return (_spanperiod(xs, px), _spanperiod(ys, py))
end

_spanperiod(xs, ::Nothing) = nothing
_spanperiod(xs, p) =
    isapprox(abs(Float64(last(xs)) - Float64(first(xs))), Float64(p); rtol = 1e-9) ?
    Float64(p) : nothing

chartspacing(s::GISpace) = _facespacing(s.grid)

# In a geographic chart a step of Δλ degrees bounds the distance along a
# parallel and Δφ the distance along a meridian. Any other chart names its own
# bound through the grid.
function _facespacing(grid)
    px, py = faceperiod(grid, 1)
    (px !== nothing && Float64(px) == 360.0 && py === nothing) || throw(ArgumentError(
        "chartspacing of a $(typeof(grid)) needs a bound in radians on the " *
        "distance between adjacent cell centres along each face axis; only a " *
        "geographic chart (faceperiod (360, nothing)) derives one, so define " *
        "`GlobalRegridding._facespacing(::$(typeof(grid)))`"))
    xs, ys = faceaxes(grid, 1)
    return (deg2rad(_axismaxstep(xs)), deg2rad(_axismaxstep(ys)))
end

# --------------------------------------------------------------------------
# Dual cells of a conforming grid
# --------------------------------------------------------------------------
#
# A point method interpolates between source sample sites, and on a conforming
# cell-centred grid the polygon it interpolates on is the dual cell of a primal
# vertex: the sites of the cells meeting there. The cells meeting at a vertex of
# the host cell are all in the host's `Vertex()` one-ring, and the host's
# `Edge()` one-ring says where one vertex's fan ends and the next begins, so the
# host's candidates come out of two ring calls and no vertex is ever numbered.

hasdualcells(s::GISpace) = GridInterface.conformity(s.grid) isa GridInterface.Conforming

# The most sample sites one dual cell holds: one node per cell meeting at a
# primal vertex. A longer fan is `WeightsDegenerate`, never a truncated stencil.
const MAX_DUAL_NODES = 12

# The most cells a host's `Vertex()` one-ring may hold.
const MAX_DUAL_RING = 12

const DualIndices = SmallVector{MAX_DUAL_NODES,Int}
const DualNodes = SmallVector{MAX_DUAL_NODES,NTuple{2,Float64}}
const DualRing = SmallVector{MAX_DUAL_RING,NTuple{2,Float64}}

"""
    GridDualCell

The [`DualCell`](@ref) a [`GISpace`](@ref) answers: fixed capacity in both its
source indices and its chart coordinates, so locating one reaches no heap.
"""
const GridDualCell = DualCell{DualIndices,DualNodes}

# The answer where no dual cell holds the point. Its status says why.
const NO_GRID_DUALCELL = DualCell(DualIndices(), DualNodes(), MeanValue)

"""
    DualState(grid)

What a [`GISpace`](@ref)'s dual sampler reads: the space's grid, whose positions
name the row, and the grid its rings and sample sites are read on. A subset of a
grid (`gridrelation` is `ParentGrid`) reads them on its parent, so a cell the
subset does not hold is a missing node — a rim — rather than a ring silently one
cell shorter.
"""
struct DualState{G,T}
    grid::G
    topo::T
end

DualState(grid) = DualState(grid, _dualtopology(grid))

_dualtopology(grid) =
    GridInterface.gridrelation(grid) isa GridInterface.ParentGrid ? parent(grid) : grid

Base.show(io::IO, st::DualState) = print(io, "DualState(", st.topo, ")")

"""
    samplerstate(s::GISpace)

A [`ChartState`](@ref) where the grid has a cell chart, a [`DualState`](@ref)
where it has dual cells, and `nothing` otherwise. Read and never written.
"""
samplerstate(s::GISpace) =
    hascellchart(s) ? ChartState(s) : hasdualcells(s) ? DualState(s.grid) : nothing

const _GIDualSampler = Sampler{<:GISpace,<:AbstractVector,<:DualState}

# The chart a dual cell is written in is the azimuthal-equidistant chart centred
# on the queried point itself, so the point is the origin and the nodes carry
# their true geodesic distances and bearings from it.
chartat(::_GIDualSampler, p) = (0.0, 0.0)

dualcellat(s::_GIDualSampler, p) = first(_locatedual(s, p))

function weightsat!(row::WeightRow, s::_GIDualSampler, p)
    empty!(row)
    cell, status = _locatedual(s, p)
    nodecount(cell) == 0 && return status
    return dualweights!(row, cell, (0.0, 0.0))
end

"""
    supportradius(::BarycentricPoint, s::GISpace)

The chart's larger axis spacing on a chart source, and twice the widest chunk
cover on a dual source: a dual cell's nodes are the sites of cells touching the
point's host cell, so the stencil reaches at most one cell width past the
point, and every cell lies inside its own chunk's cover.
"""
supportradius(::BarycentricPoint, s::GISpace) = _samplerreach(samplerstate(s), s)

_samplerreach(::ChartState, s::GISpace) = chartradius(s)
_samplerreach(::Nothing, ::GISpace) = 0.0
function _samplerreach(::DualState, s::GISpace)
    r = 0.0
    for cap in chunkextents(s)
        r = max(r, Float64(cap.radius))
    end
    return min(Float64(pi), 2 * r)
end

# A position of the space's grid as a position of its topology grid, and back.
_topoposition(grid, topo, i::Int) =
    topo === grid ? i : Int(GridInterface.parentpositions(grid)[i])
function _gridposition(grid, topo, q)
    topo === grid && return Int(q)
    i = GridInterface.localindex(grid, GridInterface.cellindex(topo, Int(q)))
    return i === nothing ? nothing : Int(i)
end

# Is `c` one of the cells in `ring`? Rings are shorter than a dozen cells and
# hold no duplicates, so the scan beats anything with a table behind it.
@inline function _inring(ring, c)
    for x in ring
        x == c && return true
    end
    return false
end

@inline _cyclicnext(k::Int, n::Int) = k == n ? 1 : k + 1

# Ring positions as bits, so a walk carries them in a register. Position `k` is
# bit `k - 1`.
@inline _isset(bits::UInt16, k::Int) = (bits >> (k - 1)) & one(UInt16) == one(UInt16)
@inline _set(bits::UInt16, k::Int) = bits | (one(UInt16) << (k - 1))

# The next set position after `k`, going round. Only called where at least one
# is set, so the walk always lands.
@inline function _nextset(bits::UInt16, k::Int, n::Int)
    j = k
    for _ in 1:n
        j = _cyclicnext(j, n)
        _isset(bits, j) && return j
    end
    return k
end

# The dual cell holding `p` and what became of the search: `p`'s host cell is
# located once, and the candidates are the dual cells of the host's own primal
# vertices — one per cell of its `Edge()` one-ring, whose nodes are the host's
# site followed by the run of its `Vertex()` one-ring that fans around that
# vertex. The candidate holding `p` supplies the stencil. The kind is always
# `MeanValue`, which on three nodes is the triangle's barycentric coordinates.
function _locatedual(s::_GIDualSampler, p)
    grid = s.state.grid
    topo = s.state.topo
    host = GridInterface.cellat(grid, p)
    if host === nothing
        # A point the topology covers but the collection does not: every dual
        # cell around it needs the host's own site, which this space has no
        # index for, so it is a rim rather than a cell built from someone else's.
        outside = topo === grid || GridInterface.cellat(topo, p) === nothing
        return (NO_GRID_DUALCELL, outside ? WeightsOutside : WeightsRim)
    end
    host = Int(host)
    hostpos = _topoposition(grid, topo, host)
    chart = TangentChart(p)
    hostsite = chartpoint(chart, GridInterface.cell_centroid(topo, hostpos))
    hostsite === nothing && return (NO_GRID_DUALCELL, WeightsDegenerate)

    fan = GridInterface.neighbors(topo, hostpos; connectivity = GridInterface.Vertex())
    m = length(fan)
    m <= MAX_DUAL_RING || return (NO_GRID_DUALCELL, WeightsDegenerate)
    sites = DualRing()
    nonlocal = zero(UInt16)
    for k in 1:m
        c = chartpoint(chart, GridInterface.cell_centroid(topo, @inbounds fan[k]))
        c === nothing && (nonlocal = _set(nonlocal, k))
        sites = SmallCollections.push(sites, c === nothing ? (0.0, 0.0) : c)
    end

    # Where one vertex's fan ends and the next begins: the cells sharing a whole
    # edge with the host, in the positions the vertex ring already put them in.
    cuts = zero(UInt16)
    edges = GridInterface.neighbors(topo, hostpos; connectivity = GridInterface.Edge())
    for k in 1:m
        _inring(edges, @inbounds fan[k]) && (cuts = _set(cuts, k))
    end
    count_ones(cuts) >= 2 || return (NO_GRID_DUALCELL, WeightsDegenerate)

    unusable = false
    for a in 1:m
        _isset(cuts, a) || continue
        b = _nextset(cuts, a, m)
        nodes = DualNodes()
        nodes = SmallCollections.push(nodes, hostsite)
        k = a
        usable = true
        while true
            if _isset(nonlocal, k) || length(nodes) == MAX_DUAL_NODES
                usable = false
                break
            end
            nodes = SmallCollections.push(nodes, @inbounds sites[k])
            k == b && break
            k = _cyclicnext(k, m)
        end
        usable || (unusable = true; continue)
        containspoint(nodes, (0.0, 0.0)) || continue

        indices = DualIndices()
        indices = SmallCollections.push(indices, host)
        k = a
        while true
            i = _gridposition(grid, topo, @inbounds fan[k])
            # A node the collection does not hold has no local index at all, so
            # this cell cannot be weighted here and nothing stands in for it.
            i === nothing && return (NO_GRID_DUALCELL, WeightsRim)
            indices = SmallCollections.push(indices, i)
            k == b && break
            k = _cyclicnext(k, m)
        end
        return (DualCell(indices, nodes, MeanValue), WeightsMapped)
    end
    return (NO_GRID_DUALCELL, unusable ? WeightsDegenerate : WeightsOutside)
end
