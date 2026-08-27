# Dual cells behind a `DGGSpace`.
#
# A point method interpolates between source sample sites, and for cell-centred
# DGGS data the polygon it interpolates on is the dual cell of a primal vertex:
# the sites of the cells meeting at that vertex. The cells meeting at a vertex
# of the host cell are all in the host's `Vertex()` one-ring, and the host's
# `Edge()` one-ring says where one vertex's fan ends and the next begins, so the
# host's candidates come out of two ring calls and no vertex is ever numbered.

# The most sample sites one dual cell holds. A dual cell of a primal vertex has
# one node per cell meeting there, so this is the largest primal-vertex valence
# supported; a fan longer than this is `WeightsDegenerate` rather than a
# truncated stencil.
const MAX_DUAL_NODES = 12

# The most cells a host's `Vertex()` one-ring may hold. The largest
# `maxneighbors` any system here declares is eleven; a longer ring is
# `WeightsDegenerate`.
const MAX_DUAL_RING = 12

const DualIndices = SmallVector{MAX_DUAL_NODES,Int}
const DualNodes = SmallVector{MAX_DUAL_NODES,NTuple{2,Float64}}
const DualRing = SmallVector{MAX_DUAL_RING,NTuple{2,Float64}}

"""
    GridDualCell

The [`GlobalRegridding.DualCell`](@ref) a [`DGGSpace`](@ref) answers: fixed
capacity in both its source indices and its chart coordinates, so locating one
reaches no heap.
"""
const GridDualCell = GR.DualCell{DualIndices,DualNodes}

# The answer where no dual cell holds the point. Its status says why.
const NO_GRID_DUALCELL = GR.DualCell(DualIndices(), DualNodes(), GR.MeanValue)

"""
    DualTopology(grid)

The complete level a [`DGGSpace`](@ref)'s dual cells read topology on.

A space over a subset reads rings and sample sites here rather than on the
subset, so a cell the subset does not hold is a missing node — a rim — instead
of a ring silently one cell shorter.
"""
struct DualTopology{G<:AbstractGrid}
    grid::G
end

Base.show(io::IO, t::DualTopology) = print(io, "DualTopology(", t.grid, ")")

_dualtopology(grid::AbstractGrid) = DualTopology(grid)
_dualtopology(grid::PartialGrid) = DualTopology(grid.complete)

# A `DGGSpace` builds dual cells from its grid's own adjacency.
GR.hasdualcells(::DGGSpace) = true

"""
    GlobalRegridding.hasdualcells(::DGGSpace{<:MultiOrderGrid}) -> false

No dual cells over mixed levels — yet.

`_locatedual` below reads the host's `Vertex()` and `Edge()` one-rings, and a
mixed-level tiling has neither: it does not conform, so a ring matched on shared
vertices is wrong across every T-junction. A point method asking this space for
a stencil gets [`GlobalRegridding.sampler`](@ref)'s own refusal rather than a
plausible fan built on the wrong topology.
"""
GR.hasdualcells(::DGGSpace{<:MultiOrderGrid}) = false

"""
    GlobalRegridding.sampler(::BarycentricPoint, space::DGGSpace{<:MultiOrderGrid})

Refuse, and say what to do instead.

The generic refusal reads `hasdualcells` and reports the type; this one names
the reason and the two working routes. Interpolating natively over mixed levels
is the dual-cell construction this file does not have yet, and interpolating on
the reference-level expansion is a different, worse function — every leaf under
a stored cell repeats one value, so the blend rebuilds the coarsening staircase
at leaf spacing. Both are said here rather than left to be discovered.
"""
GR.sampler(::GR.BarycentricPoint, space::DGGSpace{<:MultiOrderGrid}) = _nomixeddual(space)

@noinline _nomixeddual(space::DGGSpace{<:MultiOrderGrid}) = throw(ArgumentError(
    "BarycentricPoint has nothing to interpolate between on the $(ncells(space)) " *
    "stored cells of a mixed-level container: they do not tile conformingly, so " *
    "there is no ring of neighbouring sample sites to build a dual cell from, " *
    "and the native mixed-level construction has not landed yet. Use an area or " *
    "nearest-cell method, which read the stored cells as they are. To " *
    "interpolate on the leaves anyway — a different function, which rebuilds " *
    "the coarsening steps at leaf spacing — expand first, which says so: " *
    "`regrid(expand(A, DiscreteGlobalGrids.reference_level(lookup(A, Cells))); " *
    "to = ..., method = BarycentricPoint())`."))

"""
    GlobalRegridding.samplerstate(space::DGGSpace)

The complete level the space's dual cells are built on, held once per sampler.

It is read and never written, so concurrent queries share it; a query's own
working buffers are fixed-capacity and live on its stack.
"""
GR.samplerstate(space::DGGSpace) = _dualtopology(space.grid)

const _DualSampler{V} = GR.Sampler{<:DGGSpace,V,<:DualTopology}

"""
    GlobalRegridding.chartat(s::Sampler{<:DGGSpace}, p)

The origin, always.

The chart a `DGGSpace`'s dual cells are written in is the azimuthal-equidistant
chart centred on the queried point itself, so the point sits at `(0.0, 0.0)` and
the nodes carry their true geodesic distances and bearings from it.
"""
GR.chartat(::_DualSampler, p) = (0.0, 0.0)

"""
    GlobalRegridding.dualcellat(s::Sampler{<:DGGSpace}, p) -> GridDualCell

The dual cell of source sample sites containing `p`, or a cell with no nodes.

`p`'s host cell is located once, and the candidates are the dual cells of the
host's own primal vertices — one per cell of its `Edge()` one-ring, whose nodes
are the host's site followed by the run of its `Vertex()` one-ring that fans
around that vertex. The candidate holding `p` supplies the stencil; no cell
polygon is matched and no vertex is numbered.

The kind is always `MeanValue`, which on the three nodes a hexagonal source
gives is that triangle's barycentric coordinates.
"""
GR.dualcellat(s::_DualSampler, p) = first(_locatedual(s, p))

"""
    GlobalRegridding.weightsat!(row, s::Sampler{<:DGGSpace}, p)

Weight the sample sites of the dual cell holding `p`.

[`dualcellat`](@ref) locates the cell and its own kind weights the nodes, at the
origin of the chart centred on `p`. Where no cell holds `p` the row stays empty
and the status says why: `WeightsRim` where the cell exists on the level but one
of its nodes is outside this space's collection, `WeightsOutside` where the
space covers `p` nowhere, and `WeightsDegenerate` where the fan is unusable —
nonlocal on the sphere, longer than a dual cell may be, or not a simple cell.
"""
function GR.weightsat!(row::GR.WeightRow, s::_DualSampler, p)
    empty!(row)
    cell, status = _locatedual(s, p)
    GR.nodecount(cell) == 0 && return status
    return GR.dualweights!(row, cell, (0.0, 0.0))
end

"""
    GlobalRegridding.supportradius(::BarycentricPoint, space::DGGSpace)

Twice the widest chunk cover, in radians.

A dual cell's nodes are the sites of cells touching the point's host cell, so
the stencil reaches at most one cell width past the point. Every cell lies
inside its own chunk's cover, so twice the widest of those covers bounds any
cell's width and chunk discovery keeps every chunk a stencil can name.
"""
GR.supportradius(::GR.BarycentricPoint, space::DGGSpace) = _dualreach(space)

function _dualreach(space::DGGSpace)
    r = 0.0
    for cap in space.caps
        r = max(r, Float64(cap.radius))
    end
    return min(Float64(pi), 2 * r)
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

# Ring positions as bits, so a walk carries them in a register instead of a
# buffer. Position `k` is bit `k - 1`.
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

# The dual cell holding `p` and what became of the search. `dualcellat` and
# `weightsat!` are the same query: one wants the cell, the other also wants the
# reason there is none.
function _locatedual(s::_DualSampler, p)
    grid = s.space.grid
    topo = s.state.grid
    host = localindex(grid, p)
    if host === nothing
        # A point the level covers but the collection does not: every dual cell
        # around it needs the host's own site, which this space has no index
        # for, so it is a rim rather than a cell built from someone else's.
        return (NO_GRID_DUALCELL,
            localindex(topo, p) === nothing ? GR.WeightsOutside : GR.WeightsRim)
    end
    hostid = cellindex(grid, host)
    chart = GR.TangentChart(p)
    hostsite = GR.chartpoint(chart, cell_centroid(topo, hostid))
    hostsite === nothing && return (NO_GRID_DUALCELL, GR.WeightsDegenerate)

    fan = neighbors(topo, hostid, 1; connectivity = Vertex())
    m = length(fan)
    m <= MAX_DUAL_RING || return (NO_GRID_DUALCELL, GR.WeightsDegenerate)
    sites = DualRing()
    nonlocal = zero(UInt16)
    for k in 1:m
        c = GR.chartpoint(chart, cell_centroid(topo, @inbounds fan[k]))
        c === nothing && (nonlocal = _set(nonlocal, k))
        sites = SmallCollections.push(sites, c === nothing ? (0.0, 0.0) : c)
    end

    # Where one vertex's fan ends and the next begins: the cells sharing a whole
    # edge with the host, in the positions the vertex ring already put them in.
    cuts = zero(UInt16)
    edges = neighbors(topo, hostid, 1; connectivity = Edge())
    for k in 1:m
        _inring(edges, @inbounds fan[k]) && (cuts = _set(cuts, k))
    end
    count_ones(cuts) >= 2 || return (NO_GRID_DUALCELL, GR.WeightsDegenerate)

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
        GR.containspoint(nodes, (0.0, 0.0)) || continue

        indices = DualIndices()
        indices = SmallCollections.push(indices, host)
        k = a
        while true
            i = localindex(grid, @inbounds fan[k])
            # A node the collection does not hold has no local index at all, so
            # this cell cannot be weighted here and nothing stands in for it.
            i === nothing && return (NO_GRID_DUALCELL, GR.WeightsRim)
            indices = SmallCollections.push(indices, i)
            k == b && break
            k = _cyclicnext(k, m)
        end
        return (GR.DualCell(indices, nodes, GR.MeanValue), GR.WeightsMapped)
    end
    return (NO_GRID_DUALCELL, unusable ? GR.WeightsDegenerate : GR.WeightsOutside)
end
