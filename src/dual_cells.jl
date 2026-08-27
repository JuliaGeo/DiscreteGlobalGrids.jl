# Dual cells behind a `DGGSpace`.
#
# A point method interpolates between source sample sites, and for cell-centred
# DGGS data the polygon it interpolates on is the dual cell of a primal vertex:
# the sites of the cells meeting at that vertex. The cells meeting at a vertex
# of the host cell are all in the host's `Vertex()` one-ring, and the host's
# `Edge()` one-ring says where one vertex's fan ends and the next begins, so the
# host's candidates come out of two ring calls and no vertex is ever numbered.
#
# That reading needs a conforming tiling. A mixed-level container has none, and
# takes the second construction at the foot of this file: the ring is read on
# the host's own level, each member resolves to the one stored cell speaking for
# it, and the resolved sites are fanned around the host's site into triangles.

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

# One stored cell per ring member, so a fan over mixed levels is never longer
# than the ring it came from. `0` marks a direction the container covers
# nowhere.
const MixedRingIndices = SmallVector{MAX_DUAL_RING,Int}

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
    MultiOrderDualTopology(cells::MultiOrderVector)

What a mixed-level [`DGGSpace`](@ref)'s dual cells are read on: the container
itself.

A mixed-level tiling does not conform — a coarse cell's edge carries
T-junctions, and on a hexagonal hierarchy parent and child vertices never
coincide — so there is no one level to hold the topology on. Each query reads
the ring on its host cell's **own** level grid, where the ring is well defined,
and resolves every ring member to the one stored cell that speaks for it.

Held once per sampler, read and never written, so concurrent queries share it.
"""
struct MultiOrderDualTopology{ID,S<:AbstractHierarchicalGridSystem}
    cells::MultiOrderVector{ID,S}
end

Base.show(io::IO, t::MultiOrderDualTopology) =
    print(io, "MultiOrderDualTopology(", length(t.cells), " stored cells)")

_dualtopology(grid::MultiOrderGrid) = MultiOrderDualTopology(grid.cells)

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

# A mixed-level space has one chunk covering the whole sphere, so the bound
# above would report `pi` and say nothing. The stencil still reaches only one
# cell past the point, and the widest cell it can reach is the widest STORED
# one — so bound it by that instead, which is the bound that stays true once
# this space is chunked. One pass over the stored cells, once per plan.
function _dualreach(space::DGGSpace{<:MultiOrderGrid})
    mov = space.grid.cells
    sys = mov.system
    r = 0.0
    for i in eachindex(mov.cells)
        c = @inbounds mov.cells[i]
        cap = Fallbacks.cell_cap(levelgrid(sys, level(c)), c)
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

# ---------------------------------------------------------------------------
# Mixed levels.
#
# The stored cells of a `MultiOrderVector` do not tile conformingly, so the
# ring above cannot be read on them. It can be read on the host cell's OWN
# level, which is complete; each member of that ring then resolves to the one
# stored cell that speaks for it — its covering ancestor, or, where the
# container refines under it, the nearest of the stored cells beneath it. The
# resolved sites are fanned around the host's site by bearing, and the wedge
# holding the query is one triangle. Triangles only: a higher-valence cell
# needs a definition of "the cells meeting at a vertex" that survives a
# T-junction, and a uniform container never arrives here — it is a
# `PartialGrid` and keeps its own quadrilateral dual cells.
# ---------------------------------------------------------------------------

const _MixedSampler{V} = GR.Sampler{<:DGGSpace,V,<:MultiOrderDualTopology}

# The sample site of a stored cell: its centroid at its own level.
@inline function _storedsite(mov::MultiOrderVector, k::Int)
    c = @inbounds mov.cells[k]
    return cell_centroid(levelgrid(mov.system, level(c)), c)
end

# The stored cell nearest `p` among those the container holds strictly beneath
# `n`, or `0` where it holds none. Their reference-level intervals are one
# contiguous run of the container's own interval index, so the run is two
# binary searches and the scan inside it allocates nothing.
function _nearestunder(mov::MultiOrderVector, n::AbstractCellIndex, p)
    r = descendant_range(mov.system, n, mov.reference_level)
    lo = searchsortedfirst(mov.starts, first(r))
    hi = searchsortedlast(mov.starts, last(r))
    best = 0
    bestd = Inf
    for k in lo:hi
        q = _storedsite(mov, k)
        d = (q[1] - p[1])^2 + (q[2] - p[2])^2 + (q[3] - p[3])^2
        d < bestd && (bestd = d; best = k)
    end
    return best
end

"""
    GlobalRegridding.chartat(s::Sampler{<:DGGSpace,V,<:MultiOrderDualTopology}, p)

The origin, always — the mixed fan is written in the same chart centred on the
queried point that every other `DGGSpace` dual cell is.
"""
GR.chartat(::_MixedSampler, p) = (0.0, 0.0)

"""
    GlobalRegridding.dualcellat(s::Sampler{<:DGGSpace,V,<:MultiOrderDualTopology}, p)

The triangle of stored sample sites containing `p`, or a cell with no nodes.

`p`'s host is the stored cell covering it. The ring of the host's own level
resolves to one stored site per member, and the two whose bearings from the
host's site bracket the query's are the triangle's other two nodes. The kind is
`MeanValue`, which on three nodes is that triangle's barycentric coordinates.
"""
GR.dualcellat(s::_MixedSampler, p) = first(_locatemixeddual(s, p))

"""
    GlobalRegridding.weightsat!(row, s::Sampler{<:DGGSpace,V,<:MultiOrderDualTopology}, p)

Weight the stored sample sites of the triangle holding `p`.

The row is cleared on entry and left empty for every status but
`WeightsMapped`: `WeightsOutside` where the container covers `p` nowhere or the
fan does not reach it, `WeightsRim` where the wedge holding `p` needs a cell the
container holds nowhere, and `WeightsDegenerate` where the fan is unusable —
nonlocal on the sphere, longer than a ring may be, or fewer than two distinct
sites.
"""
function GR.weightsat!(row::GR.WeightRow, s::_MixedSampler, p)
    empty!(row)
    cell, status = _locatemixeddual(s, p)
    GR.nodecount(cell) == 0 && return status
    return GR.dualweights!(row, cell, (0.0, 0.0))
end

# The triangle holding `p` and what became of the search, as `_locatedual` is
# for one level.
function _locatemixeddual(s::_MixedSampler, p)
    mov = s.state.cells
    sys = mov.system
    host = localindex(mov, p)
    # A hole in the container is a real answer. No stored cell stands in for
    # one, so there is nothing to build a stencil around.
    host === nothing && return (NO_GRID_DUALCELL, GR.WeightsOutside)
    hostid = @inbounds mov.cells[host]
    hostgrid = levelgrid(sys, level(hostid))
    chart = GR.TangentChart(p)
    hostsite = GR.chartpoint(chart, cell_centroid(hostgrid, hostid))
    hostsite === nothing && return (NO_GRID_DUALCELL, GR.WeightsDegenerate)

    # The ring is read on the host's OWN level, which is complete and conforms.
    ring = neighbors(hostgrid, hostid, 1; connectivity = Vertex())
    m = length(ring)
    m <= MAX_DUAL_RING || return (NO_GRID_DUALCELL, GR.WeightsDegenerate)

    sites = DualRing()
    ids = MixedRingIndices()
    for j in 1:m
        n = @inbounds ring[j]
        k = covering_index(mov, n)
        i = k === nothing ? _nearestunder(mov, n, p) : k
        q = if i == 0
            # Nothing stored covers this direction. The member's own centroid
            # is not a site — it has no local index — but it marks where the
            # gap lies, so a wedge needing it can be told apart from one that
            # does not.
            GR.chartpoint(chart, cell_centroid(hostgrid, n))
        else
            # Several members share one coarse ancestor; that is one site, not
            # a repeated node.
            (i == host || _inring(ids, i)) && continue
            GR.chartpoint(chart, _storedsite(mov, i))
        end
        q === nothing && return (NO_GRID_DUALCELL, GR.WeightsDegenerate)
        ids = SmallCollections.push(ids, i)
        sites = SmallCollections.push(sites, q)
    end
    n = length(ids)
    n >= 2 || return (NO_GRID_DUALCELL, GR.WeightsDegenerate)

    # The wedge holding the query. The query is the chart's origin, so the
    # bearing to bracket is the one from the host's site back to it, and the
    # two sites whose bearings bracket it are neighbours in the fan. No
    # ordering of the whole ring is ever formed.
    hx, hy = hostsite
    base = atan(-hy, -hx)
    lo = hi = 0
    dlo, dhi = -1.0, 7.0
    for k in 1:n
        q = @inbounds sites[k]
        d = mod(atan(q[2] - hy, q[1] - hx) - base, 2 * Float64(pi))
        d < dhi && (dhi = d; hi = k)
        d > dlo && (dlo = d; lo = k)
    end
    lo == hi && return (NO_GRID_DUALCELL, GR.WeightsDegenerate)

    a, b = @inbounds(ids[lo]), @inbounds(ids[hi])
    # Closing the wedge across a gap would build a triangle over ground the
    # container does not hold, and weight it as if it did.
    (a == 0 || b == 0) && return (NO_GRID_DUALCELL, GR.WeightsRim)

    nodes = DualNodes()
    nodes = SmallCollections.push(nodes, hostsite)
    nodes = SmallCollections.push(nodes, @inbounds sites[lo])
    nodes = SmallCollections.push(nodes, @inbounds sites[hi])
    GR.containspoint(nodes, (0.0, 0.0)) ||
        return (NO_GRID_DUALCELL, GR.WeightsOutside)

    indices = DualIndices()
    indices = SmallCollections.push(indices, host)
    indices = SmallCollections.push(indices, a)
    indices = SmallCollections.push(indices, b)
    return (GR.DualCell(indices, nodes, GR.MeanValue), GR.WeightsMapped)
end
