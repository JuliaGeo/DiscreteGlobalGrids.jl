# Fixed capacity keeps valid dual-cell queries allocation-free.
const MAX_DUAL_NODES = 12

# This exceeds every system's declared `maxneighbors` value.
const MAX_DUAL_RING = 12

const DualIndices = SmallVector{MAX_DUAL_NODES,Int}
const DualNodes = SmallVector{MAX_DUAL_NODES,NTuple{2,Float64}}
const DualRing = SmallVector{MAX_DUAL_RING,NTuple{2,Float64}}

# `0` marks a ring direction missing from a mixed-level container.
const MixedRingIndices = SmallVector{MAX_DUAL_RING,Int}

"""
    GridDualCell

The fixed-capacity [`GlobalRegridding.DualCell`](@ref) returned by a
[`DGGSpace`](@ref). Locating one allocates no heap storage.
"""
const GridDualCell = GR.DualCell{DualIndices,DualNodes}

const NO_GRID_DUALCELL = GR.DualCell(DualIndices(), DualNodes(), GR.MeanValue)

"""
    DualTopology(grid)

Store the complete level that builds a [`DGGSpace`](@ref)'s dual cells. The
complete topology preserves missing-node detection along subset rims.
"""
struct DualTopology{G<:AbstractGrid}
    grid::G
end

Base.show(io::IO, t::DualTopology) = print(io, "DualTopology(", t.grid, ")")

_dualtopology(grid::AbstractGrid) = DualTopology(grid)
_dualtopology(grid::PartialGrid) = DualTopology(grid.complete)

GR.hasdualcells(::DGGSpace) = true

"""
    MultiOrderDualTopology(cells::MultiOrderVector)

Hold immutable mixed-level topology for concurrent sampler queries. Each query
reads its host's conforming level ring and resolves members to stored cells.
"""
struct MultiOrderDualTopology{ID,S<:AbstractHierarchicalGridSystem}
    cells::MultiOrderVector{ID,S}
end

Base.show(io::IO, t::MultiOrderDualTopology) =
    print(io, "MultiOrderDualTopology(", length(t.cells), " stored cells)")

_dualtopology(grid::MultiOrderGrid) = MultiOrderDualTopology(grid.cells)

"""
    GlobalRegridding.samplerstate(space::DGGSpace)

Create the immutable topology shared by the sampler's concurrent queries.
Each query uses its own fixed-capacity working buffers.
"""
GR.samplerstate(space::DGGSpace) = _dualtopology(space.grid)

const _DualSampler{V} = GR.Sampler{<:DGGSpace,V,<:DualTopology}

"""
    GlobalRegridding.chartat(s::Sampler{<:DGGSpace}, p)

Return the origin `(0.0, 0.0)`. A `DGGSpace` writes dual cells in the
azimuthal-equidistant chart centered on the query point, preserving each
node's geodesic distance and bearing.
"""
GR.chartat(::_DualSampler, p) = (0.0, 0.0)

"""
    GlobalRegridding.dualcellat(s::Sampler{<:DGGSpace}, p) -> GridDualCell

Return the `MeanValue` dual cell containing `p`, or an empty cell. Edge
neighbors delimit candidate fans within the host's vertex ring.
"""
GR.dualcellat(s::_DualSampler, p) = first(_locatedual(s, p))

"""
    GlobalRegridding.weightsat!(row, s::Sampler{<:DGGSpace}, p)

Weight the dual cell containing `p` at the query-centered chart origin. An
unsuccessful query leaves the row empty.

# Dual-cell weight statuses

  - `WeightsRim`: a required node lies outside the collection.
  - `WeightsOutside`: the space leaves `p` uncovered.
  - `WeightsDegenerate`: the fan exceeds capacity, crosses the chart's
    antipode, or fails to form a simple cell.
"""
function GR.weightsat!(row::GR.WeightRow, s::_DualSampler, p)
    empty!(row)
    cell, status = _locatedual(s, p)
    GR.nodecount(cell) == 0 && return status
    return GR.dualweights!(row, cell, (0.0, 0.0))
end

"""
    GlobalRegridding.supportradius(::BarycentricPoint, space::DGGSpace)

Return twice the widest chunk-cover radius. This bounds a stencil that reaches
one cell width beyond its host.
"""
GR.supportradius(::GR.BarycentricPoint, space::DGGSpace) = _dualreach(space)

function _dualreach(space::DGGSpace)
    r = 0.0
    for cap in space.caps
        r = max(r, Float64(cap.radius))
    end
    return min(Float64(pi), 2 * r)
end

# Stored-cell caps give a useful reach bound for the single-chunk mixed grid.
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

# Short duplicate-free rings favor a linear scan over a lookup table.
@inline function _inring(ring, c)
    for x in ring
        x == c && return true
    end
    return false
end

@inline _cyclicnext(k::Int, n::Int) = k == n ? 1 : k + 1

# A bit mask keeps ring positions in a register.
@inline _isset(bits::UInt16, k::Int) = (bits >> (k - 1)) & one(UInt16) == one(UInt16)
@inline _set(bits::UInt16, k::Int) = bits | (one(UInt16) << (k - 1))

# Callers guarantee at least one set bit, so the cyclic scan always lands.
@inline function _nextset(bits::UInt16, k::Int, n::Int)
    j = k
    for _ in 1:n
        j = _cyclicnext(j, n)
        _isset(bits, j) && return j
    end
    return k
end

# Share location work while preserving the status needed by `weightsat!`.
function _locatedual(s::_DualSampler, p)
    grid = s.space.grid
    topo = s.state.grid
    host = localindex(grid, p)
    if host === nothing
        # A missing host site makes this a collection rim.
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

    # Edge neighbors delimit vertex fans within the ordered vertex ring.
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
            # Every weighted node needs a local source index.
            i === nothing && return (NO_GRID_DUALCELL, GR.WeightsRim)
            indices = SmallCollections.push(indices, i)
            k == b && break
            k = _cyclicnext(k, m)
        end
        return (GR.DualCell(indices, nodes, GR.MeanValue), GR.WeightsMapped)
    end
    return (NO_GRID_DUALCELL, unusable ? GR.WeightsDegenerate : GR.WeightsOutside)
end

# --- mixed levels -----------------------------------------------------------

const _MixedSampler{V} = GR.Sampler{<:DGGSpace,V,<:MultiOrderDualTopology}

@inline function _storedsite(mov::MultiOrderVector, k::Int)
    c = @inbounds mov.cells[k]
    return cell_centroid(levelgrid(mov.system, level(c)), c)
end

# Descendant intervals make the nearest stored refinement a bounded indexed scan.
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

Return the origin used by every query-centered `DGGSpace` chart.
"""
GR.chartat(::_MixedSampler, p) = (0.0, 0.0)

"""
    GlobalRegridding.dualcellat(s::Sampler{<:DGGSpace,V,<:MultiOrderDualTopology}, p)

Return the `MeanValue` triangle of stored sample sites containing `p`, or an
empty cell. The triangle uses the host and the two resolved ring sites that
bracket the query bearing.
"""
GR.dualcellat(s::_MixedSampler, p) = first(_locatemixeddual(s, p))

"""
    GlobalRegridding.weightsat!(row, s::Sampler{<:DGGSpace,V,<:MultiOrderDualTopology}, p)

Weight the mixed-level triangle containing `p`. The function clears the row on
entry and leaves it empty unless it returns `WeightsMapped`. See
[Dual-cell weight statuses](@ref).
"""
function GR.weightsat!(row::GR.WeightRow, s::_MixedSampler, p)
    empty!(row)
    cell, status = _locatemixeddual(s, p)
    GR.nodecount(cell) == 0 && return status
    return GR.dualweights!(row, cell, (0.0, 0.0))
end

function _locatemixeddual(s::_MixedSampler, p)
    mov = s.state.cells
    sys = mov.system
    host = localindex(mov, p)
    # A container hole supplies no host site for a stencil.
    host === nothing && return (NO_GRID_DUALCELL, GR.WeightsOutside)
    hostid = @inbounds mov.cells[host]
    hostgrid = levelgrid(sys, level(hostid))
    chart = GR.TangentChart(p)
    hostsite = GR.chartpoint(chart, cell_centroid(hostgrid, hostid))
    hostsite === nothing && return (NO_GRID_DUALCELL, GR.WeightsDegenerate)

    # The host's complete level supplies a conforming ring.
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
            # A placeholder centroid preserves the angular extent of a gap.
            GR.chartpoint(chart, cell_centroid(hostgrid, n))
        else
            # Ring members sharing a coarse ancestor contribute one site.
            (i == host || _inring(ids, i)) && continue
            GR.chartpoint(chart, _storedsite(mov, i))
        end
        q === nothing && return (NO_GRID_DUALCELL, GR.WeightsDegenerate)
        ids = SmallCollections.push(ids, i)
        sites = SmallCollections.push(sites, q)
    end
    n = length(ids)
    n >= 2 || return (NO_GRID_DUALCELL, GR.WeightsDegenerate)

    # The sites bracketing the host-to-origin bearing delimit the query wedge.
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
    # A gap cannot close a valid interpolation triangle.
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
