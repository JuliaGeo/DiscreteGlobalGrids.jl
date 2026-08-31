# DGGS implementation of the `GlobalRegridding.RegridSpace` interface.

import GlobalRegridding as GR
import DimensionalData as DD

# Target cell count for automatic chunking. This affects memory use, not accuracy.
const DEFAULT_CHUNK_CELLS = 4096

# The space

"""
    DGGSpace(grid::AbstractGrid; chunklevel = nothing, chunkcells = $DEFAULT_CHUNK_CELLS)

Wrap `grid` as a `GlobalRegridding.RegridSpace`.

Chunks are non-empty ancestor subtrees at `chunklevel`. By default, the level
is chosen to keep roughly `chunkcells` cells per chunk. Grids without sorted
subtrees use one chunk. Construction computes only one covering cap per chunk.
"""
struct DGGSpace{G<:AbstractGrid,ID,C} <: GR.RegridSpace
    grid::G
    chunklevel::Int              # `< 0` means the single-chunk fallback
    chunkids::Vector{ID}         # empty in the single-chunk fallback
    ranges::Vector{UnitRange{Int}}
    starts::Vector{Int}          # `first.(ranges)`, for locating a chunk by its cells
    caps::Vector{C}
end

function DGGSpace(grid::AbstractGrid; chunklevel::Union{Nothing,Integer}=nothing,
        chunkcells::Integer=DEFAULT_CHUNK_CELLS)
    chunkcells > 0 || throw(ArgumentError("chunkcells must be positive, got $chunkcells"))
    sys = system(grid)
    lvl = level(grid)
    n = ncells(grid)
    (sys === nothing || lvl === nothing || n == 0 || !has_sorted_subtrees(sys)) &&
        return _wholechunk(grid)
    a = chunklevel === nothing ? _chunklevel(sys, lvl, n, Int(chunkcells)) : Int(chunklevel)
    (first(levels(sys)) <= a <= lvl) || throw(ArgumentError(
        "chunklevel $a is outside $(first(levels(sys))):$lvl, the levels between " *
        "the system's root and the grid's own level"))
    windows = _chunkwindows(grid, sys, lvl, a)
    windows === nothing && return _wholechunk(grid)
    ids, ranges = windows
    return DGGSpace(grid, a, ids, ranges,
        [first(r) for r in ranges], [node_extent(sys, id) for id in ids])
end

DGGSpace(space::DGGSpace) = space

function _wholechunk(grid::AbstractGrid)
    sys = system(grid)
    ID = sys === nothing ? Any : cellindextype(sys)
    n = ncells(grid)
    return DGGSpace(grid, -1, ID[], [1:n], [1],
        [Fallbacks.full_sphere_cap()])
end

_ischunked(space::DGGSpace) = space.chunklevel >= 0

Base.show(io::IO, space::DGGSpace) =
    print(io, "DGGSpace(", ncells(space.grid), " cells, ", length(space.ranges),
        _ischunked(space) ? " chunks at level $(space.chunklevel))" : " chunk)")

# Use the level grid because some systems provide their cell count there.
_levelcells(sys::AbstractHierarchicalGridSystem, l::Integer) = ncells(levelgrid(sys, l))

# Choose the ancestor level closest to `target` cells per chunk.
function _chunklevel(sys::AbstractHierarchicalGridSystem, lvl::Int, n::Int, target::Int)
    best, bestscore = lvl, Inf
    for a in first(levels(sys)):lvl
        score = abs(log(n / _levelcells(sys, a)) - log(target))
        score < bestscore && ((best, bestscore) = (a, score))
    end
    return best
end

# Return one index range per non-empty ancestor, or `nothing` when this
# cannot be determined without scanning every cell.
function _chunkwindows(grid::AbstractGrid, sys::AbstractHierarchicalGridSystem,
        lvl::Int, a::Int)
    complete = ncells(grid) == _levelcells(sys, lvl)
    (complete || grid isa PartialGrid) || return nothing
    ancestors = levelgrid(sys, a)
    ID = cellindextype(sys)
    ids = ID[]
    ranges = UnitRange{Int}[]
    for j in _ancestorindices(grid, sys, a, ancestors)
        id = cellindex(ancestors, j)
        r = descendant_range(sys, id, lvl)
        w = complete ? (Int(first(r)):Int(last(r))) : _subsetwindow(grid, r)
        isempty(w) && continue
        push!(ids, id)
        push!(ranges, w)
    end
    return ids, ranges
end

# The level-`a` ancestors that can hold any of `grid`'s cells. Every one of them
# in general — the scan is what decides which are non-empty — but a ROOTED
# `PartialGrid` holds nothing outside its root's subtree, so only that root's
# level-`a` descendants can qualify and the rest of the level need never be
# visited. That is exact, not a heuristic: the constructor checks the ancestry.
_ancestorindices(::AbstractGrid, ::AbstractHierarchicalGridSystem, ::Int,
    ancestors::AbstractGrid) = 1:ncells(ancestors)

function _ancestorindices(grid::PartialGrid, sys::AbstractHierarchicalGridSystem,
        a::Int, ancestors::AbstractGrid)
    grid.root_level >= first(levels(sys)) || return 1:ncells(ancestors)
    if grid.root_level <= a
        r = descendant_range(sys, grid.root_id, a)
        return Int(first(r)):Int(last(r))
    end
    # A root deeper than the chunk level puts the whole grid under one ancestor.
    p = globalindex(ancestors, ancestor(sys, grid.root_id, a))
    return p:p
end

# Sorted subset IDs make each ancestor's descendants a contiguous interval.
function _subsetwindow(grid::PartialGrid, r::AbstractUnitRange)
    ids = grid.ids
    isempty(ids) && return 1:0
    lo = searchsortedfirst(ids, cellindex(grid.complete, Int(first(r))))
    hi = searchsortedlast(ids, cellindex(grid.complete, Int(last(r))))
    return Int(lo):Int(hi)
end

# `RegridSpace` interface

ncells(space::DGGSpace) = ncells(space.grid)
getcell(space::DGGSpace, i::Int) = getcell(space.grid, i)

GOCore.manifold(space::DGGSpace) = GOCore.best_manifold(space.grid)

GR.nchunks(space::DGGSpace) = length(space.ranges)

GR.ownedindices(space::DGGSpace, chunk::Int) = space.ranges[chunk]

# Locate an index by binary-searching the sorted chunk starts.
function GR.chunkat(space::DGGSpace, i::Integer)
    p = Int(i)
    1 <= p <= ncells(space) || throw(BoundsError(space, p))
    k = searchsortedlast(space.starts, p)
    (1 <= k <= length(space.ranges) && p in space.ranges[k]) || throw(ArgumentError(
        "cell index $p belongs to no chunk of $space"))
    return k
end

GR.cellcentroid(space::DGGSpace, i::Int) =
    cell_centroid(space.grid, cellindex(space.grid, i))

"""
    GlobalRegridding.cellneighbors(space::DGGSpace, i::Int)

The grid's own edge one-ring of cell `i`, clipped to the collection: what a
gradient-recovering method fits a cell's neighbours from.

`Edge()` connectivity gives every system a stencil a gradient fit can use:

  - A hexagonal system answers the same six cells either way, five at a
    pentagon.
  - A quadrilateral one answers its four edge neighbours, a well-conditioned
    two-parameter fit already, where `Vertex()` would answer eight.
  - CopernicusDEM's vertex ring at a pole apex is the whole polar row, which a
    gradient fit has no use for.

On a `PartialGrid` the ring is the complete level's, clipped to members, so a
rim cell simply has fewer neighbours and its fit is one-sided.
"""
GR.cellneighbors(space::DGGSpace, i::Int) =
    neighbors(space.grid, i, 1; connectivity = Edge())

"""
    GlobalRegridding.celldiameter(space::DGGSpace) -> Float64

Twice the widest chunk cover, in radians: the same bound
`supportradius(::BarycentricPoint, ::DGGSpace)` rests on.

Every cell lies inside its own chunk's cover, so twice the widest of those
covers bounds any cell's diameter.

The bound is loose by the ratio of chunk to cell width, which costs discovery
work, never weights. The single-chunk fallback's full-sphere cover makes it
`π`, where the one source chunk leaves nothing to discover.
"""
GR.celldiameter(space::DGGSpace) = _dualreach(space)

# The sites a point method interpolates between are `cellcentroid` above, read
# through `GlobalRegridding`'s own `CentroidSites`: a pure vector that computes
# an entry on read and holds nothing, so preparing a sampler materialises
# nothing and concurrent queries share it. The collection's centroid field would
# read the same, but `CellVector` indexes every cell of a holding whose ids are
# not one subtree — 2.5e10 of them for the GLO-90 source, 187 GiB — and a
# sampler never reads more sites than its stencils name.

# Local, not global: the inverse of `cellcentroid` above, and on a `PartialGrid`
# a different numbering from the grid's own cell ids.
cellat(space::DGGSpace, p::GO.UnitSphericalPoint) = localindex(space.grid, p)

GR.celltree(space::DGGSpace) = treeify(space.grid)

GR.chunkextents(space::DGGSpace) = space.caps

# The DGG space is its own private chunk index: candidate queries descend the
# grid's existing hierarchy to `chunklevel`, with no second tree or extent
# adapter. Construct the hierarchical cursor directly because some systems use
# a different optimized tree for cell intersections.
GR.chunkindex(space::DGGSpace) = space

function GR.candidatechunks!(out::Vector{Int}, space::DGGSpace, dstcap::GR.Cap;
        radius::Real = 0.0)
    empty!(out)
    intersects = GR.DilatedIntersects(Float64(radius))
    if GR.nchunks(space) == 1
        intersects(dstcap, only(space.caps)) && push!(out, 1)
        return out
    end
    if space.chunklevel == 0 && system(space.grid) isa CopernicusDEM.CopernicusDEMSystem
        frontier = levelgrid(system(space.grid), space.chunklevel)
        _mappedfrontierchunks!(out, space, treeify(frontier), frontier, dstcap, intersects)
        sort!(out)
        unique!(out)
        return out
    end
    root = HierarchicalGridCursor(space.grid; bucket_size = 0)
    _dggcandidatechunks!(out, space, root, dstcap, intersects)
    sort!(out)
    unique!(out)
    return out
end

# CopernicusDEM has tens of thousands of level-0 roots, so the generic
# ancestor cursor would still scan a flat root fanout for every query. Its
# existing BlockCursor is the spatial hierarchy over that lattice. Traverse a
# complete grid at the requested frontier level and filter its leaf ids through
# the selected `PartialGrid` chunk ids; no source pixels are materialized.
function _mappedfrontierchunks!(out::Vector{Int}, space::DGGSpace, node, frontier,
        dstcap, intersects)
    intersects(dstcap, STI.node_extent(node)) || return out
    if STI.isleaf(node)
        for (index, cap) in STI.child_indices_extents(node)
            intersects(dstcap, cap) || continue
            id = cellindex(frontier, index)
            k = searchsortedfirst(space.chunkids, id)
            k <= length(space.chunkids) && space.chunkids[k] == id && push!(out, k)
        end
        return out
    end
    for child in STI.getchild(node)
        _mappedfrontierchunks!(out, space, child, frontier, dstcap, intersects)
    end
    return out
end

function _dggcandidatechunks!(out::Vector{Int}, space::DGGSpace,
        node::HierarchicalGridCursor, dstcap, intersects)
    intersects(dstcap, STI.node_extent(node)) || return out
    if node.level >= space.chunklevel
        inds = Engine.node_indices(node)
        isempty(inds) || push!(out, GR.chunkat(space, first(inds)))
        return out
    end
    for child in STI.getchild(node)
        _dggcandidatechunks!(out, space, child, dstcap, intersects)
    end
    return out
end

# DGGS chunks are one-dimensional index ranges.
GR.chunkranges(space::DGGSpace, chunk::Integer, ::NTuple{1,Int}) =
    (space.ranges[Int(chunk)],)

"""
    GlobalRegridding.subtree(space::DGGSpace, inds)

Return the cell tree restricted to `inds`, with leaves still addressed by the
space's local index.
In order: the whole space, a grid that can window its own tree
([`subcursor`](@ref)), an exact chunk range (the grid hierarchy in `O(1)`), and
otherwise the common packed cell-space fallback. Grids with an analytical cell
cap keep the hierarchical tree lazy; other whole spaces and small enough chunks
carry precomputed leaf caps.
"""
function GR.subtree(space::DGGSpace, inds::AbstractUnitRange{<:Integer})
    GR._iswholespace(space, inds) && return _cachedcelltree(space)
    # Dispatch-only, so it costs nothing for the grids that have no method.
    window = subcursor(space.grid, inds)
    window === nothing || return window
    cursor = _chunkcursor(space, inds)
    cursor === nothing || return _cachedchunktree(cursor, inds)
    return GR.CellSpaceRTree(space, inds)
end

# Reuse the hierarchy rooted at the ancestor for an exact chunk range.
function _chunkcursor(space::DGGSpace, inds::AbstractUnitRange{<:Integer})
    _ischunked(space) || return nothing
    isempty(inds) && return nothing
    k = searchsortedfirst(space.starts, Int(first(inds)))
    (k <= length(space.ranges) && space.ranges[k] == inds) || return nothing
    root = treeify(space.grid)
    root isa HierarchicalGridCursor && root.selection === nothing || return nothing
    complete_subtree = length(inds) ==
        length(descendant_range(root.system, space.chunkids[k], root.leaf_level))
    return typeof(root)(space.grid, root.system, root.top_level, root.leaf_level,
        root.bucket_size, space.chunklevel, space.chunkids[k],
        Int(first(inds)), Int(last(inds)), complete_subtree, nothing)
end

# Resolving `to` and `from`

"""
    GlobalRegridding._asspace(target, name)
    GlobalRegridding._asspace(target, name, src_space)

Return the [`DGGSpace`](@ref) over the cells a regridding target names. A grid
stands for itself; a [`CellLookup`](@ref), a [`CellVector`](@ref) and a
[`MultiOrderCellSet`](@ref) name the [`PartialGrid`](@ref) of their cells.

A bare system names no cells until a level is chosen. As a destination it takes
the level whose cells are closest in size to the source's, which is the only
spelling that reads `src_space`; as a source there is nothing to match against
and it is an error.
"""
GR._asspace(grid::AbstractGrid, name::AbstractString) = DGGSpace(grid)

GR._asspace(lk::AbstractCellLookup, name::AbstractString) = DGGSpace(PartialGrid(lk))

GR._asspace(cv::AbstractCellVector, name::AbstractString) = DGGSpace(PartialGrid(cv))

GR._asspace(set::MultiOrderCellSet, name::AbstractString) =
    DGGSpace(PartialGrid(CellVector(set)))

GR._asspace(sys::AbstractHierarchicalGridSystem, name::AbstractString) =
    throw(ArgumentError(
        "`$name = $(typeof(sys).name.name)()` names no cells until a level is " *
        "chosen. As a destination the level is matched to the source's cell " *
        "areas, but as a source you must name it with `levelgrid(sys, l)`."))

GR._asspace(sys::AbstractHierarchicalGridSystem, name::AbstractString,
    src_space::GR.RegridSpace) = DGGSpace(levelgrid(sys, levelfor(sys, src_space)))

# A `Cells` axis already names the cells a regrid would otherwise look for a
# raster lattice in, so a source given no `from` can point at the grid itself.
GR.dimsource(lk::AbstractCellLookup) = cellset(lk)

# Labelling the output

"""
    GlobalRegridding.destinationdims(space::DGGSpace, sampling)

Return the single [`Cells`](@ref) dimension a result over this space carries: a
[`CellLookup`](@ref) over the destination's own cells, in its local index order.

A cell holds one value however that value was measured, so the lookup is the
same whichever `sampling` the method asks for. Being the space's only axis, it
is the whole of the destination's shape, and a regrid needs no reshape to put a
result on it.
"""
GR.destinationdims(space::DGGSpace, ::DD.Lookups.Sampling) =
    (Cells(CellLookup(space.grid)),)

# The dual cells a point method interpolates on.
include("dual_cells.jl")
