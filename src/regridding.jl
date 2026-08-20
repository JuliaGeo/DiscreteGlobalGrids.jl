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

# Return one position range per non-empty ancestor, or `nothing` when this
# cannot be determined without scanning every cell.
function _chunkwindows(grid::AbstractGrid, sys::AbstractHierarchicalGridSystem,
        lvl::Int, a::Int)
    complete = ncells(grid) == _levelcells(sys, lvl)
    (complete || grid isa PartialGrid) || return nothing
    ancestors = levelgrid(sys, a)
    ID = cellindextype(sys)
    ids = ID[]
    ranges = UnitRange{Int}[]
    for j in 1:ncells(ancestors)
        id = cellindex(ancestors, j)
        r = descendant_range(sys, id, lvl)
        w = complete ? (Int(first(r)):Int(last(r))) : _subsetwindow(grid, r)
        isempty(w) && continue
        push!(ids, id)
        push!(ranges, w)
    end
    return ids, ranges
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

cellindices(space::DGGSpace, chunk::Int) = space.ranges[chunk]

# Locate a position by binary-searching the sorted chunk starts.
function GR.chunkat(space::DGGSpace, i::Integer)
    p = Int(i)
    1 <= p <= ncells(space) || throw(BoundsError(space, p))
    k = searchsortedlast(space.starts, p)
    (1 <= k <= length(space.ranges) && p in space.ranges[k]) || throw(ArgumentError(
        "cell position $p belongs to no chunk of $space"))
    return k
end

GR.cellcentroid(space::DGGSpace, i::Int) =
    cell_centroid(space.grid, cellindex(space.grid, i))

function cellat(space::DGGSpace, p::GO.UnitSphericalPoint)
    c = cellat(space.grid, p)
    c === nothing && return nothing
    return cellposition(space.grid, c)
end

GR.celltree(space::DGGSpace) = treeify(space.grid)

GR.chunktree(space::DGGSpace) = DGGChunkTree(space)

GR.chunkextents(space::DGGSpace) = space.caps

# DGGS chunks are one-dimensional position ranges.
GR.chunkranges(space::DGGSpace, chunk::Integer, ::NTuple{1,Int}) =
    (space.ranges[Int(chunk)],)

"""
    GlobalRegridding.subtree(space::DGGSpace, inds)

Return the cell tree restricted to `inds`, preserving global cell positions.
The whole space gets a cursor with decoded ids and precomputed leaf caps;
exact chunk ranges reuse the grid hierarchy in `O(1)`; other ranges use a
bounding-cap tree.
"""
function GR.subtree(space::DGGSpace, inds::AbstractUnitRange{<:Integer})
    GR._iswholespace(space, inds) && return _cachedcelltree(space)
    cursor = _chunkcursor(space, inds)
    cursor === nothing || return cursor
    return GR.CellCapTree(space, inds)
end

# Reuse the hierarchy rooted at the ancestor for an exact chunk range.
function _chunkcursor(space::DGGSpace, inds::AbstractUnitRange{<:Integer})
    _ischunked(space) || return nothing
    isempty(inds) && return nothing
    k = searchsortedfirst(space.starts, Int(first(inds)))
    (k <= length(space.ranges) && space.ranges[k] == inds) || return nothing
    root = treeify(space.grid)
    root isa HierarchicalGridCursor && root.selection === nothing || return nothing
    return typeof(root)(space.grid, root.system, root.top_level, root.leaf_level,
        root.bucket_size, space.chunklevel, space.chunkids[k],
        Int(first(inds)), Int(last(inds)), nothing)
end

"""
    DGGChunkTree(space::DGGSpace)

A one-node spatial tree whose leaves are chunk numbers. Each stored ancestor
extent covers all cells in its chunk.
"""
struct DGGChunkTree{S<:DGGSpace,E}
    space::S
    extent::E
end

DGGChunkTree(space::DGGSpace) = DGGChunkTree(space,
    isempty(space.caps) ? Fallbacks.full_sphere_cap() :
    foldl(Fallbacks.merge_caps, space.caps))

Base.show(io::IO, t::DGGChunkTree) =
    print(io, "DGGChunkTree(", length(t.space.ranges), " chunks)")

STI.isspatialtree(::Type{<:DGGChunkTree}) = true
STI.node_extent_is_expensive(::Type{<:DGGChunkTree}) = false
STI.isleaf(::DGGChunkTree) = true
STI.nchild(::DGGChunkTree) = 0
STI.getchild(::DGGChunkTree) = ()
STI.node_extent(t::DGGChunkTree) = t.extent
STI.child_indices_extents(t::DGGChunkTree) =
    zip(eachindex(t.space.ranges), t.space.caps)

GOCore.best_manifold(t::DGGChunkTree) = GOCore.best_manifold(t.space.grid)
Trees.ncells(t::DGGChunkTree) = ncells(t.space.grid)
Trees.getcell(t::DGGChunkTree, i::Int) = getcell(t.space.grid, i)
Trees.getcell(t::DGGChunkTree) = (getcell(t.space.grid, i) for i in 1:ncells(t.space.grid))

# Resolving `to`

"""
    regridgrid(x) -> AbstractGrid

Return the grid represented by a regridding target. Lookup, vector, and
multi-order targets become a [`PartialGrid`](@ref).
"""
function regridgrid end

regridgrid(grid::AbstractGrid) = grid
regridgrid(lk::CellLookup) = PartialGrid(lk)
regridgrid(cv::CellVector) = PartialGrid(cv)
regridgrid(set::MultiOrderCellSet) = PartialGrid(CellVector(set))

const RegridTarget = Union{AbstractGrid,CellLookup,CellVector,MultiOrderCellSet}

GR._asspace(target::RegridTarget, name::AbstractString) = DGGSpace(regridgrid(target))

GR._asspace(sys::AbstractHierarchicalGridSystem, name::AbstractString) =
    throw(ArgumentError(
        "`$name = $(typeof(sys).name.name)()` names no cells until a level is " *
        "chosen. As a destination the level is matched to the source's cell " *
        "areas, but as a source you must name it with `levelgrid(sys, l)`."))

# Resolving `from`

# A `Cells` axis already names the cells a regrid would otherwise look for a
# raster lattice in, so a source given no `from` can point at the grid itself.
GR.dimsource(lk::CellLookup) = cellset(lk)

# A bare system as the destination takes the level closest to the source's cells.
GR._asspace(sys::AbstractHierarchicalGridSystem, name::AbstractString,
    src_space::GR.RegridSpace) = DGGSpace(levelgrid(sys, levelfor(sys, src_space)))

# API integration

# DGGS destination plans return a cell-indexed cube. Explicit leading bounds
# keep these aliases within the plan types' declared bounds.
const _DirectToDGG =
    GR.DirectPlan{<:GR.AbstractRegriddingMethod,<:GR.AbstractMissingPolicy,<:DGGSpace}
const _ChunkedToDGG =
    GR.ChunkedPlan{<:GR.AbstractRegriddingMethod,<:GR.AbstractMissingPolicy,<:DGGSpace}

function GR.regrid(data, plan::_DirectToDGG)
    return _ascube(invoke(GR.regrid, Tuple{Any,GR.DirectPlan}, data, plan), data, plan)
end

GR.regrid(data, plan::_ChunkedToDGG) =
    _ascube(GR.LazyRegridArray(data, plan), data, plan)

# Preserve non-spatial dimensions and label the destination with `Cells`.
function _ascube(out, data, plan::GR.AbstractRegriddingPlan)
    data isa DD.AbstractDimArray || return out
    lk = CellLookup(plan.dst_space.grid)
    ds = DD.dims(data)
    sd = GR.resolvespatialdims(data, Int(ncells(plan.src_space)))
    others = Tuple(ds[i] for i in eachindex(ds) if !(i in sd))
    raw = out isa DD.AbstractDimArray ? parent(out) : out
    return DD.DimArray(raw, (Cells(lk), others...))
end
