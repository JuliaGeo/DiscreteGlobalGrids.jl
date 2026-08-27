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

A mixed-level target — a [`MultiOrderVector`](@ref) or the axis that carries
one — is expanded to its reference level first, so the destination has one cell
per leaf. As a *source* the cube resolves itself and needs no `from` at all;
see `GlobalRegridding.sourceview`.

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

# The storage container and its axis name the same cells the query-side set
# does, so all three resolve alike: expanded to the reference level.
GR._asspace(mov::MultiOrderVector, name::AbstractString) =
    DGGSpace(PartialGrid(CellVector(mov)))

GR._asspace(lk::MultiOrderLookup, name::AbstractString) =
    GR._asspace(parent(lk), name)

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

# Only reached when the axis cannot present the cube itself — a mixed-level cube
# with more than one dimension, which `expand` is not defined on.
GR.dimsource(lk::MultiOrderLookup) = cellset(lk)

"""
    GlobalRegridding.sourceview(lk::MultiOrderLookup, A, method)

Present a mixed-level cube at its [`reference_level`](@ref): `expand(A, ref)`,
whose `Cells` axis is the [`CellLookup`](@ref) over the same cells
`GlobalRegridding._asspace(lk, "from")` resolves to. A regrid reads this view,
so a mixed-level cube is a source with no `from` and no manual [`expand`](@ref).

  - Alignment: both sides come from `CellVector(mov; level = ref)`, so leaf `k`
    of the view is leaf `k` of the space — the expansion enumerates each stored
    cell's `descendant_range` in stored order, and the space's cells are those
    same ranges merged where adjacent.
  - Lazy: one stored value per multi-order cell, whatever the leaf count.
  - Only for a `method` that is `GlobalRegridding.refinementinvariant`. Others
    are refused, because every leaf under a stored cell carries one replicated
    value and interpolating between those sites rebuilds the coarsening
    staircase at leaf spacing.
  - `expand` is one-dimensional, so a cube with pass-through dimensions is
    refused outright rather than left to fail on a count.
"""
function GR.sourceview(lk::MultiOrderLookup, A::DD.AbstractDimArray, method)
    ndims(A) == 1 || _nomultidim(lk, method)
    GR.refinementinvariant(method) || _nointerpolation(lk, method)
    return expand(A, reference_level(lk))
end

GR.sourceview(::MultiOrderLookup, A, method) = nothing

@noinline _nomultidim(lk::MultiOrderLookup, method) = throw(ArgumentError(
    "a mixed-level cube presents itself refined to level " *
    "$(reference_level(lk)), and `expand` is one-dimensional, so it cannot do " *
    "that for a cube with pass-through dimensions. Regrid one slice at a time."))

@noinline _nointerpolation(lk::MultiOrderLookup, method) = throw(ArgumentError(
    "$(nameof(typeof(method))) interpolates between source sample sites, and a " *
    "mixed-level cube can only present itself refined to level " *
    "$(reference_level(lk)), where every leaf under a stored cell repeats that " *
    "cell's one value — interpolating between them rebuilds the coarsening " *
    "steps at leaf spacing. Use an area or nearest-cell method, or wait for the " *
    "native mixed-level source. To interpolate on the leaves anyway, say so: " *
    "`regrid(expand(A, DiscreteGlobalGrids.reference_level(lookup(A, Cells))); " *
    "to = ..., from = cellset(lookup(A, Cells)), method = ...)`."))

"""
    GlobalRegridding.checksource(mov::MultiOrderVector, data, space)

Refuse a `from` naming mixed-level cells against values stored one per *cell*
rather than one per leaf. The space `mov` resolves to is its reference-level
expansion, so the two counts differ whenever the container is genuinely mixed,
and the flatten step would otherwise report only a size mismatch.

A cube carrying the [`MultiOrderLookup`](@ref) passes: it expands itself.
"""
function GR.checksource(mov::MultiOrderVector, data, space::GR.RegridSpace)
    n = length(mov)
    (n == ncells(space) || _carriesmixed(data) || size(data, 1) != n) && return nothing
    throw(ArgumentError(
        "`from` names $n mixed-level cells, which stand for $(ncells(space)) " *
        "cells at reference level $(reference_level(mov)), but the source " *
        "holds $n values — one per stored cell. No `from` can pair those: the " *
        "axis itself says how the stored values spread over the cells. Put " *
        "them on it — `DimArray(values, Cells(MultiOrderLookup(mov)))` — and " *
        "regrid with no `from` at all."))
end

GR.checksource(lk::MultiOrderLookup, data, space::GR.RegridSpace) =
    GR.checksource(parent(lk), data, space)

_carriesmixed(data) = data isa DD.AbstractDimArray &&
    any(l -> l isa MultiOrderLookup, DD.lookup(data))

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
