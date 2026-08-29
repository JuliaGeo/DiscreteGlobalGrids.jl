# DGGS implementation of the `GlobalRegridding.RegridSpace` interface.

import GlobalRegridding as GR
import DimensionalData as DD

# Automatic chunking targets this cell count; accuracy remains unchanged.
const DEFAULT_CHUNK_CELLS = 4096

# The space

"""
    DGGSpace(grid::AbstractGrid; chunklevel = nothing, chunkcells = $DEFAULT_CHUNK_CELLS)

Wrap `grid` as a [`GlobalRegridding.RegridSpace`](@ref). Chunks use nonempty
ancestor subtrees at `chunklevel`; the default chooses roughly `chunkcells`
cells per chunk. Grids without sorted subtrees use one chunk. Construction
computes one covering cap per chunk.
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

# Exact ancestor windows require metadata that avoids a full cell scan.
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

# A rooted `PartialGrid` restricts candidate ancestors to its root subtree.
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

# Lazy centroid sites avoid materializing vast multi-root cell collections.

# `PartialGrid` lookup must return collection-local indices.
cellat(space::DGGSpace, p::GO.UnitSphericalPoint) = localindex(space.grid, p)

GR.celltree(space::DGGSpace) = treeify(space.grid)

GR.chunkextents(space::DGGSpace) = space.caps

# Reusing the grid hierarchy avoids a second chunk index and extent adapter.
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

# `BlockCursor` avoids scanning CopernicusDEM's flat fanout of level-0 roots.
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

Return the cell tree restricted to `inds`, with leaves addressed by local index.
Selection proceeds through these routes:

 1. reuse the cached whole-space tree;
 2. use a grid-native [`subcursor`](@ref);
 3. reuse the hierarchy for an exact chunk range;
 4. build the packed cell-space fallback.

Analytical cell caps keep hierarchy extents lazy; other small trees cache leaf
caps.
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

Resolve a regridding target into [`DGGSpace`](@ref):

  - grids stand for their own cells;
  - [`CellLookup`](@ref), [`CellVector`](@ref), and
    [`MultiOrderCellSet`](@ref) produce a [`PartialGrid`](@ref);
  - mixed-level targets expand to reference-level leaves;
  - a bare system used as a destination selects the level closest to the source
    cell size.

A bare-system source must specify its level with [`levelgrid`](@ref).
"""
GR._asspace(grid::AbstractGrid, name::AbstractString) = DGGSpace(grid)

GR._asspace(lk::AbstractCellLookup, name::AbstractString) = DGGSpace(PartialGrid(lk))

GR._asspace(cv::AbstractCellVector, name::AbstractString) = DGGSpace(PartialGrid(cv))

GR._asspace(set::MultiOrderCellSet, name::AbstractString) =
    DGGSpace(PartialGrid(CellVector(set)))

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

# Mixed-level source presentation

"""
    GlobalRegridding.sourcespacefor(mov::MultiOrderVector, method)

Return the source space through which `method` reads a mixed-level container.

  - `Points()` sampling uses one [`MultiOrderGrid`](@ref) cell per stored cell.
    Covering-ancestor lookup preserves nearest-cell values with fewer plan
    columns than reference-level expansion.
  - `Intervals` sampling uses the reference-level expansion. Descendant leaves
    provide the gap-free cover required on non-congruent hierarchies such as H3
    and IGeo7, whose child polygons may leave gaps relative to the parent.
  - A container with one stored cell per reference-level leaf uses the
    [`PartialGrid`](@ref) path for either sampling.
"""
GR.sourcespacefor(mov::MultiOrderVector, method) = _readsstored(mov, method) ?
    DGGSpace(MultiOrderGrid(mov)) : GR._asspace(mov, "from")

GR.sourcespacefor(lk::MultiOrderLookup, method) = GR.sourcespacefor(parent(lk), method)

_readsstored(mov::MultiOrderVector, method) =
    _readsstored(mov, method, GR.sourcesampling(method))

_readsstored(mov::MultiOrderVector, method, ::DD.Lookups.Points) = _expandsleaves(mov)

_readsstored(::MultiOrderVector, method, ::DD.Lookups.Intervals) = false

@noinline _readsstored(mov::MultiOrderVector, method, sampling) = throw(ArgumentError(
    "$(nameof(typeof(method))) reports unsupported source sampling $(sampling) " *
    "for a mixed-level container. Define `sourcesampling` as `Points()` to " *
    "read its $(length(mov)) stored sample sites, or `Intervals(Center())` to " *
    "read its reference-level polygon cover at level $(reference_level(mov))."))

_leafcount(mov::MultiOrderVector) = isempty(mov) ? 0 : last(mov.offsets)

_expandsleaves(mov::MultiOrderVector) = _leafcount(mov) != length(mov)

"""
    GlobalRegridding.dimsource(lk::AbstractCellLookup)

Return the exact source target named by a cell lookup.

This is usually [`cellset`](@ref). A lookup produced by expanding a
[`MultiOrderVector`](@ref) names itself because the container has
method-specific geometry that differs from the lookup's leaf cells.
"""
GR.dimsource(lk::AbstractCellLookup) = _axissource(lk, cellset(lk))

_axissource(::AbstractCellLookup, set) = set
_axissource(lk::AbstractCellLookup, ::MultiOrderVector) = lk

GR.dimsource(lk::MultiOrderLookup) = cellset(lk)

"""
    GlobalRegridding.sourceview(lk::MultiOrderLookup, A, method)

Present a mixed-level cube as the array `method` reads, matching whichever
space [`GlobalRegridding.sourcespacefor`](@ref) resolves for the same `method`.
The cube therefore needs neither `from` nor a manual [`expand`](@ref).

  - Point sampling returns `A`, with one value per stored cell, against
    `DGGSpace(MultiOrderGrid(mov))`. This route supports pass-through dimensions.
  - Area sampling returns `expand(A, ref)` against the same reference-level
    cells resolved by [`GlobalRegridding._asspace`](@ref). Both enumerate each
    stored cell's descendants in stored order, preserving leaf alignment. The
    view remains lazy over the stored values.
  - A genuine expansion requires
    [`GlobalRegridding.refinementinvariant`](@ref). Replicated leaf values change
    methods that interpolate by sample-site position.
  - Area sampling requires one-dimensional cubes because [`expand`](@ref) is
    one-dimensional.
"""
function GR.sourceview(lk::MultiOrderLookup, A::DD.AbstractDimArray, method)
    mov = parent(lk)
    _readsstored(mov, method) && return A
    ndims(A) == 1 || _nomultidim(lk, method)
    (GR.refinementinvariant(method) || !_expandsleaves(mov)) ||
        _nointerpolation(lk, method)
    return expand(A, reference_level(lk))
end

GR.sourceview(::MultiOrderLookup, A, method) = nothing

@noinline _nointerpolation(lk::MultiOrderLookup, method) = throw(ArgumentError(
    "$(nameof(typeof(method))) changes when replicated values move from stored " *
    "cells to their leaf sites, so the mixed-level source cannot refine to " *
    "level $(reference_level(lk)) implicitly. Define " *
    "`GlobalRegridding.sourcesampling(method) = Points()` to read stored sample " *
    "sites, or request leaf-site interpolation explicitly with `regrid(expand(A, " *
    "DiscreteGlobalGrids.reference_level(lookup(A, Cells))); to = ..., " *
    "method = ...)`."))

@noinline _nomultidim(lk::MultiOrderLookup, method) = throw(ArgumentError(
    "$(nameof(typeof(method))) requires the reference-level area cover at " *
    "level $(reference_level(lk)), but `expand` is one-dimensional and the " *
    "source has pass-through dimensions. Regrid one slice at a time, or use " *
    "a sample-site method such as `NearestCell` or `DirectNearest`."))

"""
    GlobalRegridding.checksource(mov::MultiOrderVector, data, space)

Validate an explicit mixed-level `from` against the source value layout.

[`GlobalRegridding.sourcespacefor`](@ref) gives `from = mov` two presentations:

  - point methods use one value per stored cell;
  - area methods use one value per reference-level leaf.

A cube carrying [`MultiOrderLookup`](@ref) selects the matching presentation.
Its explicit `from` must name the same cells and reference level because the
axis determines the value ordering even when two containers have equal counts.
"""
function GR.checksource(mov::MultiOrderVector, data, space::GR.RegridSpace)
    data isa AbstractArray || return nothing
    own = _mixedaxis(data)
    if own !== nothing
        _samecontainer(mov, own) || _fromcontradicts(mov, own)
        return nothing
    end
    n = size(data, 1)
    (n == ncells(space) || !_expandsleaves(mov)) && return nothing
    stored, leaves = length(mov), _leafcount(mov)
    n == stored && throw(ArgumentError(
        "`from` resolves $stored mixed-level cells to $(ncells(space)) leaves " *
        "at reference level $(reference_level(mov)), but the source holds " *
        "$stored values, one per stored cell. Place them on " *
        "`DimArray(values, Cells(MultiOrderLookup(mov)))` and regrid with no " *
        "`from` at all; the axis defines how each value spreads over leaves."))
    n == leaves && throw(ArgumentError(
        "the source is already expanded to $leaves values at reference level " *
        "$(reference_level(mov)), while `from` selects the $stored stored " *
        "cells. Name the leaf geometry with `from = CellVector(mov)`, or let " *
        "the expanded cube's `Cells` axis name it."))
    return nothing
end

GR.checksource(lk::MultiOrderLookup, data, space::GR.RegridSpace) =
    GR.checksource(parent(lk), data, space)

function _mixedaxis(data)
    data isa DD.AbstractDimArray || return nothing
    for l in DD.lookup(data)
        l isa MultiOrderLookup && return parent(l)
    end
    return nothing
end

_samecontainer(a::MultiOrderVector, b::MultiOrderVector) =
    a === b || (reference_level(a) == reference_level(b) && a.cells == b.cells)

@noinline _fromcontradicts(mov::MultiOrderVector, own::MultiOrderVector) =
    throw(ArgumentError(
        "`from` names a different mixed-level container than the source's own " *
        "`Cells` axis: `from` holds $(length(mov)) cells at reference level " *
        "$(reference_level(mov)), the axis holds $(length(own)) at " *
        "$(reference_level(own)). The axis fixes the value ordering, so this " *
        "pairing would combine one container's geometry with another's values. " *
        "Drop `from`, or attach the values to the intended container's axis."))

# Labelling the output

"""
    GlobalRegridding.destinationdims(space::DGGSpace, sampling)

Return one [`Cells`](@ref) dimension containing a [`CellLookup`](@ref) in the
destination's local order. Cell labels are independent of `sampling`, and the
single-axis shape requires no output reshape.
"""
GR.destinationdims(space::DGGSpace, ::DD.Lookups.Sampling) =
    (Cells(CellLookup(space.grid)),)

"""
    GlobalRegridding.destinationdims(space::DGGSpace{<:MultiOrderGrid}, sampling)

Return a [`MultiOrderLookup`](@ref) over the container's stored cells.

This method labels output only. Regridding onto mixed levels requires a separate
area-normalization policy.
"""
GR.destinationdims(space::DGGSpace{<:MultiOrderGrid}, ::DD.Lookups.Sampling) =
    (Cells(MultiOrderLookup(cellset(space.grid))),)

# The dual cells a point method interpolates on.
include("dual_cells.jl")
