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

# Routing a mixed-level source

"""
    GlobalRegridding.sourcespacefor(mov::MultiOrderVector, method)

The source space `method` reads a mixed-level container through.

  - A method that reads source **sample sites** (`GlobalRegridding.sourcesampling`
    is `Points()`) takes the stored cells as they are: `DGGSpace(MultiOrderGrid(mov))`,
    one cell per stored cell. A destination point resolves through the
    container's covering-ancestor lookup, which is the same verdict the
    reference-level expansion reaches, so nearest-cell answers are unchanged and
    the plan no longer carries a column per leaf.
  - A method that reads source **area** (`Intervals`) keeps the expansion. A
    stored cell's descendant leaves are the only gap-free cover of it on a
    non-congruent hierarchy — H3 and IGeo7, where a parent's polygon is not the
    union of its children's — so coarse polygons would leave slivers and lose
    mass.
  - A container that stores one cell per reference-level leaf expands to itself,
    so it keeps the [`PartialGrid`](@ref) path either way: same cells, same
    order, same count, and already-tested code.
"""
GR.sourcespacefor(mov::MultiOrderVector, method) = _readsstored(mov, method) ?
    DGGSpace(MultiOrderGrid(mov)) : GR._asspace(mov, "from")

GR.sourcespacefor(lk::MultiOrderLookup, method) = GR.sourcespacefor(parent(lk), method)

# The one routing decision, asked in both places it is needed: `true` means the
# stored cells ARE the source — `MultiOrderGrid` for the geometry, the cube as
# it stands for the values.
_readsstored(mov::MultiOrderVector, method) =
    _readsstored(mov, method, GR.sourcesampling(method))

_readsstored(mov::MultiOrderVector, method, ::DD.Lookups.Points) = _expandsleaves(mov)

_readsstored(::MultiOrderVector, method, ::DD.Lookups.Intervals) = false

@noinline _readsstored(mov::MultiOrderVector, method, sampling) = throw(ArgumentError(
    "$(nameof(typeof(method))) reports `sourcesampling` $(sampling), which is " *
    "neither `Points()` nor `Intervals()`, so a mixed-level container cannot " *
    "tell which of its two presentations to offer: the stored cells " *
    "themselves, which give sample sites but no gap-free polygon cover, or " *
    "the expansion to reference level $(reference_level(mov)), which gives a " *
    "gap-free cover of one cell per leaf. Declare one of the two samplings."))

# The leaves `mov` presents at its reference level, accumulated at construction.
_leafcount(mov::MultiOrderVector) = isempty(mov) ? 0 : last(mov.offsets)

# Whether presenting `mov` at its reference level names more cells than it
# stores. Exactly the question the routing turns on, and it costs one comparison.
_expandsleaves(mov::MultiOrderVector) = _leafcount(mov) != length(mov)

"""
    GlobalRegridding.dimsource(lk::AbstractCellLookup)

The cells the axis holds, as a `from` target — so a cube with a `Cells` axis is
a source with no `from` at all, the way a raster is.

Usually [`cellset`](@ref), what the collection was built *from*. The exception
is a collection built by expanding a [`MultiOrderVector`](@ref): `cellset` names
the container there, and a container resolves to a different space per method
(`GlobalRegridding.sourcespacefor`), so it is not a name for these cells. The
axis names itself instead, which is exact in every case.
"""
GR.dimsource(lk::AbstractCellLookup) = _axissource(lk, cellset(lk))

_axissource(::AbstractCellLookup, set) = set
_axissource(lk::AbstractCellLookup, ::MultiOrderVector) = lk

# The container the axis's presented view is written against, whichever
# presentation `sourceview` chose: the stored cells natively, or their
# reference-level expansion. `sourcespacefor` reads the method and resolves it
# to the matching space.
GR.dimsource(lk::MultiOrderLookup) = cellset(lk)

"""
    GlobalRegridding.sourceview(lk::MultiOrderLookup, A, method)

Present a mixed-level cube as the array `method` reads, matching whichever
space [`GlobalRegridding.sourcespacefor`](@ref) resolves for the same `method`.
Either way the cube is a source with no `from` and no manual [`expand`](@ref).

  - A method that reads sample sites takes `A` **as it stands**, one value per
    stored cell, against `DGGSpace(MultiOrderGrid(mov))`. Nothing is expanded,
    so the cube's pass-through dimensions are no obstacle either.
  - A method that reads area takes `expand(A, ref)`, on the same cells
    `GlobalRegridding._asspace(lk, "from")` resolves to. Alignment: both sides
    come from `CellVector(mov; level = ref)`, so leaf `k` of the view is leaf
    `k` of the space — the expansion enumerates each stored cell's
    `descendant_range` in stored order, and the space's cells are those same
    ranges merged where adjacent. It stays lazy: one stored value per
    multi-order cell, whatever the leaf count.
  - The expansion is offered only to a `method` that is
    `GlobalRegridding.refinementinvariant`, unless the container stores one cell
    per leaf and the expansion is therefore the identity. Otherwise every leaf
    under a stored cell carries one replicated value, and interpolating between
    those sites rebuilds the coarsening staircase at leaf spacing.
  - `expand` is one-dimensional, so a cube with pass-through dimensions can only
    take the native route; an area method over one is refused outright.
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
    "$(nameof(typeof(method))) interpolates between source sample sites, and " *
    "reports `sourcesampling` $(GR.sourcesampling(method)) — area, not points " *
    "— so a mixed-level cube can only present itself refined to level " *
    "$(reference_level(lk)), where every leaf under a stored cell repeats that " *
    "cell's one value; interpolating between them rebuilds the coarsening " *
    "steps at leaf spacing. Declare `GlobalRegridding.sourcesampling(method) = " *
    "Points()` to read the stored cells natively, or interpolate on the leaves " *
    "anyway by expanding first, which says so: `regrid(expand(A, " *
    "DiscreteGlobalGrids.reference_level(lookup(A, Cells))); to = ..., " *
    "method = ...)`."))

@noinline _nomultidim(lk::MultiOrderLookup, method) = throw(ArgumentError(
    "$(nameof(typeof(method))) reads source cell area, so a mixed-level cube " *
    "must present itself refined to level $(reference_level(lk)) — and " *
    "`expand` is one-dimensional, so it cannot do that for a cube with " *
    "pass-through dimensions. Regrid one slice at a time, or use a method that " *
    "reads sample sites (`NearestCell`, `DirectNearest`), which takes the " *
    "stored cells as they are and needs no expansion at all."))

"""
    GlobalRegridding.checksource(mov::MultiOrderVector, data, space)

Refuse a `from` naming mixed-level cells against values the space it resolved
to cannot be laid out against, and say which presentation the caller is holding.

`from = mov` names a container with two presentations, and the method picks
between them (`GlobalRegridding.sourcespacefor`). So the same spelling stands
for the stored cells under a point method and for their reference-level
expansion under an area one, and either can meet the wrong array:

  - values one per **stored cell** against the expansion — no `from` can pair
    those, because only the axis says how one stored value spreads over its
    leaves;
  - values one per **leaf** against the stored cells — the data is already
    expanded, so it should name the expansion rather than the container.

A cube carrying the [`MultiOrderLookup`](@ref) passes either way: it presents
itself to match.
"""
function GR.checksource(mov::MultiOrderVector, data, space::GR.RegridSpace)
    (data isa AbstractArray && !_carriesmixed(data)) || return nothing
    n = size(data, 1)
    (n == ncells(space) || !_expandsleaves(mov)) && return nothing
    stored, leaves = length(mov), _leafcount(mov)
    n == stored && throw(ArgumentError(
        "`from` names $stored mixed-level cells, which stand for " *
        "$(ncells(space)) cells at reference level $(reference_level(mov)) for " *
        "this method, but the source holds $stored values — one per stored " *
        "cell. No `from` can pair those: the axis itself says how the stored " *
        "values spread over the cells. Put them on it — " *
        "`DimArray(values, Cells(MultiOrderLookup(mov)))` — and regrid with no " *
        "`from` at all."))
    n == leaves && throw(ArgumentError(
        "`from` names $stored mixed-level cells, which this method reads as " *
        "the $(ncells(space)) cells stored, but the source holds $leaves " *
        "values — one per leaf at reference level $(reference_level(mov)). " *
        "The data is already expanded, so name the expansion rather than the " *
        "container: `from = CellVector(mov)`, or drop `from` and let the " *
        "expanded cube's own `Cells` axis name it."))
    return nothing
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

"""
    GlobalRegridding.destinationdims(space::DGGSpace{<:MultiOrderGrid}, sampling)

The [`MultiOrderLookup`](@ref) over the container's stored cells.

Mixed levels have no single-level [`CellLookup`](@ref) to be labelled with, so
a result written over these cells says so with the axis that carries them. This
is labelling only: routing a regrid *onto* mixed levels is a separate question
about area normalisation and does not ship here.
"""
GR.destinationdims(space::DGGSpace{<:MultiOrderGrid}, ::DD.Lookups.Sampling) =
    (Cells(MultiOrderLookup(cellset(space.grid))),)

# The dual cells a point method interpolates on.
include("dual_cells.jl")
