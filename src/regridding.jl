# ---------------------------------------------------------------------------
# The regridding face of this package: a cell collection as a
# `GlobalRegridding.RegridSpace`, the target spellings that resolve into one,
# and the cube a result comes back as.
#
# `GlobalRegridding` has no dependency on this package and no notion of what a
# DGGS is. Everything specific to one — the cursor trees, the ancestor-level
# chunking, `to = ` sugar, the `Cells` axis on the output — lives here.
#
# `cellat` and `cellindices` are `GlobalRegridding`'s bindings, imported and
# extended here exactly as `ncells` and `getcell` are
# `ConservativeRegridding.Trees`', so a session holding both surfaces sees one
# function per name rather than an ambiguity.
# ---------------------------------------------------------------------------

import GlobalRegridding as GR
import DimensionalData as DD

const _Cap = GO.UnitSpherical.SphericalCap{Float64}

# Cells per chunk the default ancestor level aims for. Chunks are the unit the
# lazy path holds accumulators and weight blocks for, so this trades weight
# reuse against residency, not accuracy.
const DEFAULT_CHUNK_CELLS = 4096

# ===========================================================================
# The space
# ===========================================================================

"""
    DGGSpace(grid::AbstractGrid; chunklevel = nothing, chunkcells = $DEFAULT_CHUNK_CELLS)

`grid` as a `GlobalRegridding.RegridSpace`, usable on either side of a
[`regrid`](@ref).

Geometry is the grid's own: cells are [`cell_polygon`](@ref)s, the cell tree is
[`treeify(grid)`](@ref treeify), point location is [`cellat`](@ref), and the
manifold is the unit sphere every grid here computes on.

**Chunks are ancestor cells.** Each chunk is one cell of `levelgrid(system(grid),
chunklevel)` that the grid holds descendants of, and its cells are that cell's
[`descendant_range`](@ref) intersected with the grid — a `UnitRange` of
positions, because positions ascend with canonical ids. Empty ancestors are not
chunks. `chunklevel` defaults to the level whose subtrees hold about `chunkcells`
of the grid's cells each; a grid with no system, no level, or no
[`has_sorted_subtrees`](@ref) (A5) is one whole chunk.

Construction reads no cell geometry beyond one covering cap per chunk.
"""
struct DGGSpace{G<:AbstractGrid,ID} <: GR.RegridSpace
    grid::G
    chunklevel::Int              # `< 0` means the single-chunk fallback
    chunkids::Vector{ID}         # empty in the single-chunk fallback
    ranges::Vector{UnitRange{Int}}
    starts::Vector{Int}          # `first.(ranges)`, for locating a chunk by its cells
    caps::Vector{_Cap}
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
    return DGGSpace{typeof(grid),eltype(ids)}(grid, a, ids, ranges,
        [first(r) for r in ranges], [node_extent(sys, id) for id in ids])
end

DGGSpace(space::DGGSpace) = space

function _wholechunk(grid::AbstractGrid)
    sys = system(grid)
    ID = sys === nothing ? Any : cellindextype(sys)
    n = ncells(grid)
    return DGGSpace{typeof(grid),ID}(grid, -1, ID[], [1:n], [1],
        [Fallbacks.full_sphere_cap()])
end

_ischunked(space::DGGSpace) = space.chunklevel >= 0

Base.show(io::IO, space::DGGSpace) =
    print(io, "DGGSpace(", ncells(space.grid), " cells, ", length(space.ranges),
        _ischunked(space) ? " chunks at level $(space.chunklevel))" : " chunk)")

# The ancestor level whose subtrees hold about `target` of this grid's cells
# each. Density is the grid's own, so a sparse subset chunks coarser than the
# complete level it is drawn from rather than into slivers.
function _chunklevel(sys::AbstractHierarchicalGridSystem, lvl::Int, n::Int, target::Int)
    best, bestscore = lvl, Inf
    for a in first(levels(sys)):lvl
        score = abs(log(n / ncells(sys, a)) - log(target))
        score < bestscore && ((best, bestscore) = (a, score))
    end
    return best
end

# One window of grid positions per ancestor cell, empties dropped. `nothing`
# means the grid's positions cannot be placed in the system's without walking
# every cell, which is the single-chunk fallback's job to absorb.
function _chunkwindows(grid::AbstractGrid, sys::AbstractHierarchicalGridSystem,
        lvl::Int, a::Int)
    complete = ncells(grid) == ncells(sys, lvl)
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

# A subset's positions ascend with its ids, so the descendants of one ancestor
# occupy one interval of them, found by two binary searches.
function _subsetwindow(grid::PartialGrid, r::AbstractUnitRange)
    ids = grid.ids
    isempty(ids) && return 1:0
    lo = searchsortedfirst(ids, cellindex(grid.complete, Int(first(r))))
    hi = searchsortedlast(ids, cellindex(grid.complete, Int(last(r))))
    return Int(lo):Int(hi)
end

# ===========================================================================
# The `RegridSpace` contract
# ===========================================================================

ncells(space::DGGSpace) = ncells(space.grid)
getcell(space::DGGSpace, i::Int) = getcell(space.grid, i)

GOCore.manifold(space::DGGSpace) = GOCore.best_manifold(space.grid)

GR.nchunks(space::DGGSpace) = length(space.ranges)

cellindices(space::DGGSpace, chunk::Int) = space.ranges[chunk]

GR.cellcentroid(space::DGGSpace, i::Int) =
    cell_centroid(space.grid, cellindex(space.grid, i))

function cellat(space::DGGSpace, p::GO.UnitSphericalPoint)
    c = cellat(space.grid, p)
    c === nothing && return nothing
    return cellposition(space.grid, c)
end

GR.celltree(space::DGGSpace) = treeify(space.grid)

GR.chunktree(space::DGGSpace) = DGGChunkTree(space)

# The caps are stored, so the generic walk of the chunk tree has nothing to add.
GR.chunkextents(space::DGGSpace) = space.caps

# A DGGS destination is one-dimensional and its chunks are position ranges, so a
# chunk is already the read.
GR.chunkranges(space::DGGSpace, chunk::Integer, ::NTuple{1,Int}) =
    (space.ranges[Int(chunk)],)

"""
    GlobalRegridding.subtree(space::DGGSpace, inds)

The cell tree restricted to `inds`, with leaf indices still `space`'s own cell
positions.

A chunk keeps the hierarchy: the node of [`treeify`](@ref) rooted at that
chunk's ancestor cell, in `O(1)`. Any other index set falls back to
`GlobalRegridding`'s bounding-cap hierarchy, which costs one cell boundary per
index.
"""
function GR.subtree(space::DGGSpace, inds::AbstractUnitRange{<:Integer})
    GR._iswholespace(space, inds) && return GR.celltree(space)
    cursor = _chunkcursor(space, inds)
    cursor === nothing || return cursor
    return GR.CellCapTree(space, inds)
end

# The tree node that a chunk's positions already are: same grid, same leaf
# level, descent restarted at the chunk's ancestor cell over the chunk's own
# position window. Leaf indices stay grid positions, which is what the weight
# builder localizes against `dst_inds`/`src_inds`.
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

The one-node `SpatialTreeInterface` tree over `space`'s chunks, with stored
extents and chunk numbers as leaf indices.

A chunk's extent is its ancestor cell's [`node_extent`](@ref), which covers
every descendant at every depth and therefore every cell
[`cellindices`](@ref) assigns to the chunk.
"""
struct DGGChunkTree{S<:DGGSpace}
    space::S
    extent::_Cap
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

# ===========================================================================
# Resolving `to`
# ===========================================================================

"""
    regridgrid(x) -> AbstractGrid

The grid behind a regridding target: a grid is itself, a [`CellLookup`](@ref) or
[`CellVector`](@ref) is its [`PartialGrid`](@ref), and a
[`MultiOrderCellSet`](@ref) is the subset it expands to at its reference level.

A bare [`AbstractHierarchicalGridSystem`](@ref) has no method here: it names no
cells until a level is chosen, which [`arealevel`](@ref) does from the source.
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
        "chosen, and the level is chosen by matching the source's cell areas. " *
        "Pass a dimensional raster to `regrid` so the source can be measured, " *
        "or name the level yourself with `levelgrid(sys, l)`."))

"""
    arealevel(sys::AbstractHierarchicalGridSystem, space; samples = 256) -> Int

The level of `sys` whose mean cell area, `4π / ncells(sys, l)`, is closest in
ratio to the median cell area of the regridding space `space`.

This is what `to = sys` resolves through: the level that neither throws away the
source's resolution nor invents one it does not have.
"""
function arealevel(sys::AbstractHierarchicalGridSystem, space::GR.RegridSpace;
        samples::Integer=256)
    target = _mediancellarea(space, Int(samples))
    best, bestscore = first(levels(sys)), Inf
    for l in levels(sys)
        area = 4 * pi / ncells(sys, l)
        score = abs(log(area) - log(target))
        score < bestscore && ((best, bestscore) = (l, score))
        # Cell areas shrink with depth, so the first level at or below the
        # target brackets it and no deeper level can score better.
        area <= target && break
    end
    return best
end

# The median area of up to `samples` of the space's cells. Positions are walked
# on a golden-ratio stride rather than a fixed one, so a raster's row-major cell
# order cannot alias the sample onto a single column of latitudes.
function _mediancellarea(space::GR.RegridSpace, samples::Int)
    n = Int(ncells(space))
    n > 0 || throw(ArgumentError("cannot match cell areas against an empty space"))
    k = clamp(samples, 1, n)
    areas = Vector{Float64}(undef, k)
    if k == n
        for i in 1:n
            areas[i] = GR.cellarea(space, i)
        end
    else
        for j in 1:k
            areas[j] = GR.cellarea(space, mod1(round(Int, j * n * 0.6180339887498949), n))
        end
    end
    sort!(areas)
    return isodd(k) ? areas[(k + 1) ÷ 2] : (areas[k ÷ 2] + areas[k ÷ 2 + 1]) / 2
end

# ===========================================================================
# The verbs
# ===========================================================================

# `to` is a keyword, so the only place a bare system can be resolved is a method
# on `plan_regrid` itself: the level comes from the source, which the generic
# `_asspace` seam never sees. Every other spelling is resolved by `_asspace`
# above, and `regrid`/`regrid!` reach this through their own plan construction.
function GR.plan_regrid(data::DD.AbstractDimArray; to, from=nothing, kwargs...)
    return invoke(GR.plan_regrid, Tuple{Any}, data;
        to=_resolvetarget(to, data, from), from, kwargs...)
end

_resolvetarget(to, data, from) = to

_resolvetarget(sys::AbstractHierarchicalGridSystem, data, from) =
    DGGSpace(levelgrid(sys, arealevel(sys,
        from === nothing ? GR._sourcespace(data) : from)))

# A plan whose DESTINATION is a DGGS. The leading parameters are spelled with
# their declared bounds rather than `<:Any`: a type variable widened past a
# struct's own bound produces a type that is not a subtype of it, and the
# method would silently never be reached.
const _DirectToDGG =
    GR.DirectPlan{<:GR.AbstractRegriddingMethod,<:GR.AbstractMissingPolicy,<:DGGSpace}
const _ChunkedToDGG =
    GR.ChunkedPlan{<:GR.AbstractRegriddingMethod,<:GR.AbstractMissingPolicy,<:DGGSpace}

function GR.regrid(data, plan::_DirectToDGG)
    return _ascube(invoke(GR.regrid, Tuple{Any,GR.DirectPlan}, data, plan), data, plan)
end

GR.regrid(data, plan::_ChunkedToDGG) =
    _ascube(GR.LazyRegridArray(data, plan), data, plan)

# The destination axis is the cells themselves, not a bare `1:n`: dimensions are
# `(Cells(lk), the source's non-spatial dimensions...)`. A source that is not a
# cube has no dimensions to carry, and comes back as the plain array it was.
function _ascube(out, data, plan::GR.AbstractRegriddingPlan)
    data isa DD.AbstractDimArray || return out
    lk = CellLookup(plan.dst_space.grid)
    ds = DD.dims(data)
    sd = GR.resolvespatialdims(data, Int(ncells(plan.src_space)))
    others = Tuple(ds[i] for i in eachindex(ds) if !(i in sd))
    raw = out isa DD.AbstractDimArray ? parent(out) : out
    return DD.DimArray(raw, (Cells(lk), others...))
end
