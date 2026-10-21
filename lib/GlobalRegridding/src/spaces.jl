# A source or destination of a regrid is a GridInterface grid
# (`GridInterface.isgrid`). Cells, trees, chunks, point location, charts and
# dual cells are read through GridInterface's verbs at the sites that need
# them; this file holds the few private helpers those sites share and the
# execution hooks a package may extend for its grids.

# --------------------------------------------------------------------------
# Chunks, read off `chunks(grid)`
# --------------------------------------------------------------------------
#
# A grid that declares no chunking (`chunks(grid) === nothing`) is one chunk
# over every cell. Nothing is held by the grid: every answer is derived when it
# is asked for, and a chunk's extent is `grid_extent` of the chunk, cheap by
# contract. The chunk list itself need not be (a subset of a level grid searches
# its positions once per chunk to build it), so each helper takes it as `ch`,
# resolved once by the caller that loops or keeps it: a [`ChunkedPlan`](@ref)
# resolves its two spaces' lists at construction.

_chunks(grid) = GridInterface.chunks(grid)

_nchunks(grid, ch = _chunks(grid)) = ch === nothing ? 1 : length(ch)

function _ownedindices(grid, chunk::Integer, ch = _chunks(grid))
    ch === nothing || return _asrange(GridInterface.parentpositions(ch[Int(chunk)]))
    Int(chunk) == 1 || throw(BoundsError(grid, Int(chunk)))
    return 1:Int(ncells(grid))
end

# A subset's positions are strictly ascending, so a contiguous one is its range:
# a full-width tile then takes the arithmetic index map, not a lookup.
_asrange(pos::AbstractUnitRange) = pos
function _asrange(pos::AbstractVector{<:Integer})
    isempty(pos) && return pos
    lo, hi = Int(first(pos)), Int(last(pos))
    return hi - lo + 1 == length(pos) ? (lo:hi) : pos
end

# One cap per chunk, in chunk order, read when asked.
function _chunkextents(grid)
    ch = _chunks(grid)
    ch === nothing && return Cap[convert(Cap, GridInterface.grid_extent(grid))]
    return Cap[convert(Cap, GridInterface.grid_extent(c)) for c in ch]
end

function _chunkextent(grid, chunk::Integer, ch = _chunks(grid))
    ch === nothing || return convert(Cap, GridInterface.grid_extent(ch[Int(chunk)]))
    Int(chunk) == 1 || throw(BoundsError(grid, Int(chunk)))
    return convert(Cap, GridInterface.grid_extent(grid))
end

function _chunkat(grid, i::Integer, ch = _chunks(grid))
    1 <= i <= Int(ncells(grid)) || throw(BoundsError(grid, Int(i)))
    ch === nothing && return 1
    return Int(GridInterface.chunkat(ch, Int(i)))
end

function _chunkat(grid, p::US.UnitSphericalPoint)
    i = GridInterface.cellat(grid, p)
    i === nothing && return nothing
    return _chunkat(grid, Int(i))
end

# The chunk whose positions are exactly `inds`, or `nothing`.
function _ownchunk(grid, inds::AbstractUnitRange{<:Integer}, ch = _chunks(grid))
    (ch === nothing || isempty(inds)) && return nothing
    chunk = ch[GridInterface.chunkat(ch, Int(first(inds)))]
    owned = GridInterface.parentpositions(chunk)
    # Strictly ascending positions with the range's ends and length are the range.
    return (length(owned) == length(inds) && Int(first(owned)) == Int(first(inds)) &&
            Int(last(owned)) == Int(last(inds))) ? chunk : nothing
end

# The chunk lists a plan holds while it builds, for the helpers a build reaches
# with a space alone (`_subtree`). A dynamic scope rather than an argument: the
# space-level builders between the plan and the tree are extension points.
struct HeldChunks{D,CD,S,CS}
    dst_space::D
    dst_chunks::CD
    src_space::S
    src_chunks::CS
end
const HELD_CHUNKS = ScopedValue{Union{Nothing,HeldChunks}}(nothing)

# Run `f()` with `held` in scope.
_holdingchunks(f, held::HeldChunks) = @with HELD_CHUNKS => held f()

# `grid`'s chunk list from the scope, else resolved here.
_heldchunks(grid) = _heldchunks(HELD_CHUNKS[], grid)
_heldchunks(::Nothing, grid) = _chunks(grid)
_heldchunks(h::HeldChunks, grid) =
    grid === h.dst_space ? h.dst_chunks : grid === h.src_space ? h.src_chunks : _chunks(grid)

# --------------------------------------------------------------------------
# Trees
# --------------------------------------------------------------------------

# The grid's tree for an intent; the sweep tree is what the dual-tree searches descend.
_celltree(grid, hint = GridInterface.ForSweeps()) = GridInterface.treeify(grid; leaf_capacity = hint)

_iswholespace(grid, inds::AbstractUnitRange{<:Integer}) =
    first(inds) == 1 && length(inds) == Int(ncells(grid))
_iswholespace(grid, inds) = false

# A tree over `inds` whose leaves emit the grid's positions, with the leaf size
# of the caller's intent: `ForSweeps` (the default) for the dual-tree searches,
# `ForQueries` for point location. The whole grid is its tree; a chunk's own
# positions take the chunk, so a typed chunk (a raster or structured tile, a
# multi-level set's window) builds its own tree in O(1); any other index set is
# the grid's `subtree`: its own tree restricted and rooted at the deepest node
# holding the set, or an R-tree over the set for a grid with no tree of its own.
function _subtree(grid, inds::AbstractUnitRange{<:Integer};
    leaf_capacity = GridInterface.ForSweeps())
    _iswholespace(grid, inds) && return _celltree(grid, leaf_capacity)
    return _chunksubtree(grid, inds, _heldchunks(grid), leaf_capacity)
end
function _chunksubtree(grid, inds, ch, leaf_capacity)
    chunk = _ownchunk(grid, inds, ch)
    chunk === nothing || return GridInterface.subtree(grid, chunk; leaf_capacity)
    return GridInterface.subtree(grid, inds; leaf_capacity)
end
_subtree(grid, inds::AbstractVector{<:Integer}; leaf_capacity = GridInterface.ForSweeps()) =
    GridInterface.subtree(grid, inds; leaf_capacity)
_subtree(grid, inds; leaf_capacity = GridInterface.ForSweeps()) =
    GridInterface.subtree(grid, collect(Int, inds); leaf_capacity)

# --------------------------------------------------------------------------
# Execution hooks: plain functions with `Any` defaults
# --------------------------------------------------------------------------

"""
    chunkranges(grid, chunk, spatialsize::NTuple{NS,Int}) -> NTuple{NS,UnitRange{Int}}

Return the rectangular array ranges that storage can read for `chunk` in one
operation, in spatial-dimension order. Flattening that block must enumerate
the chunk's positions (`parentpositions(chunks(grid)[chunk])`) in the same
order. The default reads the chunk's positions and refuses a chunk that is not
a rectangle of the array lattice; a grid whose chunks are rectangles it can
name directly, as [`RasterGrid`](@ref) does, specializes this function.
"""
function chunkranges(grid, chunk::Integer, ::NTuple{1,Int})
    inds = _ownedindices(grid, Int(chunk))
    n = length(inds)
    n == 0 && return (1:0,)
    lo, hi = Int(first(inds)), Int(last(inds))
    hi - lo + 1 == n || _notrectangular(grid, chunk)
    return (lo:hi,)
end

function chunkranges(grid, chunk::Integer, ssz::NTuple{2,Int})
    inds = _ownedindices(grid, Int(chunk))
    n = length(inds)
    n == 0 && return (1:0, 1:0)
    n1 = ssz[1]
    c1, r1 = fldmod1(Int(first(inds)), n1)
    c2, r2 = fldmod1(Int(last(inds)), n1)
    rows, cols = r1:r2, c1:c2
    length(rows) * length(cols) == n || _notrectangular(grid, chunk)
    k = 0
    for c in cols, r in rows
        k += 1
        @inbounds Int(inds[k]) == r + (c - 1) * n1 || _notrectangular(grid, chunk)
    end
    return (rows, cols)
end

chunkranges(grid, chunk::Integer, ssz::Tuple) = throw(ArgumentError(
    "a regrid flattens over one or two spatial dimensions, got $(length(ssz))"))

_notrectangular(grid, chunk) = throw(ArgumentError(
    "chunk $chunk of $(typeof(grid)) is not a rectangle of the array lattice, " *
    "so it cannot be read in one block; define " *
    "`GlobalRegridding.chunkranges(::$(typeof(grid)), chunk, spatialsize)`"))

"""
    destinationdims(grid, sampling) -> Tuple or nothing

Return the `DimensionalData` dimensions that label a result over this grid,
fastest dimension first, with lookups carrying `sampling`. The default is
`nothing`: the result keeps one flat `Cell` axis over `1:ncells(grid)`. A
transparent wrapper answers for the grid it wraps.
"""
destinationdims(::Any, ::DD.Lookups.Sampling) = nothing
destinationdims(g::GridInterface.AbstractWrappedGrid, sampling::DD.Lookups.Sampling) =
    destinationdims(parent(g), sampling)

"""
    dimsource(lookup) -> `from` target or nothing

Return the source a lookup already names, or `nothing`. A lookup that carries
its own cells is not a raster axis, so a package that supplies one extends this
and a regrid given no `from` names it instead of asking for `xdim`.
"""
dimsource(::Any) = nothing

"""
    _asgrid(x, name) -> grid
    _asgrid(x, name, src_grid) -> grid

Resolve a `to` or `from` argument into a GridInterface grid. A grid stands for
itself, and a dimensional raster or a tuple of dimensions names the
[`RasterGrid`](@ref) over its lattice — a dimension that already names cells
([`dimsource`](@ref)) is refused rather than read as one. Packages that supply
grids extend the two-argument form for their own target spellings, and the
three-argument form when the destination depends on the resolved source grid.
`name` names the keyword in error messages.
"""
function _asgrid end
