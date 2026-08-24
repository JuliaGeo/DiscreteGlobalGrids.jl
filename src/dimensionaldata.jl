# DimensionalData wrappers for `Fallbacks.CellVector`. Cell ids are computed
# from compressed position windows rather than stored as an expanded vector.

"""
    CellLookups

The DimensionalData layer: [`CellLookup`](@ref), the [`Cells`](@ref) dimension,
and the [`Covering`](@ref) selector.

A [`CellLookup`](@ref) is a one-dimensional `DimensionalData` lookup over cell
ids at a single level. It is a thin wrapper around a [`CellVector`](@ref),
which is where the compression lives: a set of **leaf position windows** —
sorted, disjoint intervals (or, where intervals are unavailable, a sorted list)
of positions in `levelgrid(sys, leaf)`. Its logical content is their
concatenation, and every operation is arithmetic over that concatenation:
`length` sums the window lengths, `lk[k]` binary-searches the cumulative
lengths and resolves one `cellindex`, [`cellposition`](@ref) runs the inverse.
Nothing is materialised.

`CellVector` provides the storage and indexing behavior; this module provides
the lookup, dimension, and selectors required by DimensionalData.
"""
module CellLookups

import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, AbstractCellVector,
    ncells, cellindex, cellposition, cellat, level, system,
    levelgrid, cellindextype, has_sorted_subtrees, descendants, query,
    neighbors, ring, neighborcount, Connectivity, Vertex, maxneighbors,
    halo, border, interior, adjacency
import ..DiscreteGlobalGrids: Helpers
import ..DiscreteGlobalGrids.Engine: PartialGrid, SubtreeIds,
    MultiOrderCoverage, MultiOrderCellSet, level_ranges
# Core collection operations delegated to `CellVector`.
import ..DiscreteGlobalGrids.Engine: CellVector, cellset, covering,
    covering_positions, windows, nwindows, RangeWindows, CellWindows, _derive,
    _windows, SubsetPositionedCell, mapneighbors, foreachneighbors,
    StorageOrder, _capacity, _ringtype

import SmallCollections
import DimensionalData as DD
import DimensionalData: Dimensions, Lookups

# ===========================================================================
# The lookup
# ===========================================================================

"""
    abstract type AbstractCellLookup{ID} <: DimensionalData.Lookups.Lookup{ID,1}

A `DimensionalData` lookup naming cells at one level — the cube face of
[`AbstractCellVector`](@ref). `Base.parent` returns that vector, and every cell
verb a cube supports is defined once here and forwarded to it.

Two lookups ship, one per backing: [`CellLookup`](@ref) over a computed
[`CellVector`](@ref), and [`ChunkedCellLookup`](@ref) over a stored
[`ChunkedCellVector`](@ref). Code that means "the cell dimension of this cube"
dispatches on this type and accepts both; naming either concrete type accepts
only cubes from one source, which is how a cube from [`dggread`](@ref) comes to
be refused by an operation that works on the identical cube built in memory.

# Required interface

`Base.parent(lk)` returns an [`AbstractCellVector`](@ref), and the lookup's
`getindex`, `length` and `eltype` agree with it. Everything else — `system`,
`level`, `cellposition`, the neighbourhood and region verbs, `PartialGrid`,
regridding and plotting — is generic over that one method.

A subtype still writes its own `Lookups.rebuild` and `Lookups.selectindices`,
because what a SUBSET of it should be is a property of the backing: a computed
window set stays compressed, and a stored axis stops being stored.
"""
abstract type AbstractCellLookup{ID} <: Lookups.Lookup{ID,1} end

"""
    CellLookup(cv::CellVector)
    CellLookup(set::MultiOrderCellSet; level = set's reference level)
    CellLookup(grid::AbstractGrid)

A `DimensionalData` lookup naming cells at one level. Pair it with
[`Cells`](@ref) to make a cube axis:

```julia
set = query(sys, MultiOrderCoverage(region); level = 9)
lk  = CellLookup(set)
A   = DimensionalData.DimArray(values, Cells(lk))
```

Semantically `lk` is the leaf id vector: `length(lk)` is the number of leaf
cells, `lk[k]` is the `k`th of them, `collect(lk)` is the vector itself.

`CellLookup` stores only a [`CellVector`](@ref). The vector represents the set
as sorted, disjoint position windows at the leaf level ([`level_ranges`](@ref)),
using O(number of windows) memory instead of O(number of leaf cells). Lookup
operations delegate to the vector's methods, including `lk[k]`,
[`cellposition`](@ref), [`cellset`](@ref), [`covering`](@ref),
and [`PartialGrid`](@ref).

`Base.parent` returns the lookup's VALUES, as `DimensionalData` requires: the
`CellVector`, which is an `AbstractVector` of the ids, is O(#windows) and
materialises nothing. [`cellset`](@ref) returns the backing — the set, or the
grid — for running a second coverage operation against without unpacking the
lookup.

Accepted inputs are:

  - a [`MultiOrderCellSet`](@ref), optionally re-expanded
    to a deeper `level` than the set's own reference level;
  - `levelgrid(sys, l)`, a whole level, which is one window;
  - a [`PartialGrid`](@ref), an arbitrary ascending subset, which is that
    subset's positions — one window when the subset is a subtree, and the
    explicit list when it is scattered.

All forms construct a [`CellVector`](@ref); an existing vector can be passed
directly.

# Selectors

```julia
A[Cells(DimensionalData.At(c))]              # a typed cell id
A[Cells(DimensionalData.Contains(8.0, 46.5))] # a lon/lat point, through `cellat`
A[Cells(Covering(polygon))]                   # a region, through `MultiOrderCoverage`
```

`At` and `Contains` resolve to one position; [`Covering`](@ref) to the
positions of every stored cell the region's coverage names, and the view it
produces carries a `CellLookup` again. Outside a cube those three are
`cellposition(cv, c)`, `cellposition(cv, lon, lat)` and
[`covering`](@ref)`(cv, polygon)`.

`At` and `Contains` are referenced as `DD.At` and `DD.Contains`. They are not
re-exported because this package exports DE9IM's unrelated [`Contains`](@ref)
geometry predicate. [`Covering`](@ref) is exported by this package.

`DD.Near` throws: cell ids ascend along a space-filling curve, so snapping to
the nearest id is not snapping to the nearest cell on the sphere, and this
lookup has no nearest-member search to offer instead. `Contains(lon, lat)`
answers the question `Near` is usually reached for.

# What the cube's own operations do to it

Indexing, concatenation, and reductions preserve the most specific valid lookup:

  - an ASCENDING subset — a range, a sorted index vector, a boolean mask, a
    selector's result, or a concatenation of disjoint ascending axes — is a
    window set again, so it is a `CellLookup`;
  - a REORDERED one — `lk[[3, 1]]`, `reverse(A; dims = Cells)` — is not, and
    falls back to an unordered `DimensionalData.Categorical` of the same ids.
    The values are right and the selectors work; the type is not `CellLookup`,
    so reversing twice restores the data but not the axis's type;
  - a reduction (`sum(A; dims = Cells)`) collapses to a single element that no
    cell id names, so the axis becomes `NoLookup`;
  - `rebuild` around anything that is not cell ids at this level throws, since
    a `CellLookup` has no free fields to put them in. Replacing an axis
    wholesale is `set(A, Cells => NoLookup())`.

!!! note "Systems without sorted subtrees (A5) are built by selection"
    [`level_ranges`](@ref) throws where [`has_sorted_subtrees`](@ref) is
    `false`, because a cell's descendants are not one interval of their level.
    The lookup is then built by selection: `descendants` names the leaves, they
    are resolved to positions and sorted, and the result is run-compressed like
    any other position list. Every method above is unchanged and every law
    still holds. [`CellVector`](@ref) documents what that costs, and the one
    consequence it inherits: on A5 a [`Covering`](@ref) selection is a superset
    of the cells that meet the region, by the same margin the refinement is.
"""
struct CellLookup{ID,C<:CellVector} <: AbstractCellLookup{ID}
    cells::C
end

CellLookup(cv::CellVector{ID}) where {ID} = CellLookup{ID,typeof(cv)}(cv)

# Forward the `level` keyword explicitly because it shadows the function name.
CellLookup(set::MultiOrderCellSet; level::Integer=set.reference_level) =
    CellLookup(CellVector(set; level=level))

CellLookup(grid::AbstractGrid) = CellLookup(CellVector(grid))

CellLookup(lk::CellLookup) = lk

# Preserve the lookup wrapper when deriving a window subset.
_derive(lk::CellLookup, w::CellWindows) = CellLookup(_derive(parent(lk), w))

windows(lk::CellLookup) = windows(parent(lk))

# --- the collection surface ------------------------------------------------

# DimensionalData derives lookup behavior from `parent`, so it must return the
# logical values (`CellVector`). Use `cellset` to access the backing set or grid.
Base.parent(lk::CellLookup) = lk.cells
Base.IndexStyle(::Type{<:AbstractCellLookup}) = Base.IndexLinear()

# Empty lookups have no lower or upper value bound.
Lookups.bounds(lk::AbstractCellLookup) =
    isempty(lk) ? (nothing, nothing) : (first(lk), last(lk))

Base.@propagate_inbounds Base.getindex(lk::AbstractCellLookup, k::Int) = parent(lk)[k]
Base.@propagate_inbounds Base.getindex(lk::AbstractCellLookup, k::CartesianIndex{1}) =
    parent(lk)[k[1]]

# Only vector indices can produce a one-dimensional window set. Shaped indices
# use Base's generic array indexing.
for f in (:getindex, :view, :dotview)
    @eval Base.$f(lk::AbstractCellLookup, ::Colon) = lk
    @eval Base.$f(lk::AbstractCellLookup, i::AbstractVector{<:Integer}) = _subset(lk, i)
end

# Route reversal through lookup indexing so the result remains a valid lookup.
Base.reverse(lk::AbstractCellLookup) = lk[lastindex(lk):-1:firstindex(lk)]

# Resolve the method ambiguity with SmallCollections vector indexing.
Base.getindex(lk::AbstractCellLookup,
    i::SmallCollections.AbstractFixedOrSmallOrPackedVector{<:Integer}) = _subset(lk, i)

# Validate boolean-mask axes against the lookup so bounds errors identify the
# indexed lookup rather than its backing vector.
#
# The body is shared but the METHOD is written per concrete lookup: `Bool <:
# Integer`, so a mask method on `AbstractCellLookup` and a concrete catch-all
# `_subset(lk, idx)` are ambiguous on exactly the `BitVector` a mask arrives as.
function _masksubset(lk::AbstractCellLookup, mask::AbstractArray{Bool})
    axes(mask) == axes(lk) || throw(BoundsError(lk, (mask,)))
    return _subset(lk, findall(mask))
end

_subset(lk::CellLookup, mask::AbstractArray{Bool}) = _masksubset(lk, mask)

# Ascending subsets remain compressed; reordered subsets use an unordered
# DimensionalData categorical lookup.
function _subset(lk::CellLookup, idx)
    sub = parent(lk)[idx]
    sub isa CellVector && return CellLookup(sub)
    return Lookups.Categorical(sub; order=Lookups.Unordered())
end

# --- what the lookup is, in this package's own vocabulary ------------------

"""
    cellset(lk::AbstractCellLookup)

Return the [`MultiOrderCellSet`](@ref) or grid used to construct the lookup.

`Base.parent(lk)` returns the logical values as an [`AbstractCellVector`](@ref).
A subset produced by indexing or selection reports its [`PartialGrid`](@ref), as
does a lookup over a stored axis, which was constructed from no set at all.
"""
cellset(lk::AbstractCellLookup) = cellset(parent(lk))

"""
    system(lk::AbstractCellLookup)

The grid system the lookup's cells are named in.
"""
system(lk::AbstractCellLookup) = system(parent(lk))

# `CellLookup` neither adds cells nor changes adjacency; it carries the same
# static degree bound as its compressed vector.
maxneighbors(lk::AbstractCellLookup, connectivity::Connectivity) =
    maxneighbors(parent(lk), connectivity)
maxneighbors(lk::AbstractCellLookup) = maxneighbors(lk, Vertex())

"""
    level(lk::AbstractCellLookup) -> Int

The one level every cell in the lookup sits at.
"""
level(lk::AbstractCellLookup) = level(parent(lk))

"""
    cellposition(lk::AbstractCellLookup, c::AbstractCellIndex) -> Union{Int,Nothing}

Position of cell `c` in the lookup, or `nothing` when the lookup does not hold
it — including when `c` is at another level. The inverse of `lk[k]`, and the
half of the bijection every selector ends at.
"""
cellposition(lk::AbstractCellLookup, c::AbstractCellIndex) = cellposition(parent(lk), c)

"""
    neighbors(lk::AbstractCellLookup, c, k = 1; connectivity = Vertex())
    ring(lk::AbstractCellLookup, c, k; connectivity = Vertex())
    neighbors(lk::AbstractCellLookup, p::Int, k = 1; connectivity = Vertex())
    ring(lk::AbstractCellLookup, p::Int, k; connectivity = Vertex())
    halo(lk::AbstractCellLookup; connectivity = Vertex(), cells = false)
    border(lk::AbstractCellLookup; connectivity = Vertex(), cells = false)
    interior(lk::AbstractCellLookup; connectivity = Vertex(), cells = false)
    adjacency(lk::AbstractCellLookup; halo = 0, connectivity = Vertex(), threaded = true)

Return the backing vector's adjacency operations. Neighbour and ring results are
clipped to lookup membership, and the lookup is a region for the four region
verbs: [`halo`](@ref) walks outside it, [`border`](@ref) and [`interior`](@ref)
split what is inside, and [`adjacency`](@ref) tables the lot.
"""
neighbors(lk::AbstractCellLookup, c::AbstractCellIndex, k::Integer=1;
    connectivity::Connectivity=Vertex()) =
    neighbors(parent(lk), c, k; connectivity)

ring(lk::AbstractCellLookup, c::AbstractCellIndex, k::Integer;
    connectivity::Connectivity=Vertex()) = ring(parent(lk), c, k; connectivity)

@inline neighbors(lk::AbstractCellLookup, p::Int, k::Integer=1;
    connectivity::Connectivity=Vertex()) =
    neighbors(parent(lk), p, k; connectivity)

ring(lk::AbstractCellLookup, p::Int, k::Integer;
    connectivity::Connectivity=Vertex()) = ring(parent(lk), p, k; connectivity)

neighborcount(lk::AbstractCellLookup, c::AbstractCellIndex;
    connectivity::Connectivity=Vertex()) = neighborcount(parent(lk), c; connectivity)

halo(lk::AbstractCellLookup; kw...) = halo(parent(lk); kw...)
border(lk::AbstractCellLookup; kw...) = border(parent(lk); kw...)
interior(lk::AbstractCellLookup; kw...) = interior(parent(lk); kw...)
adjacency(lk::AbstractCellLookup; kw...) = adjacency(parent(lk); kw...)
adjacency(lk::AbstractCellLookup, hpos::AbstractVector{<:Integer}; kw...) =
    adjacency(parent(lk), hpos; kw...)

# Positioned handles use the parent vector's positions.
neighbors(lk::AbstractCellLookup; connectivity::Connectivity=Vertex()) =
    neighbors(parent(lk); connectivity)

# Delegate neighbourhood sweeps to the parent vector.
mapneighbors(f, lk::AbstractCellLookup; kw...) = mapneighbors(f, parent(lk); kw...)
mapneighbors(f, lk::AbstractCellLookup, data::AbstractVector; kw...) =
    mapneighbors(f, parent(lk), data; kw...)
foreachneighbors(f, lk::AbstractCellLookup; kw...) =
    foreachneighbors(f, parent(lk); kw...)
foreachneighbors(f, lk::AbstractCellLookup, data::AbstractVector; kw...) =
    foreachneighbors(f, parent(lk), data; kw...)

"""
    PartialGrid(lk::AbstractCellLookup) -> PartialGrid

The lookup read as a grid: position `k` of the grid is position `k` of the
lookup, so a `Regridder` built on it lines up with a cube over the lookup's
axis without a permutation. O(1) — the ids stay lazy.
"""
PartialGrid(lk::AbstractCellLookup) = PartialGrid(parent(lk))

DGG.region(lk::AbstractCellLookup) = DGG.region(parent(lk))

# --- DimensionalData plumbing ----------------------------------------------

# The values ascend in canonical id order, which is what makes the binary
# searches below — and `searchsortedfirst` on the lookup — sound.
Lookups.order(::AbstractCellLookup) = Lookups.ForwardOrdered()
Lookups.metadata(::AbstractCellLookup) = Lookups.NoMetadata()

# DimensionalData passes concatenated values through `rebuild`. Ascending cell
# ids remain compressed; other orders become an unordered categorical lookup.
function Lookups.rebuild(lk::CellLookup; data=nothing, kw...)
    (data === nothing || data === lk || data === parent(lk)) && return lk
    return _rebuild(lk, data)
end

_rebuild(lk::CellLookup, cv::CellVector) = _derive(lk, windows(cv))

function _rebuild(lk::CellLookup, ids::AbstractVector{<:AbstractCellIndex})
    cv = parent(lk)
    positions = Vector{Int}(undef, length(ids))
    ascending = true
    for (j, c) in enumerate(ids)
        p = cellposition(cv.grid, c)
        p === nothing && throw(ArgumentError(
            "$c is not a cell of levelgrid($(system(lk)), $(level(lk))), so it " *
            "cannot join a CellLookup at that level"))
        positions[j] = p
        j > 1 && positions[j] <= positions[j-1] && (ascending = false)
    end
    ascending || return Lookups.Categorical(collect(ids); order=Lookups.Unordered())
    return _derive(lk, _windows(positions))
end

@noinline _rebuild(lk::CellLookup, data) = throw(ArgumentError(
    "a CellLookup holds cell ids at one level; it cannot be rebuilt around " *
    "$(typeof(data)). Concatenate cell axes with `vcat`/`cat`, subset them by " *
    "indexing, and replace one wholesale with `set(A, Cells => NoLookup())`."))

# A reduced cell axis no longer corresponds to a cell id.
Lookups.reducelookup(::AbstractCellLookup) = Lookups.NoLookup(Base.OneTo(1))

# Use the window membership search instead of searching all logical cell ids.
Lookups.hasselection(lk::AbstractCellLookup, sel::Lookups.At{<:AbstractCellIndex}) =
    cellposition(lk, Lookups.val(sel)) !== nothing

Lookups.hasselection(lk::AbstractCellLookup,
    sel::Lookups.Contains{<:AbstractCellIndex}) =
    cellposition(lk, Lookups.val(sel)) !== nothing

Lookups.hasselection(lk::AbstractCellLookup,
    sel::Lookups.Contains{<:Tuple{Real,Real}}) =
    cellposition(parent(lk), Lookups.val(sel)...) !== nothing

Dimensions.format(lk::AbstractCellLookup, ::Type, values, axis::AbstractRange) = lk

Base.:(==)(a::CellLookup, b::CellLookup) = parent(a) == parent(b)

function Base.show(io::IO, lk::CellLookup)
    print(io, "CellLookup(", typeof(system(lk)).name.name, ", level=", level(lk),
        ", ncells=", length(lk), ", ", nwindows(windows(lk)),
        windows(lk) isa RangeWindows ? " windows)" : " positions)")
end

Base.show(io::IO, ::MIME"text/plain", lk::CellLookup) = show(io, lk)

# ===========================================================================
# The dimension
# ===========================================================================

"""
    Cells(x)

The `DimensionalData` dimension of a cube's cell axis: `Cells(lk)` where `lk`
is a [`CellLookup`](@ref), and `Cells(selector)` when indexing.

```julia
A = DimensionalData.DimArray(values, Cells(CellLookup(set)))
A[Cells(Covering(county))]
```
"""
DD.@dim Cells "Cells"

# ---------------------------------------------------------------------------
# Positioned handles index one-dimensional cell arrays directly by storage
# position, with no membership check — the trusted-position contract. Bare
# cells keep the resolved path through the selector machinery above.
# ---------------------------------------------------------------------------

const CellsArray = DD.AbstractDimArray{T,1,<:Tuple{<:Cells}} where {T}

Base.@propagate_inbounds Base.getindex(A::CellsArray, h::SubsetPositionedCell) =
    parent(A)[h.position]
Base.@propagate_inbounds Base.setindex!(A::CellsArray, x, h::SubsetPositionedCell) =
    setindex!(parent(A), x, h.position)

# Find the first `Cells` dimension for positioned-handle indexing.
function _handle_dimnum(A::DD.AbstractDimArray)
    for (i, d) in enumerate(DD.dims(A))
        d isa Cells && return i
    end
    throw(ArgumentError(
        "cannot index with a SubsetPositionedCell: no Cells dimension in " *
        "dims $(map(DD.name, DD.dims(A)))"))
end

@inline _handle_slice(A::DD.AbstractDimArray, ::Val{D}, p::Int) where {D} =
    view(A, ntuple(i -> i == D ? p : Colon(), Val(ndims(A)))...)

# On an N-D array, a positioned handle selects the slice at its position
# along the first `Cells` dimension, as a view. The 1-D methods above keep
# the scalar fast path; the position is trusted either way.
Base.getindex(A::DD.AbstractDimArray, h::SubsetPositionedCell) =
    _handle_slice(A, Val(_handle_dimnum(A)), h.position)
Base.view(A::DD.AbstractDimArray, h::SubsetPositionedCell) =
    _handle_slice(A, Val(_handle_dimnum(A)), h.position)

# ===========================================================================
# Whole-array entry points
#
# Neighborhood operations locate the cell dimension in any `AbstractDimArray`.
# ===========================================================================

# Use the first `CellLookup` dimension unless `spatialdim` selects one
# explicitly; failures throw an informative `ArgumentError` either way.
function _cells_dimnum(A::DD.AbstractDimArray, ::Nothing)
    for (i, d) in enumerate(DD.dims(A))
        DD.lookup(d) isa AbstractCellLookup && return i
    end
    throw(ArgumentError(
        "no cell dimension found: none of the dims $(map(DD.name, DD.dims(A))) " *
        "carries a cell lookup; pass spatialdim to name one"))
end

function _cells_dimnum(A::DD.AbstractDimArray, spatialdim)
    d = DD.dims(A, spatialdim)
    d === nothing && throw(ArgumentError(
        "array has no dimension matching spatialdim = $spatialdim; its dims " *
        "are $(map(DD.name, DD.dims(A)))"))
    lk = DD.lookup(d)
    lk isa AbstractCellLookup || throw(ArgumentError(
        "dimension $(DD.name(d)) carries a $(nameof(typeof(lk))) lookup, " *
        "not a cell lookup"))
    return DD.dimnum(A, spatialdim)
end

_rebuilt(A::DD.AbstractDimArray, out::Tuple) = map(o -> DD.rebuild(A; data = o), out)
_rebuilt(A::DD.AbstractDimArray, out) = DD.rebuild(A; data = out)

"""
    Neighbors()

Pass positioned cell handles to `f(cell, neighbors)`; the handles read data
by indexing the array. This is the default for [`mapneighbors`](@ref) and
[`foreachneighbors`](@ref) on dimensional arrays. `mapneighbors` returns one
result per cell, on the cell dimension.
"""
struct Neighbors end

"""
    Values()

Pass scalar values to `f(cell, value, neighbor_values)`. On an N-D array,
[`mapneighbors`](@ref) runs the stencil independently along the cell
dimension for each index of the other dimensions, and keeps `A`'s dimensions.
"""
struct Values end

"""
    NeighborSlices()

Pass views to `f(cell, slice, neighbor_slices)`, with the cell dimension
removed from each view. [`mapneighbors`](@ref) returns one result per cell,
on the cell dimension. A one-dimensional array is refused — its per-cell
slice is a scalar, which is [`Values`](@ref).
"""
struct NeighborSlices end   # Not `Slices`: Base exports that name.

@noinline _bad_pass(pass) = throw(ArgumentError(
    "pass must be Neighbors(), Values() or NeighborSlices(), " *
    "got $(typeof(pass))"))

_need_slices(A) = ndims(A) >= 2 || throw(ArgumentError(
    "NeighborSlices() needs at least two dimensions; a one-dimensional " *
    "array's per-cell slice is its scalar — use Values()"))

# Rebuild one-result-per-cell outputs with the original cell dimension.
_rebuilt_on_cells(A, d, out::Tuple) =
    map(o -> DD.rebuild(A; data = o, dims = (d,)), out)
_rebuilt_on_cells(A, d, out) = DD.rebuild(A; data = out, dims = (d,))

"""
    mapneighbors(f, A::AbstractDimArray; spatialdim = nothing, pass = Neighbors(),
                 order = StorageOrder(), threaded = true, connectivity = Vertex())

Apply `f` to each cell and its neighbors. The result uses `A`'s wrapper and
lookups. If `f` returns a concrete tuple, each component becomes an array.

`spatialdim` accepts any selector supported by `DimensionalData.dims`.
By default, the first dimension with a [`CellLookup`](@ref) is used. An
array without one, or a selector that misses or names a non-cell dimension,
is an `ArgumentError`.

`pass` controls the callback arguments and output shape:

- [`Neighbors`](@ref): positioned handles; one result per cell, on the cell
  dimension.
- [`Values`](@ref): scalar values; the same dimensions as `A`.
- [`NeighborSlices`](@ref): views across the other dimensions; one result per
  cell, on the cell dimension.
"""
function mapneighbors(f::F, A::DD.AbstractDimArray; spatialdim = nothing,
        pass = Neighbors(), order = StorageOrder(), threaded = true,
        connectivity::Connectivity = Vertex()) where {F}
    dnum = _cells_dimnum(A, spatialdim)
    return _map_dimarray(pass, f, A, dnum, order, threaded, connectivity)
end

function _map_dimarray(::Neighbors, f::F, A, dnum, order, threaded,
        conn) where {F}
    cv = parent(DD.lookup(A, dnum))
    out = mapneighbors(f, cv; order, threaded, connectivity = conn)
    return _rebuilt_on_cells(A, DD.dims(A)[dnum], out)
end

function _map_dimarray(::Values, f::F, A, dnum, order, threaded,
        conn) where {F}
    cv = parent(DD.lookup(A, dnum))
    ndims(A) == 1 && return _rebuilt(A,
        mapneighbors(f, cv, parent(A); order, threaded, connectivity = conn))
    return _rebuilt(A, _map_slices(f, A, dnum, cv, order, threaded, conn))
end

function _map_dimarray(::NeighborSlices, f::F, A, dnum, order, threaded,
        conn) where {F}
    _need_slices(A)
    cv = parent(DD.lookup(A, dnum))
    out = _map_cell_slices(f, A, Val(dnum), cv, order, threaded, conn)
    return _rebuilt_on_cells(A, DD.dims(A)[dnum], out)
end

_map_dimarray(pass, f, A, dnum, order, threaded, conn) = _bad_pass(pass)

# Function barrier: the cell dimension's number becomes a constant, so the
# slice views are concretely typed.
function _map_cell_slices(f::F, A, ::Val{D}, cv, order, threaded,
        conn) where {F,D}
    g = (c, nbrs) -> f(c, _handle_slice(A, Val(D), cellposition(c)),
        [_handle_slice(A, Val(D), cellposition(h)) for h in nbrs])
    return mapneighbors(g, cv; order, threaded, connectivity = conn)
end

# Run a separate buffered 1-D sweep for each non-cell index, so the
# CellVector kernels own all traversal and the slices cannot interact.
function _map_slices(f::F, A, dnum::Int, cv::AbstractCellVector, order, threaded,
        connectivity::Connectivity) where {F}
    data = parent(A)
    pre = CartesianIndices(axes(data)[1:(dnum-1)])
    post = CartesianIndices(axes(data)[(dnum+1):end])
    cap = _capacity(system(cv), connectivity)
    H = SubsetPositionedCell{eltype(cv)}
    T = Base.promote_op(f, H, eltype(A), _ringtype(cap, eltype(A)))
    outs = T <: Tuple && isconcretetype(T) ?
           ntuple(j -> similar(data, fieldtype(T, j)), fieldcount(T)) :
           similar(data, T)
    buf = Vector{eltype(A)}(undef, size(data, dnum))
    for jpost in post, jpre in pre
        copyto!(buf, view(data, jpre, :, jpost))
        _slice_store!(outs,
            mapneighbors(f, cv, buf; order, threaded, connectivity),
            jpre, jpost)
    end
    return outs
end

_slice_store!(outs::Tuple, res::Tuple, jpre, jpost) =
    (map((o, r) -> copyto!(view(o, jpre, :, jpost), r), outs, res); nothing)
_slice_store!(out::AbstractArray, res::AbstractVector, jpre, jpost) =
    (copyto!(view(out, jpre, :, jpost), res); nothing)

"""
    foreachneighbors(f, A::AbstractDimArray; spatialdim = nothing, pass = Neighbors(),
                     order = StorageOrder(), threaded = false,
                     connectivity = Vertex())

Call `f` for each cell and its neighbors without collecting results.
`spatialdim` and `pass` behave as in [`mapneighbors`](@ref).
"""
function foreachneighbors(f::F, A::DD.AbstractDimArray; spatialdim = nothing,
        pass = Neighbors(), order = StorageOrder(), threaded = false,
        connectivity::Connectivity = Vertex()) where {F}
    dnum = _cells_dimnum(A, spatialdim)
    _foreach_dimarray(pass, f, A, dnum, order, threaded, connectivity)
    return nothing
end

_foreach_dimarray(::Neighbors, f::F, A, dnum, order, threaded, conn) where {F} =
    foreachneighbors(f, parent(DD.lookup(A, dnum)); order, threaded,
        connectivity = conn)

function _foreach_dimarray(::Values, f::F, A, dnum, order, threaded,
        conn) where {F}
    cv = parent(DD.lookup(A, dnum))
    data = parent(A)
    if ndims(A) == 1
        foreachneighbors(f, cv, data; order, threaded, connectivity = conn)
        return nothing
    end
    pre = CartesianIndices(axes(data)[1:(dnum-1)])
    post = CartesianIndices(axes(data)[(dnum+1):end])
    buf = Vector{eltype(A)}(undef, size(data, dnum))
    for jpost in post, jpre in pre
        copyto!(buf, view(data, jpre, :, jpost))
        foreachneighbors(f, cv, buf; order, threaded, connectivity = conn)
    end
    return nothing
end

function _foreach_dimarray(::NeighborSlices, f::F, A, dnum, order, threaded,
        conn) where {F}
    _need_slices(A)
    return _foreach_cell_slices(f, A, Val(dnum), parent(DD.lookup(A, dnum)),
        order, threaded, conn)
end

_foreach_dimarray(pass, f, A, dnum, order, threaded, conn) = _bad_pass(pass)

function _foreach_cell_slices(f::F, A, ::Val{D}, cv, order, threaded,
        conn) where {F,D}
    g = (c, nbrs) -> (f(c, _handle_slice(A, Val(D), cellposition(c)),
        [_handle_slice(A, Val(D), cellposition(h)) for h in nbrs]); nothing)
    foreachneighbors(g, cv; order, threaded, connectivity = conn)
    return nothing
end

"""
    neighbors(A::AbstractDimArray; spatialdim = nothing, connectivity = Vertex())

Iterate over each cell and its positioned neighbor handles. The cell
dimension is selected as in [`mapneighbors`](@ref); the minted positions are
that dimension's axis positions.
"""
function neighbors(A::DD.AbstractDimArray; spatialdim = nothing,
        connectivity::Connectivity = Vertex())
    dnum = _cells_dimnum(A, spatialdim)
    return neighbors(parent(DD.lookup(A, dnum)); connectivity)
end

# ===========================================================================
# Selectors
#
# Point and id selectors resolve to one position; a region selector resolves to
# all matching positions. Methods include the selector value type because DimensionalData
# reads a `Tuple`-valued selector as a pair of interval endpoints and a
# `Vector`-valued one as an elementwise map — both of which a bare
# `(::CellLookup, ::Contains)` method would be ambiguous with.
#
# ===========================================================================

"""
    Covering(target)

Select every stored cell named by [`MultiOrderCoverage`](@ref) for `target` at
the lookup's level.

```julia
A[Cells(Covering(county))]          # a GeoInterface geometry
A[Cells(Covering(extent))]          # a lon/lat Extents.Extent
A[Cells(Covering(cap))]             # a GO.UnitSpherical.SphericalCap
```

`target` is anything [`query`](@ref) accepts. The result is the intersection of
the coverage's leaf expansion with the lookup, in ascending position order. The
resulting view retains a [`CellLookup`](@ref).

Outside a cube, the equivalent selection is `covering(cv, target)`, which
returns a
[`CellVector`](@ref), or `covering_positions(cv, target)` for the positions
alone. See that docstring for what the selection costs and for the
over-covering it inherits from the coverage itself.
"""
struct Covering{T} <: Lookups.ArraySelector{T}
    val::T
end

Base.show(io::IO, sel::Covering) =
    print(io, "Covering(", typeof(sel.val).name.name, ")")

Lookups.selectindices(lk::AbstractCellLookup, sel::Covering; kw...) =
    covering_positions(parent(lk), Lookups.val(sel))

Lookups.selectindices(lk::AbstractCellLookup, sel::Covering{<:AbstractVector};
    kw...) = covering_positions(parent(lk), Lookups.val(sel))

Lookups.selectindices(lk::AbstractCellLookup,
    sel::Lookups.At{<:AbstractCellIndex}; kw...) =
    _found(lk, cellposition(lk, Lookups.val(sel)), sel)

Lookups.selectindices(lk::AbstractCellLookup,
    sel::Lookups.Contains{<:AbstractCellIndex}; kw...) =
    _found(lk, cellposition(lk, Lookups.val(sel)), sel)

# Resolve a point to a leaf-grid cell, then search its position in the windows.
Lookups.selectindices(lk::AbstractCellLookup,
    sel::Lookups.Contains{<:Tuple{Real,Real}}; kw...) =
    _found(lk, cellposition(parent(lk), Lookups.val(sel)...), sel)

Lookups.selectindices(lk::AbstractCellLookup,
    sel::Lookups.At{<:Tuple{Real,Real}}; kw...) =
    _found(lk, cellposition(parent(lk), Lookups.val(sel)...), sel)

_found(lk::AbstractCellLookup, k::Int, sel) = k
_found(lk::AbstractCellLookup, ::Nothing, sel) =
    throw(Lookups.SelectorError(lk, sel))

# `Near` on a cell axis would have to mean "nearest on the sphere". Cell ids
# ascend along a space-filling curve, so the generic order-based snap is not
# that, and there is no cheap nearest-member search over an arbitrary subset.
@noinline _no_near(lk, sel) = throw(ArgumentError(
    "Near is not defined on a cell axis: cell ids run along a space-filling " *
    "curve, so the nearest id is not the nearest cell on the sphere. Use " *
    "At(cell) for one cell, Contains(lon, lat) for the cell holding a point, " *
    "or Covering(region) for a region."))

Lookups.selectindices(lk::AbstractCellLookup, sel::Lookups.Near; kw...) =
    _no_near(lk, sel)
Lookups.selectindices(lk::AbstractCellLookup, sel::Lookups.Near{<:AbstractVector};
    kw...) = _no_near(lk, sel)

# ===========================================================================
# Selector failures
#
# The default `SelectorError` prints the lookup's full type. A cell lookup can
# say what it holds in one line, and name the mismatch — wrong level, wrong
# system — that usually explains the miss.
# ===========================================================================

"""
    show_selector_error(io, lk, sel)

Print a failed cell-axis selection: the value that matched nothing, the axis it
was applied to, and the level or system mismatch behind it, if any.

`sel` is whatever `DimensionalData` put in the error — a selector or the bare
value it was carrying.
"""
function show_selector_error(io::IO, lk, sel)
    print(io, "SelectorError: ", sel, " selects no cell of ", lk)
    lo, hi = Lookups.bounds(lk)
    lo === nothing || print(io, ", whose ids run ", lo, " to ", hi)
    _mismatch(io, lk, sel isa Lookups.Selector ? Lookups.val(sel) : sel)
    println(io)
    return nothing
end

function _mismatch(io::IO, lk, c::AbstractCellIndex)
    sys = system(lk)
    if !(typeof(c) in DGG.cellindextypes(sys))
        print(io, "\n  ", nameof(typeof(c)), " is not a ",
            nameof(typeof(sys)), " id; this axis names cells as ",
            nameof(cellindextype(sys)))
    elseif level(c) != level(lk)
        print(io, "\n  the cell is at level ", level(c), ", not the axis's ",
            "level ", level(lk))
    end
    return nothing
end

_mismatch(io::IO, lk, val) = nothing

Base.showerror(io::IO, e::Lookups.SelectorError{<:CellLookup}) =
    show_selector_error(io, e.lookup, e.selector)

end # module CellLookups
