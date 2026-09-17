"""
    CellLookups

The DimensionalData layer: [`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup), the [`Cells`](@ref DiscreteGlobalGrids.CellLookups.Cells) dimension,
and the [`Covering`](@ref) selector.

A [`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup) is a one-dimensional `DimensionalData` lookup over cell
ids at a single level. It is a thin wrapper around a [`CellVector`](@ref),
which is where the compression lives: a set of **leaf index windows** —
sorted, disjoint intervals (or, where intervals are unavailable, a sorted list)
of indices in `levelgrid(sys, leaf)`. Its logical content is their
concatenation, and every operation is arithmetic over that concatenation:
`length` sums the window lengths, `lk[k]` binary-searches the cumulative
lengths and resolves one `cellindex`, [`localindex`](@ref DiscreteGlobalGrids.localindex) runs the inverse.
Nothing is materialised.

`CellVector` provides the storage and indexing behavior; this module provides
the lookup, dimension, and selectors required by DimensionalData.

[`MultiOrderLookup`](@ref) is the mixed-level face of the same layer, over a
[`MultiOrderVector`](@ref) instead of a `CellVector`.
"""
module CellLookups

import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, AbstractCellVector,
    ncells, cellindex, localindex, globalindex, cellat, level, system,
    levelgrid, cellindextype, has_sorted_subtrees, descendants, descendant_range,
    query, neighbors, ring, neighborcount, Connectivity, Vertex, maxneighbors,
    halo, border, interior, adjacency, DE9IMPredicate, QueryPredicate
import ..DiscreteGlobalGrids: Helpers
import ..DiscreteGlobalGrids.Engine: PartialGrid, SubtreeIds,
    MultiOrderCoverage, MultiOrderCellSet, level_ranges
import ..DiscreteGlobalGrids.Engine: CellVector, cellset, covering,
    covering_indices, predicate_indices, windows, nwindows, RangeWindows, CellWindows, _derive,
    _windows, SubsetIndexedCell, mapneighbors, foreachneighbors,
    StorageOrder, _capacity, _ringtype, Neighborhood, Disc, Ring, _steps, _checkneighborhood
import ..DiscreteGlobalGrids.Engine: MultiOrderVector, reference_level,
    covering_index, aggregate, coarsen, expand

import SmallCollections
import DimensionalData as DD
import DimensionalData: Dimensions, Lookups

# --- shared cell-axis surface -----------------------------------------------

# The DimensionalData surface every cell axis shares, single-level
# (`AbstractCellLookup`) or mixed-level (`MultiOrderLookup`). A method lives here
# when all it needs is `parent(lk)` being a cell container. Private.
abstract type AbstractCellAxis{ID} <: Lookups.Lookup{ID,1} end

# --- same-level lookup ------------------------------------------------------

"""
    abstract type AbstractCellLookup{ID} <: DimensionalData.Lookups.Lookup{ID,1}

A `DimensionalData` lookup naming cells at one level — the cube face of
[`AbstractCellVector`](@ref). `Base.parent` returns that vector, and every cell
verb a cube supports is defined once here and forwarded to it.

Two lookups ship, one per backing: [`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup) over a computed
[`CellVector`](@ref), and [`ChunkedCellLookup`](@ref DiscreteGlobalGrids.ChunkedLookups.ChunkedCellLookup) over a stored
[`ChunkedCellVector`](@ref DiscreteGlobalGrids.ChunkedLookups.ChunkedCellVector). Code that means "the cell dimension of this cube"
dispatches on this type and accepts both; naming either concrete type accepts
only cubes from one source, which is how a cube from [`dggread`](@ref DiscreteGlobalGrids.dggread) comes to
be refused by an operation that works on the identical cube built in memory.

# Required interface

`Base.parent(lk)` returns an [`AbstractCellVector`](@ref), and the lookup's
`getindex`, `length` and `eltype` agree with it. Everything else — `system`,
`level`, `localindex`, the neighbourhood and region verbs, `PartialGrid`,
regridding and plotting — is generic over that one method.

A subtype still writes its own subsetting and rebuild rules (`_subset`,
`_rebuild`), because what a SUBSET of it should be is a property of the backing:
a computed window set stays compressed, and a stored axis stops being stored.

Indexing, the `At`/`Contains`/`Covering` selectors and the `DimensionalData`
plumbing are shared with the mixed-level [`MultiOrderLookup`](@ref) and are
defined one level up, on the private supertype `AbstractCellAxis`. Only the
single-level verbs — `level`, topology, `PartialGrid`, `region` — are defined on
this type.
"""
abstract type AbstractCellLookup{ID} <: AbstractCellAxis{ID} end

"""
    CellLookup(cv::CellVector)
    CellLookup(set::MultiOrderCellSet; level = set's reference level)
    CellLookup(grid::AbstractGrid)

A `DimensionalData` lookup naming cells at one level. Pair it with
[`Cells`](@ref DiscreteGlobalGrids.CellLookups.Cells) to make a cube axis:

```julia
set = query(sys, MultiOrderCoverage(region); level = 9)
lk  = CellLookup(set)
A   = DimensionalData.DimArray(values, Cells(lk))
```

Semantically `lk` is the leaf id vector: `length(lk)` is the number of leaf
cells, `lk[k]` is the `k`th of them, `collect(lk)` is the vector itself.

`CellLookup` stores only a [`CellVector`](@ref). The vector represents the set
as sorted, disjoint index windows at the leaf level ([`level_ranges`](@ref)),
using O(number of windows) memory instead of O(number of leaf cells). Lookup
operations delegate to the vector's methods, including `lk[k]`,
[`localindex`](@ref DiscreteGlobalGrids.localindex), [`cellset`](@ref), [`covering`](@ref),
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
    subset's indices — one window when the subset is a subtree, and the
    explicit list when it is scattered.

All forms construct a `CellVector`; `Base.parent(lk)` returns it. Logical
indexing and [`localindex`](@ref) remain O(log(number of windows)).

# Selectors

```julia
A[Cells(DimensionalData.At(c))]              # a typed cell id
A[Cells(DimensionalData.Contains(8.0, 46.5))] # a lon/lat point, through `cellat`
A[Cells(Covering(polygon))]                   # a region, through `MultiOrderCoverage`
```

`At` and `Contains` resolve to one index; [`Covering`](@ref) to the
indices of every stored cell the region's coverage names, and the view it
produces carries a `CellLookup` again. Outside a cube those three are
`localindex(cv, c)`, `localindex(cv, lon, lat)` and
[`covering`](@ref)`(cv, polygon)`.

`At` and `Contains` are referenced as `DD.At` and `DD.Contains`. They are not
re-exported because this package exports DE9IM's unrelated `Contains`
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
  - a reordered subset uses an unordered `DimensionalData.Categorical`;
  - a reduction (`sum(A; dims = Cells)`) collapses to a single element that no
    cell id names, so the axis becomes `NoLookup`;
  - rebuilding with incompatible data throws; use
    `set(A, Cells => NoLookup())` to replace the axis.

Systems without sorted subtrees build the backing vector by enumerating,
sorting, and compressing descendants. [`CellVector`](@ref) documents the cost
and non-congruent coverage behavior.
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

# DimensionalData requires `parent` to return the logical lookup values.
Base.parent(lk::CellLookup) = lk.cells
Base.IndexStyle(::Type{<:AbstractCellAxis}) = Base.IndexLinear()

# Empty lookups have no lower or upper value bound.
Lookups.bounds(lk::AbstractCellLookup) =
    isempty(lk) ? (nothing, nothing) : (first(lk), last(lk))

Base.@propagate_inbounds Base.getindex(lk::AbstractCellAxis, k::Int) = parent(lk)[k]
Base.@propagate_inbounds Base.getindex(lk::AbstractCellAxis, k::CartesianIndex{1}) =
    parent(lk)[k[1]]

# Base handles shaped indices; vector indices may preserve compressed windows.
for f in (:getindex, :view, :dotview)
    @eval Base.$f(lk::AbstractCellAxis, ::Colon) = lk
    @eval Base.$f(lk::AbstractCellAxis, i::AbstractVector{<:Integer}) =
        _subset(lk, _subsetindices(lk, i))
end

# Route reversal through lookup indexing so the result remains a valid lookup.
Base.reverse(lk::AbstractCellAxis) = lk[lastindex(lk):-1:firstindex(lk)]

# Resolve the method ambiguity with SmallCollections vector indexing.
Base.getindex(lk::AbstractCellAxis,
    i::SmallCollections.AbstractFixedOrSmallOrPackedVector{<:Integer}) =
    _subset(lk, _subsetindices(lk, i))

# A mask resolves to indices before the per-backing `_subset`: a mask method of
# `_subset` itself would be ambiguous with every backing's `_subset(lk, idx)`.
# The axes are checked first because `findall` discards them.
_subsetindices(lk::AbstractCellAxis, idx) = idx
function _subsetindices(lk::AbstractCellAxis, mask::AbstractArray{Bool})
    axes(mask) == axes(lk) || throw(BoundsError(lk, (mask,)))
    return findall(mask)
end

# Only ascending subsets retain compressed interval order.
function _subset(lk::CellLookup, idx)
    sub = parent(lk)[idx]
    sub isa CellVector && return CellLookup(sub)
    return Lookups.Categorical(sub; order=Lookups.Unordered())
end

# --- what the lookup is, in this package's own vocabulary ------------------

"""
    cellset(lk::AbstractCellLookup)
    cellset(lk::MultiOrderLookup)

Return the set or grid that constructed the lookup. Derived and stored-axis
lookups return their describing [`PartialGrid`](@ref). A mixed-level lookup
returns the [`MultiOrderVector`](@ref) backing it, the same collection as
`Base.parent(lk)`.
"""
cellset(lk::AbstractCellAxis) = cellset(parent(lk))

"""
    system(lk::AbstractCellLookup)
    system(lk::MultiOrderLookup)

The grid system the lookup's cells are named in.
"""
system(lk::AbstractCellAxis) = system(parent(lk))

# Wrapping preserves the backing vector's neighbor bound.
maxneighbors(lk::AbstractCellLookup, connectivity::Connectivity) =
    maxneighbors(parent(lk), connectivity)
maxneighbors(lk::AbstractCellLookup) = maxneighbors(lk, Vertex())

"""
    level(lk::AbstractCellLookup) -> Int

The one level every cell in the lookup sits at.
"""
level(lk::AbstractCellLookup) = level(parent(lk))

"""
    localindex(lk::AbstractCellLookup, c::AbstractCellIndex) -> Union{Int,Nothing}
    localindex(lk::MultiOrderLookup, c::AbstractCellIndex) -> Union{Int,Nothing}

Index of cell `c` in the lookup, or `nothing` when the lookup does not hold
it — including when `c` is at another level, or when a mixed-level lookup
stores only an ancestor or descendant of `c`. The inverse of `lk[k]`, and the
half of the bijection every selector ends at; `DimensionalData.At` resolves
through it.
"""
localindex(lk::AbstractCellAxis, c::AbstractCellIndex) = localindex(parent(lk), c)
globalindex(lk::AbstractCellLookup, c::AbstractCellIndex) = globalindex(parent(lk), c)

"""
    neighbors(lk::AbstractCellLookup, c, k = 1; connectivity = Vertex())
    ring(lk::AbstractCellLookup, c, k; connectivity = Vertex())
    neighbors(lk::AbstractCellLookup, p::Int, k = 1; connectivity = Vertex())
    ring(lk::AbstractCellLookup, p::Int, k; connectivity = Vertex())
    halo(lk::AbstractCellLookup; connectivity = Vertex(), cells = false)
    border(lk::AbstractCellLookup; connectivity = Vertex(), cells = false)
    interior(lk::AbstractCellLookup; connectivity = Vertex(), cells = false)
    adjacency(lk::AbstractCellLookup; halo = 0, connectivity = Vertex(), threaded = true)

Delegate topology to the backing vector. Neighbor and ring results are clipped
to lookup membership; region verbs treat the lookup as the same subset.
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

# Indexed handles use the parent vector's indices.
neighbors(lk::AbstractCellLookup; connectivity::Connectivity=Vertex(),
        neighborhood=Disc(1)) =
    neighbors(parent(lk); connectivity, neighborhood)

mapneighbors(f, lk::AbstractCellLookup; kw...) = mapneighbors(f, parent(lk); kw...)
mapneighbors(f, lk::AbstractCellLookup, data::AbstractVector; kw...) =
    mapneighbors(f, parent(lk), data; kw...)
foreachneighbors(f, lk::AbstractCellLookup; kw...) =
    foreachneighbors(f, parent(lk); kw...)
foreachneighbors(f, lk::AbstractCellLookup, data::AbstractVector; kw...) =
    foreachneighbors(f, parent(lk), data; kw...)

"""
    PartialGrid(lk::AbstractCellLookup) -> PartialGrid

The lookup read as a grid: index `k` of the grid is index `k` of the
lookup, so a `Regridder` built on it lines up with a cube over the lookup's
axis without a permutation. O(1) — the ids stay lazy.
"""
PartialGrid(lk::AbstractCellLookup) = PartialGrid(parent(lk))

DGG.region(lk::AbstractCellLookup) = DGG.region(parent(lk))

# --- DimensionalData plumbing ----------------------------------------------

# Canonical id order satisfies DimensionalData's ordered-lookup contract.
Lookups.order(::AbstractCellLookup) = Lookups.ForwardOrdered()
Lookups.metadata(::AbstractCellAxis) = Lookups.NoMetadata()

# Rebuild preserves compression only for ascending ids.
function Lookups.rebuild(lk::AbstractCellAxis; data=nothing, kw...)
    (data === nothing || data === lk || data === parent(lk)) && return lk
    return _rebuild(lk, data)
end

_rebuild(lk::CellLookup, cv::CellVector) = _derive(lk, windows(cv))

function _rebuild(lk::CellLookup, ids::AbstractVector{<:AbstractCellIndex})
    cv = parent(lk)
    indices = Vector{Int}(undef, length(ids))
    ascending = true
    for (j, c) in enumerate(ids)
        p = globalindex(cv.grid, c)
        p === nothing && throw(ArgumentError(
            "$c is not a cell of levelgrid($(system(lk)), $(level(lk))), so it " *
            "cannot join a CellLookup at that level"))
        indices[j] = p
        j > 1 && indices[j] <= indices[j-1] && (ascending = false)
    end
    ascending || return Lookups.Categorical(collect(ids); order=Lookups.Unordered())
    return _derive(lk, _windows(indices))
end

@noinline _rebuild(lk::CellLookup, data) = throw(ArgumentError(
    "a CellLookup holds cell ids at one level; it cannot be rebuilt around " *
    "$(typeof(data)). Concatenate cell axes with `vcat`/`cat`, subset them by " *
    "indexing, and replace one wholesale with `set(A, Cells => NoLookup())`."))

# A reduced cell axis no longer corresponds to a cell id.
Lookups.reducelookup(::AbstractCellAxis) = Lookups.NoLookup(Base.OneTo(1))

# A broadcast combines its operands' axes with `Lookups.promote_first`, whose
# `Lookup` fallback keeps a lookup only when the operands' CONCRETE types are
# identical and answers `NoLookup` otherwise. Two cell axes over the same cells
# differ in type routinely — a set-backed window run against a grid-backed one,
# a stored axis against the same cells in memory — so `a .- b` was dropping an
# axis neither operand had changed. Equal ids are the same axis whatever the
# backing; `Categorical` and `Sampled` each say so for themselves, and this is
# the cell axis saying it.
Lookups.promote_first(lk::AbstractCellLookup) = lk

# The same-type case is DimensionalData's own; it is written here too because
# the general method below would otherwise be ambiguous with it.
Lookups.promote_first(lk::L, ::L, ::L...) where {L<:AbstractCellLookup} = lk

# Unequal cell axes degrade to `NoLookup`, which is what `Sampled` does when it
# cannot promote. Under DimensionalData's default strict broadcast they never
# get this far: `comparedims` compares the values first and throws.
Lookups.promote_first(l1::AbstractCellLookup, l2::AbstractCellLookup,
    ls::AbstractCellLookup...) =
    all(==(l1), (l2, ls...)) ? l1 : Lookups.NoLookup(Base.OneTo(length(l1)))

# Use the window membership search instead of searching all logical cell ids.
Lookups.hasselection(lk::AbstractCellAxis, sel::Lookups.At{<:AbstractCellIndex}) =
    localindex(lk, Lookups.val(sel)) !== nothing

Lookups.hasselection(lk::AbstractCellLookup,
    sel::Lookups.Contains{<:AbstractCellIndex}) =
    localindex(lk, Lookups.val(sel)) !== nothing

Lookups.hasselection(lk::AbstractCellAxis,
    sel::Lookups.Contains{<:Tuple{Real,Real}}) =
    localindex(parent(lk), Lookups.val(sel)...) !== nothing

Dimensions.format(lk::AbstractCellAxis, ::Type, values, axis::AbstractRange) = lk

Base.:(==)(a::CellLookup, b::CellLookup) = parent(a) == parent(b)

function Base.show(io::IO, lk::CellLookup)
    print(io, "CellLookup(", typeof(system(lk)).name.name, ", level=", level(lk),
        ", ncells=", length(lk), ", ", nwindows(windows(lk)),
        windows(lk) isa RangeWindows ? " windows)" : " indices)")
end

Base.show(io::IO, ::MIME"text/plain", lk::CellLookup) = show(io, lk)

# --- cell dimension ---------------------------------------------------------

"""
    Cells(x)

The `DimensionalData` dimension of a cube's cell axis: `Cells(lk)` where `lk`
is a [`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup), and `Cells(selector)` when indexing.

```julia
A = DimensionalData.DimArray(values, Cells(CellLookup(set)))
A[Cells(Covering(county))]
```
"""
DD.@dim Cells "Cells"

# Indexed handles carry trusted storage positions from topology iterators.

const CellsArray = DD.AbstractDimArray{T,1,<:Tuple{<:Cells}} where {T}

Base.@propagate_inbounds Base.getindex(A::CellsArray, h::SubsetIndexedCell) =
    parent(A)[h.index]
Base.@propagate_inbounds Base.setindex!(A::CellsArray, x, h::SubsetIndexedCell) =
    setindex!(parent(A), x, h.index)

# `Base.filter` reaches for a `LogicalIndex` over the whole array whenever the
# array is not `IndexLinear` — which a cube whose data is a store is not — and
# a chunked parent cannot answer one. Resolve the mask to indices instead:
# indexing a cell axis by ascending indices is already defined, and it is what
# `filter` on a `Sampled` axis does, so the result is the surviving values on
# the axis of the cells they belong to.
Base.filter(f, A::CellsArray) = A[findall(f, parent(A))]

# Find the first `Cells` dimension for indexed-handle indexing.
function _handle_dimnum(A::DD.AbstractDimArray)
    for (i, d) in enumerate(DD.dims(A))
        d isa Cells && return i
    end
    throw(ArgumentError(
        "cannot index with a SubsetIndexedCell: no Cells dimension in " *
        "dims $(map(DD.name, DD.dims(A)))"))
end

@inline _handle_slice(A::DD.AbstractDimArray, ::Val{D}, p::Int) where {D} =
    view(A, ntuple(i -> i == D ? p : Colon(), Val(ndims(A)))...)

Base.getindex(A::DD.AbstractDimArray, h::SubsetIndexedCell) =
    _handle_slice(A, Val(_handle_dimnum(A)), h.index)
Base.view(A::DD.AbstractDimArray, h::SubsetIndexedCell) =
    _handle_slice(A, Val(_handle_dimnum(A)), h.index)

# --- whole-array entry points ----------------------------------------------

function _cells_dimnum(A::DD.AbstractDimArray, ::Nothing)
    for (i, d) in enumerate(DD.dims(A))
        DD.lookup(d) isa AbstractCellLookup && return i
    end
    throw(ArgumentError(
        "no cell dimension found: none of the dims $(map(DD.name, DD.dims(A))) " *
        "carries a cell lookup; pass the cell dimension to name one"))
end

# `DimensionalData.dims` accepts a tuple of dimensions; a cell axis is one
# dimension, so only the one-element tuple has a meaning here. The positional
# forms below take the tuple untyped: `(Cells,)` is a `Tuple{UnionAll}` to
# dispatch, which no `Tuple{Type{<:Dimension}}` signature admits.
const DimSelector = Union{Symbol, DD.Dimension, Type{<:DD.Dimension}, Tuple}

function _cells_dimnum(A::DD.AbstractDimArray, spatialdim::Tuple)
    length(spatialdim) == 1 || throw(ArgumentError(
        "a cell dimension is one dimension, not $(length(spatialdim)); " *
        "name it alone"))
    return _cells_dimnum(A, spatialdim[1])
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

# --- cell-field cube inputs -------------------------------------------------

DGG.Engine._cellknown(A::DD.AbstractDimArray, cv::CellVector, ::Type{T}) where {T} =
    _cubeknown(A, cv, T)
# Prefer the cube method over the dense-vector method for one-dimensional arrays.
DGG.Engine._cellknown(A::DD.AbstractDimArray{<:Any,1}, cv::CellVector,
    ::Type{T}) where {T} = _cubeknown(A, cv, T)

function _cubeknown(A::DD.AbstractDimArray, cv::CellVector, ::Type{T}) where {T}
    ndims(A) == 1 || throw(ArgumentError(
        "a cube `known` names one value per cell, so it is one-dimensional; " *
        "got dims $(map(DD.name, DD.dims(A)))"))
    lk = DD.lookup(A, _cells_dimnum(A, nothing))
    sub = parent(lk)
    system(sub) == system(cv) && level(sub) == level(cv) || throw(ArgumentError(
        "a cube `known` names cells of the collection being swept, so it " *
        "carries the same system and level: got $(system(sub)) level " *
        "$(level(sub)) against $(system(cv)) level $(level(cv))"))
    eltype(A) <: T || throw(ArgumentError(
        "`known` holds $(eltype(A)) where the field's element type is $T"))
    return DGG.Engine._SubsetKnown(sub, parent(A))
end

"""
    Neighbors()

Pass indexed cell handles as `f(cell, neighbors)`. This is the default for
[`mapneighbors`](@ref) and [`foreachneighbors`](@ref).
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

Pass cell-axis-free views as `f(cell, slice, neighbor_slices)`. This mode
requires at least two dimensions; use [`Values`](@ref) for scalar slices.
"""
struct NeighborSlices end   # Not `Slices`: Base exports that name.

@noinline _bad_pass(pass) = throw(ArgumentError(
    "pass must be Neighbors(), Values() or NeighborSlices(), " *
    "got $(typeof(pass))"))

_need_slices(A) = ndims(A) >= 2 || throw(ArgumentError(
    "NeighborSlices() needs at least two dimensions; a one-dimensional " *
    "array's per-cell slice is its scalar — use Values()"))

_rebuilt_on_cells(A, d, out::Tuple) =
    map(o -> DD.rebuild(A; data = o, dims = (d,)), out)
_rebuilt_on_cells(A, d, out) = DD.rebuild(A; data = out, dims = (d,))

@noinline _needs_pass(pass) = throw(ArgumentError(
    "needs cannot be combined with pass = $(typeof(pass)): a field request " *
    "already names what the callback receives; drop one of the two"))

# Field requests define callback arguments and require the default pass mode.
_checkpass(::Neighbors) = nothing
_checkpass(pass) = _needs_pass(pass)

"""
    mapneighbors(f, A::AbstractDimArray; spatialdim = nothing, pass = Neighbors(),
                 order = StorageOrder(), threaded = true, connectivity = Vertex(),
                 neighborhood = Disc(1))
    mapneighbors(f, A::AbstractDimArray; needs = (Value(a), Centroid()), ...)

Apply a neighborhood callback along a cell dimension. The default dimension is the
first cell lookup; `spatialdim` accepts a DimensionalData dimension selector.
A missing or non-cell dimension raises `ArgumentError`.
`neighborhood = Disc(k)` selects cells within `k` steps, excluding the center.
`Ring(k)` selects cells at exactly `k` steps. The default is `Disc(1)`.

[`Neighbors`](@ref) passes handles and returns one result per cell.
[`Values`](@ref) passes scalar values and preserves all input dimensions.
[`NeighborSlices`](@ref) passes views across other dimensions and returns one
result per cell; it requires at least two dimensions. Concrete tuple returns
produce one array per component, using the input wrapper and relevant lookups.

`needs` replaces the pass contract with `f(center, rings)` and produces one
result per cell. Values come from its `Value` requests, not implicitly from
`A`. Combining `needs` with a nondefault `pass` raises `ArgumentError`.

Stored arrays use chunked execution for `Values()` or `needs` with storage order.
`Index(Local())` still refers to the original cell axis. A permutation order
uses the ordinary traversal instead. See [Neighbours and stencils](@ref) for
mode examples and [Workflow execution details](@ref) for routing details.
"""
function mapneighbors(f::F, A::DD.AbstractDimArray; spatialdim = nothing,
        needs = nothing, pass = Neighbors(), order = StorageOrder(),
        threaded = true, connectivity::Connectivity = Vertex(),
        neighborhood = Disc(1)) where {F}
    dnum = _cells_dimnum(A, spatialdim)
    nb = _checkneighborhood(neighborhood)
    return _map_needs(needs, pass, f, A, dnum, order, threaded, connectivity, nb)
end

# No field request: `pass` picks the callback form, exactly as before —
# `needs = nothing` reaches those methods by dispatch, not by a branch.
_map_needs(::Nothing, pass, f::F, A, dnum, order, threaded, conn, nb) where {F} =
    _map_dimarray(pass, f, A, dnum, order, threaded, conn, nb)

# A field request is answered by the cell axis alone, whatever `A`'s
# dimensionality: the callback's values arrive through the request's `Value`
# entries, so the result is one per cell, on the cell dimension.
#
# A cube whose data is chunked on disk still takes the chunk route, under the
# rule `Values()` uses — the request may name `A` itself, or any other stored
# array, as a `Value`, and reading those along the chunk lines is the same win.
# The route is `chunks.jl`'s; a permutation `order` keeps the whole-axis path,
# because a chunked sweep visits by chunk and cannot honour one.
function _map_needs(needs, pass, f::F, A, dnum, order, threaded, conn,
        nb) where {F}
    _checkpass(pass)
    plan = _chunkedvalues(A, dnum, order, conn, nb)
    plan === nothing || return _rebuilt_on_cells(A, DD.dims(A)[dnum],
        _map_needs_chunked(f, A, dnum, plan, needs, threaded, conn, nb))
    cv = parent(DD.lookup(A, dnum))
    out = mapneighbors(f, cv; needs, order, threaded, connectivity = conn,
        neighborhood = nb)
    return _rebuilt_on_cells(A, DD.dims(A)[dnum], out)
end

function _map_needs_chunked end
function _foreach_needs_chunked end

function _map_dimarray(::Neighbors, f::F, A, dnum, order, threaded,
        conn, nb) where {F}
    cv = parent(DD.lookup(A, dnum))
    out = mapneighbors(f, cv; order, threaded, connectivity = conn,
        neighborhood = nb)
    return _rebuilt_on_cells(A, DD.dims(A)[dnum], out)
end

function _map_dimarray(::Values, f::F, A, dnum, order, threaded,
        conn, nb) where {F}
    plan = _chunkedvalues(A, dnum, order, conn, nb)
    plan === nothing || return _rebuilt(A,
        _map_values_chunked(f, A, dnum, plan, threaded, conn, nb))
    cv = parent(DD.lookup(A, dnum))
    ndims(A) == 1 && return _rebuilt(A,
        mapneighbors(f, cv, parent(A); order, threaded, connectivity = conn,
            neighborhood = nb))
    return _rebuilt(A, _map_slices(f, A, dnum, cv, order, threaded, conn, nb))
end

# Values flow through the traversal rather than being fetched by the callback,
# which is what lets a cube whose data is chunked on disk be swept along those
# chunks instead of cell by cell. `chunks.jl` owns the plan and answers these
# two; everything without a chunk grid takes the direct path above.
#
# `Neighbors()` deliberately has no such route: its callback closes over the
# ORIGINAL array and indexes it by the handles it is given, so a sweep over
# blocks would hand it indices into a block and it would read them from the
# whole cube. `foreachchunk` is the chunk-following form of that pass.
_chunkedvalues(A, dnum, order, conn, nb) = nothing
function _map_values_chunked end

function _map_dimarray(::NeighborSlices, f::F, A, dnum, order, threaded,
        conn, nb) where {F}
    _need_slices(A)
    cv = parent(DD.lookup(A, dnum))
    out = _map_cell_slices(f, A, Val(dnum), cv, order, threaded, conn, nb)
    return _rebuilt_on_cells(A, DD.dims(A)[dnum], out)
end

_map_dimarray(pass, f, A, dnum, order, threaded, conn, nb) = _bad_pass(pass)

# A value-type dimension number keeps slice views concrete.
function _map_cell_slices(f::F, A, ::Val{D}, cv, order, threaded,
        conn, nb) where {F,D}
    g = (c, nbrs) -> f(c, _handle_slice(A, Val(D), localindex(c)),
        [_handle_slice(A, Val(D), localindex(h)) for h in nbrs])
    return mapneighbors(g, cv; order, threaded, connectivity = conn,
        neighborhood = nb)
end

# Independent buffered sweeps isolate each non-cell slice.
function _map_slices(f::F, A, dnum::Int, cv::AbstractCellVector, order, threaded,
        connectivity::Connectivity, nb::Neighborhood) where {F}
    data = parent(A)
    pre = CartesianIndices(axes(data)[1:(dnum-1)])
    post = CartesianIndices(axes(data)[(dnum+1):end])
    cap = _capacity(cv.grid, nb, connectivity)
    H = SubsetIndexedCell{eltype(cv)}
    T = Base.promote_op(f, H, eltype(A), _ringtype(cap, eltype(A)))
    outs = T <: Tuple && isconcretetype(T) ?
           ntuple(j -> similar(data, fieldtype(T, j)), fieldcount(T)) :
           similar(data, T)
    buf = Vector{eltype(A)}(undef, size(data, dnum))
    for jpost in post, jpre in pre
        copyto!(buf, view(data, jpre, :, jpost))
        _slice_store!(outs,
            mapneighbors(f, cv, buf; order, threaded, connectivity,
                neighborhood = nb),
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
                     connectivity = Vertex(), neighborhood = Disc(1))
    foreachneighbors(f, A::AbstractDimArray; needs = (Value(a), Centroid()), ...)

Call `f` for each cell and its neighbourhood without collecting results.
`spatialdim`, `pass`, `needs` and `neighborhood` behave as in
[`mapneighbors`](@ref).
"""
function foreachneighbors(f::F, A::DD.AbstractDimArray; spatialdim = nothing,
        needs = nothing, pass = Neighbors(), order = StorageOrder(),
        threaded = false, connectivity::Connectivity = Vertex(),
        neighborhood = Disc(1)) where {F}
    dnum = _cells_dimnum(A, spatialdim)
    nb = _checkneighborhood(neighborhood)
    _foreach_needs(needs, pass, f, A, dnum, order, threaded, connectivity, nb)
    return nothing
end

_foreach_needs(::Nothing, pass, f::F, A, dnum, order, threaded, conn,
        nb) where {F} =
    _foreach_dimarray(pass, f, A, dnum, order, threaded, conn, nb)

function _foreach_needs(needs, pass, f::F, A, dnum, order, threaded,
        conn, nb) where {F}
    _checkpass(pass)
    plan = _chunkedvalues(A, dnum, order, conn, nb)
    plan === nothing ||
        return _foreach_needs_chunked(f, A, dnum, plan, needs, threaded, conn, nb)
    foreachneighbors(f, parent(DD.lookup(A, dnum)); needs, order, threaded,
        connectivity = conn, neighborhood = nb)
    return nothing
end

_foreach_dimarray(::Neighbors, f::F, A, dnum, order, threaded, conn,
        nb) where {F} =
    foreachneighbors(f, parent(DD.lookup(A, dnum)); order, threaded,
        connectivity = conn, neighborhood = nb)

function _foreach_dimarray(::Values, f::F, A, dnum, order, threaded,
        conn, nb) where {F}
    cv = parent(DD.lookup(A, dnum))
    data = parent(A)
    if ndims(A) == 1
        foreachneighbors(f, cv, data; order, threaded, connectivity = conn,
            neighborhood = nb)
        return nothing
    end
    pre = CartesianIndices(axes(data)[1:(dnum-1)])
    post = CartesianIndices(axes(data)[(dnum+1):end])
    buf = Vector{eltype(A)}(undef, size(data, dnum))
    for jpost in post, jpre in pre
        copyto!(buf, view(data, jpre, :, jpost))
        foreachneighbors(f, cv, buf; order, threaded, connectivity = conn,
            neighborhood = nb)
    end
    return nothing
end

function _foreach_dimarray(::NeighborSlices, f::F, A, dnum, order, threaded,
        conn, nb) where {F}
    _need_slices(A)
    return _foreach_cell_slices(f, A, Val(dnum), parent(DD.lookup(A, dnum)),
        order, threaded, conn, nb)
end

_foreach_dimarray(pass, f, A, dnum, order, threaded, conn, nb) = _bad_pass(pass)

function _foreach_cell_slices(f::F, A, ::Val{D}, cv, order, threaded,
        conn, nb) where {F,D}
    g = (c, nbrs) -> (f(c, _handle_slice(A, Val(D), localindex(c)),
        [_handle_slice(A, Val(D), localindex(h)) for h in nbrs]); nothing)
    foreachneighbors(g, cv; order, threaded, connectivity = conn,
        neighborhood = nb)
    return nothing
end

"""
    neighbors(A::AbstractDimArray; connectivity = Vertex(), neighborhood = Disc(1))
    neighbors(A::AbstractDimArray, dims; connectivity = Vertex(), neighborhood = Disc(1))

Iterate over each cell and its indexed neighbour handles; the minted indices
are the cell dimension's axis indices, and `neighborhood` is
[`mapneighbors`](@ref)'s.

`dims` names the cell dimension the way `DimensionalData.dims(A, dims)` does —
`Cells`, `:Cells`, or a `Dimension` (any `DimensionalData` dimension
selector, typed as such so a cell handle is never read as one) — and must name one that carries a cell
lookup. Without it, the first dimension carrying one is used. An array without
one, or a `dims` that misses or names a non-cell dimension, is an
`ArgumentError`.
"""
neighbors(A::DD.AbstractDimArray; connectivity::Connectivity = Vertex(),
        neighborhood = Disc(1)) =
    neighbors(parent(DD.lookup(A, _cells_dimnum(A, nothing))); connectivity,
        neighborhood = _checkneighborhood(neighborhood))

neighbors(A::DD.AbstractDimArray, dims::DimSelector;
        connectivity::Connectivity = Vertex(), neighborhood = Disc(1)) =
    neighbors(parent(DD.lookup(A, _cells_dimnum(A, dims))); connectivity,
        neighborhood = _checkneighborhood(neighborhood))

"""
    adjacency(A::AbstractDimArray; kw...) -> AdjacencyTable
    adjacency(A::AbstractDimArray, dims; kw...) -> AdjacencyTable

The one-ring [`adjacency`](@ref) table of `A`'s cell dimension, chosen as in
[`neighbors`](@ref): the first dimension carrying a cell lookup, or the one
`dims` names. The keywords are those of the [`CellVector`](@ref) method.
"""
adjacency(A::DD.AbstractDimArray; kw...) =
    adjacency(parent(DD.lookup(A, _cells_dimnum(A, nothing))); kw...)

adjacency(A::DD.AbstractDimArray, dims::DimSelector; kw...) =
    adjacency(parent(DD.lookup(A, _cells_dimnum(A, dims))); kw...)

"""
    mapneighbors(f, A::AbstractDimArray, dims; kw...)
    foreachneighbors(f, A::AbstractDimArray, dims; kw...)

The positional spelling of `spatialdim = dims`: `dims` names the cell
dimension as [`neighbors`](@ref) reads it. Every other keyword is as in the
two-argument forms.
"""
mapneighbors(f::F, A::DD.AbstractDimArray, dims::DimSelector;
        kw...) where {F} =
    mapneighbors(f, A; spatialdim = dims, kw...)

foreachneighbors(f::F, A::DD.AbstractDimArray, dims::DimSelector;
        kw...) where {F} =
    foreachneighbors(f, A; spatialdim = dims, kw...)

# Selector value types disambiguate DimensionalData's tuple and vector methods.

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
the coverage's leaf expansion with the lookup, in ascending index order. The
resulting view retains a [`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup).

Outside a cube, the equivalent selection is `covering(cv, target)`, which
returns a
[`CellVector`](@ref), or `covering_indices(cv, target)` for the indices
alone. See that docstring for what the selection costs and for the
over-covering it inherits from the coverage itself.

A [`query`](@ref) predicate — `Cells(Intersects(target))`,
`Cells(Within(target))` — is the exact selection the coverage over-covers;
see [`Cells`](@ref DiscreteGlobalGrids.CellLookups.Cells).
"""
struct Covering{T} <: Lookups.ArraySelector{T}
    val::T
end

Base.show(io::IO, sel::Covering) =
    print(io, "Covering(", typeof(sel.val).name.name, ")")

Lookups.selectindices(lk::AbstractCellAxis, sel::Covering; kw...) =
    covering_indices(parent(lk), Lookups.val(sel))

Lookups.selectindices(lk::AbstractCellAxis, sel::Covering{<:AbstractVector};
    kw...) = covering_indices(parent(lk), Lookups.val(sel))

"""
    Cells(pred::DE9IMPredicate)
    Cells(pred::CentroidCovered)

A [`query`](@ref) predicate is itself a cell selector: it selects every stored
cell that satisfies the predicate against its target, at the lookup's level.

```julia
A[Cells(Intersects(cap))]           # cells meeting a SphericalCap
A[Cells(Within(county))]            # cells lying wholly inside a polygon
A[Cells(Disjoint(extent))]          # cells clear of a lon/lat extent
A[Cells(CentroidCovered(county))]   # cells whose centroid lies in a polygon
```

The predicate and target are whatever `query` accepts — the same limits apply,
so a cap target supports `Intersects`, `Disjoint` and `Within` only. The
result is the query's answer intersected with the lookup, in ascending index
order; unlike [`Covering`](@ref) it is exact, not a coverage's over-cover. The
resulting view retains a [`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup).

Outside a cube, the equivalent selection is `predicate_indices(cv, pred)`.
"""
Lookups.selectindices(lk::AbstractCellLookup, pred::QueryPredicate; kw...) =
    predicate_indices(parent(lk), pred)

Lookups.selectindices(lk::AbstractCellAxis,
    sel::Lookups.At{<:AbstractCellIndex}; kw...) =
    _found(lk, localindex(lk, Lookups.val(sel)), sel)

Lookups.selectindices(lk::AbstractCellLookup,
    sel::Lookups.Contains{<:AbstractCellIndex}; kw...) =
    _found(lk, localindex(lk, Lookups.val(sel)), sel)

# Resolve a point to a leaf-grid cell, then search its index in the windows.
Lookups.selectindices(lk::AbstractCellAxis,
    sel::Lookups.Contains{<:Tuple{Real,Real}}; kw...) =
    _found(lk, localindex(parent(lk), Lookups.val(sel)...), sel)

Lookups.selectindices(lk::AbstractCellAxis,
    sel::Lookups.At{<:Tuple{Real,Real}}; kw...) =
    _found(lk, localindex(parent(lk), Lookups.val(sel)...), sel)

_found(lk::Lookups.Lookup, k::Int, sel) = k
_found(lk::Lookups.Lookup, ::Nothing, sel) = throw(Lookups.SelectorError(lk, sel))

# --- mixed-level lookup -----------------------------------------------------

"""
    MultiOrderLookup(mov::MultiOrderVector)
    MultiOrderLookup(set::MultiOrderCellSet)

Present a [`MultiOrderVector`](@ref) as a mixed-level DimensionalData lookup.
`Base.parent` returns the container. Pair the lookup with [`Cells`](@ref):

```julia
DimArray(vals, Cells(MultiOrderLookup(mov)))
```

Selectors use the container's interval index in O(log n): `At` tests exact
membership, `Contains` resolves the stored ancestor, and [`Covering`](@ref)
selects stored cells meeting a region. Ascending subsets retain the lookup;
reordered subsets use `Categorical`.
"""
struct MultiOrderLookup{ID,M<:MultiOrderVector} <: AbstractCellAxis{ID}
    cells::M
end

MultiOrderLookup(mov::MultiOrderVector{ID}) where {ID} =
    MultiOrderLookup{ID,typeof(mov)}(mov)

MultiOrderLookup(set::MultiOrderCellSet) = MultiOrderLookup(MultiOrderVector(set))

MultiOrderLookup(lk::MultiOrderLookup) = lk

# --- the collection surface ------------------------------------------------

# DimensionalData treats the parent as the lookup's values.
Base.parent(lk::MultiOrderLookup) = lk.cells

# Ascending subsets preserve disjoint interval order.
function _subset(lk::MultiOrderLookup, idx)
    sub = parent(lk)[idx]
    sub isa MultiOrderVector && return MultiOrderLookup(sub)
    return Lookups.Categorical(sub; order=Lookups.Unordered())
end

# --- metadata ---------------------------------------------------------------

"""
    reference_level(lk::MultiOrderLookup) -> Int

Return the backing container's interval level, which is at least as deep as
every cell on the axis. This function is public but unexported; call it as
`DiscreteGlobalGrids.reference_level(lk)`.
"""
reference_level(lk::MultiOrderLookup) = reference_level(parent(lk))

"""
    covering_index(lk::MultiOrderLookup, c::AbstractCellIndex) -> Union{Int,Nothing}

Return the index of the stored cell equal to or ancestral to `c`.
`DimensionalData.Contains` uses this lookup, including for cells deeper than
the [`reference_level`](@ref).
"""
covering_index(lk::MultiOrderLookup, c::AbstractCellIndex) =
    covering_index(parent(lk), c)

# --- DimensionalData plumbing ----------------------------------------------

# Cell-id ordering compares levels first, while selectors use interval order.
Lookups.order(::MultiOrderLookup) = Lookups.Unordered()

_rebuild(lk::MultiOrderLookup, mov::MultiOrderVector) = MultiOrderLookup(mov)

# Concatenation preserves this lookup only for strictly ascending ids.
function _rebuild(lk::MultiOrderLookup, ids::AbstractVector{<:AbstractCellIndex})
    mov = parent(lk)
    sys = system(mov)
    # Key the concatenated axis at its deepest cell level.
    ref = maximum(level, ids; init=reference_level(mov))
    starts = [first(descendant_range(sys, c, ref)) for c in ids]
    issorted(starts; lt=<=) ||
        return Lookups.Categorical(collect(ids); order=Lookups.Unordered())
    return MultiOrderLookup(MultiOrderVector(sys, ids; reference_level=ref))
end

@noinline _rebuild(lk::MultiOrderLookup, data) = throw(ArgumentError(
    "MultiOrderLookup rebuild data must be mixed-level cell ids, got " *
    "$(typeof(data)). Use `vcat` or `cat` for cell axes, indexing for subsets, " *
    "or `set(A, Cells => NoLookup())` to replace the lookup."))

# The interval index gives each selector one binary search.
Lookups.hasselection(lk::MultiOrderLookup, sel::Lookups.Contains{<:AbstractCellIndex}) =
    covering_index(lk, Lookups.val(sel)) !== nothing

Base.:(==)(a::MultiOrderLookup, b::MultiOrderLookup) = parent(a) == parent(b)

function Base.show(io::IO, lk::MultiOrderLookup)
    mov = parent(lk)
    print(io, "MultiOrderLookup(", typeof(system(lk)).name.name, ", ncells=",
        length(lk))
    isempty(lk) || print(io, ", levels=", minimum(level, mov), ":", maximum(level, mov))
    print(io, ", ref=", reference_level(lk), ")")
end

Base.show(io::IO, ::MIME"text/plain", lk::MultiOrderLookup) = show(io, lk)

# --- selectors -------------------------------------------------------------

Lookups.selectindices(lk::MultiOrderLookup,
    sel::Lookups.Contains{<:AbstractCellIndex}; kw...) =
    _found(lk, covering_index(parent(lk), Lookups.val(sel)), sel)

# --- lazy expansion values -------------------------------------------------

"""
    MultiOrderValues(values, offsets) <: AbstractVector

Present mixed-level values over their covered leaf cells. `offsets[i]` is the
cumulative leaf count for the first `i` stored cells; indexing finds the owning
stored value by binary search. Logical length counts leaves while storage
remains O(stored values).
"""
struct MultiOrderValues{T,V<:AbstractVector{T}} <: AbstractVector{T}
    values::V
    offsets::Vector{Int}

    function MultiOrderValues{T,V}(values::V, offsets::Vector{Int}) where {T,V}
        length(values) == length(offsets) || throw(ArgumentError(
            "a multi-order value vector needs one leaf count per stored value, " *
            "got $(length(values)) values against $(length(offsets)) counts"))
        return new{T,V}(values, offsets)
    end
end

MultiOrderValues(values::AbstractVector{T}, offsets::Vector{Int}) where {T} =
    MultiOrderValues{T,typeof(values)}(values, offsets)

Base.size(v::MultiOrderValues) = (isempty(v.offsets) ? 0 : @inbounds(v.offsets[end]),)
Base.IndexStyle(::Type{<:MultiOrderValues}) = Base.IndexLinear()

Base.@propagate_inbounds function Base.getindex(v::MultiOrderValues, k::Int)
    @boundscheck checkbounds(v, k)
    return @inbounds v.values[searchsortedfirst(v.offsets, k)]
end

function Base.copyto!(dest::Vector, v::MultiOrderValues)
    length(dest) >= length(v) || throw(ArgumentError(
        "destination has $(length(dest)) elements, cannot hold $(length(v))"))
    lo = 1
    @inbounds for i in eachindex(v.offsets)
        hi = v.offsets[i]
        x = v.values[i]
        for k in lo:hi
            dest[k] = x
        end
        lo = hi + 1
    end
    return dest
end

Base.collect(v::MultiOrderValues{T}) where {T} =
    copyto!(Vector{T}(undef, length(v)), v)

# --- DimArray aggregation methods ------------------------------------------

function _cell_axis(A::DD.AbstractDimArray, ::Type{L}, verb::AbstractString) where {L}
    ndims(A) == 1 || throw(ArgumentError(
        "$verb is defined on a one-dimensional array over a `Cells` axis; this " *
        "one has $(ndims(A)) dimensions, $(DD.dims(A))"))
    DD.hasdim(A, Cells) || throw(ArgumentError(
        "$verb needs a `Cells` axis; this array has $(DD.dims(A))"))
    lk = DD.lookup(A, Cells)
    lk isa L || throw(ArgumentError(
        "$verb needs a `Cells` axis whose lookup is a $(nameof(L)); got a " *
        "$(nameof(typeof(lk)))"))
    return lk
end

"""
    aggregate(f, A::AbstractDimArray, l::Integer) -> AbstractDimArray

Reduce a one-dimensional cell array to level `l`. Each level-`l` ancestor
contributes `f` of its present descendant values. The result carries a coarser
[`CellLookup`](@ref).

Build a pyramid with one call per level:

```julia
pyramid = [aggregate(sum, A, l) for l in level(lookup(A, Cells)) - 1 : -1 : 3]
```

See the `(CellVector, values)` method for what `f` sees and how partial groups
are treated.
"""
function aggregate(f, A::DD.AbstractDimArray, l::Integer)
    lk = _cell_axis(A, AbstractCellLookup, "aggregate")
    coarse, vals = aggregate(f, DGG.region(lk), parent(A), l)
    return _aggregated(A, coarse, vals)
end

# A function barrier preserves the two concrete `CellVector` window shapes.
_aggregated(A::DD.AbstractDimArray, coarse::CellVector, vals::AbstractVector) =
    DD.rebuild(A; data=vals, dims=(Cells(CellLookup(coarse)),))

"""
    coarsen(A::AbstractDimArray; atol, by = mean, minlevel = shallowest) -> AbstractDimArray

Build an adaptive one-dimensional cell array by merging complete subtrees whose
values agree within `atol`. The result carries a [`MultiOrderLookup`](@ref).

```julia
M = coarsen(A; atol = 1.0)                       # °C
M[Cells(DimensionalData.Contains(8.0, 46.5))]    # a point, on the mesh
expand(M, level(lookup(A, Cells)))               # back, within `atol`
```

See the `(CellVector, values)` method for the merge criterion, the treatment
of `missing`, and the error bound the default `by` carries.
"""
function coarsen(A::DD.AbstractDimArray; atol, kw...)
    lk = _cell_axis(A, AbstractCellLookup, "coarsen")
    mov, vals = coarsen(DGG.region(lk), parent(A); atol, kw...)
    return DD.rebuild(A; data=vals, dims=(Cells(MultiOrderLookup(mov)),))
end

"""
    expand(A::AbstractDimArray, l::Integer) -> AbstractDimArray

Present a one-dimensional mixed-level array at level `l`. The result carries a
[`CellLookup`](@ref) and lazy [`MultiOrderValues`](@ref), retaining one stored
value per mixed-level cell. `l` must cover every stored cell level.

```julia
E = expand(coarsen(A; atol), level(lookup(A, Cells)))
all(abs.(collect(parent(E)) .- parent(A)) .<= atol)   # the round-trip bound
```
"""
function expand(A::DD.AbstractDimArray, l::Integer)
    lk = _cell_axis(A, MultiOrderLookup, "expand")
    cv, data = expand(parent(lk), parent(A), l)
    return DD.rebuild(A; data=data, dims=(Cells(CellLookup(cv)),))
end

"""
    expand(mov::MultiOrderVector, values::AbstractVector, l::Integer) -> (CellVector, MultiOrderValues)

Expand a `(cells, values)` pair returned by [`coarsen`](@ref). `values` is
indexed against `mov`; the lazy result retains one value per stored cell.
"""
function expand(mov::MultiOrderVector, values::AbstractVector, l::Integer)
    length(values) == length(mov) || throw(ArgumentError(
        "values must line up with the cells they belong to: got $(length(values)) " *
        "values for $(length(mov)) cells"))
    target = Int(l)
    return CellVector(mov; level=target),
        MultiOrderValues(values, _leaf_offsets(mov, target))
end

function _leaf_offsets(mov::MultiOrderVector, l::Int)
    sys = system(mov)
    offsets = Vector{Int}(undef, length(mov))
    total = 0
    for i in eachindex(offsets)
        total += length(descendant_range(sys, @inbounds(mov[i]), l))
        @inbounds offsets[i] = total
    end
    return offsets
end

# Cell-id order cannot represent spherical nearest-neighbor distance.
@noinline _no_near(lk, sel) = throw(ArgumentError(
    "Cell-id order follows a space-filling curve: the nearest id is not the " *
    "nearest cell on the sphere. Use " *
    "At(cell) for one cell, Contains(lon, lat) for the cell holding a point, " *
    "or Covering(region) for a region."))

Lookups.selectindices(lk::AbstractCellAxis, sel::Lookups.Near; kw...) =
    _no_near(lk, sel)
Lookups.selectindices(lk::AbstractCellAxis, sel::Lookups.Near{<:AbstractVector};
    kw...) = _no_near(lk, sel)

# Cell-specific selector failures expose likely level or system mismatches.

"""
    show_selector_error(io, lk, sel)

Print a failed cell-axis selection and any level or system mismatch.
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
