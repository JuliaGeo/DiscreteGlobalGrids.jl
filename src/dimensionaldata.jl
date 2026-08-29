"""
    CellLookups

DimensionalData integration for discrete global grid cells.

  - [`CellLookup`](@ref) wraps a same-level [`CellVector`](@ref).
  - [`MultiOrderLookup`](@ref) wraps a mixed-level
    [`MultiOrderVector`](@ref).
  - [`Cells`](@ref) names the cell dimension.
  - [`Covering`](@ref) selects cells through a spatial coverage query.
"""
module CellLookups

import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, AbstractCellVector,
    ncells, cellindex, localindex, globalindex, cellat, level, system,
    levelgrid, cellindextype, has_sorted_subtrees, descendants, descendant_range,
    query, neighbors, ring, neighborcount, Connectivity, Vertex, maxneighbors,
    halo, border, interior, adjacency
import ..DiscreteGlobalGrids: Helpers
import ..DiscreteGlobalGrids.Engine: PartialGrid, SubtreeIds,
    MultiOrderCoverage, MultiOrderCellSet, level_ranges
import ..DiscreteGlobalGrids.Engine: CellVector, cellset, covering,
    covering_indices, windows, nwindows, RangeWindows, CellWindows, _derive,
    _windows, SubsetIndexedCell, mapneighbors, foreachneighbors,
    StorageOrder, _capacity, _ringtype
import ..DiscreteGlobalGrids.Engine: MultiOrderVector, reference_level,
    covering_index, aggregate, coarsen, expand

import SmallCollections
import DimensionalData as DD
import DimensionalData: Dimensions, Lookups

# --- same-level lookup ------------------------------------------------------

"""
    abstract type AbstractCellLookup{ID} <: DimensionalData.Lookups.Lookup{ID,1}

Abstract supertype for one-dimensional same-level cell lookups. `Base.parent`
returns an [`AbstractCellVector`](@ref) that supplies collection, topology, and
region operations. [`CellLookup`](@ref) uses computed windows;
[`ChunkedCellLookup`](@ref) uses a stored axis. Subtypes define rebuild and
selector behavior for their backing.
"""
abstract type AbstractCellLookup{ID} <: Lookups.Lookup{ID,1} end

"""
    CellLookup(cv::CellVector)
    CellLookup(set::MultiOrderCellSet; level = set's reference level)
    CellLookup(grid::AbstractGrid)

Create a same-level DimensionalData lookup backed by a compressed
[`CellVector`](@ref). Pair it with [`Cells`](@ref) to make a cube axis:

```julia
set = query(sys, MultiOrderCoverage(region); level = 9)
lk  = CellLookup(set)
A   = DimensionalData.DimArray(values, Cells(lk))
```

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

`DimensionalData.At` selects an exact id, `DimensionalData.Contains` selects
the cell containing a point, and [`Covering`](@ref) selects a spatial region.
`Near` is unavailable because id order is not spherical distance.

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
Base.IndexStyle(::Type{<:AbstractCellLookup}) = Base.IndexLinear()

# Empty lookups have no lower or upper value bound.
Lookups.bounds(lk::AbstractCellLookup) =
    isempty(lk) ? (nothing, nothing) : (first(lk), last(lk))

Base.@propagate_inbounds Base.getindex(lk::AbstractCellLookup, k::Int) = parent(lk)[k]
Base.@propagate_inbounds Base.getindex(lk::AbstractCellLookup, k::CartesianIndex{1}) =
    parent(lk)[k[1]]

# Base handles shaped indices; vector indices may preserve compressed windows.
for f in (:getindex, :view, :dotview)
    @eval Base.$f(lk::AbstractCellLookup, ::Colon) = lk
    @eval Base.$f(lk::AbstractCellLookup, i::AbstractVector{<:Integer}) = _subset(lk, i)
end

# Route reversal through lookup indexing so the result remains a valid lookup.
Base.reverse(lk::AbstractCellLookup) = lk[lastindex(lk):-1:firstindex(lk)]

# Resolve the method ambiguity with SmallCollections vector indexing.
Base.getindex(lk::AbstractCellLookup,
    i::SmallCollections.AbstractFixedOrSmallOrPackedVector{<:Integer}) = _subset(lk, i)

# Concrete mask methods avoid the `Bool <: Integer` catch-all ambiguity.
function _masksubset(lk::AbstractCellLookup, mask::AbstractArray{Bool})
    axes(mask) == axes(lk) || throw(BoundsError(lk, (mask,)))
    return _subset(lk, findall(mask))
end

_subset(lk::CellLookup, mask::AbstractArray{Bool}) = _masksubset(lk, mask)

# Only ascending subsets retain compressed interval order.
function _subset(lk::CellLookup, idx)
    sub = parent(lk)[idx]
    sub isa CellVector && return CellLookup(sub)
    return Lookups.Categorical(sub; order=Lookups.Unordered())
end

# --- what the lookup is, in this package's own vocabulary ------------------

"""
    cellset(lk::AbstractCellLookup)

Return the set or grid that constructed the lookup. Derived and stored-axis
lookups return their describing [`PartialGrid`](@ref).
"""
cellset(lk::AbstractCellLookup) = cellset(parent(lk))

"""
    system(lk::AbstractCellLookup)

The grid system the lookup's cells are named in.
"""
system(lk::AbstractCellLookup) = system(parent(lk))

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

Index of cell `c` in the lookup, or `nothing` when the lookup does not hold
it — including when `c` is at another level. The inverse of `lk[k]`, and the
half of the bijection every selector ends at.
"""
localindex(lk::AbstractCellLookup, c::AbstractCellIndex) = localindex(parent(lk), c)
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

neighbors(lk::AbstractCellLookup; connectivity::Connectivity=Vertex()) =
    neighbors(parent(lk); connectivity)

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
Lookups.metadata(::AbstractCellLookup) = Lookups.NoMetadata()

# Rebuild preserves compression only for ascending ids.
function Lookups.rebuild(lk::CellLookup; data=nothing, kw...)
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
Lookups.reducelookup(::AbstractCellLookup) = Lookups.NoLookup(Base.OneTo(1))

# Window membership keeps selector checks logarithmic.
Lookups.hasselection(lk::AbstractCellLookup, sel::Lookups.At{<:AbstractCellIndex}) =
    localindex(lk, Lookups.val(sel)) !== nothing

Lookups.hasselection(lk::AbstractCellLookup,
    sel::Lookups.Contains{<:AbstractCellIndex}) =
    localindex(lk, Lookups.val(sel)) !== nothing

Lookups.hasselection(lk::AbstractCellLookup,
    sel::Lookups.Contains{<:Tuple{Real,Real}}) =
    localindex(parent(lk), Lookups.val(sel)...) !== nothing

Dimensions.format(lk::AbstractCellLookup, ::Type, values, axis::AbstractRange) = lk

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

The DimensionalData dimension for a cell axis. Use `Cells(lk)` to construct an
axis and `Cells(selector)` to index it.

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
                 order = StorageOrder(), threaded = true, connectivity = Vertex())
    mapneighbors(f, A::AbstractDimArray; needs = (Value(a), Centroid()), ...)

Apply `f` to each cell and its neighbors. The result uses `A`'s wrapper and
lookups. `spatialdim` selects the cell dimension and defaults to the first cell
lookup.

`pass` controls the callback arguments and output shape:

- [`Neighbors`](@ref): indexed handles, one result per cell.
- [`Values`](@ref): scalar values; the same dimensions as `A`.
- [`NeighborSlices`](@ref): views across non-cell dimensions, one result per
  cell.

Alternatively, `needs` supplies `Cell`, `Index`, `Value`, and `Centroid` field
requests to `f(center, rings)`. Field requests produce one result per cell and
use the default pass mode.

Stored values follow chunk order when possible; see [`chunkplan`](@ref). An
explicit permutation takes precedence over chunk traversal.
"""
function mapneighbors(f::F, A::DD.AbstractDimArray; spatialdim = nothing,
        needs = nothing, pass = Neighbors(), order = StorageOrder(),
        threaded = true, connectivity::Connectivity = Vertex()) where {F}
    dnum = _cells_dimnum(A, spatialdim)
    return _map_needs(needs, pass, f, A, dnum, order, threaded, connectivity)
end

_map_needs(::Nothing, pass, f::F, A, dnum, order, threaded, conn) where {F} =
    _map_dimarray(pass, f, A, dnum, order, threaded, conn)

# Field requests produce one cell-axis result and may use chunk traversal.
function _map_needs(needs, pass, f::F, A, dnum, order, threaded, conn) where {F}
    _checkpass(pass)
    plan = _chunkedvalues(A, dnum, order, conn)
    plan === nothing || return _rebuilt_on_cells(A, DD.dims(A)[dnum],
        _map_needs_chunked(f, A, dnum, plan, needs, threaded, conn))
    cv = parent(DD.lookup(A, dnum))
    out = mapneighbors(f, cv; needs, order, threaded, connectivity = conn)
    return _rebuilt_on_cells(A, DD.dims(A)[dnum], out)
end

function _map_needs_chunked end
function _foreach_needs_chunked end

function _map_dimarray(::Neighbors, f::F, A, dnum, order, threaded,
        conn) where {F}
    cv = parent(DD.lookup(A, dnum))
    out = mapneighbors(f, cv; order, threaded, connectivity = conn)
    return _rebuilt_on_cells(A, DD.dims(A)[dnum], out)
end

function _map_dimarray(::Values, f::F, A, dnum, order, threaded,
        conn) where {F}
    plan = _chunkedvalues(A, dnum, order, conn)
    plan === nothing || return _rebuilt(A,
        _map_values_chunked(f, A, dnum, plan, threaded, conn))
    cv = parent(DD.lookup(A, dnum))
    ndims(A) == 1 && return _rebuilt(A,
        mapneighbors(f, cv, parent(A); order, threaded, connectivity = conn))
    return _rebuilt(A, _map_slices(f, A, dnum, cv, order, threaded, conn))
end

# Value traversal can follow chunk storage; indexed handles remain whole-axis indices.
_chunkedvalues(A, dnum, order, conn) = nothing
function _map_values_chunked end

function _map_dimarray(::NeighborSlices, f::F, A, dnum, order, threaded,
        conn) where {F}
    _need_slices(A)
    cv = parent(DD.lookup(A, dnum))
    out = _map_cell_slices(f, A, Val(dnum), cv, order, threaded, conn)
    return _rebuilt_on_cells(A, DD.dims(A)[dnum], out)
end

_map_dimarray(pass, f, A, dnum, order, threaded, conn) = _bad_pass(pass)

# A value-type dimension number keeps slice views concrete.
function _map_cell_slices(f::F, A, ::Val{D}, cv, order, threaded,
        conn) where {F,D}
    g = (c, nbrs) -> f(c, _handle_slice(A, Val(D), localindex(c)),
        [_handle_slice(A, Val(D), localindex(h)) for h in nbrs])
    return mapneighbors(g, cv; order, threaded, connectivity = conn)
end

# Independent buffered sweeps isolate each non-cell slice.
function _map_slices(f::F, A, dnum::Int, cv::AbstractCellVector, order, threaded,
        connectivity::Connectivity) where {F}
    data = parent(A)
    pre = CartesianIndices(axes(data)[1:(dnum-1)])
    post = CartesianIndices(axes(data)[(dnum+1):end])
    cap = _capacity(system(cv), connectivity)
    H = SubsetIndexedCell{eltype(cv)}
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
    foreachneighbors(f, A::AbstractDimArray; needs = (Value(a), Centroid()), ...)

Call `f` for side effects on each cell and its neighbors.
`spatialdim`, `pass` and `needs` behave as in [`mapneighbors`](@ref).
"""
function foreachneighbors(f::F, A::DD.AbstractDimArray; spatialdim = nothing,
        needs = nothing, pass = Neighbors(), order = StorageOrder(),
        threaded = false, connectivity::Connectivity = Vertex()) where {F}
    dnum = _cells_dimnum(A, spatialdim)
    _foreach_needs(needs, pass, f, A, dnum, order, threaded, connectivity)
    return nothing
end

_foreach_needs(::Nothing, pass, f::F, A, dnum, order, threaded, conn) where {F} =
    _foreach_dimarray(pass, f, A, dnum, order, threaded, conn)

function _foreach_needs(needs, pass, f::F, A, dnum, order, threaded,
        conn) where {F}
    _checkpass(pass)
    plan = _chunkedvalues(A, dnum, order, conn)
    plan === nothing ||
        return _foreach_needs_chunked(f, A, dnum, plan, needs, threaded, conn)
    foreachneighbors(f, parent(DD.lookup(A, dnum)); needs, order, threaded,
        connectivity = conn)
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
    g = (c, nbrs) -> (f(c, _handle_slice(A, Val(D), localindex(c)),
        [_handle_slice(A, Val(D), localindex(h)) for h in nbrs]); nothing)
    foreachneighbors(g, cv; order, threaded, connectivity = conn)
    return nothing
end

"""
    neighbors(A::AbstractDimArray; spatialdim = nothing, connectivity = Vertex())

Iterate over each cell and its indexed neighbor handles. The cell
dimension is selected as in [`mapneighbors`](@ref); the minted indices are
that dimension's axis indices.
"""
function neighbors(A::DD.AbstractDimArray; spatialdim = nothing,
        connectivity::Connectivity = Vertex())
    dnum = _cells_dimnum(A, spatialdim)
    return neighbors(parent(DD.lookup(A, dnum)); connectivity)
end

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

`target` accepts the same inputs as [`query`](@ref). The ascending intersection
retains a [`CellLookup`](@ref). Use `covering(cv, target)` for a
[`CellVector`](@ref), or `covering_indices(cv, target)` for indices.
"""
struct Covering{T} <: Lookups.ArraySelector{T}
    val::T
end

Base.show(io::IO, sel::Covering) =
    print(io, "Covering(", typeof(sel.val).name.name, ")")

Lookups.selectindices(lk::AbstractCellLookup, sel::Covering; kw...) =
    covering_indices(parent(lk), Lookups.val(sel))

Lookups.selectindices(lk::AbstractCellLookup, sel::Covering{<:AbstractVector};
    kw...) = covering_indices(parent(lk), Lookups.val(sel))

Lookups.selectindices(lk::AbstractCellLookup,
    sel::Lookups.At{<:AbstractCellIndex}; kw...) =
    _found(lk, localindex(lk, Lookups.val(sel)), sel)

Lookups.selectindices(lk::AbstractCellLookup,
    sel::Lookups.Contains{<:AbstractCellIndex}; kw...) =
    _found(lk, localindex(lk, Lookups.val(sel)), sel)

# Resolve a point to a leaf-grid cell, then search its index in the windows.
Lookups.selectindices(lk::AbstractCellLookup,
    sel::Lookups.Contains{<:Tuple{Real,Real}}; kw...) =
    _found(lk, localindex(parent(lk), Lookups.val(sel)...), sel)

Lookups.selectindices(lk::AbstractCellLookup,
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
struct MultiOrderLookup{ID,M<:MultiOrderVector} <: Lookups.Lookup{ID,1}
    cells::M
end

MultiOrderLookup(mov::MultiOrderVector{ID}) where {ID} =
    MultiOrderLookup{ID,typeof(mov)}(mov)

MultiOrderLookup(set::MultiOrderCellSet) = MultiOrderLookup(MultiOrderVector(set))

MultiOrderLookup(lk::MultiOrderLookup) = lk

# --- the collection surface ------------------------------------------------

# DimensionalData treats the parent as the lookup's values.
Base.parent(lk::MultiOrderLookup) = lk.cells
Base.IndexStyle(::Type{<:MultiOrderLookup}) = Base.IndexLinear()

"""
    cellset(lk::MultiOrderLookup)

Return the [`MultiOrderVector`](@ref) backing the axis. This is the same
collection as `Base.parent(lk)`.
"""
cellset(lk::MultiOrderLookup) = cellset(parent(lk))

Base.@propagate_inbounds Base.getindex(lk::MultiOrderLookup, k::Int) = parent(lk)[k]
Base.@propagate_inbounds Base.getindex(lk::MultiOrderLookup, k::CartesianIndex{1}) =
    parent(lk)[k[1]]

for f in (:getindex, :view, :dotview)
    @eval Base.$f(lk::MultiOrderLookup, ::Colon) = lk
    @eval Base.$f(lk::MultiOrderLookup, i::AbstractVector{<:Integer}) = _subset(lk, i)
end

Base.reverse(lk::MultiOrderLookup) = lk[lastindex(lk):-1:firstindex(lk)]

Base.getindex(lk::MultiOrderLookup,
    i::SmallCollections.AbstractFixedOrSmallOrPackedVector{<:Integer}) = _subset(lk, i)

function _subset(lk::MultiOrderLookup, mask::AbstractArray{Bool})
    axes(mask) == axes(lk) || throw(BoundsError(lk, (mask,)))
    return _subset(lk, findall(mask))
end

# Ascending subsets preserve disjoint interval order.
function _subset(lk::MultiOrderLookup, idx)
    sub = parent(lk)[idx]
    sub isa MultiOrderVector && return MultiOrderLookup(sub)
    return Lookups.Categorical(sub; order=Lookups.Unordered())
end

# --- metadata ---------------------------------------------------------------

"""
    system(lk::MultiOrderLookup)

The grid system the lookup's cells are named in.
"""
system(lk::MultiOrderLookup) = system(parent(lk))

"""
    reference_level(lk::MultiOrderLookup) -> Int

Return the backing container's interval level, which is at least as deep as
every cell on the axis. This function is public but unexported; call it as
`DiscreteGlobalGrids.reference_level(lk)`.
"""
reference_level(lk::MultiOrderLookup) = reference_level(parent(lk))

"""
    localindex(lk::MultiOrderLookup, c::AbstractCellIndex) -> Union{Int,Nothing}

Return the index of the exact cell `c`, or `nothing` when the axis stores only
an ancestor or descendant. `DimensionalData.At` uses this lookup.
"""
localindex(lk::MultiOrderLookup, c::AbstractCellIndex) = localindex(parent(lk), c)

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
Lookups.metadata(::MultiOrderLookup) = Lookups.NoMetadata()

function Lookups.rebuild(lk::MultiOrderLookup; data=nothing, kw...)
    (data === nothing || data === lk || data === parent(lk)) && return lk
    return _rebuild(lk, data)
end

_rebuild(lk::MultiOrderLookup, mov::MultiOrderVector) = MultiOrderLookup(mov)

# Concatenation preserves this lookup only for ascending, disjoint ids.
function _rebuild(lk::MultiOrderLookup, ids::AbstractVector{<:AbstractCellIndex})
    mov = parent(lk)
    sys = system(mov)
    # Key the concatenated axis at its deepest cell level.
    ref = reference_level(mov)
    for c in ids
        ref = max(ref, level(c))
    end
    ascending = true
    prev = 0
    for c in ids
        s = first(descendant_range(sys, c, ref))
        s > prev || (ascending = false)
        prev = s
    end
    ascending || return Lookups.Categorical(collect(ids); order=Lookups.Unordered())
    return MultiOrderLookup(MultiOrderVector(sys, ids; reference_level=ref))
end

@noinline _rebuild(lk::MultiOrderLookup, data) = throw(ArgumentError(
    "MultiOrderLookup rebuild data must be mixed-level cell ids, got " *
    "$(typeof(data)). Use `vcat` or `cat` for cell axes, indexing for subsets, " *
    "or `set(A, Cells => NoLookup())` to replace the lookup."))

Lookups.reducelookup(::MultiOrderLookup) = Lookups.NoLookup(Base.OneTo(1))

# The interval index gives each selector one binary search.
Lookups.hasselection(lk::MultiOrderLookup, sel::Lookups.At{<:AbstractCellIndex}) =
    localindex(lk, Lookups.val(sel)) !== nothing

Lookups.hasselection(lk::MultiOrderLookup, sel::Lookups.Contains{<:AbstractCellIndex}) =
    covering_index(lk, Lookups.val(sel)) !== nothing

Lookups.hasselection(lk::MultiOrderLookup, sel::Lookups.Contains{<:Tuple{Real,Real}}) =
    localindex(parent(lk), Lookups.val(sel)...) !== nothing

Dimensions.format(lk::MultiOrderLookup, ::Type, values, axis::AbstractRange) = lk

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

Lookups.selectindices(lk::MultiOrderLookup, sel::Covering; kw...) =
    covering_indices(parent(lk), Lookups.val(sel))

Lookups.selectindices(lk::MultiOrderLookup, sel::Covering{<:AbstractVector}; kw...) =
    covering_indices(parent(lk), Lookups.val(sel))

Lookups.selectindices(lk::MultiOrderLookup, sel::Lookups.At{<:AbstractCellIndex}; kw...) =
    _found(lk, localindex(parent(lk), Lookups.val(sel)), sel)

Lookups.selectindices(lk::MultiOrderLookup,
    sel::Lookups.Contains{<:AbstractCellIndex}; kw...) =
    _found(lk, covering_index(parent(lk), Lookups.val(sel)), sel)

Lookups.selectindices(lk::MultiOrderLookup,
    sel::Lookups.Contains{<:Tuple{Real,Real}}; kw...) =
    _found(lk, localindex(parent(lk), Lookups.val(sel)...), sel)

Lookups.selectindices(lk::MultiOrderLookup,
    sel::Lookups.At{<:Tuple{Real,Real}}; kw...) =
    _found(lk, localindex(parent(lk), Lookups.val(sel)...), sel)

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

Lookups.selectindices(lk::AbstractCellLookup, sel::Lookups.Near; kw...) =
    _no_near(lk, sel)
Lookups.selectindices(lk::AbstractCellLookup, sel::Lookups.Near{<:AbstractVector};
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
