# DimensionalData wrappers for `Fallbacks.CellVector`. Cell ids are computed
# from compressed position windows rather than stored as an expanded vector.

"""
    CellLookups

The DimensionalData layer: [`CellLookup`](@ref) and [`MultiOrderLookup`](@ref),
the [`Cells`](@ref) dimension, and the [`Covering`](@ref) selector.

A [`CellLookup`](@ref) is a one-dimensional `DimensionalData` lookup over cell
ids at a single level. It is a thin wrapper around a [`CellVector`](@ref),
which is where the compression lives: a set of **leaf position windows** —
sorted, disjoint intervals (or, where intervals are unavailable, a sorted list)
of positions in `levelgrid(sys, leaf)`. Its logical content is their
concatenation, and every operation is arithmetic over that concatenation:
`length` sums the window lengths, `lk[k]` binary-searches the cumulative
lengths and resolves one `cellindex`, [`cellposition`](@ref) runs the inverse.
Nothing is materialised.

A [`MultiOrderLookup`](@ref) is the mixed-level counterpart: a thin wrapper
around a [`MultiOrderVector`](@ref). [`coarsen`](@ref) builds one from a
`CellLookup` axis, [`expand`](@ref) presents one back at a single level, and
[`aggregate`](@ref) is the fixed-level reduction.

`CellVector` provides the storage and indexing behavior; this module provides
the lookup, dimension, and selectors required by DimensionalData.
"""
module CellLookups

import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, ncells, cellindex, cellposition, cellat, level, system,
    levelgrid, cellindextype, has_sorted_subtrees, descendants, descendant_range,
    query, neighbors, ring, halo_table, halo, neighborcount, Connectivity, Vertex,
    max_neighbors
import ..DiscreteGlobalGrids: Helpers
import ..DiscreteGlobalGrids.Fallbacks: PartialGrid, SubtreeIds,
    MultiOrderCoverage, MultiOrderCellSet, level_ranges
# Core collection operations delegated to `CellVector`.
import ..DiscreteGlobalGrids.Fallbacks: CellVector, cellset, covering,
    covering_positions, windows, nwindows, RangeWindows, CellWindows, _derive,
    _windows, SubsetPositionedCell, mapneighbors, foreachneighbors, HaloTable,
    StorageOrder
# The mixed-level container, its membership verbs, and the aggregation verbs
# whose DimArray methods close this file.
import ..DiscreteGlobalGrids.Fallbacks: MultiOrderVector, reference_level,
    covering_position, aggregate, coarsen, expand

import SmallCollections
import DimensionalData as DD
import DimensionalData: Dimensions, Lookups

# ===========================================================================
# The lookup
# ===========================================================================

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
struct CellLookup{ID,C<:CellVector} <: Lookups.Lookup{ID,1}
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
Base.IndexStyle(::Type{<:CellLookup}) = Base.IndexLinear()

# Empty lookups have no lower or upper value bound.
Lookups.bounds(lk::CellLookup) = isempty(lk) ? (nothing, nothing) : (first(lk), last(lk))

Base.@propagate_inbounds Base.getindex(lk::CellLookup, k::Int) = parent(lk)[k]
Base.@propagate_inbounds Base.getindex(lk::CellLookup, k::CartesianIndex{1}) = parent(lk)[k[1]]

# Only vector indices can produce a one-dimensional window set. Shaped indices
# use Base's generic array indexing.
for f in (:getindex, :view, :dotview)
    @eval Base.$f(lk::CellLookup, ::Colon) = lk
    @eval Base.$f(lk::CellLookup, i::AbstractVector{<:Integer}) = _subset(lk, i)
end

# Route reversal through lookup indexing so the result remains a valid lookup.
Base.reverse(lk::CellLookup) = lk[lastindex(lk):-1:firstindex(lk)]

# Resolve the method ambiguity with SmallCollections vector indexing.
Base.getindex(lk::CellLookup,
    i::SmallCollections.AbstractFixedOrSmallOrPackedVector{<:Integer}) = _subset(lk, i)

# Validate boolean-mask axes against the lookup so bounds errors identify the
# indexed lookup rather than its backing vector.
function _subset(lk::CellLookup, mask::AbstractArray{Bool})
    axes(mask) == axes(lk) || throw(BoundsError(lk, (mask,)))
    return _subset(lk, findall(mask))
end

# Ascending subsets remain compressed; reordered subsets use an unordered
# DimensionalData categorical lookup.
function _subset(lk::CellLookup, idx)
    sub = parent(lk)[idx]
    sub isa CellVector && return CellLookup(sub)
    return Lookups.Categorical(sub; order=Lookups.Unordered())
end

# --- what the lookup is, in this package's own vocabulary ------------------

"""
    cellset(lk::CellLookup)

Return the [`MultiOrderCellSet`](@ref) or grid used to construct the lookup.

`Base.parent(lk)` returns the logical values as a [`CellVector`](@ref). A subset
produced by indexing or selection reports its [`PartialGrid`](@ref).
"""
cellset(lk::CellLookup) = cellset(parent(lk))

"""
    system(lk::CellLookup)

The grid system the lookup's cells are named in.
"""
system(lk::CellLookup) = system(parent(lk))

"""
    level(lk::CellLookup) -> Int

The one level every cell in the lookup sits at.
"""
level(lk::CellLookup) = level(parent(lk))

"""
    cellposition(lk::CellLookup, c::AbstractCellIndex) -> Union{Int,Nothing}

Position of cell `c` in the lookup, or `nothing` when the lookup does not hold
it — including when `c` is at another level. The inverse of `lk[k]`, and the
half of the bijection every selector ends at.
"""
cellposition(lk::CellLookup, c::AbstractCellIndex) = cellposition(parent(lk), c)

"""
    neighbors(lk::CellLookup, c, k = 1; connectivity = Vertex())
    ring(lk::CellLookup, c, k; connectivity = Vertex())
    neighbors(lk::CellLookup, p::Int, k = 1; connectivity = Vertex()) -> Vector{Int}
    ring(lk::CellLookup, p::Int, k; connectivity = Vertex()) -> Vector{Int}
    halo_table(lk::CellLookup, k = 1; connectivity = Vertex()) -> Vector{Vector{Int}}
    halo(lk::CellLookup; connectivity = Vertex())

Return the backing vector's adjacency operations. Neighbour and ring results
are clipped to lookup membership; [`halo_table`](@ref) returns in-set positions
and [`halo`](@ref) lazily returns adjacent cells outside the lookup.
"""
neighbors(lk::CellLookup, c::AbstractCellIndex, k::Integer=1;
    connectivity::Connectivity=Vertex()) =
    neighbors(parent(lk), c, k; connectivity)

ring(lk::CellLookup, c::AbstractCellIndex, k::Integer;
    connectivity::Connectivity=Vertex()) = ring(parent(lk), c, k; connectivity)

neighbors(lk::CellLookup, p::Int, k::Integer=1;
    connectivity::Connectivity=Vertex()) =
    neighbors(parent(lk), p, k; connectivity)

ring(lk::CellLookup, p::Int, k::Integer;
    connectivity::Connectivity=Vertex()) = ring(parent(lk), p, k; connectivity)

neighborcount(lk::CellLookup, c::AbstractCellIndex;
    connectivity::Connectivity=Vertex()) = neighborcount(parent(lk), c; connectivity)

halo_table(lk::CellLookup, k::Integer=1; kw...) = halo_table(parent(lk), k; kw...)

halo(lk::CellLookup; connectivity::Connectivity=Vertex()) =
    halo(parent(lk); connectivity)

# Positioned handles use the parent vector's positions.
neighbors(lk::CellLookup; connectivity::Connectivity=Vertex()) =
    neighbors(parent(lk); connectivity)

# Delegate neighbourhood sweeps to the parent vector.
mapneighbors(f, lk::CellLookup; kw...) = mapneighbors(f, parent(lk); kw...)
mapneighbors(f, lk::CellLookup, data::AbstractVector; kw...) =
    mapneighbors(f, parent(lk), data; kw...)
foreachneighbors(f, lk::CellLookup; kw...) = foreachneighbors(f, parent(lk); kw...)
foreachneighbors(f, lk::CellLookup, data::AbstractVector; kw...) =
    foreachneighbors(f, parent(lk), data; kw...)
HaloTable(lk::CellLookup; kw...) = HaloTable(parent(lk); kw...)

"""
    PartialGrid(lk::CellLookup) -> PartialGrid

The lookup read as a grid: position `k` of the grid is position `k` of the
lookup, so a `Regridder` built on it lines up with a cube over the lookup's
axis without a permutation. O(1) — the ids stay lazy.
"""
PartialGrid(lk::CellLookup) = PartialGrid(parent(lk))

# --- DimensionalData plumbing ----------------------------------------------

# The values ascend in canonical id order, which is what makes the binary
# searches below — and `searchsortedfirst` on the lookup — sound.
Lookups.order(::CellLookup) = Lookups.ForwardOrdered()
Lookups.metadata(::CellLookup) = Lookups.NoMetadata()

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
Lookups.reducelookup(::CellLookup) = Lookups.NoLookup(Base.OneTo(1))

# Use the window membership search instead of searching all logical cell ids.
Lookups.hasselection(lk::CellLookup, sel::Lookups.At{<:AbstractCellIndex}) =
    cellposition(lk, Lookups.val(sel)) !== nothing

Lookups.hasselection(lk::CellLookup, sel::Lookups.Contains{<:AbstractCellIndex}) =
    cellposition(lk, Lookups.val(sel)) !== nothing

Lookups.hasselection(lk::CellLookup, sel::Lookups.Contains{<:Tuple{Real,Real}}) =
    cellposition(parent(lk), Lookups.val(sel)...) !== nothing

Dimensions.format(lk::CellLookup, ::Type, values, axis::AbstractRange) = lk

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
        DD.lookup(d) isa CellLookup && return i
    end
    throw(ArgumentError(
        "no cell dimension found: none of the dims $(map(DD.name, DD.dims(A))) " *
        "carries a CellLookup; pass spatialdim to name one"))
end

function _cells_dimnum(A::DD.AbstractDimArray, spatialdim)
    d = DD.dims(A, spatialdim)
    d === nothing && throw(ArgumentError(
        "array has no dimension matching spatialdim = $spatialdim; its dims " *
        "are $(map(DD.name, DD.dims(A)))"))
    lk = DD.lookup(d)
    lk isa CellLookup || throw(ArgumentError(
        "dimension $(DD.name(d)) carries a $(nameof(typeof(lk))) lookup, " *
        "not a CellLookup"))
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
function _map_slices(f::F, A, dnum::Int, cv::CellVector, order, threaded,
        connectivity::Connectivity) where {F}
    data = parent(A)
    pre = CartesianIndices(axes(data)[1:(dnum-1)])
    post = CartesianIndices(axes(data)[(dnum+1):end])
    M = max_neighbors(system(cv), connectivity)
    H = SubsetPositionedCell{eltype(cv)}
    T = Base.promote_op(f, H, eltype(A),
        SmallCollections.SmallVector{M,eltype(A)})
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

Lookups.selectindices(lk::CellLookup, sel::Covering; kw...) =
    covering_positions(parent(lk), Lookups.val(sel))

Lookups.selectindices(lk::CellLookup, sel::Covering{<:AbstractVector}; kw...) =
    covering_positions(parent(lk), Lookups.val(sel))

Lookups.selectindices(lk::CellLookup, sel::Lookups.At{<:AbstractCellIndex}; kw...) =
    _found(lk, cellposition(lk, Lookups.val(sel)), sel)

Lookups.selectindices(lk::CellLookup, sel::Lookups.Contains{<:AbstractCellIndex}; kw...) =
    _found(lk, cellposition(lk, Lookups.val(sel)), sel)

# Resolve a point to a leaf-grid cell, then search its position in the windows.
Lookups.selectindices(lk::CellLookup, sel::Lookups.Contains{<:Tuple{Real,Real}}; kw...) =
    _found(lk, cellposition(parent(lk), Lookups.val(sel)...), sel)

Lookups.selectindices(lk::CellLookup, sel::Lookups.At{<:Tuple{Real,Real}}; kw...) =
    _found(lk, cellposition(parent(lk), Lookups.val(sel)...), sel)

# Typed on `Lookup` so `MultiOrderLookup` resolves its selectors through the
# same two lines.
_found(lk::Lookups.Lookup, k::Int, sel) = k
_found(lk::Lookups.Lookup, ::Nothing, sel) = throw(Lookups.SelectorError(lk, sel))

# ===========================================================================
# The mixed-level lookup: the same layering as `CellLookup`, holding a
# `MultiOrderVector` and delegating every verb to it.
# ===========================================================================

"""
    MultiOrderLookup(mov::MultiOrderVector)
    MultiOrderLookup(set::MultiOrderCellSet)

A `DimensionalData` lookup over cells at mixed refinement levels: a thin
wrapper around the [`MultiOrderVector`](@ref) that holds them, reached as
`Base.parent`. Pair it with [`Cells`](@ref) to make the axis of an adaptively
refined mesh:

```julia
mov, vals = coarsen(cv, temperature; atol = 1.0)
A = DimensionalData.DimArray(vals, Cells(MultiOrderLookup(mov)))
```

`length(lk)` is the number of stored cells, `lk[k]` the `k`th of them, and
`collect(lk)` the id vector. Every selector resolves through the container's
interval index in O(log n):

```julia
A[Cells(DimensionalData.At(c))]               # c must be stored
A[Cells(DimensionalData.Contains(c))]         # the stored cell covering c
A[Cells(DimensionalData.Contains(8.0, 46.5))] # the stored cell a point falls in
A[Cells(Covering(polygon))]                   # every stored cell a region names
```

`At` is exact membership ([`cellposition`](@ref)); `Contains` resolves a cell —
including one deeper than the [`reference_level`](@ref) — to the stored
ancestor that holds its value ([`covering_position`](@ref)). `At` and
`Contains` are reached through `DimensionalData` because this package exports
DE9IM's geometry predicate [`Contains`](@ref).

An ascending subset stays a `MultiOrderLookup`; a reordered one falls back to
an unordered `DimensionalData.Categorical` of the same ids. A reduction
collapses the axis to `NoLookup`; `vcat`/`cat` of disjoint ascending axes
rebuild a `MultiOrderLookup`. [`coarsen`](@ref) constructs one from a
`Cells{<:CellLookup}` array; [`expand`](@ref) presents one back at a single
level.

Requires [`has_sorted_subtrees`](@ref); see [`MultiOrderVector`](@ref).
"""
struct MultiOrderLookup{ID,M<:MultiOrderVector} <: Lookups.Lookup{ID,1}
    cells::M
end

MultiOrderLookup(mov::MultiOrderVector{ID}) where {ID} =
    MultiOrderLookup{ID,typeof(mov)}(mov)

MultiOrderLookup(set::MultiOrderCellSet) = MultiOrderLookup(MultiOrderVector(set))

MultiOrderLookup(lk::MultiOrderLookup) = lk

# --- the collection surface ------------------------------------------------

# The VALUES, per DimensionalData's contract (see `CellLookup`'s `parent`).
Base.parent(lk::MultiOrderLookup) = lk.cells
Base.IndexStyle(::Type{<:MultiOrderLookup}) = Base.IndexLinear()

# No `bounds` method: the lookup is `Unordered` (see `order` below), and the
# generic `(nothing, nothing)` is the right answer.

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

# Ascending indices form sorted disjoint intervals again and stay in this
# type; anything else becomes a plain unordered id vector.
function _subset(lk::MultiOrderLookup, idx)
    sub = parent(lk)[idx]
    sub isa MultiOrderVector && return MultiOrderLookup(sub)
    return Lookups.Categorical(sub; order=Lookups.Unordered())
end

# --- what the lookup is, in this package's own vocabulary ------------------

"""
    system(lk::MultiOrderLookup)

The grid system the lookup's cells are named in.
"""
system(lk::MultiOrderLookup) = system(parent(lk))

"""
    reference_level(lk::MultiOrderLookup) -> Int

The level the backing container's intervals are stated at — no shallower than
any cell on the axis. Unexported; reach it as
`DiscreteGlobalGrids.Fallbacks.reference_level(lk)`.
"""
reference_level(lk::MultiOrderLookup) = reference_level(parent(lk))

"""
    cellposition(lk::MultiOrderLookup, c::AbstractCellIndex) -> Union{Int,Nothing}

Position of cell `c` on the axis, or `nothing` when the axis does not hold that
cell — including when `c` is a descendant or an ancestor of one it does hold.
EXACT membership, and what `At` resolves to.
"""
cellposition(lk::MultiOrderLookup, c::AbstractCellIndex) = cellposition(parent(lk), c)

"""
    covering_position(lk::MultiOrderLookup, c::AbstractCellIndex) -> Union{Int,Nothing}

Position of the stored cell that **is `c`, or is an ancestor of `c`** — the
compression verb, and what `Contains` resolves to. `c` may be deeper than
[`reference_level`](@ref).
"""
covering_position(lk::MultiOrderLookup, c::AbstractCellIndex) =
    covering_position(parent(lk), c)

# --- DimensionalData plumbing ----------------------------------------------

# `order` is DimensionalData's claim about the id VALUES under `isless`, which
# compares level first, so a mixed-level axis is genuinely unsorted: an ordered
# claim makes `cat`'s boundary check read level inversions as overlaps and
# silently return the parent array with the lookup dropped. Selectors resolve
# through the interval index and never read `order`.
Lookups.order(::MultiOrderLookup) = Lookups.Unordered()
Lookups.metadata(::MultiOrderLookup) = Lookups.NoMetadata()

function Lookups.rebuild(lk::MultiOrderLookup; data=nothing, kw...)
    (data === nothing || data === lk || data === parent(lk)) && return lk
    return _rebuild(lk, data)
end

_rebuild(lk::MultiOrderLookup, mov::MultiOrderVector) = MultiOrderLookup(mov)

# `vcat`/`cat` along `Cells` rebuild through here: an ascending disjoint id
# vector becomes a container again, anything else an unordered `Categorical`;
# overlapping ids error in the container's constructor.
function _rebuild(lk::MultiOrderLookup, ids::AbstractVector{<:AbstractCellIndex})
    mov = parent(lk)
    sys = system(mov)
    ref = reference_level(mov)
    ascending = true
    prev = 0
    for c in ids
        level(c) <= ref || throw(ArgumentError(
            "$c is at level $(level(c)), deeper than the axis's reference level " *
            "$ref, so it has no interval to be keyed by"))
        s = first(descendant_range(sys, c, ref))
        s > prev || (ascending = false)
        prev = s
    end
    ascending || return Lookups.Categorical(collect(ids); order=Lookups.Unordered())
    return MultiOrderLookup(MultiOrderVector(sys, ids; reference_level=ref))
end

@noinline _rebuild(lk::MultiOrderLookup, data) = throw(ArgumentError(
    "a MultiOrderLookup holds cell ids at mixed levels; it cannot be rebuilt " *
    "around $(typeof(data)). Concatenate cell axes with `vcat`/`cat`, subset " *
    "them by indexing, and replace one wholesale with " *
    "`set(A, Cells => NoLookup())`."))

Lookups.reducelookup(::MultiOrderLookup) = Lookups.NoLookup(Base.OneTo(1))

# One binary search each, instead of the generic scan over the values.
Lookups.hasselection(lk::MultiOrderLookup, sel::Lookups.At{<:AbstractCellIndex}) =
    cellposition(lk, Lookups.val(sel)) !== nothing

Lookups.hasselection(lk::MultiOrderLookup, sel::Lookups.Contains{<:AbstractCellIndex}) =
    covering_position(lk, Lookups.val(sel)) !== nothing

Lookups.hasselection(lk::MultiOrderLookup, sel::Lookups.Contains{<:Tuple{Real,Real}}) =
    cellposition(parent(lk), Lookups.val(sel)...) !== nothing

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
#
# One `MultiOrderVector` verb each: `At` is exact, `Contains` is covering.
# Typed on the selector's value, as in the `CellLookup` block above.

Lookups.selectindices(lk::MultiOrderLookup, sel::Covering; kw...) =
    covering_positions(parent(lk), Lookups.val(sel))

Lookups.selectindices(lk::MultiOrderLookup, sel::Covering{<:AbstractVector}; kw...) =
    covering_positions(parent(lk), Lookups.val(sel))

Lookups.selectindices(lk::MultiOrderLookup, sel::Lookups.At{<:AbstractCellIndex}; kw...) =
    _found(lk, cellposition(parent(lk), Lookups.val(sel)), sel)

# A cell the axis does not store — including one deeper than the reference
# level — resolves to its stored ancestor.
Lookups.selectindices(lk::MultiOrderLookup,
    sel::Lookups.Contains{<:AbstractCellIndex}; kw...) =
    _found(lk, covering_position(parent(lk), Lookups.val(sel)), sel)

Lookups.selectindices(lk::MultiOrderLookup,
    sel::Lookups.Contains{<:Tuple{Real,Real}}; kw...) =
    _found(lk, cellposition(parent(lk), Lookups.val(sel)...), sel)

Lookups.selectindices(lk::MultiOrderLookup,
    sel::Lookups.At{<:Tuple{Real,Real}}; kw...) =
    _found(lk, cellposition(parent(lk), Lookups.val(sel)...), sel)

# ===========================================================================
# The compression, presented at one level
# ===========================================================================

"""
    MultiOrderValues(values, offsets) <: AbstractVector

[`expand`](@ref)'s data: mixed-level values presented as the values of the
leaf cells they cover, stored as themselves. `offsets[i]` is the number of
leaf cells the first `i` stored cells cover, so `v[k]` is
`values[searchsortedfirst(offsets, k)]` and nothing is materialised:
`length(v)` is the leaf count while `Base.summarysize(v)` stays O(#stored
values).

Internal: an `expand`ed array's data, reached as `parent(A)`.
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

# A run fill, so materialising the whole expansion costs O(#leaves) rather
# than one binary search per leaf.
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

# ===========================================================================
# DimArray methods for the aggregation verbs
#
# Each is one core call with the axis rebuilt around the answer. One `Cells`
# dimension only: the cores read `(cells, values)` as parallel vectors.
# ===========================================================================

# Shared validation; the error names the verb the caller reached for.
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

Reduce a cell array to level `l`: one element per distinct level-`l` ancestor
of `A`'s cells, each carrying `f` of the values of its **present** descendants.
`A` must be one-dimensional over a [`Cells`](@ref) axis with a
[`CellLookup`](@ref); the answer carries the coarser lookup. A pyramid is this
call once per level:

```julia
pyramid = [aggregate(sum, A, l) for l in level(lookup(A, Cells)) - 1 : -1 : 3]
```

See the `(CellVector, values)` method for what `f` sees and how partial groups
are treated.
"""
function aggregate(f, A::DD.AbstractDimArray, l::Integer)
    lk = _cell_axis(A, CellLookup, "aggregate")
    coarse, vals = aggregate(f, parent(lk), parent(A), l)
    return DD.rebuild(A; data=vals, dims=(Cells(CellLookup(coarse)),))
end

"""
    coarsen(A::AbstractDimArray; atol, by = mean, minlevel = shallowest) -> AbstractDimArray

Merge each subtree of `A` whose values agree to within `atol` into the single
coarse cell that stands for them. `A` must be one-dimensional over a
[`Cells`](@ref) axis with a [`CellLookup`](@ref); the answer carries a
[`MultiOrderLookup`](@ref) over the mixed-level container.

```julia
M = coarsen(A; atol = 1.0)                       # °C
M[Cells(DimensionalData.Contains(8.0, 46.5))]    # a point, on the mesh
expand(M, level(lookup(A, Cells)))               # back, within `atol`
```

See the `(CellVector, values)` method for the merge criterion, the treatment
of `missing`, and the error bound the default `by` carries.
"""
function coarsen(A::DD.AbstractDimArray; atol, kw...)
    lk = _cell_axis(A, CellLookup, "coarsen")
    mov, vals = coarsen(parent(lk), parent(A); atol, kw...)
    return DD.rebuild(A; data=vals, dims=(Cells(MultiOrderLookup(mov)),))
end

"""
    expand(A::AbstractDimArray, l::Integer) -> AbstractDimArray

Present a mixed-level array at level `l`. `A` must be one-dimensional over a
[`Cells`](@ref) axis with a [`MultiOrderLookup`](@ref); the answer carries a
[`CellLookup`](@ref) at `l`, and its data is a lazy `MultiOrderValues` — still
one stored value per multi-order cell, so `Base.summarysize` stays O(#stored
cells) however many leaves `l` names. `l` must be no shallower than the
deepest stored cell.

```julia
E = expand(coarsen(A; atol), level(lookup(A, Cells)))
all(abs.(collect(parent(E)) .- parent(A)) .<= atol)   # the round-trip bound
```
"""
function expand(A::DD.AbstractDimArray, l::Integer)
    lk = _cell_axis(A, MultiOrderLookup, "expand")
    mov = parent(lk)
    target = Int(l)
    cv = CellVector(mov; level=target)
    data = MultiOrderValues(parent(A), _leaf_offsets(mov, target))
    return DD.rebuild(A; data=data, dims=(Cells(CellLookup(cv)),))
end

# Cumulative leaf counts at `l`, one per stored cell. The container's own
# `offsets` state these only at its reference level.
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

end # module CellLookups
