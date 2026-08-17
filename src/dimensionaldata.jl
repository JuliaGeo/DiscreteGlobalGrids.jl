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
    AbstractCellIndex, ncells, cellindex, cellposition, cellat, level, system,
    levelgrid, cellindextype, has_sorted_subtrees, descendants, query,
    neighbors, ring, halo_table, halo, neighborcount, Connectivity, Vertex,
    max_neighbors
import ..DiscreteGlobalGrids: Helpers
import ..DiscreteGlobalGrids.Fallbacks: PartialGrid, SubtreeIds,
    MultiOrderCoverage, MultiOrderCellSet, level_ranges
# Core collection operations delegated to `CellVector`.
import ..DiscreteGlobalGrids.Fallbacks: CellVector, cellset, covering,
    covering_positions, windows, nwindows, RangeWindows, CellWindows, _derive,
    _windows, SubsetPositionedCell, mapneighbors, foreachneighbors, HaloTable,
    StorageOrder

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

# The first `Cells` dimension's number, for handle indexing on any shape.
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

# On more than one dimension a handle names the whole slice at its position
# along the first `Cells` dimension, as a view; the 1-D methods above stay
# the scalar fast path, position trusted either way.
Base.getindex(A::DD.AbstractDimArray, h::SubsetPositionedCell) =
    _handle_slice(A, Val(_handle_dimnum(A)), h.position)
Base.view(A::DD.AbstractDimArray, h::SubsetPositionedCell) =
    _handle_slice(A, Val(_handle_dimnum(A)), h.position)

# ===========================================================================
# Whole-array entry points
#
# The neighborhood forms accept any `AbstractDimArray` and resolve the cell
# dimension themselves, so a downstream package never unwraps a lookup or
# slices by hand.
# ===========================================================================

# The cell dimension's number: the first `CellLookup` dim when `spatialdim` is
# `nothing`, otherwise whatever `DD.dims` accepts — informative errors both
# ways.
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

The default `pass` of [`mapneighbors`](@ref) and [`foreachneighbors`](@ref)
over an `AbstractDimArray`: `f(cell, nbrs)` receives positioned handles —
exactly the [`CellVector`](@ref) form — and reads data by indexing the array
with them. One result per cell, returned on the cell dimension.
"""
struct Neighbors end

"""
    Values()

`pass` form of [`mapneighbors`](@ref) over an `AbstractDimArray`:
`f(cell, value, values)` receives the cell's scalar and its neighbours'
scalars. On an array with other dimensions the stencil runs independently on
every cell-dimension slice — the `mapslices` semantic — and the output keeps
all of the array's dims.
"""
struct Values end

"""
    NeighborSlices()

`pass` form of [`mapneighbors`](@ref) for arrays of two or more dimensions:
`f(cell, slice, slices)` receives the array's views at the cell's and its
neighbours' positions along the cell dimension, each one dimension smaller
than the array. One result per cell, on the cell dimension. A
one-dimensional array is refused — its per-cell slice is a scalar, which is
[`Values`](@ref).
"""
struct NeighborSlices end

@noinline _bad_pass(pass) = throw(ArgumentError(
    "pass must be Neighbors(), Values() or NeighborSlices(), " *
    "got $(typeof(pass))"))

_need_slices(A) = ndims(A) >= 2 || throw(ArgumentError(
    "NeighborSlices() needs at least two dimensions; a one-dimensional " *
    "array's per-cell slice is its scalar — use Values()"))

# One result per cell lands on the cell dimension alone, lookup and all.
_rebuilt_on_cells(A, d, out::Tuple) =
    map(o -> DD.rebuild(A; data = o, dims = (d,)), out)
_rebuilt_on_cells(A, d, out) = DD.rebuild(A; data = out, dims = (d,))

"""
    mapneighbors(f, A::AbstractDimArray; spatialdim = nothing, pass = Neighbors(),
                 order = StorageOrder(), threaded = true, connectivity = Vertex())

[`mapneighbors`](@ref) over a whole array: the cell dimension is resolved,
the sweep runs on its [`CellLookup`](@ref)'s vector, and the result is `A`'s
own wrapper rebuilt with the lookups intact — a concrete-tuple-returning `f`
gets one rebuilt array per component.

`spatialdim` names the cell dimension in any form `DimensionalData.dims`
accepts; `nothing` selects the first dimension carrying a `CellLookup`. An
array without one, or a named dimension that is missing or not cell-valued,
is an `ArgumentError`.

`pass` selects what `f` sees: [`Neighbors`](@ref) (the default) hands
`f(cell, nbrs)` positioned handles and returns one value per cell on the
cell dimension; [`Values`](@ref) hands `f(cell, value, values)` scalars,
slice by independent slice on arrays with other dimensions, and returns
`A`'s full shape; [`NeighborSlices`](@ref) hands `f(cell, slice, slices)`
the views along the cell dimension and returns one value per cell.
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

# Function barrier making the cell dimension's number a constant, so the
# slice views are concretely typed.
function _map_cell_slices(f::F, A, ::Val{D}, cv, order, threaded,
        conn) where {F,D}
    g = (c, nbrs) -> f(c, _handle_slice(A, Val(D), cellposition(c)),
        [_handle_slice(A, Val(D), cellposition(h)) for h in nbrs])
    return mapneighbors(g, cv; order, threaded, connectivity = conn)
end

# The slice loop: each other-dim index gets its own buffered 1-D sweep, so
# the CellVector kernels own all traversal and the slices cannot interact.
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

[`mapneighbors`](@ref)' whole-array form without the output: same
`spatialdim` resolution, same `pass` forms, `f` called for its side effects.
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

The one-arg positioned iterator over `A`'s cell dimension, resolved as in
[`mapneighbors`](@ref); the minted positions are that dimension's axis
positions.
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

_found(lk::CellLookup, k::Int, sel) = k
_found(lk::CellLookup, ::Nothing, sel) = throw(Lookups.SelectorError(lk, sel))

end # module CellLookups
