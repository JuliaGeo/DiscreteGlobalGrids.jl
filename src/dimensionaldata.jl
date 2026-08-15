# DimensionalData wrappers for `Fallbacks.CellVector`. Cell ids are computed
# from compressed position windows rather than stored as an expanded vector.

"""
    CellLookups

The DimensionalData layer: [`CellLookup`](@ref) and [`MultiOrderLookup`](@ref),
the [`Cells`](@ref) dimension both go in, and the [`Covering`](@ref) selector.

A [`CellLookup`](@ref) is a one-dimensional `DimensionalData` lookup over cell
ids at a single level. It is a thin wrapper around a [`CellVector`](@ref),
which is where the compression lives: a set of **leaf position windows** —
sorted, disjoint intervals (or, where intervals are unavailable, a sorted list)
of positions in `levelgrid(sys, leaf)`. Its logical content is their
concatenation, and every operation is arithmetic over that concatenation:
`length` sums the window lengths, `lk[k]` binary-searches the cumulative
lengths and resolves one `cellindex`, [`cellposition`](@ref) runs the inverse.
Nothing is materialised.

A [`MultiOrderLookup`](@ref) is the same thing one layer up: a lookup over cells
at MIXED levels, a thin wrapper around a [`MultiOrderVector`](@ref), and the axis
an adaptively refined mesh's values hang off. [`coarsen`](@ref) builds one from a
`CellLookup` axis, [`expand`](@ref) presents one back at a single level, and
[`aggregate`](@ref) is the fixed-level reduction beside them.

`CellVector` provides the storage and indexing behavior; this module provides
the lookup, dimension, and selectors required by DimensionalData.
"""
module CellLookups

import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, ncells, cellindex, cellposition, cellat, level, system,
    levelgrid, cellindextype, has_sorted_subtrees, descendants, descendant_range,
    query, neighbors, ring, halo_table, halo, Connectivity, Vertex
import ..DiscreteGlobalGrids: Helpers
import ..DiscreteGlobalGrids.Fallbacks: PartialGrid, SubtreeIds,
    MultiOrderCoverage, MultiOrderCellSet, level_ranges
# Core collection operations delegated to `CellVector`.
import ..DiscreteGlobalGrids.Fallbacks: CellVector, cellset, covering,
    covering_positions, windows, nwindows, RangeWindows, CellWindows, _derive,
    _windows
# And the mixed-level half of the same story: the container, its two membership
# verbs, and the three functions born in `src/fallbacks/aggregate.jl` whose
# cube-shaped methods are defined at the bottom of this file.
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

halo_table(lk::CellLookup, k::Integer=1; connectivity::Connectivity=Vertex()) =
    halo_table(parent(lk), k; connectivity)

halo(lk::CellLookup; connectivity::Connectivity=Vertex()) =
    halo(parent(lk); connectivity)

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

# Typed on `Lookup` rather than on `CellLookup`: the mixed-level lookup below
# resolves its selectors through the same two lines, and a miss is the same
# error there for the same reason.
_found(lk::Lookups.Lookup, k::Int, sel) = k
_found(lk::Lookups.Lookup, ::Nothing, sel) = throw(Lookups.SelectorError(lk, sel))

# ===========================================================================
# The mixed-level lookup
#
# Everything above is one level; everything below is a mixture of them. The
# layering is identical — a `Lookup` that HOLDS a DimensionalData-free container
# and delegates every verb to it — because the reason for that layering
# (`parent` is the VALUES, see the long comment at `Base.parent` above) has
# nothing to do with how many levels the values sit at.
# ===========================================================================

"""
    MultiOrderLookup(mov::MultiOrderVector)
    MultiOrderLookup(set::MultiOrderCellSet)

A `DimensionalData` lookup naming cells at MIXED refinement levels, backed by
the [`MultiOrderVector`](@ref) that holds them. Pair it with [`Cells`](@ref) to
make the axis of an adaptively refined mesh:

```julia
mov, vals = coarsen(cv, temperature; atol = 1.0)
A = DimensionalData.DimArray(vals, Cells(MultiOrderLookup(mov)))
```

Semantically `lk` is the container's id vector: `length(lk)` is the number of
stored cells, `lk[k]` is the `k`th of them — at whatever level that one sits —
and `collect(lk)` is the vector itself.

# The layering

A `MultiOrderLookup` **is** a [`MultiOrderVector`](@ref) wearing a
`DimensionalData.Lookup` hat, and holds nothing else, exactly as
[`CellLookup`](@ref) is a [`CellVector`](@ref) wearing one. `Base.parent` is the
VALUES — the container — for the reason spelled out at `CellLookup`'s own
`parent`: some thirty `Lookup` methods derive their behaviour from it, and a
`parent` answering with anything else has to shadow every one of them.

The container is where the acceleration lives: a sorted interval index at one
[`reference_level`](@ref), so a point, a cell or a region resolves to a position
in O(log n) whatever the mixture of levels.

# Selectors, and the two questions they ask

```julia
A[Cells(DimensionalData.At(c))]               # c must be STORED
A[Cells(DimensionalData.Contains(c))]         # the stored cell COVERING c
A[Cells(DimensionalData.Contains(8.0, 46.5))] # the stored cell a point falls in
A[Cells(Covering(polygon))]                   # every stored cell a region names
```

`At` is exact membership ([`cellposition`](@ref)) and `Contains` is the
compression verb ([`covering_position`](@ref)): a leaf resolves to the ancestor
that holds its value. That is `Contains`' interval reading, which is
`DimensionalData`'s own for an interval lookup — and it is why a cell DEEPER
than the reference level answers too, through its reference-level ancestor.
Outside a cube the four are `cellposition(mov, c)`, `covering_position(mov, c)`,
`cellposition(mov, lon, lat)` and [`covering`](@ref)`(mov, polygon)`.

`At` and `Contains` are `DimensionalData`'s and are reached through it, for the
reason [`CellLookup`](@ref) gives: this package exports DE9IM's
[`Contains`](@ref), a predicate about two geometries, and the two names must
never collide in a caller's namespace.

# What the cube's own operations do to it

The same fork [`CellLookup`](@ref) takes, for the same reason — an ascending
subset is a sorted disjoint interval list again and stays in this type; a
reordered one is not, and falls back to an unordered
`DimensionalData.Categorical` of the same ids. A reduction collapses the axis to
`NoLookup`, and `vcat`/`cat` of two disjoint ascending axes rebuild into a
`MultiOrderLookup` again.

# The verbs that build one and take one apart

[`coarsen`](@ref) is the constructor from data — a `Cells{<:CellLookup}` array
in, a `Cells{<:MultiOrderLookup}` array out — and [`expand`](@ref) is the
inverse presentation, which names one level's cells while still storing one
value per multi-order cell.

Requires [`has_sorted_subtrees`](@ref); see [`MultiOrderVector`](@ref) for why
A5 has no container to wrap.
"""
struct MultiOrderLookup{ID,M<:MultiOrderVector} <: Lookups.Lookup{ID,1}
    cells::M
end

MultiOrderLookup(mov::MultiOrderVector{ID}) where {ID} =
    MultiOrderLookup{ID,typeof(mov)}(mov)

MultiOrderLookup(set::MultiOrderCellSet) = MultiOrderLookup(MultiOrderVector(set))

MultiOrderLookup(lk::MultiOrderLookup) = lk

# --- the collection surface ------------------------------------------------

# The VALUES, per DimensionalData's contract. See the long comment at
# `CellLookup`'s `parent` for what depends on that.
Base.parent(lk::MultiOrderLookup) = lk.cells
Base.IndexStyle(::Type{<:MultiOrderLookup}) = Base.IndexLinear()

# No `bounds` method: the lookup is `Unordered` (see `order` below), and
# DimensionalData's generic answers `(nothing, nothing)` for exactly that
# reason — a `(first, last)` pair here would be an inverted-`isless` pair
# masquerading as an interval.

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

# The container's fork, wearing a lookup's clothes: ascending indices are a
# sorted disjoint interval list again, anything else is a plain id vector and
# takes the ordinary lookup that can hold one.
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

# Unordered, and deliberately so. The container IS ordered in its own key — the
# reference-level interval starts ascend strictly, and every verb here searches
# them — but `order` is DimensionalData's claim about the VALUES under `isless`,
# and `isless` on an id compares its level first, so a mixed-level axis is
# unsorted in exactly the sense the claim is read in. The reading is not
# hypothetical: `cat`'s boundary check walks `last(l1) < first(l2)` on the ids
# whenever the lookup claims an order, and a coarse cell following a deep one
# then reads as an overlap — DimensionalData warns and silently hands back the
# PARENT array, axis gone, at whichever split indices the levels happen to
# invert. Under `Unordered` the same check asks the right question instead —
# are the two id sets disjoint — and the join reaches `rebuild` below, where
# the interval index decides.
#
# What the claim costs is only DimensionalData's value-ordered machinery
# (`Near`, `a .. b`), which was never part of this axis's contract. `At`,
# `Contains` and [`Covering`](@ref) resolve through the interval index and do
# not read `order`. (`CellLookup` says `ForwardOrdered` and means it: on a
# single level the ids do ascend as values.)
Lookups.order(::MultiOrderLookup) = Lookups.Unordered()
Lookups.metadata(::MultiOrderLookup) = Lookups.NoMetadata()

function Lookups.rebuild(lk::MultiOrderLookup; data=nothing, kw...)
    (data === nothing || data === lk || data === parent(lk)) && return lk
    return _rebuild(lk, data)
end

_rebuild(lk::MultiOrderLookup, mov::MultiOrderVector) = MultiOrderLookup(mov)

# `vcat`/`cat` along a `Cells` dimension hand the joined id vector back through
# here. An ascending disjoint union of subtrees is a container again; anything
# else is not, and takes the same honest fallback `getindex` takes. Overlap —
# an ancestor beside its own descendant — is neither, and the container's own
# constructor is what says so.
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

# The two membership questions, asked directly rather than left to a scan over
# the values: both are one binary search over the interval index.
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
# One `MultiOrderVector` verb each, and the pairing is the whole semantics:
# `At` is exact, `Contains` is covering. Typed on the selector's VALUE as well
# as on the lookup, for the reason the `CellLookup` block above gives.

Lookups.selectindices(lk::MultiOrderLookup, sel::Covering; kw...) =
    covering_positions(parent(lk), Lookups.val(sel))

Lookups.selectindices(lk::MultiOrderLookup, sel::Covering{<:AbstractVector}; kw...) =
    covering_positions(parent(lk), Lookups.val(sel))

Lookups.selectindices(lk::MultiOrderLookup, sel::Lookups.At{<:AbstractCellIndex}; kw...) =
    _found(lk, cellposition(parent(lk), Lookups.val(sel)), sel)

# The compression verb: a cell the axis does not store resolves to the stored
# ancestor that speaks for it — including a cell DEEPER than the reference
# level, which the container keys through its reference-level ancestor.
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

[`expand`](@ref)'s data: the values of a mixed-level array, presented as the
values of the leaf cells they cover, and stored as themselves.

`offsets[i]` is the number of leaf cells the first `i` stored cells name at the
expansion level, so `v[k]` is `values[searchsortedfirst(offsets, k)]` — one
binary search, nothing materialised. `length(v)` is `offsets[end]`, which is the
leaf count; `Base.summarysize(v)` is O(length(values)), which is the stored
count. Expanding three levels deeper multiplies the first and does not move the
second.

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

# Sequential materialisation is a run fill, not `length(v)` binary searches:
# that is what makes reading the whole expansion once — plotting it, writing it
# — cost the leaves rather than the leaves times a log.
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
# The cube-shaped face of the aggregation verbs
#
# Each is one core call with an axis rebuilt around the answer. One dimension
# only in v1: the cores read `(cells, values)` as parallel vectors, and a cube
# with a second dimension would have to say which slice of it the tolerance and
# the reducer speak about — a decision with no obvious default, so it is refused
# rather than guessed.
# ===========================================================================

# The one gate all three go through, so that the restriction is stated once and
# the error names the verb the caller reached for.
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

Reduce a cell array to level `l`: one element per distinct level-`l` ancestor of
`A`'s cells, each carrying `f` of the values of its **present** descendants.

The cube-shaped [`aggregate`](@ref): `A` carries a [`Cells`](@ref) axis whose
lookup is a [`CellLookup`](@ref), and the answer carries the coarser one. A
pyramid is this call once per level:

```julia
pyramid = [aggregate(sum, A, l) for l in level(lookup(A, Cells)) - 1 : -1 : 3]
```

See the `(CellVector, values)` method for what `f` sees, how partial groups are
treated, and why policy about `missing` stays with the caller.

!!! note "One dimension, in v1"
    `A` must be one-dimensional over `Cells`; anything else is an
    `ArgumentError`. The core reads cells and values as parallel vectors, and a
    second dimension would have to name the slice the reduction speaks about.
"""
function aggregate(f, A::DD.AbstractDimArray, l::Integer)
    lk = _cell_axis(A, CellLookup, "aggregate")
    coarse, vals = aggregate(f, parent(lk), parent(A), l)
    return DD.rebuild(A; data=vals, dims=(Cells(CellLookup(coarse)),))
end

"""
    coarsen(A::AbstractDimArray; atol, by = mean, minlevel = shallowest) -> AbstractDimArray

Merge each subtree of `A` whose values agree to within `atol` into the single
coarse cell that stands for them: a cell array in, an adaptively refined mesh
out.

The cube-shaped [`coarsen`](@ref): `A` carries a [`Cells`](@ref) axis whose
lookup is a [`CellLookup`](@ref), and the answer carries a
[`MultiOrderLookup`](@ref) over the mixed-level container.

```julia
M = coarsen(A; atol = 1.0)                       # °C
length(M) < length(A)                            # what the tolerance bought
M[Cells(DimensionalData.Contains(8.0, 46.5))]    # a point, on the mesh
expand(M, level(lookup(A, Cells)))               # and back, within `atol`
```

See the `(CellVector, values)` method for the merge criterion, the three
readings of `missing`, and the error bound the default `by` buys.

!!! note "One dimension, in v1"
    `A` must be one-dimensional over `Cells`; anything else is an
    `ArgumentError`, for the reason [`aggregate`](@ref) gives.
"""
function coarsen(A::DD.AbstractDimArray; atol, kw...)
    lk = _cell_axis(A, CellLookup, "coarsen")
    mov, vals = coarsen(parent(lk), parent(A); atol, kw...)
    return DD.rebuild(A; data=vals, dims=(Cells(MultiOrderLookup(mov)),))
end

"""
    expand(A::AbstractDimArray, l::Integer) -> AbstractDimArray

Present a mixed-level array at level `l`: the answer's axis names
`CellVector(mov; level = l)`'s cells, and its data is still one value per
multi-order cell.

The inverse of [`coarsen`](@ref) as a *presentation*, and the compression this
whole file exists for. `A` carries a [`Cells`](@ref) axis whose lookup is a
[`MultiOrderLookup`](@ref); the answer carries a [`CellLookup`](@ref) at `l`,
and its data is a lazy `MultiOrderValues` — `Base.summarysize` of it is
O(#stored cells) however many leaves `l` names, so re-expanding three levels
deeper multiplies the axis's length and moves nothing.

```julia
E = expand(coarsen(A; atol), level(lookup(A, Cells)))
all(abs.(collect(parent(E)) .- parent(A)) .<= atol)   # the round-trip bound
```

`l` must be no shallower than the deepest stored cell, which is what the
container's own expansion requires.

!!! note "One dimension, in v1"
    `A` must be one-dimensional over `Cells`; anything else is an
    `ArgumentError`, for the reason [`aggregate`](@ref) gives.
"""
function expand(A::DD.AbstractDimArray, l::Integer)
    lk = _cell_axis(A, MultiOrderLookup, "expand")
    mov = parent(lk)
    target = Int(l)
    cv = CellVector(mov; level=target)
    data = MultiOrderValues(parent(A), _leaf_offsets(mov, target))
    return DD.rebuild(A; data=data, dims=(Cells(CellLookup(cv)),))
end

# Cumulative leaf counts at `l`, one entry per STORED cell. The container's own
# `offsets` are the same thing at its reference level and cannot be reused for
# another; neither can the expanded ranges, which merge adjacent subtrees and so
# lose exactly the per-cell boundaries this needs.
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
