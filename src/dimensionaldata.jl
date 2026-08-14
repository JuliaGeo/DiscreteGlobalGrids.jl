# A DimensionalData cell axis. Every cell in it sits at one leaf level, but the
# stored form is the multi-order set the region came back as, never the expanded
# id vector: memory is O(#coverage entries), the ids are computed on demand.
#
# All of that arithmetic is `Fallbacks.CellVector`, which is DimensionalData-free
# and lives in `src/fallbacks/cell_vector.jl`. This file is the cube-shaped face
# of it: a `Lookup` that HOLDS one, the `Cells` dimension it goes in, and the
# selectors. Every method below either delegates to a `CellVector` verb or is
# DimensionalData plumbing that has no counterpart outside a cube.

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

The split is deliberate. `CellVector` is usable — and is used — with no cube
library in sight; this module is what makes one an axis.
"""
module CellLookups

import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, ncells, cellindex, cellposition, cellat, level, system,
    levelgrid, cellindextype, has_sorted_subtrees, descendants, query,
    neighbors, ring, halo_table, Connectivity, Vertex
import ..DiscreteGlobalGrids: Helpers
import ..DiscreteGlobalGrids.Fallbacks: PartialGrid, SubtreeIds,
    MultiOrderCoverage, MultiOrderCellSet, level_ranges
# The core type and its verbs. Everything this file does to a cell axis, it does
# by calling one of these.
import ..DiscreteGlobalGrids.Fallbacks: CellVector, cellset, covering,
    covering_positions, windows, nwindows, RangeWindows, CellWindows, _derive,
    _windows

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

A `DimensionalData` lookup naming cells at ONE level, backed by the compact set
they came from. Pair it with [`Cells`](@ref) to make a cube axis:

```julia
set = query(sys, MultiOrderCoverage(region); level = 9)
lk  = CellLookup(set)
A   = DimensionalData.DimArray(values, Cells(lk))
```

Semantically `lk` is the leaf id vector: `length(lk)` is the number of leaf
cells, `lk[k]` is the `k`th of them, `collect(lk)` is the vector itself.

# The layering

A `CellLookup` **is** a [`CellVector`](@ref) wearing a `DimensionalData.Lookup`
hat, and holds nothing else. The `CellVector` is where the compression lives —
the set's expansion to sorted, disjoint position windows at the leaf level
([`level_ranges`](@ref)), so memory is O(#entries in the set) rather than
O(#leaf cells). On a Switzerland-sized region at IGEO7 level 9 that is some
hundreds of words standing for tens of thousands of cells, and the same words
for the level-12 re-expansion, which names twenty million.

That type is DimensionalData-free on purpose: regridding, chunking and plain
`Array` code use it directly, and this module exists so that a *cube* can too.
Every method here delegates to one of its verbs — `lk[k]`,
[`cellposition`](@ref), [`cellset`](@ref), [`covering`](@ref),
[`PartialGrid`](@ref) — and adds only what a `Lookup` owes DimensionalData.

`Base.parent` returns the lookup's VALUES, as `DimensionalData` requires: the
`CellVector`, which is an `AbstractVector` of the ids, is O(#windows) and
materialises nothing. [`cellset`](@ref) returns the backing — the set, or the
grid — for running a second coverage operation against without unpacking the
lookup.

# The three ways in

  - a [`MultiOrderCellSet`](@ref), the compressed form, optionally re-expanded
    to a deeper `level` than the set's own reference level;
  - `levelgrid(sys, l)`, a whole level, which is one window;
  - a [`PartialGrid`](@ref), an arbitrary ascending subset, which is that
    subset's positions — one window when the subset is a subtree, and the
    explicit list when it is scattered.

The last two are the degenerate cases of the first, and answer every method
below identically. All three build the [`CellVector`](@ref) first; a caller who
already has one hands it over directly.

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

`At` and `Contains` are `DimensionalData`'s selectors and are reached through
it (`DD.At`, `DD.Contains`) rather than re-exported here. That is not an
oversight: this package exports DE9IM's [`Contains`](@ref), a *predicate about
two geometries*, and a `using` of both packages would otherwise leave one name
meaning two unrelated things. Only [`Covering`](@ref), which `DimensionalData`
has no spelling for, is exported.

# What the cube's own operations do to it

Indexing, `vcat`/`cat` and reductions all reach the lookup, and each answers
with the most specific thing that is still true:

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
    The lookup is then built by SELECTION: `descendants` names the leaves, they
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

# The keyword shadows the `level` function, so it is forwarded by name to the
# core constructor, which takes the same care one call in.
CellLookup(set::MultiOrderCellSet; level::Integer=set.reference_level) =
    CellLookup(CellVector(set; level=level))

CellLookup(grid::AbstractGrid) = CellLookup(CellVector(grid))

CellLookup(lk::CellLookup) = lk

# A subset of an existing lookup is a subset of its cell vector, wrapped again.
_derive(lk::CellLookup, w::CellWindows) = CellLookup(_derive(parent(lk), w))

windows(lk::CellLookup) = windows(parent(lk))

# --- the collection surface ------------------------------------------------

# `parent` is the VALUES — the lazy `CellVector` — and nothing else. That is
# DimensionalData's contract, and it is load-bearing rather than decorative:
# some thirty `Lookup` methods derive their behaviour from `parent(l)`, and a
# `parent` that answered with the backing instead had to be shadowed by a
# hand-written override for each of them. Every one that was missed became a
# crash inside DimensionalData's own machinery, or — for `Where`, which filters
# `parent(lookup)` and returns the surviving indices as positions — a plausible
# subset of the wrong collection, with no error at all. The backing is reached
# through [`cellset`](@ref), which is this package's own name and cannot be
# mistaken for the values by code that has never heard of it.
#
# So the list below is short by design. It is the methods where the lazy form
# beats the generic one, not a re-implementation of the generic ones.
Base.parent(lk::CellLookup) = lk.cells
Base.IndexStyle(::Type{<:CellLookup}) = Base.IndexLinear()

# `first`/`last` on an empty lookup are a `BoundsError`, so the one caller that
# asks about an empty axis rather than indexing it answers `nothing` instead.
Lookups.bounds(lk::CellLookup) = isempty(lk) ? (nothing, nothing) : (first(lk), last(lk))

Base.@propagate_inbounds Base.getindex(lk::CellLookup, k::Int) = parent(lk)[k]
Base.@propagate_inbounds Base.getindex(lk::CellLookup, k::CartesianIndex{1}) = parent(lk)[k[1]]

# `AbstractVector`, not `AbstractArray`, for the reason the core file gives at
# the same signature: an index with a shape wants an answer with that shape, and
# a window set has none to give. A matrix index falls through to Base's generic
# here too, so both faces answer it the same way.
for f in (:getindex, :view, :dotview)
    @eval Base.$f(lk::CellLookup, ::Colon) = lk
    @eval Base.$f(lk::CellLookup, i::AbstractVector{<:Integer}) = _subset(lk, i)
end

# DimensionalData reverses only the lookups it knows, and a `Lookup` it does not
# know falls through to `Base.reverse` on the AbstractArray — which answers with
# a bare `Vector` where a lookup is expected, and everything downstream of the
# reversed dimension then recurses on it. Routing through `getindex` reverses
# into the type the rest of this file understands.
Base.reverse(lk::CellLookup) = lk[lastindex(lk):-1:firstindex(lk)]

# SmallCollections' own `getindex(::AbstractVector, ::AbstractFixedOrSmall...)`
# is neither more nor less specific than the line above, and a neighbour list is
# exactly one of those vectors, so the tie is broken towards the same subset
# rather than left as an ambiguity for whoever indexes an axis by a halo.
Base.getindex(lk::CellLookup,
    i::SmallCollections.AbstractFixedOrSmallOrPackedVector{<:Integer}) = _subset(lk, i)

# A mask names a position by its index, so one of the wrong length is a bounds
# error rather than a shorter answer. The core checks this too, and would catch
# it a moment later through the delegation below; the check is repeated here so
# that the error names the axis the caller indexed rather than the cell vector
# behind it.
function _subset(lk::CellLookup, mask::AbstractArray{Bool})
    axes(mask) == axes(lk) || throw(BoundsError(lk, (mask,)))
    return _subset(lk, findall(mask))
end

# The fork is the cell vector's: an ascending index set is windows again and
# comes back as a `CellVector`, anything else is a plain id vector. This file
# only has to say what a *lookup* wears in each case — its own type, or the
# ordinary DimensionalData lookup that can hold an unordered list.
function _subset(lk::CellLookup, idx)
    sub = parent(lk)[idx]
    sub isa CellVector && return CellLookup(sub)
    return Lookups.Categorical(sub; order=Lookups.Unordered())
end

# --- what the lookup is, in this package's own vocabulary ------------------

"""
    cellset(lk::CellLookup)

The thing the lookup was built from — a [`MultiOrderCellSet`](@ref), or the grid
— for running a second coverage operation against without unpacking the lookup.

`Base.parent(lk)` is deliberately NOT this: it is the lookup's VALUES, the
[`CellVector`](@ref), because that is what `DimensionalData` derives some thirty
`Lookup` methods from. A subset produced by indexing or by a selector reports
the [`PartialGrid`](@ref) describing it.
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

Adjacency on the axis: the system's neighbourhood clipped to the cells the
lookup holds, in ids or in positions along the axis. Every one of them is the
[`CellVector`](@ref)'s, because a `CellLookup` is one wearing a `Lookup` hat —
so a stencil over a cube and a stencil over a plain vector are the same call
with the same answer.
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

# A `CellLookup` has no properties to vary — its values are its windows — but
# it does have to survive being rebuilt around new VALUES, because that is how
# DimensionalData concatenates: `vcat`/`cat` along a `Cells` dimension hand the
# joined id vector back through `rebuild`. An ascending disjoint union of window
# sets is still a window set, so the join stays in this type; anything else is
# not a set of windows and takes the same honest fallback `getindex` takes.
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

# Reducing over the cell axis collapses it to one element, and no single cell
# id names that element — the same answer `Categorical` gives, for the same
# reason.
Lookups.reducelookup(::CellLookup) = Lookups.NoLookup(Base.OneTo(1))

# Asked directly rather than left to the generic search over the values: the
# window search is the exact membership test, and it is O(log #windows) where
# the generic is O(log #cells) of `cellindex` calls.
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
# Three questions, one answer each. Point and id resolve to a position;
# a region resolves to the positions its coverage names. Every method is typed
# on the selector's VALUE as well as on the lookup, because DimensionalData
# reads a `Tuple`-valued selector as a pair of interval endpoints and a
# `Vector`-valued one as an elementwise map — both of which a bare
# `(::CellLookup, ::Contains)` method would be ambiguous with.
#
# Each of them is one line, because each of them is one `CellVector` verb.
# ===========================================================================

"""
    Covering(target)

The region selector for a [`Cells`](@ref) dimension: every stored cell that
[`MultiOrderCoverage`](@ref) names for `target`, at the lookup's own level.

```julia
A[Cells(Covering(county))]          # a GeoInterface geometry
A[Cells(Covering(extent))]          # a lon/lat Extents.Extent
A[Cells(Covering(cap))]             # a GO.UnitSpherical.SphericalCap
```

`target` is anything [`query`](@ref) accepts. The result is the intersection of
the coverage's leaf expansion with the lookup, in ascending position order, and
the view it produces carries a [`CellLookup`](@ref) again — so what a subset
*stores* is windows, never an id vector.

This is [`covering`](@ref) under a `DimensionalData` hat: outside a cube the
same selection is `covering(cv, target)`, which answers with a
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

# A point names a cell before it names a position: `cellat` on the leaf grid,
# then the same window search every other selector ends in — which is exactly
# what the three-argument `cellposition` on a `CellVector` is.
Lookups.selectindices(lk::CellLookup, sel::Lookups.Contains{<:Tuple{Real,Real}}; kw...) =
    _found(lk, cellposition(parent(lk), Lookups.val(sel)...), sel)

Lookups.selectindices(lk::CellLookup, sel::Lookups.At{<:Tuple{Real,Real}}; kw...) =
    _found(lk, cellposition(parent(lk), Lookups.val(sel)...), sel)

_found(lk::CellLookup, k::Int, sel) = k
_found(lk::CellLookup, ::Nothing, sel) = throw(Lookups.SelectorError(lk, sel))

end # module CellLookups
