# A DimensionalData cell axis. Every cell in it sits at one leaf level, but the
# stored form is the multi-order set the region came back as, never the expanded
# id vector: memory is O(#coverage entries), the ids are computed on demand.

"""
    CellLookups

The DimensionalData layer: [`CellLookup`](@ref), the [`Cells`](@ref) dimension,
and the [`Covering`](@ref) selector.

A [`CellLookup`](@ref) is a one-dimensional `DimensionalData` lookup over cell
ids at a single level. What it stores is a set of **leaf position windows** —
sorted, disjoint intervals (or, where intervals are unavailable, a sorted list)
of positions in `levelgrid(sys, leaf)`. Its logical content is their
concatenation, and every operation is arithmetic over that concatenation:
`length` sums the window lengths, `lk[k]` binary-searches the cumulative
lengths and resolves one `cellindex`, [`cellposition`](@ref) runs the inverse.
Nothing is materialised.
"""
module CellLookups

import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, ncells, cellindex, cellposition, cellat, level, system,
    levelgrid, cellindextype, has_sorted_subtrees, descendants, query
import ..DiscreteGlobalGrids: Helpers
import ..DiscreteGlobalGrids.Fallbacks: PartialGrid, SubtreeIds,
    MultiOrderCoverage, MultiOrderCellSet, level_ranges

import SmallCollections
import DimensionalData as DD
import DimensionalData: Dimensions, Lookups

# ===========================================================================
# Windows: the stored form
#
# Two shapes, chosen by the system's `has_sorted_subtrees` trait and by how
# well an explicit position list compresses. Both answer the same two
# questions, which is the whole of what the lookup needs:
#
#   `leafposition(w, k)`   concatenation position -> leaf grid position
#   `windowposition(w, p)` leaf grid position -> concatenation position or `nothing`
# ===========================================================================

# Ranges are held as three parallel `Int` vectors rather than as a
# `Vector{UnitRange{Int}}` plus offsets: same three words per window, and both
# searches become a plain `searchsortedfirst` over a stored vector instead of a
# `last.(ranges)` allocation on every lookup.
struct RangeWindows
    starts::Vector{Int}
    stops::Vector{Int}
    offsets::Vector{Int}     # offsets[j] = cells in windows 1:j
end

struct PositionWindows
    positions::Vector{Int}   # sorted, strictly ascending
end

const CellWindows = Union{RangeWindows,PositionWindows}

Base.length(w::RangeWindows) = isempty(w.offsets) ? 0 : @inbounds w.offsets[end]
Base.length(w::PositionWindows) = length(w.positions)

nwindows(w::RangeWindows) = length(w.starts)
nwindows(w::PositionWindows) = length(w.positions)

@inline function leafposition(w::RangeWindows, k::Int)
    j = searchsortedfirst(w.offsets, k)
    base = j == 1 ? 0 : @inbounds w.offsets[j-1]
    return @inbounds w.starts[j] + (k - base) - 1
end

@inline leafposition(w::PositionWindows, k::Int) = @inbounds w.positions[k]

@inline function windowposition(w::RangeWindows, p::Int)
    j = searchsortedfirst(w.stops, p)
    j <= length(w.stops) || return nothing
    @inbounds start = w.starts[j]
    p >= start || return nothing
    base = j == 1 ? 0 : @inbounds w.offsets[j-1]
    return base + (p - start) + 1
end

@inline function windowposition(w::PositionWindows, p::Int)
    j = searchsortedfirst(w.positions, p)
    (j <= length(w.positions) && @inbounds(w.positions[j]) == p) || return nothing
    return j
end

leafpositions(w::CellWindows) = (leafposition(w, k) for k in 1:length(w))

# Two windowings are equal when they imply the same leaf positions, whatever
# shape each of them chose to store: the compression heuristic below is an
# implementation detail and must not be observable through `==`.
Base.:(==)(a::RangeWindows, b::RangeWindows) =
    a.starts == b.starts && a.stops == b.stops
Base.:(==)(a::CellWindows, b::CellWindows) =
    length(a) == length(b) && all(((x, y),) -> x == y, zip(leafpositions(a), leafpositions(b)))

_empty_windows() = RangeWindows(Int[], Int[], Int[])

function _range_windows(ranges)
    starts, stops, offsets = Int[], Int[], Int[]
    total = 0
    for r in ranges
        isempty(r) && continue
        push!(starts, first(r))
        push!(stops, last(r))
        total += length(r)
        push!(offsets, total)
    end
    return RangeWindows(starts, stops, offsets)
end

# Run-compress a sorted position list, then keep whichever shape is smaller.
# A run costs three words and a bare position one, so ranges win exactly when
# they are at most a third of the positions; the factor is stated here rather
# than tuned, because both shapes answer identically and only memory moves.
function _windows(positions::Vector{Int})
    isempty(positions) && return _empty_windows()
    runs = 1
    for i in 2:length(positions)
        @inbounds positions[i] > positions[i-1] || throw(ArgumentError(
            "leaf positions must be strictly ascending, got $(positions[i-1]) " *
            "then $(positions[i]) at index $i"))
        @inbounds positions[i] == positions[i-1] + 1 || (runs += 1)
    end
    3 * runs <= length(positions) || return PositionWindows(positions)
    starts = Vector{Int}(undef, runs)
    stops = Vector{Int}(undef, runs)
    offsets = Vector{Int}(undef, runs)
    j = 1
    @inbounds starts[1] = positions[1]
    for i in 2:length(positions)
        @inbounds if positions[i] != positions[i-1] + 1
            stops[j] = positions[i-1]
            offsets[j] = i - 1
            j += 1
            starts[j] = positions[i]
        end
    end
    @inbounds stops[runs] = positions[end]
    @inbounds offsets[runs] = length(positions)
    return RangeWindows(starts, stops, offsets)
end

# ===========================================================================
# The lazy id vector
# ===========================================================================

"""
    CellIds(windows, grid) <: AbstractVector

The ids a windowing names, resolved one at a time through `cellindex(grid, ·)`.
This is what a [`CellLookup`](@ref) holds and what [`PartialGrid`](@ref) reads
when a lookup is handed to the regridder, so neither ever materialises.
"""
struct CellIds{W<:CellWindows,G<:AbstractGrid,ID} <: AbstractVector{ID}
    windows::W
    grid::G
end

CellIds(windows::CellWindows, grid::AbstractGrid) =
    CellIds{typeof(windows),typeof(grid),cellindextype(system(grid))}(windows, grid)

Base.size(v::CellIds) = (length(v.windows),)
Base.IndexStyle(::Type{<:CellIds}) = Base.IndexLinear()

Base.@propagate_inbounds function Base.getindex(v::CellIds, k::Int)
    @boundscheck checkbounds(v, k)
    return cellindex(v.grid, leafposition(v.windows, k))
end

# A complete level grid's positions ascend in canonical id order, so windows —
# which are ascending positions by construction — name ascending ids. The O(n)
# verification `PartialGrid` runs on an arbitrary vector has nothing to find.
Helpers.strictly_increasing(::CellIds) = true

# ===========================================================================
# The lookup
# ===========================================================================

"""
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
cells, `lk[k]` is the `k`th of them, `collect(lk)` is the vector itself. What
is *stored* is the set's expansion to sorted, disjoint position windows at the
leaf level ([`level_ranges`](@ref)), so memory is O(#entries in the set) rather
than O(#leaf cells) — on a Switzerland-sized region at IGEO7 level 9 that is
some hundreds of words standing for tens of thousands of cells. `lk[k]`
binary-searches the windows' cumulative lengths and resolves one `cellindex`;
`cellposition(lk, c)` runs the inverse.

`Base.parent` returns the lookup's VALUES, as `DimensionalData` requires: the
lazy cell id vector, which is O(#windows) and materialises nothing.
[`cellset`](@ref) returns the backing — the set, or the grid — for running a
second coverage operation against without unpacking the lookup.

# The three ways in

  - a [`MultiOrderCellSet`](@ref), the compressed form, optionally re-expanded
    to a deeper `level` than the set's own reference level;
  - `levelgrid(sys, l)`, a whole level, which is one window;
  - a [`PartialGrid`](@ref), an arbitrary ascending subset, which is that
    subset's positions — one window when the subset is a subtree, and the
    explicit list when it is scattered.

The last two are the degenerate cases of the first, and answer every method
below identically.

# Selectors

```julia
A[Cells(DimensionalData.At(c))]              # a typed cell id
A[Cells(DimensionalData.Contains(8.0, 46.5))] # a lon/lat point, through `cellat`
A[Cells(Covering(polygon))]                   # a region, through `MultiOrderCoverage`
```

`At` and `Contains` resolve to one position; [`Covering`](@ref) to the
positions of every stored cell the region's coverage names, and the view it
produces carries a `CellLookup` again.

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
    still holds.

    What that costs is the construction, which walks the leaves — O(#leaf cells
    in the subset) in time and in transient memory, against O(#entries) for the
    windowed path. The *stored* form is whatever the compression finds, which
    on a coverage of a connected region is usually a handful of windows and in
    the worst case is one position per cell. Either way it is bounded by the
    subset, never by the globe.

    The alternative was to refuse: it would make A5 the one system on which a
    cube cannot carry a cell axis, for a property that
    `PartialGrid(sys, cell, level)` already degrades gracefully around.

    One consequence is inherited rather than introduced. An A5 cell's
    descendants need not lie inside its own footprint, so expanding a coverage
    names leaves the target does not touch — most visibly inside a hole. A
    [`Covering`](@ref) selection on A5 is therefore a superset of the cells that
    meet the region, by the same margin the refinement itself is; see
    [`MultiOrderCoverage`](@ref).
"""
struct CellLookup{ID,B,V<:CellIds} <: Lookups.Lookup{ID,1}
    backing::B
    ids::V
    level::Int
end

function CellLookup(backing, ids::CellIds{<:Any,<:Any,ID}, l::Integer) where {ID}
    return CellLookup{ID,typeof(backing),typeof(ids)}(backing, ids, Int(l))
end

# The keyword shadows the `level` function, so the work is a positional `Int`
# one call in — the `_multi_order` pattern.
CellLookup(set::MultiOrderCellSet; level::Integer=set.reference_level) =
    _celllookup(set, Int(level))

function _celllookup(set::MultiOrderCellSet, l::Int)
    sys = system(set)
    grid = levelgrid(sys, l)
    windows = has_sorted_subtrees(sys) ?
              _range_windows(level_ranges(set, l)) :
              _windows(_selection_positions(set, grid, l))
    return CellLookup(set, CellIds(windows, grid), l)
end

# The selection-mode expansion: `descendants` names the leaves as a list, and
# their positions are sorted afterwards because the concatenation of two
# siblings' subtrees need not be ascending where subtrees are not sorted.
function _selection_positions(set::MultiOrderCellSet, grid::AbstractGrid, l::Int)
    sys = system(set)
    out = Int[]
    for c in set
        level(c) <= l || throw(ArgumentError(
            "cannot expand to level $l: the set contains a level-$(level(c)) cell"))
        for d in descendants(sys, c, l)
            p = cellposition(grid, d)
            p === nothing && throw(ArgumentError(
                "descendant $d of $c is not a cell of levelgrid($sys, $l)"))
            push!(out, p)
        end
    end
    return sort!(out)
end

function CellLookup(grid::AbstractGrid)
    sys = system(grid)
    l = level(grid)
    complete = levelgrid(sys, l)
    return CellLookup(grid, CellIds(_grid_windows(grid, complete), complete), l)
end

CellLookup(lk::CellLookup) = lk

# A grid holding every cell of its level is positions `1:ncells` by the
# completeness contract, whatever wrapper it wears; only a proper subset has to
# be walked.
function _grid_windows(grid::AbstractGrid, complete::AbstractGrid)
    n = ncells(grid)
    n == ncells(complete) && return _range_windows((1:n,))
    return _windows(_grid_positions(grid, complete))
end

# A rooted subset over sorted subtrees is one window and knows it, so the walk
# below is skipped rather than performed and compressed.
function _grid_windows(grid::PartialGrid{<:Any,<:SubtreeIds}, complete::AbstractGrid)
    ids = grid.ids
    ids.n == 0 && return _empty_windows()
    return _range_windows((ids.first:(ids.first+ids.n-1),))
end

function _grid_positions(grid::AbstractGrid, complete::AbstractGrid)
    out = Vector{Int}(undef, ncells(grid))
    for i in eachindex(out)
        c = cellindex(grid, i)
        p = cellposition(complete, c)
        p === nothing && throw(ArgumentError(
            "$c is not a cell of levelgrid($(system(grid)), $(level(grid)))"))
        out[i] = p
    end
    return out
end

# A subset of an existing lookup keeps the same leaf grid and level; its backing
# becomes the `PartialGrid` those windows describe, which is O(1) to build
# because the ids stay lazy.
function _derive(lk::CellLookup, windows::CellWindows)
    ids = CellIds(windows, lk.ids.grid)
    return CellLookup(PartialGrid(system(lk), lk.level, ids), ids, lk.level)
end

windows(lk::CellLookup) = lk.ids.windows

# --- the collection surface ------------------------------------------------

# `parent` is the VALUES — the lazy id vector — and nothing else. That is
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
Base.parent(lk::CellLookup) = lk.ids
Base.IndexStyle(::Type{<:CellLookup}) = Base.IndexLinear()

# `first`/`last` on an empty lookup are a `BoundsError`, so the one caller that
# asks about an empty axis rather than indexing it answers `nothing` instead.
Lookups.bounds(lk::CellLookup) = isempty(lk) ? (nothing, nothing) : (first(lk), last(lk))

Base.@propagate_inbounds Base.getindex(lk::CellLookup, k::Int) = lk.ids[k]
Base.@propagate_inbounds Base.getindex(lk::CellLookup, k::CartesianIndex{1}) = lk.ids[k[1]]

for f in (:getindex, :view, :dotview)
    @eval Base.$f(lk::CellLookup, ::Colon) = lk
    @eval Base.$f(lk::CellLookup, i::AbstractArray{<:Integer}) = _subset(lk, i)
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

# Ascending indices keep the windowed form; anything else — a permutation, a
# repeat, a reversal — is not a set of leaf windows and is answered with the
# ordinary DimensionalData lookup that can hold it. That case materialises, and
# only that case.
#
# A `Bool` array is a mask rather than a list of ones, and `Bool <: Integer`, so
# the branch is here rather than in a second signature: dispatching on
# `AbstractArray{Bool}` would re-open the ambiguity the method above closes.
_subset(lk::CellLookup, mask::AbstractArray{Bool}) = _subset(lk, findall(mask))

function _subset(lk::CellLookup, idx::AbstractArray{<:Integer})
    n = length(lk)
    positions = Vector{Int}(undef, length(idx))
    ascending = true
    for (j, k) in enumerate(idx)
        1 <= k <= n || throw(BoundsError(lk, k))
        positions[j] = leafposition(windows(lk), Int(k))
        j > 1 && positions[j] <= positions[j-1] && (ascending = false)
    end
    ascending || return Lookups.Categorical([lk[Int(k)] for k in idx];
        order=Lookups.Unordered())
    return _derive(lk, _windows(positions))
end

# --- what the lookup is, in this package's own vocabulary ------------------

"""
    cellset(lk::CellLookup)

The thing the lookup was built from — a [`MultiOrderCellSet`](@ref), or the grid
— for running a second coverage operation against without unpacking the lookup.

`Base.parent(lk)` is deliberately NOT this: it is the lookup's VALUES, the lazy
cell id vector, because that is what `DimensionalData` derives some thirty
`Lookup` methods from. A subset produced by indexing or by a selector reports
the [`PartialGrid`](@ref) describing it.
"""
cellset(lk::CellLookup) = lk.backing

"""
    system(lk::CellLookup)

The grid system the lookup's cells are named in.
"""
system(lk::CellLookup) = system(lk.ids.grid)

"""
    level(lk::CellLookup) -> Int

The one level every cell in the lookup sits at.
"""
level(lk::CellLookup) = lk.level

"""
    cellposition(lk::CellLookup, c::AbstractCellIndex) -> Union{Int,Nothing}

Position of cell `c` in the lookup, or `nothing` when the lookup does not hold
it — including when `c` is at another level. The inverse of `lk[k]`, and the
half of the bijection every selector ends at.
"""
function cellposition(lk::CellLookup, c::AbstractCellIndex)
    p = cellposition(lk.ids.grid, c)
    p === nothing && return nothing
    return windowposition(windows(lk), p)
end

"""
    PartialGrid(lk::CellLookup) -> PartialGrid

The lookup read as a grid: position `k` of the grid is position `k` of the
lookup, so a `Regridder` built on it lines up with a cube over the lookup's
axis without a permutation. O(1) — the ids stay lazy.
"""
PartialGrid(lk::CellLookup) = PartialGrid(system(lk), lk.level, lk.ids)

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
    (data === nothing || data === lk || data === lk.ids) && return lk
    return _rebuild(lk, data)
end

_rebuild(lk::CellLookup, ids::CellIds) = _derive(lk, ids.windows)

function _rebuild(lk::CellLookup, ids::AbstractVector{<:AbstractCellIndex})
    grid = lk.ids.grid
    positions = Vector{Int}(undef, length(ids))
    ascending = true
    for (j, c) in enumerate(ids)
        p = cellposition(grid, c)
        p === nothing && throw(ArgumentError(
            "$c is not a cell of levelgrid($(system(lk)), $(lk.level)), so it " *
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
    _at_point(lk, Lookups.val(sel)...) !== nothing

Dimensions.format(lk::CellLookup, ::Type, values, axis::AbstractRange) = lk

Base.:(==)(a::CellLookup, b::CellLookup) =
    system(a) == system(b) && a.level == b.level && windows(a) == windows(b)

function Base.show(io::IO, lk::CellLookup)
    print(io, "CellLookup(", typeof(system(lk)).name.name, ", level=", lk.level,
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

!!! note "What it stores and what it spends are different questions"
    Selecting walks the coverage's leaves one at a time and builds the vector of
    matching positions, because that vector is what indexes the cube's data.
    Both are O(#leaf cells the coverage names), transiently, however compact the
    two ends are: on a level-12 axis that is hundreds of megabytes passing
    through. The storage claim is about the lookup, not about the selection —
    select at the level you mean to read at.

Coverage means *covering*: the selection is a superset of the cells that meet
`target`, by whatever margin the system's refinement is non-congruent. See
[`MultiOrderCoverage`](@ref) for the size of that margin per system.
"""
struct Covering{T} <: Lookups.ArraySelector{T}
    val::T
end

Base.show(io::IO, sel::Covering) =
    print(io, "Covering(", typeof(sel.val).name.name, ")")

Lookups.selectindices(lk::CellLookup, sel::Covering; kw...) =
    _covering(lk, Lookups.val(sel))

Lookups.selectindices(lk::CellLookup, sel::Covering{<:AbstractVector}; kw...) =
    _covering(lk, Lookups.val(sel))

function _covering(lk::CellLookup, target)
    sys = system(lk)
    set = query(sys, MultiOrderCoverage(target); level=lk.level)
    out = Int[]
    w = windows(lk)
    _each_leaf_position(sys, set, lk.ids.grid, lk.level) do p
        k = windowposition(w, p)
        k === nothing || push!(out, k)
    end
    return issorted(out) ? out : sort!(out)
end

# The one place the two expansions are chosen between: ranges where subtrees are
# sorted, `descendants` where they are not. Both visit each leaf once.
function _each_leaf_position(f, sys, set::MultiOrderCellSet, grid::AbstractGrid, l::Int)
    if has_sorted_subtrees(sys)
        for r in level_ranges(set, l), p in r
            f(p)
        end
    else
        for p in _selection_positions(set, grid, l)
            f(p)
        end
    end
    return nothing
end

Lookups.selectindices(lk::CellLookup, sel::Lookups.At{<:AbstractCellIndex}; kw...) =
    _found(lk, cellposition(lk, Lookups.val(sel)), sel)

Lookups.selectindices(lk::CellLookup, sel::Lookups.Contains{<:AbstractCellIndex}; kw...) =
    _found(lk, cellposition(lk, Lookups.val(sel)), sel)

Lookups.selectindices(lk::CellLookup, sel::Lookups.Contains{<:Tuple{Real,Real}}; kw...) =
    _found(lk, _at_point(lk, Lookups.val(sel)...), sel)

Lookups.selectindices(lk::CellLookup, sel::Lookups.At{<:Tuple{Real,Real}}; kw...) =
    _found(lk, _at_point(lk, Lookups.val(sel)...), sel)

# A point names a cell before it names a position: `cellat` on the leaf grid,
# then the same window search every other selector ends in.
function _at_point(lk::CellLookup, lon::Real, lat::Real)
    c = cellat(lk.ids.grid, Float64(lon), Float64(lat))
    c === nothing && return nothing
    return cellposition(lk, c)
end

_found(lk::CellLookup, k::Int, sel) = k
_found(lk::CellLookup, ::Nothing, sel) = throw(Lookups.SelectorError(lk, sel))

end # module CellLookups
