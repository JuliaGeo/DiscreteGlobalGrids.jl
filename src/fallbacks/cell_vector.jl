# ---------------------------------------------------------------------------
# `CellVector` — the compressed cell collection
#
# A strictly ascending run of cells at ONE level of one system, stored as the
# leaf position WINDOWS they occupy rather than as their ids. It is an
# `AbstractVector` of the ids: the compression is the storage, never the
# semantics.
#
# This is the substrate the DimensionalData layer (`src/dimensionaldata.jl`)
# wraps, and it is deliberately DimensionalData-free so that everything else —
# regridding, chunking, a bare `PartialGrid` workflow — gets the same
# compression without importing a cube library. `CellLookup` is this type
# wearing a `DimensionalData.Lookup` hat and nothing more; every verb it
# answers is one of the verbs below.
#
# The two shapes below answer the same two questions, which is the whole of
# what the vector needs:
#
#   `leafposition(w, k)`   concatenation position -> leaf grid position
#   `windowposition(w, p)` leaf grid position -> concatenation position or `nothing`
# ---------------------------------------------------------------------------

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

# Every window set in this file is built in CANONICAL form: runs are maximal, so
# two windowings name the same leaf positions if and only if their runs match
# one for one. That is what lets the `RangeWindows` pair below decide equality
# on three vectors instead of on every position, and it is why the set
# operations normalise (`_windows_from_intervals`) rather than emitting whatever
# split their merge happened to produce.
Base.:(==)(a::RangeWindows, b::RangeWindows) =
    a.starts == b.starts && a.stops == b.stops

# Across the two shapes the runs cannot be compared directly, so the fallback
# says the same thing the slow way: the compression heuristic is an
# implementation detail and must not be observable through `==`.
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

# --- the two shapes, read as intervals -------------------------------------
#
# The set operations are interval arithmetic, and interval arithmetic does not
# care which shape stored the interval. A bare position is a one-cell interval,
# so both shapes answer this and the merges below are O(#windows) on either.

intervals(w::RangeWindows) =
    [(@inbounds(w.starts[j]), @inbounds(w.stops[j])) for j in eachindex(w.starts)]

intervals(w::PositionWindows) = [(p, p) for p in w.positions]

# Sorted, disjoint, maximal intervals -> windows, in the canonical form the
# equality above assumes. The heuristic is `_windows`', restated over intervals
# so that a set operation never has to expand to positions to choose a shape.
function _windows_from_intervals(ivs::Vector{Tuple{Int,Int}})
    isempty(ivs) && return _empty_windows()
    merged = Tuple{Int,Int}[]
    for (lo, hi) in ivs
        if !isempty(merged) && lo <= last(merged)[2] + 1
            merged[end] = (last(merged)[1], max(last(merged)[2], hi))
        else
            push!(merged, (lo, hi))
        end
    end
    total = sum(hi - lo + 1 for (lo, hi) in merged)
    if 3 * length(merged) <= total
        n = length(merged)
        starts, stops, offsets = Vector{Int}(undef, n), Vector{Int}(undef, n), Vector{Int}(undef, n)
        acc = 0
        for (j, (lo, hi)) in enumerate(merged)
            starts[j] = lo
            stops[j] = hi
            acc += hi - lo + 1
            offsets[j] = acc
        end
        return RangeWindows(starts, stops, offsets)
    end
    positions = Vector{Int}(undef, total)
    k = 0
    for (lo, hi) in merged, p in lo:hi
        positions[k+=1] = p
    end
    return PositionWindows(positions)
end

# ===========================================================================
# The vector
# ===========================================================================

"""
    CellVector(set::MultiOrderCellSet; level = set's reference level)
    CellVector(grid::AbstractGrid)
    CellVector(sys, level, ids::AbstractVector)

An immutable, lazy `AbstractVector` of **strictly ascending cell ids at one
level** of one hierarchical system, stored as the leaf position *windows* they
occupy rather than as the ids themselves.

Semantically `cv` **is** the id vector: `length(cv)` is the number of cells,
`cv[k]` is the `k`th of them, `collect(cv)` is the vector itself. What is
*stored* is sorted, disjoint intervals (or, where intervals are unavailable, a
sorted position list) in `levelgrid(system(cv), level(cv))`, so memory is
O(#windows) rather than O(#cells) — on a Switzerland-sized region at IGEO7
level 9 that is 666 windows standing for 60,861 cells, and the *same* 666
windows for the level-12 re-expansion of the same coverage, which names
20,875,323. `cv[k]` binary-searches the windows' cumulative lengths and
resolves one `cellindex`; [`cellposition`](@ref) runs the inverse. Nothing is
materialised.

This is the compression itself, with no cube library attached.
[`CellLookup`](@ref) is this type wearing a `DimensionalData.Lookup` hat, and
every method it answers delegates to one of the verbs here — so regridding,
chunking and plain-array code get the same laziness without importing
`DimensionalData`.

# The three ways in

  - a [`MultiOrderCellSet`](@ref), the compressed coverage, optionally
    re-expanded to a deeper `level` than the set's own reference level;
  - an [`AbstractGrid`](@ref) — `levelgrid(sys, l)`, a whole level, which is one
    window; or a [`PartialGrid`](@ref), which is that subset's positions (one
    window when the subset is a subtree, the explicit list when it is
    scattered);
  - `sys, level, ids` — an explicit **strictly ascending** id vector, validated
    and run-compressed.

`CellVector(cv)` is the identity. The last two forms are the degenerate cases
of the first and answer every method below identically.

# The verbs

```julia
cv[k]                          # the kth cell id
cellposition(cv, c)            # its inverse: Int, or `nothing`
c in cv                        # membership, O(log #windows)
cellposition(cv, lon, lat)     # the position of the cell a point falls in
cellat(cv, lon, lat)           # that cell's id instead
covering(cv, polygon)          # the sub-vector a region's coverage names
intersect(cv, other)           # two vectors at the same level, O(#windows)
PartialGrid(cv)                # read as a grid, O(1) — the regridding handshake
cellset(cv)                    # what it was built from
```

# Memory, and where it is spent

`Base.summarysize(cv)` is O(#windows) plus whatever [`cellset`](@ref) points
at, and re-expanding one coverage to a deeper level does not move it — that
invariance is the reason the type exists. Reading is O(log #windows) and
allocation-free.

Construction is a different question. Where [`has_sorted_subtrees`](@ref)
holds, [`level_ranges`](@ref) is already the answer and building is
O(#entries). Where it does not (A5), a cell's descendants are not one interval
of their level, so the vector is built by SELECTION: `descendants` names the
leaves, they are resolved to positions and sorted, and the result is
run-compressed like any other position list. Every method above is unchanged
and every law still holds; what it costs is the *construction*, which walks the
leaves. The stored form is still whatever the compression finds — usually a
handful of windows on a connected region, one position per cell in the worst
case, and bounded by the subset either way.

One consequence there is inherited rather than introduced: an A5 cell's
descendants need not lie inside its own footprint, so expanding a coverage
names leaves the target does not touch — most visibly inside a hole. A
[`covering`](@ref) subset on A5 is therefore a superset of the cells that meet
the region, by the same margin the refinement itself is; see
[`MultiOrderCoverage`](@ref).
"""
struct CellVector{ID,W<:CellWindows,G<:AbstractGrid,B} <: AbstractVector{ID}
    windows::W
    grid::G                  # `levelgrid(system, level)` — every `cv[k]` reads it
    backing::B               # what it was built from, or `nothing` when derived
    level::Int
end

function CellVector(windows::CellWindows, grid::AbstractGrid, backing, l::Integer)
    ID = cellindextype(system(grid))
    return CellVector{ID,typeof(windows),typeof(grid),typeof(backing)}(
        windows, grid, backing, Int(l))
end

# The keyword shadows the `level` function, so the work is a positional `Int`
# one call in — the `_multi_order` pattern.
CellVector(set::MultiOrderCellSet; level::Integer=set.reference_level) =
    _cellvector(set, Int(level))

function _cellvector(set::MultiOrderCellSet, l::Int)
    sys = system(set)
    grid = levelgrid(sys, l)
    w = has_sorted_subtrees(sys) ?
        _range_windows(level_ranges(set, l)) :
        _windows(_selection_positions(set, grid, l))
    return CellVector(w, grid, set, l)
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

function CellVector(grid::AbstractGrid)
    sys = system(grid)
    l = level(grid)
    complete = levelgrid(sys, l)
    return CellVector(_grid_windows(grid, complete), complete, grid, l)
end

CellVector(sys::AbstractHierarchicalGridSystem, l::Integer, ids::AbstractVector) =
    CellVector(PartialGrid(sys, l, ids))

CellVector(cv::CellVector) = cv

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

# The round trip. `PartialGrid(cv)` is O(1) and so is coming back: the grid's
# ids ARE a `CellVector`, and its windows are already the answer.
_grid_windows(grid::PartialGrid{<:Any,<:CellVector}, complete::AbstractGrid) =
    grid.ids.windows

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

# A subset of an existing vector keeps the same leaf grid and level, and has no
# backing of its own: `cellset` answers such a vector with the `PartialGrid`
# describing it, which is O(1) to build because the ids stay lazy.
_derive(cv::CellVector, w::CellWindows) = CellVector(w, cv.grid, nothing, cv.level)

windows(cv::CellVector) = cv.windows

# --- the collection surface ------------------------------------------------

Base.size(cv::CellVector) = (length(cv.windows),)
Base.IndexStyle(::Type{<:CellVector}) = Base.IndexLinear()

Base.@propagate_inbounds function Base.getindex(cv::CellVector, k::Int)
    @boundscheck checkbounds(cv, k)
    return cellindex(cv.grid, leafposition(cv.windows, k))
end

# Immutable, so the whole-vector slice is the vector rather than a copy of it.
Base.getindex(cv::CellVector, ::Colon) = cv

# Ascending indices keep the windowed form; anything else — a permutation, a
# repeat, a reversal — is not a set of leaf windows and is answered with the
# ordinary `Vector` that can hold it. That case materialises, and only that
# case. [`CellLookup`](@ref) takes the same fork one level up.
#
# `AbstractVector`, not `AbstractArray`: indexing an array by a
# higher-dimensional index returns something of the INDEX's shape, and a window
# set has no shape to give back. Narrowing here lets `cv[matrix]` fall through
# to Base's generic, which answers with a matrix of ids — so the answer's shape
# no longer depends on whether the index happened to be ascending.
Base.getindex(cv::CellVector, idx::AbstractVector{<:Integer}) = _subset(cv, idx)

# SmallCollections' own `getindex(::AbstractVector, ::AbstractFixedOrSmall...)`
# is neither more nor less specific than the line above, and a neighbour list is
# exactly one of those vectors, so the tie is broken towards the same subset
# rather than left as an ambiguity for whoever indexes a cell vector by a halo.
Base.getindex(cv::CellVector,
    i::SmallCollections.AbstractFixedOrSmallOrPackedVector{<:Integer}) = _subset(cv, i)

# A `Bool` array is a mask rather than a list of ones, and `Bool <: Integer`, so
# the branch is here rather than in a second `getindex` signature: dispatching
# on `AbstractArray{Bool}` would re-open the ambiguity the method above closes.
#
# A mask names a position by its INDEX, so a mask of the wrong length is a
# bounds error rather than a shorter answer. `findall` cannot tell — it reports
# the `true`s it was given and nothing about the ones it was not — so the axes
# are checked before it runs, which is what Base's own logical indexing does.
# The signature stays `AbstractArray` so that a mask of the wrong RANK lands
# here and is caught by the same check.
function _subset(cv::CellVector, mask::AbstractArray{Bool})
    axes(mask) == axes(cv) || throw(BoundsError(cv, (mask,)))
    return _subset(cv, findall(mask))
end

function _subset(cv::CellVector, idx::AbstractArray{<:Integer})
    n = length(cv)
    positions = Vector{Int}(undef, length(idx))
    ascending = true
    for (j, k) in enumerate(idx)
        1 <= k <= n || throw(BoundsError(cv, (idx,)))
        positions[j] = leafposition(cv.windows, Int(k))
        j > 1 && positions[j] <= positions[j-1] && (ascending = false)
    end
    ascending || return [cv[Int(k)] for k in idx]
    return _derive(cv, _windows(positions))
end

# The ids ascend in canonical order by construction — the windows are ascending
# positions of a complete level grid — so the O(n) verification `PartialGrid`
# runs on an arbitrary vector has nothing to find.
Helpers.strictly_increasing(::CellVector) = true

# Same system, same level and the same windows name the same ids, so equality is
# O(#windows) instead of O(#cells) — and it agrees with the elementwise answer
# `AbstractVector` would give, because an id is self-describing about its level.
Base.:(==)(a::CellVector, b::CellVector) =
    system(a) == system(b) && a.level == b.level && a.windows == b.windows

# --- what the vector is, in this package's own vocabulary ------------------

"""
    system(cv::CellVector)

The grid system the vector's cells are named in.
"""
system(cv::CellVector) = system(cv.grid)

"""
    level(cv::CellVector) -> Int

The one level every cell in the vector sits at.
"""
level(cv::CellVector) = cv.level

"""
    cellset(cv::CellVector)
    cellset(lk::CellLookup)

The thing the collection was built from — a [`MultiOrderCellSet`](@ref), or the
grid — for running a second coverage operation against without unpacking it.

A collection *derived* from another one, by indexing or by [`covering`](@ref),
has no such origin and reports the [`PartialGrid`](@ref) describing it instead.

For a [`CellLookup`](@ref), `Base.parent` is deliberately NOT this: it is the
lookup's VALUES, the [`CellVector`](@ref), because that is what
`DimensionalData` derives some thirty `Lookup` methods from.
"""
cellset(cv::CellVector) = _origin(cv, cv.backing)

# Dispatched on the field's VALUE rather than on the fourth type parameter: a
# `CellVector{<:Any,<:Any,<:Any,Nothing}` signature is a subtype of the general
# one but is not *more specific* than it — the two slots left open have declared
# bounds the wildcard widens — so the general method wins and every derived
# vector answers `nothing`. Passing the field in leaves the choice to ordinary
# dispatch on `Nothing`, and the type parameter still constant-folds it.
_origin(::CellVector, backing) = backing
_origin(cv::CellVector, ::Nothing) = PartialGrid(cv)

"""
    cellposition(cv::CellVector, c::AbstractCellIndex) -> Union{Int,Nothing}
    cellposition(cv::CellVector, p::GO.UnitSphericalPoint) -> Union{Int,Nothing}
    cellposition(cv::CellVector, lon::Real, lat::Real) -> Union{Int,Nothing}

Position of a cell in the vector, or `nothing` when the vector does not hold it
— including when `c` is at another level. The inverse of `cv[k]`, and the half
of the bijection every selection ends at.

The point forms name the cell through [`cellat`](@ref) first and then answer
with its position, which is the one-call form of "where in my data does this
point land". Degrees for the `(lon, lat)` method, as everywhere else.
"""
function cellposition(cv::CellVector, c::AbstractCellIndex)
    p = cellposition(cv.grid, c)
    p === nothing && return nothing
    return windowposition(cv.windows, p)
end

function cellposition(cv::CellVector, p::GO.UnitSphericalPoint)
    c = cellat(cv.grid, p)
    c === nothing && return nothing
    return cellposition(cv, c)
end

cellposition(cv::CellVector, lon::Real, lat::Real) =
    cellposition(cv, unit_point(lon, lat))

"""
    cellat(cv::CellVector, p::GO.UnitSphericalPoint) -> Union{AbstractCellIndex,Nothing}
    cellat(cv::CellVector, lon::Real, lat::Real)

The cell **of `cv`** containing the point, or `nothing` when the point falls
outside the cells the vector holds. Same contract as [`cellat`](@ref) on a
grid, restricted to the subset: a point inside the level grid but outside this
vector answers `nothing` rather than naming a cell that is not here.

[`cellposition`](@ref) is the same question answered as a position.
"""
function cellat(cv::CellVector, p::GO.UnitSphericalPoint)
    c = cellat(cv.grid, p)
    c === nothing && return nothing
    return cellposition(cv, c) === nothing ? nothing : c
end

cellat(cv::CellVector, lon::Real, lat::Real) = cellat(cv, unit_point(lon, lat))

# Membership is the window search, which is the exact test and O(log #windows)
# where the generic scan over the values is O(#cells) of `cellindex` calls.
Base.in(c::AbstractCellIndex, cv::CellVector) = cellposition(cv, c) !== nothing

"""
    PartialGrid(cv::CellVector) -> PartialGrid

The vector read as a grid: position `k` of the grid is position `k` of the
vector, so a `Regridder` built on it lines up with data laid out against the
vector without a permutation. **O(1)** — the ids stay lazy, and this is the
handshake that lets a regridder consume a compressed coverage directly.

The grid keeps the windows and drops the backing: [`cellset`](@ref)'s origin is
provenance the grid never reads, and a regridder should not hold a coverage set
alive behind it. `CellVector(PartialGrid(cv))` comes back in O(1) either way.
"""
PartialGrid(cv::CellVector) = PartialGrid(system(cv), cv.level, _bare(cv, cv.backing))

# Dispatched on the backing's value, for the reason `_origin` gives above.
_bare(cv::CellVector, backing) = CellVector(cv.windows, cv.grid, nothing, cv.level)
_bare(cv::CellVector, ::Nothing) = cv

# --- selecting a region ----------------------------------------------------

"""
    covering(cv::CellVector, target) -> CellVector

The cells of `cv` that a [`MultiOrderCoverage`](@ref) of `target` names, at
`cv`'s own level — the region selector, as a `CellVector` again, so what a
subset *stores* is windows rather than an id vector.

`target` is anything [`query`](@ref) accepts: a GeoInterface geometry, an
`Extents.Extent` in lon/lat degrees, or a `GO.UnitSpherical.SphericalCap`.

```julia
cv    = CellVector(query(sys, MultiOrderCoverage(canton); level = 9))
basin = covering(cv, watershed)          # a CellVector again
data[covering_positions(cv, watershed)]  # the same selection as positions
```

[`covering_positions`](@ref) is the position-space form, for indexing a data
array laid out against `cv`. This is what the `DimensionalData` selector
[`Covering`](@ref) is spelled as outside `DimensionalData`.

!!! note "What it stores and what it spends are different questions"
    Selecting walks the coverage's leaves one at a time, because that walk is
    what decides which of them `cv` holds. It is O(#leaf cells the coverage
    names), transiently, however compact the two ends are: on a level-12 vector
    that is hundreds of megabytes passing through. The storage claim is about
    the result, not about the selection — select at the level you mean to read
    at.

Coverage means *covering*: the result is a superset of the cells that meet
`target`, by whatever margin the system's refinement is non-congruent. See
[`MultiOrderCoverage`](@ref) for the size of that margin per system.
"""
covering(cv::CellVector, target) =
    _derive(cv, _windows(_covering_leafpositions(cv, target)))

"""
    covering_positions(cv::CellVector, target) -> Vector{Int}

The positions in `cv` of the cells [`covering`](@ref) selects, ascending — for
indexing a data array laid out against `cv` without building the sub-vector.

`covering(cv, target)` and `cv[covering_positions(cv, target)]` name the same
cells; this form is the one a cube's `getindex` needs, and is what the
[`Covering`](@ref) selector resolves to.
"""
function covering_positions(cv::CellVector, target)
    out = Int[]
    w = cv.windows
    _each_leaf_position(cv, target) do p
        k = windowposition(w, p)
        k === nothing || push!(out, k)
    end
    return issorted(out) ? out : sort!(out)
end

function _covering_leafpositions(cv::CellVector, target)
    out = Int[]
    w = cv.windows
    _each_leaf_position(cv, target) do p
        windowposition(w, p) === nothing || push!(out, p)
    end
    return issorted(out) ? out : sort!(out)
end

# The one place the two expansions are chosen between: ranges where subtrees are
# sorted, `descendants` where they are not. Both visit each leaf once.
function _each_leaf_position(f, cv::CellVector, target)
    sys = system(cv)
    set = query(sys, MultiOrderCoverage(target); level=cv.level)
    if has_sorted_subtrees(sys)
        for r in level_ranges(set, cv.level), p in r
            f(p)
        end
    else
        for p in _selection_positions(set, cv.grid, cv.level)
            f(p)
        end
    end
    return nothing
end

# --- set arithmetic over the windows ---------------------------------------
#
# `intersect` and `issubset` are Base's, with Base's meaning: the first keeps
# the left operand's order, which for an ascending vector IS ascending, and the
# second is a pure set question. Both are answered here in O(#windows) instead
# of O(#cells).
#
# They differ on operands from different systems or levels, and the asymmetry is
# Base's rather than this file's: `issubset` is a question every pair of sets can
# answer — the empty vector is a subset of anything, and across mismatched
# domains nothing else is — while `intersect` has to RETURN a cell vector, and
# there is no system or level to build one in. So the first answers and the
# second throws.
#
# `union` is deliberately absent. Base's `union` maintains first-appearance
# order, which for two ascending vectors is NOT ascending, so an ascending
# answer would be a different function wearing Base's name. Merge two coverages
# by querying their union, or by `CellVector(sys, l, sort(vcat(a, b)))`.

function _same_space(a::CellVector, b::CellVector, verb::AbstractString)
    system(a) == system(b) || throw(ArgumentError(
        "cannot $verb cell vectors from $(typeof(system(a))) and $(typeof(system(b)))"))
    a.level == b.level || throw(ArgumentError(
        "cannot $verb cell vectors at levels $(a.level) and $(b.level)"))
    return nothing
end

function Base.intersect(a::CellVector, b::CellVector)
    _same_space(a, b, "intersect")
    return _derive(a, _windows_from_intervals(_intersect_intervals(a.windows, b.windows)))
end

function _intersect_intervals(x::CellWindows, y::CellWindows)
    A, B = intervals(x), intervals(y)
    out = Tuple{Int,Int}[]
    i = j = 1
    while i <= length(A) && j <= length(B)
        lo = max(A[i][1], B[j][1])
        hi = min(A[i][2], B[j][2])
        lo <= hi && push!(out, (lo, hi))
        A[i][2] < B[j][2] ? (i += 1) : (j += 1)
    end
    return out
end

function Base.issubset(a::CellVector, b::CellVector)
    system(a) == system(b) && a.level == b.level || return isempty(a)
    B = intervals(b.windows)
    j = 1
    for (lo, hi) in intervals(a.windows)
        while j <= length(B) && B[j][2] < lo
            j += 1
        end
        (j <= length(B) && B[j][1] <= lo && hi <= B[j][2]) || return false
    end
    return true
end

# --- show ------------------------------------------------------------------

Base.show(io::IO, cv::CellVector) =
    print(io, "CellVector(", typeof(system(cv)).name.name, ", level=", cv.level,
        ", ncells=", length(cv), ", ", nwindows(cv.windows),
        cv.windows isa RangeWindows ? " windows)" : " positions)")

Base.show(io::IO, ::MIME"text/plain", cv::CellVector) = show(io, cv)
