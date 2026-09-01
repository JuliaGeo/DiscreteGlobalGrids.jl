# `CellVector` is a strictly ascending collection of cells from one system and
# level. It stores leaf-grid index windows and resolves ids on demand.
# DimensionalData's `CellLookup` delegates to this dependency-free collection.
#
# Both window representations support these mappings:
#
#   `leafindex(w, k)`    concatenation index -> leaf grid index
#   `windowindex(w, p)`  leaf grid index -> concatenation index or `nothing`
# ---------------------------------------------------------------------------

# Parallel vectors allow both searches to use `searchsortedfirst` without
# deriving range endpoints for each lookup.
struct RangeWindows
    starts::Vector{Int}
    stops::Vector{Int}
    offsets::Vector{Int}     # offsets[j] = cells in windows 1:j
end

struct IndexWindows
    indices::Vector{Int}   # sorted, strictly ascending
end

const CellWindows = Union{RangeWindows,IndexWindows}

Base.length(w::RangeWindows) = isempty(w.offsets) ? 0 : @inbounds w.offsets[end]
Base.length(w::IndexWindows) = length(w.indices)

nwindows(w::RangeWindows) = length(w.starts)
nwindows(w::IndexWindows) = length(w.indices)

@inline function leafindex(w::RangeWindows, k::Int)
    j = searchsortedfirst(w.offsets, k)
    base = j == 1 ? 0 : @inbounds w.offsets[j-1]
    return @inbounds w.starts[j] + (k - base) - 1
end

@inline leafindex(w::IndexWindows, k::Int) = @inbounds w.indices[k]

@inline function windowindex(w::RangeWindows, p::Int)
    j = searchsortedfirst(w.stops, p)
    j <= length(w.stops) || return nothing
    @inbounds start = w.starts[j]
    p >= start || return nothing
    base = j == 1 ? 0 : @inbounds w.offsets[j-1]
    return base + (p - start) + 1
end

@inline function windowindex(w::IndexWindows, p::Int)
    j = searchsortedfirst(w.indices, p)
    (j <= length(w.indices) && @inbounds(w.indices[j]) == p) || return nothing
    return j
end

leafindices(w::CellWindows) = (leafindex(w, k) for k in 1:length(w))

# Map a complete-level index interval to the contiguous logical window that a
# `CellVector` stores from it.  Tree descent asks this question for every child.
# Answering it from the run table avoids binary-searching decoded cell IDs where
# every comparison would itself decode a `CellVector` element.
@inline function subset_window_bounds(w::RangeWindows, lo::Int, hi::Int)
    lo <= hi || return (1, 0)
    first_run = searchsortedfirst(w.stops, lo)
    first_run <= length(w.stops) || return (length(w) + 1, length(w))
    @inbounds first_global = max(lo, w.starts[first_run])
    if first_global > hi
        @inbounds base = first_run == 1 ? 0 : w.offsets[first_run - 1]
        return (base + 1, base)
    end

    last_run = searchsortedlast(w.starts, hi)
    last_run >= first_run || return (1, 0)
    @inbounds last_global = min(hi, w.stops[last_run])

    @inbounds first_base = first_run == 1 ? 0 : w.offsets[first_run - 1]
    @inbounds last_base = last_run == 1 ? 0 : w.offsets[last_run - 1]
    @inbounds first_local = first_base + first_global - w.starts[first_run] + 1
    @inbounds last_local = last_base + last_global - w.starts[last_run] + 1
    return (first_local, last_local)
end

@inline function subset_window_bounds(w::IndexWindows, lo::Int, hi::Int)
    lo <= hi || return (1, 0)
    return (searchsortedfirst(w.indices, lo), searchsortedlast(w.indices, hi))
end

# Classify how much of the leaf block `lo:hi` is stored. Both representations
# inspect the first entry that reaches `lo`:
#
#   * nothing reaches `lo`, or the first thing that does starts past `hi` — the
#     block is empty of stored cells;
#   * that one run covers the whole block — the block is stored entire, and it
#     takes only one run to say so because runs are maximal and disjoint;
#   * anything else — the block is partly stored, which is all the walk needs.
#
# For `IndexWindows`, strict ascent makes the block complete exactly when
# the entry `hi - lo` slots later is `hi`.
@inline function span_windows(w::RangeWindows, lo::Int, hi::Int)
    j = searchsortedfirst(w.stops, lo)
    j <= length(w.stops) || return _SPAN_NONE
    @inbounds start = w.starts[j]
    start > hi && return _SPAN_NONE
    (start <= lo && @inbounds(w.stops[j]) >= hi) && return _SPAN_ALL
    return _SPAN_SOME
end

@inline function span_windows(w::IndexWindows, lo::Int, hi::Int)
    j = searchsortedfirst(w.indices, lo)
    j <= length(w.indices) || return _SPAN_NONE
    @inbounds p = w.indices[j]
    p > hi && return _SPAN_NONE
    if p == lo
        k = j + (hi - lo)
        (k <= length(w.indices) && @inbounds(w.indices[k]) == hi) &&
            return _SPAN_ALL
    end
    return _SPAN_SOME
end

# Runs are stored maximally, so equal `RangeWindows` have identical boundaries.
# Set operations normalize their output to preserve this invariant.
Base.:(==)(a::RangeWindows, b::RangeWindows) =
    a.starts == b.starts && a.stops == b.stops

# Compare logical indices when the window representations differ.
Base.:(==)(a::CellWindows, b::CellWindows) =
    length(a) == length(b) && all(((x, y),) -> x == y, zip(leafindices(a), leafindices(b)))

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

# Run-compress sorted indices and keep the smaller representation. A range
# uses three integers; an explicit index uses one.
function _windows(indices::Vector{Int})
    isempty(indices) && return _empty_windows()
    runs = 1
    for i in 2:length(indices)
        @inbounds indices[i] > indices[i-1] || throw(ArgumentError(
            "leaf indices must be strictly ascending, got $(indices[i-1]) " *
            "then $(indices[i]) at index $i"))
        @inbounds indices[i] == indices[i-1] + 1 || (runs += 1)
    end
    3 * runs <= length(indices) || return IndexWindows(indices)
    starts = Vector{Int}(undef, runs)
    stops = Vector{Int}(undef, runs)
    offsets = Vector{Int}(undef, runs)
    j = 1
    @inbounds starts[1] = indices[1]
    for i in 2:length(indices)
        @inbounds if indices[i] != indices[i-1] + 1
            stops[j] = indices[i-1]
            offsets[j] = i - 1
            j += 1
            starts[j] = indices[i]
        end
    end
    @inbounds stops[runs] = indices[end]
    @inbounds offsets[runs] = length(indices)
    return RangeWindows(starts, stops, offsets)
end

# --- the two shapes, read as intervals -------------------------------------
#
# Set operations convert both representations to intervals. An explicit
# index becomes a one-cell interval.

intervals(w::RangeWindows) =
    [(@inbounds(w.starts[j]), @inbounds(w.stops[j])) for j in eachindex(w.starts)]

intervals(w::IndexWindows) = [(p, p) for p in w.indices]

# Convert sorted intervals to maximal windows without expanding ranges merely
# to choose a representation.
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
    indices = Vector{Int}(undef, total)
    k = 0
    for (lo, hi) in merged, p in lo:hi
        indices[k+=1] = p
    end
    return IndexWindows(indices)
end

# ===========================================================================
# The vector
# ===========================================================================

"""
    CellVector(set::MultiOrderCellSet; level = set's reference level)
    CellVector(grid::AbstractGrid)
    CellVector(sys, level, ids::AbstractVector)

An immutable, lazy `AbstractVector` of strictly ascending cell ids from one
level of one hierarchical system. It stores their leaf-grid index windows
instead of the ids.

Semantically `cv` **is** the id vector: `length(cv)` is the number of cells,
`cv[k]` is the `k`th of them, `collect(cv)` is the vector itself. What is
*stored* is sorted, disjoint intervals or a sorted index list in
`levelgrid(system(cv), level(cv))`, so memory is O(number of windows) rather
than O(number of cells). `cv[k]` searches the windows' cumulative lengths and
resolves one `cellindex`; [`localindex`](@ref) runs the inverse. Nothing is
materialised.

[`CellLookup`](@ref) provides the DimensionalData wrapper and delegates its
collection operations to this type.

# The three ways in

  - a [`MultiOrderCellSet`](@ref), the compressed coverage, optionally
    re-expanded to a deeper `level` than the set's own reference level;
  - an [`AbstractGrid`](@ref) — `levelgrid(sys, l)`, a whole level, which is one
    window; or a [`PartialGrid`](@ref), which is that subset's indices (one
    window when the subset is a subtree, the explicit list when it is
    scattered);
  - `sys, level, ids` — an explicit **strictly ascending** id vector, validated
    and run-compressed.

`CellVector(cv)` returns `cv` unchanged.

# The verbs

```julia
cv[k]                          # the kth cell id
localindex(cv, c)              # its inverse: Int, or `nothing`
c in cv                        # membership, O(log #windows)
localindex(cv, lon, lat)       # the index of the cell a point falls in
cellat(cv, lon, lat)           # that cell's id instead
covering(cv, polygon)          # the sub-vector a region's coverage names
intersect(cv, other)           # two vectors at the same level, O(#windows)
PartialGrid(cv)                # read as a grid, O(1) — the regridding handshake
cellset(cv)                    # what it was built from
```

`Base.summarysize(cv)` is O(number of windows) plus the object referenced by
[`cellset`](@ref). Indexing is O(log(number of windows)) and allocation-free.

Where [`has_sorted_subtrees`](@ref) holds, [`level_ranges`](@ref) constructs the
windows in
O(#entries). Where it does not (A5), a cell's descendants are not one interval
of their level, so the vector is built by SELECTION: `descendants` names the
leaves, they are resolved to indices and sorted, and the result is
run-compressed like any other index list. Every method above is unchanged
as any other index list. This construction visits every leaf. The stored
form ranges from a small set of windows to one index per cell.

For A5, descendants need not lie inside their parent's footprint, so expanding a coverage
names leaves the target does not touch — most visibly inside a hole. A
[`covering`](@ref) subset on A5 is therefore a superset of the cells that meet
the region, by the same margin the refinement itself is; see
[`MultiOrderCoverage`](@ref).
"""
struct CellVector{ID,W<:CellWindows,G<:AbstractGrid,B} <: AbstractCellVector{ID}
    windows::W
    grid::G                  # `levelgrid(system, level)` — every `cv[k]` reads it
    backing::B               # what it was built from, or `nothing` when derived
    level::Int
end

@inline subset_window_bounds(cv::CellVector, lo::Int, hi::Int) =
    subset_window_bounds(cv.windows, lo, hi)

# A compressed `CellVector` already stores complete-grid index runs.  Consult
# those runs directly during cursor descent: generic `searchsorted*` over the
# logical vector makes every comparison pay a second binary search plus cell-ID
# decoding.
@inline function _partial_child_window(
        ids::CellVector, complete, range, first_index::Int, last_index::Int)
    lo, hi = subset_window_bounds(ids, Int(first(range)), Int(last(range)))
    lo = max(lo, first_index)
    hi = min(hi, last_index)
    return hi < lo ? (first_index, first_index - 1) : (lo, hi)
end

function CellVector(windows::CellWindows, grid::AbstractGrid, backing, l::Integer)
    ID = cellindextype(system(grid))
    return CellVector{ID,typeof(windows),typeof(grid),typeof(backing)}(
        windows, grid, backing, Int(l))
end

# Convert the shadowing `level` keyword to a positional `Int` internally.
CellVector(set::MultiOrderCellSet; level::Integer=set.reference_level) =
    _cellvector(set, Int(level))

function _cellvector(set::MultiOrderCellSet, l::Int)
    sys = system(set)
    grid = levelgrid(sys, l)
    w = has_sorted_subtrees(sys) ?
        _range_windows(level_ranges(set, l)) :
        _windows(_selection_indices(set, grid, l))
    return CellVector(w, grid, set, l)
end

# Selection mode sorts leaf indices because sibling descendant lists need not
# concatenate in canonical order.
function _selection_indices(set::MultiOrderCellSet, grid::AbstractGrid, l::Int)
    sys = system(set)
    out = Int[]
    for c in set
        level(c) <= l || throw(ArgumentError(
            "cannot expand to level $l: the set contains a level-$(level(c)) cell"))
        for d in descendants(sys, c, l)
            p = globalindex(grid, d)
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

# A complete level occupies indices `1:ncells`; proper subsets are resolved
# cell by cell.
function _grid_windows(grid::AbstractGrid, complete::AbstractGrid)
    n = ncells(grid)
    n == ncells(complete) && return _range_windows((1:n,))
    return _windows(_grid_indices(grid, complete))
end

# A rooted, sorted subtree maps directly to one index window.
function _grid_windows(grid::PartialGrid{<:Any,<:SubtreeIds}, complete::AbstractGrid)
    ids = grid.ids
    ids.n == 0 && return _empty_windows()
    return _range_windows((ids.first:(ids.first+ids.n-1),))
end

# Reuse the windows when converting a `CellVector`-backed `PartialGrid` back to
# a vector.
_grid_windows(grid::PartialGrid{<:Any,<:CellVector}, complete::AbstractGrid) =
    grid.ids.windows

function _grid_indices(grid::AbstractGrid, complete::AbstractGrid)
    out = Vector{Int}(undef, ncells(grid))
    for i in eachindex(out)
        c = cellindex(grid, i)
        p = globalindex(complete, c)
        p === nothing && throw(ArgumentError(
            "$c is not a cell of levelgrid($(system(grid)), $(level(grid)))"))
        out[i] = p
    end
    return out
end

# Derived subsets keep the leaf grid and level but replace provenance with a
# lazily constructed `PartialGrid`.
_derive(cv::CellVector, w::CellWindows) = CellVector(w, cv.grid, nothing, cv.level)

# A `CellVector` IS the region container, so `region` is the identity on it.
DGG.region(cv::CellVector) = cv

windows(cv::CellVector) = cv.windows

# --- the collection surface ------------------------------------------------

Base.size(cv::CellVector) = (length(cv.windows),)
Base.IndexStyle(::Type{<:CellVector}) = Base.IndexLinear()

Base.@propagate_inbounds function Base.getindex(cv::CellVector, k::Int)
    @boundscheck checkbounds(cv, k)
    return cellindex(cv.grid, leafindex(cv.windows, k))
end

# Immutable, so the whole-vector slice is the vector rather than a copy of it.
Base.getindex(cv::CellVector, ::Colon) = cv

# Ascending unique indices remain windowed; permutations, repetitions, and
# reversals return an ordinary materialized vector.
#
# `AbstractVector`, not `AbstractArray`: indexing an array by a
# higher-dimensional index returns something of the INDEX's shape, and a window
# set has no shape to give back. Narrowing here lets `cv[matrix]` fall through
# to Base's generic, which answers with a matrix of ids — so the answer's shape
# no longer depends on whether the index happened to be ascending.
Base.getindex(cv::CellVector, idx::AbstractVector{<:Integer}) = _subset(cv, idx)

# Resolve the method ambiguity with SmallCollections vector indexing.
Base.getindex(cv::CellVector,
    i::SmallCollections.AbstractFixedOrSmallOrPackedVector{<:Integer}) = _subset(cv, i)

# Handle boolean masks inside `_subset` because `Bool <: Integer`. Validate axes
# before `findall` so masks with the wrong length or rank throw `BoundsError`.
function _subset(cv::CellVector, mask::AbstractArray{Bool})
    axes(mask) == axes(cv) || throw(BoundsError(cv, (mask,)))
    return _subset(cv, findall(mask))
end

function _subset(cv::CellVector, idx::AbstractArray{<:Integer})
    n = length(cv)
    indices = Vector{Int}(undef, length(idx))
    ascending = true
    for (j, k) in enumerate(idx)
        1 <= k <= n || throw(BoundsError(cv, (idx,)))
        indices[j] = leafindex(cv.windows, Int(k))
        j > 1 && indices[j] <= indices[j-1] && (ascending = false)
    end
    ascending || return [cv[Int(k)] for k in idx]
    return _derive(cv, _windows(indices))
end

# Ascending leaf-grid indices imply ascending canonical ids.
Helpers.strictly_increasing(::CellVector) = true

# Equal systems, levels, and logical windows imply elementwise equality.
Base.:(==)(a::CellVector, b::CellVector) =
    system(a) == system(b) && a.level == b.level && a.windows == b.windows

# --- what the vector is, in this package's own vocabulary ------------------

"""
    system(cv::CellVector)

The grid system the vector's cells are named in.
"""
system(cv::CellVector) = system(cv.grid)

# A subset can only lose neighbours from the complete system's one-ring, so
# the system-wide bound remains valid for every compressed collection shape.
maxneighbors(cv::CellVector, connectivity::Connectivity) =
    maxneighbors(cv.grid, connectivity)
maxneighbors(cv::CellVector) = maxneighbors(cv, Vertex())

"""
    level(cv::CellVector) -> Int

The one level every cell in the vector sits at.
"""
level(cv::CellVector) = cv.level

"""
    cellset(cv::CellVector)
    cellset(lk::CellLookup)

Return the [`MultiOrderCellSet`](@ref) or grid used to build the collection.

A collection *derived* from another one, by indexing or by [`covering`](@ref),
has no such origin and reports the [`PartialGrid`](@ref) describing it instead.

For [`CellLookup`](@ref), `Base.parent` returns the logical values as a
[`CellVector`](@ref).
"""
cellset(cv::CellVector) = _origin(cv, cv.backing)

# Dispatch on the backing value to distinguish `Nothing` without ambiguous
# partially specified `CellVector` parameter bounds.
_origin(::CellVector, backing) = backing
_origin(cv::CellVector, ::Nothing) = PartialGrid(cv)

"""
    localindex(cv::CellVector, c::AbstractCellIndex) -> Union{Int,Nothing}
    localindex(cv::CellVector, p::GO.UnitSphericalPoint) -> Union{Int,Nothing}
    localindex(cv::CellVector, lon::Real, lat::Real) -> Union{Int,Nothing}

Local index of a cell in the vector, or `nothing` when the vector does not hold
it — including when `c` is at another level. The inverse of `cv[k]`, and the
half of the bijection every selection ends at.

The point forms name the cell through [`cellat`](@ref) first and then answer
with its local index, which is the one-call form of "where in my data does
this point land". Degrees for the `(lon, lat)` method, as everywhere else.
"""
# A subset's storage is carved out of a complete level, so its global index is
# that level's — the number two different subsets of one level agree on.
globalindex(cv::CellVector, c::AbstractCellIndex) = globalindex(cv.grid, c)

function localindex(cv::CellVector, c::AbstractCellIndex)
    p = globalindex(cv.grid, c)
    p === nothing && return nothing
    return windowindex(cv.windows, p)
end

function localindex(cv::CellVector, p::GO.UnitSphericalPoint)
    c = cellat(cv.grid, p)
    c === nothing && return nothing
    return localindex(cv, c)
end

localindex(cv::CellVector, lon::Real, lat::Real) =
    localindex(cv, unit_point(lon, lat))

"""
    cellat(cv::CellVector, p::GO.UnitSphericalPoint) -> Union{AbstractCellIndex,Nothing}
    cellat(cv::CellVector, lon::Real, lat::Real)

The cell **of `cv`** containing the point, or `nothing` when the point falls
outside the cells the vector holds. Same contract as [`cellat`](@ref) on a
grid, restricted to the subset: a point inside the level grid but outside this
vector answers `nothing` rather than naming a cell that is not here.

[`localindex`](@ref) is the same question answered as an index.
"""
function cellat(cv::CellVector, p::GO.UnitSphericalPoint)
    c = cellat(cv.grid, p)
    c === nothing && return nothing
    return localindex(cv, c) === nothing ? nothing : c
end

cellat(cv::CellVector, lon::Real, lat::Real) = cellat(cv, unit_point(lon, lat))

# Membership uses the O(log(number of windows)) inverse index search.
Base.in(c::AbstractCellIndex, cv::CellVector) = localindex(cv, c) !== nothing

# Classify a complete-grid index block directly from the windows.
subset_span(cv::CellVector, lo::Int, hi::Int) = span_windows(cv.windows, lo, hi)

"""
    PartialGrid(cv::CellVector) -> PartialGrid

Return a grid whose index `k` is `cv[k]`. Construction is O(1) and keeps the
ids lazy, so data indexed by the vector needs no permutation.

A vector built from a rooted `PartialGrid` returns that grid, preserving its
root and bucket size. Other vectors return an unrooted grid because their
windows do not identify an ancestor.
"""
PartialGrid(cv::CellVector) = _partial_grid(cv, cv.backing)

# Preserve a rooted backing grid. Rebuild unrooted grids around the compressed
# vector so membership continues to use its windows.
_partial_grid(cv::CellVector, pg::PartialGrid) =
    _is_rooted(pg) ? pg : PartialGrid(system(cv), cv.level, _bare(cv, pg))
_partial_grid(cv::CellVector, backing) =
    PartialGrid(system(cv), cv.level, _bare(cv, backing))

# Delegate membership to the compressed vector to avoid decoding one id per
# binary-search probe. Spell the system parameter with its declared bound to
# keep this method more specific than the general `PartialGrid` method.
localindex(grid::PartialGrid{<:AbstractHierarchicalGridSystem,<:CellVector},
    c::AbstractCellIndex) = localindex(grid.ids, c)

# Preserve an already bare vector; otherwise drop its backing.
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
data[covering_indices(cv, watershed)]    # the same selection as indices
```

[`covering_indices`](@ref) is the index-space form, for indexing a data
array laid out against `cv`. This is what the `DimensionalData` selector
[`Covering`](@ref) is spelled as outside `DimensionalData`.

Selection visits each leaf named by the coverage, even though the result is
stored compactly. Select at the level being read to avoid unnecessary expansion.

Coverage means *covering*: the result is a superset of the cells that meet
`target`, by whatever margin the system's refinement is non-congruent. See
[`MultiOrderCoverage`](@ref) for the size of that margin per system.
"""
covering(cv::CellVector, target) =
    _derive(cv, _windows(_covering_leafindices(cv, target)))

"""
    covering_indices(cv::CellVector, target) -> Vector{Int}

The indices in `cv` of the cells [`covering`](@ref) selects, ascending — for
indexing a data array laid out against `cv` without building the sub-vector.

`covering(cv, target)` and `cv[covering_indices(cv, target)]` name the same
cells; this form is the one a cube's `getindex` needs, and is what the
[`Covering`](@ref) selector resolves to.
"""
function covering_indices(cv::CellVector, target)
    out = Int[]
    w = cv.windows
    _each_leaf_index(cv, target) do p
        k = windowindex(w, p)
        k === nothing || push!(out, k)
    end
    return issorted(out) ? out : sort!(out)
end

"""
    predicate_indices(cv::CellVector, pred::DE9IMPredicate) -> Vector{Int}

The indices in `cv` of the cells that satisfy `pred` at `cv`'s level, ascending
— `query(system(cv), pred; level = level(cv))` intersected with `cv`, answered
in index space so a data array laid out against `cv` can be indexed by it.

`pred` is any predicate [`query`](@ref) implements, over any target it accepts;
`cv[predicate_indices(cv, Within(cap))]` names the cells of `cv` lying wholly
inside `cap`. This is what a predicate used as a [`Cells`](@ref) selector
resolves to. Unlike [`covering_indices`](@ref), the answer is exact: it
inherits no over-covering from a coverage.
"""
function predicate_indices(cv::CellVector, pred::DE9IM.DE9IMPredicate)
    out = Int[]
    for c in query(cv.grid, pred)
        k = localindex(cv, c)
        k === nothing || push!(out, k)
    end
    return issorted(out) ? out : sort!(out)
end

function _covering_leafindices(cv::CellVector, target)
    out = Int[]
    w = cv.windows
    _each_leaf_index(cv, target) do p
        windowindex(w, p) === nothing || push!(out, p)
    end
    return issorted(out) ? out : sort!(out)
end

# Use ranges for sorted subtrees and explicit descendants otherwise.
function _each_leaf_index(f, cv::CellVector, target)
    sys = system(cv)
    set = query(sys, MultiOrderCoverage(target); level=cv.level)
    if has_sorted_subtrees(sys)
        for r in level_ranges(set, cv.level), p in r
            f(p)
        end
    else
        for p in _selection_indices(set, cv.grid, cv.level)
            f(p)
        end
    end
    return nothing
end

# --- set arithmetic over the windows ---------------------------------------
#
# `intersect` preserves the left operand's ascending order; `issubset` performs
# a set comparison. Both run in O(number of windows).
#
# Across different systems or levels, only an empty left vector is a subset.
# `intersect` throws because no result system and level can be selected.
#
# `union` is not specialized because Base preserves first-appearance order,
# which need not be ascending for two ascending operands.

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
        cv.windows isa RangeWindows ? " windows)" : " indices)")

Base.show(io::IO, ::MIME"text/plain", cv::CellVector) = show(io, cv)
