# Parallel endpoints support binary search in both directions.
struct RangeWindows
    starts::Vector{Int}
    stops::Vector{Int}
    offsets::Vector{Int}
end

struct IndexWindows
    indices::Vector{Int}
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

# Run-table lookup keeps tree descent from repeatedly decoding cell ids.
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

# Maximal runs make logical equality equivalent to endpoint equality.
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

# Choose ranges only when their three arrays use no more integers than indices.
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

# --- interval form ----------------------------------------------------------

intervals(w::RangeWindows) =
    [(@inbounds(w.starts[j]), @inbounds(w.stops[j])) for j in eachindex(w.starts)]

intervals(w::IndexWindows) = [(p, p) for p in w.indices]

# Preserve compact ranges while normalizing sorted intervals.
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

# --- vector -----------------------------------------------------------------

"""
    CellVector(set::MultiOrderCellSet; level = set's reference level)
    CellVector(grid::AbstractGrid)
    CellVector(sys, level, ids::AbstractVector)

`CellVector` is an immutable, lazy vector of strictly ascending cell ids from
one level of one hierarchical system. It stores sorted level-grid index runs
or indices and resolves ids on demand. Memory is O(number of windows), and
indexing is O(log(number of windows)). [`CellLookup`](@ref) provides its
DimensionalData wrapper.

Construct a `CellVector` from:

  - a [`MultiOrderCellSet`](@ref), optionally expanded to a deeper level;
  - an [`AbstractGrid`](@ref), including a [`PartialGrid`](@ref); or
  - an explicit strictly ascending cell-id vector with `sys` and `level`.

`CellVector(cv)` returns `cv`.

```julia
cv[k]
localindex(cv, c)
covering(cv, polygon)
PartialGrid(cv)
```

Systems with [`has_sorted_subtrees`](@ref) construct windows from
[`level_ranges`](@ref) in O(entries). Other systems enumerate, sort, and
compress descendant indices. On A5, non-congruent refinement can make an
expanded [`covering`](@ref) selection over-cover the target; see
[`MultiOrderCoverage`](@ref).
"""
struct CellVector{ID,W<:CellWindows,G<:AbstractGrid,B} <: AbstractCellVector{ID}
    windows::W
    grid::G
    backing::B
    level::Int
end

@inline subset_window_bounds(cv::CellVector, lo::Int, hi::Int) =
    subset_window_bounds(cv.windows, lo, hi)

# Direct run lookup avoids nested searches and id decoding during cursor descent.
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

# Sorting restores canonical order for non-contiguous descendant lists.
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

# Whole levels use one run; subsets resolve their complete-grid indices.
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

# Reuse the original windows when round-tripping through `PartialGrid`.
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

# Derived subsets retain their level grid and drop source provenance.
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

# Immutability makes the whole-vector slice safe to return directly.
Base.getindex(cv::CellVector, ::Colon) = cv

# Vector-only specialization lets Base preserve higher-dimensional index shapes.
Base.getindex(cv::CellVector, idx::AbstractVector{<:Integer}) = _subset(cv, idx)

# Resolve the method ambiguity with SmallCollections vector indexing.
Base.getindex(cv::CellVector,
    i::SmallCollections.AbstractFixedOrSmallOrPackedVector{<:Integer}) = _subset(cv, i)

# Validate boolean axes before `findall`, which discards the mask shape.
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

# Subsetting can only reduce the complete grid's neighbor count.
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

Return the collection that constructed `cv`: a
[`MultiOrderCellSet`](@ref), expanded [`MultiOrderVector`](@ref), or grid.
Derived collections return their describing [`PartialGrid`](@ref).
For [`CellLookup`](@ref), `Base.parent` returns the logical [`CellVector`](@ref).
"""
cellset(cv::CellVector) = _origin(cv, cv.backing)

_origin(::CellVector, backing) = backing
_origin(cv::CellVector, ::Nothing) = PartialGrid(cv)

"""
    localindex(cv::CellVector, c::AbstractCellIndex) -> Union{Int,Nothing}
    localindex(cv::CellVector, p::GO.UnitSphericalPoint) -> Union{Int,Nothing}
    localindex(cv::CellVector, lon::Real, lat::Real) -> Union{Int,Nothing}

Return the local index of cell `c`, or `nothing` when `cv` does not contain it.
This is the inverse of `cv[k]`. The point forms first call [`cellat`](@ref) and
return the matching data index. Longitude and latitude are in degrees.
"""
# Global indices come from the shared complete level.
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

Return the cell in `cv` containing the point, or `nothing` when the point lies
outside the subset. [`localindex`](@ref) returns the corresponding index.
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

Return an O(1) grid view whose index `k` is `cv[k]`. A vector built from a
rooted `PartialGrid` preserves that grid's root and bucket size. Other vectors
return an unrooted grid.
"""
PartialGrid(cv::CellVector) = _partial_grid(cv, cv.backing)

# Preserve rooted metadata while keeping window-based membership for other grids.
_partial_grid(cv::CellVector, pg::PartialGrid) =
    _is_rooted(pg) ? pg : PartialGrid(system(cv), cv.level, _bare(cv, pg))
_partial_grid(cv::CellVector, backing) =
    PartialGrid(system(cv), cv.level, _bare(cv, backing))

# Delegate membership to the compressed windows without decoding probe ids.
localindex(grid::PartialGrid{<:AbstractHierarchicalGridSystem,<:CellVector},
    c::AbstractCellIndex) = localindex(grid.ids, c)

# Preserve an already bare vector; otherwise drop its backing.
_bare(cv::CellVector, backing) = CellVector(cv.windows, cv.grid, nothing, cv.level)
_bare(cv::CellVector, ::Nothing) = cv

# --- selecting a region ----------------------------------------------------

"""
    covering(cv::CellVector, target) -> CellVector

Return the cells of `cv` named by a [`MultiOrderCoverage`](@ref) of `target` at
`cv`'s level. The result remains a windowed `CellVector`.

`target` is anything [`query`](@ref) accepts: a GeoInterface geometry, an
`Extents.Extent` in lon/lat degrees, or a `GO.UnitSpherical.SphericalCap`.

```julia
cv    = CellVector(query(sys, MultiOrderCoverage(canton); level = 9))
basin = covering(cv, watershed)          # a CellVector again
data[covering_indices(cv, watershed)]    # the same selection as indices
```

[`covering_indices`](@ref) returns the same selection in index space and backs
the DimensionalData [`Covering`](@ref) selector.

Selection visits each leaf named by the coverage, even though the result is
stored compactly. Selection at the target level limits expansion work.

Non-congruent refinement can make the result over-cover `target`; see
[`MultiOrderCoverage`](@ref).
"""
covering(cv::CellVector, target) =
    _derive(cv, _windows(_covering_leafindices(cv, target)))

"""
    covering_indices(cv::CellVector, target) -> Vector{Int}

Return the ascending indices selected by [`covering`](@ref), suitable for a
parallel data vector. The [`Covering`](@ref) selector resolves to this form.
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

# --- set arithmetic ---------------------------------------------------------

# Base's first-appearance union order may violate the ascending invariant.

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
