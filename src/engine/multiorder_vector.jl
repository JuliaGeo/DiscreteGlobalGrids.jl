# ---------------------------------------------------------------------------
# `MultiOrderVector`: cells at mixed refinement levels, pairwise disjoint as
# subtrees, sorted by their descendant-range intervals at a reference level.
# The storage form of a MOC (`MultiOrderCellSet` is the query form). Every
# query is a binary search over the flat interval index — the layout astronomy
# MOC libraries use; there is no tree. Everything requires
# `has_sorted_subtrees(sys)`: without descendant ranges a subtree is not one
# interval of a level. A5 throws.
# ---------------------------------------------------------------------------

"""
    MultiOrderVector(set::MultiOrderCellSet)
    MultiOrderVector(sys, cells::AbstractVector; reference_level = deepest cell level)

An immutable `AbstractVector` of cells at mixed refinement levels of one
hierarchical system, pairwise disjoint as subtrees (no member is an ancestor
of another, none repeats), ordered by the position intervals their subtrees
occupy at a **reference level** — depth-first order, the order a sorted
single-level axis has and coarse-ancestor store chunking assumes
([`dggwrite`](@ref)). Semantically it is that id vector; what is
stored beside the ids is the interval index (`starts`, `stops`, cumulative
`offsets`) that answers every query below in O(log n).

Construction from a [`MultiOrderCellSet`](@ref) is O(n) and keeps the set's
reference level; from `sys, cells` the vector is sorted and validated, with
`reference_level` defaulting to the deepest level present (a deeper one is
legal). The sort discards the input order, so a parallel value vector must be
permuted with it — `cellposition(mov, c)` gives each cell's stored position.
The reference level is a property of the container, not of the system: two
containers of the same cells at different reference levels compare equal.

```julia
mov[k]                          # the kth cell id
cellposition(mov, c)            # exact membership: Int, or `nothing`
c in mov                        # the same as a Bool
covering_position(mov, c)       # the stored cell that is `c` or an ancestor of it
cellat(mov, lon, lat)           # the stored cell a point falls in
covering(mov, polygon)          # the sub-container a region names
CellVector(mov; level = l)      # the leaf-level expansion, as windows
union(a, b); intersect(a, b); setdiff(a, b); complement(a)
issubset(a, b); isdisjoint(a, b); issetequal(a, b); symdiff(a, b)
```

Membership and point queries are one binary search; the set operations are
O(n + m) interval arithmetic and return **normalized minimal** containers —
the coarsest cells tiling the result — so `union(a, a)` may be coarser than a
non-minimal `a`. `==` compares the stored cells, so those two are unequal even
though they hold the same leaves; `issetequal` is the leaf-set comparison. The
system must have [`has_sorted_subtrees`](@ref).
"""
struct MultiOrderVector{ID,S<:AbstractHierarchicalGridSystem} <: AbstractVector{ID}
    system::S
    cells::Vector{ID}
    starts::Vector{Int}      # first(descendant_range(sys, c, reference_level))
    stops::Vector{Int}       # ... and its inclusive last
    offsets::Vector{Int}     # offsets[i] = leaves at the reference level in cells 1:i
    reference_level::Int

    # Validation lives in the inner constructor; `checked = false` is for
    # callers that have just established the invariants themselves.
    function MultiOrderVector{ID,S}(system::S, cells::Vector{ID}, starts::Vector{Int},
        stops::Vector{Int}, offsets::Vector{Int}, ref::Int,
        checked::Bool) where {ID,S<:AbstractHierarchicalGridSystem}
        checked && _check_multiorder(system, cells, starts, stops, ref)
        return new{ID,S}(system, cells, starts, stops, offsets, ref)
    end
end

# `offsets` is derived, never passed in, so it cannot disagree with the
# intervals.
function _multiorder_vector(sys::S, cells::Vector{ID}, starts::Vector{Int},
    stops::Vector{Int}, ref::Int, checked::Bool) where {ID,S<:AbstractHierarchicalGridSystem}
    length(cells) == length(starts) == length(stops) || throw(ArgumentError(
        "a multi-order vector needs one interval per cell, got $(length(cells)) " *
        "cells against $(length(starts)) starts and $(length(stops)) stops"))
    offsets = Vector{Int}(undef, length(cells))
    total = 0
    for i in eachindex(offsets)
        total += stops[i] - starts[i] + 1
        offsets[i] = total
    end
    return MultiOrderVector{ID,S}(sys, cells, starts, stops, offsets, ref, checked)
end

# Checked before any `descendant_range` call, which would otherwise fail with
# a less specific error.
function _check_reference(sys::AbstractHierarchicalGridSystem, cells, ref::Int)
    has_sorted_subtrees(sys) || throw(ArgumentError(
        "$(typeof(sys)) has no descendant ranges, so a cell's subtree is not one " *
        "interval of a level and a multi-order vector cannot key its cells"))
    ref in levels(sys) || throw(ArgumentError(
        "reference level $ref is outside $(typeof(sys))'s levels $(levels(sys))"))
    for c in cells
        level(c) <= ref || throw(ArgumentError(
            "reference level $ref is shallower than the level-$(level(c)) cell $c, " *
            "which therefore has no interval at it"))
    end
    return nothing
end

function _check_multiorder(sys::AbstractHierarchicalGridSystem, cells, starts, stops,
    ref::Int)
    _check_reference(sys, cells, ref)
    n = ncells(levelgrid(sys, ref))
    for i in eachindex(starts)
        starts[i] <= stops[i] || throw(ArgumentError(
            "the interval of $(cells[i]) at level $ref is empty: " *
            "$(starts[i]):$(stops[i])"))
        # Catches an id outside the system here rather than at the first
        # geometry call that decodes it.
        1 <= starts[i] && stops[i] <= n || throw(ArgumentError(
            "$(cells[i]) is not a cell of $(typeof(sys)): its level-$ref interval " *
            "$(starts[i]):$(stops[i]) leaves the level's positions 1:$n"))
        i == 1 && continue
        # Overlapping reference-level intervals mean one cell contains the
        # other, or repeats it.
        starts[i] > stops[i-1] || throw(ArgumentError(
            "multi-order cells must be pairwise disjoint: $(cells[i-1]) holds " *
            "positions $(starts[i-1]):$(stops[i-1]) at level $ref and $(cells[i]) " *
            "holds $(starts[i]):$(stops[i]), so one descends from the other or " *
            "the two repeat"))
    end
    return nothing
end

# --- the two ways in -------------------------------------------------------

function MultiOrderVector(set::MultiOrderCellSet)
    sys = system(set)
    ref = set.reference_level
    _check_reference(sys, set.cells, ref)
    # `set.keys` ARE the reference-level interval starts, already sorted, so
    # only the stops are new work.
    stops = [last(descendant_range(sys, c, ref)) for c in set.cells]
    return _multiorder_vector(sys, copy(set.cells), copy(set.keys), stops, ref, true)
end

# The keyword shadows the `reference_level` accessor, so the work takes a
# positional `Int`.
MultiOrderVector(sys::AbstractHierarchicalGridSystem, cells::AbstractVector;
    reference_level::Integer=_default_reference_level(sys, cells)) =
    _multiorder_from_cells(sys, cells, Int(reference_level))

_default_reference_level(sys::AbstractHierarchicalGridSystem, cells) =
    isempty(cells) ? first(levels(sys)) : maximum(level, cells)

function _multiorder_from_cells(sys::AbstractHierarchicalGridSystem, cells, ref::Int)
    # An abstract input eltype would make every downstream cell call dynamic.
    E = eltype(cells)
    ID = isconcretetype(E) && E <: AbstractCellIndex ? E : cellindextype(sys)
    cs = collect(ID, cells)
    _check_reference(sys, cs, ref)
    ranges = [descendant_range(sys, c, ref) for c in cs]
    perm = sortperm(ranges; by=first)
    return _multiorder_vector(sys, cs[perm], [first(ranges[p]) for p in perm],
        [last(ranges[p]) for p in perm], ref, true)
end

MultiOrderVector(mov::MultiOrderVector) = mov

# A subset of a valid container keeps its system, its reference level and its
# disjointness; only the intervals it keeps are new.
_derive(mov::MultiOrderVector, cells::Vector, starts::Vector{Int}, stops::Vector{Int}) =
    _multiorder_vector(mov.system, cells, starts, stops, mov.reference_level, false)

# Re-keying to a DEEPER reference level is a change of units: the intervals
# scale, their order and their disjointness do not.
function _rekey(mov::MultiOrderVector, ref::Int)
    ref == mov.reference_level && return mov
    ref >= mov.reference_level || throw(ArgumentError(
        "cannot re-key a multi-order vector from reference level " *
        "$(mov.reference_level) to the shallower level $ref: a stored cell may " *
        "be deeper than $ref"))
    ranges = [descendant_range(mov.system, c, ref) for c in mov.cells]
    return _multiorder_vector(mov.system, mov.cells, first.(ranges), last.(ranges), ref, false)
end

# --- the collection surface ------------------------------------------------

Base.size(mov::MultiOrderVector) = (length(mov.cells),)
Base.IndexStyle(::Type{<:MultiOrderVector}) = Base.IndexLinear()

Base.@propagate_inbounds function Base.getindex(mov::MultiOrderVector, k::Int)
    @boundscheck checkbounds(mov, k)
    return @inbounds mov.cells[k]
end

# Immutable, so the whole-vector slice is the vector rather than a copy of it.
Base.getindex(mov::MultiOrderVector, ::Colon) = mov

# Ascending indices keep the interval index; any other index is answered with
# a plain `Vector`. `AbstractVector` rather than `AbstractArray`, as in
# `CellVector`, so a shaped index falls through to Base's generic.
Base.getindex(mov::MultiOrderVector, idx::AbstractVector{<:Integer}) = _subset(mov, idx)

# Ambiguity tie-break against SmallCollections' own `getindex`, as in
# `cell_vector.jl`.
Base.getindex(mov::MultiOrderVector,
    i::SmallCollections.AbstractFixedOrSmallOrPackedVector{<:Integer}) = _subset(mov, i)

# A mask of the wrong shape is a bounds error; `findall` alone cannot check it.
function _subset(mov::MultiOrderVector, mask::AbstractArray{Bool})
    axes(mask) == axes(mov) || throw(BoundsError(mov, (mask,)))
    return _subset(mov, findall(mask))
end

function _subset(mov::MultiOrderVector, idx::AbstractArray{<:Integer})
    n = length(mov)
    ascending = true
    prev = 0
    for k in idx
        1 <= k <= n || throw(BoundsError(mov, (idx,)))
        k <= prev && (ascending = false)
        prev = Int(k)
    end
    ascending || return [mov.cells[Int(k)] for k in idx]
    ks = Int[Int(k) for k in idx]
    return _derive(mov, mov.cells[ks], mov.starts[ks], mov.stops[ks])
end

"""
    system(mov::MultiOrderVector)

The grid system the container's cells are named in.
"""
system(mov::MultiOrderVector) = mov.system

"""
    reference_level(mov::MultiOrderVector) -> Int

The level the container's intervals are stated at, no shallower than any
stored cell. A property of the container, not of the system, and not part of
equality. Public but unexported: `DiscreteGlobalGrids.reference_level(mov)`.
"""
reference_level(mov::MultiOrderVector) = mov.reference_level

# The container's own intervals, in the shape the set operations merge.
intervals(mov::MultiOrderVector) =
    [(@inbounds(mov.starts[i]), @inbounds(mov.stops[i])) for i in eachindex(mov.starts)]

# --- membership, and the covering ancestor ---------------------------------

"""
    cellposition(mov::MultiOrderVector, c::AbstractCellIndex) -> Union{Int,Nothing}
    cellposition(mov::MultiOrderVector, p::GO.UnitSphericalPoint) -> Union{Int,Nothing}
    cellposition(mov::MultiOrderVector, lon::Real, lat::Real) -> Union{Int,Nothing}

Position of `c` in the container, or `nothing` when it does not hold that
exact cell — a stored ancestor or descendant of `c` does not count. One binary
search, O(log n). The point forms locate the point's covering stored cell, as
[`cellat`](@ref) does; lon/lat in degrees.

[`covering_position`](@ref) answers which stored cell speaks for `c`.
"""
function cellposition(mov::MultiOrderVector, c::AbstractCellIndex)
    lc = level(c)
    lc <= mov.reference_level || return nothing
    r = descendant_range(mov.system, c, mov.reference_level)
    j = searchsortedfirst(mov.starts, first(r))
    j <= length(mov.starts) || return nothing
    # Same interval and same level identify the cell uniquely.
    (@inbounds(mov.starts[j]) == first(r) && @inbounds(mov.stops[j]) == last(r) &&
     level(@inbounds mov.cells[j]) == lc) || return nothing
    return j
end

"""
    covering_position(mov::MultiOrderVector, c::AbstractCellIndex) -> Union{Int,Nothing}

Position of the stored cell that is `c` or an ancestor of `c`, or `nothing`
when none covers it — resolves a leaf to the cell holding its value. `c` may
be deeper than the container's [`reference_level`](@ref); it is keyed through
its reference-level ancestor. O(log n). [`cellposition`](@ref) is exact
membership.
"""
function covering_position(mov::MultiOrderVector, c::AbstractCellIndex)
    ref = mov.reference_level
    keyed = level(c) > ref ? ancestor(mov.system, c, ref) : c
    r = descendant_range(mov.system, keyed, ref)
    j = searchsortedlast(mov.starts, first(r))
    j >= 1 || return nothing
    return last(r) <= @inbounds(mov.stops[j]) ? j : nothing
end

# The interval search is O(log n); Base's generic scan would be O(n).
Base.in(c::AbstractCellIndex, mov::MultiOrderVector) = cellposition(mov, c) !== nothing

"""
    cellat(mov::MultiOrderVector, p::GO.UnitSphericalPoint) -> Union{AbstractCellIndex,Nothing}
    cellat(mov::MultiOrderVector, lon::Real, lat::Real)

The stored cell containing the point, at whatever level it sits, or `nothing`
when none covers it: one point location at the reference level plus one
binary search. [`cellposition`](@ref)'s point forms answer with the position
instead.
"""
function cellat(mov::MultiOrderVector, p::GO.UnitSphericalPoint)
    k = cellposition(mov, p)
    return k === nothing ? nothing : @inbounds mov.cells[k]
end

cellat(mov::MultiOrderVector, lon::Real, lat::Real) = cellat(mov, unit_point(lon, lat))

function cellposition(mov::MultiOrderVector, p::GO.UnitSphericalPoint)
    c = cellat(levelgrid(mov.system, mov.reference_level), p)
    c === nothing && return nothing
    return covering_position(mov, c)
end

cellposition(mov::MultiOrderVector, lon::Real, lat::Real) =
    cellposition(mov, unit_point(lon, lat))

# --- the bridge into everything `CellVector` knows -------------------------

"""
    CellVector(mov::MultiOrderVector; level = reference_level(mov))

The container expanded to single-level ids: the leaf cells its subtrees name
at `level`, as a windowed [`CellVector`](@ref). `level` must be no shallower
than the deepest stored cell. The windows are the stored cells'
`descendant_range`s merged where adjacent, agreeing with
[`level_ranges`](@ref) window for window. O(n) in the number of stored cells,
whatever `level` is. [`cellset`](@ref) on the result returns `mov`.
"""
CellVector(mov::MultiOrderVector; level::Integer=reference_level(mov)) =
    _cellvector(mov, Int(level))

function _cellvector(mov::MultiOrderVector, l::Int)
    sys = mov.system
    for c in mov.cells
        level(c) <= l || throw(ArgumentError(
            "cannot expand to level $l: the container holds a level-$(level(c)) cell"))
    end
    return CellVector(_range_windows(_expanded_ranges(mov, l)), levelgrid(sys, l), mov, l)
end

# Adjacent subtrees merge into one window, matching `level_ranges`.
function _expanded_ranges(mov::MultiOrderVector, l::Int)
    out = UnitRange{Int}[]
    for c in mov.cells
        r = descendant_range(mov.system, c, l)
        if !isempty(out) && first(r) == last(out[end]) + 1
            out[end] = first(out[end]):last(r)
        else
            push!(out, r)
        end
    end
    return out
end

# --- selecting a region ----------------------------------------------------

"""
    covering(mov::MultiOrderVector, target) -> MultiOrderVector

The cells of `mov` whose subtrees meet a [`MultiOrderCoverage`](@ref) of
`target`, kept whole, never clipped. `target` is anything [`query`](@ref)
accepts: a GeoInterface geometry, an `Extents.Extent` in lon/lat degrees, or
a `GO.UnitSpherical.SphericalCap`. O(log n) per coverage interval.
[`covering_positions`](@ref) is the position-space form.
"""
function covering(mov::MultiOrderVector, target)
    ks = covering_positions(mov, target)
    return _derive(mov, mov.cells[ks], mov.starts[ks], mov.stops[ks])
end

"""
    covering_positions(mov::MultiOrderVector, target) -> Vector{Int}

The positions in `mov` of the cells [`covering`](@ref) selects, ascending — for
indexing a data array laid out against the container without building the
sub-container.
"""
function covering_positions(mov::MultiOrderVector, target)
    out = Int[]
    isempty(mov) && return out
    ref = mov.reference_level
    set = query(mov.system, MultiOrderCoverage(target); level=ref)
    seen = 0
    for r in level_ranges(set, ref)
        lo, hi = first(r), last(r)
        j = searchsortedfirst(mov.stops, lo)
        while j <= length(mov.starts) && @inbounds(mov.starts[j]) <= hi
            # A coarse cell can overlap several coverage intervals; name it once.
            j > seen && (push!(out, j); seen = j)
            j += 1
        end
    end
    return out
end

# --- set arithmetic over the intervals -------------------------------------
#
# Operands are rekeyed to the deeper of the two reference levels, merged as
# sorted intervals, and decomposed by `_cells_from_intervals` into the
# coarsest cells that tile the result — hence normalized minimal answers.

function _same_space(a::MultiOrderVector, b::MultiOrderVector, verb::AbstractString)
    system(a) == system(b) || throw(ArgumentError(
        "cannot $verb multi-order vectors from $(typeof(system(a))) and " *
        "$(typeof(system(b)))"))
    return nothing
end

function Base.union(a::MultiOrderVector, b::MultiOrderVector)
    _same_space(a, b, "union")
    ref = max(a.reference_level, b.reference_level)
    return _cells_from_intervals(a.system,
        _union_intervals(intervals(_rekey(a, ref)), intervals(_rekey(b, ref))), ref)
end

function Base.intersect(a::MultiOrderVector, b::MultiOrderVector)
    _same_space(a, b, "intersect")
    ref = max(a.reference_level, b.reference_level)
    return _cells_from_intervals(a.system,
        _intersect_intervals(intervals(_rekey(a, ref)), intervals(_rekey(b, ref))), ref)
end

function Base.setdiff(a::MultiOrderVector, b::MultiOrderVector)
    _same_space(a, b, "setdiff")
    ref = max(a.reference_level, b.reference_level)
    return _cells_from_intervals(a.system,
        _setdiff_intervals(intervals(_rekey(a, ref)), intervals(_rekey(b, ref))), ref)
end

# Base's n-ary forms go through `union!` on a `Set`, losing order and
# normalization; fold the binary forms instead.
Base.union(a::MultiOrderVector, b::MultiOrderVector, rest::MultiOrderVector...) =
    foldl(union, rest; init=union(a, b))
Base.intersect(a::MultiOrderVector, b::MultiOrderVector, rest::MultiOrderVector...) =
    foldl(intersect, rest; init=intersect(a, b))
Base.setdiff(a::MultiOrderVector, b::MultiOrderVector, rest::MultiOrderVector...) =
    foldl(setdiff, rest; init=setdiff(a, b))

# The predicates Base would otherwise answer by scanning the stored ids, which
# reads a parent and its children as disjoint. Across systems only an empty
# left operand is a subset, as in `CellVector`.
function Base.issubset(a::MultiOrderVector, b::MultiOrderVector)
    system(a) == system(b) || return isempty(a)
    ref = max(a.reference_level, b.reference_level)
    return isempty(_setdiff_intervals(intervals(_rekey(a, ref)),
        intervals(_rekey(b, ref))))
end

function Base.isdisjoint(a::MultiOrderVector, b::MultiOrderVector)
    system(a) == system(b) || return true
    ref = max(a.reference_level, b.reference_level)
    return isempty(_intersect_intervals(intervals(_rekey(a, ref)),
        intervals(_rekey(b, ref))))
end

Base.issetequal(a::MultiOrderVector, b::MultiOrderVector) =
    issubset(a, b) && issubset(b, a)

Base.symdiff(a::MultiOrderVector, b::MultiOrderVector) =
    union(setdiff(a, b), setdiff(b, a))

Base.symdiff(a::MultiOrderVector, b::MultiOrderVector, rest::MultiOrderVector...) =
    foldl(symdiff, rest; init=symdiff(a, b))

"""
    complement(mov::MultiOrderVector) -> MultiOrderVector

Everything `mov` does not hold: the coarsest cells tiling the whole sphere
minus its subtrees, at its own [`reference_level`](@ref).
`complement(complement(mov))` is `mov` normalized.
"""
function complement(mov::MultiOrderVector)
    ref = mov.reference_level
    n = ncells(levelgrid(mov.system, ref))
    return _cells_from_intervals(mov.system,
        _setdiff_intervals([(1, n)], intervals(mov)), ref)
end

# Sorted, disjoint, maximal intervals from a sorted list that may touch or
# overlap. Adjacency counts as overlap, or a parent spanning two touching
# intervals is never found.
function _merged_intervals(ivs::Vector{Tuple{Int,Int}})
    out = Tuple{Int,Int}[]
    for (lo, hi) in ivs
        if !isempty(out) && lo <= last(out)[2] + 1
            out[end] = (last(out)[1], max(last(out)[2], hi))
        else
            push!(out, (lo, hi))
        end
    end
    return out
end

# Both operands are sorted by start, so union is an O(n + m) merge; ties keep
# `A` first.
function _union_intervals(A::Vector{Tuple{Int,Int}}, B::Vector{Tuple{Int,Int}})
    both = Vector{Tuple{Int,Int}}(undef, length(A) + length(B))
    i = j = 1
    for k in eachindex(both)
        take_a = if i > length(A)
            false
        elseif j > length(B)
            true
        else
            @inbounds A[i][1] <= B[j][1]
        end
        if take_a
            @inbounds both[k] = A[i]
            i += 1
        else
            @inbounds both[k] = B[j]
            j += 1
        end
    end
    return _merged_intervals(both)
end

function _intersect_intervals(A::Vector{Tuple{Int,Int}}, B::Vector{Tuple{Int,Int}})
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

function _setdiff_intervals(A::Vector{Tuple{Int,Int}}, B::Vector{Tuple{Int,Int}})
    out = Tuple{Int,Int}[]
    j = 1
    for (lo, hi) in A
        cur = lo
        while j <= length(B) && B[j][2] < cur
            j += 1
        end
        k = j
        while cur <= hi && k <= length(B) && B[k][1] <= hi
            blo, bhi = B[k]
            blo > cur && push!(out, (cur, blo - 1))
            cur = max(cur, bhi + 1)
            k += 1
        end
        cur <= hi && push!(out, (cur, hi))
    end
    return out
end

# Intervals back to cells: at each position, climb ancestors while they start
# there and fit inside the interval, emit the highest, resume past it.
# O(cells × depth); the emitted cells are sorted and disjoint, so the
# container is built unchecked.
function _cells_from_intervals(sys::AbstractHierarchicalGridSystem,
    ivs::Vector{Tuple{Int,Int}}, ref::Int)
    grid = levelgrid(sys, ref)
    top = first(levels(sys))
    cells = cellindextype(sys)[]
    starts, stops = Int[], Int[]
    for (lo, hi) in _merged_intervals(ivs)
        pos = lo
        while pos <= hi
            best = cellindex(grid, pos)
            first_pos, last_pos = pos, pos
            for l in (ref-1):-1:top
                a = ancestor(sys, best, l)
                r = descendant_range(sys, a, ref)
                (first(r) == pos && last(r) <= hi) || break
                best, first_pos, last_pos = a, first(r), last(r)
            end
            push!(cells, best)
            push!(starts, first_pos)
            push!(stops, last_pos)
            pos = last_pos + 1
        end
    end
    return _multiorder_vector(sys, cells, starts, stops, ref, false)
end

# --- equality and geometry -------------------------------------------------

# Same system and same cells. The reference level is the unit the intervals
# are stated in, not part of what is stored, so containers keyed at different
# levels compare equal without rekeying — and intervals alone cannot be
# compared, since an interval names a cell only together with its level.
function Base.:(==)(a::MultiOrderVector, b::MultiOrderVector)
    system(a) == system(b) || return false
    length(a) == length(b) || return false
    return a.cells == b.cells
end

"""
    cell_polygons(mov::MultiOrderVector) -> Vector{<:GI.Polygon}

Every cell as a unit-sphere polygon, in container order, each read from its
own level's grid.
"""
cell_polygons(mov::MultiOrderVector) =
    [cell_polygon(levelgrid(mov.system, level(c)), c) for c in mov.cells]

# --- show ------------------------------------------------------------------

function Base.show(io::IO, mov::MultiOrderVector)
    print(io, "MultiOrderVector(", typeof(mov.system).name.name, ", ",
        length(mov.cells), " cells")
    isempty(mov.cells) || print(io, ", levels ",
        minimum(level, mov.cells), ":", maximum(level, mov.cells))
    print(io, ", ref ", mov.reference_level, ")")
end

Base.show(io::IO, ::MIME"text/plain", mov::MultiOrderVector) = show(io, mov)
