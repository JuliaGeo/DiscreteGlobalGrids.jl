# ---------------------------------------------------------------------------
# `MultiOrderVector` — the mixed-level cell container: cells at different
# refinement levels, pairwise-disjoint subtrees, sorted and indexed by their
# descendant-range intervals at a reference level. The storage form of a MOC,
# where `MultiOrderCellSet` is the query form.
#
# The acceleration is the flat sorted interval index and nothing else — three
# `Int` vectors beside the ids. Point, cell and region queries are one binary
# search over `starts` plus a bound test against `stops`, which is how
# astronomy MOCs are stored and searched; there is no tree to build, to balance
# or to keep in step with the cells. `treeify` over one is a future extension,
# not a missing piece.
#
# Everything here is interval arithmetic, so everything here requires
# `has_sorted_subtrees(sys)`: without descendant ranges a subtree is not one
# interval of its level and there is nothing to sort by. A5 throws, in the same
# style `level_ranges` does.
#
# Contract: docs/design/moc-storage.md §1.
# ---------------------------------------------------------------------------

"""
    MultiOrderVector(set::MultiOrderCellSet)
    MultiOrderVector(sys, cells::AbstractVector; reference_level = deepest cell level)

An immutable `AbstractVector` of **cells at mixed refinement levels** of one
hierarchical system, pairwise disjoint as subtrees — no member is an ancestor
of another, and none repeats — held in the order of the position INTERVALS
their subtrees occupy at one **reference level**.

Semantically `mov` **is** that id vector: `length(mov)` is the number of cells,
`mov[k]` is the `k`th of them, `collect(mov)` is the vector itself. What is
*stored* beside the ids is the interval index — `starts`, `stops` and the
cumulative interval lengths `offsets` — which is what answers every query below
in O(log n).

This is the axis of an adaptively refined mesh: one value per cell, cells as
coarse as the data allows. [`MultiOrderCellSet`](@ref) is the same mixed-level
shape as the *result of a query*; this is the shape that data hangs off, and
the two convert in one call.

# The two ways in

  - a [`MultiOrderCellSet`](@ref) — a coverage, read as storage. Its cells are
    already sorted by descendant-range start and already disjoint, so this is
    O(n) and keeps the set's own reference level;
  - `sys, cells` — an explicit cell vector in any order, sorted and validated
    here. The `reference_level` keyword defaults to the deepest level present,
    which is the shallowest level that can key every cell; a deeper one is
    legal and changes nothing but the integers.

The reference level is a property of the CONTAINER, not of the system: it is
the depth the intervals speak about, and two containers with different
reference levels can still be equal.

# The verbs

```julia
mov[k]                          # the kth cell id
cellposition(mov, c)            # EXACT membership: Int, or `nothing`
c in mov                        # the same question as a Bool
covering_position(mov, c)       # the stored cell that IS `c` or an ANCESTOR of it
cellat(mov, lon, lat)           # the stored cell a point falls in
covering(mov, polygon)          # the sub-container a region names
CellVector(mov; level = l)      # the leaf-level expansion, as windows
union(a, b); intersect(a, b); setdiff(a, b); complement(a)
```

[`cellposition`](@ref) and [`covering_position`](@ref) are the two halves that
make this a *storage* type rather than a set: the first asks whether a cell is
stored, the second asks which stored cell speaks for it. Compression reads the
second — a leaf resolves to the ancestor that holds its value.

# Complexity

Construction is O(n) from a set and O(n log n) from a loose vector (the sort).
`cellposition`, `covering_position` and the point forms are one binary search,
O(log n), allocation-free. [`covering`](@ref) is O(log n) per interval of the
region's coverage plus the coverage query itself. The set operations are
O(n + m) interval arithmetic followed by an O(cells × depth) decomposition into
the coarsest cells that tile the result.

# What the set operations return

**Normalized minimal** containers: the coarsest cells whose subtrees tile the
result, which is what astronomy MOC libraries do too. A consequence worth
stating rather than discovering — `union(a, a)` may be *coarser* than `a` if
`a` was not itself minimal, because a complete sibling family in `a` comes back
as its parent. The two name the same leaves, and `==` says so.
"""
struct MultiOrderVector{ID,S<:AbstractHierarchicalGridSystem} <: AbstractVector{ID}
    system::S
    cells::Vector{ID}
    starts::Vector{Int}      # first(descendant_range(sys, c, reference_level))
    stops::Vector{Int}       # ... and its inclusive last
    offsets::Vector{Int}     # offsets[i] = leaves at the reference level in cells 1:i
    reference_level::Int

    # The invariants live here so that no route in can dodge them, and
    # `checked = false` is the escape hatch for the routes that have just
    # established them — a subset of a valid container, a rekeying, the
    # decomposition of disjoint intervals.
    function MultiOrderVector{ID,S}(system::S, cells::Vector{ID}, starts::Vector{Int},
        stops::Vector{Int}, offsets::Vector{Int}, ref::Int,
        checked::Bool) where {ID,S<:AbstractHierarchicalGridSystem}
        checked && _check_multiorder(system, cells, starts, stops, ref)
        return new{ID,S}(system, cells, starts, stops, offsets, ref)
    end
end

# The one assembly point: `offsets` is derived, never passed in, so it cannot
# disagree with the intervals it summarises.
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

# Split from the disjointness test because the constructors need it FIRST: both
# resolve `descendant_range`, which is a `MethodError` without the trait and an
# `ArgumentError` above the cell's own level.
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
    for i in eachindex(starts)
        starts[i] <= stops[i] || throw(ArgumentError(
            "the interval of $(cells[i]) at level $ref is empty: " *
            "$(starts[i]):$(stops[i])"))
        i == 1 && continue
        # Strictly ascending starts with no overlap IS "no member is an ancestor
        # of another, and none repeats": two cells share positions at the
        # reference level exactly when one contains the other.
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

# The keyword shadows the `reference_level` accessor, so the work is a
# positional `Int` one call in — the `_multi_order` pattern.
MultiOrderVector(sys::AbstractHierarchicalGridSystem, cells::AbstractVector;
    reference_level::Integer=_default_reference_level(sys, cells)) =
    _multiorder_from_cells(sys, cells, Int(reference_level))

_default_reference_level(sys::AbstractHierarchicalGridSystem, cells) =
    isempty(cells) ? first(levels(sys)) : maximum(level, cells)

function _multiorder_from_cells(sys::AbstractHierarchicalGridSystem, cells, ref::Int)
    ID = eltype(cells) <: AbstractCellIndex ? eltype(cells) : cellindextype(sys)
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

# Ascending indices keep the interval index; anything else — a permutation, a
# repeat, a reversal — is not a sorted disjoint interval list and is answered
# with the ordinary `Vector` that can hold it. Same fork, and the same
# rationale, as `CellVector`: `AbstractVector` rather than `AbstractArray`, so
# that a shaped index falls through to Base's generic and comes back with the
# INDEX's shape whatever its values.
Base.getindex(mov::MultiOrderVector, idx::AbstractVector{<:Integer}) = _subset(mov, idx)

# The SmallCollections tie-break `cell_vector.jl` documents: that package's own
# `getindex(::AbstractVector, ::AbstractFixedOrSmall...)` is neither more nor
# less specific than the line above.
Base.getindex(mov::MultiOrderVector,
    i::SmallCollections.AbstractFixedOrSmallOrPackedVector{<:Integer}) = _subset(mov, i)

# A mask names a position by its INDEX, so a mask of the wrong length or rank is
# a bounds error rather than a shorter answer — `findall` cannot tell.
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

The level the container's intervals are stated at: no shallower than any stored
cell, and the level [`covering_position`](@ref) locates points at.

It is a property of the container rather than of the system — two containers of
the same cells at different reference levels hold the same cells and compare
equal. Unexported; reach it as
`DiscreteGlobalGrids.Fallbacks.reference_level(mov)`.
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

Position of a cell in the container, or `nothing` when it does not hold that
cell — including when `c` is a descendant or an ancestor of one it does hold.
This is EXACT membership and the inverse of `mov[k]`: one binary search over
the interval starts, O(log n).

[`covering_position`](@ref) is the other question, and the one compression
asks: which stored cell speaks for `c`.

The point forms are [`cellat`](@ref)'s answer read as a position, so they *are*
the covering question — a point falls inside a stored cell rather than naming
one. Degrees for the `(lon, lat)` method, as everywhere else.
"""
function cellposition(mov::MultiOrderVector, c::AbstractCellIndex)
    lc = level(c)
    lc <= mov.reference_level || return nothing
    r = descendant_range(mov.system, c, mov.reference_level)
    j = searchsortedfirst(mov.starts, first(r))
    j <= length(mov.starts) || return nothing
    # Same interval and same level is the same cell: positions of a level are a
    # bijection onto its cells, so nothing else can occupy both.
    (@inbounds(mov.starts[j]) == first(r) && @inbounds(mov.stops[j]) == last(r) &&
     level(@inbounds mov.cells[j]) == lc) || return nothing
    return j
end

"""
    covering_position(mov::MultiOrderVector, c::AbstractCellIndex) -> Union{Int,Nothing}

The position of the stored cell that **is `c`, or is an ancestor of `c`** — the
compression verb, which resolves a leaf to the cell holding its value.
`nothing` when no stored cell covers `c`.

`c` may be deeper than the container's [`reference_level`](@ref): it is keyed
through its reference-level ancestor, which has the same covering cell because
every stored cell is at or above that level.

O(log n): `searchsortedlast` over the interval starts, then one bound test
against that interval's stop. [`cellposition`](@ref) is the exact question.
"""
function covering_position(mov::MultiOrderVector, c::AbstractCellIndex)
    ref = mov.reference_level
    keyed = level(c) > ref ? ancestor(mov.system, c, ref) : c
    r = descendant_range(mov.system, keyed, ref)
    j = searchsortedlast(mov.starts, first(r))
    j >= 1 || return nothing
    return last(r) <= @inbounds(mov.stops[j]) ? j : nothing
end

# Membership is the interval search, which is exact and O(log n) where the
# generic scan over the values is O(n) id comparisons.
Base.in(c::AbstractCellIndex, mov::MultiOrderVector) = cellposition(mov, c) !== nothing

"""
    cellat(mov::MultiOrderVector, p::GO.UnitSphericalPoint) -> Union{AbstractCellIndex,Nothing}
    cellat(mov::MultiOrderVector, lon::Real, lat::Real)

The cell **of `mov`** containing the point, at whatever level that cell sits —
or `nothing` when no stored cell covers it. The point is located once in
`levelgrid(system(mov), reference_level(mov))` and the interval index does the
rest, so this is one point location plus one binary search whatever the mixture
of levels.

[`cellposition`](@ref)'s point forms are the same question answered as a
position.
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

The container expanded to single-level ids: the leaf cells its subtrees name at
`level`, as the windowed [`CellVector`](@ref) — and from there a `PartialGrid`,
a regridder, a halo table or a cube axis.

`level` must be no shallower than the deepest stored cell. The windows are the
stored cells' `descendant_range`s at that level, merged where adjacent, which
is exactly what [`level_ranges`](@ref) produces for the [`MultiOrderCellSet`](@ref)
the container came from — so the two agree window for window, at every level.

O(n) whatever `level` is: expanding three levels deeper multiplies the cells
named and moves nothing that is stored. [`cellset`](@ref) on the result answers
with the container it came from, which is the provenance a re-expansion at
another level needs.
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

# `level_ranges`' merge, over the container's cells: adjacent subtrees are one
# window, which is what makes the expansion canonical and comparable.
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

The cells of `mov` that a [`MultiOrderCoverage`](@ref) of `target` names, as a
container again. A stored cell is kept **whole** whenever its subtree meets the
coverage — never clipped to it — because a mixed-level cell is one datum and
half of one is not a cell.

`target` is anything [`query`](@ref) accepts: a GeoInterface geometry, an
`Extents.Extent` in lon/lat degrees, or a `GO.UnitSpherical.SphericalCap`.

[`covering_positions`](@ref) is the position-space form, for indexing a data
array laid out against `mov`. The coverage is run at
[`reference_level`](@ref)`(mov)` and matched against the interval index, so the
selection costs O(log n) per interval of the coverage rather than a walk over
its leaves.
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
            # One stored cell can overlap several coverage intervals — a coarse
            # cell spanning a hole in the target does — and must be named once.
            j > seen && (push!(out, j); seen = j)
            j += 1
        end
    end
    return out
end

# --- set arithmetic over the intervals -------------------------------------
#
# The operands are rekeyed to the deeper of the two reference levels, so that
# both speak in the same units, and the arithmetic is then the ordinary sorted
# interval merge. What is not ordinary is the way back: an interval is not a
# cell, and the answer has to be cells again. `_cells_from_intervals` decomposes
# each interval into the COARSEST cells that tile it, which is what makes the
# result normalized and minimal — and what makes `union(a, a)` coarser than a
# non-minimal `a`, since a complete sibling family comes back as its parent.

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

# Base folds its own n-ary forms through `union!`/`intersect!`, which build a
# `Set` and lose both the order and the normalization. Folding the binary form
# keeps the answer a container.
Base.union(a::MultiOrderVector, b::MultiOrderVector, rest::MultiOrderVector...) =
    foldl(union, rest; init=union(a, b))
Base.intersect(a::MultiOrderVector, b::MultiOrderVector, rest::MultiOrderVector...) =
    foldl(intersect, rest; init=intersect(a, b))
Base.setdiff(a::MultiOrderVector, b::MultiOrderVector, rest::MultiOrderVector...) =
    foldl(setdiff, rest; init=setdiff(a, b))

"""
    complement(mov::MultiOrderVector) -> MultiOrderVector

Everything `mov` does not hold: the coarsest cells tiling the whole sphere
minus its subtrees, at its own [`reference_level`](@ref).

```julia
ocean = complement(land)
isempty(intersect(land, ocean))          # true, at every level
```

`complement(complement(mov))` is `mov` NORMALIZED — the same leaves as the same
coarsest cells — which equals `mov` itself whenever `mov` was already minimal.
"""
function complement(mov::MultiOrderVector)
    ref = mov.reference_level
    n = ncells(levelgrid(mov.system, ref))
    return _cells_from_intervals(mov.system,
        _setdiff_intervals([(1, n)], intervals(mov)), ref)
end

# Sorted, disjoint, maximal intervals from a sorted list that may touch or
# overlap. Adjacency counts as overlap: two neighbouring intervals must become
# one before decomposition, or the shared parent is never found.
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

# Both operands are already sorted by start, so the union is their MERGE — the
# O(n + m) the type's complexity note promises, where sorting the concatenation
# would have been O((n + m) log(n + m)) over data that was never out of order.
# Ties keep `A` first, which is the stable sort's answer too.
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

# The way back from intervals to cells, and the reason the set operations
# normalize. Walking an interval left to right, the cell at the current position
# is climbed as far as its ancestors keep starting there and keep fitting — the
# coarsest cell that tiles the interval from here — and the walk resumes past
# it. O(cells × depth), and the emitted cells are sorted and disjoint by
# construction, so the container is built unchecked.
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

# Same system and the same cells — the CELLS, which is the question `mov`
# semantically is its id vector answers, and the one `cellposition` answers too:
# an interval names a cell only together with the level it is read at, which is
# why that function tests the level as well and why this one cannot test the
# intervals alone.
#
# The reference level does not enter, and that is the point rather than an
# omission: it is the unit the intervals are stated in, not part of what is
# stored, so two containers of the same cells keyed at different levels compare
# equal here without either being re-keyed first.
function Base.:(==)(a::MultiOrderVector, b::MultiOrderVector)
    system(a) == system(b) || return false
    length(a) == length(b) || return false
    return a.cells == b.cells
end

"""
    cell_polygons(mov::MultiOrderVector) -> Vector{<:GI.Polygon}

Every cell of the container as a unit-sphere polygon, in its own order — what a
plot of an adaptively refined mesh needs, in one call. Each cell is read from
its own level's grid, so the mixture of levels costs the caller nothing.
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
