# ---------------------------------------------------------------------------
# Multi-order coverage
#
# The other shape of a spatial query: instead of every cell at one level, the
# COARSEST cells that describe a region — a breadth-first walk from the roots
# that emits a cell whole once it is entirely inside the target, and recurses
# only where the target's boundary crosses it.
#
# The result is sorted in **space-filling-curve order**: each cell owns a
# disjoint position interval at a reference depth (`descendant_range`), and
# ordering by that interval's start is depth-first curve order. That is what
# generalises HEALPix's MOC/NUNIQ, and it buys sibling compaction, binary-search
# membership, and lazy expansion to any level as sorted disconnected ranges —
# the handshake the lookup layer consumes.
# ---------------------------------------------------------------------------

"""
    MultiOrderCoverage(target)

A multi-order coverage query: hand it to [`query`](@ref) with a system and a
maximum depth to get a [`MultiOrderCellSet`](@ref).

```julia
set = query(sys, MultiOrderCoverage(polygon); level = 8)
```

`target` takes the same forms as any other query target — a GeoInterface
geometry, an `Extents.Extent` in lon/lat degrees, or a
`GO.UnitSpherical.SphericalCap`. `Base.parent` unwraps it, as it does for a
DE9IM predicate.

The traversal emits a cell as soon as the cell lies entirely inside the target,
and recurses into the children of a cell the target's boundary crosses, down to
the requested level; cells still crossed at that level are emitted too, so the
set **covers** the target rather than being covered by it.

!!! note "Coverage is a statement about cells, not about the union of their descendants"
    A cell is emitted because the *cell* is inside the target. In a system
    where children overhang their parent, expanding that cell to a deeper level
    yields descendants that may poke outside — see [`level_ranges`](@ref).
"""
struct MultiOrderCoverage{T}
    target::T
end

Base.parent(coverage::MultiOrderCoverage) = coverage.target

Base.show(io::IO, coverage::MultiOrderCoverage) =
    print(io, "MultiOrderCoverage(", typeof(coverage.target).name.name, ")")

"""
    MultiOrderCellSet

A set of cells at **mixed levels**, in space-filling-curve order — the result of
a [`MultiOrderCoverage`](@ref) query.

Iterating it yields the typed cell ids, coarsest-first within each branch and
in curve order overall (`length`, `getindex`, `eltype` and `collect` all work).
[`level_ranges(set, l)`](@ref level_ranges) expands it to one level as sorted,
disjoint position ranges.

Order is by the start of each cell's `descendant_range` at the set's reference
level, which is exactly depth-first curve order and makes sibling intervals
adjacent. A system without [`has_sorted_subtrees`](@ref) has no such intervals,
and falls back to `(level, id)` order.
"""
struct MultiOrderCellSet{S<:AbstractHierarchicalGridSystem,ID}
    system::S
    cells::Vector{ID}
    keys::Vector{Int}
    reference_level::Int
end

"""
    MultiOrderCellSet(sys, coverage::MultiOrderCoverage; level)

Run a [`MultiOrderCoverage`](@ref) against `sys`, recursing no deeper than
`level`. Equivalent to `query(sys, coverage; level)`.
"""
MultiOrderCellSet(sys::AbstractHierarchicalGridSystem, coverage::MultiOrderCoverage;
    level::Integer) = _multi_order(sys, coverage.target, Int(level))

"""
    query(sys, coverage::MultiOrderCoverage; level) -> MultiOrderCellSet

The multi-order form of [`query`](@ref): the coarsest cells covering the
target, down to `level`.
"""
query(sys::AbstractHierarchicalGridSystem, coverage::MultiOrderCoverage; level::Integer) =
    _multi_order(sys, coverage.target, Int(level))

# The keyword `level` shadows the `level` function, so the whole traversal
# lives here, where the maximum depth is a plain positional `Int`.
function _multi_order(sys::AbstractHierarchicalGridSystem, target_value, maxlevel::Int)
    maxlevel in levels(sys) || throw(ArgumentError(
        "level $maxlevel is outside $(typeof(sys))'s levels $(levels(sys))"))
    target = _query_target(target_value)
    cells = cellindextype(sys)[]
    # One level grid per level, built once rather than per visited cell: the
    # traversal touches every level from the roots down, and `levelgrid` is
    # cheap but not free.
    top = first(levels(sys))
    grids = [levelgrid(sys, l) for l in top:maxlevel]
    for c in rootcells(sys)
        _coverage_visit!(cells, sys, target, c, maxlevel, grids, top)
    end
    return _sorted_cell_set(sys, cells, maxlevel)
end

function _coverage_visit!(cells, sys, target, c, maxlevel::Int, grids, top::Int)
    # The subtree prune, under the covering law: a cell whose node extent misses
    # the target has no descendant that can meet it.
    intersects_cap(target.cap, node_extent(sys, c)) || return nothing
    lc = level(c)
    grid = grids[lc-top+1]
    # Intersection first — it has the sandwich fast path, and a cell that does
    # not even meet the target is not worth a containment predicate.
    _matches(DE9IM.Intersects(nothing), target, grid, c) || return nothing
    if _matches(DE9IM.Within(nothing), target, grid, c)
        push!(cells, c)                      # entirely inside: emit whole
        return nothing
    end
    if lc >= maxlevel
        push!(cells, c)                      # still crossed at the deepest level
        return nothing
    end
    for child in children(sys, c)
        _coverage_visit!(cells, sys, target, child, maxlevel, grids, top)
    end
    return nothing
end

function _sorted_cell_set(sys::AbstractHierarchicalGridSystem, cells::Vector{ID},
        reference_level::Int) where {ID}
    if has_sorted_subtrees(sys)
        keys = [first(descendant_range(sys, c, reference_level)) for c in cells]
        perm = sortperm(keys)
        return MultiOrderCellSet{typeof(sys),ID}(sys, cells[perm], keys[perm],
            reference_level)
    end
    # No curve intervals to order by; `(level, id)` is the documented fallback,
    # and the keys become the cells' own positions within their level, which is
    # still a total order but not a curve order.
    perm = sortperm(cells; by=c -> (level(c), c))
    ordered = cells[perm]
    keys = [something(cellposition(levelgrid(sys, level(c)), c), 0) for c in ordered]
    return MultiOrderCellSet{typeof(sys),ID}(sys, ordered, keys, reference_level)
end

# --- the collection surface ------------------------------------------------

Base.length(set::MultiOrderCellSet) = length(set.cells)
Base.eltype(::Type{MultiOrderCellSet{S,ID}}) where {S,ID} = ID
Base.eltype(set::MultiOrderCellSet) = eltype(typeof(set))
Base.isempty(set::MultiOrderCellSet) = isempty(set.cells)
Base.iterate(set::MultiOrderCellSet, state...) = iterate(set.cells, state...)
Base.getindex(set::MultiOrderCellSet, i::Int) = set.cells[i]
Base.firstindex(::MultiOrderCellSet) = 1
Base.lastindex(set::MultiOrderCellSet) = length(set.cells)
Base.collect(set::MultiOrderCellSet) = copy(set.cells)

"""
    system(set::MultiOrderCellSet)

The system the set's cells are named in.
"""
system(set::MultiOrderCellSet) = set.system

"""
    curve_keys(set::MultiOrderCellSet) -> Vector{Int}

The sort key of each cell, ascending: the start of its position interval at the
set's reference level. Sibling intervals are adjacent, which is what makes
compaction and binary-search membership cheap.
"""
curve_keys(set::MultiOrderCellSet) = set.keys

function Base.show(io::IO, set::MultiOrderCellSet)
    print(io, "MultiOrderCellSet(", typeof(set.system).name.name, ", ",
        length(set.cells), " cells")
    isempty(set.cells) || print(io, ", levels ",
        minimum(level, set.cells), ":", maximum(level, set.cells))
    print(io, ")")
end

Base.show(io::IO, ::MIME"text/plain", set::MultiOrderCellSet) = show(io, set)

"""
    level_ranges(set::MultiOrderCellSet, l::Integer) -> Vector{UnitRange{Int}}

The set expanded to level `l`, as **sorted, disjoint position ranges** in
`levelgrid(sys, l)` — the form a lookup layer slices data arrays with.

Adjacent ranges are merged, so a set whose cells happen to be a compacted
sibling group comes back as one range rather than as its parts.

Requires [`has_sorted_subtrees`](@ref) (there are no position intervals
otherwise) and `l` at least as deep as every cell in the set: expanding to a
coarser level would have to replace a cell by an ancestor, which covers more
than the set does.
"""
function level_ranges(set::MultiOrderCellSet, l::Integer)
    has_sorted_subtrees(set.system) || throw(ArgumentError(
        "$(typeof(set.system)) has no descendant ranges, so a multi-order set " *
        "cannot be expanded to position ranges"))
    target = Int(l)
    out = UnitRange{Int}[]
    for c in set.cells
        level(c) <= target || throw(ArgumentError(
            "cannot expand to level $target: the set contains a level-$(level(c)) cell"))
        r = descendant_range(set.system, c, target)
        if !isempty(out) && first(r) == last(out[end]) + 1
            out[end] = first(out[end]):last(r)
        else
            push!(out, r)
        end
    end
    return out
end

"""
    cellindices(set::MultiOrderCellSet, l::Integer) -> Vector{<:AbstractCellIndex}

The set expanded to level `l` as typed ids, ascending — [`level_ranges`](@ref)
resolved through `cellindex`. O(cells at `l`), so reach for the ranges instead
wherever the positions are what is wanted.
"""
function cellindices(set::MultiOrderCellSet, l::Integer)
    grid = levelgrid(set.system, Int(l))
    out = cellindextype(set.system)[]
    for r in level_ranges(set, l), i in r
        push!(out, cellindex(grid, i))
    end
    return out
end
