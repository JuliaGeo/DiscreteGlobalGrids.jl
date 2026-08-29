# Multi-order queries return the coarsest cells covering a region, ordered by
# descendant-range start at a reference level — the generalisation of HEALPix's
# MOC/NUNIQ.

"""
    MultiOrderCoverage(target)

A multi-order coverage query for [`query`](@ref), bounded by a maximum depth or
by a cell budget.

```julia
set = query(sys, MultiOrderCoverage(polygon); level = 8)      # accuracy first
set = query(sys, MultiOrderCoverage(polygon); maxcells = 10)  # cardinality first
```

`target` accepts a GeoInterface geometry, lon/lat `Extents.Extent`, or
`GO.UnitSpherical.SphericalCap`.

# The two modes

`level` is accuracy first. Traversal is depth first: it emits a cell contained
by the target and recurses through boundary crossings, to the requested level;
cells still crossed there are emitted too, so the set covers the target rather
than being covered by it. The cardinality is whatever the outline needs.

`maxcells` is cardinality first — "ten cells that cover California, or a
hundred". Refinement is coarsest first over the crossing cells, a level at a
time; a cell whose replacement would not fit is kept whole and the search moves
on. A seed already larger than the budget is the one set returned over it. The
depth is then whatever the budget bought, and varies from branch to branch.

Neither mode approximates the other: a `level` set is the exact answer at a
fixed depth, a `maxcells` set the best a fixed cardinality can say, with the
deepest level it reached as its reference level. The keywords are mutually
exclusive; [`query`](@ref) states the rules.

[`iscontained`](@ref) reports which emissions were *proven* to fit — a cell
emitted at the deepest level is never asked; [`coarsest_contained`](@ref) is the
accessor that uses it.

!!! warning "Coverage is a statement about the leaves, not about the drawn cells"
    The guarantee is at the deepest level: every cell there that meets the
    target is a member of the set or a descendant of one, and no member
    descends from another. [`level_ranges`](@ref) therefore contains the
    single-level `Intersects` query, and equals it where refinement is
    congruent.

    The union of the emitted polygons is not that region. Replacing a subtree
    by its root swaps the subtree's footprint for the root's, and the two agree
    only under congruent refinement: HEALPix, S2 and ISEA4R tile; on a
    state-sized target the slivers cover under 2% of it on IGEO7, 15% on H3 and
    30% on A5 — the bounds [`query`](@ref) states and its suite asserts. Draw a
    set as *which cells were chosen*; expand it before computing with it as a
    region.

    The same non-congruence runs the other way: a member's descendants may lie
    outside the target, inside a hole in it for instance, so the expansion
    over-covers exactly where the refinement does.
"""
struct MultiOrderCoverage{T}
    target::T
end

Base.parent(coverage::MultiOrderCoverage) = coverage.target

Base.show(io::IO, coverage::MultiOrderCoverage) =
    print(io, "MultiOrderCoverage(", typeof(coverage.target).name.name, ")")

"""
    MultiOrderCellSet

A mixed-level cell set returned by [`MultiOrderCoverage`](@ref). Iteration yields
typed ids in descendant-range order at the reference level; systems without
sorted subtrees use `(level, id)` order. [`level_ranges`](@ref level_ranges)
expands it to sorted, disjoint ranges at one level.

[`iscontained`](@ref) reports which cells were *proven* to fit inside the
target — not which ones do; that docstring draws the line.
[`cell_polygon`](@ref) and [`cell_polygons`](@ref) read mixed-level geometry
without the caller resolving a level grid per cell.

The REFERENCE LEVEL is the depth the set speaks about: the `level` the query was
given, or — in `maxcells` mode, where no depth was asked for — the deepest level
the budget reached. It is the default expansion level for
[`level_ranges`](@ref), [`cellindices`](@ref) and `CellLookup`, and the level at
which the covering guarantee is stated.

!!! note "Expansion needs sorted subtrees"
    `level_ranges` throws where [`has_sorted_subtrees`](@ref) is `false` (A5),
    because a cell's descendants do not occupy one interval of their level.
    `descendants(sys, c, l)` still names them, as a list rather than ranges.
"""
struct MultiOrderCellSet{S<:AbstractHierarchicalGridSystem,ID}
    system::S
    cells::Vector{ID}
    keys::Vector{Int}
    contained::BitVector
    reference_level::Int
end

"""
    MultiOrderCellSet(sys, coverage::MultiOrderCoverage; level)
    MultiOrderCellSet(sys, coverage::MultiOrderCoverage; maxcells, maxlevel = deepest)

Run a [`MultiOrderCoverage`](@ref) against `sys`, recursing no deeper than
`level`, or refining until `maxcells` cells are spent. Equivalent to
`query(sys, coverage; ...)`, which documents both modes.
"""
MultiOrderCellSet(sys::AbstractHierarchicalGridSystem, coverage::MultiOrderCoverage;
    level::Union{Integer,Nothing}=nothing, maxcells::Union{Integer,Nothing}=nothing,
    maxlevel::Union{Integer,Nothing}=nothing) =
    _multi_order_query(sys, coverage.target, level, maxcells, maxlevel)

"""
    query(sys, coverage::MultiOrderCoverage; level) -> MultiOrderCellSet
    query(sys, coverage::MultiOrderCoverage; maxcells, maxlevel = deepest) -> MultiOrderCellSet

The multi-order form of [`query`](@ref): the coarsest cells covering the target.
The two keywords are the two modes of [`MultiOrderCoverage`](@ref), and exactly
one of them must be given.

`level` is ACCURACY FIRST: refine everything the target's boundary crosses down
to `level`, and let the cell count fall where it may.

```julia
set = query(sys, MultiOrderCoverage(california); level = 7)   # thousands of cells
```

`maxcells` is CARDINALITY FIRST: refine the crossing cells coarsest first,
keeping whole any cell whose replacement would not fit.

```julia
set = query(sys, MultiOrderCoverage(california); maxcells = 10)   # ten cells
```

`maxlevel` bounds how deep the budget may descend, and defaults to the system's
deepest level, so that the budget alone decides. It is worth setting on a target
much smaller than a root cell, where refinement replaces one crossing cell by one
crossing cell for level after level and the budget never binds.

# Edge cases, all of them by design

  - Giving both keywords, or neither, is an `ArgumentError`: they are modes, not
    a bound and a hint. `maxlevel` belongs with `maxcells`; in `level` mode the
    level IS the bound.
  - `maxcells < 1` is an `ArgumentError`. A covering of a target the system
    meets at all needs at least one cell.
  - **A seed larger than the budget is returned whole**, over budget, rather
    than throwing. The seed is the set of coarsest cells that meet the target,
    and it is the smallest covering this traversal can express; there is nothing
    to refine away. A target spanning most of the sphere with `maxcells = 3`
    lands here. `length(set) <= maxcells` holds in every other case.
  - A target smaller than one cell of the seed comes back as that one crossing
    cell when the budget is 1, and as a chain of single crossing cells descending
    towards it as the budget grows.

# What "cover" means here, and where it is exact

Every point of the target lies inside one of the emitted cells. The seed is
the coarsest cells meeting the target; refinement replaces a crossing cell by
its meeting children. Non-congruent refinement descends through missing cells
too, as `level` mode does, so a crossing cell with no meeting child is dropped
for the cells that do cover its share — and its reserved slot pays for what the
walk still owes: cells found and not afforded, or the cell itself kept whole
where its share was never certified.

It is EXACT, at every budget, on the three systems whose four children tile
their parent: HEALPix, S2 and ISEA4R. Where children do not tile their parent it
degrades in the way [`MultiOrderCoverage`](@ref)'s warning already describes, and
for the same reason — replacing a cell by its children swaps one footprint for
another. Measured on a state-sized outline as the fraction of the target lying
in no emitted cell: under 2% on IGEO7, 3% on its authalic wrap, 15% on H3 and
30% on A5.

The LEAF statement `level` mode makes — every reference-level cell meeting the
target is a member or the descendant of one — is a law here on those same three
systems only. Elsewhere the budget makes the same overhang descent but has no
fixed depth to carry it to: where it stops paying, the search goes on only for
the dropped cells whose share is still unproven. On the same outline the leaf
statement misses under 1% of the target on IGEO7, 2% on its authalic wrap and
on H3, and 18% on A5. `test/systems/crosssystem/multiorder_budget.jl` asserts
both statements' bounds per system, at three budgets and on four targets.

What a budget does NOT buy is a tight picture of the target: at ten cells the
set over-covers California by a wide margin, and [`iscontained`](@ref) says so.

# Composition

A budget set is a `MultiOrderCellSet` like any other. It sorts in curve order,
answers [`coarsest_contained`](@ref), [`cell_polygons`](@ref) and
[`level_ranges`](@ref), and backs a `CellLookup` at its own reference level or
at any deeper one — so ten cells chosen for the picture can still name every
leaf under them for the data.
"""
query(sys::AbstractHierarchicalGridSystem, coverage::MultiOrderCoverage;
    level::Union{Integer,Nothing}=nothing, maxcells::Union{Integer,Nothing}=nothing,
    maxlevel::Union{Integer,Nothing}=nothing) =
    _multi_order_query(sys, coverage.target, level, maxcells, maxlevel)

# One place where the two modes are told apart, so that both entry points give
# the same errors for the same keyword combinations.
function _multi_order_query(sys::AbstractHierarchicalGridSystem, target_value,
        level_kw, maxcells, maxlevel)
    if level_kw !== nothing
        maxcells === nothing || throw(ArgumentError(
            "`level` and `maxcells` are the two modes of a multi-order coverage and " *
            "cannot both be given: `level` refines to a fixed depth, `maxcells` " *
            "spends a cell budget"))
        maxlevel === nothing || throw(ArgumentError(
            "`maxlevel` bounds the budget traversal and belongs with `maxcells`; " *
            "in `level` mode the level is the bound"))
        return _multi_order(sys, target_value, Int(level_kw))
    end
    maxcells === nothing && throw(ArgumentError(
        "a multi-order coverage needs one of `level` (refine to a fixed depth) or " *
        "`maxcells` (spend a cell budget)"))
    cap = maxlevel === nothing ? last(levels(sys)) : Int(maxlevel)
    return _multi_order_budget(sys, target_value, Int(maxcells), cap)
end

# The keyword `level` shadows the `level` function, so the whole traversal
# lives here, where the maximum depth is a plain positional `Int`.
function _multi_order(sys::AbstractHierarchicalGridSystem, target_value, maxlevel::Int)
    maxlevel in levels(sys) || throw(ArgumentError(
        "level $maxlevel is outside $(typeof(sys))'s levels $(levels(sys))"))
    target = _query_target(target_value)
    cells = cellindextype(sys)[]
    # Parallel to `cells`: `true` exactly where `Within` was asked and held.
    # Emissions at `maxlevel` are never asked, so `false` there means unproven,
    # not outside. `iscontained` documents the asymmetry.
    contained = BitVector()
    # One level grid per level, built once rather than per visited cell: the
    # traversal touches every level from the roots down, and `levelgrid` is
    # cheap but not free.
    top = first(levels(sys))
    grids = [levelgrid(sys, l) for l in top:maxlevel]
    for c in rootcells(sys)
        _coverage_visit!(cells, contained, sys, target, c, maxlevel, grids, top)
    end
    return _sorted_cell_set(sys, cells, contained, maxlevel)
end

function _coverage_visit!(cells, contained, sys, target, c, maxlevel::Int, grids, top::Int)
    # Only `node_extent` may prune descendants; child geometry can overhang its
    # parent. Exact cell geometry determines emission, not descent.
    #
    # Both prunes read that extent: the target's cap first, at one distance,
    # then the boundary-arc proof, which is what keeps the traversal
    # output-sensitive when the cap is loose or the whole sphere.
    extent = node_extent(sys, c)
    Extents.intersects(target.cap, extent) || return nothing
    _subtree_outside(target, extent) && return nothing
    lc = level(c)
    grid = grids[lc-top+1]
    meets = _matches(DE9IM.Intersects(nothing), target, grid, c)
    if lc >= maxlevel
        # Emitted to cover, and flagged unproven WITHOUT asking `Within`: many
        # of these cells do fit, and asking costs one ~48 KB call per boundary
        # cell to label cells the traversal is done with. See `iscontained`.
        meets && (push!(cells, c); push!(contained, false))
        return nothing
    end
    # Containment is asked only of a cell already known to meet the target: it
    # is the expensive predicate, and it has no fast path of its own.
    if meets && _matches(DE9IM.Within(nothing), target, grid, c)
        push!(cells, c)                      # entirely inside: emit whole
        push!(contained, true)
        return nothing
    end
    for child in children(sys, c)
        _coverage_visit!(cells, contained, sys, target, child, maxlevel, grids, top)
    end
    return nothing
end


function _sorted_cell_set(sys::AbstractHierarchicalGridSystem, cells::Vector{ID},
        contained::BitVector, reference_level::Int) where {ID}
    if has_sorted_subtrees(sys)
        keys = [first(descendant_range(sys, c, reference_level)) for c in cells]
        perm = sortperm(keys)
        return MultiOrderCellSet{typeof(sys),ID}(sys, cells[perm], keys[perm],
            contained[perm], reference_level)
    end
    # No curve intervals to order by; `(level, id)` is the documented fallback,
    # and the keys become the cells' own indices within their level, which is
    # still a total order but not a curve order.
    perm = sortperm(cells; by=c -> (level(c), c))
    ordered = cells[perm]
    keys = [something(globalindex(levelgrid(sys, level(c)), c), 0) for c in ordered]
    return MultiOrderCellSet{typeof(sys),ID}(sys, ordered, keys, contained[perm],
        reference_level)
end
