# Multi-order queries return the coarsest cells covering a region, ordered by
# descendant-range start at a reference level — the generalisation of HEALPix's
# MOC/NUNIQ.

"""
    MultiOrderCoverage(target)

Request cells at several levels through [`query`](@ref). `target` accepts a
GeoInterface geometry, a longitude/latitude extent in degrees, or a spherical cap.
Supply exactly one mode to `query`:

- `level=l`: refine boundary crossings to level `l`. Expanding the result at
  that level includes every intersecting cell. With congruent refinement, the
  expansion equals the fixed-level intersection result.
- `maxcells=n`: refine within a cell budget, optionally limited by `maxlevel`.
  If the initial seed exceeds `n`, it is returned over budget. On noncongruent
  systems, this mode does not guarantee full polygon or reference-level coverage.

HEALPix, S2, and ISEA4R have congruent refinement. IGeo7, H3, and A5 do not:
a parent's polygon can differ from the union of its descendants.
Use [`iscontained`](@ref) for proven containment and [`level_ranges`](@ref)
for represented leaves where sorted subtrees are available.
See [Multi-order coverage](@ref) for the mode comparison and examples.
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

The REFERENCE LEVEL is the depth the set speaks about — the `level` the query
was given, or in `maxcells` mode the deepest level the budget reached. It is
the default expansion level for [`level_ranges`](@ref), [`cellindices`](@ref)
and `CellLookup`. Coverage at that level depends on the query mode and
refinement traits; see [`MultiOrderCoverage`](@ref).

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
`level`, or refining within `maxcells` (subject to the seed exception). Equivalent to
`query(sys, coverage; ...)`, which documents both modes.
"""
MultiOrderCellSet(sys::AbstractHierarchicalGridSystem, coverage::MultiOrderCoverage;
    level::Union{Integer,Nothing}=nothing, maxcells::Union{Integer,Nothing}=nothing,
    maxlevel::Union{Integer,Nothing}=nothing) =
    _multi_order_query(sys, coverage.target, level, maxcells, maxlevel)

"""
    query(sys, coverage::MultiOrderCoverage; level) -> MultiOrderCellSet
    query(sys, coverage::MultiOrderCoverage; maxcells, maxlevel = deepest) -> MultiOrderCellSet

Return a mixed-level representation of the target. Supply exactly one of
`level` and `maxcells`; invalid combinations raise `ArgumentError`.

`level` fixes the finest level. Its expanded membership includes every
intersecting cell at that level and equals that set on congruent systems.
The emitted polygons need not cover the same region as their descendants.

`maxcells` is a refinement budget, not an unconditional maximum. The initial
coarsest-cell seed is returned whole if it already exceeds the budget.
Otherwise `length(result) <= maxcells`. `maxcells < 1` raises `ArgumentError`.
Set `maxlevel` with this mode to limit depth, especially for small targets.

On HEALPix, S2, and ISEA4R, the budget result covers the target and represents
all intersecting leaves at its reference level. IGeo7, H3, and A5 have no such
budget-mode guarantee: both polygon and leaf coverage can miss part of the target.
Increasing a budget does not turn empirical error measurements into a contract.

The reference level is the requested `level`, or the deepest level reached in
budget mode. It is the default for expansion and `CellLookup` construction.
[`level_ranges`](@ref) requires sorted subtrees; A5 does not support that operation.
[`cell_polygons`](@ref) displays chosen members, and [`iscontained`](@ref)
reports cells proven to fit inside the target.

See [Multi-order coverage](@ref) for a worked comparison of the two modes.
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
    # One grid per level, built once: `levelgrid` is cheap but not free.
    top = first(levels(sys))
    grids = [levelgrid(sys, l) for l in top:maxlevel]
    for c in rootcells(sys)
        _coverage_visit!(cells, contained, sys, target, c, maxlevel, grids, top)
    end
    return _sorted_cell_set(sys, cells, contained, maxlevel)
end

function _coverage_visit!(cells, contained, sys, target, c, maxlevel::Int, grids, top::Int)
    # Only `node_extent` may prune descent — child geometry can overhang its
    # parent; exact cell geometry decides emission only. The one-distance cap
    # test runs before the boundary-arc proof, keeping the traversal
    # output-sensitive on a loose or whole-sphere cap.
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
    # `(level, id)` is the documented fallback; the keys are the cells' own
    # within-level indices — a total order, not a curve order.
    perm = sortperm(cells; by=c -> (level(c), c))
    ordered = cells[perm]
    keys = [something(globalindex(levelgrid(sys, level(c)), c), 0) for c in ordered]
    return MultiOrderCellSet{typeof(sys),ID}(sys, ordered, keys, contained[perm],
        reference_level)
end
