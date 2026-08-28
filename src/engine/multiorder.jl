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
time, and a cell whose replacement would not fit is kept whole rather than
ending the search. The one set larger than the budget is a seed that already
was. The depth is then whatever the budget bought, and varies from branch to
branch.

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

`maxcells` is CARDINALITY FIRST: refine the crossing cells coarsest first, and
keep whole any cell whose replacement would not fit rather than ending the
search.

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

Every point of the target lies inside one of the emitted cells. That is the
plain reading of "ten cells that cover California", and it is the guarantee this
mode is built around. The seed is the coarsest cells that meet the target, and
those tile the sphere between them; refinement then replaces a crossing cell
only by the children that meet the target. Where refinement is not congruent it
descends through cells that MISS as well, as `level` mode does — except out of a
cell kept whole, which is never opened — so that a cell meeting the target under
a parent that does not is still found and admitted; a crossing cell none of whose
children meet is then dropped rather than kept — except where the budget refused
an admission, which ends the descent everywhere at once, and then the dropped
cells that no member of the covering lies against are kept after all, as many as
the budget still has room for.

It is EXACT, at every budget, on the three systems whose four children tile
their parent: HEALPix, S2 and ISEA4R. Where children do not tile their parent it
degrades in the way [`MultiOrderCoverage`](@ref)'s warning already describes, and
for the same reason — replacing a cell by its children swaps one footprint for
another. Measured on a state-sized outline as the fraction of the target lying
in no emitted cell: under 2% on IGEO7, 3% on its authalic wrap, 15% on H3 and
30% on A5. Those are not decorative figures — they are the bounds
`test/systems/crosssystem/multiorder_budget.jl` asserts, per system, so this
paragraph cannot quietly stop being true.

The LEAF statement the `level` mode makes — every cell of the reference level
that meets the target is a member or the descendant of one — is a law here on
those same three systems only. The `level` mode earns it everywhere by
descending into cells that miss the target, because a child can overhang its
parent, and carrying that descent all the way to `level`. A budget has no fixed
depth to carry it to: it makes the same descent, but each cell it finds costs a
cell of the budget, and the first one it cannot pay for ends the search, so a
branch stalled for want of cells stops there. On the same outline the leaf
statement misses under 1% of the target on IGEO7, 2% on its authalic wrap and
on H3, and 18% on A5. Both statements are pinned per system, at three budgets
and on four targets, in `test/systems/crosssystem/multiorder_budget.jl`.

What a budget does NOT buy is a tight picture of the target: at ten cells the
set over-covers California by a wide margin, and it says so through
[`iscontained`](@ref) rather than by pretending otherwise.

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

# ---------------------------------------------------------------------------
# The budget traversal
#
# Same target preparation, same prunes and the same two exact predicates as the
# depth-first walk above. What changes is the SCHEDULE and what stops it.
#
# The covering is maintained explicitly as {contained cells} ∪ {crossing
# members}. Contained cells are never refined: a cell proven to lie inside the
# target is already the tightest thing that can be said about its own footprint,
# and splitting it spends budget to say the same thing in more words. Only
# crossing members are candidates, and they are visited COARSEST FIRST, ties
# broken by index within the level — the curve order, so the schedule is a
# function of the inputs alone and of nothing else.
#
# PHANTOMS are the second stream, and they carry no cells into the set. Where
# refinement is not congruent a cell can meet the target while its parent misses
# it, so a descent that follows only cells that meet loses whole lineages of the
# covering. A phantom is a cell that fails the exact test but passes both
# `node_extent` prunes — the descent `_coverage_visit!` makes — and every cell
# this traversal opens contributes its missing children to the stream. A phantom
# costs nothing against the budget and exists only to reach the meeting cells
# beneath it, which are admitted one at a time as the budget allows. Where
# `has_congruent_refinement` holds there are none: children inside their parent
# means a meeting cell always has a meeting parent, and the member stream
# already reaches everything.
#
# A cell KEPT WHOLE is not opened: no phantom descends from it, because a member
# admitted under a member would be emitted as its own ancestor's descendant.
#
# The priority queue is spelled here as two vectors and a pass counter rather
# than as a heap: each pass sorts what it holds on `(level, key)` and refines in
# that order, so the schedule stays total and a function of the inputs alone.
# Every cell of either stream advances exactly one level per pass, which is the
# order a heap on `(level, key)` would pop. Members are refined before phantoms
# within a pass; that order decides only which admission meets the budget first.
#
# WHEN A REPLACEMENT DOES NOT FIT the cell is kept whole and the traversal moves
# on to the next candidate rather than stopping. What that gains is the tail of
# the budget — a coarse crossing cell with seven intersecting children can
# overrun a budget that a finer one with two children still fits into, and
# stopping on the first miss would strand those cells unspent. The queue is
# walked once and never revisited: the covering grows, except where a cell is
# set aside for its overhang and hands its slot back, and a replacement refused
# before that is not offered the slot afterwards.
#
# PENDING CELLS are the only way the set ever shrinks. A member none of whose
# children meets the target is held pending rather than kept or refined: it is
# many times the target's size, and the phantom stream reaches the members that
# cover what it covers. It leaves the covering, and its slot with it, when it
# is set aside; the end phase then drops it for good or hands it back, and is
# the one place the covering can GROW — a freed slot spent on an entrant, then
# another freed slot handing a pending cell back.
#
# WHY THE DROP IS SOUND, at its strongest: `node_extent` bounds a cell's whole
# subtree, so where no child even reaches, nothing under the cell meets the
# target, `level` mode emits nothing there either, and the drop reproduces that
# mode exactly. Where some child does reach, the phantom stream carries the
# descent on and what it admits covers this cell's share of the target. Both
# arguments need a phantom stream, so a cell is pended only on a system where
# one runs: on a congruent system a childless member is kept whole instead. Where the budget
# refused an admission the stream stopped everywhere, including over lineages
# that had asked for nothing, and the end phase hands back the pending cells
# that were left uncovered by it. `_budget_resolve!` has the last word.
# ---------------------------------------------------------------------------

function _multi_order_budget(sys::AbstractHierarchicalGridSystem, target_value,
        maxcells::Int, maxlevel::Int)
    maxcells >= 1 || throw(ArgumentError(
        "maxcells must be at least 1, got $maxcells"))
    maxlevel in levels(sys) || throw(ArgumentError(
        "maxlevel $maxlevel is outside $(typeof(sys))'s levels $(levels(sys))"))
    target = _query_target(target_value)
    ID = cellindextype(sys)
    top = first(levels(sys))
    grids = [levelgrid(sys, l) for l in top:maxlevel]

    # Children inside their parent means the member stream reaches every cell
    # that meets the target, and a phantom could only add work.
    congruent = has_congruent_refinement(sys)

    contained = ID[]        # proven `Within`: never refined, never requeued
    stalled = ID[]          # crossing, and kept whole: the budget said no
    current = ID[]          # crossing members, at the level being refined
    phantoms = ID[]         # non-members that passed the prunes, same level
    pending = ID[]          # members no child of which meets: settled at the end
    entrant_refused = false # a meeting cell the budget could not pay for
    for c in rootcells(sys)
        # One prune per cell, shared: admission needs it and so does the phantom
        # test, and `node_extent` is the expensive half of both.
        _budget_reaches(sys, target, c) || continue
        _budget_admit!(contained, current, sys, target, c, grids, top, maxlevel) && continue
        congruent || push!(phantoms, c)
    end
    # The running size of the covering, maintained rather than recomputed: a
    # replacement is committed exactly when the size it would leave behind fits.
    # Neither phantoms nor pending cells are in it — a pending cell is one the
    # end phase has to be given a reason to keep.
    total = length(contained) + length(current)
    schedule_key(c) = (level(c), _budget_key(grids[level(c)-top+1], c))

    kids_in = ID[]
    kids_out = ID[]
    kids_miss = ID[]
    for _ in top:(maxlevel-1)
        (isempty(current) && isempty(phantoms)) && break
        sort!(current; by=schedule_key)
        sort!(phantoms; by=schedule_key)
        next = ID[]
        next_phantoms = ID[]
        # A cell that misses is worth carrying only while the phantom stream is
        # running: `congruent` switches it off outright, and a refusal has
        # already stopped it, so collecting either way fills a vector that gets
        # thrown away. Only the phantom pass below can flip this, so it is read
        # once for the whole member pass.
        phantoms_live = !(congruent || entrant_refused)
        for c in current
            # Nothing at `maxlevel` has children to reach for. Both streams
            # advance one level per pass, so nothing reaches this today; it is
            # what keeps `grids` in range should a cell ever enter deeper.
            if level(c) >= maxlevel
                push!(stalled, c)
                continue
            end
            empty!(kids_in)
            empty!(kids_out)
            empty!(kids_miss)
            for child in children(sys, c)
                _budget_reaches(sys, target, child) || continue
                _budget_admit!(kids_in, kids_out, sys, target, child, grids,
                    top, maxlevel) && continue
                phantoms_live && push!(kids_miss, child)
            end
            k = length(kids_in) + length(kids_out)
            if k == 0 && !congruent
                # No child meets: the target sits in this cell's overhang
                # annulus (non-congruent refinement). Refining it would drop the
                # target, keeping it buys a cell many times the target's size,
                # so it leaves the covering here and the end phase confirms it.
                push!(pending, c)
                total -= 1
                append!(next_phantoms, kids_miss)
            elseif k == 0 || total + k - 1 > maxcells
                # Kept whole, and with it everything under it: opening the cell
                # would put members below a member.
                #
                # `k == 0` lands here on a CONGRUENT system, where the drop has
                # no argument behind it: children that tile their parent cannot
                # lose the target between them, so a childless member is a
                # boundary-polyline sliver rather than an overhang, and no
                # phantom stream runs to find what covers it. Dropping it would
                # delete a piece of the covering against this mode's one
                # unqualified promise, that it is exact where refinement is.
                push!(stalled, c)
            else
                append!(contained, kids_in)
                append!(next, kids_out)
                total += k - 1
                append!(next_phantoms, kids_miss)
            end
        end
        for p in phantoms
            entrant_refused && break
            # Dead defence, like its twin in the member loop above: both streams
            # advance one level per pass, so nothing reaches `maxlevel` here
            # today. It is what keeps `grids` in range should one ever enter
            # deeper.
            level(p) >= maxlevel && continue
            for child in children(sys, p)
                _budget_reaches(sys, target, child) || continue
                empty!(kids_in)
                empty!(kids_out)
                if !_budget_admit!(kids_in, kids_out, sys, target, child, grids,
                        top, maxlevel)
                    push!(next_phantoms, child)
                    continue
                end
                # A meeting cell whose parent misses replaces nothing, so it
                # costs a whole cell rather than the difference.
                if total + 1 > maxcells
                    entrant_refused = true
                    break
                end
                append!(contained, kids_in)
                append!(next, kids_out)
                total += 1
            end
        end
        current = next
        # The first refusal ends the search: the budget only shrinks from here,
        # so no later phantom could be paid for either.
        phantoms = entrant_refused ? empty!(next_phantoms) : next_phantoms
    end
    append!(stalled, current)       # members only: what `maxlevel` left unrefined
    _budget_resolve!(stalled, contained, sys, pending, grids, top, maxcells,
        entrant_refused)

    cells = vcat(contained, stalled)
    flags = falses(length(cells))
    flags[1:length(contained)] .= true
    # No depth was asked for, so the set speaks about the deepest level it
    # reached. An empty set reached nothing and reports the top.
    reference = isempty(cells) ? top : maximum(level, cells)
    return _sorted_cell_set(sys, cells, flags, reference)
end

# Can anything in this cell's SUBTREE belong to the covering? Both prunes read
# `node_extent`, which bounds the subtree, so a cell that fails them shelters
# nothing that meets the target — which is what makes a phantom worth carrying
# and a pruned cell worth forgetting. `_coverage_visit!` descends on exactly
# this test.
function _budget_reaches(sys, target, c)
    extent = node_extent(sys, c)
    Extents.intersects(target.cap, extent) || return false
    return !_subtree_outside(target, extent)
end

# Classify one candidate cell, or reject it: the `Intersects` question, then the
# `Within` one. Both are the predicates `_coverage_visit!` uses, in the same
# order and for the same reasons — including the `maxlevel` arm, where `Within`
# is not asked because no decision depends on the answer and the call is the
# expensive one. `iscontained` documents what that leaves unproven.
#
# CALLERS RUN `_budget_reaches` FIRST. A rejected cell is a phantom candidate,
# and deciding that needs the same `node_extent` prune admission would run, so
# the prune sits at the call site and is paid once instead of twice.
function _budget_admit!(contained, crossing, sys, target, c, grids, top::Int,
        maxlevel::Int)
    lc = level(c)
    grid = grids[lc-top+1]
    _matches(DE9IM.Intersects(nothing), target, grid, c) || return false
    if lc < maxlevel && _matches(DE9IM.Within(nothing), target, grid, c)
        push!(contained, c)
    else
        push!(crossing, c)
    end
    return true
end

# The last word on the pending cells, and the two cases where one comes back. A
# pending cell is a member none of whose children meets the target, so it is
# many times the target's size and the phantom stream reaches the members that
# cover what it covers. It left the covering the moment it was set aside — and
# the budget went with it, which is why `room` below is what the covering has
# actually spent and not what it was allowed.
#
# THE SEARCH STOPPED SHORT where an entrant was refused: from there on the
# stream was buying nothing, and a part of the target can be left with no member
# at all — an offshore island whose one replacement cell the budget never
# reached. So a pending cell the rest of the covering says nothing about is kept
# whole. Two tests decide "says nothing about", and only the first is a proof:
#
#   * no member is a DESCENDANT of it. An emitted cell under an emitted cell is
#     a malformed set, not a loose covering, so this is exact and structural.
#   * no member's CELL CAP touches its own. That is a proximity test on bounding
#     caps — necessary for a member to cover any of this cell, nowhere near
#     sufficient — so the rule is a heuristic and the union bounds in
#     `test/systems/crosssystem/multiorder_budget.jl` are what gate the
#     approximation. It reads caps rather than node extents because the question
#     is whether a member's own footprint lies near this one, not whether their
#     subtrees might meet.
#
# Deepest first, ties by index within the level, so that where two of them
# overlap the smaller is the one kept and the order is total.
#
# NOTHING IS KEPT ONCE THE BUDGET IS SPENT. A refusal happens at
# `total == maxcells` exactly, so the only slots this rule has to spend are the
# ones the pending cells themselves handed back, and `room` counts them. The
# rule that hands a cell back must not be the one thing that breaks the bound,
# which is what `room` is for and why it decrements per keep.
#
# A COVERING WITH NOTHING IN IT is the other case. A budget of one, on a target
# smaller than a cell, can spend itself before the descent finds anything to
# keep, and the answer to a target the system meets is a cell, not an empty set.
# The SHALLOWEST pending cell is the one kept: every pending cell meets the
# target, but only its own footprint is offered as the covering, and the coarser
# the cell the more of the target that footprint can hold. It is reachable only
# with `contained` and `stalled` both empty, so `room` is the whole budget there
# and one cell is within it.
function _budget_resolve!(stalled, contained, sys, pending, grids, top::Int,
        maxcells::Int, entrant_refused::Bool)
    isempty(pending) && return nothing
    room = maxcells - (length(contained) + length(stalled))
    if entrant_refused && room > 0
        cellcap(c) = cell_cap(grids[level(c)-top+1], c)
        caps = [cellcap(m) for m in contained]
        append!(caps, cellcap(m) for m in stalled)
        sort!(pending; by=c -> (-level(c), _budget_key(grids[level(c)-top+1], c)))
        for c in pending
            room > 0 || break
            (_ancestor_of_any(sys, c, contained) ||
             _ancestor_of_any(sys, c, stalled)) && continue
            cap = cellcap(c)
            any(k -> Extents.intersects(cap, k), caps) && continue
            push!(stalled, c)
            push!(caps, cap)
            room -= 1
        end
    end
    (isempty(contained) && isempty(stalled)) && push!(stalled, argmin(level, pending))
    return nothing
end

# Is `c` an ancestor of anything in `members`? Emitting both would put a cell
# under its own ancestor, which every consumer of a set — `level_ranges` most
# visibly — reads as two disjoint claims about the same leaves.
function _ancestor_of_any(sys, c, members)
    lc = level(c)
    return any(m -> level(m) > lc && ancestor(sys, m, lc) == c, members)
end

# The tie-break inside one level. `globalindex` is the level's own order, which
# is curve order on every system here and is a bijection, so no two cells of a
# level ever tie and the schedule has nothing left to decide. A missing index
# would break that bijection and silently alias two cells onto one key, which is
# a determinism bug wearing a plausible answer — so it is an error, not a zero.
function _budget_key(grid, c)::Int
    pos = globalindex(grid, c)
    pos === nothing && throw(ArgumentError(
        "$(typeof(grid)) has no index for the cell $c it just produced. The " *
        "budget schedule orders each level by `globalindex` and needs it total: " *
        "without it two cells share a key and the traversal stops being " *
        "deterministic"))
    return pos
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

# --- the collection surface ------------------------------------------------

Base.length(set::MultiOrderCellSet) = length(set.cells)
Base.eltype(::Type{MultiOrderCellSet{S,ID}}) where {S,ID} = ID
Base.eltype(set::MultiOrderCellSet) = eltype(typeof(set))
Base.isempty(set::MultiOrderCellSet) = isempty(set.cells)
Base.iterate(set::MultiOrderCellSet, state...) = iterate(set.cells, state...)
Base.getindex(set::MultiOrderCellSet, i::Int) = set.cells[i]
Base.firstindex(::MultiOrderCellSet) = 1
Base.lastindex(set::MultiOrderCellSet) = length(set.cells)
Base.eachindex(set::MultiOrderCellSet) = Base.OneTo(length(set.cells))
Base.collect(set::MultiOrderCellSet) = copy(set.cells)

"""
    system(set::MultiOrderCellSet)

The system the set's cells are named in.
"""
system(set::MultiOrderCellSet) = set.system

# --- which cells fit inside the target -------------------------------------

"""
    iscontained(set::MultiOrderCellSet, i::Integer) -> Bool

Whether the set's `i`th cell was **proven** to lie inside the coverage target.
`true` means `Within` was asked of it and held. `false` means one of two things,
told apart by the cell's level:

  - above the traversal's maximum depth, the cell was asked and the target's
    boundary crosses it. The flag is exact there, in both directions.
  - at the maximum depth, the cell was never asked. The traversal ran out of
    depth and emitted it to cover; it may or may not fit.

In `level` mode the maximum depth *is* the reference level, so the blind spot
sits exactly there. In `maxcells` mode they part: the reference level is the
deepest level the budget reached, the maximum depth is the `maxlevel` cap, and a
budget that stopped short of the cap — the ordinary case — carries an exact flag
on every member.

The asymmetry is the contract. `Within` allocates ~48 KB per call against ~600 B
for `Intersects`, and a deep coverage ends on thousands of reference-level
cells: labelling them exactly costs hundreds of megabytes and changes no member
of the set. The direction that matters survives — `true` implies inside — and
code needing exact containment at the reference level asks `Within` of the few
cells it cares about.

`argmin(level, set)` is therefore not "the coarsest cell inside the target":
every emission can be unproven. [`coarsest_contained`](@ref) reads this flag.
"""
iscontained(set::MultiOrderCellSet, i::Integer) = set.contained[i]

"""
    coarsest_contained(set::MultiOrderCellSet) -> cell id or `nothing`

The shallowest cell of `set` **proven** inside the coverage target, or `nothing`
when no cell of it was — see [`iscontained`](@ref) for what "proven" leaves
out. Maximum-depth cells are never tested, so a set of nothing but those answers
`nothing` even where some fit; a target smaller than one cell is the clearest
way there, not the only one.

A budget set answers `nothing` for a second, more ordinary reason: at ten cells
over a state nothing has been refined far enough to fit inside it, and the
accessor says so rather than hand back the shallowest crossing cell. Raise the
budget and the answer appears.

```julia
set = query(sys, MultiOrderCoverage(tile); level = 12)
cell = coarsest_contained(set)          # `nothing`, or a cell that fits in `tile`
```

Ties go to the first such cell in the set's own order.
"""
function coarsest_contained(set::MultiOrderCellSet)
    best = nothing
    for i in eachindex(set)
        set.contained[i] || continue
        (best === nothing || level(set.cells[i]) < level(best)) && (best = set.cells[i])
    end
    return best
end

# --- geometry, without a level grid per cell -------------------------------

# A set is not a grid — no indices, and its cells are at different levels —
# but it does know which level grid each cell belongs to. `levelgrid` is O(1),
# so nothing here is worth caching.

"""
    cell_boundary(set::MultiOrderCellSet, c) -> Vector{UnitSphericalPoint}
    cell_centroid(set::MultiOrderCellSet, c) -> UnitSphericalPoint
    cell_polygon(set::MultiOrderCellSet, c) -> GI.Polygon

Geometry of one member of a mixed-level set, read from
`levelgrid(system(set), level(c))`. Same values as that grid gives, without the
caller resolving it.
"""
cell_boundary(set::MultiOrderCellSet, c::AbstractCellIndex) =
    cell_boundary(levelgrid(set.system, level(c)), c)

cell_centroid(set::MultiOrderCellSet, c::AbstractCellIndex) =
    cell_centroid(levelgrid(set.system, level(c)), c)

cell_polygon(set::MultiOrderCellSet, c::AbstractCellIndex) =
    cell_polygon(levelgrid(set.system, level(c)), c)

"""
    cell_polygons(set::MultiOrderCellSet) -> Vector{<:GI.Polygon}

Every cell of the set as a unit-sphere polygon, in the set's own order: what a
plot of a coverage needs, in one call.

```julia
poly(GO.transform(GO.GeographicFromUnitSphere(), cell_polygons(set)))
```
"""
cell_polygons(set::MultiOrderCellSet) =
    [cell_polygon(levelgrid(set.system, level(c)), c) for c in set.cells]

"""
    curve_keys(set::MultiOrderCellSet) -> Vector{Int}

Return stored-cell sort keys. For sorted-subtree systems, each key is the start
of the cell's reference-level descendant range. Otherwise it is the cell's
index within its own level and is not comparable across levels.
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

Expand the set to sorted, disjoint index ranges in `levelgrid(sys, l)`,
merging adjacent ranges. Requires sorted subtrees and `l` no shallower than any
cell in the set.

!!! warning "Two things the expansion is not"
    Not universal: it throws where [`has_sorted_subtrees`](@ref) is `false`
    (A5). Branch on the trait, or expand with `descendants(sys, c, l)`.

    Not a covering. A cell is in the set because the *cell* is inside the
    target; under non-congruent refinement its descendants need not be, so the
    expansion can name leaves the target does not touch — most visibly inside a
    hole. See [`MultiOrderCoverage`](@ref).
"""
function level_ranges(set::MultiOrderCellSet, l::Integer)
    has_sorted_subtrees(set.system) || throw(ArgumentError(
        "$(typeof(set.system)) has no descendant ranges, so a multi-order set " *
        "cannot be expanded to index ranges"))
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
wherever the indices are what is wanted.
"""
function cellindices(set::MultiOrderCellSet, l::Integer)
    grid = levelgrid(set.system, Int(l))
    out = cellindextype(set.system)[]
    for r in level_ranges(set, l), i in r
        push!(out, cellindex(grid, i))
    end
    return out
end
