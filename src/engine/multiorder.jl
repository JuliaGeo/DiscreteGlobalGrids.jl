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

Every point of the target lies inside one of the emitted cells. That is the
plain reading of "ten cells that cover California", and it is the guarantee this
mode is built around. The seed is the coarsest cells that meet the target, and
those tile the sphere between them; refinement then replaces a crossing cell by
the children that meet the target. Non-congruent refinement gets the same
descent `level` mode has — through cells that miss the target as well — so a
cell meeting the target under a parent that misses it is still found, and a
crossing cell none of whose children meet is dropped in favour of the cells
that do cover its share. When the budget cuts that descent short, the dropped
cells' reserved slots pay for what the walk still owes: the cells it found
and could not afford, or a dropped cell kept whole after all where its share
was never certified.

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
depth to carry it to: it makes the same descent, and where it stops paying, the
search goes on only for the dropped cells whose share is still unproven — the
rest of the frontier ends there. On the same outline the leaf
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
# PHANTOMS carry the descent through cells that miss the target. Under
# non-congruent refinement a cell can meet the target while its parent misses
# it; a phantom is a cell that fails the exact test but passes both
# `node_extent` prunes — the descent `_coverage_visit!` makes — and every
# opened cell contributes its missing children to the stream. Phantoms are
# free: only the meeting cells found beneath them are charged, one whole cell
# each, as ENTRANTS. On a congruent system every meeting cell has a meeting
# parent, so the member stream reaches everything and no phantoms run. A cell
# kept whole is never opened: a member admitted beneath it would be emitted
# under its own ancestor.
#
# The queue is two vectors and a pass counter: each pass sorts both streams on
# `(level, key)` and advances every cell one level — the order a heap on that
# key would pop. Members refine before phantoms within a pass.
#
# PENDING cells are members none of whose children meet: the target lies in
# their overhang annulus, and the phantom stream reaches the cells that cover
# it. A pended cell leaves the covering but keeps its budget slot as a CLAIM
# until one of three settlements takes it:
#   * REDEEMED — at the budget's margin, by an entrant that provably stands
#     in for it (see `_claim_match`);
#   * COMPLETED — its neighbourhood searched to the end, every meeting cell
#     admitted: the share is covered, the slot released to the stream;
#   * RESOLVED — by `_budget_resolve!` at the end: dropped where the search
#     covered its share, handed back on its own slot where the evidence
#     stands on its cell.
# Only non-congruent systems pend.
#
# A REPLACEMENT that does not fit leaves its cell kept whole and the walk
# moves on: a coarse cell with seven meeting children can overrun a budget
# that a finer one with two still fits into. An ENTRANT costs a fresh slot,
# or — at the margin — a claim it provably stands in for; one that can pay
# neither is REFUSED and recorded: a meeting cell the covering still owes,
# which `_budget_resolve!` pays for out of the slots the claims held. After a
# refusal, every branch no unsettled claim keeps alive dies where it stands,
# recorded as ABANDONED — territory the walk never revisits, which is what
# lets later refusals stay local and settled claims release soundly. Refused
# and abandoned cells together bound everything the walk left uncovered:
# they are the EVIDENCE the end phase settles claims against.
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

    # Congruent systems need no phantoms: every meeting cell has a meeting parent.
    congruent = has_congruent_refinement(sys)

    contained = ID[]        # proven `Within`: never refined, never requeued
    stalled = ID[]          # crossing, and kept whole: the budget said no
    current = ID[]          # crossing members, at the level being refined
    phantoms = ID[]         # non-members that passed the prunes, same level
    pending = ID[]          # members no child of which meets: settled at the end
    settled = Bool[]        # per pending cell: no longer holds a budget slot
    covered = Bool[]        # per pending cell: share proven covered
    refused = ID[]          # meeting cells the budget could not pay for
    abandoned = ID[]        # branches a refusal killed: unsearched territory
    entrant_refused = false
    for c in rootcells(sys)
        # The prune is shared: `_budget_admit!` assumes it ran, and rejection
        # feeds the phantom stream.
        _budget_reaches(sys, target, c) || continue
        _budget_admit!(contained, current, sys, target, c, grids, top, maxlevel) && continue
        congruent || push!(phantoms, c)
    end
    # The budget's running claim: contained + crossing + stalled + unsettled
    # pending. A pended cell holds its claim until it is settled — redeemed,
    # released on completion, or resolved at the end — so the end phase always
    # has room for what is still uncovered. Phantoms are outside it.
    total = length(contained) + length(current)
    schedule_key(c) = (level(c), _budget_key(grids[level(c)-top+1], c))
    cellcap(c) = cell_cap(grids[level(c)-top+1], c)
    # What a redeeming entrant must contain whole (see `_claim_match`).
    pieces = congruent ? nothing : _target_pieces(target_value, target)

    kids_in = ID[]
    kids_out = ID[]
    kids_miss = ID[]
    for _ in top:(maxlevel-1)
        (isempty(current) && isempty(phantoms)) && break
        sort!(current; by=schedule_key)
        sort!(phantoms; by=schedule_key)
        next = ID[]
        next_phantoms = ID[]
        for c in current
            # Unreachable while both streams enter at the roots and advance one
            # level per pass; keeps `grids` in range if entry ever deepens.
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
                congruent || push!(kids_miss, child)
            end
            k = length(kids_in) + length(kids_out)
            if k == 0 && !congruent
                # No child meets: the target lies in this cell's overhang
                # annulus, and the phantom stream reaches the cells that cover
                # it. The claim is held only where handing the cell back could
                # ever be the right answer — commensurate with a piece it
                # guards; a cell that dwarfs every piece it touches frees its
                # slot to the stream, and the refused entrants stand in for
                # its share at the end.
                push!(pending, c)
                held = _claim_worth_holding(cellcap(c),
                    ncells(grids[level(c)-top+1]), pieces)
                push!(settled, !held)
                push!(covered, false)
                held || (total -= 1)
                append!(next_phantoms, kids_miss)
            elseif k == 0 || total + k - 1 > maxcells
                # Kept whole; the subtree stays closed so no member lands under
                # a member. Congruent `k == 0` lands here too: children that
                # tile their parent cannot lose the target between them — the
                # cell meets it on a boundary sliver — and dropping it would
                # break the exactness this mode promises where refinement is.
                push!(stalled, c)
            else
                append!(contained, kids_in)
                append!(next, kids_out)
                total += k - 1
                append!(next_phantoms, kids_miss)
            end
        end
        # The unsettled claims, deepest first: what marginal entrants redeem,
        # and — after a refusal — the only territory still worth searching.
        live = [i for i in eachindex(pending) if !covered[i]]
        sort!(live; by=i -> (-level(pending[i]), schedule_key(pending[i])[2]))
        livecaps = [cellcap(pending[i]) for i in live]
        for p in phantoms
            # After a refusal, a branch matters only where a claim's share is
            # still unproven; the rest die where they stand, recorded.
            if entrant_refused
                ext = node_extent(sys, p)
                if !any(k -> !covered[live[k]] &&
                            Extents.intersects(livecaps[k], ext), eachindex(live))
                    push!(abandoned, p)
                    continue
                end
            end
            # Same `grids`-range guard as the member loop.
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
                # costs a whole cell: a fresh slot, or — at the margin — the
                # slot of a claim it provably stands in for.
                if total + 1 <= maxcells
                    total += 1
                elseif (m = _claim_match(grids[level(child)-top+1], child,
                        pieces, live, livecaps, settled, total, maxcells)) != 0
                    settled[m] = true
                    covered[m] = true
                else
                    entrant_refused = true
                    push!(refused, child)
                    continue
                end
                append!(contained, kids_in)
                append!(next, kids_out)
            end
        end
        # A claim's share is proven covered by COMPLETION when nothing can
        # change its covering any more: no live phantom still reaches it, and
        # no refusal was recorded against it — nor can one be later, since
        # any branch able to reach it would still be live now. Everything
        # found there was admitted, so the share is covered; a slot still
        # held is surplus, released for the stream to spend.
        if !congruent && !isempty(pending)
            exts = [node_extent(sys, p) for p in next_phantoms]
            aexts = [node_extent(sys, a)
                     for a in Iterators.flatten((abandoned, refused))]
            for i in eachindex(pending)
                covered[i] && continue
                cap = cellcap(pending[i])
                any(x -> Extents.intersects(cap, x), exts) && continue
                any(x -> Extents.intersects(cap, x), aexts) && continue
                covered[i] = true
                if !settled[i]
                    settled[i] = true
                    total -= 1
                end
            end
        end
        current = next
        phantoms = next_phantoms
    end
    append!(stalled, current)       # members only: what `maxlevel` left unrefined
    _budget_resolve!(stalled, contained, sys, pending, settled, covered,
        grids, top, maxcells, abandoned, refused, pieces)

    cells = vcat(contained, stalled)
    flags = falses(length(cells))
    flags[1:length(contained)] .= true
    # The reference level is the deepest level reached; an empty set reports the top.
    reference = isempty(cells) ? top : maximum(level, cells)
    return _sorted_cell_set(sys, cells, flags, reference)
end

# `node_extent` bounds the whole subtree, so a cell failing these prunes
# shelters nothing that meets the target. The same test `_coverage_visit!`
# descends on.
function _budget_reaches(sys, target, c)
    extent = node_extent(sys, c)
    Extents.intersects(target.cap, extent) || return false
    return !_subtree_outside(target, extent)
end

# Classify one candidate: `Intersects`, then `Within` — the predicates
# `_coverage_visit!` uses, in the same order; at `maxlevel`, `Within` is skipped
# (see `iscontained`). Callers run `_budget_reaches` first: a rejected cell is a
# phantom candidate, and both decisions share that one prune.
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

# One (bounding cap, target) pair per connected piece of the target — a
# MultiPolygon's parts, each prepared on its own; any other target is one
# piece, itself.
function _target_pieces(target_value, target)
    (GI.isgeometry(target_value) &&
     GI.trait(target_value) isa GI.MultiPolygonTrait) ||
        return [(target.cap, target)]
    return [(points_cap([query_point(p) for p in GI.getpoint(g)]),
             _query_target(g)) for g in GI.getpolygon(target_value)]
end

# A held claim promises the cell may come back whole at the end. That is only
# ever the right answer for a cell commensurate with a piece it touches: a
# giant handed back over a speck is the annulus pathology, so its claim is
# worthless as insurance and the slot serves the stream better. The yardstick
# is the piece cap's spherical area against the cell's: the cell must fit the
# cap, measured on the annulus and refusal-isolation laws — a cap already
# overshoots its piece, so 1 carries the slack a covering cell needs.
const _CLAIM_PROPORTION = 1.0

_cap_fraction(cap) = (1 - cos(min(Float64(cap.radius), Float64(pi)))) / 2

function _claim_worth_holding(cellcap, ncells_level, pieces)
    yardstick = 0.0
    for (cap, _) in pieces
        Extents.intersects(cellcap, cap) || continue
        yardstick = max(yardstick, _cap_fraction(cap))
    end
    return 1 / ncells_level <= _CLAIM_PROPORTION * yardstick
end

# Whether one whole piece of the target lies inside the cell — the only
# containment the engine answers exactly. A piece kind without an exact test
# certifies nothing.
_piece_in_cell(piece::GeometryTarget, grid, c) =
    _matches(DE9IM.Contains(nothing), piece, grid, c)
_piece_in_cell(piece, grid, c) = false

# The deepest unsettled claim the entrant provably stands in for, or 0: the
# entrant's CELL must contain every piece of the target the claim's cap
# touches — a claim's share is inside its pieces, so the entrant covers it.
# A cap test is not enough: a cap swallows more than its cell does, and a
# claim redeemed on a sliver loses the rest of its share. `live` comes
# deepest first, so the first hit is the smallest cell. Redemption is a
# transfer, not a purchase, and is refused while the covering itself is over
# budget (a seed larger than `maxcells`).
function _claim_match(grid, child, pieces, live, livecaps, settled, total::Int,
        maxcells::Int)
    total <= maxcells || return 0
    holds = zeros(UInt8, length(pieces))    # 0 unknown, 1 contains, 2 does not
    for k in eachindex(live)
        settled[live[k]] && continue
        ok = false
        for (j, (cap, piece)) in enumerate(pieces)
            Extents.intersects(livecaps[k], cap) || continue
            if holds[j] == 0x00
                holds[j] = _piece_in_cell(piece, grid, child) ? 0x01 : 0x02
            end
            ok = holds[j] == 0x01
            ok || break
        end
        ok && return live[k]
    end
    return 0
end

# Hands back what the walk still owes, finest cells first, inside `room`.
#
# Everything the covering leaves uncovered lies inside a refused entrant or
# under a culled branch: those are the evidence. The candidates are the
# refused cells themselves — meeting cells the walk already proved belong in
# the covering, each its own evidence — and the pending cells not proven
# covered, which stand for evidence ON THEIR CELL: not their cap — a coarse
# cell's cap reaches territory its footprint never held, and a giant handed
# back on far-away evidence is the annulus pathology.
#
# One pass, deepest first: emitting a refused cell settles its evidence, so a
# coarser pending cell comes back only for evidence nothing finer explained,
# and of two nested candidates the smaller comes back and blocks the coarser
# (emitting both a cell and its ancestor would claim the same leaves twice).
#
# STARVED pieces come before everything: a piece of the target no member
# stands on has only this pass to be covered at all, and room spent first on
# refused slivers of a well-covered piece is the multipart pathology — a
# budget exhausted on the mainland losing the island whole.
#
# A covering is never empty: with nothing else kept, the shallowest unproven
# cell — the one whose footprint holds the most target — is the answer.
function _budget_resolve!(stalled, contained, sys, pending, settled, covered,
        grids, top::Int, maxcells::Int, abandoned, refused, pieces)
    unresolved = [pending[i] for i in eachindex(pending) if !covered[i]]
    (isempty(unresolved) && isempty(refused)) && return nothing
    room = maxcells - (length(contained) + length(stalled))
    ev = vcat(refused, abandoned)
    evcaps = [cell_cap(grids[level(a)-top+1], a) for a in ev]
    members = Iterators.flatten((contained, stalled))
    evlive = [j > length(refused) ||
              !_refused_covered(sys, ev[j], evcaps[j], members, grids, top)
              for j in eachindex(ev)]
    if any(evlive) && room > 0
        # Held claims only: a freed cell was freed because handing it back
        # could never be proportionate (see `_claim_worth_holding`).
        cands = [(pending[i], 0) for i in eachindex(pending) if !settled[i]]
        append!(cands, [(refused[j], j) for j in eachindex(refused) if evlive[j]])
        # A NEGLECTED piece — far fewer members than the dominant one — is a
        # component the budget passed over, and losing it whole is the
        # multipart pathology; its candidates come first, coarsest first,
        # covering the most with the least. Everything else is tightening,
        # finest first.
        mcount = zeros(Int, length(pieces))
        for m in members
            for j in eachindex(pieces)
                _cell_meets_cap(pieces[j][1], grids[level(m)-top+1], m, false) &&
                    (mcount[j] += 1)
            end
        end
        mmax = maximum(mcount; init=0)
        rescues(c) = any(j -> 2 * mcount[j] < mmax && Extents.intersects(
                cell_cap(grids[level(c)-top+1], c), pieces[j][1]),
            eachindex(pieces))
        sort!(cands; by=((c, _),) -> begin
            r = rescues(c)
            (!r, r ? Int(level(c)) : -Int(level(c)),
             _budget_key(grids[level(c)-top+1], c))
        end)
        for (c, j) in cands
            room > 0 || break
            (_ancestor_of_any(sys, c, contained) ||
             _ancestor_of_any(sys, c, stalled) ||
             _descends_from_any(sys, c, contained) ||
             _descends_from_any(sys, c, stalled)) && continue
            if j == 0
                grid = grids[level(c)-top+1]
                any(k -> evlive[k] && _cell_meets_cap(evcaps[k], grid, c, false),
                    eachindex(ev)) || continue
            end
            push!(stalled, c)
            j == 0 || (evlive[j] = false)
            room -= 1
        end
    end
    (isempty(contained) && isempty(stalled) && !isempty(unresolved)) &&
        push!(stalled, argmin(level, unresolved))
    return nothing
end

# Emitting both a cell and its ancestor claims the same leaves twice; consumers
# read a set as disjoint subtrees.
function _ancestor_of_any(sys, c, members)
    lc = level(c)
    return any(m -> level(m) > lc && ancestor(sys, m, lc) == c, members)
end

# Whether the cap lies entirely inside the cell: centre in the cell, and every
# boundary arc at least the radius away. The converse of `_cell_meets_cap`,
# with the same exactness for convex caps.
function _cell_contains_cap(cap, grid, c)
    ring, n = open_ring(cell_boundary(grid, c))
    point_in_cell(ring, cap.point) === true || return false
    threshold = cos(min(Float64(pi), Float64(cap.radius)))
    for i in 1:n
        a = ring[i]
        b = ring[i == n ? 1 : i+1]
        arc = BoundaryArc(USPoint(a[1], a[2], a[3]), USPoint(b[1], b[2], b[3]))
        _arc_cos_distance(arc, cap.point) > threshold && return false
    end
    return true
end

# A refused cell some single member's footprint holds whole is already
# covered: not a candidate, and not evidence either.
function _refused_covered(sys, e, ecap, members, grids, top)
    for m in members
        mcap = cell_cap(grids[level(m)-top+1], m)
        US.spherical_distance(mcap.point, ecap.point) <= mcap.radius || continue
        _cell_contains_cap(ecap, grids[level(m)-top+1], m) && return true
    end
    return false
end

function _descends_from_any(sys, c, members)
    lc = level(c)
    return any(m -> level(m) < lc && ancestor(sys, c, level(m)) == m, members)
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
