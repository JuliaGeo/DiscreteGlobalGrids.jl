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

A multi-order coverage query: hand it to [`query`](@ref) with a system and
either a maximum depth or a cell budget to get a [`MultiOrderCellSet`](@ref).

```julia
set = query(sys, MultiOrderCoverage(polygon); level = 8)      # accuracy first
set = query(sys, MultiOrderCoverage(polygon); maxcells = 10)  # cardinality first
```

`target` takes the same forms as any other query target — a GeoInterface
geometry, an `Extents.Extent` in lon/lat degrees, or a
`GO.UnitSpherical.SphericalCap`. `Base.parent` unwraps it, as it does for a
DE9IM predicate.

# The two modes

`level` is the ACCURACY-FIRST mode, and the older one. The traversal is depth
first: it emits a cell as soon as the cell lies entirely inside the target, and
recurses into the children of a cell the target's boundary crosses, down to the
requested level; cells still crossed at that level are emitted too, so the set
**covers** the target rather than being covered by it. How many cells that takes
is whatever the outline needs — a coastline at a fine level is tens of thousands
of them.

`maxcells` is the CARDINALITY-FIRST mode. It answers "give me ten cells that
cover California, or a hundred", and it never returns more than the budget
(with one documented exception, below). Refinement is breadth first over the
cells the boundary crosses, coarsest level first, and it stops when the next
replacement would not fit. The depth is then whatever the budget bought, and it
differs from branch to branch: the two modes trade the same currency in
opposite directions.

Neither mode is the other's approximation. A `level` set is the exact answer at
a fixed depth; a `maxcells` set is the best a fixed number of cells can say, and
its own reference level is the deepest level it reached. [`query`](@ref)
documents the keyword rules; the two are mutually exclusive.

[`is_contained`](@ref) reports which emissions were *proven* to fit — the first
kind, and only those; a cell emitted at the deepest level is never asked. The
shallowest cell of a set fits inside the target only when that flag says so, and
[`coarsest_contained`](@ref) is the accessor that asks.

!!! warning "Coverage is a statement about the LEAVES, not about the drawn cells"
    The set is a statement about the deepest level: every cell of that level
    which meets the target is the set's own member or the descendant of one, and
    no member is the descendant of another. That is the guarantee — and it is
    the one a lookup layer needs, because it makes the expansion
    ([`level_ranges`](@ref)) a superset of the single-level `Intersects` query,
    equal to it wherever the refinement is congruent.

    It is **not** a guarantee about the union of the emitted cells' polygons.
    Replacing a subtree by its root replaces the subtree's footprint by the
    root's, and those two agree only where the refinement is congruent. On the
    six shipped systems:

      - HEALPix, S2 and ISEA4R refine congruently — four children exactly tile
        their parent — and the emitted polygons do tile the target.
      - IGEO7 and H3 are aperture 7: the seven children are a rotated rosette
        that matches the parent in area but not in footprint. Roughly 3% of a
        state-sized target falls in slivers between a mixed-level set's cells.
      - A5's four children cover their parent's area without covering its
        footprint at all, and the figure is nearer 17%.

    Draw a multi-order set as a picture of *which cells were chosen*, and expand
    it before computing with it as a region. In the other direction the same
    non-congruence means a member's descendants can lie outside the target —
    inside a hole in it, for instance — so the expansion over-covers exactly
    where the refinement does.
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

The REFERENCE LEVEL is the depth the set speaks about: the `level` the query was
given, or — in `maxcells` mode, where no depth was asked for — the deepest level
the budget reached. It is the default expansion level for
[`level_ranges`](@ref), [`cellindices`](@ref) and `CellLookup`, and the level at
which the covering guarantee is stated.

!!! note "Expansion needs sorted subtrees"
    `level_ranges` is the compressed form of the set and exists only where
    `has_sorted_subtrees(sys)` holds; on A5 it throws, because a cell's
    descendants are scattered through their level rather than occupying one
    interval of it. `descendants(sys, c, l)` still names them, so the set can
    always be expanded — just not to a short list of ranges.

[`is_contained`](@ref) says which stored cells were proven to fit inside the
target — not which ones do, a distinction that docstring spells out;
[`cell_polygon`](@ref) and [`cell_polygons`](@ref) read the geometry of a
mixed-level member without the caller resolving a level grid per cell.
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
stop when the next replacement would not fit in the budget.

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
only by the children that meet the target, and a cell that no child of the
replacement would cover is not replaced at all.

It is EXACT, at every budget, on the three systems whose four children tile
their parent: HEALPix, S2 and ISEA4R. Where children do not tile their parent it
degrades in the way [`MultiOrderCoverage`](@ref)'s warning already describes, and
for the same reason — replacing a cell by its children swaps one footprint for
another. Measured on a state-sized outline as the fraction of the target lying
in no emitted cell: under 2% on IGEO7 and its authalic wrap, under 11% on H3,
under 25% on A5.

The LEAF statement the `level` mode makes — every cell of the reference level
that meets the target is a member or the descendant of one — is a law here on
those same three systems only. The `level` mode earns it everywhere by
descending into cells that miss the target, because a child can overhang its
parent, and carrying that descent all the way to `level`. A budget has no fixed
depth to carry it to; stopping early is the whole point, and a branch stopped
early at a cell that misses the target is a branch whose overhang is not
followed. Both statements are pinned per system, at three budgets and on four
targets, in `test/systems/crosssystem/multiorder_budget.jl`.

What a budget does NOT buy is a tight picture of the target: at ten cells the
set over-covers California by a wide margin, and it says so through
[`is_contained`](@ref) rather than by pretending otherwise.

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
    # Emissions at `maxlevel` are never asked, so `false` there records that
    # nothing was proven, not that the cell sticks out. `is_contained` documents
    # the asymmetry and why it is the contract.
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
    # The ONLY sound subtree prune is the covering law: a cell whose node extent
    # misses the target has no descendant that can meet it.
    #
    # A cell's own geometry is emphatically NOT a prune. Children overhang their
    # parents wherever the refinement is not congruent — under aperture 7 they
    # poke out past the parent's edges, which is the whole reason `node_extent`
    # exists — so a cell disjoint from the target can still have a child that
    # meets it, and descending only into cells that meet the target drops that
    # child from the coverage silently. The exact test below therefore decides
    # what is EMITTED, never what is descended into.
    #
    # Both prunes read the node extent and nothing else: the target's own cap
    # first, because it is one distance, then the boundary-arc proof, which is
    # what keeps the traversal output-sensitive when that cap is the whole
    # sphere (a target wider than a hemisphere) or merely much bigger than the
    # target inside it (any long or thin one).
    extent = node_extent(sys, c)
    intersects_cap(target.cap, extent) || return nothing
    _subtree_outside(target, extent) && return nothing
    lc = level(c)
    grid = grids[lc-top+1]
    meets = _matches(DE9IM.Intersects(nothing), target, grid, c)
    if lc >= maxlevel
        # Emitted so that the set covers the target, and flagged unproven
        # WITHOUT asking `Within`. Many of these cells do fit; asking would cost
        # one ~48 KB predicate call per boundary cell, thousands of them at a
        # deep level, to label cells the traversal is finished with. The flag is
        # a record of proof, not of geometry — see `is_contained`.
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
# The covering is maintained explicitly as {contained cells} ∪ {crossing cells}.
# Contained cells are never refined: a cell proven to lie inside the target is
# already the tightest thing that can be said about its own footprint, and
# splitting it spends budget to say the same thing in more words. Only crossing
# cells are candidates, and they are visited COARSEST FIRST, ties broken by
# position within the level — the curve order, so the schedule is a function of
# the inputs alone and of nothing else.
#
# Coarsest-first plus "children are one level deeper" means the queue is at all
# times a level and its successor, so the priority queue is spelled here as two
# vectors and a level counter rather than as a heap. That is not a shortcut: it
# is the same order a heap on `(level, key)` would pop, and it makes the
# breadth-first shape of the traversal visible.
#
# WHEN A REPLACEMENT DOES NOT FIT the cell is set aside and the traversal moves
# on to the next candidate rather than stopping. Nothing is lost by trying: the
# set only ever grows, so a replacement rejected once is rejected for good, and
# a single pass over the queue is enough. What is gained is the tail of the
# budget — a coarse crossing cell with seven intersecting children can overrun a
# budget that a finer one with two children still fits into, and stopping on the
# first miss would strand those cells unspent.
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

    contained = ID[]        # proven `Within`: never refined, never requeued
    stalled = ID[]          # crossing, and kept whole: the budget said no
    current = ID[]          # crossing, at the level being refined
    for c in rootcells(sys)
        _budget_admit!(contained, current, sys, target, c, grids, top, maxlevel)
    end
    # The running size of {contained} ∪ {crossing}, maintained rather than
    # recomputed: a replacement is committed exactly when the size it would
    # leave behind fits.
    total = length(contained) + length(current)

    kids_in = ID[]
    kids_out = ID[]
    for _ in top:(maxlevel-1)
        isempty(current) && break
        sort!(current; by=c -> (level(c), _budget_key(grids[level(c)-top+1], c)))
        next = ID[]
        for c in current
            empty!(kids_in)
            empty!(kids_out)
            for child in children(sys, c)
                _budget_admit!(kids_in, kids_out, sys, target, child, grids, top, maxlevel)
            end
            k = length(kids_in) + length(kids_out)
            # `k == 0` is a cell that meets the target and has no child that
            # does, which non-congruent refinement allows. Replacing it by
            # nothing would shrink the covering, so it is not a replacement.
            if k == 0 || total + k - 1 > maxcells
                push!(stalled, c)
                continue
            end
            append!(contained, kids_in)
            append!(next, kids_out)
            total += k - 1
        end
        current = next
    end
    append!(stalled, current)       # whatever `maxlevel` left unrefinable

    cells = vcat(contained, stalled)
    flags = falses(length(cells))
    flags[1:length(contained)] .= true
    # No depth was asked for, so the set speaks about the deepest level it
    # reached. An empty set reached nothing and reports the top.
    reference = isempty(cells) ? top : maximum(level, cells)
    return _sorted_cell_set(sys, cells, flags, reference)
end

# Classify one candidate cell, or reject it. Both prunes and both predicates are
# the ones `_coverage_visit!` uses, in the same order and for the same reasons —
# including the `maxlevel` arm, where `Within` is not asked because no decision
# depends on the answer and the call is the expensive one. `is_contained`
# documents what that leaves unproven.
function _budget_admit!(contained, crossing, sys, target, c, grids, top::Int,
        maxlevel::Int)
    extent = node_extent(sys, c)
    intersects_cap(target.cap, extent) || return false
    _subtree_outside(target, extent) && return false
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

# The tie-break inside one level. `cellposition` is the level's own order, which
# is curve order on every system here and is a bijection, so no two cells of a
# level ever tie and the schedule has nothing left to decide.
_budget_key(grid, c) = something(cellposition(grid, c), 0)

function _sorted_cell_set(sys::AbstractHierarchicalGridSystem, cells::Vector{ID},
        contained::BitVector, reference_level::Int) where {ID}
    if has_sorted_subtrees(sys)
        keys = [first(descendant_range(sys, c, reference_level)) for c in cells]
        perm = sortperm(keys)
        return MultiOrderCellSet{typeof(sys),ID}(sys, cells[perm], keys[perm],
            contained[perm], reference_level)
    end
    # No curve intervals to order by; `(level, id)` is the documented fallback,
    # and the keys become the cells' own positions within their level, which is
    # still a total order but not a curve order.
    perm = sortperm(cells; by=c -> (level(c), c))
    ordered = cells[perm]
    keys = [something(cellposition(levelgrid(sys, level(c)), c), 0) for c in ordered]
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
    is_contained(set::MultiOrderCellSet, i::Integer) -> Bool

Whether the set's `i`th cell was **proven** to lie inside the coverage target.
`true` means the traversal asked `Within` of that cell and it held. `false`
means one of two different things, told apart by the cell's level:

  - above the traversal's maximum depth, the cell *was* asked and the target's
    boundary crosses it — that is why the traversal descended into it. There
    the flag is exact in both directions.
  - at the maximum depth, the cell was **never asked**. The traversal ran out
    of depth and emitted it so that the set covers; it may fit inside the target
    or it may not.

"Maximum depth" is the `level` keyword in accuracy-first mode, where the set's
reference level and the maximum depth are the same number, so the blind spot is
exactly the reference level. In `maxcells` mode they are not: the reference
level is the deepest level the budget reached, the maximum depth is the
`maxlevel` cap it was allowed, and every cell the budget stopped short of the
cap *was* asked. A budget set whose refinement never reached its cap — which is
the ordinary case — therefore carries an exact flag on every member, and a
`true` at its deepest level is a real containment rather than a gap in the
record.

That asymmetry is the contract, not an oversight. `Within` costs on the order of
48 KB of allocation per call against 600 bytes for `Intersects`, and a deep
coverage finishes on thousands of reference-level cells; asking each of them
once more would cost hundreds of megabytes and change no cell of the result.
What the flag gives up is a label, and what it keeps is the direction that
matters: `true` implies inside. Code that needs exact containment at the
reference level asks `Within` itself, of the few cells it cares about.

`argmin(level, set)` is therefore *not* "the coarsest cell inside the target" —
every emission can be unproven. [`coarsest_contained`](@ref) reads this flag.
"""
is_contained(set::MultiOrderCellSet, i::Integer) = set.contained[i]

"""
    coarsest_contained(set::MultiOrderCellSet) -> cell id or `nothing`

The shallowest cell of `set` **proven** to lie inside the coverage target, or
`nothing` when no cell of it was — see [`is_contained`](@ref) for what "proven"
leaves out. Cells at the traversal's maximum depth are never tested, so a set
made only of them answers `nothing` even where some of them do fit; a target
smaller than one such cell is the clearest way to land there, not the only one.

A budget set answers `nothing` for a second and more ordinary reason: at ten
cells over a state, nothing has been refined far enough to fit inside it yet,
and the accessor says so instead of handing back the shallowest crossing cell.
Raise the budget and the answer appears.

```julia
set = query(sys, MultiOrderCoverage(tile); level = 12)
cell = coarsest_contained(set)          # `nothing`, or a cell that fits in `tile`
```

Ties are broken by the set's own order, so the answer is the first shallowest
cell in curve order and does not depend on how the traversal was scheduled.
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

# A `MultiOrderCellSet` is not a grid — it holds no positions and its cells are
# at different levels — but it does know which level grid each of its cells
# belongs to, and that is the only thing a caller was missing. `levelgrid` is
# O(1), so resolving it per cell costs nothing worth caching.

"""
    cell_boundary(set::MultiOrderCellSet, c) -> Vector{UnitSphericalPoint}
    cell_centroid(set::MultiOrderCellSet, c) -> UnitSphericalPoint
    cell_polygon(set::MultiOrderCellSet, c) -> GI.Polygon

The geometry of one member of a mixed-level set, read from
`levelgrid(system(set), level(c))`. Same values as asking that grid directly;
the set spares the caller from resolving it per cell.
"""
cell_boundary(set::MultiOrderCellSet, c::AbstractCellIndex) =
    cell_boundary(levelgrid(set.system, level(c)), c)

cell_centroid(set::MultiOrderCellSet, c::AbstractCellIndex) =
    cell_centroid(levelgrid(set.system, level(c)), c)

cell_polygon(set::MultiOrderCellSet, c::AbstractCellIndex) =
    cell_polygon(levelgrid(set.system, level(c)), c)

"""
    cell_polygons(set::MultiOrderCellSet) -> Vector{<:GI.Polygon}

Every cell of the set as a unit-sphere polygon, in the set's own order — what a
plot of a coverage needs, in one call:

```julia
poly(GO.transform(GO.GeographicFromUnitSphere(), cell_polygons(set)))
```
"""
cell_polygons(set::MultiOrderCellSet) =
    [cell_polygon(levelgrid(set.system, level(c)), c) for c in set.cells]

"""
    curve_keys(set::MultiOrderCellSet) -> Vector{Int}

The sort key of each cell, in the order the set stores them.

For a system with [`has_sorted_subtrees`](@ref) these are curve keys proper: the
start of each cell's position interval at the set's reference level, ascending,
with sibling intervals adjacent — which is what makes compaction and
binary-search membership cheap.

Without sorted subtrees there are no position intervals to key on, and the set
falls back to ordering cells by `(level, id)`. The keys are then each cell's own
position within its own level, so they ascend only *within* a level and restart
at the next one; they are reported for inspection, not to be compared across
levels.
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

!!! warning "Two things the expansion is not"
    It is not available everywhere. `has_sorted_subtrees(A5System())` is
    `false` — an A5 cell's descendants are scattered through their level rather
    than forming one interval of it — and this throws there. Generic code either
    branches on the trait or expands with `descendants(sys, c, l)`, which is
    always available and gives a list rather than ranges.

    It is not a covering of the target. A cell is in the set because the *cell*
    is inside the target; where the refinement is not congruent its descendants
    are not, so the expansion can name leaves the target does not touch — most
    visibly inside a hole. See [`MultiOrderCoverage`](@ref).
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
