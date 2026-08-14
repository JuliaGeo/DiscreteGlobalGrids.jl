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

  - above the set's reference level, the cell *was* asked and the target's
    boundary crosses it — that is why the traversal descended into it. There
    the flag is exact in both directions.
  - at the reference level, the cell was **never asked**. The traversal ran out
    of depth and emitted it so that the set covers; it may fit inside the target
    or it may not.

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
`nothing` when no cell above the set's reference level was — see
[`is_contained`](@ref) for what "proven" leaves out. Reference-level cells are
never tested, so a set made only of them answers `nothing` even where some of
them do fit; a target smaller than one such cell is the clearest way to land
there, not the only one.

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
