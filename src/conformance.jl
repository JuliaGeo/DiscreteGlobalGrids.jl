# ---------------------------------------------------------------------------
# T3 — Conformance harness.
#
# The executable form of the interface contracts: `test_grid_interface(grid)`
# and `test_hierarchical_system(sys)`, property suites a third-party
# implementor runs with two calls. Cursor-free by design — it exercises
# interface primitives only, so it does not depend on the T2 fallbacks.
#
# The structure is two layers, and the reason is that a conformance harness has
# two different callers. Every law is implemented once as a function returning
# a `Vector{String}` of *problems* (empty means conforming); the `@testset`
# layer wraps those in labelled test sets so a failure names the violated
# contract, and the `check_*` predicates wrap the same functions in a `Bool` so
# a caller can assert that a deliberately broken implementation is *caught*
# without capturing test-set internals. The harness's own test suite does
# exactly that, which is how the harness is tested rather than merely run.
#
# Sampling is seeded (`rng` kwarg) and every law is checked on the same sampled
# cells, so two runs of the same call examine the same cells and a reported
# failure is reproducible by re-running it.
# ---------------------------------------------------------------------------

module Conformance

using Test

# `Random` is reached through `Test` rather than declared as a second package
# dependency: `Test` depends on it, so `Test.Random` is the same module object
# that `import Random` would bind, and shipping the harness inside the package
# costs exactly one new dependency instead of two.
const Random = Test.Random

import GeometryOps as GO
const USph = GO.UnitSpherical

import ..DiscreteGlobalGrids as DGG
using ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, Connectivity, Vertex, Edge

export test_grid_interface, test_hierarchical_system

"""Seed of the default RNG, so that an unparameterised run is reproducible."""
const DEFAULT_SEED = 20260813

"""How far a boundary vertex may sit off the unit sphere. "A few `eps`", generously."""
const DEFAULT_UNIT_ATOL = 1e-9

"""Angular slack (radians) allowed before a point is *outside* a node extent."""
const DEFAULT_CAP_ATOL = 1e-12

_default_rng() = Random.MersenneTwister(DEFAULT_SEED)

# ===========================================================================
# Small vector helpers
#
# Written out rather than taken from LinearAlgebra: the interface layer's only
# geometric dependency is the unit sphere itself, and three-component dot and
# cross products are less code than the import they would replace.
# ===========================================================================

_dot(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
_cross(a, b) = (a[2] * b[3] - a[3] * b[2],
                a[3] * b[1] - a[1] * b[3],
                a[1] * b[2] - a[2] * b[1])
_norm(a) = sqrt(_dot(a, a))
_chord(a, b) = sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2 + (a[3] - b[3])^2)

"""
    spherical_signed_area(ring) -> Float64

The signed area, in steradians, enclosed by the implicitly closed spherical
polygon `ring`.

Positive exactly when the ring winds **counter-clockwise seen from outside the
sphere**, which is the orientation [`cell_boundary`](@ref) requires. Computed
by summing the signed solid angles of the fan triangulation from the first
vertex with the Van Oosterom–Strackee formula, whose `atan2` form is stable for
the sliver triangles a densified boundary produces and whose sign is carried by
the numerator, so the sum is correct for non-convex rings too — the apex may be
any vertex, and a re-entrant ring's negative-area triangles cancel exactly.

!!! warning "Only defined modulo 4π"
    Each term lands in `(-2π, 2π)`, so a ring enclosing **more** than half the
    sphere reports its area minus 4π — negative, and indistinguishable from a
    clockwise ring of the complementary area. Winding is therefore only
    decidable this way for rings enclosing less than 2π steradians. That covers
    every cell of every grid with three or more cells, which is why
    [`boundary_problems`](@ref) uses it, but a deliberately hemispheric cell
    would need a winding-based test instead.
"""
function spherical_signed_area(ring)
    n = length(ring)
    n < 3 && return 0.0
    i0 = firstindex(ring)
    a = ring[i0]
    total = 0.0
    for i in (i0 + 1):(i0 + n - 2)
        b = ring[i]
        c = ring[i + 1]
        num = _dot(a, _cross(b, c))
        den = 1 + _dot(a, b) + _dot(b, c) + _dot(c, a)
        total += 2 * atan(num, den)
    end
    return total
end

# ===========================================================================
# Sampling
# ===========================================================================

"""
    sample_positions(rng, npositions, n_samples) -> Vector{Int}

`n_samples` distinct positions drawn from `1:npositions`, sorted; all of them
when the grid is small enough to check exhaustively.

Sampling is with replacement then deduplicated rather than a `randperm`,
because a level grid may have more cells than there is memory for a
permutation of them.
"""
function sample_positions(rng, npositions::Integer, n_samples::Integer)
    n = Int(npositions)
    n <= 0 && return Int[]
    n <= n_samples && return collect(1:n)
    return sort!(unique(rand(rng, 1:n, Int(n_samples))))
end

"""
    sample_cells(rng, grid, n_samples)

The typed ids at [`sample_positions`](@ref), paired with those positions.
"""
function sample_cells(rng, grid, n_samples::Integer)
    positions = sample_positions(rng, DGG.ncells(grid), n_samples)
    return positions, [DGG.cellindex(grid, i) for i in positions]
end

"""
    sample_levels(rng, levels, n_levels) -> Vector{Int}

`n_levels` levels drawn from `levels`, sorted and always including the first
and last: the coarsest and deepest levels are where level-dependent arithmetic
breaks, so they are never sampled away.
"""
function sample_levels(rng, levels, n_levels::Integer)
    ls = collect(Int, levels)
    length(ls) <= n_levels && return ls
    keep = Set{Int}((first(ls), last(ls)))
    while length(keep) < n_levels
        push!(keep, rand(rng, ls))
    end
    return sort!(collect(keep))
end

# ===========================================================================
# Base-interface laws, as problem collectors
# ===========================================================================

"""
    boundary_problems(pts; unit_atol) -> Vector{String}

Everything wrong with one [`cell_boundary`](@ref) ring: too few vertices,
vertices off the unit sphere, an explicitly repeated closing vertex (rings are
*implicitly* closed), or clockwise winding.

Winding is read off [`spherical_signed_area`](@ref) and so assumes a cell
encloses less than half the sphere; see that function's warning. A ring whose
area is within `area_atol` of zero is reported as degenerate rather than as
reversed, because at zero the sign carries no information.
"""
function boundary_problems(pts; unit_atol::Real = DEFAULT_UNIT_ATOL,
        area_atol::Real = 1e-12)
    problems = String[]
    n = length(pts)
    if n < 3
        push!(problems, "ring has $n vertices; a spherical polygon needs at least 3")
        return problems
    end

    worst = 0.0
    for p in pts
        worst = max(worst, abs(_norm(p) - 1))
    end
    if worst > unit_atol
        push!(problems, "ring vertex is off the unit sphere by $worst (tolerance $unit_atol)")
    end

    if _chord(first(pts), last(pts)) <= unit_atol
        push!(problems, "ring repeats its first vertex at the end; cell_boundary rings are implicitly closed")
    end

    area = spherical_signed_area(pts)
    if abs(area) <= area_atol
        push!(problems, "ring encloses no area ($area sr); its vertices are degenerate or collinear")
    elseif area < 0
        push!(problems,
            "spherical signed area is $area sr; the ring must wind counter-clockwise seen " *
            "from outside (a negative area also results from a ring enclosing more than " *
            "half the sphere, which this check cannot represent)")
    end
    return problems
end

"""
    centroid_problems(centroid, pts; unit_atol) -> Vector{String}

Everything wrong with one [`cell_centroid`](@ref): off the unit sphere, outside
its own cell's boundary ring, or *on* that ring — the contract requires a
strictly interior point, since a boundary point has no well-defined owning cell.

A ring whose containment is numerically indeterminate (`spherical_ring_contains`
returning `nothing`, which happens for near-hemispherical or vertex-symmetric
rings) is not reported as a problem: the harness never fails an implementation
on a predicate that declined to answer.
"""
function centroid_problems(centroid, pts; unit_atol::Real = DEFAULT_UNIT_ATOL)
    problems = String[]
    off = abs(_norm(centroid) - 1)
    if off > unit_atol
        push!(problems, "centroid is off the unit sphere by $off (tolerance $unit_atol)")
    end

    n = length(pts)
    n < 3 && return problems

    # `cell_boundary` may return any `AbstractVector`, and the ring predicates
    # index from 1, so anything else is normalised before it is handed over.
    v = firstindex(pts) == 1 ? pts : collect(pts)

    inside = USph.spherical_ring_contains(v, n, centroid)
    if inside === false
        push!(problems, "centroid lies outside its own cell boundary")
    end
    for j in 1:n
        if USph.point_on_spherical_arc(centroid, v[j], v[mod1(j + 1, n)])
            push!(problems,
                "centroid lies on the cell boundary; the contract requires a strictly interior point")
            break
        end
    end
    return problems
end

"""
    grid_interface_problems(grid; n_samples, rng, unit_atol) -> Vector{String}

Every base-interface law violated by `grid`, over a sampled set of positions.
Empty means conforming. This is the `Bool`-free core that both
[`test_grid_interface`](@ref) and [`check_grid_interface`](@ref) are written
against.
"""
function grid_interface_problems(grid;
        n_samples::Integer = 16,
        rng::Random.AbstractRNG = _default_rng(),
        unit_atol::Real = DEFAULT_UNIT_ATOL)
    problems = String[]
    n = DGG.ncells(grid)
    positions, cells = sample_cells(rng, grid, n_samples)

    if !allunique(cells)
        push!(problems, "cellindex is not injective: distinct positions returned the same id")
    end

    for (i, c) in zip(positions, cells)
        pos = DGG.cellposition(grid, c)
        pos == i || push!(problems, "cellposition(grid, cellindex(grid, $i)) == $pos, not $i")

        pts = DGG.cell_boundary(grid, c)
        for p in boundary_problems(pts; unit_atol)
            push!(problems, "cell $c: $p")
        end

        centroid = DGG.cell_centroid(grid, c)
        for p in centroid_problems(centroid, pts; unit_atol)
            push!(problems, "cell $c: $p")
        end

        if DGG.cell_boundary(grid, c) != pts || DGG.cell_centroid(grid, c) != centroid
            push!(problems, "cell $c: repeated calls returned different geometry")
        end
    end
    return problems
end

"""
    check_grid_interface(grid; kwargs...) -> Bool

Whether `grid` satisfies the base-interface laws — [`test_grid_interface`](@ref)
with the assertions collapsed into a single boolean, for callers that need to
assert a *failure* (a conformance harness's own tests) rather than a success.
"""
check_grid_interface(grid; kwargs...) = isempty(grid_interface_problems(grid; kwargs...))

# ===========================================================================
# Hierarchical laws, as problem collectors
# ===========================================================================

"""
    hierarchy_problems(sys, c) -> Vector{String}

The [`parent`](@ref)/[`children`](@ref) laws at one cell: the two are inverses,
children are distinct, non-empty, ascending, and exactly one level deeper.
"""
function hierarchy_problems(sys, c)
    problems = String[]
    l = DGG.level(c)
    levelrange = DGG.levels(sys)

    if l > first(levelrange)
        p = Base.parent(sys, c)
        if DGG.level(p) != l - 1
            push!(problems, "parent of a level-$l cell is at level $(DGG.level(p)), not $(l - 1)")
        end
        kids = collect(DGG.children(sys, p))
        c in kids || push!(problems, "children(parent($c)) does not contain $c")
    end

    if l < DGG.max_level(sys)
        kids = collect(DGG.children(sys, c))
        if isempty(kids)
            push!(problems, "$c has no children below max_level")
        end
        allunique(kids) || push!(problems, "children($c) contains duplicates")
        issorted(kids) || push!(problems, "children($c) is not in ascending canonical order")
        for k in kids
            if DGG.level(k) != l + 1
                push!(problems, "child $k of a level-$l cell is at level $(DGG.level(k)), not $(l + 1)")
            end
            back = Base.parent(sys, k)
            back == c || push!(problems, "parent(child $k) == $back, not $c")
        end
    end
    return problems
end

"""
    ring_sample_points(pts, densify) -> Vector

The vertices of a boundary ring, optionally with `densify` interior points per
edge interpolated along the great-circle arc.

Vertices alone are a sound proxy for the whole boundary when the region they
are being tested against is geodesically convex (see the covering law), which
for a `SphericalCap` means an angular radius of at most 90°. `densify` is what
the harness reaches for when that does not hold.
"""
function ring_sample_points(pts, densify::Integer)
    verts = collect(pts)
    densify <= 0 && return verts
    n = length(verts)
    out = similar(verts, 0)
    for j in 1:n
        a = verts[j]
        b = verts[mod1(j + 1, n)]
        push!(out, a)
        for t in 1:densify
            push!(out, USph.slerp(a, b, t / (densify + 1)))
        end
    end
    return out
end

"""
    cap_overshoot(cap, p) -> Float64

How far outside `cap` the point `p` lies, in radians; non-positive when
contained. The signed form of `GO.UnitSpherical._contains`, so that a covering
violation can be reported with its magnitude rather than as a bare `false`.
"""
cap_overshoot(cap, p) = USph.spherical_distance(cap.point, p) - cap.radius

"""
    node_extent_problems(cap) -> Vector{String}

Whether a [`node_extent`](@ref) is a well-formed `SphericalCap` at all: a
unit-norm centre and a finite, non-negative radius.
"""
function node_extent_problems(cap; unit_atol::Real = DEFAULT_UNIT_ATOL)
    problems = String[]
    cap isa USph.SphericalCap ||
        return push!(problems, "node_extent returned a $(typeof(cap)), not a SphericalCap")
    off = abs(_norm(cap.point) - 1)
    off > unit_atol && push!(problems, "node extent centre is off the unit sphere by $off")
    isfinite(cap.radius) || push!(problems, "node extent radius is $(cap.radius)")
    cap.radius < 0 && push!(problems, "node extent radius is negative ($(cap.radius))")
    return problems
end

"""
    covering_law_problems(sys; levels, n_samples, rng, descent_depth,
                          branch_samples, densify, arc_samples, atol) -> Vector{String}

Violations of **the covering law**: `node_extent(sys, c)` must contain the
geometry of every descendant of `c` at every depth.

For each sampled cell the harness walks `descent_depth` levels of its subtree,
following `branch_samples` children at each step, and asserts that every
boundary point of every cell it reaches lies inside the node extent of *every*
ancestor on the path — not merely its immediate parent, which is the weaker
statement a naive recursion would check and which a system can satisfy while
still violating the law.

Depth alone is not the whole law — an extent that covers three levels and fails
at four still violates it — so alongside the bounded bushy walk each sampled
cell also follows one random chain of children all the way to
[`max_level`](@ref), which costs `O(max_level)` cells. Pass `deep_chain =
false` to skip it, or `cells` (a `level => cells` mapping) to check exactly the
cells another law sampled.

Boundary **vertices** are the default proxy for the boundary, which is exact
when the ancestor extents are geodesically convex (radius ≤ 90°, so the cap
contains the great-circle arc between any two points it contains). When a
sampled extent is larger than that the harness silently switches to
`arc_samples` interpolated points per edge, so the check stays sound;
[`test_hierarchical_system`](@ref) additionally *asserts* convexity, so the
weakening is reported rather than hidden. Pass `densify` to force a point
count either way.
"""
function covering_law_problems(sys;
        levels = DGG.levels(sys),
        n_samples::Integer = 8,
        rng::Random.AbstractRNG = _default_rng(),
        descent_depth::Integer = 3,
        branch_samples::Integer = 2,
        deep_chain::Bool = true,
        cells = nothing,
        densify::Union{Nothing,Integer} = nothing,
        arc_samples::Integer = 2,
        atol::Real = DEFAULT_CAP_ATOL)
    problems = String[]
    cache = Dict{Int,Any}()
    getgrid(l) = get!(() -> DGG.levelgrid(sys, l), cache, Int(l))

    for l in levels
        grid = getgrid(l)
        sampled = cells === nothing ? last(sample_cells(rng, grid, n_samples)) : cells[l]
        for c in sampled
            root = [(c, DGG.node_extent(sys, c))]
            _covering_descend!(problems, sys, c, root, getgrid, rng,
                Int(descent_depth), Int(branch_samples), densify, Int(arc_samples), atol)
            # The bushy walk above is bounded in depth, but the law is not: an
            # extent that covers three levels and fails at four is still a
            # violation, and `cap_inflation` is documented as a limiting
            # quantity. One chain per cell reaches max_level for O(depth) cells.
            deep_chain && _covering_chain!(problems, sys, c, root, getgrid, rng,
                densify, Int(arc_samples), atol)
        end
    end
    return problems
end

"Check one cell's boundary against every ancestor extent on its path."
function _covering_check!(problems, sys, c, ancestors, getgrid, densify, arc_samples, atol)
    l = DGG.level(c)
    pts = DGG.cell_boundary(getgrid(l), c)

    # Vertices suffice against convex (≤ 90°) caps; densify against anything wider.
    d = densify === nothing ?
        (all(a -> last(a).radius <= π / 2 + 1e-9, ancestors) ? 0 : arc_samples) :
        Int(densify)
    samples = ring_sample_points(pts, d)

    for (ac, acap) in ancestors
        worst = -Inf
        for p in samples
            worst = max(worst, cap_overshoot(acap, p))
        end
        if worst > atol
            push!(problems,
                "covering law: the boundary of $c escapes node_extent($ac) by $worst rad " *
                "(ancestor at level $(DGG.level(ac)), descendant at level $l)")
        end
    end
    return problems
end

function _covering_descend!(problems, sys, c, ancestors, getgrid, rng,
        depth, branch_samples, densify, arc_samples, atol)
    _covering_check!(problems, sys, c, ancestors, getgrid, densify, arc_samples, atol)

    (depth <= 0 || DGG.level(c) >= DGG.max_level(sys)) && return problems

    kids = collect(DGG.children(sys, c))
    isempty(kids) && return problems
    chosen = length(kids) <= branch_samples ? kids :
        kids[sort!(unique(rand(rng, 1:length(kids), branch_samples)))]
    for k in chosen
        _covering_descend!(problems, sys, k, vcat(ancestors, [(k, DGG.node_extent(sys, k))]),
            getgrid, rng, depth - 1, branch_samples, densify, arc_samples, atol)
    end
    return problems
end

"Follow a single random child per level to `max_level`, checking all the way down."
function _covering_chain!(problems, sys, c, ancestors, getgrid, rng, densify, arc_samples, atol)
    cur = c
    chain = ancestors
    while DGG.level(cur) < DGG.max_level(sys)
        kids = collect(DGG.children(sys, cur))
        isempty(kids) && break
        cur = kids[rand(rng, 1:length(kids))]
        chain = vcat(chain, [(cur, DGG.node_extent(sys, cur))])
        _covering_check!(problems, sys, cur, chain, getgrid, densify, arc_samples, atol)
    end
    return problems
end

"""
    check_covering_law(sys; kwargs...) -> Bool

Whether `sys` obeys the covering law over a sampled walk of its subtrees. The
boolean face of [`covering_law_problems`](@ref); a `false` here is the single
most consequential conformance failure there is, because a node extent that
under-covers makes every pruned traversal silently drop cells.
"""
check_covering_law(sys; kwargs...) = isempty(covering_law_problems(sys; kwargs...))

"""
    neighbor_problems(grid, c; connectivity, sys) -> Vector{String}

The [`neighbors`](@ref) laws at one cell: determinism across calls, exclusion of
`c` itself, distinctness, a container typed at the system's cell index type, a
count within [`max_neighbors`](@ref), membership in the grid, a common level,
and symmetry — `c′ ∈ neighbors(c)` implies `c ∈ neighbors(c′)`.

Partial coverage is handled by the contract rather than by an exemption: a
neighbour beyond the grid's coverage must be **absent** from the result, so a
returned cell that has no [`cellposition`](@ref) is itself a violation, and
every cell that *is* returned is in the grid and can be asked for its own
neighbours in turn. Symmetry is therefore total over whatever the grid returns.
"""
function neighbor_problems(grid, c; connectivity::Connectivity = Vertex(), sys = DGG.system(grid))
    problems = String[]
    ns = collect(DGG.neighbors(grid, c, 1; connectivity))

    collect(DGG.neighbors(grid, c, 1; connectivity)) == ns ||
        push!(problems, "neighbors($c) is not deterministic across calls")
    c in ns && push!(problems, "neighbors($c) contains $c itself")
    allunique(ns) || push!(problems, "neighbors($c) contains duplicates")

    if sys !== nothing
        T = DGG.cellindextype(sys)
        eltype(ns) === T ||
            push!(problems, "neighbors($c) has eltype $(eltype(ns)), not the system's $T")
        bound = DGG.max_neighbors(sys, connectivity)
        length(ns) <= bound ||
            push!(problems, "neighbors($c) returned $(length(ns)) cells, over max_neighbors of $bound")
    end

    lc = DGG.level(c)
    for nb in ns
        DGG.level(nb) == lc ||
            push!(problems, "neighbour $nb of $c is at level $(DGG.level(nb)), not $lc")
        if DGG.cellposition(grid, nb) === nothing
            push!(problems, "neighbour $nb of $c is not a cell of the grid")
            continue
        end
        back = collect(DGG.neighbors(grid, nb, 1; connectivity))
        c in back || push!(problems, "neighbours are not symmetric: $nb ∈ neighbors($c) but $c ∉ neighbors($nb)")
    end

    isempty(collect(DGG.neighbors(grid, c, 0; connectivity))) ||
        push!(problems, "neighbors($c, 0) is not empty")
    return problems
end

"""
    check_neighbors(grid; n_samples, rng, connectivity) -> Bool

Whether `grid`'s neighbour relation satisfies [`neighbor_problems`](@ref)'s laws
over sampled cells.
"""
function check_neighbors(grid;
        n_samples::Integer = 16,
        rng::Random.AbstractRNG = _default_rng(),
        connectivity::Connectivity = Vertex())
    _, cells = sample_cells(rng, grid, n_samples)
    return all(c -> isempty(neighbor_problems(grid, c; connectivity)), cells)
end

"""
    descendants_at(sys, c, l) -> Vector

Every level-`l` descendant of `c`, computed by recursing on [`children`](@ref)
alone. This is the *definition* the harness holds
[`descendant_range`](@ref) and `descendants` against, so it deliberately
duplicates work a system can do in closed form: an oracle that shared the
system's shortcut would not be one.
"""
function descendants_at(sys, c, l::Integer)
    lc = DGG.level(c)
    lc == l && return [c]
    out = Vector{DGG.cellindextype(sys)}()
    for k in DGG.children(sys, c)
        append!(out, descendants_at(sys, k, l))
    end
    return out
end

"""
    descendant_range_problems(sys, c, l, grid) -> Vector{String}

Violations of the two-sided [`descendant_range`](@ref) contract at one cell:
the positions of `c`'s actual level-`l` descendants in `grid` must *exactly*
fill the returned range — every descendant in it, and every position in it a
descendant.
"""
function descendant_range_problems(sys, c, l::Integer, grid)
    problems = String[]
    r = DGG.descendant_range(sys, c, l)
    r isa UnitRange{Int} ||
        push!(problems, "descendant_range($c, $l) returned a $(typeof(r)), not a UnitRange{Int}")

    actual = descendants_at(sys, c, l)
    positions = Int[]
    for d in actual
        pos = DGG.cellposition(grid, d)
        if pos === nothing
            push!(problems, "level-$l descendant $d of $c has no position in its own level grid")
        else
            push!(positions, pos)
        end
    end

    # A range wildly wider than the descendant set is itself the violation, and
    # materialising it to say so is how a harness runs a machine out of memory.
    if length(r) > 8 * max(length(positions), 1) + 64
        push!(problems,
            "descendant_range($c, $l) spans $(length(r)) positions for $(length(positions)) descendants")
        return problems
    end

    missed = setdiff(positions, r)
    extra = setdiff(collect(r), positions)
    isempty(missed) ||
        push!(problems, "descendant_range($c, $l) = $r omits descendant positions $missed")
    isempty(extra) ||
        push!(problems, "descendant_range($c, $l) = $r contains $(length(extra)) positions that are not descendants of $c")
    return problems
end

"""
    check_descendant_ranges(sys; levels, n_samples, rng, depth) -> Bool

Whether [`descendant_range`](@ref) exactly matches the descendants reached by
[`children`](@ref) recursion, over sampled cells and depths.
"""
function check_descendant_ranges(sys;
        levels = DGG.levels(sys),
        n_samples::Integer = 8,
        rng::Random.AbstractRNG = _default_rng(),
        depth::Integer = 2)
    DGG.has_sorted_subtrees(sys) || return true
    maxl = DGG.max_level(sys)
    for l in levels
        grid = DGG.levelgrid(sys, l)
        _, cells = sample_cells(rng, grid, n_samples)
        for c in cells, d in 1:depth
            l + d > maxl && continue
            isempty(descendant_range_problems(sys, c, l + d, DGG.levelgrid(sys, l + d))) || return false
        end
    end
    return true
end

"""
    hierarchical_system_problems(sys; levels, n_samples, rng, kwargs...) -> Vector{String}

Every hierarchical law violated by `sys`: [`hierarchy_problems`](@ref) over
sampled cells at sampled levels, plus the covering law. The boolean face is
[`check_hierarchical_system`](@ref).
"""
function hierarchical_system_problems(sys;
        levels = DGG.levels(sys),
        n_samples::Integer = 8,
        rng::Random.AbstractRNG = _default_rng(),
        descent_depth::Integer = 3,
        branch_samples::Integer = 2)
    problems = String[]
    sampled = Dict{Int,Any}()
    for l in levels
        grid = DGG.levelgrid(sys, l)
        _, cells = sample_cells(rng, grid, n_samples)
        sampled[Int(l)] = cells
        for c in cells
            append!(problems, hierarchy_problems(sys, c))
        end
    end
    # The same cells, so a covering failure names a cell the other laws saw too.
    append!(problems, covering_law_problems(sys; levels, rng, cells = sampled,
        descent_depth, branch_samples))
    return problems
end

"""
    check_hierarchical_system(sys; kwargs...) -> Bool

Whether `sys` satisfies the hierarchical laws.
"""
check_hierarchical_system(sys; kwargs...) = isempty(hierarchical_system_problems(sys; kwargs...))

# ===========================================================================
# The two entry points
# ===========================================================================

"""
    test_grid_interface(grid; n_samples = 32, rng = MersenneTwister(20260813),
                        unit_atol = 1e-9, label = <grid type>)

Property-test `grid` against the base-interface contracts, as a labelled
`Test.@testset`: the [`cellindex`](@ref)/[`cellposition`](@ref) bijection over
sampled positions (including `nothing` for a cell that is not in the grid),
boundary rings of unit-norm points that are implicitly closed and wind
counter-clockwise seen from outside, centroids strictly inside their own cell,
`cellat(cell_centroid(grid, c)) == c` where [`cellat`](@ref) is implemented, and
determinism of repeated calls.

Each law is its own nested test set, so a failure names the contract it
violated rather than a line number. Every law examines the same sampled cells,
which are drawn once per call from a freshly seeded RNG — so two identical
calls examine identical cells, and a failure is reproduced by re-running it.
Passing an `rng` explores a different sample; note that an RNG is mutable and
advances as it is drawn from, so reproducing a specific run means constructing
a *new* generator from the same seed, not reusing the one it consumed. Up to
`n_samples` distinct positions are drawn (sampling is with replacement and then
deduplicated, so a large request may yield somewhat fewer), or every position
when the grid is smaller than that.

Optional methods are skipped, not failed: a grid that does not implement
[`cellat`](@ref) reports that test set as skipped and conforms regardless.

Out of scope by design, because this harness is **cursor-free** and so cannot
depend on the generic fallback layer: [`treeify`](@ref), [`query`](@ref),
[`cell_polygon`](@ref), [`cell_area`](@ref), [`cell_extent`](@ref) and
[`getcell`](@ref). Those are derived from the four required primitives this
suite does check, and belong to the fallback layer's own tests.

A third-party grid implementor calls this and [`test_hierarchical_system`](@ref)
to validate an implementation.

See also [`check_grid_interface`](@ref), the same laws as a `Bool`.
"""
function test_grid_interface(grid;
        n_samples::Integer = 32,
        rng::Random.AbstractRNG = _default_rng(),
        unit_atol::Real = DEFAULT_UNIT_ATOL,
        label::AbstractString = string(nameof(typeof(grid))))
    n = DGG.ncells(grid)
    positions, cells = sample_cells(rng, grid, n_samples)
    sys = DGG.system(grid)

    @testset "grid interface: $label" begin
        @testset "ncells" begin
            @test n isa Integer
            @test n >= 0
            @test DGG.ncells(grid) == n   # stable over the grid's lifetime
        end

        @testset "provenance (system/level)" begin
            @test sys === nothing || sys isa AbstractHierarchicalGridSystem
            l = DGG.level(grid)
            @test l === nothing || l isa Integer
            # A grid made by a system reports both; a standalone grid reports neither.
            @test (sys === nothing) == (l === nothing)
            if sys !== nothing
                @test l in DGG.levels(sys)
                @test all(c -> DGG.level(c) == l, cells)
                @test all(c -> c isa DGG.cellindextype(sys), cells)
                # A system grid's dense order IS the system's canonical id order
                # (`isless`), which is what every `searchsorted` fast path
                # assumes. A standalone grid may choose its own order, so this
                # is asserted only for grids that have a system.
                @test issorted(cells)
            end
        end

        @testset "cellindex/cellposition bijection" begin
            @test all(c -> c isa AbstractCellIndex, cells)
            @test allunique(cells)
            for (i, c) in zip(positions, cells)
                @test DGG.cellposition(grid, c) == i
            end
        end

        @testset "cellindex bounds" begin
            @test_throws BoundsError DGG.cellindex(grid, 0)
            @test_throws BoundsError DGG.cellindex(grid, n + 1)
        end

        @testset "cellposition of a cell outside the grid" begin
            # A cell one level down is, by the level contract, not in this grid.
            other = _foreign_cell(grid, sys, cells)
            if other === nothing
                @test_skip "no foreign cell constructible from the interface alone"
            else
                @test DGG.cellposition(grid, other) === nothing
            end
        end

        @testset "cell_boundary rings" begin
            for c in cells
                @test boundary_problems(DGG.cell_boundary(grid, c); unit_atol) == String[]
            end
        end

        @testset "cell_centroid" begin
            for c in cells
                @test centroid_problems(DGG.cell_centroid(grid, c),
                    DGG.cell_boundary(grid, c); unit_atol) == String[]
            end
        end

        @testset "cellat(cell_centroid(c)) == c" begin
            if isempty(cells) || !applicable(DGG.cellat, grid, DGG.cell_centroid(grid, first(cells)))
                @test_skip "cellat is not implemented for this grid"
            else
                for c in cells
                    centroid = DGG.cell_centroid(grid, c)
                    @test DGG.cellat(grid, centroid) == c
                    # Interior points that are not the centroid. Without these a
                    # nearest-centroid lookup conforms, and nearest-centroid is
                    # wrong for every tessellation that is not a Voronoi one.
                    for v in DGG.cell_boundary(grid, c)
                        @test DGG.cellat(grid, USph.slerp(centroid, v, 0.5)) == c
                    end
                end
            end
        end

        @testset "determinism of repeated calls" begin
            for (i, c) in zip(positions, cells)
                @test DGG.cellindex(grid, i) == c
                @test DGG.cell_boundary(grid, c) == DGG.cell_boundary(grid, c)
                @test DGG.cell_centroid(grid, c) == DGG.cell_centroid(grid, c)
                @test DGG.cellposition(grid, c) == DGG.cellposition(grid, c)
            end
        end
    end
end

"""
    _foreign_cell(grid, sys, cells) -> Union{AbstractCellIndex,Nothing}

A cell that must *not* be in `grid`, or `nothing` if none can be built from the
interface alone.

A same-level cell that a partial grid does not cover is preferred, because that
is the case exercising a real membership test; most implementations answer the
other case, a cell at `level(grid) + 1`, with a level comparison alone. The
child is the fallback, and it is portable to complete grids.
"""
function _foreign_cell(grid, sys, cells)
    (sys === nothing || isempty(cells)) && return nothing
    l = DGG.level(grid)
    l === nothing && return nothing

    full = DGG.levelgrid(sys, l)
    if DGG.ncells(full) > DGG.ncells(grid)
        for i in 1:min(DGG.ncells(full), 4096)
            c = DGG.cellindex(full, i)
            DGG.cellposition(grid, c) === nothing && return c
        end
    end

    l >= DGG.max_level(sys) && return nothing
    kids = collect(DGG.children(sys, first(cells)))
    return isempty(kids) ? nothing : first(kids)
end

"""
    test_hierarchical_system(sys; levels = levels(sys), n_levels = 4,
                             n_samples = 8, rng = MersenneTwister(20260813),
                             descent_depth = 3, branch_samples = 2,
                             connectivities = (Vertex(), Edge()),
                             label = <system type>)

Property-test `sys` against the hierarchical contracts, as a labelled
`Test.@testset`.

The laws, each its own nested test set:

  - **`levels`/`rootcells`** — a non-empty level range whose last element is
    [`max_level`](@ref); roots all at the first level, ascending and distinct,
    and exactly the cells of the coarsest level grid.
  - **`levelgrid` consistency** — `system`/`level` agree with the system,
    `cellindex`/`cellposition` round-trip, and cell counts increase with level.
  - **`parent`/`children`** — mutual inverses; children distinct, non-empty,
    ascending and one level deeper; `parent` throws an `ArgumentError` on a
    root and `children` throws one at `max_level`.
  - **the covering law** — see [`covering_law_problems`](@ref): sampled
    subtrees walked `descent_depth` levels down, every descendant's boundary
    contained in *every* ancestor's [`node_extent`](@ref). Node extents are
    also asserted to be geodesically convex, which is what makes sampling
    boundary vertices a sound proxy for the whole boundary.
  - **neighbours** — determinism, symmetry, `k = 0` empty, `ring(grid, c, 0) ==
    [c]`, `Edge()` a subset of `Vertex()`, all under both connectivities.
  - **`descendant_range`** — when [`has_sorted_subtrees`](@ref) is `true`, the
    positions of the actual descendants exactly fill the returned range.

`levels` selects the levels to test; when it has more than `n_levels` entries a
seeded sample of that many is used, always including the coarsest and deepest.
Derived methods that a system has not implemented (`ancestor`, `descendants`,
[`neighbors`](@ref), [`ring`](@ref)) are skipped rather than failed.

See also [`check_hierarchical_system`](@ref) and [`check_covering_law`](@ref),
the same laws as `Bool`s.
"""
function test_hierarchical_system(sys;
        levels = DGG.levels(sys),
        n_levels::Integer = 4,
        n_samples::Integer = 8,
        rng::Random.AbstractRNG = _default_rng(),
        descent_depth::Integer = 3,
        branch_samples::Integer = 2,
        require_convex_extents::Bool = true,
        neighbor_k::Integer = 2,
        connectivities = (Vertex(), Edge()),
        label::AbstractString = string(nameof(typeof(sys))))
    levelrange = DGG.levels(sys)
    tested = sample_levels(rng, levels, n_levels)
    maxl = DGG.max_level(sys)
    grids = Dict{Int,Any}(l => DGG.levelgrid(sys, l) for l in tested)
    samples = Dict{Int,Any}(l => sample_cells(rng, grids[l], n_samples) for l in tested)

    @testset "hierarchical system: $label" begin
        @testset "levels and traits" begin
            @test levelrange isa AbstractUnitRange{Int}
            @test !isempty(levelrange)
            @test maxl == last(levelrange)
            @test DGG.cellindextype(sys) <: AbstractCellIndex
            @test DGG.cap_inflation(sys) >= 1
            @test DGG.has_sorted_subtrees(sys) isa Bool
            for conn in connectivities
                @test DGG.max_neighbors(sys, conn) >= 1
            end
            @test DGG.max_neighbors(sys) == DGG.max_neighbors(sys, Vertex())
            if Vertex() in connectivities && Edge() in connectivities
                # Edge() is a restriction of Vertex(), so its bound cannot be larger.
                @test DGG.max_neighbors(sys, Edge()) <= DGG.max_neighbors(sys, Vertex())
            end
            @test_throws ArgumentError DGG.levelgrid(sys, first(levelrange) - 1)
            @test_throws ArgumentError DGG.levelgrid(sys, maxl + 1)
        end

        @testset "rootcells" begin
            roots = collect(DGG.rootcells(sys))
            @test !isempty(roots)
            @test allunique(roots)
            @test issorted(roots)
            @test all(c -> DGG.level(c) == first(levelrange), roots)
            root_grid = DGG.levelgrid(sys, first(levelrange))
            @test DGG.ncells(root_grid) == length(roots)
            @test Set(roots) == Set(DGG.cellindex(root_grid, i) for i in 1:DGG.ncells(root_grid))
        end

        @testset "levelgrid consistency" begin
            for l in tested
                grid = grids[l]
                @test DGG.system(grid) === sys
                @test DGG.level(grid) == l
                @test DGG.ncells(grid) > 0
                positions, cells = samples[l]
                for (i, c) in zip(positions, cells)
                    @test DGG.level(c) == l
                    @test DGG.cellposition(grid, c) == i
                end
            end
            # Refinement adds cells: counts strictly increase with level.
            for (a, b) in zip(tested, Iterators.drop(tested, 1))
                @test DGG.ncells(grids[b]) > DGG.ncells(grids[a])
            end
        end

        @testset "parent/children inverses" begin
            for l in tested, c in last(samples[l])
                @test hierarchy_problems(sys, c) == String[]
            end
        end

        @testset "parent throws on roots, children throws at max_level" begin
            @test_throws ArgumentError Base.parent(sys, first(DGG.rootcells(sys)))
            deepest = DGG.levelgrid(sys, maxl)
            @test_throws ArgumentError DGG.children(sys, DGG.cellindex(deepest, 1))
        end

        @testset "ancestor/descendants" begin
            probe = last(samples[first(tested)])
            if isempty(probe) || !applicable(DGG.ancestor, sys, first(probe), 0)
                @test_skip "ancestor is not implemented for this system"
            else
                for l in tested, c in last(samples[l])
                    @test DGG.ancestor(sys, c, l) == c
                    # Every depth, not just one: a closed-form "drop level(c) - l
                    # digits" override goes wrong at depth 2, never at depth 1.
                    walked = c
                    for j in (l - 1):-1:first(levelrange)
                        walked = Base.parent(sys, walked)
                        @test DGG.ancestor(sys, c, j) == walked
                    end
                    l < maxl && @test_throws ArgumentError DGG.ancestor(sys, c, l + 1)
                end
            end
            if isempty(probe) || !applicable(DGG.descendants, sys, first(probe), 0)
                @test_skip "descendants is not implemented for this system"
            else
                for l in tested, c in last(samples[l])
                    @test collect(DGG.descendants(sys, c, l)) == [c]
                    # Depth 2 as well as 1: at depth 1 the oracle is just
                    # `children`, so depth 1 alone tests nothing new.
                    for d in 1:2
                        l + d > maxl && continue
                        @test collect(DGG.descendants(sys, c, l + d)) ==
                              sort!(descendants_at(sys, c, l + d))
                    end
                    l > first(levelrange) && @test_throws ArgumentError DGG.descendants(sys, c, l - 1)
                end
            end
        end

        @testset "node_extent well-formed" begin
            for l in tested, c in last(samples[l])
                @test node_extent_problems(DGG.node_extent(sys, c)) == String[]
            end
        end

        @testset "node_extent convexity (vertex-sampling proxy)" begin
            # A cap of radius ≤ 90° is geodesically convex, so containing a
            # boundary's vertices implies containing its arcs. `node_extent`
            # explicitly permits wider extents ("owes the full law, not the
            # proxy"), and the covering check densifies for them, so a system
            # that means it opts out with `require_convex_extents = false`.
            if !require_convex_extents
                @test_skip "convex node extents not required; the covering check densifies"
            else
                for l in tested, c in last(samples[l])
                    @test DGG.node_extent(sys, c).radius <= π / 2
                end
            end
        end

        @testset "covering law" begin
            # The same sampled cells every other law used, so a covering failure
            # is attributable to a cell the rest of the report also names.
            @test covering_law_problems(sys; levels = tested, rng,
                cells = Dict(l => last(samples[l]) for l in tested),
                descent_depth, branch_samples) == String[]
        end

        @testset "neighbors" begin
            probe = last(samples[first(tested)])
            if isempty(probe) || !applicable(DGG.neighbors, grids[first(tested)], first(probe), 1)
                @test_skip "neighbors is not implemented for this system"
            else
                for conn in connectivities, l in tested
                    grid = grids[l]
                    for c in last(samples[l])
                        @test neighbor_problems(grid, c; connectivity = conn, sys) == String[]
                    end
                end
                # Edge() is the opt-in restriction of the Vertex() default.
                if Vertex() in connectivities && Edge() in connectivities
                    for l in tested, c in last(samples[l])
                        vs = collect(DGG.neighbors(grids[l], c, 1; connectivity = Vertex()))
                        es = collect(DGG.neighbors(grids[l], c, 1; connectivity = Edge()))
                        @test issubset(Set(es), Set(vs))
                    end
                end
            end
        end

        @testset "ring" begin
            probe = last(samples[first(tested)])
            if isempty(probe) || !applicable(DGG.ring, grids[first(tested)], first(probe), 0)
                @test_skip "ring is not implemented for this system"
            else
                for conn in connectivities, l in tested
                    grid = grids[l]
                    for c in last(samples[l])
                        @test collect(DGG.ring(grid, c, 0; connectivity = conn)) == [c]
                        # `neighbors(c, k)` is the union of the shells 1..k, and
                        # the shells are disjoint. At k = 1 this collapses to an
                        # identity every implementation gets right, so it is the
                        # k ≥ 2 case that carries the weight.
                        shells = Set{DGG.cellindextype(sys)}()
                        for k in 1:neighbor_k
                            shell = Set(DGG.ring(grid, c, k; connectivity = conn))
                            @test isempty(intersect(shell, shells))
                            union!(shells, shell)
                            @test Set(DGG.neighbors(grid, c, k; connectivity = conn)) == shells
                        end
                    end
                end
            end
        end

        @testset "descendant_range" begin
            if !DGG.has_sorted_subtrees(sys)
                @test_skip "has_sorted_subtrees(sys) is false; descendant_range is not offered"
            else
                for l in tested, c in last(samples[l])
                    if l > first(levelrange)
                        @test_throws ArgumentError DGG.descendant_range(sys, c, l - 1)
                    end
                    @test DGG.descendant_range(sys, c, l) ==
                          DGG.cellposition(grids[l], c):DGG.cellposition(grids[l], c)
                    for d in 1:2
                        l + d > maxl && continue
                        target = get!(() -> DGG.levelgrid(sys, l + d), grids, l + d)
                        @test descendant_range_problems(sys, c, l + d, target) == String[]
                    end
                    # Sibling ranges partition the parent's range, in order.
                    if l < maxl
                        kids = collect(DGG.children(sys, c))
                        target = get!(() -> DGG.levelgrid(sys, l + 1), grids, l + 1)
                        ranges = [DGG.descendant_range(sys, k, l + 1) for k in kids]
                        @test reduce(vcat, collect.(ranges); init = Int[]) ==
                              collect(DGG.descendant_range(sys, c, l + 1))
                    end
                end
            end
        end
    end
end

end # module Conformance
