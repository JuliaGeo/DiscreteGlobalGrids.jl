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

module DiscreteGlobalGridsConformanceTesting

using Test
import Random

import GeometryOps as GO
const USph = GO.UnitSpherical

import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, Connectivity, Vertex, Edge

export test_grid_interface, test_hierarchical_system

"""Seed of the default RNG, so that an unparameterised run is reproducible."""
const DEFAULT_SEED = 20260813

"""How far a boundary vertex may sit off the unit sphere. "A few `eps`", generously."""
const DEFAULT_UNIT_ATOL = 1e-9

"""Angular slack (radians) allowed before a point is *outside* a node extent."""
const DEFAULT_CAP_ATOL = 1e-12

"""
Angular slack (radians) on the 90° geodesic-convexity threshold for a
`SphericalCap`.

Defined once and used by both places that ask the question, because they must
agree: [`covering_law_problems`](@ref) decides whether sampling boundary
*vertices* is enough on it, and [`test_hierarchical_system`](@ref) asserts
convexity on it. A radius inside the band would otherwise be simultaneously
"convex enough to skip densification" and "not convex", which reports a
conformance failure while quietly weakening the check that failure is about.
"""
const CONVEX_RADIUS_SLACK = 1e-9

"""Below this tangent-plane projection length a neighbour names no direction."""
const DEFAULT_PROJ_ATOL = 1e-9

"""Azimuths (radians) closer than this are treated as the same direction."""
const DEFAULT_ANGLE_ATOL = 1e-9

_default_rng() = Random.MersenneTwister(DEFAULT_SEED)

"""
    has_nonfallback_method(f, args...) -> Bool

Whether dispatch for `f(args...)` reaches a method written **for this kind of
grid or system** rather than the interface-wide generic — the harness's test
for "is this primitive implemented here?".

# Strategy: specificity, not provenance

The test is **specificity**: the matched method counts as an implementation
when one of its arguments is a type *strictly narrower* than `AbstractGrid` or
`AbstractHierarchicalGridSystem`. Whoever owns the module is not consulted.

`applicable` alone answered this question only before the fallback substrate
landed; since then every generic in the interface is applicable to every grid,
because the substrate supplies a method for
`AbstractGrid`/`AbstractHierarchicalGridSystem`. So applicability is uniformly
`true` and carries no information. The `applicable` call is kept as the cheap
first clause — it rules out a signature that genuinely has no method at all,
and it keeps `which` off the error path for that case.

This used to test **provenance**: `DGG.Fallbacks` as a sentinel module, a
method owned by it meaning "not implemented". That is right for the generics
and wrong for everything else the substrate ships, because `Fallbacks` also
defines the package's own concrete grid types and their specialised methods.
`PartialGrid` and `AuthalicGrid` were therefore reported as implementing
nothing at all, and `cellat`, `neighbors` and `ring` were skipped on them —
precisely the wrapped and subset paths most worth testing, and the ones with
no other harness coverage. Specificity restores them with no special-case
list: a method on `PartialGrid` is a method about `PartialGrid`, whichever
module it was typed in.

A consequence worth stating, unchanged by the move: this is a test of what the
method is *written for*, not of whether the answer is correct. A grid that
implements a primitive badly is tested and fails, which is the point; a grid
that leaves it to the interface-wide generic is skipped, because that generic
has its own tests and this harness is deliberately cursor-free.

# Deliberate non-guards

Two checks in [`test_hierarchical_system`](@ref) run **unconditionally**, and
that is not an oversight:

  - [`node_extent`](@ref) — the generic default is parameterised by the
    system's own [`cap_inflation`](@ref), so even a fully generic extent is the
    system's claim about its own geometry. The covering law is therefore always
    the system's to answer for.
  - [`cellposition`](@ref) — the `cellindex`/`cellposition` bijection is a law
    the grid owes regardless of which module computes it. A generic
    `cellposition` that disagrees with the grid's own `cellindex` is a real
    conformance failure, not an unimplemented method.
"""
function has_nonfallback_method(f, args...)
    applicable(f, args...) || return false
    m = try
        which(f, Base.typesof(args...))
    catch err
        # `which` throws `ArgumentError` on an ambiguous match. An ambiguity is
        # not a specificity answer either, so it is not evidence that anything
        # was implemented: report "not implemented" and let the caller skip
        # rather than let the harness die inside its own guard.
        (err isa ArgumentError || err isa MethodError) || rethrow()
        return false
    end
    return has_specific_subject(m)
end

"""
    has_specific_subject(m::Method) -> Bool

Whether `m` takes a grid or a system argument **strictly narrower** than
`AbstractGrid` / `AbstractHierarchicalGridSystem` — the specificity test behind
[`has_nonfallback_method`](@ref).

Each parameter of the signature is rewrapped in the method's own `where`
clauses before the comparison, so a parametric subject type
(`PartialGrid{S,V,ID,G}`, `AuthalicGrid{G}`) is compared as the closed type
`PartialGrid{S,V,ID,G} where {S,V,ID,G}` rather than as an open body whose free
type variables make subtyping undecidable.

"Strictly narrower" is spelled `T <: base && !(base <: T)` rather than
`T !== base` on purpose: a method written `f(g::G) where {G<:AbstractGrid}`
has a parameter that is *equal* to `AbstractGrid` without being `===` to it,
and that method is the interface-wide generic, not an implementation.
"""
function has_specific_subject(m::Method)
    sig = m.sig
    body = Base.unwrap_unionall(sig)
    body isa DataType || return false
    params = body.parameters
    # Parameter 1 is `typeof(f)`; the subject can be in any of the rest.
    for i in 2:length(params)
        T = Base.rewrap_unionall(params[i], sig)
        T isa Type || continue
        (_strictly_narrower(T, DGG.AbstractGrid) ||
         _strictly_narrower(T, DGG.AbstractHierarchicalGridSystem)) && return true
    end
    return false
end

_strictly_narrower(T, base) = T <: base && !(base <: T)

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

Vertices alone are a *sound* proxy for the whole boundary when the region they
are being tested against is geodesically convex (see the covering law), which
for a `SphericalCap` means an angular radius of at most 90°: a convex region
containing two points contains the arc between them, so containing every vertex
is containing every arc.

Against a non-convex region no finite sample is sound, and `densify` does not
make it so — it **strengthens** the check without closing it, because an arc can
still bulge outside the region between two interpolated samples that are both
inside. `densify` is what the harness reaches for when convexity does not hold,
and what it buys is a smaller escape, not a proof.
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
when the ancestor extents are geodesically convex (radius ≤ 90° up to
[`CONVEX_RADIUS_SLACK`](@ref), so the cap contains the great-circle arc between
any two points it contains). When a sampled extent is larger than that the
harness silently switches to `arc_samples` interpolated points per edge.

That switch **strengthens** the check; it does not restore soundness. Against a
cap wider than 90° no finite sample of a boundary is sound at all — an arc
between two sampled points that are both inside the cap can still bulge outside
it — so densification shrinks the escape rather than closing it. What keeps the
harness honest is that [`test_hierarchical_system`](@ref) additionally *asserts*
convexity: a system with wider extents has to opt out of that assertion
explicitly, so the weakening is declared rather than discovered. Pass `densify`
to force a point count either way.
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

    # Vertices suffice against convex (≤ 90°) caps; densify against anything
    # wider. The threshold is `CONVEX_RADIUS_SLACK` so that it is the *same*
    # question `test_hierarchical_system`'s convexity assertion asks — a cap in
    # a band between the two would be densified and reported, or asserted
    # convex and sampled at its vertices only.
    d = densify === nothing ?
        (all(a -> last(a).radius <= π / 2 + CONVEX_RADIUS_SLACK, ancestors) ? 0 : arc_samples) :
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
    neighbor_problems(grid, c; connectivity, sys, two_hop) -> Vector{String}

The [`neighbors`](@ref) laws at one cell: determinism across calls, exclusion of
`c` itself, distinctness, a container typed at the system's cell index type, a
count within [`max_neighbors`](@ref), membership in the grid, a common level,
and symmetry — `c′ ∈ neighbors(c)` implies `c ∈ neighbors(c′)`, checked in both
directions (see `two_hop` below).

This function is about *membership*; the ordering half of the contract lives in
[`neighbor_order_problems`](@ref). Nothing here requires or rewards an
id-sorted result — the order is rotational, and a harness that sorted would be
asserting the opposite of the contract.

Partial coverage is handled by the contract rather than by an exemption: a
neighbour beyond the grid's coverage must be **absent** from the result, so a
returned cell that has no [`cellposition`](@ref) is itself a violation, and
every cell that *is* returned is in the grid and can be asked for its own
neighbours in turn. Symmetry is therefore total over whatever the grid returns.

`two_hop` (default `true`) adds the closure described below; passing `false`
leaves only the one-directional sweep, which is what the harness's own tests use
to show that the closure catches strictly more.
"""
function neighbor_problems(grid, c; connectivity::Connectivity = Vertex(),
        sys = DGG.system(grid), two_hop::Bool = true)
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

    # ---- Two-hop symmetry closure -----------------------------------------
    #
    # The sweep above is one-directional: it walks out from `c` and asks whether
    # each cell that `c` names lists `c` back, so it only fires when the sampled
    # cell is the VICTIM of an omission. A cell that forgets one of its own
    # neighbours passes it — from the forgetful cell everything it still lists
    # is symmetric — and the violation is invisible unless the sampler happens
    # to draw the forgotten cell too. On a grid sampled at 8 cells out of
    # millions, it usually does not.
    #
    # The closure fixes that from `c` alone. A cell `x` that is truly adjacent
    # to `c` but missing from `neighbors(c)` cannot be far away: the
    # neighbourhood of `c` is a cycle, so `x` shares a boundary feature with at
    # least one of `c`'s SURVIVING neighbours and is therefore reachable in two
    # hops. Gathering `X = ⋃ neighbors(nb)` over `nb ∈ neighbors(c)` and asking,
    # for every `x ∈ X` that `c` does not list, whether `x` lists `c` is
    # therefore a COMPLETE test for a single omission at `c` — the exhaustive
    # both-directions sweep of `test/systems/H3/runtests.jl`, made affordable on
    # a sampled grid.
    if two_hop
        # Seeded with `c` and everything `c` already lists, so the loop only
        # ever visits the cells that are candidates for the omission, once each.
        checked = Set(ns)
        push!(checked, c)
        for nb in ns
            DGG.cellposition(grid, nb) === nothing && continue
            for x in DGG.neighbors(grid, nb, 1; connectivity)
                x in checked && continue
                push!(checked, x)
                DGG.cellposition(grid, x) === nothing && continue
                if c in collect(DGG.neighbors(grid, x, 1; connectivity))
                    push!(problems,
                        "neighbours are not symmetric: $c ∈ neighbors($x) but $x ∉ neighbors($c) " *
                        "(two-hop closure; $c is the cell that forgot a neighbour)")
                end
            end
        end
    end

    isempty(collect(DGG.neighbors(grid, c, 0; connectivity))) ||
        push!(problems, "neighbors($c, 0) is not empty")
    return problems
end

# ===========================================================================
# Rotational order
#
# `neighbors` and `ring` are ordered, and the order is part of the contract:
# each ring runs counter-clockwise seen from OUTSIDE the sphere, and the disc is
# its rings concatenated outward. These are the checks for that, and they are
# the reason nothing else in this harness ever sorts a neighbour list.
# ===========================================================================

"""
    tangent_frame(p) -> (east, north)

A right-handed orthonormal basis of the tangent plane at the unit vector `p`,
oriented so that `east × north == p`.

That orientation is the whole point: with the outward normal completing the
basis, an azimuth `atan(d ⋅ north, d ⋅ east)` **increases counter-clockwise seen
from outside the sphere**, which is the direction the [`neighbors`](@ref)
contract is written in.

The seed axis is the one `p` leans on least, so the Gram–Schmidt step stays
well-conditioned everywhere — including at the poles, where a frame built from
lon/lat east/north degenerates.
"""
function tangent_frame(p)
    ax = abs(p[1]) <= abs(p[2]) ?
        (abs(p[1]) <= abs(p[3]) ? (1.0, 0.0, 0.0) : (0.0, 0.0, 1.0)) :
        (abs(p[2]) <= abs(p[3]) ? (0.0, 1.0, 0.0) : (0.0, 0.0, 1.0))
    s = _dot(ax, p)
    e = (ax[1] - s * p[1], ax[2] - s * p[2], ax[3] - s * p[3])
    n = _norm(e)
    east = (e[1] / n, e[2] / n, e[3] / n)
    return east, _cross(p, east)   # p × east, so that east × north == p
end

"""
    tangent_offset(p, q) -> NTuple{3,Float64}

`q` projected onto the tangent plane at `p`: the direction from `p` to `q`, as a
vector. Its length shrinks to zero as `q` approaches `p` or its antipode, which
is exactly when `q` names no direction at all.
"""
function tangent_offset(p, q)
    r = _dot(q, p)
    return (q[1] - r * p[1], q[2] - r * p[2], q[3] - r * p[3])
end

"""
    winding_problems(grid, c, shell; label, proj_atol, ang_atol) -> Vector{String}

Whether `shell` — one ring of cells about `c`, in the order the grid returned
them — is a single **counter-clockwise cycle seen from outside the sphere**.
This is the rotational half of the [`neighbors`](@ref)/[`ring`](@ref) contract,
and the property that sorting by id destroys.

# Method

Project each member's [`cell_centroid`](@ref) onto the tangent plane at `c`'s
centroid and take its azimuth in a right-handed frame ([`tangent_frame`](@ref)).
A single counter-clockwise cycle is exactly a cyclic sequence of azimuths that
increases and **wraps around once**. A clockwise cycle wraps `n - 1` times; a
sequence merely sorted by id wraps some arbitrary number of times. Counting
wraparounds is invariant under where the ring starts, which is what lets this
check hold every system to the rotational law while leaving each of them free to
document its own ring-1 start.

The equivalent local statement — and the geometric definition the count is a
robust rendering of — is that for consecutive tangent offsets `a`, `b` about the
centroid `p`, counter-clockwise seen from outside means the signed volume
`dot(cross(a, b), p) > 0`, which is the sign of the wrapped azimuth increment.
Asserting that pair by pair is *not* robust: a ring with an angular gap wider
than 180° has one honest increment past `π` that the sign test reads as
negative. So the wraparound count carries the law, and the signed volume is used
to name the offending step in the report.

# Robustness

  - A member that projects to nothing (`proj_atol`) names no direction — a
    coarse grid really does have a ring member at the subject's antipode — and
    is dropped rather than given a meaningless azimuth.
  - Two consecutive members whose azimuths agree to within `ang_atol` are not
    counted as a wraparound in either direction.
  - A shell left with fewer than three usable directions is not checked: two
    points have no winding, and the contract's order is vacuous there.
"""
function winding_problems(grid, c, shell;
        label::AbstractString = "ring",
        proj_atol::Real = DEFAULT_PROJ_ATOL,
        ang_atol::Real = DEFAULT_ANGLE_ATOL)
    problems = String[]
    members = collect(shell)
    length(members) < 3 && return problems

    p = DGG.cell_centroid(grid, c)
    east, north = tangent_frame(p)
    azimuths = Float64[]
    offsets = NTuple{3,Float64}[]
    kept = similar(members, 0)
    for m in members
        d = tangent_offset(p, DGG.cell_centroid(grid, m))
        _norm(d) <= proj_atol && continue
        push!(azimuths, atan(_dot(d, north), _dot(d, east)))
        push!(offsets, d)
        push!(kept, m)
    end
    n = length(azimuths)
    n < 3 && return problems

    wraps = 0
    for i in 1:n
        j = i == n ? 1 : i + 1
        azimuths[j] < azimuths[i] - ang_atol && (wraps += 1)
    end
    wraps == 1 && return problems

    # Name a step that turns the wrong way, if there is one. There need not be:
    # a sequence that winds twice turns counter-clockwise at every step.
    culprit = ""
    for i in 1:n
        j = i == n ? 1 : i + 1
        v = _dot(_cross(offsets[i], offsets[j]), p)
        if v <= 0 && abs(azimuths[j] - azimuths[i]) > ang_atol
            culprit = "; the step $(kept[i]) → $(kept[j]) turns clockwise " *
                      "(signed volume $v ≤ 0)"
            break
        end
    end
    push!(problems,
        "$label about $c is not a single counter-clockwise cycle: its members' azimuths " *
        "about cell_centroid($c) wrap $wraps times, not once (a clockwise ring wraps " *
        "$(n - 1) times, an id-sorted one arbitrarily)$culprit")
    return problems
end

"""
    neighbor_order_problems(grid, c; connectivity, k, require_rotational_rings,
                            proj_atol, ang_atol) -> Vector{String}

The **order** laws of [`neighbors`](@ref)/[`ring`](@ref) at one cell, for every
`j` in `1:k`:

  - **concatenation** — `neighbors(grid, c, j) == vcat(ring(grid, c, 1), …,
    ring(grid, c, j))`, element for element. Not as sets: rings are never
    interleaved, and the result is never sorted by id.
  - **tail block** — `ring(grid, c, j)` is the last `length(ring)` elements of
    `neighbors(grid, c, j)`, element for element. Implied by the concatenation
    law, and checked separately because it is the form callers rely on and the
    form an independent-walk implementation breaks first.
  - **rotational winding** — each ring is a single counter-clockwise cycle seen
    from outside the sphere ([`winding_problems`](@ref)).

At `j = 1` the concatenation law collapses to an identity every implementation
gets right, so it is `j ≥ 2` that carries the weight there; conversely the
winding law is sharpest at `j = 1`, where the ring is the cell's own
neighbourhood and unambiguously a cycle.

Ring 1's winding is **always** checked. `require_rotational_rings` (default
`true`) gates rings `2:k`, because a system whose outer shells are genuinely
irregular — a lattice whose distance-`k` shell is not a simple cycle at all —
can conform to the contract's spirit without the outer rings reading as one
turn. Opting out is a documented claim about that system's shells, not a way to
silence a failure: ring 1 is not negotiable.
"""
function neighbor_order_problems(grid, c;
        connectivity::Connectivity = Vertex(),
        k::Integer = 2,
        require_rotational_rings::Bool = true,
        proj_atol::Real = DEFAULT_PROJ_ATOL,
        ang_atol::Real = DEFAULT_ANGLE_ATOL)
    problems = String[]
    rings = [collect(DGG.ring(grid, c, j; connectivity)) for j in 1:Int(k)]
    for j in 1:Int(k)
        disc = collect(DGG.neighbors(grid, c, j; connectivity))
        shell = rings[j]

        concatenated = reduce(vcat, view(rings, 1:j); init = similar(disc, 0))
        disc == concatenated ||
            push!(problems,
                "neighbors($c, $j) is not its rings concatenated outward: got $disc, " *
                "but vcat(ring 1..$j) is $concatenated")

        if length(disc) < length(shell)
            push!(problems,
                "neighbors($c, $j) has $(length(disc)) cells, fewer than ring($c, $j)'s " *
                "$(length(shell)); ring cannot be its tail block")
        else
            tail = disc[(length(disc) - length(shell) + 1):end]
            tail == shell ||
                push!(problems,
                    "ring($c, $j) is not the tail block of neighbors($c, $j): the tail is " *
                    "$tail, the ring is $shell")
        end

        if j == 1 || require_rotational_rings
            append!(problems, winding_problems(grid, c, shell;
                label = "ring $j", proj_atol, ang_atol))
        end
    end
    return problems
end

"""
    check_neighbor_order(grid; n_samples, rng, connectivity, k,
                         require_rotational_rings) -> Bool

Whether `grid` satisfies [`neighbor_order_problems`](@ref)' laws over sampled
cells — the boolean face, for a caller asserting that a deliberately misordered
implementation is caught.
"""
function check_neighbor_order(grid;
        n_samples::Integer = 16,
        rng::Random.AbstractRNG = _default_rng(),
        connectivity::Connectivity = Vertex(),
        k::Integer = 2,
        require_rotational_rings::Bool = true)
    _, cells = sample_cells(rng, grid, n_samples)
    return all(c -> isempty(neighbor_order_problems(grid, c;
        connectivity, k, require_rotational_rings)), cells)
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
    if !(r isa UnitRange{Int})
        # Reported, not thrown: the checks below index and difference `r`, so
        # carrying on with something that is not a range turns a conformance
        # report into an uncaught error from inside the harness.
        push!(problems, "descendant_range($c, $l) returned a $(typeof(r)), not a UnitRange{Int}")
        return problems
    end

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

!!! note "The `cellat` probe assumes star-shaped cells"
    `cellat(grid, cell_centroid(grid, c)) == c` is a law of the interface. The
    probe that goes with it is not: to catch a nearest-centroid lookup — which
    conforms at the centroid alone, and is wrong for every tessellation that is
    not a Voronoi one — the suite also walks `slerp(centroid, v, 0.5)` out
    towards each boundary vertex `v` and asserts those points belong to `c`.
    That is only guaranteed for a cell that is **star-shaped about its reported
    centroid**: the half-way point of the geodesic from the centroid to a vertex
    is inside the cell exactly when every such geodesic stays inside it.

    A cell that is not star-shaped about its centroid — a strongly re-entrant
    cell, or one whose centroid is a computed interior point rather than a
    kernel point — can fail this probe while satisfying `cellat` perfectly. If
    that is your system: **report it** rather than working around it. The probe
    is a proxy for "does `cellat` really locate, or does it just match
    centroids?", not a law of the interface, and it should be replaced for such
    a system by interior points drawn some other way (a triangulation's
    incentres, say). Do not move the centroid to make the probe pass; the
    centroid contract is a strictly interior point, not a kernel point, and
    changing it to satisfy a proxy trades a real law for a convenient one.

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
            if isempty(cells) ||
                    !has_nonfallback_method(DGG.cellat, grid,
                        DGG.cell_centroid(grid, first(cells)))
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
                             atol = 1e-12, unit_atol = 1e-9,
                             require_convex_extents = true,
                             neighbor_k = 2, require_rotational_rings = true,
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
  - **neighbours** — determinism, symmetry in both directions (including the
    two-hop closure of [`neighbor_problems`](@ref)), `k = 0` empty, ring 1 a
    single counter-clockwise cycle, `Edge()` a subset of `Vertex()`, all under
    both connectivities.
  - **`ring`** — `ring(grid, c, 0) == [c]`, shells disjoint, and the *order*
    laws of [`neighbor_order_problems`](@ref): the disc is its rings
    concatenated outward, each ring is the tail block of its disc, and each ring
    winds counter-clockwise.
  - **`descendant_range`** — when [`has_sorted_subtrees`](@ref) is `true`, the
    positions of the actual descendants exactly fill the returned range.

`levels` selects the levels to test; when it has more than `n_levels` entries a
seeded sample of that many is used, always including the coarsest and deepest.
Derived methods that a system has not implemented (`ancestor`, `descendants`,
[`neighbors`](@ref), [`ring`](@ref)) are skipped rather than failed.

# The keyword arguments that carry semantics

  - `atol` — angular slack, in radians, before a boundary point counts as
    *outside* an ancestor's [`node_extent`](@ref); forwarded to
    [`covering_law_problems`](@ref). A system whose subtree caps are exact
    rather than inflated puts cell corners *on* the cap rim, where the covering
    test is a floating-point coin toss — HEALPix is the case in this package —
    and such a system must loosen this rather than inflate its extents to hide
    the rounding.
  - `unit_atol` — how far a [`node_extent`](@ref) centre may sit off the unit
    sphere; forwarded to [`node_extent_problems`](@ref).
  - `require_convex_extents` — whether to *assert* that every sampled node
    extent has radius ≤ 90° (up to [`CONVEX_RADIUS_SLACK`](@ref)), which is what
    makes sampling a boundary's vertices a sound proxy for the whole boundary.
    `node_extent` explicitly permits wider extents, so a system that means it
    passes `false`; the covering check then densifies, which strengthens it
    without making it sound (see [`covering_law_problems`](@ref)). Opting out is
    a claim about the system, and it belongs in that system's documentation.
  - `neighbor_k` — how many rings out to test the order and shell laws. The
    default of 2 is the smallest value that tests anything: at `k = 1` the
    concatenation law is an identity.
  - `require_rotational_rings` — whether rings `2:neighbor_k` must each read as
    a single counter-clockwise cycle. Ring 1 is checked unconditionally. A
    system whose outer shells are genuinely irregular can pass `false`, with the
    reason documented on its side; it does not exempt ring 1.

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
        atol::Real = DEFAULT_CAP_ATOL,
        unit_atol::Real = DEFAULT_UNIT_ATOL,
        require_convex_extents::Bool = true,
        neighbor_k::Integer = 2,
        require_rotational_rings::Bool = true,
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
            if isempty(probe) || !has_nonfallback_method(DGG.ancestor, sys, first(probe), 0)
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
            if isempty(probe) || !has_nonfallback_method(DGG.descendants, sys, first(probe), 0)
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
                @test node_extent_problems(DGG.node_extent(sys, c); unit_atol) == String[]
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
                    # The same threshold the covering check's densify decision
                    # uses, to the bit: see `CONVEX_RADIUS_SLACK`.
                    @test DGG.node_extent(sys, c).radius <= π / 2 + CONVEX_RADIUS_SLACK
                end
            end
        end

        @testset "covering law" begin
            # The same sampled cells every other law used, so a covering failure
            # is attributable to a cell the rest of the report also names.
            @test covering_law_problems(sys; levels = tested, rng,
                cells = Dict(l => last(samples[l]) for l in tested),
                descent_depth, branch_samples, atol) == String[]
        end

        @testset "neighbors" begin
            probe = last(samples[first(tested)])
            if isempty(probe) ||
                    !has_nonfallback_method(DGG.neighbors, grids[first(tested)], first(probe), 1)
                @test_skip "neighbors is not implemented for this system"
            else
                for conn in connectivities, l in tested
                    grid = grids[l]
                    for c in last(samples[l])
                        @test neighbor_problems(grid, c; connectivity = conn, sys) == String[]
                        # Ring 1 is `neighbors(c, 1)`, so its winding is testable
                        # here without `ring` — and it is the one ring whose
                        # rotational order is never negotiable.
                        @test winding_problems(grid, c,
                            DGG.neighbors(grid, c, 1; connectivity = conn);
                            label = "neighbors(c, 1)") == String[]
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
            if isempty(probe) ||
                    !has_nonfallback_method(DGG.ring, grids[first(tested)], first(probe), 0)
                @test_skip "ring is not implemented for this system"
            else
                # The order laws relate `ring` to `neighbors`, so they are the
                # composite's to answer for only when the system owns both.
                ordered = has_nonfallback_method(DGG.neighbors, grids[first(tested)],
                    first(probe), 1)
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
                        # ...and the same relation as SEQUENCES, which is the
                        # actual contract: the disc is its rings concatenated
                        # outward, each ring is the tail block of its disc, and
                        # each ring winds counter-clockwise. The set laws above
                        # are what an implementation that sorts by id still
                        # satisfies.
                        if ordered
                            @test neighbor_order_problems(grid, c; connectivity = conn,
                                k = neighbor_k, require_rotational_rings) == String[]
                        end
                    end
                end
                ordered || @test_skip "neighbors is a fallback; the ring/disc order laws are its own"
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

end # module DiscreteGlobalGridsConformanceTesting
