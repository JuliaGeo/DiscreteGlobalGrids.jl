# Conformance checks for the public grid and hierarchical-system interfaces.
# Each law has a problem collector, a labelled `@testset` wrapper, and a
# boolean `check_*` wrapper. Seeded sampling makes failures reproducible, and
# all laws in one run examine the same sampled cells.

module DiscreteGlobalGridsConformanceTesting

using Test
import Random

import GeometryOps as GO
const USph = GO.UnitSpherical

import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, Connectivity, Vertex, Edge

export test_grid_interface, test_hierarchical_system, test_generic_fallbacks

"""Seed used for reproducible default sampling."""
const DEFAULT_SEED = 20260813

"""Maximum Euclidean deviation of a boundary vertex from the unit sphere."""
const DEFAULT_UNIT_ATOL = 1e-9

"""Angular slack (radians) allowed before a point is *outside* a node extent."""
const DEFAULT_CAP_ATOL = 1e-12

"""
Angular slack (radians) on the 90° geodesic-convexity threshold for a
`SphericalCap`.

[`covering_law_problems`](@ref) and [`test_hierarchical_system`](@ref) use this
same tolerance when deciding whether a cap is geodesically convex.
"""
const CONVEX_RADIUS_SLACK = 1e-9

"""Below this tangent-plane projection length a neighbour names no direction."""
const DEFAULT_PROJ_ATOL = 1e-9

"""Azimuths (radians) closer than this are treated as the same direction."""
const DEFAULT_ANGLE_ATOL = 1e-9

_default_rng() = Random.MersenneTwister(DEFAULT_SEED)

"""The subjects a specialized grid or system method dispatches on."""
const DISPATCH_SUBJECTS = (DGG.AbstractGrid, DGG.AbstractHierarchicalGridSystem)

"""
    has_nonfallback_method(f, args...) -> Bool

Whether dispatch for `f(args...)` reaches a method specialized for a grid or
system type rather than an interface-wide generic.

# Method selection

The matched method counts as specialized when any grid or system parameter is
strictly narrower than `AbstractGrid` or
`AbstractHierarchicalGridSystem`. Module ownership is irrelevant because a
module may define both generic fallbacks and specialized methods.

This function reports only whether a specialized method exists. The
conformance checks separately validate the method's result. `node_extent` and
`cellposition` are checked unconditionally because their contracts also apply
when an interface-wide implementation provides them.
"""
function has_nonfallback_method(f, args...)
    applicable(f, args...) || return false
    m = try
        which(f, Base.typesof(args...))
    catch err
        # Ambiguous dispatch does not identify a specialized implementation.
        (err isa ArgumentError || err isa MethodError) || rethrow()
        return false
    end
    return has_specific_subject(m)
end

"""
    has_specific_subject(m::Method, bases = (AbstractGrid, AbstractHierarchicalGridSystem)) -> Bool

Whether `m` takes an argument **strictly narrower** than one of `bases` — the
specificity test behind [`has_nonfallback_method`](@ref), whose subjects are a
grid and a system.

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
function has_specific_subject(m::Method, bases = DISPATCH_SUBJECTS)
    sig = m.sig
    body = Base.unwrap_unionall(sig)
    body isa DataType || return false
    params = body.parameters
    # Parameter 1 is `typeof(f)`; the subject can be in any of the rest.
    for i in 2:length(params)
        T = Base.rewrap_unionall(params[i], sig)
        T isa Type || continue
        any(base -> _strictly_narrower(T, base), bases) && return true
    end
    return false
end

_strictly_narrower(T, base) = T <: base && !(base <: T)

"""
    dispatches_generically(f, args...) -> Bool

Whether `f(args...)` is answered by an interface-wide method: one specialized
for no grid, system **or cell index** type.

Stricter than `!has_nonfallback_method(f, args...)`, which is a statement about
gate labels and reads a method dispatching on a system's cell id type alone as
generic. This is the assertion [`test_generic_fallbacks`](@ref) makes before it
runs a law, so "the fallback was tested" is checked rather than assumed.
"""
function dispatches_generically(f, args...)
    applicable(f, args...) || return false
    m = try
        which(f, Base.typesof(args...))
    catch err
        (err isa ArgumentError || err isa MethodError) || rethrow()
        return false
    end
    return !has_specific_subject(m, (DISPATCH_SUBJECTS..., DGG.AbstractCellIndex))
end

# ===========================================================================
# Small vector helpers
#
# Local three-component operations avoid adding LinearAlgebra as a dependency.
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
            # The sampled branches are depth-limited. One additional chain
            # reaches `max_level` so deeper extent failures remain observable.
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

    # Vertices suffice for convex caps with radius at most 90°. Wider caps use
    # densified arcs. The system-level convexity check uses the same threshold.
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

`require_nonempty` (default `false`) adds the one law here that is not
relational: on a grid of more than one cell, a cell has at least one neighbour.
It is off by default because a system answering for its own cells knows its own
degrees. [`test_generic_fallbacks`](@ref) turns it on, because an empty answer
is exactly how the geometric fallback fails — a tessellation whose neighbouring
cells do not compute identical corner coordinates has no adjacency at all under
it, and every relational law above then passes vacuously.
"""
function neighbor_problems(grid, c; connectivity::Connectivity = Vertex(),
        sys = DGG.system(grid), two_hop::Bool = true,
        require_nonempty::Bool = false)
    problems = String[]
    ns = collect(DGG.neighbors(grid, c, 1; connectivity))

    if require_nonempty && isempty(ns) && DGG.ncells(grid) > 1
        push!(problems,
            "neighbors($c) is empty on a grid of $(DGG.ncells(grid)) cells; every cell " *
            "of a tessellation has at least one $(nameof(typeof(connectivity)))-neighbour")
    end

    collect(DGG.neighbors(grid, c, 1; connectivity)) == ns ||
        push!(problems, "neighbors($c) is not deterministic across calls")
    c in ns && push!(problems, "neighbors($c) contains $c itself")
    allunique(ns) || push!(problems, "neighbors($c) contains duplicates")

    if sys !== nothing
        T = DGG.cellindextype(sys)
        eltype(ns) === T ||
            push!(problems, "neighbors($c) has eltype $(eltype(ns)), not the system's $T")
        bound = DGG.max_neighbors(sys, connectivity)
        # `nothing` means the system declares no static bound, so there is no
        # ceiling to check against.
        bound === nothing || length(ns) <= bound ||
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

    # A direct symmetry check only detects omissions when the omitted cell is
    # sampled. The two-hop closure also examines cells reached through each
    # reported neighbour, exposing a cell that lists `c` when `c` omits it.
    if two_hop
        # Skip the subject and its reported neighbours; only omission
        # candidates need a reverse-adjacency check.
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
# Each ring runs counter-clockwise as seen from outside the sphere, and a
# neighbourhood is its rings concatenated outward. Sorting would destroy this
# rotational order.
# ===========================================================================

"""
    tangent_frame(p) -> (east, north)

A right-handed orthonormal basis of the tangent plane at the unit vector `p`,
oriented so that `east × north == p`.

With the outward normal completing the basis, an azimuth
`atan(d ⋅ north, d ⋅ east)` **increases counter-clockwise seen
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

# Azimuths of `points` about `p`, in one right-handed frame, with the ones that
# name no direction dropped. Returns the azimuths, their tangent offsets, and
# which inputs survived.
function _azimuths_about(p, points, proj_atol::Real)
    east, north = tangent_frame(p)
    azimuths = Float64[]
    offsets = NTuple{3,Float64}[]
    kept = Int[]
    for (i, q) in enumerate(points)
        d = tangent_offset(p, q)
        _norm(d) <= proj_atol && continue
        push!(azimuths, atan(_dot(d, north), _dot(d, east)))
        push!(offsets, d)
        push!(kept, i)
    end
    return azimuths, offsets, kept
end

# How many times a cyclic azimuth sequence steps backwards. One for a single
# counter-clockwise turn, `n - 1` for a clockwise one.
function _wrap_count(azimuths::Vector{Float64}, ang_atol::Real)
    n = length(azimuths)
    wraps = 0
    for i in 1:n
        j = i == n ? 1 : i + 1
        azimuths[j] < azimuths[i] - ang_atol && (wraps += 1)
    end
    return wraps
end

"""
    winding_problems(grid, c, shell; label, proj_atol, ang_atol) -> Vector{String}

Whether `shell` — one ring of cells about `c`, in the order the grid returned
them — is a single **counter-clockwise cycle seen from outside the sphere**.
This is the rotational half of the [`neighbors`](@ref)/[`ring`](@ref) contract,
and the property that sorting by id destroys.

**Which counter-clockwise.** The package has one handedness and
[`cell_boundary`](@ref) fixes it: a boundary ring winds counter-clockwise seen
from outside, and a neighbour ring turns the same way. So the check measures
both in the SAME frame — `c`'s own boundary vertices and the shell's centroids —
and holds them to the same wrap count. A shell that turns against its own cell's
boundary is reported here rather than left to the reader to reconcile against
[`boundary_problems`](@ref).

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
  - The boundary comparison is skipped when the boundary itself does not wrap
    once about the centroid, which is a defect [`boundary_problems`](@ref) owns.
"""
function winding_problems(grid, c, shell;
        label::AbstractString = "ring",
        proj_atol::Real = DEFAULT_PROJ_ATOL,
        ang_atol::Real = DEFAULT_ANGLE_ATOL)
    problems = String[]
    members = collect(shell)
    length(members) < 3 && return problems

    p = DGG.cell_centroid(grid, c)
    azimuths, offsets, keptix = _azimuths_about(p,
        (DGG.cell_centroid(grid, m) for m in members), proj_atol)
    n = length(azimuths)
    n < 3 && return problems
    kept = members[keptix]

    wraps = _wrap_count(azimuths, ang_atol)
    if wraps != 1
        # Name a step that turns the wrong way, if there is one. There need not
        # be: a sequence that winds twice turns counter-clockwise at every step.
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

    # The sense the ring just passed is the sense the cell's own boundary winds
    # in, measured the same way.
    bnd, _, _ = _azimuths_about(p, DGG.cell_boundary(grid, c), proj_atol)
    length(bnd) < 3 && return problems
    bwraps = _wrap_count(bnd, ang_atol)
    bwraps == length(bnd) - 1 &&
        push!(problems,
            "$label about $c turns the opposite way from cell_boundary($c): measured in " *
            "one frame about cell_centroid($c), the ring is a counter-clockwise cycle " *
            "and the boundary is a clockwise one. Both must wind counter-clockwise seen " *
            "from outside the sphere")
    return problems
end

"""
    ring_cycle_problems(grid, c, shell) -> Vector{String}

Whether `shell` walks *around* `c` rather than jumping about it: consecutive
members, cyclically, must be neighbours of one another. A complete one-ring
closes, so it has no break; a ring clipped by coverage is an arc, so it has one.
Anything with two or more breaks is not a rotational order at all — an id-sorted
ring of degree `n` typically has `n` of them.

This is the adjacency half of the order contract, and it catches what winding
alone cannot: a sequence whose centroids happen to sweep counter-clockwise while
stepping across the cell rather than around it.

Vertex connectivity only. Under [`Edge()`](@ref Edge) the ring members are
generally not edge-adjacent to one another — the four edge-neighbours of a
lattice cell touch only at corners — so the law does not apply.
"""
function ring_cycle_problems(grid, c, shell)
    problems = String[]
    members = collect(shell)
    length(members) < 3 && return problems
    breaks = 0
    culprit = ""
    for i in eachindex(members)
        j = i == lastindex(members) ? firstindex(members) : i + 1
        a, b = members[i], members[j]
        if !(b in DGG.neighbors(grid, a, 1; connectivity = Vertex()))
            breaks += 1
            isempty(culprit) && (culprit = "$a → $b")
        end
    end
    breaks <= 1 && return problems
    push!(problems,
        "ring 1 about $c is not a walk around the cell: $breaks of its " *
        "$(length(members)) consecutive pairs are not neighbours of each other " *
        "(first at $culprit). A closed one-ring has none and a coverage-clipped " *
        "arc has one; more means the members are in some order other than " *
        "rotational")
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
    from outside the sphere, turning the same way the cell's own
    [`cell_boundary`](@ref) winds ([`winding_problems`](@ref)).
  - **adjacency** — ring 1 walks around `c`: consecutive members are neighbours
    of one another, with at most the one break a coverage-clipped arc has
    ([`ring_cycle_problems`](@ref), `Vertex()` only).
  - **idiom agreement** — the position forms answer the id forms read through
    [`cellposition`](@ref), element for element, so the two idioms cannot
    present the same neighbourhood in two orders.

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

None of these laws says where a ring *starts*: the contract guarantees the
direction everywhere and leaves the phase to each system, so every check here is
invariant under rotating a ring. A system that wants its start pinned owes its
own oracle vector.
"""
function neighbor_order_problems(grid, c;
        connectivity::Connectivity = Vertex(),
        k::Integer = 2,
        require_rotational_rings::Bool = true,
        proj_atol::Real = DEFAULT_PROJ_ATOL,
        ang_atol::Real = DEFAULT_ANGLE_ATOL)
    problems = String[]
    rings = [collect(DGG.ring(grid, c, j; connectivity)) for j in 1:Int(k)]
    connectivity isa Vertex && !isempty(rings) &&
        append!(problems, ring_cycle_problems(grid, c, first(rings)))
    append!(problems, _position_form_problems(grid, c, rings; connectivity))
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

# The position form is the id form read through `cellposition`, element for
# element. Written once, here, so no system can present one neighbourhood in two
# orders — the drift that makes an oriented stencil silently change meaning when
# a caller moves from ids to indices.
function _position_form_problems(grid, c, rings; connectivity::Connectivity)
    problems = String[]
    p = DGG.cellposition(grid, c)
    p isa Int || return problems
    for (j, shell) in enumerate(rings)
        want = Int[]
        for x in shell
            q = DGG.cellposition(grid, x)
            q isa Int && push!(want, q)
        end
        got = DGG.ring(grid, p, j; connectivity)
        collect(got) == want || push!(problems,
            "ring(grid, $p, $j) is not ring(grid, $c, $j) read through cellposition: " *
            "got $got, the ids map to $want. The position and id forms are one order")
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
# Skip reporting
# ===========================================================================

"""
    skip!(skips, reason) -> nothing

Record one skipped law group as a `Broken` result and keep its `reason` for the
end-of-run report. `Test` never displays a `@test_skip` message, and a summary
cannot tell a skip from a known failure, so an unreported skip is a law the
caller believes was checked.
"""
function skip!(skips::Vector{String}, reason::AbstractString)
    push!(skips, reason)
    @test_skip reason
    return nothing
end

"""
    unimplemented_skip(what, subject, signature) -> String

The reason string for a law group that `subject` does not open: what was
skipped, the method whose absence skipped it, and where the law is checked
instead.
"""
unimplemented_skip(what::AbstractString, subject, signature::AbstractString) =
    "skipped: $what — $(subject_name(subject)) does not implement `$signature`. " *
    "A generic fallback answers it; its laws run under test_generic_fallbacks."

"""The name a message calls a system or grid by; a wrapper names what it wraps."""
subject_name(x) = string(nameof(typeof(x)))

# The `@info` is the only place a skip states itself. It goes inside the test
# set so it is emitted whether or not a later law fails.
function report_skips(skips::Vector{String}, label::AbstractString)
    isempty(skips) && return nothing
    @info "conformance: $(length(skips)) law group(s) skipped for $label\n" *
          join(("  " * s for s in skips), "\n")
    return nothing
end

# ===========================================================================
# Forcing the generic fallbacks
# ===========================================================================

"""
    GenericFallbackSystem(sys) <: AbstractHierarchicalGridSystem

`sys` with its optional fast paths hidden: the same cells, ids, geometry,
hierarchy and traits, presented to dispatch as a system type that no
specialization was written for.

The required surface is forwarded through `sys`'s own
[`levelgrid`](@ref) — so a system that answers the grid contract with a grid
type of its own is wrapped as faithfully as one that implements the five
system-level primitives — and so are the traits and the two covering methods
([`node_extent`](@ref), [`cap_inflation`](@ref)) and
[`descendant_range`](@ref). Nothing else is: [`cellat`](@ref),
[`neighbors`](@ref), [`ring`](@ref), `one_ring`, [`ancestor`](@ref) and
[`descendants`](@ref) are left to the interface-wide implementations, and the
traversal engines to the generic ones.

Hiding by *type* rather than by `invoke` is what makes this total. A system's
fast paths hang on `HierarchicalLevelGrid{TheSystem}` or on `TheSystem`, and a
wrapped system appears as `HierarchicalLevelGrid{GenericFallbackSystem{TheSystem}}`
— outside every one of those signatures, including ones written for a family
supertype, and including any this package has not seen. There is no argument
form that reaches a specialization from here, so a fallback law cannot silently
be answered by the specialized path.

[`test_generic_fallbacks`](@ref) additionally asserts
[`dispatches_generically`](@ref) for each verb before running its laws.
"""
struct GenericFallbackSystem{S<:AbstractHierarchicalGridSystem} <: AbstractHierarchicalGridSystem
    system::S
end

# The wrapped system's own complete-level grid: the entry point every consumer
# uses, so both implementor styles forward identically.
_inner(w::GenericFallbackSystem, l::Integer) = DGG.levelgrid(w.system, Int(l))

DGG.cellindextype(w::GenericFallbackSystem) = DGG.cellindextype(w.system)
DGG.levels(w::GenericFallbackSystem) = DGG.levels(w.system)
DGG.rootcells(w::GenericFallbackSystem) = DGG.rootcells(w.system)
Base.parent(w::GenericFallbackSystem, c::AbstractCellIndex) = Base.parent(w.system, c)
DGG.children(w::GenericFallbackSystem, c::AbstractCellIndex) = DGG.children(w.system, c)

DGG.node_extent(w::GenericFallbackSystem, c::AbstractCellIndex) =
    DGG.node_extent(w.system, c)
DGG.cap_inflation(w::GenericFallbackSystem) = DGG.cap_inflation(w.system)
DGG.has_sorted_subtrees(w::GenericFallbackSystem) = DGG.has_sorted_subtrees(w.system)
DGG.descendant_range(w::GenericFallbackSystem, c::AbstractCellIndex, l::Integer) =
    DGG.descendant_range(w.system, c, l)
DGG.max_neighbors(w::GenericFallbackSystem, conn::Connectivity) =
    DGG.max_neighbors(w.system, conn)

DGG.ncells(w::GenericFallbackSystem, l::Integer) = DGG.ncells(_inner(w, l))
DGG.cellindex(w::GenericFallbackSystem, l::Integer, i::Int) = DGG.cellindex(_inner(w, l), i)
DGG.cellposition(w::GenericFallbackSystem, c::AbstractCellIndex) =
    DGG.cellposition(_inner(w, DGG.level(c)), c)
DGG.cell_boundary(w::GenericFallbackSystem, c::AbstractCellIndex) =
    DGG.cell_boundary(_inner(w, DGG.level(c)), c)
DGG.cell_centroid(w::GenericFallbackSystem, c::AbstractCellIndex) =
    DGG.cell_centroid(_inner(w, DGG.level(c)), c)

Base.show(io::IO, w::GenericFallbackSystem) =
    print(io, "GenericFallbackSystem(", w.system, ")")

subject_name(w::GenericFallbackSystem) = subject_name(w.system)

# The coarsest five levels: enough for the geometric fallbacks, and above the
# cell size at which the harness's own area and closure tolerances give out.
function _coarse_levels(sys)
    ls = DGG.levels(sys)
    return first(ls):min(last(ls), first(ls) + 4)
end

# ===========================================================================
# The two entry points
# ===========================================================================

"""
    test_grid_interface(grid; n_samples = 32, rng = MersenneTwister(20260813),
                        unit_atol = 1e-9, fallback_laws = false,
                        label = <grid type>)

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
[`cellat`](@ref) reports that test set as skipped and conforms regardless. Every
skip states its reason in an `@info` at the end of the run, because a `Broken`
count alone cannot be told from a known failure.

# Keyword arguments

  - `n_samples` — how many distinct positions to draw (see above).
  - `rng` — the sampling generator; a fresh seeded one per call by default, so
    two identical calls examine identical cells.
  - `unit_atol` — how far a boundary vertex or a centroid may sit off the unit
    sphere; forwarded to [`boundary_problems`](@ref) and
    [`centroid_problems`](@ref).
  - `fallback_laws` — whether to run the [`cellat`](@ref) law group even when a
    generic fallback is what answers. `false` (the default) tests what this grid
    implements; `true` tests whatever dispatch selects. Forcing the fallback
    *path* as well is [`test_generic_fallbacks`](@ref)' job — set here alone,
    this only stops the gate from skipping.
  - `label` — the name in the test-set header and in the skip report.

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
        fallback_laws::Bool = false,
        label::AbstractString = string(nameof(typeof(grid))))
    n = DGG.ncells(grid)
    positions, cells = sample_cells(rng, grid, n_samples)
    sys = DGG.system(grid)
    skips = String[]

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
            # System provenance and level are either both present or both absent.
            @test (sys === nothing) == (l === nothing)
            if sys !== nothing
                @test l in DGG.levels(sys)
                @test all(c -> DGG.level(c) == l, cells)
                @test all(c -> c isa DGG.cellindextype(sys), cells)
                # System grids use canonical cell-id order; standalone grids
                # may define a different dense order.
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
                skip!(skips, "skipped: cellposition of a cell outside the grid — no such " *
                    "cell is constructible from the interface alone (no system provenance, " *
                    "or the grid is a complete deepest level)")
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
            if isempty(cells)
                skip!(skips, "skipped: cellat(cell_centroid(c)) == c — the grid has no cells")
            elseif !(fallback_laws || has_nonfallback_method(DGG.cellat, grid,
                        DGG.cell_centroid(grid, first(cells))))
                skip!(skips, unimplemented_skip("cellat(cell_centroid(c)) == c",
                    sys === nothing ? grid : sys,
                    "cellat(::$(typeof(grid)), ::UnitSphericalPoint)"))
            else
                for c in cells
                    centroid = DGG.cell_centroid(grid, c)
                    @test DGG.cellat(grid, centroid) == c
                    # Midpoints between the centroid and each vertex verify
                    # interior lookup away from the centroid itself.
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

        report_skips(skips, label)
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
                             fallback_laws = false, label = <system type>)

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

Derived methods that a system has not implemented (`ancestor`, `descendants`,
[`neighbors`](@ref), [`ring`](@ref)) are skipped rather than failed — the
generic fallbacks answer them, and [`test_generic_fallbacks`](@ref) is what
holds those answers to the same laws. Every skip states its reason in an
`@info` at the end of the run.

# Keyword arguments

  - `levels` — which levels to test. Defaults to all of them.
  - `n_levels` — how many of `levels` to keep when it holds more: a seeded
    sample of that many, always including the coarsest and the deepest, where
    level-dependent arithmetic breaks.
  - `n_samples` — how many cells to draw per tested level. Every law in one run
    examines the same drawn cells, so a covering failure names a cell the
    hierarchy laws saw too.
  - `rng` — the sampling generator; a fresh seeded one per call by default, so
    two identical calls examine identical cells and levels. An RNG advances as
    it is drawn from: reproducing a run means a *new* generator on the same
    seed.
  - `descent_depth`, `branch_samples` — the shape of the covering law's subtree
    walk: how many levels down, and how many children followed at each step.
    Forwarded to [`covering_law_problems`](@ref), which also follows one chain
    to [`max_level`](@ref) regardless of the depth.
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
  - `connectivities` — the adjacencies to test the neighbour and ring laws
    under. Both of them by default; a system with one meaningful adjacency
    passes just that one, and the `Edge() ⊆ Vertex()` law then does not apply.
  - `fallback_laws` — whether to run the gated law groups (`ancestor`,
    `descendants`, [`neighbors`](@ref), [`ring`](@ref) and the order laws) even
    when a generic fallback is what answers. `false` (the default) tests what
    this system implements. Forcing the fallback *path* as well is
    [`test_generic_fallbacks`](@ref)' job — set here alone, this only stops the
    gates from skipping.
  - `label` — the name in the test-set header and in the skip report.

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
        fallback_laws::Bool = false,
        label::AbstractString = string(nameof(typeof(sys))))
    levelrange = DGG.levels(sys)
    tested = sample_levels(rng, levels, n_levels)
    maxl = DGG.max_level(sys)
    grids = Dict{Int,Any}(l => DGG.levelgrid(sys, l) for l in tested)
    samples = Dict{Int,Any}(l => sample_cells(rng, grids[l], n_samples) for l in tested)
    skips = String[]

    @testset "hierarchical system: $label" begin
        @testset "levels and traits" begin
            @test levelrange isa AbstractUnitRange{Int}
            @test !isempty(levelrange)
            @test maxl == last(levelrange)
            @test DGG.cellindextype(sys) <: AbstractCellIndex
            @test DGG.cap_inflation(sys) >= 1
            @test DGG.has_sorted_subtrees(sys) isa Bool
            # The trait obliges the method. Declaring one without the other is
            # caught here rather than deep inside a pruned traversal, and the
            # interface-wide method that reports it does not count as one.
            @test !DGG.has_sorted_subtrees(sys) ||
                  has_nonfallback_method(DGG.descendant_range, sys,
                      first(DGG.rootcells(sys)), first(levelrange))
            # A declared bound must be a usable capacity; `nothing` is the
            # legal "no static bound declared" answer, which costs the system
            # the fixed-capacity neighbour buffers and nothing else — said out
            # loud below rather than passed in silence.
            for conn in connectivities
                mn = DGG.max_neighbors(sys, conn)
                @test mn === nothing || mn >= 1
                mn === nothing && push!(skips,
                    "note: max_neighbors — $label declares no bound for $conn; subset neighbour buffers will be heap Vectors")
            end
            @test DGG.max_neighbors(sys) == DGG.max_neighbors(sys, Vertex())
            if Vertex() in connectivities && Edge() in connectivities
                # Edge() is a restriction of Vertex(), so its bound cannot be larger.
                mv = DGG.max_neighbors(sys, Vertex())
                me = DGG.max_neighbors(sys, Edge())
                @test mv === nothing || me === nothing || me <= mv
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
            if isempty(probe) || !(fallback_laws ||
                    has_nonfallback_method(DGG.ancestor, sys, first(probe), 0))
                skip!(skips, unimplemented_skip("ancestor", sys,
                    "ancestor(::$(subject_name(sys)), ::$(DGG.cellindextype(sys)), ::Integer)"))
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
            if isempty(probe) || !(fallback_laws ||
                    has_nonfallback_method(DGG.descendants, sys, first(probe), 0))
                skip!(skips, unimplemented_skip("descendants", sys,
                    "descendants(::$(subject_name(sys)), ::$(DGG.cellindextype(sys)), ::Integer)"))
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
            # For caps with radius at most 90°, containing a boundary's vertices
            # also contains its arcs. Wider caps require densified boundary checks.
            if !require_convex_extents
                skip!(skips, "skipped: node_extent convexity — require_convex_extents = false. " *
                    "The covering check densifies boundaries instead, which strengthens it " *
                    "without making vertex sampling sound")
            else
                for l in tested, c in last(samples[l])
                    # Match the threshold used to select boundary densification.
                    @test DGG.node_extent(sys, c).radius <= π / 2 + CONVEX_RADIUS_SLACK
                end
            end
        end

        @testset "covering law" begin
            # Reuse the samples from the other laws for consistent reports.
            @test covering_law_problems(sys; levels = tested, rng,
                cells = Dict(l => last(samples[l]) for l in tested),
                descent_depth, branch_samples, atol) == String[]
        end

        @testset "neighbors" begin
            probe = last(samples[first(tested)])
            if isempty(probe) || !(fallback_laws ||
                    has_nonfallback_method(DGG.neighbors, grids[first(tested)], first(probe), 1))
                skip!(skips, unimplemented_skip("neighbors", sys,
                    "neighbors(::$(typeof(grids[first(tested)])), " *
                    "::$(DGG.cellindextype(sys)), ::Integer)"))
            else
                for conn in connectivities, l in tested
                    grid = grids[l]
                    for c in last(samples[l])
                        @test neighbor_problems(grid, c; connectivity = conn, sys,
                            require_nonempty = fallback_laws) == String[]
                        # Ring 1 is `neighbors(c, 1)`, so its winding can be
                        # checked even when `ring` is not implemented.
                        @test winding_problems(grid, c,
                            DGG.neighbors(grid, c, 1; connectivity = conn);
                            label = "neighbors(c, 1)") == String[]
                    end
                end
                # Edge neighbours are a subset of vertex neighbours.
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
            if isempty(probe) || !(fallback_laws ||
                    has_nonfallback_method(DGG.ring, grids[first(tested)], first(probe), 0))
                skip!(skips, unimplemented_skip("ring", sys,
                    "ring(::$(typeof(grids[first(tested)])), " *
                    "::$(DGG.cellindextype(sys)), ::Integer)"))
            else
                # The order laws relate `ring` to `neighbors`, so they need both
                # to come from the same implementation — either both specialized
                # or, in fallback mode, both generic.
                ordered = fallback_laws || has_nonfallback_method(DGG.neighbors,
                    grids[first(tested)], first(probe), 1)
                for conn in connectivities, l in tested
                    grid = grids[l]
                    for c in last(samples[l])
                        @test collect(DGG.ring(grid, c, 0; connectivity = conn)) == [c]
                        # `neighbors(c, k)` is the union of disjoint shells 1:k.
                        shells = Set{DGG.cellindextype(sys)}()
                        for k in 1:neighbor_k
                            shell = Set(DGG.ring(grid, c, k; connectivity = conn))
                            @test isempty(intersect(shell, shells))
                            union!(shells, shell)
                            @test Set(DGG.neighbors(grid, c, k; connectivity = conn)) == shells
                        end
                        # The sequence checks also verify concatenation, tail
                        # blocks, and counter-clockwise winding.
                        if ordered
                            @test neighbor_order_problems(grid, c; connectivity = conn,
                                k = neighbor_k, require_rotational_rings) == String[]
                        end
                    end
                end
                ordered || skip!(skips, unimplemented_skip("the ring/disc order laws", sys,
                    "neighbors(::$(typeof(grids[first(tested)])), " *
                    "::$(DGG.cellindextype(sys)), ::Integer)"))
            end
        end

        @testset "descendant_range" begin
            if !DGG.has_sorted_subtrees(sys)
                skip!(skips, "skipped: descendant_range — $(subject_name(sys)) declares " *
                    "has_sorted_subtrees = false, so it offers no range to check")
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

        report_skips(skips, label)
    end
end

"""
    test_generic_fallbacks(sys; levels = <coarsest five>, n_levels = 3, n_samples = 4,
                           rng = MersenneTwister(20260813), label = <system type>,
                           kwargs...)

Property-test the **generic fallbacks** against `sys`'s geometry: the same laws
[`test_grid_interface`](@ref) and [`test_hierarchical_system`](@ref) check, run
against the implementations that answer when a system implements no fast path.

Default mode tests what a system implemented; this tests what the package will
actually answer for anyone who reaches it through a path the specializations do
not cover, and what a minimal implementation of the interface gets. The two
questions are different, and a system can pass one and fail the other: a
tessellation whose neighbouring cells do not share vertex coordinates, or whose
centroid is far enough off-centre to reorder a ring, breaks the geometric
fallbacks while its own closed-form `neighbors` stays right.

`sys` is wrapped in a [`GenericFallbackSystem`](@ref), which hides every
specialization by type, and each forced verb is asserted
[`dispatches_generically`](@ref) before its laws run — so a specialization that
found its way through the wrapper fails the run instead of quietly answering
for the fallback. The laws themselves run with `fallback_laws = true`, which is
what stops the gates from skipping a group because the answer comes from a
fallback.

This is expensive. The generic adjacency is a spatial-tree query and a
boundary-vertex comparison per cell, against a system's own `O(1)` arithmetic,
so it belongs in a certification run rather than in the loop a system is
developed in. Extra keywords (`atol`, `neighbor_k`, `connectivities`, …) are
forwarded to [`test_hierarchical_system`](@ref).

`levels` defaults to the **coarsest five**, not to all of them, because these
laws are about geometry that does not change with depth while two of the
harness's own tolerances do: a cell smaller than `1e-12` sr reads as a
degenerate boundary ring, and one whose diameter is under `unit_atol` reads as
explicitly closed. Deep levels are id arithmetic, which is what default mode
samples the deepest level for.
"""
function test_generic_fallbacks(sys;
        levels = _coarse_levels(sys),
        n_levels::Integer = 3,
        n_samples::Integer = 4,
        rng::Random.AbstractRNG = _default_rng(),
        label::AbstractString = string(nameof(typeof(sys))),
        kwargs...)
    wrapped = GenericFallbackSystem(sys)
    tested = sample_levels(rng, levels, n_levels)
    grid = DGG.levelgrid(wrapped, first(tested))
    probe = DGG.cellindex(grid, 1)

    @testset "generic fallbacks: $label" begin
        @testset "the wrapper reaches no specialization" begin
            @test dispatches_generically(DGG.cellat, grid, DGG.cell_centroid(grid, probe))
            @test dispatches_generically(DGG.neighbors, grid, probe, 1)
            @test dispatches_generically(DGG.ring, grid, probe, 1)
            @test dispatches_generically(DGG.ancestor, wrapped, probe, 0)
            @test dispatches_generically(DGG.descendants, wrapped, probe, 0)
            # ...while the contract the fallbacks consume is still the system's.
            @test DGG.cell_boundary(grid, probe) ==
                  DGG.cell_boundary(DGG.levelgrid(sys, first(tested)), probe)
            @test DGG.node_extent(wrapped, probe) == DGG.node_extent(sys, probe)
        end

        for l in tested
            test_grid_interface(DGG.levelgrid(wrapped, l); n_samples, rng,
                fallback_laws = true, label = "$label level $l (generic fallbacks)")
        end
        test_hierarchical_system(wrapped; levels = tested, n_levels = length(tested),
            n_samples, rng, fallback_laws = true,
            label = "$label (generic fallbacks)", kwargs...)
    end
end

end # module DiscreteGlobalGridsConformanceTesting
