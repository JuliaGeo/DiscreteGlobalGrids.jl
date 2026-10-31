# ===========================================================================
# Oracle helpers the system suites share: the problem collectors that used to
# live in the conformance library this package used to ship (deleted in favour of
# GridInterface's `GridInterfaceConformance`, which states the interface laws
# for every GridInterface provider). What stays here is this package's own,
# stricter-than-the-interface vocabulary — rotational order of rings, the
# two-hop neighbour closure, the id-form/index-form agreement — which the
# system suites still assert cell by cell.
#
# Included by each suite that needs it:  include("../../conformance_helpers.jl")
# ===========================================================================
module ConformanceHelpers

import Random
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: Connectivity, Vertex, Edge

export spherical_signed_area, sample_indices, sample_cells, sample_levels,
    boundary_problems, neighbor_problems, winding_problems, ring_cycle_problems,
    neighbor_order_problems

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
    sample_indices(rng, nindices, n_samples) -> Vector{Int}

`n_samples` distinct indices drawn from `1:nindices`, sorted; all of them
when the grid is small enough to check exhaustively.

Sampling is with replacement then deduplicated rather than a `randperm`,
because a level grid may have more cells than there is memory for a
permutation of them.
"""
function sample_indices(rng, nindices::Integer, n_samples::Integer)
    n = Int(nindices)
    n <= 0 && return Int[]
    n <= n_samples && return collect(1:n)
    return sort!(unique(rand(rng, 1:n, Int(n_samples))))
end

"""
    sample_cells(rng, grid, n_samples)

The typed ids at [`sample_indices`](@ref), paired with those indices.
"""
function sample_cells(rng, grid, n_samples::Integer)
    indices = sample_indices(rng, DGG.ncells(grid), n_samples)
    return indices, [DGG.cellindex(grid, i) for i in indices]
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
    neighbor_problems(grid, c; connectivity, sys, two_hop) -> Vector{String}

The [`neighbors`](@ref) laws at one cell: determinism across calls, exclusion of
`c` itself, distinctness, a container typed at the system's cell index type, a
count within [`maxneighbors`](@ref), membership in the grid, a common level,
and symmetry — `c′ ∈ neighbors(c)` implies `c ∈ neighbors(c′)`, checked in both
directions (see `two_hop` below).

This function is about *membership*; the ordering half of the contract lives in
[`neighbor_order_problems`](@ref). Nothing here requires or rewards an
id-sorted result — the order is rotational, and a harness that sorted would be
asserting the opposite of the contract.

Partial coverage is handled by the contract rather than by an exemption: a
neighbour beyond the grid's coverage must be **absent** from the result, so a
returned cell that has no [`localindex`](@ref) is itself a violation, and
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
        bound = DGG.maxneighbors(sys, connectivity)
        # `nothing` means the system declares no static bound, so there is no
        # ceiling to check against.
        bound === nothing || length(ns) <= bound ||
            push!(problems, "neighbors($c) returned $(length(ns)) cells, over maxneighbors of $bound")
    end

    lc = DGG.level(c)
    for nb in ns
        DGG.level(nb) == lc ||
            push!(problems, "neighbour $nb of $c is at level $(DGG.level(nb)), not $lc")
        if DGG.localindex(grid, nb) === nothing
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
            DGG.localindex(grid, nb) === nothing && continue
            for x in DGG.neighbors(grid, nb, 1; connectivity)
                x in checked && continue
                push!(checked, x)
                DGG.localindex(grid, x) === nothing && continue
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

# One wrapper per literal `k`, so inference sees the same shape a caller's
# kernel does: the grid and the cell arrive as arguments while `k` is written
# into the source. Constant-propagating `k` is what collapses the branches, and

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
  - **idiom agreement** — the index forms answer the id forms read through
    [`localindex`](@ref), element for element, so the two idioms cannot
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
    append!(problems, _index_form_problems(grid, c, rings; connectivity))
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

# The index form is the id form read through `localindex`, element for
# element. Written once, here, so no system can present one neighbourhood in two
# orders — the drift that makes an oriented stencil silently change meaning when
# a caller moves from ids to indices.
function _index_form_problems(grid, c, rings; connectivity::Connectivity)
    problems = String[]
    p = DGG.localindex(grid, c)
    p isa Int || return problems
    for (j, shell) in enumerate(rings)
        want = Int[]
        for x in shell
            q = DGG.localindex(grid, x)
            q isa Int && push!(want, q)
        end
        got = DGG.ring(grid, p, j; connectivity)
        collect(got) == want || push!(problems,
            "ring(grid, $p, $j) is not ring(grid, $c, $j) read through localindex: " *
            "got $got, the ids map to $want. The index and id forms are one order")
    end
    return problems
end

# ===========================================================================
# Running GridInterface's conformance checks as tests, with the gaps this
# package knows about recorded as broken rather than red.
#
# `GridInterfaceConformance.test_grid`/`test_system` emit one `@test` per law.
# A few laws fail on shipped systems for reasons outside a test run's reach
# (an allocating pruning cap on the projected systems, point location in the
# sliver between libh3's bent edge and the straight one, near-edge location on
# IGeo7, …); each system's suite names those as `known`, and they are emitted
# with `@test_broken` so a fix shows up as a broken test starting to pass.
# ===========================================================================

using Test
import GridInterfaceConformance
const GIC = GridInterfaceConformance

export conformance_grid, conformance_system

"Whether `law` (a report law string) is one of the `known` gaps (substrings)."
_known(law::AbstractString, known) = any(k -> occursin(k, law), known)

"""
    emit_conformance(rep; known = String[])

`GridInterfaceConformance.emit_tests`, with a failing law whose name contains
one of `known` emitted as `@test_broken`. Every failure is printed with its
detail, so the log says what failed, not just that something did.
"""
function emit_conformance(rep; known=String[])
    @testset "$(rep.name)" begin
        for fam in GIC.FAMILIES
            rs = filter(r -> r.family === fam, rep.results)
            isempty(rs) && continue
            @testset "$fam" begin
                for r in rs
                    if r.status === :pass
                        @test true
                    elseif r.status === :fail
                        status = _known(r.law, known) ? "known gap" : "FAIL"
                        println("  [", status, "] ", rep.name, " / ", fam, " / ", r.law,
                            isempty(r.detail) ? "" : ": " * r.detail)
                        if _known(r.law, known)
                            @test_broken false
                        else
                            @test false
                        end
                    elseif r.status === :skip
                        @test_skip "$fam: $(r.detail)"
                    end
                end
            end
        end
    end
    return rep
end

"""
    conformance_grid(grid; known = String[], name, nsamples = 64, kw...)

`GridInterfaceConformance.test_grid` with the known gaps recorded as broken.
"""
conformance_grid(grid; known=String[], nsamples::Int=64, kw...) =
    emit_conformance(GIC.check_grid(grid; nsamples, kw...); known)

"""
    conformance_system(sys; levels = (2, 3), known = String[], name, kw...)

`GridInterfaceConformance.test_system` with the known gaps recorded as broken:
the hierarchy laws, one level grid and a sample multi-level set at
`first(levels)`, then the grid laws at each further level.
"""
function conformance_system(sys; levels=(2, 3), known=String[],
        name::AbstractString=string(nameof(typeof(sys))), kw...)
    ls = [l for l in levels if l in DGG.levels(sys)]
    isempty(ls) && (ls = [DGG.levels(sys)[min(end, 3)]])
    rep = emit_conformance(GIC.check_system(sys; grid_level=first(ls), name, kw...); known)
    for l in ls[2:end]
        emit_conformance(GIC.check_grid(DGG.levelgrid(sys, l); name="$name level $l", kw...); known)
    end
    return rep
end

end # module ConformanceHelpers
