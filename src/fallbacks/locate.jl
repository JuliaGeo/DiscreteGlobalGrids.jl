# ---------------------------------------------------------------------------
# Location and topology: `cellat`, geometric `neighbors`, `ring`
#
# Both are tree descents against the cell that is being asked about — the
# generic answers, correct for any grid, and the ones a system with closed-form
# arithmetic overrides. "Slow but correct" is the deal the design makes here:
# nothing in this file assumes a hierarchy, a projection, or a lattice.
# ---------------------------------------------------------------------------

"""
    cellat(grid, p::UnitSphericalPoint) -> Union{AbstractCellIndex,Nothing}
    cellat(grid, lon::Real, lat::Real)

The cell containing a point, or `nothing` outside the grid's coverage. The
`(lon, lat)` method takes **degrees** and converts; the unit-sphere method is
the primitive.

Generic implementation: descend [`treeify(grid)`](@ref treeify) to the cells
whose extents contain the point, then test point-in-cell exactly, on the
sphere. Candidates are visited in ascending canonical id order, so a point on a
shared boundary — which both cells legitimately contain — resolves to the first
of them, deterministically and never by floating-point luck.
"""
function cellat(grid::AbstractGrid, p::GO.UnitSphericalPoint)
    tree = treeify(grid)
    positions = STI.query(tree, cap -> cap_contains(cap, p))
    isempty(positions) && return nothing
    candidates = sort!([cellindex(grid, i) for i in positions])
    # A ring the exact predicate cannot decide (degenerate edges, a
    # near-hemispheric cell) is kept as a last resort rather than silently
    # dropped: `nothing` from the whole scan would claim the point is outside
    # the coverage, which is a different — and wrong — answer.
    undecided = nothing
    for c in candidates
        verdict = point_in_cell(cell_boundary(grid, c), p)
        verdict === true && return c
        verdict === nothing && undecided === nothing && (undecided = c)
    end
    return undecided
end

cellat(grid::AbstractGrid, lon::Real, lat::Real) = cellat(grid, unit_point(lon, lat))

# ===========================================================================
# Geometric adjacency
# ===========================================================================

"""
    adjacent_cells(grid, c, connectivity = Vertex(), tree = treeify(grid))

The cells of `grid` sharing at least a vertex (or, under [`Edge()`](@ref Edge),
at least two) with cell `c`, excluding `c` itself.

This is a **set-producing primitive**, not a neighbour answer: the order it
returns is ascending by canonical id, chosen only so that a tree query's
arrival order cannot leak into the result. The rotational order the interface
promises is *not* imposed here — [`neighbors`](@ref) and [`ring`](@ref) apply it
where the shells are assembled, because a shell's winding is defined about the
subject cell `c` and this function is called once per frontier cell rather than
once per shell. Callers wanting the contract order must go through `neighbors`
or `ring`.

It assumes a **conforming** tessellation: cells that meet do so at coincident
vertices, matched here to within a tolerance scaled to the cell's own size. That
covers every grid this package describes; a grid with T-junctions needs its own
`neighbors` method, which is the fast path every system writes anyway.

Candidate generation is sound: two cells that share a boundary point share a
point of both bounding caps, so no touching cell is ever pruned.
"""
function adjacent_cells(grid::AbstractGrid, c::AbstractCellIndex,
        connectivity::Connectivity=Vertex(), tree=treeify(grid))
    boundary = cell_boundary(grid, c)
    cap = cell_cap(grid, c)
    tol = _match_tolerance(boundary)
    needed = connectivity isa Edge ? 2 : 1
    out = typeof(c)[]
    for i in STI.query(tree, other -> intersects_cap(cap, other))
        d = cellindex(grid, i)
        d == c && continue
        _shared_vertices(boundary, cell_boundary(grid, d), tol) >= needed && push!(out, d)
    end
    return sort!(out)
end

# How close two vertices must be to count as the same corner: a thousandth of
# the cell's own shortest edge. Scaled from the ring rather than from its
# bounding cap, because `cell_cap` degrades to the full sphere for a cell wider
# than a hemisphere and a tolerance derived from that would be a radian wide.
function _match_tolerance(boundary)
    ring, n = open_ring(boundary)
    shortest = Inf
    for i in 1:n
        a = ring[i]
        b = ring[i == n ? 1 : i+1]
        d2 = (a[1] - b[1])^2 + (a[2] - b[2])^2 + (a[3] - b[3])^2
        d2 > 0 && (shortest = min(shortest, d2))
    end
    isfinite(shortest) || return 1e-9        # an all-degenerate ring
    return max(1e-12, 1e-3 * sqrt(shortest))
end

function _shared_vertices(a, b, tol::Float64)
    tol2 = tol * tol
    shared = 0
    for p in a
        for q in b
            dx = p[1] - q[1]
            dy = p[2] - q[2]
            dz = p[3] - q[3]
            if dx * dx + dy * dy + dz * dz <= tol2
                shared += 1
                break
            end
        end
        shared >= 2 && return shared      # nothing above 2 changes any answer
    end
    return shared
end

"""
    neighbors(grid, c, k = 1; connectivity = Vertex())

All cells within `k` adjacency steps of `c`, excluding `c`, as the rings
1, 2, …, `k` **concatenated outward** — see the [`neighbors`](@ref) interface
docstring, which is authoritative for the order.

The fallback has no lattice to read a direction off, so it realises the
rotational contract by measurement: each shell is ordered by azimuth about
`cell_centroid(grid, c)`, counter-clockwise seen from outside the sphere, with
the **first ring-1 neighbour** as the zero direction so that every shell starts
on the same spoke. Ring 1's own start is its smallest canonical id, which is the
only arbitrary choice in the scheme and the one that makes it reproducible.
Exact azimuth ties break by ascending canonical id. Nothing here sorts the
result by id: that would interleave the shells and destroy the tail-block law
with [`ring`](@ref).

The walk itself is a breadth-first sweep over [`adjacent_cells`](@ref), which is
geometric and therefore slow; it exists so that every grid has correct
neighbours, and every system overrides it. Cells outside the grid's coverage are
simply never produced, which is the "absent, not padded" rule holding by
construction.
"""
function neighbors(grid::AbstractGrid, c::AbstractCellIndex, k::Integer=1;
        connectivity::Connectivity=Vertex())
    shells = _adjacency_shells(grid, c, Int(k), connectivity)
    isempty(shells) && return typeof(c)[]
    return reduce(vcat, shells)
end

# The breadth-first walk both `neighbors` and `ring` are reads of: shell `j` of
# the result is the cells at adjacency distance exactly `j`, already wound. One
# tree, built once — for a grid with no hierarchy that build is O(ncells), so
# doing it per step (or twice, for `ring`'s two discs) is the difference between
# one sweep and several.
#
# The winding happens HERE rather than in the two callers, because that is what
# makes `ring(grid, c, k)` the literal tail block of `neighbors(grid, c, k)`:
# both read the same wound shells, so they cannot disagree element for element.
function _adjacency_shells(grid::AbstractGrid, c::AbstractCellIndex, steps::Int,
        connectivity::Connectivity)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    T = typeof(c)
    shells = Vector{T}[]
    steps == 0 && return shells
    tree = treeify(grid)
    seen = Set{T}((c,))
    frontier = T[c]
    centre = cell_centroid(grid, c)
    # Fixed once, from ring 1, and reused by every outer shell: one spoke for
    # the whole disc is what makes position `j` of a ring mean a direction.
    frame = nothing
    for _ in 1:steps
        next = T[]
        for x in frontier
            for y in adjacent_cells(grid, x, connectivity, tree)
                y in seen && continue
                push!(seen, y)
                push!(next, y)
            end
        end
        if frame === nothing && !isempty(next)
            frame = _ring_frame(grid, centre, next)
        end
        frame === nothing || _wind!(next, grid, centre, frame)
        push!(shells, next)
        isempty(next) && break
        frontier = next
    end
    return shells
end

"""
    ring(grid, c, k; connectivity = Vertex())

The cells at adjacency distance **exactly** `k`. `ring(grid, c, 0)` is `[c]`.

Shares [`neighbors`](@ref)' breadth-first walk — one shell of it, rather than
the difference of two discs — and therefore its order exactly: counter-clockwise
seen from outside the sphere, from the same spoke every other shell starts on.
The result is the tail block of `neighbors(grid, c, k)` element for element,
which is the law the interface docstring states and the reason this reads a
shell instead of walking its own.
"""
function ring(grid::AbstractGrid, c::AbstractCellIndex, k::Integer;
        connectivity::Connectivity=Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    steps == 0 && return typeof(c)[c]
    shells = _adjacency_shells(grid, c, steps, connectivity)
    # A walk that ran out of cells before reaching `steps` has an empty shell
    # there: the ring is genuinely empty, not missing.
    steps <= length(shells) || return typeof(c)[]
    return shells[steps]
end

# ===========================================================================
# Rotational shell order
#
# The fallback's own azimuth machinery. It deliberately duplicates what the
# HEALPix and H3 systems do with their own lattices: this layer may not depend
# on a system module, and a system that has a lattice direction should be using
# it rather than measuring one.
# ===========================================================================

# The frame every shell's azimuth is measured in, built once from ring 1:
# the tangent basis at `centre` whose zero direction points at the ring-1
# neighbour with the smallest canonical id.
#
# `zero` is the anchor's own measured azimuth rather than the literal `0.0`.
# In exact arithmetic they are the same number; in floating point the anchor's
# tangential component can round to a hair below zero, and `mod(-eps, 2pi)` is
# `2pi` — which would sort the spoke's own cell to the END of its ring.
# Subtracting the measurement cancels that exactly.
function _ring_frame(grid::AbstractGrid, centre, shell::AbstractVector)
    anchor = cell_centroid(grid, minimum(shell))
    e1, e2 = _tangent_basis(centre, anchor)
    return (e1, e2, _azimuth(centre, e1, e2, anchor))
end

# Order one shell counter-clockwise about `centre`, from the frame's spoke.
# Exact ties go to the smaller canonical id, so the result is total.
function _wind!(shell::AbstractVector, grid::AbstractGrid, centre, frame)
    length(shell) <= 1 && return shell
    e1, e2, zero = frame
    turn = 2 * Float64(pi)
    sort!(shell; by=d -> (mod(_azimuth(centre, e1, e2, cell_centroid(grid, d)) - zero,
            turn), d))
    return shell
end

# A tangent basis at `centre` with `e1` pointing at `toward` (projected into the
# tangent plane) and `e2 = centre x e1`. Then `e1 x e2 == centre`, i.e. the pair
# is right-handed SEEN FROM OUTSIDE the sphere, which is what makes increasing
# `atan(u.e2, u.e1)` run counter-clockwise viewed from above the plane rather
# than from below it.
function _tangent_basis(centre, toward)
    u = (toward[1] - centre[1], toward[2] - centre[2], toward[3] - centre[3])
    radial = u[1] * centre[1] + u[2] * centre[2] + u[3] * centre[3]
    t = (u[1] - radial * centre[1], u[2] - radial * centre[2], u[3] - radial * centre[3])
    n = sqrt(t[1]^2 + t[2]^2 + t[3]^2)
    # A neighbour whose centroid coincides with the subject's, or sits exactly
    # antipodal to it, has no tangential direction. Neither can happen in a real
    # tessellation, but an arbitrary-but-valid basis keeps this total instead of
    # quietly returning NaNs and corrupting an order.
    if n <= eps(Float64)
        t = abs(centre[3]) < 0.9 ? (0.0, 0.0, 1.0) : (1.0, 0.0, 0.0)
        radial = t[1] * centre[1] + t[2] * centre[2] + t[3] * centre[3]
        t = (t[1] - radial * centre[1], t[2] - radial * centre[2], t[3] - radial * centre[3])
        n = sqrt(t[1]^2 + t[2]^2 + t[3]^2)
    end
    e1 = (t[1] / n, t[2] / n, t[3] / n)
    e2 = (centre[2] * e1[3] - centre[3] * e1[2],
        centre[3] * e1[1] - centre[1] * e1[3],
        centre[1] * e1[2] - centre[2] * e1[1])
    return e1, e2
end

# The azimuth of `p` about `centre`, in the `(e1, e2)` frame: `p`'s offset
# projected into the tangent plane, then read as an angle.
function _azimuth(centre, e1, e2, p)
    u = (p[1] - centre[1], p[2] - centre[2], p[3] - centre[3])
    radial = u[1] * centre[1] + u[2] * centre[2] + u[3] * centre[3]
    t = (u[1] - radial * centre[1], u[2] - radial * centre[2], u[3] - radial * centre[3])
    return atan(t[1] * e2[1] + t[2] * e2[2] + t[3] * e2[3],
        t[1] * e1[1] + t[2] * e1[2] + t[3] * e1[3])
end
