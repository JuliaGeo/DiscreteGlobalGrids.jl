# Generic location and adjacency operations for grids without closed-form
# hierarchy or lattice methods.

"""
    cellat(grid, p::UnitSphericalPoint) -> Union{AbstractCellIndex,Nothing}
    cellat(grid, lon::Real, lat::Real)

Return the cell containing a point, or `nothing` outside grid coverage. The
`(lon, lat)` overload accepts degrees. Candidates are pruned by the spatial tree,
tested exactly, and ordered by canonical id for deterministic boundary ties.
"""
function cellat(grid::AbstractGrid, p::GO.UnitSphericalPoint)
    tree = treeify(grid)
    positions = STI.query(tree, cap -> cap_contains(cap, p))
    isempty(positions) && return nothing
    candidates = sort!([cellindex(grid, i) for i in positions])
    # Preserve one undecidable candidate rather than report it outside coverage.
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

Return cells sharing at least one vertex, or two under [`Edge()`](@ref Edge),
with `c`, excluding `c`. Results are sorted by canonical id; the public
rotational ordering is applied by [`one_ring`](@ref).

The fallback assumes a conforming tessellation with coincident shared vertices.
Cap intersection safely prunes candidates because touching cells have
intersecting caps.
"""
function adjacent_cells(grid::AbstractGrid, c::AbstractCellIndex,
        connectivity::Connectivity=Vertex(), tree=treeify(grid))
    boundary = cell_boundary(grid, c)
    cap = cell_cap(grid, c)
    tol = _match_tolerance(boundary)
    needed = connectivity isa Edge ? 2 : 1
    out = typeof(c)[]
    for i in STI.query(tree, Base.Fix1(Extents.intersects, cap))
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
    one_ring(grid, c, connectivity, tree = treeify(grid))

The geometric one-ring: [`adjacent_cells`](@ref) wound counter-clockwise about
`c`'s centroid, starting at the smallest-id neighbour. `tree` is the traversal
pruner, taken as an argument so a shell walk builds one tree rather than one per
cell.

This is the generic method of the `one_ring` hook; a system with native adjacency
overrides the three-argument form.
"""
function one_ring(grid::AbstractGrid, c::AbstractCellIndex,
        connectivity::Connectivity=Vertex(), tree=treeify(grid))
    shell = adjacent_cells(grid, c, connectivity, tree)
    length(shell) <= 1 && return shell
    centre = cell_centroid(grid, c)
    return _wind!(shell, grid, centre, _ring_frame(grid, centre, first(shell)))
end

"""
    checked_steps(k) -> Int

`k` as an `Int`, or an `ArgumentError` naming it. The one place the neighbourhood
family's non-negativity precondition is written.
"""
function checked_steps(k::Integer)
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    return steps
end

"""
    adjacency_shells(walk, grid, c, steps) -> Vector{Vector{eltype}}
    adjacency_shells(grid, c, steps, connectivity)

The shared breadth-first shell walk: entry `j` holds the cells at adjacency
distance exactly `j`, counter-clockwise, so that `ring(grid, c, k)` is the exact
tail of `neighbors(grid, c, k)` and neither can drift from the other.

`walk(x)` supplies one cell's ordered one-ring; the four-argument form uses
`one_ring`. Ring 1 arrives already wound and keeps the system's own start; rings
`2:steps` come out of the frontier in whatever order it reached them and are
wound about `c`'s centroid from the spoke through ring 1's first cell.

A walk that exhausts the component stops early, so a shell past the end is
absent rather than empty.
"""
adjacency_shells(walk::W, grid::AbstractGrid, c::AbstractCellIndex,
    steps::Int) where {W} = adjacency_shells(walk, grid, c, steps, Unordered())

# A declared turn is carried outward; anything else is measured.
adjacency_shells(walk::W, grid::AbstractGrid, c::AbstractCellIndex, steps::Int,
    ::Union{CounterClockwise,Clockwise}) where {W} =
    _shells_wound(walk, grid, c, steps)

adjacency_shells(walk::W, grid::AbstractGrid, c::AbstractCellIndex, steps::Int,
    ::Winding) where {W} = _shells_azimuth(walk, grid, c, steps)

adjacency_shells(grid::AbstractGrid, c::AbstractCellIndex, steps::Int,
    connectivity::Connectivity) =
    adjacency_shells(x -> one_ring(grid, x, connectivity), grid, c, steps,
        winding(grid, connectivity))

# The geometric walk hoists the spatial tree out of the per-cell one-ring.
function _geometric_shells(grid::AbstractGrid, c::AbstractCellIndex, steps::Int,
        connectivity::Connectivity)
    tree = treeify(grid)
    return adjacency_shells(x -> one_ring(grid, x, connectivity, tree), grid, c,
        steps, winding(grid, connectivity))
end

# --- the measured walk: no declared turn, so every outer shell is sorted -----

function _shells_azimuth(walk::W, grid::AbstractGrid, c::AbstractCellIndex,
        steps::Int) where {W}
    T = typeof(c)
    shells = Vector{T}[]
    steps == 0 && return shells
    seen = Set{T}((c,))
    frontier = T[c]
    # Fixed once, from ring 1, and reused by every outer shell: one spoke for
    # the whole disc is what makes position `j` of a ring mean a direction.
    centre = nothing
    frame = nothing
    for j in 1:steps
        next = T[]
        for x in frontier
            for y in walk(x)
                y in seen && continue
                push!(seen, y)
                push!(next, y)
            end
        end
        if j > 1 && frame === nothing
            centre = cell_centroid(grid, c)
            frame = _ring_frame(grid, centre, first(first(shells)))
        end
        frame === nothing || _wind!(next, grid, centre, frame)
        push!(shells, next)
        isempty(next) && break
        frontier = next
    end
    return shells
end

# --- the carried walk: the declared turn supplies the order ------------------
#
# Each frontier cell's one-ring is already a turn, so the cells it adds to the
# next shell come off it in that turn's order once the arc is started at the
# right place: just past the inward neighbours, which are contiguous in the
# turn. Concatenating over a frontier that is itself in turn order gives the
# whole shell in turn order, with no centroid and no sort.
#
# What that does NOT give is the shell's PHASE — which of its cells the ring
# starts at. `_pin_phase!` supplies it, and is the only geometry here.

function _shells_wound(walk::W, grid::AbstractGrid, c::AbstractCellIndex,
        steps::Int) where {W}
    T = typeof(c)
    shells = Vector{T}[]
    steps == 0 && return shells
    r1 = collect(T, walk(c))
    push!(shells, r1)
    (steps == 1 || isempty(r1)) && return shells
    centre = cell_centroid(grid, c)
    frame = _ring_frame(grid, centre, first(r1))
    prev2 = T[c]
    prev1 = r1
    for _ in 2:steps
        next = T[]
        for f in prev1
            R = walk(f)
            m = length(R)
            m == 0 && continue
            # Start the arc after the last inward neighbour. The inward set is
            # contiguous in the turn, so one index is the whole boundary.
            a = 0
            for i in 1:m
                _member(prev2, R[i]) && (a = i)
            end
            for t in 1:m
                y = @inbounds R[mod1(a + t, m)]
                (_member(prev2, y) || _member(prev1, y) || _member(next, y)) &&
                    continue
                push!(next, y)
            end
        end
        _pin_phase!(next, grid, centre, frame)
        push!(shells, next)
        isempty(next) && break
        prev2, prev1 = prev1, next
    end
    return shells
end

# Linear membership over a shell. Shells are `O(k)` and hold no duplicates, so
# this beats hashing at every size the walk reaches.
@inline function _member(v, y)
    for x in v
        x == y && return true
    end
    return false
end

"""
    _pin_phase!(shell, grid, centre, frame) -> shell

Rotate `shell` — already in turn order — so it starts where [`_wind!`](@ref)
would have put it: at the smallest phase about `frame`'s zero spoke.

Phase increases along a turn and wraps exactly once, so the starting cell is the
wrap point of a rotated sorted sequence and a binary search finds it in
`O(log n)` centroids rather than one per member. That search is only valid
because the sequence really is a rotation of a sorted one, which is what a
declared [`winding`](@ref) asserts and `test_grid_interface` checks.
"""
function _pin_phase!(shell::AbstractVector, grid, centre, frame)
    n = length(shell)
    n <= 1 && return shell
    e1, e2, zero = frame
    ph(i) = _phase(_azimuth(centre, e1, e2,
        cell_centroid(grid, @inbounds shell[i])) - zero)
    # Against a FIXED reference, so the search costs one centroid per step
    # rather than two: phase rises to the end of the turn and restarts, so
    # `phase < phase[1]` is false along the tail and true from the wrap on, and
    # the first index where it holds is the start. No wrap means index 1.
    first_phase = ph(1)
    first_phase <= ph(n) && return shell
    lo, hi = 1, n
    while lo < hi
        mid = (lo + hi) >> 1
        if ph(mid) >= first_phase
            lo = mid + 1
        else
            hi = mid
        end
    end
    lo == 1 && return shell
    rot = similar(shell)
    @inbounds for i in 1:n
        rot[i] = shell[mod1(lo + i - 1, n)]
    end
    copyto!(shell, rot)
    return shell
end

"""
    neighbors(grid, c, k = 1; connectivity = Vertex())

Return cells within `k` adjacency steps, excluding `c`, as outward-concatenated
rings in the counter-clockwise order [`neighbors`](@ref) states, starting at the
smallest-id ring-1 neighbour. Cells outside grid coverage are omitted.

This walk measures distance INSIDE the grid it is given, so on a subset it
answers induced-subgraph distance rather than the clipped-to-membership law
[`neighbors`](@ref) states — a hole here lengthens the path around itself. The
two agree at `k == 1` and part company from `k == 2`. The law belongs to the
named subset types ([`PartialGrid`](@ref), [`CellVector`](@ref), `CellLookup`),
which override this method; a new subset grid owes its own override.
"""
function neighbors(grid::AbstractGrid, c::AbstractCellIndex, k::Integer=1;
        connectivity::Connectivity=Vertex())
    steps = checked_steps(k)
    steps == 0 && return typeof(c)[]
    steps == 1 && return one_ring(grid, c, connectivity)
    shells = _geometric_shells(grid, c, steps, connectivity)
    isempty(shells) && return typeof(c)[]
    return reduce(vcat, shells)
end

"""
    ring(grid, c, k; connectivity = Vertex())

Return cells at adjacency distance exactly `k`; `k == 0` returns `[c]`. Ordering
matches the corresponding tail block of [`neighbors`](@ref).
"""
function ring(grid::AbstractGrid, c::AbstractCellIndex, k::Integer;
        connectivity::Connectivity=Vertex())
    steps = checked_steps(k)
    steps == 0 && return typeof(c)[c]
    shells = _geometric_shells(grid, c, steps, connectivity)
    # A walk that ran out of cells before reaching `steps` has an empty shell
    # there: the ring is genuinely empty, not missing.
    steps <= length(shells) || return typeof(c)[]
    return shells[steps]
end

# Rotational shell ordering by measured azimuth.

"""
    _ring_frame(grid, centre, anchor) -> (e1, e2, zero)

Tangent basis at `centre` whose zero azimuth points at the cell `anchor`.
Storing the anchor's own azimuth keeps `mod(-eps, 2π)` from winding it to the
end of the order.

An internal extension point, re-exported as
`DiscreteGlobalGrids._ring_frame`: a system walking its own adjacency shells
pairs it with `_wind!` to get the package's ring ordering rather than a private
one.
"""
function _ring_frame(grid, centre, anchor::AbstractCellIndex)
    p = cell_centroid(grid, anchor)
    e1, e2 = _tangent_basis(centre, p)
    return (e1, e2, _azimuth(centre, e1, e2, p))
end

"""
    _wind!(shell, grid, centre, frame) -> shell

Sort `shell` in place counter-clockwise about `centre` from `frame`'s zero
spoke, `frame` being a `_ring_frame` triple. Exact azimuth ties go to the
smaller canonical id, so the order is total.

A cell lying ON the spoke begins the ring rather than ends it: `mod` sends an
azimuth a hair below the spoke to nearly a full turn, and which side of an
exactly-aligned cell falls on is decided by the last bit of an `atan`. Phases
within [`SPOKE_ATOL`](@ref) of a full turn are therefore read as zero. This is
the only tolerance in the order contract, and it is far below the angular
separation of two distinct cells at any level.

An internal extension point, re-exported as `DiscreteGlobalGrids._wind!`.
"""
function _wind!(shell::AbstractVector, grid, centre, frame)
    length(shell) <= 1 && return shell
    e1, e2, zero = frame
    # One centroid per member, not one per comparison: an outer shell has O(k)
    # members and `sort!(by = ...)` would project each of them O(log k) times.
    keyed = [(_phase(_azimuth(centre, e1, e2, cell_centroid(grid, d)) - zero), d)
             for d in shell]
    sort!(keyed)
    for i in eachindex(shell, keyed)
        @inbounds shell[i] = keyed[i][2]
    end
    return shell
end

"""
    SPOKE_ATOL

How close to a full turn a phase must be to count as sitting on the ring's
starting spoke: `1e-9` radians, orders of magnitude below the angular separation
of two distinct cells and orders of magnitude above the `atan` noise it exists
to absorb.
"""
const SPOKE_ATOL = 1e-9

@inline function _phase(a::Float64)
    turn = 2 * Float64(pi)
    p = mod(a, turn)
    return p >= turn - SPOKE_ATOL ? 0.0 : p
end

# Right-handed tangent basis at `centre`, with `e1` toward the anchor and
# `e2 = centre × e1`; increasing azimuth is counter-clockwise from outside.
function _tangent_basis(centre, toward)
    u = (toward[1] - centre[1], toward[2] - centre[2], toward[3] - centre[3])
    radial = u[1] * centre[1] + u[2] * centre[2] + u[3] * centre[3]
    t = (u[1] - radial * centre[1], u[2] - radial * centre[2], u[3] - radial * centre[3])
    n = sqrt(t[1]^2 + t[2]^2 + t[3]^2)
    # Use an arbitrary valid tangent for coincident or antipodal centroids.
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
