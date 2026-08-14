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
with `c`, excluding `c`. Results are sorted by canonical id; public rotational
ordering is applied by [`neighbors`](@ref) and [`ring`](@ref).

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

Return cells within `k` adjacency steps, excluding `c`, as outward-concatenated
rings. A breadth-first traversal orders each ring counter-clockwise by measured
azimuth from the smallest-id ring-1 neighbor; exact ties use canonical id.
Cells outside grid coverage are omitted.

This walk measures distance INSIDE the grid it is given, so on a subset it
answers induced-subgraph distance rather than the clipped-to-membership law
[`neighbors`](@ref) states — a hole here lengthens the path around itself. The
two agree at `k == 1` and part company from `k == 2`. The law belongs to the
named subset types ([`PartialGrid`](@ref), [`CellVector`](@ref), `CellLookup`),
which override this method; a new subset grid owes its own override.
"""
function neighbors(grid::AbstractGrid, c::AbstractCellIndex, k::Integer=1;
        connectivity::Connectivity=Vertex())
    shells = _adjacency_shells(grid, c, Int(k), connectivity)
    isempty(shells) && return typeof(c)[]
    return reduce(vcat, shells)
end

# Shared breadth-first traversal producing already wound distance shells. It
# builds one tree and makes `ring(..., k)` the exact tail of `neighbors(..., k)`.
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

Return cells at adjacency distance exactly `k`; `k == 0` returns `[c]`. Ordering
matches the corresponding tail block of [`neighbors`](@ref).
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

# Rotational shell ordering by measured azimuth.

# Tangent frame anchored at the smallest-id ring-1 neighbor. Subtract its
# measured azimuth to prevent `mod(-eps, 2π)` from moving the anchor to the end.
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
