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
        ring = cell_boundary(grid, c)
        verdict = US.spherical_ring_encloses(ring, length(ring), p)
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
at least two) with cell `c`, ascending by id and excluding `c` itself.

This is the geometric fallback [`neighbors`](@ref) is built on, and it assumes
a **conforming** tessellation: cells that meet do so at coincident vertices,
matched here to within a tolerance scaled to the cell's own size. That covers
every grid this package describes; a grid with T-junctions needs its own
`neighbors` method, which is the fast path every system writes anyway.

Candidate generation is sound: two cells that share a boundary point share a
point of both bounding caps, so no touching cell is ever pruned.
"""
function adjacent_cells(grid::AbstractGrid, c::AbstractCellIndex,
        connectivity::Connectivity=Vertex(), tree=treeify(grid))
    boundary = cell_boundary(grid, c)
    cap = cell_cap(grid, c)
    tol = 1e-4 * Float64(cap.radius) + 1e-12
    needed = connectivity isa Edge ? 2 : 1
    out = typeof(c)[]
    for i in STI.query(tree, other -> intersects_cap(cap, other))
        d = cellindex(grid, i)
        d == c && continue
        _shared_vertices(boundary, cell_boundary(grid, d), tol) >= needed && push!(out, d)
    end
    return sort!(out)
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

All cells within `k` adjacency steps of `c`, excluding `c`, ascending by
canonical id — see the interface docstring for the contract.

The generic implementation is a breadth-first walk over
[`adjacent_cells`](@ref), which is geometric and therefore slow; it exists so
that every grid has correct neighbours, and every system overrides it. Cells
outside the grid's coverage are simply never produced, which is the "absent,
not padded" rule holding by construction.
"""
function neighbors(grid::AbstractGrid, c::AbstractCellIndex, k::Integer=1;
        connectivity::Connectivity=Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    T = typeof(c)
    out = T[]
    steps == 0 && return out
    tree = treeify(grid)
    seen = Set{T}((c,))
    frontier = T[c]
    for _ in 1:steps
        next = T[]
        for x in frontier
            for y in adjacent_cells(grid, x, connectivity, tree)
                y in seen && continue
                push!(seen, y)
                push!(next, y)
                push!(out, y)
            end
        end
        isempty(next) && break
        frontier = next
    end
    return sort!(out)
end

"""
    ring(grid, c, k; connectivity = Vertex())

The cells at adjacency distance **exactly** `k`. `ring(grid, c, 0)` is `[c]`.

Derived from [`neighbors`](@ref) as the difference of two discs, and
overridable by any system that can walk a shell directly.
"""
function ring(grid::AbstractGrid, c::AbstractCellIndex, k::Integer;
        connectivity::Connectivity=Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    steps == 0 && return typeof(c)[c]
    outer = neighbors(grid, c, steps; connectivity)
    inner = neighbors(grid, c, steps - 1; connectivity)
    return sort!(setdiff(outer, inner))
end
