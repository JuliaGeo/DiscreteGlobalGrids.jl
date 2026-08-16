# A5 supplies unordered adjacency sets. Shells are therefore sorted by azimuth
# counter-clockwise in a tangent frame. Ring 1 starts at the smallest id; outer
# rings use the same spoke. `edge_only=true` implements `Edge()` connectivity.

# Container capacity for one-step results under either connectivity.
const MAX_NEIGHBORS = 11

_edge_only(::Vertex) = false
_edge_only(::Edge) = true

# Return one-step adjacency as raw ids after validating the cell and grid level.
function _native_neighbors(grid::LevelGrid, c::A5Cell, connectivity::Connectivity)
    level(c) == grid.level || throw(ArgumentError(
        "A5 cell $c is at resolution $(level(c)), not this grid's $(grid.level)"))
    isvalid(c) || throw(ArgumentError("A5 cell $c is not a valid cell"))
    return A5Native._get_global_cell_neighbors(c.id; edge_only=_edge_only(connectivity))
end

"""
    neighbors(grid::LevelGrid, c::A5Cell, k = 1; connectivity = Vertex())

Cells within `k` grid steps, excluding `c`, as counter-clockwise shells
concatenated outward. Ring 1 starts at the smallest [`A5Cell`](@ref) id and all
outer rings use the same spoke. [`ring`](@ref) is the final shell.

[`Vertex()`](@ref Vertex) includes corner-only neighbours; [`Edge()`](@ref Edge)
does not. See [`max_neighbors`](@ref).

`k <= 1` returns a `SmallVector{11,A5Cell}` — sized by the `Vertex()` bound
under both connectivities. `k >= 2` returns a `Vector{A5Cell}`.

Throws `ArgumentError` unless `c` is valid at the grid resolution.
"""
function neighbors(grid::LevelGrid, c::A5Cell, k::Integer=1;
        connectivity::Connectivity=Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    steps == 0 && return SmallVector{MAX_NEIGHBORS,A5Cell}()
    shells = _shells(grid, c, steps, connectivity)
    if steps == 1
        out = SmallVector{MAX_NEIGHBORS,A5Cell}()
        for x in @inbounds shells[1]
            out = SmallCollections.push(out, x)
        end
        return out
    end
    return reduce(vcat, shells)
end

"""
    ring(grid::LevelGrid, c::A5Cell, k; connectivity = Vertex())

The cells at grid distance **exactly** `k` from `c`, counter-clockwise seen from
outside. `ring(grid, c, 0)` is `[c]`.

The result is the final shell returned by [`neighbors`](@ref)`(grid, c, k)`.
"""
function ring(grid::LevelGrid, c::A5Cell, k::Integer;
        connectivity::Connectivity=Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    steps == 0 && return A5Cell[c]
    shells = _shells(grid, c, steps, connectivity)
    # No shell exists after the traversal exhausts the connected component.
    steps <= length(shells) || return A5Cell[]
    return @inbounds shells[steps]
end

# ===========================================================================
# Breadth-first shells
# ===========================================================================

# Shell `j` contains cells at distance `j`, already wound for both callers.
function _shells(grid::LevelGrid, c::A5Cell, steps::Int, connectivity::Connectivity)
    shells = Vector{Vector{A5Cell}}()
    seen = Set{UInt64}((c.id,))
    frontier = UInt64[c.id]
    centre = cell_centroid(grid, c)
    # Use one ring-1 reference direction for every shell.
    frame = nothing
    for _ in 1:steps
        next = UInt64[]
        for x in frontier
            for y in _native_neighbors(grid, A5Cell(x), connectivity)
                y in seen && continue
                push!(seen, y)
                push!(next, y)
            end
        end
        shell = [A5Cell(y) for y in next]
        if frame === nothing && !isempty(shell)
            frame = _spoke_frame(grid, centre, shell)
        end
        frame === nothing || _wind!(shell, grid, centre, frame)
        push!(shells, shell)
        isempty(shell) && break
        frontier = next
    end
    return shells
end

# ===========================================================================
# Rotational shell order
#
# ===========================================================================

# Tangent frame whose zero direction points at the smallest ring-1 id.
#
# Store the anchor azimuth so a small negative rounding error cannot wrap it to
# the end of the ordered shell.
function _spoke_frame(grid::LevelGrid, centre, shell::AbstractVector{A5Cell})
    anchor = cell_centroid(grid, minimum(shell))
    e1, e2 = _tangent_basis(centre, anchor)
    return (e1, e2, _azimuth(centre, e1, e2, anchor))
end

# Order one shell counter-clockwise about `centre`, from the frame's spoke.
# Exact ties go to the smaller canonical id, so the result is total.
#
# Compute each projected-centroid key once before sorting.
function _wind!(shell::AbstractVector{A5Cell}, grid::LevelGrid, centre, frame)
    length(shell) <= 1 && return shell
    e1, e2, spoke = frame
    turn = 2 * Float64(pi)
    keyed = [(mod(_azimuth(centre, e1, e2, cell_centroid(grid, d)) - spoke, turn), d)
             for d in shell]
    sort!(keyed)
    for i in eachindex(shell, keyed)
        @inbounds shell[i] = keyed[i][2]
    end
    return shell
end

# Right-handed tangent basis with `e1` toward a neighbour and
# `e2 = centre × e1`. A neighbour seed remains defined at polar cell centres.
function _tangent_basis(centre, toward)
    u = (toward[1] - centre[1], toward[2] - centre[2], toward[3] - centre[3])
    radial = u[1] * centre[1] + u[2] * centre[2] + u[3] * centre[3]
    t = (u[1] - radial * centre[1], u[2] - radial * centre[2], u[3] - radial * centre[3])
    n = sqrt(t[1]^2 + t[2]^2 + t[3]^2)
    # Use a deterministic tangent direction when the reference has no tangent
    # component.
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
