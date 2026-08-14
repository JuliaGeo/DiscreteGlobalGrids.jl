# A5 supplies unordered adjacency sets. Shells are therefore sorted by azimuth
# counter-clockwise in a tangent frame. Ring 1 starts at the smallest id; outer
# rings use the same spoke. `edge_only=true` implements `Edge()` connectivity.

# The `Vertex()` bound, which is also the container capacity for both
# connectivities so that `neighbors` has one concrete return type at k <= 1.
const MAX_NEIGHBORS = 11

_edge_only(::Vertex) = false
_edge_only(::Edge) = true

# a5's adjacency for one cell, as raw ids. Guarded, because `deserialize` raises
# a `BoundsError` out of `ORIGINS` for an id that names no face and the lattice
# walk would otherwise answer for a cell that does not exist.
function _native_neighbors(grid::LevelGrid, c::A5Cell, connectivity::Connectivity)
    level(c) == grid.level || throw(ArgumentError(
        "A5 cell $c is at resolution $(level(c)), not this grid's $(grid.level)"))
    isvalid(c) || throw(ArgumentError("A5 cell $c is not a valid cell"))
    return A5Native._get_global_cell_neighbors(c.id; edge_only=_edge_only(connectivity))
end

# The native adjacency walk allocates internally through a `Set{UInt64}` even
# for `k == 1`: 3.0 KB at res 2, plateauing at 6.4 KB from res 6 down.
"""
    neighbors(grid::LevelGrid, c::A5Cell, k = 1; connectivity = Vertex())

Cells within `k` grid steps, excluding `c`, as counter-clockwise shells
concatenated outward. Ring 1 starts at the smallest [`A5Cell`](@ref) id and all
outer rings use the same spoke. [`ring`](@ref) is the final shell.

[`Vertex()`](@ref Vertex) includes corner-only neighbours; [`Edge()`](@ref Edge)
does not. See [`max_neighbors`](@ref).

`k <= 1` returns a `SmallVector{11,A5Cell}` — sized by the `Vertex()` bound
under both connectivities. `k >= 2` returns a `Vector{A5Cell}`.

Throws `ArgumentError` unless `c` is valid at the grid resolution. The native
adjacency walk allocates internally.
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

Identical, element for element, to the last `length` entries of
[`neighbors`](@ref)`(grid, c, k)` — the two share one breadth-first walk and one
winding step, so the disc really is its shells concatenated rather than merely
agreeing with them as a set.
"""
function ring(grid::LevelGrid, c::A5Cell, k::Integer;
        connectivity::Connectivity=Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    steps == 0 && return A5Cell[c]
    shells = _shells(grid, c, steps, connectivity)
    # A walk that ran out of cells before reaching `steps` has an empty shell
    # there: the ring is genuinely empty, not missing.
    steps <= length(shells) || return A5Cell[]
    return @inbounds shells[steps]
end

# ===========================================================================
# The breadth-first walk both entry points are reads of
# ===========================================================================

# Shell `j` contains cells at distance `j`, already wound for both callers.
function _shells(grid::LevelGrid, c::A5Cell, steps::Int, connectivity::Connectivity)
    shells = Vector{Vector{A5Cell}}()
    seen = Set{UInt64}((c.id,))
    frontier = UInt64[c.id]
    centre = cell_centroid(grid, c)
    # Fixed once, from ring 1, and reused by every outer shell: one spoke for
    # the whole disc is what makes position `j` of a ring mean a direction.
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
# Written out here rather than borrowed from `Fallbacks`: a system module may
# not reach into the fallback substrate, and the arithmetic is six lines.
# ===========================================================================

# The frame every shell's azimuth is measured in, built once from ring 1: the
# tangent basis at `centre` whose zero direction points at the ring-1 neighbour
# with the smallest canonical id.
#
# Store the measured anchor azimuth so a negative rounding error cannot rotate
# the anchor to the end via `mod(-eps, 2pi)`.
function _spoke_frame(grid::LevelGrid, centre, shell::AbstractVector{A5Cell})
    anchor = cell_centroid(grid, minimum(shell))
    e1, e2 = _tangent_basis(centre, anchor)
    return (e1, e2, _azimuth(centre, e1, e2, anchor))
end

# Order one shell counter-clockwise about `centre`, from the frame's spoke.
# Exact ties go to the smaller canonical id, so the result is total.
#
# Cache keys to avoid recomputing projected centroids during comparisons.
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
    # A neighbour whose centroid coincides with the subject's, or sits exactly
    # antipodal to it, names no tangential direction. Neither can happen in a
    # real tessellation, but an arbitrary-but-valid basis keeps this total.
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
