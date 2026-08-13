# ---------------------------------------------------------------------------
# Adjacency
#
# THE ORDER IS ROTATIONAL, NOT SORTED. `neighbors(grid, c, k)` is rings
# 1, 2, ..., k concatenated outward, and each ring runs counter-clockwise seen
# from outside the sphere. So `ring(grid, c, k)` is exactly the tail block of
# `neighbors(grid, c, k)`, and a stencil consumer can read a weight vector
# straight off the result.
#
# WHERE THE SET COMES FROM. `A5Native._get_global_cell_neighbors` — a5's own
# adjacency, a within-quintant lattice step unioned with the quintant-seam,
# face-seam, apex and corner special cases. It is exact and it is fast; what it
# is not is ordered. It ends in `sort!(collect(::Set))`, i.e. ASCENDING BY ID,
# which is precisely the order the interface forbids. Nothing about a5's
# internal walk order is usable either: the set is assembled from two
# independently computed families, so there is no cycle to recover from it.
#
# WHERE THE ORDER COMES FROM, therefore: geometry. Each shell is sorted by
# azimuth about `cell_centroid(grid, c)` in a right-handed tangent frame, which
# is counter-clockwise seen from outside BY CONSTRUCTION — a sorted sequence of
# angles in `[0, 2pi)` wraps exactly once, which is the winding law's
# definition. This is the same scheme the generic geometric fallback uses, and
# the suite re-derives the winding with an independently written frame rather
# than restating this one.
#
# THE START. Ring 1 begins at the neighbour with the SMALLEST `A5Cell` id, and
# every outer shell begins on that same spoke — the member whose azimuth,
# measured counter-clockwise from the ring-1 anchor's direction, is smallest.
# The anchor is an integer comparison, so the start is exactly reproducible
# across platforms even though the ordering after it is floating point. The
# choice is arbitrary in the way the interface says a start may be (a lattice
# direction would be better, but A5's lattice direction changes meaning at every
# quintant seam and would not be uniform either); what matters is that it is
# stated, and `test/systems/A5/runtests.jl` pins one cell's literal sequence so
# that a silent rotation cannot pass.
#
# CONNECTIVITY IS REAL HERE. Unlike the icosahedral hex systems, A5's `Vertex()`
# and `Edge()` name different sets — see `max_neighbors` in `system.jl` for the
# measurement. `edge_only = true` selects von Neumann; the a5 default, `false`,
# is Moore, and is the interface's default too.
# ---------------------------------------------------------------------------

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

# ALLOCATION — KNOWN, TRACKED, NOT FIXED HERE. A `k = 1` call allocates 3.0 KB
# at res 2 and 6.4 KB from res 6 down (it plateaus there, since the deep levels
# all take the same walk), nearly all of it inside
# `A5Native._get_global_cell_neighbors`: it
# accumulates candidates in a `Set{UInt64}` with intermediate vectors and hands
# back `sort!(collect(set))`. No amount of work on this side recovers that — the
# set is built and thrown away before this module is handed anything, so the
# `SmallVector` in `neighbors` is copying an allocation that already happened.
#
# THE FIX, when someone takes it: the degree is bounded by 11 (see
# `max_neighbors`), so the walk in `native.jl` needs no `Set` at all. Give it a
# fixed-capacity scratch buffer — a `SmallVector{11,UInt64}` of candidates with a
# linear dedup pass, which beats hashing at that size anyway — and let it fill a
# caller-supplied destination instead of returning a fresh vector. That makes
# `k = 1` allocation-free, matching `children`.
#
# Deferred deliberately rather than overlooked: `native.jl` is a verbatim
# carry-over of the pre-redesign arithmetic, kept diff-clean against its source
# so the port stays auditable. Changing its hot loop is a separate, testable
# piece of work and wants its own before/after numbers.
#
# (Keep this note ABOVE the docstring. A comment between a docstring and the
# function it documents silently DETACHES it — the docstring stops registering
# for the method and `@doc` falls back to the interface's generic one.)
"""
    neighbors(grid::LevelGrid, c::A5Cell, k = 1; connectivity = Vertex())

The cells within `k` grid steps of `c`, excluding `c`, in **rotational order**:
rings `1..k` concatenated outward, each ring counter-clockwise seen from outside
the sphere, every ring starting on the same spoke.

This makes `ring(grid, c, k)` the tail block of `neighbors(grid, c, k)`, and it
is what lets a stencil index into the result positionally. Ring 1 starts at the
neighbour with the **smallest [`A5Cell`](@ref) id**; see this file's header for
why the order is measured rather than read off a lattice.

`connectivity` genuinely changes the answer on A5 — the default
[`Vertex()`](@ref Vertex) adds the corner-only neighbours that
[`Edge()`](@ref Edge) leaves out, 1 to 3 of them below level 1 and 8 of them at
level 1 — which is unlike the two icosahedral hexagonal systems in this package.
See [`max_neighbors`](@ref) for the degrees.

# Container

`k <= 1` returns a `SmallVector{11,A5Cell}` — sized by the `Vertex()` bound
under both connectivities, so the type does not depend on a keyword. `k >= 2`
returns a `Vector{A5Cell}`, since a disc outgrows any static bound. The return
type is therefore a two-way union across `k`, with the boundary at `k = 1` /
`k = 2`.

Throws an `ArgumentError` for a cell that is not a valid cell of this grid's
resolution — there is no `nothing` in this contract to return, and answering
for a neighbouring cell instead would be worse.

Allocates: a few kilobytes even at `k = 1`, inside a5's own adjacency walk
rather than here. See the note above this docstring for the numbers, the cause
and the fix.
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

# Shell `j` of the result is the cells at adjacency distance exactly `j`,
# already wound. The winding happens HERE rather than in the two callers,
# because that is what makes `ring(grid, c, k)` the literal tail block of
# `neighbors(grid, c, k)`: both read the same shells, so they cannot disagree
# element for element.
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
# The frame's third entry is the anchor's own measured azimuth rather than the
# literal `0.0`. In exact arithmetic they are the same number; in floating point
# the anchor's tangential component can round to a hair below zero, and
# `mod(-eps, 2pi)` is `2pi` — which would sort the spoke's own cell to the END
# of its ring. Subtracting the measurement cancels that exactly.
function _spoke_frame(grid::LevelGrid, centre, shell::AbstractVector{A5Cell})
    anchor = cell_centroid(grid, minimum(shell))
    e1, e2 = _tangent_basis(centre, anchor)
    return (e1, e2, _azimuth(centre, e1, e2, anchor))
end

# Order one shell counter-clockwise about `centre`, from the frame's spoke.
# Exact ties go to the smaller canonical id, so the result is total.
#
# The keys are computed once per cell and sorted alongside them, rather than
# handed to `sort!` as a `by` function: `cell_centroid` is a native projection
# round trip, and a `by` key is recomputed at every comparison, which would
# spend O(n log n) projections to order n cells.
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

# A tangent basis at `centre` with `e1` pointing at `toward` (projected into the
# tangent plane) and `e2 = centre x e1`. Then `e1 x e2 == centre`, i.e. the pair
# is right-handed SEEN FROM OUTSIDE the sphere, which is what makes increasing
# `atan(u.e2, u.e1)` run counter-clockwise viewed from above the plane rather
# than from below it.
#
# Seeding from a neighbour's direction rather than from local east is not a
# stylistic choice: two of A5's twelve res-0 cells are centred exactly on a
# geographic pole, where east and north do not exist, and a frame built from
# them would emit NaNs and corrupt an order rather than fail.
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
