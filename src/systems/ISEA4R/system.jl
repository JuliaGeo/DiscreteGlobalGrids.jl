# ISEA4R uses 0-based Morton identifiers over ten diamonds. Parent/child and
# subtree operations are radix-4 arithmetic; grid positions are identifier + 1.

# ===========================================================================
# Types
# ===========================================================================

"""
    ISEA4RSystem() <: AbstractHierarchicalGridSystem

Ten equal-area ISEA rhombus charts refined by aperture-4 subdivision. Level
`l in 0:29` contains `10*4^l` cells of solid angle `4π/(10*4^l)`.

Canonical [`LevelIndex`](@ref) values store
`diamond*4^level + morton(ix,iy)`. Morton order gives contiguous subtrees.
Vertex connectivity has at most 9 neighbors and edge connectivity at most 4.
Chart edges are curved and [`cell_boundary`](@ref) densifies them.
"""
struct ISEA4RSystem <: DGG.AbstractHierarchicalGridSystem end

# Grid descriptor for all `10 * 4^l` cells in Morton order.
const LevelGrid = DGG.HierarchicalLevelGrid{ISEA4RSystem}

# Convert a validated level to its diamond side and cell count.
@inline _nside(level::Integer) = Int64(1) << Int(level)
@inline _ncells(level::Integer) = 10 * (Int64(1) << (2 * Int(level)))

"The deepest level at which `10 * 4^level` still fits a signed 64-bit integer."
const MAX_LEVEL = 29

"The number of diamonds — the ten rhombi of `diamonds.jl`, and the root count."
const NDIAMONDS = 10

# ===========================================================================
# System interface
# ===========================================================================

DGG.cellindextype(::ISEA4RSystem) = DGG.LevelIndex
DGG.levels(::ISEA4RSystem) = 0:MAX_LEVEL
DGG.has_sorted_subtrees(::ISEA4RSystem) = true

"""
    max_neighbors(ISEA4RSystem(), connectivity) -> Int

`9` under `Vertex()`, `4` under `Edge()`.

Interior cells have eight vertex neighbours. At vertices 0 and 11, five
diamond corners meet and a cell can have nine. Other valence-3 icosahedron
vertices give seven. At level zero every diamond has six vertex neighbours.
"""
DGG.max_neighbors(::ISEA4RSystem, ::DGG.Vertex) = 9
DGG.max_neighbors(::ISEA4RSystem, ::DGG.Edge) = 4

# Row-major codecs remain internal; only the canonical Morton index is exposed.
DGG.cellindextypes(::ISEA4RSystem) = (DGG.LevelIndex,)

"""
    rootcells(ISEA4RSystem())

The ten level-0 cells: one per diamond, `LevelIndex(0, 0:9)`. At level 0 the
Morton code is empty, so the id *is* the diamond number.
"""
DGG.rootcells(::ISEA4RSystem) = [DGG.LevelIndex(0, d) for d in 0:(NDIAMONDS - 1)]

"""
    parent(ISEA4RSystem(), c) -> LevelIndex

The Morton parent: `index ÷ 4`, one level up. Throws an `ArgumentError` on a
level-0 cell, which has no parent.
"""
function Base.parent(::ISEA4RSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    l > 0 || throw(ArgumentError(
        "level-0 ISEA4R cell $c is a root and has no parent"))
    return DGG.LevelIndex(l - 1, c.index >> 2)
end

"""
    children(ISEA4RSystem(), c)

The four Morton children `4*index .+ (0:3)`, one level down, ascending.

ISEA4R refinement is a uniform quadtree, so every cell has exactly four
children. Icosahedron vertices affect adjacency but not subdivision.

In `(ix, iy)` terms the four are `(2ix, 2iy)`, `(2ix+1, 2iy)`, `(2ix, 2iy+1)`,
`(2ix+1, 2iy+1)` in that order, since `morton(2ix + a, 2iy + b) == 4*morton(ix,
iy) + a + 2b`. Throws an `ArgumentError` at `max_level`.
"""
function DGG.children(sys::ISEA4RSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    l < DGG.max_level(sys) || throw(ArgumentError(
        "ISEA4R cell $c is at max_level $(DGG.max_level(sys)) and has no children"))
    base = c.index << 2
    return [DGG.LevelIndex(l + 1, base + k) for k in 0:3]
end

"""
    ancestor(ISEA4RSystem(), c, l) -> LevelIndex

The ancestor at level `l`, in one shift: `index >> 2Δ`.

Sound because the id is `diamond * 4^level + morton(ix, iy)` and the Morton code
is positional — dropping the low `2Δ` bits drops `Δ` bits from each of `ix` and
`iy`, which is exactly `Δ` steps up the quadtree, and the diamond term divides
through untouched.
"""
function DGG.ancestor(sys::ISEA4RSystem, c::DGG.LevelIndex, l::Integer)
    target = Int(l)
    lc = DGG.level(c)
    target <= lc || throw(ArgumentError(
        "ancestor level $target is deeper than the cell's own level $lc"))
    target >= 0 || throw(ArgumentError(
        "ancestor level $target is above the root level 0"))
    return DGG.LevelIndex(target, c.index >> (2 * (lc - target)))
end

"""
    descendant_range(ISEA4RSystem(), c, l) -> UnitRange{Int}

The contiguous **positions** in `levelgrid(sys, l)` occupied by `c`'s level-`l`
descendants: `index * 4^Δ` through `(index + 1) * 4^Δ - 1` in 0-based Morton
ids, shifted into 1-based positions.

Exact and hole-free in both directions, which is what
`has_sorted_subtrees(sys) == true` asserts: Morton order is depth-first curve
order by construction, every position in the range names a real cell (ISEA4R
has no id gaps — ten dense diamonds, no pentagons), and sibling ranges tile the
parent's exactly.
"""
function DGG.descendant_range(sys::ISEA4RSystem, c::DGG.LevelIndex, l::Integer)
    target = Int(l)
    lc = DGG.level(c)
    target >= lc || throw(ArgumentError(
        "descendant level $target is above the cell's own level $lc"))
    target <= DGG.max_level(sys) || throw(ArgumentError(
        "descendant level $target is past max_level $(DGG.max_level(sys))"))
    shift = 2 * (target - lc)
    lo = c.index << shift
    hi = ((c.index + 1) << shift) - 1
    return Int(lo + 1):Int(hi + 1)
end

"""
    descendants(ISEA4RSystem(), c, l)

Every level-`l` descendant of `c`, ascending.

Morton ids are dense and subtree-contiguous, so this is
[`descendant_range`](@ref) read off as consecutive ids, with no `children`
recursion and no sort.
"""
function DGG.descendants(sys::ISEA4RSystem, c::DGG.LevelIndex, l::Integer)
    r = DGG.descendant_range(sys, c, l)      # validates `l` both ways
    target = Int(l)
    return [DGG.LevelIndex(target, i - 1) for i in r]
end

# ===========================================================================
# The level grid: size, and positions <-> ids
# ===========================================================================

DGG.ncells(::ISEA4RSystem, l::Integer) = Int(_ncells(l))

# The grid bounds-checks `i`, so this is the bijection and nothing else.
DGG.cellindex(::ISEA4RSystem, l::Integer, i::Int) = DGG.LevelIndex(l, i - 1)

"""
    cellposition(ISEA4RSystem(), c) -> Union{Int,Nothing}

Return `index + 1` for an in-range id, or `nothing` otherwise. The grid must
reject cells from another level first.
"""
function DGG.cellposition(::ISEA4RSystem, c::DGG.LevelIndex)
    0 <= c.index < _ncells(DGG.level(c)) || return nothing
    return Int(c.index + 1)
end

# Validate a Morton id before decoding it into diamond coordinates.
@inline function _checked_index(c::DGG.LevelIndex)
    l = DGG.level(c)
    0 <= c.index < _ncells(l) || throw(ArgumentError(
        "ISEA4R id $(c.index) is out of range 0:$(_ncells(l) - 1) at level $l"))
    return c.index
end

# The grid-level form the topology entry points use, which additionally pins the
# cell to the grid it was handed to.
@inline function _checked_index(g::LevelGrid, c::DGG.LevelIndex)
    DGG.level(c) == g.level || throw(ArgumentError(
        "cell $c is at level $(DGG.level(c)), not the grid's level $(g.level)"))
    return _checked_index(c)
end

# ===========================================================================
# Geometry
# ===========================================================================

# Eight great-circle segments approximate each curved chart edge. Shared points
# are bit-identical within a diamond; cross-diamond incidence requires tolerance.
#
# Hashed coordinate incidence is valid only within one diamond. Cross-diamond
# boundaries require a tolerance, and boundary coordinates do not contain
# negative zero.
const BOUNDARY_SEGMENTS = 8

"""
    _perimeter_points(ix, iy, diamond, nside, nseg) -> Vector{UnitSphericalPoint}

The cell's boundary walked in `nseg` equal chart steps per edge, starting at the
`(x+, y+)` corner and running in the [`cell_corners`](@ref) order
`(x+,y+) → (x-,y+) → (x-,y-) → (x+,y-)`. `nseg == 1` reproduces
[`cell_corners`](@ref) exactly.

Implicitly closed — each edge contributes its start vertex and its interior
points, never its end vertex, so the next edge's start is not duplicated.
"""
function _perimeter_points(ix::Integer, iy::Integer, diamond::Integer,
        nside::Integer, nseg::Integer)
    n = nside
    x0 = Int64(ix)
    y0 = Int64(iy)
    pts = Vector{GO.UnitSphericalPoint{Float64}}(undef, 4 * nseg)
    k = 0
    for i in 0:(nseg - 1)          # (x+,y+) -> (x-,y+), along y = (iy+1)/n
        t = i / nseg
        pts[k += 1] = xyd_to_point((x0 + 1 - t) / n, (y0 + 1) / n, diamond)
    end
    for i in 0:(nseg - 1)          # (x-,y+) -> (x-,y-), along x = ix/n
        t = i / nseg
        pts[k += 1] = xyd_to_point(x0 / n, (y0 + 1 - t) / n, diamond)
    end
    for i in 0:(nseg - 1)          # (x-,y-) -> (x+,y-), along y = iy/n
        t = i / nseg
        pts[k += 1] = xyd_to_point((x0 + t) / n, y0 / n, diamond)
    end
    for i in 0:(nseg - 1)          # (x+,y-) -> (x+,y+), along x = (ix+1)/n
        t = i / nseg
        pts[k += 1] = xyd_to_point((x0 + 1) / n, (y0 + t) / n, diamond)
    end
    return pts
end

"""
    cell_boundary(grid, c) -> Vector{UnitSphericalPoint}

Return the implicitly closed boundary counterclockwise from `(x+,y+)`, seen from
outside the sphere. Each curved chart edge uses `BOUNDARY_SEGMENTS`
great-circle segments. Use
[`cell_area`](@ref) for the exact equal-area value.
"""
function DGG.cell_boundary(::ISEA4RSystem, c::DGG.LevelIndex)
    nside = _nside(DGG.level(c))
    ix, iy, d = morton_to_xyd(_checked_index(c), nside)
    return _perimeter_points(ix, iy, d, nside, BOUNDARY_SEGMENTS)
end

"""
    cell_area(grid, c) -> Float64

Exact cell area in steradians: `4π/(10*4^level)`. This `O(1)` value is
independent of the approximate boundary polygon.
"""
DGG.cell_area(g::LevelGrid, c::DGG.LevelIndex) =
    (_checked_index(g, c); 4 * Float64(π) / _ncells(g.level))

"""
    cell_centroid(grid, c) -> UnitSphericalPoint

The cell centre: the chart evaluated at the lattice cell's midpoint
`((ix + 0.5)/nside, (iy + 0.5)/nside)`.

This is the centre *by definition* — the chart is exactly equal-area, so the
midpoint of the chart rectangle is the canonical centre — and it is strictly
interior, as [`cell_centroid`](@ref) requires. It is not the spherical centroid
of the published 4-gon; no equal-area DGGS claims that of its cell centres.
"""
function DGG.cell_centroid(::ISEA4RSystem, c::DGG.LevelIndex)
    nside = _nside(DGG.level(c))
    ix, iy, d = morton_to_xyd(_checked_index(c), nside)
    return cell_center(ix, iy, d, nside)
end

# ===========================================================================
# node_extent — the subtree cap
# ===========================================================================

# How many chart samples per edge the subtree cap is built from. See
# `_subtree_cap` for why this number, and not the boundary's, sets the bound.
const CAP_EDGE_SEGMENTS = 8

"""
    _subtree_cap(ix, iy, diamond, nside) -> SphericalCap

Bounding cap for a cell and its subtree. Children tile the parent chart
rectangle exactly, so a cap covering that rectangle covers all descendants.
The radius is the sampled perimeter maximum plus half the largest sample gap
and one outward ULP. Sampling only the perimeter suffices because the distance
from the centre is maximised there — in fact at one of the four corners, all of
which are samples.
"""
function _subtree_cap(ix::Integer, iy::Integer, diamond::Integer, nside::Integer)
    center = cell_center(ix, iy, diamond, nside)
    pts = _perimeter_points(ix, iy, diamond, nside, CAP_EDGE_SEGMENTS)
    rmax = 0.0
    gap = 0.0
    prev = pts[end]
    for p in pts
        rmax = max(rmax, US.spherical_distance(center, p))
        gap = max(gap, US.spherical_distance(prev, p))
        prev = p
    end
    radius = min(Float64(π), rmax + gap / 2)
    return SphericalCap(center, nextfloat(radius))
end

"""
    node_extent(ISEA4RSystem(), c) -> SphericalCap

Return the uninflated subtree cap. Exact chart nesting makes the cell rectangle
a bound for every descendant; [`_subtree_cap`](@ref) supplies a conservative
sampled radius. All returned caps are geodesically convex: the widest in the
system is a level-0 diamond's, 62.3°, against the 90° bound.
"""
function DGG.node_extent(::ISEA4RSystem, c::DGG.LevelIndex)
    l = DGG.level(c)
    nside = _nside(l)
    ix, iy, d = morton_to_xyd(c.index, nside)
    return _subtree_cap(ix, iy, d, nside)
end

# ===========================================================================
# Location
# ===========================================================================

"""
    cellat(grid, p::UnitSphericalPoint) -> LevelIndex

Return the cell containing `p` via the analytic chart inverse. Complete grids
never return `nothing`. Boundary ties use Snyder's face choice and the
higher-side lattice cell, deterministically per floating-point platform.
"""
DGG.cellat(g::LevelGrid, p::GO.UnitSphericalPoint) =
    DGG.LevelIndex(g.level, point_to_morton(p, _nside(g.level)))

# ===========================================================================
# Topology
# ===========================================================================

"""
    _one_ring(grid, c, connectivity) -> SmallVector{9,LevelIndex}

The immediate neighbours of `c` in **counter-clockwise rotational order seen
from outside the sphere, starting at the `(+1, 0)` chart direction**, from
[`lattice_neighbors`](@ref).

The subject's id is validated once, by `_checked_index`; the neighbours are then
encoded through [`xyd_to_morton_unchecked`](@ref), because
[`lattice_neighbors`](@ref) derives them from that validated cell by table
lookup and they cannot be off the lattice. Re-deriving `ispow2(nside)` and two
range tests per neighbour is the whole cost of this function otherwise, and this
is the inner loop of every breadth-first shell walk.
"""
function _one_ring(g::LevelGrid, c::DGG.LevelIndex, connectivity::DGG.Connectivity)
    nside = _nside(g.level)
    ix, iy, d = morton_to_xyd(_checked_index(g, c), nside)
    out = SmallVector{9,DGG.LevelIndex}()
    for (jx, jy, jd) in lattice_neighbors(ix, iy, d, nside, connectivity)
        out = SmallCollections.push(out,
            DGG.LevelIndex(g.level, xyd_to_morton_unchecked(jx, jy, jd, nside)))
    end
    return out
end

"""
    neighbors(grid, c, k = 1; connectivity = Vertex())

Cells within `k` lattice steps, excluding `c`, with rings concatenated outward.
Each ring is counterclockwise seen from outside the sphere, starting at chart
direction `(+1,0)` — a chart direction, not a compass one, so the same slot
points a different way on each diamond. Vertex connectivity has degree 7–9 at
corner cells and 8 in the interior; edge connectivity uses
the four axis offsets. Missing corner slots are omitted and multi-cell vertex
slots follow fan order. Outer rings use the same starting azimuth.

At level zero each diamond has six vertex-neighbors; the 7–9 counts apply to
corner cells at finer levels.

`k == 0` returns an empty container; `k == 1` returns a
`SmallCollections.SmallVector` sized by [`max_neighbors`](@ref).
"""
function DGG.neighbors(g::LevelGrid, c::DGG.LevelIndex, k::Integer = 1;
        connectivity::DGG.Connectivity = DGG.Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    steps == 0 && return SmallVector{9,DGG.LevelIndex}()
    steps == 1 && return _one_ring(g, c, connectivity)
    shells = _shells(g, c, steps, connectivity)
    isempty(shells) && return DGG.LevelIndex[]
    return reduce(vcat, shells)
end

"""
    ring(grid, c, k; connectivity = Vertex())

Cells at lattice distance exactly `k`, counterclockwise from the first ring-1
direction. `k == 0` returns `[c]`; outer-ring azimuth ties use canonical order.
"""
function DGG.ring(g::LevelGrid, c::DGG.LevelIndex, k::Integer;
        connectivity::DGG.Connectivity = DGG.Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    steps == 0 && return DGG.LevelIndex[c]
    shells = _shells(g, c, steps, connectivity)
    steps <= length(shells) || return DGG.LevelIndex[]
    return shells[steps]
end

# Breadth-first expansion over the lattice one-ring; shell `j` is the set at
# distance exactly `j`, each returned in CCW rotational order. Shared by
# `neighbors` and `ring` so the two cannot disagree about what a shell is or
# what order it is in.
function _shells(g::LevelGrid, c::DGG.LevelIndex, steps::Int,
        connectivity::DGG.Connectivity)
    shells = Vector{DGG.LevelIndex}[]
    steps == 0 && return shells
    seen = Set{DGG.LevelIndex}((c,))
    frontier = DGG.LevelIndex[c]
    for j in 1:steps
        next = DGG.LevelIndex[]
        for x in frontier
            for y in _one_ring(g, x, connectivity)
                y in seen && continue
                push!(seen, y)
                push!(next, y)
            end
        end
        # Shell 1 is already the lattice cycle, in order. Outer shells come out
        # of the breadth-first walk in whatever order the frontier happened to
        # reach them, so they are wound geometrically, from the spoke shell 1
        # starts on.
        if j > 1 && !isempty(first(shells))
            _sort_ccw!(next, g, c, first(first(shells)))
        end
        push!(shells, next)
        isempty(next) && break
        frontier = next
    end
    return shells
end

# Order a shell counter-clockwise about `c`'s centre, from the azimuth of the
# first ring-1 neighbour. `e1 × e2 == centre` makes `(e1, e2)` right-handed SEEN
# FROM OUTSIDE, which is what puts increasing `atan(u·e2, u·e1)`
# counter-clockwise from outside rather than from inside.
function _sort_ccw!(cells::Vector{DGG.LevelIndex}, g::LevelGrid,
        c::DGG.LevelIndex, reference::DGG.LevelIndex)
    length(cells) <= 1 && return cells
    centre = DGG.cell_centroid(g, c)
    e1, e2 = _tangent_basis(centre, DGG.cell_centroid(g, reference))
    # Family-wide observation, not a defect: an outer-ring cell whose azimuth
    # falls within floating-point noise BELOW the starting spoke sorts to the
    # END of the ring, because `mod(-ε, 2π)` is `~2π`. That is a start rotation
    # only — the single-CCW-cycle law the contract states is cyclic and holds
    # either way — and HEALPix and IGeo7 share the construction and the
    # behaviour. Nothing here depends on which side of the spoke such a cell
    # lands, so it is left alone rather than nudged by a tolerance that would
    # itself need a documented width.
    key(x) = begin
        p = DGG.cell_centroid(g, x)
        (mod(_azimuth(centre, e1, e2, p), 2 * Float64(π)), x)
    end
    sort!(cells; by = key)
    return cells
end

# A right-handed-from-outside tangent basis at `centre`, with `e1` pointing at
# `toward` (projected into the tangent plane).
function _tangent_basis(centre, toward)
    u = (toward[1] - centre[1], toward[2] - centre[2], toward[3] - centre[3])
    dot = u[1] * centre[1] + u[2] * centre[2] + u[3] * centre[3]
    t = (u[1] - dot * centre[1], u[2] - dot * centre[2], u[3] - dot * centre[3])
    n = sqrt(t[1]^2 + t[2]^2 + t[3]^2)
    # A degenerate reference (a neighbour whose centroid projects onto the cell
    # centre) cannot happen for a real cell, but a fallback keeps this total.
    n <= eps(Float64) && (t = abs(centre[3]) < 0.9 ? (0.0, 0.0, 1.0) : (1.0, 0.0, 0.0);
                          n = 1.0)
    e1 = (t[1] / n, t[2] / n, t[3] / n)
    e2 = (centre[2] * e1[3] - centre[3] * e1[2],
          centre[3] * e1[1] - centre[1] * e1[3],
          centre[1] * e1[2] - centre[2] * e1[1])
    return e1, e2
end

function _azimuth(centre, e1, e2, p)
    u = (p[1] - centre[1], p[2] - centre[2], p[3] - centre[3])
    return atan(u[1] * e2[1] + u[2] * e2[2] + u[3] * e2[3],
                u[1] * e1[1] + u[2] * e1[2] + u[3] * e1[3])
end
