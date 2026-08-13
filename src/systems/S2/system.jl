# ---------------------------------------------------------------------------
# The S2 system on the grid interface
#
# S2 is the interface's second DENSE-ORDINAL system, and structurally the
# closest relative of nested HEALPix: at level `l` the cells are exactly
# `0:6*4^l - 1`, cell `p`'s children are `4p:4p+3`, and its parent is `p ÷ 4`.
# Every hierarchy method below is that arithmetic, and nothing here needs a
# lookup table.
#
# ## The id, and why `÷ 4` is legal
#
# The canonical id is the interface's own `LevelIndex(level, index)`, where
#
#     index  ==  the 0-BASED SCAFFOLD ORDINAL `face * 4^level + hilbert_position`
#     position in `levelgrid(sys, level)`  ==  index + 1
#
# The Hilbert position is built two bits per level, most significant first
# (`chart.jl`), so dropping its low two bits steps exactly one level up. The
# face term `face * 4^level` divides through by 4 untouched, so the SAME
# statement holds of the whole ordinal — which is what makes `÷ 4` the parent
# map, `4p + k` the children, and subtrees contiguous runs of positions
# (`has_sorted_subtrees == true`).
#
# That claim is checked by exhaustion over levels 0-6 in
# `test/systems/S2/runtests.jl` rather than inferred from the Hilbert tables.
# It is the one property the whole hierarchy rests on, and the tables it would
# follow from are transcribed constants.
#
# ## What is a fast path here
#
#   * `cellat` is closed-form (`point_to_xyf`) — no tree descent.
#   * `descendant_range` is `[p*4^Δ, (p+1)*4^Δ)` shifted into position space.
#   * `node_extent` is the EXACT four-corner cap, not the inflated default.
#   * `neighbors` / `ring` walk the lattice and the seam table, not the geometry.
#   * `ancestor` drops `2Δ` bits in one shift.
#
# ## What is deliberately NOT here
#
# The native 64-bit `s2_cellid` (face bits, Hilbert bits, lsb sentinel) as an
# alternate scheme via `reindex`. It is one codec away — `xyf_to_hilbert`
# already transcribes `S2CellId::FromFaceIJ`'s tables — but this repository has
# no s2geometry fixtures, so shipping the conversion would mean publishing an
# interoperability claim that nothing checks. The scaffold ordinal is canonical
# either way, so adding the scheme later is additive.
# ---------------------------------------------------------------------------

# ===========================================================================
# Types
# ===========================================================================

"""
    S2System() <: AbstractHierarchicalGridSystem

The [S2](https://s2geometry.io) discrete global grid system on the unit sphere:
six cube-face charts, each refined by aperture-4 quadrant subdivision, with
cells ordered along a Hilbert curve within each face.

Level `l` has `6 * 4^l` cells. S2 is **not** equal-area — the quadratic
`ST → UV` projection narrows the within-level area spread to about 2.08×, and
what it buys instead is that every cell edge is an exact great-circle arc, so
[`cell_boundary`](@ref) is a four-vertex ring that *is* the cell rather than a
densification of it.

# Ids

The canonical cell index is [`LevelIndex`](@ref), whose `index` field holds the
**0-based scaffold ordinal** `face * 4^level + hilbert_position`; a cell's
position in `levelgrid(sys, l)` is `index + 1`. The native 64-bit `s2_cellid`
is not offered as an alternate scheme — see the note in `system.jl`.

# Levels

`0:30`. The bound is s2geometry's own leaf level, and it is also where the
`Int64` cell count runs out: level 30 has `6 * 4^30 = 6917529027641081856`
cells and level 31 would overflow.

# Traits

`has_sorted_subtrees` is `true` — the Hilbert scaffold ordinal nests exactly,
so a subtree occupies a contiguous run of positions and
[`descendant_range`](@ref) is exact and hole-free. [`max_neighbors`](@ref) is 8
under `Vertex()` and 4 under `Edge()`. [`node_extent`](@ref) is overridden with
the cell's own exact four-corner cap — S2 children tile their parent exactly,
so nothing needs inflating and [`cap_inflation`](@ref) is never consulted.
"""
struct S2System <: DGG.AbstractHierarchicalGridSystem end

"""
    S2Grid(level) <: AbstractGrid

The complete S2 grid at refinement `level`: all `6 * 4^level` cells in scaffold
ordinal (face-major, Hilbert-within-face) order. Built by
[`levelgrid`](@ref); a lightweight descriptor, not a materialised cell list.
"""
struct S2Grid <: DGG.AbstractGrid
    level::Int
end

"""
    MAX_LEVEL

The deepest S2 level, 30 — s2geometry's leaf level, and the last one whose cell
count `6 * 4^level` fits in a signed 64-bit integer.
"""
const MAX_LEVEL = 30

# The one place `nside` is derived from a level, and the one place the cell
# count is. Both guard nothing: callers reach them only after `levels(sys)` has
# been checked, and a shift of 62 is what a bad level would produce.
@inline _nside(level::Integer) = Int64(1) << Int(level)
@inline _ncells(level::Integer) = 6 * (Int64(1) << (2 * Int(level)))

# ===========================================================================
# System interface
# ===========================================================================

DGG.cellindextype(::S2System) = DGG.LevelIndex
DGG.levels(::S2System) = 0:MAX_LEVEL
DGG.has_sorted_subtrees(::S2System) = true

DGG.max_neighbors(::S2System, ::DGG.Vertex) = 8
DGG.max_neighbors(::S2System, ::DGG.Edge) = 4

function DGG.levelgrid(sys::S2System, l::Integer)
    lvl = Int(l)
    lvl in DGG.levels(sys) || throw(ArgumentError(
        "level $lvl is outside $(DGG.levels(sys)) for $(nameof(typeof(sys)))"))
    return S2Grid(lvl)
end

"""
    rootcells(S2System())

The six cube faces, as level-0 cells `LevelIndex(0, 0:5)` in face order — the
s2geometry face numbering `+x, +y, +z, -x, -y, -z` (see `FACE_NORMAL`).
"""
DGG.rootcells(::S2System) = [DGG.LevelIndex(0, i) for i in 0:5]

"""
    parent(S2System(), c) -> LevelIndex

The scaffold-ordinal parent: `index ÷ 4`, one level up. Throws an
`ArgumentError` on a level-0 cell, which has no parent.

Sound because the ordinal is `face * 4^level + hilbert_position` and the Hilbert
position is positional two bits per level: `÷ 4` drops one level's worth of
Hilbert bits and divides the face term through exactly.
"""
function Base.parent(::S2System, c::DGG.LevelIndex)
    l = DGG.level(c)
    l > 0 || throw(ArgumentError(
        "level-0 S2 cell $c is a root and has no parent"))
    return DGG.LevelIndex(l - 1, c.index >> 2)
end

"""
    children(S2System(), c)

The four children `4*index .+ (0:3)`, one level down, ascending.

Always exactly four: S2 refinement is a uniform quadtree over the six face
charts with no exceptional cells, so unlike the icosahedral systems this count
never varies. Throws an `ArgumentError` at `max_level`.
"""
function DGG.children(sys::S2System, c::DGG.LevelIndex)
    l = DGG.level(c)
    l < DGG.max_level(sys) || throw(ArgumentError(
        "S2 cell $c is at max_level $(DGG.max_level(sys)) and has no children"))
    base = c.index << 2
    return [DGG.LevelIndex(l + 1, base + k) for k in 0:3]
end

"""
    ancestor(S2System(), c, l) -> LevelIndex

The ancestor at level `l`, in one shift: `index >> 2Δ` — [`parent`](@ref)
applied `Δ` times, with the same soundness argument.
"""
function DGG.ancestor(sys::S2System, c::DGG.LevelIndex, l::Integer)
    target = Int(l)
    lc = DGG.level(c)
    target <= lc || throw(ArgumentError(
        "ancestor level $target is deeper than the cell's own level $lc"))
    target >= 0 || throw(ArgumentError(
        "ancestor level $target is above the root level 0"))
    return DGG.LevelIndex(target, c.index >> (2 * (lc - target)))
end

"""
    descendant_range(S2System(), c, l) -> UnitRange{Int}

The contiguous **positions** in `levelgrid(sys, l)` occupied by `c`'s level-`l`
descendants: `index * 4^Δ` through `(index + 1) * 4^Δ - 1` in 0-based scaffold
ordinals, shifted into 1-based positions.

Exact and hole-free in both directions, which is what
`has_sorted_subtrees(sys) == true` asserts: the Hilbert scaffold ordinal is
depth-first curve order by construction, every position in the range names a
real cell (S2 has no id gaps), and sibling ranges tile the parent's exactly.
"""
function DGG.descendant_range(sys::S2System, c::DGG.LevelIndex, l::Integer)
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
    descendants(S2System(), c, l)

Every level-`l` descendant of `c`, ascending.

Scaffold ordinals are dense and subtree-contiguous, so this is
[`descendant_range`](@ref) read off as consecutive ids, with no `children`
recursion and no sort.
"""
function DGG.descendants(sys::S2System, c::DGG.LevelIndex, l::Integer)
    r = DGG.descendant_range(sys, c, l)      # validates `l` both ways
    target = Int(l)
    return [DGG.LevelIndex(target, i - 1) for i in r]
end

# ===========================================================================
# Grid interface
# ===========================================================================

DGG.system(::S2Grid) = S2System()
DGG.level(g::S2Grid) = g.level
DGG.ncells(g::S2Grid) = Int(_ncells(g.level))

function DGG.cellindex(g::S2Grid, i::Int)
    1 <= i <= DGG.ncells(g) || throw(BoundsError(g, i))
    return DGG.LevelIndex(g.level, i - 1)
end

"""
    cellposition(grid, c) -> Union{Int,Nothing}

Closed form: `index + 1` for a cell at the grid's own level and in range, and
`nothing` otherwise (a different level, or an ordinal no cell has). Replaces
the fallback's linear scan.
"""
function DGG.cellposition(g::S2Grid, c::DGG.LevelIndex)
    DGG.level(c) == g.level || return nothing
    0 <= c.index < _ncells(g.level) || return nothing
    return Int(c.index + 1)
end

# The id guard every geometry entry point needs: `hilbert_to_xyf` validates the
# ordinal's range for us, but not that the cell belongs to THIS grid's level.
@inline function _checked_index(g::S2Grid, c::DGG.LevelIndex)
    DGG.level(c) == g.level || throw(ArgumentError(
        "cell $c is at level $(DGG.level(c)), not the grid's level $(g.level)"))
    0 <= c.index < _ncells(g.level) || throw(ArgumentError(
        "scaffold ordinal $(c.index) is out of range 0:$(_ncells(g.level) - 1) at level $(g.level)"))
    return c.index
end

# ===========================================================================
# Geometry
# ===========================================================================

"""
    cell_boundary(grid, c) -> Vector{UnitSphericalPoint}

The cell's boundary ring: **four vertices**, counter-clockwise seen from outside
the sphere, implicitly closed (the first vertex is not repeated).

The ring starts at the `(s+, t+)` corner and runs `(s+, t+)`, `(s-, t+)`,
`(s-, t-)`, `(s+, t-)` — `cell_corners`' lattice order, the same slot order as
HEALPix's `pixel_corners`.

**Four vertices is exact, not an approximation.** An S2 cell edge is a chart
line `u = const` or `v = const` on a cube face; lifted by the gnomonic map it is
the intersection of a plane through the origin with the sphere, i.e. a
great-circle arc. So the 4-gon *is* the cell and there is nothing to densify —
the one structural difference from HEALPix, whose chart-line edges bulge off the
geodesic and whose ring therefore carries 32 vertices.

Consequently `cell_area` needs no override: the generic spherical polygon area
of this ring is the cell's true area.
"""
function DGG.cell_boundary(g::S2Grid, c::DGG.LevelIndex)
    nside = _nside(g.level)
    ix, iy, face = hilbert_to_xyf(_checked_index(g, c), nside)
    corners = cell_corners(ix, iy, face, nside)
    return GO.UnitSphericalPoint{Float64}[corners[1], corners[2], corners[3], corners[4]]
end

"""
    cell_centroid(grid, c) -> UnitSphericalPoint

The S2 cell centre: the chart evaluated at the lattice cell's midpoint
`((ix + 0.5)/nside, (iy + 0.5)/nside)`.

This is S2's own definition of a centre — `S2CellId::ToPoint` is precisely the
ST-space midpoint pushed through `ST → UV → XYZ` and normalised — and it is
strictly interior, as [`cell_centroid`](@ref) requires. It is **not** the area
centroid of the spherical quadrilateral, and S2 does not claim it is: the chart
is not equal-area, so the two differ by a fraction of a cell.
"""
function DGG.cell_centroid(g::S2Grid, c::DGG.LevelIndex)
    nside = _nside(g.level)
    ix, iy, face = hilbert_to_xyf(_checked_index(g, c), nside)
    return cell_center(ix, iy, face, nside)
end

# ===========================================================================
# node_extent — the exact four-corner cap
# ===========================================================================

"""
    _cell_cap(ix, iy, face, nside) -> SphericalCap

The bounding cap of an S2 cell — *and therefore of its whole subtree*.

# Why the four corners are the whole answer

Two facts compose, and each is exact rather than sampled:

 1. **Children tile their parent exactly.** A child's chart rectangle is a
    quadrant of the parent's, and the chart is a homeomorphism of the closed
    unit square onto the face's spherical patch. So every descendant boundary
    point at every depth lies in the closed spherical quadrilateral of this
    cell. (This is what the aperture-7 icosahedral systems cannot say, and why
    the generic `node_extent` has to inflate.)
 2. **The cap through the four corners contains the cell.** The cell's edges
    are great-circle arcs (see [`cell_boundary`](@ref)) and the gnomonic map
    carries a convex planar rectangle to a geodesically convex patch, so the
    cell *is* the geodesic hull of its four corners. `rmax` is at most a
    level-0 face's `acos(1/√3) ≈ 0.9553 < π/2`, and a cap of radius below `π/2`
    is itself geodesically convex — so containing the four corners, it contains
    their hull, which is the cell.

The `π/2` bound is the load-bearing half of step 2 and not a stray remark: it
is what makes the cap convex, and hence what turns "contains the corners" into
"contains the cell". Note that it is a bound on the corner distance from the
cell's OWN centre; it says nothing about how far apart two points of the cell
may be (on a level-0 face, up to 1.91 rad), and it is the former the cap
argument needs.

The measured margin agrees: sweeping levels 0-3 with a 17×17 chart sampling of
each cell plus a 6-levels-deeper descendant lattice, the worst overshoot of
`rmax` is exactly `0.0` — corners of descendants land on the parent's corners
*bit-identically*, because `(ix + 1)/nside` and `(2ix + 2)/(2 nside)` are the
same exactly-representable dyadic and both divisions are exact.

`nextfloat` is one ulp of insurance on top of that, so a descendant corner sits
strictly inside rather than exactly on the rim — which keeps the conformance
harness's covering check off a floating-point coin toss without loosening its
`atol`.
"""
function _cell_cap(ix::Integer, iy::Integer, face::Integer, nside::Integer)
    centre = cell_center(ix, iy, face, nside)
    rmax = 0.0
    for p in cell_corners(ix, iy, face, nside)
        rmax = max(rmax, US.spherical_distance(centre, p))
    end
    return SphericalCap(centre, nextfloat(rmax))
end

"""
    node_extent(S2System(), c) -> SphericalCap

The cell's own exact four-corner cap — an override of the generic inflated
default, and the reason [`cap_inflation`](@ref) is never consulted for this
system.

S2 children tile their parent exactly and an S2 cell is the geodesic convex hull
of its four corners, so the cap through those corners already bounds the whole
subtree at every depth. See `_cell_cap` for the two-step argument and for the
measured margin.

The widest extent in the system is a level-0 face at 0.9553 rad ≈ 54.7°, well
inside the 90° that makes a cap geodesically convex. That one bound does double
duty: it is what makes the cap contain the cell in the first place (`_cell_cap`
step 2), and it is what makes the conformance harness's vertex-sampling proxy
for the covering law sound — so `require_convex_extents` needs no opting out.
"""
function DGG.node_extent(::S2System, c::DGG.LevelIndex)
    nside = _nside(DGG.level(c))
    ix, iy, face = hilbert_to_xyf(c.index, nside)
    return _cell_cap(ix, iy, face, nside)
end

# ===========================================================================
# Location
# ===========================================================================

"""
    cellat(grid, p::UnitSphericalPoint) -> LevelIndex

The cell containing `p`, in closed form — [`point_to_xyf`](@ref), the chart's
analytic inverse, with no tree descent and no point-in-polygon test.

Never `nothing`: a complete S2 level grid covers the sphere.

**Ties.** A point exactly on a cell boundary is legitimately contained by every
cell meeting there, and the tie is broken in two stages, both deterministic and
self-consistent (the returned cell's own centroid maps back to it):

  - the **face** is the axis of `p`'s largest-magnitude component, ties broken
    toward the lower axis in `x, y, z` order and a zero component counting as
    positive ([`xyz_to_face`](@ref)) — a pure comparison of three `Float64`s,
    so this half of the rule is bit-identical on every platform;
  - the **lattice cell** is chosen by `floor`, which puts a point on a cut line
    on the higher side of it.
"""
function DGG.cellat(g::S2Grid, p::GO.UnitSphericalPoint)
    nside = _nside(g.level)
    ix, iy, face = point_to_xyf(p, nside)
    return DGG.LevelIndex(g.level, xyf_to_hilbert(ix, iy, face, nside))
end

# ===========================================================================
# Topology
# ===========================================================================

"""
    _one_ring(grid, c, connectivity) -> SmallVector{8,LevelIndex}

The immediate neighbours of `c` in **counter-clockwise rotational order seen
from outside the sphere, starting at the `+s` lattice direction**
([`NEIGHBOR_OFFSETS`](@ref)), with cube-corner steps and repeats dropped.
"""
function _one_ring(g::S2Grid, c::DGG.LevelIndex, connectivity::DGG.Connectivity)
    _checked_index(g, c)
    out = SmallVector{8,DGG.LevelIndex}()
    for h in lattice_neighbors(c.index, g.level, connectivity)
        out = SmallCollections.push(out, DGG.LevelIndex(g.level, h))
    end
    return out
end

"""
    neighbors(grid, c, k = 1; connectivity = Vertex())

The cells within `k` lattice steps of `c`, excluding `c`, in **rotational
order**: the rings `1:k` concatenated outward, each ring counter-clockwise seen
from outside the sphere.

So `ring(grid, c, k)` is exactly the trailing block of
`neighbors(grid, c, k)`, and `neighbors(grid, c, k)` is
`vcat(ring(grid, c, 1), ..., ring(grid, c, k))`.

# Connectivity

`Vertex()` (the default) is the 3×3 lattice neighbourhood: **8 cells**, and
**7** at the 24 cells sitting in a face corner, where three cube faces meet and
the diagonal step has no fourth cell to name (`neighbors.jl`). `Edge()` keeps
the four axis steps, which share a whole cell edge; the four diagonal steps
share a single corner. An S2 cell is an axis-aligned rectangle in its chart, so
that labelling is the obvious one — and unlike HEALPix, where the pixel is a
diamond rotated 45° against the lattice and the obvious labelling is exactly
backwards.

At **level 0** a cell is a whole cube face, and both connectivities give 4: the
four faces sharing a cube edge. The three faces meeting at each cube corner are
pairwise edge-adjacent already, so `Vertex()` adds nothing.

# Order

**Counter-clockwise seen from outside, starting at the `+s` lattice
direction** — the cycle `+s, +s+t, +t, -s+t, -s, -s-t, -t, +s-t` under
`Vertex()` and its restriction `+s, +t, -s, -t` under `Edge()`. `+s` is the
direction of increasing `ix`, i.e. of increasing chart coordinate `s`, and the
gnomonic chart is orientation-preserving, so counter-clockwise in the lattice
plane is counter-clockwise on the sphere with no reversal (contrast HEALPix,
whose reference compass tuple runs the other way). Absent neighbours drop out of
the cycle rather than leaving a hole, so a ring of seven is still in rotational
order.

Rings beyond the first are ordered by azimuth about the cell centre, from the
spoke through the first ring-1 neighbour; see [`ring`](@ref).

`k == 0` returns an empty container; `k == 1` returns a
`SmallCollections.SmallVector` sized by `max_neighbors`.
"""
function DGG.neighbors(g::S2Grid, c::DGG.LevelIndex, k::Integer = 1;
        connectivity::DGG.Connectivity = DGG.Vertex())
    steps = Int(k)
    steps >= 0 || throw(ArgumentError("k must be non-negative, got $steps"))
    steps == 0 && return SmallVector{8,DGG.LevelIndex}()
    steps == 1 && return _one_ring(g, c, connectivity)
    shells = _shells(g, c, steps, connectivity)
    isempty(shells) && return DGG.LevelIndex[]
    return reduce(vcat, shells)
end

"""
    ring(grid, c, k; connectivity = Vertex())

The cells at lattice distance **exactly** `k`, counter-clockwise seen from
outside the sphere. `ring(grid, c, 0)` is `[c]`.

`k == 1` is the lattice cycle (see [`neighbors`](@ref)). For `k >= 2` there is
no lattice cycle to read off — an outer ring crosses cube seams and may wrap a
cube corner — so the shell is ordered **by azimuth about the cell centre**,
measured counter-clockwise from the direction of the **first ring-1
neighbour**. That is the extension the interface's `neighbors` docstring
recommends for exactly this case, and it makes every ring start on the same
spoke, so slot `j` of ring 2 points the same way as slot 1 of ring 1 does.
Ties in azimuth break by canonical id, so the order is total and deterministic.
"""
function DGG.ring(g::S2Grid, c::DGG.LevelIndex, k::Integer;
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
function _shells(g::S2Grid, c::DGG.LevelIndex, steps::Int,
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
        # reach them, so they are wound geometrically — from the spoke through
        # shell 1's first entry, which is what keeps every ring on one start.
        j > 1 && _sort_ccw!(next, g, c, shells[1])
        push!(shells, next)
        isempty(next) && break
        frontier = next
    end
    return shells
end

# Order a shell counter-clockwise about `c`'s centre, from the azimuth of the
# first ring-1 neighbour. `e1 x e2 == centre` makes `(e1, e2)` right-handed SEEN
# FROM OUTSIDE, which is what puts increasing `atan(u.e2, u.e1)` counter-
# clockwise from outside rather than from inside.
function _sort_ccw!(cells::Vector{DGG.LevelIndex}, g::S2Grid, c::DGG.LevelIndex,
        ring1::Vector{DGG.LevelIndex})
    length(cells) <= 1 && return cells
    centre = DGG.cell_centroid(g, c)
    e1, e2 = _tangent_basis(centre)
    ref = isempty(ring1) ? 0.0 :
        _azimuth(centre, e1, e2, DGG.cell_centroid(g, first(ring1)))
    key(x) = begin
        p = DGG.cell_centroid(g, x)
        (mod(_azimuth(centre, e1, e2, p) - ref, 2 * Float64(π)), x)
    end
    sort!(cells; by = key)
    return cells
end

# A right-handed-from-outside tangent basis at `centre`. The seed axis is the
# one `centre` leans on least, so the Gram-Schmidt step stays well-conditioned
# everywhere, poles included.
function _tangent_basis(centre)
    ax = abs(centre[1]) <= abs(centre[2]) ?
        (abs(centre[1]) <= abs(centre[3]) ? (1.0, 0.0, 0.0) : (0.0, 0.0, 1.0)) :
        (abs(centre[2]) <= abs(centre[3]) ? (0.0, 1.0, 0.0) : (0.0, 0.0, 1.0))
    s = _dot3(ax, centre)
    t = (ax[1] - s * centre[1], ax[2] - s * centre[2], ax[3] - s * centre[3])
    n = sqrt(_dot3(t, t))
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
