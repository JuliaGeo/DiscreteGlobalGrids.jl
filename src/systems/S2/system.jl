# ---------------------------------------------------------------------------
# S2 system interface. At level `l`, scaffold ordinals are `0:6*4^l-1` and
# complete-grid position is ordinal plus one. Two Hilbert bits per level make
# parent `p ÷ 4`, children `4p:4p+3`, and subtrees contiguous.
#
#   * `cellat` is closed-form (`point_to_xyf`) — no tree descent.
#   * `descendant_range` is `[p*4^Δ, (p+1)*4^Δ)` shifted into position space.
#   * `node_extent` is the EXACT four-corner cap, not the inflated default.
#   * `neighbors` / `ring` walk the lattice and the seam table, not the geometry.
#   * `ancestor` drops `2Δ` bits in one shift.
#
# Native `s2_cellid` reindexing is unavailable because compatibility has not
# been verified against s2geometry fixtures.
# ---------------------------------------------------------------------------

# ===========================================================================
# Types
# ===========================================================================

"""
    S2System() <: AbstractHierarchicalGridSystem

The [S2](https://s2geometry.io) discrete global grid system on the unit sphere:
six cube-face charts, each refined by aperture-4 quadrant subdivision, with
cells ordered along a Hilbert curve within each face.

Level `l` has `6 * 4^l` cells. S2 is not equal-area; the quadratic `ST → UV`
projection gives an approximately 2.08× area spread. Every cell edge is an
exact great-circle arc, so [`cell_boundary`](@ref) is the exact four-vertex
cell boundary.

# Ids

The canonical [`LevelIndex`](@ref) stores the 0-based scaffold ordinal
`face * 4^level + hilbert_position`; complete-grid position is `index + 1`.
Native 64-bit `s2_cellid` is not an alternate scheme.

# Levels

`0:30`. Level 30 has `6 * 4^30 = 6917529027641081856` cells; level 31 would
overflow `Int64`.

# Traits

`has_sorted_subtrees` is `true`. [`max_neighbors`](@ref) is 8 for `Vertex()` and
4 for `Edge()`. [`node_extent`](@ref) uses the exact four-corner subtree cap.
"""
struct S2System <: DGG.AbstractHierarchicalGridSystem end

# `levelgrid(S2System(), l)` is the package's `HierarchicalLevelGrid`: all
# `6 * 4^l` cells in scaffold ordinal (face-major, Hilbert-within-face) order.
# S2's fast paths hang off this alias, and the five primitives it forwards to
# are the `(sys, ...)` methods further down.
const LevelGrid = DGG.HierarchicalLevelGrid{S2System}

"""
    MAX_LEVEL

The deepest S2 level, 30 — s2geometry's leaf level, and the last one whose cell
count `6 * 4^level` fits in a signed 64-bit integer.
"""
const MAX_LEVEL = 30

# Callers validate levels before these unchecked conversions.
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

"""
    rootcells(S2System())

The six cube faces, as level-0 cells `LevelIndex(0, 0:5)` in face order — the
s2geometry face numbering `+x, +y, +z, -x, -y, -z` (see `FACE_NORMAL`).
"""
DGG.rootcells(::S2System) = [DGG.LevelIndex(0, i) for i in 0:5]

"""
    parent(S2System(), c) -> LevelIndex

The scaffold-ordinal parent `index ÷ 4`, one level up. Throws `ArgumentError`
for a level-0 cell.
"""
function Base.parent(::S2System, c::DGG.LevelIndex)
    l = DGG.level(c)
    l > 0 || throw(ArgumentError(
        "level-0 S2 cell $c is a root and has no parent"))
    return DGG.LevelIndex(l - 1, c.index >> 2)
end

"""
    children(S2System(), c)

The four children `4*index .+ (0:3)`, ascending at the next level. Throws
`ArgumentError` at `max_level`.
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

The ancestor at level `l`: `index >> 2Δ`, equivalent to applying
[`parent`](@ref) `Δ` times.
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

The range is exact and hole-free; sibling ranges partition the parent's range.
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

Reads the dense, subtree-contiguous [`descendant_range`](@ref) as consecutive
ids.
"""
function DGG.descendants(sys::S2System, c::DGG.LevelIndex, l::Integer)
    r = DGG.descendant_range(sys, c, l)      # validates `l` both ways
    target = Int(l)
    return [DGG.LevelIndex(target, i - 1) for i in r]
end

# ===========================================================================
# The level grid: size, and positions <-> ids
# ===========================================================================

DGG.ncells(::S2System, l::Integer) = Int(_ncells(l))

# The grid bounds-checks `i`, so this is the bijection and nothing else.
DGG.cellindex(::S2System, l::Integer, i::Int) = DGG.LevelIndex(l, i - 1)

"""
    cellposition(S2System(), c) -> Union{Int,Nothing}

Returns `index + 1` for an in-range ordinal, otherwise `nothing`. The grid has
already rejected a cell from another level.
"""
function DGG.cellposition(::S2System, c::DGG.LevelIndex)
    0 <= c.index < _ncells(DGG.level(c)) || return nothing
    return Int(c.index + 1)
end

# Validate the ordinal against its own level. `hilbert_to_xyf` would otherwise
# un-Hilbert an ordinal no cell has, yielding the geometry of a cell that does
# not exist.
@inline function _checked_index(c::DGG.LevelIndex)
    l = DGG.level(c)
    0 <= c.index < _ncells(l) || throw(ArgumentError(
        "scaffold ordinal $(c.index) is out of range 0:$(_ncells(l) - 1) at level $l"))
    return c.index
end

# The grid-level form, which additionally pins the cell to this grid's level.
@inline function _checked_index(g::LevelGrid, c::DGG.LevelIndex)
    DGG.level(c) == g.level || throw(ArgumentError(
        "cell $c is at level $(DGG.level(c)), not the grid's level $(g.level)"))
    return _checked_index(c)
end

# ===========================================================================
# Geometry
# ===========================================================================

"""
    cell_boundary(grid, c) -> Vector{UnitSphericalPoint}

The exact four-vertex boundary ring, counter-clockwise seen from outside and
implicitly closed.

The ring starts at the `(s+, t+)` corner and runs `(s+, t+)`, `(s-, t+)`,
`(s-, t-)`, `(s+, t-)`, following `cell_corners`.

Chart lines map to great-circle arcs, so the four-corner polygon is the cell and
requires no densification. Generic `cell_area` therefore returns its true area.
"""
function DGG.cell_boundary(::S2System, c::DGG.LevelIndex)
    nside = _nside(DGG.level(c))
    ix, iy, face = hilbert_to_xyf(_checked_index(c), nside)
    corners = cell_corners(ix, iy, face, nside)
    return GO.UnitSphericalPoint{Float64}[corners[1], corners[2], corners[3], corners[4]]
end

"""
    cell_centroid(grid, c) -> UnitSphericalPoint

The S2 cell centre: the chart evaluated at the lattice cell's midpoint
`((ix + 0.5)/nside, (iy + 0.5)/nside)`.

This matches `S2CellId::ToPoint` and is strictly interior. It is not the
spherical quadrilateral's area centroid.
"""
function DGG.cell_centroid(::S2System, c::DGG.LevelIndex)
    nside = _nside(DGG.level(c))
    ix, iy, face = hilbert_to_xyf(_checked_index(c), nside)
    return cell_center(ix, iy, face, nside)
end

# ===========================================================================
# node_extent — the exact four-corner cap
# ===========================================================================

# Bounding cap for an S2 cell and its subtree. Children tile their parent, and
# each geodesically convex cell is the hull of its four corners. The maximum
# corner radius is `acos(1/√3) < π/2`, so the cap is convex. `nextfloat` includes
# boundary corners despite rounding.
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

The cell's four-corner cap, used without [`cap_inflation`](@ref).

Children tile their parent, and each cell is the geodesic convex hull of its
corners, so the cap covers the whole subtree.

The maximum radius is `acos(1/√3) ≈ 0.9553` rad, below `π/2`, so every extent is
geodesically convex.
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

The cell containing `p`, computed by [`point_to_xyf`](@ref)'s analytic inverse.

Never `nothing`: a complete S2 level grid covers the sphere.

**Ties.** Boundary points use deterministic face and lattice rules, each
self-consistent — the returned cell's own centroid maps back to it:

  - the **face** is the axis of `p`'s largest-magnitude component, ties broken
    toward the lower axis in `x, y, z` order and zero counting as positive
    ([`xyz_to_face`](@ref));
  - the **lattice cell** is chosen by `floor`, which puts a point on a cut line
    on the higher side of it.
"""
function DGG.cellat(g::LevelGrid, p::GO.UnitSphericalPoint)
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
function _one_ring(g::LevelGrid, c::DGG.LevelIndex, connectivity::DGG.Connectivity)
    _checked_index(g, c)
    out = SmallVector{8,DGG.LevelIndex}()
    for h in lattice_neighbors(c.index, g.level, connectivity)
        out = SmallCollections.push(out, DGG.LevelIndex(g.level, h))
    end
    return out
end

"""
    neighbors(grid, c, k = 1; connectivity = Vertex())

The cells within `k` lattice steps of `c`, excluding `c`. Rings `1:k` are
concatenated outward, each counter-clockwise seen from outside.

So `ring(grid, c, k)` is exactly the trailing block of
`neighbors(grid, c, k)`, and `neighbors(grid, c, k)` is
`vcat(ring(grid, c, 1), ..., ring(grid, c, k))`.

# Connectivity

`Vertex()` uses the 3×3 lattice neighbourhood: 8 cells, or 7 at the 24 face
corners where no fourth diagonal cell exists. `Edge()` keeps the four axis
steps; diagonal steps share only a vertex.

At level 0 both connectivities give the four faces sharing a cube edge; the
three faces at each corner are already pairwise edge-adjacent.

# Order

**Counter-clockwise seen from outside, starting at the `+s` lattice
direction** — the cycle `+s, +s+t, +t, -s+t, -s, -s-t, -t, +s-t` under
`Vertex()` and its restriction `+s, +t, -s, -t` under `Edge()`. `+s` increases
`ix`. The orientation-preserving chart carries this order to the sphere.
Missing neighbours are omitted without leaving an ordering gap.

Rings beyond the first are ordered by azimuth about the cell centre, from the
spoke through the first ring-1 neighbour; see [`ring`](@ref).

`k == 0` returns an empty container; `k == 1` returns a
`SmallCollections.SmallVector` sized by `max_neighbors`.
"""
function DGG.neighbors(g::LevelGrid, c::DGG.LevelIndex, k::Integer = 1;
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

`k == 1` uses the lattice cycle. For `k >= 2`, cells are sorted by azimuth about
the centre, counter-clockwise from the first ring-1 neighbour. Azimuth ties use
canonical id.
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
function _sort_ccw!(cells::Vector{DGG.LevelIndex}, g::LevelGrid, c::DGG.LevelIndex,
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
