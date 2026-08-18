# ---------------------------------------------------------------------------
# S2 system interface. At level `l`, scaffold ordinals are `0:6*4^l-1` and
# complete-grid position is ordinal plus one. Two Hilbert bits per level make
# parent `p ÷ 4`, children `4p:4p+3`, and subtrees contiguous — so the hierarchy,
# the level-grid arithmetic, and the subtree engines are the quad-face family's,
# and this file writes only what the cube-face chart decides:
#
#   * `cellat` is closed-form (`point_to_xyf`) — no tree descent.
#   * `node_extent` is the exact four-corner cap.
#   * `neighbors` / `ring` walk the lattice and the seam table, not the geometry.
#
# Native `s2_cellid` reindexing is unavailable because compatibility has not
# been verified against s2geometry fixtures.
# ---------------------------------------------------------------------------

# ===========================================================================
# Types
# ===========================================================================

"""
    S2System() <: AbstractQuadFaceGridSystem

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

`has_sorted_subtrees` is `true`. [`maxneighbors`](@ref) is 8 for `Vertex()` and
4 for `Edge()`. [`node_extent`](@ref) uses the exact four-corner subtree cap.
[`subtree_border`](@ref) is an `O(border)` walk over the subtree's square block.
"""
struct S2System <: DGG.AbstractQuadFaceGridSystem end

# Grid descriptor for all `6 * 4^l` cells in face-major Hilbert order.
const LevelGrid = DGG.HierarchicalLevelGrid{S2System}

"""
    MAX_LEVEL

The deepest S2 level, 30 — s2geometry's leaf level, and the last one whose cell
count `6 * 4^level` fits in a signed 64-bit integer.
"""
const MAX_LEVEL = 30

# ===========================================================================
# System interface
# ===========================================================================

# The six cube faces, which `rootcells` names `LevelIndex(0, 0:5)` in the
# s2geometry face order `+x, +y, +z, -x, -y, -z` (see `FACE_NORMAL`).
DGG.nbasefaces(::S2System) = 6
DGG.systemname(::S2System) = "S2"
DGG.idname(::S2System) = "scaffold ordinal"

DGG.levels(::S2System) = 0:MAX_LEVEL

DGG.maxneighbors(::S2System, ::DGG.Vertex) = 8
DGG.maxneighbors(::S2System, ::DGG.Edge) = 4

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
function DGG.cell_boundary(sys::S2System, c::DGG.LevelIndex)
    nside = DGG.nside(DGG.level(c))
    ix, iy, face = hilbert_to_xyf(DGG.checked_id(sys, c), nside)
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
function DGG.cell_centroid(sys::S2System, c::DGG.LevelIndex)
    nside = DGG.nside(DGG.level(c))
    ix, iy, face = hilbert_to_xyf(DGG.checked_id(sys, c), nside)
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
    nside = DGG.nside(DGG.level(c))
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
    nside = DGG.nside(g.level)
    ix, iy, face = point_to_xyf(p, nside)
    return DGG.LevelIndex(g.level, xyf_to_hilbert(ix, iy, face, nside))
end

# ===========================================================================
# Topology
# ===========================================================================

"""
    one_ring(grid, c, connectivity) -> SmallVector{8,LevelIndex}

The immediate neighbours of `c` in **counter-clockwise rotational order seen
from outside the sphere, starting at the `+s` lattice direction**
([`NEIGHBOR_OFFSETS`](@ref)), with cube-corner steps and repeats dropped.
"""
function DGG.one_ring(g::LevelGrid, c::DGG.LevelIndex, connectivity::DGG.Connectivity)
    DGG.checked_id(g, c)
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
`SmallCollections.SmallVector` sized by `maxneighbors`.
"""
function DGG.neighbors(g::LevelGrid, c::DGG.LevelIndex, k::Integer = 1;
        connectivity::DGG.Connectivity = DGG.Vertex())
    steps = DGG.checked_steps(k)
    steps == 0 && return SmallVector{8,DGG.LevelIndex}()
    steps == 1 && return DGG.one_ring(g, c, connectivity)
    shells = DGG.adjacency_shells(g, c, steps, connectivity)
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
    steps = DGG.checked_steps(k)
    steps == 0 && return DGG.LevelIndex[c]
    steps == 1 && return DGG.one_ring(g, c, connectivity)
    shells = DGG.adjacency_shells(g, c, steps, connectivity)
    steps <= length(shells) || return DGG.LevelIndex[]
    return shells[steps]
end
