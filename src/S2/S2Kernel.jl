# ---------------------------------------------------------------------------
# S2 operations-kernel wiring — dense enumeration + geometry
#
# This file connects `S2DGGS` to the package's geometry generics
# (`src/core/kernel.jl`) and to the user-facing `cell_polygon(::AbstractDGGS,
# level, id)` of `src/core/interface.jl`. Every method here goes through the
# *chart kernel* (`chart.jl`) — the same `cell_corners` / `cell_center` closed
# forms the dense face grid (`face_grid.jl`) evaluates, reached by way of
# `hilbert_to_xyf`. That sharing is the point: the DGGS-level cell and the
# face-grid cell for the same ordinal are LITERALLY the same four
# `Float64`-triples, not two evaluations that agree to a tolerance, and
# `test/S2/test_s2_kernel.jl` asserts it bitwise against
# `Trees.getcell(treeify(S2FaceGrid(2^level; ordering = HilbertOrder())), id + 1)`.
#
# ## The id these methods take, and why the hierarchy is NOT wired
#
# `id` here is the **scaffold ordinal** `face * 4^level + hilbert_position`
# recorded by the `S2DGGS` registry entry (`src/core/systems/s2.jl`) — a
# 0-based dense index over the `nside = 2^level` lattice, and exactly
# `HilbertOrder`'s data position minus one. It is **not** the native 64-bit
# `s2_cellid`, whose face/Hilbert/sentinel bit encoding this package has not
# ported (`S2CellId` orientation / `range_min` / `range_max` are the registry's
# recorded prerequisite).
#
# `canonical_index_name(S2DGGS()) === :s2_cellid`, so declaring
# `has_ordinal_ids` over the scaffold ordinals would misrepresent the canonical
# id space: the whole hierarchy / ordinal / pruning group would then answer in
# a coordinate system that is not the system's own, silently, with no error to
# catch it. So the hierarchy group stays unwired — `cell_children`,
# `cell_parent`, `cell_descendants`, `cell_to_ordinal`, `descendant_range` and
# `root_ids` all still throw `NotPortedError`, and
# `test/S2/test_s2_kernel.jl` pins that boundary.
#
# Recorded as now-possible, because it is one trait away: over the scaffold
# ordinals the entire radix-4 group is already exact and free. The Hilbert
# position drops its low two bits per level up (`xyf_to_hilbert`), so an
# ordinal's children really are `4p:4p + 3`, its ancestor at `parent_level` is
# `p ÷ 4^Δ`, and its `leaf_level` descendants really are the contiguous interval
# `[p * 4^Δ, (p + 1) * 4^Δ)` — the same arithmetic `has_ordinal_ids(HEALPixDGGS())
# = true` buys HEALPix, with `root_count = 6` and `radix = 4` already wired.
# What is missing is not code but a *decision*: whether `S2DGGS` ids mean
# scaffold ordinals or native `s2_cellid`s. Wiring the hierarchy before that
# decision would pin the wrong answer into a public interface, so it waits.
#
# ## Caps
#
# `cell_cap` is the exact four-corner cap rather than the generic
# 1.2-inflated center-to-vertex cap; the argument is at the method.
# `subtree_cap` is deliberately left alone: its generic fallback needs
# `cell_descendants` above level equality, so it answers exactly for
# `level == leaf_level` and throws `NotPortedError` otherwise — which is the
# honest state of a system with no ported hierarchy.
# ---------------------------------------------------------------------------

import ..DiscreteGlobalGrids as DGG
import ConservativeRegridding: Trees

# --------------------------------------------------------------------------
# Dense geometry enumeration
# --------------------------------------------------------------------------

# These methods deliberately do not opt S2 into `has_ordinal_ids`: they
# enumerate the scaffold ids already accepted by the geometry methods, but make
# no claim that those ids are native `s2_cellid`s or that the hierarchy kernel
# is wired. This is the minimal bridge from `1:num_cells` to cell polygons.
DGG.num_cells(system::DGG.S2DGGS, level::Integer) =
    Int64(DGG.leaf_count(system, level))

function DGG.ordinal_to_cell(system::DGG.S2DGGS, level::Integer, ordinal::Integer)
    total = DGG.num_cells(system, level)
    1 <= ordinal <= total || throw(DGG.OrdinalRangeError(
        DGG.system_name(system), Int(level), Int(ordinal), Int(total)))
    return DGG.cell_id_type(system)(ordinal - 1)
end

# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------

# `level` is the lattice level, so `nside = 2^level` and `hilbert_to_xyf`
# un-Hilberts the scaffold ordinal into the `(ix, iy, face)` coordinates the
# chart is written against. `hilbert_to_xyf` is also where id validity is
# decided: it throws `ArgumentError` for an id outside `0:6 * 4^level - 1`.
# Both calls are in `chart.jl`, which is s2geometry-free by construction — this
# is the single evaluation the face grid shares.
_cell_corners(level::Integer, id) =
    cell_corners(hilbert_to_xyf(id, 2^Int(level))..., 2^Int(level))

"""
    cell_boundary(::S2DGGS, level, id; closed=false)

The four cell corners in ring order, from the chart kernel
([`cell_corners`](@ref)).

`id` is the scaffold ordinal `face * 4^level + hilbert_position` (0-based, and
`HilbertOrder`'s data position minus one), **not** a native 64-bit `s2_cellid`;
`level` is the lattice level, `0 <= level <= 30`, so `nside = 2^level`. An id
outside `0:6 * 4^level - 1` is an `ArgumentError` from
[`hilbert_to_xyf`](@ref).

The ring is **counter-clockwise as seen from outside the sphere**, which is a
hard contract rather than a convention: the convex-clip kernel that computes
spherical intersections clips a clockwise ring to EMPTY, so a reversed ring
yields silent zero areas instead of an error. The order is `chart.jl`'s lattice
order — `(x+, y+), (x-, y+), (x-, y-), (x+, y-)` — identical to what the face
grid emits for the same cell, corner for corner and bit for bit.

Unlike HEALPix's, S2 cell edges **are** great-circle arcs (a chart line
`u = const` is a central plane section, see [`stf_to_point`](@ref)), so this
4-vertex ring is the cell *exactly* — there is no curvature between the corners
for the caps to absorb and nothing to densify.
"""
function DGG.cell_boundary(::DGG.S2DGGS, level::Integer, id; closed::Bool=false)
    corners = _cell_corners(level, id)
    points = Vector{eltype(corners)}(undef, closed ? 5 : 4)
    @inbounds for i in 1:4
        points[i] = corners[i]
    end
    closed && (@inbounds points[5] = corners[1])
    return points
end

"""
    cell_center(::S2DGGS, level, id)

Native cell center: the chart evaluated at the lattice cell's midpoint
([`cell_center`](@ref)), which is S2's own definition of a center —
`S2CellId::ToPoint` is precisely the ST-space midpoint pushed through
`ST → UV → XYZ` and normalised. Preferred over the kernel's normalized
boundary mean: exact, one chart evaluation instead of four, and identical to
the face grid's center for the same cell.

(It is *not* the centroid of the spherical quadrilateral, and S2 does not claim
it is — the chart is not equal-area.) `id` is the scaffold ordinal; see
[`cell_boundary`](@ref).
"""
DGG.cell_center(::DGG.S2DGGS, level::Integer, id) =
    cell_center(hilbert_to_xyf(id, 2^Int(level))..., 2^Int(level))

# Exact cell cap over the *chart* corners, replacing the generic
# `CELL_CAP_INFLATION` formula (max center-to-vertex distance × 1.2). The
# inflation exists for hierarchies where a parent does not geographically
# contain its children; S2 cells partition their parent exactly — children are
# sub-rectangles of the parent's chart rectangle and the edges are geodesics —
# so there is no overhang to leave room for, and the four-corner cap is both
# sound and tighter.
#
# Soundness is the FourCornerCap proof already recorded at
# `cap_policy(::S2FaceSystem)` in `face_grid.jl`, read on a 1x1 block: block
# edges are great-circle arcs, every boundary point is strictly within a
# quarter turn of the corner-mean centre, so the farthest boundary point from
# that centre is one of the four corners — the midpoints and the `1.0001` slack
# `circle_from_four_corners` applies are pure insurance for S2 (for HEALPix,
# whose edges bulge off the geodesic, they are load-bearing).
#
# Corner ORDER: `circle_from_four_corners` documents its argument as
# `(BL, TL, BR, TR)` and internally reorders to the CCW walk
# `(BL, BR, TR, TL)` before taking slerp midpoints of consecutive pairs. The
# chart ring is `(TR, TL, BL, BR)` in those terms, so it is permuted into the
# documented slots here — `(corners[3], corners[2], corners[4], corners[1])` —
# and consequently ALL FOUR sampled midpoints are true cell-edge midpoints.
# Two further consequences worth having: the cap is then bitwise the
# `FourCornerCap` node extent the face-grid layer would compute for the same
# cell as a 1x1 block (same function, same four points, same argument order),
# and the cap is strictly the tightest of the two orderings. This deliberately
# differs from `HealpixKernel.jl`, which passes its ring through unpermuted and
# documents the resulting two diagonal midpoints: those caps are bitwise-pinned
# by existing tests and cannot move, whereas nothing pins these yet.
DGG.cell_cap(::DGG.S2DGGS, level::Integer, id) =
    (c = _cell_corners(level, id); Trees.circle_from_four_corners((c[3], c[2], c[4], c[1]), ()))

# --------------------------------------------------------------------------
# The user-facing connection
# --------------------------------------------------------------------------

"""
    cell_polygon(::S2DGGS, level, id) -> GI.Polygon

Cell `(level, id)` as a closed 4-gon on the **unit sphere**, i.e.
[`cell_polygon_unitsphere`](@ref) — the interface-level entry point, which
before this wiring threw `NotPortedError` for S2.

`id` is the scaffold ordinal `face * 4^level + hilbert_position` (see
[`cell_boundary`](@ref)), not a native `s2_cellid`. The ring is
counter-clockwise as seen from outside the sphere and, because S2 cell edges
are geodesics, it *is* the cell rather than an approximation of it.

```julia
cell_polygon(S2DGGS(), 3, 17)                    # this method
S2.cell_polygon(hilbert_to_xyf(17, 8)..., 8)     # the same polygon, chart-side
```

Dense geometry enumeration (`num_cells` / `ordinal_to_cell`) is wired over the
same scaffold ids. The hierarchy group (`cell_children`, `cell_parent`,
`descendant_range`, ...) still throws `NotPortedError` pending the native
`s2_cellid` port — see the file header.
"""
DGG.cell_polygon(system::DGG.S2DGGS, level::Integer, id::Integer) =
    DGG.cell_polygon_unitsphere(system, level, id)
