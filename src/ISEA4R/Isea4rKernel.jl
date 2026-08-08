# ---------------------------------------------------------------------------
# ISEA4R operations-kernel wiring — dense enumeration + geometry
#
# This file connects `ISEA4RDGGS` to the package's geometry generics
# (`src/core/kernel.jl`) and to the user-facing `cell_polygon(::AbstractDGGS,
# level, id)` of `src/core/interface.jl`. Every method here goes through the
# *chart kernel* (`chart.jl`) — the same `cell_corners` / `cell_center` closed
# forms the dense diamond grid (`face_grid.jl`) evaluates, reached by way of
# `morton_to_xyd`. That sharing is the point: the DGGS-level cell and the
# face-grid cell for the same ordinal are LITERALLY the same four
# `Float64`-triples, not two evaluations that agree to a tolerance, and
# `test/ISEA4R/test_isea4r_kernel.jl` asserts it bitwise against
# `Trees.getcell(treeify(Isea4rFaceGrid(2^r; ordering = MortonOrder())), id + 1)`.
#
# ## The id these methods take
#
# `id` is the **`isea4r_ordinal`** `diamond * 4^r + morton_position` — 0-based,
# and exactly `MortonOrder`'s data position minus one. Unlike S2's scaffold
# ordinal this *is* the system's canonical index
# (`canonical_index_name(ISEA4RDGGS()) === :isea4r_ordinal`), with `position`
# pinned to the Z-order Morton code by milestone 2 (`ix` in the even bit
# positions, `iy` in the odd ones) in this package's own ten-diamond layout.
#
# **That layout has no external oracle behind it: no DGGAL,
# SphericalSpatialTrees or other external identifier compatibility is claimed
# and none may be inferred** — the caveat is the same one `diamonds.jl`, the
# registry entry and `docs/design/isea4r_diamond_layout.md` carry, and it
# applies to these ids verbatim.
# `level` is the refinement level `r`, `nside = 2^r`, `0 <= r <= 29`
# (the bound where `10 * nside^2` still fits `Int64`, matching
# `max_nside(Isea4rFaceSystem())`).
#
# ## What is NOT wired, and what it would cost
#
# The hierarchy / pruning group — `cell_children`, `cell_parent`,
# `cell_descendants`, `cell_to_ordinal`, `descendant_range`, `root_ids` — still
# throws `NotPortedError`, and
# `test/ISEA4R/test_isea4r_kernel.jl` pins that boundary. This is a scope line,
# not an obstacle: this milestone is geometry.
#
# Recorded as now-possible, because unlike S2 there is no open question of *which
# id space is canonical* — this one is. The radix-4 arithmetic is already exact
# over these ordinals: the Morton code drops its low two bits per level up, and
# the within-diamond lattice nesting is bit-exact (`fl(ix/n) === fl(2ix/2n)`, see
# `xyd_to_point`), so a cell's children really are `4p:4p + 3`, its ancestor is
# `p ÷ 4^Δ`, and its `leaf_level` descendants really are the contiguous interval
# `[p * 4^Δ, (p + 1) * 4^Δ)`. `root_count = 10` and `radix = 4` are wired, so
# `has_ordinal_ids(::ISEA4RDGGS) = true` — one line — would derive the whole
# group from `src/core/kernel.jl` exactly as it does for HEALPix. What that line
# additionally *buys* is `DGGSGrid` / `DGGSPartialGrid` / `subtree_grid` support;
# what it *costs* is the kernel-test battery those grids deserve (ordinal round
# trips, the two-sided `descendant_range` contract, subtree-cap containment) and
# a decision on `has_exact_subtree_cap`. Deferred as a unit, deliberately.
#
# ## Caps
#
# `cell_cap` is the exact four-corner cap rather than the generic
# 1.2-inflated center-to-vertex cap; the argument is at the method.
# `subtree_cap` is left alone: its generic fallback needs `cell_descendants`
# above level equality, so it answers exactly for `level == leaf_level` and
# throws `NotPortedError` otherwise — the honest state of a system whose
# hierarchy is not wired.
# ---------------------------------------------------------------------------

import ..DiscreteGlobalGrids as DGG
import ConservativeRegridding: Trees

# --------------------------------------------------------------------------
# Dense geometry enumeration
# --------------------------------------------------------------------------

# The canonical ISEA4R ids are dense Morton ordinals. Wire only the count and
# ordinal-to-id direction needed to enumerate geometry; the hierarchy trait and
# all parent/child operations remain deliberately untouched.
DGG.num_cells(system::DGG.ISEA4RDGGS, level::Integer) =
    Int64(DGG.leaf_count(system, level))

function DGG.ordinal_to_cell(system::DGG.ISEA4RDGGS, level::Integer, ordinal::Integer)
    total = DGG.num_cells(system, level)
    1 <= ordinal <= total || throw(DGG.OrdinalRangeError(
        DGG.system_name(system), Int(level), Int(ordinal), Int(total)))
    return DGG.cell_id_type(system)(ordinal - 1)
end

# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------

# `level` is the refinement level `r`, so `nside = 2^r` and `morton_to_xyd`
# de-interleaves the ordinal into the `(ix, iy, diamond)` coordinates the chart
# is written against. `morton_to_xyd` is also where id validity is decided: it
# throws `ArgumentError` for an id outside `0:10 * 4^r - 1`. Both calls are in
# `chart.jl` — this is the single evaluation the face grid shares.
_cell_corners(level::Integer, id) =
    cell_corners(morton_to_xyd(id, 2^Int(level))..., 2^Int(level))

"""
    cell_boundary(::ISEA4RDGGS, level, id; closed=false)

The four cell corners in ring order, from the chart kernel
([`cell_corners`](@ref)).

`id` is the `isea4r_ordinal` `diamond * 4^level + morton_position` (0-based,
and `MortonOrder`'s data position minus one); `level` is the refinement level
`r`, `0 <= r <= 29`, so `nside = 2^r`. An id outside `0:10 * 4^r - 1` is an
`ArgumentError` from [`morton_to_xyd`](@ref). The ten-diamond layout these ids
are written against is this package's own convention with no external oracle
behind it — see `docs/design/isea4r_diamond_layout.md` before reading any
DGGAL / SST / other external compatibility into them.

The ring is **counterclockwise as seen from outside the sphere**, which is a
hard contract rather than a convention: the convex-clip kernel that computes
spherical intersections clips a clockwise ring to EMPTY, so a reversed ring
yields silent zero areas instead of an error. The order is `chart.jl`'s lattice
order — `(x+, y+), (x-, y+), (x-, y-), (x+, y-)` — identical to what the
diamond grid emits for the same cell, corner for corner and bit for bit (same
chart functions, same seam-ownership rule).

Snyder ISEA maps the chart's straight lines to curves on the sphere, so — as
with HEALPix, and unlike S2 — the 4-vertex ring *approximates* the cell rather
than being it; densify through [`xyd_to_point`](@ref) along the edges if more
is needed, and note that the caps below bound the curvature separately.
"""
function DGG.cell_boundary(::DGG.ISEA4RDGGS, level::Integer, id; closed::Bool=false)
    corners = _cell_corners(level, id)
    points = Vector{eltype(corners)}(undef, closed ? 5 : 4)
    @inbounds for i in 1:4
        points[i] = corners[i]
    end
    closed && (@inbounds points[5] = corners[1])
    return points
end

"""
    cell_center(::ISEA4RDGGS, level, id)

Native cell center: the chart evaluated at the lattice cell's midpoint
([`cell_center`](@ref)). The rhombus chart is *exactly* equal-area
([`xyd_to_point`](@ref)), so the chart midpoint is the canonical center — the
same argument HEALPix's `pixel_center` rests on. Preferred over the kernel's
normalized boundary mean: exact, one chart evaluation instead of four, and
identical to the diamond grid's center for the same cell.

(It is not the spherical centroid of the 4-gon; no equal-area DGGS claims that
of its cell centers.) `id` is the `isea4r_ordinal`; see [`cell_boundary`](@ref).
"""
DGG.cell_center(::DGG.ISEA4RDGGS, level::Integer, id) =
    cell_center(morton_to_xyd(id, 2^Int(level))..., 2^Int(level))

# Exact cell cap over the *chart* corners, replacing the generic
# `CELL_CAP_INFLATION` formula (max center-to-vertex distance × 1.2). The
# inflation exists for hierarchies where a parent does not geographically
# contain its children; an ISEA4R cell is the union of its four children —
# their chart rectangles quarter its own, bit-exactly (`fl(ix/n) ===
# fl(2ix/2n)`, see `xyd_to_point`) — so there is no overhang to leave room for,
# and the four-corner cap is both sound and tighter.
#
# Soundness is the pre-registered measurement already recorded at
# `cap_policy(::Isea4rFaceSystem)` in `face_grid.jl`, read on its 1x1 case: a
# dense 17x17 sampling of every block's chart rectangle — every leaf block
# among them, on every diamond, at `nside ∈ (1, 2, 3, 4, 5, 8, 16)` — showed a
# worst PRE-SLACK overhang of exactly `0.0`, i.e. the farthest sampled point of
# a cell is always one of its own four corners, seam-straddling cells included.
# The measurement ships as a standing test in `test/ISEA4R/test_face_grid.jl`;
# `test/ISEA4R/test_isea4r_kernel.jl` re-runs the containment check against
# these caps specifically. If the standing assert ever fires, this method must
# go with `cap_policy`'s override.
#
# Corner ORDER: `circle_from_four_corners` documents its argument as
# `(BL, TL, BR, TR)` and internally reorders to the CCW walk
# `(BL, BR, TR, TL)` before taking slerp midpoints of consecutive pairs. The
# chart ring is `(TR, TL, BL, BR)` in those terms, so it is permuted into the
# documented slots here — `(corners[3], corners[2], corners[4], corners[1])` —
# and consequently ALL FOUR sampled midpoints are true cell-edge midpoints,
# which is what the curvature of the Snyder edges actually needs bounded. The
# cap is then also bitwise the `FourCornerCap` node extent the face-grid layer
# would compute for the same cell as a 1x1 block. This deliberately differs
# from `HealpixKernel.jl`, which passes its ring through unpermuted and
# documents the resulting two diagonal midpoints: those caps are bitwise-pinned
# by existing tests and cannot move, whereas nothing pins these yet.
DGG.cell_cap(::DGG.ISEA4RDGGS, level::Integer, id) =
    (c = _cell_corners(level, id); Trees.circle_from_four_corners((c[3], c[2], c[4], c[1]), ()))

# --------------------------------------------------------------------------
# The user-facing connection
# --------------------------------------------------------------------------

"""
    cell_polygon(::ISEA4RDGGS, level, id) -> GI.Polygon

Cell `(level, id)` as a closed 4-gon on the **unit sphere**, i.e.
[`cell_polygon_unitsphere`](@ref) — the interface-level entry point, which
before this wiring threw `NotPortedError` for ISEA4R.

`id` is the `isea4r_ordinal` `diamond * 4^level + morton_position` (see
[`cell_boundary`](@ref)), in this package's own ten-diamond layout, which no
external fixture pins. The ring is counterclockwise as seen from outside the
sphere and, Snyder edges not being great circles, approximates the cell to
4-gon accuracy.

```julia
cell_polygon(ISEA4RDGGS(), 3, 17)                  # this method
ISEA4R.cell_polygon(morton_to_xyd(17, 8)..., 8)    # the same polygon, chart-side
```

Dense geometry enumeration (`num_cells` / `ordinal_to_cell`) is wired over the
canonical Morton ordinals. The hierarchy group (`cell_children`, `cell_parent`,
`descendant_range`, ...) still throws `NotPortedError`. See the file header.
"""
DGG.cell_polygon(system::DGG.ISEA4RDGGS, level::Integer, id::Integer) =
    DGG.cell_polygon_unitsphere(system, level, id)
