# ---------------------------------------------------------------------------
# ISEA9R operations-kernel wiring — geometry only
#
# This file connects `ISEA9RDGGS` to the package's geometry generics
# (`src/core/kernel.jl`) and to the user-facing `cell_polygon(::AbstractDGGS,
# level, id)` of `src/core/interface.jl`. Every method here goes through the
# *chart kernel* — `ISEA4R`'s `cell_corners` / `cell_center` closed forms, which
# is also what the dense diamond grid (`face_grid.jl`) evaluates, reached by way
# of this module's base-9 `morton_to_xyd`. That sharing is the point: the
# DGGS-level cell and the face-grid cell for the same ordinal are LITERALLY the
# same four `Float64`-triples, not two evaluations that agree to a tolerance,
# and `test/ISEA9R/test_isea9r_kernel.jl` asserts it bitwise against
# `Trees.getcell(treeify(Isea9rFaceGrid(3^r; ordering = MortonOrder())), id + 1)`.
#
# It goes one step further than the ISEA4R sibling can: because the chart is
# *the same object*, the ISEA9R cell at `(level, id)` is bitwise the ISEA4R cell
# at the same lattice position whenever `3^level` is a resolution both admit.
# `test/ISEA9R/test_delegation.jl` pins that too.
#
# ## The id these methods take
#
# `id` is the **`isea9r_ordinal`** `diamond * 9^r + morton9_position` — 0-based,
# and exactly `MortonOrder`'s data position minus one. This *is* the system's
# canonical index (`canonical_index_name(ISEA9RDGGS()) === :isea9r_ordinal`),
# with `position` pinned to the base-9 Z-order Morton code (base-3 digit `k` of
# `ix` in the low slot of base-9 digit `k`, of `iy` in the high slot) in this
# package's own ten-diamond layout.
#
# **That layout has no external oracle behind it, and this ordinal is NOT
# DGGAL's `{LevelChar}{RootRhombus}-{HexIndex}`**: the two differ by a
# root→diamond permutation, a transpose of the in-rhombus square, and a 0.05°
# icosahedron orientation offset. No DGGRID, DGGAL or SphericalSpatialTrees
# identifier or geometry compatibility is claimed and none may be inferred —
# see `docs/design/isea9r_layout.md`, which records the three deltas and the
# fixture dump that would settle them. `level` is the refinement level `r`,
# `nside = 3^r`, `0 <= r <= 18` (the bound where `10 * nside^2` still fits
# `Int64` inside the borrowed chart's stated domain, matching
# `max_nside(Isea9rFaceSystem())`).
#
# ## What is NOT wired, and what it would cost
#
# The hierarchy / ordinal / pruning group — `cell_children`, `cell_parent`,
# `cell_descendants`, `cell_to_ordinal`, `ordinal_to_cell`, `descendant_range`,
# `num_cells`, `root_ids` — still throws `NotPortedError`, and
# `test/ISEA9R/test_isea9r_kernel.jl` pins that boundary. This is a scope line,
# not an obstacle, and it is deliberately the SAME line the ISEA4R sibling
# holds, so the two systems stay one decision rather than two.
#
# Recorded as now-possible: the radix-9 arithmetic is already exact over these
# ordinals. The base-9 Morton code drops its low base-9 digit per level up, and
# the within-diamond lattice nesting is bit-exact (`fl(ix/n) === fl(3ix/3n)`,
# see `ISEA4R.xyd_to_point` and `MortonOrder`), so a cell's children really are
# `9p:9p + 8`, its ancestor is `p ÷ 9^Δ`, and its `leaf_level` descendants
# really are the contiguous interval `[p * 9^Δ, (p + 1) * 9^Δ)` — which is
# exactly what `supports_prefix_ranges(ISEA9RDGGS()) == true` asserts at the
# interface level. `root_count = 10` and `radix = 9` are wired, so
# `has_ordinal_ids(::ISEA9RDGGS) = true` — one line — would derive the whole
# group from `src/core/kernel.jl` exactly as it does for HEALPix. What that line
# additionally *buys* is `DGGSGrid` / `DGGSPartialGrid` / `subtree_grid`
# support; what it *costs* is the kernel-test battery those grids deserve
# (ordinal round trips, the two-sided `descendant_range` contract, subtree-cap
# containment) and a decision on `has_exact_subtree_cap`. Deferred as a unit,
# deliberately, and jointly with ISEA4R's.
#
# ## Caps
#
# `cell_cap` is the exact four-corner cap rather than the generic 1.2-inflated
# center-to-vertex cap; the argument is at the method. `subtree_cap` is left
# alone: its generic fallback needs `cell_descendants` above level equality, so
# it answers exactly for `level == leaf_level` and throws `NotPortedError`
# otherwise — the honest state of a system whose hierarchy is not wired.
# ---------------------------------------------------------------------------

import ..DiscreteGlobalGrids as DGG
import ConservativeRegridding: Trees

# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------

# `level` is the refinement level `r`, so `nside = 3^r` and `morton_to_xyd`
# de-interleaves the ordinal into the `(ix, iy, diamond)` coordinates the chart
# is written against. `morton_to_xyd` is also where id validity is decided: it
# throws `ArgumentError` for an id outside `0:10 * 9^r - 1`. The de-interleave
# is this module's (`chart.jl`); `cell_corners` is `ISEA4R`'s, imported there —
# this is the single evaluation the face grid shares.
_cell_corners(level::Integer, id) =
    cell_corners(morton_to_xyd(id, 3^Int(level))..., 3^Int(level))

"""
    cell_boundary(::ISEA9RDGGS, level, id; closed=false)

The four cell corners in ring order, from the chart kernel
([`ISEA4R.cell_corners`](@ref)).

`id` is the `isea9r_ordinal` `diamond * 9^level + morton9_position` (0-based,
and `MortonOrder`'s data position minus one); `level` is the refinement level
`r`, `0 <= r <= 18`, so `nside = 3^r`. An id outside `0:10 * 9^r - 1` is an
`ArgumentError` from [`morton_to_xyd`](@ref). The ten-diamond layout these ids
are written against is this package's own convention with no external oracle
behind it, and the ordinal is *not* DGGAL's zone id — see
`docs/design/isea9r_layout.md` before reading any DGGRID / DGGAL / SST
compatibility into them.

The ring is **counterclockwise as seen from outside the sphere**, which is a
hard contract rather than a convention: the convex-clip kernel that computes
spherical intersections clips a clockwise ring to EMPTY, so a reversed ring
yields silent zero areas instead of an error. The order is the chart's lattice
order — `(x+, y+), (x-, y+), (x-, y-), (x+, y-)` — identical to what the diamond
grid emits for the same cell, corner for corner and bit for bit (same chart
functions, same seam-ownership rule), and identical to what `ISEA4RDGGS` emits
for the same lattice cell at the same `nside`.

Snyder ISEA maps the chart's straight lines to curves on the sphere, so — as
with HEALPix, and unlike S2 — the 4-vertex ring *approximates* the cell rather
than being it; densify through [`ISEA4R.xyd_to_point`](@ref) along the edges if
more is needed, and note that the caps below bound the curvature separately.
"""
function DGG.cell_boundary(::DGG.ISEA9RDGGS, level::Integer, id; closed::Bool=false)
    corners = _cell_corners(level, id)
    points = Vector{eltype(corners)}(undef, closed ? 5 : 4)
    @inbounds for i in 1:4
        points[i] = corners[i]
    end
    closed && (@inbounds points[5] = corners[1])
    return points
end

"""
    cell_center(::ISEA9RDGGS, level, id)

Native cell center: the chart evaluated at the lattice cell's midpoint
([`ISEA4R.cell_center`](@ref)). The rhombus chart is *exactly* equal-area
([`ISEA4R.xyd_to_point`](@ref)), so the chart midpoint is the canonical center —
the same argument HEALPix's `pixel_center` rests on, and the same one OGC
21-038r1 leans on when it calls ISEA9R zones equal-area. Preferred over the
kernel's normalized boundary mean: exact, one chart evaluation instead of four,
and identical to the diamond grid's center for the same cell.

(It is not the spherical centroid of the 4-gon; no equal-area DGGS claims that
of its cell centers.) `id` is the `isea9r_ordinal`; see [`cell_boundary`](@ref).
"""
DGG.cell_center(::DGG.ISEA9RDGGS, level::Integer, id) =
    cell_center(morton_to_xyd(id, 3^Int(level))..., 3^Int(level))

# Exact cell cap over the *chart* corners, replacing the generic
# `CELL_CAP_INFLATION` formula (max center-to-vertex distance × 1.2). The
# inflation exists for hierarchies where a parent does not geographically
# contain its children; an ISEA9R cell is the union of its nine children —
# their chart rectangles ninth its own, bit-exactly (`fl(ix/n) === fl(3ix/3n)`,
# see `MortonOrder`) — so there is no overhang to leave room for, and the
# four-corner cap is both sound and tighter.
#
# Soundness is the pre-registered measurement recorded at
# `cap_policy(::Isea9rFaceSystem)` in `face_grid.jl`, read on its 1x1 case: a
# dense 17x17 sampling of every block's chart rectangle — every leaf block among
# them, on every diamond, at `nside ∈ (3, 9, 27)` — showed a worst PRE-SLACK
# overhang of exactly `0.0`, i.e. the farthest sampled point of a cell is always
# one of its own four corners, seam-straddling cells included. The measurement
# ships as a standing test in `test/ISEA9R/test_face_grid.jl`;
# `test/ISEA9R/test_isea9r_kernel.jl` re-runs the containment check against
# these caps specifically. If the standing assert ever fires, this method must
# go with `cap_policy`'s override.
#
# Corner ORDER: `circle_from_four_corners` documents its argument as
# `(BL, TL, BR, TR)` and internally reorders to the CCW walk `(BL, BR, TR, TL)`
# before taking slerp midpoints of consecutive pairs. The chart ring is
# `(TR, TL, BL, BR)` in those terms, so it is permuted into the documented slots
# here — `(corners[3], corners[2], corners[4], corners[1])` — and consequently
# ALL FOUR sampled midpoints are true cell-edge midpoints, which is what the
# curvature of the Snyder edges actually needs bounded. The cap is then also
# bitwise the `FourCornerCap` node extent the face-grid layer would compute for
# the same cell as a 1x1 block. Same permutation, same reasoning, as
# `Isea4rKernel.jl`.
DGG.cell_cap(::DGG.ISEA9RDGGS, level::Integer, id) =
    (c = _cell_corners(level, id); Trees.circle_from_four_corners((c[3], c[2], c[4], c[1]), ()))

# --------------------------------------------------------------------------
# The user-facing connection
# --------------------------------------------------------------------------

"""
    cell_polygon(::ISEA9RDGGS, level, id) -> GI.Polygon

Cell `(level, id)` as a closed 4-gon on the **unit sphere**, i.e.
[`cell_polygon_unitsphere`](@ref) — the interface-level entry point, which
before this wiring threw `NotPortedError` for ISEA9R.

`id` is the `isea9r_ordinal` `diamond * 9^level + morton9_position` (see
[`cell_boundary`](@ref)), in this package's own ten-diamond layout, which no
external fixture pins. The ring is counterclockwise as seen from outside the
sphere and, Snyder edges not being great circles, approximates the cell to
4-gon accuracy.

```julia
cell_polygon(ISEA9RDGGS(), 2, 17)                  # this method
ISEA9R.cell_polygon(morton_to_xyd(17, 9)..., 9)    # the same polygon, chart-side
```

Only geometry is wired for ISEA9R: the hierarchy group (`cell_children`,
`cell_parent`, `descendant_range`, `num_cells`, ...) still throws
`NotPortedError` — deferred rather than blocked, since the radix-9 arithmetic
over these ordinals is exact. See the file header.
"""
DGG.cell_polygon(system::DGG.ISEA9RDGGS, level::Integer, id::Integer) =
    DGG.cell_polygon_unitsphere(system, level, id)
