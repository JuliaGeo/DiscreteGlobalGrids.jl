# ---------------------------------------------------------------------------
# HEALPix operations-kernel wiring
#
# Nested HEALPix is the reference *ordinal-id* system of `core/kernel.jl`: the
# canonical ids at `level` are exactly `0:12 * 4^level - 1` and pixel `p`'s
# children are `4p:4p + 3`, so `has_ordinal_ids` alone derives the whole
# hierarchy / ordinal / pruning group from `root_count = 12` and `radix = 4`.
# Only geometry is wired here, and it goes through the *chart kernel*
# (`chart.jl`) — the same `pixel_corners` / `pixel_center` closed forms the
# dense face grid (`face_grid.jl`) evaluates, reached by way of
# `nested_to_xyf`. That sharing is the point: an id-hierarchy cell and the
# face-grid cell for the same pixel are now LITERALLY the same four
# `Float64`-triples, so a `Regridder` built on `DGGSGrid(HEALPixDGGS(), level)`
# and one built on `HealpixFaceGrid(2^level; ordering = NestedOrder())` agree
# bitwise rather than to a tolerance.
#
# (Until milestone 3 this file called `Healpix.boundariesRing` / `pix2vecNest`
# instead, which was byte-identical to the pre-migration per-system HEALPix
# tree. The chart forms agree with those to 1.1e-16 per coordinate on corners
# and 8.9e-16 on centers — including the corner ORDER, which is the same ring
# `boundariesRing` emits — and `test/HEALPix/test_healpix_kernel.jl` keeps that
# cross-check running against Healpix.jl for as long as the dependency is here.)
#
# Caps are exact rather than inflated: a HEALPix parent pixel *is* the
# geographic union of its children, so the pixel's own 4-corner cap bounds its
# entire subtree and `subtree_cap` overrides the generic union/inflation logic
# with it (`test/HEALPix/test_healpix_kernel.jl` verifies the containment).
#
# `subtree_border` is overridden for the *same* structural fact, stated on the
# lattice rather than on the sphere: a subtree here is a square block of the
# leaf lattice, so its rim is that block's perimeter and can be emitted
# directly in ascending id order, without a single neighbor query. The
# derivation sits above the wiring below.
# ---------------------------------------------------------------------------

import ..DiscreteGlobalGrids as DGG
import ConservativeRegridding: Trees
import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeoInterface as GI
# For `cell_neighbors` only: the neighbor tables live in `HealpixLookups`
# (`nested_neighbors`, a healpix_base.cc transcription over Healpix.jl's xyf
# codec), so this one operation rides Healpix.jl rather than the chart.
import Healpix
using SmallCollections: SmallVector

# --------------------------------------------------------------------------
# Id model
#
# Nothing else from the hierarchy/ordinal/pruning group is wired: the kernel
# defaults are the nested-HEALPix arithmetic verbatim (`cell_children` =
# `4p:4p+3`, `descendant_range` = `[p * 4^Δ, (p + 1) * 4^Δ)`, `cell_to_ordinal`
# = `p + 1`, `num_cells` = `12 * 4^level`).
# --------------------------------------------------------------------------

DGG.has_ordinal_ids(::DGG.HEALPixDGGS) = true

# HEALPix parents geographically contain their children and the 4-corner cap is
# O(1), so the exact parent cap beats a stored-id union cap at every internal
# node — which is what the old per-system HEALPix node extents did.
DGG.has_exact_subtree_cap(::DGG.HEALPixDGGS) = true

# The same containment, stated as the trait: nested HEALPix is a congruent
# refinement. This is the fact the `subtree_cap` override above rests on, and
# what makes the densified pixel outline wired below (`subtree_polygon_
# unitsphere`) an exact subtree bound for any traversal that wants polygon-
# level subtree classification.
DGG.has_congruent_geometry(::DGG.HEALPixDGGS) = true

# --------------------------------------------------------------------------
# Neighbors
#
# The 3×3 lattice-block neighborhood (edge + corner neighbors), the HEALPix
# stencil convention: 8 per pixel, 7 at the 24 pixels sitting on a degree-3
# base-tiling vertex. `nested_neighbors` answers in compass order with -1 for
# the missing diagonal; the kernel contract is ascending ids, existing
# neighbors only.
# --------------------------------------------------------------------------

DGG.max_neighbors(::DGG.HEALPixDGGS) = 8

function DGG.cell_neighbors(system::DGG.HEALPixDGGS, level::Integer, id)
    # The ordinal id guard of the kernel's hierarchy defaults: pix2xyfNest
    # happily un-Mortons an id no pixel has, into neighbors of a cell that
    # does not exist.
    total = DGG.num_cells(system, level)
    0 <= id < total || throw(ArgumentError(
        "$(DGG.system_name(system)) level-$(Int(level)) cell id $id is out of range 0:$(total - 1)"))
    compass = HealpixLookups.nested_neighbors(Healpix.Resolution(2^Int(level)), Int(id))
    out = SmallVector{8,Int64}()
    for neighbor in compass
        neighbor < 0 && continue
        out = DGG._insert_sorted(out, Int64(neighbor))
    end
    return out
end

# --------------------------------------------------------------------------
# Subtree border
#
# The subtree rim without a single neighbor query — the same operation IGEO7
# reads off its Z7 digits (`IGeo7/grid.jl`), here read off the leaf lattice.
#
# WHY THE SUBTREE IS A SQUARE BLOCK. A nested id is
# `face * nside^2 + morton(ix, iy)` (`xyf_to_nested`, chart.jl), and refining
# Δ levels scales the lattice by `s = 2^Δ` *on the same face*: pixel
# `(ix, iy, face)` at `level` covers exactly `[ix*s, (ix+1)*s) x
# [iy*s, (iy+1)*s)` of the `leaf_level` lattice. Bit-interleaving is
# positional, so that scaling is a shift of the Morton code:
#
#     morton(ix*s + dx, iy*s + dy) = morton(ix, iy) * s^2 + morton(dx, dy)
#
# — the high `2*level` bits stay the parent's code, the low `2Δ` bits are
# free. Two consequences, and the whole method is built on them:
#
#   * the subtree is the contiguous id range `[id * 4^Δ, (id+1) * 4^Δ)`, which
#     is exactly what `descendant_range` already returns; and
#   * a descendant's OFFSET within that range *is* `morton(dx, dy)`. So
#     ascending id order over the subtree is Morton order over the block, and
#     the operation reduces to "emit the perimeter of an `s x s` square in
#     Morton order".
#
# WHY THE RIM IS THAT PERIMETER. `cell_neighbors` here is the 3x3 lattice
# block, so a descendant at offset `(dx, dy)` has every neighbor inside the
# subtree iff `dx ± 1` and `dy ± 1` all stay in `0:s-1` — iff `0 < dx < s-1`
# and `0 < dy < s-1`. The rim is therefore `dx ∈ {0, s-1} || dy ∈ {0, s-1}`,
# of size `4s - 4` against the `s^2` cells the subtree holds.
#
# Three things could break that argument. None does:
#
#   * *Corner vs edge adjacency.* The 8-neighborhood gives the IDENTICAL rim
#     to a 4-neighborhood: every perimeter cell already has an EDGE neighbor
#     outside (step outward across the side it sits on), and every interior
#     cell has all 8 inside. The thicker HEALPix neighborhood costs nothing
#     here — the rim is one cell wide either way — so this override is not
#     quietly answering a different question than the aperture-7 ones do.
#   * *Face seams.* A lattice step off a face edge wraps through
#     `NB_FACEARRAY` (`HealpixLookups.nested_neighbors`), which could in
#     principle deposit the neighbor back inside the block. It cannot,
#     structurally: no entry of any NON-CENTER row of that table maps a face
#     to itself, so a wrapped neighbor always lands on a different base face —
#     and a subtree block never leaves its own face. A seam-touching block is
#     rim-exposed for precisely the same reason an interior one is, with no
#     special case.
#   * *The 24 degenerate pixels.* The pixels with only 7 neighbors sit on the
#     8 degree-3 vertices of the base tiling, which in lattice terms are face
#     CORNERS (`ix, iy ∈ {0, nside-1}`). A face corner is a corner of whatever
#     block contains it, and a block corner has 5 of its 8 neighbors outside —
#     only its own inward 2x2 quadrant is inside — so losing the one missing
#     diagonal still leaves 4. The margin, not luck, is why the missing
#     neighbor cannot un-rim it.
#
# WHY THE WALK RECURSES. Walking the perimeter geometrically (along the
# bottom, up the right, back along the top) is NOT monotone in Morton order,
# so it would need a sort. `_rim_walk!` instead recurses in Morton-quadrant
# order: the four quadrants of an `s x s` square occupy consecutive id blocks
# at offsets `q * (s/2)^2` for `q = qx + 2*qy`, so visiting `q = 0, 1, 2, 3`
# and pruning the quadrants that inherit none of the parent square's exposed
# sides emits in strictly ascending id order by construction. The pruning is
# what keeps it `O(rim)`: a quadrant touching no exposed side is interior, and
# its whole `(s/2)^2` subtree is skipped rather than filtered.
#
# COST, and what is actually won. Θ(rim) = Θ(2^Δ) time, ONE allocation, O(Δ)
# stack. The generic fallback is *also* Θ(2^Δ) in cells — its level-by-level
# expansion prunes to the rim as it goes — so the win here is the constant, not
# the exponent: the fallback runs a `cell_neighbors` sweep and a `cell_parent`
# ascent for every child of every rim cell at every intermediate level (about
# `Σ_d 4 * (4 * 2^d - 4)` children examined, each with up to 8 neighbor
# lookups), where this walk does one integer addition per cell it emits and
# never touches the neighbor tables at all. Measured 50x / 123x / 149x faster
# at Δ = 4 / 6 / 8, and the gap widens with Δ because the fallback's
# intermediate levels accumulate.
# --------------------------------------------------------------------------

# Which sides of a square are exposed to the outside of the subtree, one bit
# per side of the face-local lattice. A quadrant inherits only the sides it
# actually shares with its parent square — the left half can be exposed at
# `dx = 0` but never at `dx = s-1` — which is the pruning rule in one line.
const _RIM_XMIN = 0x1
const _RIM_XMAX = 0x2
const _RIM_YMIN = 0x4
const _RIM_YMAX = 0x8

"""
    _rim_walk!(out, k, base, code, sz, mask) -> k

Append the perimeter of the `sz x sz` sub-square whose lower-left corner has
within-subtree Morton code `code`, restricted to the sides named by `mask`
(`_RIM_XMIN` etc.), to `out` starting at `out[k + 1]`; return the new fill
mark. Cells are emitted as `base + code'`, `base` being the subtree's first
id — `descendant_range`'s `lo`.

Ascending in the id, because Morton quadrants are visited in id order; `O(cells
emitted)`, because a quadrant inheriting no exposed side is skipped whole
rather than descended and filtered. Never called with `mask == 0`, which is
exactly that skipped case.
"""
function _rim_walk!(out::Vector{Int64}, k::Int, base::Int64, code::Int64,
        sz::Int64, mask::UInt8)
    if sz == 1
        @inbounds out[k += 1] = base + code
        return k
    end
    half = sz >> 1
    quarter = half * half
    for q in Int64(0):Int64(3)
        # `q`'s low bit picks the x half, its high bit the y half (the Morton
        # convention `nested_to_xyf` inverts: even bits -> ix, odd -> iy), so a
        # quadrant keeps the parent's XMIN side only if it is the low-x half,
        # and so on. `m == 0` means every side of this quadrant is interior.
        m = ((q & 1) == 0 ? (mask & _RIM_XMIN) : (mask & _RIM_XMAX)) |
            ((q >> 1) == 0 ? (mask & _RIM_YMIN) : (mask & _RIM_YMAX))
        m == 0x0 && continue
        k = _rim_walk!(out, k, base, code + q * quarter, half, m)
    end
    return k
end

# `descendant_range` runs FIRST, and it is the only guard this method needs: it
# raises the kernel's `level <= leaf_level` `ArgumentError` (the same one the
# generic fallback raises, for the reason `cell_descendants` gives), then the
# ordinal id range check, and only then hands back `lo`. Nothing below
# re-derives either.
#
# In particular the depth-0 answer is `[lo]` — the validated range's own single
# id — never an unchecked `[id]`. That is the same routing the generic fallback
# takes through `cell_descendants` (which for HEALPix *is* `descendant_range`),
# and it is what makes this branch reject a pixel that does not exist instead of
# confidently returning it.
function DGG.subtree_border(system::DGG.HEALPixDGGS, level::Integer, id, leaf_level::Integer)
    lo, _ = DGG.descendant_range(system, level, id, leaf_level)
    delta = Int(leaf_level) - Int(level)
    delta == 0 && return Int64[lo]           # `cell_id_type(HEALPixDGGS()) === Int64`
    s = Int64(1) << delta
    out = Vector{Int64}(undef, 4 * s - 4)
    k = _rim_walk!(out, 0, Int64(lo), Int64(0), s,
        _RIM_XMIN | _RIM_XMAX | _RIM_YMIN | _RIM_YMAX)
    # Deliberately NOT an `@assert`. Assertions are elided under
    # `--check-bounds=no`-style flags, and `@assert` in this package is for
    # load-time table construction (`ISEA4R/diamonds.jl`), never for a returned
    # result. This is the one invariant here whose violation would be *silent*:
    # a short walk leaves `undef` slots, i.e. arbitrary `Int64`s handed back as
    # cell ids and then indexed with. It costs one comparison per call against
    # Θ(2^Δ) work, so it stays on unconditionally.
    #
    # (The alternative the IGEO7 side takes — `sizehint!` plus `push!`, making a
    # wrong census cost a reallocation rather than a result — is not needed
    # here: `4s - 4` is exact for every Δ, so the exact allocation is free and
    # this check is the cheap way to keep it honest.)
    k == length(out) || error(
        "HEALPix subtree_border filled $k of $(length(out)) rim slots for \
         level-$(Int(level)) pixel $id at leaf level $(Int(leaf_level))")
    return out
end

# --------------------------------------------------------------------------
# Geometry
# --------------------------------------------------------------------------

# Canonical ids are 0-based NESTED (EOPF convention) and `level` is a nested
# refinement level, so `nside = 2^level` and `nested_to_xyf` un-Mortons the id
# into the `(ix, iy, face)` lattice coordinates the chart is written against.
# Both calls are in `chart.jl`, which is Healpix.jl-free by construction — this
# is the single evaluation the face grid shares.
_pixel_corners(level::Integer, pixel::Integer) =
    pixel_corners(nested_to_xyf(pixel, 2^Int(level))..., 2^Int(level))

"""
    cell_boundary(::HEALPixDGGS, level, id; closed=false)

The four pixel corners in ring order, from the chart kernel
([`pixel_corners`](@ref)).

The ring is **counter-clockwise as seen from outside the sphere**, which is a
hard contract rather than a convention: the convex-clip kernel that computes
spherical intersections clips a clockwise ring to EMPTY, so a reversed ring
yields silent zero areas instead of an error. The order is the one
`Healpix.boundariesRing` emits (`(north, west, south, east)`), so order-
sensitive consumers written against the pre-chart implementation see the same
ring.

HEALPix cell edges are curved (constant-latitude / constant-longitude chart
arcs, not great circles), so the 4-vertex ring is the same spherical-polygon
approximation the old per-system HEALPix tree handed to
ConservativeRegridding; the caps below bound the curvature separately.
"""
function DGG.cell_boundary(::DGG.HEALPixDGGS, level::Integer, id; closed::Bool=false)
    corners = _pixel_corners(level, id)
    points = Vector{eltype(corners)}(undef, closed ? 5 : 4)
    @inbounds for i in 1:4
        points[i] = corners[i]
    end
    closed && (@inbounds points[5] = corners[1])
    return points
end

"""
    cell_center(::HEALPixDGGS, level, id)

Native pixel center: the chart evaluated at the lattice cell's midpoint
([`pixel_center`](@ref)), which is the HEALPix center by definition because the
chart is equal-area. Preferred over the kernel's normalized boundary mean —
exact, one chart evaluation instead of four, and identical to the face grid's
center for the same pixel. (Agrees with `Healpix.pix2vecNest` to 8.9e-16 per
coordinate.)
"""
DGG.cell_center(::DGG.HEALPixDGGS, level::Integer, id) =
    pixel_center(nested_to_xyf(id, 2^Int(level))..., 2^Int(level))

# Exact pixel cap over the *chart* corners: `circle_from_four_corners` centers
# on the corner mean and takes the max distance to the corners and their
# great-circle edge midpoints (plus 1.0001 slack), which covers the curved
# pixel boundary — no cell-level inflation needed, so this replaces the generic
# `CELL_CAP_INFLATION` formula. Same inputs and same formula as the face grid's
# `node_extent` on a 1x1 block, so leaf extents on both paths are algebraically
# identical (not asserted bitwise).
#
# Corner ORDER, and why it is sound to pass the ring straight through:
# `pixel_corners` yields the ring `(N, W, S, E)`, whereas
# `circle_from_four_corners` is written for `(BL, TL, BR, TR)` and internally
# reorders to walk its own perimeter. Under that mismatch only two of the four
# midpoints it samples — W–N and S–E — are true pixel-edge midpoints; the other
# two land on the pixel's diagonals instead. The cap is nevertheless a valid
# bound: a HEALPix pixel is mirror-symmetric across its W–E diagonal, so each
# unsampled edge is the reflection of a sampled one and has the same distance
# from the (symmetric) corner-mean center, and the 1.0001 slack absorbs the
# ulp-level asymmetry. Do NOT "fix" the order: the caps are bitwise-pinned by
# the cap tests, and reordering would perturb every recorded value.
DGG.cell_cap(::DGG.HEALPixDGGS, level::Integer, id) =
    Trees.circle_from_four_corners(_pixel_corners(level, id), ())

# A parent pixel contains its descendants geographically, so its own exact cap
# is the subtree cap — O(1) instead of the generic union over up to
# `SUBTREE_CAP_EXACT_LIMIT` leaves, and tighter than the inflated fallback.
DGG.subtree_cap(system::DGG.HEALPixDGGS, level::Integer, id, leaf_level::Integer) =
    DGG.cell_cap(system, level, id)

# Where the subtree outline below stops being worth building: the ring costs
# `4 * 2^Δ` chart evaluations and any predicate over it scales with as many
# segments. 8 is the old planar descent's `MAX_QUERY_DENSIFY_DELTA` cutoff:
# a 1,025-vertex ring bounding a 65,536-leaf subtree.
const MAX_SUBTREE_POLYGON_DELTA = 8

"""
    subtree_polygon_unitsphere(::HEALPixDGGS, level, id, leaf_level)

The exact spherical outline of pixel `(level, id)`'s `leaf_level` subtree:
the pixel's perimeter walked in **leaf-lattice** steps, every vertex the
chart evaluated at a leaf corner (`xyf_to_point`, the same closed form
`pixel_corners` evaluates). Descendant corners on the parent's perimeter are
those very lattice points, so the ring's great-circle edges are *identical* —
bit-for-bit endpoints — to the boundary edges of the perimeter leaf polygons,
and the ring bounds the union of descendant [`cell_polygon_unitsphere`](@ref)
4-gons exactly. (The pixel's own 4 corners would not: HEALPix edges are
chart arcs, not great circles, so leaf corners along an edge bulge off the
corner-to-corner chord — the spherical restatement of the `step = 2^Δ`
densification the old planar query descent needed.)

The ring is counter-clockwise seen from outside the sphere (the
`pixel_corners` contract) and has `4 * 2^Δ` vertices, so past
`MAX_SUBTREE_POLYGON_DELTA` levels of separation the method answers
`nothing` rather than build an outline that big — a traversal descends one
level and asks again, with Δ back under the cap.
"""
function DGG.subtree_polygon_unitsphere(system::DGG.HEALPixDGGS, level::Integer, id,
        leaf_level::Integer)
    level <= leaf_level || throw(ArgumentError("expected level <= leaf_level"))
    delta = Int(leaf_level) - Int(level)
    delta <= MAX_SUBTREE_POLYGON_DELTA || return nothing
    # The ordinal id guard of the kernel's hierarchy defaults: `nested_to_xyf`
    # happily un-Mortons an id no pixel has.
    total = DGG.num_cells(system, level)
    0 <= id < total || throw(ArgumentError(
        "$(DGG.system_name(system)) level-$(Int(level)) cell id $id is out of range 0:$(total - 1)"))
    s = 1 << delta                       # leaf cells per pixel side
    n = 1 << Int(leaf_level)             # nside of the leaf lattice
    ix, iy, face = nested_to_xyf(id, 1 << Int(level))
    x0 = ix * s
    y0 = iy * s
    points = Vector{GO.UnitSphericalPoint{Float64}}(undef, 4 * s + 1)
    k = 0
    # Perimeter of the pixel's leaf-lattice square, starting at the "north"
    # corner `((ix+1)/nside, (iy+1)/nside)` and running in `pixel_corners`'
    # N → W → S → E order (CCW from outside the sphere).
    for i in 0:(s - 1)
        points[k += 1] = xyf_to_point((x0 + s - i) / n, (y0 + s) / n, face)
    end
    for i in 0:(s - 1)
        points[k += 1] = xyf_to_point(x0 / n, (y0 + s - i) / n, face)
    end
    for i in 0:(s - 1)
        points[k += 1] = xyf_to_point((x0 + i) / n, y0 / n, face)
    end
    for i in 0:(s - 1)
        points[k += 1] = xyf_to_point((x0 + s) / n, (y0 + i) / n, face)
    end
    points[k += 1] = points[1]
    return GI.Polygon([GI.LinearRing(points)])
end

"""
    cell_polygon(::HEALPixDGGS, level, id) -> GI.Polygon

Pixel `(level, id)` as a closed 4-gon on the **unit sphere**, i.e.
[`cell_polygon_unitsphere`](@ref) — the interface-level entry point of
`src/core/interface.jl`, whose fallback throws `NotPortedError`.

`id` is the 0-based nested (EOPF) pixel id. The ring is counter-clockwise as
seen from outside the sphere, and — HEALPix pixel edges being constant-latitude
/ constant-longitude chart arcs rather than geodesics — approximates the pixel
to 4-gon accuracy; see [`cell_boundary`](@ref).

Wired here rather than left to the fallback purely for uniformity across the
kernel-wired systems: S2 and ISEA4R answer this generic over their own ordinals
(`S2Kernel.jl`, `Isea4rKernel.jl`), and a system whose geometry group is fully
wired should not be the one that still throws at the interface level.
"""
DGG.cell_polygon(system::DGG.HEALPixDGGS, level::Integer, id::Integer) =
    DGG.cell_polygon_unitsphere(system, level, id)

# --------------------------------------------------------------------------
# Lookup convenience
# --------------------------------------------------------------------------

"""
    DGGSPartialGrid(l::HealpixLookup; kwargs...)

The lookup's stored cells as a generic partial grid. Cell-id validity at
`l.level` is trusted from the lookup (which range-checks ids against
`12 * 4^level` on construction); the generic constructor only re-checks
ordering and element type. `kwargs` reach `DGGSPartialGrid`'s `bucket_size` /
`root_level` / `root_id`.
"""
DGG.DGGSPartialGrid(l::HealpixLookups.HealpixLookup; kwargs...) =
    DGG.DGGSPartialGrid(DGG.HEALPixDGGS(), l.level, l.data; kwargs...)

# What the generic lookup operations (`neighbor_indices`, `stencil`, `zonal`)
# ask of a lookup: which system, which level.
DGG.dggs_system(::HealpixLookups.HealpixLookup) = DGG.HEALPixDGGS()
DGG.dggs_level(l::HealpixLookups.HealpixLookup) = l.level

# Treeifying a lookup directly is the shortest path from a `DimensionalData`
# dimension to a `Regridder`, and it needs nothing from this file: the method
# is generic over `AbstractDGGSLookup` (`core/lookups.jl`), routing a stored id
# vector through the constructor above and a `DGGSGlobeIds` through the dense
# `DGGSGrid` instead. All that is per-system is the manifold, which is what
# makes the one-argument `treeify(l)` resolve at all.
GOCore.best_manifold(::HealpixLookups.HealpixLookup) = GO.Spherical()
