# ---------------------------------------------------------------------------
# S2 cube-face charts — pure closed forms, no s2geometry
#
# S2 (Google's spherical geometry library) is six continuous charts
# `[0, 1]² → S²`, one per face of a cube inscribed in the sphere. Each chart is
# the composition of three steps: the quadratic `ST → UV` transform, the linear
# `UV → XYZ` lift onto the cube face, and radial normalisation onto the sphere.
# Unlike HEALPix's, these charts are *not* equal-area — S2 trades that away for
# geodesic cell edges: a line `u = const` on a cube face spans a plane through
# the origin, so its spherical image is a great-circle arc and a four-corner
# ring *is* the cell exactly, with no bulge between the corners.
#
# That is the whole reason a *chart* layer is worth having separately from the
# id layer — cell geometry, vertex geometry, and refinement all fall out of
# evaluating one function on the `(s, t)` lattice, with no Hilbert arithmetic in
# between, and lattice points shared between neighbouring cells come out
# bit-identical (⇒ the tessellation is exact, not merely consistent to
# rounding).
#
# Provenance: everything here is transcribed from the *published conventions*
# of s2geometry — the face frames and `S2_QUADRATIC_PROJECTION` of `s2coords.h`,
# and the `kIJtoPos` / `kPosToIJ` / `kPosToOrientation` Hilbert tables of
# `s2cell_id.h` — as documented at
#
#   https://s2geometry.io/devguide/s2cell_hierarchy.html
#   https://github.com/google/s2geometry
#   https://pkg.go.dev/github.com/golang/geo/s2
#
# No s2geometry code is vendored and no s2geometry library is called: the point
# is a self-contained kernel, and the tests here are correspondingly oracle-free
# (analytic invariants and internal consistency rather than a reference
# implementation — see `test/S2/test_chart.jl`, and contrast
# `test/HEALPix/test_chart.jl`, which does have Healpix.jl to check against).
#
# ## Index conventions
#
# Read this block before calling anything here.
#
# - `ix`, `iy` — 0-based face-local lattice coordinates in `0:nside-1`. `ix`
#   runs along the `s`/`u` axis and `iy` along the `t`/`v` axis; continuous
#   chart coordinates are `s = ix/nside`, `t = iy/nside`.
# - `face` — 0-based, `0:5`, in s2geometry's numbering: face `f` has its normal
#   on axis `f % 3`, positive for `f < 3` and negative for `f >= 3`.
# - Row-major id — 0-BASED: `face * nside² + iy * nside + ix`.
# - Hilbert id — 0-BASED: `face * nside² + hilbert_position`, which is exactly
#   the scaffold ordinal shape the `S2DGGS` registry entry records
#   (`face * 4^level + hilbert_position`, see `src/core/systems/s2.jl`).
#
# Both id spaces are 0-based and data positions are `id + 1` for BOTH orderings.
# That is a deliberate simplification relative to `HEALPix/chart.jl`, whose
# mixed 1-based-RING / 0-based-NESTED convention exists only because RING
# doubles as the position in a HEALPix FITS field. S2 has no such external
# layout to match, so there is no reason to carry two conventions.
#
# Argument order is `(ix, iy, face, nside)` throughout and ids come first in the
# inverse direction (`rowmajor_to_xyf(q, nside)`), matching `HEALPix/chart.jl`
# slot for slot so the two files diff against each other cleanly.
#
# Naming: the S2 literature says "cell", not "pixel", so the geometry helpers
# here are `cell_corners` / `cell_center` (and `cell_polygon` in `face_grid.jl`)
# — slot-for-slot siblings of HEALPix's `pixel_corners` / `pixel_center` /
# `pixel_polygon`. These are functions *in the `S2` namespace*, not methods of
# the top-level `cell_polygon(::AbstractDGGS, ...)` or of the kernel generic
# `cell_center`; that shadowing is established package precedent (see the export
# comment block in `src/DiscreteGlobalGrids.jl`).
# ---------------------------------------------------------------------------

import GeometryOps as GO

# ---------------------------------------------------------------------------
# Face frames
#
# The six faces of the inscribed cube, each carried by a right-handed
# orthonormal frame `(û, v̂, ŵ)` with `û × v̂ == ŵ` and `ŵ` the outward face
# normal. `FACE_NORMAL[f + 1]` is a signed coordinate axis, each of the six
# appearing exactly once; s2geometry's face numbering puts the normal on axis
# `f % 3`, positive for `f < 3`.
# ---------------------------------------------------------------------------

"""
    FACE_U_AXIS[face + 1] -> NTuple{3, Int}

Unit vector `û` of face `face` (0-based, `0:5`): the direction the chart's `u`
coordinate advances in. Together with [`FACE_V_AXIS`](@ref) and
[`FACE_NORMAL`](@ref) it forms a right-handed orthonormal frame, `û × v̂ == ŵ`.
"""
const FACE_U_AXIS = ((0, 1, 0), (-1, 0, 0), (-1, 0, 0), (0, 0, -1), (0, 0, -1), (0, 1, 0))

"""
    FACE_V_AXIS[face + 1] -> NTuple{3, Int}

Unit vector `v̂` of face `face` (0-based, `0:5`): the direction the chart's `v`
coordinate advances in. See [`FACE_U_AXIS`](@ref).
"""
const FACE_V_AXIS = ((0, 0, 1), (0, 0, 1), (0, -1, 0), (0, -1, 0), (1, 0, 0), (1, 0, 0))

"""
    FACE_NORMAL[face + 1] -> NTuple{3, Int}

Outward normal `ŵ` of face `face` (0-based, `0:5`) — the centre of the face,
and a signed coordinate axis. The six entries are the six signed axes, each
exactly once: `+x, +y, +z, -x, -y, -z`, which is s2geometry's face numbering.
"""
const FACE_NORMAL = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (-1, 0, 0), (0, -1, 0), (0, 0, -1))

"""
    face_uv_to_xyz(face, u, v) -> NTuple{3, Float64}

Lift face-local chart coordinates `(u, v) ∈ [-1, 1]²` onto the plane of cube
face `face` (0-based, `0:5`), returning an *unnormalised* point — the frame sum
`ŵ + u û + v v̂` of [`FACE_U_AXIS`](@ref) / [`FACE_V_AXIS`](@ref) /
[`FACE_NORMAL`](@ref).

Written as an explicit switch rather than as that sum for a reason that is
numerical, not stylistic: the switch places the literal `±1.0` on the fixed
axis and plain negations (`-u`, `-v`) on the others, so the coordinate two
faces share along a seam comes out *bit-identically* on both sides — IEEE
negation is exact, whereas `0 * x + 1 * y` style accumulation is not obliged to
be. The cross-face seam identities the chart tests pin (`==`, not `≈`) rest on
this.
"""
function face_uv_to_xyz(face::Integer, u::Float64, v::Float64)
    face == 0 && return (1.0, u, v)
    face == 1 && return (-u, 1.0, v)
    face == 2 && return (-u, -v, 1.0)
    face == 3 && return (-1.0, -v, -u)
    face == 4 && return (v, -1.0, -u)
    return (v, u, -1.0)
end

# ---------------------------------------------------------------------------
# ST <-> UV: the quadratic transform
#
# THE DECISION, recorded: this is s2geometry's `S2_QUADRATIC_PROJECTION`, and it
# is what makes the `(s, t)` lattice at `nside = 2^level` land on the *exact*
# canonical S2 cell boundaries. A face-grid cell is then the same spherical
# quadrilateral an `s2_cellid` cell is — the property the eventual native-id
# wiring needs. The scaffold-ordinal geometry wiring (`S2Kernel.jl`) already
# rests on it in the weaker sense that matters today: it and the face grid share
# `stf_to_point` and therefore clip literally the same polygons, exactly as the
# HEALPix face grid and the nested kernel share `xyf_to_point`.
#
# Rejected alternative: the linear map `u = 2s - 1` (`S2_LINEAR_PROJECTION`).
# Simpler and monotone too, but its cells are NOT canonical S2 cells, and its
# cell areas vary by a factor of ~5.2 within a level against the quadratic
# map's ~2.08. (The `S2_TAN_PROJECTION` is more uniform still but costs a
# `tan`/`atan` per coordinate and is likewise not the canonical convention.)
# ---------------------------------------------------------------------------

"""
    st_to_uv(s) -> Float64

s2geometry's `S2_QUADRATIC_PROJECTION`, mapping a chart coordinate
`s ∈ [0, 1]` to the cube-face coordinate `u ∈ [-1, 1]`:

```
u = (4s² - 1) / 3              for s >= 1/2
u = (1 - 4(1 - s)²) / 3        otherwise
```

Strictly increasing, odd about `s = 1/2` (`st_to_uv(1 - s) == -st_to_uv(s)`,
and *exactly* so in floating point — both branches consume the same `1 - s`),
and C¹ across the branch seam (both one-sided derivatives are `4/3`). The three
anchors are exact: `st_to_uv(0) == -1.0`, `st_to_uv(1/2) == 0.0`,
`st_to_uv(1) == 1.0`.

Choosing this transform rather than the linear `u = 2s - 1` is what puts the
`(s, t)` lattice at `nside = 2^level` on canonical S2 cell boundaries; see the
section comment above for the full decision record. It does *not* make the grid
equal-area — S2 is `is_equal_area == false` and stays that way — it only
narrows the within-level cell-area spread from ~5.2× to ~2.08×.

See [`uv_to_st`](@ref) for the inverse.
"""
st_to_uv(s::Real) = (sf = Float64(s); sf >= 0.5 ? (4 * sf^2 - 1) / 3 : (1 - 4 * (1 - sf)^2) / 3)

"""
    uv_to_st(u) -> Float64

Inverse of [`st_to_uv`](@ref): cube-face coordinate `u ∈ [-1, 1]` back to chart
coordinate `s ∈ [0, 1]`.

```
s = sqrt(1 + 3u) / 2           for u >= 0
s = 1 - sqrt(1 - 3u) / 2       otherwise
```

Shipped for the round-trip tests, and because point-to-cell lookup — still
unported — is `xyz → face/uv → st → lattice` and needs exactly this step.
"""
uv_to_st(u::Real) = (uf = Float64(u); uf >= 0 ? 0.5 * sqrt(1 + 3 * uf) : 1 - 0.5 * sqrt(1 - 3 * uf))

"""
    stf_to_point(s, t, face) -> GO.UnitSphericalPoint

Evaluate the S2 chart of `face` (0-based, `0:5`) at continuous face coordinates
`(s, t) ∈ [0, 1]²`, returning the corresponding point on the unit sphere:
quadratic [`st_to_uv`](@ref) per axis, then [`face_uv_to_xyz`](@ref), then
radial normalisation.

Defined for *any* real `s, t` in the unit square — nothing here is quantised to
a lattice or restricted to power-of-two `nside`, which is exactly what makes it
usable as the chart underlying refinement of arbitrary depth (and what lets the
face grid exist at `nside = 3`, `5`, `7`, which have no S2 id space at all).

The final step is the gnomonic projection of the cube face onto the sphere, and
it is why S2 cell edges are exact great-circle arcs: a chart line `u = const`
lifts to a straight line on the face plane, which together with the origin spans
a plane, whose intersection with the sphere is a great circle. Consequently the
four-corner ring of a cell is the cell *exactly*, with nothing to densify.

Unlike HEALPix's `xyf_to_point`, the normalisation here needs no accuracy
special-case: `n = sqrt(x² + y² + z²)` runs over `[1, √3]` with one component
always `±1`, so there is no cancellation to guard against (contrast the
`have_sintheta` branch, which exists because `sqrt((1-z)(1+z))` loses its
significant digits as `z → ±1`).
"""
function stf_to_point(s::Real, t::Real, face::Integer)
    x, y, z = face_uv_to_xyz(Int(face), st_to_uv(s), st_to_uv(t))
    n = sqrt(x * x + y * y + z * z)
    return GO.UnitSphericalPoint(x / n, y / n, z / n)
end

# ---------------------------------------------------------------------------
# Cell geometry
#
# Why the corner order below is counter-clockwise seen from outside, on all six
# faces, analytically (the tests then check it by exhaustion):
#
#  1. The emission order `(s+, t+), (s-, t+), (s-, t-), (s+, t-)` is CCW in the
#     `(s, t)` plane — its shoelace sum is `+2` on the unit square.
#  2. `st_to_uv` is strictly increasing, so the image is still CCW in `(u, v)`.
#  3. Every face frame `(û, v̂, ŵ)` is in SO(3), and the CCW-seen-from-outside
#     measure is SO(3)-invariant. So verifying one face verifies all six.
#  4. Face 0 at `nside = 1`: the unnormalised corners are `(1, 1, 1)`,
#     `(1, -1, 1)`, `(1, -1, -1)`, `(1, 1, -1)`, whose `Σ pᵢ × pᵢ₊₁` is
#     `(8, 0, 0)`; the outward direction `Σ pᵢ` is `(4, 0, 0)`, and the dot
#     product is `32 > 0`. CCW.
#  5. It extends to arbitrary sub-cells because the gnomonic map is an
#     orientation-preserving diffeomorphism from the face square onto the face's
#     spherical patch: any positively-oriented rectangle maps to a
#     positively-oriented spherical quadrilateral.
# ---------------------------------------------------------------------------

"""
    cell_corners(ix, iy, face, nside) -> NTuple{4, GO.UnitSphericalPoint}

The four corners of cell `(ix, iy)` on `face` at resolution `nside`, as lattice
points evaluated through [`stf_to_point`](@ref).

The corners are emitted **counter-clockwise as seen from outside the sphere**,
in the order `((ix+1)/n, (iy+1)/n)`, `(ix/n, (iy+1)/n)`, `(ix/n, iy/n)`,
`((ix+1)/n, iy/n)` — the same lattice order as HEALPix's `pixel_corners`. CCW
is a hard contract, not a stylistic preference: the convex-clip kernel used for
spherical intersection clips a clockwise ring to EMPTY, so a reversed ring
silently produces zero intersection area instead of an error. The winding holds
on every face for a structural reason (the emission order is CCW in `(s, t)`,
`st_to_uv` is increasing, and each face frame is in SO(3), under which the
measure is invariant); the file comment above spells the argument out.

Because every corner is a lattice point evaluated by the same function, two
cells sharing a lattice corner get *bit-identical* points and the tessellation
is exact. And unlike HEALPix — where the corners describe the pixel only to
4-gon accuracy, because pixel edges follow constant-`z` / constant-φ chart
lines — S2 cell edges **are** great circles ([`stf_to_point`](@ref)), so the
4-gon *is* the cell, exactly. There is nothing to densify.
"""
function cell_corners(ix::Integer, iy::Integer, face::Integer, nside::Integer)
    n = nside
    return (
        stf_to_point((ix + 1) / n, (iy + 1) / n, face),
        stf_to_point( ix      / n, (iy + 1) / n, face),
        stf_to_point( ix      / n,  iy      / n, face),
        stf_to_point((ix + 1) / n,  iy      / n, face),
    )
end

"""
    cell_center(ix, iy, face, nside) -> GO.UnitSphericalPoint

The center of cell `(ix, iy)` on `face` at resolution `nside`: the chart
evaluated at the lattice cell's midpoint `((ix + 0.5)/nside, (iy + 0.5)/nside)`.

This is S2's own definition of a cell center, not a convention borrowed from
HEALPix: `S2CellId::ToPoint` is precisely the ST-space midpoint pushed through
`ST → UV → XYZ` and normalised. (It is *not* the centroid of the spherical
quadrilateral, and S2 does not claim it is — the chart is not equal-area.)
"""
cell_center(ix::Integer, iy::Integer, face::Integer, nside::Integer) =
    stf_to_point((ix + 0.5) / nside, (iy + 0.5) / nside, face)

# ---------------------------------------------------------------------------
# Row-major order
#
# The plain lexicographic layout of the `nside × nside` lattice, `ix` fastest,
# then `iy`, then `face`. All closed-form arithmetic, valid for ANY
# `nside >= 1` — the power-of-two restriction belongs to the *Hilbert* index
# only, exactly as RING/NESTED split in `HEALPix/chart.jl`.
# ---------------------------------------------------------------------------

"""
    xyf_to_rowmajor(ix, iy, face, nside) -> Int64

0-based row-major id of cell `(ix, iy)` on `face` at resolution `nside`:
`face * nside² + iy * nside + ix`. Valid for any `nside >= 1`.

Inputs are assumed lattice-valid (`0 <= ix, iy < nside`, `0 <= face <= 5`) and
are not checked — garbage in, garbage out, matching `xyf_to_ring` and keeping
the cursor hot path branch-free. [`rowmajor_to_xyf`](@ref) *is* checked, being
the direction a user-supplied id enters through.
"""
xyf_to_rowmajor(ix::Integer, iy::Integer, face::Integer, nside::Integer) =
    Int64(face) * Int64(nside)^2 + Int64(iy) * Int64(nside) + Int64(ix)

"""
    rowmajor_to_xyf(q, nside) -> (ix, iy, face)

Inverse of [`xyf_to_rowmajor`](@ref): face-local lattice coordinates and 0-based
face of the 0-based row-major id `q`. Valid for any `nside >= 1`; throws an
`ArgumentError` for an id outside `0:6nside²-1`.
"""
function rowmajor_to_xyf(q::Integer, nside::Integer)
    npface = Int64(nside)^2
    qi = Int64(q)
    0 <= qi < 6npface || throw(ArgumentError(
        "row-major id $qi out of range for nside=$nside (expected 0:$(6npface - 1))"))
    face, r = divrem(qi, npface)
    iy, ix = divrem(r, Int64(nside))
    return (Int(ix), Int(iy), Int(face))
end

# ---------------------------------------------------------------------------
# Hilbert order
#
# The Hilbert curve is to S2 what the Morton curve is to HEALPix: the canonical
# within-face ordering, and the one the registry's scaffold ordinal is written
# against (`face * 4^level + hilbert_position`, see `src/core/systems/s2.jl`).
# So [`HilbertOrder`](@ref) is to S2 exactly what `NestedOrder` is to HEALPix,
# and — like the Morton code — it needs a `2^k × 2^k` face and says so loudly.
#
# The tables below are s2geometry's `kIJtoPos` / `kPosToIJ` /
# `kPosToOrientation`, indexed by a two-bit orientation state built from
# `SWAP_MASK` and `INVERT_MASK`. At each level the quadrant `ij = 2i + j`
# (`i` = the `ix`/s-axis bit, high; `j` = the `iy`/t-axis bit, low) is looked up
# to a position along the curve, the position contributes two bits, and the
# orientation is advanced by XOR. Odd faces start swapped, which is what chains
# the six per-face curves into a single closed loop over the whole sphere.
#
# A plain Morton ordering is deliberately NOT shipped alongside: S2 has no
# Morton convention to be compatible with, `RowMajorOrder` already covers every
# `nside`, and the [`AbstractS2Ordering`](@ref) contract makes a Morton ordering
# a ~15-line user addition if anyone wants one.
# ---------------------------------------------------------------------------

"""
    SWAP_MASK

Orientation bit meaning "the `i` and `j` axes are exchanged" in the Hilbert
state machine (s2geometry's `kSwapMask`). See [`IJ_TO_POS`](@ref).
"""
const SWAP_MASK = 1

"""
    INVERT_MASK

Orientation bit meaning "both axes are reflected" in the Hilbert state machine
(s2geometry's `kInvertMask`). See [`IJ_TO_POS`](@ref).
"""
const INVERT_MASK = 2

"""
    IJ_TO_POS[orientation + 1][ij + 1] -> Int

Quadrant `ij = 2i + j` (`i` the `ix`/s-axis bit, `j` the `iy`/t-axis bit) to its
0-based position along the Hilbert curve, under the current `orientation` in
`0:3` (`SWAP_MASK | INVERT_MASK`). s2geometry's `kIJtoPos`.
"""
const IJ_TO_POS = ((0, 1, 3, 2), (0, 3, 1, 2), (2, 3, 1, 0), (2, 1, 3, 0))

"""
    POS_TO_IJ[orientation + 1][pos + 1] -> Int

Inverse of [`IJ_TO_POS`](@ref): position along the curve back to the quadrant
`ij = 2i + j`. s2geometry's `kPosToIJ`.
"""
const POS_TO_IJ = ((0, 1, 3, 2), (0, 2, 3, 1), (3, 2, 0, 1), (3, 1, 0, 2))

"""
    POS_TO_ORIENTATION[pos + 1] -> Int

Orientation delta XORed into the state after descending into the sub-cell at
0-based curve position `pos`. s2geometry's `kPosToOrientation`.
"""
const POS_TO_ORIENTATION = (SWAP_MASK, 0, 0, SWAP_MASK | INVERT_MASK)

"""
    xyf_to_hilbert(ix, iy, face, nside) -> Int64

0-based Hilbert id of cell `(ix, iy)` on `face` at resolution `nside`, which
must be a power of two — throws `ArgumentError` otherwise.

The id is `face * nside² + hilbert_position`, i.e. the scaffold ordinal shape
the `S2DGGS` registry entry records (`face * 4^level + hilbert_position`).
Within a face the position is built two bits per level, most significant first,
by walking [`IJ_TO_POS`](@ref) and advancing the orientation state through
[`POS_TO_ORIENTATION`](@ref); odd faces start at `SWAP_MASK`.

Dropping the low two bits of the position steps exactly one level up, so the
Hilbert order nests across resolutions the same way the Morton order does —
which is what makes cross-resolution refinement observable as contiguous blocks
(the face-grid tests check this directly, `4j-3:4j`).

!!! note "Alignment with native `s2_cellid` is intended, not yet verified"
    The tables and the odd-face initial swap transcribe s2geometry's `kIJtoPos`
    / `kPosToIJ` / `kPosToOrientation` and `S2CellId::FromFaceIJ`, so agreement
    with native S2 cell ids is intended *by construction*. It is **not**
    oracle-verified: this repository carries no s2geometry fixtures, and the
    tests here pin bijectivity, Hilbert locality, prefix nesting across levels
    and hand-computed level-1/2 tables instead. Claims about native ids wait for
    the id-hierarchy milestone.
"""
function xyf_to_hilbert(ix::Integer, iy::Integer, face::Integer, nside::Integer)
    ispow2(nside) || throw(ArgumentError(
        "the Hilbert index is only defined for nside = 2^k, got nside=$nside; \
         use `xyf_to_rowmajor` for arbitrary nside"))
    (0 <= ix < nside && 0 <= iy < nside) || throw(ArgumentError(
        "lattice coordinates ($ix, $iy) out of range for nside=$nside (expected 0:$(nside - 1))"))
    0 <= face <= 5 || throw(ArgumentError("face $face out of range (expected 0:5)"))
    k = trailing_zeros(Int64(nside))                # nside == 2^k
    # Odd faces start swapped; that is what joins the six per-face curves into
    # one continuous loop over the sphere.
    orientation = isodd(face) ? SWAP_MASK : 0
    pos = Int64(0)
    for b in (k - 1):-1:0
        ij = 2 * ((Int(ix) >> b) & 1) + ((Int(iy) >> b) & 1)
        p = IJ_TO_POS[orientation + 1][ij + 1]
        pos = (pos << 2) | p
        orientation ⊻= POS_TO_ORIENTATION[p + 1]
    end
    return Int64(face) * Int64(nside)^2 + pos
end

"""
    hilbert_to_xyf(h, nside) -> (ix, iy, face)

Inverse of [`xyf_to_hilbert`](@ref): face-local lattice coordinates and 0-based
face of the 0-based Hilbert id `h`. `nside` must be a power of two — throws
`ArgumentError` otherwise, and for an id outside `0:6nside²-1`.

The mirror image of the forward walk: split off the face, then consume the
position two bits at a time through [`POS_TO_IJ`](@ref), rebuilding `ix` and
`iy` one bit per level while advancing the same orientation state.

!!! note "Alignment with native `s2_cellid` is intended, not yet verified"
    See [`xyf_to_hilbert`](@ref) — the same caveat applies verbatim.
"""
function hilbert_to_xyf(h::Integer, nside::Integer)
    ispow2(nside) || throw(ArgumentError(
        "the Hilbert index is only defined for nside = 2^k, got nside=$nside; \
         use `rowmajor_to_xyf` for arbitrary nside"))
    npface = Int64(nside)^2
    hid = Int64(h)
    0 <= hid < 6npface || throw(ArgumentError(
        "Hilbert id $hid out of range for nside=$nside (expected 0:$(6npface - 1))"))
    face, pos = divrem(hid, npface)
    k = trailing_zeros(Int64(nside))
    orientation = isodd(face) ? SWAP_MASK : 0
    ix = 0
    iy = 0
    for b in (k - 1):-1:0
        p = Int((pos >> (2b)) & 3)
        ij = POS_TO_IJ[orientation + 1][p + 1]
        ix = (ix << 1) | ((ij >> 1) & 1)
        iy = (iy << 1) | (ij & 1)
        orientation ⊻= POS_TO_ORIENTATION[p + 1]
    end
    return (ix, iy, Int(face))
end
