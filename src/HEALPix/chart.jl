# ---------------------------------------------------------------------------
# HEALPix face charts — pure closed forms, no Healpix.jl
#
# HEALPix (Górski et al. 2005, DOI:10.1086/427976) is twelve continuous charts
# `[0, 1]² → S²`, one per base face. Every chart is an equal-area map: a
# rectangle of area `A` in `(x, y)` always covers solid angle `A * 4π/12`,
# whatever the face and wherever on the face it sits. That is the whole reason
# a *chart* layer is worth having separately from the pixel layer — pixel
# geometry, vertex geometry, and refinement all fall out of evaluating one
# function on the `(x, y)` lattice, with no ring/nest arithmetic in between,
# and lattice points shared between neighbouring pixels come out bit-identical
# (⇒ the tessellation is exact, not merely consistent to rounding).
#
# Provenance: `xyf_to_point`, the pixel-corner order, and the four RING
# helpers are ports of the closed forms in ConservativeRegridding.jl's
# RingGrids extension (`ext/ConservativeRegriddingRingGridsExt/healpix.jl`),
# which in turn transcribes HEALPix's own `xyf2loc` / `xyf2ring` / `pix2xyf`.
# The NESTED maps below are not in that extension and are written here from
# the standard Morton (bit-interleave) convention. Nothing in this file may
# depend on Healpix.jl: the point is a self-contained kernel that the tests can
# then cross-validate *against* Healpix.jl rather than inherit from it.
#
# ## Index conventions
#
# These are fixed by what each index has to interoperate with, so they are not
# uniform — read this block before calling anything here.
#
# - `ix`, `iy` — 0-based face-local lattice coordinates in `0:nside-1`. `x`
#   grows toward the east (increasing φ), `y` toward the west; both grow
#   northward (increasing z). Continuous chart coordinates are `x = ix/nside`.
# - `face` — 0-based, `0:11`. Faces `0:3` are the north polar caps, `4:7` the
#   equatorial belt, `8:11` the south polar caps.
# - RING index — 1-BASED. It doubles as the position of the pixel in a
#   ring-ordered data vector, which is how RingGrids and Healpix.jl both index
#   HEALPix fields, so making it 0-based here would just move an off-by-one to
#   every call site.
# - NESTED id — 0-BASED. This is the EOPF canonical id convention already used
#   by `HealpixKernel.jl` and the `12 * 4^level` id space of `HealpixLookup`
#   (Healpix.jl's own `nest2ring`/`ring2nest` take 1-based pixel numbers, so
#   converting to them means `+ 1`).
#
# Argument order is `(ix, iy, face, nside)` throughout — deliberately the
# reference extension's order, so that a line-by-line diff against
# `ConservativeRegridding.jl` stays trivial as either side evolves. Ids come
# first in the inverse direction (`ring_to_xyf(ipix, nside)`), again matching
# the reference.
# ---------------------------------------------------------------------------

import GeometryOps as GO

# Górski face constants, indexed by 0-based face number `f` as `JRLL[f + 1]`.
# `JRLL` is the face's row in the 3-row base tiling (2 = north cap, 3 =
# equatorial belt, 4 = south cap) and `JPLL` its column in units of 45° of
# longitude; together they place each unit square on the sphere.
const JRLL = (2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4)
const JPLL = (1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7)

"""
    xyf_to_point(x, y, face) -> GO.UnitSphericalPoint

Evaluate the HEALPix chart of `face` (0-based, `0:11`) at continuous face
coordinates `(x, y) ∈ [0, 1]²`, returning the corresponding point on the unit
sphere.

This is HEALPix's `xyf2loc` (Górski et al. 2005, DOI:10.1086/427976 §4),
ported from ConservativeRegridding.jl's RingGrids extension. It is defined for
*any* real `x, y` in the unit square — nothing here is quantised to a lattice
or restricted to power-of-two `nside`, which is exactly what makes it usable
as the chart underlying refinement of arbitrary depth.

The two branches are the two regimes of the HEALPix projection: within the
equatorial belt `|z| ≤ 2/3` the map is Lambert cylindrical (`z` linear in the
chart coordinate, φ linear too), while in the polar caps it is the "collapsing
meridian" regime where each ring shrinks toward the pole and φ is rescaled by
the ring's own width `nr`.

The `have_sintheta` branch matters for accuracy, not correctness: near a pole
`z → ±1`, so recovering `sinθ` as `sqrt((1-z)(1+z))` loses most of its
significant digits to cancellation. In the cap branches `tmp = nr²/3` is
already `1 - |z|` computed *without* cancellation, so `sqrt(tmp * (2 - tmp))`
is the accurate form and is used whenever `|z| > 0.99`.
"""
function xyf_to_point(x::Real, y::Real, face::Integer)
    xf = Float64(x)
    yf = Float64(y)
    # `jr` is the continuous ring coordinate, 0 at the north pole and 4 at the
    # south pole; the base face's row `JRLL` sets the offset.
    jr = JRLL[face + 1] - xf - yf
    sintheta = 0.0
    have_sintheta = false
    if jr < 1                                       # north polar cap
        nr = jr
        tmp = nr * nr / 3
        z = 1 - tmp
        if z > 0.99
            sintheta = sqrt(tmp * (2 - tmp))
            have_sintheta = true
        end
    elseif jr > 3                                   # south polar cap
        nr = 4 - jr
        tmp = nr * nr / 3
        z = tmp - 1
        if z < -0.99
            sintheta = sqrt(tmp * (2 - tmp))
            have_sintheta = true
        end
    else                                            # equatorial belt
        nr = 1.0
        z = (2 - jr) * (2 / 3)
    end
    # `t` is longitude in units of 45°, wrapped into `[0, 8)`. In the caps it
    # is divided by the ring width `nr`, which is what fans the shrinking rings
    # back out over the full 2π.
    t = JPLL[face + 1] * nr + xf - yf
    t < 0 && (t += 8)
    t >= 8 && (t -= 8)
    # Exactly at a pole `nr == 0` and φ is undefined; 0 is the conventional
    # choice and `sinθ == 0` makes it irrelevant to the returned point.
    phi = nr < 1e-15 ? 0.0 : (Float64(π) / 4 * t) / nr
    st = have_sintheta ? sintheta : sqrt((1 - z) * (1 + z))
    return GO.UnitSphericalPoint(st * cos(phi), st * sin(phi), z)
end

"""
    pixel_corners(ix, iy, face, nside) -> NTuple{4, GO.UnitSphericalPoint}

The four corners of pixel `(ix, iy)` on `face` at resolution `nside`, as
lattice points evaluated through [`xyf_to_point`](@ref).

The corners are emitted **counter-clockwise as seen from outside the sphere**,
in the order `(north, west, south, east)` — that is,
`((ix+1)/n, (iy+1)/n)`, `(ix/n, (iy+1)/n)`, `(ix/n, iy/n)`, `((ix+1)/n, iy/n)`.
CCW is a hard contract, not a stylistic preference: the convex-clip kernel
used for spherical intersection clips a clockwise ring to EMPTY, so a reversed
ring silently produces zero intersection area instead of an error.

Because every corner is a lattice point evaluated by the same function, two
pixels sharing a lattice corner get *bit-identical* points and the
tessellation is exact. Note that HEALPix pixel edges are not great circles
(they follow constant-`z` and constant-φ chart lines), so the four corners
describe the pixel only to the accuracy of a 4-gon; densify through
[`xyf_to_point`](@ref) along the edges if more is needed.
"""
function pixel_corners(ix::Integer, iy::Integer, face::Integer, nside::Integer)
    n = nside
    return (
        xyf_to_point((ix + 1) / n, (iy + 1) / n, face),
        xyf_to_point( ix      / n, (iy + 1) / n, face),
        xyf_to_point( ix      / n,  iy      / n, face),
        xyf_to_point((ix + 1) / n,  iy      / n, face),
    )
end

"""
    pixel_center(ix, iy, face, nside) -> GO.UnitSphericalPoint

The center of pixel `(ix, iy)` on `face` at resolution `nside`: the chart
evaluated at the lattice cell's midpoint `((ix + 0.5)/nside, (iy + 0.5)/nside)`.

This is the HEALPix pixel center by definition (the chart is equal-area, so
the chart midpoint is the canonical center), and agrees with Healpix.jl's
`pix2vecNest`/`pix2vecRing` to round-off.
"""
pixel_center(ix::Integer, iy::Integer, face::Integer, nside::Integer) =
    xyf_to_point((ix + 0.5) / nside, (iy + 0.5) / nside, face)

# ---------------------------------------------------------------------------
# RING order
#
# The ring layout numbers pixels north→south along iso-latitude rings, and
# west→east within a ring. Ring `jr` (1-based, `1:4nside-1`) holds `4jr` pixels
# in the north cap, `4nside` through the belt, mirrored in the south. All of
# it is closed-form arithmetic, valid for ANY `nside >= 1` — the power-of-two
# restriction belongs to the *nested* index only.
# ---------------------------------------------------------------------------

"""
    ring_nlon(jr, nside) -> Int

Number of pixels on global ring `jr` (1-based, north→south): `4jr` while the
cap is still growing, `4nside` across the equatorial belt, mirrored in the
south cap.
"""
ring_nlon(jr::Integer, nside::Integer) =
    jr < nside ? 4jr : (jr <= 3nside ? 4nside : 4 * (4nside - jr))

"""
    ring_first(jr, nside) -> Int

1-based RING index of the first (westernmost) pixel of global ring `jr`.

Each branch is the closed-form prefix sum of [`ring_nlon`](@ref): the north
cap accumulates `4 + 8 + … = 2jr(jr-1)`, the belt adds whole `4nside` rings on
top of the `ncap = 2nside(nside-1)` cap pixels, and the south cap is counted
backwards from the total `12nside²`.
"""
function ring_first(jr::Integer, nside::Integer)
    ncap = 2nside * (nside - 1)
    if jr < nside
        return 2jr * (jr - 1) + 1
    elseif jr <= 3nside
        return ncap + (jr - nside) * 4nside + 1
    else
        js = 4nside - jr                            # mirrored cap ring index
        return 12nside^2 - 2js * (js + 1) + 1
    end
end

"""
    xyf_to_ring(ix, iy, face, nside) -> Int

RING index (1-based, i.e. the position in a ring-ordered data vector) of pixel
`(ix, iy)` on `face` at resolution `nside`. Valid for any `nside >= 1`.

Port of HEALPix's `xyf2ring`. The ring number `jr` follows straight from the
face row `JRLL`; `kshift` is the half-pixel stagger that alternate belt rings
carry (the belt rings are offset by half a pixel from one another, which is
what makes the pixel centers lie on a proper lattice rather than a grid), and
`jp` is the pixel's 1-based position within its ring.

Inputs are assumed lattice-valid (`0 <= ix, iy < nside`, `0 <= face <= 11`) and
are not checked — garbage in, garbage out. This matches the reference extension
and keeps the cursor hot path branch-free.
"""
function xyf_to_ring(ix::Integer, iy::Integer, face::Integer, nside::Integer)
    jr = JRLL[face + 1] * nside - ix - iy - 1
    nr = ring_nlon(jr, nside) >> 2                  # pixels per quadrant of this ring
    kshift = (jr < nside || jr > 3nside) ? 0 : ((jr + nside) & 1)
    jp = (JPLL[face + 1] * nr + ix - iy + 1 + kshift) ÷ 2
    jp < 1 && (jp += 4nside)                        # wrap across the φ = 0 seam
    return ring_first(jr, nside) - 1 + jp
end

"""
    ring_to_xyf(ipix, nside) -> (ix, iy, face)

Inverse of [`xyf_to_ring`](@ref): the face-local lattice coordinates and
0-based face of the pixel at 1-based RING index `ipix`. Valid for any
`nside >= 1`.

Port of HEALPix's `pix2xyf`. Each branch first recovers `(iring, iphi)` — ring
number and 1-based position within the ring — then converts to `(ix, iy, face)`
through the same `JRLL`/`JPLL` relations `xyf_to_ring` uses. In the cap
branches the ring number comes from inverting the triangular number
`2jr(jr-1)` with `isqrt`; the belt branch has to *derive* the face from the
two candidate faces `ifp`/`ifm` that a belt ring straddles (`ifp == ifm` means
the pixel is on an equatorial face, otherwise it belongs to the north or south
cap face depending on which candidate is smaller).
"""
function ring_to_xyf(ipix::Integer, nside::Integer)
    ncap = 2nside * (nside - 1)
    npix = 12nside^2
    n2 = 2nside
    1 <= ipix <= npix || throw(ArgumentError(
        "ring index $ipix out of range for nside=$nside (expected 1:$npix)"))
    if ipix <= ncap                                 # north polar cap
        jr = (1 + isqrt(2ipix - 1)) >> 1
        iphi = ipix - 2jr * (jr - 1)
        kshift = 0
        nr = jr
        iring = jr
        face = (iphi - 1) ÷ nr
    elseif ipix <= npix - ncap                      # equatorial belt
        ip = ipix - 1 - ncap
        tmp = ip ÷ (4nside)
        iring = tmp + nside
        iphi = ip - tmp * 4nside + 1
        kshift = (iring + nside) & 1
        nr = nside
        ire = iring - nside + 1
        irm = n2 + 2 - ire
        ifm = (iphi - ire ÷ 2 + nside - 1) ÷ nside  # candidate face from the north
        ifp = (iphi - irm ÷ 2 + nside - 1) ÷ nside  # candidate face from the south
        face = ifp == ifm ? (ifp | 4) : (ifp < ifm ? ifp : ifm + 8)
    else                                            # south polar cap
        ip = npix - ipix + 1
        jr2 = (1 + isqrt(2ip - 1)) >> 1
        iphi = 4jr2 + 1 - (ip - 2jr2 * (jr2 - 1))
        kshift = 0
        nr = jr2
        iring = 4nside - jr2
        face = 8 + (iphi - 1) ÷ nr
    end
    irt = iring - JRLL[face + 1] * nside + 1
    ipt = 2iphi - JPLL[face + 1] * nr - kshift - 1
    ipt >= n2 && (ipt -= 8nside)                    # unwrap the φ = 0 seam
    return (Int((ipt - irt) >> 1), Int((-ipt - irt) >> 1), Int(face))
end

# ---------------------------------------------------------------------------
# NESTED order
#
# The nested id is a Morton (Z-order) code within a face: it exists precisely
# so that a pixel's four children at the next level are `4p:4p+3`, which needs
# the face to be a `2^k × 2^k` quadtree. Hence — unlike the RING maps above —
# these two functions are only defined for `nside == 2^k` and say so loudly.
# ---------------------------------------------------------------------------

"""
    xyf_to_nested(ix, iy, face, nside) -> Int64

0-based NESTED (EOPF canonical) id of pixel `(ix, iy)` on `face` at resolution
`nside`, which must be a power of two — throws `ArgumentError` otherwise.

The id is `face * nside² + morton(ix, iy)`, where the Morton code places the
bits of `ix` in the even positions and the bits of `iy` in the odd positions.
That is what makes the quadtree work: dropping the low two bits of the code
drops one bit of `ix` and one of `iy`, i.e. steps exactly one level up, so
`p ÷ 4` is the parent and `4p:4p+3` the children — the id arithmetic
`HealpixKernel.jl` gets for free from `has_ordinal_ids`.

The interleave is an explicit bit loop rather than the magic-mask "spread"
trick; both are allocation-free, and at `nside ≤ 2^29` the loop runs at most
29 iterations while staying obviously correct.
"""
function xyf_to_nested(ix::Integer, iy::Integer, face::Integer, nside::Integer)
    ispow2(nside) || throw(ArgumentError(
        "the NESTED index is only defined for nside = 2^k, got nside=$nside; \
         use `xyf_to_ring` for arbitrary nside"))
    (0 <= ix < nside && 0 <= iy < nside) || throw(ArgumentError(
        "lattice coordinates ($ix, $iy) out of range for nside=$nside (expected 0:$(nside - 1))"))
    0 <= face <= 11 || throw(ArgumentError("face $face out of range (expected 0:11)"))
    x = Int64(ix)
    y = Int64(iy)
    morton = Int64(0)
    shift = 0
    while (x | y) != 0
        morton |= ((x & 1) << shift) | ((y & 1) << (shift + 1))
        x >>= 1
        y >>= 1
        shift += 2
    end
    return Int64(face) * Int64(nside)^2 + morton
end

"""
    nested_to_xyf(p, nside) -> (ix, iy, face)

Inverse of [`xyf_to_nested`](@ref): face-local lattice coordinates and 0-based
face of the 0-based NESTED id `p`. `nside` must be a power of two — throws
`ArgumentError` otherwise.

De-interleaves the within-face Morton code: even bits rebuild `ix`, odd bits
rebuild `iy`.
"""
function nested_to_xyf(p::Integer, nside::Integer)
    ispow2(nside) || throw(ArgumentError(
        "the NESTED index is only defined for nside = 2^k, got nside=$nside; \
         use `ring_to_xyf` for arbitrary nside"))
    npface = Int64(nside)^2
    pid = Int64(p)
    0 <= pid < 12npface || throw(ArgumentError(
        "nested id $pid out of range for nside=$nside (expected 0:$(12npface - 1))"))
    face = pid ÷ npface
    code = pid - face * npface
    ix = Int64(0)
    iy = Int64(0)
    shift = 0
    while code != 0
        ix |= (code & 1) << shift
        iy |= ((code >> 1) & 1) << shift
        code >>= 2
        shift += 1
    end
    return (Int(ix), Int(iy), Int(face))
end

"""
    nested_to_ring(p, nside) -> Int

0-based NESTED id → 1-based RING index, via the shared `(ix, iy, face)` chart
coordinates. Equivalent to Healpix.jl's `nest2ring(res, p + 1)`.
"""
nested_to_ring(p::Integer, nside::Integer) =
    xyf_to_ring(nested_to_xyf(p, nside)..., nside)

"""
    ring_to_nested(ipix, nside) -> Int64

1-based RING index → 0-based NESTED id, via the shared `(ix, iy, face)` chart
coordinates. Equivalent to Healpix.jl's `ring2nest(res, ipix) - 1`.
"""
ring_to_nested(ipix::Integer, nside::Integer) =
    xyf_to_nested(ring_to_xyf(ipix, nside)..., nside)
