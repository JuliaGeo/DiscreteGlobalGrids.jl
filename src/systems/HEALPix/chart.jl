# HEALPix (Górski et al. 2005, DOI:10.1086/427976) comprises twelve continuous charts
# `[0, 1]² → S²`, one per base face. Every chart is an equal-area map: a
# rectangle of area `A` in `(x, y)` always covers solid angle `A * 4π/12`,
# whatever the face or location. Shared lattice points are evaluated by the
# same function and are bit-identical.
#
# ## Index conventions
#
# - `ix`, `iy` — 0-based face-local lattice coordinates in `0:nside-1`. `x`
#   grows toward the east (increasing φ), `y` toward the west; both grow
#   northward (increasing z). Continuous chart coordinates are `x = ix/nside`.
# - `face` — 0-based, `0:11`. Faces `0:3` are the north polar caps, `4:7` the
#   equatorial belt, `8:11` the south polar caps.
# - RING indices are 1-based; NESTED ids are 0-based.
# Argument order is `(ix, iy, face, nside)`; inverse codecs take the id first.

# Górski face constants, indexed by 0-based face number `f` as `JRLL[f + 1]`.
# `JRLL` is the face's row in the 3-row base tiling (2 = north cap, 3 =
# equatorial belt, 4 = south cap) and `JPLL` its column in units of 45° of
# longitude; together they place each unit square on the sphere.
const JRLL = (2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4)
const JPLL = (1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7)

"""
    xyf_to_point(x, y, face) -> GO.UnitSphericalPoint

Evaluate the HEALPix chart for 0-based `face` at continuous
`(x, y) ∈ [0, 1]²`.

The equatorial branch is Lambert cylindrical; polar rings shrink their
longitude span toward the pole.

For `|z| > 0.99`, `sinθ` uses `sqrt(tmp * (2 - tmp))` to avoid cancellation in
`sqrt((1-z)(1+z))`.
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

Corners are counter-clockwise seen from outside in `(north, west, south, east)` order:
`((ix+1)/n, (iy+1)/n)`, `(ix/n, (iy+1)/n)`, `(ix/n, iy/n)`, `((ix+1)/n, iy/n)`.
The winding is contractual: the spherical convex-clip kernel clips a clockwise
ring to empty, so a reversed ring silently yields zero intersection area rather
than an error. Shared lattice corners are bit-identical, and carry no negative
zero, so they are safe as `Set`/`Dict` keys — those compare with `isequal`,
under which `-0.0` and `0.0` differ. Pixel edges are chart curves rather than
great circles, so four corners alone are only a polygonal approximation.
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
# it is closed-form arithmetic valid for any `nside >= 1`. Only the nested
# index requires a power-of-two side length.
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
        return 12 * nside * nside - 2js * (js + 1) + 1
    end
end

"""
    xyf_to_ring(ix, iy, face, nside) -> Int

RING index (1-based, i.e. the index in a ring-ordered data vector) of pixel
`(ix, iy)` on `face` at resolution `nside`. Valid for any `nside >= 1`.

Implements HEALPix `xyf2ring`; `kshift` is the alternating half-pixel belt
stagger and `jp` is the 1-based index within the ring.

Inputs must satisfy `0 <= ix, iy < nside` and `0 <= face <= 11`; unlike the
NESTED codecs, they are not checked.
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

Implements HEALPix `pix2xyf` by recovering `(iring, iphi)` and applying the
`JRLL`/`JPLL` relations. Polar-ring triangular counts are inverted with `isqrt`.
"""
function ring_to_xyf(ipix::Integer, nside::Integer)
    ncap = 2nside * (nside - 1)
    npix = 12 * nside * nside
    n2 = 2nside
    1 <= ipix <= npix || _throw_bad_ring(ipix, nside, npix)
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
# the face to be a `2^k × 2^k` quadtree. These functions therefore require
# `nside == 2^k`.
# ---------------------------------------------------------------------------

"""
    xyf_to_nested(ix, iy, face, nside) -> Int64

0-based NESTED (EOPF canonical) id of pixel `(ix, iy)` on `face` at resolution
`nside`, which must be a power of two — throws `ArgumentError` otherwise.

The id is `face * nside² + morton(ix, iy)`, with `ix` bits in even positions and
`iy` bits in odd positions. Thus `p ÷ 4` is the parent and `4p:4p+3` are the
children.
"""
function xyf_to_nested(ix::Integer, iy::Integer, face::Integer, nside::Integer)
    ispow2(nside) || _throw_not_pow2(nside, "xyf_to_ring")
    (0 <= ix < nside && 0 <= iy < nside) || _throw_bad_xy(ix, iy, nside)
    0 <= face <= 11 || _throw_bad_face(face)
    n = Int64(nside)
    return Int64(face) * (n * n) + DGG.morton_encode(ix, iy)
end

# The argument-error paths build interpolated strings, which is far more code
# than the arithmetic they guard. Kept behind `@noinline` so the hot path stays
# small enough for the compiler to inline it into callers.
@noinline _throw_not_pow2(nside, alt) = throw(ArgumentError(
    "the NESTED index is only defined for nside = 2^k, got nside=$nside; \
     use `$alt` for arbitrary nside"))
@noinline _throw_bad_xy(ix, iy, nside) = throw(ArgumentError(
    "lattice coordinates ($ix, $iy) out of range for nside=$nside (expected 0:$(nside - 1))"))
@noinline _throw_bad_face(face) = throw(ArgumentError("face $face out of range (expected 0:11)"))
@noinline _throw_bad_nested(pid, nside, npface) = throw(ArgumentError(
    "nested id $pid out of range for nside=$nside (expected 0:$(12npface - 1))"))
@noinline _throw_bad_ring(ipix, nside, npix) = throw(ArgumentError(
    "ring index $ipix out of range for nside=$nside (expected 1:$npix)"))

"""
    nested_to_xyf(p, nside) -> (ix, iy, face)

Inverse of [`xyf_to_nested`](@ref): face-local lattice coordinates and 0-based
face of the 0-based NESTED id `p`. `nside` must be a power of two — throws
`ArgumentError` otherwise.

De-interleaves the within-face Morton code: even bits rebuild `ix`, odd bits
rebuild `iy`.
"""
function nested_to_xyf(p::Integer, nside::Integer)
    ispow2(nside) || _throw_not_pow2(nside, "ring_to_xyf")
    n = Int64(nside)
    npface = n * n
    pid = Int64(p)
    0 <= pid < 12npface || _throw_bad_nested(pid, nside, npface)
    # `nside` is a power of two, so the face quotient and within-face remainder
    # are a shift and a mask rather than a 64-bit division.
    shift = 2 * trailing_zeros(n)
    face = pid >> shift
    code = pid & (npface - 1)
    ix, iy = DGG.morton_decode(code)
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

# Closed-form inverse chart using HEALPix `loc2xyf`. Its angular convention
# differs from `xyf_to_point`: there `t` is longitude in
# units of 45 degrees (`[0, 8)`); here `tt` is longitude in units of 90 degrees
# (`[0, 4)`).

"""
    point_to_xyf(p, nside, order) -> (ix, iy, face)

Return the face-local lattice coordinates and 0-based face containing `p` for
`nside == 2^order`. This O(1) inverse quantises the HEALPix chart without tree
descent or polygon tests. Above `|z| = 0.99`, it derives polar ring width from
the equatorial radius to avoid cancellation in `1 - |z|`. Boundary ownership
uses the deterministic higher-side `floor` convention and may differ from
other HEALPix implementations.
"""
function point_to_xyf(p, nside::Integer, order::Integer)
    z = Float64(p[3])
    za = abs(z)
    # sinθ from the equatorial radius rather than `sqrt(1 - z^2)`: exact for a
    # unit-norm point, and free of the cancellation that form suffers near a pole.
    st = sqrt(Float64(p[1])^2 + Float64(p[2])^2)
    phi = atan(Float64(p[2]), Float64(p[1]))
    # longitude in units of 90 degrees, [0, 4). `atan` returns `(-π, π]`, so the
    # scaled value is in `(-2, 2]` and the wrap is a single compare-and-add --
    # bit-identical to `mod(·, 4.0)` over that range, without the `fmod` call.
    # The `+ 0.0` reproduces `mod`'s normalisation of `-0.0` to `+0.0`.
    tt = phi * (2 / Float64(π))
    tt = tt < 0 ? tt + 4.0 : tt + 0.0
    n = Int64(nside)

    if za <= 2 / 3                              # equatorial belt
        temp1 = n * (0.5 + tt)
        temp2 = n * (z * 0.75)
        jp = floor(Int64, temp1 - temp2)        # index of the ascending edge line
        jm = floor(Int64, temp1 + temp2)        # index of the descending edge line
        ifp = jp >> order                       # face candidate from the north
        ifm = jm >> order                       # face candidate from the south
        face = ifp == ifm ? (ifp | 4) : (ifp < ifm ? ifp : ifm + 8)
        ix = jm & (n - 1)
        iy = n - (jp & (n - 1)) - 1
    else                                        # polar caps
        ntt = min(3, floor(Int, tt))
        tp = tt - ntt
        tmp = za < 0.99 ? n * sqrt(3 * (1 - za)) : n * st / sqrt((1 + za) / 3)
        jp = min(floor(Int64, tp * tmp), n - 1)
        jm = min(floor(Int64, (1 - tp) * tmp), n - 1)
        if z >= 0
            face = Int64(ntt)
            ix = n - jm - 1
            iy = n - jp - 1
        else
            face = Int64(ntt) + 8
            ix = jp
            iy = jm
        end
    end
    return (Int(ix), Int(iy), Int(face))
end

"""
    point_to_nested(p, level) -> Int64

The 0-based NESTED id of the pixel containing unit-sphere point `p` at
refinement `level` (`nside = 2^level`). [`point_to_xyf`](@ref) composed with
[`xyf_to_nested`](@ref); equivalent to `Healpix.ang2pixNest(res, θ, φ) - 1`.
"""
function point_to_nested(p, level::Integer)
    nside = 1 << Int(level)
    return xyf_to_nested(point_to_xyf(p, nside, Int(level))..., nside)
end
