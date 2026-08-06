# ---------------------------------------------------------------------------
# ISEA9R's chart IS ISEA4R's — the delegation record, and the base-9 index maps
#
# This file sits in the slot `src/ISEA4R/chart.jl`, `src/S2/chart.jl` and
# `src/HEALPix/chart.jl` occupy for their systems, and it answers the same
# question — *where does this system's `[0, 1]² → S²` come from?* — with an
# import rather than a definition. ISEA9R is the ten-diamond rhombus chart at
# `nside = 3^level`; ISEA4R is the same chart at `nside = 2^level`. There is one
# chart, and it already ships.
#
# ## Why delegation is the correct shape here, not a smell
#
# The rhombus chart is **aperture-agnostic by construction**: `xyd_to_point(x,
# y, diamond)` takes *continuous* coordinates and quantises nothing, so a
# resolution is only ever a divisor applied by the caller. `docs/design/
# isea4r_diamond_layout.md` §7 posed exactly this reuse as the cheap branch of
# the ISEA9R question, and `docs/design/isea9r_layout.md` records the
# primary-source resolution of that question (OGC 21-038r1 Annex B.2: *"The ten
# root rhombuses are formed by combining two icosahedron triangles at their
# base"*; DGGAL `src/dggrs/RI9R.ec` `countZones(level) = 10 * 9^level`). The
# ten-diamond layout is ISEA9R's layout, so re-deriving the `DIAMONDS` table
# here would produce a second copy of one table, with two chances to drift and
# an `@assert` that only guards one of them.
#
# Hoisting the shared chart into a common submodule (an `ISEARhombic`, or into
# `ISEA` itself) is a live option and is **deliberately deferred**: the layout
# is not oracle-locked and `ISEA` is charter-bound to oracle-locked machinery
# (see the `ISEA4R` module docstring), so the hoist is a design decision of its
# own rather than a side effect of adding an aperture. Until it happens, a plain
# cross-submodule import is the honest expression of "one chart, two apertures"
# — and it is *checkable*: `test/ISEA9R/test_delegation.jl` asserts the function
# objects are `===` identical and that the two systems' polygons at a common
# `nside = 3^k` are bitwise equal.
#
# What is imported is the whole geometry surface, unchanged:
#
#   `xyd_to_point`     the chart itself, seam-ownership rule included
#   `cell_corners`     the four CCW corners of a lattice cell
#   `cell_center`      the chart at the lattice cell's midpoint
#   `xyd_to_rowmajor` / `rowmajor_to_xyd`
#                      the row-major index maps, which carry no aperture at all
#                      (`diamond * nside² + iy * nside + ix`, any `nside >= 1`)
#   `DIAMONDS`         the ten-diamond layout table
#
# What is *added* here is the one thing aperture 9 needs and aperture 4 cannot
# supply: the **base-9 Morton** index maps below.
#
# ## Index conventions
#
# Identical to `src/ISEA4R/chart.jl`, which is the point:
#
# - `ix`, `iy` — 0-based diamond-local lattice coordinates in `0:nside-1`.
#   Continuous chart coordinates are `x = ix/nside`, `y = iy/nside`.
# - `diamond` — 0-based, `0:9`, the `ISEA4R.DIAMONDS` numbering. Diamonds `0:4`
#   are northern (apex at base 0), `5:9` southern (apex at base 11).
# - Row-major id — 0-BASED: `diamond * nside² + iy * nside + ix`.
# - Morton id — 0-BASED: `diamond * nside² + morton9(ix, iy)`, which is exactly
#   the ordinal shape the `ISEA9RDGGS` registry entry records
#   (`diamond * 9^level + position`, see `src/core/systems/isea4r_isea9r.jl`).
#
# Both id spaces are 0-based and data positions are `id + 1` for BOTH orderings.
# Argument order is `(ix, iy, diamond, nside)` throughout and ids come first in
# the inverse direction, matching the sibling charts slot for slot.
#
# ## What is NOT claimed
#
# The ten-diamond numbering is this package's own convention with no external
# oracle behind it, and the *ordering* below is a package-canonical choice that
# DGGAL does not share (DGGAL's within-rhombus index is row-major over a
# transposed square under a non-identity root permutation — see
# `docs/design/isea9r_layout.md` §6). **No DGGRID / DGGAL / SST identifier or
# geometry compatibility is claimed and none may be inferred.**
# ---------------------------------------------------------------------------

# The chart, imported rather than redefined — see the file header. `import`
# rather than `using` so the delegation is spelled out name by name and shows up
# in a grep for `ISEA4R` from this directory.
import ..ISEA4R: xyd_to_point, cell_corners, cell_center, xyd_to_rowmajor,
    rowmajor_to_xyd, DIAMONDS

# ---------------------------------------------------------------------------
# Powers of three
#
# `Base.ispow2` has no base-3 twin, and the aperture-9 orderings need one in
# three places (both codecs and `validate_ordering`). Written as a division loop
# rather than as `3^round(Int, log(3, n)) == n`, which is a floating-point
# question about an integer.
# ---------------------------------------------------------------------------

"""
    ispow3(n) -> Bool

Whether `n` is `3^k` for some integer `k >= 0`. The base-3 counterpart of
`Base.ispow2`, which is what the aperture-9 index maps test `nside` against
(`ISEA4R`'s aperture-4 ones call `ispow2` in the same places).

`ispow3(0)` and `ispow3(-3)` are `false`: a resolution is a positive count.
"""
function ispow3(n::Integer)
    n >= 1 || return false
    m = Int64(n)
    while m % 3 == 0
        m ÷= 3
    end
    return m == 1
end

# ---------------------------------------------------------------------------
# Base-9 Morton order
#
# The aperture-9 analogue of `ISEA4R`'s aperture-4 Morton code, and the reason
# it is a *radix* question rather than a bit question: a Z-order code for
# aperture `a²` interleaves one base-`a` DIGIT of `ix` with one of `iy` per
# level. For `a = 2` the digits are bits and the interleave is a shift-and-mask
# loop (`ISEA4R.xyd_to_morton`); for `a = 3` they are base-3 digits and the
# interleave is the arithmetic below. Everything else about the two is the same,
# including the property that makes them worth having:
#
#   dropping the low base-9 digit of the code drops one base-3 digit of `ix` and
#   one of `iy`, i.e. steps exactly one level up,
#
# so `p ÷ 9` is the parent and `9p:9p+8` the children — the radix-9 prefix
# arithmetic the `ISEA9RDGGS` registry entry records as
# `diamond * 9^level + position`, and the reason
# `supports_prefix_ranges(ISEA9RDGGS())` is `true` against THIS ordinal.
#
# Digit convention: base-9 digit `k` of the code is `ix_k + 3 * iy_k`, i.e. `ix`
# in the low base-3 slot and `iy` in the high one — the exact analogue of
# `ISEA4R`'s "`ix` in the even bit positions, `iy` in the odd ones", since a
# base-9 digit is a base-3 digit pair.
# ---------------------------------------------------------------------------

"""
    xyd_to_morton(ix, iy, diamond, nside) -> Int64

0-based base-9 Morton id of cell `(ix, iy)` on `diamond` at resolution `nside`,
which must be a power of three — throws `ArgumentError` otherwise.

The id is `diamond * nside² + morton9(ix, iy)`, where `morton9` interleaves the
*base-3 digits* of `ix` and `iy`: base-9 digit `k` of the code is
`ix_k + 3 * iy_k`. That is what makes the aperture-9 quadtree work — dropping
the low base-9 digit drops one base-3 digit of each coordinate, i.e. steps
exactly one level up, so `p ÷ 9` is the parent and `9p:9p+8` the children, which
is the radix-9 prefix arithmetic the `ISEA9RDGGS` registry entry records as
`diamond * 9^level + position`.

The digit loop is the base-3 twin of [`ISEA4R.xyd_to_morton`](@ref)'s bit loop
(there, `ix` goes to the even bit positions and `iy` to the odd ones — the same
statement, since a base-9 digit is a base-3 digit pair). At `nside <= 3^18` it
runs at most 18 iterations.

This ordering is **package-canonical and is not DGGAL's**: DGGAL indexes within
a root rhombus in row-major order over a transposed square, under a non-identity
root→diamond permutation. See `docs/design/isea9r_layout.md`.
"""
function xyd_to_morton(ix::Integer, iy::Integer, diamond::Integer, nside::Integer)
    ispow3(nside) || throw(ArgumentError(
        "the base-9 Morton index is only defined for nside = 3^k, got nside=$nside; \
         use `xyd_to_rowmajor` for arbitrary nside"))
    (0 <= ix < nside && 0 <= iy < nside) || throw(ArgumentError(
        "lattice coordinates ($ix, $iy) out of range for nside=$nside (expected 0:$(nside - 1))"))
    0 <= diamond <= 9 || throw(ArgumentError("diamond $diamond out of range (expected 0:9)"))
    x = Int64(ix)
    y = Int64(iy)
    morton = Int64(0)
    scale = Int64(1)                                  # 9^k, the digit's place value
    while (x | y) != 0
        morton += ((x % 3) + 3 * (y % 3)) * scale
        x ÷= 3
        y ÷= 3
        scale *= 9
    end
    return Int64(diamond) * Int64(nside)^2 + morton
end

"""
    morton_to_xyd(p, nside) -> (ix, iy, diamond)

Inverse of [`xyd_to_morton`](@ref): diamond-local lattice coordinates and
0-based diamond of the 0-based base-9 Morton id `p`. `nside` must be a power of
three — throws `ArgumentError` otherwise, as it does for an id outside
`0:10nside²-1`.

De-interleaves the within-diamond base-9 code: the low base-3 digit of each
base-9 digit rebuilds `ix`, the high one rebuilds `iy`.
"""
function morton_to_xyd(p::Integer, nside::Integer)
    ispow3(nside) || throw(ArgumentError(
        "the base-9 Morton index is only defined for nside = 3^k, got nside=$nside; \
         use `rowmajor_to_xyd` for arbitrary nside"))
    npd = Int64(nside)^2
    pid = Int64(p)
    0 <= pid < 10npd || throw(ArgumentError(
        "base-9 Morton id $pid out of range for nside=$nside (expected 0:$(10npd - 1))"))
    diamond = pid ÷ npd
    code = pid - diamond * npd
    ix = Int64(0)
    iy = Int64(0)
    scale = Int64(1)                                  # 3^k, the digit's place value
    while code != 0
        digit = code % 9
        ix += (digit % 3) * scale
        iy += (digit ÷ 3) * scale
        code ÷= 9
        scale *= 3
    end
    return (Int(ix), Int(iy), Int(diamond))
end
