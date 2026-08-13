# ---------------------------------------------------------------------------
# ISEA4R diamond charts — the rhombus chart over ISEA's Snyder machinery
#
# ISEA4R is ten continuous charts `[0, 1]² → S²`, one per icosahedron diamond
# (`diamonds.jl`). Every chart is an equal-area map: a rectangle of area `A` in
# `(x, y)` always covers solid angle `A · 4π/10`, whatever the diamond and
# wherever in it the rectangle sits. That is exact, not asymptotic — Snyder ISEA
# is exactly equal-area per face and the two affine halves have
# `|det| == 2π/5 == 4π/10` to the last bit (see [`Diamond`](@ref)).
#
# That is the whole reason a *chart* layer is worth having separately from any
# id layer: cell geometry, vertex geometry, and refinement all fall out of
# evaluating one function on the `(x, y)` lattice, with no ordinal arithmetic in
# between, and lattice points shared between neighbouring cells come out
# bit-identical (⇒ the tessellation is exact, not merely consistent to
# rounding). It is also what lets the grid exist at `nside = 3`, `5`, `7`, which
# have no aperture-4 id space at all.
#
# Provenance: the projection is entirely `ISEA`'s — `snyder_inv_xyz` and the
# face frames, unchanged. What is added here is the layout (`diamonds.jl`) and
# the piecewise-affine assembly of the two face triangles into one square chart.
# The file is otherwise a deliberate sibling of `src/systems/HEALPix/chart.jl`,
# section for section.
#
# Ported from the pre-redesign `src/ISEA4R/chart.jl`; the closed-form INVERSE
# (`point_to_xyd` / `point_to_morton`) at the foot of the file is new, and is
# what makes `cellat` a chart evaluation rather than a tree descent.
#
# ## Index conventions
#
# Read this block before calling anything here.
#
# - `ix`, `iy` — 0-based diamond-local lattice coordinates in `0:nside-1`.
#   Continuous chart coordinates are `x = ix/nside`, `y = iy/nside`.
# - `diamond` — 0-based, `0:9`, the [`DIAMONDS`](@ref) numbering. Diamonds `0:4`
#   are northern (apex at base 0), `5:9` southern (apex at base 11).
# - Row-major id — 0-BASED: `diamond * nside² + iy * nside + ix`.
# - Morton id — 0-BASED: `diamond * nside² + morton(ix, iy)`, which is exactly
#   the canonical `LevelIndex.index` of `ISEA4RSystem` — `diamond * 4^level +
#   position` with `position` pinned to the Z-order Morton code.
#
# Both id spaces are 0-based and data positions are `id + 1` for BOTH orderings.
# (HEALPix's mixed 1-based-RING / 0-based-NESTED convention exists only because
# RING doubles as the position in a HEALPix FITS field — ISEA4R has no such
# external layout to match.)
#
# Argument order is `(ix, iy, diamond, nside)` throughout and ids come first in
# the inverse direction (`rowmajor_to_xyd(q, nside)`), matching the sibling
# chart slot for slot so the files diff against each other cleanly.
# ---------------------------------------------------------------------------

"""
    xyd_to_point(x, y, diamond) -> GO.UnitSphericalPoint

Evaluate the ISEA4R chart of `diamond` (0-based, `0:9`) at continuous chart
coordinates `(x, y) ∈ [0, 1]²`, returning the corresponding point on the unit
sphere (grid frame — the standard ISEA icosahedron placement, `ISEA`'s
identity `Orientation`).

Defined for *any* real `x, y` in the square — nothing here is quantised to a
lattice or restricted to power-of-two `nside`, which is what makes it usable as
the chart underlying refinement of arbitrary depth.

The two branches are the diamond's two icosahedron faces. Their planar
triangles are glued along the seam `y == x`, and each half is a plain affine
map of the chart square's half-triangle onto its face's Snyder plane
([`Diamond`](@ref)); the Snyder inverse then lifts the planar point to the
sphere. Because the map is affine into an equal-area projection, the chart is
**exactly equal-area**: a chart rectangle of area `A` covers solid angle
`A · 4π/10`, everywhere, on every diamond.

# Seam ownership

The half `y >= x` — including the seam itself — is evaluated through the
*upper* face. That is a decision, not an accident, and it is what makes the
tessellation exact:

* Lattice points with `iy >= ix` (every seam point among them) go through the
  upper branch, so any two cells sharing such a point evaluate it through the
  identical branch of the identical pure function and get *bit-identical*
  results. The two halves agree on the seam only to `~5e-15` rad (they are two
  independent developments of one icosahedron edge — see the seam assertions in
  [`_make_diamonds`](@ref)), so a rule that let one cell take the upper
  development and its neighbour the lower one would leave a gap of that size
  between two polygons that ought to share an edge.
* The continuous predicate `Float64(iy/nside) >= Float64(ix/nside)` coincides
  with the integer predicate `iy >= ix` for all `0 <= ix, iy <= nside <= 2^29`:
  correctly-rounded division is monotone, and it is injective on
  `{0..n}/n` because the spacing `1/n >= 2^-29` dwarfs the ulp. So the branch a
  lattice point takes is decided exactly, and the same point evaluated at
  resolution `n` and at `2n` takes the same branch too — `fl(ix/n)` and
  `fl(2ix/2n)` are the same `Float64`, so cross-resolution nesting is
  bit-exact as well.

# Cell edges are not great circles

Snyder ISEA maps the chart's straight lines to curves on the sphere, so the
four corners of a cell describe it only to the accuracy of a 4-gon — as with
HEALPix's `xyf_to_point`, and unlike S2's `stf_to_point`, whose cell edges *are*
great circles. Densify through `xyd_to_point` along the edges if more is needed.
"""
function xyd_to_point(x::Real, y::Real, diamond::Integer)
    dm = @inbounds DIAMONDS[diamond+1]
    xf = Float64(x)
    yf = Float64(y)
    if yf >= xf                                     # upper half OWNS the seam
        p = ISEA.snyder_inv_xyz(dm.upper, dm.cP0 + xf * dm.aP + yf * dm.bP)
    else
        p = ISEA.snyder_inv_xyz(dm.lower, dm.cQ0 + xf * dm.aQ + yf * dm.bQ)
    end
    return GO.UnitSphericalPoint(p[1], p[2], p[3])
end

# ---------------------------------------------------------------------------
# Cell geometry
#
# Why the corner order below is counterclockwise seen from outside, on all ten
# diamonds (the tests then check it by exhaustion):
#
#  1. The emission order `(x+, y+), (x-, y+), (x-, y-), (x+, y-)` is CCW in the
#     `(x, y)` plane — its shoelace sum is `+2` on the unit square.
#  2. Both affine halves are orientation-preserving: `imag(conj(a) * b)` is
#     `+2π/5`, asserted exactly at build time for all twenty half-maps.
#  3. `snyder_inv_xyz` is orientation-preserving from a face's planar frame onto
#     the sphere seen from outside — the frame is `(u, w = c × u)` with `c` the
#     outward face centre, so `(u, w)` is right-handed about `c`, and the radial
#     lift is a diffeomorphism onto the face's spherical patch.
#
# So a positively-oriented chart rectangle maps to a positively-oriented
# spherical quadrilateral, on either half and hence across the seam too.
# ---------------------------------------------------------------------------

"""
    cell_corners(ix, iy, diamond, nside) -> NTuple{4, GO.UnitSphericalPoint}

The four corners of cell `(ix, iy)` on `diamond` at resolution `nside`, as
lattice points evaluated through [`xyd_to_point`](@ref).

The corners are emitted **counterclockwise as seen from outside the sphere**, in
the order `((ix+1)/n, (iy+1)/n)`, `(ix/n, (iy+1)/n)`, `(ix/n, iy/n)`,
`((ix+1)/n, iy/n)` — the same lattice order as HEALPix's `pixel_corners` and
S2's `cell_corners`. CCW is a hard contract, not a stylistic preference: the
convex-clip kernel used for spherical intersection clips a clockwise ring to
EMPTY, so a reversed ring silently produces zero intersection area instead of an
error. The winding holds on every diamond for the structural reason spelled out
in the file comment above.

Because every corner is a lattice point evaluated by the same function under the
same seam-ownership rule, two cells sharing a lattice corner get *bit-identical*
points and the tessellation is exact — within a diamond, across the seam, and
across resolutions `n` / `2n` alike (see [`xyd_to_point`](@ref)). Note that
Snyder cell edges are NOT great circles, so the four corners describe the cell
only to the accuracy of a 4-gon; densify through [`xyd_to_point`](@ref) along
the edges if more is needed.
"""
function cell_corners(ix::Integer, iy::Integer, diamond::Integer, nside::Integer)
    n = nside
    return (
        xyd_to_point((ix + 1) / n, (iy + 1) / n, diamond),
        xyd_to_point(ix / n, (iy + 1) / n, diamond),
        xyd_to_point(ix / n, iy / n, diamond),
        xyd_to_point((ix + 1) / n, iy / n, diamond),
    )
end

"""
    cell_center(ix, iy, diamond, nside) -> GO.UnitSphericalPoint

The center of cell `(ix, iy)` on `diamond` at resolution `nside`: the chart
evaluated at the lattice cell's midpoint `((ix + 0.5)/nside, (iy + 0.5)/nside)`.

The chart is exactly equal-area ([`xyd_to_point`](@ref)), so the chart midpoint
*is* the canonical center — the same argument HEALPix's `pixel_center` rests on.
(It is not the spherical centroid of the 4-gon; no equal-area DGGS claims that
of its cell centers.)

This is `ISEA4R.cell_center`, a function in the `ISEA4R` namespace; the
interface generic is [`cell_centroid`](@ref), which `system.jl` implements by
calling this.
"""
cell_center(ix::Integer, iy::Integer, diamond::Integer, nside::Integer) =
    xyd_to_point((ix + 0.5) / nside, (iy + 0.5) / nside, diamond)

# ---------------------------------------------------------------------------
# Row-major order
#
# The plain lexicographic layout of the `nside × nside` lattice, `ix` fastest,
# then `iy`, then `diamond`. All closed-form arithmetic, valid for ANY
# `nside >= 1` — the power-of-two restriction belongs to the *Morton* index
# only, exactly as row-major/Hilbert split in `S2/chart.jl` and RING/NESTED in
# `HEALPix/chart.jl`.
#
# Row-major is NOT an id scheme of `ISEA4RSystem` — the canonical id is Morton,
# because only Morton makes `diamond * 4^level + position` a radix-4 prefix
# hierarchy. It is kept here because it is the natural dense layout of an
# `nside`-by-`nside` diamond at arbitrary (non-power-of-two) `nside`, which is
# what the chart itself supports, and because `xyd_to_rowmajor` is the cheapest
# way to enumerate a diamond in scan order in a test.
# ---------------------------------------------------------------------------

"""
    xyd_to_rowmajor(ix, iy, diamond, nside) -> Int64

0-based row-major id of cell `(ix, iy)` on `diamond` at resolution `nside`:
`diamond * nside² + iy * nside + ix`. Valid for any `nside >= 1`.

Inputs are assumed lattice-valid (`0 <= ix, iy < nside`, `0 <= diamond <= 9`)
and are not checked — garbage in, garbage out, matching the sibling charts and
keeping the cursor hot path branch-free. [`rowmajor_to_xyd`](@ref) *is* checked,
being the direction a user-supplied id enters through.
"""
xyd_to_rowmajor(ix::Integer, iy::Integer, diamond::Integer, nside::Integer) =
    Int64(diamond) * Int64(nside)^2 + Int64(iy) * Int64(nside) + Int64(ix)

"""
    rowmajor_to_xyd(q, nside) -> (ix, iy, diamond)

Inverse of [`xyd_to_rowmajor`](@ref): diamond-local lattice coordinates and
0-based diamond of the 0-based row-major id `q`. Valid for any `nside >= 1`;
throws an `ArgumentError` for an id outside `0:10nside²-1`.
"""
function rowmajor_to_xyd(q::Integer, nside::Integer)
    npd = Int64(nside)^2
    qi = Int64(q)
    0 <= qi < 10npd || throw(ArgumentError(
        "row-major id $qi out of range for nside=$nside (expected 0:$(10npd - 1))"))
    diamond, r = divrem(qi, npd)
    iy, ix = divrem(r, Int64(nside))
    return (Int(ix), Int(iy), Int(diamond))
end

# ---------------------------------------------------------------------------
# Morton order
#
# The Morton (Z-order) code is to ISEA4R what it is to HEALPix: the within-face
# ordering under which a cell's four children at the next level are contiguous,
# which is what makes `diamond * 4^r + position` a radix-4 prefix hierarchy. It
# needs a `2^k × 2^k` diamond and says so loudly.
#
# The two functions below are line-for-line ports of
# `src/systems/HEALPix/chart.jl`'s `xyf_to_nested` / `nested_to_xyf` with
# `12 → 10` and `0:11 → 0:9`. That is the whole difference: the code is a
# property of the lattice, not of the sphere.
# ---------------------------------------------------------------------------

"""
    xyd_to_morton(ix, iy, diamond, nside) -> Int64

0-based Morton id of cell `(ix, iy)` on `diamond` at resolution `nside`, which
must be a power of two — throws `ArgumentError` otherwise.

The id is `diamond * nside² + morton(ix, iy)`, where the Morton code places the
bits of `ix` in the even positions and the bits of `iy` in the odd positions.
That is what makes the quadtree work: dropping the low two bits of the code
drops one bit of `ix` and one of `iy`, i.e. steps exactly one level up, so
`p ÷ 4` is the parent and `4p:4p+3` the children — the radix-4 prefix arithmetic
[`ISEA4RSystem`](@ref)'s canonical id is `diamond * 4^level + position`.

The interleave is an explicit bit loop rather than the magic-mask "spread"
trick; both are allocation-free, and at `nside ≤ 2^29` the loop runs at most 29
iterations while staying obviously correct.

This is the CHECKED entry point, which is what a user-supplied coordinate comes
in through. [`xyd_to_morton_unchecked`](@ref) is the same interleave without the
three guards, for callers that produced the coordinates themselves.
"""
function xyd_to_morton(ix::Integer, iy::Integer, diamond::Integer, nside::Integer)
    ispow2(nside) || throw(ArgumentError(
        "the Morton index is only defined for nside = 2^k, got nside=$nside; \
         use `xyd_to_rowmajor` for arbitrary nside"))
    (0 <= ix < nside && 0 <= iy < nside) || throw(ArgumentError(
        "lattice coordinates ($ix, $iy) out of range for nside=$nside (expected 0:$(nside - 1))"))
    0 <= diamond <= 9 || throw(ArgumentError("diamond $diamond out of range (expected 0:9)"))
    return xyd_to_morton_unchecked(ix, iy, diamond, nside)
end

"""
    xyd_to_morton_unchecked(ix, iy, diamond, nside) -> Int64

[`xyd_to_morton`](@ref) with the three validity guards dropped: `nside` a power
of two, `(ix, iy)` on the lattice, `diamond` in `0:9`.

For callers that produced the coordinates themselves and know they are valid —
`lattice_neighbors`' output, which is derived from a validated cell by table
lookup, and the rim walk in `border.jl`, which enumerates a block it computed.
Both run per neighbour and per rim cell respectively, where re-deriving what
the caller already knows is pure overhead. Garbage in, garbage out: an
out-of-range coordinate yields the id of a cell that does not exist rather than
an error, which is exactly why the checked form is the one a user-supplied
coordinate goes through.
"""
@inline function xyd_to_morton_unchecked(ix::Integer, iy::Integer, diamond::Integer,
        nside::Integer)
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
    return Int64(diamond) * Int64(nside)^2 + morton
end

"""
    morton_to_xyd(p, nside) -> (ix, iy, diamond)

Inverse of [`xyd_to_morton`](@ref): diamond-local lattice coordinates and
0-based diamond of the 0-based Morton id `p`. `nside` must be a power of two —
throws `ArgumentError` otherwise.

De-interleaves the within-diamond Morton code: even bits rebuild `ix`, odd bits
rebuild `iy`.
"""
function morton_to_xyd(p::Integer, nside::Integer)
    ispow2(nside) || throw(ArgumentError(
        "the Morton index is only defined for nside = 2^k, got nside=$nside; \
         use `rowmajor_to_xyd` for arbitrary nside"))
    npd = Int64(nside)^2
    pid = Int64(p)
    0 <= pid < 10npd || throw(ArgumentError(
        "Morton id $pid out of range for nside=$nside (expected 0:$(10npd - 1))"))
    diamond = pid ÷ npd
    code = pid - diamond * npd
    ix = Int64(0)
    iy = Int64(0)
    shift = 0
    while code != 0
        ix |= (code & 1) << shift
        iy |= ((code >> 1) & 1) << shift
        code >>= 2
        shift += 1
    end
    return (Int(ix), Int(iy), Int(diamond))
end

# ---------------------------------------------------------------------------
# The chart's analytic inverse
#
# `xyd_to_point` is a composition of two invertible maps — a plain affine of the
# chart square's half-triangle onto a face's Snyder plane, then `snyder_inv_xyz`
# — and `ISEA` ships the inverse of the second (`snyder_fwd`, which also names
# the face). So the whole chart inverts in closed form, with no search and no
# point-in-polygon test, which is what makes `cellat` O(1) for this system.
#
# The one thing that has to be derived here is the FACE-TO-DIAMOND direction:
# `snyder_fwd` answers in `ISEA`'s twenty-face vocabulary and the chart is
# written in the ten-diamond one.
#
# Why the face `snyder_fwd` picks is the right half of the right diamond:
# `snyder_fwd` assigns `p` to the face whose centre it is closest to, and on a
# regular icosahedron the spherical Voronoi cell of a face centre IS that face's
# spherical triangle — the bisector between two adjacent face centres is exactly
# their shared edge, and non-adjacent centres are strictly farther. So the face
# determines the diamond and which of its two affine halves to invert, and the
# recovered `(x, y)` lands in that half of the unit square by construction.
# ---------------------------------------------------------------------------

"""
    FACE_TO_DIAMOND

For each of the twenty `ISEA` faces (indexed by `face + 1`), the pair
`(diamond, is_upper)` naming the diamond it belongs to and which of that
diamond's two halves it is. Derived from [`DIAMONDS`](@ref); the ten `(upper,
lower)` pairs use every face exactly once, which `_make_diamonds` asserts, so
this is a bijection `0:19 → (0:9) × Bool`.
"""
const FACE_TO_DIAMOND = let table = Vector{Tuple{Int,Bool}}(undef, 20)
    for d in 0:9
        table[DIAMONDS[d+1].upper+1] = (d, true)
        table[DIAMONDS[d+1].lower+1] = (d, false)
    end
    ntuple(i -> table[i], 20)
end

# Solve `u = x·a + y·b` for the real coefficients, `u`, `a`, `b` read as plane
# vectors: `imag(conj(a) * b)` is the determinant, which the diamond build
# asserts is exactly `2π/5` on every half of every diamond, so this never
# divides by anything near zero.
@inline function _affine_solve(u::ComplexF64, a::ComplexF64, b::ComplexF64)
    det = imag(conj(a) * b)
    return (-imag(conj(b) * u) / det, imag(conj(a) * u) / det)
end

"""
    point_to_xy(p) -> (x, y, diamond)

The exact inverse of [`xyd_to_point`](@ref): the continuous chart coordinates
`(x, y) ∈ [0, 1]²` and 0-based `diamond` of the unit-sphere point `p` (grid
frame).

`x` and `y` are clamped to `[0, 1]`. The clamp is round-off insurance, not a
projection of far-away points: `p` always lies in some face triangle, and the
affine inverse of that triangle's own map is in the square by construction —
but a point exactly on a diamond edge can come back as `-1e-17` or
`1 + 1e-17`, and every caller wants that to read as the edge.

`round_trip: xyd_to_point(point_to_xy(p)...) ≈ p` to `~1e-15` rad, and
`point_to_xy(xyd_to_point(x, y, d)) == (x, y, d)` to the same order (exactly on
the diamond corners, which are `ISEA.VERTICES` on both sides).
"""
function point_to_xy(p)
    f, w = ISEA.snyder_fwd((Float64(p[1]), Float64(p[2]), Float64(p[3])))
    d, is_upper = @inbounds FACE_TO_DIAMOND[f+1]
    dm = @inbounds DIAMONDS[d+1]
    x, y = is_upper ? _affine_solve(w - dm.cP0, dm.aP, dm.bP) :
           _affine_solve(w - dm.cQ0, dm.aQ, dm.bQ)
    return (clamp(x, 0.0, 1.0), clamp(y, 0.0, 1.0), d)
end

"""
    point_to_xyd(p, nside) -> (ix, iy, diamond)

The lattice cell of resolution `nside` containing the unit-sphere point `p`:
[`point_to_xy`](@ref) floored onto the `nside × nside` diamond lattice.

# Ties

A point on a shared cell boundary is contained by every cell meeting there, so
which one comes back is a tie, and it is resolved deterministically by the
arithmetic rather than by a rule applied afterwards:

  - **within a diamond**, `floor(x * nside)` puts the point on the HIGHER side
    of each chart cut line, so a point on the line `x = k/nside` belongs to the
    cell with `ix == k`;
  - **on a diamond edge**, the clamp in [`point_to_xy`](@ref) keeps `x` and `y`
    in `[0, 1]` and the `min(..., nside - 1)` below keeps the floor on the
    lattice, so the point lands in the diamond `snyder_fwd` named — i.e. the
    face whose centre it is closest to, lowest face index on an exact tie.

Both rules are decided by comparisons of `Float64`s that this port computes, so
the answer is deterministic per platform and self-consistent (the returned
cell's own centre maps back to it), which is what [`cellat`](@ref) requires.
"""
function point_to_xyd(p, nside::Integer)
    x, y, d = point_to_xy(p)
    n = Int(nside)
    ix = min(floor(Int, x * n), n - 1)
    iy = min(floor(Int, y * n), n - 1)
    return (ix, iy, d)
end

"""
    point_to_morton(p, nside) -> Int64

The 0-based Morton id `diamond * nside² + morton(ix, iy)` of the cell
containing `p`; [`point_to_xyd`](@ref) composed with [`xyd_to_morton`](@ref).
`nside` must be a power of two.
"""
point_to_morton(p, nside::Integer) =
    xyd_to_morton(point_to_xyd(p, nside)..., nside)
