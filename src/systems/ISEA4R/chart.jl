# Ten equal-area charts `[0,1]^2 -> S^2`, one per diamond. Coordinates and
# row-major/Morton identifiers are 0-based; data positions are identifier + 1.
# Morton order is canonical for power-of-two `nside`.

"""
    xyd_to_point(x, y, diamond) -> GO.UnitSphericalPoint

Evaluate `diamond`'s chart at `(x,y) in [0,1]^2` in the standard ISEA frame.
The two affine Snyder faces meet at `y == x`; the upper face owns the seam so
shared lattice points are bit-identical. For `nside <= 2^29` the branch test on
lattice coordinates agrees with the integer `iy >= ix`, and `fl(ix/n)` equals
`fl(2ix/2n)`, so the tiling is bit-exact across levels as well. The chart is
exactly equal-area: chart area `A` maps to solid angle `A*4π/10`.

Chart lines map to spherical curves, not great circles. Densify cell edges by
sampling this function when four corners are insufficient.
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

# Cell corners follow the chart's orientation-preserving counterclockwise order.

"""
    cell_corners(ix, iy, diamond, nside) -> NTuple{4, GO.UnitSphericalPoint}

Return four corners counterclockwise as seen from outside the sphere — a hard
contract, since the convex-clip kernel reduces a clockwise ring to empty rather
than erroring. Shared lattice corners are bit-identical under the chart's seam
rule. Snyder edges are curved; sample [`xyd_to_point`](@ref) for densification.
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

Return the chart midpoint of cell `(ix,iy)` on `diamond`. This canonical center
is not necessarily the spherical centroid of its four-corner polygon.
"""
cell_center(ix::Integer, iy::Integer, diamond::Integer, nside::Integer) =
    xyd_to_point((ix + 0.5) / nside, (iy + 0.5) / nside, diamond)

# ---------------------------------------------------------------------------
# Row-major order is `ix` fastest, then `iy`, then diamond. It supports any
# positive `nside`; canonical hierarchical identifiers use Morton order.
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
# Morton order interleaves `(ix,iy)` bits and requires power-of-two `nside`.
# ---------------------------------------------------------------------------

"""
    xyd_to_morton(ix, iy, diamond, nside) -> Int64

Return `diamond*nside^2 + morton(ix,iy)`, with `ix` bits in even positions and
`iy` bits in odd positions. `nside` must be a power of two. This checked entry
point validates coordinates and diamond; internal callers may use
[`xyd_to_morton_unchecked`](@ref).
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

Unchecked form of [`xyd_to_morton`](@ref). The caller must provide a
power-of-two `nside`, valid coordinates, and `diamond in 0:9`.
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
# The analytic inverse applies Snyder forward projection, maps the selected face
# to its diamond half, and solves the affine chart without spatial search.
# `snyder_fwd` names the nearest face centre, and on a regular icosahedron a face
# centre's spherical Voronoi cell IS that face's triangle, so the half it names
# is the one holding `p` and the affine solve lands inside the unit square.
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

Return continuous `(x,y)` coordinates and 0-based diamond for grid-frame point
`p`. Coordinates are clamped to `[0,1]` to absorb edge roundoff.
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

Return the lattice cell containing `p`. Within a diamond, boundary coordinates
belong to the higher-index cell; diamond-edge ties use the face selected by
Snyder projection. Results are deterministic per floating-point platform.
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
