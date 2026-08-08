module S2ChartTestSuite

# Tests for `src/S2/chart.jl`: the pure closed-form S2 cube-face charts and the
# row-major / Hilbert index maps over the `nside × nside` lattice.
#
# There is deliberately **no S2 oracle here**. `chart.jl` carries no s2geometry
# dependency and this repository vendors no s2geometry fixtures, so — unlike
# `test/HEALPix/test_chart.jl`, which has Healpix.jl as ground truth at every
# power-of-two `nside` — the ground truth below is analytic invariants plus
# internal consistency:
#
#   * the face frames are orthonormal and right-handed, and the evaluation
#     switch reproduces the frame sum exactly;
#   * `st_to_uv` hits its three anchors exactly, is monotone, odd about 1/2 and
#     C¹ across its branch seam, and `uv_to_st` inverts it;
#   * chart points are unit vectors, iso-`u` / iso-`v` sections are great-circle
#     arcs (a 3×3 determinant of three sampled points vanishes), corner rings
#     wind CCW seen from outside, and the twelve cube-edge seams agree
#     bit-identically face to face;
#   * the tessellation closes: distinct cell corners come out at `6 nside² + 2`,
#     the Euler count for `6 nside²` quads;
#   * both index maps are bijections, and the Hilbert one is local, nests across
#     resolutions, and reproduces hand-computed level-1/2 tables.
#
# Alignment with native `s2_cellid` values is intended by construction (the
# tables transcribe s2geometry's) but is NOT asserted anywhere here; that lands
# with the id-hierarchy milestone.

using Test
using Printf
import GeometryOps as GO

using DiscreteGlobalGrids.S2: FACE_U_AXIS, FACE_V_AXIS, FACE_NORMAL,
    face_uv_to_xyz, st_to_uv, uv_to_st, stf_to_point,
    cell_corners, cell_center,
    xyf_to_rowmajor, rowmajor_to_xyf,
    SWAP_MASK, INVERT_MASK, IJ_TO_POS, POS_TO_IJ, POS_TO_ORIENTATION,
    xyf_to_hilbert, hilbert_to_xyf

# Componentwise sup-norm between a chart point and a plain 3-tuple/vector.
maxdev(p, q) = max(abs(p[1] - q[1]), abs(p[2] - q[2]), abs(p[3] - q[3]))

cross3(a, b) = (a[2] * b[3] - a[3] * b[2],
                a[3] * b[1] - a[1] * b[3],
                a[1] * b[2] - a[2] * b[1])

# Scalar triple product; zero exactly when the three points are coplanar with
# the origin, i.e. lie on one great circle.
det3(a, b, c) = a[1] * (b[2] * c[3] - b[3] * c[2]) -
                a[2] * (b[1] * c[3] - b[3] * c[1]) +
                a[3] * (b[1] * c[2] - b[2] * c[1])

# Signed area of the 4-gon, seen from *outside* the sphere: `Σ pᵢ × pᵢ₊₁` is
# `2 * area * n̂` for a ring wound CCW about `n̂`, so dotting it with an
# outward radial direction is positive exactly when the ring is CCW.
function ccw_measure(corners)
    acc = (0.0, 0.0, 0.0)
    for i in 1:4
        acc = acc .+ cross3(corners[i], corners[i % 4 + 1])
    end
    outward = reduce((a, b) -> a .+ Tuple(b), corners; init=(0.0, 0.0, 0.0))
    return sum(acc .* outward)
end

# Recorded so the numbers land in the test log (and in the milestone report).
# `-Inf` rather than `0.0` as the neutral element: some of the quantities below
# are minima of positive measures, and a `0.0` floor would report a genuinely
# tiny one as zero.
const MEASURED = Dict{String,Float64}()
record!(key, value) = (MEASURED[key] = max(get(MEASURED, key, -Inf), value))

# Quantised key, for the multiplicity / Euler counts across face seams. Cell
# corners agree bit-identically *within* a face and across the six
# same-orientation seams, but a reversed seam pairs lattice fraction `k/n`
# against `1 - k/n`, which is a different Float64 from `(n-k)/n` at
# non-dyadic `n`. So the global counts are taken at 1e-10, well below any real
# separation and well above the ~1e-16 the seams actually differ by.
qkey(p) = (round(Int, p[1] * 1e10), round(Int, p[2] * 1e10), round(Int, p[3] * 1e10))

# ---------------------------------------------------------------------------
# 1. Face frames
# ---------------------------------------------------------------------------

@testset "face frames" begin
    @test length(FACE_U_AXIS) == 6 && length(FACE_V_AXIS) == 6 && length(FACE_NORMAL) == 6

    # Right-handed orthonormal on every face. This is what lets the CCW
    # argument for `cell_corners` be made on face 0 alone: the six frames are
    # in SO(3), under which the CCW-from-outside measure is invariant.
    for f in 0:5
        u, v, w = FACE_U_AXIS[f + 1], FACE_V_AXIS[f + 1], FACE_NORMAL[f + 1]
        @test sum(abs2, u) == 1 && sum(abs2, v) == 1 && sum(abs2, w) == 1
        @test sum(u .* v) == 0 && sum(v .* w) == 0 && sum(w .* u) == 0
        @test cross3(u, v) == w
    end

    # The six normals are the six signed coordinate axes, each exactly once —
    # s2geometry's numbering, with face `f`'s normal on axis `f % 3`, positive
    # for `f < 3`.
    @test Set(FACE_NORMAL) == Set(((1, 0, 0), (0, 1, 0), (0, 0, 1),
                                   (-1, 0, 0), (0, -1, 0), (0, 0, -1)))
    for f in 0:5
        @test FACE_NORMAL[f + 1][f % 3 + 1] == (f < 3 ? 1 : -1)
    end

    # The evaluation switch is the frame sum — exactly, not to a tolerance.
    # The switch exists so that shared-edge coordinates come out as literal
    # `±1.0` and plain IEEE negations; this pins that it costs no accuracy.
    for f in 0:5, a in 0:8, b in 0:8
        u = -1 + 2a / 8
        v = -1 + 2b / 8
        expected = FACE_NORMAL[f + 1] .+ u .* FACE_U_AXIS[f + 1] .+ v .* FACE_V_AXIS[f + 1]
        @test face_uv_to_xyz(f, u, v) == expected
    end
end

# ---------------------------------------------------------------------------
# 2. The quadratic ST <-> UV transform
# ---------------------------------------------------------------------------

@testset "st_to_uv / uv_to_st" begin
    # Anchors, exactly: these are what put the lattice endpoints on the cube
    # edges and the face centre on the axis.
    @test st_to_uv(0) == -1.0
    @test st_to_uv(0.5) == 0.0
    @test st_to_uv(1) == 1.0
    @test uv_to_st(-1) == 0.0
    @test uv_to_st(0) == 0.5
    @test uv_to_st(1) == 1.0

    # Strictly increasing (so the CCW orientation of the `(s, t)` emission
    # order survives into `(u, v)`), and inside `[-1, 1]`.
    prev = -Inf
    for k in 0:1000
        u = st_to_uv(k / 1000)
        @test u > prev
        @test -1.0 <= u <= 1.0
        prev = u
    end

    # Round trip, on a grid that deliberately includes non-dyadic fractions.
    worst_rt = -Inf
    for n in (64, 97), k in 0:n
        s = k / n
        worst_rt = max(worst_rt, abs(uv_to_st(st_to_uv(s)) - s))
    end
    record!("uv_to_st(st_to_uv(s)) - s", worst_rt)
    @test worst_rt <= 1e-15
    worst_rt2 = -Inf
    for k in 0:200
        u = -1 + 2k / 200
        worst_rt2 = max(worst_rt2, abs(st_to_uv(uv_to_st(u)) - u))
    end
    record!("st_to_uv(uv_to_st(u)) - u", worst_rt2)
    @test worst_rt2 <= 1e-15

    # Odd about s = 1/2, and *exactly* so in floating point: both branches
    # consume the same `1 - s`, and IEEE gives `fl(1 - w) == -fl(w - 1)` and
    # `fl(-x/3) == -fl(x/3)`. The seam tests below lean on this.
    for k in 0:64
        s = k / 64
        @test st_to_uv(1 - s) == -st_to_uv(s)
    end

    # C¹ across the branch seam at s = 1/2: both one-sided derivatives are 4/3.
    for h in (1e-6, 1e-7, 1e-8)
        @test (st_to_uv(0.5 + h) - st_to_uv(0.5)) / h ≈ 4 / 3 atol = 1e-5
        @test (st_to_uv(0.5) - st_to_uv(0.5 - h)) / h ≈ 4 / 3 atol = 1e-5
    end
end

# ---------------------------------------------------------------------------
# 3. stf_to_point
# ---------------------------------------------------------------------------

@testset "stf_to_point — norms, centers, cube corners" begin
    worst_norm = -Inf
    for face in 0:5, i in 0:20, j in 0:20
        p = stf_to_point(i / 20, j / 20, face)
        worst_norm = max(worst_norm, abs(sum(abs2, p) - 1))
    end
    record!("unit-norm residual (|p|^2 - 1)", worst_norm)
    @test worst_norm <= 1e-15

    # The face centre is the face normal, exactly. (`==` and not `===`: the
    # switch produces `-0.0` on the off-axes of faces 1-5, and `-0.0 == 0.0`
    # while `-0.0 === 0.0` is false. The point is the same point.)
    for face in 0:5
        n = FACE_NORMAL[face + 1]
        @test stf_to_point(0.5, 0.5, face) ==
            GO.UnitSphericalPoint(Float64(n[1]), Float64(n[2]), Float64(n[3]))
    end

    # The 24 face-square corner evaluations collapse to the 8 cube corners,
    # each shared by exactly 3 faces — and they collapse *bitwise*, which is
    # the strongest form of the "charts agree at the seams" claim.
    counts = Dict{NTuple{3,Float64},Int}()
    for face in 0:5, s in (0.0, 1.0), t in (0.0, 1.0)
        p = stf_to_point(s, t, face)
        key = (p[1], p[2], p[3])
        counts[key] = get(counts, key, 0) + 1
    end
    @test length(counts) == 8
    @test all(==(3), values(counts))
    # ... and they are the ±(1,1,1)/√3 cube corners.
    r = 1 / sqrt(3)
    for key in keys(counts)
        @test all(c -> abs(c) == r, key)
    end
    @test Set(keys(counts)) ==
        Set((a * r, b * r, c * r) for a in (-1.0, 1.0), b in (-1.0, 1.0), c in (-1.0, 1.0))

    # `cell_center` is the chart at the cell midpoint, by definition.
    for nside in (1, 3, 4), face in 0:5, ix in 0:(nside - 1), iy in 0:(nside - 1)
        @test cell_center(ix, iy, face, nside) ===
            stf_to_point((ix + 0.5) / nside, (iy + 0.5) / nside, face)
    end
end

# ---------------------------------------------------------------------------
# 4. Cell edges are great circles
#
# The property that distinguishes S2 from HEALPix at this layer: a chart line
# `u = const` is a central plane section, so its image is a great-circle arc and
# the four-corner ring IS the cell. Three points of a great circle are coplanar
# with the origin, so their scalar triple product vanishes.
# ---------------------------------------------------------------------------

@testset "iso-u and iso-v sections are great-circle arcs" begin
    worst = -Inf
    for face in 0:5, k in 0:6
        fixed = k / 6
        for m in 0:5
            a, b = m / 6, (m + 1) / 6
            mid = (a + b) / 2
            worst = max(worst, abs(det3(stf_to_point(fixed, a, face),
                                        stf_to_point(fixed, mid, face),
                                        stf_to_point(fixed, b, face))))
            worst = max(worst, abs(det3(stf_to_point(a, fixed, face),
                                        stf_to_point(mid, fixed, face),
                                        stf_to_point(b, fixed, face))))
        end
    end
    record!("|det[p(a), p(mid), p(b)]| on iso-u/iso-v sections", worst)
    @test worst <= 1e-15
end

# ---------------------------------------------------------------------------
# 5. Winding
# ---------------------------------------------------------------------------

@testset "cell_corners wind CCW from outside (nside = $nside)" for nside in (1, 2, 3, 5)
    worst = Inf
    for face in 0:5, ix in 0:(nside - 1), iy in 0:(nside - 1)
        worst = min(worst, ccw_measure(cell_corners(ix, iy, face, nside)))
    end
    # CCW as seen from outside the sphere — the convex-clip kernel clips a CW
    # ring to EMPTY, so this is a correctness contract, not cosmetics.
    @test worst > 0
    record!("min CCW measure (nside=$nside)", worst)
end

# ---------------------------------------------------------------------------
# 6. Seams
#
# The twelve cube edges, as pairs of chart evaluations that must land on the
# same physical point. Six are "same-orientation" (the shared parameter runs
# the same way on both faces) and six are "reversed" (`t ↔ 1 - t`, with the
# roles of the two axes swapped).
#
# All twelve are asserted **bit-identical** (`==`), at every `nside`, dyadic or
# not. That is stronger than one might expect of the reversed rows, and it is
# not luck:
#
#   * the switch in `face_uv_to_xyz` places literal `±1.0` and exact IEEE
#     negations on the shared coordinate, so no arithmetic reorders across a
#     seam;
#   * `st_to_uv(1 - s) == -st_to_uv(s)` exactly for *any* Float64 `s ∈ [0, 1]`
#     — for `s >= 1/2`, `1 - s` is exact by Sterbenz and `1 - (1 - s) == s`;
#     for `s < 1/2` both sides consume the identical `fl(1 - s)`, and IEEE
#     gives `fl(1 - w) == -fl(w - 1)` and `fl(-x/3) == -fl(x/3)`;
#   * the normalisation divides by the same `n = sqrt(...)` on both sides,
#     the summands being a permutation-with-signs of each other.
#
# `==` and never `===`: the seam midpoints produce `±0.0` on the off-axes, and
# `-0.0 === 0.0` is false while the points are the same point.
# ---------------------------------------------------------------------------

# (label, left evaluation, right evaluation) as functions of the shared
# lattice fraction.
const SAME_SEAMS = (
    ("face 0 s=1 | face 1 s=0", a -> stf_to_point(1.0, a, 0), a -> stf_to_point(0.0, a, 1)),
    ("face 0 t=0 | face 5 t=1", a -> stf_to_point(a, 0.0, 0), a -> stf_to_point(a, 1.0, 5)),
    ("face 1 t=1 | face 2 t=0", a -> stf_to_point(a, 1.0, 1), a -> stf_to_point(a, 0.0, 2)),
    ("face 2 s=1 | face 3 s=0", a -> stf_to_point(1.0, a, 2), a -> stf_to_point(0.0, a, 3)),
    ("face 3 t=1 | face 4 t=0", a -> stf_to_point(a, 1.0, 3), a -> stf_to_point(a, 0.0, 4)),
    ("face 4 s=1 | face 5 s=0", a -> stf_to_point(1.0, a, 4), a -> stf_to_point(0.0, a, 5)),
)

const REVERSED_SEAMS = (
    ("face 0 s=0 | face 4 t=1", a -> stf_to_point(0.0, a, 0), a -> stf_to_point(1 - a, 1.0, 4)),
    ("face 0 t=1 | face 2 s=0", a -> stf_to_point(a, 1.0, 0), a -> stf_to_point(0.0, 1 - a, 2)),
    ("face 1 s=1 | face 3 t=0", a -> stf_to_point(1.0, a, 1), a -> stf_to_point(1 - a, 0.0, 3)),
    ("face 1 t=0 | face 5 s=1", a -> stf_to_point(a, 0.0, 1), a -> stf_to_point(1.0, 1 - a, 5)),
    ("face 2 t=1 | face 4 s=0", a -> stf_to_point(a, 1.0, 2), a -> stf_to_point(0.0, 1 - a, 4)),
    ("face 3 s=1 | face 5 t=0", a -> stf_to_point(1.0, a, 3), a -> stf_to_point(1 - a, 0.0, 5)),
)

@testset "cube-edge seams are bit-identical (nside = $nside)" for nside in (2, 3, 4, 5, 7, 8)
    for (kind, rows) in (("same", SAME_SEAMS), ("reversed", REVERSED_SEAMS))
        for (label, left, right) in rows
            worst = -Inf
            for k in 0:nside
                a = k / nside
                p, q = left(a), right(a)
                @test p == q
                worst = max(worst, maxdev(p, q))
            end
            record!("seam sup-deviation, $kind", worst)
        end
    end
end

@testset "seam multiplicity and Euler count (nside = $nside)" for nside in (1, 2, 3, 4, 5, 7)
    # Every boundary-lattice point of every face, each counted once per face
    # per point (`4 nside` per face, `24 nside` total). Quantised, because a
    # reversed seam pairs `1 - k/n` against `(n-k)/n`, which differ in the last
    # bits at non-dyadic `n`.
    mult = Dict{NTuple{3,Int},Int}()
    total = 0
    for face in 0:5
        perimeter = NTuple{2,Float64}[]
        for k in 0:(nside - 1); push!(perimeter, (k / nside, 0.0)); end
        for k in 0:(nside - 1); push!(perimeter, (1.0, k / nside)); end
        for k in nside:-1:1;    push!(perimeter, (k / nside, 1.0)); end
        for k in nside:-1:1;    push!(perimeter, (0.0, k / nside)); end
        @test length(perimeter) == 4nside
        for (s, t) in perimeter
            total += 1
            key = qkey(stf_to_point(s, t, face))
            mult[key] = get(mult, key, 0) + 1
        end
    end
    @test total == 24nside
    # `12(nside - 1)` edge-interior points shared by 2 faces, plus the 8 cube
    # corners shared by 3 — which is `24 nside` evaluations over `12 nside - 4`
    # distinct points.
    @test length(mult) == 12nside - 4
    @test count(==(2), values(mult)) == 12 * (nside - 1)
    @test count(==(3), values(mult)) == 8

    # And the whole tessellation closes: `V - E + F = 2` with `F = 6 nside²`
    # quads and `E = 2F` gives `V = 6 nside² + 2`. Hitting that count is a
    # global statement — the charts agree across every face seam and no cell is
    # missing or doubled.
    corners = Set{NTuple{3,Int}}()
    for face in 0:5, ix in 0:(nside - 1), iy in 0:(nside - 1)
        for p in cell_corners(ix, iy, face, nside)
            push!(corners, qkey(p))
        end
    end
    @test length(corners) == 6nside^2 + 2

    # Within a face, shared lattice corners are bit-identical (not merely
    # quantised-equal), so the tessellation is exact there.
    for face in 0:5, ix in 0:(nside - 2), iy in 0:(nside - 2)
        c = cell_corners(ix, iy, face, nside)
        east = cell_corners(ix + 1, iy, face, nside)
        north = cell_corners(ix, iy + 1, face, nside)
        @test c[1] === east[2] && c[4] === east[3]
        @test c[1] === north[4] && c[2] === north[3]
    end
end

# ---------------------------------------------------------------------------
# 7. Row-major maps
# ---------------------------------------------------------------------------

@testset "row-major maps (nside = $nside)" for nside in (1, 2, 3, 4, 5, 7)
    ncell = 6nside^2

    seen = Set{NTuple{3,Int}}()
    for q in 0:(ncell - 1)
        ix, iy, face = rowmajor_to_xyf(q, nside)
        @test 0 <= ix < nside && 0 <= iy < nside && 0 <= face <= 5
        push!(seen, (ix, iy, face))
        @test xyf_to_rowmajor(ix, iy, face, nside) == q
    end
    @test length(seen) == ncell                       # onto the full lattice

    for face in 0:5, ix in 0:(nside - 1), iy in 0:(nside - 1)
        q = xyf_to_rowmajor(ix, iy, face, nside)
        @test 0 <= q < ncell
        @test rowmajor_to_xyf(q, nside) == (ix, iy, face)
        # The formula itself, spelled out.
        @test q == face * nside^2 + iy * nside + ix
    end

    # Face `f` owns exactly `[f nside², (f+1) nside²)`.
    for face in 0:5
        @test xyf_to_rowmajor(0, 0, face, nside) == face * nside^2
        @test xyf_to_rowmajor(nside - 1, nside - 1, face, nside) == (face + 1) * nside^2 - 1
    end

    # The inverse is the checked direction (ids arrive from users); the forward
    # one is GIGO, like `xyf_to_ring`.
    @test_throws ArgumentError rowmajor_to_xyf(-1, nside)
    @test_throws ArgumentError rowmajor_to_xyf(ncell, nside)
end

# ---------------------------------------------------------------------------
# 8. Hilbert maps
# ---------------------------------------------------------------------------

@testset "Hilbert maps reject what they cannot index" begin
    for nside in (3, 5, 6, 7)
        @test_throws ArgumentError xyf_to_hilbert(0, 0, 0, nside)
        @test_throws ArgumentError hilbert_to_xyf(0, nside)
    end
    @test_throws ArgumentError xyf_to_hilbert(4, 0, 0, 4)
    @test_throws ArgumentError xyf_to_hilbert(0, -1, 0, 4)
    @test_throws ArgumentError xyf_to_hilbert(0, 0, 6, 4)
    @test_throws ArgumentError xyf_to_hilbert(0, 0, -1, 4)
    @test_throws ArgumentError hilbert_to_xyf(-1, 4)
    @test_throws ArgumentError hilbert_to_xyf(6 * 16, 4)
end

@testset "Hilbert tables" begin
    # The transcribed s2geometry state machine, pinned as data. `IJ_TO_POS` and
    # `POS_TO_IJ` are per-orientation inverse permutations of `0:3`.
    @test SWAP_MASK == 1 && INVERT_MASK == 2
    @test length(IJ_TO_POS) == 4 && length(POS_TO_IJ) == 4
    for o in 0:3
        @test sort(collect(IJ_TO_POS[o + 1])) == [0, 1, 2, 3]
        for ij in 0:3
            @test POS_TO_IJ[o + 1][IJ_TO_POS[o + 1][ij + 1] + 1] == ij
        end
    end
    @test POS_TO_ORIENTATION == (SWAP_MASK, 0, 0, SWAP_MASK | INVERT_MASK)

    # Hand-computed level-1 and level-2 openings of the curve. An even face
    # starts unswapped, so the first four cells walk `+t` then `+s` then `-t`;
    # an odd face starts at `SWAP_MASK`, which transposes that.
    @test [hilbert_to_xyf(p, 2)[1:2] for p in 0:3] ==
        [(0, 0), (0, 1), (1, 1), (1, 0)]                     # face 0, even
    @test [hilbert_to_xyf(1 * 4 + p, 2)[1:2] for p in 0:3] ==
        [(0, 0), (1, 0), (1, 1), (0, 1)]                     # face 1, odd
    # At nside = 4 the first four positions are the level-2 opening *inside*
    # the first quadrant, whose orientation has already been advanced by
    # `POS_TO_ORIENTATION[1] == SWAP_MASK`.
    @test [hilbert_to_xyf(p, 4)[1:2] for p in 0:3] ==
        [(0, 0), (1, 0), (1, 1), (0, 1)]                     # face 0, even
    # And the forward direction agrees with all of it.
    for (n, f, table) in ((2, 0, [(0, 0), (0, 1), (1, 1), (1, 0)]),
                          (2, 1, [(0, 0), (1, 0), (1, 1), (0, 1)]),
                          (4, 0, [(0, 0), (1, 0), (1, 1), (0, 1)]))
        for (p, (ix, iy)) in enumerate(table)
            @test xyf_to_hilbert(ix, iy, f, n) == f * n^2 + (p - 1)
        end
    end
end

@testset "Hilbert bijectivity and block structure (nside = $nside)" for nside in (1, 2, 4, 8, 16)
    ncell = 6nside^2

    seen = Set{NTuple{3,Int}}()
    for h in 0:(ncell - 1)
        ix, iy, face = hilbert_to_xyf(h, nside)
        @test 0 <= ix < nside && 0 <= iy < nside && 0 <= face <= 5
        push!(seen, (ix, iy, face))
        @test xyf_to_hilbert(ix, iy, face, nside) == h
    end
    @test length(seen) == ncell                       # onto the full lattice

    for face in 0:5, ix in 0:(nside - 1), iy in 0:(nside - 1)
        h = xyf_to_hilbert(ix, iy, face, nside)
        # Face `f` occupies exactly `[f nside², (f+1) nside²)` — the block
        # structure the scaffold ordinal `face * 4^level + position` names.
        @test face * nside^2 <= h < (face + 1) * nside^2
        @test hilbert_to_xyf(h, nside) == (ix, iy, face)
    end
end

@testset "Hilbert locality and prefix nesting" begin
    # Locality: consecutive within-face positions are lattice-adjacent. This is
    # the defining property of a Hilbert curve, and the one a transcription
    # error in the tables would break immediately.
    for nside in (2, 4, 8), face in 0:5
        for p in 0:(nside^2 - 2)
            a = hilbert_to_xyf(face * nside^2 + p, nside)
            b = hilbert_to_xyf(face * nside^2 + p + 1, nside)
            @test a[3] == face && b[3] == face
            @test abs(a[1] - b[1]) + abs(a[2] - b[2]) == 1
        end
    end

    # Prefix nesting: dropping the low two bits of the within-face position
    # halves both lattice coordinates, i.e. steps exactly one level up. This is
    # what makes cross-resolution refinement contiguous (the face-grid suite
    # observes it as rows `4j-3:4j` of a coarse column).
    for nside in (1, 2, 4, 8), face in 0:5
        for ix in 0:(2nside - 1), iy in 0:(2nside - 1)
            fine = xyf_to_hilbert(ix, iy, face, 2nside) - face * (2nside)^2
            coarse = xyf_to_hilbert(ix >> 1, iy >> 1, face, nside) - face * nside^2
            @test (fine >> 2) == coarse
        end
    end
end

@printf("[S2 chart] measured:\n")
for key in sort!(collect(keys(MEASURED)))
    @printf("[S2 chart]   %-52s %+.3e\n", key, MEASURED[key])
end

end # module S2ChartTestSuite
