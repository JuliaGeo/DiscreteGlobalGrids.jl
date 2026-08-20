# S2 system tests:
#
#   1. ORACLE. There is no external S2 implementation or fixture in this
#      repository:
#      `src/systems/S2/` carries no s2geometry dependency and no s2geometry
#      fixtures are vendored, so unlike HEALPix (which has Healpix.jl) the
#      ground truth here is analytic invariants, internal consistency, and
#      hand-computed tables: the face frames, the
#      quadratic transform's anchors and exact oddness, great-circle sections,
#      CCW winding, the twelve cube-edge seams bit-identical, the Euler count,
#      and both index codecs' bijectivity, locality and prefix nesting.
#
#   2. CONTRACT. The two conformance suites from
#      `DiscreteGlobalGridsConformanceTesting`, with default kwargs.
#
#   3. STRUCTURAL. The hierarchy over the scaffold ordinal, the cube-edge seam topology checked
#      GEOMETRICALLY against shared boundary points, the rotational neighbour
#      order and its documented start, the exact four-corner `node_extent`, and
#      the chart inverse behind `cellat`.
#
# `test/runtests.jl` includes this file; it also runs standalone:
#     julia --project=test --startup-file=no test/systems/S2/runtests.jl

module S2SystemTests

using Test
using Printf
using Random

using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
const S2 = DiscreteGlobalGrids.S2

using DiscreteGlobalGridsConformanceTesting

import GeometryOps as GO
const US = GO.UnitSpherical

using DiscreteGlobalGrids.S2: FACE_U_AXIS, FACE_V_AXIS, FACE_NORMAL,
    face_uv_to_xyz, xyz_to_face, st_to_uv, uv_to_st, stf_to_point, point_to_xyf,
    cell_corners, cell_center,
    xyf_to_rowmajor, rowmajor_to_xyf,
    SWAP_MASK, INVERT_MASK, IJ_TO_POS, POS_TO_IJ, POS_TO_ORIENTATION,
    xyf_to_hilbert, hilbert_to_xyf,
    SEAM, NEIGHBOR_OFFSETS, wrap_xyf,
    S2System, MAX_LEVEL

const SYS = S2System()

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

cross3(a, b) = (a[2] * b[3] - a[3] * b[2],
                a[3] * b[1] - a[1] * b[3],
                a[1] * b[2] - a[2] * b[1])

# Scalar triple product; zero exactly when the three points are coplanar with
# the origin, i.e. lie on one great circle.
det3(a, b, c) = a[1] * (b[2] * c[3] - b[3] * c[2]) -
                a[2] * (b[1] * c[3] - b[3] * c[1]) +
                a[3] * (b[1] * c[2] - b[2] * c[1])

# Signed area of the ring, seen from *outside* the sphere: `Σ pᵢ × pᵢ₊₁` is
# `2 * area * n̂` for a ring wound CCW about `n̂`, so dotting it with an outward
# radial direction is positive exactly when the ring is CCW.
function ccw_measure(corners)
    acc = (0.0, 0.0, 0.0)
    n = length(corners)
    for i in 1:n
        acc = acc .+ cross3(corners[i], corners[i % n + 1])
    end
    outward = reduce((a, b) -> a .+ Tuple(b), corners; init = (0.0, 0.0, 0.0))
    return sum(acc .* outward)
end

"""
Signed-zero-normalised coordinate triple, for identifying a shared boundary
vertex through a `Set` or a `Dict`.

Cell corners agree BIT-IDENTICALLY across a cube seam except for the sign of
zero: `face_uv_to_xyz` puts exact IEEE negations on the shared coordinate, and
`-(0.0) === -0.0`. The two points compare `==` and are the same point on the
sphere, but `Set`/`Dict` use `isequal`, under which `-0.0` and `0.0` differ. So
every set-based incidence test in this file goes through here — the alternative,
a tolerance, would weaken a check that is otherwise exact. Adding `0.0` maps
`-0.0` to `0.0` and leaves every other Float64 alone.
"""
vkey(p) = (p[1] + 0.0, p[2] + 0.0, p[3] + 0.0)

"""
A deterministic cell sample: ALL `6 * 4^level` cells when the level fits under
`n`, a seeded draw of at most `n` of them when it does not. The draw is
`unique`d before sorting, so it can return fewer than `n` ids.
"""
function cell_sample(level::Integer, n::Integer = 64)
    ncell = 6 * Int64(4)^level
    ncell <= n && return collect(Int64, 0:(ncell - 1))
    rng = MersenneTwister(20_260_813 + level)
    return sort!(unique(rand(rng, Int64(0):(ncell - 1), n)))
end

# Retain the largest measurement for each key in the test log.
const MEASURED = Dict{String,Float64}()
record!(key, value) = (MEASURED[key] = max(get(MEASURED, key, -Inf), value))

@testset "S2 system (T11)" begin

# =========================================================================
# 1. Chart: frames, transform, geometry, seams
#
# Analytic invariants and internal consistency provide the oracle because no
# external reference implementation is available here.
# =========================================================================

@testset "face frames" begin
    @test length(FACE_U_AXIS) == 6 && length(FACE_V_AXIS) == 6 && length(FACE_NORMAL) == 6

    # Right-handed orthonormal on every face. This is what lets the CCW
    # argument for `cell_corners` be made on face 0 alone: the six frames are in
    # SO(3), under which the CCW-from-outside measure is invariant.
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

    # The evaluation switch is the frame sum — exactly, not to a tolerance. The
    # switch exists so that shared-edge coordinates come out as literal `±1.0`
    # and plain IEEE negations; this pins that it costs no accuracy.
    for f in 0:5, a in 0:8, b in 0:8
        u = -1 + 2a / 8
        v = -1 + 2b / 8
        expected = FACE_NORMAL[f + 1] .+ u .* FACE_U_AXIS[f + 1] .+ v .* FACE_V_AXIS[f + 1]
        @test face_uv_to_xyz(f, u, v) == expected
    end
end

@testset "st_to_uv / uv_to_st" begin
    # Anchors, exactly: these are what put the lattice endpoints on the cube
    # edges and the face centre on the axis.
    @test st_to_uv(0) == -1.0
    @test st_to_uv(0.5) == 0.0
    @test st_to_uv(1) == 1.0
    @test uv_to_st(-1) == 0.0
    @test uv_to_st(0) == 0.5
    @test uv_to_st(1) == 1.0

    # Strictly increasing (so the CCW orientation of the `(s, t)` emission order
    # survives into `(u, v)`), and inside `[-1, 1]`.
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

    # Odd about s = 1/2, and *exactly* so in floating point. This is not a
    # curiosity: `neighbors.jl` derives the cube-seam correspondence in `(u, v)`
    # and applies it to the `(s, t)` lattice, which is legal only because a sign
    # flip in `u` is exactly the reflection `s ↦ 1 - s`.
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

@testset "stf_to_point — norms, centres, cube corners" begin
    worst_norm = -Inf
    for face in 0:5, i in 0:20, j in 0:20
        p = stf_to_point(i / 20, j / 20, face)
        worst_norm = max(worst_norm, abs(sum(abs2, p) - 1))
    end
    record!("unit-norm residual (|p|^2 - 1)", worst_norm)
    @test worst_norm <= 1e-15

    # The face centre is the face normal, exactly. (`==` and not `===`: the
    # switch produces `-0.0` on the off-axes of faces 1-5.)
    for face in 0:5
        n = FACE_NORMAL[face + 1]
        @test stf_to_point(0.5, 0.5, face) ==
            GO.UnitSphericalPoint(Float64(n[1]), Float64(n[2]), Float64(n[3]))
    end

    # The 24 face-square corner evaluations collapse to the 8 cube corners, each
    # shared by exactly 3 faces — and they collapse *bitwise* (modulo the sign of
    # zero, which `vkey` normalises), the strongest form of "the charts agree at
    # the seams".
    counts = Dict{NTuple{3,Float64},Int}()
    for face in 0:5, s in (0.0, 1.0), t in (0.0, 1.0)
        k = vkey(stf_to_point(s, t, face))
        counts[k] = get(counts, k, 0) + 1
    end
    @test length(counts) == 8
    @test all(==(3), values(counts))
    r = 1 / sqrt(3)
    @test Set(keys(counts)) ==
        Set((a * r, b * r, c * r) for a in (-1.0, 1.0), b in (-1.0, 1.0), c in (-1.0, 1.0))

    # `cell_center` is the chart at the cell midpoint, by definition.
    for nside in (1, 3, 4), face in 0:5, ix in 0:(nside - 1), iy in 0:(nside - 1)
        @test cell_center(ix, iy, face, nside) ===
            stf_to_point((ix + 0.5) / nside, (iy + 0.5) / nside, face)
    end
end

@testset "iso-u and iso-v sections are great-circle arcs" begin
    # The property that distinguishes S2 from HEALPix at this layer, and the
    # reason `cell_boundary` needs no densification: a chart line `u = const` is
    # a central plane section, so its image is a great-circle arc.
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

@testset "cell_corners wind CCW from outside (nside = $nside)" for nside in (1, 2, 3, 5)
    worst = Inf
    for face in 0:5, ix in 0:(nside - 1), iy in 0:(nside - 1)
        worst = min(worst, ccw_measure(cell_corners(ix, iy, face, nside)))
    end
    # A clockwise ring clips to EMPTY in the convex-clip kernel, silently, so
    # this is a correctness contract rather than cosmetics.
    @test worst > 0
    record!("min CCW measure, negated (< 0 ⇒ all CCW, nside=$nside)", -worst)
end

# The twelve cube edges, as pairs of chart evaluations that must land on the
# same physical point. Six are "same-orientation" (the shared parameter runs the
# same way on both faces) and six are "reversed" (`t ↔ 1 - t`, with the roles of
# the two axes swapped).
SAME_SEAMS = (
    ("face 0 s=1 | face 1 s=0", a -> stf_to_point(1.0, a, 0), a -> stf_to_point(0.0, a, 1)),
    ("face 0 t=0 | face 5 t=1", a -> stf_to_point(a, 0.0, 0), a -> stf_to_point(a, 1.0, 5)),
    ("face 1 t=1 | face 2 t=0", a -> stf_to_point(a, 1.0, 1), a -> stf_to_point(a, 0.0, 2)),
    ("face 2 s=1 | face 3 s=0", a -> stf_to_point(1.0, a, 2), a -> stf_to_point(0.0, a, 3)),
    ("face 3 t=1 | face 4 t=0", a -> stf_to_point(a, 1.0, 3), a -> stf_to_point(a, 0.0, 4)),
    ("face 4 s=1 | face 5 s=0", a -> stf_to_point(1.0, a, 4), a -> stf_to_point(0.0, a, 5)),
)

REVERSED_SEAMS = (
    ("face 0 s=0 | face 4 t=1", a -> stf_to_point(0.0, a, 0), a -> stf_to_point(1 - a, 1.0, 4)),
    ("face 0 t=1 | face 2 s=0", a -> stf_to_point(a, 1.0, 0), a -> stf_to_point(0.0, 1 - a, 2)),
    ("face 1 s=1 | face 3 t=0", a -> stf_to_point(1.0, a, 1), a -> stf_to_point(1 - a, 0.0, 3)),
    ("face 1 t=0 | face 5 s=1", a -> stf_to_point(a, 0.0, 1), a -> stf_to_point(1.0, 1 - a, 5)),
    ("face 2 t=1 | face 4 s=0", a -> stf_to_point(a, 1.0, 2), a -> stf_to_point(0.0, 1 - a, 4)),
    ("face 3 s=1 | face 5 t=0", a -> stf_to_point(1.0, a, 3), a -> stf_to_point(1 - a, 0.0, 5)),
)

@testset "cube-edge seams are bit-identical (nside = $nside)" for nside in (2, 3, 4, 5, 7, 8)
    # All twelve are `==` at every `nside`, dyadic or not: the switch in
    # `face_uv_to_xyz` puts literal `±1.0` and exact IEEE negations on the shared
    # coordinate, `st_to_uv` is exactly odd, and both sides divide by the same
    # normalisation. `==` and never `===`: seam midpoints produce `±0.0` on the
    # off-axes and are the same point regardless.
    for rows in (SAME_SEAMS, REVERSED_SEAMS), (label, left, right) in rows
        for k in 0:nside
            a = k / nside
            @test left(a) == right(a)
        end
    end
end

@testset "seam multiplicity and Euler count (nside = $nside)" for nside in (1, 2, 4, 8)
    # Every boundary-lattice point of every face, counted once per face. At
    # power-of-two `nside` the reversed seams pair `1 - k/n` against `(n-k)/n`,
    # which are the same Float64. Exact keys therefore suffice because the id
    # space contains only dyadic `nside` values.
    mult = Dict{NTuple{3,Float64},Int}()
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
            k = vkey(stf_to_point(s, t, face))
            mult[k] = get(mult, k, 0) + 1
        end
    end
    @test total == 24nside
    @test length(mult) == 12nside - 4
    @test count(==(2), values(mult)) == 12 * (nside - 1)
    @test count(==(3), values(mult)) == 8

    # And the whole tessellation closes: `V - E + F = 2` with `F = 6 nside²`
    # quads and `E = 2F` gives `V = 6 nside² + 2`. A global statement — the
    # charts agree across every seam and no cell is missing or doubled.
    corners = Set{NTuple{3,Float64}}()
    for face in 0:5, ix in 0:(nside - 1), iy in 0:(nside - 1)
        for p in cell_corners(ix, iy, face, nside)
            push!(corners, vkey(p))
        end
    end
    @test length(corners) == 6nside^2 + 2

    # Within a face, shared lattice corners are bit-identical (`===`, signed
    # zero included), so the tessellation is exact there with nothing normalised.
    for face in 0:5, ix in 0:(nside - 2), iy in 0:(nside - 2)
        c = cell_corners(ix, iy, face, nside)
        east = cell_corners(ix + 1, iy, face, nside)
        north = cell_corners(ix, iy + 1, face, nside)
        @test c[1] === east[2] && c[4] === east[3]
        @test c[1] === north[4] && c[2] === north[3]
    end
end

# =========================================================================
# 2. Codecs: row-major, Hilbert, and the chart inverse
# =========================================================================

@testset "row-major maps (nside = $nside)" for nside in (1, 2, 3, 4, 5, 7)
    ncell = 6nside^2
    seen = Set{NTuple{3,Int}}()
    for q in 0:(ncell - 1)
        ix, iy, face = rowmajor_to_xyf(q, nside)
        @test 0 <= ix < nside && 0 <= iy < nside && 0 <= face <= 5
        push!(seen, (ix, iy, face))
        @test xyf_to_rowmajor(ix, iy, face, nside) == q
    end
    @test length(seen) == ncell
    for face in 0:5, ix in 0:(nside - 1), iy in 0:(nside - 1)
        @test xyf_to_rowmajor(ix, iy, face, nside) == face * nside^2 + iy * nside + ix
    end
    # The inverse is the checked direction (ids arrive from users); the forward
    # one is garbage-in-garbage-out by documented design.
    @test_throws ArgumentError rowmajor_to_xyf(-1, nside)
    @test_throws ArgumentError rowmajor_to_xyf(ncell, nside)
end

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
    # The transcribed s2geometry state machine, pinned as data.
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
    @test [hilbert_to_xyf(p, 4)[1:2] for p in 0:3] ==
        [(0, 0), (1, 0), (1, 1), (0, 1)]                     # face 0, even
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
    @test length(seen) == ncell
    for face in 0:5, ix in 0:(nside - 1), iy in 0:(nside - 1)
        h = xyf_to_hilbert(ix, iy, face, nside)
        # Face `f` occupies exactly `[f nside², (f+1) nside²)` — the block
        # structure the scaffold ordinal `face * 4^level + position` names.
        @test face * nside^2 <= h < (face + 1) * nside^2
        @test hilbert_to_xyf(h, nside) == (ix, iy, face)
    end
end

@testset "the Hilbert order is LOCAL and NESTS (the trait's whole basis)" begin
    # Locality: consecutive within-face positions are lattice-adjacent. The
    # defining property of a Hilbert curve, and the one a transcription error in
    # the tables would break immediately.
    for nside in (2, 4, 8), face in 0:5
        for p in 0:(nside^2 - 2)
            a = hilbert_to_xyf(face * nside^2 + p, nside)
            b = hilbert_to_xyf(face * nside^2 + p + 1, nside)
            @test a[3] == face && b[3] == face
            @test abs(a[1] - b[1]) + abs(a[2] - b[2]) == 1
        end
    end

    # NESTING — verified by EXHAUSTION over levels 0-6, because it is the single
    # claim `has_sorted_subtrees(S2System()) == true` and the whole `÷ 4` /
    # `4p + k` hierarchy rest on, and the Hilbert tables it would otherwise
    # follow from are transcribed constants nothing here can check.
    #
    # Three statements, each stronger than the last:
    #   (a) dropping the low two bits of the within-face POSITION halves both
    #       lattice coordinates, i.e. steps exactly one level up;
    #   (b) the same holds of the whole scaffold ORDINAL, because the face term
    #       `face * 4^level` divides through by 4 untouched;
    #   (c) so the four children of ordinal `p` are exactly `4p:4p+3`.
    for level in 0:6
        n = 1 << level
        n2 = 2n
        for face in 0:5, ix in 0:(n2 - 1), iy in 0:(n2 - 1)
            fine = xyf_to_hilbert(ix, iy, face, n2)
            coarse = xyf_to_hilbert(ix >> 1, iy >> 1, face, n)
            @test (fine - face * n2^2) >> 2 == coarse - face * n^2      # (a)
            @test fine ÷ 4 == coarse                                    # (b)
        end
        for p in 0:(6 * 4^level - 1)                                    # (c)
            ix, iy, face = hilbert_to_xyf(p, n)
            kids = sort([xyf_to_hilbert(2ix + a, 2iy + b, face, n2) for a in 0:1 for b in 0:1])
            @test kids == collect(4p:(4p + 3))
        end
    end
end

@testset "the chart inverse: point_to_xyf" begin
    # Round trip from cell centres — exhaustive at levels 0-3, sampled at 6, 12
    # and 30 (the deepest level, where the codec is 60 bits wide).
    for level in (0, 1, 2, 3, 6, 12, MAX_LEVEL)
        nside = Int64(1) << level
        for h in cell_sample(level, 96)
            ix, iy, face = hilbert_to_xyf(h, nside)
            @test point_to_xyf(cell_center(ix, iy, face, nside), nside) == (ix, iy, face)
        end
    end

    # Interior points that are not the centre: the chart at a 5x5 sub-lattice of
    # the cell, offset off the border so no sample lands on a cut line.
    for level in (0, 2, 4)
        nside = Int64(1) << level
        for h in cell_sample(level, 48)
            ix, iy, face = hilbert_to_xyf(h, nside)
            for a in 1:5, b in 1:5
                p = stf_to_point((ix + a / 6) / nside, (iy + b / 6) / nside, face)
                @test point_to_xyf(p, nside) == (ix, iy, face)
            end
        end
    end

    # The face rule, and its documented tie-break: the largest-magnitude
    # component, ties toward the lower axis, a zero component counting positive.
    @test xyz_to_face(GO.UnitSphericalPoint(1.0, 0.0, 0.0)) == 0
    @test xyz_to_face(GO.UnitSphericalPoint(0.0, 1.0, 0.0)) == 1
    @test xyz_to_face(GO.UnitSphericalPoint(0.0, 0.0, 1.0)) == 2
    @test xyz_to_face(GO.UnitSphericalPoint(-1.0, 0.0, 0.0)) == 3
    @test xyz_to_face(GO.UnitSphericalPoint(0.0, -1.0, 0.0)) == 4
    @test xyz_to_face(GO.UnitSphericalPoint(0.0, 0.0, -1.0)) == 5
    r = 1 / sqrt(3)
    @test xyz_to_face(GO.UnitSphericalPoint(r, r, r)) == 0        # three-way tie -> x
    @test xyz_to_face(GO.UnitSphericalPoint(-r, r, r)) == 3       # ... signed
    s = 1 / sqrt(2)
    @test xyz_to_face(GO.UnitSphericalPoint(0.0, s, s)) == 1      # y/z tie -> y
end

# =========================================================================
# 3. Ids and hierarchy over the scaffold ordinal
#
# Hierarchy and adjacency over the scaffold ordinal.
# =========================================================================

@testset "system traits and level grids" begin
    @test cellindextype(SYS) === LevelIndex
    @test cellindextypes(SYS) == (LevelIndex,)          # no alternate scheme yet
    @test levels(SYS) == 0:30
    @test maxlevel(SYS) == 30
    @test has_sorted_subtrees(SYS)
    @test maxneighbors(SYS, Vertex()) == 8
    @test maxneighbors(SYS, Edge()) == 4
    @test maxneighbors(SYS) == 8

    for l in (0, 1, 5, 30)
        g = levelgrid(SYS, l)
        @test g isa DGG.HierarchicalLevelGrid{S2System}
        @test system(g) === SYS
        @test level(g) == l
        @test ncells(g) == 6 * 4^l
    end
    # Level 30 is the last one whose cell count fits in Int64; level 31 would
    # overflow, which is why the range stops there.
    @test ncells(levelgrid(SYS, 30)) == 6917529027641081856
    @test_throws ArgumentError levelgrid(SYS, -1)
    @test_throws ArgumentError levelgrid(SYS, 31)

    roots = rootcells(SYS)
    @test roots == [LevelIndex(0, i) for i in 0:5]
    @test issorted(roots) && allunique(roots)
end

@testset "cellindex / cellposition is the identity up to the base offset" begin
    for l in (0, 1, 3, 12, 30)
        g = levelgrid(SYS, l)
        n = ncells(g)
        for i in (1, 2, n ÷ 2, n)
            c = cellindex(g, i)
            @test c === LevelIndex(l, i - 1)
            @test cellposition(g, c) == i
            @test rawid(c) == i - 1
        end
        @test_throws BoundsError cellindex(g, 0)
        @test_throws BoundsError cellindex(g, n + 1)
        # Not in the grid: out of range, or the right ordinal at a wrong level.
        @test cellposition(g, LevelIndex(l, -1)) === nothing
        @test cellposition(g, LevelIndex(l, n)) === nothing
        l < 30 && @test cellposition(g, LevelIndex(l + 1, 0)) === nothing
        l > 0 && @test cellposition(g, LevelIndex(l - 1, 0)) === nothing
    end
end

@testset "parent / children / ancestor / descendants, exhaustively" begin
    for l in 0:4
        g = levelgrid(SYS, l)
        for i in 1:ncells(g)
            c = cellindex(g, i)
            kids = children(SYS, c)
            @test kids == [LevelIndex(l + 1, 4 * rawid(c) + k) for k in 0:3]
            @test issorted(kids) && allunique(kids)
            for k in kids
                @test parent(SYS, k) == c
                @test ancestor(SYS, k, l) == c
            end
            @test ancestor(SYS, c, l) === c
            @test descendants(SYS, c, l) == [c]
            @test descendants(SYS, c, l + 1) == kids
            @test descendants(SYS, c, l + 2) ==
                  [LevelIndex(l + 2, j) for j in (16 * rawid(c)):(16 * rawid(c) + 15)]
        end
    end
    @test_throws ArgumentError parent(SYS, LevelIndex(0, 0))
    @test_throws ArgumentError children(SYS, LevelIndex(30, 0))
    @test_throws ArgumentError ancestor(SYS, LevelIndex(2, 0), 3)
    @test_throws ArgumentError descendants(SYS, LevelIndex(2, 0), 1)
    @test_throws ArgumentError descendant_range(SYS, LevelIndex(2, 0), 1)
    @test_throws ArgumentError descendant_range(SYS, LevelIndex(2, 0), 31)

    # The ancestor shortcut at real depth: `>> 2Δ` must agree with `parent`
    # applied Δ times, which is where a closed form goes wrong at depth ≥ 2.
    for l in (8, 20, 30), h in cell_sample(min(l, 6), 24)
        c = LevelIndex(l, h)
        walked = c
        for j in (l - 1):-1:0
            walked = parent(SYS, walked)
            @test ancestor(SYS, c, j) == walked
        end
    end
end

@testset "descendant_range is exact and hole-free" begin
    for l in 0:3
        g = levelgrid(SYS, l)
        for i in 1:ncells(g)
            c = cellindex(g, i)
            @test descendant_range(SYS, c, l) == i:i
            for d in 1:3
                target = levelgrid(SYS, l + d)
                r = descendant_range(SYS, c, l + d)
                actual = descendants(SYS, c, l + d)
                @test length(r) == 4^d
                # Two-sided: every descendant's POSITION is in the range, and
                # every position in the range is a descendant.
                @test [cellposition(target, x) for x in actual] == collect(r)
            end
            # Sibling ranges partition the parent's, in order.
            if l < 30
                kids = children(SYS, c)
                @test reduce(vcat, [collect(descendant_range(SYS, k, l + 1)) for k in kids]) ==
                      collect(descendant_range(SYS, c, l + 1))
            end
        end
    end
    # And at the extreme: the six roots' level-30 ranges tile the whole level.
    total = 0
    for c in rootcells(SYS)
        r = descendant_range(SYS, c, 30)
        @test length(r) == 4^30
        total += length(r)
    end
    @test total == ncells(levelgrid(SYS, 30))
end

# =========================================================================
# 4. Geometry
# =========================================================================

@testset "cell_boundary is the cell, exactly" begin
    for l in 0:3
        g = levelgrid(SYS, l)
        nside = Int64(1) << l
        worst = Inf
        for i in 1:ncells(g)
            c = cellindex(g, i)
            ring = cell_boundary(g, c)
            # FOUR vertices — not a densification. S2 cell edges are geodesics.
            @test length(ring) == 4
            ix, iy, face = hilbert_to_xyf(rawid(c), nside)
            chart = cell_corners(ix, iy, face, nside)
            # Bitwise the chart's own corners: one evaluation, not two that
            # agree to a tolerance.
            @test all(k -> ring[k] === chart[k], 1:4)
            @test ring[1] != ring[end]              # implicitly closed
            worst = min(worst, ccw_measure(ring))
        end
        @test worst > 0
    end
    # The geometry entry points guard the id: a cell from another level is an
    # error, not the silent geometry of a cell that does not exist.
    g = levelgrid(SYS, 2)
    @test_throws ArgumentError cell_boundary(g, LevelIndex(3, 0))
    @test_throws ArgumentError cell_centroid(g, LevelIndex(1, 0))
    @test_throws ArgumentError cell_boundary(g, LevelIndex(2, 96))
    @test_throws ArgumentError cell_boundary(g, LevelIndex(2, -1))
end

@testset "cell_centroid is the S2 cell centre, and interior" begin
    for l in 0:3
        g = levelgrid(SYS, l)
        nside = Int64(1) << l
        for i in 1:ncells(g)
            c = cellindex(g, i)
            ix, iy, face = hilbert_to_xyf(rawid(c), nside)
            @test cell_centroid(g, c) === cell_center(ix, iy, face, nside)
            @test sum(abs2, cell_centroid(g, c)) ≈ 1.0 atol = 1e-15
        end
    end
end

@testset "the cells partition the sphere" begin
    # The generic `cell_area` is the spherical area of the four-corner ring, and
    # for S2 that ring IS the cell — so the level total is 4π to rounding, with
    # no densification error to absorb. (Contrast HEALPix, whose densified rings
    # sum to 4π only to ~3e-3.)
    for l in 0:3
        g = levelgrid(SYS, l)
        areas = [cell_area(g, cellindex(g, i)) for i in 1:ncells(g)]
        @test all(>(0), areas)
        rel = abs(sum(areas) - 4π) / (4π)
        record!("|Σ cell_area - 4π| / 4π (level $l)", rel)
        @test rel < 1e-14
    end
    # S2 is NOT equal-area, and this is how far from it: the within-level spread
    # of the quadratic projection, which converges to about 2.08.
    g = levelgrid(SYS, 5)
    areas = [cell_area(g, cellindex(g, i)) for i in 1:ncells(g)]
    spread = maximum(areas) / minimum(areas)
    record!("within-level cell-area spread (level 5)", spread)
    @test 1.5 < spread < 2.1
end

# =========================================================================
# 5. node_extent — the exact four-corner cap
# =========================================================================

@testset "node_extent covers the subtree" begin
    # Measure the covering law: every descendant boundary
    # point of every sampled cell, several levels down, against the ancestor's
    # own cap. The margin is expected to be tiny and NEGATIVE — descendant
    # corners land on the ancestor's corners bit-identically, and `nextfloat` on
    # the radius is what keeps them strictly inside instead of on the border.
    worst = -Inf
    for l in 0:3
        g = levelgrid(SYS, l)
        for i in 1:ncells(g)
            c = cellindex(g, i)
            cap = node_extent(SYS, c)
            for d in 1:4
                target = levelgrid(SYS, l + d)
                r = descendant_range(SYS, c, l + d)
                step = max(1, length(r) ÷ 8)
                for pos in first(r):step:last(r)
                    for p in cell_boundary(target, cellindex(target, pos))
                        worst = max(worst, US.spherical_distance(cap.point, p) - cap.radius)
                    end
                end
            end
        end
    end
    record!("worst descendant overshoot of node_extent (rad)", worst)
    @test worst < 0

    # ... and against a dense sampling of the cell's own chart rectangle, which
    # is where every descendant at EVERY depth lives (children tile their parent
    # exactly, so this is the whole subtree and not a proxy for it).
    worst_chart = -Inf
    for l in 0:2
        g = levelgrid(SYS, l)
        nside = Int64(1) << l
        for i in 1:ncells(g)
            c = cellindex(g, i)
            cap = node_extent(SYS, c)
            ix, iy, face = hilbert_to_xyf(rawid(c), nside)
            for a in 0:16, b in 0:16
                p = stf_to_point((ix + a / 16) / nside, (iy + b / 16) / nside, face)
                worst_chart = max(worst_chart, US.spherical_distance(cap.point, p) - cap.radius)
            end
        end
    end
    record!("worst chart-rectangle overshoot of node_extent (rad)", worst_chart)
    @test worst_chart < 0
end

@testset "node_extent is convex, and tighter than the inflated default" begin
    # Radius ≤ 90° is what makes sampling a boundary's VERTICES a sound proxy
    # for the whole boundary, which is how the conformance harness checks the
    # covering law with `require_convex_extents` at its default.
    biggest = -Inf
    for l in 0:3
        g = levelgrid(SYS, l)
        for i in 1:ncells(g)
            cap = node_extent(SYS, cellindex(g, i))
            @test cap.radius <= π / 2
            biggest = max(biggest, cap.radius)
        end
    end
    record!("widest node_extent radius (rad)", biggest)
    # A level-0 face: acos(1/√3), the cube-corner half-angle.
    @test node_extent(SYS, LevelIndex(0, 0)).radius ≈ acos(1 / sqrt(3)) atol = 1e-15

    # Strictly tighter than the generic `cap_inflation`-scaled default the
    # override replaces. `cap_inflation` itself is never consulted, and stays at
    # the interface default.
    @test DGG.cap_inflation(SYS) == 1.2
    for l in 0:2
        g = levelgrid(SYS, l)
        for i in 1:ncells(g)
            c = cellindex(g, i)
            cap = node_extent(SYS, c)
            centre = cell_centroid(g, c)
            r = maximum(US.spherical_distance(centre, p) for p in cell_boundary(g, c))
            @test cap.radius < r * DGG.cap_inflation(SYS)
        end
    end
end

# =========================================================================
# 6. Topology: the cube-edge seam, geometrically
# =========================================================================

@testset "the seam table is the geometry (level $level)" for level in 0:3
    # The neighbour relation computed from `SEAM`'s integer
    # arithmetic must equal the one computed from SHARED BOUNDARY POINTS, over
    # a whole level, seams and cube corners included.
    #
    #   * `Edge()` neighbours share TWO boundary points — a whole cell edge,
    #     since an S2 ring is four corners and its edges are geodesics;
    #   * the corner-only neighbours, `setdiff(Vertex(), Edge())`, share exactly
    #     ONE.
    g = levelgrid(SYS, level)
    n = ncells(g)

    incidence = Dict{NTuple{3,Float64},Vector{LevelIndex}}()
    for i in 1:n, p in cell_boundary(g, cellindex(g, i))
        push!(get!(() -> LevelIndex[], incidence, vkey(p)), cellindex(g, i))
    end

    for i in 1:n
        c = cellindex(g, i)
        # The geometric answer: count how many boundary points each other cell
        # shares with `c`.
        shared = Dict{LevelIndex,Int}()
        for p in cell_boundary(g, c), o in incidence[vkey(p)]
            o == c && continue
            shared[o] = get(shared, o, 0) + 1
        end
        geo_v = Set(k for (k, v) in shared if v >= 1)
        geo_e = Set(k for (k, v) in shared if v >= 2)

        vs = collect(neighbors(g, c, 1))
        es = collect(neighbors(g, c, 1; connectivity = Edge()))
        @test Set(vs) == geo_v
        @test Set(es) == geo_e
        @test issubset(Set(es), Set(vs))

        # The counts, spelled out rather than inferred from the sets.
        ring_c = Set(vkey.(cell_boundary(g, c)))
        for x in es
            @test length(intersect(ring_c, Set(vkey.(cell_boundary(g, x))))) == 2
        end
        for x in setdiff(vs, es)
            @test length(intersect(ring_c, Set(vkey.(cell_boundary(g, x))))) == 1
        end
    end
end

@testset "neighbour counts: 8 in a face, 7 in a cube corner, 4 at level 0" begin
    # Level 0 is its own case: a cell IS a face, and the four faces sharing a
    # cube edge are the whole neighbourhood — the three faces meeting at a cube
    # corner are pairwise edge-adjacent already, so `Vertex()` adds nothing.
    g0 = levelgrid(SYS, 0)
    for i in 1:6
        c = cellindex(g0, i)
        @test length(neighbors(g0, c, 1)) == 4
        @test length(neighbors(g0, c, 1; connectivity = Edge())) == 4
    end

    for level in 1:3
        g = levelgrid(SYS, level)
        nside = Int64(1) << level
        sevens = 0
        for i in 1:ncells(g)
            c = cellindex(g, i)
            ix, iy, face = hilbert_to_xyf(rawid(c), nside)
            corner = (ix == 0 || ix == nside - 1) && (iy == 0 || iy == nside - 1)
            @test length(neighbors(g, c, 1)) == (corner ? 7 : 8)
            @test length(neighbors(g, c, 1; connectivity = Edge())) == 4
            corner && (sevens += 1)
        end
        # Four cells per face sit in a face corner, at every level ≥ 1.
        @test sevens == 24
    end

    # `wrap_xyf` is what makes the corner case a 7 rather than an 8: the only
    # step it refuses is the one that would cross a cube CORNER.
    @test wrap_xyf(1, 1, 0, 4) == (1, 1, 0)
    @test wrap_xyf(-1, -1, 0, 4) === nothing
    @test wrap_xyf(4, 4, 0, 4) === nothing
    @test wrap_xyf(-1, 4, 0, 4) === nothing
    @test wrap_xyf(4, 2, 0, 4) !== nothing
    # And the table is a genuine involution on cube edges: stepping out of a
    # face and back in returns the cell you left.
    for level in 1:3
        nside = Int64(1) << level
        for face in 0:5, k in 0:(nside - 1)
            for (out, back) in (((nside, k), 1), ((-1, k), 3), ((k, nside), 2), ((k, -1), 4))
                w = wrap_xyf(out[1], out[2], face, nside)
                @test w !== nothing
                # The cell we landed on is adjacent to `(k, ·)` on `face`, so it
                # must list the cell we came from among its own neighbours.
                src = xyf_to_hilbert(clamp(out[1], 0, nside - 1), clamp(out[2], 0, nside - 1),
                                     face, nside)
                dst = xyf_to_hilbert(w[1], w[2], w[3], nside)
                gg = levelgrid(SYS, level)
                @test LevelIndex(level, src) in neighbors(gg, LevelIndex(level, dst), 1)
            end
        end
    end
end

@testset "ring 1 winds counter-clockwise seen from outside" begin
    # MEASURED, not argued: the azimuths of the neighbour centres about the cell
    # centre, in the order returned, must wrap exactly once. A clockwise cycle
    # wraps `n - 1` times and an id-sorted list wraps arbitrarily. This is the
    # check that caught both hex systems' "natural" orders running the wrong
    # way; S2's lattice order needs no reversal, and this is why we know.
    for level in 0:3, conn in (Vertex(), Edge())
        g = levelgrid(SYS, level)
        for i in 1:ncells(g)
            c = cellindex(g, i)
            @test DiscreteGlobalGridsConformanceTesting.winding_problems(
                g, c, neighbors(g, c, 1; connectivity = conn);
                label = "ring 1 ($conn)") == String[]
        end
    end
end

@testset "ring 1 starts where the docstring says" begin
    # ORACLE PIN on the documented START of the rotational order.
    #
    # The neighbouring testsets check the CYCLE (CCW winding, `Edge()`
    # preserving it), and the conformance harness's winding law is
    # start-invariant by construction — every rotation of a CCW cycle is still a
    # CCW cycle. So the START is documented and otherwise unpinned, and a
    # rotation would pass everything while silently rotating the weight vector
    # of any stencil that has baked "slot 1 is +s" into itself.
    #
    # Values produced by the implementation and checked against the documented
    # rule: the cycle `+s, +s+t, +t, -s+t, -s, -s-t, -t, +s-t` of
    # `NEIGHBOR_OFFSETS`, read off the lattice.
    @test NEIGHBOR_OFFSETS ==
          ((1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1))

    g = levelgrid(SYS, 3)

    # (a) A FACE-INTERIOR cell: all eight offsets stay on face 1, so this pins
    #     the cycle itself with no seam arithmetic in the way.
    c = cellindex(g, 100)
    @test rawid(c) == 99
    @test hilbert_to_xyf(99, 8) == (5, 4, 1)
    @test [rawid(x) for x in neighbors(g, c, 1)] ==
          [100, 103, 98, 97, 96, 95, 92, 91]
    @test [rawid(x) for x in neighbors(g, c, 1; connectivity = Edge())] ==
          [100, 98, 96, 92]
    # ... and the ids really are the lattice offsets, in that order.
    @test [rawid(x) for x in neighbors(g, c, 1)] ==
          [xyf_to_hilbert(5 + dx, 4 + dy, 1, 8) for (dx, dy) in NEIGHBOR_OFFSETS]

    # (b) A CUBE-CORNER cell: the `-s-t` diagonal has no cell to name, so the
    #     cycle is seven long and the two seam crossings sit where the cycle
    #     puts them.
    corner = cellindex(g, 1)
    @test rawid(corner) == 0
    @test hilbert_to_xyf(0, 8) == (0, 0, 0)
    @test [rawid(x) for x in neighbors(g, corner, 1)] ==
          [3, 2, 1, 297, 298, 383, 382]
    @test [rawid(x) for x in neighbors(g, corner, 1; connectivity = Edge())] ==
          [3, 1, 298, 383]
end

@testset "ring is the tail block of neighbors" begin
    # The composition contract: `neighbors(k)` is rings 1..k concatenated
    # outward, so `ring(k)` is exactly the trailing block.
    for level in (1, 2, 3), conn in (Vertex(), Edge())
        g = levelgrid(SYS, level)
        for h in cell_sample(level, 24)
            c = LevelIndex(level, h)
            acc = LevelIndex[]
            for k in 1:3
                shell = collect(ring(g, c, k; connectivity = conn))
                append!(acc, shell)
                disc = collect(neighbors(g, c, k; connectivity = conn))
                @test disc == acc                                    # concatenated outward
                @test disc[(end - length(shell) + 1):end] == shell    # tail block
                @test allunique(disc)
                @test !(c in disc)
            end
            @test collect(ring(g, c, 0; connectivity = conn)) == [c]
            @test isempty(neighbors(g, c, 0; connectivity = conn))
            @test_throws ArgumentError neighbors(g, c, -1)
            @test_throws ArgumentError ring(g, c, -1)
        end
    end
end

@testset "ring shells are disjoint and outer rings wind CCW too" begin
    for level in (2, 3), conn in (Vertex(), Edge())
        g = levelgrid(SYS, level)
        for h in cell_sample(level, 16)
            c = LevelIndex(level, h)
            seen = Set{LevelIndex}()
            for k in 1:3
                shell = ring(g, c, k; connectivity = conn)
                @test isempty(intersect(Set(shell), seen))
                union!(seen, shell)
                @test Set(neighbors(g, c, k; connectivity = conn)) == seen
                @test DiscreteGlobalGridsConformanceTesting.winding_problems(
                    g, c, shell; label = "ring $k ($conn)") == String[]
            end
            @test !(c in seen)
        end
    end
end

@testset "neighbours are symmetric over a whole level" begin
    # Exhaustive both-directions symmetry, which the conformance harness can
    # only sample. Cheap at these sizes, and it is where a seam table typically
    # breaks: an asymmetric entry names a cell that does not name you back.
    for level in 0:3, conn in (Vertex(), Edge())
        g = levelgrid(SYS, level)
        adjacency = Dict{LevelIndex,Set{LevelIndex}}()
        for i in 1:ncells(g)
            c = cellindex(g, i)
            adjacency[c] = Set(neighbors(g, c, 1; connectivity = conn))
        end
        for (c, ns) in adjacency
            @test !(c in ns)
            for x in ns
                @test c in adjacency[x]
            end
        end
    end
end

# =========================================================================
# 7. Location
# =========================================================================

@testset "cellat inverts cell_centroid everywhere" begin
    for level in 0:4
        g = levelgrid(SYS, level)
        for i in 1:ncells(g)
            c = cellindex(g, i)
            @test cellat(g, cell_centroid(g, c)) == c
        end
    end
    # Deep levels, sampled: the closed form is 60 bits wide at level 30.
    for level in (10, 20, MAX_LEVEL)
        g = levelgrid(SYS, level)
        for h in cell_sample(min(level, 5), 48)
            c = LevelIndex(level, h)
            @test cellat(g, cell_centroid(g, c)) == c
        end
    end
end

@testset "cellat locates interior points, not just centroids" begin
    # A nearest-centroid lookup would pass the test above and be wrong for every
    # tessellation that is not a Voronoi one. These are interior points that are
    # NOT the centroid, walked out towards each corner.
    for level in 0:3
        g = levelgrid(SYS, level)
        for i in 1:ncells(g)
            c = cellindex(g, i)
            centroid = cell_centroid(g, c)
            for v in cell_boundary(g, c)
                @test cellat(g, US.slerp(centroid, v, 0.5)) == c
            end
        end
    end
    # A complete S2 level grid covers the sphere, so `cellat` is never `nothing`.
    rng = MersenneTwister(20_260_813)
    g = levelgrid(SYS, 4)
    for _ in 1:512
        x, y, z = randn(rng), randn(rng), randn(rng)
        n = sqrt(x^2 + y^2 + z^2)
        c = cellat(g, GO.UnitSphericalPoint(x / n, y / n, z / n))
        @test c !== nothing
        @test cellposition(g, c) !== nothing
    end
    # The lon/lat wrapper agrees with the unit-sphere primitive.
    for (lon, lat) in ((0.0, 0.0), (-73.9, 40.7), (139.7, 35.7), (0.0, 90.0), (0.0, -90.0))
        @test cellat(g, lon, lat) == cellat(g, US.UnitSphereFromGeographic()((lon, lat)))
    end
end

# =========================================================================
# 8. Contract: the conformance suites, default kwargs
# =========================================================================

@testset "conformance" begin
    for l in (0, 1, 3)
        test_grid_interface(levelgrid(SYS, l); label = "S2 level $l")
    end
    test_hierarchical_system(SYS)
end

end # @testset "S2 system (T11)"

@printf("[S2] measured:\n")
for key in sort!(collect(keys(MEASURED)))
    @printf("[S2]   %-56s %+.3e\n", key, MEASURED[key])
end

end # module S2SystemTests
