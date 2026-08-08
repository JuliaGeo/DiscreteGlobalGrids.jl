module HealpixChartTestSuite

# Tests for `src/HEALPix/chart.jl`: the pure closed-form HEALPix face charts
# and the RING / NESTED index maps over the `nside × nside` lattice.
#
# `chart.jl` deliberately carries no Healpix.jl dependency, so Healpix.jl is
# free to be the *ground truth* here rather than the implementation: at
# power-of-two `nside` every index map and every geometric output is checked
# against `nest2ring` / `ring2nest` / `pix2vecNest` / `boundariesRing`.
#
# Healpix.jl's `Resolution` only admits `nside = 2^k` (verified below), so at
# non-power-of-two `nside` — which the RING maps and the charts support and
# which is the whole reason this layer exists separately from the nested
# kernel — the ground truth is instead a set of analytic invariants taken
# straight from Górski et al. 2005 (DOI:10.1086/427976): the ring-length and
# ring-latitude closed forms, exact corner sharing on the lattice, and the
# Euler characteristic of the resulting tessellation.

using Test
using Printf
import Healpix
import GeometryOps as GO

using DiscreteGlobalGrids.HEALPix: JRLL, JPLL,
    xyf_to_point, pixel_corners, pixel_center,
    ring_nlon, ring_first, xyf_to_ring, ring_to_xyf,
    xyf_to_nested, nested_to_xyf, nested_to_ring, ring_to_nested

# Componentwise sup-norm between a chart point and a plain 3-tuple/vector.
maxdev(p, q) = max(abs(p[1] - q[1]), abs(p[2] - q[2]), abs(p[3] - q[3]))

cross3(a, b) = (a[2] * b[3] - a[3] * b[2],
                a[3] * b[1] - a[1] * b[3],
                a[1] * b[2] - a[2] * b[1])

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
const MEASURED = Dict{String,Float64}()
record!(key, value) = (MEASURED[key] = max(get(MEASURED, key, 0.0), value))

# ---------------------------------------------------------------------------

@testset "chart constants" begin
    @test length(JRLL) == 12 && length(JPLL) == 12
    # Rows 2 / 3 / 4 of the base tiling: north caps, equatorial belt, south caps.
    @test JRLL == (2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4)
    # Columns in units of 45°: the two cap rows are offset by 45° from the belt,
    # which is what interlocks the 12 faces into a rhombic dodecahedral tiling.
    @test all(isodd, JPLL[1:4]) && all(iseven, JPLL[5:8]) && all(isodd, JPLL[9:12])
    @test JPLL[1:4] == JPLL[9:12]                     # caps share their longitudes
end

@testset "xyf_to_point — poles, norms, branch boundaries" begin
    worst_norm = 0.0
    for face in 0:11, i in 0:20, j in 0:20
        p = xyf_to_point(i / 20, j / 20, face)
        worst_norm = max(worst_norm, abs(sum(abs2, p) - 1))
    end
    record!("unit-norm residual (|p|^2 - 1)", worst_norm)
    @test worst_norm <= 1e-15

    # The north pole is the (x, y) = (1, 1) corner of all four north cap faces
    # and it must come out *exactly*, not to a tolerance: `have_sintheta` exists
    # precisely so this corner does not drift.
    for face in 0:3
        @test xyf_to_point(1, 1, face) == GO.UnitSphericalPoint(0.0, 0.0, 1.0)
    end
    for face in 8:11
        @test xyf_to_point(0, 0, face) == GO.UnitSphericalPoint(0.0, 0.0, -1.0)
    end

    # The cap/belt seam is `jr == 1` (north) and `jr == 3` (south), i.e.
    # `z = ±2/3` — the HEALPix transition latitude.
    for face in 0:11
        @test xyf_to_point(1, 1, face)[3] ≈ (face < 4 ? 1.0 : face < 8 ? 2 / 3 : 0.0) atol = 1e-15
        @test xyf_to_point(0, 0, face)[3] ≈ (face < 4 ? 0.0 : face < 8 ? -2 / 3 : -1.0) atol = 1e-15
    end

    # Equator: `jr == 2`, i.e. the chart line `x + y == JRLL - 2`. (Not bitwise
    # zero — `jr` is formed as `JRLL - x - y` in floating point, so it lands
    # within an ulp of 2 rather than on it.)
    for face in 0:11
        s = JRLL[face + 1] - 2
        0 <= s <= 2 || continue
        for t in 0:10
            x = clamp(s * t / 10, 0, 1)
            y = s - x
            0 <= y <= 1 || continue
            @test xyf_to_point(x, y, face)[3] ≈ 0.0 atol = 1e-15
        end
    end

    # Cross-face continuity at the cap/belt seam: the north corner of belt face
    # 4 and the west corner of cap face 0 are the same physical vertex, and the
    # charts must agree there (this is the seam a discontinuous port breaks).
    @test maxdev(xyf_to_point(1, 1, 4), xyf_to_point(0, 1, 0)) <= 1e-15
    @test maxdev(xyf_to_point(1, 1, 4), xyf_to_point(1, 0, 3)) <= 1e-15
end

# ---------------------------------------------------------------------------
# Power-of-two `nside`: Healpix.jl is ground truth for everything.
# ---------------------------------------------------------------------------

@testset "power-of-two nside vs Healpix.jl (nside = $nside)" for (nside, stride) in
        ((1, 1), (2, 1), (4, 1), (8, 1), (16, 5))

    res = Healpix.Resolution(nside)
    npix = 12 * nside^2
    pixels = 0:stride:(npix - 1)

    # --- NESTED bijectivity -------------------------------------------------
    seen = Set{NTuple{3,Int}}()
    for p in 0:(npix - 1)
        ix, iy, face = nested_to_xyf(p, nside)
        @test 0 <= ix < nside && 0 <= iy < nside && 0 <= face <= 11
        push!(seen, (ix, iy, face))
        @test xyf_to_nested(ix, iy, face, nside) == p
    end
    @test length(seen) == npix                        # onto the full lattice

    # --- NESTED <-> RING against Healpix.jl ---------------------------------
    # Healpix.jl numbers pixels 1-based in both schemes; our nested ids are
    # 0-based (EOPF) and our ring indices 1-based, hence the asymmetric `+ 1`.
    for p in 0:(npix - 1)
        ix, iy, face = nested_to_xyf(p, nside)
        @test xyf_to_ring(ix, iy, face, nside) == Healpix.nest2ring(res, p + 1)
        @test nested_to_ring(p, nside) == Healpix.nest2ring(res, p + 1)
    end
    for ipix in 1:npix
        @test ring_to_nested(ipix, nside) == Healpix.ring2nest(res, ipix) - 1
        ix, iy, face = ring_to_xyf(ipix, nside)
        @test xyf_to_ring(ix, iy, face, nside) == ipix          # round trip
    end

    # --- Centers vs pix2vecNest ---------------------------------------------
    worst_center = 0.0
    for p in pixels
        ix, iy, face = nested_to_xyf(p, nside)
        worst_center = max(worst_center,
            maxdev(pixel_center(ix, iy, face, nside), Healpix.pix2vecNest(res, p + 1)))
    end
    record!("centers vs Healpix.pix2vecNest", worst_center)
    @test worst_center <= 1e-12

    # --- Corners vs boundariesRing ------------------------------------------
    # Corner *order* is not guaranteed to agree between the two, so each of our
    # corners is matched to its nearest Healpix corner and the resulting cyclic
    # offsets are collected: a single offset across every pixel of every face is
    # the real assertion (a per-pixel nearest match alone would also accept a
    # scrambled ring).
    worst_corner = 0.0
    offsets = Set{Int}()
    worst_ccw = Inf
    for p in pixels
        ix, iy, face = nested_to_xyf(p, nside)
        mine = pixel_corners(ix, iy, face, nside)
        b = Healpix.boundariesRing(res, Healpix.nest2ring(res, p + 1), 1, Float64)
        theirs = ntuple(k -> (b[k, 1], b[k, 2], b[k, 3]), 4)
        for i in 1:4
            devs = ntuple(k -> maxdev(mine[i], theirs[k]), 4)
            k = argmin(devs)
            push!(offsets, mod(k - i, 4))
            worst_corner = max(worst_corner, devs[k])
        end
        worst_ccw = min(worst_ccw, ccw_measure(mine))
    end
    record!("corners vs Healpix.boundariesRing", worst_corner)
    @test length(offsets) == 1                        # one consistent ring rotation
    @test worst_corner <= 1e-9
    # CCW as seen from outside the sphere — the convex-clip kernel clips a CW
    # ring to EMPTY, so this is a correctness contract, not cosmetics.
    @test worst_ccw > 0

    # --- Exact corner sharing on the lattice --------------------------------
    for face in 0:11, ix in 0:(nside - 2), iy in 0:(nside - 2)
        c = pixel_corners(ix, iy, face, nside)
        east = pixel_corners(ix + 1, iy, face, nside)
        north = pixel_corners(ix, iy + 1, face, nside)
        @test c[1] === east[2] && c[4] === east[3]    # bit-identical, not ≈
        @test c[1] === north[4] && c[2] === north[3]
    end
end

@testset "nested id space matches the 12 * 4^level convention" begin
    # The nested id is `face * nside² + morton(ix, iy)`, so face `f` owns
    # exactly `[f * nside², (f+1) * nside²)` — this is the block structure
    # `HealpixKernel.jl`'s `has_ordinal_ids` derives the whole hierarchy from.
    for level in 0:5
        nside = 2^level
        for face in 0:11
            @test xyf_to_nested(0, 0, face, nside) == face * nside^2
            @test xyf_to_nested(nside - 1, nside - 1, face, nside) == (face + 1) * nside^2 - 1
        end
    end
    # One level up = drop the low two Morton bits = halve both coordinates.
    for level in 1:6
        nside = 2^level
        for face in (0, 5, 11), ix in (0, 1, nside ÷ 2, nside - 1), iy in (0, nside - 1)
            p = xyf_to_nested(ix, iy, face, nside)
            @test p ÷ 4 == xyf_to_nested(ix ÷ 2, iy ÷ 2, face, nside ÷ 2)
        end
    end
end

# ---------------------------------------------------------------------------
# Non-power-of-two `nside`: analytic invariants.
# ---------------------------------------------------------------------------

@testset "Healpix.jl rejects non-power-of-two nside" begin
    # Documents *why* the block below validates against closed forms instead of
    # against Healpix.jl. (Healpix.jl signals this by calling `DomainError()`
    # with no arguments, so what actually escapes is a `MethodError` — hence
    # the deliberately loose `Exception`.)
    for nside in (3, 5, 7)
        @test_throws Exception Healpix.Resolution(nside)
    end
    @test Healpix.Resolution(4) isa Healpix.Resolution
end

@testset "arbitrary nside — RING maps and chart invariants (nside = $nside)" for nside in (3, 5, 7)
    npix = 12 * nside^2

    # --- The nested index must refuse to exist here -------------------------
    @test_throws ArgumentError xyf_to_nested(0, 0, 0, nside)
    @test_throws ArgumentError nested_to_xyf(0, nside)
    @test_throws ArgumentError nested_to_ring(0, nside)
    @test_throws ArgumentError ring_to_nested(1, nside)

    # --- Ring structure closed forms ---------------------------------------
    # `ring_first` must be the exact prefix sum of `ring_nlon`, and the rings
    # must exhaust the sphere: `Σ nlon == 12 nside²`.
    total = 0
    for jr in 1:(4nside - 1)
        @test ring_first(jr, nside) == total + 1
        nlon = ring_nlon(jr, nside)
        # Górski's ring lengths: 4jr growing through the north cap, 4nside
        # across the belt, mirrored south.
        expected = jr < nside ? 4jr : (jr <= 3nside ? 4nside : 4 * (4nside - jr))
        @test nlon == expected
        @test nlon == ring_nlon(4nside - jr, nside)   # north/south mirror symmetry
        total += nlon
    end
    @test total == npix

    # --- RING bijectivity over the whole lattice ----------------------------
    perface = [Set{NTuple{2,Int}}() for _ in 0:11]
    for ipix in 1:npix
        ix, iy, face = ring_to_xyf(ipix, nside)
        @test 0 <= ix < nside && 0 <= iy < nside && 0 <= face <= 11
        @test xyf_to_ring(ix, iy, face, nside) == ipix
        push!(perface[face + 1], (ix, iy))
    end
    @test all(s -> length(s) == nside^2, perface)     # each face fully covered
    for face in 0:11, ix in 0:(nside - 1), iy in 0:(nside - 1)
        @test (ix, iy) in perface[face + 1]
        @test ring_to_xyf(xyf_to_ring(ix, iy, face, nside), nside) == (ix, iy, face)
    end
    @test_throws ArgumentError ring_to_xyf(0, nside)
    @test_throws ArgumentError ring_to_xyf(npix + 1, nside)

    # --- Ring geometry: constant z, Górski latitudes, increasing φ ----------
    worst_z = 0.0
    worst_norm = 0.0
    for jr in 1:(4nside - 1)
        first = ring_first(jr, nside)
        nlon = ring_nlon(jr, nside)
        # Górski et al. 2005 eqs. (4)/(5), evaluated at the ring's pixel centers
        # (`jr` here is the integer ring index, `jr / nside` the chart's
        # continuous ring coordinate).
        zref = if jr < nside
            1 - jr^2 / (3 * nside^2)
        elseif jr <= 3nside
            (2nside - jr) * 2 / (3 * nside)
        else
            (4nside - jr)^2 / (3 * nside^2) - 1
        end
        phis = Float64[]
        for ipix in first:(first + nlon - 1)
            ix, iy, face = ring_to_xyf(ipix, nside)
            c = pixel_center(ix, iy, face, nside)
            worst_z = max(worst_z, abs(c[3] - zref))
            worst_norm = max(worst_norm, abs(sum(abs2, c) - 1))
            push!(phis, mod(atan(c[2], c[1]), 2π))
        end
        # A ring is an iso-latitude circle traversed west→east in RING order.
        @test issorted(phis; lt = <)
    end
    record!("ring-center z vs Górski closed form", worst_z)
    record!("unit-norm residual (|p|^2 - 1)", worst_norm)
    @test worst_z <= 1e-14
    @test worst_norm <= 1e-15

    # --- Exact tessellation: shared lattice corners are bit-identical -------
    for face in 0:11, ix in 0:(nside - 2), iy in 0:(nside - 2)
        c = pixel_corners(ix, iy, face, nside)
        east = pixel_corners(ix + 1, iy, face, nside)
        north = pixel_corners(ix, iy + 1, face, nside)
        @test c[1] === east[2] && c[4] === east[3]
        @test c[1] === north[4] && c[2] === north[3]
    end

    # --- Winding and Euler characteristic -----------------------------------
    worst_ccw = Inf
    corners = Set{NTuple{3,Int}}()
    for face in 0:11, ix in 0:(nside - 1), iy in 0:(nside - 1)
        c = pixel_corners(ix, iy, face, nside)
        worst_ccw = min(worst_ccw, ccw_measure(c))
        for p in c
            push!(corners, (round(Int, p[1] * 1e10), round(Int, p[2] * 1e10), round(Int, p[3] * 1e10)))
        end
    end
    @test worst_ccw > 0
    # V - E + F = 2 with F = 12 nside² quads and E = 2F ⇒ V = 12 nside² + 2.
    # Distinct corners hitting that count is a global statement: the charts
    # agree across every face seam, and no pixel is missing or doubled.
    @test length(corners) == npix + 2
end

@printf("[HEALPix chart] measured maxima:\n")
for key in sort!(collect(keys(MEASURED)))
    @printf("[HEALPix chart]   %-42s %.3e\n", key, MEASURED[key])
end

end # module HealpixChartTestSuite
