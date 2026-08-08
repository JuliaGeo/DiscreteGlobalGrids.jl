# test_engine.jl — Eisenstein digit engine (src/IGeo7/engine.jl).
#
# Assertable vectors come from spec/aperture7-indexing-spec.md section 7 (the
# convention-independent machinery, exercised with that spec's canonical digit
# map delta(d) = UNITS[d] and uniform chirality c = 3+omega) plus the fitted
# IGEO7 convention constants (chi alternation, sigma) from
# spec/igeo7-geometry-diagnosis.md §4.

using Test

# --- helpers built on the engine primitives, in the SPEC-CANONICAL convention
# (uniform chi = c, delta(d) = UNITS[d]); these mirror spec section 2.2 / 2.4.
function _spec_horner(ds)
    a = Int64(0); b = Int64(0)
    for d in ds
        (a, b) = IGeo7.mul_c(a, b)
        if d != 0
            u = IGeo7.UNITS[d]
            a += u[1]; b += u[2]
        end
    end
    return (a, b)
end

function _spec_decode(a::Int64, b::Int64, r::Int)
    ds = zeros(Int, r)
    for k in r:-1:1
        j = IGeo7.RES_TO_J_C[IGeo7.res_c(a, b)+1]
        if j >= 0
            u = IGeo7.UNITS[j+1]
            a -= u[1]; b -= u[2]
            ds[k] = j + 1                      # spec-canonical: digit = j+1
        end
        (a, b) = IGeo7.div_c(a, b)
    end
    return ds, (a, b)
end

@testset "engine" begin

    @testset "Eisenstein basics" begin
        @test IGeo7.OMEGA ≈ cis(2pi / 3)
        @test IGeo7.UNITS == ((1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1))
        for (j, u) in enumerate(IGeo7.UNITS)
            @test IGeo7.norm_eis(u[1], u[2]) == 1
            @test abs(mod(rad2deg(angle(IGeo7.ecpx(u...))) - 60.0 * (j - 1) + 180, 360) - 180) < 1e-12
        end
        @test IGeo7.norm_eis(3, 1) == 7      # N(c) = 7
        @test IGeo7.norm_eis(2, -1) == 7     # N(cbar) = 7
        @test IGeo7.ALPHA_DEG ≈ 19.106605350869096 atol = 1e-12
        @test abs(rad2deg(angle(IGeo7.ecpx(3, 1))) - IGeo7.ALPHA_DEG) < 1e-12
        @test abs(rad2deg(angle(IGeo7.ecpx(2, -1))) + IGeo7.ALPHA_DEG) < 1e-12
    end

    @testset "multiplication / exact division" begin
        @test IGeo7.mul_c(1, 0) == (3, 1)
        @test IGeo7.mul_cbar(1, 0) == (2, -1)
        for a in -50:50, b in -50:50
            # matrix form vs complex arithmetic
            (ca, cb) = IGeo7.mul_c(a, b)
            @test abs(IGeo7.ecpx(ca, cb) - IGeo7.ecpx(3, 1) * IGeo7.ecpx(a, b)) < 1e-9
            (da, db) = IGeo7.mul_cbar(a, b)
            @test abs(IGeo7.ecpx(da, db) - IGeo7.ecpx(2, -1) * IGeo7.ecpx(a, b)) < 1e-9
            # c * cbar == 7
            @test IGeo7.mul_cbar(ca, cb) == (7a, 7b)
            @test IGeo7.mul_c(da, db) == (7a, 7b)
            # exact rounded division inverts multiplication
            @test IGeo7.div_c(ca, cb) == (a, b)
            @test IGeo7.div_cbar(da, db) == (a, b)
            # residue of a multiple is 0
            @test IGeo7.res_c(ca, cb) == 0
            @test IGeo7.res_cbar(da, db) == 0
        end
    end

    @testset "residue tables" begin
        # the seven digits are the seven cosets (spec S2.3)
        @test sort([IGeo7.res_c(0, 0); [IGeo7.res_c(u...) for u in IGeo7.UNITS]]) == collect(0:6)
        @test sort([IGeo7.res_cbar(0, 0); [IGeo7.res_cbar(u...) for u in IGeo7.UNITS]]) == collect(0:6)
        # residue -> unit index (-1 = the zero digit)
        @test IGeo7.RES_TO_J_C == (-1, 0, 4, 5, 2, 1, 3)
        @test IGeo7.RES_TO_J_CBAR == (-1, 0, 2, 1, 4, 5, 3)
        # spec S2.3 tables, in the spec-canonical map digit = j+1
        @test [j < 0 ? 0 : j + 1 for j in IGeo7.RES_TO_J_C] == [0, 1, 5, 6, 3, 2, 4]
        @test [j < 0 ? 0 : j + 1 for j in IGeo7.RES_TO_J_CBAR] == [0, 1, 3, 2, 5, 6, 4]
    end

    @testset "unit rotations" begin
        for a in -12:12, b in -12:12
            @test IGeo7.unitmul(Int64(a), Int64(b), 0) == (a, b)
            @test IGeo7.unitmul(Int64(a), Int64(b), 6) == (a, b)
            @test IGeo7.unitmul(Int64(a), Int64(b), -1) ==
                  IGeo7.unitmul(Int64(a), Int64(b), 5)
            for g in -3:3
                (ra, rb) = IGeo7.unitmul(Int64(a), Int64(b), g)
                @test abs(IGeo7.ecpx(ra, rb) -
                          IGeo7.ecpx(a, b) * cis(deg2rad(60.0 * g))) < 1e-12
                @test IGeo7.norm_eis(ra, rb) == IGeo7.norm_eis(a, b)
            end
            # group law
            for g1 in -2:2, g2 in -2:2
                z = IGeo7.unitmul(Int64(a), Int64(b), g1)
                @test IGeo7.unitmul(z[1], z[2], g2) ==
                      IGeo7.unitmul(Int64(a), Int64(b), g1 + g2)
            end
        end
    end

    @testset "spec S7 worked examples (spec-canonical convention)" begin
        cases = ([1, 3] => (3, 2), [2, 0, 5] => (2, 7), [6, 6, 6] => (6, -6),
            [1, 1, 1] => (12, 6), [4, 2] => (-2, 0), [0, 3, 1] => (0, 2),
            [5, 0, 0, 2] => (0, -18), [0, 0, 4, 3, 3] => (-9, -2))
        for (ds, z) in cases
            @test _spec_horner(ds) == z
            got, leftover = _spec_decode(Int64(z[1]), Int64(z[2]), length(ds))
            @test got == ds
            @test leftover == (0, 0)
        end
        # the decode trace of (2,0,5) spelled out in spec S2.6
        @test IGeo7.res_c(2, 7) == 2
        @test IGeo7.div_c(2 + 1, 7 + 1) == (2, 3)
        @test IGeo7.res_c(2, 3) == 0
        @test IGeo7.div_c(2, 3) == (1, 1)
        @test IGeo7.res_c(1, 1) == 5
    end

    @testset "powers of c and growth bound" begin
        z = (Int64(1), Int64(0))
        pows = [z]
        for r in 1:20
            z = IGeo7.mul_c(z...)
            push!(pows, z)
        end
        @test pows[3] == (8, 5)
        @test pows[4] == (19, 18)
        @test pows[5] == (39, 55)
        @test pows[6] == (62, 149)
        @test pows[21] == (323103824, 122884025)
        for r in 0:20
            @test IGeo7.norm_eis(pows[r+1]...) == Int64(7)^r
        end
    end

    @testset "fitted convention tables" begin
        # [fitted; see spec/igeo7-geometry-diagnosis.md §3 — the chirality
        # alternates cbar-first: cbar at odd levels, c at even]
        @test !IGeo7.chi_is_c(1) && IGeo7.chi_is_c(2)
        @test all(IGeo7.chi_is_c(k) == iseven(k) for k in 1:20)
        # [fitted, A-gauge; see spec/igeo7-geometry-diagnosis.md §4]
        # digit -> unit index (dev angle 60j from the reference edge)
        @test IGeo7.SIGMA_J == (5, 3, 4, 1, 0, 2)
        @test IGeo7.DIGIT_OF_J == (5, 4, 6, 2, 3, 1)
        for d in 1:6
            @test IGeo7.DIGIT_OF_J[IGeo7.SIGMA_J[d]+1] == d
            @test IGeo7.sigma(d) == IGeo7.UNITS[IGeo7.SIGMA_J[d]+1]
        end
        @test IGeo7.sigma(0) == (0, 0)
        # GBT complement pairs (1,6) (2,5) (3,4) are antipodal
        for (d, e) in ((1, 6), (2, 5), (3, 4))
            @test IGeo7.sigma(d) .+ IGeo7.sigma(e) == (0, 0)
        end
        # digit angles: 5->0, 4->60, 6->120, 2->180, 3->240, 1->300 — the
        # published GBT cycle 4,6,2,3,1,5 CCW seen from outside [z7 spec 3.4]
        for (d, ang) in ((5, 0), (4, 60), (6, 120), (2, 180), (3, 240), (1, 300))
            @test mod(60 * IGeo7.SIGMA_J[d], 360) == ang
        end
        # residue -> digit under the fitted sigma
        @test IGeo7.RES_TO_DIGIT_C == (0, 5, 3, 1, 6, 4, 2)
        @test IGeo7.RES_TO_DIGIT_CBAR == (0, 5, 6, 4, 3, 1, 2)
        for k in 1:20, d in 0:6
            u = IGeo7.sigma(d)
            @test IGeo7.digit_of_res(IGeo7.reschi(u[1], u[2], k), k) == d
        end
    end

    @testset "chi products P_r" begin
        # P_r = prod_{k=1..r} chi_k with alternating chirality: c*cbar = 7, so
        # P_{2m} = 7^m and P_{2m+1} = 7^m * cbar — argP alternates 0 / -alpha.
        @test IGeo7.P_R[1] == (1, 0)
        @test IGeo7.P_R[2] == (2, -1)
        @test IGeo7.P_R[3] == (7, 0)
        for m in 0:10
            @test IGeo7.P_R[2m+1] == (Int64(7)^m, 0)
            2m + 1 <= 20 && @test IGeo7.P_R[2m+2] == (2 * Int64(7)^m, -Int64(7)^m)
        end
        for r in 0:20
            @test IGeo7.norm_eis(IGeo7.P_R[r+1]...) == Int64(7)^r
            @test IGeo7.ARGP_DEG[r+1] ≈
                  (iseven(r) ? 0.0 : 360.0 - IGeo7.ALPHA_DEG) atol = 1e-12
        end
        # each level's chi divides the running product exactly
        for r in 1:20
            @test IGeo7.mulchi(IGeo7.P_R[r]..., r) == IGeo7.P_R[r+1]
            @test IGeo7.divchi(IGeo7.P_R[r+1]..., r) == IGeo7.P_R[r]
        end
    end

    @testset "Horner encode / digit decode (fitted convention)" begin
        @test IGeo7.horner(Int[]) == (0, 0)
        @test IGeo7.horner([0, 0, 0, 0]) == (0, 0)
        @test IGeo7.horner([5]) == (1, 0)          # sigma(5) = unit j=0
        @test IGeo7.horner([1]) == (0, -1)         # sigma(1) = unit j=5
        # level 2 uses c: horner([5,5]) = c*(1,0) + (1,0) = (4,1)
        @test IGeo7.horner([5, 5]) == (4, 1)
        # round-trip over random digit strings at every resolution
        rng_state = 12345
        rand16() = (rng_state = (1103515245 * rng_state + 12345) % 2147483648)
        for r in 1:20, _ in 1:20
            ds = [rand16() % 7 for _ in 1:r]
            (a, b) = IGeo7.horner(ds)
            got = zeros(Int, r)
            la, lb = IGeo7.digit_decode!(got, a, b, r)
            @test got == ds
            @test (la, lb) == (0, 0)
            # streaming form used by grid.jl: one exact step per level
            aa, bb = a, b
            for k in r:-1:1
                d, aa, bb = IGeo7.decode_step(aa, bb, k)
                @test d == ds[k]
            end
            @test (aa, bb) == (0, 0)
        end
        # growth bound of spec S2.5: |a|,|b| <= (2/sqrt3)(7^(r/2)-1)/(sqrt7-1)
        for r in 1:8
            bound = (2 / sqrt(3)) * (7.0^(r / 2) - 1) / (sqrt(7) - 1)
            for _ in 1:200
                ds = [rand16() % 7 for _ in 1:r]
                (a, b) = IGeo7.horner(ds)
                @test max(abs(a), abs(b)) <= bound + 1e-9
            end
        end
        # a nonzero leftover marks a lattice point outside this base's subtree
        la, lb = IGeo7.digit_decode!(zeros(Int, 2), Int64(5), Int64(0), 2)
        @test (la, lb) != (0, 0)
    end

    @testset "hex rounding" begin
        @test IGeo7.hex_round(3.0 + 1e-9, 2.0 + 1e-9) == (3, 2)
        @test IGeo7.hex_round(3.0 - 1e-9, 2.0 - 1e-9) == (3, 2)
        @test IGeo7.hex_round(0.0, 0.0) == (0, 0)
        for (a, b) in IGeo7.UNITS
            @test IGeo7.hex_round(Float64(a), Float64(b)) == (a, b)
        end
        # idempotence on exact lattice points, and the cube-sum invariant
        for a in -20:20, b in -20:20
            @test IGeo7.hex_round(Float64(a), Float64(b)) == (a, b)
        end
        # spec S3.4: P/b0 = (19+omega)/49 times c^2 lands on cell (3,2)
        w = (19 + IGeo7.OMEGA) / 49 * IGeo7.ecpx(8, 5)
        x, y = real(w), imag(w)
        @test IGeo7.hex_round(x + y / sqrt(3), 2y / sqrt(3)) == (3, 2)
        # rounding really returns the nearest lattice point (brute force)
        rs = 987654321
        nextr() = ((rs = (1103515245 * rs + 12345) % 2147483648) / 2147483648)
        for _ in 1:2000
            x = 12 * (nextr() - 0.5); y = 12 * (nextr() - 0.5)
            at = x + y / sqrt(3); bt = 2y / sqrt(3)
            (ra, rb) = IGeo7.hex_round(at, bt)
            p = complex(x, y)
            dbest = abs(p - IGeo7.ecpx(ra, rb))
            for da in -3:3, db in -3:3
                @test abs(p - IGeo7.ecpx(ra + da, rb + db)) >= dbest - 1e-12
            end
        end
    end

    @testset "allocation-free scalar paths" begin
        IGeo7.mul_c(3, 4); IGeo7.div_c(3, 4); IGeo7.decode_step(Int64(9), Int64(4), 3)
        IGeo7.hex_round(1.2, 3.4); IGeo7.unitmul(Int64(3), Int64(4), 2)
        IGeo7.mulchi(Int64(3), Int64(4), 5)
        @test @allocated(IGeo7.mul_c(3, 4)) == 0
        @test @allocated(IGeo7.decode_step(Int64(9), Int64(4), 3)) == 0
        @test @allocated(IGeo7.hex_round(1.2, 3.4)) == 0
        @test @allocated(IGeo7.unitmul(Int64(3), Int64(4), 2)) == 0
        @test @allocated(IGeo7.mulchi(Int64(3), Int64(4), 5)) == 0
        # Horner over a tuple of digits allocates nothing
        ds = (1, 0, 4, 3, 3)
        IGeo7.horner(ds)
        @test @allocated(IGeo7.horner(ds)) == 0
    end

end
