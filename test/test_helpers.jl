module HelpersTestSuite

using Test
using DiscreteGlobalGrids
using DiscreteGlobalGrids.Helpers

@testset "shared helpers" begin
    empty = empty_small_list(Val(3), 0)
    @test isempty(empty)
    @test isbitstype(typeof(empty))

    full = small_push(small_push(small_push(empty, 2), 1), 2)
    @test collect(full) == [2, 1, 2]
    @test collect(small_sort(full)) == [1, 2, 2]
    @test_throws BoundsError small_push(full, 3)
    @test_throws BoundsError full[4]

    @test strictly_increasing(Int[])
    @test strictly_increasing([1, 2, 3])
    @test !strictly_increasing([1, 1, 2])
    @test !strictly_increasing([2, 1])

    @test sorted_index([1, 3, 5], 3) == 2
    @test sorted_index([1, 3, 5], 2) == 0
    @test to_uint64_id("ff") == 0xff
    @test to_uint64_id("0xFF") == 0xff
end

# ---------------------------------------------------------------------------
# Authalic latitude
#
# Two independent reference sets, deliberately kept separate:
#
#   * BIGFLOAT_* — the closed form ξ = asin(q/q_p) (Snyder 1987 eq. 3-11/3-12)
#     evaluated at `setprecision(BigFloat, 300)` and rounded once to Float64,
#     with the inverse obtained by 400-step bisection on the same function.
#     This is the correctly rounded answer, so the series is checked against it
#     for *bitwise* equality, not a tolerance.
#
#   * PROJ_* — PROJ 9's `+proj=cea` (Lambert cylindrical equal area). PROJ has
#     no public authalic-latitude operator, but ellipsoidal `cea` with
#     `lat_ts=0` computes `y(φ) = a·q(φ)/2` through its own `pj_qsfn`, so
#     `ξ = asin(y(φ)/y(90°))` and its inverse follow with no series anywhere.
#     Generated with (Proj.jl 1.9, PROJ 9.x):
#
#         fwd = Proj.Transformation("+proj=longlat +ellps=WGS84 +no_defs",
#                   "+proj=cea +ellps=WGS84 +lat_ts=0 +lon_0=0 +units=m +no_defs")
#         yp  = fwd(0.0, 90.0)[2]
#         beta(d) = rad2deg(asin(fwd(0.0, d)[2] / yp))
#         phi(d)  = Proj.inv(fwd)(0.0, yp * sind(d))[2]
#
#     PROJ's own `asin(q/q_p)` is ill-conditioned near the poles (condition
#     number 1/cos ξ), so these agree only to ~1e-12 degrees there — which is
#     exactly why the bitwise reference above is BigFloat and not PROJ.
# ---------------------------------------------------------------------------

const BIGFLOAT_GEODETIC_TO_AUTHALIC = (
    (0.0, 0.0), (1.0, 0.9955309909569205), (5.0, 4.977763096123115),
    (10.0, 9.956198098935735), (15.0, 14.935956949386629),
    (22.5, 22.40940211374809), (30.0, 29.888997034459564),
    (40.0, 39.87369373453434), (45.0, 44.87170287343394),
    (50.0, 49.87361022076773), (60.0, 59.888785569885165),
    (70.0, 69.91741174086397), (75.0, 74.93574548414333),
    (80.0, 79.95604114354843), (85.0, 84.97767958186002),
    (89.0, 88.995513957862), (89.9, 89.89955130506617), (90.0, 90.0),
    (-30.0, -29.888997034459564), (-45.0, -44.87170287343394),
    (-60.0, -59.888785569885165), (-89.0, -88.995513957862), (-90.0, -90.0),
)

const BIGFLOAT_AUTHALIC_TO_GEODETIC = (
    (0.0, 0.0), (1.0, 1.0044890625471037), (5.0, 5.022335225436584),
    (10.0, 10.04398667355641), (15.0, 15.064291967518724),
    (22.5, 22.59088523889082), (30.0, 30.11125171864826),
    (40.0, 40.1264043490986), (45.0, 45.12829693352109),
    (50.0, 50.126291349769936), (60.0, 60.11096559341369),
    (70.0, 70.08230543334226), (75.0, 75.06400584025931),
    (80.0, 80.04377430172227), (85.0, 85.02222222460503),
    (89.0, 89.00446601553388), (89.9, 89.90044669066359), (90.0, 90.0),
    (-30.0, -30.11125171864826), (-45.0, -45.12829693352109),
    (-60.0, -60.11096559341369), (-89.0, -89.00446601553388), (-90.0, -90.0),
)

const PROJ_GEODETIC_TO_AUTHALIC = (
    (0.0, 0.0), (1.0, 0.99553099095692066), (5.0, 4.9777630961231152),
    (15.0, 14.935956949386627), (30.0, 29.888997034459557),
    (45.0, 44.871702873433939), (60.0, 59.888785569885165),
    (75.0, 74.935745484143339), (89.0, 88.99551395786284), (90.0, 90.0),
    (-45.0, -44.871702873433939), (-89.0, -88.99551395786284),
)

const PROJ_AUTHALIC_TO_GEODETIC = (
    (0.0, 0.0), (1.0, 1.0044890625471037), (5.0, 5.0223352254365841),
    (15.0, 15.064291967518722), (30.0, 30.111251718648251),
    (45.0, 45.128296933521099), (60.0, 60.110965593413681),
    (75.0, 75.064005840259327), (89.0, 89.004466015533623), (90.0, 90.0),
    (-45.0, -45.128296933521099), (-89.0, -89.004466015533623),
)

@testset "authalic latitude" begin
    t = WGS84_AUTHALIC

    @testset "ellipsoid parameterization" begin
        # `f` and `1/f` reach the third flattening by the same `f/(2 − f)`, so
        # they agree bitwise.
        from_invf = AuthalicTransform(; inverse_flattening=298.257223563)
        from_f = AuthalicTransform(; flattening=1 / 298.257223563)
        @test from_invf.fwd == t.fwd
        @test from_invf.inv == t.inv
        @test from_f.fwd == t.fwd
        @test from_f.inv == t.inv

        # The `e²` route reconstructs `n` as `e²/(1 + √(1 − e²))²` instead,
        # which is a different (equally cancellation-free) rounding of the
        # same quantity — so it lands within a few ulps, not bitwise. What
        # must hold is that the transforms are interchangeable in use.
        from_e2 = AuthalicTransform(; eccentricity_squared=t.eccentricity_squared)
        @test all(isapprox.(from_e2.fwd, t.fwd; rtol=1e-14))
        @test all(isapprox.(from_e2.inv, t.inv; rtol=1e-14))
        for φ in range(-90.0, 90.0; length=181)
            @test geodetic_to_authalicd(from_e2, φ) ≈ geodetic_to_authalicd(t, φ) atol = 1e-13
        end

        # e² = f(2 − f), reproduced from the defining 1/f.
        @test t.eccentricity_squared ≈ (1 / 298.257223563) * (2 - 1 / 298.257223563)
        @test t.semimajor_axis == 6378137.0
        @test isbitstype(typeof(t))
        @test occursin("AuthalicTransform{Float64}", sprint(show, t))

        Shape = DiscreteGlobalGrids.Helpers.EllipsoidShapeError
        @test_throws Shape AuthalicTransform()
        @test_throws Shape AuthalicTransform(; flattening=0.1, eccentricity_squared=0.1)
        @test occursin("exactly one", sprint(showerror, Shape(0)))
        @test_throws DomainError AuthalicTransform(; eccentricity_squared=1.0)
        @test_throws DomainError AuthalicTransform(; eccentricity_squared=-0.1)
        @test_throws DomainError AuthalicTransform(; flattening=1.5)
        # Prolate (f < 0) would need atan in place of atanh; rejected, not NaN.
        @test_throws DomainError AuthalicTransform(; flattening=-0.1)

        # Reprecision to the same type is a no-op, not a rebuild.
        @test AuthalicTransform{Float64}(t) === t
    end

    @testset "authalic radius" begin
        # Snyder (3-13). This is the correctly rounded value: PROJ's cea gives
        # sqrt(a·y(90°)) = 6.371007180918474e6 and so does a 300-bit
        # evaluation of the closed form.
        @test authalic_radius(t) === 6.371007180918474e6
        @test authalic_radius(t) === t.authalic_radius
        @test authalic_radius(6378137.0, t.eccentricity_squared) === authalic_radius(t)

        # The repo's own constant, used for grid areas, is ONE ULP above the
        # correctly rounded value (9.3e-10 m — 0.9 nm — high). Pinned here so
        # the discrepancy is recorded rather than rediscovered; `R_AUTHALIC`
        # is a documented contract constant and is deliberately not changed.
        repo = DiscreteGlobalGrids.ISEA.R_AUTHALIC
        @test repo != authalic_radius(t)
        @test reinterpret(Int64, repo) - reinterpret(Int64, authalic_radius(t)) == 1
        @test isapprox(repo, authalic_radius(t); rtol=2eps(Float64))

        # Sphere: the authalic radius is the sphere's own radius, exactly.
        sphere = AuthalicTransform(; eccentricity_squared=0.0)
        @test authalic_radius(sphere) === 6378137.0
        @test authalic_radius(1.0, 0.0) === 1.0
        @test_throws DomainError authalic_radius(1.0, 1.0)
    end

    @testset "reference values — 300-bit BigFloat (bitwise)" begin
        for (φ, ξ) in BIGFLOAT_GEODETIC_TO_AUTHALIC
            @test geodetic_to_authalicd(t, φ) === ξ
        end
        for (ξ, φ) in BIGFLOAT_AUTHALIC_TO_GEODETIC
            @test authalic_to_geodeticd(t, ξ) === φ
        end
        # The radian entry points must land on the same values.
        for (φ, ξ) in BIGFLOAT_GEODETIC_TO_AUTHALIC
            @test rad2deg(geodetic_to_authalic(t, deg2rad(φ))) ≈ ξ atol = 1e-13
        end
        for (ξ, φ) in BIGFLOAT_AUTHALIC_TO_GEODETIC
            @test rad2deg(authalic_to_geodetic(t, deg2rad(ξ))) ≈ φ atol = 1e-13
        end
    end

    @testset "reference values — PROJ +proj=cea" begin
        # 1e-11 deg ≈ 1.1 mm: loose enough for PROJ's asin conditioning near
        # the poles, tight enough that a wrong series or ellipsoid fails.
        for (φ, ξ) in PROJ_GEODETIC_TO_AUTHALIC
            @test geodetic_to_authalicd(t, φ) ≈ ξ atol = 1e-11
        end
        for (ξ, φ) in PROJ_AUTHALIC_TO_GEODETIC
            @test authalic_to_geodeticd(t, ξ) ≈ φ atol = 1e-11
        end
    end

    @testset "round trip" begin
        worst_deg = 0.0
        worst_rad = 0.0
        for φ in range(-90.0, 90.0; length=4001)
            worst_deg = max(worst_deg,
                abs(authalic_to_geodeticd(t, geodetic_to_authalicd(t, φ)) - φ))
            r = deg2rad(φ)
            worst_rad = max(worst_rad,
                abs(authalic_to_geodetic(t, geodetic_to_authalic(t, r)) - r))
        end
        # Claimed: 1 ulp of the largest latitude in the sweep, i.e. eps(90.0)
        # in degrees and eps(π/2) in radians — about 1.4 nm on the authalic
        # sphere either way.
        @test worst_deg <= eps(90.0)
        @test worst_rad <= eps(Float64(π) / 2)

        # ... and the other way round, authalic → geodetic → authalic.
        worst = 0.0
        for ξ in range(-90.0, 90.0; length=4001)
            worst = max(worst, abs(geodetic_to_authalicd(t, authalic_to_geodeticd(t, ξ)) - ξ))
        end
        @test worst <= eps(90.0)
    end

    @testset "closed form agreement" begin
        # Independent in-package check: the series against asin(q/q_p) built
        # from `authalic_q` (Snyder 3-12). Restricted to |φ| ≤ 80° because
        # asin loses ~3 digits at the poles — see `authalic_q`'s docstring.
        e2 = t.eccentricity_squared
        qp = authalic_q(e2, 1.0)
        @test qp ≈ 2 * (authalic_radius(t) / t.semimajor_axis)^2
        @test authalic_q(e2, 0.0) == 0.0
        @test authalic_q(0.0, 0.5) === 1.0          # sphere: q = 2 sin φ

        worst = 0.0
        for φd in range(-80.0, 80.0; length=1601)
            φ = deg2rad(φd)
            worst = max(worst, abs(asin(authalic_q(e2, sin(φ)) / qp) - geodetic_to_authalic(t, φ)))
        end
        @test worst < 1e-14                          # rad; ~6e-8 m
    end

    @testset "edge cases" begin
        # Poles and equator are fixed points to the last bit, both directions,
        # both conventions. `sincospi` makes this structural, not incidental.
        for f in (geodetic_to_authalicd, authalic_to_geodeticd)
            @test f(t, 90.0) === 90.0
            @test f(t, -90.0) === -90.0
            @test f(t, 0.0) === 0.0
        end
        for f in (geodetic_to_authalic, authalic_to_geodetic)
            @test f(t, π / 2) === π / 2
            @test f(t, -π / 2) === -π / 2
            @test f(t, 0.0) === 0.0
        end

        # Antisymmetry: ξ(−φ) = −ξ(φ), exactly.
        for φ in (1.0, 17.5, 45.0, 63.25, 88.0)
            @test geodetic_to_authalicd(t, -φ) === -geodetic_to_authalicd(t, φ)
            @test authalic_to_geodeticd(t, -φ) === -authalic_to_geodeticd(t, φ)
        end

        # The sphere is the identity, with no NaN from the e² = 0 division in
        # `q` and no ±0.0 leaking out of the all-zero coefficients.
        sphere = AuthalicTransform(; eccentricity_squared=0.0)
        @test all(iszero, sphere.fwd)
        @test all(iszero, sphere.inv)
        for φ in (-90.0, -45.0, 0.0, 33.3, 90.0)
            @test geodetic_to_authalicd(sphere, φ) === φ
            @test authalic_to_geodeticd(sphere, φ) === φ
            @test !isnan(geodetic_to_authalicd(sphere, φ))
        end

        # A far more oblate body than WGS84 still round-trips (the n-series is
        # in f/(2−f), so this is a much harder case).
        oblate = AuthalicTransform(; flattening=1 / 30)
        @test geodetic_to_authalicd(oblate, 90.0) === 90.0
        worst = 0.0
        for φ in range(-90.0, 90.0; length=721)
            worst = max(worst, abs(authalic_to_geodeticd(oblate, geodetic_to_authalicd(oblate, φ)) - φ))
        end
        @test worst < 1e-9
    end

    @testset "Float32" begin
        t32 = AuthalicTransform{Float32}(WGS84_AUTHALIC)
        @test t32 isa AuthalicTransform{Float32}
        @test eltype(t32.fwd) === Float32
        @test isbitstype(typeof(t32))

        # Coefficients are computed in Float64 and narrowed once, so they are
        # the correctly rounded Float32 versions of the Float64 ones.
        @test t32.fwd == map(Float32, WGS84_AUTHALIC.fwd)
        @test t32.inv == map(Float32, WGS84_AUTHALIC.inv)
        @test t32.authalic_radius === Float32(WGS84_AUTHALIC.authalic_radius)

        @test geodetic_to_authalicd(t32, 45.0f0) isa Float32
        @test geodetic_to_authalicd(t32, 45.0f0) ≈ 44.871702873433939f0 atol = 1f-4
        @test geodetic_to_authalicd(t32, 90.0f0) === 90.0f0
        @test geodetic_to_authalicd(t32, -90.0f0) === -90.0f0
        @test geodetic_to_authalicd(t32, 0.0f0) === 0.0f0
        @test authalic_to_geodeticd(t32, 90.0f0) === 90.0f0

        worst = 0.0f0
        for φ in range(-90.0f0, 90.0f0; length=1001)
            worst = max(worst, abs(authalic_to_geodeticd(t32, geodetic_to_authalicd(t32, φ)) - φ))
        end
        @test worst < 1.0f-4
    end

    @testset "type stability and allocation" begin
        for f in (geodetic_to_authalic, authalic_to_geodetic,
            geodetic_to_authalicd, authalic_to_geodeticd)
            @test @inferred(f(t, 0.7)) isa Float64
            @test @inferred(f(AuthalicTransform{Float32}(t), 0.7f0)) isa Float32
            @test @allocated(f(t, 0.7)) == 0
        end
        @test @inferred(authalic_radius(t)) isa Float64
        @test @inferred(AuthalicTransform(; flattening=0.003)) isa AuthalicTransform{Float64}

        # An Int latitude must not break inference or change the result type.
        @test geodetic_to_authalicd(t, 45) === geodetic_to_authalicd(t, 45.0)
        @test @inferred(geodetic_to_authalicd(t, 45)) isa Float64
    end
end

end # module HelpersTestSuite
