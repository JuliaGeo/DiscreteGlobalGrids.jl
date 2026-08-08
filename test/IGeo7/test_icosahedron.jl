# test_icosahedron.jl — locks the standard ISEA icosahedron placement of
# src/ISEA/icosahedron.jl (spec/design.md sections 2, 5, 12), plus the res-0
# lock on the IGEO7 public API built on top of it.
#
# Provenance of every expected value asserted here:
#   * vertex layout (closed form)   -> spec/isea-projection-spec.md §4.1–4.2 /
#                                      spec/z7-paper-spec.md §1.2
#   * res-0 centers                 -> test/IGeo7/vectors/dggrid_res0_centers.csv
#                                      (dggrid CLI output; its own print noise
#                                      is <= 5.913e-9 deg, so the comparison
#                                      tolerance is 1e-8, not 1e-9)
#   * authalic R, deleted digits    -> spec/interface-contract.md
#   * adjacency / antipodal pairs   -> spec/aperture7-indexing-spec.md S5.1
#   * per-base reference edge       -> fitted; see
#                                      spec/igeo7-geometry-diagnosis.md §4

using Test

const _ICO_VEC = joinpath(@__DIR__, "vectors")

"rows of dggrid_res0_centers.csv as (base, lon, lat)"
function _load_dggrid_centers(path)
    rows = Tuple{Int,Float64,Float64}[]
    for line in eachline(path)
        isempty(strip(line)) && continue
        p = split(line, ',')
        push!(rows, (parse(Int, p[1]), parse(Float64, p[2]), parse(Float64, p[3])))
    end
    return rows
end

@testset "icosahedron" begin

    @testset "sphere helpers" begin
        @test ISEA.lonlat_to_xyz(0.0, 0.0) == (1.0, 0.0, 0.0)
        @test all(abs.(ISEA.lonlat_to_xyz(90.0, 0.0) .- (0.0, 1.0, 0.0)) .< 1e-16)
        @test ISEA.lonlat_to_xyz(0.0, 90.0)[3] == 1.0
        for (lon, lat) in ((0.0, 0.0), (37.5, -12.25), (-179.0, 88.0), (123.456, -45.678))
            v = ISEA.lonlat_to_xyz(lon, lat)
            @test abs(ISEA.vnorm(v) - 1) < 1e-15
            lo, la = ISEA.xyz_to_lonlat(v)
            @test abs(lo - lon) < 1e-12 && abs(la - lat) < 1e-12
        end
        a = ISEA.lonlat_to_xyz(10.0, 20.0)
        b = ISEA.lonlat_to_xyz(-30.0, 44.0)
        @test ISEA.angdist(a, a) < 1e-15
        @test abs(ISEA.angdist(a, b) - acosd(ISEA.vdot(a, b))) < 1e-9
        @test abs(ISEA.vdot(ISEA.vcross(a, b), a)) < 1e-15
    end

    @testset "exact constants" begin
        # [contract] authalic radius; res-0 cell area = 4 pi R^2 / 12
        @test ISEA.R_AUTHALIC == 6371007.180918475
        @test 4pi * ISEA.R_AUTHALIC^2 / 12 ≈ 4.250546847700739e13 rtol = 1e-14
        # cos of the icosahedron edge arc
        @test ISEA.ADJ_DOT ≈ 1 / sqrt(5) atol = 1e-16
        @test abs(acosd(ISEA.ADJ_DOT) - 63.43494882292201) < 1e-10
        # the standard ISEA placement constants [isea spec 4.1]
        @test ISEA.ISEA_LON0 == 11.25
        @test ISEA.ISEA_LAT_HI === atand((1 + sqrt(5.0)) / 2)
        @test abs(ISEA.ISEA_LAT_HI - 58.282525588538995) < 1e-13
    end

    @testset "vertex table vs dggrid_res0_centers.csv" begin
        rows = _load_dggrid_centers(joinpath(_ICO_VEC, "dggrid_res0_centers.csv"))
        @test length(rows) == 12
        @test sort(first.(rows)) == collect(0:11)
        # TOLERANCE: the CSV carries DGGRID's own numerical noise, up to
        # 5.913e-9 deg (unstructured, in both lon and lat). 1e-9 is NOT
        # achievable against that file; 1e-8 is the honest tolerance. The
        # lower guard (worst > 1e-10) documents that the residual is the
        # CSV's, not ours: our vertices ARE the exact closed form.
        worst = 0.0
        for (b, lon, lat) in rows
            d = ISEA.angdist(ISEA.vertex(b), ISEA.lonlat_to_xyz(lon, lat))
            @test d < 1e-8
            worst = max(worst, d)
        end
        @test worst < 6e-9
        @test worst > 1e-10
    end

    @testset "vertex table closed form" begin
        V = ISEA.VERTICES
        @test length(V) == 12
        @test all(abs(ISEA.vnorm(v) - 1) < 1e-15 for v in V)
        # vertex 0 at (11.25 E, atand(phi) N) exactly [isea spec 4.1]
        lon0, lat0 = ISEA.xyz_to_lonlat(ISEA.vertex(0))
        @test abs(lon0 - 11.25) < 1e-12
        @test abs(lat0 - ISEA.ISEA_LAT_HI) < 1e-12
        # latitude classes per base [isea spec 4.2 / T1]:
        # 0,1 -> +LAT_HI; 2,5 -> +(90-LAT_HI); 3,4,6,10 -> 0;
        # 7,9 -> -(90-LAT_HI); 8,11 -> -LAT_HI
        hi, mid = ISEA.ISEA_LAT_HI, 90.0 - ISEA.ISEA_LAT_HI
        expected_lat = (hi, hi, mid, 0.0, 0.0, mid, 0.0, -mid, -hi, -mid, 0.0, -hi)
        for b in 0:11
            _, lat = ISEA.xyz_to_lonlat(ISEA.vertex(b))
            @test abs(lat - expected_lat[b+1]) < 1e-12
        end
        # every pair is adjacent, antipodal, or at the complementary arc
        for i in 0:11, j in 0:11
            i < j || continue
            d = ISEA.vdot(ISEA.vertex(i), ISEA.vertex(j))
            @test min(abs(d - ISEA.ADJ_DOT), abs(d + ISEA.ADJ_DOT),
                abs(d + 1)) < 1e-14
        end
        @test_throws ArgumentError ISEA.vertex(12)
        @test_throws ArgumentError ISEA.vertex(-1)
    end

    @testset "neighbor rings" begin
        # [spec/aperture7-indexing-spec.md S5.1]
        expected = ((1, 2, 3, 4, 5), (0, 2, 5, 6, 10), (0, 1, 3, 6, 7),
            (0, 2, 4, 7, 8), (0, 3, 5, 8, 9), (0, 1, 4, 9, 10),
            (1, 2, 7, 10, 11), (2, 3, 6, 8, 11), (3, 4, 7, 9, 11),
            (4, 5, 8, 10, 11), (1, 5, 6, 9, 11), (6, 7, 8, 9, 10))
        for b in 0:11
            @test ISEA.NEIGHBORS[b+1] == expected[b+1]
            @test Tuple(sort(collect(ISEA.NBRS_CCW[b+1]))) == expected[b+1]
        end
        # antipodal pairs [spec S5.1]
        for (i, j) in ((0, 11), (1, 8), (2, 9), (3, 10), (4, 6), (5, 7))
            @test ISEA.vdot(ISEA.vertex(i), ISEA.vertex(j)) ≈ -1 atol = 1e-15
        end
        # adjacency is symmetric
        for b in 0:11, n in ISEA.NEIGHBORS[b+1]
            @test b in ISEA.NEIGHBORS[n+1]
        end
    end

    @testset "per-base reference edge and CCW rings (fitted)" begin
        # [fitted; see spec/igeo7-geometry-diagnosis.md §4]
        @test ISEA.REFERENCE_EDGE == (1, 10, 6, 7, 8, 9, 11, 11, 11, 11, 11, 8)
        for b in 0:11
            @test ISEA.REFERENCE_EDGE[b+1] in ISEA.NEIGHBORS[b+1]
            @test ISEA.NBRS_CCW[b+1][1] == ISEA.REFERENCE_EDGE[b+1]
        end
        # the ring really is counterclockwise (seen from outside) on a perfect
        # 72-degree comb: measure azimuths in a locally built tangent frame
        for b in 0:11
            Vb = ISEA.vertex(b)
            tangent(n) = ISEA.vnormalize(ISEA.vsub(ISEA.vertex(n),
                ISEA.vscale(Vb, ISEA.vdot(ISEA.vertex(n), Vb))))
            ut = tangent(ISEA.REFERENCE_EDGE[b+1])
            wt = ISEA.vcross(Vb, ut)          # ut x wt = Vb: CCW from outside
            for (k, n) in enumerate(ISEA.NBRS_CCW[b+1])
                t = tangent(n)
                az = mod(atand(ISEA.vdot(t, wt), ISEA.vdot(t, ut)), 360.0)
                @test abs(mod(az - 72.0 * (k - 1) + 180, 360) - 180) < 1e-8
            end
        end
    end

    @testset "deleted digit per base" begin
        # [contract: pentagon_chains.csv, complete through depth 6]
        # Single source of truth: `Z7_DELETED_DIGIT` in src/z7.jl (the table is
        # a digit-alphabet fact, not spherical geometry) — see PROVENANCE.md.
        @test IGeo7.Z7_DELETED_DIGIT == (2, 2, 2, 2, 2, 2, 5, 5, 5, 5, 5, 5)
        for b in 0:11
            @test IGeo7.z7_deleted_digit(b) == (b <= 5 ? 2 : 5)
        end
        @test_throws InvalidZ7Error IGeo7.z7_deleted_digit(12)
        @test_throws InvalidZ7Error IGeo7.z7_deleted_digit(-1)
    end

    @testset "cross-file constant consistency" begin
        # Numeric facts deliberately restated in two layers: `NBASE` is shared
        # geometry (`ISEA`), `Z7_NUM_BASES` is the pure integer Z7 layer, and
        # neither may depend on the other. Pinned equal here so the duplication
        # cannot drift — see PROVENANCE.md.
        @test ISEA.NBASE == IGeo7.Z7_NUM_BASES
        @test IGeo7.MAX_DIGITS == IGeo7.Z7_MAX_RESOLUTION
        @test IGeo7.MAX_RESOLUTION + 1 == IGeo7.Z7_MAX_RESOLUTION
    end

    @testset "nearest vertex" begin
        for b in 0:11
            Vb = ISEA.vertex(b)
            @test ISEA.nearest_vertex(Vb) == b
            @test ISEA.nearest_vertex(ISEA.xyz_to_lonlat(Vb)...) == b
            # a point most of the way from the vertex toward each neighbor
            # still belongs to the vertex
            for n in ISEA.NEIGHBORS[b+1]
                p = ISEA.vnormalize(ISEA.vadd(ISEA.vscale(Vb, 0.8),
                    ISEA.vscale(ISEA.vertex(n), 0.2)))
                @test ISEA.nearest_vertex(p) == b
            end
        end
    end

    @testset "orientation" begin
        O = ISEA.ORIENT_IDENTITY
        @test O.identity
        @test O.R == (1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
        p = ISEA.lonlat_to_xyz(41.0, -7.5)
        @test ISEA.to_grid(O, p) === p
        @test ISEA.from_grid(O, p) === p
        # a non-identity rotation: 90 degrees about +z, world -> grid
        Rz = ISEA.Orientation((0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0))
        @test !Rz.identity
        q = ISEA.to_grid(Rz, p)
        @test all(abs.(q .- ISEA.lonlat_to_xyz(41.0 - 90.0, -7.5)) .< 1e-15)
        @test all(abs.(ISEA.from_grid(Rz, q) .- p) .< 1e-15)
        # the constructor detects the identity matrix
        @test ISEA.Orientation((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)).identity

        # ...and it is the *only* way in. `identity` is a cached `R == I`
        # predicate, not a degree of freedom: `to_grid`/`from_grid` skip the
        # rotation entirely when it is set, so `Orientation(R, true)` beside a
        # non-identity `R` would return every point unrotated — a silently
        # misplaced grid, not an error. The two-argument form is gone.
        @test_throws MethodError ISEA.Orientation(Rz.R, true)
        @test_throws MethodError ISEA.Orientation(ISEA.ORIENT_IDENTITY.R, false)
        @test ISEA.Orientation(Rz.R).identity == Rz.identity == false
        # The default entry points still see the identity, so the res-0 lock
        # below (and every geometric default) is unmoved.
        @test ISEA.ORIENT_IDENTITY.identity
        @test ISEA.ORIENT_IDENTITY.R == ISEA._IDENTITY_R
    end

    # -----------------------------------------------------------------------
    # Public-API res-0 lock: the identity orientation IS the standard ISEA
    # placement, so cell_center / lonlat_to_z7 reproduce the dggrid CLI's
    # res-0 rows directly (tolerance = the CSV's own noise, see above).
    # -----------------------------------------------------------------------
    @testset "public API res-0 lock" begin
        rows = _load_dggrid_centers(joinpath(_ICO_VEC, "dggrid_res0_centers.csv"))

        @testset "cell_center" begin
            worst = 0.0
            for (b, lon, lat) in rows
                z = z7_from_string(lpad(b, 2, '0'))
                @test z7_resolution(z) == 0
                @test z7_base_cell(z) == b
                @test z7_is_pentagon(z)
                clon, clat = cell_center(z)
                d = ISEA.angdist(ISEA.lonlat_to_xyz(clon, clat),
                    ISEA.lonlat_to_xyz(lon, lat))
                @test d < 1e-8
                worst = max(worst, d)
            end
            @test worst < 6e-9
            # base 0 lands exactly on the published closed form
            @test cell_center(z7_from_string("00"))[1] ≈ 11.25 atol = 1e-12
            @test cell_center(z7_from_string("00"))[2] ≈
                  ISEA.ISEA_LAT_HI atol = 1e-12
        end

        @testset "pentagon centers stay pentagons" begin
            # A base cell's published center is the icosahedron vertex, so at
            # EVERY resolution it must decode to that base's all-0-digit
            # pentagon. This exercises decode + collapse + chart well past
            # the res-0 data.
            for (b, lon, lat) in rows, r in (0, 1, 2, 3, 5, 8)
                z = lonlat_to_z7(lon, lat, r)
                @test z == z7_from_string(lpad(b, 2, '0') * "0"^r)
                @test z7_is_pentagon(z)
                @test z7_base_cell(z) == b
                @test z7_resolution(z) == r
            end
            # lonlat_to_cell agrees with lonlat_to_z7
            for (b, lon, lat) in rows
                @test lonlat_to_cell(lon, lat, 5) == lonlat_to_z7(lon, lat, 5)
            end
            # ... and the encode side round-trips back: the res-5 pentagon's
            # own center decodes to itself
            for (b, _, _) in rows
                z = z7_from_string(lpad(b, 2, '0') * "00000")
                clon, clat = cell_center(z)
                @test lonlat_to_z7(clon, clat, 5) == z
            end
        end
    end

    @testset "allocation-free scalar paths" begin
        p = ISEA.lonlat_to_xyz(12.0, 34.0)
        ISEA.nearest_vertex(p)
        ISEA.lonlat_to_xyz(12.0, 34.0)
        ISEA.xyz_to_lonlat(p)
        @test @allocated(ISEA.nearest_vertex(p)) == 0
        @test @allocated(ISEA.lonlat_to_xyz(12.0, 34.0)) == 0
        @test @allocated(ISEA.xyz_to_lonlat(p)) == 0
    end

end
