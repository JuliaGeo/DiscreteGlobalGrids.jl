# test/test_chart.jl — unit tests for src/ISEA/snyder.jl
# (the per-face Snyder ISEA plane and the dev-frame slot maps).
#
# Everything here is derived from analytic properties of the published
# construction (spec/isea-projection-spec.md §3, §5–§7) plus that spec's own
# §9 forward/inverse vectors — the chart is self-validated against published
# numbers. The composed-grid behavior (encode/decode against the dggrid CLI
# dumps) lives in test_grid.jl.

module TestChart

using Test
using Random

# The chart is shared ISEA machinery: `snyder.jl` plus the `icosahedron.jl`
# vertex tables and vector helpers it builds on (design.md Section 2).
using DiscreteGlobalGrids.ISEA

const RNG = MersenneTwister(20260805)

@testset verbose = true "snyder.jl — per-face Snyder ISEA plane" begin

    # -------------------------------------------------------------------
    @testset "1. constants and closed forms" begin
        @test R_EA ≈ sqrt(4pi / (15 * sqrt(3.0))) rtol = 0
        @test L_PLANE ≈ sqrt(3.0) * R_EA rtol = 0
        @test R_EA2 == R_EA * R_EA
        @test COS_LG ≈ cotd(60.0) * cotd(36.0) atol = 1e-16
        @test TAN_LG ≈ 3 - sqrt(5.0) atol = 1e-15          # tan g
        @test acosd(COS_LG) ≈ 37.37736814064969 atol = 1e-12
        @test SNY_COTT == sqrt(3.0)                        # cot 30
        @test SNY_SECTOR ≈ 2pi / 3 rtol = 0
        @test DEV_CONE_DEG == 300.0
        # the equal-area face: 20 triangles of area 4pi/20, circumradius R_EA,
        # planar area (3*sqrt(3)/4)*R_EA^2 per triangle
        @test 20 * (3 * sqrt(3.0) / 4) * R_EA2 ≈ 4pi rtol = 1e-15
        # g is the vertex -> face-center arc of the icosahedron
        f = FACES[1]
        @test acosd(COS_LG) ≈ angdist(VERTICES[f.verts[1]+1], f.c) atol = 1e-12
    end

    # -------------------------------------------------------------------
    @testset "2. faces: triples, frames, corners" begin
        @test length(FACES) == 20
        # every base is a corner of exactly 5 faces; triples are ascending
        for b in 0:11
            @test count(t -> b in t, FACE_TRIPLES) == 5
        end
        for (fi, tri) in enumerate(FACE_TRIPLES)
            @test issorted(collect(tri))
            f = FACES[fi]
            @test f.verts == tri
            # orthonormal right-handed frame at the face center
            @test abs(vnorm(f.u) - 1) < 1e-15
            @test abs(vnorm(f.w) - 1) < 1e-15
            @test abs(vdot(f.u, f.c)) < 1e-15
            @test abs(vdot(f.w, f.c)) < 1e-15
            @test all(abs.(vcross(f.c, f.u) .- f.w) .< 1e-14)
            # corners at radius R_EA on the rays 0/±120°, aligned with verts
            for (i, v) in enumerate(tri)
                @test abs(abs(f.corner[i]) - R_EA) < 1e-15
                @test abs(angdist(VERTICES[v+1], f.c) - acosd(COS_LG)) < 1e-12
            end
            @test abs(f.corner[findfirst(==(minimum(tri)), tri)] -
                      complex(R_EA, 0.0)) < 1e-15
        end
    end

    # -------------------------------------------------------------------
    @testset "3. spec §9 vectors: T6 forward" begin
        # (lon, lat, f, x, y) [spec/isea-projection-spec.md 9.2]
        T6 = [
            (0.0, 0.0, 6, -0.34773547074696676, 0.20898643445399323),
            (11.25, 0.0, 6, -0.34773547074696670, -5.8e-16),
            (90.0, 0.0, 8, -0.16130676939722316, -2.7e-16),
            (-90.0, 0.0, 10, 0.08065338469861187, -0.13969576010039336),
            (180.0, 0.0, 7, -0.34773547074696676, 0.20898643445399323),
            (0.0, 90.0, 0, 0.17386773537348396, -0.30114775146381334),
            (11.25, 58.2825255885389, 6, 0.69547094149393129, -3.1e-16),
            (101.25, 69.0948425521107, 0, 5e-17, -2e-17),
            (56.25, 35.2643896827547, 2, 5.2e-16, -5.7e-16),
            (-77.0, 38.9, 1, -0.25992037346515245, -0.49557483299448002),
            (139.7, 35.7, 3, -0.05195418494398664, 0.08639975918728288),
            (2.35, 48.85, 6, 0.50877992213758572, 0.09725187947879274),
            (151.2, -33.87, 15, -0.07110028496677201, -0.01323046773125910),
            (-58.4, -34.6, 17, -0.17591654329840159, 0.33267017266645826),
            (37.6, 55.75, 2, 0.42847866031377657, -0.07800805690131862),
            (-157.86, 21.31, 7, 0.00748063799522557, -0.17424567915016106),
            (18.42, -33.92, 12, -0.20130655269856126, 0.15398866257092678),
            (-43.2, -22.9, 17, 0.13031000184752703, 0.21068632492196290),
            (103.85, 1.29, 9, 0.17305878492015556, 0.25215196404683909),
        ]
        worst_xy = 0.0
        nface = 0
        for (lon, lat, f, x, y) in T6
            (fg, w) = snyder_fwd(lonlat_to_xyz(lon, lat))
            fg == f || (nface += 1)
            worst_xy = max(worst_xy, abs(w - complex(x, y)))
        end
        @test nface == 0
        @test worst_xy < 1e-13
    end

    # -------------------------------------------------------------------
    @testset "4. spec §9 vectors: T7 inverse + invariants" begin
        # (f, x, y, lon, lat) [spec 9.3]
        T7 = [
            (0, 0.0, 0.0, 101.2500000000000142, 69.0948425521106913),
            (6, 0.3, 0.1, 3.5525584486344801, 37.5079522022028300),
            (6, -0.3, -0.2, 22.0644755061567679, 2.9262858193115817),
            (8, 0.69547094149393307, 0.0, 42.9674744114610334, 0.0),
            (8, -0.34773547074696654, 0.60229550292762746, 101.25, -31.7174744114610050),
        ]
        worst_inv = 0.0
        for (f, x, y, lon, lat) in T7
            p = snyder_inv_xyz(f, complex(x, y))
            worst_inv = max(worst_inv, angdist(p, lonlat_to_xyz(lon, lat)))
        end
        @test worst_inv < 1e-11

        # invariants: face centers map to 0, vertices to radius R_EA
        for f in 0:19
            (fg, w) = snyder_fwd(FACES[f+1].c)
            @test fg == f
            @test abs(w) < 1e-13
        end
        worst_v = 0.0
        for f in 0:19, v in FACES[f+1].verts
            (fg, w) = snyder_fwd(VERTICES[v+1])
            fg == f || continue                      # lowest-index tie rule
            worst_v = max(worst_v, abs(abs(w) - R_EA))
        end
        @test worst_v < 1e-12
    end

    # -------------------------------------------------------------------
    @testset "5. forward / inverse round trip" begin
        worst_rt = 0.0
        for lon in -175.0:5.0:175.0, lat in -85.0:4.25:85.0
            (f, w) = snyder_fwd(lonlat_to_xyz(lon, lat))
            p = snyder_inv_xyz(f, w)
            worst_rt = max(worst_rt, angdist(p, lonlat_to_xyz(lon, lat)))
        end
        # random points too
        for _ in 1:20000
            lon = 360.0 * rand(RNG) - 180.0
            lat = asind(2.0 * rand(RNG) - 1.0)
            (f, w) = snyder_fwd(lonlat_to_xyz(lon, lat))
            p = snyder_inv_xyz(f, w)
            worst_rt = max(worst_rt, angdist(p, lonlat_to_xyz(lon, lat)))
        end
        @test worst_rt < 1e-11
    end

    # -------------------------------------------------------------------
    @testset "6. equal area: Jacobian is identically 1" begin
        # d(x,y)/d(lon_rad, lat_rad) must equal cos(lat) (the spherical area
        # element) at interior face points. Finite differences.
        h = 1e-5
        worst = 0.0
        ntested = 0
        for lon in -170.0:11.0:170.0, lat in -80.0:8.0:80.0
            # skip points whose 4-point stencil straddles a face boundary
            fs = [snyder_fwd(lonlat_to_xyz(lon + dl, lat + dt))[1]
                  for (dl, dt) in ((-rad2deg(h), 0.0), (rad2deg(h), 0.0),
                (0.0, -rad2deg(h)), (0.0, rad2deg(h)))]
            allequal(fs) || continue
            ux = (snyder_fwd(lonlat_to_xyz(lon + rad2deg(h), lat))[2] -
                  snyder_fwd(lonlat_to_xyz(lon - rad2deg(h), lat))[2]) / (2h)
            uy = (snyder_fwd(lonlat_to_xyz(lon, lat + rad2deg(h)))[2] -
                  snyder_fwd(lonlat_to_xyz(lon, lat - rad2deg(h)))[2]) / (2h)
            det = real(ux) * imag(uy) - imag(ux) * real(uy)
            worst = max(worst, abs(det / cosd(lat) - 1))
            ntested += 1
        end
        @test ntested > 400
        @test worst < 1e-6
    end

    # -------------------------------------------------------------------
    @testset "7. dev-frame slot maps" begin
        # rigid, invertible, consistent at slot boundaries
        nbad = 0
        for b in 0:11, j in 0:4, _ in 1:20
            u = (0.9 * L_PLANE * rand(RNG)) * cis(deg2rad(60.0 * j + 60.0 * rand(RNG)))
            s = DEV_SLOTS[b+1][j+1]
            w = s.cb + s.rot * u
            abs(face_to_dev(b, s.f, w) - u) < 1e-13 || (nbad += 1)
            DEV_SLOT_OF_FACE[b+1][s.f+1] == j || (nbad += 1)
            dev_slot_index(u) == j || (nbad += 1)
        end
        @test nbad == 0
        # the slot-0 ray points at the fitted reference neighbor: dev position
        # on the ray at distance L maps to the neighbor's vertex
        for b in 0:11
            s = DEV_SLOTS[b+1][1]
            w = s.cb + s.rot * complex(L_PLANE, 0.0)
            p = snyder_inv_xyz(s.f, w)
            @test angdist(p, VERTICES[REFERENCE_EDGE[b+1]+1]) < 1e-12
        end
        # dev origin is the base vertex, through every slot
        for b in 0:11, j in 0:4
            s = DEV_SLOTS[b+1][j+1]
            @test angdist(snyder_inv_xyz(s.f, s.cb), VERTICES[b+1]) < 1e-12
        end
        # dev_to_xyz round-trips through the containing face
        nbad2 = 0
        for b in 0:11, _ in 1:50
            u = (0.95 * L_PLANE * rand(RNG)) * cis(deg2rad(299.9 * rand(RNG)))
            p = dev_to_xyz(b, u)
            f, w = snyder_fwd(p)
            # p may sit on a face this base does not corner (fringe); only
            # check the inversion when it does
            DEV_SLOT_OF_FACE[b+1][f+1] < 0 && continue
            abs(face_to_dev(b, f, w) - u) < 1e-12 || (nbad2 += 1)
        end
        @test nbad2 == 0
        # the cut guard band: an FP shadow of dev angle 0 reads as slot 0
        @test dev_slot_index(complex(0.4, -1e-12)) == 0
        @test dev_slot_index(0.4 * cis(deg2rad(299.999))) == 4
        @test dev_angle_deg(complex(0.4, -1e-12)) < 0.0
        @test dev_angle_deg(0.4 * cis(deg2rad(299.999))) ≈ 299.999 atol = 1e-9
    end

    # -------------------------------------------------------------------
    @testset "8. allocation-free scalar paths" begin
        p = lonlat_to_xyz(3.0, 40.0)
        u = 0.3 * cis(0.7)
        function alloc_probe(p, u)
            a1 = @allocated snyder_fwd(p)
            a2 = @allocated snyder_inv_xyz(4, u)
            a3 = @allocated dev_to_xyz(2, u)
            a4 = @allocated face_to_dev(0, 1, u)
            a5 = @allocated dev_slot_index(u)
            a6 = @allocated dev_angle_deg(u)
            return (a1, a2, a3, a4, a5, a6)
        end
        alloc_probe(p, u)          # compile
        @test alloc_probe(p, u) == (0, 0, 0, 0, 0, 0)
    end
end

end # module
