# test/test_grid.jl — the primary geometry suite for src/IGeo7/grid.jl:
# encode (cell_center / cell_boundary / cell_area) and decode (lonlat_to_z7)
# of the single grid law (spec/igeo7-geometry-diagnosis.md §6–§7,
# spec/design.md section 12).
#
# Oracle suites (dggrid CLI black-box facts):
#   * dggrid_true_res{1..5}_centers.txt — ALL 196,080 dggrid cells: exact z7
#     from center decode, centers reproduced within 1e-8 deg (the data's own
#     print noise is <= 5.913e-9 deg; res 5 pins the alternating-chirality
#     continuation, diagnosis §8.1)
#   * dggrid_res0_centers.csv — res-0 lock (also in test_icosahedron.jl)
#
# No dggrid boundary dumps exist (diagnosis §8.3), so the boundary and
# round-trip behavior at res 0:19 is property-tested: center-in-own-boundary,
# ring closure, CCW winding, closed-form areas, pentagon counts, orientation
# threading, adversarial points, allocation-free scalar paths.

module TestGrid

using Test
using Random
using Printf

# `grid.jl` composes the shared ISEA geometry with the IGEO7 Z7/engine layers
# (design.md Section 2), so this suite needs both namespaces.
using DiscreteGlobalGrids.ISEA
using DiscreteGlobalGrids.IGeo7
# not exported by IGeo7 (an internal digit-alphabet accessor), but the random
# id generator and the cut/rim regression families below need it
using DiscreteGlobalGrids.IGeo7: z7_deleted_digit

const VEC_DIR = normpath(joinpath(@__DIR__, "vectors"))

const RNG = MersenneTwister(20260805)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

"angular distance in degrees between two (lon, lat) pairs"
lldist(a::NTuple{2,Float64}, b::NTuple{2,Float64}) =
    angdist(lonlat_to_xyz(a[1], a[2]), lonlat_to_xyz(b[1], b[2]))

"random valid Z7 id at resolution `res` (deleted pentagon digits skipped)"
function random_valid_id(rng, res::Int)
    base = rand(rng, 0:11)
    z = (UInt64(base) << 60) | UInt64(0x0fffffffffffffff)   # res-0 pentagon
    deleted = z7_deleted_digit(base)
    pentagon = true
    for _ in 1:res
        d = rand(rng, 0:6)
        while pentagon && d == deleted
            d = rand(rng, 0:6)
        end
        z = z7_child(z, d)
        pentagon &= (d == 0)
    end
    return z
end

# Both star tests subtract the reference point before the cross product:
# p·(a×b) == p·((a−p)×(b−p)) exactly, but the subtracted form evaluates it
# without the O(1) cancellation — at res 19 the signal is ~4e-17 while the
# unsubtracted form carries ~1e-16 of FP noise and flips signs at random.

"strict spherical point-in-ring: all geodesic edge planes on one side"
function inside_corner_ring(ring::Vector{NTuple{3,Float64}}, p::NTuple{3,Float64})
    smin, smax = Inf, -Inf
    n = length(ring)
    for i in 1:n
        a = vsub(ring[i], p)
        b = vsub(ring[mod1(i + 1, n)], p)
        s = vdot(vcross(a, b), p)
        smin = min(smin, s)
        smax = max(smax, s)
    end
    return smin > 0.0 || smax < 0.0
end

"ring winding seen from outside: +1 CCW, -1 CW, 0 mixed (star test from `c`)"
function ring_winding(ring::Vector{NTuple{3,Float64}}, c::NTuple{3,Float64})
    n = length(ring)
    pos = 0
    neg = 0
    for i in 1:n
        a = vsub(ring[i], c)
        b = vsub(ring[mod1(i + 1, n)], c)
        s = vdot(vcross(a, b), c)
        s > 0 ? (pos += 1) : (neg += 1)
    end
    pos == n && return 1
    neg == n && return -1
    return 0
end

"rows (z7, lon, lat) of dggrid_true_res{r}_centers.txt"
function load_true_centers(r::Int)
    rows = Tuple{UInt64,Float64,Float64}[]
    for line in eachline(joinpath(VEC_DIR, "dggrid_true_res$(r)_centers.txt"))
        s = strip(line)
        isempty(s) && continue
        p = split(s, ',')
        push!(rows, (z7_from_string(String(p[1])), parse(Float64, p[2]),
            parse(Float64, p[3])))
    end
    return rows
end

# ===========================================================================

@testset verbose = true "grid.jl — encode/decode of the composed grid" begin

    # -------------------------------------------------------------------
    @testset "1. dggrid dumps res 1..5: exhaustive encode + decode" begin
        # ALL 196,080 cells of the actual dggrid binary; res 5 pins the
        # alternating-chirality continuation (diagnosis §8.1).
        total = 0
        for r in 1:5
            rows = load_true_centers(r)
            total += length(rows)
            worst = 0.0
            ndecode = 0
            for (z, lon, lat) in rows
                c = cell_center(z)
                worst = max(worst, lldist(c, (lon, lat)))
                lonlat_to_z7(lon, lat, r) == z || (ndecode += 1)
            end
            @printf("    res %d: %6d cells   center max err %.3e deg   decode mismatches %d\n",
                r, length(rows), worst, ndecode)
            @test worst <= 1e-8
            @test ndecode == 0
        end
        @test total == 196080
    end

    # -------------------------------------------------------------------
    @testset "2. round-trip z7 -> center -> decode at every res" begin
        nbad = 0
        ntried = 0
        for res in 0:19, _ in 1:60
            z = random_valid_id(RNG, res)
            lon, lat = cell_center(z)
            ntried += 1
            lonlat_to_z7(lon, lat, res) == z || (nbad += 1)
        end
        @test ntried == 1200
        @test nbad == 0
    end

    # -------------------------------------------------------------------
    @testset "3. boundary structure: counts, closure, winding, containment" begin
        nbad_ring = 0
        nbad_wind = 0
        for res in 0:19, _ in 1:12
            z = random_valid_id(RNG, res)
            c = lonlat_to_xyz(cell_center(z)...)
            ring = cell_boundary_cartesian(z; closed_ring=false)
            inside_corner_ring(ring, c) || (nbad_ring += 1)
            # boundary winding is CCW seen from outside under the identity
            # orientation (diagnosis §6)
            ring_winding(ring, c) == 1 || (nbad_wind += 1)
        end
        @test nbad_ring == 0
        @test nbad_wind == 0
        for res in (0, 1, 3, 7, 19)
            zp = z7_from_string("07" * "0"^res)              # pentagon
            @test length(cell_boundary(zp; closed_ring=false)) == 5
            ringp = cell_boundary(zp)                        # default closed
            @test length(ringp) == 6
            @test ringp[1] == ringp[end]
            if res >= 1
                zh = z7_from_string("07" * "1" * "0"^(res - 1))  # hexagon
                @test length(cell_boundary(zh; closed_ring=false)) == 6
                ringh = cell_boundary(zh; closed_ring=true)
                @test length(ringh) == 7
                @test ringh[1] == ringh[end]
            end
            cart = cell_boundary_cartesian(zp; closed_ring=false)
            @test length(cart) == 5
            @test all(abs(vnorm(p) - 1) < 1e-14 for p in cart)
            cartc = cell_boundary_cartesian(zp; closed_ring=true)
            @test length(cartc) == 6 && cartc[1] == cartc[end]
            # lonlat and cartesian describe the same ring
            @test all(angdist(lonlat_to_xyz(ll[1], ll[2]), p) < 1e-12
                      for (ll, p) in zip(cell_boundary(zp; closed_ring=false), cart))
        end
    end

    # -------------------------------------------------------------------
    @testset "4. pentagons: counts, vertices, res-0 rings" begin
        # the 12 pentagons per res sit exactly at the icosahedron vertices
        for res in (0, 1, 5, 19), b in 0:11
            z = z7_from_string(string(b, pad=2) * "0"^res)
            @test is_pentagon(z)
            @test angdist(lonlat_to_xyz(cell_center(z)...), vertex(b)) < 1e-12
        end
        # res-0 rings are the five adjacent Snyder face centers
        for b in 0:11
            z = (UInt64(b) << 60) | UInt64(0x0fffffffffffffff)
            ring = cell_boundary_cartesian(z; closed_ring=false)
            fcs = [FACES[f+1].c for f in 0:19 if b in FACES[f+1].verts]
            @test length(ring) == 5 && length(fcs) == 5
            w = maximum(minimum(angdist(v, q) for q in ring) for v in fcs)
            @test w < 1e-12
        end
    end

    # -------------------------------------------------------------------
    @testset "5. areas: closed form" begin
        # analytic closed form [contract]: hexagons 4*pi*R^2/(10*7^r),
        # pentagons exactly 5/6 of that; 12 pentagons + (10*7^r - 10)
        # hexagons tile the authalic sphere. The Snyder chart is exactly
        # equal-area, so the closed form IS the geometry.
        for res in 0:19
            zp = z7_from_string("03" * "0"^res)
            hex_a = 4pi * R_AUTHALIC^2 / (10 * 7.0^res)
            @test cell_area(zp) == 5 * hex_a / 6
            if res >= 1
                zh = z7_from_string("03" * "1" * "0"^(res - 1))
                @test cell_area(zh) == hex_a
                @test isapprox(12 * cell_area(zp) + (10 * 7.0^res - 10) * cell_area(zh),
                    4pi * R_AUTHALIC^2; rtol=1e-12)
            else
                @test isapprox(12 * cell_area(zp), 4pi * R_AUTHALIC^2; rtol=1e-12)
            end
        end
    end

    # -------------------------------------------------------------------
    @testset "6. orientation parameter" begin
        # a non-trivial rigid rotation (90 deg about x then 30 deg about z)
        c30, s30 = cosd(30.0), sind(30.0)
        Rz = ((c30, -s30, 0.0), (s30, c30, 0.0), (0.0, 0.0, 1.0))
        Rx = ((1.0, 0.0, 0.0), (0.0, 0.0, -1.0), (0.0, 1.0, 0.0))
        # R = Rz * Rx, row-major
        mat = ntuple(9) do i
            r, c = divrem(i - 1, 3)
            sum(Rz[r+1][k] * Rx[k][c+1] for k in 1:3)
        end
        O = Orientation(mat)
        @test !O.identity
        nbad = 0
        for res in (0, 1, 4, 9, 15, 19), _ in 1:8
            z = random_valid_id(RNG, res)
            # encode consistency: rotated center == R' * (grid-frame center)
            cg = lonlat_to_xyz(cell_center(z)...)
            co = cell_center(z; orientation=O)
            angdist(lonlat_to_xyz(co[1], co[2]), from_grid(O, cg)) < 1e-12 || (nbad += 1)
            # decode round-trip in the rotated frame
            lonlat_to_z7(co[1], co[2], res; orientation=O) == z || (nbad += 1)
            # boundary rotates with the frame
            bg = cell_boundary_cartesian(z; closed_ring=false)
            bo = cell_boundary_cartesian(z; closed_ring=false, orientation=O)
            all(angdist(from_grid(O, p), q) < 1e-12 for (p, q) in zip(bg, bo)) || (nbad += 1)
        end
        @test nbad == 0
        # aliases and composition
        @test lonlat_to_cell(12.3, 45.6, 7) == lonlat_to_z7(12.3, 45.6, 7)
        @test lonlat_to_index(12.3, 45.6, 7) ==
              cell_to_index(lonlat_to_z7(12.3, 45.6, 7))
    end

    # -------------------------------------------------------------------
    @testset "7. adversarial points: vertices, edge midpoints, meridians" begin
        # every point must decode and round-trip strictly (decode -> center
        # -> decode is a fixed point)
        nbad = 0
        for res in (0, 1, 5, 10, 19)
            pts = NTuple{2,Float64}[]
            for b in 0:11
                push!(pts, xyz_to_lonlat(vertex(b)))
                for n in NEIGHBORS[b+1]
                    push!(pts, xyz_to_lonlat(vnormalize(vadd(vertex(b), vertex(n)))))
                end
            end
            for lon in (0.0, 11.25, -168.75, 36.0, 180.0), lat in -85.0:17.0:85.0
                push!(pts, (lon, lat))
            end
            for (lon, lat) in pts
                z = lonlat_to_z7(lon, lat, res)
                c = cell_center(z)
                lonlat_to_z7(c[1], c[2], res) == z || (nbad += 1)
            end
        end
        @test nbad == 0
        # cut/rim regression families: cells whose subtree hugs the cone cut
        # at every res (design Section 7 guard bands) — every digit pattern
        # d,0,0,... and d,5,0,5,0,... for every base
        nbad2 = 0
        for base in 0:11, d in 1:6, res in (1, 2, 5, 10, 19)
            d == z7_deleted_digit(base) && continue
            z = z7_from_string(string(base, pad=2) * string(d) * "0"^(res - 1))
            lonlat_to_z7(cell_center(z)..., res) == z || (nbad2 += 1)
            # after the first nonzero digit the deleted digit is legal again,
            # so the "d5050..." family is valid for every base
            z2 = z7_from_string(string(base, pad=2) * string(d) * repeat("50", 10)[1:res-1])
            lonlat_to_z7(cell_center(z2)..., res) == z2 || (nbad2 += 1)
        end
        @test nbad2 == 0
    end

    # -------------------------------------------------------------------
    @testset "8. argument validation" begin
        z20 = random_valid_id(RNG, 20)
        @test z7_resolution(z20) == 20
        @test_throws InvalidZ7Error cell_center(z20)
        @test_throws InvalidZ7Error cell_boundary(z20)
        @test_throws InvalidZ7Error cell_boundary_cartesian(z20)
        @test_throws InvalidZ7Error cell_area(z20)
        @test_throws InvalidZ7Error lonlat_to_z7(0.0, 0.0, 20)
        @test_throws InvalidZ7Error lonlat_to_z7(0.0, 0.0, -1)
        @test_throws InvalidZ7Error lonlat_to_cell(0.0, 0.0, 20)
        # structurally invalid ids: bad base, deleted-pentagon digit chain
        @test_throws InvalidZ7Error cell_center(0xffffffffffffffff)
        bad = z7_from_string("013") ⊻ (UInt64(1) << 57)   # digit 3 -> 2 (deleted)
        @test !is_valid_z7(bad)
        @test_throws InvalidZ7Error cell_center(bad)
        @test_throws InvalidZ7Error cell_boundary(bad)
        @test_throws InvalidZ7Error cell_area(bad)
    end

    # -------------------------------------------------------------------
    @testset "9. allocation-free, type-stable scalar paths" begin
        z = z7_from_string("0803425160")
        @test (@inferred cell_center(z)) isa NTuple{2,Float64}
        @test (@inferred lonlat_to_z7(12.3, 45.6, 10)) isa UInt64
        @test (@inferred cell_area(z)) isa Float64
        f_dec(lon, lat, r) = @allocated lonlat_to_z7(lon, lat, r)
        f_cen(zz) = @allocated cell_center(zz)
        f_area(zz) = @allocated cell_area(zz)
        f_dec(12.3, 45.6, 10); f_cen(z); f_area(z)      # warmup
        @test f_dec(12.3, 45.6, 10) == 0
        @test f_dec(-91.0, -37.0, 19) == 0
        @test f_cen(z) == 0
        @test f_cen(random_valid_id(RNG, 19)) == 0
        @test f_area(z) == 0
        # a fringe-family decode and a near-vertex decode stay allocation-free
        zf = z7_from_string("013")
        lonf, latf = cell_center(zf)
        @test f_dec(lonf, latf, 1) == 0
        lonv, latv = xyz_to_lonlat(vertex(4))
        @test f_dec(lonv, latv, 10) == 0
        # quick timing report (informational; the batched-median perf
        # methodology and recorded numbers live in spec/design.md §6 / §11)
        function bench(f, n)
            t0 = time_ns()
            s = 0.0
            for _ in 1:n
                s += f()
            end
            return (time_ns() - t0) / n, s
        end
        lon0, lat0 = 12.3, 45.6
        tdec, _ = bench(() -> Float64(lonlat_to_z7(lon0, lat0, 10) % 7), 100_000)
        tcen, _ = bench(() -> cell_center(z)[1], 100_000)
        @printf("    lonlat_to_z7 res 10: %.0f ns/call   cell_center res 10: %.0f ns/call\n",
            tdec, tcen)
    end
end

end # module TestGrid
