# Tests for the ellipsoid wrapper in `src/fallbacks/authalic_grid.jl`:
#
#   1. Bounds. `authalic_shift` and `authalic_stretch` are derived
#      analytically from the series coefficients, so they are checked against a
#      dense sweep of the transform they claim to bound — and against the
#      Lipschitz property itself, measured on point pairs rather than restated
#      from the derivation.
#   2. FORWARDING. Everything that is not geometry must come through the wrapper
#      unchanged, ids and indices above all.
#   3. Covering law. The conformance harness is the oracle, but the
#      necessity of the inflation is shown separately: against the TIGHTEST
#      sound base cap the warp really does push descendants out, so a wrapper
#      that forwarded `node_extent` would under-cover.
#   4. REGISTRATION. A wrapped grid regridded
#      against a geodetic lon/lat grid lands where the data says it does, and an
#      unwrapped one misses by the authalic shift.
#
# Direct collector calls pin the warped geometry in addition to the conformance
# harness's specialized-method checks.

module AuthalicWrapperTests

using Test
using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
const FB = DGG.Fallbacks
const H = DGG.Helpers
import GeometryOps as GO
import GeometryOpsCore as GOCore
import ConservativeRegridding as CR
import GeometryOps: SpatialTreeInterface as STI
using DiscreteGlobalGridsConformanceTesting

const US = GO.UnitSpherical
const CT = DiscreteGlobalGridsConformanceTesting
const WGS84 = H.WGS84_AUTHALIC
const BASE = HEALPixSystem()
const SYS = AuthalicSystem(BASE)

# A stub whose node extent becomes non-convex if inflated. The test only uses
# its `node_extent` method.
struct WideExtentStub <: DGG.AbstractHierarchicalGridSystem
    radius::Float64
end

DGG.node_extent(sys::WideExtentStub, ::DGG.AbstractCellIndex) =
    US.SphericalCap(FB.USPoint(0.0, 0.0, 1.0), sys.radius)

# A deterministic spread of cells: no RNG, so a failure names the same cell on
# every run and on every machine.
spread(grid, n::Int) = [cellindex(grid, i) for i in 1:max(1, ncells(grid) ÷ n):ncells(grid)]

@testset "authalic wrapper" begin

# ===========================================================================
# 1. The bounds
# ===========================================================================

@testset "authalic_shift / authalic_stretch bound the transform" begin
    shift = FB.authalic_shift(WGS84)
    stretch = FB.authalic_stretch(WGS84)

    # The familiar WGS84 figures: 0.1284 deg of latitude, 0.45% of stretch.
    @test rad2deg(shift) ≈ 0.128463 atol = 1e-6
    @test stretch - 1 ≈ 4.48998e-3 rtol = 1e-5

    # Swept maxima of the two quantities the derivation bounds: the latitude
    # displacement, and both singular values of dPhi = diag(cos φ/cos ξ, φ').
    worst_shift = 0.0
    worst_north = 0.0
    worst_east = 0.0
    steps = 20_000
    for k in 0:steps
        ξ = -pi / 2 + pi * k / steps
        φ = H.authalic_to_geodetic(WGS84, ξ)
        worst_shift = max(worst_shift, abs(φ - ξ))
        h = 1e-6
        worst_north = max(worst_north,
            (H.authalic_to_geodetic(WGS84, ξ + h) - H.authalic_to_geodetic(WGS84, ξ - h)) / 2h)
        abs(cos(ξ)) > 1e-6 && (worst_east = max(worst_east, cos(φ) / cos(ξ)))
    end
    @test worst_shift <= shift
    @test worst_shift > 0.99 * shift         # and the bound is not slack
    @test worst_north <= stretch
    @test worst_north > 0.999 * stretch      # attained at the equator, where cos(2jξ) = 1
    @test worst_east <= stretch

    # The Lipschitz property itself, measured: no pair of points may have its
    # separation stretched past the constant. Deterministic point set, densest
    # where the stretch is (meridional pairs about the equator).
    worst_ratio = 0.0
    for lat0 in -88:4:88, dlat in (0.05, 0.5, 5.0, 30.0), dlon in (0.0, 10.0, 90.0)
        lat1 = lat0 + dlat
        lat1 > 90 && continue
        p = FB.unit_point(0.0, Float64(lat0))
        q = FB.unit_point(dlon, Float64(lat1))
        d0 = US.spherical_distance(p, q)
        d1 = US.spherical_distance(FB.geodetic_point(WGS84, p), FB.geodetic_point(WGS84, q))
        worst_ratio = max(worst_ratio, d1 / d0)
    end
    @test worst_ratio <= stretch
    @test worst_ratio > 1                     # the warp really does stretch

    # The degenerate (spherical) transform is the identity, exactly.
    sphere = H.AuthalicTransform(GO.Spherical(; radius=1.0))
    @test FB.authalic_shift(sphere) == 0
    @test FB.authalic_stretch(sphere) == 1
end

@testset "the point warp" begin
    for lon in (-180.0, -75.0, 0.0, 12.5, 179.9), lat in (-90.0, -45.0, -0.0, 0.0, 23.5, 45.0, 90.0)
        p = FB.unit_point(lon, lat)
        w = FB.geodetic_point(WGS84, p)
        wlon, wlat = FB.lonlat(w)
        @test abs(sqrt(w[1]^2 + w[2]^2 + w[3]^2) - 1) < 1e-15
        # Longitude is untouched; latitude is the series, and the poles and the
        # equator are fixed points.
        abs(lat) == 90 || @test wlon ≈ lon atol = 1e-12
        @test wlat ≈ H.authalic_to_geodeticd(WGS84, lat) atol = 1e-12
        # ... and the input warp is its inverse.
        back = FB.authalic_point(WGS84, w)
        @test US.spherical_distance(back, p) < 1e-14
    end
    # The poles come back bit-identical rather than reconstructed.
    @test FB.geodetic_point(WGS84, FB.USPoint(0.0, 0.0, 1.0)) === FB.USPoint(0.0, 0.0, 1.0)
    # And so does everything, under the degenerate transform.
    sphere = H.AuthalicTransform(GO.Spherical(; radius=1.0))
    p = FB.unit_point(31.0, 47.0)
    @test FB.geodetic_point(sphere, p) === p
    @test FB.authalic_point(sphere, p) === p
end

# ===========================================================================
# 2. Forwarding
# ===========================================================================

@testset "ids, indices and hierarchy forward unchanged" begin
    for l in 0:4
        g = levelgrid(SYS, l)
        b = levelgrid(BASE, l)
        @test g isa AuthalicGrid
        @test ncells(g) == ncells(b)
        @test level(g) == l
        @test parent(g) === b
        for i in (1, 2, ncells(g) ÷ 3, ncells(g))
            c = cellindex(g, i)
            @test c === cellindex(b, i)
            @test globalindex(g, c) == i
        end
        # A cell that is not in the grid is a `nothing`, never an error.
        l < maxlevel(SYS) && @test globalindex(g, first(children(SYS, cellindex(g, 1)))) === nothing
        @test_throws BoundsError cellindex(g, 0)
        @test_throws BoundsError cellindex(g, ncells(g) + 1)
    end

    # `system(levelgrid(sys, l)) === sys` is what the cursor and the conformance
    # harness both assume, and it holds by value equality of an immutable pair.
    @test system(levelgrid(SYS, 3)) === SYS
    @test AuthalicSystem(BASE) === SYS
    @test AuthalicSystem(BASE, GO.Geodesic()) === SYS   # WGS84 by another name
    @test parent(SYS) === BASE
    @test system(AuthalicGrid(levelgrid(BASE, 2))) === SYS

    # The hierarchical surface is the base system's, verbatim.
    @test cellindextype(SYS) === cellindextype(BASE)
    @test cellindextypes(SYS) === cellindextypes(BASE)
    @test levels(SYS) == levels(BASE)
    @test maxlevel(SYS) == maxlevel(BASE)
    @test collect(rootcells(SYS)) == collect(rootcells(BASE))
    @test has_sorted_subtrees(SYS) == has_sorted_subtrees(BASE)
    @test has_congruent_refinement(SYS) == has_congruent_refinement(BASE)
    @test DGG.cap_inflation(SYS) == DGG.cap_inflation(BASE)
    @test maxneighbors(SYS) == maxneighbors(BASE)
    @test maxneighbors(SYS, Edge()) == maxneighbors(BASE, Edge())
    g = levelgrid(SYS, 3)
    for c in spread(g, 12)
        @test collect(children(SYS, c)) == collect(children(BASE, c))
        @test parent(SYS, c) == parent(BASE, c)
        @test ancestor(SYS, c, 1) == ancestor(BASE, c, 1)
        @test collect(descendants(SYS, c, 5)) == collect(descendants(BASE, c, 5))
        @test descendant_range(SYS, c, 5) == descendant_range(BASE, c, 5)
        @test collect(border(subtree(SYS, c, 5); cells = true)) ==
              collect(border(subtree(BASE, c, 5); cells = true))
        @test collect(neighbors(g, c)) == collect(neighbors(levelgrid(BASE, 3), c))
        @test collect(ring(g, c, 2)) == collect(ring(levelgrid(BASE, 3), c, 2))

        # The rotational order is measured about the wrapped grid's warped
        # centroids, so run the order collectors directly on this geometry.
        @test CT.winding_problems(g, c, collect(ring(g, c, 1)); label="warped ring 1") == String[]
        @test CT.neighbor_order_problems(g, c; k=2) == String[]
        @test CT.neighbor_problems(g, c) == String[]
    end
    @test_throws ArgumentError parent(SYS, first(rootcells(SYS)))
    @test_throws ArgumentError levelgrid(SYS, -1)
end

@testset "a second wrap is refused" begin
    g = levelgrid(SYS, 2)
    @test_throws ArgumentError AuthalicGrid(g)
    @test_throws ArgumentError AuthalicGrid(g, GO.Geodesic())
    @test_throws ArgumentError AuthalicSystem(SYS)
    # The message says how to get back rather than only what went wrong.
    @test occursin("parent(grid)", sprint(showerror, try
        AuthalicGrid(g)
    catch err
        err
    end))

    # A subset composes the other way round: the wrapper goes on the system.
    ids = [cellindex(levelgrid(BASE, 3), i) for i in 1:20]
    @test_throws ArgumentError AuthalicGrid(PartialGrid(BASE, 3, ids))
    @test occursin("PartialGrid(AuthalicSystem(sys), level, ids)",
        sprint(showerror, try
            AuthalicGrid(PartialGrid(BASE, 3, ids))
        catch err
            err
        end))

    # An unnamed frame is refused where it always is.
    @test_throws ArgumentError AuthalicGrid(levelgrid(BASE, 2), GO.Planar())
    @test_throws ArgumentError AuthalicSystem(BASE, GO.AutoManifold())
end

# ===========================================================================
# 3. Geometry
# ===========================================================================

@testset "boundary and centroid are warped, and only in latitude" begin
    for l in (0, 3, 6)
        g = levelgrid(SYS, l)
        b = levelgrid(BASE, l)
        for c in spread(g, 8)
            wb = cell_boundary(g, c)
            bb = cell_boundary(b, c)
            @test length(wb) == length(bb)
            for (p, q) in zip(bb, wb)
                plon, plat = FB.lonlat(p)
                qlon, qlat = FB.lonlat(q)
                abs(plat) == 90 || @test qlon ≈ plon atol = 1e-12
                @test qlat ≈ H.authalic_to_geodeticd(WGS84, plat) atol = 1e-12
                @test abs(sqrt(q[1]^2 + q[2]^2 + q[3]^2) - 1) < 1e-14
            end
            wc = cell_centroid(g, c)
            @test wc == FB.geodetic_point(WGS84, cell_centroid(b, c))
            # Repeated calls agree — the harness checks this too, and it is what
            # a lazily warped container would break.
            @test cell_boundary(g, c) == wb
        end
    end

    # A grid wrapped on a sphere is the grid, coordinate for coordinate.
    flat = AuthalicGrid(levelgrid(BASE, 3), GO.Spherical(; radius=1.0))
    for c in spread(flat, 8)
        @test cell_boundary(flat, c) == cell_boundary(levelgrid(BASE, 3), c)
        @test cell_centroid(flat, c) === cell_centroid(levelgrid(BASE, 3), c)
    end

    # The warp is a homeomorphism of the sphere and adjacent cells share their
    # corners, so the warped cells still tile the sphere exactly.
    g = levelgrid(SYS, 3)
    @test sum(cell_area(g, cellindex(g, i)) for i in 1:ncells(g)) ≈ 4pi rtol = 1e-12
    # ... but they are no longer equal-area: that is what the warp costs, and
    # what makes the base grid the one to read areas off.
    areas = [cell_area(g, cellindex(g, i)) for i in 1:ncells(g)]
    @test maximum(areas) / minimum(areas) > 1.005
end

@testset "cellat inverse-warps its input" begin
    for l in (2, 5)
        g = levelgrid(SYS, l)
        for c in spread(g, 24)
            centroid = cell_centroid(g, c)
            @test cellat(g, centroid) == c
            # Interior points that are not the centroid, so that a lookup that
            # merely matched centroids could not pass.
            for v in cell_boundary(g, c)
                @test cellat(g, US.slerp(centroid, v, 0.5)) == c
            end
            # The degree method takes geodetic longitude and latitude.
            lon, lat = FB.lonlat(centroid)
            @test cellat(g, lon, lat) == c
        end
    end

    # Feeding it the unwarped (authalic) latitude asks a different question, and
    # the answer differs once the cells are smaller than the shift — level 12
    # pixels are 0.014° across, so a 0.128° miss is a dozen cells wide. The
    # shift vanishes at the equator and at the poles, so the cells that sit
    # there are not evidence either way and are skipped rather than asserted on.
    fine = levelgrid(SYS, 12)
    unwarped = levelgrid(BASE, 12)
    displaced = 0
    for c in spread(fine, 16)
        geodetic = FB.lonlat(cell_centroid(fine, c))
        authalic = FB.lonlat(cell_centroid(unwarped, c))
        @test cellat(fine, geodetic...) == c
        abs(geodetic[2] - authalic[2]) < 0.05 && continue
        displaced += 1
        @test cellat(fine, authalic...) != c
    end
    @test displaced > 0
end

# ===========================================================================
# 4. The covering law
# ===========================================================================

@testset "node_extent inflates, and has to" begin
    stretch = FB.authalic_stretch(WGS84)
    for l in 0:3
        g = levelgrid(SYS, l)
        for c in spread(g, 6)
            basecap = node_extent(BASE, c)
            cap = node_extent(SYS, c)
            @test cap.radius ≈ stretch * basecap.radius rtol = 1e-12
            @test cap.point == FB.geodetic_point(WGS84, basecap.point)
            # Convex, so the harness's vertex sampling is a sound proxy.
            @test cap.radius <= pi / 2
        end
    end

    # NECESSITY. Against the tightest cap that is sound in the base frame — the
    # exact maximum distance from the cell's centroid to any descendant vertex —
    # the warp pushes descendants OUT for some cells: a wrapper that forwarded
    # `node_extent` unchanged would under-cover, which is the silent kind of
    # wrong. The same measurement shows the bound is not merely necessary but
    # sufficient: `stretch` covers every one of them.
    pushed_out = 0
    checked = 0
    worst = 0.0
    for l in (2, 3)
        g = levelgrid(BASE, l)
        wg = levelgrid(SYS, l)
        for c in spread(g, 8)
            centre = cell_centroid(g, c)
            wcentre = cell_centroid(wg, c)
            deep = levelgrid(BASE, l + 2)
            wdeep = levelgrid(SYS, l + 2)
            tightest = 0.0
            warped = 0.0
            for d in descendants(BASE, c, l + 2)
                for p in cell_boundary(deep, d)
                    tightest = max(tightest, US.spherical_distance(centre, p))
                end
                for p in cell_boundary(wdeep, d)
                    warped = max(warped, US.spherical_distance(wcentre, p))
                end
            end
            checked += 1
            warped > tightest && (pushed_out += 1)
            worst = max(worst, warped / tightest)
            @test warped <= stretch * tightest
        end
    end
    @test checked > 0
    @test pushed_out > 0
    @test worst > 1

    # The convexity threshold, from both sides. A base extent wide enough that
    # `stretch` would carry it past 90° loses the argument that containing the
    # warped vertices contains the arcs between them, and the only sound answer
    # left is the whole sphere. No system in this package produces one — the
    # widest is 0.907 rad — so the branch is pinned here on a stub rather than
    # left to a system that might one day grow into it.
    threshold = (pi / 2) / stretch
    c0 = DGG.LevelIndex(0, 0)
    narrow = node_extent(AuthalicSystem(WideExtentStub(0.999 * threshold)), c0)
    @test narrow.radius <= pi / 2
    @test narrow.radius < FB.full_sphere_cap().radius
    # ... and everything above it, including a base extent that is already the
    # full sphere. Note that `1.5` would NOT be here: 1.5 · L is 1.507, still
    # inside the quadrant, which is the arithmetic this threshold is about.
    for r in (1.001 * threshold, 1.57, pi / 2, FB.full_sphere_cap().radius)
        wide = node_extent(AuthalicSystem(WideExtentStub(r)), c0)
        @test wide.radius == FB.full_sphere_cap().radius
        @test wide.point == FB.full_sphere_cap().point
    end

    # SOUNDNESS, directly: every descendant vertex, several levels down, inside
    # the extent of every ancestor on the path. (The conformance suite below
    # does this too, on its own sampled cells — this one is deterministic.)
    for l in (0, 2)
        g = levelgrid(SYS, l)
        for c in spread(g, 4)
            chain = [(c, node_extent(SYS, c))]
            cur = c
            for _ in 1:3
                cur = children(SYS, cur)[end]
                push!(chain, (cur, node_extent(SYS, cur)))
                deep = levelgrid(SYS, level(cur))
                for p in cell_boundary(deep, cur), (_, cap) in chain
                    @test US.spherical_distance(cap.point, p) <= cap.radius + 1e-12
                end
            end
        end
    end
end

# ===========================================================================
# 5. Conformance
# ===========================================================================

@testset "conformance suites" begin
    for l in (0, 3)
        test_grid_interface(levelgrid(SYS, l); label="AuthalicGrid(HEALPix $l)")
    end
    test_grid_interface(AuthalicGrid(levelgrid(IGeo7System(), 2));
        label="AuthalicGrid(IGeo7 2)")
    test_hierarchical_system(SYS; label="AuthalicSystem(HEALPix)")
    test_hierarchical_system(AuthalicSystem(IGeo7System()); label="AuthalicSystem(IGeo7)")
end

# ===========================================================================
# 6. Trees, subsets and regridding
# ===========================================================================

@testset "PartialGrid over an AuthalicSystem" begin
    complete = levelgrid(SYS, 4)
    ids = [cellindex(complete, i) for i in 1:7:ncells(complete)]
    pg = PartialGrid(SYS, 4, ids)
    @test system(pg) === SYS
    @test ncells(pg) == length(ids)
    @test cell_boundary(pg, ids[3]) == cell_boundary(complete, ids[3])

    # The tree descends wrapped node extents in the subset's local-index space.
    # Both the tree query and fallback `cellat` must retain each stored cell's
    # centroid through every pruning step.
    tree = treeify(pg)
    for c in ids[1:13:end]
        hits = STI.query(tree, cap -> FB.cap_contains(cap, cell_centroid(pg, c)))
        @test localindex(pg, c) in hits
        @test cellat(pg, cell_centroid(pg, c)) == c
    end

    # The subtree form, which is where `descendant_range` has to survive the wrap.
    root = cellindex(levelgrid(SYS, 1), 5)
    sub = subtree(SYS, root, 4)
    @test ncells(sub) == length(descendants(SYS, root, 4))
    @test system(sub) === SYS
    subtree_tree = treeify(sub)
    for i in (1, 7, ncells(sub))
        hits = STI.query(subtree_tree,
            cap -> FB.cap_contains(cap, cell_centroid(sub, cellindex(sub, i))))
        @test i in hits
    end
end

@testset "regridding against a geodetic lon/lat grid" begin
    # A 5-degree lon/lat destination, as unit-sphere corner points: the frame
    # every terrestrial dataset is in, and the reason this wrapper exists.
    to_sphere = US.UnitSphereFromGeographic()
    nlon, nlat = 72, 36
    lons = collect(range(0, 360; length=nlon + 1))
    lats = collect(range(-90, 90; length=nlat + 1))
    dst = [to_sphere((x, y)) for x in lons, y in lats]
    manifold = GO.Spherical(; radius=1.0)

    src = levelgrid(SYS, 6)
    base_src = levelgrid(BASE, 6)
    @test GOCore.best_manifold(src) === manifold      # NOT the authalic sphere
    regridder = CR.Regridder(manifold, dst, src)
    base_regridder = CR.Regridder(manifold, dst, base_src)

    # Conservation: the destination tiles the sphere, so every source cell's
    # column of intersection areas sums to its own area, and both budgets are 4pi.
    A = regridder.intersections
    @test size(A) == (nlon * nlat, ncells(src))
    @test maximum(abs.(vec(sum(A; dims=1)) .- regridder.src_areas) ./ regridder.src_areas) < 1e-10
    @test maximum(abs.(vec(sum(A; dims=2)) .- regridder.dst_areas)) /
          maximum(regridder.dst_areas) < 1e-10
    @test sum(regridder.src_areas) ≈ 4pi rtol = 1e-12

    # REGISTRATION. `z` is the sine of GEODETIC latitude on the wrapped grid, so
    # sampling it at wrapped cell centroids is a geodetically registered field,
    # and a lon/lat cell's exact average of it is (sin(lat0) + sin(lat1))/2 — an
    # analytic destination truth that owes nothing to either regridder.
    analytic = vec([(sind(lats[j]) + sind(lats[j+1])) / 2 for i in 1:nlon, j in 1:nlat])
    values = [cell_centroid(src, cellindex(src, i))[3] for i in 1:ncells(src)]

    wrapped = zeros(length(regridder.dst_areas))
    CR.regrid!(wrapped, regridder, values)
    unwrapped = zeros(length(base_regridder.dst_areas))
    CR.regrid!(unwrapped, base_regridder, values)

    rms(x) = sqrt(sum(abs2, x) / length(x))
    err_wrapped = rms(wrapped .- analytic)
    err_unwrapped = rms(unwrapped .- analytic)
    # The wrapped source lands where the data says it does, to the source grid's
    # own discretisation; the unwrapped one is off by the authalic shift, whose
    # scale is `cos(lat) * 2.24e-3` in `z` and shows up as a rms five times
    # larger. That gap IS the 14 km misregistration the wrapper removes.
    @test err_wrapped < err_unwrapped / 4
    @test err_wrapped < 5e-4
    @test 5e-4 < err_unwrapped < 5e-3

    # Conservative means conservative: the area-weighted mean survives exactly.
    @test sum(values .* regridder.src_areas) / sum(regridder.src_areas) ≈
          sum(wrapped .* regridder.dst_areas) / sum(regridder.dst_areas) atol = 1e-14
end

end # @testset "authalic wrapper"

end # module AuthalicWrapperTests

# ---------------------------------------------------------------------------
# The authalic latitude math itself (`src/Helpers/authalic.jl`).
#
# These tests cover the authalic-latitude series used by `AuthalicSystem` and
# compare it with two independent reference sets.
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
#
# Neither set is regenerable from this repo. Do not "fix" a failure by editing
# a tuple.
# ---------------------------------------------------------------------------

module HelpersAuthalicTests

using Test
using DiscreteGlobalGrids
using DiscreteGlobalGrids.Helpers
import DiscreteGlobalGrids as DGG

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

        Shape = DGG.Helpers.EllipsoidShapeError
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
        repo = DGG.ISEA.R_AUTHALIC
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

end # module HelpersAuthalicTests
