# ---------------------------------------------------------------------------
# T9 — the ellipsoid wrapper (`src/fallbacks/authalic_grid.jl`).
#
# Four kinds of test here, and the distinction matters when one fails:
#
#   1. THE BOUNDS. `authalic_shift` and `authalic_stretch` are derived
#      analytically from the series coefficients, so they are checked against a
#      dense sweep of the transform they claim to bound — and against the
#      Lipschitz property itself, measured on point pairs rather than restated
#      from the derivation.
#   2. FORWARDING. Everything that is not geometry must come through the wrapper
#      unchanged, ids and positions above all.
#   3. THE COVERING LAW. The conformance harness is the oracle, but the
#      necessity of the inflation is shown separately: against the TIGHTEST
#      sound base cap the warp really does push descendants out, so a wrapper
#      that forwarded `node_extent` would under-cover.
#   4. REGISTRATION. The point of the whole type: a wrapped grid regridded
#      against a geodetic lon/lat grid lands where the data says it does, and an
#      unwrapped one misses by the authalic shift.
#
# The conformance harness skips `cellat`, `neighbors` and `ring` on this grid —
# they dispatch into `DGG.Fallbacks`, which is its "the system did not implement
# this" sentinel (the same blind spot it documents for `PartialGrid`) — so those
# three are tested directly below.
# ---------------------------------------------------------------------------

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
const WGS84 = H.WGS84_AUTHALIC
const BASE = HEALPixSystem()
const SYS = AuthalicSystem(BASE)

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

@testset "ids, positions and hierarchy forward unchanged" begin
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
            @test cellposition(g, c) == i
        end
        # A cell that is not in the grid is a `nothing`, never an error.
        l < max_level(SYS) && @test cellposition(g, first(children(SYS, cellindex(g, 1)))) === nothing
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
    @test max_level(SYS) == max_level(BASE)
    @test collect(rootcells(SYS)) == collect(rootcells(BASE))
    @test has_sorted_subtrees(SYS) == has_sorted_subtrees(BASE)
    @test cap_inflation(SYS) == cap_inflation(BASE)
    @test max_neighbors(SYS) == max_neighbors(BASE)
    @test max_neighbors(SYS, Edge()) == max_neighbors(BASE, Edge())
    g = levelgrid(SYS, 3)
    for c in spread(g, 12)
        @test collect(children(SYS, c)) == collect(children(BASE, c))
        @test parent(SYS, c) == parent(BASE, c)
        @test ancestor(SYS, c, 1) == ancestor(BASE, c, 1)
        @test collect(descendants(SYS, c, 5)) == collect(descendants(BASE, c, 5))
        @test descendant_range(SYS, c, 5) == descendant_range(BASE, c, 5)
        @test collect(subtree_border(SYS, c, 5)) == collect(subtree_border(BASE, c, 5))
        @test collect(neighbors(g, c)) == collect(neighbors(levelgrid(BASE, 3), c))
        @test collect(ring(g, c, 2)) == collect(ring(levelgrid(BASE, 3), c, 2))
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
            # The degree method takes GEODETIC lon/lat, which is the whole point.
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

    # The tree descends the WRAPPED node extents over the subset's own position
    # space, so a stored cell's own centroid must survive every prune on the way
    # down to it. Asserted on the tree rather than through `cellat`, because the
    # generic point-in-cell test that `cellat` finishes with misjudges the
    # centroid of ~9% of HEALPix's densified rings — a pre-existing weakness of
    # the fallback locate path (HEALPix's own `cellat` never reaches it), which
    # hits the base grid and the wrapped grid identically and has nothing to say
    # about the warp.
    tree = treeify(pg)
    for c in ids[1:13:end]
        hits = STI.query(tree, cap -> FB.cap_contains(cap, cell_centroid(pg, c)))
        @test cellposition(pg, c) in hits
    end

    # The subtree form, which is where `descendant_range` has to survive the wrap.
    root = cellindex(levelgrid(SYS, 1), 5)
    sub = PartialGrid(SYS, root, 4)
    @test ncells(sub) == length(descendants(SYS, root, 4))
    @test system(sub) === SYS
    subtree = treeify(sub)
    for i in (1, 7, ncells(sub))
        hits = STI.query(subtree, cap -> FB.cap_contains(cap, cell_centroid(sub, cellindex(sub, i))))
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
