# ---------------------------------------------------------------------------
# Copernicus DEM system tests.
#
# Three kinds of test, and the distinction matters when one fails:
#
#   1. ORACLE. The lattice itself — which tiles are how many columns wide, and
#      where their pixel centres sit — is checked against 79 real AWS Open Data
#      COGs, measured once with `ArchGDAL` and committed to `fixtures.jl`. No
#      network is touched from `test/`; the fixtures ARE the oracle. They are the
#      only external evidence for the one rule neither primary source states: a
#      tile's band is chosen by its EQUATOR-WARD edge, so `S50` is 1x and full
#      width while `N50` is 1.5x.
#
#      HOW MUCH OF EACH SWEEP. The tile lattice is 64 800 cells, which is small
#      enough to walk exhaustively in seconds even at GLO-30, so the prefix-sum,
#      pole-row and 4-pi testsets below iterate all of it rather than sampling —
#      they say so where they do. Level 1 is 620 524 800 000 cells at GLO-30 and
#      is always sampled, from a seeded `MersenneTwister`, or swept along one
#      structural line (a raster row, a tile's pixels, the 180 tile rows).
#
#   2. STRUCTURAL. Ids, prefix sums, raster order, and the pole degeneracies:
#      this package's own design, which no oracle has an opinion about.
#
#   3. GEOMETRY. That the boxes partition the sphere, that the published rings
#      are convex and how far they sit from their boxes, that `node_extent`
#      covers the subtree, and that `cellat` inverts `cell_box` on the edges
#      where floating point makes the two disagree.
#
#   4. CONTRACT. The two conformance suites, in section (k): in full on a scaled
#      twin of the shipped lattices, and — because the harness materialises a
#      cell's whole child list, and reduces one range per child quadratically —
#      at the interface level on GLO-30 and GLO-90 themselves, with their
#      hierarchical runs behind `DGG_COPDEM_FULL=1`. That section carries the
#      seed the harness's absolute `area_atol` forces on this system's pole rows,
#      and the assertion that the seed still does what it claims.
#
# Every testset names the mutant it kills. A test that kills no mutant no other
# test already kills does not belong here.
#
# `test/runtests.jl` includes this file; it also runs standalone:
#     julia --project=test --startup-file=no test/systems/CopernicusDEM/runtests.jl
# ---------------------------------------------------------------------------

module CopernicusDEMSystemTests

using Test
using Random

using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
const CD = DiscreteGlobalGrids.CopernicusDEM

using DiscreteGlobalGridsConformanceTesting
# ... and by name as well, for `sample_positions` and `boundary_problems`, which
# section (k) calls directly: the seed-away assertion there has to reproduce the
# harness's own draw and the harness's own verdict, not an imitation of either.
import DiscreteGlobalGridsConformanceTesting as CT

import GeometryOps as GO
import GeoInterface as GI
const US = GO.UnitSpherical
using GeometryOps.UnitSpherical: spherical_orient

# The manifold every grid in this package computes on. Named once: a bare vector
# of `UnitSphericalPoint`s carries no manifold, and `best_manifold` would guess
# the WGS84 sphere for it, which is a factor of R^2 in every area.
const MANIFOLD = GO.Spherical(; radius = 1.0)

const GLO30 = DGG.CopernicusDEMSystem(30)
const GLO90 = DGG.CopernicusDEMSystem(90)
# The scaled twin: the same code, the same band table, the same pole clamps, at
# 1/120 the pixel count. Task 5 runs the conformance harness on it; here it is a
# third independent `N` for every law that is stated in `N`.
const TWIN = CD.CopernicusDEMSystem{30}()
const ALL_SYSTEMS = (GLO30, GLO90, TWIN)

include("fixtures.jl")

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

# Tile rows that straddle every band edge in both hemispheres, plus the equator
# and both pole rows. This is the latitude axis of every geometry sweep below:
# the band edges are where the column count changes, and the pole rows are where
# the cells degenerate.
const PROBE_LATS = (89, 88, 85, 84, 80, 79, 70, 69, 60, 59, 50, 49, 0, -1,
                    -50, -51, -60, -61, -70, -71, -80, -81, -85, -86, -89, -90)
# The antimeridian tile from both sides, and the prime meridian.
const PROBE_LONS = (-180, 0, 179)

"""
A sweep's failures, capped. These loops run to 64 800 iterations and the
assertion is "none of them failed", not "each of them passed" — one `@test` per
cell would be 200 000 test results for no extra information. The cap keeps a
failure message readable while still naming the first few offenders.
"""
note!(bad::Vector{String}, msg) = (length(bad) < 5 && push!(bad, string(msg)); bad)

# Does this ring turn right anywhere, i.e. is it non-convex? Copied VERBATIM
# from `test/systems/crosssystem/regridding_conservation.jl:120-132` — copied
# rather than `include`d, because the two suites are separate modules and neither
# should be able to break the other by editing a shared helper. Consecutive
# repeated vertices are skipped: a zero-length edge has no turn to measure, and
# `spherical_orient` goes through `robust_cross_product`, which returns an
# arbitrary perpendicular for two identical points, so a duplicated vertex would
# otherwise read as a random reflex turn.
function has_reflex_vertex(poly)
    pts = collect(GI.getpoint(GI.getexterior(poly)))
    while length(pts) > 1 && pts[end] == pts[1]
        pop!(pts)
    end
    keep = [i for i in eachindex(pts) if i == 1 || pts[i] != pts[i - 1]]
    pts = pts[keep]
    n = length(pts)
    n < 4 && return false
    return any(1:n) do i
        spherical_orient(pts[i], pts[mod1(i + 1, n)], pts[mod1(i + 2, n)]) < 0
    end
end

_walk(z, k) = (for _ in 1:abs(k); z = k > 0 ? nextfloat(z) : prevfloat(z); end; z)

"""
The `UnitSphericalPoint` on the prime meridian whose decoded latitude is EXACTLY
`lat`, or `nothing` when no such point exists.

`GeographicFromUnitSphere ∘ UnitSphereFromGeographic` is not the identity on
latitude — it is `asind ∘ cosd(90 - ·)` in the third coordinate, which moves a
value by up to 10 ulps and, being expansive near the equator, SKIPS values
entirely. So the naive `TO_SPHERE((0, lat))` usually decodes to a neighbouring
`Float64`, and testing `cellat` on an exact box edge means searching a few ulps
of `z` around the naive image for one that decodes back exactly. 40 steps each
way is ample; a third of all tile-row south edges have no such point at all.
"""
function exact_latitude_point(lon, lat)
    p0 = CD.TO_SPHERE((lon, lat))
    z0 = p0[3]
    for k in 0:40
        for z in (k == 0 ? (z0,) : (_walk(z0, k), _walk(z0, -k)))
            abs(z) > 1 && continue
            q = GO.UnitSphericalPoint(p0[1], p0[2], z)
            CD.FROM_SPHERE(q)[2] === lat && return q
        end
    end
    return nothing
end

# The seed section (k)'s conformance calls run on, and the sample size they run
# it at. Both are named here rather than left to the harness's defaults because
# `sampled_tile_lats` below has to reproduce the harness's draw EXACTLY — see
# section (k) for what the reproduction is for.
const CONFORMANCE_SEED = 20260815
const GI_SAMPLES = 32

"""
The tile-row latitudes `test_grid_interface(levelgrid(sys, l); n_samples, rng =
MersenneTwister(seed))` will sample: the harness's own draw, reproduced.

`sample_positions` is the whole of the harness's sampling — `sample_cells` is it
composed with `cellindex` — and `test_grid_interface` calls it once, before any
law runs, so a caller can predict the cells from the seed alone. Calling the
harness's own function rather than reimplementing `rand` is the point: an
upstream change to how positions are drawn moves this too, and section (k)'s
assertion goes red instead of silently testing different cells.
"""
function sampled_tile_lats(sys, l, seed, n_samples)
    g = levelgrid(sys, l)
    positions = CT.sample_positions(MersenneTwister(seed), ncells(g), n_samples)
    return [CD.tilecorner(sys, cellindex(g, p))[1] for p in positions]
end

@testset "CopernicusDEM system" begin

# =========================================================================
# (a) The band table, against the measured tiles
# =========================================================================

# KILLS: choosing the band from the tile LABEL rather than from
# `min(|lat_s|, |lat_s + 1|)` — the `S50`/`S60`/`S70`/`S80`/`S85` rows are the
# only ones that separate the two — and closing the band intervals on the wrong
# side.
@testset "the band is the equator-ward edge" begin
    for (sys, N, fixtures) in ((GLO30, 3600, GLO30_TILES), (GLO90, 1200, GLO90_TILES))
        for f in fixtures
            @test CD.ncols_at(sys, f.lat_s) == f.ncols
            t = CD.tilecell(sys, f.lat_s, f.lon_w)
            @test length(children(sys, t)) == f.ncols * N
        end
    end

    # The six half-open band intervals, read off at the pair of rows that
    # straddles each edge. The southern member of each pair is the one the tile
    # LABEL gets wrong: `S50` spans -50 to -49, so its equator-ward edge is 49
    # and it is a 1x tile, while `N50` spans 50 to 51 and is 1.5x.
    for (lat_s, cols30, cols90) in (( 49, 3600, 1200), ( 50, 2400,  800),
                                    ( 59, 2400,  800), ( 60, 1800,  600),
                                    ( 69, 1800,  600), ( 70, 1200,  400),
                                    ( 79, 1200,  400), ( 80,  720,  240),
                                    ( 84,  720,  240), ( 85,  360,  120),
                                    (-50, 3600, 1200), (-51, 2400,  800),
                                    (-60, 2400,  800), (-61, 1800,  600),
                                    (-70, 1800,  600), (-71, 1200,  400),
                                    (-80, 1200,  400), (-81,  720,  240),
                                    (-85,  720,  240), (-86,  360,  120))
        @test CD.ncols_at(GLO30, lat_s) == cols30
        @test CD.ncols_at(GLO90, lat_s) == cols90
        # and the twin, which is the same table at 1/120 the width
        @test CD.ncols_at(TWIN, lat_s) == cols30 ÷ 120
    end
    # The two rows that touch a pole are both in the 10x band.
    @test CD.ncols_at(GLO30, 89) == 360
    @test CD.ncols_at(GLO30, -90) == 360
end

# =========================================================================
# (b) Registration: the measured geotransforms
# =========================================================================

# KILLS: dropping the half-pixel outset, or applying it with the wrong sign in y
# — the origin is half a pixel WEST and half a pixel NORTH of the first pixel
# centre, and a sign flip in y lands half a pixel south, which is 1e-4 degrees at
# GLO-30 and would sail past any `rtol`.
@testset "registration reproduces the measured geotransforms" begin
    for f in GEOTRANSFORMS
        sys = CD.CopernicusDEMSystem{f.N}()
        t = CD.tilecell(sys, f.lat_s, f.lon_w)
        west, _, _, north = CD.cell_box(sys, t)
        # The raster SPACINGS are read off an interior pixel row, not off the
        # tile box's own height. `cell_box` extends the `lat_s = -90` bottom row
        # down to -90 to close the half-pixel gap ring at the pole (and clamps
        # the `lat_s = 89` top row up to +90), so those two tile rows are N + 1/2
        # pixels tall BY THIS SYSTEM'S OWN CONVENTION and `(north - south) / N`
        # is half a pixel out: measured 3.86e-8 deg at GLO-30 and 3.47e-7 at
        # GLO-90. The COG has no such row — its south edge is -89.99986111111112
        # — so the fixture's Δlat is the ordinary row height, and row 1 is
        # unaffected by either correction in both rows.
        pw, pe, ps, pn = CD.cell_box(sys, CD.pixelcell(sys, t, 1, 0))
        gt = (west, pe - pw, 0.0, north, 0.0, -(pn - ps))
        expected = (f.origin_x, f.dlon_arcsec / 3600, 0.0, f.origin_y, 0.0, -1 / f.N)
        @test all(abs.(gt .- expected) .<= 1e-12)

        # The first pixel CENTRE is the integer degree pair the file name carries,
        # exactly — this is what "pixel-is-point" means and it is not an
        # approximation, so `==` and not `isapprox`.
        cw, ce, cs, cn = CD.cell_box(sys, CD.pixelcell(sys, t, 0, 0))
        @test (cw + ce) / 2 == Float64(f.lon_w)
        @test (cs + cn) / 2 == Float64(f.lat_s + 1)
    end
end

# =========================================================================
# (c) Ids, prefix sums, and the descendant windows
# =========================================================================

# KILLS: an off-by-one in the prefix sums — which corrupts every tile after the
# first, and which no single-tile test can see — and a row/column transposition
# in `tilebase`, which survives any test that only looks at square rows.
@testset "ids, prefix sums and descendant windows" begin
    for sys in ALL_SYSTEMS
        N = CD.lat_intervals(sys)
        g0 = levelgrid(sys, 0)
        g1 = levelgrid(sys, 1)

        # EXHAUSTIVE over all 64 800 tiles, in level-0 order: their level-1
        # windows must abut with no gap and no overlap, start at 1, and end at
        # `ncells(sys, 1)`. That is the whole partition, stated once.
        bad = String[]
        prev_stop = 0
        for r in 0:179, q in 0:359
            t = LevelIndex(0, CD.tileordinal(r, q))
            nc = CD.ncols(sys, r)
            base = CD.tilebase(sys, r, q)
            w = descendant_range(sys, t, 1)
            w == Int(base + 1):Int(base + nc * Int64(N)) ||
                note!(bad, "window $(w) != tilebase form at (r=$r, q=$q)")
            first(w) == prev_stop + 1 ||
                note!(bad, "gap or overlap before (r=$r, q=$q): $(first(w)) after $prev_stop")
            prev_stop = last(w)
        end
        @test bad == String[]
        @test prev_stop == ncells(sys, 1)
        # `prev_stop` IS `sum(ncols(sys, r) * N * 360 for r in 0:179)` — the loop
        # above walked that sum one window at a time — so asserting the closed
        # form separately restates the loop and is not a test.

        # positions -> ids -> positions, at both levels. Only this direction: the
        # other one, `cellindex(g, cellposition(g, c)) == c`, is the line above
        # composed with itself on a `c` this loop built with `cellindex`, and
        # cannot fail unless the line above already has.
        rng = MersenneTwister(20260815)
        for g in (g0, g1)
            n = ncells(g)
            for i in unique([1, 2, n - 1, n, rand(rng, 1:n, 32)...])
                c = cellindex(g, i)
                @test cellposition(g, c) == i
            end
        end
    end
end

# =========================================================================
# (d) Within-tile order
# =========================================================================

# KILLS: column-major ordering, which would silently transpose every DEM raster
# read through this system, and a registration drift that breaks the
# bit-identical shared edges the quad tessellation depends on.
@testset "within-tile order is AWS raster order" begin
    for sys in ALL_SYSTEMS
        N = CD.lat_intervals(sys)
        for lat_s in (0, 50, 60, 70, 80, 85, -1, -51, -86, 89, -90),
            lon_w in (-180, 6, 179)

            t = CD.tilecell(sys, lat_s, lon_w)
            r, q, _, _ = CD.decode(sys, t)
            nc = Int(CD.ncols(sys, r))
            pg = PartialGrid(sys, t, 1)
            @test ncells(pg) == nc * N
            _, _, tile_s, tile_n = CD.cell_box(sys, t)

            for j in (0, 1, N ÷ 2, N - 2, N - 1), i in (0, 1, nc ÷ 2, nc - 2, nc - 1)
                c = CD.pixelcell(sys, t, j, i)
                @test CD.decode(sys, c) == (r, q, j, i)
                # north row first, west to east: position `j * ncols + i + 1`.
                k = j * nc + i + 1
                @test cellindex(pg, k) == c
                @test cellposition(pg, c) == k

                _, east, south, north = CD.cell_box(sys, c)
                if i < nc - 1
                    # `===`, not `==`: the east edge of column `i` and the west
                    # edge of column `i+1` are the SAME Float64, which is what
                    # makes the two quads share a geodesic rather than nearly
                    # share one.
                    @test east === CD.cell_box(sys, CD.pixelcell(sys, t, j, i + 1))[1]
                end
                if j < N - 1
                    @test south === CD.cell_box(sys, CD.pixelcell(sys, t, j + 1, i))[4]
                end
                j == 0 && @test north === tile_n
                j == N - 1 && @test south === tile_s
            end
        end
    end
end

# =========================================================================
# (e) The pole rows
# =========================================================================

# This testset carries the pole rows ON ITS OWN. The conformance harness cannot:
# `boundary_problems` judges a ring degenerate against an ABSOLUTE `1e-12`
# steradian floor that no caller can reach, and a GLO-30 `lat_s = 89` top-row
# pixel is a legitimate 1.4e-16 sr. Task 5 samples those rows away; everything
# anyone would want to assert about them lives here instead, which is why the
# level-0 sweep below is exhaustive rather than sampled.
#
# KILLS: emitting a duplicated pole vertex (a degenerate quad that
# `has_reflex_vertex` reads as a random reflex turn), building the pole from
# `UnitSphereFromGeographic((lon, ±90))` rather than from the `NORTH_POLE` /
# `SOUTH_POLE` literals, and forgetting the clamp/extend, so the top row runs
# past latitude 90.
#
# That middle mutant is killed by ONE character — the `===` in the pole-vertex
# count below — and only because the level-0 sweep is exhaustive. The mutant's
# point is NOT ~6e-17 off the pole: GeometryOps goes through `sincosd`,
# `cosd(±90)` is exactly `0.0`, and `TO_SPHERE((lon, ±90))` is exactly
# `(±0.0, ±0.0, ±1.0)` — four bit patterns whose only difference from the literal
# is the SIGN BIT, set in x wherever `cosd(lon) < 0` and in y wherever
# `sind(lon) < 0`. Numerically the two points never differ, so `==` cannot see it
# (`-0.0 == 0.0`) and `===` can. Measured: 91 of the 360 integer longitudes —
# exactly `0:90`, where both sign bits are clear — reproduce the literal
# bit-for-bit, so a probe restricted to that quadrant would let the mutant live.
# Sweeping all of `lon_w in -180:179` is what makes the kill certain.
@testset "pole cells are triangles" begin
    for sys in ALL_SYSTEMS
        N = CD.lat_intervals(sys)
        g0 = levelgrid(sys, 0)
        g1 = levelgrid(sys, 1)
        bad = String[]

        function check_triangle!(bad, g, c, pole, label)
            ring = cell_boundary(sys, c)
            length(ring) == 3 ||
                return note!(bad, "$label: $(length(ring)) vertices, not 3")
            # In a 3-ring every pair of vertices is consecutive, so `allunique`
            # already is the consecutive-duplicate check.
            allunique(ring) || note!(bad, "$label: repeated vertex")
            all(p -> abs(hypot(p[1], p[2], p[3]) - 1) < 1e-15, ring) ||
                note!(bad, "$label: not unit norm")
            count(p -> p === pole, ring) == 1 ||
                note!(bad, "$label: pole vertex is not the exact literal")
            GO.area(MANIFOLD, cell_polygon(g, c)) > 0 ||
                note!(bad, "$label: signed area is not positive")
            return bad
        end

        # EXHAUSTIVE: all 360 `lat_s = 89` tiles and all 360 `lat_s = -90` tiles.
        for (lat_s, pole) in ((89, CD.NORTH_POLE), (-90, CD.SOUTH_POLE)),
            lon_w in -180:179

            check_triangle!(bad, g0, CD.tilecell(sys, lat_s, lon_w), pole,
                            "tile ($lat_s, $lon_w)")
        end

        # Level 1: the pole-most raster row of a seeded sample of pole tiles.
        rng = MersenneTwister(8990)
        for (lat_s, pole, j) in ((89, CD.NORTH_POLE, 0), (-90, CD.SOUTH_POLE, N - 1))
            nc = Int(CD.ncols_at(sys, lat_s))
            for lon_w in unique(rand(rng, -180:179, 16)), i in (0, nc ÷ 2, nc - 1)
                c = CD.pixelcell(sys, CD.tilecell(sys, lat_s, lon_w), j, i)
                check_triangle!(bad, g1, c, pole, "pixel ($lat_s, $lon_w, $j, $i)")
            end
        end

        # ... and a dense sweep down one tile of each pole row: EXACTLY one
        # raster row of each degenerates, the rest are honest quads, and no row
        # escapes [-90, 90].
        for (lat_s, jpole) in ((89, 0), (-90, N - 1))
            t = CD.tilecell(sys, lat_s, 0)
            for j in 0:(N - 1)
                c = CD.pixelcell(sys, t, j, 0)
                ring = cell_boundary(sys, c)
                _, _, south, north = CD.cell_box(sys, c)
                length(ring) == (j == jpole ? 3 : 4) ||
                    note!(bad, "row $j of tile $lat_s has $(length(ring)) vertices")
                (-90.0 <= south && north <= 90.0) ||
                    note!(bad, "row $j of tile $lat_s spans ($south, $north)")
                GO.area(MANIFOLD, cell_polygon(g1, c)) > 0 ||
                    note!(bad, "row $j of tile $lat_s is not counter-clockwise")
            end
        end
        @test bad == String[]

        # The clamp and the extension, stated exactly.
        @test CD.cell_box(sys, CD.pixelcell(sys, CD.tilecell(sys, 89, 0), 0, 0))[4] === 90.0
        @test CD.cell_box(sys, CD.pixelcell(sys, CD.tilecell(sys, -90, 0), N - 1, 0))[3] === -90.0

        # And the poles themselves locate, at both levels, into a cell whose cap
        # contains them.
        for pole in (CD.NORTH_POLE, CD.SOUTH_POLE), g in (g0, g1)
            c = cellat(g, pole)
            @test DGG.level(c) == g.level
            cap = node_extent(sys, c)
            @test US.spherical_distance(cap.point, pole) <= cap.radius
        end
    end
end

# =========================================================================
# (f) The boxes partition the sphere
# =========================================================================

# KILLS: using the pixel Δlon where the tile Δlon belongs (a factor of `ncols`),
# swapping `sin φ_N` and `sin φ_S`, and a missing pole clamp — all three of which
# show up as a total that is not 4π.
@testset "the boxes partition the sphere" begin
    for sys in ALL_SYSTEMS
        g0 = levelgrid(sys, 0)
        # MATERIALISE, then let Julia's pairwise `sum` reduce it. This is not
        # decoration: the closed form does not telescope (890 of the 64 800
        # GLO-90 tiles have `east - west !== 1.0`, by up to 1.4e-14), so a
        # generator argument — summed strictly left to right — lands 2.9e-12
        # (GLO-30) / 1.4e-12 (GLO-90) from 4π and fails this `rtol`, while the
        # pairwise reduction of the same numbers lands 3.6e-15 out. Do not
        # loosen the tolerance to accommodate a generator; the areas are fine.
        areas = [cell_area(g0, cellindex(g0, i)) for i in 1:ncells(g0)]
        @test sum(areas) ≈ 4π rtol = 1e-14

        # Latitude coverage is exactly [-90, 90], by `==`: tile row 0's north
        # edge is the clamp and row 179's south edge is the extension.
        @test CD.cell_box(sys, CD.tilecell(sys, 89, 0))[4] === 90.0
        @test CD.cell_box(sys, CD.tilecell(sys, -90, 0))[3] === -90.0
    end

    # One tile per band, its pixels summing to the tile. GLO-90 and the twin
    # only, and the real argument is the ARITHMETIC, not the clock: `cell_area`
    # is one closed form, generic in `N` and in the band's column count, so a
    # third lattice re-runs the same expression on bigger numbers and tests no
    # line the other two leave untested. The cost is only what settles it —
    # GLO-30's six tiles are 36.3 million pixels, measured in a fresh process
    # (GLO-30 first, so it pays every bit of compilation) at 0.82 s and 0.27 GiB
    # against 0.11 s and 0.033 GiB for GLO-90.
    for sys in (GLO90, TWIN)
        g0 = levelgrid(sys, 0)
        g1 = levelgrid(sys, 1)
        for lat_s in (0, 50, 60, 70, 80, 85)
            t = CD.tilecell(sys, lat_s, 6)
            pixels = [cell_area(g1, cellindex(g1, k))
                      for k in descendant_range(sys, t, 1)]
            @test sum(pixels) ≈ cell_area(g0, t) rtol = 1e-12
        end
    end
end

# =========================================================================
# (g) The published rings: convex, and how far from the box
# =========================================================================

# KILLS: any future densification of the edges — the change that would silently
# put this system in the same broken-as-a-regridding-DESTINATION class as
# HEALPix, ISEA4R and A5 — and a corner order that is clockwise or crossed.
#
# The thresholds are MEASURED worst cases, not round numbers, and they are
# two-scoped on purpose. The gap between the box (what `cell_area` returns) and
# the published geodesic quad (what `GO.area` measures) is strongly
# latitude-dependent, and the `lat_s = ±90` rows are half-pixel slivers and pole
# triangles whose gap is four orders above every other row: one flat threshold
# over a sample containing them would fail on geometry that is correct.
#
# Run on the two shipped lattices, which is what the thresholds were measured on.
# The twin is deliberately absent: its pixels are 1/30 of a degree tall and from
# 1/30 of a degree wide in the 1x band to 1/3 in the 10x band (its bands are 30,
# 20, 15, 10, 6 and 3 columns), so they bow like small tiles rather than like
# pixels (measured 5.6e-6), and a threshold that admitted them would admit a real
# defect at GLO-30.
@testset "rings are convex, and how far they are from the box" begin
    for sys in (GLO30, GLO90)
        N = CD.lat_intervals(sys)
        g0 = levelgrid(sys, 0)
        g1 = levelgrid(sys, 1)
        reflex = String[]
        worst_tile = 0.0
        worst_pixel = 0.0
        worst_pole_pixel = 0.0

        gap(g, c) = abs(GO.area(MANIFOLD, cell_polygon(g, c)) - cell_area(g, c)) /
                    cell_area(g, c)

        for lat_s in PROBE_LATS, lon_w in PROBE_LONS
            t = CD.tilecell(sys, lat_s, lon_w)
            nc = Int(CD.ncols_at(sys, lat_s))
            has_reflex_vertex(cell_polygon(g0, t)) &&
                note!(reflex, "tile ($lat_s, $lon_w)")
            worst_tile = max(worst_tile, gap(g0, t))
            for (j, i) in ((0, 0), (0, nc - 1), (1, 1), (N ÷ 2, nc ÷ 2),
                           (N - 1, 0), (N - 1, nc - 1))
                c = CD.pixelcell(sys, t, j, i)
                has_reflex_vertex(cell_polygon(g1, c)) &&
                    note!(reflex, "pixel ($lat_s, $lon_w, $j, $i)")
                if lat_s == 89 || lat_s == -90
                    worst_pole_pixel = max(worst_pole_pixel, gap(g1, c))
                else
                    worst_pixel = max(worst_pixel, gap(g1, c))
                end
            end
        end

        @info "$sys ring vs box" worst_tile worst_pixel worst_pole_pixel
        @test reflex == String[]
        @test worst_tile < 1e-4              # measured 5.1e-5, at the ±90 rows
        @test worst_pixel < 1e-8             # measured 8.1e-10 (GLO-30), 3.8e-9 (GLO-90)
        @test worst_pole_pixel < 1e-4        # measured 1.4e-5 (GLO-30), 1.2e-6 (GLO-90)
    end
end

# =========================================================================
# (h) `node_extent` covers the subtree
# =========================================================================

# KILLS: a pad that is too small — the level-0 to level-1 bow really is 1.9e-5
# rad (120 m), so a cap built from corner distances alone under-covers — and a
# centroid computed from the wrong box.
@testset "node_extent covers the subtree" begin
    for sys in ALL_SYSTEMS
        worst_perimeter = -Inf
        max_radius = 0.0
        for lat_s in PROBE_LATS, lon_w in PROBE_LONS
            t = CD.tilecell(sys, lat_s, lon_w)
            cap = node_extent(sys, t)
            max_radius = max(max_radius, cap.radius)
            west, east, south, north = CD.cell_box(sys, t)
            # The tile's own box perimeter, re-sampled 256 points to an edge
            # straight from `cell_box`. This is a TEST-ONLY densification of the
            # continuous truth the cap has to bound; the published ring stays
            # undensified, which is the whole point of the system.
            for k in 0:255
                f = k / 256
                for p in (CD.TO_SPHERE((west + f * (east - west), south)),
                          CD.TO_SPHERE((west + f * (east - west), north)),
                          CD.TO_SPHERE((west, south + f * (north - south))),
                          CD.TO_SPHERE((east, south + f * (north - south))))
                    worst_perimeter = max(worst_perimeter,
                                          US.spherical_distance(cap.point, p) - cap.radius)
                end
            end
        end
        @info "$sys node_extent vs the densely sampled box" worst_perimeter max_radius
        @test worst_perimeter < 0            # strictly inside: the slack is never consumed
        @test max_radius <= π / 2            # convex, as `require_convex_extents` asserts
    end

    # The covering law itself, over EVERY child of every probe tile. Only on the
    # twin: a GLO-30 tile has 12 960 000 children and a GLO-90 tile 1 440 000,
    # and the perimeter probe above is the sharper statement anyway — a child's
    # box vertices all lie on or inside the parent's box, which is the region
    # that probe bounds.
    worst_child = -Inf
    checked = 0
    for lat_s in PROBE_LATS, lon_w in PROBE_LONS
        t = CD.tilecell(TWIN, lat_s, lon_w)
        cap = node_extent(TWIN, t)
        for child in children(TWIN, t), p in cell_boundary(TWIN, child)
            worst_child = max(worst_child, US.spherical_distance(cap.point, p) - cap.radius)
            checked += 1
        end
    end
    @info "twin node_extent vs every child vertex" worst_child checked
    @test checked > 0
    @test worst_child < 0
end

# =========================================================================
# (i) `cellat` on the exact box edges
# =========================================================================

# KILLS: dropping the two repair lines in `cellat` and trusting
# `floor(Int, lat - Δlat/2)`, which is not the inverse of the
# `Float64(lat_s) + Δlat/2` that `cell_box` builds. Measured BEFORE the repair
# landed: the exact-edge probe failed at GLO-30 `lat_s ∈ {4, 16, 64}`, GLO-90
# `{4}` and twin `{1, 16}`; the round-tripped sweep failed at GLO-30 `{16, 64}`
# and GLO-90 `{-4, 4}`. Neither probe catches everything the other does — the
# round-tripped sweep never sees the twin's failures — so both are here.
@testset "cellat agrees with cell_box on south edges" begin
    # PROBE 1: a point whose decoded latitude is EXACTLY the tile's south edge.
    # Not every edge is reachable — `asind` skips values — so the reachable count
    # is pinned rather than assumed: if an upstream change to the coordinate
    # conversion moves it, this goes red instead of quietly testing less.
    for (sys, reachable) in ((GLO30, 120), (GLO90, 120), (TWIN, 125))
        g0 = levelgrid(sys, 0)
        g1 = levelgrid(sys, 1)
        found = 0
        level0_fail = Int[]
        level1_fail = Int[]
        for lat_s in -90:89
            t = CD.tilecell(sys, lat_s, 0)
            south = CD.cell_box(sys, t)[3]
            p = exact_latitude_point(0.0, south)
            p === nothing && continue
            found += 1
            cellat(g0, p) == t || push!(level0_fail, lat_s)
            parent(sys, cellat(g1, p)) == t || push!(level1_fail, lat_s)
        end
        @test found == reachable
        @test level0_fail == Int[]
        @test level1_fail == Int[]
    end

    # PROBE 2: the same 180 south edges with no skips, judged by an independent
    # linear scan of `cell_box` over the latitude that actually SURVIVED the
    # round trip. Complete coverage, at the cost of not always landing on an edge.
    for sys in ALL_SYSTEMS
        g0 = levelgrid(sys, 0)
        bad = String[]
        for lat_s in -90:89
            south = CD.cell_box(sys, CD.tilecell(sys, lat_s, 0))[3]
            p = CD.TO_SPHERE((0.0, south))
            lat = CD.FROM_SPHERE(p)[2]
            owner = nothing
            for ls in -90:89
                b = CD.cell_box(sys, CD.tilecell(sys, ls, 0))
                b[3] <= lat < b[4] && (owner = ls)
            end
            owner === nothing && (note!(bad, "no tile owns $lat"); continue)
            cellat(g0, p) == CD.tilecell(sys, owner, 0) ||
                note!(bad, "south edge of $lat_s decoded to $lat, owned by $owner, " *
                           "but cellat said $(CD.tilecorner(sys, cellat(g0, p)))")
        end
        @test bad == String[]
    end

    # The ANTIMERIDIAN branch, through the real `cellat`. Everything below this
    # is arithmetic that only mirrors `cellat`'s longitude path; these two are
    # the only assertions in the file that drive the `floor(s) >= 180` branch
    # through `cellat` itself. A point a ten-thousandth of a degree west of the
    # antimeridian is inside the E179 tile, but adding `Δlon/2` (1/7200 at
    # GLO-30, 1/2400 at GLO-90, both larger than 1e-4) carries `s` over 180, and
    # the branch sends it round to the W180 tile — the same tile, approached from
    # the east.
    for sys in (GLO30, GLO90)
        g0 = levelgrid(sys, 0)
        @test CD.tilecorner(sys, cellat(g0, CD.TO_SPHERE((179.9999, 0.0)))) == (-1, -180)
    end

    # The mirrored LONGITUDE repair is deliberately ABSENT from `cellat`, and
    # this is the measurement that justifies its absence rather than a repair
    # that never fires: `west` is built as `lon_w - Δlon/2`, and `cellat`'s
    # longitude path — normalise into [-180, 180), add `Δlon/2` back, `floor`,
    # and send a `floor` of 180 round to -180 — recovers `lon_w` for every one of
    # the 3 x 64 800 tiles, with no repair term anywhere.
    #
    # Those last two steps are DIFFERENT lines and this sweep weighs them
    # differently. Deleting the `[-180, 180)` normalisation costs the sweep
    # nothing (measured: 0 misses on all three systems) — every `west` fed in is
    # already in range, and `-180 - Δlon/2` floors to -180 unaided; that line
    # earns its place against the `lon = 180` that `atand` returns at a pole, not
    # against these. Deleting `floor(s) >= 180 -> lon_w = -180` costs exactly 180
    # misses per system, every one of them at `lon_w = -180`: normalisation sends
    # that tile's west edge to `180 - Δlon/2`, and adding the half-pixel back
    # lands on 180.0 exactly. That is the W180 tile reached from the east, which
    # is the branch the two `cellat` assertions above exercise directly.
    #
    # And this is why the sweep is arithmetic rather than a `cellat` probe on the
    # west edge itself: the exact west edge is no more reachable through a
    # `UnitSphericalPoint` than the exact south edge was. Measured, feeding
    # `TO_SPHERE((west, mid))` to `cellat` lands on the WEST neighbour for
    # 10 869 / 10 188 / 11 807 of the 64 800 tiles (GLO-30 / GLO-90 / twin) —
    # a property of the round trip, not of the inversion.
    for sys in ALL_SYSTEMS
        misses = 0
        for lat_s in -90:89
            half_dlon = (1 / CD.ncols_at(sys, lat_s)) / 2
            for lon_w in -180:179
                west = CD.cell_box(sys, CD.tilecell(sys, lat_s, lon_w))[1]
                s = (west - 360 * floor((west + 180) / 360)) + half_dlon
                recovered = floor(s)
                recovered >= 180 && (recovered = -180.0)
                recovered == lon_w || (misses += 1)
            end
        end
        @test misses == 0
    end
end

# =========================================================================
# (j) Cross-resolution nesting: one lattice inside another
# =========================================================================

# What research §5.4 measures is that the two products' pixel CENTRES coincide —
# GLO-90 post `(j, i)` is GLO-30 post `(3j, 3i)`, to 8.9e-16 degrees — and that
# the column counts divide exactly. That, and not a cell-box tiling, is what
# `refine`/`coarsen` implement and what this testset pins: the k x k index block,
# its ascending order, the round trip, and the co-location of the block's
# north-west post with the coarse post.
#
# The cell BOXES are deliberately not asserted to tile the coarse box, because
# they do not. Both products are pixel-is-point, so each outsets its box by half
# of ITS OWN pixel, and the block's box is the coarse box translated south-east,
# on each axis, by `Δ_coarse * (1 - 1/k) / 2` — half a coarse pixel less half a
# fine one, which at k = 3 is one whole GLO-30 pixel. That is a FRACTION of a
# coarse pixel, so the arcsecond figure is per-band (in longitude 1.0" in
# [0,50), 1.5" in [50,60), 10" in [85,90); in latitude 1.0" everywhere), and it
# is the fraction, not the arcseconds, that `worst_shift` asserts below.
#
# The block's areas therefore sum to that translated box (asserted, via
# `worst_union`) and miss the coarse cell's own area by about `tan(φ) * 1"`:
# worst 6.9e-5 at k = 3 and 4.1e-3 at k = 40. In the pole tile rows, where the
# two systems' pole clamps differ by half a pixel each, the miss is 1.78 (+90)
# and 0.40 (-90) at k = 3, and 2.90 / 0.54 at k = 40. Every gap quoted here and
# in `refine`'s docstring is the relative gap `|block - coarse| / coarse`, which
# is the one convention both use. Those gaps are @info-logged, not asserted:
# they are properties of the registration, not tolerances to tighten. See
# `refine`'s docstring for why no uniform, tile-local index scheme fixes them.
#
# KILLS: a half-pixel registration shift between the two products — the
# co-location probe sits ten orders inside the 2.8e-4 degrees the smallest such
# shift moves a post, and measured 2.8e-3 against a `half_dlon = 0` mutant in
# `cell_box`; hardcoding `k = 3`, which the k = 40 pair catches; an off-by-one
# in the k-scaling, which breaks the round trip and the tile-corner extremes
# together; a block anchored anywhere but `(k*j, k*i)` — the "centred" block
# `refine`'s docstring rules out reads a shift of 0 rather than `(1 - 1/k) / 2`,
# and an off-by-one fine column reads `1/k` off, both of which `worst_shift`
# catches eleven orders above its own noise; and a j/i transposition in the
# block base, which lands in a different tile row entirely.
@testset "one lattice nests k-fold inside another" begin
    # One tile per band per hemisphere, plus both pole rows.
    nest_lats = (89, 85, 80, 70, 60, 50, 0, -1, -51, -61, -71, -81, -86, -90)

    # The shipped pair at k = 3, and the twin inside GLO-90 at k = 40 — the case
    # that says the code is written in `k` and not in 3. `k = 40` is also even,
    # so its block has no centre column: the block is an INDEX block, and only
    # the post lattice nests.
    for (coarse, fine) in ((GLO90, GLO30), (TWIN, GLO90))
        k = CD.nesting_factor(coarse, fine)
        Nc = CD.lat_intervals(coarse)
        gc = levelgrid(coarse, 1)
        gf = levelgrid(fine, 1)
        rng = MersenneTwister(20260815)
        bad = String[]
        worst_post = 0.0        # coarse post vs the block's north-west post, degrees
        worst_shift = 0.0       # block box vs the coarse box, in coarse pixels
        worst_union = 0.0       # block area sum vs the block's own box
        worst_coarse = 0.0      # block area sum vs the coarse cell, off the pole rows
        gap_lat50 = 0.0         # the same, in the lat_s = 50 tile row alone
        gap_equator = 0.0       # the same, in the lat_s = 0 tile row alone
        worst_pole_n = 0.0      # the same, in the +90 tile row
        worst_pole_s = 0.0      # the same, in the -90 tile row

        for lat_s in nest_lats
            CD.ncols_at(fine, lat_s) == k * CD.ncols_at(coarse, lat_s) ||
                note!(bad, "row $lat_s: the column counts do not scale by $k")
            nc = CD.ncols_at(coarse, lat_s)
            for lon_w in PROBE_LONS
                tc = CD.tilecell(coarse, lat_s, lon_w)
                tf = CD.tilecell(fine, lat_s, lon_w)
                # Level 0: the tile with the same lower-left corner, both ways.
                CD.refine(coarse, fine, tc) == [tf] ||
                    note!(bad, "tile ($lat_s, $lon_w): refine is not the same corner")
                CD.coarsen(fine, coarse, tf) == tc ||
                    note!(bad, "tile ($lat_s, $lon_w): coarsen is not the same corner")

                for j in unique([0, Nc - 1, rand(rng, 0:(Nc - 1))]),
                    i in unique([0, nc - 1, rand(rng, 0:(nc - 1))])

                    label = "($lat_s, $lon_w) pixel ($j, $i)"
                    p = CD.pixelcell(coarse, tc, j, i)
                    fs = CD.refine(coarse, fine, p)
                    length(fs) == k * k ||
                        note!(bad, "$label: $(length(fs)) cells, not $(k * k)")
                    # Ascending, but NOT contiguous: `k` runs of `k` ids, each run
                    # `ncols(fine, r)` past the last.
                    (issorted(fs) && allunique(fs)) ||
                        note!(bad, "$label: not ascending and distinct")
                    all(f -> CD.coarsen(fine, coarse, f) == p, fs) ||
                        note!(bad, "$label: coarsen does not invert refine")
                    # `parent(fine, f) == tf` for every `f` is NOT re-asserted:
                    # `coarsen` reads its tile straight out of `decode(fine, f)`,
                    # so a block that left its tile could not round-trip back to
                    # `p` above.

                    # The block tiles its own box bitwise — inside a row the east
                    # edge of one cell IS the west edge of the next, and between
                    # rows the south edge of one IS the north edge of the one
                    # below — so the union of the block is just the box spanned
                    # by its north-west and south-east cells, and the areas sum
                    # to it. That shared-edge identity is testset (d)'s, on the
                    # fine lattice itself and for every system; re-checking it on
                    # these blocks tests `cell_box`, not `refine`.
                    boxes = [CD.cell_box(fine, f) for f in fs]
                    west, north = boxes[1][1], boxes[1][4]
                    east, south = boxes[end][2], boxes[end][3]
                    union_area = deg2rad(east - west) * (sind(north) - sind(south))
                    total = sum([cell_area(gf, f) for f in fs])
                    worst_union = max(worst_union, abs(total - union_area) / union_area)

                    cw, ce, cs, cn = CD.cell_box(coarse, p)

                    # The shift the docstring documents: the block's box is the
                    # coarse box translated south-east by `Δ_coarse * (1-1/k)/2`
                    # on each axis. Measured as a FRACTION of the coarse pixel,
                    # which makes the assertion band-free and k-free.
                    #
                    # The pole rows: the longitude form holds there unchanged —
                    # the clamps touch only north and south — and in fact reads
                    # tighter (2.3e-12 vs 2.3e-11), the [85,90) and [-90,-89)
                    # columns being 10" wide. The latitude form does NOT: at
                    # `lat_s = 89, j = 0` both systems clamp their north edge to
                    # exactly 90.0, so the block is not shifted at all there and
                    # the ratio reads 0 against a target of 0.333 (k = 3) or
                    # 0.4875 (k = 40); at `lat_s = -90, j = Nc - 1` the coarse
                    # cell is one and a half pixels tall, so the same ratio comes
                    # out at 2/3 of the target. Those two pixels, and only those
                    # two, are skipped in latitude.
                    target = (1 - 1 / k) / 2
                    worst_shift = max(worst_shift,
                        abs((boxes[1][1] - cw) / (ce - cw) - target))
                    clamped = (lat_s == 89 && j == 0) ||
                              (lat_s == -90 && j == Nc - 1)
                    clamped || (worst_shift = max(worst_shift,
                        abs((cn - boxes[1][4]) / (cn - cs) - target)))

                    gap = abs(total - cell_area(gc, p)) / cell_area(gc, p)
                    if lat_s == 89
                        worst_pole_n = max(worst_pole_n, gap)
                    elseif lat_s == -90
                        worst_pole_s = max(worst_pole_s, gap)
                    else
                        worst_coarse = max(worst_coarse, gap)
                        # The two rows `refine`'s docstring quotes `tan(φ)·1"` at.
                        lat_s == 50 && (gap_lat50 = max(gap_lat50, gap))
                        lat_s == 0 && (gap_equator = max(gap_equator, gap))
                        # The co-location itself, on the box midpoints — which are
                        # the posts, by testset (b). Skipped in the ±90 rows only
                        # because `cell_box` clamps those to the pole by half of
                        # the system's OWN pixel, which moves the midpoint off the
                        # post by a different amount on each side; the posts there
                        # are `lat_n - j/N` on both sides, as everywhere else.
                        worst_post = max(worst_post,
                            abs((cw + ce) / 2 - (boxes[1][1] + boxes[1][2]) / 2),
                            abs((cs + cn) / 2 - (boxes[1][3] + boxes[1][4]) / 2))
                    end
                end
            end
        end

        @info("$coarse inside $fine (k = $k)",
              worst_post, worst_shift, worst_union,
              worst_coarse, gap_lat50, gap_equator, worst_pole_n, worst_pole_s)
        @test bad == String[]
        # Measured 2.8e-14 degrees, i.e. 3 nm on the ground, against the 4.2e-4
        # degrees (1.5") a half-pixel registration slip would move a post.
        @test worst_post < 1e-12
        # The documented south-east shift, as a fraction of a coarse pixel.
        # Measured 2.3e-11 (k = 3) and 4.4e-13 (k = 40): cancellation, not
        # registration. `boxes[1][1] - cw` differences two values near ±180,
        # whose ulp is 2.8e-14, and the shift itself is only 2.8e-4 degrees at
        # k = 3 in the widest band. That noise floor sits eleven orders below
        # what a wrongly anchored block reads — 0.333 (k = 3) / 0.4875 (k = 40)
        # for a block centred on the coarse box, `1/k` for an off-by-one fine
        # column, i.e. 0.025 even at k = 40.
        @test worst_shift < 1e-10
        # Measured 3.6e-16 (k = 3) and, at k = 40, 9.7e-16 standalone against
        # 8.2e-15 under `Pkg.test`'s `--check-bounds=yes`, which costs `sum` its
        # SIMD path and so changes the order of the 1600 additions. The areas are
        # materialised so Julia's pairwise `sum` reduces them, as `cell_area`'s
        # docstring asks: a generator over the same terms lands at 2.6e-14.
        @test worst_union < 1e-12

        # The block decomposition IS the fine tile's raster cut into k x k
        # blocks: the first coarse pixel's block starts at the fine tile's first
        # child, and the last one's ends at its last. The count identity
        # `k^2 * length(children(coarse, tc)) == length(ch)` is not asserted —
        # both sides are `ncols * N` scaled by k, so it restates the per-row
        # `ncols_at(fine) == k * ncols_at(coarse)` check at the top of the loop.
        for lat_s in (0, 50, 89, -90)
            tc = CD.tilecell(coarse, lat_s, 7)
            tf = CD.tilecell(fine, lat_s, 7)
            nc = CD.ncols_at(coarse, lat_s)
            ch = children(fine, tf)
            @test first(CD.refine(coarse, fine, CD.pixelcell(coarse, tc, 0, 0))) == first(ch)
            @test last(CD.refine(coarse, fine, CD.pixelcell(coarse, tc, Nc - 1, nc - 1))) ==
                  last(ch)
        end
    end

    # Not an integer refinement, in both directions and for a lattice that is
    # perfectly valid on its own: `N = 1800` is a legal twin (30 | 1800) but
    # 1800/1200 is 1.5, so no GLO-90 cell is a whole number of its cells.
    t30 = CD.tilecell(GLO30, 0, 0)
    @test_throws ArgumentError CD.refine(GLO30, GLO90, t30)
    @test_throws ArgumentError CD.coarsen(GLO90, GLO30, CD.tilecell(GLO90, 0, 0))
    @test_throws ArgumentError CD.refine(GLO90, CD.CopernicusDEMSystem{1800}(),
                                         CD.tilecell(GLO90, 0, 0))
    # and it says which two lattices it is talking about
    msg = sprint(showerror, try
        CD.refine(GLO30, GLO90, t30)
    catch err
        err
    end)
    @test occursin("3600", msg) && occursin("1200", msg)
end

# =========================================================================
# (k) Contract: the conformance suites
# =========================================================================

# WHY THE TWIN. `test_hierarchical_system` calls `collect(children(sys, c))` for
# every sampled cell, materialises `descendants_at` by recursion, and builds
# `reduce(vcat, collect.(ranges))` over every sibling range (harness `:1610`). A
# GLO-30 tile has up to 12 960 000 children and a GLO-90 tile 1 440 000, and the
# last of those three is quadratic in that count, so the default `n_samples = 8`
# is not affordable on either shipped lattice at any budget: ONE sampled tile of
# GLO-90 measured 143.9 s and 3.36 TiB of allocation, and the note above the
# opt-in runs below says why. `CopernicusDEMSystem{30}()` is the SAME CODE at a
# different `N` — same band table, same six reduction factors, same tile lattice,
# same pole clamps, 900 children per tile — so every law the harness states is
# checked on the code that ships. The real lattices are then checked where they
# differ from the twin, which is arithmetic: by
# `test_grid_interface` at both levels below, by the opt-in hierarchical runs
# below that, and by the fixture, prefix-sum and geometry testsets above.
@testset "conformance (scaled twin)" begin
    test_grid_interface(levelgrid(TWIN, 0); label = "CopernicusDEM twin level 0")
    test_grid_interface(levelgrid(TWIN, 1); label = "CopernicusDEM twin level 1")
    test_hierarchical_system(TWIN; label = "CopernicusDEM twin")

    # And the property that lets the twin run unseeded where the shipped pair
    # cannot: its pole-most rings clear the harness's absolute degeneracy floor.
    # Put to `boundary_problems` itself rather than re-derived from an area, so
    # this is the harness's own verdict on the harness's own threshold. Measured
    # margins: 2.5e-10 sr at +90 and 2.2e-9 at -90, against a floor of 1e-12.
    # Drop the twin's `N` far enough and this goes red HERE, deterministically,
    # instead of turning the run above into a one-draw-in-many flake.
    N = CD.lat_intervals(TWIN)
    for (lat_s, j) in ((89, 0), (-90, N - 1))
        c = CD.pixelcell(TWIN, CD.tilecell(TWIN, lat_s, 0), j, 0)
        @test CT.boundary_problems(cell_boundary(TWIN, c)) == String[]
    end
end

# THE SEED, AND WHAT IT IS FOR. `boundary_problems` (harness `:287-311`) calls a
# ring degenerate when `abs(spherical_signed_area(pts)) <= area_atol`, and
# `area_atol` is an ABSOLUTE steradian tolerance defaulting to `1e-12` that no
# caller can reach: both call sites (`:385`, `:1231`) invoke
# `boundary_problems(pts; unit_atol)`, and neither `test_grid_interface` nor
# `test_hierarchical_system` forwards it. This system's pole-most level-1 rings
# are legitimately smaller than that — a GLO-30 `lat_s = 89` top-row pixel is
# 1/360 degrees wide and 1/7200 tall, i.e. 1.4e-16 sr — so a level-1 draw that
# lands in a ±90 tile row reports a conformance failure on correct geometry.
# Measured, counting raster rows from the pole inward until the ring area clears
# `1e-12`:
#
#   system   N89 top rows   S90 bottom rows   smallest ring area (sr)
#   GLO-30   878            877               1.4e-16 / 1.3e-15
#   GLO-90    33             32               3.8e-15 / 3.5e-14
#   twin       0              0               2.5e-10 / 2.2e-09
#
# That is 3.7e-4 of GLO-30's level 1, so the default draw of 32 fails about one
# seed in eighty. The lever the harness does give a caller is `rng`, and the
# sample is a pure function of it — so the level-1 calls below pass an explicit
# `MersenneTwister` and ASSERT what it buys, by reproducing the harness's own
# `sample_positions` at the same `n_samples` and checking that no sampled cell
# lies in a ±90 tile row at all. (That is the wider property: it implies the
# narrower "no cell in a sub-`1e-12` row", and unlike it, it does not have to be
# restated when the threshold or the geometry moves.) A seed that stops working
# after an upstream change to the sampling then goes red instead of quietly
# re-introducing the flake. The pole rows themselves are NOT sampled away from
# the suite — testset (e) carries them, exhaustively at level 0.
#
# UPSTREAM: `boundary_problems` should judge degeneracy against the ring's own
# scale — `max(area_atol, unit_atol^2)`, say — or take a floor threaded through
# `grid_interface_problems`/`test_grid_interface`, because an absolute steradian
# floor is wrong for any system with legitimately tiny cells and will bite the
# next one too. Not fixed here: the harness is shared by six other systems.
@testset "conformance (GLO-30 and GLO-90)" begin
    for sys in (GLO90, GLO30)
        test_grid_interface(levelgrid(sys, 0); label = "$sys level 0")

        # The seed-away, asserted before it is leaned on. `n_samples` is passed
        # explicitly rather than left to the harness's default so that the draw
        # reproduced here and the draw the call makes are the same draw.
        @test !any(in((89, -90)),
                   sampled_tile_lats(sys, 1, CONFORMANCE_SEED, GI_SAMPLES))
        test_grid_interface(levelgrid(sys, 1); n_samples = GI_SAMPLES,
                            rng = MersenneTwister(CONFORMANCE_SEED),
                            label = "$sys level 1")
    end

    # THE HIERARCHICAL RUNS ARE OPT-IN, and the cost is the harness's, not this
    # system's. `test_hierarchical_system`'s sibling-partition check (`:1610`)
    # evaluates `reduce(vcat, collect.(ranges); init = Int[])` over one range per
    # child; the `init` keyword takes it off `Base`'s `_typed_vcat` fast path and
    # onto `foldl`, which is QUADRATIC in the number of children. Measured in
    # isolation: 2.3 s at 100 000 ranges, 31.6 s at 300 000.
    #
    # So even at `n_samples = 1` this is not cheap. GLO-90 at the seed below
    # draws the tile at (-54, 158), 960 000 children, and MEASURED 143.9 s and
    # 3.36 TiB allocated (55% of it in GC) for 51 passes. GLO-30 draws the same
    # tile — 8 640 000 children, 9x — and the quadratic makes that about 81x the
    # work; it was not run to completion here, and no wall clock is quoted for it
    # because none was observed. Neither belongs in a routine `Pkg.test()`.
    #
    # What a default run therefore loses on the shipped lattices is the covering
    # law under a real descent and the sibling partition. Both are covered:
    # in full on the twin above, and — for the descendant windows specifically —
    # by testset (c), which walks all 64 800 of them EXHAUSTIVELY on all three
    # systems, which is the stronger statement the harness samples one cell of.
    if get(ENV, "DGG_COPDEM_FULL", "0") == "1"
        for sys in (GLO90, GLO30)
            test_hierarchical_system(sys; n_samples = 1,
                                     rng = MersenneTwister(CONFORMANCE_SEED),
                                     label = "CopernicusDEM $sys")
        end
    end
end

end # @testset "CopernicusDEM system"

end # module CopernicusDEMSystemTests
