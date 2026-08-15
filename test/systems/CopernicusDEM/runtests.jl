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
# Every testset names the mutant it kills. A test that kills no mutant no other
# test already kills does not belong here.
#
# The CONTRACT layer — the two conformance suites — is deliberately absent: it
# arrives with Task 5, along with the sampling workaround the harness's absolute
# `area_atol` forces on this system's pole rows.
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

end # @testset "CopernicusDEM system"

end # module CopernicusDEMSystemTests
