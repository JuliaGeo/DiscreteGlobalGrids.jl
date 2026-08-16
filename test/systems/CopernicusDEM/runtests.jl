# ---------------------------------------------------------------------------
# Copernicus DEM system tests.
#
# Four kinds of test, and the distinction matters when one fails:
#
#   1. ORACLE. The lattice itself — which tiles are how many columns wide, and
#      where their pixel centres sit — is checked against 79 real AWS Open Data
#      COGs read with `ArchGDAL` and committed to `fixtures.jl`. No network is
#      touched from `test/`; the fixtures ARE the oracle. They are the only
#      external evidence for the one rule neither primary source states: a
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
#      hierarchical runs behind `DGG_COPDEM_FULL` (`90`, `30` or `all`; off by
#      default). That section carries the seed the harness's absolute
#      `area_atol` forces on this system's pole rows, and the assertion that the
#      seed still does what it claims.
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
# section (l) calls directly: the seed-away assertion there has to reproduce the
# harness's own draw and the harness's own verdict, not an imitation of either.
import DiscreteGlobalGridsConformanceTesting as CT

import GeometryOps as GO
import GeoInterface as GI
# Section (k) checks the block cursor by building a regridder with it and with
# the generic cursor, and demanding the two matrices agree bit for bit.
import ConservativeRegridding as CR
const US = GO.UnitSpherical
using GeometryOps.UnitSpherical: spherical_orient

# The manifold every grid in this package computes on. Named once: a bare vector
# of `UnitSphericalPoint`s carries no manifold, and `best_manifold` would guess
# the WGS84 sphere for it, which is a factor of R^2 in every area.
const MANIFOLD = GO.Spherical(; radius = 1.0)

const GLO30 = DGG.CopernicusDEMSystem(30)
const GLO90 = DGG.CopernicusDEMSystem(90)
# The scaled twin: the same code, the same band table, the same pole clamps, at
# 1/120 the pixel count. Section (l) runs the conformance harness on it; the
# testsets before that use it as a third independent `N` for every law stated in `N`.
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
# first, and which no single-tile test can see — a row/column transposition in
# `tilebase`, which survives any test that only looks at square rows, and a
# `children`/`decode` pair that stops inverting each other at the SHIPPED `N`s.
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

        # children -> parent, AT REAL `N`. `children` numbers a tile's pixels
        # from `tilebase`, `parent` reads the tile back out of an id through
        # `decode`, and both go through the same `N`-dependent prefix sums — so
        # a disagreement between them (a base off by a row, a decode that
        # searches on the wrong side) shows up here and — outside section (k)'s
        # OPT-IN hierarchical runs — in nothing else that starts from a real
        # tile's `children`; testset (i) closes the same loop from `cellat`,
        # which shares the decode but not the base arithmetic. Every
        # band-edge and pole row, at the three probe longitudes in turn, three
        # children each: the range is LAZY, so `ch[k]` is O(1) and the
        # 12 960 000-element child list is never built.
        bad = String[]
        for (i, lat_s) in enumerate(PROBE_LATS)
            lon_w = PROBE_LONS[mod1(i, length(PROBE_LONS))]
            t = CD.tilecell(sys, lat_s, lon_w)
            ch = children(sys, t)
            nch = length(ch)
            for k in (1, (nch + 1) ÷ 2, nch)
                p = parent(sys, ch[k])
                p == t || note!(bad, "child $k/$nch of (lat_s=$lat_s, lon_w=$lon_w) -> $p")
            end
        end
        @test bad == String[]
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
# pixel is a legitimate 1.4e-16 sr. Section (l) samples those rows away, so
# everything asserted about them is asserted here, which is why the level-0
# sweep below is exhaustive rather than sampled.
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
# (`-0.0 == 0.0`) and `===` can. Only the 91 longitudes `0:90`, where both sign
# bits are clear, reproduce the literal bit for bit, so a probe restricted to that
# quadrant would let the mutant live; sweeping all of `lon_w in -180:179` is what
# makes the kill certain.
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
        # decoration: the closed form does not telescope, because not every tile
        # has `east - west === 1.0`, so a generator argument — summed strictly
        # left to right — lands orders further from 4π than the pairwise
        # reduction of the same numbers and fails this `rtol`. Do not loosen the
        # tolerance to accommodate a generator; the areas are fine.
        areas = [cell_area(g0, cellindex(g0, i)) for i in 1:ncells(g0)]
        @test sum(areas) ≈ 4π rtol = 1e-14

        # Latitude coverage is exactly [-90, 90], by `==`: tile row 0's north
        # edge is the clamp and row 179's south edge is the extension.
        @test CD.cell_box(sys, CD.tilecell(sys, 89, 0))[4] === 90.0
        @test CD.cell_box(sys, CD.tilecell(sys, -90, 0))[3] === -90.0
    end

    # One tile per band, its pixels summing to the tile. GLO-90 and the twin
    # only: `cell_area` is one closed form, generic in `N` and in the band's
    # column count, so a third lattice re-runs the same expression on bigger
    # numbers and tests no line the other two leave untested. GLO-30's six
    # tiles would be 36.3 million pixels for that.
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
# The thresholds are two-scoped on purpose. The gap between the box (what
# `cell_area` returns) and the published geodesic quad (what `GO.area` measures)
# is strongly latitude-dependent, and the `lat_s = ±90` rows are half-pixel
# slivers and pole triangles whose gap is four orders above every other row: one
# flat threshold over a sample containing them would fail on geometry that is
# correct.
#
# Run on the two shipped lattices, whose pixel sizes the thresholds are for. The
# twin is deliberately absent: its pixels are 1/30 of a degree tall and from 1/30
# of a degree wide in the 1x band to 1/3 in the 10x band (its bands are 30, 20,
# 15, 10, 6 and 3 columns), so they bow like small tiles rather than like pixels,
# and a threshold loose enough to admit them would admit a real defect at GLO-30.
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
        @test worst_tile < 1e-4              # the ±90 tile rows set this one
        @test worst_pixel < 1e-8             # pixels outside the ±90 tile rows
        @test worst_pole_pixel < 1e-4        # slivers and pole triangles
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
# `Float64(lat_s) + Δlat/2` that `cell_box` builds. The two probes below catch
# that mutant in different tile rows — the round-tripped sweep does not reach the
# twin's — so both are here.
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
    # differently. Deleting the `[-180, 180)` normalisation costs this sweep
    # nothing — every `west` fed in is already in range, and `-180 - Δlon/2`
    # floors to -180 unaided; that line
    # earns its place against the `lon = 180` that `atand` returns at a pole, not
    # against these. Deleting `floor(s) >= 180 -> lon_w = -180` costs exactly 180
    # misses per system, every one of them at `lon_w = -180`: normalisation sends
    # that tile's west edge to `180 - Δlon/2`, and adding the half-pixel back
    # lands on 180.0 exactly. That is the W180 tile reached from the east, which
    # is the branch the two `cellat` assertions above exercise directly.
    #
    # And this is why the sweep is arithmetic rather than a `cellat` probe on the
    # west edge itself: the exact west edge is no more reachable through a
    # `UnitSphericalPoint` than the exact south edge was. Feeding
    # `TO_SPHERE((west, mid))` to `cellat` lands on the WEST neighbour for
    # thousands of the 64 800 tiles on every lattice — a property of the round
    # trip, not of the inversion.
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

# The two products' pixel CENTRES coincide — GLO-90 post `(j, i)` is GLO-30 post
# `(3j, 3i)` — and their column counts divide exactly. That, and not a cell-box
# tiling, is what `refine`/`coarsen` implement and what this testset pins: the
# k x k index block, its ascending order, the round trip, and the co-location of
# the block's north-west post with the coarse post.
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
# `worst_union`) and miss the coarse cell's own area by about `tan(φ) * 1"`. In
# the pole tile rows, where the two systems' pole clamps differ by half a pixel
# each, that miss is of order one instead. Every gap quoted here and in
# `refine`'s docstring is the relative gap `|block - coarse| / coarse`, which is
# the one convention both use. These gaps are @info-logged, not asserted: they
# are properties of the registration, not tolerances to tighten. See `refine`'s
# docstring for why no uniform, tile-local index scheme fixes them.
#
# KILLS: a half-pixel registration shift between the two products — the
# co-location probe sits eight orders inside the 2.8e-4 degrees the smallest
# such shift moves a post; hardcoding `k = 3`, which the k = 40 pair catches;
# an off-by-one
# in the k-scaling, which breaks the round trip and the tile-corner extremes
# together; a block anchored anywhere but `(k*j, k*i)` — the "centred" block
# `refine`'s docstring rules out reads a shift of 0 rather than `(1 - 1/k) / 2`,
# and an off-by-one fine column reads `1/k` off, both far above the 1e-10
# `worst_shift` asserts; and a j/i transposition in the block base, which lands
# in a different tile row entirely.
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
                    # tighter, the [85,90) and [-90,-89) columns being 10" wide.
                    # The latitude form does NOT: at
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
        # A half-pixel registration slip would move a post by 4.2e-4 degrees
        # (1.5"), eight orders above this threshold.
        @test worst_post < 1e-12
        # The documented south-east shift, as a fraction of a coarse pixel. What
        # this threshold has to clear is cancellation rather than registration:
        # `boxes[1][1] - cw` differences two values near ±180, whose ulp is
        # 2.8e-14, while the shift itself is only 2.8e-4 degrees at k = 3 in the
        # widest band. What it has to catch is far above it — a block centred on
        # the coarse box reads 0.333 (k = 3) or 0.4875 (k = 40), an off-by-one
        # fine column reads `1/k`, which is 0.025 even at k = 40.
        @test worst_shift < 1e-10
        # This one has to hold under `Pkg.test`'s `--check-bounds=yes`, which
        # costs `sum` its SIMD path and so reorders the 1600 additions. The areas
        # are materialised so Julia's pairwise `sum` reduces them, as `cell_area`'s
        # docstring asks; a generator over the same terms lands orders worse.
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
# (k) The block cursor: an interior tree over a two-level lattice
# =========================================================================

# WHAT IS AT STAKE. `treeify` hands ConservativeRegridding a cursor, and the dual
# tree search trusts three things about it: a node's cap COVERS every cell
# beneath it (or pairs are silently dropped), a leaf's indices are GRID
# POSITIONS (or the matrix is transposed cell-wise into nonsense), and the
# leaves partition the grid exactly once (or a cell is regridded twice or not at
# all). `CD.BlockCursor` asserts all three by arithmetic on the band tables
# rather than by sampling, so all three are checked here by construction rather
# than by spot check — and then the whole thing is checked at once by demanding
# the intersection matrix be BIT-IDENTICAL to the generic cursor's.
@testset "the block cursor is a tree over the lattice" begin
    STI = GO.SpatialTreeInterface

    # ---- which grids get it -------------------------------------------------
    # KILLS: a `_block_cursor` that accepts a grid whose ids are not one
    # rectangle held as one contiguous run. `position = id - origin` is false
    # there, and every leaf index would be wrong by a growing offset.
    tile90 = CD.tilecell(GLO90, 50, 6)
    rect = PartialGrid(GLO90, tile90, 1)
    ncols90 = CD.ncols_at(GLO90, 50)
    first_id = cellindex(rect, 1).index
    rows = PartialGrid(GLO90, 1, [LevelIndex(1, k)
                                  for k in first_id:(first_id + 4 * ncols90 - 1)])
    midrow = PartialGrid(GLO90, 1, [LevelIndex(1, k)
                                    for k in (first_id + 3):(first_id + 500)])
    scattered = PartialGrid(GLO90, 1, [LevelIndex(1, first_id + 2k) for k in 0:99])
    @test treeify(levelgrid(GLO90, 0)) isa CD.BlockCursor
    @test treeify(levelgrid(GLO90, 1)) isa CD.BlockCursor
    @test treeify(rect) isa CD.BlockCursor
    @test treeify(rows) isa CD.BlockCursor
    # A window that starts mid-row is contiguous in id but is NOT a rectangle,
    # and a scattered list is neither: both fall back, and stay correct.
    @test treeify(midrow) isa DGG.HierarchicalGridCursor
    @test treeify(scattered) isa DGG.HierarchicalGridCursor

    # A level-1 run of WHOLE tiles is a rectangle too — one tile row segment, or
    # whole tile rows, the same two shapes the level-0 rule admits — and it gets
    # the tile-node machinery a level-1 level grid already uses.
    # KILLS: a dispatch that reads "more than one tile" as "not a rectangle" and
    # sends the window down the generic cursor's 64 800-child root; and one that
    # reads it as a rectangle without checking that the END tile is whole, which
    # would claim positions the grid does not hold.
    ntwin = Int(CD.lat_intervals(TWIN))
    nc_twin = Int(CD.ncols_at(TWIN, 50))
    two_lo = CD.pixelcell(TWIN, CD.tilecell(TWIN, 50, 6), 0, 0).index
    two_hi = CD.pixelcell(TWIN, CD.tilecell(TWIN, 50, 7), ntwin - 1, nc_twin - 1).index
    two_tiles = PartialGrid(TWIN, 1, [LevelIndex(1, k) for k in two_lo:two_hi])
    # Whole tile rows, at level 1: the pole-most two TWIN rows are 3 columns to a
    # tile, so this is 64 800 pixels rather than millions.
    pole_lo = CD.pixelcell(TWIN, CD.tilecell(TWIN, 89, -180), 0, 0).index
    pole_hi = CD.pixelcell(TWIN, CD.tilecell(TWIN, 88, 179), ntwin - 1,
                           Int(CD.ncols_at(TWIN, 88)) - 1).index
    pole_rows = PartialGrid(TWIN, 1, [LevelIndex(1, k) for k in pole_lo:pole_hi])
    part_end = PartialGrid(TWIN, 1, [LevelIndex(1, k) for k in two_lo:(two_hi-2)])
    @test treeify(two_tiles) isa CD.BlockCursor
    @test treeify(pole_rows) isa CD.BlockCursor
    @test treeify(part_end) isa DGG.HierarchicalGridCursor

    # An id no cell has. `decode` throws on one and `PartialGrid` does not
    # range-check its ids, while `treeify` is documented to answer for every
    # grid. KILLS: dropping the range check in `_block_cursor`, which turns a
    # fallback into an `ArgumentError` out of `treeify`.
    beyond = ncells(TWIN, 1)
    @test treeify(PartialGrid(TWIN, 1, [LevelIndex(1, beyond + k) for k in 0:3])) isa
          DGG.HierarchicalGridCursor

    # ---- the leaves partition the positions, exactly once each --------------
    # KILLS: an off-by-one in the near-equal split (`_part`), which drops or
    # repeats whole blocks; and a `_position` that forgets the band's `ncols` or
    # the grid's own first id.
    # The counts are kept rather than short-circuited: "some position was wrong"
    # and "position 7 was yielded twice by two different leaves" are different
    # bugs, and the failure message should say which — with the grid and the
    # strategy that produced it, since every case below runs the same assertion.
    function leaf_positions(tree, n)
        seen = falses(n)
        stack = [tree]
        nodes = 0
        dup = 0
        oob = 0
        while !isempty(stack)
            node = pop!(stack)
            nodes += 1
            if STI.isleaf(node)
                for (i, _) in STI.child_indices_extents(node)
                    if !(1 <= i <= n)
                        oob += 1
                    elseif seen[i]
                        dup += 1
                    else
                        seen[i] = true
                    end
                end
            else
                append!(stack, collect(STI.getchild(node)))
            end
        end
        return (nodes = nodes, seen = count(seen), dup = dup, oob = oob)
    end

    twin_tile = PartialGrid(TWIN, CD.tilecell(TWIN, 50, 6), 1)
    for (label, grid) in (("twin tile", twin_tile),
                          ("4 GLO-90 rows", rows),
                          ("GLO-90 tiles", levelgrid(GLO90, 0)),
                          ("two twin tiles", two_tiles),
                          ("two twin tile rows", pole_rows))
        for strategy in (CD.Blocked{3}(), CD.Bisected())
            r = leaf_positions(CD.BlockCursor(grid; strategy), ncells(grid))
            @test (label, string(typeof(strategy)), r.seen, r.dup, r.oob, r.nodes > 1) ==
                  (label, string(typeof(strategy)), ncells(grid), 0, 0, true)
        end
    end

    # ---- the covering law, at every node ------------------------------------
    # KILLS: a `_node_box` that takes one band's half-pixel longitude offset for
    # a block spanning several — the tile lattice's west edges step at latitude
    # 50/60/70/80/85 — or that drops the +90 clamp and the -90 extension. The
    # level-0 globe below crosses every band boundary and both poles.
    function covering_slack(tree, grid)
        worst = -Inf
        stack = [tree]
        while !isempty(stack)
            node = pop!(stack)
            cap = STI.node_extent(node)
            if STI.isleaf(node)
                for (i, _) in STI.child_indices_extents(node),
                    p in cell_boundary(grid, cellindex(grid, i))

                    worst = max(worst, US.spherical_distance(cap.point, p) - cap.radius)
                end
            else
                append!(stack, collect(STI.getchild(node)))
            end
        end
        return worst
    end

    # Both strategies: the caps are derived from a node's rectangle, and the two
    # strategies cut different rectangles out of the same grid — `Blocked{3}`'s
    # edge blocks are the ceil-divided ones, which is where an off-by-one in
    # `_part` would put a cell outside its own node's box.
    for (label, grid) in (("twin tile", twin_tile), ("GLO-90 tiles", levelgrid(GLO90, 0))),
        strategy in (CD.Blocked{3}(), CD.Bisected())

        slack = covering_slack(CD.BlockCursor(grid; strategy), grid)
        @info "block cursor covering slack, $label" strategy slack
        # every leaf vertex strictly inside every ancestor cap
        @test (label, string(typeof(strategy)), slack < 0) ==
              (label, string(typeof(strategy)), true)
    end

    # ---- and the caps against the BOX, not just the corners it was built from --
    # The walk above bottoms out at leaf cells whose ring vertices ARE the box
    # corners the cap radius came from, so its leaf arm can only re-derive
    # `_leaf_pad`. This samples the node's own box perimeter instead — 65 points
    # to an edge, a TEST-ONLY densification of the continuous truth a cap has to
    # bound, as in section (h) — at both scales, in three bands and both pole
    # rows, with the two TILE blocks straddling the band edges at latitude 50 and
    # 85 where the column count steps.
    # KILLS: a `_box_cap` centred anywhere but the box midpoint, and a
    # `_node_box` tile branch that takes one band's half-pixel west offset for a
    # block spanning several.
    g0_90, g1_90 = levelgrid(GLO90, 0), levelgrid(GLO90, 1)
    n90 = Int(CD.lat_intervals(GLO90))
    worst_box = -Inf
    for (lat_s, lon_w, pixels) in ((50, 6, true), (89, 0, true), (-90, 100, true),
                                   (86, -180, true), (51, 0, false), (86, 100, false))
        r, q = 89 - lat_s, lon_w + 180
        nc = Int(CD.ncols(GLO90, r))
        node = pixels ?
               CD.BlockCursor(g1_90, GLO90, CD.Bisected(), 1, Int64(-1), r, r, q, q,
                              0, min(n90, 64) - 1, 0, min(nc, 64) - 1, true) :
               CD.BlockCursor(g0_90, GLO90, CD.Bisected(), 0, Int64(-1),
                              r, min(r + 3, 179), q, min(q + 2, 359), 0, 0, 0, 0, false)
        cap = STI.node_extent(node)
        west, east, south, north = CD._node_box(node)
        for k in 0:64
            f = k / 64
            for p in (CD.TO_SPHERE((west + f * (east - west), south)),
                      CD.TO_SPHERE((west + f * (east - west), north)),
                      CD.TO_SPHERE((west, south + f * (north - south))),
                      CD.TO_SPHERE((east, south + f * (north - south))))
                worst_box = max(worst_box,
                                US.spherical_distance(cap.point, p) - cap.radius)
            end
        end
    end
    @info "block cursor caps vs the densely sampled node box" worst_box
    @test worst_box < 0

    # ---- what `_leaf_pad` is, which is not a bound on the bow ---------------
    # Its docstring says the CORNER sets the cap radius and a bowed edge never
    # reaches it, so the pad is belt-and-braces against rounding rather than the
    # thing that keeps a ring inside. This is that statement as an assertion:
    # rebuild each leaf's cap with pad ZERO and sample its ring's edges, at both
    # scales and in both pole rows.
    # KILLS: a rationale that has drifted from the code. If `cell_boundary` ever
    # densified its rings, or a band's Δλ grew enough for the bow to matter, this
    # goes red and `_leaf_pad`'s docstring needs rewriting.
    worst_vertex = -Inf
    worst_interior = -Inf
    for lat_s in (89, 86, 50, 0, -60, -90), lon_w in (-180, 6)
        r, q = 89 - lat_s, lon_w + 180
        nc = Int(CD.ncols(GLO90, r))
        for (lvl, j, i) in ((1, 0, 0), (1, n90 ÷ 2, nc ÷ 2), (1, n90 - 1, nc - 1),
                            (0, 0, 0))
            grid = lvl == 0 ? g0_90 : g1_90
            node = CD.BlockCursor(grid, GLO90, CD.Bisected(), lvl, Int64(-1),
                                  r, r, q, q, j, j, i, i, lvl == 1)
            west, east, south, north = CD._node_box(node)
            bare = CD._box_cap(west, east, south, north, 0.0)
            c = lvl == 0 ? LevelIndex(0, CD.tileordinal(r, q)) :
                CD.pixelcell(GLO90, LevelIndex(0, CD.tileordinal(r, q)), j, i)
            ring = cell_boundary(grid, c)
            m = length(ring)
            for k in 1:m, t in range(0, 1; length = 65)
                a, b = ring[k], ring[mod1(k + 1, m)]
                x = (1 - t) * a[1] + t * b[1]
                y = (1 - t) * a[2] + t * b[2]
                z = (1 - t) * a[3] + t * b[3]
                nrm = sqrt(x * x + y * y + z * z)
                p = US.UnitSphericalPoint(x / nrm, y / nrm, z / nrm)
                s = US.spherical_distance(bare.point, p) - bare.radius
                # A vertex IS a corner the radius was built from, so its slack is
                # the rounding; the interior of an edge is where the bow lives.
                if t == 0 || t == 1
                    worst_vertex = max(worst_vertex, s)
                else
                    worst_interior = max(worst_interior, s)
                end
            end
        end
    end
    @info "leaf ring vs its UNPADDED corner cap" worst_vertex worst_interior
    @test worst_interior < 0        # the bow never leaves the unpadded cap at all
    @test worst_vertex < 1e-15      # and the corners leave it by float rounding

    # ---- the level-1 TILE node path ----------------------------------------
    # A level-1 grid over more than one tile descends TILE nodes before it ever
    # reaches a raster, and that is a separate arm in three places: `isleaf` (a
    # tile block of a level-1 grid is never a leaf, however few tiles it holds),
    # `_childspace` (a one-tile block descends straight into that tile's raster,
    # spending no level on itself), and `_node_box` (the tile-rectangle branch).
    # KILLS: `isleaf` reading a one-tile level-1 node as a leaf, which would
    # yield one position where a whole raster belongs; and a `_childspace` that
    # descends into the wrong tile's column count.
    g1twin = levelgrid(TWIN, 1)
    root = treeify(g1twin)
    @test root isa CD.BlockCursor
    @test !STI.isleaf(root)
    holds(nd, r, q, j, i) = nd.inpixels ?
                            (nd.r0 == r && nd.q0 == q && nd.j0 <= j <= nd.j1 && nd.i0 <= i <= nd.i1) :
                            (nd.r0 <= r <= nd.r1 && nd.q0 <= q <= nd.q1)
    worst_ancestor = -Inf
    worst_globe = -Inf
    reached = 0
    for (lat_s, lon_w) in ((-90, 0), (89, 179), (50, 6), (49, 6), (86, -180))
        tile = CD.tilecell(TWIN, lat_s, lon_w)
        r, q, _, _ = CD.decode(TWIN, tile)
        nc = Int(CD.ncols(TWIN, r))
        for (j, i) in ((0, 0), (ntwin - 1, nc - 1))
            c = CD.pixelcell(TWIN, tile, j, i)
            pos = cellposition(g1twin, c)
            ring = cell_boundary(g1twin, c)
            node = root
            depth = 0
            while !STI.isleaf(node) && depth < 90
                cap = STI.node_extent(node)
                for p in ring
                    s = US.spherical_distance(cap.point, p) - cap.radius
                    # The root spans the globe, so `_box_cap` takes its POLAR
                    # branch and clamps the radius at π — and the pole opposite
                    # that cap's centre then sits exactly ON it. A cap of radius
                    # π IS the whole sphere, so 0.0 is the right answer there and
                    # only there; every narrower cap must contain the ring
                    # strictly, which is what the pad is for.
                    if cap.radius >= Float64(π)
                        worst_globe = max(worst_globe, s)
                    else
                        worst_ancestor = max(worst_ancestor, s)
                    end
                end
                k = findfirst(k -> holds(STI.getchild(node, k), r, q, j, i),
                              1:STI.nchild(node))
                k === nothing && break
                node = STI.getchild(node, k)
                depth += 1
            end
            STI.isleaf(node) &&
                any(idx == pos for (idx, _) in STI.child_indices_extents(node)) &&
                (reached += 1)
        end
    end
    @info "level-1 tile descent" worst_ancestor worst_globe
    @test reached == 10          # every target found, in its own leaf
    @test worst_ancestor < 0     # and strictly inside every cap on the way down
    @test worst_globe <= 0       # the whole-sphere root included, on its rim

    # ---- the index space ----------------------------------------------------
    # `Trees.getcell(tree, i)` and `child_indices_extents`'s `i` are the same
    # space, and it is the grid's, at every node — the contract the sparse
    # matrix's row and column numbers are read out of.
    root = CD.BlockCursor(twin_tile)
    @test DGG.ncells(root) == ncells(twin_tile)
    for i in (1, 2, ncells(twin_tile) ÷ 3, ncells(twin_tile))
        @test getcell(root, i) == cell_polygon(twin_tile, cellindex(twin_tile, i))
    end

    # ---- and the whole thing at once ----------------------------------------
    # KILLS: any cap that under-covers (pairs vanish), any position that is
    # wrong (entries land in the wrong column), any leaf double-count (entries
    # are added twice). The generic `HierarchicalGridCursor` is the oracle: it
    # descends the SYSTEM's hierarchy and knows nothing about rectangles, so
    # agreeing with it to the last bit is agreeing about the geometry, not about
    # a shared implementation. Both strategies, both scales.
    # Two destinations, because they descend differently: HEALPix bisects a
    # 4-fold hierarchy and IGEO7 a 7-fold one, so the pairs the dual search puts
    # this cursor's nodes against are not the same pairs. And a multi-tile
    # level-1 source, which is the only case whose root is a level-1 TILE node.
    for (label, src, dst) in
        (("twin tile -> HEALPix 5", twin_tile, levelgrid(DGG.HEALPixSystem(), 5)),
         ("twin tile -> IGEO7 4", twin_tile, levelgrid(DGG.IGeo7System(), 4)),
         ("two twin tiles -> HEALPix 5", two_tiles, levelgrid(DGG.HEALPixSystem(), 5)),
         ("GLO-90 tiles -> HEALPix 2", levelgrid(GLO90, 0),
          levelgrid(DGG.HEALPixSystem(), 2)))

        reference = CR.Regridder(MANIFOLD, dst,
            DGG.HierarchicalGridCursor(src)).intersections
        @test length(reference.nzval) > 0
        for strategy in (CD.Blocked{3}(), CD.Bisected())
            blocked = CR.Regridder(MANIFOLD, dst,
                CD.BlockCursor(src; strategy)).intersections
            @test (label, string(typeof(strategy)), blocked == reference) ==
                  (label, string(typeof(strategy)), true)
        end
    end
end

# =========================================================================
# (l) Contract: the conformance suites
# =========================================================================

# WHY THE TWIN. Bare line numbers in this section are into the conformance
# harness, `lib/DiscreteGlobalGridsConformanceTesting/src/DiscreteGlobalGridsConformanceTesting.jl`.
# `test_hierarchical_system` calls `collect(children(sys, c))` for every sampled
# cell, materialises `descendants_at` by recursion, and builds
# `reduce(vcat, collect.(ranges))` over every sibling range (`:1610`). A GLO-30
# tile has up to 12 960 000 children and a GLO-90 tile 1 440 000, and the last of
# those three is quadratic in that count, so the default `n_samples = 8` is not
# affordable on either shipped lattice at any budget: ONE sampled tile of GLO-90
# allocates over 3 TiB, by the closed form the opt-in note below gives.
# `CopernicusDEMSystem{30}()` is the SAME CODE at a different `N` — same band
# table, same six reduction factors, same tile lattice, same pole clamps, 900
# children per tile — so every law the harness states is checked on the code
# that ships. The real lattices are then checked where they differ from the
# twin, which is arithmetic: by `test_grid_interface` at both levels below; by
# the fixture, prefix-sum and geometry testsets above; and — only when
# `DGG_COPDEM_FULL` asks for them, which nothing routine does — by the OPT-IN
# hierarchical runs below that, which a default `Pkg.test()` never reaches and
# which therefore defend nothing on their own.
@testset "conformance (scaled twin)" begin
    test_grid_interface(levelgrid(TWIN, 0); label = "CopernicusDEM twin level 0")
    test_grid_interface(levelgrid(TWIN, 1); label = "CopernicusDEM twin level 1")
    test_hierarchical_system(TWIN; label = "CopernicusDEM twin")

    # And the property that lets the twin run unseeded where the shipped pair
    # cannot: its pole-most rings clear the harness's absolute degeneracy floor.
    # Put to `boundary_problems` itself rather than re-derived from an area, so
    # this is the harness's own verdict on the harness's own threshold.
    # Drop the twin's `N` far enough and this goes red HERE, deterministically,
    # instead of turning the run above into a one-draw-in-many flake.
    N = CD.lat_intervals(TWIN)
    for (lat_s, j) in ((89, 0), (-90, N - 1))
        c = CD.pixelcell(TWIN, CD.tilecell(TWIN, lat_s, 0), j, 0)
        @test CT.boundary_problems(cell_boundary(TWIN, c)) == String[]
    end
end

# THE SEED, AND WHAT IT IS FOR. `boundary_problems` (`:287-311`) calls a
# ring degenerate when `abs(spherical_signed_area(pts)) <= area_atol`, and
# `area_atol` is an ABSOLUTE steradian tolerance defaulting to `1e-12` that no
# caller can reach: both call sites (`:385`, `:1231`) invoke
# `boundary_problems(pts; unit_atol)`, and neither `test_grid_interface` nor
# `test_hierarchical_system` forwards it. This system's pole-most level-1 rings
# are legitimately smaller than that — a GLO-30 `lat_s = 89` top-row pixel is
# 1/360 degrees wide and 1/7200 tall, i.e. 1.4e-16 sr — so a level-1 draw that
# lands in a ±90 tile row reports a conformance failure on correct geometry.
# Counting raster rows from the pole inward until the ring area clears `1e-12`:
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
# ONE PIECE OF EVIDENCE, NOT TWO. The loop below asserts the seed-away on both
# shipped lattices, but at this seed the two draws land on the SAME tile rows —
# the positions differ, the latitudes do not. `rand(rng, 1:n, 32)` turns each
# raw word of the stream into `floor(u * n)` for a `u` that does not depend on
# `n`, and GLO-30's row prefix sums are exactly 9x GLO-90's (3x the columns, 3x
# the rows), so both lattices cut `1:n` into rows at the same fractions of `n`
# and the row sequence is invariant in `N`. So this is one piece of evidence
# read twice, not two: a seed safe for GLO-90 is safe for GLO-30 for that
# reason and not for an independent one, and the second assertion earns its
# place as a guard on the invariance, not as a second sample.
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

        # ...and that the harness really did make THAT draw. The line above
        # predicts the sampled cells from the seed, which is only sound while
        # `test_grid_interface` consumes its `rng` in exactly one
        # `sample_positions(rng, ncells(g), n_samples)` before any law runs
        # (`:1178`). That is a fact about upstream, not a promise, so it
        # is asserted here rather than assumed: `ref` is advanced by the
        # reproduction alone, `r` by the harness alone, and the two generators
        # must end in the same state. One extra draw upstream, a different
        # `n_samples`, or sampling moved after a law that itself draws, and this
        # goes red — instead of the seed-away above quietly guarding cells the
        # suite does not actually visit.
        r = MersenneTwister(CONFORMANCE_SEED)
        ref = MersenneTwister(CONFORMANCE_SEED)
        CT.sample_positions(ref, ncells(levelgrid(sys, 1)), GI_SAMPLES)
        test_grid_interface(levelgrid(sys, 1); n_samples = GI_SAMPLES,
                            rng = r, label = "$sys level 1")
        @test r == ref
    end

    # THE HIERARCHICAL RUNS ARE OPT-IN, and the cost is the harness's, not this
    # system's. `test_hierarchical_system`'s sibling-partition check (`:1610`)
    # evaluates `reduce(vcat, collect.(ranges); init = Int[])` over one range per
    # child; the `init` keyword takes it off `Base`'s `_typed_vcat` fast path and
    # onto `foldl`, which rebuilds the accumulator at every step. The children
    # are at `max_level`, so each of those ranges holds ONE position and the
    # foldl writes 1 + 2 + ... + n `Int`s for n children — `8 * n^2 / 2` bytes,
    # a closed form that needs no machine to state. At n = 960 000 it is
    # 3.353 TiB; the per-range `collect`s add the rest.
    #
    # So even at `n_samples = 1` this is not cheap. GLO-90 at the seed below
    # draws the tile at (-54, 158), 960 000 children, and allocates 3.36 TiB —
    # over half the run in GC — for 51 passes and 2 broken, in minutes rather
    # than seconds. GLO-30 draws the same tile with 8 640 000 children, 9x as
    # many, so the quadratic makes it about 81x the work; that is a ratio and
    # nothing more, since no GLO-30 run is on record as having finished. Neither
    # belongs in a routine `Pkg.test()`.
    #
    # SO THE GATE IS BY VALUE:
    #
    #   unset, or "0"   neither. The default, and what CI runs.
    #   "90"            GLO-90 only: 51 pass, 2 broken; minutes, not seconds.
    #   "1"             the same as "90".
    #   "30"            GLO-30 only. Not known to have completed; see above.
    #   "all"           both.
    #
    # No run here skips anything the harness offers this system. `neighbors` and
    # `ring` are both closed form, so the harness's whole neighbour and ring
    # families run — the shell, concatenation, tail-block and winding laws
    # included. A default `Pkg.test()` reports 17 broken; none of them is a
    # `@test_broken` written here.
    #
    # WHAT A DEFAULT RUN THEREFORE NEVER RUNS ON THE SHIPPED LATTICES, listed in
    # full rather than by the two laws that are easiest to defend: the
    # parent/child inverses (`hierarchy_problems`), the `rootcells` laws,
    # `ancestor` and `descendants`, `node_extent` well-formedness and its
    # convexity proxy, the covering law under a real descent, the sibling
    # partition, the `ArgumentError` guards on `parent`/`children`/`levelgrid`,
    # and the trait block. None of the harness's own versions of those laws runs
    # at `N = 1200` or `N = 3600` unless the gate is set.
    #
    # The defence is not that those laws are minor — it is that they are not
    # about `N`. Of the interface, only `children`, `descendant_range`,
    # `cell_box` and `cellat` have `N` in their method bodies; everything in the
    # list above reaches the lattice through the shared `tables` / `decode` /
    # `pixelcell` helpers and is otherwise `N`-free, so the twin above runs THE
    # SAME CODE, not an analogue of it. What that leaves is the four parametric
    # bodies and those three helpers, which do have to be right at the shipped
    # `N`s — and that is what testsets (a), (c), (d), (h), (i) and (j) do at
    # real `N`: the band table and child counts, the prefix sums exhaustively
    # and the children -> parent round trip, the raster order, `node_extent`
    # over real subtrees, `cellat` on real box edges, and the cross-lattice
    # nesting.
    gate = get(ENV, "DGG_COPDEM_FULL", "0")
    gated = gate in ("0", "")   ? () :
            gate in ("1", "90") ? (GLO90,) :
            gate == "30"        ? (GLO30,) :
            gate == "all"       ? (GLO90, GLO30) :
            throw(ArgumentError("DGG_COPDEM_FULL=$gate is not one of 0, 1, 90, 30, all"))
    for sys in gated
        test_hierarchical_system(sys; n_samples = 1,
                                 rng = MersenneTwister(CONFORMANCE_SEED),
                                 label = "CopernicusDEM $sys")
    end
end

# =========================================================================
# (m) `neighbors` and `ring`: the closed form, and the bounds it forces
# =========================================================================

# WHY THE GENERIC WALK IS NOT THE ORACLE HERE. The walk of
# `src/fallbacks/locate.jl` decides adjacency by counting COINCIDENT RING
# VERTICES within a tolerance. That model cannot see two cells that share a
# boundary SEGMENT without sharing an endpoint, and this lattice makes those at
# every band boundary, where the two sides carry different half-pixel
# registrations. So the walk is a comparator below, not a definition, and the
# oracle is an independent reimplementation of the geometry itself.
fallback_neighbors(g, c, k = 1; connectivity = DGG.Vertex()) =
    invoke(DGG.neighbors, Tuple{DGG.AbstractGrid,DGG.AbstractCellIndex,Integer},
           g, c, k; connectivity)

# Bytes one warmed-up call allocates. Used to assert that no cell reaches a
# spatial tree, without comparing any answers.
_neighbor_bytes(g, c) = (DGG.neighbors(g, c, 1); @allocated DGG.neighbors(g, c, 1))

# -------------------------------------------------------------------------
# THE ORACLE: adjacency from first principles, in exact `Rational{Int}`.
#
# A cell is a longitude interval on one of `180N` latitude rows — `180` tile
# rows at level 0. Consecutive rows abut and rows two apart are a whole row
# apart, so two cells meet exactly where their intervals meet on a shared
# parallel, plus the one point a pole row's triangles share at their apex.
#
# Written against the REGISTRATION `cell_box` documents — a pixel centre plus
# or minus half a column, a tile `ncols` of those — and not against any helper
# `neighbors` uses; and it SCANS a window of candidate columns where the
# implementation solves for one. The id codec is the only thing shared, because
# naming a cell is testset (c)'s subject rather than this one's.
# -------------------------------------------------------------------------

oracle_row(sys, level, J) = level == 1 ? fld(J, CD.lat_intervals(sys)) : J
oracle_nrows(sys, level) = level == 1 ? CD.NROWS * CD.lat_intervals(sys) : CD.NROWS
oracle_width(sys, level, row) = level == 1 ? 360 * CD.ncols(sys, row) : 360

# Cell `k`'s longitude interval, in degrees east of -180, exactly.
function oracle_lon(sys, level, row, k)
    nc = CD.ncols(sys, row)
    return level == 1 ? ((2k - 1)//(2nc), (2k + 1)//(2nc)) :
                        ((2nc * k - 1)//(2nc), (2nc * (k + 1) - 1)//(2nc))
end

# AND ITS REGISTRATION IS `cell_box`'s, checked rather than assumed. The formula
# above writes the registration out a second time, so a wrong one in BOTH would
# survive every comparison this oracle makes. These five cells tie it back to
# `cell_box` — which testset (b) pins against the measured geotransforms — one
# kind at a time: a pixel interior to its tile, the pixel straddling +-180, a
# band-boundary tile, and the two pole cells whose north and south edges are
# corrected. `cell_box` reaches a longitude in three roundings where this is one
# exact rational, so a PIXEL edge can land one Float64 step off; tiles and pole
# cells are bit-exact, and one step is the whole tolerance. Latitudes are not in
# scope here: the oracle decides adjacency on longitude and row index alone.
@testset "the oracle is registered as cell_box is" begin
    N = CD.lat_intervals(TWIN)
    # Representable `Float64` steps apart, which is 0 for a bit-exact match.
    ulps(a, b) = a == b ? 0 :
        signbit(a) == signbit(b) ?
            Int(abs(reinterpret(Int64, a) - reinterpret(Int64, b))) : typemax(Int)
    bad = String[]
    for (label, level, x) in
            (("tile interior pixel", 1, CD.pixelcell(TWIN, CD.tilecell(TWIN, 0, 0), 15, 15)),
             ("antimeridian pixel", 1, CD.pixelcell(TWIN, CD.tilecell(TWIN, 0, -180), 15, 0)),
             ("band-boundary tile", 0, CD.tilecell(TWIN, 50, 0)),
             ("lat_s = 89 sliver pixel", 1, CD.pixelcell(TWIN, CD.tilecell(TWIN, 89, 0), 0, 0)),
             ("lat_s = -90 extension pixel", 1,
              CD.pixelcell(TWIN, CD.tilecell(TWIN, -90, 0), N - 1, 0)))
        r, q, _, i = CD.decode(TWIN, x)
        K = level == 0 ? q : q * CD.ncols(TWIN, r) + i
        a1, a2 = oracle_lon(TWIN, level, r, K)
        west, east = CD.cell_box(TWIN, x)[1:2]
        for (side, got, want) in (("west", west, Float64(-180 + a1)),
                                  ("east", east, Float64(-180 + a2)))
            d = ulps(got, want)
            d <= 1 || note!(bad, "$label $side: $got is $d ulps from $want")
        end
    end
    @test bad == String[]
end

# How two closed arcs of the 360-degree circle meet.
function oracle_meet(a1, a2, b1, b2)
    verdict = :none
    for shift in (-360//1, 0//1, 360//1)
        lo = max(a1, b1 + shift)
        hi = min(a2, b2 + shift)
        hi > lo && return :segment
        hi == lo && (verdict = :point)
    end
    return verdict
end

function oracle_id(sys, level, J, K)
    level == 0 && return DGG.LevelIndex(0, CD.tileordinal(J, K))
    N = CD.lat_intervals(sys)
    r = fld(J, N)
    nc = CD.ncols(sys, r)
    return CD.pixelcell(sys, DGG.LevelIndex(0, CD.tileordinal(r, fld(K, nc))),
                        J - r * N, mod(K, nc))
end

function oracle_neighbors(sys, level, J, K; edge::Bool)
    rows = oracle_nrows(sys, level)
    row = oracle_row(sys, level, J)
    m = oracle_width(sys, level, row)
    a1, a2 = oracle_lon(sys, level, row, K)
    # The two laterals share a whole meridian edge, at every latitude, both
    # levels, and on either side of a pole apex.
    out = Tuple{Int,Int}[(J, mod(K - 1, m)), (J, mod(K + 1, m))]
    if !edge && (J == 0 || J == rows - 1)
        # A pole row's cells are triangles meeting at the exact +-90 point, so
        # the whole row shares that one point and nothing more.
        for l in 0:(m - 1)
            (l == K || l == mod(K - 1, m) || l == mod(K + 1, m)) || push!(out, (J, l))
        end
    end
    for J2 in (J - 1, J + 1)
        0 <= J2 < rows || continue
        row2 = oracle_row(sys, level, J2)
        m2 = oracle_width(sys, level, row2)
        centre = round(Int, (K + 1//2) * m2 // m)
        for l in (centre - 4):(centre + 4)          # scan; do not solve
            l2 = mod(l, m2)
            met = oracle_meet(a1, a2, oracle_lon(sys, level, row2, l2)...)
            (met === :segment || (met === :point && !edge)) && push!(out, (J2, l2))
        end
    end
    return Set(oracle_id(sys, level, p[1], p[2]) for p in unique!(out))
end

# One cell, three laws: the set the geometry says, no duplicates, and — by the
# conformance harness's own measurement, not an imitation of it — that ring 1 is
# a single counter-clockwise cycle about the cell's centroid.
function oracle_problem(sys, level, J, K)
    g = levelgrid(sys, level)
    c = oracle_id(sys, level, J, K)
    for (conn, edge) in ((DGG.Vertex(), false), (DGG.Edge(), true))
        got = DGG.neighbors(g, c, 1; connectivity = conn)
        allunique(got) || return "($J, $K) under $conn: a repeated neighbour"
        want = oracle_neighbors(sys, level, J, K; edge)
        Set(got) == want ||
            return "($J, $K) under $conn: $(length(got)) analytic against $(length(want)) by geometry"
        wound = CT.winding_problems(g, c, got; label = "neighbors($c, 1) under $conn")
        isempty(wound) || return first(wound)
    end
    return nothing
end

# `b in neighbors(a)` iff `a in neighbors(b)`. Pure geometry is symmetric, so
# what this kills is a case gate that fires on one side of a seam and not the
# other. Quadratic in the neighbour count, which is why the pole rings get four
# columns rather than all of them: the ring is one set shared by the whole row,
# and the link a gate could break there is the one to the row below, which each
# of the four exercises.
function symmetry_problem(sys, level, J, K)
    g = levelgrid(sys, level)
    c = oracle_id(sys, level, J, K)
    for conn in (DGG.Vertex(), DGG.Edge())
        for nb in DGG.neighbors(g, c, 1; connectivity = conn)
            c in DGG.neighbors(g, nb, 1; connectivity = conn) ||
                return "($J, $K) under $conn: $nb does not name $c back"
        end
    end
    return nothing
end

@testset "neighbors is closed form over the whole lattice" begin
    N = CD.lat_intervals(TWIN)
    g1 = levelgrid(TWIN, 1)
    g0 = levelgrid(TWIN, 0)

    # THE ANCHOR. Every other check here compares two implementations, so a
    # mutant that moved BOTH would survive all of them. These ids are written
    # out instead of computed: tile `lat_s = 0, lon_w = 0` is tile row 89, tile
    # column 180, thirty columns wide, and the subject is its raster
    # `(j, i) = (15, 15)`. Killed here: any wrong offset, a row and column
    # stride transposed, and the winding started elsewhere or run backwards.
    c = CD.pixelcell(TWIN, CD.tilecell(TWIN, 0, 0), 15, 15)
    @test c == DGG.LevelIndex(1, 21384465)
    @test DGG.neighbors(g1, c, 1) == DGG.LevelIndex.(1,
        [21384434, 21384464, 21384494, 21384495, 21384496, 21384466, 21384436, 21384435])
    @test DGG.neighbors(g1, c, 1; connectivity = DGG.Edge()) ==
          DGG.LevelIndex.(1, [21384435, 21384464, 21384495, 21384466])
    # The same ids as raster positions, so the literals above are legible as the
    # documented NW, W, SW, S, SE, E, NE, N rather than as eight numbers.
    @test [CD.decode(TWIN, x)[3:4] for x in DGG.neighbors(g1, c, 1)] ==
          [(14, 14), (15, 14), (16, 14), (16, 15), (16, 16), (15, 16), (14, 16), (14, 15)]
    @test [CD.decode(TWIN, x)[3:4]
           for x in DGG.neighbors(g1, c, 1; connectivity = DGG.Edge())] ==
          [(14, 15), (15, 14), (16, 15), (15, 16)]

    # NOTHING REACHES A TREE. A tile edge, the antimeridian, a band boundary and
    # a level-0 tile each cost what the tile interior costs. A method that
    # delegated any of them to the generic walk would allocate orders of
    # magnitude more, because the walk builds and queries a spatial tree over
    # the level — so this separates the two structurally, without timing either.
    for x in (c,
              CD.pixelcell(TWIN, CD.tilecell(TWIN, 0, 0), 15, 0),        # west tile edge
              CD.pixelcell(TWIN, CD.tilecell(TWIN, 0, -180), 0, 0),      # antimeridian corner
              CD.pixelcell(TWIN, CD.tilecell(TWIN, 50, 0), N - 1, 5))    # a band boundary
        @test _neighbor_bytes(g1, x) < 1024
    end
    @test _neighbor_bytes(g0, CD.tilecell(TWIN, 50, 0)) < 1024

    # WHY SWEEPING THE TWIN SWEEPS THE SHIPPED LATTICES. Across a boundary the
    # two rows' breakpoints interleave on a pattern fixed by the REDUCED ratio
    # of their column counts — shift the narrow side by `p` columns and the wide
    # side by `q` and both move the same longitude, so the alignment repeats
    # with period `p`, and `p` divides a tile's column count. Every alignment
    # the lattice has therefore appears in one boundary row's columns. And
    # `ncols = 2N / factor`, so that ratio is a ratio of two band factors and is
    # the SAME at every `N`: asserted here rather than argued, so a change to
    # the band table takes this with it. Both terms odd — the only way two
    # corners can coincide across a boundary — happens at +-80 alone.
    reduced(sys) = [(CD.ncols(sys, r), CD.ncols(sys, r + 1)) .÷
                    gcd(CD.ncols(sys, r), CD.ncols(sys, r + 1))
                    for r in 0:(CD.NROWS - 2) if CD.ncols(sys, r) != CD.ncols(sys, r + 1)]
    @test reduced(TWIN) == reduced(GLO30) == reduced(GLO90)
    @test reduced(TWIN) == [(1, 2), (3, 5), (2, 3), (3, 4), (2, 3),
                            (3, 2), (4, 3), (3, 2), (5, 3), (2, 1)]

    # THE BAND BOUNDARIES, every column of both facing rows at all ten of them,
    # in both hemispheres. Exhaustive by the paragraph above.
    bad = String[]
    for r in 0:(CD.NROWS - 2)
        a, b = CD.ncols(TWIN, r), CD.ncols(TWIN, r + 1)
        a == b && continue
        for (J, m) in ((r * N + N - 1, 360a), ((r + 1) * N, 360b)), K in 0:(m - 1)
            p = oracle_problem(TWIN, 1, J, K)
            p === nothing || note!(bad, p)
            p = symmetry_problem(TWIN, 1, J, K)
            p === nothing || note!(bad, p)
        end
    end
    @test bad == String[]

    # BOTH POLE ROWS IN FULL, and the rows facing them. Row 0 is a half-pixel
    # sliver of triangles, row `180N - 1` is one and a half pixels tall, and the
    # oracle gets the same interval machinery for them as for anything else —
    # the slivers are not special-cased on either side of this comparison.
    bad = String[]
    for J in (0, 1, 180N - 2, 180N - 1), K in 0:(360 * CD.ncols(TWIN, fld(J, N)) - 1)
        p = oracle_problem(TWIN, 1, J, K)
        p === nothing || note!(bad, p)
    end
    @test bad == String[]
    bad = String[]
    for J in (0, 180N - 1), K in (0, 1, 180 * CD.ncols(TWIN, fld(J, N)), 360 * CD.ncols(TWIN, fld(J, N)) - 1)
        p = symmetry_problem(TWIN, 1, J, K)
        p === nothing || note!(bad, p)
    end
    @test bad == String[]

    # THE TILE SEAMS AND THE ANTIMERIDIAN, at every probe latitude. Column 0 of
    # a row is the cell that STRADDLES +-180, and a tile's east neighbour is 359
    # tile columns away across it; neither is an offset from anything. Four
    # corners and the edge midpoints of each probe tile.
    bad = String[]
    for lat_s in PROBE_LATS, lon_w in PROBE_LONS
        row = CD._row(lat_s)
        nc = CD.ncols(TWIN, row)
        q = CD._col(lon_w)
        for j in (0, N ÷ 2, N - 1), i in (0, 1, nc ÷ 2, nc - 1)
            p = oracle_problem(TWIN, 1, row * N + j, q * nc + i)
            p === nothing || note!(bad, p)
        end
    end
    @test bad == String[]

    # LEVEL 0, all 64 800 tiles. The tile lattice is small enough to walk
    # exhaustively, and it is its own case: a tile's box is offset half a PIXEL,
    # so the offset MOVES with the band and tiles across a boundary have no
    # corner in common at all.
    bad = String[]
    for r in 0:(CD.NROWS - 1), q in 0:(CD.NCOLS_TILES - 1)
        p = oracle_problem(TWIN, 0, r, q)
        p === nothing || note!(bad, p)
        r in (0, CD.NROWS - 1) && continue          # the pole rings, four tiles each below
        p = symmetry_problem(TWIN, 0, r, q)
        p === nothing || note!(bad, p)
    end
    for r in (0, CD.NROWS - 1), q in (0, 1, 180, 359)
        p = symmetry_problem(TWIN, 0, r, q)
        p === nothing || note!(bad, p)
    end
    @test bad == String[]

    # GLO-90 SPOT CHECK: the same code where the row stride is 40 times wider
    # and the column counts are 40 times larger, so the integers `_facing`
    # multiplies are too. One cell of each kind, against the same oracle.
    N90 = CD.lat_intervals(GLO90)
    bad = String[]
    for (J, K) in ((89N90 + 600, 180 * CD.ncols(GLO90, 89) + 600),        # tile interior
                   (89N90, 180 * CD.ncols(GLO90, 89)),                    # tile NW corner
                   (CD._row(85) * N90 + N90 - 1, 180 * CD.ncols(GLO90, CD._row(85)) + 1),
                   (CD._row(84) * N90, 180 * CD.ncols(GLO90, CD._row(84)) + 3),
                   (0, 5), (180N90 - 1, 7))                               # both pole rows
        p = oracle_problem(GLO90, 1, J, K)
        p === nothing || note!(bad, p)
    end
    @test bad == String[]

    # `ring` AND `neighbors` ARE ONE ANSWER, which is the interface's law:
    # ring 1 IS `neighbors(c, 1)`, the disc is the rings concatenated outward,
    # the shells are disjoint, and ring 2 is wound too.
    # The level-1 pole cell is in the list because its ring 1 is the whole apex
    # row, so ring 2 is the first shell here that the winding has to order with
    # no lattice order to inherit.
    for (g, x) in ((g1, c),
                   (g1, CD.pixelcell(TWIN, CD.tilecell(TWIN, 50, 0), N - 1, 7)),
                   (g1, CD.pixelcell(TWIN, CD.tilecell(TWIN, 89, 0), 0, 0)),
                   (g0, CD.tilecell(TWIN, 50, 0)),
                   (g0, CD.tilecell(TWIN, 89, 0)))
        for conn in (DGG.Vertex(), DGG.Edge())
            r1 = DGG.ring(g, x, 1; connectivity = conn)
            r2 = DGG.ring(g, x, 2; connectivity = conn)
            @test DGG.ring(g, x, 0; connectivity = conn) == [x]
            @test r1 == DGG.neighbors(g, x, 1; connectivity = conn)
            @test DGG.neighbors(g, x, 2; connectivity = conn) == vcat(r1, r2)
            @test isempty(intersect(Set(r1), Set(r2)))
            @test CT.winding_problems(g, x, r2; label = "ring 2") == String[]
        end
        @test DGG.neighbors(g, x, 0) == DGG.LevelIndex[]
    end
    @test_throws ArgumentError DGG.neighbors(g1, c, -1)
    @test_throws ArgumentError DGG.ring(g1, c, -1)
    # A pixel handed to the tile grid is an error, not an answer about tiles.
    @test_throws ArgumentError DGG.neighbors(g0, c, 1)
    @test_throws ArgumentError DGG.neighbors(g1, DGG.LevelIndex(1, DGG.ncells(TWIN, 1)), 1)
end

@testset "where the generic walk parts company" begin
    N = CD.lat_intervals(TWIN)
    g1 = levelgrid(TWIN, 1)
    g0 = levelgrid(TWIN, 0)

    # THE WALK IS ALWAYS A SUBSET, and away from a band boundary it is the whole
    # set. Both halves matter: the first says the closed form never drops a
    # neighbour the vertex-matching model can see, the second says the two
    # disagree in exactly one place and not diffusely.
    inside = String[]
    outside = String[]
    for lat_s in PROBE_LATS, lon_w in PROBE_LONS
        row = CD._row(lat_s)
        nc = CD.ncols(TWIN, row)
        q = CD._col(lon_w)
        for j in (0, N ÷ 2, N - 1), i in (0, nc ÷ 2, nc - 1),
                conn in (DGG.Vertex(), DGG.Edge())
            x = CD.pixelcell(TWIN, DGG.LevelIndex(0, CD.tileordinal(row, q)), j, i)
            got = Set(DGG.neighbors(g1, x, 1; connectivity = conn))
            walk = Set(fallback_neighbors(g1, x, 1; connectivity = conn))
            issubset(walk, got) || note!(outside, "($lat_s, $lon_w, $j, $i) under $conn")
            J = row * N + j
            straddles = (J > 0 && CD.ncols(TWIN, fld(J - 1, N)) != nc) ||
                        (J < 180N - 1 && CD.ncols(TWIN, fld(J + 1, N)) != nc)
            (straddles || got == walk) ||
                note!(inside, "($lat_s, $lon_w, $j, $i) under $conn")
        end
    end
    @test outside == String[]
    @test inside == String[]

    # THE MECHANISM, at a level-0 band boundary tile. `N50_00_E000_00` and the
    # tile below it share more than 0.99 degrees of the parallel half a pixel
    # north of latitude 50 — a segment, not a point — and NO corner, because a
    # tile's box is offset half a pixel and the pixel is 1/20 of a degree above
    # the boundary against 1/30 below. The walk counts coincident ring vertices
    # within a thousandth of the cell's shortest edge; the nearest pair of
    # corners is more than eight times that tolerance apart, so it counts none
    # and reports the five neighbours inside the tile's own band. The lattice
    # has seven.
    tile = CD.tilecell(TWIN, 50, 0)
    below = CD.tilecell(TWIN, 49, 0)
    ring_a = DGG.cell_boundary(TWIN, tile)
    ring_b = DGG.cell_boundary(TWIN, below)
    tol = DGG.Fallbacks._match_tolerance(ring_a)
    gap = minimum(sqrt(sum(abs2, p .- q)) for p in ring_a, q in ring_b)
    @test DGG.Fallbacks._shared_vertices(ring_a, ring_b, tol) == 0
    @test gap > 8 * tol
    west_a, east_a = CD.cell_box(TWIN, tile)[1:2]
    west_b, east_b = CD.cell_box(TWIN, below)[1:2]
    @test min(east_a, east_b) - max(west_a, west_b) > 0.99
    # THE WALKED COUNTS HERE AND IN THE TABLE BELOW PIN A WRONG ANSWER ON
    # PURPOSE: teaching `adjacent_cells` that a shared segment counts would turn
    # them into the closed form's 7, 7, 7 and 8, so red on a walked count after
    # such a fix means DELETE THESE PINS, not that anything regressed.
    @test length(fallback_neighbors(g0, tile)) == 5
    @test length(DGG.neighbors(g0, tile)) == 7
    @test issubset(Set(fallback_neighbors(g0, tile)), Set(DGG.neighbors(g0, tile)))
    # ... and it is the tolerance that decides, given that model: raise it past
    # the corner offset and the same two rings match.
    @test DGG.Fallbacks._shared_vertices(ring_a, ring_b, 1e-3) == 2

    # THE SAME AT LEVEL 1, on both sides of the same boundary, and at +-80 where
    # some corners DO coincide (5:3, the one ratio with two odd terms) so the
    # walk finds part of the answer rather than none of it.
    for (lat_s, j, i, walked, closed) in ((50, N - 1, 5, 5, 7), (49, 0, 5, 5, 7),
                                          (80, N - 1, 2, 7, 8))
        x = CD.pixelcell(TWIN, CD.tilecell(TWIN, lat_s, 0), j, i)
        @test length(fallback_neighbors(g1, x)) == walked
        @test length(DGG.neighbors(g1, x)) == closed
    end
end

@testset "max_neighbors is attained" begin
    # WHAT SETS THE `Vertex()` BOUND, and it is not the lattice interior. Raster
    # row 0 of a `lat_s = 89` tile and row `N - 1` of a `lat_s = -90` tile are
    # triangles sharing one apex — the exact +-90 vertex — with every other cell
    # of their ring, so each is adjacent to the whole ring: `360 * ncols - 1` of
    # them and the three below. WHAT SETS THE `Edge()` BOUND is a band boundary,
    # where a cell on the wide side faces up to three on the narrow side; the
    # 2:1 ratio at +-85 is where three is reached. Both are ATTAINED on all
    # three lattices, which is what makes `==` the assertion — a bound merely
    # respected would survive being raised, and these would not.
    for sys in ALL_SYSTEMS
        N = CD.lat_intervals(sys)
        g = levelgrid(sys, 1)
        for (lat_s, j) in ((89, 0), (-90, N - 1))
            c = CD.pixelcell(sys, CD.tilecell(sys, lat_s, 0), j, 0)
            @test length(DGG.neighbors(g, c)) == DGG.max_neighbors(sys, DGG.Vertex())
            # An apex is one point, so `Edge()` sees only the two laterals and
            # the one overlap equatorward: three, inside the von Neumann four.
            @test length(DGG.neighbors(g, c; connectivity = DGG.Edge())) == 3
        end
        for (lat_s, j) in ((85, N - 1), (-86, 0))       # the wide side of +-85
            t = CD.tilecell(sys, lat_s, 0)
            nc = CD.ncols_at(sys, lat_s)
            @test maximum(length(DGG.neighbors(g, CD.pixelcell(sys, t, j, i);
                                               connectivity = DGG.Edge()))
                          for i in 0:(nc - 1)) == DGG.max_neighbors(sys, DGG.Edge())
        end
        @test DGG.max_neighbors(sys, DGG.Vertex()) == 36 * N + 2
        @test DGG.max_neighbors(sys, DGG.Edge()) == 6
    end

    N = CD.lat_intervals(TWIN)
    g1 = levelgrid(TWIN, 1)
    g0 = levelgrid(TWIN, 0)

    # BAND BOUNDARIES, column by column, because the count there depends on the
    # column: how many of the far side's cells a cell spans is a function of
    # where its edges fall between theirs. Both extremes are pinned, so a mutant
    # that dropped the touch cases or double-counted an overlap moves one of
    # them. Six is the smallest `Vertex()` count any of the ten produces.
    span = Dict{Int,NTuple{4,Int}}()
    for r in 0:(CD.NROWS - 2)
        a, b = CD.ncols(TWIN, r), CD.ncols(TWIN, r + 1)
        a == b && continue
        v = Int[]
        e = Int[]
        for (t, j, n) in ((DGG.LevelIndex(0, CD.tileordinal(r, 180)), N - 1, a),
                          (DGG.LevelIndex(0, CD.tileordinal(r + 1, 180)), 0, b))
            for i in 0:(n - 1)
                x = CD.pixelcell(TWIN, t, j, i)
                push!(v, length(DGG.neighbors(g1, x)))
                push!(e, length(DGG.neighbors(g1, x; connectivity = DGG.Edge())))
            end
        end
        span[CD._lat_s(r)] = (minimum(v), maximum(v), minimum(e), maximum(e))
    end
    @test sort(collect(keys(span))) ==
          [-85, -80, -70, -60, -50, 50, 60, 70, 80, 85]
    @test all(==((6, 8, 4, 6)), values(span))

    # AND NOTHING OFF A POLE ROW COMES NEAR THE `Vertex()` BOUND. The tile edges
    # and corners of every probe latitude, on the antimeridian from both sides
    # and on the prime meridian: the sweep that goes red if a tile edge, a band
    # boundary or the seam could produce a ninth `Vertex()` neighbour or a
    # seventh `Edge()` one.
    bad = String[]
    worst_v = 0
    worst_e = 0
    for lat_s in PROBE_LATS, lon_w in PROBE_LONS
        t = CD.tilecell(TWIN, lat_s, lon_w)
        nc = CD.ncols_at(TWIN, lat_s)
        for j in (0, N - 1), i in (0, nc - 1)
            (lat_s == 89 && j == 0) && continue          # the north pole row
            (lat_s == -90 && j == N - 1) && continue     # the south pole row
            x = CD.pixelcell(TWIN, t, j, i)
            v = length(DGG.neighbors(g1, x))
            e = length(DGG.neighbors(g1, x; connectivity = DGG.Edge()))
            worst_v = max(worst_v, v)
            worst_e = max(worst_e, e)
            (v <= 8 && e <= 6) ||
                note!(bad, "($lat_s, $lon_w) pixel ($j, $i): $v vertex, $e edge")
        end
    end
    @test bad == String[]
    @test worst_v == 8
    @test worst_e == 6

    # LEVEL 0 is the same geometry one level up — 360 tiles to a pole ring
    # rather than `360 * ncols` pixels, and seven neighbours across a band
    # boundary rather than eight, because a tile's half-pixel offset moves with
    # the band and no corner survives the crossing. Swept over every tile row,
    # so this is a claim about the lattice and not about three of its rows.
    counts0 = [(length(DGG.neighbors(g0, DGG.LevelIndex(0, CD.tileordinal(r, 180)))),
                length(DGG.neighbors(g0, DGG.LevelIndex(0, CD.tileordinal(r, 180));
                                     connectivity = DGG.Edge())))
               for r in 0:(CD.NROWS - 1)]
    @test sort(unique(first.(counts0))) == [7, 8, 362]
    @test sort(unique(last.(counts0))) == [3, 4, 5]
    @test maximum(first.(counts0)) < DGG.max_neighbors(TWIN, DGG.Vertex())
    @test maximum(last.(counts0)) < DGG.max_neighbors(TWIN, DGG.Edge())
end

end # @testset "CopernicusDEM system"

end # module CopernicusDEMSystemTests
