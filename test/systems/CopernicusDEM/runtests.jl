# ---------------------------------------------------------------------------
# Copernicus DEM tests use committed COG fixtures and no network.
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
# Section (l) calls these conformance helpers directly.
import DiscreteGlobalGridsConformanceTesting as CT

import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeoInterface as GI
import ConservativeRegridding as CR
# Section (k) checks the tree the regridder gets for one source chunk.
import GlobalRegridding as GR
import Extents
const US = GO.UnitSpherical
using GeometryOps.UnitSpherical: spherical_orient

# A bare point vector would select the WGS84 sphere instead of the unit sphere.
const MANIFOLD = GO.Spherical(; radius = 1.0)

const GLO30 = DGG.CopernicusDEMSystem(30)
const GLO90 = DGG.CopernicusDEMSystem(90)
# Scaled conformance twin.
const TWIN = CD.CopernicusDEMSystem{30}()
const ALL_SYSTEMS = (GLO30, GLO90, TWIN)

include("fixtures.jl")

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

# Tile rows spanning all band edges, the equator, and both poles.
const PROBE_LATS = (89, 88, 85, 84, 80, 79, 70, 69, 60, 59, 50, 49, 0, -1,
                    -50, -51, -60, -61, -70, -71, -80, -81, -85, -86, -89, -90)
# The antimeridian tile from both sides, and the prime meridian.
const PROBE_LONS = (-180, 0, 179)

"A sweep's failures, capped so a failure message names the first few offenders."
note!(bad::Vector{String}, msg) = (length(bad) < 5 && push!(bad, string(msg)); bad)

# Detect reflex turns, ignoring consecutive duplicate vertices.
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
`lat`, or `nothing` when no such point exists. The round trip through the sphere
moves a latitude by up to 10 ulps and skips some values entirely, so this
searches a few ulps of `z` around the naive image.
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

# Section (l) reproduces this exact seeded draw.
const CONFORMANCE_SEED = 20260815
const GI_SAMPLES = 32

"""
The tile-row latitudes `test_grid_interface` will sample: the harness's own
draw, reproduced by calling the harness's own `sample_indices` — so an
upstream change to the sampling moves this too and section (l)'s assertion
goes red instead of silently testing different cells.
"""
function sampled_tile_lats(sys, l, seed, n_samples)
    g = levelgrid(sys, l)
    indices = CT.sample_indices(MersenneTwister(seed), ncells(g), n_samples)
    return [CD.tilecorner(sys, cellindex(g, p))[1] for p in indices]
end

@testset "CopernicusDEM system" begin

# =========================================================================
# (a) The band table, against the measured tiles
# =========================================================================

# Southern boundary rows distinguish equatorward-edge selection from tile labels.
@testset "the band is the equator-ward edge" begin
    for (sys, N, fixtures) in ((GLO30, 3600, GLO30_TILES), (GLO90, 1200, GLO90_TILES))
        for f in fixtures
            @test CD.ncols_at(sys, f.lat_s) == f.ncols
            t = CD.tilecell(sys, f.lat_s, f.lon_w)
            @test length(children(sys, t)) == f.ncols * N
        end
    end

    # Probe both sides of every half-open band boundary.
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
        # The twin uses the same scaled table.
        @test CD.ncols_at(TWIN, lat_s) == cols30 ÷ 120
    end
    # The two rows that touch a pole are both in the 10x band.
    @test CD.ncols_at(GLO30, 89) == 360
    @test CD.ncols_at(GLO30, -90) == 360
end

# =========================================================================
# (b) Registration: the measured geotransforms
# =========================================================================

# Pixel-is-point origins lie half a pixel west and north of the first centre.
@testset "registration reproduces the measured geotransforms" begin
    for f in GEOTRANSFORMS
        sys = CD.CopernicusDEMSystem{f.N}()
        t = CD.tilecell(sys, f.lat_s, f.lon_w)
        west, _, _, north = CD.cell_box(sys, t)
        # Use interior pixels because pole tile boxes include the clamp/extension.
        pw, pe, ps, pn = CD.cell_box(sys, CD.pixelcell(sys, t, 1, 0))
        gt = (west, pe - pw, 0.0, north, 0.0, -(pn - ps))
        expected = (f.origin_x, f.dlon_arcsec / 3600, 0.0, f.origin_y, 0.0, -1 / f.N)
        @test all(abs.(gt .- expected) .<= 1e-12)

        # The first pixel centre equals the filename's integer-degree corner.
        cw, ce, cs, cn = CD.cell_box(sys, CD.pixelcell(sys, t, 0, 0))
        @test (cw + ce) / 2 == Float64(f.lon_w)
        @test (cs + cn) / 2 == Float64(f.lat_s + 1)
    end
end

# =========================================================================
# (c) Ids, prefix sums, and the descendant windows
# =========================================================================

# Check prefix sums, nonsquare indexing, and shipped-resolution round trips.
@testset "ids, prefix sums and descendant windows" begin
    for sys in ALL_SYSTEMS
        N = CD.lat_intervals(sys)
        g0 = levelgrid(sys, 0)
        g1 = levelgrid(sys, 1)

        # All tile windows partition level 1 without gaps or overlaps.
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

        # Indices round-trip through ids at both levels.
        rng = MersenneTwister(20260815)
        for g in (g0, g1)
            n = ncells(g)
            for i in unique([1, 2, n - 1, n, rand(rng, 1:n, 32)...])
                c = cellindex(g, i)
                @test globalindex(g, c) == i
            end
        end

        # Sample lazy children at shipped resolutions and round-trip to parents.
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

# Pin raster order and bit-identical shared edges.
@testset "within-tile order is AWS raster order" begin
    for sys in ALL_SYSTEMS
        N = CD.lat_intervals(sys)
        for lat_s in (0, 50, 60, 70, 80, 85, -1, -51, -86, 89, -90),
            lon_w in (-180, 6, 179)

            t = CD.tilecell(sys, lat_s, lon_w)
            r, q, _, _ = CD.decode(sys, t)
            nc = Int(CD.ncols(sys, r))
            pg = subtree(sys, t, 1)
            @test ncells(pg) == nc * N
            _, _, tile_s, tile_n = CD.cell_box(sys, t)

            for j in (0, 1, N ÷ 2, N - 2, N - 1), i in (0, 1, nc ÷ 2, nc - 2, nc - 1)
                c = CD.pixelcell(sys, t, j, i)
                @test CD.decode(sys, c) == (r, q, j, i)
                # North row first, then west to east.
                k = j * nc + i + 1
                @test cellindex(pg, k) == c
                @test localindex(pg, c) == k

                _, east, south, north = CD.cell_box(sys, c)
                if i < nc - 1
                    # Shared edges must be bit-identical.
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

# Pole rows are tested here because section (l) excludes them from level-1 sampling.
#
# Sweep all level-0 longitudes to catch duplicate vertices and signed-zero poles.
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
            # Every pair in a triangle is consecutive.
            allunique(ring) || note!(bad, "$label: repeated vertex")
            all(p -> abs(hypot(p[1], p[2], p[3]) - 1) < 1e-15, ring) ||
                note!(bad, "$label: not unit norm")
            count(p -> p === pole, ring) == 1 ||
                note!(bad, "$label: pole vertex is not the exact literal")
            GO.area(MANIFOLD, cell_polygon(g, c)) > 0 ||
                note!(bad, "$label: signed area is not positive")
            return bad
        end

        # Exhaust all level-0 pole tiles.
        for (lat_s, pole) in ((89, CD.NORTH_POLE), (-90, CD.SOUTH_POLE)),
            lon_w in -180:179

            check_triangle!(bad, g0, CD.tilecell(sys, lat_s, lon_w), pole,
                            "tile ($lat_s, $lon_w)")
        end

        # Sample the polemost level-1 raster rows.
        rng = MersenneTwister(8990)
        for (lat_s, pole, j) in ((89, CD.NORTH_POLE, 0), (-90, CD.SOUTH_POLE, N - 1))
            nc = Int(CD.ncols_at(sys, lat_s))
            for lon_w in unique(rand(rng, -180:179, 16)), i in (0, nc ÷ 2, nc - 1)
                c = CD.pixelcell(sys, CD.tilecell(sys, lat_s, lon_w), j, i)
                check_triangle!(bad, g1, c, pole, "pixel ($lat_s, $lon_w, $j, $i)")
            end
        end

        # Exactly one raster row per pole tile degenerates; all boxes stay bounded.
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

        # Pole locations must fall within their cell caps.
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

# Global tile areas must sum to 4π.
@testset "the boxes partition the sphere" begin
    for sys in ALL_SYSTEMS
        g0 = levelgrid(sys, 0)
        # Materialise for pairwise summation.
        areas = [cell_area(g0, cellindex(g0, i)) for i in 1:ncells(g0)]
        @test sum(areas) ≈ 4π rtol = 1e-14

        # Pole clamp and extension give exact latitude coverage.
        @test CD.cell_box(sys, CD.tilecell(sys, 89, 0))[4] === 90.0
        @test CD.cell_box(sys, CD.tilecell(sys, -90, 0))[3] === -90.0
    end

    # Pixel areas sum to their tile in every band on affordable lattices.
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

# Pin undensified, counter-clockwise convex rings. Pole rows use a separate bow bound.
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

# The ring moved from a heap `Vector` to inline `Helpers.SmallList` storage, which
# is a representation change and nothing else: the vertices, their order, and both
# pole degeneracies have to be bit-for-bit what the heap version emitted.
@testset "rings are inline, and identical to the heap version" begin
    # Verbatim copy of the pre-inline body. It is the oracle; do not "fix" it.
    function heap_cell_boundary(sys, c)
        west, east, south, north = CD.cell_box(sys, c)
        north == 90.0 && return [CD.TO_SPHERE((west, south)), CD.TO_SPHERE((east, south)),
                                 CD.NORTH_POLE]
        south == -90.0 && return [CD.SOUTH_POLE, CD.TO_SPHERE((east, north)),
                                  CD.TO_SPHERE((west, north))]
        return [CD.TO_SPHERE((west, south)), CD.TO_SPHERE((east, south)),
                CD.TO_SPHERE((east, north)), CD.TO_SPHERE((west, north))]
    end

    differed = String[]
    for sys in ALL_SYSTEMS
        N = CD.lat_intervals(sys)
        for lat_s in PROBE_LATS, lon_w in PROBE_LONS
            t = CD.tilecell(sys, lat_s, lon_w)
            nc = Int(CD.ncols_at(sys, lat_s))
            cells = [CD.pixelcell(sys, t, j, i)
                     for (j, i) in ((0, 0), (0, nc - 1), (1, 1), (N ÷ 2, nc ÷ 2),
                                    (N - 1, 0), (N - 1, nc - 1))]
            pushfirst!(cells, t)
            for c in cells
                got = cell_boundary(sys, c)
                want = heap_cell_boundary(sys, c)
                # `===` on an isbits point compares bits: a signed zero differs.
                (length(got) == length(want) &&
                 all(got[k] === want[k] for k in eachindex(want))) ||
                    note!(differed, "$sys $c: $(collect(got)) != $want")
            end
        end
    end
    @test differed == String[]

    # Inline end to end: the ring, its closed form, and the published polygon.
    g1 = levelgrid(GLO90, 1)
    quad = CD.pixelcell(GLO90, CD.tilecell(GLO90, 0, 0), 3, 3)
    npole = CD.pixelcell(GLO90, CD.tilecell(GLO90, 89, 0), 0, 0)
    spole = CD.pixelcell(GLO90, CD.tilecell(GLO90, -90, 0),
                         CD.lat_intervals(GLO90) - 1, 0)
    for c in (quad, npole, spole)
        ring = cell_boundary(GLO90, c)
        @test ring isa DGG.Helpers.SmallList
        @test isbits(ring)
        @test isbits(DGG.Fallbacks.closed_ring(ring))
        @test isbits(cell_polygon(g1, c))
    end
    # The fourth slot is capacity, not a vertex: a pole cell still reports three.
    @test length(cell_boundary(GLO90, quad)) == 4
    @test length(cell_boundary(GLO90, npole)) == 3
    @test length(cell_boundary(GLO90, spole)) == 3
end

# =========================================================================
# (h) `node_extent` covers the subtree
# =========================================================================

# Caps must cover descendant edge bow as well as corners.
@testset "node_extent covers the subtree" begin
    for sys in ALL_SYSTEMS
        worst_perimeter = -Inf
        max_radius = 0.0
        for lat_s in PROBE_LATS, lon_w in PROBE_LONS
            t = CD.tilecell(sys, lat_s, lon_w)
            cap = node_extent(sys, t)
            max_radius = max(max_radius, cap.radius)
            west, east, south, north = CD.cell_box(sys, t)
            # Densify only the test's continuous reference edge.
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

    # Check every child of each twin probe tile.
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

# Probe Float64 edge cases in the inverse of `cell_box`.
@testset "cellat agrees with cell_box on south edges" begin
    # Probe reachable exact south edges.
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

    # Probe every row after the sphere round trip, using a linear box oracle.
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

    # Half-pixel shifting across +180 wraps into the W180 tile.
    for sys in (GLO30, GLO90)
        g0 = levelgrid(sys, 0)
        @test CD.tilecorner(sys, cellat(g0, CD.TO_SPHERE((179.9999, 0.0)))) == (-1, -180)
    end

    # Longitude inversion needs no mirrored repair on any tile.
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

# Pin block order, round trips, and post co-location. Cell boxes do not nest.
@testset "one lattice nests k-fold inside another" begin
    # One tile per band per hemisphere, plus both pole rows.
    nest_lats = (89, 85, 80, 70, 60, 50, 0, -1, -51, -61, -71, -81, -86, -90)

    # Exercise both shipped `k = 3` and scaled `k = 40` nesting.
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
                # Level 0 preserves the tile corner.
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
                    # Fine pixels are ascending in `k` separated runs.
                    (issorted(fs) && allunique(fs)) ||
                        note!(bad, "$label: not ascending and distinct")
                    all(f -> CD.coarsen(fine, coarse, f) == p, fs) ||
                        note!(bad, "$label: coarsen does not invert refine")

                    # Shared edges make the block union its corner-spanned box.
                    boxes = [CD.cell_box(fine, f) for f in fs]
                    west, north = boxes[1][1], boxes[1][4]
                    east, south = boxes[end][2], boxes[end][3]
                    union_area = deg2rad(east - west) * (sind(north) - sind(south))
                    total = sum([cell_area(gf, f) for f in fs])
                    worst_union = max(worst_union, abs(total - union_area) / union_area)

                    cw, ce, cs, cn = CD.cell_box(coarse, p)

                    # Measure the southeast shift in coarse-pixel units; skip clamped edges.
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
                        lat_s == 50 && (gap_lat50 = max(gap_lat50, gap))
                        lat_s == 0 && (gap_equator = max(gap_equator, gap))
                        # Box midpoints are posts except in clamped pole rows.
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
        # This threshold is well below a half-pixel registration shift.
        @test worst_post < 1e-12
        # Allow antimeridian cancellation while catching index shifts.
        @test worst_shift < 1e-10
        # Materialise areas for pairwise summation.
        @test worst_union < 1e-12

        # The first and last blocks span the fine tile's full raster.
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

    # A valid `N = 1800` lattice does not nest with `N = 1200`.
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
# (j.1) `cellat` on a holding of pixels
# =========================================================================

# A holding is a `PartialGrid` over the lattice, and locating a point in one
# must give the complete lattice's pixel wherever that pixel is held.
_cellat_bytes(g, p) = (cellat(g, p); @allocated cellat(g, p))

@testset "cellat on a holding of pixels" begin
    searched(g, p) = invoke(cellat, Tuple{DGG.AbstractGrid,GO.UnitSphericalPoint}, g, p)
    complete = levelgrid(GLO90, 1)
    tile = CD.tilecell(GLO90, 46, 10)
    nc = Int(CD.ncols_at(GLO90, 46))
    first_id = CD.pixelcell(GLO90, tile, 0, 0).index
    # Four whole raster rows of one tile — the shape a tile holding is built of,
    # and small enough to search cell by cell for the comparison.
    held = [LevelIndex(1, k) for k in first_id:(first_id + 4 * nc - 1)]
    pg = PartialGrid(GLO90, 1, held)

    rng = MersenneTwister(20260825)
    probes = [cell_centroid(complete, rand(rng, held)) for _ in 1:300]
    @test all(p -> cellat(pg, p) === searched(pg, p), probes)
    @test all(p -> cellat(pg, p) === cellat(complete, p), probes)

    # A pixel of the same tile the holding does not reach is outside it, even
    # though the complete lattice names it.
    outside = LevelIndex(1, first_id + 8 * nc)
    @test cellat(complete, cell_centroid(complete, outside)) === outside
    @test cellat(pg, cell_centroid(complete, outside)) === nothing

    # Locating is the lattice arithmetic and a membership search, with no tree
    # query, no candidate list and no boundary polygon behind it.
    @test _cellat_bytes(pg, cell_centroid(complete, held[2 * nc])) == 0 skip =
        VERSION < v"1.12"
end

# =========================================================================
# (k) The block cursor: an interior tree over a two-level lattice
# =========================================================================

# Check cap coverage, index coverage, leaf partitioning, and intersection parity.
@testset "the block cursor is a tree over the lattice" begin
    STI = GO.SpatialTreeInterface

    # ---- which grids get it -------------------------------------------------
    # The complete lattice is one rectangle and gets the block cursor; a holding
    # of tiles, in any arrangement, gets the tiled raster tree.
    tile90 = CD.tilecell(GLO90, 50, 6)
    rect = subtree(GLO90, tile90, 1)
    ncols90 = CD.ncols_at(GLO90, 50)
    first_id = cellindex(rect, 1).index
    rows = PartialGrid(GLO90, 1, [LevelIndex(1, k)
                                  for k in first_id:(first_id + 4 * ncols90 - 1)])
    midrow = PartialGrid(GLO90, 1, [LevelIndex(1, k)
                                    for k in (first_id + 3):(first_id + 500)])
    scattered = PartialGrid(GLO90, 1, [LevelIndex(1, first_id + 2k) for k in 0:99])
    # `treeify` hands the block cursor back memoized; `BlockCursor` is the bare one.
    @test treeify(levelgrid(GLO90, 0)) isa CD.MemoBlockCursor
    @test treeify(levelgrid(GLO90, 1)) isa CD.MemoBlockCursor
    @test treeify(rect) isa DGG.TiledRasterCursor
    @test treeify(rows) isa DGG.TiledRasterCursor
    # Mid-row and scattered holdings are rectangles too, one per run.
    @test treeify(midrow) isa DGG.TiledRasterCursor
    @test treeify(scattered) isa DGG.TiledRasterCursor

    ntwin = Int(CD.lat_intervals(TWIN))
    nc_twin = Int(CD.ncols_at(TWIN, 50))
    two_lo = CD.pixelcell(TWIN, CD.tilecell(TWIN, 50, 6), 0, 0).index
    two_hi = CD.pixelcell(TWIN, CD.tilecell(TWIN, 50, 7), ntwin - 1, nc_twin - 1).index
    two_tiles = PartialGrid(TWIN, 1, [LevelIndex(1, k) for k in two_lo:two_hi])
    # Whole polemost twin tile rows remain affordable at level 1.
    pole_lo = CD.pixelcell(TWIN, CD.tilecell(TWIN, 89, -180), 0, 0).index
    pole_hi = CD.pixelcell(TWIN, CD.tilecell(TWIN, 88, 179), ntwin - 1,
                           Int(CD.ncols_at(TWIN, 88)) - 1).index
    pole_rows = PartialGrid(TWIN, 1, [LevelIndex(1, k) for k in pole_lo:pole_hi])
    part_end = PartialGrid(TWIN, 1, [LevelIndex(1, k) for k in two_lo:(two_hi-2)])
    # Two tiles either side of a latitude row are two id runs, never one
    # rectangle, and are the shape that used to fall to the generic cursor.
    crossrow = PartialGrid(TWIN, 1, sort!(reduce(vcat,
        [collect(children(TWIN, CD.tilecell(TWIN, lat, lon)))
         for (lat, lon) in ((50, 6), (50, 7), (51, 6), (51, 7))])))
    @test treeify(two_tiles) isa DGG.TiledRasterCursor
    @test treeify(pole_rows) isa DGG.TiledRasterCursor
    @test treeify(part_end) isa DGG.TiledRasterCursor
    @test treeify(crossrow) isa DGG.TiledRasterCursor

    # An id this lattice does not name has no rectangle, and keeps the generic
    # cursor rather than throwing where `PartialGrid` promises not to.
    beyond = ncells(TWIN, 1)
    @test treeify(PartialGrid(TWIN, 1, [LevelIndex(1, beyond + k) for k in 0:3])) isa
          DGG.HierarchicalGridCursor

    # ---- the leaves partition the indices, exactly once each ----------------
    # Leaves must cover each grid index exactly once.
    function leaf_indices(tree, n)
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

    twin_tile = subtree(TWIN, CD.tilecell(TWIN, 50, 6), 1)
    for (label, grid) in (("twin tiles", levelgrid(TWIN, 0)),
                          ("GLO-90 tiles", levelgrid(GLO90, 0)))
        for strategy in (CD.Blocked{3}(), CD.Bisected())
            r = leaf_indices(CD.BlockCursor(grid; strategy), ncells(grid))
            @test (label, string(typeof(strategy)), r.seen, r.dup, r.oob, r.nodes > 1) ==
                  (label, string(typeof(strategy)), ncells(grid), 0, 0, true)
        end
    end

    # A raster window of the complete level-1 lattice covers exactly its own
    # index range, which is the shape `subcursor` cuts for a chunk.
    g1twin = levelgrid(TWIN, 1)
    twin_tilecell = CD.tilecell(TWIN, 50, 6)
    twin_range = descendant_range(TWIN, twin_tilecell, 1)
    tr, tq, _, _ = CD.decode(TWIN, twin_tilecell)
    nc_tr = Int(CD.ncols(TWIN, tr))
    function leaf_index_list(tree)
        out = Int[]
        stack = Any[tree]
        while !isempty(stack)
            node = pop!(stack)
            if STI.isleaf(node)
                append!(out, first(e) for e in STI.child_indices_extents(node))
            else
                append!(stack, collect(STI.getchild(node)))
            end
        end
        return sort!(out)
    end
    for strategy in (CD.Blocked{3}(), CD.Bisected())
        window = CD.BlockCursor(g1twin, TWIN, strategy, 1, Int64(-1), tr, tr, tq, tq,
                                0, ntwin - 1, 0, nc_tr - 1, true)
        @test (string(typeof(strategy)), leaf_index_list(window)) ==
              (string(typeof(strategy)), collect(twin_range))
    end
    @test leaf_index_list(DGG.subcursor(g1twin, first(twin_range):last(twin_range))) ==
          collect(twin_range)

    # ---- the covering law, at every node ------------------------------------
    # The global root crosses every band offset and both pole corrections.
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

    # Exercise bisection and uneven `Blocked{3}` edge blocks, at both scales.
    for strategy in (CD.Blocked{3}(), CD.Bisected())
        for (label, tree, grid) in
            (("twin tile raster",
              CD.BlockCursor(g1twin, TWIN, strategy, 1, Int64(-1), tr, tr, tq, tq,
                             0, ntwin - 1, 0, nc_tr - 1, true), g1twin),
             ("GLO-90 tiles", CD.BlockCursor(levelgrid(GLO90, 0); strategy),
              levelgrid(GLO90, 0)))

            slack = covering_slack(tree, grid)
            @info "block cursor covering slack, $label" strategy slack
            # Every leaf vertex lies strictly inside every ancestor cap.
            @test (label, string(typeof(strategy)), slack < 0) ==
                  (label, string(typeof(strategy)), true)
        end
    end

    # ---- and the caps against the BOX, not just the corners it was built from --
    # Densify node-box perimeters at both scales, including band-straddling blocks.
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

    # ---- the corner cap alone contains the ring -----------------------------
    # With zero pad, corner caps must still contain undensified ring edges.
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
                # Vertices test rounding; edge interiors test bow.
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
    # Multi-tile level-1 grids descend through tile nodes into each raster.
    g1twin = levelgrid(TWIN, 1)
    root = treeify(g1twin)
    @test root isa CD.MemoBlockCursor
    @test !STI.isleaf(root)
    holds(t, r, q, j, i) = let nd = t.node
        nd.inpixels ?
            (nd.r0 == r && nd.q0 == q && nd.j0 <= j <= nd.j1 && nd.i0 <= i <= nd.i1) :
            (nd.r0 <= r <= nd.r1 && nd.q0 <= q <= nd.q1)
    end
    worst_ancestor = -Inf
    worst_globe = -Inf
    reached = 0
    for (lat_s, lon_w) in ((-90, 0), (89, 179), (50, 6), (49, 6), (86, -180))
        tile = CD.tilecell(TWIN, lat_s, lon_w)
        r, q, _, _ = CD.decode(TWIN, tile)
        nc = Int(CD.ncols(TWIN, r))
        for (j, i) in ((0, 0), (ntwin - 1, nc - 1))
            c = CD.pixelcell(TWIN, tile, j, i)
            gi = globalindex(g1twin, c)
            ring = cell_boundary(g1twin, c)
            node = root
            depth = 0
            while !STI.isleaf(node) && depth < 90
                cap = STI.node_extent(node)
                for p in ring
                    s = US.spherical_distance(cap.point, p) - cap.radius
                    # Only a radius-π cap may have zero boundary slack.
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
                any(idx == gi for (idx, _) in STI.child_indices_extents(node)) &&
                (reached += 1)
        end
    end
    @info "level-1 tile descent" worst_ancestor worst_globe
    @test reached == 10          # every target found, in its own leaf
    @test worst_ancestor < 0     # and strictly inside every cap on the way down
    @test worst_globe <= 0       # the whole-sphere root included, on its border

    # ---- the index space ----------------------------------------------------
    # Tree and leaf indices both use the grid's own index space.
    root = CD.BlockCursor(levelgrid(TWIN, 0))
    @test DGG.ncells(root) == ncells(TWIN, 0)
    for i in (1, 2, ncells(TWIN, 0) ÷ 3, ncells(TWIN, 0))
        @test getcell(root, i) ==
              cell_polygon(levelgrid(TWIN, 0), cellindex(levelgrid(TWIN, 0), i))
    end

    # ---- and the whole thing at once ----------------------------------------
    # Intersection matrices must match the generic cursor for both destinations.
    for (label, src, dst) in
        (("GLO-90 tiles -> HEALPix 2", levelgrid(GLO90, 0),
          levelgrid(DGG.HEALPixSystem(), 2)),)

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

    # A holding's tree answers the same intersections as the generic cursor.
    for (label, src, dst) in
        (("twin tile -> HEALPix 5", twin_tile, levelgrid(DGG.HEALPixSystem(), 5)),
         ("twin tile -> IGEO7 4", twin_tile, levelgrid(DGG.IGeo7System(), 4)),
         ("two twin tiles -> HEALPix 5", two_tiles, levelgrid(DGG.HEALPixSystem(), 5)))

        reference = CR.Regridder(MANIFOLD, dst,
            DGG.HierarchicalGridCursor(src)).intersections
        @test length(reference.nzval) > 0
        tiled = CR.Regridder(MANIFOLD, dst, treeify(src)).intersections
        @test (label, tiled == reference) == (label, true)
    end
end

# =========================================================================
# (k2) The chunked source: `GR.subtree` windows the cursor
# =========================================================================

# `DGGSpace` chunks a Copernicus source by tile and asks for one tree per chunk.
# `subcursor` answers with the index rectangle rather than the bounding-cap
# fallback, whose row-major bisection prunes nothing on longitude.
@testset "a source chunk keeps the block cursor" begin
    STI = GO.SpatialTreeInterface

    # Every leaf index under `tree`, ascending.
    function leafindices(tree)
        out = Int[]
        stack = Any[tree]
        while !isempty(stack)
            node = pop!(stack)
            if STI.isleaf(node)
                append!(out, first(e) for e in STI.child_indices_extents(node))
            else
                append!(stack, (STI.getchild(node, k) for k in 1:STI.nchild(node)))
            end
        end
        return sort!(out)
    end

    # ---- the complete level grid: 97 M cells, none materialized here --------
    g1twin = levelgrid(TWIN, 1)
    dense = DGG.DGGSpace(g1twin; chunklevel = 0)
    @test GR.nchunks(dense) == 64_800
    # A pole tile, a band edge on both sides, the equator, and the south pole.
    for lat_s in (89, 50, 49, 0, -90), lon_w in (-180, 7)
        tile = CD.tilecell(TWIN, lat_s, lon_w)
        k = GR.chunkat(dense, localindex(g1twin, CD.pixelcell(TWIN, tile, 0, 0)))
        inds = GR.ownedindices(dense, k)
        @test length(inds) == Int(CD.lat_intervals(TWIN)) * Int(CD.ncols_at(TWIN, lat_s))
        tree = GR.subtree(dense, inds)
        @test tree isa CD.MemoBlockCursor
        # The window is the chunk exactly — no cell of another tile leaks in.
        @test leafindices(tree) == collect(inds)
    end

    # Any window of whole raster rows is a rectangle too; a mid-row window is
    # not one, and takes the bounding-cap fallback.
    r = GR.ownedindices(dense, GR.chunkat(dense,
        localindex(g1twin, CD.pixelcell(TWIN, CD.tilecell(TWIN, 12, 40), 0, 0))))
    nc = Int(CD.ncols_at(TWIN, 12))
    rows = (first(r) + nc):(first(r) + 3nc - 1)
    @test GR.subtree(dense, rows) isa CD.MemoBlockCursor
    @test leafindices(GR.subtree(dense, rows)) == collect(rows)
    @test GR.subtree(dense, (first(r) + 1):(first(r) + nc)) isa GR.CellSpaceRTree
    # And a run spanning two tiles is not one rectangle either.
    @test GR.subtree(dense, (last(r) - 1):(last(r) + 1)) isa GR.CellSpaceRTree

    # ---- a sparse holding: a scattered `PartialGrid`, one tile per chunk ----
    tiles = [CD.tilecell(TWIN, 89, 10),     # polemost band, 3 columns
             CD.tilecell(TWIN, 50, -3),     # band edge, west of the prime meridian
             CD.tilecell(TWIN, 12, 40)]     # full-width equatorial tile
    ids = sort!(reduce(vcat, [collect(children(TWIN, t)) for t in tiles]))
    holding = PartialGrid(TWIN, 1, ids)
    @test treeify(holding) isa DGG.TiledRasterCursor
    src = DGG.DGGSpace(holding; chunklevel = 0)
    @test GR.nchunks(src) == length(tiles)

    # One destination covering the three tiles, at a level fine enough that most
    # of its cells take weight from the source.
    tilebox(t) = let (lat_s, lon_w) = CD.tilecorner(TWIN, t)
        Extents.Extent(X = (Float64(lon_w), lon_w + 1.0), Y = (Float64(lat_s), lat_s + 1.0))
    end
    sys7 = DGG.IGeo7System()
    dstgrid = PartialGrid(sys7, 7, sort!(unique(reduce(vcat,
        [query(sys7, Intersects(tilebox(t)); level = 7) for t in tiles]))))
    dsttree = GR.subtree(DGG.DGGSpace(dstgrid), 1:ncells(dstgrid))

    op = CR.DefaultIntersectionOperator(MANIFOLD)
    for k in 1:GR.nchunks(src)
        inds = GR.ownedindices(src, k)
        tree = GR.subtree(src, inds)
        @test tree isa DGG.TiledRasterCursor
        @test leafindices(tree) == collect(inds)
        # The weights are the fallback's, entry for entry.
        fast = CR.intersection_areas(MANIFOLD, GOCore.False(), dsttree, tree;
            intersection_operator = op)
        slow = CR.intersection_areas(MANIFOLD, GOCore.False(), dsttree,
            GR.CellSpaceRTree(src, inds); intersection_operator = op)
        @test length(fast.nzval) > 0
        @test (k, fast == slow) == (k, true)
    end
end

# =========================================================================
# (k2b) The tiled raster tree over a holding of tiles
# =========================================================================

# A holding is a collection of tiles in no particular arrangement. Its tree
# packs the tiles by their caps and bisects each tile's raster beneath, and its
# leaf index is the tile's offset in the grid plus the pixel's row-major index
# within the tile.
@testset "the tiled raster tree over a holding" begin
    STI = GO.SpatialTreeInterface
    capholds(cap, p) = US.spherical_distance(cap.point, p) <= cap.radius

    holding(spec) = PartialGrid(TWIN, 1, sort!(reduce(vcat,
        [collect(children(TWIN, CD.tilecell(TWIN, lat, lon))) for (lat, lon) in spec])))

    onetile = [(50, 6)]
    onerow = [(50, 6), (50, 7)]
    square = [(50, 6), (50, 7), (51, 6), (51, 7)]
    ragged = [(50, 6), (50, 7), (51, 6)]                 # an L, not a rectangle
    scattered = [(89, 10), (50, -3), (12, 40), (-90, 100)]

    # ---- the index law -----------------------------------------------------
    # A leaf's index is the grid's own, for every pixel of every shape.
    for spec in (onetile, onerow, square, ragged, scattered)
        g = holding(spec)
        n = ncells(g)
        seen = falses(n)
        dup = oob = badindex = badcap = 0
        stack = Any[treeify(g)]
        while !isempty(stack)
            node = pop!(stack)
            if STI.isleaf(node)
                for (i, cap) in STI.child_indices_extents(node)
                    if !(1 <= i <= n)
                        oob += 1
                        continue
                    end
                    seen[i] ? (dup += 1) : (seen[i] = true)
                    c = cellindex(g, i)
                    localindex(g, c) == i || (badindex += 1)
                    all(capholds(cap, p) for p in cell_boundary(g, c)) ||
                        (badcap += 1)
                end
            else
                append!(stack, collect(STI.getchild(node)))
            end
        end
        @test (length(spec), count(seen), dup, oob, badindex, badcap) ==
              (length(spec), n, 0, 0, 0, 0)
    end

    # And it is the stated law: the tile's offset plus the pixel's row-major
    # index within the tile, for a random pixel of a random tile.
    rng = MersenneTwister(20260825)
    for spec in (square, ragged)
        g = holding(spec)
        for _ in 1:200
            lat, lon = spec[rand(rng, eachindex(spec))]
            tile = CD.tilecell(TWIN, lat, lon)
            r, q, _, _ = CD.decode(TWIN, tile)
            nc = Int(CD.ncols(TWIN, r))
            j, i = rand(rng, 0:(Int(CD.lat_intervals(TWIN)) - 1)), rand(rng, 0:(nc - 1))
            offset = localindex(g, CD.pixelcell(TWIN, tile, 0, 0)) - 1
            @test localindex(g, CD.pixelcell(TWIN, tile, j, i)) == offset + j * nc + i + 1
        end
    end

    # ---- the extents nest at the tile layer --------------------------------
    # Every packed node's cap is a merge of its children's, so a tile's cap is
    # inside every cap above it.
    let g = holding(scattered)
        tree = treeify(g)
        worst = -Inf
        stack = Any[tree]
        while !isempty(stack)
            node = pop!(stack)
            node.inraster && continue
            cap = STI.node_extent(node)
            for k in 1:STI.nchild(node)
                child = STI.getchild(node, k)
                child.inraster && continue
                cc = STI.node_extent(child)
                worst = max(worst,
                    US.spherical_distance(cap.point, cc.point) + cc.radius - cap.radius)
                push!(stack, child)
            end
        end
        @test worst <= 0
    end

    # ---- the same answers the block cursor gave -----------------------------
    # A holding's leaf caps are the complete lattice's, so a point query over
    # the holding names exactly what the per-tile block cursors name.
    complete = levelgrid(TWIN, 1)
    function blockcursor_hits(g, spec, p)
        out = Int[]
        for (lat, lon) in spec
            r = descendant_range(TWIN, CD.tilecell(TWIN, lat, lon), 1)
            window = DGG.subcursor(complete, first(r):last(r))
            for i in STI.query(window, cap -> capholds(cap, p))
                push!(out, localindex(g, cellindex(complete, i)))
            end
        end
        return sort!(out)
    end
    for spec in (onetile, onerow, square)
        g = holding(spec)
        tree = treeify(g)
        probes = GO.UnitSphericalPoint{Float64}[]
        for _ in 1:60
            c = cellindex(g, rand(rng, 1:ncells(g)))
            push!(probes, cell_centroid(g, c))
            append!(probes, cell_boundary(g, c))       # edges and corners
        end
        bad = 0
        for p in probes
            sort!(collect(STI.query(tree, cap -> capholds(cap, p)))) ==
                blockcursor_hits(g, spec, p) || (bad += 1)
        end
        @test (length(spec), length(probes) >= 300, bad) == (length(spec), true, 0)
    end

    # ---- the interface the regridder reads ---------------------------------
    let g = holding(square), tree = treeify(g)
        @test STI.isspatialtree(typeof(tree))
        # Tile caps are stored and raster caps memoized, so a repeat ask is a load.
        @test !STI.node_extent_is_expensive(typeof(tree))
        @test !STI.isleaf(tree)
        @test STI.nchild(tree) > 1
        @test [STI.getchild(tree, k) for k in 1:STI.nchild(tree)] ==
              collect(STI.getchild(tree))
        @test_throws BoundsError STI.getchild(tree, STI.nchild(tree) + 1)
        @test_throws ArgumentError STI.child_indices_extents(tree)
        @test DGG.ncells(tree) == ncells(g)
        @test getcell(tree, 7) == cell_polygon(g, cellindex(g, 7))
        @test CR.Trees.split_weight(tree) == ncells(g)
        @test GOCore.best_manifold(tree) == GOCore.best_manifold(g)
        @test treeify(tree) === tree
        # A node's weight is the pixels beneath it, and the children partition them.
        @test sum(CR.Trees.split_weight(c) for c in STI.getchild(tree)) ==
              CR.Trees.split_weight(tree)
        # The memo answers what the hooks answer, whoever asks and in what order.
        node = STI.getchild(STI.getchild(tree, 1), 1)
        want = STI.node_extent(node)
        @test STI.node_extent(node) === want
        @test all(==(0), fetch.([Threads.@spawn(count(_ -> STI.node_extent(node) !== want,
                                                     1:64)) for _ in 1:4]))
    end

    # ---- the hot path infers ------------------------------------------------
    # A search descends, reads extents and reaches leaves without a dynamic
    # dispatch at any node, so every accessor and every hook it calls has a
    # concrete return type. Construction may still return a union.
    let g = holding(square), tree = treeify(g)
        packed = STI.getchild(tree, 1)        # one tile, still a packed node
        raster = STI.getchild(packed, 1)      # a rectangle inside that tile
        leaf = raster
        while !STI.isleaf(leaf)
            leaf = STI.getchild(leaf, 1)
        end
        for node in (tree, packed, raster, leaf)
            @test @inferred(STI.isleaf(node)) isa Bool
            @test @inferred(STI.nchild(node)) isa Int
            @test @inferred(CR.Trees.split_weight(node)) isa Int
            @test @inferred(STI.node_extent(node)) isa US.SphericalCap
        end
        for node in (tree, packed, raster)
            @test @inferred(STI.getchild(node, 1)) isa DGG.TiledRasterCursor
        end
        @test @inferred(STI.child_indices_extents(leaf)) isa DGG.Engine.LeafCells

        # And the hooks the tree calls at every node.
        tile = first(DGG.raster_tiles(g, 1:ncells(g)))
        @test @inferred(DGG.raster_shape(g, tile)) isa Tuple{Int,Int}
        @test @inferred(DGG.raster_localindex(g, tile, 1, 2)) isa Int
        @test @inferred(DGG.raster_cap(g, tile, 0, 1, 0, 1)) isa US.SphericalCap
        @test @inferred(DGG.Engine.rect_part(0, 7, 2, 1)) isa Tuple{Int,Int}
        @test @inferred(DGG.Engine.bisect_parts(4, 7)) isa Tuple{Int,Int}
        cap = DGG.raster_cap(g, tile, 0, 0, 0, 0)
        @test @inferred(DGG.Engine.leaf_cells(k -> (k, cap), 3)) isa
              DGG.Engine.LeafCells
    end

    # An empty holding is a tree with no cells rather than an error.
    let empty_tree = treeify(PartialGrid(TWIN, 1, LevelIndex[]))
        @test empty_tree isa DGG.TiledRasterCursor
        @test STI.isleaf(empty_tree)
        @test isempty(STI.child_indices_extents(empty_tree))
    end

    # ---- a level-0 holding: the cell is the tile ---------------------------
    let g = PartialGrid(TWIN, 0, sort!([CD.tilecell(TWIN, lat, lon)
                                        for (lat, lon) in scattered]))
        tree = treeify(g)
        @test tree isa DGG.TiledRasterCursor
        hits = Int[]
        stack = Any[tree]
        while !isempty(stack)
            node = pop!(stack)
            STI.isleaf(node) ?
                append!(hits, first(e) for e in STI.child_indices_extents(node)) :
                append!(stack, collect(STI.getchild(node)))
        end
        @test sort!(hits) == collect(1:ncells(g))
    end

    # ---- the weights a regrid builds are the generic cursor's ---------------
    for spec in (square, ragged)
        g = holding(spec)
        dst = levelgrid(DGG.HEALPixSystem(), 5)
        reference = CR.Regridder(MANIFOLD, dst,
            DGG.HierarchicalGridCursor(g)).intersections
        tiled = CR.Regridder(MANIFOLD, dst, treeify(g)).intersections
        @test length(reference.nzval) > 0
        @test (length(spec), tiled == reference) == (length(spec), true)
    end
end

# =========================================================================
# (k3) The extent memo: the same caps, no false hits, no shared state
# =========================================================================

# `treeify`/`subcursor` hand back a `MemoBlockCursor`, which answers
# `node_extent` from a per-task direct-mapped table instead of re-deriving the
# box and its cap. The table must be invisible: every answer is the bare
# cursor's, bit for bit, whoever asks and in whatever order.
@testset "the node-extent memo" begin
    STI = GO.SpatialTreeInterface

    # The block cursor addresses the complete lattice and the windows
    # `subcursor` cuts out of it; one tile's raster is such a window.
    g1twin = levelgrid(TWIN, 1)
    tilerange = descendant_range(TWIN, CD.tilecell(TWIN, 50, 6), 1)
    memoroot = DGG.subcursor(g1twin, first(tilerange):last(tilerange))
    @test memoroot isa CD.MemoBlockCursor
    bareroot = memoroot.node

    # ---- the interface the bare cursor implements, all of it ----------------
    @test STI.isspatialtree(typeof(memoroot))
    # A hit is a compare and a load, so the search must not cache child extents.
    @test !STI.node_extent_is_expensive(typeof(memoroot))
    @test STI.isleaf(memoroot) == STI.isleaf(bareroot)
    @test STI.nchild(memoroot) == STI.nchild(bareroot)
    @test all(c isa CD.MemoBlockCursor for c in STI.getchild(memoroot))
    @test [c.node for c in STI.getchild(memoroot)] ==
          [STI.getchild(bareroot, k) for k in 1:STI.nchild(bareroot)]
    @test STI.getchild(memoroot, 2).node == STI.getchild(bareroot, 2)
    @test DGG.ncells(memoroot) == DGG.ncells(bareroot)
    @test getcell(memoroot, first(tilerange)) == getcell(bareroot, first(tilerange))
    @test CR.Trees.split_weight(memoroot) == CR.Trees.split_weight(bareroot)
    @test GOCore.best_manifold(memoroot) == GOCore.best_manifold(bareroot)
    @test treeify(memoroot) === memoroot
    @test sprint(show, memoroot) == "Memo" * sprint(show, bareroot)

    # ---- (a) every answer is the unmemoized one, bit for bit ---------------
    # `===` on an immutable cap compares the bits, so this is stricter than `==`.
    function sameextents(a, b)
        STI.node_extent(a) === STI.node_extent(b) || return false
        STI.isleaf(a) == STI.isleaf(b) || return false
        STI.isleaf(a) && return isequal(collect(STI.child_indices_extents(a)),
                                        collect(STI.child_indices_extents(b)))
        STI.nchild(a) == STI.nchild(b) || return false
        return all(sameextents(STI.getchild(a, k), STI.getchild(b, k))
                   for k in 1:STI.nchild(a))
    end
    @test sameextents(memoroot, bareroot)
    # A second walk is served from the table and answers the same.
    @test sameextents(memoroot, bareroot)

    # The root's slot really does hold the root's cap after the first ask.
    rootkey = CD._nodekey(bareroot)
    slot = CD._nodeslot(rootkey, CD._MEMO_EXTENT_SLOTS)
    rootcap = STI.node_extent(memoroot)
    thememo = CD._taskmemo(bareroot)
    @test thememo.keys[slot] == rootkey
    @test thememo.vals[slot] === rootcap === STI.node_extent(bareroot)

    # ---- (b) a slot collision is a miss, never a false hit ------------------
    # The whole level-1 lattice, breadth-first, past the point where the table
    # is full: with more live nodes than slots, collisions are the common case.
    globe = treeify(levelgrid(TWIN, 1))
    nodes = CD.MemoBlockCursor[globe]
    k = 1
    while k <= length(nodes) && length(nodes) < 4 * CD._MEMO_EXTENT_SLOTS
        n = nodes[k]
        k += 1
        STI.isleaf(n) || append!(nodes, (STI.getchild(n, c) for c in 1:STI.nchild(n)))
    end
    @test length(nodes) > CD._MEMO_EXTENT_SLOTS   # more nodes than slots: real pressure

    # Two distinct nodes, with distinct caps, that land in one slot.
    seen = Dict{Int,CD.MemoBlockCursor}()
    collide = nothing
    for n in nodes
        s = CD._nodeslot(CD._nodekey(n.node), CD._MEMO_EXTENT_SLOTS)
        prev = get(seen, s, nothing)
        if prev !== nothing && CD._nodekey(prev.node) != CD._nodekey(n.node) &&
           STI.node_extent(prev.node) !== STI.node_extent(n.node)
            collide = (prev, n)
            break
        end
        seen[s] = n
    end
    @test collide !== nothing
    a, b = collide
    @test CD._nodeslot(CD._nodekey(a.node), CD._MEMO_EXTENT_SLOTS) ==
          CD._nodeslot(CD._nodekey(b.node), CD._MEMO_EXTENT_SLOTS)
    # Alternating evicts on every ask; each ask still gets its own node's cap.
    for _ in 1:4
        @test STI.node_extent(a) === STI.node_extent(a.node)
        @test STI.node_extent(b) === STI.node_extent(b.node)
    end
    @test STI.node_extent(a) !== STI.node_extent(b)

    # One tile rectangle is one key, but the cap depends on the lattice it sits
    # in: the level-0 pad is a whole tile, the level-1 pad one pixel. The table
    # is cleared when a task turns to another lattice, so neither answers for
    # the other.
    tr, tq, _, _ = CD.decode(TWIN, CD.tilecell(TWIN, 50, 6))
    tilenode(g, l) = CD.MemoBlockCursor(CD.BlockCursor(g, TWIN, CD.Bisected(), l,
        Int64(-1), tr, tr, tq, tq, 0, 0, 0, 0, false))
    lvl0 = tilenode(levelgrid(TWIN, 0), 0)
    lvl1 = tilenode(levelgrid(TWIN, 1), 1)
    @test CD._nodekey(lvl0.node) == CD._nodekey(lvl1.node)
    @test STI.node_extent(lvl0.node) !== STI.node_extent(lvl1.node)
    for _ in 1:4
        @test STI.node_extent(lvl0) === STI.node_extent(lvl0.node)
        @test STI.node_extent(lvl1) === STI.node_extent(lvl1.node)
    end
    # And turning back to the tile grid still answers for the tile grid.
    @test STI.node_extent(memoroot) === STI.node_extent(bareroot)

    # ---- (c) concurrent readers of one shared grid --------------------------
    # The table is task-local, so 8 tasks walking the same nodes share no slot.
    bare = [n.node for n in nodes]
    want = [STI.node_extent(n) for n in bare]
    tasks = map(1:8) do _
        Threads.@spawn begin
            bad = 0
            for (k, n) in enumerate(nodes)
                STI.node_extent(n) === want[k] || (bad += 1)
                k % 32 == 0 && yield()
            end
            bad
        end
    end
    @test all(==(0), fetch.(tasks))
    # Tasks that alternate between two lattices keep their own reset straight.
    mixed = map(1:8) do t
        Threads.@spawn begin
            bad = 0
            for _ in 1:8
                STI.node_extent(lvl0) === STI.node_extent(lvl0.node) || (bad += 1)
                yield()
                STI.node_extent(lvl1) === STI.node_extent(lvl1.node) || (bad += 1)
                yield()
            end
            bad
        end
    end
    @test all(==(0), fetch.(mixed))
end

# =========================================================================
# (k4) The leaf cells: the same pairs, in the same order, for nothing
# =========================================================================

# One warmed, type-stable pass over a leaf's cells. Kept out of the testset so
# the node's type is concrete at the call, which is what `@allocated` needs to
# mean anything: an inferred call to a heap-free build is the whole point.
function _leafcell_sum(node)
    s = 0
    for (i, _) in GO.SpatialTreeInterface.child_indices_extents(node)
        s += i
    end
    return s
end

# The `Ref` is allocated before the measurement and forces the call to happen.
function leafcell_bytes(node)
    warm = _leafcell_sum(node)
    total = Ref(warm)
    bytes = @allocated (total[] += _leafcell_sum(node))
    return (bytes, total[] - warm)
end

# `child_indices_extents` hands back a `LeafCells` instead of a fresh
# vector per call. It must behave as the vector did — indexable, iterable,
# collectable, same length and eltype — and yield the same pairs, bit for bit.
@testset "the leaf cells" begin
    STI = GO.SpatialTreeInterface
    Cap = US.SphericalCap{Float64}

    # The materializing implementation `LeafCells` replaced, verbatim.
    function materialized(c::CD.BlockCursor)
        entries = Tuple{Int,Cap}[]
        if c.inpixels
            for j in c.j0:c.j1, i in c.i0:c.i1
                leaf = CD.BlockCursor(c.grid, c.sys, c.strategy, c.level, c.origin,
                    c.r0, c.r0, c.q0, c.q0, j, j, i, i, true)
                push!(entries, (CD._index(c, c.r0, c.q0, j, i), STI.node_extent(leaf)))
            end
        else
            for r in c.r0:c.r1, q in c.q0:c.q1
                leaf = CD.BlockCursor(c.grid, c.sys, c.strategy, c.level, c.origin,
                    r, r, q, q, 0, 0, 0, 0, false)
                push!(entries, (CD._index(c, r, q, 0, 0), STI.node_extent(leaf)))
            end
        end
        return entries
    end

    function someleaves(root, cap)
        out = typeof(root)[]
        stack = [root]
        while !isempty(stack) && length(out) < cap
            n = pop!(stack)
            if STI.isleaf(n)
                push!(out, n)
            else
                for k in 1:STI.nchild(n)
                    push!(stack, STI.getchild(n, k))
                end
            end
        end
        return out
    end

    # Tile leaves, pixel leaves, pole rows, and both split strategies.
    for (label, grid) in (("twin tiles", levelgrid(TWIN, 0)),
                          ("twin pixels", levelgrid(TWIN, 1)),
                          ("GLO-90 tiles", levelgrid(GLO90, 0)),
                          ("GLO-90 pixels", levelgrid(GLO90, 1))),
        strategy in (CD.Blocked{3}(), CD.Bisected())

        root = CD.BlockCursor(grid; strategy)
        ls = someleaves(root, 400)
        @test (label, isempty(ls)) == (label, false)
        bad = 0
        for l in ls
            e = STI.child_indices_extents(l)
            want = materialized(l)
            # `===` on a tuple of `Int` and an immutable cap compares the bits.
            (e isa AbstractVector && eltype(e) == Tuple{Int,Cap} &&
             size(e) == (length(want),) && length(e) == length(want) &&
             all(e[k] === want[k] for k in eachindex(want)) &&
             collect(e) == want && all(x === want[k] for (k, x) in enumerate(e))) ||
                (bad += 1)
        end
        @test (label, string(typeof(strategy)), bad) ==
              (label, string(typeof(strategy)), 0)
    end

    # The whole value is `isbits`, which is why it never reaches the heap.
    @test isbitstype(CD.LeafCells)
    @test Base.IndexStyle(CD.LeafCells) === IndexLinear()

    # An interior node is still an error, as it was for the vector.
    globe = CD.BlockCursor(levelgrid(GLO90, 0))
    @test !STI.isleaf(globe)
    @test_throws ArgumentError STI.child_indices_extents(globe)

    # ---- it allocates nothing: construction and a full pass -----------------
    # This is the point of `LeafCells`. 71.6% of a production regrid's
    # allocation was the per-leaf vector it replaces.
    for grid in (levelgrid(GLO90, 0), levelgrid(GLO90, 1))
        leaf = someleaves(CD.BlockCursor(grid), 1)[1]
        bytes, sum1 = leafcell_bytes(leaf)
        @test sum1 == sum(i for (i, _) in materialized(leaf))
        @test bytes == 0
        # A memoized cursor forwards the entries, so it allocates nothing either.
        mbytes, sum2 = leafcell_bytes(CD.MemoBlockCursor(leaf))
        @test sum2 == sum1
        @test mbytes == 0
    end

    # ---- two leaves' entries live at once, which is how the search reads them
    # `dual_depth_first_search` binds both leaves' entries and loops over them
    # nested, so a shared per-task buffer would serve one leaf's cells from the
    # other's. A self-join is exactly that shape: both sides are `BlockCursor`.
    let root = CD.BlockCursor(levelgrid(TWIN, 0))
        ls = someleaves(root, 3)
        a, b = ls[1], ls[end]
        @test a != b
        ea, eb = STI.child_indices_extents(a), STI.child_indices_extents(b)
        wa, wb = materialized(a), materialized(b)
        # Holding both is safe: neither call disturbed the other's answers.
        @test all(ea[k] === wa[k] for k in eachindex(wa))
        @test all(eb[k] === wb[k] for k in eachindex(wb))

        pairs = Tuple{Int,Int}[]
        STI.dual_depth_first_search((_, _) -> true, a, b) do i1, i2
            push!(pairs, (i1, i2))
        end
        @test pairs == [(i, j) for (i, _) in wa for (j, _) in wb]
        # And the search's pruning still agrees with the materialized caps.
        overlaps(x, y) =
            US.spherical_distance(x.point, y.point) <= x.radius + y.radius
        near = Tuple{Int,Int}[]
        STI.dual_depth_first_search(overlaps, a, b) do i1, i2
            push!(near, (i1, i2))
        end
        @test near == [(i, j) for (i, ca) in wa for (j, cb) in wb if overlaps(ca, cb)]
    end
end

# =========================================================================
# (l) Contract: the conformance suites
# =========================================================================

# The hierarchical harness is quadratic in child count, so it runs on the scaled twin.
@testset "conformance (scaled twin)" begin
    test_grid_interface(levelgrid(TWIN, 0); label = "CopernicusDEM twin level 0")
    test_grid_interface(levelgrid(TWIN, 1); label = "CopernicusDEM twin level 1")
    test_hierarchical_system(TWIN; label = "CopernicusDEM twin")

    # Twin pole rings clear the harness's absolute degeneracy floor.
    N = CD.lat_intervals(TWIN)
    for (lat_s, j) in ((89, 0), (-90, N - 1))
        c = CD.pixelcell(TWIN, CD.tilecell(TWIN, lat_s, 0), j, 0)
        @test CT.boundary_problems(cell_boundary(TWIN, c)) == String[]
    end
end

# The harness's absolute `1e-12` degeneracy floor rejects valid level-1 pole rings.
# Use a reproducible sample without pole rows; testset (e) covers them separately.
@testset "conformance (GLO-30 and GLO-90)" begin
    for sys in (GLO90, GLO30)
        test_grid_interface(levelgrid(sys, 0); label = "$sys level 0")

        # Reproduce and validate the harness sample before using it.
        @test !any(in((89, -90)),
                   sampled_tile_lats(sys, 1, CONFORMANCE_SEED, GI_SAMPLES))

        # Matching RNG states confirm the harness consumed the reproduced draw.
        r = MersenneTwister(CONFORMANCE_SEED)
        ref = MersenneTwister(CONFORMANCE_SEED)
        CT.sample_indices(ref, ncells(levelgrid(sys, 1)), GI_SAMPLES)
        test_grid_interface(levelgrid(sys, 1); n_samples = GI_SAMPLES,
                            rng = r, label = "$sys level 1")
        @test r == ref
    end

    # Hierarchical runs are opt-in because sibling partitioning is quadratic:
    #
    #   unset, or "0"   neither. The default, and what CI runs.
    #   "90" / "1"      GLO-90 only.
    #   "30"            GLO-30 only; not known to have completed.
    #   "all"           both.
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

# The generic vertex-matching walk misses endpoint-free shared segments at band
# boundaries. Reached through the generic `one_ring` method, since dispatching
# `neighbors` on the level grid now lands on this system's own hook.
fallback_neighbors(g, c; connectivity = DGG.Vertex()) =
    invoke(DGG.one_ring,
           Tuple{DGG.AbstractGrid,DGG.AbstractCellIndex,DGG.Connectivity},
           g, c, connectivity)

# Warmed-call allocations detect accidental spatial-tree fallback.
_neighbor_bytes(g, c) = (DGG.neighbors(g, c, 1); @allocated DGG.neighbors(g, c, 1))

# -------------------------------------------------------------------------
# Exact rational adjacency oracle, independent of implementation helpers.
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

# Anchor the oracle to measured `cell_box` registration within one Float64 step.
@testset "the oracle is registered as cell_box is" begin
    N = CD.lat_intervals(TWIN)
    # Distance in representable Float64 steps.
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
    # Laterals share a meridian edge at both levels.
    out = Tuple{Int,Int}[(J, mod(K - 1, m)), (J, mod(K + 1, m))]
    if !edge && (J == 0 || J == rows - 1)
        # Every pole-row triangle shares the apex point.
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

# Check the exact set, uniqueness, and counter-clockwise ring order.
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

# Check symmetric adjacency; sample pole rings because this is quadratic.
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

    # Literal ids anchor the documented offset, stride, and winding.
    c = CD.pixelcell(TWIN, CD.tilecell(TWIN, 0, 0), 15, 15)
    @test c == DGG.LevelIndex(1, 21384465)
    @test DGG.neighbors(g1, c, 1) == DGG.LevelIndex.(1,
        [21384434, 21384464, 21384494, 21384495, 21384496, 21384466, 21384436, 21384435])
    @test DGG.neighbors(g1, c, 1; connectivity = DGG.Edge()) ==
          DGG.LevelIndex.(1, [21384435, 21384464, 21384495, 21384466])
    # Decode literals as the documented `NW, W, SW, S, SE, E, NE, N`.
    @test [CD.decode(TWIN, x)[3:4] for x in DGG.neighbors(g1, c, 1)] ==
          [(14, 14), (15, 14), (16, 14), (16, 15), (16, 16), (15, 16), (14, 16), (14, 15)]
    @test [CD.decode(TWIN, x)[3:4]
           for x in DGG.neighbors(g1, c, 1; connectivity = DGG.Edge())] ==
          [(14, 15), (15, 14), (16, 15), (15, 16)]

    # ORACLE PIN on the OUTER rings' start. The winding law is start-invariant
    # by construction, so nothing else here would notice ring 2 rotating.
    # Ring 2 begins on the SAME spoke ring 1 does — the `NW` direction.
    # Hand-checked at this cell: its ring-1 compass bearings run NW(320.96)
    # W(270.28) SW(219.94) S(180) SE(140.06) E(89.72) NE(39.04) N(360),
    # decreasing through one turn, which is counter-clockwise seen from outside
    # and the sense its four boundary corners wind in (219.72, 140.28, 39.27,
    # 320.73); ring 2's first entry sits at 301.97, the first bearing
    # counter-clockwise past 320.96.
    let d = DGG.cellindex(g0, 20000)
        @test d == DGG.LevelIndex(0, 19999)
        @test DGG.ring(g0, d, 2) == DGG.LevelIndex.(0,
            [19637, 19997, 20357, 20717, 20718, 20719, 20720, 20721,
                20361, 20001, 19641, 19281, 19280, 19279, 19278, 19277])
    end

    # Closed-form adjacency must not build a spatial tree.
    for x in (c,
              CD.pixelcell(TWIN, CD.tilecell(TWIN, 0, 0), 15, 0),        # west tile edge
              CD.pixelcell(TWIN, CD.tilecell(TWIN, 0, -180), 0, 0),      # antimeridian corner
              CD.pixelcell(TWIN, CD.tilecell(TWIN, 50, 0), N - 1, 5))    # a band boundary
        @test _neighbor_bytes(g1, x) < 1024
    end
    @test _neighbor_bytes(g0, CD.tilecell(TWIN, 50, 0)) < 1024

    # Boundary alignment periods divide tile widths and are invariant across `N`.
    reduced(sys) = [(CD.ncols(sys, r), CD.ncols(sys, r + 1)) .÷
                    gcd(CD.ncols(sys, r), CD.ncols(sys, r + 1))
                    for r in 0:(CD.NROWS - 2) if CD.ncols(sys, r) != CD.ncols(sys, r + 1)]
    @test reduced(TWIN) == reduced(GLO30) == reduced(GLO90)
    @test reduced(TWIN) == [(1, 2), (3, 5), (2, 3), (3, 4), (2, 3),
                            (3, 2), (4, 3), (3, 2), (5, 3), (2, 1)]

    # Exhaust both facing rows at every twin band boundary.
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

    # Exhaust both twin pole rows and their facing rows.
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

    # Probe tile seams and antimeridian corners and midpoints.
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

    # Exhaust all level-0 tiles; band offsets remove shared crossing corners.
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

    # Spot-check GLO-90 against the same oracle.
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

    # Rings are disjoint ordered blocks whose concatenation is `neighbors`.
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

    # The generic walk is exact away from band boundaries and otherwise a subset.
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
            walk = Set(fallback_neighbors(g1, x; connectivity = conn))
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

    # At N50, an endpoint-free shared segment makes the walk report five, not seven.
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
    # These counts pin the walk's wrong answers on purpose: an upstream fix that
    # counts shared segments turns them into the closed form's 7, 7, 7 and 8 —
    # red here after such a fix means delete these pins, not a regression.
    @test length(fallback_neighbors(g0, tile)) == 5
    @test length(DGG.neighbors(g0, tile)) == 7
    @test issubset(Set(fallback_neighbors(g0, tile)), Set(DGG.neighbors(g0, tile)))
    # Raising the vertex tolerance past the offset makes these rings match.
    @test DGG.Fallbacks._shared_vertices(ring_a, ring_b, 1e-3) == 2

    # Repeat at level 1 and at ±80, where some corners coincide.
    for (lat_s, j, i, walked, closed) in ((50, N - 1, 5, 5, 7), (49, 0, 5, 5, 7),
                                          (80, N - 1, 2, 7, 8))
        x = CD.pixelcell(TWIN, CD.tilecell(TWIN, lat_s, 0), j, i)
        @test length(fallback_neighbors(g1, x)) == walked
        @test length(DGG.neighbors(g1, x)) == closed
    end
end

@testset "maxneighbors is attained" begin
    # Both advertised neighbour bounds must be attained on every lattice.
    for sys in ALL_SYSTEMS
        N = CD.lat_intervals(sys)
        g = levelgrid(sys, 1)
        for (lat_s, j) in ((89, 0), (-90, N - 1))
            c = CD.pixelcell(sys, CD.tilecell(sys, lat_s, 0), j, 0)
            @test length(DGG.neighbors(g, c)) == DGG.maxneighbors(sys, DGG.Vertex())
            # Under `Edge()`, a pole cell has two laterals and one equatorward overlap.
            @test length(DGG.neighbors(g, c; connectivity = DGG.Edge())) == 3
        end
        for (lat_s, j) in ((85, N - 1), (-86, 0))       # the wide side of +-85
            t = CD.tilecell(sys, lat_s, 0)
            nc = CD.ncols_at(sys, lat_s)
            @test maximum(length(DGG.neighbors(g, CD.pixelcell(sys, t, j, i);
                                               connectivity = DGG.Edge()))
                          for i in 0:(nc - 1)) == DGG.maxneighbors(sys, DGG.Edge())
        end
        @test DGG.maxneighbors(sys, DGG.Vertex()) == 36 * N + 2
        @test DGG.maxneighbors(sys, DGG.Edge()) == 6
    end

    N = CD.lat_intervals(TWIN)
    g1 = levelgrid(TWIN, 1)
    g0 = levelgrid(TWIN, 0)

    # Sweep column-dependent neighbour counts at every band boundary.
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

    # Off-pole twin cells stay below nine `Vertex()` and seven `Edge()` neighbours.
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

    # At level 0, band-boundary tiles have seven neighbours because no corner crosses.
    counts0 = [(length(DGG.neighbors(g0, DGG.LevelIndex(0, CD.tileordinal(r, 180)))),
                length(DGG.neighbors(g0, DGG.LevelIndex(0, CD.tileordinal(r, 180));
                                     connectivity = DGG.Edge())))
               for r in 0:(CD.NROWS - 1)]
    @test sort(unique(first.(counts0))) == [7, 8, 362]
    @test sort(unique(last.(counts0))) == [3, 4, 5]
    @test maximum(first.(counts0)) < DGG.maxneighbors(TWIN, DGG.Vertex())
    @test maximum(last.(counts0)) < DGG.maxneighbors(TWIN, DGG.Edge())
end

end # @testset "CopernicusDEM system"

end # module CopernicusDEMSystemTests
