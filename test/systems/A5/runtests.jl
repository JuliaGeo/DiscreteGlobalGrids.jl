# ---------------------------------------------------------------------------
# T10 — the A5 system suite.
#
# Three kinds of test, deliberately kept apart:
#
#   * ORACLE tests, which check the wiring against `A5Native` — the ported
#     upstream arithmetic — and against the sealed vectors the pre-redesign
#     `test/A5/` suite pinned (the Tokyo index, its centre, the antimeridian
#     ring, the res-30 children). `A5Native` is A5's definition here, so
#     anywhere this package could have drifted is checked against a direct call
#     rather than against a value someone typed in.
#   * CONFORMANCE tests, the two property suites from
#     `DiscreteGlobalGridsConformanceTesting`, which check the interface
#     contracts every system owes.
#   * SELECTION-MODE tests. A5 is the first real system with
#     `has_sorted_subtrees == false`, so `treeify`, `query` and
#     `MultiOrderCoverage` take the generic substrate's fallback paths here —
#     the materialised-selection cursor and the `(level, position)` set order —
#     and this suite is where those paths meet a real grid rather than a mock.
#
# Two measurement batteries print their numbers rather than only asserting
# them, because both are claims about A5's geometry that the docstrings quote:
# NEIGHBOR-VALIDATION for the `max_neighbors` bounds and the Moore/von Neumann
# split, CAP-VALIDATION for the raised `cap_inflation`.
# ---------------------------------------------------------------------------

module A5TestSuite

using Test
using Printf
using Random
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGridsConformanceTesting
using SmallCollections: SmallVector
import GeometryOps as GO
import GeoInterface as GI
import Extents
import GeometryOps: SpatialTreeInterface as STI
import ConservativeRegridding: Trees

const A5 = DGG.A5
const A5N = A5.A5Native
const S = A5.A5System()
const SD = GO.UnitSpherical.spherical_distance

# Complete levels, built by expanding the hierarchy one level at a time — never
# by the dense-ordinal arithmetic the ordinal testset checks against them.
const ROOTS = DGG.rootcells(S)
const RES1 = reduce(vcat, [collect(DGG.children(S, r)) for r in ROOTS])
const RES2 = reduce(vcat, [collect(DGG.children(S, c)) for c in RES1])
const RES3 = reduce(vcat, [collect(DGG.children(S, c)) for c in RES2])
const RES4 = reduce(vcat, [collect(DGG.children(S, c)) for c in RES3])

const CONNECTIVITIES = (DGG.Vertex(), DGG.Edge())
conn_name(::DGG.Vertex) = "Vertex"
conn_name(::DGG.Edge) = "Edge"

"A deterministic evenly spaced sweep of `n` positions of level `l`."
function ordinal_sample(l::Int, n::Int)
    grid = DGG.levelgrid(S, l)
    total = DGG.ncells(grid)
    step = max(1, (total - 1) ÷ max(1, n - 1))
    return [DGG.cellindex(grid, i) for i in 1:step:total]
end

# The lon/lat probe grid for `cellat`, avoiding the poles' coordinate
# degeneracy but keeping high latitudes.
const PROBE_LONLAT = [(lon, lat) for lon in -180.0:24.0:170.0, lat in -78.0:12.0:78.0]

# ---------------------------------------------------------------------------
# Rotational-order machinery: is a sequence of cells counter-clockwise about a
# centre, seen from OUTSIDE the sphere?
#
# Build a right-handed tangent basis at the centre (`east x north` is the
# outward normal), take each cell's azimuth in it, and sum the signed
# increments. A counter-clockwise cycle winds exactly +1 turn; a clockwise one
# -1; a sequence that is merely deterministic and not rotational winds 0.
#
# This is written independently of `neighbors.jl`'s own frame — a different
# seed axis, chosen as the one the centre leans on least rather than as a
# neighbour's direction — so that it tests the order rather than restating it.
# The seed matters here: two of A5's twelve res-0 cells are centred exactly on
# a pole, where a lon/lat east/north frame does not exist at all.
# ---------------------------------------------------------------------------
dot3(a, b) = a[1] * b[1] + a[2] * b[2] + a[3] * b[3]
cross3(a, b) = (a[2] * b[3] - a[3] * b[2],
                a[3] * b[1] - a[1] * b[3],
                a[1] * b[2] - a[2] * b[1])

function tangent_frame(p)
    ax = abs(p[1]) <= abs(p[2]) ?
         (abs(p[1]) <= abs(p[3]) ? (1.0, 0.0, 0.0) : (0.0, 0.0, 1.0)) :
         (abs(p[2]) <= abs(p[3]) ? (0.0, 1.0, 0.0) : (0.0, 0.0, 1.0))
    s = dot3(ax, p)
    e = (ax[1] - s * p[1], ax[2] - s * p[2], ax[3] - s * p[3])
    n = sqrt(dot3(e, e))
    east = (e[1] / n, e[2] / n, e[3] / n)
    return east, cross3(p, east)
end

"Signed winding of `cells` about `c`, in turns: +1 is counter-clockwise."
function winding_turns(grid, c, cells)
    p = DGG.cell_centroid(grid, c)
    east, north = tangent_frame(p)
    phis = Float64[]
    for m in cells
        q = DGG.cell_centroid(grid, m)
        r = dot3(q, p)
        d = (q[1] - r * p[1], q[2] - r * p[2], q[3] - r * p[3])
        sqrt(dot3(d, d)) <= 1e-9 && continue
        push!(phis, atan(dot3(d, north), dot3(d, east)))
    end
    length(phis) < 3 && return NaN
    total = 0.0
    for i in eachindex(phis)
        j = i == lastindex(phis) ? firstindex(phis) : i + 1
        d = phis[j] - phis[i]
        while d > pi
            d -= 2pi
        end
        while d <= -pi
            d += 2pi
        end
        total += d
    end
    return total / 2pi
end

# The five (or three) CORNERS of a cell, undensified — the vertices two
# adjacent cells actually share, as opposed to the interpolated points the
# published ring carries between them.
corners(c) = [DGG.Fallbacks.unit_point(p[1], p[2])
              for p in A5N.cell_boundary(DGG.rawid(c); closed_ring=false, segments=1)]

chord(p, q) = sqrt((p[1] - q[1])^2 + (p[2] - q[2])^2 + (p[3] - q[3])^2)

function shortest_edge(ring)
    n = length(ring)
    m = Inf
    for i in 1:n
        d = chord(ring[i], ring[i == n ? 1 : i+1])
        d > 0 && (m = min(m, d))
    end
    return m
end

"How many corners of `a` coincide with a corner of `b`."
shared_corners(a, b, tol) = count(p -> any(q -> chord(p, q) <= tol, b), a)

ring_points(polygon) = collect(GI.getpoint(GI.getexterior(polygon)))

@testset "A5" begin

    # =======================================================================
    @testset "native layer (ported oracles)" begin
        roots = A5N.res0_cells()
        @test length(roots) == 12
        @test issorted(roots)
        @test A5N.num_cells(-1) == 1
        @test A5N.num_cells(0) == 12
        @test A5N.num_cells(1) == 60
        @test A5N.num_cells(2) == 240
        # Res 30 exists only for the 42 quintants that fit 64 bits, so there is
        # no complete grid there and the native count refuses rather than guessing.
        @test A5N.MAX_RESOLUTION == 30
        @test A5N.MAX_GRID_RESOLUTION == 29
        @test_throws ArgumentError A5N.num_cells(30)

        # The canonical A5 documentation example, sealed.
        tokyo = A5N.lonlat_to_cell(139.7623824402441, 35.677369792795794, 4)
        @test tokyo == parse(UInt64, "8708000000000000"; base=16)
        @test A5N.get_resolution(tokyo) == 4
        lon, lat = A5N.cell_to_lonlat(tokyo)
        @test isapprox(lon, 141.11036610726177; atol=1e-10)
        @test isapprox(lat, 34.33639263740687; atol=1e-10)

        # A cell that straddles the antimeridian still reports one contiguous
        # longitude run: the native ring re-centres its longitudes.
        antimeridian = parse(UInt64, "eb60000000000000"; base=16)
        lons = first.(A5N.cell_boundary(antimeridian; segments=1))
        @test maximum(lons) - minimum(lons) < 180

        # Res 30 is reachable for a supported quintant and refused for an
        # unsupported one — the encoding edge `A5Cell` documents.
        supported = A5N.serialize(A5.NativeCell(A5N.ORIGINS[1], 0, 0, 28))
        @test length(A5N.cell_to_children(supported, 30)) == 16
        @test all(==(30) ∘ A5N.get_resolution, A5N.cell_to_children(supported, 30))
        unsupported_origin = A5N.ORIGINS[10]
        unsupported = A5N.serialize(
            A5.NativeCell(unsupported_origin, unsupported_origin.first_quintant, 0, 28))
        @test_throws ArgumentError A5N.cell_to_children(unsupported, 30)

        for resolution in 0:8
            sampled = A5N.lonlat_to_cell(12.5, -34.25, resolution)
            @test A5N.get_resolution(sampled) == resolution
            @test A5N.cell_to_parent(sampled, resolution) == sampled
            resolution > 0 &&
                @test sampled in A5N.cell_to_children(A5N.cell_to_parent(sampled), resolution)
        end
    end

    # =======================================================================
    @testset "A5Cell encoding" begin
        for l in (0, 1, 2, 5, 12, 29), c in ordinal_sample(l, 6)
            @test DGG.level(c) == A5N.get_resolution(DGG.rawid(c))
            @test DGG.level(c) == l
            @test DGG.rawid(c) isa UInt64
            @test A5.A5Cell(DGG.rawid(c)) == c
            @test isvalid(c)
        end

        a = A5.A5Cell(0x0080000000000000)      # res 2, quintant 0, S 0
        b = A5.A5Cell(0x0180000000000000)      # res 2, quintant 0, S 1
        @test a < b
        @test isless(a, b)
        @test a == A5.A5Cell(0x0080000000000000)
        @test hash(a) == hash(A5.A5Cell(0x0080000000000000))
        @test occursin("res 2", sprint(show, a))

        # `level` is a READ of the encoding, so it answers for the two ids that
        # name no cell of any grid. `isvalid` is the predicate that says so.
        world = A5.A5Cell(A5N.WORLD_CELL)
        @test DGG.level(world) == -1
        @test !isvalid(world)
        res30 = A5.A5Cell(A5N.serialize(A5.NativeCell(A5N.ORIGINS[1], 0, UInt64(0), 30)))
        @test DGG.level(res30) == 30
        @test !isvalid(res30)
        # ... and for the three ways an id can look plausible and name nothing.
        @test !isvalid(A5.A5Cell(UInt64(60) << 58 | UInt64(1) << 56))   # quintant 60
        @test !isvalid(A5.A5Cell(UInt64(12) << 58 | UInt64(1) << 57))   # face 12 at res 0
        # Junk in the padding below the resolution marker. Most such bits move
        # the marker and so change the reported level, but not all: the marker
        # walk only ever examines odd bit positions, so an even bit below it
        # (here bit 6) is invisible to `level`, is DISCARDED by `deserialize`,
        # and would otherwise decode to — and be given the position of — a
        # perfectly good cell that it is not. That is what the round-trip half
        # of `isvalid` is for.
        junk = A5.A5Cell(DGG.rawid(RES3[1]) | (UInt64(1) << 6))
        @test DGG.level(junk) == 3
        @test A5N.deserialize(DGG.rawid(junk)) == A5N.deserialize(DGG.rawid(RES3[1]))
        @test junk != RES3[1]
        @test !isvalid(junk)
        @test isvalid(RES3[1])
    end

    # =======================================================================
    @testset "dense order: cellindex/cellposition" begin
        # The hierarchy-built complete levels, against the ordinal arithmetic.
        for (l, cells) in ((0, ROOTS), (1, RES1), (2, RES2), (3, RES3), (4, RES4))
            grid = DGG.levelgrid(S, l)
            @test DGG.ncells(grid) == length(cells)
            @test issorted(cells)
            @test allunique(cells)
            @test [DGG.cellindex(grid, i) for i in 1:length(cells)] == cells
            @test [DGG.cellposition(grid, c) for c in cells] == collect(1:length(cells))
            @test_throws BoundsError DGG.cellindex(grid, 0)
            @test_throws BoundsError DGG.cellindex(grid, length(cells) + 1)
        end
        @test DGG.ncells(DGG.levelgrid(S, 29)) == 60 * Int64(4)^28

        # Deep levels, where no complete level can be held: a sweep of positions
        # must come back ascending, at the right resolution, and round-trip.
        for l in (5, 9, 15, 22, 29)
            grid = DGG.levelgrid(S, l)
            total = DGG.ncells(grid)
            step = (total - 1) ÷ 17
            positions = sort!(unique!(vcat([1, 2, total - 1, total],
                [1 + k * step for k in 1:16])))
            cells = [DGG.cellindex(grid, i) for i in positions]
            @test all(c -> DGG.level(c) == l, cells)
            @test issorted(cells)
            @test [DGG.cellposition(grid, c) for c in cells] == positions
        end

        # Everything that is not a cell of the grid answers `nothing`, never an
        # error — a different level, the world cell, a res-30 id, a malformed
        # id, and an id of another system entirely.
        g2 = DGG.levelgrid(S, 2)
        @test DGG.cellposition(g2, RES3[1]) === nothing
        @test DGG.cellposition(g2, ROOTS[1]) === nothing
        @test DGG.cellposition(g2, A5.A5Cell(A5N.WORLD_CELL)) === nothing
        @test DGG.cellposition(g2, A5.A5Cell(UInt64(60) << 58 | UInt64(1) << 55)) === nothing
        @test DGG.cellposition(g2, A5.A5Cell(DGG.rawid(RES2[1]) | UInt64(1))) === nothing
        @test DGG.cellposition(g2, DGG.LevelIndex(2, 1)) === nothing
    end

    # =======================================================================
    @testset "hierarchy: 12 -> 5 -> 4" begin
        @test length(ROOTS) == 12
        @test all(c -> DGG.level(c) == 0, ROOTS)
        @test length(RES1) == 60
        @test length(RES2) == 240
        @test length(RES3) == 960

        for c in vcat(ROOTS, RES1, RES2[1:17:end], ordinal_sample(9, 8), ordinal_sample(28, 4))
            l = DGG.level(c)
            kids = DGG.children(S, c)
            @test length(kids) == (l == 0 ? 5 : 4)
            @test issorted(kids)
            @test allunique(kids)
            @test typeof(kids) === SmallVector{5,A5.A5Cell}
            for k in kids
                @test DGG.level(k) == l + 1
                @test parent(S, k) == c
                @test DGG.rawid(parent(S, k)) == A5N.cell_to_parent(DGG.rawid(k), l)
            end
            if l > 0
                p = parent(S, c)
                @test c in DGG.children(S, p)
                for j in 0:(l-1)
                    @test DGG.rawid(DGG.ancestor(S, c, j)) == A5N.cell_to_parent(DGG.rawid(c), j)
                end
            end
            @test DGG.ancestor(S, c, l) == c
            @test collect(DGG.descendants(S, c, l)) == [c]
            for d in 1:2
                l + d > 29 && continue
                # The oracle is the level-by-level expansion through `children`,
                # not the native multi-level call the implementation uses.
                expected = [c]
                for _ in 1:d
                    expected = reduce(vcat, [collect(DGG.children(S, x)) for x in expected])
                end
                @test DGG.descendants(S, c, l + d) == sort(expected)
            end
        end

        # A complete level really is its parents' subtrees, concatenated and sorted.
        @test sort(reduce(vcat, [DGG.descendants(S, r, 2) for r in ROOTS])) == RES2
        @test sort(reduce(vcat, [DGG.descendants(S, c, 3) for c in RES1])) == RES3
        @test all(c -> DGG.ancestor(S, c, 0) in ROOTS, RES3)

        # Error contract at the two ends of the hierarchy.
        @test_throws ArgumentError parent(S, ROOTS[1])
        @test_throws ArgumentError DGG.children(S, DGG.cellindex(DGG.levelgrid(S, 29), 1))
        @test_throws ArgumentError DGG.ancestor(S, RES2[1], 3)
        @test_throws ArgumentError DGG.descendants(S, RES2[1], 1)
        @test_throws ArgumentError DGG.descendants(S, RES2[1], 30)
        @test_throws ArgumentError DGG.levelgrid(S, -1)
        @test_throws ArgumentError DGG.levelgrid(S, 30)
        # Garbage in is an error, not a confidently enumerated subtree of cells
        # that do not exist.
        @test_throws ArgumentError DGG.children(S, A5.A5Cell(A5N.WORLD_CELL))
        @test_throws ArgumentError DGG.descendants(S, A5.A5Cell(UInt64(60) << 58 | UInt64(1) << 55), 3)
    end

    # =======================================================================
    @testset "boundary, centroid and area" begin
        # Ring shape: a root is a pentagon, a quintant a triangle, everything
        # below a pentagon again — each subdivided `2^(6-level)` times per edge
        # down to level 6, once above it.
        for (l, c, corner_count) in ((0, ROOTS[1], 5), (1, RES1[1], 3), (2, RES2[1], 5),
                (3, RES3[1], 5), (7, ordinal_sample(7, 1)[1], 5))
            grid = DGG.levelgrid(S, l)
            ring = DGG.cell_boundary(grid, c)
            @test ring isa Vector{GO.UnitSphericalPoint{Float64}}
            @test length(ring) == corner_count * (l <= 6 ? 2^(6 - l) : 1)
            @test length(corners(c)) == corner_count
            @test all(p -> abs(sqrt(sum(abs2, p)) - 1) < 1e-12, ring)
            @test ring[1] != ring[end]                       # implicitly closed
            @test DGG.cell_area(grid, c) > 0

            centroid = DGG.cell_centroid(grid, c)
            lon, lat = A5N.cell_to_lonlat(DGG.rawid(c))
            @test abs(sqrt(sum(abs2, centroid)) - 1) < 1e-12
            @test centroid[3] ≈ sin(deg2rad(lat)) atol = 1e-14
            @test atan(centroid[2], centroid[1]) ≈ deg2rad(lon) atol = 1e-12

            polygon = DGG.cell_polygon(grid, c)
            @test GI.trait(polygon) isa GI.PolygonTrait
            @test ring_points(polygon) == vcat(ring, ring[1:1])
        end

        # Counter-clockwise seen from OUTSIDE, over whole levels rather than a
        # sample: the winding the base interface requires, and the one the
        # native ring gets by reversing itself on the way out.
        for (l, cells) in ((0, ROOTS), (1, RES1), (2, RES2))
            grid = DGG.levelgrid(S, l)
            @test all(c -> DGG.cell_area(grid, c) > 0, cells)
        end

        # A whole level partitions the sphere. A5 is equal-area on the
        # ELLIPSOID and reports geodetic coordinates, so on the unit sphere the
        # individual areas carry the authalic -> geodetic latitude conversion
        # and spread about 1% peak to peak; the sum is still 4pi.
        for (l, cells) in ((0, ROOTS), (1, RES1))
            grid = DGG.levelgrid(S, l)
            areas = [DGG.cell_area(grid, c) for c in cells]
            @test sum(areas) ≈ 4pi rtol = 1e-7
            @test maximum(areas) - minimum(areas) < 0.015 * maximum(areas)
        end
        # From level 2 the subdivided rings' spherical areas drift a fraction of
        # a per mille in the sum; that residual is the native geometry's.
        res2_areas = [DGG.cell_area(DGG.levelgrid(S, 2), c) for c in RES2]
        @test sum(res2_areas) ≈ 4pi rtol = 1e-3
        @test maximum(res2_areas) - minimum(res2_areas) < 2e-2 * maximum(res2_areas)

        # The two res-0 cells centred exactly on a pole, which is why nothing in
        # this system may build a tangent frame from local east/north.
        polar = [c for c in ROOTS if abs(abs(DGG.cell_centroid(DGG.levelgrid(S, 0), c)[3]) - 1) < 1e-12]
        @test length(polar) == 2
        ext = DGG.cell_extent(DGG.levelgrid(S, 0), polar[1])
        @test ext.X == (-180.0, 180.0)
        @test 90.0 in ext.Y || -90.0 in ext.Y
    end

    # =======================================================================
    @testset "cellat vs lonlat_to_cell" begin
        for l in (0, 1, 4, 9, 20)
            grid = DGG.levelgrid(S, l)
            for (lon, lat) in PROBE_LONLAT
                expected = A5.A5Cell(A5N.lonlat_to_cell(lon, lat, l))
                @test DGG.cellat(grid, lon, lat) == expected
                # The unit-sphere primitive agrees with the lon/lat wrapper.
                @test DGG.cellat(grid, DGG.Fallbacks.unit_point(lon, lat)) == expected
            end
        end
        # Round trip, and the interior probe that a nearest-centroid lookup
        # would fail: half way out to every boundary vertex is still this cell.
        for l in (0, 1, 3, 6, 12, 29)
            grid = DGG.levelgrid(S, l)
            for c in ordinal_sample(l, 12)
                centroid = DGG.cell_centroid(grid, c)
                @test DGG.cellat(grid, centroid) == c
                @test all(v -> DGG.cellat(grid, GO.UnitSpherical.slerp(centroid, v, 0.5)) == c,
                    DGG.cell_boundary(grid, c))
            end
        end
    end

    # =======================================================================
    # NEIGHBOR-VALIDATION. Both halves of the `max_neighbors` docstring are
    # claims about A5's geometry, and both are measured here rather than read
    # off the cell shape:
    #
    #   * the DEGREES, which are constant per regime rather than merely
    #     bounded — that is the stronger statement, and the one that makes the
    #     bound safe to extrapolate past the sampled levels;
    #   * the MOORE/VON NEUMANN SPLIT, against the ring corners themselves.
    #     "A5 cells are pentagons, so Vertex() == Edge()" is exactly the
    #     inherited assumption this battery refutes.
    # =======================================================================
    @testset "NEIGHBOR-VALIDATION: the degree bounds are measured" begin
        @test DGG.max_neighbors(S, DGG.Vertex()) == 11
        @test DGG.max_neighbors(S, DGG.Edge()) == 5
        @test DGG.max_neighbors(S) == 11

        println("\n  A5 NEIGHBOR-VALIDATION — degrees per level and connectivity")
        groups = (("res 0 (all 12)", 0, ROOTS), ("res 1 (all 60)", 1, RES1),
            ("res 2 (all 240)", 2, RES2), ("res 3 (all 960)", 3, RES3),
            ("res 4 (all 3840)", 4, RES4),
            ("res 9 (200 sample)", 9, ordinal_sample(9, 200)),
            ("res 15 (200 sample)", 15, ordinal_sample(15, 200)),
            ("res 22 (200 sample)", 22, ordinal_sample(22, 200)),
            ("res 29 (200 sample)", 29, ordinal_sample(29, 200)))
        worst = Dict(DGG.Vertex() => 0, DGG.Edge() => 0)
        for (label, l, cells) in groups
            grid = DGG.levelgrid(S, l)
            observed = Dict(conn =>
                sort(unique(length(DGG.neighbors(grid, c, 1; connectivity=conn)) for c in cells))
                             for conn in CONNECTIVITIES)
            @printf("  %-20s (%5d cells) Vertex: %-10s Edge: %s\n", label, length(cells),
                join(observed[DGG.Vertex()], ", "), join(observed[DGG.Edge()], ", "))
            for conn in CONNECTIVITIES
                @test maximum(observed[conn]) <= DGG.max_neighbors(S, conn)
                worst[conn] = max(worst[conn], maximum(observed[conn]))
            end
            # Constant per regime, not merely bounded: 5 at res 0 (a
            # dodecahedron face has 5 face-adjacent faces and shares a vertex
            # with no other), 3 edge / 11 vertex at res 1 (a quintant is a
            # triangle whose corners are the face centre and two dodecahedron
            # vertices), 5 edge / 6-8 vertex from res 2 down.
            @test observed[DGG.Edge()] == [l == 1 ? 3 : 5]
            @test observed[DGG.Vertex()] ⊆ (l == 0 ? [5] : l == 1 ? [11] : [6, 7, 8])
        end
        @printf("  worst observed — Vertex %d (bound %d), Edge %d (bound %d)\n\n",
            worst[DGG.Vertex()], DGG.max_neighbors(S, DGG.Vertex()),
            worst[DGG.Edge()], DGG.max_neighbors(S, DGG.Edge()))
        # Both bounds are attained, so neither has slack to trim.
        @test worst[DGG.Vertex()] == DGG.max_neighbors(S, DGG.Vertex())
        @test worst[DGG.Edge()] == DGG.max_neighbors(S, DGG.Edge())
    end

    @testset "Vertex() really is Moore and Edge() really is von Neumann" begin
        # THE ASSUMPTION UNDER TEST: `Vertex()` docs say hexagonal and
        # pentagonal grids coincide. That holds for the icosahedral hex systems,
        # where three cells meet at every vertex. A5 tiles in the manner of a
        # Cairo pentagon tiling — four cells meet at some corners — so it does
        # not, and the two connectivities name different sets.
        #
        # Checked against the ring CORNERS, not the densified boundary: an edge
        # neighbour shares two corners, a corner-only neighbour exactly one.
        # `probe` is the cells examined; `level_cells` is the WHOLE level, which
        # the completeness half needs — "nothing else touches this cell" is only
        # a statement if the search covers every cell there is.
        for (l, probe, level_cells) in ((1, RES1, RES1), (2, RES2, RES2),
                (3, RES3[1:29:end], RES3))
            grid = DGG.levelgrid(S, l)
            rings = Dict(x => corners(x) for x in level_cells)
            for c in probe
                mine = rings[c]
                tol = 0.05 * shortest_edge(mine)
                vs = Set(DGG.neighbors(grid, c, 1; connectivity=DGG.Vertex()))
                es = Set(DGG.neighbors(grid, c, 1; connectivity=DGG.Edge()))
                @test issubset(es, vs)
                @test es != vs                       # they genuinely differ
                @test all(x -> shared_corners(mine, rings[x], tol) == 2, es)
                @test all(x -> shared_corners(mine, rings[x], tol) == 1, setdiff(vs, es))
                # ... and nothing else in the level touches this cell at all, so
                # `Vertex()` is the complete Moore neighbourhood rather than
                # merely a bigger set than `Edge()`.
                touching = Set(x for x in level_cells
                               if x != c && shared_corners(mine, rings[x], tol) >= 1)
                @test touching == vs
            end
        end
        # Res 0 is the one level where the two coincide, and it is a fact about
        # the dodecahedron rather than about the wiring.
        grid = DGG.levelgrid(S, 0)
        @test all(ROOTS) do c
            collect(DGG.neighbors(grid, c, 1; connectivity=DGG.Vertex())) ==
            collect(DGG.neighbors(grid, c, 1; connectivity=DGG.Edge()))
        end
    end

    # =======================================================================
    @testset "neighbors: membership, symmetry and containers" begin
        for (l, cells) in ((0, ROOTS), (1, RES1), (2, RES2), (3, RES3), (4, RES4))
            grid = DGG.levelgrid(S, l)
            stored = Set(cells)
            for conn in CONNECTIVITIES
                table = Dict(c => collect(DGG.neighbors(grid, c, 1; connectivity=conn))
                             for c in cells)
                @test all(c -> allunique(table[c]), cells)
                @test all(c -> !(c in table[c]), cells)
                @test all(c -> all(n -> DGG.level(n) == l, table[c]), cells)
                @test all(c -> all(n -> n in stored, table[c]), cells)
                # With every cell of the level in hand, "b in N(a) implies
                # a in N(b)" over all a IS symmetry.
                @test all(c -> all(n -> c in table[n], table[c]), cells)
                @test all(c -> Set(table[c]) == Set(DGG.neighbors(grid, c, 1; connectivity=conn)),
                    cells)
            end
        end

        # Deep levels, where no complete level can be held: the same checks on a
        # patch — sampled cells closed under one neighbour step — so both
        # directions are testable for every pair inside the patch.
        for l in (9, 15, 22, 29), conn in CONNECTIVITIES
            grid = DGG.levelgrid(S, l)
            seeds = ordinal_sample(l, 40)
            patch = Set(seeds)
            for c in seeds, n in DGG.neighbors(grid, c, 1; connectivity=conn)
                push!(patch, n)
            end
            cells = collect(patch)
            table = Dict(c => collect(DGG.neighbors(grid, c, 1; connectivity=conn)) for c in cells)
            @test all(c -> allunique(table[c]) && !(c in table[c]), cells)
            @test all(c -> all(n -> DGG.level(n) == l, table[c]), cells)
            # Validity without a level to enumerate: the position round trip,
            # which only closes for an id that really is a cell of `l`.
            @test all(c -> all(n ->
                    DGG.cellindex(grid, DGG.cellposition(grid, n)) == n, table[c]), cells)
            @test all(c -> all(n -> !(n in patch) || (c in table[n]), table[c]), cells)
        end

        # k = 0 is empty and has the same concrete type as k = 1, which is what
        # keeps a neighbour sweep type-stable across the `k` boundary.
        grid = DGG.levelgrid(S, 5)
        c = DGG.cellindex(grid, 1000)
        @test isempty(DGG.neighbors(grid, c, 0))
        @test typeof(DGG.neighbors(grid, c, 0)) === SmallVector{11,A5.A5Cell}
        @test typeof(DGG.neighbors(grid, c, 1)) === SmallVector{11,A5.A5Cell}
        @test typeof(DGG.neighbors(grid, c)) === SmallVector{11,A5.A5Cell}
        @test typeof(DGG.neighbors(grid, c, 1; connectivity=DGG.Edge())) ===
              SmallVector{11,A5.A5Cell}
        @test typeof(DGG.neighbors(grid, c, 2)) === Vector{A5.A5Cell}
        @test typeof(DGG.ring(grid, c, 1)) === Vector{A5.A5Cell}
        @test typeof(DGG.children(S, c)) === SmallVector{5,A5.A5Cell}
        DGG.children(S, c)
        @test @allocated(DGG.children(S, c)) == 0

        # An id that is not a cell of this grid's resolution is an error, not a
        # confident answer about some other cell.
        @test_throws ArgumentError DGG.neighbors(grid, ROOTS[1])
        @test_throws ArgumentError DGG.neighbors(grid, A5.A5Cell(A5N.WORLD_CELL))
        @test_throws ArgumentError DGG.neighbors(grid, c, -1)
        @test_throws ArgumentError DGG.ring(grid, c, -1)
    end

    # =======================================================================
    @testset "rotational order: winding, shells and the disc" begin
        for l in (0, 1, 2, 3, 5, 9)
            grid = DGG.levelgrid(S, l)
            cells = l <= 1 ? [DGG.cellindex(grid, i) for i in 1:DGG.ncells(grid)] :
                    ordinal_sample(l, 40)
            for conn in CONNECTIVITIES, c in cells
                @test DGG.ring(grid, c, 0; connectivity=conn) == [c]
                shells = Vector{A5.A5Cell}[]
                seen = Set{A5.A5Cell}()
                for k in 1:3
                    shell = DGG.ring(grid, c, k; connectivity=conn)
                    @test allunique(shell)
                    @test isempty(intersect(Set(shell), seen))
                    union!(seen, shell)
                    push!(shells, collect(shell))
                    # THE ORDER IS ROTATIONAL: counter-clockwise seen from
                    # outside, exactly one turn. This is the property the
                    # native `sort!(collect(::Set))` destroyed.
                    length(shell) >= 3 && @test winding_turns(grid, c, shell) ≈ 1.0 atol = 0.02
                    # The disc is the shells CONCATENATED OUTWARD, element for
                    # element — not merely their union — so the ring is the
                    # tail block of the disc.
                    disc = collect(DGG.neighbors(grid, c, k; connectivity=conn))
                    @test disc == reduce(vcat, shells)
                    @test Set(disc) == seen
                    @test disc[end-length(shell)+1:end] == collect(shell)
                end
            end
        end
    end

    # =======================================================================
    @testset "ring 1 starts where the docstring says" begin
        # ORACLE PIN on the documented START of the rotational order.
        #
        # The conformance harness's winding law is start-invariant by design —
        # it asserts a ring is one CCW cycle, and every rotation of a cycle is
        # still one cycle. So without this literal, nothing would catch A5's
        # ring 1 quietly beginning at a different neighbour. That is worth
        # catching: the value of a rotational order is that slot `j` names a
        # fixed direction, so a stencil that bakes "slot 1" into a weight
        # vector is silently rotated if the start moves.
        #
        # Produced by the implementation and checked against the documented
        # rule — the smallest-id neighbour first, then counter-clockwise — not
        # transcribed from the native walk, which sorts by id throughout and is
        # exactly what this module replaces.
        grid = DGG.levelgrid(S, 2)
        c = DGG.cellindex(grid, 100)
        @test string(DGG.rawid(c); base=16) == "6380000000000000"
        @test [string(DGG.rawid(x); base=16) for x in DGG.neighbors(grid, c, 1)] ==
              ["b80000000000000", "c80000000000000", "6180000000000000",
               "5380000000000000", "5180000000000000", "6a80000000000000"]
        @test [string(DGG.rawid(x); base=16)
               for x in DGG.neighbors(grid, c, 1; connectivity=DGG.Edge())] ==
              ["b80000000000000", "c80000000000000", "6180000000000000",
               "5180000000000000", "6a80000000000000"]

        # The res-1 regime, where the two connectivities are furthest apart.
        g1 = DGG.levelgrid(S, 1)
        c1 = DGG.cellindex(g1, 7)
        @test string(DGG.rawid(c1); base=16) == "1900000000000000"
        @test [string(DGG.rawid(x); base=16) for x in DGG.neighbors(g1, c1, 1)] ==
              ["100000000000000", "ed00000000000000", "dd00000000000000",
               "e100000000000000", "d900000000000000", "d500000000000000",
               "1d00000000000000", "2100000000000000", "2500000000000000",
               "1500000000000000", "1100000000000000"]
        @test [string(DGG.rawid(x); base=16)
               for x in DGG.neighbors(g1, c1, 1; connectivity=DGG.Edge())] ==
              ["1500000000000000", "dd00000000000000", "1d00000000000000"]

        # The rule itself, over whole levels: ring 1 begins at the smallest id
        # of the neighbourhood, under both connectivities and every regime.
        for (l, cells) in ((0, ROOTS), (1, RES1), (2, RES2), (3, RES3[1:7:end]))
            grid = DGG.levelgrid(S, l)
            for conn in CONNECTIVITIES, c in cells
                shell = DGG.neighbors(grid, c, 1; connectivity=conn)
                @test first(shell) == minimum(shell)
            end
        end
    end

    # =======================================================================
    @testset "subtree_border (generic fallback) vs brute force" begin
        # A5 keeps the generic fallback — an A5 subtree is four Hilbert children
        # that cover their parent's area but not its footprint, so there is no
        # digit predicate to read a rim off. What is checked is that the
        # fallback and the definition agree on this hierarchy.
        function brute_border(root, target, conn)
            grid = DGG.levelgrid(S, target)
            l = DGG.level(root)
            return [d for d in DGG.descendants(S, root, target)
                    if any(n -> DGG.ancestor(S, n, l) != root,
                           DGG.neighbors(grid, d, 1; connectivity=conn))]
        end
        for (root, depths) in ((ROOTS[1], 0:3), (RES1[7], 0:3), (RES2[33], 0:2))
            l = DGG.level(root)
            for d in depths, conn in CONNECTIVITIES
                border = DGG.subtree_border(S, root, l + d; connectivity=conn)
                @test border == brute_border(root, l + d, conn)
                @test issorted(border)
                @test border ⊆ DGG.descendants(S, root, l + d)
                interior = DGG.subtree_interior(S, root, l + d; connectivity=conn)
                @test sort(vcat(border, interior)) == DGG.descendants(S, root, l + d)
                @test isempty(intersect(Set(border), Set(interior)))
            end
        end
        @test DGG.subtree_border(S, RES1[3], 1) == [RES1[3]]
        @test isempty(DGG.subtree_interior(S, RES1[3], 1))
        # Every quintant of a face is on its rim: each has a neighbour in
        # another face.
        @test DGG.subtree_border(S, ROOTS[1], 1) == collect(DGG.children(S, ROOTS[1]))
        @test_throws ArgumentError DGG.subtree_border(S, RES2[1], 1)
    end

    # =======================================================================
    # CAP-VALIDATION. `cap_inflation` is a claim about A5's refinement geometry
    # and this is its evidence: the union ratio, i.e. the farthest descendant
    # vertex from a cell's own cap centre over that cell's own radius. A5's
    # pentagon is not a rep-4 tile — the four Hilbert children cover their
    # parent's area but not its footprint — so the ratio is large and
    # essentially level-independent, which is why the factor is 1.75 rather
    # than the shared 1.2.
    # =======================================================================
    @testset "CAP-VALIDATION: subtree union ratios" begin
        inflation = DGG.cap_inflation(S)
        @test inflation == 1.75

        function union_ratios(cells, l, deltas)
            ratios = zeros(length(deltas))
            of_cap = 0.0
            grid = DGG.levelgrid(S, l)
            for c in cells
                cap = DGG.Fallbacks.cell_cap(grid, c)
                raw = maximum(SD(cap.point, p) for p in DGG.cell_boundary(grid, c))
                for (k, delta) in enumerate(deltas)
                    leaf = l + delta
                    leafgrid = DGG.levelgrid(S, leaf)
                    worst = 0.0
                    for d in DGG.descendants(S, c, leaf), p in DGG.cell_boundary(leafgrid, d)
                        distance = SD(cap.point, p)
                        distance > worst && (worst = distance)
                    end
                    ratios[k] = max(ratios[k], worst / raw)
                    of_cap = max(of_cap, worst / DGG.node_extent(S, c).radius)
                end
            end
            return ratios, of_cap
        end

        groups = (
            ("res 0 (all 12)", ROOTS, 0, 1:5),
            ("res 1 (all 60)", RES1, 1, 1:5),
            ("res 2 (all 240)", RES2, 2, 1:5),
            ("res 3 (96 sample)", RES3[1:10:end], 3, 1:6),
            ("res 5 (200 sample)", ordinal_sample(5, 200), 5, 1:6),
            ("res 8 (60 sample)", ordinal_sample(8, 60), 8, 1:6),
            ("res 2 (deep probe)", RES2[1:20:end], 2, 1:8),
        )
        println("\n  A5 CAP-VALIDATION — union ratio (max descendant-vertex distance / " *
                "the cell's own radius before the $(inflation) inflation)")
        worst_ratio = 0.0
        worst_of_cap = 0.0
        worst_extrapolated = 0.0
        convergence_bad = 0
        for (label, cells, l, deltas) in groups
            ratios, of_cap = union_ratios(cells, l, deltas)
            increments = diff(ratios)
            @printf("  %-22s (%4d cells) deltas %d:%d\n    ratios     %s\n    increments %s\n",
                label, length(cells), first(deltas), last(deltas),
                join((@sprintf("%9.5f", r) for r in ratios), " "),
                join((@sprintf("%9.6f", d) for d in increments), " "))
            # Convergence envelope. A5 refines by aperture 4 below the
            # quintants, so the overhang shrinks by ~1/2 per delta — but it can
            # approach its limit from either side (res 0 converges downwards),
            # so the envelope is on the increments' magnitude, not their sign.
            # The first three deltas are exempt: delta 1 spans the 5-way
            # quintant cut, and the lattice's drift only settles into geometric
            # decay once the subtree is a few levels deep.
            for k in 4:length(increments)
                abs(increments[k]) <= max(0.55 * abs(increments[k-1]), 1e-5) ||
                    (convergence_bad += 1)
            end
            # Geometric tail beyond the last measured delta at ratio 0.55:
            # sup <= r_last + |inc_last| * 0.55 / (1 - 0.55).
            worst_extrapolated = max(worst_extrapolated,
                ratios[end] + abs(increments[end]) * 1.23)
            worst_ratio = max(worst_ratio, maximum(ratios))
            worst_of_cap = max(worst_of_cap, of_cap)
        end
        @printf("  worst union ratio %.5f | worst fraction of the wired cap radius %.5f | extrapolated supremum %.5f\n\n",
            worst_ratio, worst_of_cap, worst_extrapolated)

        @test convergence_bad == 0
        # A5 needs far more than the shared 1.2 default — that is the finding,
        # and the reason `cap_inflation` is raised. 1.45363 measured / 1.47078
        # extrapolated on this group set; the gates leave the measurement room
        # to drift a few percent without pretending the budget is looser.
        @test worst_ratio > 1.2 * 1.15
        @test worst_ratio <= 1.50
        @test worst_extrapolated <= inflation * 0.85
        @test worst_of_cap < 0.90                # every descendant inside the cap
    end

    # =======================================================================
    @testset "node_extent covers the subtree" begin
        # The covering law itself, checked directly and all the way down: a
        # cell's node extent contains every vertex of every descendant on a
        # chain to max_level. This is the property `cap_inflation` exists to
        # buy, and the one a sampled ratio can only be evidence for.
        for c in vcat(ROOTS, RES1[1:11:end], RES2[1:53:end])
            cap = DGG.node_extent(S, c)
            @test cap isa GO.UnitSpherical.SphericalCap
            @test cap.radius <= pi / 2        # convex, so vertices imply arcs
            current = c
            while DGG.level(current) < 29
                kids = DGG.children(S, current)
                current = kids[1+(DGG.level(current)%length(kids))]
                grid = DGG.levelgrid(S, DGG.level(current))
                for p in DGG.cell_boundary(grid, current)
                    @test DGG.Fallbacks.cap_contains(cap, p)
                end
            end
        end
    end

    # =======================================================================
    @testset "traits" begin
        @test DGG.cellindextype(S) === A5.A5Cell
        @test DGG.levels(S) === 0:29
        @test DGG.max_level(S) == 29
        @test !DGG.has_sorted_subtrees(S)
        # No `descendant_range` is offered, and asking is a MethodError rather
        # than a wrong answer.
        @test_throws MethodError DGG.descendant_range(S, ROOTS[1], 2)
        @test DGG.cellindextypes(S) === (A5.A5Cell,)
        @test DGG.cap_inflation(S) == 1.75
        for l in 0:29
            g = DGG.levelgrid(S, l)
            @test DGG.system(g) === S
            @test DGG.level(g) == l
        end
        @test occursin("A5Grid", sprint(show, DGG.levelgrid(S, 3)))
        @test sprint(show, S) == "A5System()"
    end

    # =======================================================================
    # SELECTION MODE. A5 is the first real system with
    # `has_sorted_subtrees == false`, so everything below runs on the generic
    # substrate's fallback paths: a `HierarchicalGridCursor` that materialises
    # the positions each node owns, and a `MultiOrderCellSet` ordered by
    # `(level, position)` rather than by curve interval.
    # =======================================================================
    @testset "treeify takes the selection-mode cursor" begin
        for l in (0, 1, 2, 3)
            grid = DGG.levelgrid(S, l)
            tree = DGG.treeify(grid)
            @test tree isa DGG.HierarchicalGridCursor
            @test tree.selection isa Vector{Int}      # the selection path, not a window
            @test Trees.ncells(tree) == DGG.ncells(grid)
            # Full leaf coverage, exactly once each, and leaf index i is
            # position i of the grid.
            @test sort(STI.depth_first_search(Returns(true), tree)) == collect(1:DGG.ncells(grid))
            @test all(i -> ring_points(Trees.getcell(tree, i)) ==
                           ring_points(DGG.cell_polygon(grid, DGG.cellindex(grid, i))),
                1:min(DGG.ncells(grid), 40))
            @test_throws BoundsError Trees.getcell(tree, DGG.ncells(grid) + 1)
        end

        # A partial grid — the shape the selection cursor actually exists for.
        for (ids, bucket) in ((RES3, 0), (RES3[1:13:end], 0), (RES3[1:13:end], 7),
                (DGG.descendants(S, ROOTS[3], 4), 0))
            grid = DGG.PartialGrid(S, DGG.level(first(ids)), ids; bucket_size=bucket)
            tree = DGG.treeify(grid)
            @test tree.selection isa Vector{Int}
            @test Trees.ncells(tree) == length(ids)
            @test sort(STI.depth_first_search(Returns(true), tree)) == collect(1:length(ids))
        end

        # A subtree-rooted chunk answers in chunk-local indices, and its root is
        # the chunk's own cell rather than the synthetic whole-sphere node.
        chunk = DGG.PartialGrid(S, RES1[20], 4)
        chunk_tree = DGG.treeify(chunk)
        @test DGG.Fallbacks.node_cell(chunk_tree) == RES1[20]
        @test Trees.ncells(chunk_tree) == 64
        @test sort(STI.depth_first_search(Returns(true), chunk_tree)) == collect(1:64)
        for child in STI.getchild(chunk_tree)
            @test DGG.level(DGG.Fallbacks.node_cell(child)) == 2
            @test parent(S, DGG.Fallbacks.node_cell(child)) == RES1[20]
        end

        # An empty grid is a leaf with no entries, not an error.
        empty_tree = DGG.treeify(DGG.PartialGrid(S, 3, A5.A5Cell[]))
        @test STI.isleaf(empty_tree)
        @test isempty(STI.child_indices_extents(empty_tree))
        @test Trees.ncells(empty_tree) == 0
    end

    @testset "query through the selection cursor" begin
        cap = GO.UnitSpherical.SphericalCap(GO.UnitSphericalPoint(0.0, 0.0, 1.0), 0.5)
        for l in (1, 2, 3)
            grid = DGG.levelgrid(S, l)
            hits = DGG.query(grid, DGG.Intersects(cap))
            @test issorted(hits)
            @test allunique(hits)
            @test eltype(hits) === A5.A5Cell
            # Sandwich: every cell with a boundary vertex strictly inside the
            # cap must be in the answer, and every cell in the answer must have
            # its own cap meet the target. Pruning never appears in the answer.
            inside = [c for c in [DGG.cellindex(grid, i) for i in 1:DGG.ncells(grid)]
                      if any(p -> SD(cap.point, p) < cap.radius, DGG.cell_boundary(grid, c))]
            @test issubset(Set(inside), Set(hits))
            @test all(c -> DGG.Fallbacks.intersects_cap(cap, DGG.Fallbacks.cell_cap(grid, c)), hits)
            # The system-level form answers the same without a grid in hand.
            @test DGG.query(S, DGG.Intersects(cap); level=l) == hits
        end
        # A lon/lat extent target, lifted to the sphere at the call boundary.
        ext = Extents.Extent(X=(-10.0, 10.0), Y=(40.0, 60.0))
        hits = DGG.query(DGG.levelgrid(S, 3), DGG.Intersects(ext))
        @test !isempty(hits)
        @test issorted(hits)
    end

    @testset "MultiOrderCellSet takes the (level, position) fallback sort" begin
        cap = GO.UnitSpherical.SphericalCap(GO.UnitSphericalPoint(0.0, 0.0, 1.0), 0.5)
        set = DGG.query(S, DGG.MultiOrderCoverage(cap); level=4)
        @test set isa DGG.MultiOrderCellSet
        @test !isempty(set)
        cells = collect(set)
        @test allunique(cells)
        @test all(c -> DGG.level(c) <= 4, cells)
        @test length(unique(DGG.level.(cells))) > 1     # genuinely multi-order
        # WITHOUT sorted subtrees the order is `(level, id)` — levels grouped,
        # ascending inside each — rather than curve order.
        @test issorted(cells; by=c -> (DGG.level(c), c))
        # ... and the reported keys are each cell's position within its own
        # level, which is what the fallback documents.
        @test DGG.Fallbacks.curve_keys(set) ==
              [DGG.cellposition(DGG.levelgrid(S, DGG.level(c)), c) for c in cells]
        # There are no position intervals to expand to, and asking says so.
        @test_throws ArgumentError DGG.level_ranges(set, 4)

        # COVERAGE, which is the property the set is for: every level-4 cell
        # that meets the target has itself or an ancestor in the set.
        grid4 = DGG.levelgrid(S, 4)
        owners = Set(cells)
        @test all(DGG.query(grid4, DGG.Intersects(cap))) do c
            any(l -> DGG.ancestor(S, c, l) in owners, 0:4)
        end
    end

    # =======================================================================
    @testset "conformance" begin
        for l in (0, 1, 3, 9)
            test_grid_interface(DGG.levelgrid(S, l); label="A5Grid(res $l)")
        end
        test_hierarchical_system(S)
        # Levels the seeded sampler will not have drawn, so the hierarchy laws
        # meet all three regimes rather than whichever four came up.
        test_hierarchical_system(S; levels=0:4, n_levels=5, label="A5System (regimes 0-4)")
    end
end

end # module A5TestSuite
