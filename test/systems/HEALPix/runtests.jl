# ---------------------------------------------------------------------------
# T6 — HEALPix system tests.
#
# Three kinds of test, and the distinction matters when one fails:
#
#   1. ORACLE. Everything the port could have transcribed wrong is checked
#      against Healpix.jl, an independent implementation of the same standard:
#      the nested/ring/xyf codecs both ways, point location, and pixel centres.
#      The implementation itself never calls Healpix.jl — `src/systems/HEALPix/`
#      is pure Julia — which is exactly what makes it an oracle rather than a
#      tautology. Adapted from the pre-redesign `test/HEALPix/` suites.
#
#   2. CONTRACT. The two conformance suites from
#      `DiscreteGlobalGridsConformanceTesting`, with default kwargs.
#
#   3. STRUCTURAL. The things the oracle has no opinion about because they are
#      this package's own design: the exact subtree cap's covering margin, the
#      Morton rim walk, the `Edge()` restriction, and the 0-based/1-based
#      conventions.
#
# Run directly:
#     julia --project=test --startup-file=no test/systems/HEALPix/runtests.jl
# (T7 wires it into `test/runtests.jl`.)
# ---------------------------------------------------------------------------

module HEALPixSystemTests

using Test
using Random

using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
const HP = DiscreteGlobalGrids.HEALPix

using DiscreteGlobalGridsConformanceTesting

import GeometryOps as GO
const US = GO.UnitSpherical

# The independent implementation. Note the capitalisation: `Healpix` is the
# registered package, `HP` is this package's own submodule.
import Healpix

const SYS = HP.HEALPixSystem()

"A deterministic pixel sample: all of them when small, a seeded draw when not."
function pixel_sample(level::Integer, n::Integer = 64)
    npix = 12 * 4^level
    npix <= n && return collect(0:(npix - 1))
    rng = MersenneTwister(911_017 + level)
    return sort!(unique(rand(rng, 0:(npix - 1), n)))
end

@testset "HEALPix system (T6)" begin

# =========================================================================
# 1. Oracle: the codecs against Healpix.jl
# =========================================================================

@testset "codecs vs Healpix.jl" begin
    for level in 0:6
        nside = 1 << level
        res = Healpix.Resolution(nside)
        @testset "level $level (nside $nside)" begin
            for p in pixel_sample(level)
                # nested -> xyf, and back
                mine = HP.nested_to_xyf(p, nside)
                @test mine == Healpix.pix2xyfNest(res, p + 1)
                @test HP.xyf_to_nested(mine..., nside) == p

                # nested <-> ring. Healpix.jl numbers pixels 1-based in BOTH
                # schemes; ours is 0-based nested and 1-based ring, so exactly
                # one `+ 1` appears here and it is on the nested side.
                r = HP.nested_to_ring(p, nside)
                @test r == Healpix.nest2ring(res, p + 1)
                @test HP.ring_to_nested(r, nside) == p

                # ring -> xyf, and back
                rxyf = HP.ring_to_xyf(r, nside)
                @test rxyf == Healpix.pix2xyfRing(res, r)
                @test HP.xyf_to_ring(rxyf..., nside) == r
                @test HP.xyf_to_ring(mine..., nside) == r
            end
        end
    end
end

@testset "codecs at a deep level" begin
    # The Int64 codec well past anything exhaustive, to catch a shift that only
    # overflows when the Morton code is wide.
    for level in (12, 20, 29)
        nside = 1 << level
        npix = 12 * Int64(4)^level
        rng = MersenneTwister(4242 + level)
        for _ in 1:64
            p = rand(rng, 0:(npix - 1))
            ix, iy, f = HP.nested_to_xyf(p, nside)
            @test 0 <= ix < nside && 0 <= iy < nside && 0 <= f <= 11
            @test HP.xyf_to_nested(ix, iy, f, nside) == p
        end
    end
    # Healpix.jl agrees at level 12 too (its own Resolution tops out at 29).
    let level = 12, nside = 1 << level, res = Healpix.Resolution(nside)
        rng = MersenneTwister(99)
        for _ in 1:64
            p = rand(rng, 0:(12 * 4^level - 1))
            @test HP.nested_to_xyf(p, nside) == Healpix.pix2xyfNest(res, p + 1)
            @test HP.nested_to_ring(p, nside) == Healpix.nest2ring(res, p + 1)
        end
    end
end

# =========================================================================
# 2. Oracle: point location and centres
# =========================================================================

# Points deliberately placed ON a cell boundary: the poles, the antimeridian,
# the ±2/3 belt boundaries where `point_to_xyf` switches branch, and the 0.99
# cancellation-guard latitude.
BELT = asind(2 / 3)
GUARD = asind(0.99)
DEGENERATE_POINTS = [
    (0.0, 90.0), (0.0, -90.0), (180.0, 0.0), (-180.0, 0.0), (0.0, 0.0),
    (45.0, BELT), (45.0, -BELT), (123.0, GUARD), (-77.0, -GUARD),
    (359.9, 89.99), (0.0, BELT + 1e-9), (0.0, BELT - 1e-9)]

@testset "cellat vs Healpix.vec2pixNest (interior points)" begin
    # `vec2pixNest`, not `ang2pixNest`, is the sharp oracle: it takes the SAME
    # Cartesian point this port does, so the two implementations see
    # bit-identical input and no difference can be blamed on the (lon, lat)
    # conversion. Points in the interior of a cell must agree exactly.
    rng = MersenneTwister(20260813)
    sample = [(360 * rand(rng) - 180, asind(2 * rand(rng) - 1)) for _ in 1:400]
    for level in 0:8
        res = Healpix.Resolution(1 << level)
        grid = levelgrid(SYS, level)
        @testset "level $level" begin
            for (lon, lat) in sample
                p = US.UnitSphereFromGeographic()((lon, lat))
                expected = Healpix.vec2pixNest(res, p[1], p[2], p[3]) - 1
                @test HP.point_to_nested(p, level) == expected
                @test cellat(grid, p) == LevelIndex(level, expected)
            end
        end
    end
end

@testset "cellat on cell boundaries: deterministic, incident, self-consistent" begin
    # A point exactly on a shared boundary is legitimately contained by every
    # cell that meets there, so "which one" is a TIE, and the interface asks
    # only that the answer be deterministic and documented — not that it match
    # another implementation. Healpix.jl breaks two of these ties the other way:
    # (1, 0, 0) and (45°, -asin(2/3)) are exact 4-way lattice corners, and both
    # libraries are self-consistent (each maps its own chosen pixel's centre
    # back to that pixel). So this asserts the properties that are actually
    # contractual, and does NOT assert equality with Healpix.jl.
    for level in 0:8
        grid = levelgrid(SYS, level)
        for (lon, lat) in DEGENERATE_POINTS
            p = US.UnitSphereFromGeographic()((lon, lat))
            c = cellat(grid, p)
            @test c isa LevelIndex
            @test DGG.level(c) == level
            @test cellat(grid, p) == c                    # deterministic
            @test cellposition(grid, c) !== nothing       # a real cell of the grid
            # incident: the point is inside the cell's own extent, and the cell
            # is the one this port's centre round-trips to.
            cap = node_extent(SYS, c)
            @test US.spherical_distance(cap.point, p) <= cap.radius
            @test cellat(grid, cell_centroid(grid, c)) == c
        end
    end
end

@testset "cellat vs Healpix.ang2pixNest (lon/lat entry point)" begin
    # The (lon, lat) wrapper against Healpix.jl's own (θ, φ) entry point. Only
    # the random sample: a point placed exactly on a pixel boundary is a tie
    # that the two paths' different roundings of `z` may break differently, and
    # that is a property of the input conversion, not of the port. Ties are
    # covered above, where both sides get the identical Cartesian point.
    rng = MersenneTwister(31337)
    for level in 0:8
        res = Healpix.Resolution(1 << level)
        grid = levelgrid(SYS, level)
        for _ in 1:200
            lon = 360 * rand(rng) - 180
            lat = asind(2 * rand(rng) - 1)
            expected = Healpix.ang2pixNest(res, deg2rad(90 - lat), deg2rad(mod(lon, 360))) - 1
            @test cellat(grid, lon, lat) == LevelIndex(level, expected)
        end
    end
end

@testset "cell_centroid vs Healpix.pix2vecNest" begin
    worst = 0.0
    for level in 0:6
        nside = 1 << level
        res = Healpix.Resolution(nside)
        grid = levelgrid(SYS, level)
        for p in pixel_sample(level)
            mine = cell_centroid(grid, LevelIndex(level, p))
            theirs = Healpix.pix2vecNest(res, p + 1)
            worst = max(worst, maximum(abs.(Tuple(mine) .- Tuple(theirs))))
        end
    end
    @info "cell_centroid vs Healpix.pix2vecNest: worst coordinate deviation" worst
    @test worst < 1e-14
end

@testset "cell corners vs Healpix.boundariesRing" begin
    # The four CORNERS of the ring, in the same (north, west, south, east)
    # order `boundariesRing` emits. `cell_boundary` densifies each edge into
    # `BOUNDARY_SEGMENTS` segments, so the corners are every 8th vertex.
    step = HP.BOUNDARY_SEGMENTS
    worst = 0.0
    for level in 2:6
        nside = 1 << level
        res = Healpix.Resolution(nside)
        grid = levelgrid(SYS, level)
        for p in pixel_sample(level, 32)
            mine = cell_boundary(grid, LevelIndex(level, p))
            @test length(mine) == 4 * step
            # (4, 3) matrix, one row per corner
            theirs = Healpix.boundariesRing(res, Healpix.nest2ring(res, p + 1), 1, Float64)
            for k in 1:4
                corner = mine[1 + (k - 1) * step]
                worst = max(worst, maximum(abs.(Tuple(corner) .- Tuple(theirs[k, :]))))
            end
        end
    end
    @info "cell corners vs Healpix.boundariesRing: worst coordinate deviation" worst
    @test worst < 1e-14
end

@testset "densified edge points lie on the true chart edge" begin
    # The densification is not an interpolation of the corners: every vertex is
    # the chart evaluated at a lattice-fraction point, so it is ON the true
    # pixel edge, not on the chord. Checked by evaluating the chart directly.
    for level in (0, 2, 5)
        nside = 1 << level
        grid = levelgrid(SYS, level)
        for p in pixel_sample(level, 16)
            ix, iy, f = HP.nested_to_xyf(p, nside)
            ring = cell_boundary(grid, LevelIndex(level, p))
            @test ring == HP._perimeter_points(ix, iy, f, nside, HP.BOUNDARY_SEGMENTS)
            @test all(q -> abs(hypot(q[1], q[2], q[3]) - 1) < 1e-15, ring)
        end
    end
end

@testset "shared edges are bit-identical between neighbours" begin
    # What makes the tessellation exact rather than merely consistent: a
    # densification point of one pixel is the SAME Float64 triple in its
    # neighbour's ring. Checked within a face, where the two pixels share the
    # chart and the argument arithmetic must agree bitwise.
    level = 4
    nside = 1 << level
    grid = levelgrid(SYS, level)
    shared = 0
    for ix in 1:6, iy in 1:6
        a = LevelIndex(level, HP.xyf_to_nested(ix, iy, 0, nside))
        b = LevelIndex(level, HP.xyf_to_nested(ix - 1, iy, 0, nside))
        ra = Set(cell_boundary(grid, a))
        rb = Set(cell_boundary(grid, b))
        common = intersect(ra, rb)
        # the shared edge contributes its BOUNDARY_SEGMENTS start points plus
        # the far corner, which is the next edge's start in the other ring
        @test length(common) >= HP.BOUNDARY_SEGMENTS
        shared += length(common)
    end
    @test shared > 0
end

# =========================================================================
# 3. Structural: id conventions
# =========================================================================

@testset "0-based nested ids, 1-based positions" begin
    for level in (0, 1, 4)
        grid = levelgrid(SYS, level)
        @test ncells(grid) == 12 * 4^level
        @test cellindex(grid, 1) == LevelIndex(level, 0)          # position 1 -> id 0
        @test cellindex(grid, ncells(grid)) == LevelIndex(level, 12 * 4^level - 1)
        @test rawid(cellindex(grid, 1)) == 0
        for i in (1, 7, ncells(grid))
            @test cellposition(grid, cellindex(grid, i)) == i
        end
        # A cell from another level is simply not in this grid.
        @test cellposition(grid, LevelIndex(level + 1, 0)) === nothing
        # An id no pixel has is not in it either.
        @test cellposition(grid, LevelIndex(level, 12 * 4^level)) === nothing
        @test cellposition(grid, LevelIndex(level, -1)) === nothing
        @test_throws BoundsError cellindex(grid, 0)
        @test_throws BoundsError cellindex(grid, ncells(grid) + 1)
    end
end

@testset "ring reindex round-trip" begin
    @test cellindextype(SYS) === LevelIndex
    @test cellindextypes(SYS) == (LevelIndex, HP.HEALPixRingIndex)
    for level in 0:5
        nside = 1 << level
        res = Healpix.Resolution(nside)
        grid = levelgrid(SYS, level)
        for p in pixel_sample(level, 32)
            c = LevelIndex(level, p)
            r = reindex(HP.HEALPixRingIndex, SYS, c)
            @test DGG.level(r) == level
            # 1-BASED on the ring side: equal to Healpix.jl's own ring number.
            @test rawid(r) == Healpix.nest2ring(res, p + 1)
            @test reindex(LevelIndex, SYS, r) == c
            @test reindex(LevelIndex, SYS, c) === c
            # A ring-named cell still finds its position in the nested grid.
            @test cellposition(grid, r) == cellposition(grid, c)
        end
    end
end

# =========================================================================
# 4. Structural: hierarchy arithmetic
# =========================================================================

@testset "hierarchy arithmetic" begin
    @test length(rootcells(SYS)) == 12
    @test issorted(rootcells(SYS))
    @test levels(SYS) == 0:29
    @test max_level(SYS) == 29
    @test has_sorted_subtrees(SYS)
    @test max_neighbors(SYS) == 8
    @test max_neighbors(SYS, Vertex()) == 8
    @test max_neighbors(SYS, Edge()) == 4

    @test_throws ArgumentError parent(SYS, LevelIndex(0, 0))
    @test_throws ArgumentError children(SYS, LevelIndex(29, 0))
    @test_throws ArgumentError levelgrid(SYS, -1)
    @test_throws ArgumentError levelgrid(SYS, 30)

    for level in 0:4, p in pixel_sample(level, 24)
        c = LevelIndex(level, p)
        kids = children(SYS, c)
        @test length(kids) == 4                    # always 4: no pentagons here
        @test kids == [LevelIndex(level + 1, 4p + k) for k in 0:3]
        @test all(k -> parent(SYS, k) == c, kids)
        # `ancestor` in one shift agrees with walking `parent`
        walked = c
        for l in (level - 1):-1:0
            walked = parent(SYS, walked)
            @test ancestor(SYS, c, l) == walked
        end
        @test ancestor(SYS, c, level) == c
        @test_throws ArgumentError ancestor(SYS, c, level + 1)
    end
end

@testset "descendant_range is exact and hole-free" begin
    for level in 0:3, p in pixel_sample(level, 16)
        c = LevelIndex(level, p)
        @test descendant_range(SYS, c, level) == (p + 1):(p + 1)
        @test_throws ArgumentError descendant_range(SYS, c, level - 1)
        for d in 1:3
            target = level + d
            r = descendant_range(SYS, c, target)
            grid = levelgrid(SYS, target)
            @test length(r) == 4^d
            # Every position in the range is a descendant, and every descendant
            # is in the range — the two-sided contract.
            actual = [cellposition(grid, x) for x in descendants(SYS, c, target)]
            @test sort(actual) == collect(r)
            @test all(x -> ancestor(SYS, cellindex(grid, x), level) == c, r)
        end
        # Sibling ranges tile the parent's, in order.
        if level < 3
            kid_ranges = [descendant_range(SYS, k, level + 1) for k in children(SYS, c)]
            @test reduce(vcat, collect.(kid_ranges)) ==
                  collect(descendant_range(SYS, c, level + 1))
        end
    end
end

# =========================================================================
# 5. Structural: the exact subtree cap
# =========================================================================

@testset "node_extent covers the true pixel, with margin" begin
    # The conformance suite samples DESCENDANT VERTICES; this samples the
    # continuous truth those vertices lie on, 32x finer than the cap is built
    # from. A negative worst-overshoot means the cap is not merely passing the
    # sampled law but bounding the real region.
    worst_overshoot = -Inf
    max_radius = 0.0
    for level in 0:6
        nside = 1 << level
        for p in pixel_sample(level, 512)
            ix, iy, f = HP.nested_to_xyf(p, nside)
            cap = node_extent(SYS, LevelIndex(level, p))
            max_radius = max(max_radius, cap.radius)
            for q in HP._perimeter_points(ix, iy, f, nside, 256)
                worst_overshoot = max(worst_overshoot,
                    US.spherical_distance(cap.point, q) - cap.radius)
            end
        end
    end
    @info "node_extent vs densely sampled true perimeter" worst_overshoot max_radius
    @test worst_overshoot < 0            # strictly inside: slack never consumed
    @test max_radius <= π / 2            # geodesically convex, as the harness asserts
end

@testset "node_extent covers deep descendants" begin
    # The covering law stated directly: descendants many levels down, boundary
    # points in the ANCESTOR's extent. The conformance suite does this on a
    # sample; here it is exhaustive over a whole subtree.
    for (level, p, depth) in ((0, 0, 5), (0, 4, 5), (0, 11, 5), (2, 37, 4), (3, 500, 3))
        c = LevelIndex(level, p)
        cap = node_extent(SYS, c)
        target = level + depth
        grid = levelgrid(SYS, target)
        worst = -Inf
        for pos in descendant_range(SYS, c, target)
            d = cellindex(grid, pos)
            for q in cell_boundary(grid, d)
                worst = max(worst, US.spherical_distance(cap.point, q) - cap.radius)
            end
            # a descendant's own extent is contained in its ancestor's, too
            dcap = node_extent(SYS, d)
            @test US.spherical_distance(cap.point, dcap.point) + dcap.radius <=
                  cap.radius + 1e-12
        end
        @test worst < 0
    end
end

# =========================================================================
# 6. Structural: neighbours
# =========================================================================

@testset "neighbours" begin
    for level in 0:4
        grid = levelgrid(SYS, level)
        npix = 12 * 4^level
        counts = Dict{Int,Int}()
        for p in 0:(npix - 1)
            c = LevelIndex(level, p)
            vs = collect(neighbors(grid, c, 1))
            es = collect(neighbors(grid, c, 1; connectivity = Edge()))

            @test eltype(vs) === LevelIndex
            @test issorted(vs) && allunique(vs)
            @test issorted(es) && allunique(es)
            @test !(c in vs)
            @test issubset(Set(es), Set(vs))        # Edge() restricts Vertex()
            @test length(vs) <= 8
            @test length(es) <= 4
            counts[length(vs)] = get(counts, length(vs), 0) + 1

            # symmetry, both connectivities, exhaustively over the level
            for n in vs
                @test c in collect(neighbors(grid, n, 1))
            end
            for n in es
                @test c in collect(neighbors(grid, n, 1; connectivity = Edge()))
            end
            @test isempty(collect(neighbors(grid, c, 0)))
            @test collect(ring(grid, c, 0)) == [c]
        end
        # The 24 pixels on a degree-3 base-tiling vertex have seven neighbours.
        #
        # Level 0 is its own case, and measured rather than assumed: with
        # nside == 1 a base pixel IS a face, every lattice offset wraps through
        # the face tables, and all twelve come out with exactly six neighbours
        # (the `E` and `W` diagonals are absent for every face). So the base
        # tiling is 6-regular, and the "24 pixels with 7" rule starts at
        # level 1.
        if level == 0
            @test counts == Dict(6 => 12)
        else
            @test counts[7] == 24
            @test counts[8] == npix - 24
            @test length(counts) == 2
        end
    end
end

@testset "ring shells and neighbour discs agree" begin
    for level in (2, 3), conn in (Vertex(), Edge())
        grid = levelgrid(SYS, level)
        for p in pixel_sample(level, 24)
            c = LevelIndex(level, p)
            seen = Set{LevelIndex}()
            for k in 1:3
                shell = Set(ring(grid, c, k; connectivity = conn))
                @test isempty(intersect(shell, seen))
                union!(seen, shell)
                @test Set(neighbors(grid, c, k; connectivity = conn)) == seen
            end
            @test !(c in seen)
        end
    end
end

# =========================================================================
# 7. Structural: the Morton rim walk
# =========================================================================

@testset "subtree_border" begin
    for level in 0:2, p in pixel_sample(level, 12)
        c = LevelIndex(level, p)
        @test HP.subtree_border(SYS, c, level) == [c]
        for depth in 1:5
            target = level + depth
            rim = HP.subtree_border(SYS, c, target)
            s = 1 << depth
            @test length(rim) == 4s - 4
            @test issorted(rim) && allunique(rim)

            # Against the definition: a descendant is on the rim iff it has a
            # neighbour outside the subtree. Computed here by neighbour query,
            # which is precisely what the rim walk exists not to do.
            grid = levelgrid(SYS, target)
            inside = Set(descendant_range(SYS, c, target))
            brute = LevelIndex[]
            for pos in inside
                d = cellindex(grid, pos)
                any(n -> !(cellposition(grid, n) in inside),
                    neighbors(grid, d, 1)) && push!(brute, d)
            end
            @test rim == sort!(brute)
            # and the same answer under Edge() connectivity
            brute_e = LevelIndex[]
            for pos in inside
                d = cellindex(grid, pos)
                any(n -> !(cellposition(grid, n) in inside),
                    neighbors(grid, d, 1; connectivity = Edge())) && push!(brute_e, d)
            end
            @test rim == sort!(brute_e)
        end
    end
    @test_throws ArgumentError HP.subtree_border(SYS, LevelIndex(3, 0), 2)
end

# =========================================================================
# 8. Geometry sanity
# =========================================================================

@testset "cell_area is exactly equal-area" begin
    # The override: the closed form, identical for every pixel of a level, and
    # summing to the whole sphere.
    for level in 0:5
        grid = levelgrid(SYS, level)
        expected = 4π / (12 * 4^level)
        for p in pixel_sample(level, 64)
            @test cell_area(grid, LevelIndex(level, p)) == expected
        end
        # Exact in closed form; the running sum carries O(ncells * eps) of
        # accumulation error, which is why the tolerance is not 0.
        @test ncells(grid) * cell_area(grid, cellindex(grid, 1)) ≈ 4π rtol = 1e-15
        total = sum(cell_area(grid, cellindex(grid, i)) for i in 1:ncells(grid))
        @test total ≈ 4π rtol = 1e-10
    end
end

@testset "boundary polygon area approaches the exact area" begin
    # The densified ring is an approximation of the equal-area truth, and this
    # pins how good it is. A POSITIVE signed area is also the orientation
    # check: counter-clockwise seen from outside, as `cell_boundary` requires.
    # The bound is the `BOUNDARY_SEGMENTS = 8` column of the table in
    # `system.jl`; if that constant changes, this number moves with it.
    for level in 0:5
        grid = levelgrid(SYS, level)
        exact = 4π / (12 * 4^level)
        worst = 0.0
        for p in pixel_sample(level, 64)
            poly = cell_polygon(grid, LevelIndex(level, p))
            a = GO.area(GO.Spherical(radius = 1.0), poly)
            @test a > 0                       # counter-clockwise from outside
            worst = max(worst, abs(a - exact) / exact)
        end
        @test worst < 3e-3
    end
end

@testset "the full sphere is partitioned by the polygons" begin
    # Total DENSIFIED POLYGON area over a whole level, against 4π: the cells
    # tile the sphere with no gaps or overlaps, to the densification's accuracy.
    for level in (1, 2, 3)
        grid = levelgrid(SYS, level)
        total = sum(GO.area(GO.Spherical(radius = 1.0),
                            cell_polygon(grid, cellindex(grid, i)))
                    for i in 1:ncells(grid))
        @test isapprox(total, 4π; rtol = 3e-3)
    end
end

@testset "cellat inverts cell_centroid everywhere" begin
    for level in 0:5
        grid = levelgrid(SYS, level)
        for p in pixel_sample(level, 256)
            c = LevelIndex(level, p)
            @test cellat(grid, cell_centroid(grid, c)) == c
        end
    end
end

# =========================================================================
# 9. Contract: the conformance suites, default kwargs
# =========================================================================

@testset "conformance" begin
    for l in (0, 1, 3)
        test_grid_interface(levelgrid(SYS, l); label = "HEALPixGrid(level=$l)")
    end
    test_hierarchical_system(SYS)
end

end # @testset "HEALPix system (T6)"

end # module HEALPixSystemTests
