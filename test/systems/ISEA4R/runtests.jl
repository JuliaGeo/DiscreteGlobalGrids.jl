# ---------------------------------------------------------------------------
# T12 — ISEA4R system tests.
#
# Three kinds of test, and the distinction matters when one fails:
#
#   1. ORACLE. ISEA4R has no independent implementation to check against — the
#      ten-diamond layout is this package's own convention with no external
#      fixture behind it (see the `ISEA4R` module docstring). What plays the
#      oracle's role instead is the pre-redesign `test/ISEA4R/` suite, whose
#      vectors are restated here INDEPENDENTLY of the source: the pinned
#      layout literals are written out a second time, the structural facts (the
#      ten pairs partition the twenty faces, every seam is a real icosahedron
#      edge, the northern/southern apexes are bases 0 and 11) are re-derived
#      rather than restated, and the chart's own claims — equal-areaness, seam
#      ownership, bit-exact nesting — are measured against closed forms written
#      out here.
#
#      HOW MUCH OF EACH LEVEL. Levels 0-3 (10, 40, 160 and 640 cells) are walked
#      EXHAUSTIVELY almost everywhere below; `cell_sample` is used where a
#      deeper level is wanted and takes a seeded, reproducible draw above its
#      cap. Where a testset iterates a whole level it says so.
#
#   2. CONTRACT. The two conformance suites from
#      `DiscreteGlobalGridsConformanceTesting`, with default kwargs.
#
#   3. STRUCTURAL. The things nothing else has an opinion about because they
#      are this port's own design: the seam topology and its 9-neighbour corner
#      cells, the rotational start, the exact subtree cap, and the 0-based id /
#      1-based position convention.
#
# `test/runtests.jl` includes this file (T13 flips that line on); it also runs
# standalone:
#     julia --project=test --startup-file=no test/systems/ISEA4R/runtests.jl
# ---------------------------------------------------------------------------

module ISEA4RSystemTests

using Test
using Printf
using Random

using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
const I4 = DiscreteGlobalGrids.ISEA4R
const ISEA = DiscreteGlobalGrids.ISEA

using DiscreteGlobalGridsConformanceTesting

import GeometryOps as GO
const US = GO.UnitSpherical

const SYS = I4.ISEA4RSystem()

# Recorded so the numbers land in the test log (and in the milestone report).
const MEASURED = Dict{String,Float64}()
record!(key, value) = (MEASURED[key] = max(get(MEASURED, key, -Inf), value))

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

sd(a, b) = US.spherical_distance(a, b)

"`ISEA.VERTICES[v + 1]` as a `UnitSphericalPoint`."
vpoint(v) = (w = ISEA.VERTICES[v + 1]; GO.UnitSphericalPoint(w[1], w[2], w[3]))

"""
A deterministic cell sample: ALL `10 * 4^level` ids when the level fits under
`n`, a seeded draw of at most `n` of them when it does not. At the default
`n = 64` that is every cell at levels 0 (10) and 1 (40), and a draw from level 2
(160) up.
"""
function cell_sample(level::Integer, n::Integer = 64)
    ncell = 10 * 4^level
    ncell <= n && return collect(0:(ncell - 1))
    rng = MersenneTwister(4_120_017 + level)
    return sort!(unique(rand(rng, 0:(ncell - 1), n)))
end

cross3(a, b) = (a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])

"Signed area of a ring seen from OUTSIDE the sphere; positive iff CCW."
function ccw_measure(corners)
    acc = (0.0, 0.0, 0.0)
    n = length(corners)
    for i in 1:n
        acc = acc .+ cross3(corners[i], corners[i % n + 1])
    end
    outward = reduce((a, b) -> a .+ Tuple(b), corners; init = (0.0, 0.0, 0.0))
    return sum(acc .* outward)
end

"Area of a spherical triangle by L'Huilier's formula (unsigned)."
function lhuilier(a, b, c)
    A = sd(b, c); B = sd(a, c); C = sd(a, b)
    s = (A + B + C) / 2
    t = tan(s / 2) * tan((s - A) / 2) * tan((s - B) / 2) * tan((s - C) / 2)
    return 4 * atan(sqrt(max(t, 0.0)))
end

"""
Area of a closed spherical ring, as a fan of L'Huilier triangles from the ring's
own centroid. The centroid apex is what makes this correct for the densified,
gently curved rings below — a fan from a *vertex* folds over on them, and
`GO.area` on a densified ring of hundreds of nearly collinear vertices
disagrees with both by several percent.
"""
function ring_area(pts)
    c = reduce((a, p) -> a .+ (p[1], p[2], p[3]), pts; init = (0.0, 0.0, 0.0))
    c = c ./ sqrt(sum(abs2, c))
    o = GO.UnitSphericalPoint(c[1], c[2], c[3])
    return sum(lhuilier(o, pts[i], pts[mod1(i + 1, length(pts))]) for i in eachindex(pts))
end

"The upper-half affine evaluation, written out — the seam-ownership reference."
function upper_half_point(x, y, diamond)
    dm = I4.DIAMONDS[diamond + 1]
    p = ISEA.snyder_inv_xyz(dm.upper, dm.cP0 + x * dm.aP + y * dm.bP)
    return GO.UnitSphericalPoint(p[1], p[2], p[3])
end

"The lower-half affine evaluation, written out."
function lower_half_point(x, y, diamond)
    dm = I4.DIAMONDS[diamond + 1]
    p = ISEA.snyder_inv_xyz(dm.lower, dm.cQ0 + x * dm.aQ + y * dm.bQ)
    return GO.UnitSphericalPoint(p[1], p[2], p[3])
end

"""
A hashable key for a boundary point, with `-0.0` normalised to `0.0`.

`Set` and `Dict` compare with `isequal`, under which `-0.0` and `0.0` are
DISTINCT even though `==` says they are equal. Two closed forms that agree on
every bit but the sign of a zero coordinate therefore vanish from a hashed
incidence test, silently and only at the handful of points where a coordinate is
exactly zero. That is a real trap in this family — the sibling S2 port hit it on
its cube seams — so every hashed comparison below goes through this, and
`signed zeros cannot break a hashed incidence test` pins the two facts that make
ISEA4R immune to it today. `x + 0.0` is the identity on every `Float64` except
`-0.0`, which it maps to `0.0`.
"""
vkey(p) = (p[1] + 0.0, p[2] + 0.0, p[3] + 0.0)

"""
How many boundary points two cells share, to a tolerance.

A tolerance rather than a `Set` intersection, and that is load-bearing: shared
lattice points are BIT-identical only within a diamond. Across a diamond rim the
two sides are two independent developments of the same icosahedron edge and
agree to ~4e-15 rad but to NO bits at all, so an exact set intersection would
report a seam neighbour as sharing nothing. The measured worst mismatch across a
rim is recorded by `cross-diamond borders line up point for point` below, three
orders under this tolerance.
"""
function shared_boundary_points(g, a, b; tol = 1e-12)
    ra = cell_boundary(g, a)
    rb = cell_boundary(g, b)
    return count(p -> any(q -> sd(p, q) <= tol, rb), ra)
end

"The `(ix, iy, diamond)` of a cell id at a level."
xyd_of(c) = I4.morton_to_xyd(rawid(c), 1 << level(c))

# --------------------------------------------------------------------------
# The layout pin, restated independently of `diamonds.jl`'s own literals.
#
# `diamonds.jl` derives the table and asserts it against a pin it carries
# itself, which catches upstream `ISEA` drift but not a coordinated edit of both
# — so the same numbers are written out a second time here.
# --------------------------------------------------------------------------

# (v00, v10, v11, v01) — chart (0,0), (1,0), (1,1), (0,1); seam = (v00, v11).
const EXPECTED_VERTS = (
    (1, 6, 2, 0), (2, 7, 3, 0), (3, 8, 4, 0), (4, 9, 5, 0), (5, 10, 1, 0),
    (10, 11, 6, 1), (6, 11, 7, 2), (7, 11, 8, 3), (8, 11, 9, 4), (9, 11, 10, 5),
)
const EXPECTED_FACES = (
    (1, 4), (5, 11), (6, 12), (2, 8), (0, 3),
    (7, 13), (10, 16), (17, 19), (14, 18), (9, 15),
)

@testset "ISEA4R system (T12)" begin

# =========================================================================
# 1. Oracle: the ten-diamond layout table
# =========================================================================

@testset "the ten-diamond table is the pinned one" begin
    @test length(I4.DIAMONDS) == 10
    for d in 0:9
        dm = I4.DIAMONDS[d + 1]
        @test dm.verts == EXPECTED_VERTS[d + 1]
        @test (dm.upper, dm.lower) == EXPECTED_FACES[d + 1]
    end

    # Structural facts, re-derived rather than restated.
    faces = sort!(vcat([dm.upper for dm in I4.DIAMONDS], [dm.lower for dm in I4.DIAMONDS]))
    @test faces == collect(0:19)                       # the twenty faces, once each
    for d in 0:9
        dm = I4.DIAMONDS[d + 1]
        v00, v10, v11, v01 = dm.verts
        @test length(unique(dm.verts)) == 4
        # Each half's three corners are the face's three corners.
        @test issetequal((v00, v01, v11), ISEA.FACE_TRIPLES[dm.upper + 1])
        @test issetequal((v00, v10, v11), ISEA.FACE_TRIPLES[dm.lower + 1])
        # The seam is a real icosahedron edge.
        @test abs(ISEA.vdot(ISEA.VERTICES[v00 + 1], ISEA.VERTICES[v11 + 1]) - ISEA.ADJ_DOT) < 1e-9
        # Northern diamonds apex on base 0, southern ones on base 11.
        @test (d < 5) == (v01 == 0)
        @test (d >= 5) == (v10 == 11)
        # Equal-area, exactly, on both halves.
        @test imag(conj(dm.aP) * dm.bP) == 2pi / 5
        @test imag(conj(dm.aQ) * dm.bQ) == 2pi / 5
    end
end

@testset "chart corners land on the icosahedron vertices" begin
    worst = 0.0
    for d in 0:9
        verts = I4.DIAMONDS[d + 1].verts
        for ((x, y), v) in (((0.0, 0.0), verts[1]), ((1.0, 0.0), verts[2]),
                            ((1.0, 1.0), verts[3]), ((0.0, 1.0), verts[4]))
            worst = max(worst, sd(I4.xyd_to_point(x, y, d), vpoint(v)))
        end
    end
    @test worst < 1e-14
    record!("chart corner vs ISEA.VERTICES (rad)", worst)
end

@testset "the two halves agree on the seam" begin
    worst = 0.0
    for d in 0:9, t in range(0.0, 1.0; length = 257)
        worst = max(worst, sd(upper_half_point(t, t, d), lower_half_point(t, t, d)))
    end
    # The floor is Snyder's own round-off: two independent developments of one
    # icosahedron edge can agree only to the accuracy of `snyder_inv_xyz`. That
    # is exactly why the chart picks ONE of them for the seam.
    @test worst < 1e-13
    record!("upper/lower half disagreement on the seam (rad)", worst)
end

@testset "evaluations of one icosahedron vertex cluster to round-off" begin
    by_vertex = Dict{Int,Vector{GO.UnitSphericalPoint{Float64}}}()
    for d in 0:9
        verts = I4.DIAMONDS[d + 1].verts
        for ((x, y), v) in (((0.0, 0.0), verts[1]), ((1.0, 0.0), verts[2]),
                            ((1.0, 1.0), verts[3]), ((0.0, 1.0), verts[4]))
            push!(get!(by_vertex, v, GO.UnitSphericalPoint{Float64}[]),
                  I4.xyd_to_point(x, y, d))
        end
    end
    # All twelve vertices are hit; bases 0 and 11 five times (they are apexes of
    # five diamonds each), the other ten three times. This is the SAME count the
    # nine-neighbour corner cells come from — see the seam-topology testsets.
    @test sort!(collect(keys(by_vertex))) == collect(0:11)
    @test length(by_vertex[0]) == 5 && length(by_vertex[11]) == 5
    @test all(v -> length(by_vertex[v]) == 3, 1:10)

    worst = 0.0
    for (_v, ps) in by_vertex, a in ps, b in ps
        worst = max(worst, sd(a, b))
    end
    @test worst < 1e-14
    record!("icosahedron-vertex cluster spread (rad)", worst)
end

@testset "cross-diamond borders line up point for point" begin
    nside = 8
    border = Dict{Int,Vector{GO.UnitSphericalPoint{Float64}}}()
    for d in 0:9
        ps = GO.UnitSphericalPoint{Float64}[]
        for i in 0:nside, j in 0:nside
            (i == 0 || i == nside || j == 0 || j == nside) || continue
            push!(ps, I4.xyd_to_point(i / nside, j / nside, d))
        end
        border[d] = ps
    end
    worst = 0.0
    for d in 0:9, p in border[d]
        best = Inf
        for e in 0:9
            e == d && continue
            for q in border[e]
                best = min(best, sd(p, q))
            end
        end
        worst = max(worst, best)
    end
    @test worst < 1e-13
    record!("cross-diamond border mismatch (rad, nside=8)", worst)
end

@testset "seam ownership: y >= x is the upper half (nside = $nside)" for
        nside in (1, 2, 4, 8, 16)
    # The integer predicate and the Float64 predicate agree on the whole
    # lattice, so the branch a lattice point takes is exactly decidable.
    for ix in 0:nside, iy in 0:nside
        @test ((iy / nside) >= (ix / nside)) == (iy >= ix)
    end
    # ... and everything on or above the seam evaluates through the upper half.
    # `===`, not `≈`: bit-identical shared corners are what makes the
    # tessellation exact.
    for d in 0:9, ix in 0:nside, iy in 0:nside
        p = I4.xyd_to_point(ix / nside, iy / nside, d)
        if iy >= ix
            @test p === upper_half_point(ix / nside, iy / nside, d)
        else
            @test p === lower_half_point(ix / nside, iy / nside, d)
        end
    end
    # Cross-resolution nesting is bit-exact — the reason `node_extent` can be
    # the cell's own cap. See the `node_extent` section.
    if 2nside <= 32
        for ix in 0:nside
            @test (2ix) / (2nside) === ix / nside
        end
        for d in 0:9, ix in 0:nside, iy in 0:nside
            @test I4.xyd_to_point((2ix) / (2nside), (2iy) / (2nside), d) ===
                  I4.xyd_to_point(ix / nside, iy / nside, d)
        end
    end
end

@testset "cell_corners winds CCW from outside (nside = $nside)" for nside in (3, 4, 5)
    worst = Inf
    for d in 0:9, ix in 0:(nside - 1), iy in 0:(nside - 1)
        worst = min(worst, ccw_measure(I4.cell_corners(ix, iy, d, nside)))
    end
    @test worst > 0
end

# =========================================================================
# 2. Oracle: the chart's index maps
# =========================================================================

@testset "row-major and Morton index maps (nside = $nside)" for nside in (1, 2, 3, 4, 5, 8)
    ncell = 10 * nside^2
    for d in 0:9, ix in 0:(nside - 1), iy in 0:(nside - 1)
        q = I4.xyd_to_rowmajor(ix, iy, d, nside)
        @test 0 <= q < ncell
        @test q == d * nside^2 + iy * nside + ix
        @test I4.rowmajor_to_xyd(q, nside) == (ix, iy, d)
        if ispow2(nside)
            p = I4.xyd_to_morton(ix, iy, d, nside)
            @test 0 <= p < ncell
            @test I4.morton_to_xyd(p, nside) == (ix, iy, d)
        end
    end
    @test_throws ArgumentError I4.rowmajor_to_xyd(-1, nside)
    @test_throws ArgumentError I4.rowmajor_to_xyd(ncell, nside)
    if ispow2(nside)
        @test_throws ArgumentError I4.morton_to_xyd(-1, nside)
        @test_throws ArgumentError I4.morton_to_xyd(ncell, nside)
        @test_throws ArgumentError I4.xyd_to_morton(nside, 0, 0, nside)
        @test_throws ArgumentError I4.xyd_to_morton(0, 0, 10, nside)
    else
        @test_throws ArgumentError I4.xyd_to_morton(0, 0, 0, nside)
        @test_throws ArgumentError I4.morton_to_xyd(0, nside)
    end
end

@testset "the Morton code is a plain bit interleave" begin
    # Written out independently: `ix` into the even bit positions, `iy` into the
    # odd ones. Same convention as HEALPix's NESTED id, 12 → 10.
    interleave(ix, iy) = sum(((ix >> b) & 1) << (2b) | ((iy >> b) & 1) << (2b + 1) for b in 0:15)
    for lvl in 0:4
        nside = 2^lvl
        for d in 0:9, ix in 0:(nside - 1), iy in 0:(nside - 1)
            @test I4.xyd_to_morton(ix, iy, d, nside) == d * 4^lvl + interleave(ix, iy)
        end
    end
    for ix in 0:7, iy in 0:7
        @test I4.xyd_to_morton(ix, iy, 0, 8) ÷ 4 == I4.xyd_to_morton(ix ÷ 2, iy ÷ 2, 0, 4)
    end
end

# =========================================================================
# 3. Ids: 0-based ids, 1-based positions
# =========================================================================

@testset "0-based Morton ids, 1-based positions" begin
    @test DGG.cellindextype(SYS) === LevelIndex
    @test cellindextypes(SYS) == (LevelIndex,)
    @test levels(SYS) == 0:29
    @test max_level(SYS) == 29
    @test rootcells(SYS) == [LevelIndex(0, d) for d in 0:9]
    @test max_neighbors(SYS, Vertex()) == 9
    @test max_neighbors(SYS, Edge()) == 4
    @test max_neighbors(SYS) == 9
    @test has_sorted_subtrees(SYS)

    for lvl in 0:4
        g = levelgrid(SYS, lvl)
        @test system(g) === SYS
        @test level(g) == lvl
        @test ncells(g) == 10 * 4^lvl
        @test cellindex(g, 1) == LevelIndex(lvl, 0)
        @test cellindex(g, ncells(g)) == LevelIndex(lvl, ncells(g) - 1)
        @test_throws BoundsError cellindex(g, 0)
        @test_throws BoundsError cellindex(g, ncells(g) + 1)
        for p in cell_sample(lvl)
            @test cellposition(g, LevelIndex(lvl, p)) == p + 1
            @test cellindex(g, p + 1) == LevelIndex(lvl, p)
        end
        # `cellposition` returns `nothing` on any miss and never throws — an
        # id below the range, above it, or at another level entirely.
        @test cellposition(g, LevelIndex(lvl, -1)) === nothing
        @test cellposition(g, LevelIndex(lvl, ncells(g))) === nothing
        @test cellposition(g, LevelIndex(lvl + 1, 0)) === nothing
        lvl > 0 && @test cellposition(g, LevelIndex(lvl - 1, 0)) === nothing
    end

    @test_throws ArgumentError levelgrid(SYS, -1)
    @test_throws ArgumentError levelgrid(SYS, 30)
end

@testset "the id is diamond * 4^level + morton" begin
    for lvl in 0:4
        nside = 1 << lvl
        for d in 0:9, ix in 0:(nside - 1), iy in 0:(nside - 1)
            id = I4.xyd_to_morton(ix, iy, d, nside)
            @test id ÷ (nside^2) == d                      # the diamond is the high part
            @test xyd_of(LevelIndex(lvl, id)) == (ix, iy, d)
        end
    end
    # At level 0 the Morton code is empty, so the id IS the diamond number.
    @test [rawid(c) for c in rootcells(SYS)] == collect(0:9)
end

# =========================================================================
# 4. Hierarchy
# =========================================================================

@testset "hierarchy arithmetic" begin
    for lvl in 0:3
        for p in cell_sample(lvl, 40)
            c = LevelIndex(lvl, p)
            kids = children(SYS, c)
            @test length(kids) == 4
            @test kids == [LevelIndex(lvl + 1, 4p + k) for k in 0:3]
            @test issorted(kids)
            @test all(k -> parent(SYS, k) == c, kids)
            # The lattice reading of the same statement: the four children are
            # the four quadrants of the parent's chart rectangle, in the order
            # (2ix, 2iy), (2ix+1, 2iy), (2ix, 2iy+1), (2ix+1, 2iy+1).
            ix, iy, d = xyd_of(c)
            @test [xyd_of(k) for k in kids] ==
                  [(2ix, 2iy, d), (2ix + 1, 2iy, d), (2ix, 2iy + 1, d), (2ix + 1, 2iy + 1, d)]
            # `ancestor` is the same shift, several levels at once.
            for k in kids, kk in children(SYS, k)
                @test ancestor(SYS, kk, lvl) == c
                @test ancestor(SYS, kk, lvl + 1) == k
                @test ancestor(SYS, kk, level(kk)) == kk
            end
        end
    end
    @test_throws ArgumentError parent(SYS, LevelIndex(0, 0))
    @test_throws ArgumentError children(SYS, LevelIndex(29, 0))
    @test_throws ArgumentError ancestor(SYS, LevelIndex(2, 0), 3)
    @test_throws ArgumentError descendants(SYS, LevelIndex(3, 0), 2)
end

@testset "descendant_range is exact and hole-free" begin
    # `has_sorted_subtrees` is the claim that a subtree is a CONTIGUOUS run of
    # positions, in both directions. This is the empirical verification of the
    # trait rather than a restatement of it: the range is compared against the
    # descendants reached by walking `children`, so a Morton ordering that only
    # LOOKED contiguous would fail here.
    for lvl in 0:2, depth in 1:3
        g = levelgrid(SYS, lvl + depth)
        for p in cell_sample(lvl, 20)
            c = LevelIndex(lvl, p)
            r = descendant_range(SYS, c, lvl + depth)
            @test r isa UnitRange{Int}
            @test length(r) == 4^depth

            # Forward: the walked subtree is exactly the range's positions.
            walked = [c]
            for _ in 1:depth
                walked = reduce(vcat, (children(SYS, x) for x in walked))
            end
            @test sort!([cellposition(g, x) for x in walked]) == collect(r)
            # Backward: every position in the range is a descendant.
            @test all(i -> ancestor(SYS, cellindex(g, i), lvl) == c, r)
            # `descendants` agrees with the range and is ascending.
            @test descendants(SYS, c, lvl + depth) == [cellindex(g, i) for i in r]
        end
        @test descendant_range(SYS, LevelIndex(lvl, 0), lvl) == 1:1
    end
    # Sibling ranges tile the parent's, in order.
    for lvl in 1:3, p in cell_sample(lvl - 1, 20)
        c = LevelIndex(lvl - 1, p)
        rs = [descendant_range(SYS, k, lvl) for k in children(SYS, c)]
        @test reduce(vcat, [collect(r) for r in rs]) == collect(descendant_range(SYS, c, lvl))
    end
    @test_throws ArgumentError descendant_range(SYS, LevelIndex(3, 0), 2)
    @test_throws ArgumentError descendant_range(SYS, LevelIndex(3, 0), 30)
end

# =========================================================================
# 5. node_extent — the exact subtree cap
# =========================================================================

@testset "the four-corner cap already bounds the whole cell" begin
    # THE measurement `node_extent` rests on: the farthest point of a cell from
    # its own centre is one of its four corners, seam-straddling cells included,
    # so a cap through the corners needs no inflation for the cell itself. A
    # dense 17x17 sampling of every cell's chart rectangle, on every diamond, at
    # every `nside` below. The pre-redesign face-grid layer recorded a worst
    # overhang of exactly 0.0; this re-runs it rather than trusting the record.
    worst = -Inf
    for nside in (1, 2, 3, 4, 5, 8, 16)
        for d in 0:9, ix in 0:(nside - 1), iy in 0:(nside - 1)
            ctr = I4.cell_center(ix, iy, d, nside)
            rc = maximum(sd(ctr, p) for p in I4.cell_corners(ix, iy, d, nside))
            for i in 0:16, j in 0:16
                worst = max(worst, sd(ctr, I4.xyd_to_point((ix + i / 16) / nside,
                                                           (iy + j / 16) / nside, d)) - rc)
            end
        end
    end
    @test worst <= 0.0
    record!("interior overhang past the four-corner cap (rad)", worst)
end

@testset "node_extent covers the true cell, with margin" begin
    # The conformance suite samples DESCENDANT VERTICES; this samples the cell's
    # own boundary at 32x the shipped densification, which is the part a
    # vertex-only check can miss on a curved edge.
    for lvl in 0:4
        worst = -Inf
        for p in cell_sample(lvl, 128)
            c = LevelIndex(lvl, p)
            cap = node_extent(SYS, c)
            ix, iy, d = xyd_of(c)
            for q in I4._perimeter_points(ix, iy, d, 1 << lvl, 32 * I4.CAP_EDGE_SEGMENTS)
                worst = max(worst, sd(cap.point, q) - cap.radius)
            end
        end
        @test worst < 0                                   # strictly inside
        record!("node_extent slack over the dense perimeter (rad)", worst)
    end
end

@testset "node_extent covers deep descendants" begin
    for lvl in 0:2, depth in 1:3
        gd = levelgrid(SYS, lvl + depth)
        for p in cell_sample(lvl, 20)
            c = LevelIndex(lvl, p)
            cap = node_extent(SYS, c)
            for q in descendants(SYS, c, lvl + depth), v in cell_boundary(gd, q)
                @test sd(cap.point, v) <= cap.radius
            end
        end
    end
end

@testset "node_extent is geodesically convex, and beats the inflated default" begin
    for lvl in 0:4
        rmax = 0.0
        for p in cell_sample(lvl, 128)
            cap = node_extent(SYS, LevelIndex(lvl, p))
            rmax = max(rmax, cap.radius)
        end
        @test rmax < pi / 2                    # what makes vertex sampling sound
        record!("node_extent radius at level $lvl (deg)", rad2deg(rmax))
    end
    # `cap_inflation` is never consulted; it keeps its default and the override
    # is tighter than what it would have produced.
    @test cap_inflation(SYS) == 1.2
    for lvl in 1:3, p in cell_sample(lvl, 32)
        c = LevelIndex(lvl, p)
        g = levelgrid(SYS, lvl)
        ctr = cell_centroid(g, c)
        inflated = 1.2 * maximum(sd(ctr, v) for v in cell_boundary(g, c))
        @test node_extent(SYS, c).radius < inflated
    end
end

# =========================================================================
# 6. Topology: the seam tables
# =========================================================================

@testset "the seam tables are an involution over the twenty rim edges" begin
    # Restated independently of `topology.jl`'s own build asserts: a diamond's
    # four rim edges are icosahedron edges, the ten seams are the other ten, and
    # forty rim slots pair up into twenty joins.
    pairs = Set{Tuple{Int,Int}}()
    for d in 0:9, e in 1:4
        a, b = I4._edge_pair(d, e)
        @test abs(ISEA.vdot(ISEA.VERTICES[a + 1], ISEA.VERTICES[b + 1]) - ISEA.ADJ_DOT) < 1e-9
        push!(pairs, minmax(a, b))
        d2, e2 = I4.EDGE_NEIGHBORS[d + 1][e]
        @test I4._edge_pair(d2, e2) == (b, a)             # reversed, always
        @test I4.EDGE_NEIGHBORS[d2 + 1][e2] == (d, e)     # involutive
        @test (d2, e2) != (d, e)
    end
    @test length(pairs) == 20                             # twenty rim edges
    seams = Set(minmax(dm.verts[1], dm.verts[3]) for dm in I4.DIAMONDS)
    @test length(seams) == 10                             # ten seams
    @test isempty(intersect(pairs, seams))                # and they are disjoint
    @test length(union(pairs, seams)) == 30               # the icosahedron's edges
end

@testset "the vertex fans close, and the valences are 5, 5 and ten times 3" begin
    counts = zeros(Int, 12)
    for d in 0:9, c in 1:4
        v = I4.DIAMONDS[d + 1].verts[c]
        counts[v + 1] += 1
        fan = I4.CORNER_FANS[d + 1][c]
        # Every fan member is a diamond-corner at the SAME icosahedron vertex,
        # they are distinct, and they do not include the subject.
        @test all(((d2, c2),) -> I4.DIAMONDS[d2 + 1].verts[c2] == v, fan)
        @test allunique(fan)
        @test !((d, c) in fan)
        @test length(fan) == (v in (0, 11) ? 4 : 2)
        # The ends of the fan are the two corners the AXIS offsets reach.
        @test fan[1][1] == I4.EDGE_NEIGHBORS[d + 1][mod1(c - 1, 4)][1]
        @test fan[end][1] == I4.EDGE_NEIGHBORS[d + 1][c][1]
    end
    @test counts == [5, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 5]
end

@testset "neighbour counts over whole levels" begin
    # Exhaustive, every cell of levels 0-3. The distribution IS the topology:
    # 9 at the twenty corner cells on vertices 0 and 11, 7 at the thirty on the
    # valence-3 vertices, 8 everywhere else. Level 0 is its own case — a diamond
    # is one cell, so every cell is all four corners at once and the base tiling
    # comes out 6-regular.
    expected = Dict(
        0 => Dict(6 => 10),
        1 => Dict(7 => 30, 9 => 10),
        2 => Dict(7 => 30, 8 => 120, 9 => 10),
        3 => Dict(7 => 30, 8 => 600, 9 => 10),
    )
    for lvl in 0:3
        g = levelgrid(SYS, lvl)
        counts = Dict{Int,Int}()
        for p in 0:(ncells(g) - 1)
            k = length(neighbors(g, LevelIndex(lvl, p), 1))
            counts[k] = get(counts, k, 0) + 1
        end
        @test counts == expected[lvl]
        # `Edge()` is uniformly four: every cell has four rim segments, and a
        # rim segment always has exactly one cell on the other side.
        @test all(p -> length(neighbors(g, LevelIndex(lvl, p), 1; connectivity = Edge())) == 4,
                  0:(ncells(g) - 1))
    end
end

@testset "neighbours are symmetric and Edge() is Vertex() restricted" begin
    # Exhaustive over levels 0-3, both directions, both connectivities.
    for lvl in 0:3, conn in (Vertex(), Edge())
        g = levelgrid(SYS, lvl)
        adj = [Set(neighbors(g, LevelIndex(lvl, p), 1; connectivity = conn))
               for p in 0:(ncells(g) - 1)]
        for p in 0:(ncells(g) - 1)
            c = LevelIndex(lvl, p)
            @test !(c in adj[p + 1])
            @test length(adj[p + 1]) == length(collect(neighbors(g, c, 1; connectivity = conn)))
            for q in adj[p + 1]
                @test level(q) == lvl
                @test cellposition(g, q) !== nothing
                @test c in adj[rawid(q) + 1]
            end
        end
    end
    for lvl in 0:3
        g = levelgrid(SYS, lvl)
        for p in 0:(ncells(g) - 1)
            c = LevelIndex(lvl, p)
            vs = collect(neighbors(g, c, 1))
            es = collect(neighbors(g, c, 1; connectivity = Edge()))
            @test es == filter(in(Set(es)), vs)   # a subset, order preserved
        end
    end
end

@testset "Edge() neighbours share a whole edge, corner-only ones share a point" begin
    # THE test that decides which four of the eight `Edge()` keeps, and the one
    # every other property in this section is satisfied by just as happily under
    # the WRONG four. An ISEA4R cell is an axis-aligned square in its chart, so
    # the axis offsets are the edge-sharing ones and the diagonals are the
    # corner-only ones — the OPPOSITE of HEALPix, whose pixel is rotated 45°
    # against its lattice. Only geometry can tell the two apart.
    #
    # Exhaustive over levels 1-3, so the rim and vertex cells are all in it: a
    # tolerance-based count (see `shared_boundary_points`) rather than a set
    # intersection, because bit-identity holds only within a diamond.
    for lvl in 1:3
        g = levelgrid(SYS, lvl)
        for p in 0:(ncells(g) - 1)
            c = LevelIndex(lvl, p)
            vs = collect(neighbors(g, c, 1))
            es = collect(neighbors(g, c, 1; connectivity = Edge()))
            for n in es
                # A whole densified edge: `BOUNDARY_SEGMENTS` start points plus
                # the far corner, which the other ring supplies.
                @test shared_boundary_points(g, c, n) == I4.BOUNDARY_SEGMENTS + 1
            end
            for n in setdiff(vs, es)
                @test shared_boundary_points(g, c, n) == 1     # a single corner
            end
        end
    end
end

@testset "the vertex cells: nine neighbours, and every one of them at the vertex" begin
    # The degenerate case spelled out. At icosahedron vertices 0 and 11 five
    # diamonds meet, so the corner cell there has FOUR cells sharing the vertex
    # with it — two reached by axis offsets and two only by the diagonal — and
    # nine neighbours in total.
    for lvl in 1:4
        g = levelgrid(SYS, lvl)
        nside = 1 << lvl
        found = 0
        for d in 0:9, c in 1:4
            v = I4.DIAMONDS[d + 1].verts[c]
            v in (0, 11) || continue
            found += 1
            ix, iy = I4._corner_cell(c, nside)
            cell = LevelIndex(lvl, I4.xyd_to_morton(ix, iy, d, nside))
            ns = collect(neighbors(g, cell, 1))
            @test length(ns) == 9
            # Exactly four of the nine touch the vertex — the other five touch
            # only an edge or another corner of the cell.
            vp = vpoint(v)
            atvertex = [n for n in ns
                        if minimum(sd(vp, q) for q in cell_boundary(g, n)) < 1e-12]
            @test length(atvertex) == 4
            # ... and those four are corner cells of the OTHER four diamonds.
            @test length(unique(last.(xyd_of.(atvertex)))) == 4
            @test !(d in last.(xyd_of.(atvertex)))
        end
        @test found == 10       # five northern v01 corners, five southern v10
    end
end

# =========================================================================
# 7. Rotational order
# =========================================================================

@testset "neighbour order is CCW seen from outside" begin
    # The rotational-order contract, measured rather than argued: sum the signed
    # azimuth steps of the neighbour centres around the cell centre. A closed
    # loop winding counter-clockwise seen from outside totals +2π; a clockwise
    # one totals -2π. `e1 × e2 == centre` is what makes the frame right-handed
    # SEEN FROM OUTSIDE; building it the other way would flip every sign below
    # and the test would pass happily on a clockwise implementation.
    function basis(c)
        a = abs(c[3]) < 0.9 ? (0.0, 0.0, 1.0) : (1.0, 0.0, 0.0)
        e1 = (a[2] * c[3] - a[3] * c[2], a[3] * c[1] - a[1] * c[3], a[1] * c[2] - a[2] * c[1])
        n = sqrt(sum(e1 .^ 2)); e1 = e1 ./ n
        e2 = (c[2] * e1[3] - c[3] * e1[2], c[3] * e1[1] - c[1] * e1[3], c[1] * e1[2] - c[2] * e1[1])
        return e1, e2
    end
    azim(c, e1, e2, p) = atan(sum((p .- c) .* e2), sum((p .- c) .* e1))
    wrap(d) = (d = mod(d, 2pi); d > pi ? d - 2pi : d)
    function winding(grid, c, cells)
        ctr = Tuple(cell_centroid(grid, c))
        e1, e2 = basis(ctr)
        as = [azim(ctr, e1, e2, Tuple(cell_centroid(grid, x))) for x in cells]
        return sum(wrap(as[mod1(i + 1, length(as))] - as[i]) for i in eachindex(as))
    end

    # Exhaustive over levels 0-3, so every rim, corner and vertex cell is wound.
    for lvl in 0:3, conn in (Vertex(), Edge())
        g = levelgrid(SYS, lvl)
        for p in 0:(ncells(g) - 1)
            c = LevelIndex(lvl, p)
            ns = collect(neighbors(g, c, 1; connectivity = conn))
            length(ns) < 3 && continue
            @test winding(g, c, ns) ≈ 2pi atol = 1e-6
        end
    end
    # Outer rings are wound geometrically; they must come out CCW too.
    for lvl in (3, 4), k in 2:3
        g = levelgrid(SYS, lvl)
        for p in cell_sample(lvl, 24)
            shell = collect(ring(g, LevelIndex(lvl, p), k))
            length(shell) < 3 && continue
            @test winding(g, LevelIndex(lvl, p), shell) ≈ 2pi atol = 1e-6
        end
    end
end

@testset "ring 1 follows the documented chart offset cycle" begin
    # In a diamond's interior the whole cycle is readable off the lattice, so
    # check it there against the offsets written out independently.
    offsets = ((1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1))
    for lvl in 2:4
        nside = 1 << lvl
        g = levelgrid(SYS, lvl)
        checked = 0
        for p in cell_sample(lvl, 96)
            ix, iy, d = I4.morton_to_xyd(p, nside)
            (1 <= ix <= nside - 2 && 1 <= iy <= nside - 2) || continue
            c = LevelIndex(lvl, p)
            expected = [LevelIndex(lvl, I4.xyd_to_morton(ix + dx, iy + dy, d, nside))
                        for (dx, dy) in offsets]
            @test collect(neighbors(g, c, 1)) == expected
            @test collect(neighbors(g, c, 1; connectivity = Edge())) == expected[[1, 3, 5, 7]]
            checked += 1
        end
        @test checked > 0        # the interior filter did not empty the loop
    end
end

@testset "ring 1 starts where the docstring says" begin
    # ORACLE PIN on the documented START of the rotational order.
    #
    # The neighbouring testsets check the CYCLE (CCW winding, the offset tuple,
    # Edge() preserving the order), and the conformance harness's winding law is
    # start-invariant by construction — every rotation of a CCW cycle is still a
    # CCW cycle. So the START itself is documented and otherwise unpinned, and a
    # rotation would pass everything while silently rotating the weight vector
    # of any stencil that has baked "slot 1 is the +x chart direction" into
    # itself.
    #
    # Three cells at level 3, one of each kind the topology has. Values produced
    # by the implementation and checked against the documented rule (the chart
    # offset cycle from `(+1, 0)`), then frozen as literals.
    g = levelgrid(SYS, 3)

    # (a) a diamond-interior cell: the plain eight.
    interior = LevelIndex(3, I4.xyd_to_morton(3, 5, 0, 8))
    @test rawid(interior) == 39
    @test [rawid(x) for x in neighbors(g, interior, 1)] ==
          [50, 56, 45, 44, 38, 36, 37, 48]
    @test [rawid(x) for x in neighbors(g, interior, 1; connectivity = Edge())] ==
          [50, 45, 38, 37]

    # (b) a corner cell on a VALENCE-3 icosahedron vertex: seven, the NE
    #     diagonal dropping out of the cycle with no hole left behind.
    corner3 = LevelIndex(3, I4.xyd_to_morton(7, 7, 0, 8))
    @test rawid(corner3) == 63
    @test [rawid(x) for x in neighbors(g, corner3, 1)] ==
          [426, 64, 66, 62, 60, 61, 424]
    @test [rawid(x) for x in neighbors(g, corner3, 1; connectivity = Edge())] ==
          [426, 64, 62, 61]

    # (c) a corner cell on icosahedron VERTEX 0: nine, the NW diagonal yielding
    #     the two fan cells (ids 170 and 234, on diamonds 2 and 3) between its
    #     two axis neighbours.
    vertex0 = LevelIndex(3, I4.xyd_to_morton(0, 7, 0, 8))
    @test rawid(vertex0) == 42
    @test [rawid(x) for x in neighbors(g, vertex0, 1)] ==
          [43, 104, 106, 170, 234, 298, 299, 40, 41]
    @test [rawid(x) for x in neighbors(g, vertex0, 1; connectivity = Edge())] ==
          [43, 106, 298, 40]

    # (d) level 0, where a diamond is one cell and the base tiling is 6-regular.
    @test [rawid(x) for x in neighbors(levelgrid(SYS, 0), LevelIndex(0, 0), 1)] ==
          [6, 1, 2, 3, 4, 5]
    @test [rawid(x) for x in neighbors(levelgrid(SYS, 0), LevelIndex(0, 0), 1;
                                       connectivity = Edge())] == [6, 1, 4, 5]
end

@testset "ring is the tail block of neighbors" begin
    # The composition contract: neighbors(k) is rings 1..k concatenated outward,
    # so ring(k) is exactly the trailing block.
    for lvl in (2, 3, 4), conn in (Vertex(), Edge())
        g = levelgrid(SYS, lvl)
        for p in cell_sample(lvl, 24)
            c = LevelIndex(lvl, p)
            acc = LevelIndex[]
            seen = Set{LevelIndex}()
            for k in 1:3
                shell = collect(ring(g, c, k; connectivity = conn))
                @test isempty(intersect(Set(shell), seen))
                union!(seen, shell)
                append!(acc, shell)
                disc = collect(neighbors(g, c, k; connectivity = conn))
                @test disc == acc                                    # concatenated outward
                @test disc[(end - length(shell) + 1):end] == shell    # tail block
                @test allunique(disc)
                @test !(c in disc)
                @test Set(disc) == seen
            end
        end
    end
    for lvl in (2, 3)
        g = levelgrid(SYS, lvl)
        @test ring(g, LevelIndex(lvl, 0), 0) == [LevelIndex(lvl, 0)]
        @test isempty(neighbors(g, LevelIndex(lvl, 0), 0))
        @test_throws ArgumentError neighbors(g, LevelIndex(lvl, 0), -1)
        @test_throws ArgumentError ring(g, LevelIndex(lvl, 0), -1)
    end
end

# =========================================================================
# 8. The subtree rim
# =========================================================================

@testset "subtree_border is the lattice block's boundary ring" begin
    # The automaton against the definition, brute-forced independently: a
    # level-`l` descendant is on the rim iff one of its neighbours is not a
    # descendant. Exhaustive over the level-0 and level-1 subtrees at depths
    # 1-3, so blocks against a diamond rim and blocks at icosahedron vertices
    # are both in it.
    function brute_border(c, l, conn)
        g = levelgrid(SYS, l)
        kids = descendants(SYS, c, l)
        inside = Set(kids)
        return [k for k in kids
                if any(n -> !(n in inside), neighbors(g, k, 1; connectivity = conn))]
    end
    for lvl in 0:1, depth in 1:3
        for p in cell_sample(lvl, 12)
            c = LevelIndex(lvl, p)
            l = lvl + depth
            for conn in (Vertex(), Edge())
                border = subtree_border(SYS, c, l; connectivity = conn)
                @test border == brute_border(c, l, conn)      # order included
                @test issorted(border)
                @test eltype(border) === LevelIndex
                @test length(border) == 4 * (1 << depth) - 4
            end
            # `Vertex()` and `Edge()` give the same rim; the docstring says so.
            @test subtree_border(SYS, c, l) ==
                  subtree_border(SYS, c, l; connectivity = Edge())
            # Border and interior partition the subtree.
            interior = subtree_interior(SYS, c, l)
            kids = descendants(SYS, c, l)
            @test issorted(interior)
            @test isempty(intersect(Set(interior), Set(subtree_border(SYS, c, l))))
            @test sort!(vcat(collect(interior), collect(subtree_border(SYS, c, l)))) == kids
        end
    end
    # A depth-0 subtree is the cell itself, and it is all rim.
    @test subtree_border(SYS, LevelIndex(2, 7), 2) == [LevelIndex(2, 7)]
    @test isempty(subtree_interior(SYS, LevelIndex(2, 7), 2))
    @test_throws ArgumentError subtree_border(SYS, LevelIndex(2, 7), 1)
    @test_throws ArgumentError subtree_border(SYS, LevelIndex(2, 7), 30)
end

@testset "the generic substrate reaches ISEA4R" begin
    # A smoke test of the fallback layer over this system, so that the
    # `has_sorted_subtrees` fast paths — the subtree grid backed by
    # `descendant_range`, and the cursor the tree descends — are exercised
    # somewhere. The substrate has its own suite; what is checked here is that
    # this system plugs into it and that the answers are the same ones a
    # brute-force sweep of the level gives.
    g = levelgrid(SYS, 3)
    @test treeify(g) !== nothing

    # A cap query. The tree only ever prunes, so every cell the brute-force
    # sweep finds must be in the answer.
    target = US.SphericalCap(cell_centroid(g, LevelIndex(3, 100)), 0.15)
    hits = query(g, Intersects(target))
    @test eltype(hits) === LevelIndex
    @test issorted(hits)
    @test allunique(hits)
    @test LevelIndex(3, 100) in hits
    brute = [LevelIndex(3, p) for p in 0:(ncells(g) - 1)
             if sd(target.point, cell_centroid(g, LevelIndex(3, p))) <= target.radius ||
                any(v -> sd(target.point, v) <= target.radius,
                    cell_boundary(g, LevelIndex(3, p)))]
    @test issubset(Set(brute), Set(hits))
    @test length(hits) < ncells(g)                 # and it really did prune

    # The subtree grid: `PartialGrid(sys, cell, level)` takes the
    # `descendant_range` path when `has_sorted_subtrees` is true, which is the
    # trait's payoff.
    sub = PartialGrid(SYS, LevelIndex(0, 0), 3)
    @test ncells(sub) == 4^3
    @test system(sub) === SYS
    @test level(sub) == 3
    @test all(i -> cellposition(sub, cellindex(sub, i)) == i, 1:ncells(sub))
    @test Set(cellindex(sub, i) for i in 1:ncells(sub)) ==
          Set(descendants(SYS, LevelIndex(0, 0), 3))
    # A cell of another diamond is simply not in it — `nothing`, not an error.
    @test cellposition(sub, LevelIndex(3, 700)) === nothing
    # Coverage: a neighbour outside the subtree is ABSENT, never padded.
    corner = LevelIndex(3, I4.xyd_to_morton(0, 0, 0, 8))
    @test length(neighbors(sub, corner, 1)) < length(neighbors(g, corner, 1))
    @test all(n -> cellposition(sub, n) !== nothing, neighbors(sub, corner, 1))
    @test Set(neighbors(sub, corner, 1)) ==
          Set(n for n in neighbors(g, corner, 1) if cellposition(sub, n) !== nothing)
end

# =========================================================================
# 9. Geometry
# =========================================================================

@testset "cell_area is exactly equal-area" begin
    for lvl in 0:4
        g = levelgrid(SYS, lvl)
        exact = 4pi / (10 * 4^lvl)
        for p in cell_sample(lvl, 64)
            @test cell_area(g, LevelIndex(lvl, p)) == exact
        end
        # The whole level sums to the sphere.
        @test sum(cell_area(g, LevelIndex(lvl, p)) for p in 0:(ncells(g) - 1)) ≈ 4pi
    end
end

@testset "the densified boundary polygon converges to the exact area" begin
    # The `BOUNDARY_SEGMENTS` schedule, measured: the relative error is flat
    # across levels and falls as segments^-2, which is why the shipped count
    # does not depend on the level. LEVEL 0 IS EXACT and is not evidence about
    # any other level — a diamond's rim edges are icosahedron edges, hence
    # great circles, so a level-0 cell IS its 4-gon.
    for nseg in (2, 4, 8, 16)
        worst = 0.0
        for lvl in 0:3
            nside = 1 << lvl
            exact = 4pi / (10 * 4.0^lvl)
            for d in 0:9, ix in 0:(nside - 1), iy in 0:(nside - 1)
                a = ring_area(I4._perimeter_points(ix, iy, d, nside, nseg))
                worst = max(worst, abs(a - exact) / exact)
            end
        end
        record!("densified ring rel. area error at nseg=$nseg", worst)
        @test worst < 3.0e-2 / (nseg / 2)^2 * 1.2      # O(nseg^-2), with headroom
    end
    # Level 0 is exact to round-off.
    for d in 0:9
        @test ring_area(I4._perimeter_points(0, 0, d, 1, 8)) ≈ 4pi / 10 rtol = 1e-13
    end
    # The shipped ring is within 0.2% of the closed form everywhere.
    for lvl in 0:3
        g = levelgrid(SYS, lvl)
        for p in cell_sample(lvl, 64)
            c = LevelIndex(lvl, p)
            @test ring_area(cell_boundary(g, c)) ≈ cell_area(g, c) rtol = 2e-3
        end
    end
end

@testset "cell_boundary is a closed CCW ring of unit points" begin
    for lvl in 0:3
        g = levelgrid(SYS, lvl)
        for p in cell_sample(lvl, 48)
            c = LevelIndex(lvl, p)
            ring_c = cell_boundary(g, c)
            @test length(ring_c) == 4 * I4.BOUNDARY_SEGMENTS
            @test allunique(ring_c)                      # implicitly closed
            @test all(q -> abs(sqrt(sum(abs2, q)) - 1) < 1e-12, ring_c)
            @test ccw_measure(ring_c) > 0
            # The four corners are vertices 1, 9, 17, 25 of the 32.
            ix, iy, d = xyd_of(c)
            @test ring_c[1:I4.BOUNDARY_SEGMENTS:end] ==
                  collect(I4.cell_corners(ix, iy, d, 1 << lvl))
        end
    end
    @test_throws ArgumentError cell_boundary(levelgrid(SYS, 2), LevelIndex(3, 0))
    @test_throws ArgumentError cell_boundary(levelgrid(SYS, 2), LevelIndex(2, 160))
    @test_throws ArgumentError cell_centroid(levelgrid(SYS, 2), LevelIndex(2, -1))
end

@testset "shared edges are bit-identical within a diamond" begin
    # What makes the tessellation exact rather than merely consistent to
    # rounding — and the reason the densification count is a power of two. Only
    # WITHIN a diamond: across a rim the two sides are two developments of one
    # icosahedron edge, which `cross-diamond borders line up point for point`
    # measures instead.
    for lvl in 2:3
        nside = 1 << lvl
        g = levelgrid(SYS, lvl)
        for d in 0:9, ix in 0:(nside - 2), iy in 0:(nside - 2)
            c = LevelIndex(lvl, I4.xyd_to_morton(ix, iy, d, nside))
            east = LevelIndex(lvl, I4.xyd_to_morton(ix + 1, iy, d, nside))
            north = LevelIndex(lvl, I4.xyd_to_morton(ix, iy + 1, d, nside))
            # Hashed containers are legitimate HERE and nowhere near a seam;
            # `vkey` guards the signed-zero trap regardless (see its docstring).
            rc = Set(vkey.(cell_boundary(g, c)))
            @test length(intersect(rc, Set(vkey.(cell_boundary(g, east))))) ==
                  I4.BOUNDARY_SEGMENTS + 1
            @test length(intersect(rc, Set(vkey.(cell_boundary(g, north))))) ==
                  I4.BOUNDARY_SEGMENTS + 1
        end
    end
end

@testset "signed zeros cannot break a hashed incidence test" begin
    # The trap `vkey` exists for, pinned as two facts about this system.
    #
    # (1) No boundary point anywhere carries a NEGATIVE zero coordinate, so a
    #     hashed comparison of two bit-identical points can never split on the
    #     sign of a zero. Exactly-zero coordinates do occur (32 of them over
    #     levels 0-3, at the chart corners that land on an icosahedron vertex
    #     with a zero component), which is why the fact is worth pinning rather
    #     than assuming.
    #
    # (2) Across a diamond rim there is NO bit-identity to lose in the first
    #     place: the two sides are independent developments of one icosahedron
    #     edge and agree to ~4e-15 rad but on no coordinate exactly. That is why
    #     every cross-seam incidence test above is written with a tolerance, and
    #     the reason a future "optimisation" of one of them into a `Set`
    #     intersection would report seam neighbours as sharing nothing.
    negzeros = 0
    exactzeros = 0
    for lvl in 0:3
        g = levelgrid(SYS, lvl)
        for p in 0:(ncells(g) - 1), q in cell_boundary(g, LevelIndex(lvl, p)), k in 1:3
            q[k] == 0.0 || continue
            exactzeros += 1
            signbit(q[k]) && (negzeros += 1)
        end
    end
    @test exactzeros > 0            # the case really does arise
    @test negzeros == 0             # ... and never with the sign bit set

    coincident = 0
    bitequal = 0
    lvl = 2
    g = levelgrid(SYS, lvl)
    for p in 0:(ncells(g) - 1)
        c = LevelIndex(lvl, p)
        _, _, d = xyd_of(c)
        for n in neighbors(g, c, 1; connectivity = Edge())
            last(xyd_of(n)) == d && continue          # only cross-diamond joins
            for a in cell_boundary(g, c), b in cell_boundary(g, n)
                sd(a, b) <= 1e-12 || continue
                coincident += 1
                a == b && (bitequal += 1)
            end
        end
    end
    @test coincident > 0
    @test bitequal == 0
end

@testset "cell_centroid is the chart midpoint, strictly inside" begin
    for lvl in 0:3
        g = levelgrid(SYS, lvl)
        nside = 1 << lvl
        for p in cell_sample(lvl, 48)
            c = LevelIndex(lvl, p)
            ix, iy, d = xyd_of(c)
            @test cell_centroid(g, c) === I4.cell_center(ix, iy, d, nside)
            @test cell_centroid(g, c) ===
                  I4.xyd_to_point((ix + 0.5) / nside, (iy + 0.5) / nside, d)
            # Strictly interior: closer to the centre than any boundary point.
            ctr = cell_centroid(g, c)
            @test minimum(sd(ctr, q) for q in cell_boundary(g, c)) > 0
        end
    end
end

# =========================================================================
# 10. Location
# =========================================================================

@testset "cellat inverts cell_centroid over whole levels" begin
    for lvl in 0:4
        g = levelgrid(SYS, lvl)
        @test all(0:(ncells(g) - 1)) do p
            cellat(g, cell_centroid(g, LevelIndex(lvl, p))) == LevelIndex(lvl, p)
        end
    end
end

@testset "cellat is the chart inverse, and round-trips the chart" begin
    for lvl in 0:4
        nside = 1 << lvl
        g = levelgrid(SYS, lvl)
        for p in cell_sample(lvl, 64)
            ix, iy, d = I4.morton_to_xyd(p, nside)
            # An interior chart point of the cell lands in the cell.
            for (fx, fy) in ((0.25, 0.25), (0.5, 0.75), (0.75, 0.5), (0.1, 0.9))
                q = I4.xyd_to_point((ix + fx) / nside, (iy + fy) / nside, d)
                @test cellat(g, q) == LevelIndex(lvl, p)
            end
        end
    end
    # The continuous inverse round-trips the chart itself.
    worst = 0.0
    for d in 0:9, i in 0:8, j in 0:8
        x = (i + 0.5) / 9; y = (j + 0.5) / 9
        p = I4.xyd_to_point(x, y, d)
        xx, yy, dd = I4.point_to_xy(p)
        @test dd == d
        worst = max(worst, max(abs(xx - x), abs(yy - y)))
    end
    @test worst < 1e-12
    record!("point_to_xy chart round-trip (chart units)", worst)
    # ... and inverts the sphere direction too.
    worst2 = 0.0
    for d in 0:9, i in 0:8, j in 0:8
        x = (i + 0.5) / 9; y = (j + 0.5) / 9
        p = I4.xyd_to_point(x, y, d)
        worst2 = max(worst2, sd(p, I4.xyd_to_point(I4.point_to_xy(p)...)))
    end
    @test worst2 < 1e-13
    record!("cellat/chart sphere round-trip (rad)", worst2)
end

@testset "cellat on cell boundaries: deterministic, incident, self-consistent" begin
    # A point exactly on a shared boundary belongs to ONE cell; the contract is
    # that the answer is one of the cells genuinely incident to the point, the
    # same on every call, and self-consistent.
    lvl = 3
    nside = 1 << lvl
    g = levelgrid(SYS, lvl)
    for d in 0:9, ix in 0:2:(nside - 1), iy in 0:2:(nside - 1)
        for (x, y) in ((ix / nside, iy / nside),                 # a lattice corner
                       ((ix + 0.5) / nside, iy / nside),         # an edge midpoint
                       (ix / nside, (iy + 0.5) / nside))
            q = I4.xyd_to_point(x, y, d)
            c = cellat(g, q)
            @test cellat(g, q) == c                              # deterministic
            @test cellposition(g, c) !== nothing                 # a real cell
            # Incident: the point is on that cell's own boundary (or inside it).
            @test minimum(sd(q, v) for v in cell_boundary(g, c)) < 1e-9
        end
    end

    # THE TWELVE ICOSAHEDRON VERTICES, explicitly. The even-lattice sweep above
    # steps `0:2:nside-1`, so it probes chart corner (0,0) but never (0,1) or
    # (1,0) — and those are exactly the corners that carry vertices 0 and 11,
    # the two degenerate five-diamond points. A vertex is the hardest tie in the
    # system: three or five cells are incident and `snyder_fwd`'s face choice is
    # a five-way argmax on equal dot products. Taken from `ISEA.VERTICES`
    # directly rather than through the chart, so this probes the point itself
    # and not a chart evaluation of it.
    for v in 0:11
        q = vpoint(v)
        c = cellat(g, q)
        @test cellat(g, q) == c                                  # deterministic
        @test cellposition(g, c) !== nothing
        @test minimum(sd(q, w) for w in cell_boundary(g, c)) < 1e-9
        # The cell named really is one of the cells at that vertex: its chart
        # corner slot on its own diamond carries this vertex.
        ix, iy, d = xyd_of(c)
        slot = findfirst(==(v), I4.DIAMONDS[d + 1].verts)
        @test slot !== nothing
        @test (ix, iy) == I4._corner_cell(slot, nside)
    end

    # Random points anywhere: `cellat` never returns `nothing` on a complete
    # level grid, and the cell it names contains the point's own cell centre
    # neighbourhood — checked by the round trip through the chart.
    rng = MersenneTwister(4120)
    for _ in 1:2000
        v = randn(rng, 3)
        q = GO.UnitSphericalPoint((v ./ sqrt(sum(abs2, v)))...)
        c = cellat(g, q)
        @test c isa LevelIndex
        @test cellposition(g, c) !== nothing
        ix, iy, d = xyd_of(c)
        x, y, dd = I4.point_to_xy(q)
        @test dd == d
        @test ix <= x * nside <= ix + 1
        @test iy <= y * nside <= iy + 1
    end
end

@testset "the ten diamonds partition the sphere" begin
    # Every random point lands in exactly one cell, and the level's polygons sum
    # to the whole sphere.
    for lvl in 0:2
        g = levelgrid(SYS, lvl)
        total = sum(ring_area(cell_boundary(g, LevelIndex(lvl, p))) for p in 0:(ncells(g) - 1))
        @test total ≈ 4pi rtol = 2e-3
    end
end

# =========================================================================
# 11. Contract: the conformance suites, default kwargs
# =========================================================================

@testset "conformance" begin
    for l in (0, 1, 3)
        test_grid_interface(levelgrid(SYS, l); label = "ISEA4RGrid(level=$l)")
    end
    test_hierarchical_system(SYS)
end

end # @testset "ISEA4R system (T12)"

@printf("[ISEA4R] measured:\n")
for key in sort!(collect(keys(MEASURED)))
    @printf("[ISEA4R]   %-52s %+.3e\n", key, MEASURED[key])
end

end # module ISEA4RSystemTests
