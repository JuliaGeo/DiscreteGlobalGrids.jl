module Isea4rDiamondTestSuite

# Tests for `src/ISEA4R/diamonds.jl` and `src/ISEA4R/chart.jl`: the ten-diamond
# layout table and the piecewise-affine rhombus chart built on it.
#
# Three things are being pinned down here.
#
# 1. *The layout is what it says it is.* `diamonds.jl` derives the table from
#    `ISEA`'s icosahedron tables and asserts it against literals it carries
#    itself; that catches upstream drift but not a coordinated edit of both. So
#    the same numbers are restated independently below, and the structural
#    facts (the ten pairs partition the twenty faces, every seam is a real
#    icosahedron edge, northern/southern apexes are 0/11) are re-derived rather
#    than restated.
#
# 2. *The chart is one continuous, exactly equal-area map of the square.* The
#    two affine halves agree on the seam, the four chart corners land on their
#    `ISEA.VERTICES` entries, evaluations of the same icosahedron vertex from
#    different diamonds cluster to round-off, and a densified cell boundary
#    encloses exactly `4π/(10 nside²)`.
#
# 3. *Seam ownership is exactly decidable.* Every lattice point with `iy >= ix`
#    — the seam included — goes through the upper half, `===` and not `≈`,
#    because bit-identical shared corners are what makes the tessellation exact.
#    The two supporting facts (`fl(ix/n) >= fl(iy/n)` iff `ix >= iy`, and
#    `fl(2ix/2n) === fl(ix/n)`) are checked by exhaustion at every `nside` used.
#
# Areas here are computed by an in-test L'Huilier fan from the ring's centroid
# rather than through `GO.area`: on a *densified* ring (hundreds of nearly
# collinear vertices) the two disagree by several percent, and the fan is the
# one that matches an independent Jacobian integral of the chart. On the plain
# 4-gons the two agree bitwise, and the 4-gon areas the `Regridder` computes are
# `test_face_grid.jl`'s business.

using Test
using Printf
import GeometryOps as GO
import GeoInterface as GI

using DiscreteGlobalGrids
import DiscreteGlobalGrids.ISEA as ISEA
import DiscreteGlobalGrids.ISEA4R as ISEA4R
using DiscreteGlobalGrids.ISEA4R: DIAMONDS, Diamond, xyd_to_point, cell_corners,
    cell_center, xyd_to_rowmajor, rowmajor_to_xyd, xyd_to_morton, morton_to_xyd

const US = GO.UnitSpherical

# Recorded so the numbers land in the test log (and in the milestone report).
const MEASURED = Dict{String,Float64}()
record!(key, value) = (MEASURED[key] = max(get(MEASURED, key, -Inf), value))

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

sd(a, b) = US.spherical_distance(a, b)

"`ISEA.VERTICES[v + 1]` as a `UnitSphericalPoint`."
vpoint(v) = (w = ISEA.VERTICES[v+1]; GO.UnitSphericalPoint(w[1], w[2], w[3]))

cross3(a, b) = (a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1])

# Signed area of the ring seen from *outside* the sphere; positive ⇔ CCW.
function ccw_measure(corners)
    acc = (0.0, 0.0, 0.0)
    n = length(corners)
    for i in 1:n
        acc = acc .+ cross3(corners[i], corners[i%n+1])
    end
    outward = reduce((a, b) -> a .+ Tuple(b), corners; init=(0.0, 0.0, 0.0))
    return sum(acc .* outward)
end

"Area of a spherical triangle by L'Huilier's formula (unsigned)."
function lhuilier(a, b, c)
    A = sd(b, c)
    B = sd(a, c)
    C = sd(a, b)
    s = (A + B + C) / 2
    t = tan(s / 2) * tan((s - A) / 2) * tan((s - B) / 2) * tan((s - C) / 2)
    return 4 * atan(sqrt(max(t, 0.0)))
end

"""
Area of a closed spherical ring, as a fan of L'Huilier triangles from the ring's
own centroid. The centroid apex is what makes this correct for the densified,
gently curved rings below — a fan from a *vertex* folds over on them.
"""
function ring_area(pts)
    c = reduce((a, p) -> a .+ (p[1], p[2], p[3]), pts; init=(0.0, 0.0, 0.0))
    c = c ./ sqrt(sum(abs2, c))
    o = GO.UnitSphericalPoint(c[1], c[2], c[3])
    return sum(lhuilier(o, pts[i], pts[mod1(i + 1, length(pts))]) for i in eachindex(pts))
end

"The boundary of cell `(ix, iy)` of `diamond`, `k` chart points per edge."
function densified_cell(ix, iy, diamond, nside, k)
    pts = GO.UnitSphericalPoint{Float64}[]
    corners = ((ix + 1, iy + 1), (ix, iy + 1), (ix, iy), (ix + 1, iy))
    for i in 1:4
        a = corners[i]
        b = corners[i%4+1]
        for t in 0:(k - 1)
            push!(pts, xyd_to_point((a[1] + (b[1] - a[1]) * t / k) / nside,
                (a[2] + (b[2] - a[2]) * t / k) / nside, diamond))
        end
    end
    return pts
end

"The upper-half affine evaluation, written out — the seam-ownership reference."
function upper_half_point(x, y, diamond)
    dm = DIAMONDS[diamond+1]
    p = ISEA.snyder_inv_xyz(dm.upper, dm.cP0 + x * dm.aP + y * dm.bP)
    return GO.UnitSphericalPoint(p[1], p[2], p[3])
end

"The lower-half affine evaluation, written out."
function lower_half_point(x, y, diamond)
    dm = DIAMONDS[diamond+1]
    p = ISEA.snyder_inv_xyz(dm.lower, dm.cQ0 + x * dm.aQ + y * dm.bQ)
    return GO.UnitSphericalPoint(p[1], p[2], p[3])
end

# --------------------------------------------------------------------------
# 0. The layout table
#
# `diamonds.jl` carries its own pin and asserts against it at load; these
# literals are a second, independent copy, so a coordinated edit of the source
# pin still shows up here.
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

@testset "the ten-diamond table is the pinned one" begin
    @test length(DIAMONDS) == 10
    @test all(d -> d isa Diamond, DIAMONDS)
    for d in 0:9
        dm = DIAMONDS[d+1]
        @test dm.verts == EXPECTED_VERTS[d+1]
        @test (dm.upper, dm.lower) == EXPECTED_FACES[d+1]
        # The corner ids really are corners of the faces they are claimed on.
        v00, v10, v11, v01 = dm.verts
        @test issetequal(ISEA.FACE_TRIPLES[dm.upper+1], (v00, v11, v01))
        @test issetequal(ISEA.FACE_TRIPLES[dm.lower+1], (v00, v11, v10))
        # The seam is a real icosahedron edge (adjacent vertices dot to 1/√5).
        @test isapprox(ISEA.vdot(ISEA.VERTICES[v00+1], ISEA.VERTICES[v11+1]),
            ISEA.ADJ_DOT; atol=1e-12)
    end

    # The ten pairs partition the twenty faces.
    faces = sort!(vcat([dm.upper for dm in DIAMONDS], [dm.lower for dm in DIAMONDS]))
    @test faces == collect(0:19)
    @test sort!([dm.upper for dm in DIAMONDS[1:5]]) == [0, 1, 2, 5, 6]
    @test sort!([dm.lower for dm in DIAMONDS[1:5]]) == [3, 4, 8, 11, 12]

    # The anchor convention: northern diamonds hang off base 0, southern ones
    # off its antipode base 11.
    @test all(d -> DIAMONDS[d+1].verts[4] == 0, 0:4)
    @test all(d -> DIAMONDS[d+1].verts[2] == 11, 5:9)
    @test isapprox(ISEA.vdot(ISEA.VERTICES[1], ISEA.VERTICES[12]), -1.0; atol=1e-12)

    # Both affine halves are positively oriented and carry exactly 4π/10 — the
    # ten diamonds account for 4π with nothing left over. Exact, not `isapprox`:
    # the snapped coefficients make it exact (see `_snap_edge`).
    for dm in DIAMONDS
        @test imag(conj(dm.aP) * dm.bP) == 2π / 5
        @test imag(conj(dm.aQ) * dm.bQ) == 2π / 5
        @test isapprox(abs(dm.aP), ISEA.L_PLANE; atol=1e-15)
        @test isapprox(abs(dm.bQ), ISEA.L_PLANE; atol=1e-15)
    end
    @test isapprox(10 * (2π / 5), 4π; rtol=1e-15)
end

# --------------------------------------------------------------------------
# 1. The chart
# --------------------------------------------------------------------------

@testset "chart corners land on the icosahedron vertices" begin
    worst = 0.0
    for d in 0:9
        verts = DIAMONDS[d+1].verts
        for ((x, y), v) in (((0.0, 0.0), verts[1]), ((1.0, 0.0), verts[2]),
            ((1.0, 1.0), verts[3]), ((0.0, 1.0), verts[4]))
            worst = max(worst, sd(xyd_to_point(x, y, d), vpoint(v)))
        end
    end
    @test worst < 1e-14
    record!("chart corner vs ISEA.VERTICES (rad)", worst)
end

@testset "the two halves agree on the seam" begin
    worst = 0.0
    for d in 0:9, t in range(0.0, 1.0; length=257)
        worst = max(worst, sd(upper_half_point(t, t, d), lower_half_point(t, t, d)))
    end
    # The floor here is Snyder's own round-off: the two planar segments are two
    # independent developments of the same icosahedron edge, so they can agree
    # only to the accuracy of `snyder_inv_xyz`. That is exactly why the chart
    # picks ONE of them for the seam (see the ownership testset below) instead
    # of letting neighbouring cells take different ones.
    @test worst < 1e-13
    record!("upper/lower half disagreement on the seam (rad)", worst)
end

@testset "evaluations of one icosahedron vertex cluster to round-off" begin
    by_vertex = Dict{Int,Vector{GO.UnitSphericalPoint{Float64}}}()
    for d in 0:9
        verts = DIAMONDS[d+1].verts
        for ((x, y), v) in (((0.0, 0.0), verts[1]), ((1.0, 0.0), verts[2]),
            ((1.0, 1.0), verts[3]), ((0.0, 1.0), verts[4]))
            push!(get!(by_vertex, v, GO.UnitSphericalPoint{Float64}[]), xyd_to_point(x, y, d))
        end
    end
    # All twelve vertices are hit; bases 0 and 11 five times (they are apexes of
    # five diamonds each), the other ten three times.
    @test sort!(collect(keys(by_vertex))) == collect(0:11)
    @test length(by_vertex[0]) == 5 && length(by_vertex[11]) == 5
    @test all(v -> length(by_vertex[v]) == 3, 1:10)

    worst = 0.0
    for (_v, ps) in by_vertex, a in ps, b in ps
        worst = max(worst, sd(a, b))
    end
    # No snapping to `ISEA.VERTICES` anywhere: the measured spread is already at
    # the Newton floor, snapping would break the pure-chart contract (the chart
    # would stop being one function of `(x, y)`), and the total conservation
    # leakage this leaves is ~1e-13 sr against a 4π·1e-10 budget.
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
            push!(ps, xyd_to_point(i / nside, j / nside, d))
        end
        border[d] = ps
    end
    # Every border lattice point of every diamond has a partner on some other
    # diamond's border: the ten charts tile the sphere with matching lattices.
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
        nside in (1, 2, 3, 4, 5, 7, 8, 16)

    # The integer predicate and the Float64 predicate agree on the whole
    # lattice, so the branch a lattice point takes is exactly decidable.
    for ix in 0:nside, iy in 0:nside
        @test ((iy / nside) >= (ix / nside)) == (iy >= ix)
    end

    # ... and every point on the seam, plus everything above it, evaluates
    # through the upper half. `===`, not `≈`: bit-identical shared corners are
    # what makes the tessellation exact.
    for ix in 0:nside
        @test xyd_to_point(ix / nside, ix / nside, 0) === upper_half_point(ix / nside, ix / nside, 0)
    end
    for d in 0:9, ix in 0:nside, iy in 0:nside
        p = xyd_to_point(ix / nside, iy / nside, d)
        if iy >= ix
            @test p === upper_half_point(ix / nside, iy / nside, d)
        else
            @test p === lower_half_point(ix / nside, iy / nside, d)
        end
    end

    # Cross-resolution nesting is bit-exact: the coordinate of lattice point
    # `ix` at `nside` and of `2ix` at `2nside` is the same `Float64`, so the
    # refined cell's corners are the parent's corners, bitwise.
    if 2nside <= 32
        for ix in 0:nside
            @test (2ix) / (2nside) === ix / nside
        end
        for d in 0:9, ix in 0:nside, iy in 0:nside
            @test xyd_to_point((2ix) / (2nside), (2iy) / (2nside), d) ===
                  xyd_to_point(ix / nside, iy / nside, d)
        end
    end
end

@testset "the chart is exactly equal-area (nside = $nside)" for nside in (1, 2, 3, 4, 5)
    # A chart rectangle of area A covers solid angle A·4π/10, so every cell is
    # `4π/(10 nside²)` — exactly, not on average. Checked on the *densified*
    # boundary (64 chart points per edge): the plain 4-gon is a chord
    # approximation of a curved cell and deviates by up to ~15%, which is
    # `test_face_grid.jl`'s subject, not this one's.
    target = 4π / (10 * nside^2)
    worst = 0.0
    for d in 0:9, ix in 0:(nside - 1), iy in 0:(nside - 1)
        a = ring_area(densified_cell(ix, iy, d, nside, 64))
        worst = max(worst, abs(a - target) / target)
    end
    # The residual is the geodesic-chord discretisation of the curved boundary,
    # and it falls as `k^-2`; it is not a property of the chart.
    @test worst < 1e-4
    record!("densified cell area vs 4π/(10n²) (relative)", worst)
end

@testset "cell_center and cell_corners are chart evaluations" begin
    for nside in (1, 3, 4), d in 0:9, ix in 0:(nside - 1), iy in 0:(nside - 1)
        @test cell_center(ix, iy, d, nside) ===
              xyd_to_point((ix + 0.5) / nside, (iy + 0.5) / nside, d)
        @test cell_corners(ix, iy, d, nside) === (
            xyd_to_point((ix + 1) / nside, (iy + 1) / nside, d),
            xyd_to_point(ix / nside, (iy + 1) / nside, d),
            xyd_to_point(ix / nside, iy / nside, d),
            xyd_to_point((ix + 1) / nside, iy / nside, d),
        )
    end
    # Neighbouring cells share bit-identical corners — the exactness claim.
    nside = 5
    for d in 0:9, ix in 0:(nside - 2), iy in 0:(nside - 1)
        left = cell_corners(ix, iy, d, nside)
        right = cell_corners(ix + 1, iy, d, nside)
        @test left[1] === right[2]      # (ix+1, iy+1)
        @test left[4] === right[3]      # (ix+1, iy)
    end
end

# --------------------------------------------------------------------------
# 2. CCW discipline at the chart level
#
# The convex-clip kernel clips a clockwise ring to EMPTY, so a reversed ring
# yields silent zero intersection areas rather than an error.
# --------------------------------------------------------------------------

@testset "cell_corners winds CCW from outside (nside = $nside)" for nside in (3, 4, 5)
    worst = Inf
    for d in 0:9, ix in 0:(nside - 1), iy in 0:(nside - 1)
        worst = min(worst, ccw_measure(cell_corners(ix, iy, d, nside)))
    end
    @test worst > 0
    record!("min CCW measure, chart level (nside=$nside)", worst)
end

# --------------------------------------------------------------------------
# 3. The chart's index maps
#
# The orderings built on these are `test_face_grid.jl`'s subject; what is
# checked here is the closed forms themselves, including the paths that throw.
# --------------------------------------------------------------------------

@testset "row-major and Morton index maps (nside = $nside)" for nside in (1, 2, 3, 4, 5, 8)
    ncell = 10 * nside^2
    for d in 0:9, ix in 0:(nside - 1), iy in 0:(nside - 1)
        q = xyd_to_rowmajor(ix, iy, d, nside)
        @test 0 <= q < ncell
        @test q == d * nside^2 + iy * nside + ix
        @test rowmajor_to_xyd(q, nside) == (ix, iy, d)
        if ispow2(nside)
            p = xyd_to_morton(ix, iy, d, nside)
            @test 0 <= p < ncell
            @test morton_to_xyd(p, nside) == (ix, iy, d)
        end
    end
    @test_throws ArgumentError rowmajor_to_xyd(-1, nside)
    @test_throws ArgumentError rowmajor_to_xyd(ncell, nside)
    if ispow2(nside)
        @test_throws ArgumentError morton_to_xyd(-1, nside)
        @test_throws ArgumentError morton_to_xyd(ncell, nside)
        @test_throws ArgumentError xyd_to_morton(nside, 0, 0, nside)
        @test_throws ArgumentError xyd_to_morton(0, 0, 10, nside)
    else
        @test_throws ArgumentError xyd_to_morton(0, 0, 0, nside)
        @test_throws ArgumentError morton_to_xyd(0, nside)
    end
end

@testset "the Morton code is a plain bit interleave" begin
    # Written out independently: `ix` into the even bit positions, `iy` into the
    # odd ones. Same convention as HEALPix's NESTED id, 12 → 10.
    interleave(ix, iy) = sum(((ix >> b) & 1) << (2b) | ((iy >> b) & 1) << (2b + 1) for b in 0:15)
    for level in 0:4
        nside = 2^level
        for d in 0:9, ix in 0:(nside - 1), iy in 0:(nside - 1)
            @test xyd_to_morton(ix, iy, d, nside) == d * 4^level + interleave(ix, iy)
        end
    end
    # The quadtree property the radix-4 scaffold ordinal rests on: dropping the
    # low two bits of a level-`r` code is the parent's level-`r-1` code.
    for ix in 0:7, iy in 0:7
        @test xyd_to_morton(ix, iy, 0, 8) ÷ 4 == xyd_to_morton(ix ÷ 2, iy ÷ 2, 0, 4)
    end
end

@printf("[ISEA4R diamonds] measured:\n")
for key in sort!(collect(keys(MEASURED)))
    @printf("[ISEA4R diamonds]   %-52s %+.3e\n", key, MEASURED[key])
end

end # module Isea4rDiamondTestSuite
