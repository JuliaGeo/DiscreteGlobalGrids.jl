module Isea9rDelegationTestSuite

# ISEA9R has no chart of its own: `src/ISEA9R/chart.jl` *imports* `ISEA4R`'s
# (`xyd_to_point`, `cell_corners`, `cell_center`, the row-major index maps and
# the `DIAMONDS` table) because the rhombus chart takes continuous `(x, y)` and
# quantises nothing, so it carries no aperture. That is a load-bearing claim —
# it is the reason no diamond table, no seam rule and no cap argument was
# re-derived for aperture 9 — and this file is where it stops being a claim.
#
# Four things are pinned, in rising order of strength:
#
# 1. *The delegation is literally an import.* The two modules' chart bindings
#    are `===` the same function objects. A copy-paste fork would pass every
#    numeric test in `test/ISEA9R/` and fail here.
#
# 2. *The lattice properties survive.* Shared lattice points come out
#    bit-identical between adjacent cells, across the diamond seam, and across
#    resolutions `n` → `3n` — the aperture-9 reading of the `fl(ix/n) ===
#    fl(3ix/3n)` argument the `MortonOrder` docstring and
#    `docs/design/isea4r_diamond_layout.md` §3.3 rest on. `test/ISEA4R/` proves
#    these at `2^k`; the exponent is not incidental to the argument, so they are
#    re-proved at `3^k` here rather than assumed.
#
# 3. *The two systems are the same grid at a common resolution.* Every polygon
#    of `Isea9rFaceGrid(n)` is bitwise a polygon of `Isea4rFaceGrid(n)`, under
#    the permutation the two orderings define.
#
# 4. *...and therefore so is every Regridder built on them.* Against a common
#    destination the two intersection matrices are EXACTLY equal under that
#    permutation — not `isapprox`. This is the ISEA4R/ISEA9R analogue of the
#    HEALPix ring↔nested exact-permutation test, one step stronger because it
#    crosses two systems rather than two orderings of one.
#
# If a future milestone hoists the shared chart into a common submodule (an
# `ISEARhombic`, or `ISEA` itself — deliberately deferred, see
# `docs/design/isea9r_layout.md` §5), this file is the regression net for that
# move: (1) changes shape, (2)-(4) must not move at all.

using Test
using Printf
import GeometryOps as GO
import GeoInterface as GI
import ConservativeRegridding as CR
import ConservativeRegridding: Trees

using DiscreteGlobalGrids
using DiscreteGlobalGrids.HEALPix
import DiscreteGlobalGrids as DGG
import DiscreteGlobalGrids.ISEA9R as ISEA9R
import DiscreteGlobalGrids.ISEA4R as ISEA4R
using DiscreteGlobalGrids.ISEA9R: Isea9rFaceSpace, Isea9rFaceGrid

const US = GO.UnitSpherical
const UNIT = GO.Spherical(radius=1.0)

const MEASURED = Dict{String,Float64}()
record!(key, value) = (MEASURED[key] = max(get(MEASURED, key, -Inf), value))

identical(a, b) = all(k -> a[k] === b[k], 1:3)
ring_points(poly) = collect(GI.getpoint(GI.getexterior(poly)))
open_ring(poly) = ring_points(poly)[1:4]

# --------------------------------------------------------------------------
# 1. The delegation is an import, not a copy
# --------------------------------------------------------------------------

@testset "ISEA9R's chart bindings ARE ISEA4R's" begin
    # The chart and everything derived from it.
    @test ISEA9R.xyd_to_point === ISEA4R.xyd_to_point
    @test ISEA9R.cell_corners === ISEA4R.cell_corners
    @test ISEA9R.cell_center === ISEA4R.cell_center
    # The row-major index maps, which carry no aperture at all.
    @test ISEA9R.xyd_to_rowmajor === ISEA4R.xyd_to_rowmajor
    @test ISEA9R.rowmajor_to_xyd === ISEA4R.rowmajor_to_xyd
    # ...and the layout table itself: one table, not two copies that could drift.
    @test ISEA9R.DIAMONDS === ISEA4R.DIAMONDS
    @test length(ISEA9R.DIAMONDS) == 10

    # The base-9 Morton codec is the one thing ISEA9R adds, so it must NOT be
    # the sibling's radix-4 one.
    @test ISEA9R.xyd_to_morton !== ISEA4R.xyd_to_morton
    @test ISEA9R.morton_to_xyd !== ISEA4R.morton_to_xyd

    # `face_chart` is the shared layer's contract method: both systems route it
    # to the same function, so the two charts are one chart by construction.
    for diamond in 0:9, (x, y) in ((0.0, 0.0), (0.25, 0.75), (1 / 3, 1 / 3), (1.0, 1.0))
        @test identical(DGG.face_chart(ISEA9R.Isea9rFaceSystem(), x, y, diamond),
                        DGG.face_chart(ISEA4R.Isea4rFaceSystem(), x, y, diamond))
    end

    # ...but the systems themselves are distinct, and so are their bounds and
    # ordering families — sharing geometry is not sharing identity.
    @test ISEA9R.Isea9rFaceSystem() !== ISEA4R.Isea4rFaceSystem()
    @test DGG.max_nside(ISEA9R.Isea9rFaceSystem()) == 3^18
    @test DGG.max_nside(ISEA4R.Isea4rFaceSystem()) == 2^29
    @test DGG.ordering_family(ISEA9R.Isea9rFaceSystem()) === ISEA9R.AbstractIsea9rOrdering
    @test DGG.nfaces(ISEA9R.Isea9rFaceSystem()) == DGG.nfaces(ISEA4R.Isea4rFaceSystem()) == 10
    # Both opted into the O(1) block cap — each on its own measurement, at its
    # own block shapes (see the two `cap_policy` overrides).
    @test DGG.cap_policy(ISEA9R.Isea9rFaceSystem()) isa DGG.FourCornerCap
    @test DGG.cap_policy(ISEA4R.Isea4rFaceSystem()) isa DGG.FourCornerCap
end

# --------------------------------------------------------------------------
# 2. Lattice bit-exactness at 3^k
#
# The `2^k` versions of these live in `test/ISEA4R/test_diamonds.jl`. The
# exponent is not incidental — the whole argument is about which `Float64` a
# division lands on — so they are re-run here at the aperture-9 resolutions.
# --------------------------------------------------------------------------

@testset "shared lattice points are bit-identical between neighbours (nside = $nside)" for
        nside in (3, 9, 27)

    corner(ix, iy, diamond) = ISEA4R.cell_corners(ix, iy, diamond, nside)
    for diamond in 0:9, ix in 0:(nside - 1), iy in 0:(nside - 1)
        here = corner(ix, iy, diamond)
        # East neighbour shares the two `x = (ix+1)/n` corners. Ring order is
        # `(TR, TL, BL, BR)`, so `here[1] === east[2]` and `here[4] === east[3]`.
        if ix + 1 < nside
            east = corner(ix + 1, iy, diamond)
            @test identical(here[1], east[2])
            @test identical(here[4], east[3])
        end
        # North neighbour shares the two `y = (iy+1)/n` corners.
        if iy + 1 < nside
            north = corner(ix, iy + 1, diamond)
            @test identical(here[1], north[4])
            @test identical(here[2], north[3])
        end
    end
end

# The seam (`y == x`, the glued icosahedron edge) is the case that could break:
# the two halves of the chart are two independent developments of one edge and
# agree only to ~5e-15 rad, so a cell that took the lower development where its
# neighbour took the upper one would leave a gap. The rule — the half `y >= x`,
# seam included, evaluates through the upper face — makes the choice a function
# of the lattice point alone. Diagonal cells `(k, k)` straddle the seam and are
# where it is tested.
@testset "the seam-ownership rule is exact on the 3^k lattice (nside = $nside)" for
        nside in (3, 9, 27)

    # The integer predicate and the floating-point one decide identically...
    for i in 0:nside, j in 0:nside
        @test (Float64(j / nside) >= Float64(i / nside)) == (j >= i)
    end
    # ...so seam-straddling cells share bit-identical corners with both
    # neighbours, exactly as interior cells do.
    for diamond in 0:9, k in 0:(nside - 2)
        here = ISEA4R.cell_corners(k, k, diamond, nside)
        east = ISEA4R.cell_corners(k + 1, k, diamond, nside)
        north = ISEA4R.cell_corners(k, k + 1, diamond, nside)
        @test identical(here[1], east[2]) && identical(here[4], east[3])
        @test identical(here[1], north[4]) && identical(here[2], north[3])
    end
end

# The claim `MortonOrder` and `supports_prefix_ranges(ISEA9RDGGS())` rest on:
# `fl(ix/3^k) === fl(3ix/3^(k+1))`. Both sides are correctly-rounded quotients
# of exactly-representable integers with an identical real value, so they are
# the identical `Float64` — which makes a level-L cell's corners a bit-identical
# SUBSET of the level-(L+1) lattice, and the aperture-9 nesting exact rather
# than merely consistent to rounding.
@testset "cross-level corner sharing: level L corners are level L+1 lattice points" begin
    for level in 0:3
        n = 3^level
        fine = 3n
        # (a) the float identity itself, exhaustively over the coarse lattice.
        for ix in 0:n
            @test Float64(ix) / Float64(n) === Float64(3ix) / Float64(3n)
        end
        # (b) and therefore, corner for corner, on the sphere.
        for diamond in 0:9, ix in 0:(n - 1), iy in 0:(n - 1)
            coarse = ISEA4R.cell_corners(ix, iy, diamond, n)
            # The coarse cell's four corners are four points of the fine
            # lattice: its own bottom-left, and the corners of the 3x3 block.
            block = ISEA4R.cell_corners(3ix, 3iy, diamond, fine)          # BL sub-cell
            far = ISEA4R.cell_corners(3ix + 2, 3iy + 2, diamond, fine)    # TR sub-cell
            @test identical(coarse[3], block[3])                          # BL
            @test identical(coarse[1], far[1])                            # TR
            # The other two corners come from the off-diagonal sub-cells.
            @test identical(coarse[2], ISEA4R.cell_corners(3ix, 3iy + 2, diamond, fine)[2])
            @test identical(coarse[4], ISEA4R.cell_corners(3ix + 2, 3iy, diamond, fine)[4])
        end
    end
end

# The same statement one level up, through the ordinal: the nine children of
# `isea9r_ordinal` p tile p's chart rectangle, and the union of their corners
# contains p's four corners bitwise. This is the geometric content of
# `descendant_range` being exact over this ordinal.
@testset "the nine base-9 Morton children tile their parent's rectangle" begin
    for level in 0:1
        n = 3^level
        fine = 3n
        for p in 0:(10 * n^2 - 1)
            ix, iy, diamond = ISEA9R.morton_to_xyd(p, n)
            parent = ISEA4R.cell_corners(ix, iy, diamond, n)
            kid_corners = Set{NTuple{3,Float64}}()
            for k in 0:8
                cx, cy, cd = ISEA9R.morton_to_xyd(9p + k, fine)
                @test cd == diamond
                @test cx ÷ 3 == ix && cy ÷ 3 == iy
                for c in ISEA4R.cell_corners(cx, cy, cd, fine)
                    push!(kid_corners, (c[1], c[2], c[3]))
                end
            end
            for c in parent
                @test (c[1], c[2], c[3]) in kid_corners
            end
        end
    end
end

# --------------------------------------------------------------------------
# 3. THE CROSS-SYSTEM GATE: one chart, one lattice, one grid
#
# `Isea4rFaceGrid(n)` and `Isea9rFaceGrid(n)` are the same chart at the same
# lattice, so their polygon SETS must be bitwise identical up to the permutation
# their orderings define. Anything less would mean the aperture leaked into the
# geometry.
# --------------------------------------------------------------------------

"""
σ[i] = the ISEA4R row-major data position of the cell that ISEA9R's `ordering`
puts at data position `i`. For `ISEA9R.RowMajorOrder` this is the identity (both
systems compute the same row-major index — the arithmetic has no aperture in
it); for `ISEA9R.MortonOrder` it is a genuine permutation.
"""
function sigma_to_isea4r(ordering, nside)
    space = Isea9rFaceSpace(nside)
    n4 = ISEA4R.Isea4rFaceSpace(nside)
    return [DGG.data_index(ISEA4R.RowMajorOrder(), n4,
                DGG.lattice_index(ordering, space, i)...) for i in 1:(10nside^2)]
end

@testset "Isea4rFaceGrid(n) and Isea9rFaceGrid(n) are bitwise the same grid (nside = $nside)" for
        nside in (3, 9)

    ncell = 10nside^2
    i4 = treeify(ISEA4R.Isea4rFaceGrid(nside; ordering=ISEA4R.RowMajorOrder()))

    for ordering in (ISEA9R.RowMajorOrder(), ISEA9R.MortonOrder())
        i9 = treeify(Isea9rFaceGrid(nside; ordering))
        sigma = sigma_to_isea4r(ordering, nside)
        @test sort(sigma) == collect(1:ncell)                  # a permutation
        if ordering isa ISEA9R.RowMajorOrder || nside <= 3
            # Row-major over `10 nside²` cells is one expression with no
            # aperture in it, so the two systems agree slot for slot.
            #
            # Base-9 Morton coincides with it at `nside <= 3`, and that is
            # arithmetic rather than luck: a ONE-digit base-9 code is
            # `ix + 3 * iy`, which is exactly the row-major index `iy * 3 + ix`.
            # So the two ISEA9R orderings are the same ordering at levels 0 and
            # 1 and first diverge at level 2 (`nside = 9`), where the second
            # base-9 digit puts `ix`'s high digit below `iy`'s low one. Worth
            # pinning: it is the reason the permutation test below has to run at
            # `nside = 9` to be testing anything at all.
            @test sigma == collect(1:ncell)
        else
            @test sigma != collect(1:ncell)                     # genuinely permuted
        end

        # Bitwise, corner for corner — `==` on the ring vectors would already be
        # strict, but `identical` also refuses `0.0` for `-0.0`.
        for i in 1:ncell
            a = ring_points(Trees.getcell(i9, i))
            b = ring_points(Trees.getcell(i4, sigma[i]))
            @test length(a) == length(b) == 5
            @test all(k -> identical(a[k], b[k]), 1:5)
        end
    end
end

# ...and the same statement through the clipper: Regridders built on the two
# grids against a COMMON destination must be exactly equal under the same
# permutation. Exact, not `isapprox`: both sides hand the clipper literally the
# same polygons and assembly with `threaded = false` is deterministic, so any
# drift here would mean the two systems' geometry had forked.
@testset "Regridders against a common grid are exactly equal under σ (nside = $nside)" for
        nside in (3, 9)

    common = treeify(HealpixFaceGrid(4; ordering=RingOrder()))
    i4 = treeify(ISEA4R.Isea4rFaceGrid(nside; ordering=ISEA4R.RowMajorOrder()))
    R4 = CR.Regridder(UNIT, i4, common; threaded=false, normalize=false)

    for ordering in (ISEA9R.RowMajorOrder(), ISEA9R.MortonOrder())
        i9 = treeify(Isea9rFaceGrid(nside; ordering))
        R9 = CR.Regridder(UNIT, i9, common; threaded=false, normalize=false)
        sigma = sigma_to_isea4r(ordering, nside)

        @test size(R9.intersections) == size(R4.intersections) == (10nside^2, 12 * 16)
        # Rows are destination cells, so the ISEA permutation acts on rows only.
        @test R9.intersections == R4.intersections[sigma, :]
        @test R9.dst_areas == R4.dst_areas[sigma]
        @test R9.src_areas == R4.src_areas                     # the common grid
        # Conservation, on both.
        @test isapprox(sum(R9.intersections), 4π; rtol=1e-10)
        record!("|sum(intersections) - 4π|, cross-system gate nside=$nside",
            abs(sum(R9.intersections) - 4π))
    end
end

# The delegation is not restricted to `3^k`: at ANY nside both systems evaluate
# the same chart, and the row-major orderings agree slot for slot. `nside = 4`
# is the interesting case — it is a resolution where ISEA4R has a Morton id
# space and ISEA9R does not, so the grids can only be compared through
# row-major, and they are still identical.
@testset "the two systems agree at non-power-of-three nside too (nside = $nside)" for
        nside in (2, 4, 5)

    i4 = treeify(ISEA4R.Isea4rFaceGrid(nside; ordering=ISEA4R.RowMajorOrder()))
    i9 = treeify(Isea9rFaceGrid(nside; ordering=ISEA9R.RowMajorOrder()))
    for j in 1:(10nside^2)
        a = ring_points(Trees.getcell(i9, j))
        b = ring_points(Trees.getcell(i4, j))
        @test all(k -> identical(a[k], b[k]), 1:5)
    end
end

# --------------------------------------------------------------------------
# 4. The DGGS-level statement of the same thing
#
# At level 0 the two ordinals coincide (ten diamonds, one cell each) and the two
# systems' `cell_polygon` answers are bitwise identical. Above level 0 they are
# different id spaces over the same geometry, which the lattice comparison below
# states precisely.
# --------------------------------------------------------------------------

@testset "cell_polygon agrees across the two registry systems" begin
    for id in 0:9
        a = ring_points(cell_polygon(ISEA9RDGGS(), 0, id))
        b = ring_points(cell_polygon(ISEA4RDGGS(), 0, id))
        @test all(k -> identical(a[k], b[k]), 1:5)
    end
    # Level 2 of ISEA9R is nside 9; level 2 of ISEA4R is nside 4. Different
    # lattices, so the comparison has to go through the lattice, not the id:
    # ISEA9R level 2 cell `(ix, iy, d)` is ISEA4R's *nside-9* cell of the same
    # lattice coordinates — which ISEA4R can only address through its
    # `RowMajorOrder` grid, not through its ordinal (9 is not `2^k`).
    root4 = treeify(ISEA4R.Isea4rFaceGrid(9; ordering=ISEA4R.RowMajorOrder()))
    space4 = ISEA4R.Isea4rFaceSpace(9)
    for id in 0:809
        ix, iy, d = ISEA9R.morton_to_xyd(id, 9)
        a = ring_points(cell_polygon(ISEA9RDGGS(), 2, id))
        b = ring_points(Trees.getcell(root4,
            DGG.data_index(ISEA4R.RowMajorOrder(), space4, ix, iy, d)))
        @test all(k -> identical(a[k], b[k]), 1:5)
    end
end

@printf("[ISEA9R delegation] measured:\n")
for key in sort!(collect(keys(MEASURED)))
    @printf("[ISEA9R delegation]   %-52s %+.3e\n", key, MEASURED[key])
end

end # module Isea9rDelegationTestSuite
