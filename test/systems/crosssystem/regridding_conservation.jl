# ---------------------------------------------------------------------------
# T17 — conservative regridding, in BOTH directions, on every system.
#
# `examples/regridding.jl` and the T7 work checked one direction only: a DGGS
# as the regridder's SOURCE, against a lon/lat destination. That direction is
# conservative to `1e-13` on all seven. Nothing checked the other one, and the
# other one is broken — which is how a package whose headline use case is
# regridding shipped a silently non-conservative destination direction (gap
# inventory entry 3). This file states the law for both directions so that the
# untested half can never be the untested half again.
#
# THE LAW. Build a regridder between a global lon/lat mesh and a DGGS grid.
# Both tile the same sphere, so the intersection matrix `A` must satisfy, in
# whichever order the two grids were given:
#
#     sum(A; dims = 2) == dst_areas      (every destination cell fully covered)
#     sum(A; dims = 1) == src_areas      (every source cell fully consumed)
#     regrid!(out, r, ones) == ones      (the two above, as a user sees them)
#
# and it must satisfy them with the DGGS on EITHER side. `A` is a different
# matrix in the two cases, not a transpose: each is assembled by clipping the
# source cells against the destination cells, and the clipper is not symmetric.
#
# WHY SOME ARMS ARE `@test_broken`, AND HOW THEY DECIDE.
#
# `ConservativeRegridding.DefaultIntersectionOperator` on a `Spherical`
# manifold is `GO.intersection(ConvexConvexSutherlandHodgman(m), p1, p2)`, and
# Sutherland-Hodgman clips `p1` (the SUBJECT) against the half-space of every
# edge of `p2` (the CLIP WINDOW). That is `p1 ∩ p2` only when the CLIP ring is
# convex — a precondition `ConvexConvexSutherlandHodgman`'s own docstring
# states. For a non-convex clip ring it returns the strictly smaller
# `p1 ∩ (⋂ half-spaces of p2)`, and when that comes out empty it decides
# "the subject contains the clip" from the clip ring's FIRST VERTEX ALONE and
# hands back the clip polygon WHOLE.
#
# ConservativeRegridding always passes the SOURCE cell as `p1` and the
# DESTINATION cell as `p2` (`src/regridder/intersection_areas.jl:147-150`), so
# the destination is always the clip window. Hence the asymmetry: a grid whose
# rings are non-convex is clipped correctly as a SOURCE and incorrectly as a
# DESTINATION.
#
# Which grids have non-convex rings is not a property of the system but of the
# system AND the level, so this file does not carry a list — it MEASURES the
# rings (`has_reflex_vertex` below) and expects conservation exactly of the
# grids whose rings are all convex. Today that classification comes out:
#
#     IGeo7, S2, Authalic(IGeo7)  convex at every level swept   -> conservative
#     H3                          convex at L0 and L2, NOT at
#                                 L1 (21 cells in 120) or L3    -> L1 is broken
#     HEALPix, ISEA4R, A5         non-convex at every level     -> broken
#
# HEALPix and ISEA4R and A5 are the systems whose cells are curvilinear, so
# `cell_boundary` densifies each chart edge into eight great-circle segments
# and a straight chart edge is a CURVE on the sphere: two of a cell's four
# sides bow inward and every densified vertex on them is a reflex vertex. That
# is correct geometry, correctly densified — the rings still tile the sphere
# exactly, which is why the source direction is exact and why clipping these
# same cells the other way round (as the subject) gives the right answer to
# `1e-16`. Nothing here is fixable in this package.
#
# UPSTREAM. The defect and its fix are in
# `GeometryOps/src/methods/clipping/sutherland_hodgman.jl`:
#
#   :306-315  clips the subject against every clip edge's great circle,
#             i.e. `subject ∩ (⋂ half-spaces)`, valid only for a convex clip;
#   :317-324  infers "subject contains clip" from `clip_points[1]` alone and
#             returns the WHOLE clip polygon, which turns a grazing sliver into
#             a full-cell credit.
#
# The fix — promote whichever of the two rings is convex into the clip slot,
# since intersection is symmetric and the SUBJECT may be non-convex; and
# require EVERY clip vertex to be inside the subject before declaring
# containment — takes the seven cases below to `4e-13` (ISEA4R), `5e-12`
# (HEALPix) and `1e-8` (A5), all of them passing. When that lands, every
# `@test_broken` here starts reporting "Unexpectedly Pass" and the branch
# below can be deleted along with this paragraph.
#
# TOLERANCE. `1e-10`: three orders above the largest residual any conservative
# arm shows (`2.5e-13`, IGeo7 as the source), and eight orders below the
# smallest defect (`2.2e-2`, A5 level 3). Nothing in between, so the number is
# not load-bearing.
# ---------------------------------------------------------------------------

module RegriddingConservationTests

using Test
import DiscreteGlobalGrids as DGG
import ConservativeRegridding as CR
import GeoInterface as GI
import GeometryOps as GO
using GeometryOps.UnitSpherical: spherical_orient

const TOL = 1e-10

# The manifold every grid in this package computes on. Named once: a bare
# matrix of `UnitSphericalPoint`s carries no manifold, and `best_manifold`
# guesses the WGS84 sphere for it, which is a factor of R^2 in every area.
const MANIFOLD = GO.Spherical(; radius = 1.0)
const TO_SPHERE = GO.UnitSpherical.UnitSphereFromGeographic()

# A 5-degree lon/lat mesh as a corner matrix. Neighbouring quads share their
# two corner POINTS exactly, so the great-circle quads they span partition the
# sphere with no gaps and no overlaps — which is what makes the sums below
# theorems rather than measurements.
const MESH = [TO_SPHERE((x, y)) for x in range(0, 360; length = 73),
                                    y in range(-90, 90; length = 37)]
const MESH_CELLS = 72 * 36

# Does this ring turn right anywhere, i.e. is it non-convex? Repeated vertices
# are skipped: `spherical_orient` goes through `robust_cross_product`, which
# returns an arbitrary perpendicular for two identical points, so a duplicated
# vertex would otherwise read as a random reflex turn.
function has_reflex_vertex(poly)
    pts = collect(GI.getpoint(GI.getexterior(poly)))
    while length(pts) > 1 && pts[end] == pts[1]
        pop!(pts)
    end
    unique!(pts)
    n = length(pts)
    n < 4 && return false
    return any(1:n) do i
        spherical_orient(pts[i], pts[mod1(i + 1, n)], pts[mod1(i + 2, n)]) < 0
    end
end

all_rings_convex(grid) = !any(has_reflex_vertex(DGG.cell_polygon(grid, DGG.cellindex(grid, i)))
                              for i in 1:DGG.ncells(grid))

"Largest relative disagreement between the matrix's marginals and the cell areas."
function conservation_errors(r)
    A = r.intersections
    row = maximum(abs.(vec(sum(A; dims = 2)) .- r.dst_areas) ./ r.dst_areas)
    col = maximum(abs.(vec(sum(A; dims = 1)) .- r.src_areas) ./ r.src_areas)
    return row, col
end

"`regrid!` of a field of ones; conservative iff it comes back as ones."
function regrid_ones(r)
    out = zeros(length(r.dst_areas))
    CR.regrid!(out, r, ones(length(r.src_areas)))
    return out
end

base(sys) = sys isa DGG.AuthalicSystem ? parent(sys) : sys
label(sys) = sys isa DGG.AuthalicSystem ?
             "Authalic($(nameof(typeof(parent(sys)))))" : string(nameof(typeof(sys)))

# Coarse on purpose — the defect is a property of a cell's shape, not of how
# many there are, and it shows at every level (see the header's table).
demo_level(sys) = base(sys) isa DGG.H3System ? 2 :
                  base(sys) isa DGG.IGeo7System ? 2 : 3

# Every registered system at its demo level, plus the authalic wrap, plus H3 at
# level 1 — the case that stops this file from reading as "three systems are
# fine and three are not". H3's hexagons are convex at level 0 and level 2 and
# a minority of them are NOT at level 1, and level 1 fails: it is ring
# convexity that decides, not the system.
cases() = [(sys, demo_level(sys)) for sys in DGG.systems()] ∪
          [(DGG.AuthalicSystem(DGG.IGeo7System()), 2), (DGG.H3System(), 1)]

@testset "regridding conserves in both directions" begin
    for (sys, l) in cases()
        @testset "$(label(sys)) level $l" begin
            grid = DGG.levelgrid(sys, l)
            convex = all_rings_convex(grid)

            # ---- the DGGS as the regridder's SOURCE. Correct on every system:
            # the DGGS cell lands in the subject slot, where Sutherland-Hodgman
            # accepts a non-convex ring.
            forward = CR.Regridder(MANIFOLD, MESH, grid)
            @test size(forward.intersections) == (MESH_CELLS, DGG.ncells(grid))
            row, col = conservation_errors(forward)
            @test col <= TOL
            @test row <= TOL
            @test all(v -> isapprox(v, 1.0; atol = TOL), regrid_ones(forward))

            # ---- the DGGS as the regridder's DESTINATION. THE SAME TWO GRIDS,
            # the other way round. Correct exactly when the DGGS rings are
            # convex, because the destination cell is the clip window.
            reverse = CR.Regridder(MANIFOLD, grid, MESH)
            @test size(reverse.intersections) == (DGG.ncells(grid), MESH_CELLS)
            rrow, rcol = conservation_errors(reverse)
            ones_back = regrid_ones(reverse)
            if convex
                @test rrow <= TOL
                @test rcol <= TOL
                @test all(v -> isapprox(v, 1.0; atol = TOL), ones_back)
            else
                # UPSTREAM DEFECT — GeometryOps
                # `src/methods/clipping/sutherland_hodgman.jl:306-324` clips
                # against a non-convex clip window out of contract. See this
                # file's header for the mechanism and the fix. Delete this
                # branch, not the assertions, when the fix lands.
                @test_broken rrow <= TOL
                @test_broken rcol <= TOL
                @test_broken all(v -> isapprox(v, 1.0; atol = TOL), ones_back)
            end

            # Both directions agree about how much sphere there is, whatever
            # the clipper did to the individual weights — so a failure above is
            # never a disagreement about the grids themselves.
            @test sum(forward.src_areas) ≈ 4pi rtol = 1e-12
            @test sum(reverse.dst_areas) ≈ 4pi rtol = 1e-12
            @test sum(forward.dst_areas) ≈ sum(reverse.src_areas) rtol = 1e-14
        end
    end
end

end # module RegriddingConservationTests
