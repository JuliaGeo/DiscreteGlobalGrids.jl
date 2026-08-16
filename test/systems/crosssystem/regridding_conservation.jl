# Test conservative regridding in both directions between each DGGS and a
# global longitude/latitude mesh. Marginal intersection areas are compared with
# independently computed source and destination cell areas.
#
# Both grids tile the same sphere, so the intersection matrix `A` must satisfy
# `sum(A; dims = 2) == dst_areas`, `sum(A; dims = 1) == src_areas`, and
# `regrid!(out, r, ones) == ones`, with the DGGS on either side. `A` is a
# different matrix in the two cases, not a transpose, because the clipper is
# not symmetric. The DGGS-as-source direction is conservative to `1e-13`.
#
# WHY SOME ARMS ARE `@test_broken`. On a `Spherical` manifold the default
# intersection operator is Sutherland-Hodgman, which clips `p1` (the subject)
# against the half-space of every edge of `p2` (the clip window). That is
# `p1 ∩ p2` only when the CLIP ring is convex, a precondition
# `ConvexConvexSutherlandHodgman`'s own docstring states. For a non-convex clip
# ring it returns the strictly smaller `p1 ∩ (⋂ half-spaces of p2)`, and when
# that is empty it infers "the subject contains the clip" from the clip ring's
# first vertex alone and returns the clip polygon whole. ConservativeRegridding
# always passes the source as `p1` and the destination as `p2`
# (`src/regridder/intersection_areas.jl:147-150`), so a grid with non-convex
# rings is clipped correctly as a source and incorrectly as a destination.
#
# Convexity is a property of the system AND the level, so this file measures the
# rings (`has_reflex_vertex` below) rather than carrying a list, and expects
# conservation exactly of the grids whose rings are all convex. That
# classification currently comes out:
#
#     IGeo7, S2, Authalic(IGeo7)  convex at every level swept   -> conservative
#     H3                          convex at L0 and L2, NOT at
#                                 L1 (21 of 120 sampled, and
#                                 150 of all 842) or L3         -> L1 is broken
#     HEALPix, ISEA4R, A5         non-convex at every level     -> broken
#
# The non-convex systems are the curvilinear ones: `cell_boundary` densifies
# each chart edge into eight great-circle segments, and a straight chart edge is
# a curve on the sphere, so two of a cell's four sides bow inward and every
# densified vertex on them is reflex. That is correct geometry — the rings still
# tile the sphere exactly, which is why the source direction is exact and why
# clipping the same cells as the subject is right to `1e-16`. It is not fixable
# in this package.
#
# UPSTREAM, in `GeometryOps/src/methods/clipping/sutherland_hodgman.jl`:
# `:306-315` clips the subject against every clip edge's great circle, valid
# only for a convex clip; `:317-324` infers containment from `clip_points[1]`
# alone and returns the whole clip polygon, turning a grazing sliver into a
# full-cell credit. The fix — promote whichever ring is convex into the clip
# slot, and require every clip vertex to be inside the subject before declaring
# containment — takes the seven cases below to `4e-13` (ISEA4R), `5e-12`
# (HEALPix) and `1e-8` (A5), all passing. When it lands, these `@test_broken`
# arms start reporting "Unexpectedly Pass".
#
# TOLERANCE `1e-10`: three orders above the largest residual any conservative
# arm shows (`2.5e-13`, IGeo7 as the source), and eight orders below the
# smallest defect (`2.2e-2`, A5 level 3). Nothing lies in between, so the number
# is not load-bearing.
#
# ISEA3H/4H expose a documented finite approximation to their Snyder edges.
# Their analytic areas are exact, but those polygons are not yet an
# implementation-gating conservative-regridding surface. The broken arms below
# pin that limitation, so a future crack-free canonical edge construction turns
# into an Unexpectedly Pass instead of silently changing the contract.

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

const MESH = [TO_SPHERE((x, y)) for x in range(0, 360; length = 73),
                                    y in range(-90, 90; length = 37)]
const MESH_CELLS = 72 * 36

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
approximate_boundary(sys) = base(sys) isa Union{DGG.ISEA3HSystem,DGG.ISEA4HSystem}

demo_level(sys) = base(sys) isa DGG.H3System ? 2 :
                  base(sys) isa DGG.IGeo7System ? 2 : 3

cases() = [(sys, demo_level(sys)) for sys in DGG.systems()] ∪
          [(DGG.AuthalicSystem(DGG.IGeo7System()), 2), (DGG.H3System(), 1)]

@testset "regridding conserves in both directions" begin
    for (sys, l) in cases()
        @testset "$(label(sys)) level $l" begin
            grid = DGG.levelgrid(sys, l)
            convex = all_rings_convex(grid)
            approximate = approximate_boundary(sys)

            # ---- the DGGS as the regridder's SOURCE. Correct on every system:
            # the DGGS cell lands in the subject slot, where Sutherland-Hodgman
            # accepts a non-convex ring.
            forward = CR.Regridder(MANIFOLD, MESH, grid)
            @test size(forward.intersections) == (MESH_CELLS, DGG.ncells(grid))
            row, col = conservation_errors(forward)
            @test col <= TOL
            if approximate
                @test_broken row <= TOL
                @test_broken all(v -> isapprox(v, 1.0; atol = TOL), regrid_ones(forward))
            else
                @test row <= TOL
                @test all(v -> isapprox(v, 1.0; atol = TOL), regrid_ones(forward))
            end

            # Use the same grids with the DGGS as the destination. This
            # direction conserves when the DGGS ring is a convex clip window.
            reverse = CR.Regridder(MANIFOLD, grid, MESH)
            @test size(reverse.intersections) == (DGG.ncells(grid), MESH_CELLS)
            rrow, rcol = conservation_errors(reverse)
            ones_back = regrid_ones(reverse)
            if convex && !approximate
                @test rrow <= TOL
                @test rcol <= TOL
                @test all(v -> isapprox(v, 1.0; atol = TOL), ones_back)
            else
                @test_broken rrow <= TOL
                @test_broken rcol <= TOL
                @test_broken all(v -> isapprox(v, 1.0; atol = TOL), ones_back)
            end

            # Both directions agree about how much sphere there is, whatever the
            # clipper did to individual weights, so a failure above is never a
            # disagreement about the grids themselves. The systems with
            # approximate edges cannot close to 4pi and are pinned broken.
            if approximate
                @test_broken sum(forward.src_areas) ≈ 4pi rtol = 1e-12
                @test_broken sum(reverse.dst_areas) ≈ 4pi rtol = 1e-12
            else
                @test sum(forward.src_areas) ≈ 4pi rtol = 1e-12
                @test sum(reverse.dst_areas) ≈ 4pi rtol = 1e-12
            end
            @test sum(forward.dst_areas) ≈ sum(reverse.src_areas) rtol = 1e-14
        end
    end
end

end # module RegriddingConservationTests
