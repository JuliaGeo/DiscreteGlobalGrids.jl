# Test conservative regridding in both directions between each DGGS and a
# global longitude/latitude mesh. Marginal intersection areas are compared with
# independently computed source and destination cell areas.

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

demo_level(sys) = base(sys) isa DGG.H3System ? 2 :
                  base(sys) isa DGG.IGeo7System ? 2 : 3

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

            # Use the same grids with the DGGS as the destination. This
            # direction conserves when the DGGS ring is a convex clip window.
            reverse = CR.Regridder(MANIFOLD, grid, MESH)
            @test size(reverse.intersections) == (DGG.ncells(grid), MESH_CELLS)
            rrow, rcol = conservation_errors(reverse)
            ones_back = regrid_ones(reverse)
            if convex
                @test rrow <= TOL
                @test rcol <= TOL
                @test all(v -> isapprox(v, 1.0; atol = TOL), ones_back)
            else
                @test_broken rrow <= TOL
                @test_broken rcol <= TOL
                @test_broken all(v -> isapprox(v, 1.0; atol = TOL), ones_back)
            end

            @test sum(forward.src_areas) ≈ 4pi rtol = 1e-12
            @test sum(reverse.dst_areas) ≈ 4pi rtol = 1e-12
            @test sum(forward.dst_areas) ≈ sum(reverse.src_areas) rtol = 1e-14
        end
    end
end

end # module RegriddingConservationTests
