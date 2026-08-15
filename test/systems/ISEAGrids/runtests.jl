module ISEAGridsTests

using Test
using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
import GeometryOps as GO
using DiscreteGlobalGridsConformanceTesting

const G = DGG.ISEAGrids
const ISEA = DGG.ISEA
const ORACLES = normpath(joinpath(@__DIR__, "..", "..", "oracles"))

function center_rows(family, level)
    file = joinpath(ORACLES, family, "dggrid-9.0b",
        "level-$(lpad(level, 2, '0'))-centers.txt")
    return [(s = split(line, ',');
        (String(s[1]), parse(Float64, s[2]), parse(Float64, s[3]))) for line in eachline(file)]
end

angular_degrees(a, b) = ISEA.angdist(Tuple(a), Tuple(b))

function i4h_cell(text)
    level = length(text) - 2
    root = parse(Int, text[1:2])
    path = level == 0 ? 0 : parse(Int, text[3:end], base=4)
    return DGG.LevelIndex(level, G._i4encode(root, path, level))
end

function i4h_text(c)
    root, path = G._i4decode(c.index, level(c))
    suffix = level(c) == 0 ? "" : lpad(string(path, base=4), level(c), '0')
    return lpad(root, 2, '0') * suffix
end

@testset "ISEA literature grids" begin
    @testset "counts, ids, and primary prefix trees" begin
        for (sys, A, maxlevel) in ((G.ISEA3HSystem(), 3, 5), (G.ISEA4HSystem(), 4, 4))
            for level in 0:maxlevel
                grid = levelgrid(sys, level)
                @test ncells(grid) == 10A^level + 2
                @test all(cellposition(grid, cellindex(grid, i)) == i for i in 1:ncells(grid))
                @test count(G.is_pentagon, (cellindex(grid, i) for i in 1:ncells(grid))) == 12
            end
            for root in rootcells(sys), child in children(sys, root)
                @test parent(sys, child) == root
            end
        end

        @test G.z3_string(G.Z3Cell("071201")) == "071201"
        @test_throws ArgumentError G.Z3Cell("001") # north polar root has only zero chain
        @test length(children(G.ISEA3HSystem(), G.Z3Cell("000"))) == 1
        # The sealed labels are 00, 01..10, 11 (decimal roots): one polar
        # chain, ten full root blocks, one polar chain.  This pin protects
        # against the older prose error that treated roots 10 and 11 as the
        # two poles.
        @test map(c -> length(children(G.ISEA3HSystem(), c)),
            rootcells(G.ISEA3HSystem())) == [1; fill(3, 10); 1]
        @test map(c -> length(children(G.ISEA4HSystem(), c)),
            rootcells(G.ISEA4HSystem())) == [1; fill(4, 10); 1]

        # Central-place geometry is not nested; the relation above is the
        # explicitly chosen indexing tree.  Still, every descendant block is
        # contiguous in canonical root-major prefix order.
        sys = G.ISEA3HSystem(); c = G.Z3Cell("0210")
        r = descendant_range(sys, c, 5)
        @test length(r) == 3^3
        @test all(ancestor(sys, cellindex(sys, 5, i), 2) == c for i in r)
    end

    @testset "DGGRID black-box center oracles" begin
        for level in 0:5
            sys = G.ISEA3HSystem()
            for (text, lon, lat) in center_rows("ISEA3H", level)
                c = G.Z3Cell(text)
                expected = ISEA.lonlat_to_xyz(lon, lat)
                @test angular_degrees(cell_centroid(sys, c), expected) < 1.2e-6
            end
        end
        for level in 0:4
            sys = G.ISEA4HSystem()
            for (text, lon, lat) in center_rows("ISEA4H", level)
                c = i4h_cell(text)
                expected = ISEA.lonlat_to_xyz(lon, lat)
                @test angular_degrees(cell_centroid(sys, c), expected) < 1.2e-6
            end
        end
    end

    @testset "hex boundaries and analytic nominal areas" begin
        for sys in (G.ISEA3HSystem(), G.ISEA4HSystem())
            for level in 0:3
                grid = levelgrid(sys, level)
                for i in 1:ncells(grid)
                    c = cellindex(grid, i)
                    ring = cell_boundary(grid, c)
                    @test length(ring) == 8 * (G.is_pentagon(c) ? 5 : 6)
                    @test all(abs(sum(abs2, p) - 1) < 2e-14 for p in ring)
                    @test cellat(grid, cell_centroid(grid, c)) == c
                end
            end
        end
        @test isapprox(sum(G.equal_area_steradians(G.Z3Cell(text))
            for (text, _, _) in center_rows("ISEA3H", 3)), 4pi; atol=5e-14)
        for sys in (G.ISEA3HSystem(), G.ISEA4HSystem()), level in 0:4
            grid = levelgrid(sys, level)
            @test sum(cell_area(grid, cellindex(grid, i)) for i in 1:ncells(grid)) ≈ 4pi
            polygon_areas = [GO.area(GO.Spherical(radius=1.0),
                cell_polygon(grid, cellindex(grid, i))) for i in 1:ncells(grid)]
            analytic = [cell_area(grid, cellindex(grid, i)) for i in 1:ncells(grid)]
            # Eight Snyder-face samples are a finite polygon approximation;
            # pin both global closure and the observed per-cell error envelope.
            @test sum(polygon_areas) ≈ 4pi rtol=5e-4
            @test maximum(abs.(polygon_areas .- analytic) ./ analytic) < 0.04
        end

        # Exact, exhaustive one-ring comparison with the sealed black-box
        # oracle at every acquired level.
        for (family, sys, maxlevel, parsecell, showcell) in (
                ("ISEA3H", G.ISEA3HSystem(), 5, G.Z3Cell, G.z3_string),
                ("ISEA4H", G.ISEA4HSystem(), 4, i4h_cell, i4h_text))
            for l in 0:maxlevel
                grid = levelgrid(sys, l)
                file = joinpath(ORACLES, family, "dggrid-9.0b",
                    "level-$(lpad(l, 2, '0'))-neighbors.txt")
                for line in eachline(file)
                    fields = split(line)
                    got = Set(showcell.(neighbors(grid, parsecell(fields[1]))))
                    @test got == Set(String.(fields[2:end]))
                end
            end
        end

        # The locators are O(level) bounded-candidate searches. Exercise the
        # advertised deep levels, well beyond the acquired oracle depths.
        for sys in (G.ISEA3HSystem(), G.ISEA4HSystem()), l in (10, 20, last(levels(sys)))
            grid = levelgrid(sys, l); total = ncells(grid)
            for i in unique(round.(Int, range(1, total; length=97)))
                c = cellindex(grid, i)
                @test cellat(grid, cell_centroid(grid, c)) == c
                @test all(c in neighbors(grid, n) for n in neighbors(grid, c))
            end
        end
    end

    @testset "ISEA4T exact face/path subdivision" begin
        sys = G.ISEA4TSystem()
        for level in 0:5
            grid = levelgrid(sys, level)
            @test ncells(grid) == 20 * 4^level
            for i in 1:ncells(grid)
                c = cellindex(grid, i)
                @test cellposition(grid, c) == i
                @test length(cell_boundary(grid, c)) == 24
                @test cellat(grid, cell_centroid(grid, c)) == c
                level > 0 && @test c in children(sys, parent(sys, c))
            end
        end
        @test all(cell_area(levelgrid(sys, 3), DGG.LevelIndex(3, i)) == pi / (5 * 4^3)
            for i in 0:(20 * 4^3 - 1))

        # DGGRID exposes only SEQNUM for this family.  Without inventing a
        # face/path crosswalk, compare the complete center set at levels 0:2.
        for level in 0:2
            grid = levelgrid(sys, level)
            ours = [cell_centroid(grid, cellindex(grid, i)) for i in 1:ncells(grid)]
            for (_, lon, lat) in center_rows("ISEA4T", level)
                p = ISEA.lonlat_to_xyz(lon, lat)
                @test minimum(angular_degrees(q, p) for q in ours) < 1.2e-6
            end
        end

        # The advertised ceiling is a numerical geometry ceiling, not merely
        # an integer-code limit. Probe the deepest supported levels directly.
        for level in (9, 12)
            grid = levelgrid(sys, level)
            step = max(1, ncells(grid) ÷ 97)
            for i in 1:step:ncells(grid)
                c = cellindex(grid, i)
                @test cellat(grid, cell_centroid(grid, c)) == c
                edge = neighbors(grid, c; connectivity=Edge())
                @test length(edge) == 3
                @test all(c in neighbors(grid, n; connectivity=Edge()) for n in edge)
            end
        end
        @test last(levels(sys)) == 12


        for level in 0:4
            grid = levelgrid(sys, level)
            for i in 1:ncells(grid)
                c = cellindex(grid, i)
                edge = neighbors(grid, c; connectivity=Edge())
                @test length(edge) == 3
                @test all(c in neighbors(grid, n; connectivity=Edge()) for n in edge)
            end
        end
        for level in 0:2
            grid = levelgrid(sys, level)
            @test all(9 <= length(neighbors(grid, cellindex(grid, i);
                connectivity=Vertex())) <= 12 for i in 1:ncells(grid))
        end
    end


    # `cap_inflation` is the one constant on these two systems that is a
    # correctness bug when it is too small and a silent tax on every tree descent
    # when it is too large. It was 4.5 and 3.5, from a triangle-inequality bound
    # that ignores the rotation of the central-place digit directions; the true
    # planar suprema are 2 and sqrt(3) (see `hex.jl`). Nothing reproduced either
    # number, which is how the oversize survived.
    #
    # This is that measurement: the smallest factor that satisfies the covering
    # law, exhaustively over every ancestor at the low levels. It kills the mutant
    # in both directions — a factor lowered under the geometry, and a refinement
    # change that moves the geometry out from under the factor.
    @testset "cap_inflation covers the subtree, and is not wildly over it" begin
        function worst_factor(sys, levels, depth)
            worst = 0.0
            for la in levels
                ga = levelgrid(sys, la)
                for i in 1:ncells(ga)
                    c = cellindex(ga, i)
                    cap = DGG.Fallbacks.cell_cap(ga, c)
                    for j in 1:depth
                        gd = levelgrid(sys, la + j)
                        for d in DGG.descendants(sys, c, la + j)
                            capd = DGG.Fallbacks.cell_cap(gd, d)
                            f = (GO.UnitSpherical.spherical_distance(cap.point,
                                     capd.point) + capd.radius) / cap.radius
                            worst = max(worst, f)
                        end
                    end
                end
            end
            return worst
        end

        for (sys, levels, depth, expected) in
            ((G.ISEA3HSystem(), 0:3, 5, 1.9953), (G.ISEA4HSystem(), 0:2, 4, 1.6788))

            worst = worst_factor(sys, levels, depth)
            # The geometry has not moved.
            @test worst ≈ expected atol = 1e-3
            # The declared factor covers it, with margin but not a factor of two.
            @test worst <= DGG.cap_inflation(sys)
            @test DGG.cap_inflation(sys) <= 1.3 * worst
            # And every extent it produces is a convex cap, which is what lets
            # the conformance call below keep `require_convex_extents` on.
            g0 = levelgrid(sys, 0)
            @test maximum(DGG.node_extent(sys, cellindex(g0, i)).radius
                          for i in 1:ncells(g0)) < π / 2
        end
    end

    @testset "package interface conformance" begin
        test_hierarchical_system(G.ISEA3HSystem(); levels=0:8, n_levels=5,
            n_samples=16)
        test_hierarchical_system(G.ISEA4HSystem(); levels=0:8, n_levels=5,
            n_samples=16)
        test_hierarchical_system(G.ISEA4TSystem(); levels=0:8, n_levels=5,
            n_samples=16)
    end
end

end # module
