# test/core/test_globe_trees.jl — the tree layer over a globe-complete lookup
# (src/core/lookups.jl; docs/design/full_globe_lookups.md §3).
#
# The claim of the step: a lookup whose ids are a `DGGSGlobeIds` treeifies to
# the *dense* `DGGSGrid` cursor rather than to the `DGGSPartialGrid` cursor its
# id vector would otherwise take — decided by the type of that vector, never by
# its length. Three things follow, each checked below for all four kernel-wired
# systems:
#
#   * the tree is a whole-globe tree at any level, built in O(1). Down the
#     partial path a res-15 H3 globe is an id vector of some 4.6 PB, so a
#     regression is an `OutOfMemoryError`, not a slow test.
#   * alignment survives: leaf index `i` is lookup position `i`, which is what
#     makes a `ConservativeRegridding.Regridder` line up with a `DimArray` over
#     the dimension column for column. The partial cursor gets that from
#     `grid.ids === l.data`; the dense cursor gets it from leaf index = ordinal,
#     which a `DGGSGlobeIds` position *is* by definition. Same guarantee, a
#     different argument — hence its own test rather than the partial lookups'.
#   * `DGGSPartialGrid(l)` is guarded, not forbidden. A caller may explicitly
#     want the partial cursor over a globe, so the constructor stays open and
#     its O(N) check is the one short-circuited by type.
#
# One module rather than four per-system additions, because the implementation
# is one generic method over `AbstractDGGSLookup`. `num_cells` is deliberately
# unexported from the package (the system submodules own that name), so it is
# qualified here and never enters the shared test namespace.
module GlobeTreeTests

using Test
import SparseArrays

using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
import ConservativeRegridding as CR
import GeoInterface as GI
import GeometryOps as GO
import GeometryOps: SpatialTreeInterface as STI

using DiscreteGlobalGrids.A5.A5Lookups: A5Lookup
using DiscreteGlobalGrids.H3.H3Lookups: H3Lookup
using DiscreteGlobalGrids.HEALPix.HealpixLookups: HealpixLookup
using DiscreteGlobalGrids.IGeo7.IGeo7Lookups: IGeo7Lookup

# One entry per kernel-wired system: its lookup type, the system, a level whose
# globe is small enough to compare cell for cell, and one whose explicit id
# vector could not exist. The levels are `test/core/test_globe_ids.jl`'s.
const SYSTEMS = (
    (H3Lookup, H3DGGS(), 0, 15),
    (A5Lookup, A5DGGS(), 1, 25),
    (IGeo7Lookup, IGEO7DGGS(), 1, 19),
    (HealpixLookup, HEALPixDGGS(), 1, 29),
)

ring_points(polygon) = collect(GI.getpoint(GI.getexterior(polygon)))

all_leaf_indices(tree) = STI.depth_first_search(Returns(true), tree)

@testset "a globe lookup treeifies to the dense grid" begin
    for (Lookup, system, level, huge) in SYSTEMS
        l = Lookup(DGGSGlobeIds(system, level))
        total = Int(DGG.num_cells(system, level))
        tree = treeify(l)
        @test tree isa DGGSCursor
        @test tree.grid isa DGGSGrid            # the whole point: not partial
        @test tree.grid.system === system
        @test tree.grid.level == level
        @test ncells(tree) == total == length(l)
        # Indistinguishable from the grid's own tree, which is what makes
        # "globe-complete dimension" and "whole grid" the same thing to every
        # consumer downstream.
        @test tree == treeify(DGGSGrid(system, level))
        # The manifold-explicit form is the one `Regridder` itself calls.
        @test DGG.Trees.treeify(GO.Spherical(), l) == tree

        # Alignment, cell for cell over the whole globe: the traversal emits
        # every position exactly once, in order, and leaf `i` is the cell the
        # lookup holds at `i`.
        @test all_leaf_indices(tree) == collect(1:total)
        for i in (1, 2, total ÷ 2, total)
            @test ring_points(getcell(tree, i)) ==
                  DGG.cell_boundary(system, level, l[i]; closed=true)
        end

        # ...and the same at a level whose id vector is unrepresentable, where
        # only the dense path can answer at all.
        huge_tree = treeify(Lookup(DGGSGlobeIds(system, huge)))
        @test huge_tree.grid isa DGGSGrid
        @test huge_tree.grid.level == huge
        @test ncells(huge_tree) == DGG.num_cells(system, huge)
        @test ring_points(getcell(huge_tree, 1)) == DGG.cell_boundary(
            system, huge, DGG.ordinal_to_cell(system, huge, 1); closed=true)
    end
end

# The other branch of the same method, which the globe branch must not have
# taken over: a slice is an ordinary `Vector`, hence an ordinary partial tree.
@testset "a sliced globe lookup treeifies to the partial cursor" begin
    for (Lookup, system, _, huge) in SYSTEMS
        sliced = Lookup(DGGSGlobeIds(system, huge))[1:8]
        tree = treeify(sliced)
        @test tree.grid isa DGGSPartialGrid
        @test tree.grid.ids === sliced.data     # passed through, not copied
        @test ncells(tree) == 8
        @test all_leaf_indices(tree) == collect(1:8)
    end
end

# ---------------------------------------------------------------------------
# `ConservativeRegridding.Regridder` over a globe lookup — the consumer the
# alignment guarantee exists for. A self-regrid cannot see a permutation, since
# both of its axes would carry it, so the claim is made against a destination
# grid holding one *named* cell: the column that cell's whole area lands in is
# its position in the lookup. Off by even one position and the column would
# hold a shared-edge sliver or nothing at all, never an area.
# ---------------------------------------------------------------------------
@testset "a Regridder over a globe lookup maps leaf i to position i" begin
    for (Lookup, system, level, _) in SYSTEMS
        l = Lookup(DGGSGlobeIds(system, level))
        total = length(l)

        self = CR.Regridder(l, l; threaded=false, normalize=false)
        @test size(self.intersections) == (total, total)
        @test SparseArrays.nnz(self.intersections) >= total
        # The diagonal only, and loosely: A5's coarse rings are not geodesically
        # convex in either of its frames, so the spherical clipper reports real
        # off-diagonal overlap there and misses ~2e-6 of each cell's own area
        # (the A5 Regridder testset owns that fact and asserts the same way).
        # Every system's cell still regrids onto itself.
        @test all(i -> isapprox(self.intersections[i, i], self.dst_areas[i]; rtol=1e-5),
            1:total)

        for k in (1, 2, total ÷ 2, total)
            probe = CR.Regridder(DGGSPartialGrid(system, level, [l[k]]), l;
                threaded=false, normalize=false)
            @test size(probe.intersections) == (1, total)
            # The whole cell, up to the clipper: A5's non-convex coarse rings
            # cost a couple of parts per million here, orders of magnitude less
            # than the gap between "the cell itself" and the shared-edge sliver
            # or structural zero a misaligned column would hold.
            @test isapprox(probe.intersections[1, k], only(probe.dst_areas); rtol=1e-5)
            # ...and the position really is the ordinal, which is the argument
            # the dense branch rests on: a position in a `DGGSGlobeIds` *is*
            # `cell_to_ordinal` of the cell there, hence the dense cursor's own
            # leaf index.
            @test DGG.cell_to_ordinal(system, level, l[k]) == k
        end
    end
end

@testset "DGGSPartialGrid over a globe lookup stays constructible and O(1)" begin
    for (Lookup, system, level, huge) in SYSTEMS
        # `treeify` routes a globe to the dense grid, but the partial cursor is
        # still a legitimate thing to ask for, so the constructor is guarded
        # rather than closed: strict ascent is answered from the type
        # (`globe_ids.jl`) instead of walking `num_cells` ids, and every other
        # check reads only the two endpoints. A regression here does not fail
        # this test, it hangs it.
        l = Lookup(DGGSGlobeIds(system, huge))
        grid = DGGSPartialGrid(l)
        @test grid.ids === l.data               # passed through, never collected
        @test grid.level == huge
        @test length(grid.ids) == DGG.num_cells(system, huge)
        @test grid.root_level == -1
        @test DGGSPartialGrid(l; bucket_size=8).bucket_size == 8

        # And it means what it says: over a globe small enough to walk, the
        # partial cursor's leaves are the same cells in the same order as the
        # dense cursor's — the two alignment arguments meeting on one grid.
        small = Lookup(DGGSGlobeIds(system, level))
        partial_tree = treeify(DGGSPartialGrid(small))
        @test partial_tree.grid isa DGGSPartialGrid
        @test ncells(partial_tree) == length(small)
        @test all_leaf_indices(partial_tree) == collect(1:length(small))
        for i in (1, length(small))
            @test ring_points(getcell(partial_tree, i)) ==
                  ring_points(getcell(treeify(small), i))
        end
    end
end

end # module GlobeTreeTests
