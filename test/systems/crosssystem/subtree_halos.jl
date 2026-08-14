# ---------------------------------------------------------------------------
# The OUTSIDE face of a subtree boundary.
#
# `subtree_iterators.jl` walks the inside of that boundary — the descendants
# with a neighbour that is not one. This file walks the outside: the level-`l`
# cells that are NOT descendants but have a neighbour that is. Same boundary,
# opposite side, so the two files share a fixture vocabulary and nothing else.
#
# Written against the generic interface only: no system module is imported, and
# every law runs against every system in `systems()`.
# ---------------------------------------------------------------------------

module SubtreeHaloTests

using Test
using DiscreteGlobalGrids
using DiscreteGlobalGrids: systems, levelgrid, level, levels, max_level, ncells,
    cellindex, cellposition, neighbors, ancestor, descendants, descendant_range,
    subtree_border, has_sorted_subtrees, Vertex, Edge, Connectivity,
    SubtreeHaloIterator, subtree_halo, halo
import DiscreteGlobalGrids as DGG

# ---------------------------------------------------------------------------
# Depth zero: a cell's own one-ring
# ---------------------------------------------------------------------------

@testset "depth zero is the cell's own one-ring" begin
    for sys in systems()
        grid = levelgrid(sys, 1)
        c = cellindex(grid, 1)
        for conn in (Vertex(), Edge())
            expected = sort!(collect(neighbors(grid, c, 1; connectivity = conn)))
            it = SubtreeHaloIterator(sys, c, 1; connectivity = conn)
            @test collect(it) == expected
            @test subtree_halo(sys, c, 1; connectivity = conn) == expected
            @test eltype(it) == DGG.cellindextype(sys)
        end
    end
end

@testset "level validation" begin
    for sys in systems()
        grid = levelgrid(sys, 1)
        c = cellindex(grid, 1)
        @test_throws ArgumentError SubtreeHaloIterator(sys, c, 0)
        @test_throws ArgumentError SubtreeHaloIterator(sys, c, max_level(sys) + 1)
    end
end

end # module
