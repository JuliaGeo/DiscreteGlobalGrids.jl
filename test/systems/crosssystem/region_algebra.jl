# Region algebra: growth by rings, `union`/`vcat`, bulk level expansion, and
# compaction to a multi-order set. One hex system, one quad system, and A5 —
# whose expansion is the whole reason the level walk may not assume that a
# cell's descendants are contiguous or ascending.

module RegionAlgebraTests

using Test
import DiscreteGlobalGrids as DGG

# (system, root level, leaf level). The leaf level is shallow on purpose: every
# oracle here materializes the subtree it compares against.
const SWEEP = [
    (DGG.IGeo7System(), 2, 4),
    (DGG.HEALPixSystem(), 2, 4),
    (DGG.A5System(), 1, 3),
]

sysname(sys) = string(nameof(typeof(sys)))

# A handful of cells of `l` that are not descendants of `c`, taken from the far
# end of the level so no accidental sibling group is complete.
function elsewhere(sys, c, l, n)
    grid = DGG.levelgrid(sys, l)
    out = [DGG.cellindex(grid, DGG.ncells(grid) - 2k) for k in n:-1:1]
    return sort!(out)
end

@testset "region algebra" begin

    # THE A5 MUTANT. `descendant_range` does not exist on A5 and the native
    # child enumeration is rotated by the level-0 segment, so an expansion that
    # reaches for ranges, or concatenates children in enumeration order, is
    # either a `MethodError` or not ascending. Only the `descendants`-and-sort
    # path answers, and it answers exactly.
    @testset "A5 expansion is the descendant set, exactly" begin
        sys = DGG.A5System()
        cells = [DGG.cellindex(DGG.levelgrid(sys, 0), i) for i in (1, 3, 8, 11)]
        cv = DGG.CellVector(sys, 0, cells)

        @test !DGG.has_sorted_subtrees(sys)
        @test_throws MethodError DGG.descendant_range(sys, first(cells), 1)

        @test collect(DGG.expand(cv, 1)) ==
              sort!(reduce(vcat, [collect(DGG.children(sys, c)) for c in cells]))
        @test collect(DGG.expand(cv, 2)) ==
              sort!(reduce(vcat, [DGG.descendants(sys, c, 2) for c in cells]))
        @test DGG.expand(cv, 0) === cv
        @test_throws ArgumentError DGG.expand(DGG.expand(cv, 2), 1)
    end

    @testset "$(sysname(sys))" for (sys, rootlevel, leaflevel) in SWEEP
        root = DGG.cellindex(DGG.levelgrid(sys, rootlevel), 3)
        pg = DGG.PartialGrid(sys, root, leaflevel)
        cv = DGG.CellVector(pg)
        strays = elsewhere(sys, root, leaflevel, 3)

        # GROWTH IS THE HALO MACHINERY. One ring is the subtree plus exactly the
        # cells `subtree_halo` names, and rings nest.
        @testset "grow" begin
            @test DGG.grow(cv, 0) === cv
            @test Set(DGG.grow(pg, 1)) ==
                  union(Set(DGG.descendants(sys, root, leaflevel)),
                        Set(DGG.subtree_halo(sys, root, leaflevel)))
            @test issubset(DGG.grow(pg, 1), DGG.grow(pg, 2))
            @test length(DGG.grow(pg, 2)) > length(DGG.grow(pg, 1))
        end

        # UNION DEDUPES AND RE-SORTS. Overlapping operands answer the set union
        # once each, ascending, which is the container's invariant rather than
        # Base's first-appearance order.
        @testset "union and vcat" begin
            left = DGG.CellVector(sys, leaflevel, collect(cv)[1:end-2])
            right = DGG.CellVector(sys, leaflevel,
                sort!(vcat(collect(cv)[end-3:end], strays)))
            u = union(left, right)
            @test u isa DGG.CellVector
            @test Set(u) == union(Set(left), Set(right))
            @test issorted(u) && allunique(u)
            @test length(u) < length(left) + length(right)

            # Ordered concatenation survives as a vector only where it ascends.
            @test vcat(left, DGG.CellVector(sys, leaflevel, strays)) isa DGG.CellVector
            @test vcat(right, left) == vcat(collect(right), collect(left))
        end

        # CROSS-LEVEL UNION IS A MULTI-ORDER SET, with no member descending from
        # another: the fine cells under the coarse one are dropped.
        @testset "cross-level union" begin
            fine = DGG.CellVector(sys, leaflevel,
                sort!(vcat(collect(cv)[1:2], strays)))
            set = union(DGG.CellVector(sys, rootlevel, [root]), fine)
            @test set isa DGG.MultiOrderCellSet
            @test Set(set) == Set(vcat(root, strays))
        end

        # COMPACTION IS LOSSLESS. A complete subtree collapses to its root, the
        # strays stay where they are, and expanding the set back names the same
        # cells in the same order.
        @testset "compact round-trips" begin
            mixed = union(cv, DGG.CellVector(sys, leaflevel, strays))
            set = DGG.compact(mixed)
            @test set isa DGG.MultiOrderCellSet
            @test root in collect(set)
            @test minimum(DGG.level, set) == rootlevel
            @test length(set) < length(mixed)
            @test DGG.expand(set, leaflevel) == mixed
        end
    end
end

end # module
