# test/IGeo7/test_gbt_neighbors.jl — the GBT digit-arithmetic edge neighbors
# (`IGeo7._cell_neighbors`, src/IGeo7/gbt_neighbors.jl).
#
# The tables in gbt_neighbors.jl are ported from IGEO7.jl, so the thing worth
# testing is not that they encode adjacency — they do, upstream — but that
# they encode adjacency *in this module's digit conventions*, which were fitted
# independently. That is a claim about two implementations agreeing, so the
# test is a differential one: `_cell_neighbors` against
# `_cell_neighbors_geometric`, the clean-room lattice implementation that
# test_neighbors.jl separately pins to the module's oracle-validated cell
# boundaries. Adjacency itself is ground-truthed there; here it is only
# transferred.
#
# Exhaustive where exhaustive is affordable (every cell through res 5 —
# 196,092 of them, and the geometric reference costs ~2us each), sampled above
# it, and hand-aimed at the two places the port could plausibly diverge: the
# pentagon chains and the base-cell crossings.
#
# Big loops accumulate a failure counter and assert once (the idiom of
# test_indexing.jl) so the suite's test count stays readable.

module IGeo7GBTNeighborTests

using Test
using Random

using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
using DiscreteGlobalGrids: Helpers, IGeo7

"All cells of resolution `res` in canonical (ascending id) order."
level_cells(res) = (IGeo7.index_to_cell(i, res) for i in 1:IGeo7.num_cells(res))

"The pentagon of `base` at resolution `res` (its all-zero-digit descendant)."
function pentagon(base::Integer, res::Integer)
    z = (UInt64(base) << 60) | IGeo7.Z7_PAD_MASK
    for _ in 1:res
        z = IGeo7.z7_child(z, 0)
    end
    return z
end

@testset "GBT neighbors" begin

    @testset "container and contract" begin
        z = IGeo7.index_to_cell(1000, 6)
        @test IGeo7._cell_neighbors(z) isa Helpers.SmallList{6,UInt64}
        @test length(IGeo7._cell_neighbors(z)) == 6
        @test issorted(IGeo7._cell_neighbors(z))
        @test allunique(IGeo7._cell_neighbors(z))
        # allocation-free on the scalar path
        IGeo7._cell_neighbors(z)
        @test (@allocated IGeo7._cell_neighbors(z)) == 0

        # invalid ids and res-20 ids keep the module's error contract
        @test_throws IGeo7.InvalidZ7Error IGeo7._cell_neighbors(UInt64(0xffffffffffffffff))
        res20 = IGeo7.index_to_cell(1, 19)
        res20 = IGeo7.z7_child(res20, 1)
        @test IGeo7.z7_resolution(res20) == 20
        @test_throws IGeo7.InvalidZ7Error IGeo7._cell_neighbors(res20)
    end

    # The ported EXCLUDE_NEIGHBOURS and the independently fitted
    # Z7_DELETED_DIGIT are the same table; pinned rather than aliased so the
    # port and the fit stay separate evidence.
    @testset "ported exclusion table matches the fitted deleted digit" begin
        @test IGeo7.EXCLUDE_NEIGHBOURS == IGeo7.Z7_DELETED_DIGIT
    end

    @testset "agrees with the geometric reference, res 0:5 exhaustive" begin
        for res in 0:5
            bad = 0
            first_bad = nothing
            for z in level_cells(res)
                if IGeo7._cell_neighbors(z) != IGeo7._cell_neighbors_geometric(z)
                    bad += 1
                    first_bad === nothing && (first_bad = z)
                end
            end
            @test bad == 0
            bad == 0 || @info "res $res first mismatch" cell = IGeo7.z7_to_string(first_bad)
        end
    end

    @testset "agrees with the geometric reference, sampled res 6:19" begin
        rng = Random.MersenneTwister(20260813)
        for res in (6, 8, 10, 12, 15, 19)
            n = IGeo7.num_cells(res)
            bad = 0
            for _ in 1:2000
                z = IGeo7.index_to_cell(rand(rng, 1:n), res)
                IGeo7._cell_neighbors(z) == IGeo7._cell_neighbors_geometric(z) || (bad += 1)
            end
            @test bad == 0
        end
    end

    @testset "pentagons have five neighbors, all bases, res 0:8" begin
        bad = 0
        for res in 0:8, base in 0:11
            z = pentagon(base, res)
            nbs = IGeo7._cell_neighbors(z)
            length(nbs) == 5 || (bad += 1)
            nbs == IGeo7._cell_neighbors_geometric(z) || (bad += 1)
        end
        @test bad == 0
    end

    # A cell whose neighbors leave its base cell is where the carry, the
    # adjacency table and the frame rotations all have to line up. Take the
    # rim of each base's subtree so every crossing is exercised.
    @testset "base-cell crossings" begin
        bad = 0
        crossings = 0
        for res in 1:4, base in 0:11
            root = (UInt64(base) << 60) | IGeo7.Z7_PAD_MASK
            for z in IGeo7.border_descendants(root, res)
                nbs = IGeo7._cell_neighbors(z)
                crossings += count(n -> IGeo7.z7_base_cell(n) != base, nbs)
                nbs == IGeo7._cell_neighbors_geometric(z) || (bad += 1)
            end
        end
        @test bad == 0
        @test crossings > 0   # the fixture actually crosses
    end

    @testset "symmetry: adjacency is mutual" begin
        bad = 0
        for res in 0:4, z in level_cells(res)
            for n in IGeo7._cell_neighbors(z)
                z in IGeo7._cell_neighbors(n) || (bad += 1)
            end
        end
        @test bad == 0
    end

    # The kernel wiring and everything above it (halo tables, stencils) now
    # run on the GBT path; check the re-seating into `SmallVector` preserved
    # the contract.
    @testset "kernel wiring" begin
        S = IGEO7DGGS()
        bad = 0
        for res in 0:3, z in level_cells(res)
            collect(DGG.cell_neighbors(S, res, z)) == collect(IGeo7._cell_neighbors(z)) ||
                (bad += 1)
        end
        @test bad == 0
    end
end

end # module IGeo7GBTNeighborTests
