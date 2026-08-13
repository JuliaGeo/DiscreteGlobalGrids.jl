# `DGGSSubtreeIds` and the generic neighbor-stepper layer
# (`src/core/subtree_ids.jl`). Everything here is system-agnostic; the IGEO7
# twiddle that specializes `neighbor_stepper` is tested in
# `test/IGeo7/test_tile_neighbors.jl`.

using Test
using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
import SmallCollections

@testset "DGGSSubtreeIds is the descendant vector, lazily" begin
    for (S, rl, root, leaf) in ((IGEO7DGGS(), 5, 0x0c4d9fffffffffff, 9),
                                (IGEO7DGGS(), 2, DGG.IGeo7.z7_from_string("0313"), 6),
                                (HEALPixDGGS(), 3, 100, 6))
        t = DGGSSubtreeIds(S, rl, root, leaf)
        want = DGG.cell_descendants(S, rl, root, leaf)
        @test length(t) == length(want)
        @test collect(t) == want
        @test issorted(t)
        # the point of the type: it does not hold them
        @test Base.summarysize(t) < 200
    end
end

@testset "membership and position are exact" begin
    S, rl, root, leaf = IGEO7DGGS(), 5, 0x0c4d9fffffffffff, 9
    t = DGGSSubtreeIds(S, rl, root, leaf)
    ids = collect(t)

    @test all(id -> id in t, ids)
    @test all(i -> subtree_position(t, ids[i]) == i, eachindex(ids))

    # a sibling subtree's cells are outside, and report position 0
    sib = DGG.cell_descendants(S, rl, 0x0c4dbfffffffffff, leaf)
    @test !any(id -> id in t, sib)
    @test all(id -> subtree_position(t, id) == 0, sib)
end

@testset "rim and interior partition the tile" begin
    S = IGEO7DGGS()
    t = DGGSSubtreeIds(S, 5, 0x0c4d9fffffffffff, 9)
    border = subtree_border_positions(t)
    interior = subtree_interior_positions(t, border)

    @test issorted(border)
    @test issorted(interior)
    @test isempty(intersect(Set(border), Set(interior)))
    @test length(border) + length(interior) == length(t)
    @test sort(vcat(border, interior)) == collect(1:length(t))

    # the rim is `subtree_border`, positioned
    @test [t[i] for i in border] == DGG.subtree_border(S, 5, 0x0c4d9fffffffffff, 9)

    # and it is genuinely the truncated-neighborhood set
    s = GenericNeighborStepper(t)
    @test Set(findall(i -> length(step_neighbors(s, i)) < DGG.max_neighbors(S),
        1:length(t))) == Set(border)
end

@testset "generic stepper agrees with cell_neighbors" begin
    S = HEALPixDGGS()
    t = DGGSSubtreeIds(S, 3, 100, 6)
    s = neighbor_stepper(t)
    for i in 1:length(t)
        want = sort!([subtree_position(t, nb) for nb in DGG.cell_neighbors(S, t.level, t[i])])
        filter!(!iszero, want)
        @test collect(step_neighbors(s, i)) == want
    end
end

@testset "table stepper is a faithful materialization" begin
    S = IGEO7DGGS()
    t = DGGSSubtreeIds(S, 5, 0x0c4d9fffffffffff, 8)
    computed = neighbor_stepper(t)
    table = neighbor_table(t)
    @test all(i -> collect(step_neighbors(table, i)) == collect(step_neighbors(computed, i)),
        1:length(t))
end

@testset "subtree_stencil sweeps every cell exactly once" begin
    S = IGEO7DGGS()
    t = DGGSSubtreeIds(S, 5, 0x0c4d9fffffffffff, 8)
    data = collect(1.0:length(t))

    # count the neighbors each cell saw
    counts = subtree_stencil((c, nb) -> length(nb), data, t)
    s = neighbor_stepper(t)
    @test counts == [length(step_neighbors(s, i)) for i in 1:length(t)]

    # a real reduction, and the border/interior split does not change the answer
    mean_nb = subtree_stencil((c, nb) -> isempty(nb) ? c : sum(nb) / length(nb), data, t)
    reference = [let v = step_neighbors(s, i)
            isempty(v) ? data[i] : sum(data[j] for j in v) / length(v)
        end for i in 1:length(t)]
    @test mean_nb ≈ reference

    # the table stepper is a drop-in
    @test subtree_stencil((c, nb) -> length(nb), data, t; stepper = neighbor_table(t)) == counts

    @test_throws DimensionMismatch subtree_stencil((c, nb) -> 1.0, data[1:end-1], t)
end

@testset "construction errors" begin
    S = IGEO7DGGS()
    @test_throws ArgumentError DGGSSubtreeIds(S, 5, 0x0c4d9fffffffffff, 3)   # leaf < root
    @test_throws ArgumentError DGGSSubtreeIds(S, -1, 0x0c4d9fffffffffff, 5)
    # A5 has no descendant ranges, so a subtree is not an id interval
    @test !DGG.has_descendant_ranges(A5DGGS())
    @test_throws ArgumentError DGGSSubtreeIds(A5DGGS(), 1, 0, 3)
end
