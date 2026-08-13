# `DGGSSubtreeIds` and the generic neighbor-stepper layer
# (`src/core/subtree_ids.jl`). Everything here is system-agnostic. No system
# currently specializes `neighbor_stepper`, so the generic stepper below is
# also the one every wired system actually takes.

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

# `edge_cells` / `interior_cells` are the same two sets without the arrays, so
# the contract is that they are indistinguishable from the materialized forms —
# and, for the interior, that its two independent access paths (the merge walk
# in `iterate`, the binary search in `getindex`) never disagree.
@testset "edge_cells and interior_cells are lazy views of the same sets" begin
    S = IGEO7DGGS()
    for (rl, root, lvl) in ((5, 0x0c4d9fffffffffff, 9),   # hexagon-rooted
                            (0, 0x2fffffffffffffff, 4),   # pentagon-rooted
                            (3, 0x0c4fffffffffffff, 3))   # depth 0: all rim
        t = DGGSSubtreeIds(S, rl, root, lvl)
        e, v = edge_cells(t), interior_cells(t)

        @test e isa AbstractVector{Int}
        @test v isa AbstractVector{Int}
        @test collect(e) == subtree_border_positions(t)
        @test collect(v) == subtree_interior_positions(t)
        @test length(e) + length(v) == length(t)
        @test issorted(e) && issorted(v)
        @test sort(vcat(collect(e), collect(v))) == collect(1:length(t))

        # iterate (merge) and getindex (binary search) are separate code paths
        @test [v[k] for k in eachindex(v)] == collect(v)
        @test [e[k] for k in eachindex(e)] == collect(e)

        @test_throws BoundsError v[length(v)+1]
        @test_throws BoundsError e[length(e)+1]
        length(v) > 0 && @test_throws BoundsError v[0]

        # ids are recoverable through the tile
        @test [t[i] for i in e] == DGG.subtree_border(S, rl, root, lvl)
    end
end

# The interior is the whole point of the laziness: it is ~97% of a tile, so
# materializing it is close to writing down `1:n`.
@testset "interior_cells does not materialize the interior" begin
    S = IGEO7DGGS()
    t = DGGSSubtreeIds(S, 1, 0x33ffffffffffffff, 7)
    e = edge_cells(t)
    interior_cells(t, e)                        # warm
    @test length(t) > 100_000                   # the fixture is big enough to matter
    @test (@allocated interior_cells(t, e)) == 0
    # and eager interior storage is what it avoids
    @test length(interior_cells(t, e)) > 0.9 * length(t)
end

# The complement logic is index arithmetic over an ascending set; check it
# against a brute-force complement on shapes the tile fixtures do not reach
# (empty rim, full rim, rim at both ends, singletons).
@testset "interior complement arithmetic, brute force" begin
    S = IGEO7DGGS()
    t = DGGSSubtreeIds(S, 5, 0x0c4d9fffffffffff, 7)
    n = length(t)
    for edge in (Int[], [1], [n], [1, n], collect(1:n), collect(1:2:n),
                 collect(2:2:n), [1, 2, 3, n-1, n], collect(1:n-1), [n÷2])
        v = interior_cells(t, edge)
        expected = setdiff(1:n, edge)
        @test length(v) == length(expected)
        @test collect(v) == expected
        @test [v[k] for k in eachindex(v)] == expected
    end
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
