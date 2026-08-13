# The Z7 tile neighbor stepper (`src/IGeo7/tile_neighbors.jl`): edge adjacency
# inside a subtree from the digits alone.
#
# The whole file is one claim — the integer twiddle agrees with the geometric
# `cell_neighbors` on every cell of every tile it accepts — plus the boundary
# of where it stops applying (pentagon roots) and the properties that make it
# usable as a stencil (ascending, rim self-describing).

using Test
using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG

const S = IGEO7DGGS()

# The reference: `cell_neighbors` resolved through the tile, which is the
# generic stepper. Spelled out rather than reused so a bug in
# `GenericNeighborStepper` cannot make both sides agree.
function reference_neighbors(t::DGGSSubtreeIds, i::Int)
    out = Int[]
    for nb in DGG.cell_neighbors(S, t.level, t[i])
        pos = subtree_position(t, nb)
        pos == 0 || push!(out, pos)
    end
    return sort!(out)
end

agrees(t, stepper, i) = collect(step_neighbors(stepper, i)) == reference_neighbors(t, i)

@testset "Z7 tile neighbors == geometric neighbors" begin
    # The tutorial tile, exhaustively at a depth small enough to be quick and
    # sampled at the depth a real workflow uses.
    root = 0x0c4d9fffffffffff
    for leaf in 5:10
        t = DGGSSubtreeIds(S, 5, root, leaf)
        s = neighbor_stepper(t)
        @test s isa DGG.IGeo7.Z7TileNeighborStepper
        @test all(i -> agrees(t, s, i), 1:length(t))
    end

    t = DGGSSubtreeIds(S, 5, root, 12)
    s = neighbor_stepper(t)
    @test all(i -> agrees(t, s, i), 1:1237:length(t))
end

@testset "across base cells and root-level parities" begin
    # The chirality alternates on absolute level parity, so a tile rooted at an
    # odd level and one rooted at an even level exercise different arithmetic.
    for base in (0, 3, 7, 11), root_level in 1:5
        z = DGG.IGeo7.z7_from_string(string(base; base = 10, pad = 2) *
                                     "1" * "3"^(root_level - 1))
        DGG.IGeo7.z7_is_pentagon(z) && continue
        for leaf in root_level:(root_level + 3)
            t = DGGSSubtreeIds(S, root_level, z, leaf)
            s = neighbor_stepper(t)
            @test all(i -> agrees(t, s, i), 1:length(t))
        end
    end
end

@testset "pentagon roots fall back rather than lie" begin
    # A pentagon subtree is `p(d) = (5*7^d + 1)/6` cells, not `7^d`, so the
    # suffix is not a base-7 numeral and the construction does not apply.
    for base in 0:2
        pent = DGG.IGeo7.z7_from_string(string(base; base = 10, pad = 2) * "000")
        @test DGG.IGeo7.z7_is_pentagon(pent)
        t = DGGSSubtreeIds(S, 3, pent, 6)
        @test length(t) == DGG.subtree_leaf_count(S, 3, pent, 6)
        @test length(t) != 7^3                      # the reason for the guard
        s = neighbor_stepper(t)
        @test s isa GenericNeighborStepper          # not the twiddle
        @test all(i -> agrees(t, s, i), 1:length(t))
    end
end

@testset "stepper contract" begin
    t = DGGSSubtreeIds(S, 5, 0x0c4d9fffffffffff, 9)
    s = neighbor_stepper(t)
    border = Set(subtree_border_positions(t))

    for i in 1:length(t)
        v = step_neighbors(s, i)
        @test issorted(v)                           # ascending
        @test all(>(0), v)                          # off-tile absent, not zero
        @test allunique(v)
        # the rim is exactly the cells with a truncated neighborhood
        @test (length(v) < 6) == (i in border)
    end
end

@testset "table stepper matches the computed one" begin
    t = DGGSSubtreeIds(S, 5, 0x0c4d9fffffffffff, 9)
    computed = neighbor_stepper(t)
    table = neighbor_table(t)
    @test table isa TableNeighborStepper
    @test all(i -> collect(step_neighbors(table, i)) == collect(step_neighbors(computed, i)),
        1:length(t))

    # and it can be built over an explicitly chosen source
    generic = GenericNeighborStepper(t)
    @test all(i -> collect(step_neighbors(neighbor_table(t, generic), i)) ==
                   collect(step_neighbors(computed, i)), 1:1:length(t))
end

@testset "no allocation in the hot path" begin
    t = DGGSSubtreeIds(S, 5, 0x0c4d9fffffffffff, 9)
    s = neighbor_stepper(t)
    sweep(s, n) = (acc = 0; for i in 1:n, j in step_neighbors(s, i); acc += j end; acc)
    sweep(s, length(t))
    @test @allocated(sweep(s, length(t))) == 0
end
