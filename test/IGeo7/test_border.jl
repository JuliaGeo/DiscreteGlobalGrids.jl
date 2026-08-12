# test/IGeo7/test_border.jl — the subtree rim for IGEO7: the native digit
# automaton (`IGeo7.border_descendants`, src/IGeo7/grid.jl) and the kernel
# wiring (`DGG.subtree_border(::IGEO7DGGS, ...)`, src/IGeo7/IGeo7Kernel.jl).
#
# Ground truth is the definition itself, computed the slow way: enumerate the
# *whole* subtree and keep the cells with an edge neighbor outside it
# (`brute_border` below). That reaches outside this file for its truth twice
# over — `cell_neighbors` is validated geometrically in test_neighbors.jl,
# against shared boundary corners rather than against any hierarchy — so the
# automaton is never checked against itself. It costs `7^depth` neighbor
# sweeps, which caps it at depth 6.
#
# Past that depth the reference is the *generic* kernel fallback
# (src/core/kernel.jl), reached with `invoke` since IGEO7 overrides it. It is
# an independent implementation: it decides membership by asking for actual
# neighbors at every level, where the native path reads the Z7 digits and never
# looks at the grid. The two share only the pruning premise — that a cell whose
# neighbors are all inside the subtree has no descendant touching the outside —
# and the brute-force depths are what pin that premise down.
#
# Fixtures deliberately cover the cases the digit rule could plausibly get
# wrong: pentagon roots in both hemispheres (the two `Z7_DELETED_DIGIT`
# values), a root whose rim crosses into another base cell, and roots at even
# *and* odd resolutions — the transition table is level-parity dependent, so a
# suite that only ever tested res-3 roots would pass on half a rule.
#
# Big sweeps accumulate a failure counter and assert once (the idiom of
# test_indexing.jl and test_neighbors.jl) so the suite's test count stays
# readable.

module IGeo7BorderTests

using Test
using Random

using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
using DiscreteGlobalGrids: Helpers, ISEA, IGeo7

const S = IGEO7DGGS()

"Descendants of `parent` at `res` with an edge neighbor outside the subtree,
by exhaustive enumeration of the subtree — `O(7^depth)`, the definition."
function brute_border(parent::UInt64, res::Integer)
    return UInt64[c for c in IGeo7.cell_to_children(parent, res)
                  if any(n -> !IGeo7.z7_is_descendant(n, parent),
        DGG.cell_neighbors(S, res, c))]
end

"The generic kernel fallback, which IGEO7 overrides — reached past the override."
generic_border(parent::UInt64, level::Integer, res::Integer) =
    invoke(DGG.subtree_border, Tuple{DGG.AbstractDGGS,Integer,Any,Integer},
        S, level, parent, res)

# Roots chosen for the failure modes the rule could have: hexagons at three
# different resolutions (the table is level-parity dependent), a hexagon whose
# rim reaches a neighboring base cell, and a pentagon from each half of
# `Z7_DELETED_DIGIT` (2 for bases 0-5, 5 for bases 6-11).
const ROOTS = (
    ("hexagon, res 3", IGeo7.lonlat_to_z7(10.0, 45.0, 3)),
    ("hexagon, res 3, southern", IGeo7.lonlat_to_z7(-120.0, -30.0, 3)),
    ("hexagon, res 2 (even parity)", IGeo7.lonlat_to_z7(-120.0, -30.0, 2)),
    ("hexagon, res 4", IGeo7.lonlat_to_z7(77.0, 28.0, 4)),
    ("hexagon, rim crosses a base cell", IGeo7.z7_from_string("00555")),
    ("pentagon, base 0 (deleted digit 2)",
        IGeo7.cell_to_children(IGeo7.res0_cells()[1], 3)[1]),
    ("pentagon, base 6 (deleted digit 5)",
        IGeo7.cell_to_children(IGeo7.res0_cells()[7], 3)[1]),
)

@testset "IGeo7 subtree border" begin

    @testset "matches the definition (brute force)" begin
        bad = 0
        for (_, root) in ROOTS, depth in 1:5
            level = IGeo7.get_resolution(root)
            native = IGeo7.border_descendants(root, level + depth)
            native == brute_border(root, level + depth) || (bad += 1)
        end
        @test bad == 0

        # One depth deeper for a single root: 117_649 descendants swept.
        root = ROOTS[1][2]
        @test IGeo7.border_descendants(root, 9) == brute_border(root, 9)
    end

    @testset "matches the generic kernel fallback past brute-force range" begin
        hexagon = ROOTS[1][2]
        for depth in 6:9
            @test IGeo7.border_descendants(hexagon, 3 + depth) ==
                  generic_border(hexagon, 3, 3 + depth)
        end
        pentagon = ROOTS[end][2]
        @test IGeo7.border_descendants(pentagon, 10) == generic_border(pentagon, 3, 10)
    end

    @testset "result shape" begin
        bad = 0
        for (_, root) in ROOTS, depth in 0:5
            level = IGeo7.get_resolution(root)
            rim = IGeo7.border_descendants(root, level + depth)
            issorted(rim) || (bad += 1)
            allunique(rim) || (bad += 1)
            all(c -> IGeo7.get_resolution(c) == level + depth, rim) || (bad += 1)
            all(c -> IGeo7.z7_is_descendant(c, root), rim) || (bad += 1)
            rim isa Vector{UInt64} || (bad += 1)
        end
        @test bad == 0

        # The subtree of depth 0 is the cell itself, and its whole neighborhood
        # lies outside it.
        for (_, root) in ROOTS
            @test IGeo7.border_descendants(root, IGeo7.get_resolution(root)) == [root]
        end
    end

    # The census `B(d) = 3 B(d-1) + (6 hexagon | 5 pentagon)` of grid.jl's
    # block comment, in closed form. Independent of the enumeration above:
    # `3^d`, not `7^d`, is the whole point of the operation.
    @testset "rim size" begin
        bad = 0
        for (_, root) in ROOTS, depth in 1:6
            level = IGeo7.get_resolution(root)
            expected = IGeo7.is_pentagon(root) ? (5 * (3^depth - 1)) ÷ 2 : 3^(depth + 1) - 3
            length(IGeo7.border_descendants(root, level + depth)) == expected || (bad += 1)
        end
        @test bad == 0

        # res 3 -> res 12: 59_046 rim cells out of 40_353_607 descendants, which
        # is the size the operation exists to avoid touching.
        @test length(IGeo7.border_descendants(ROOTS[1][2], 12)) == 3^10 - 3
        @test DGG.subtree_leaf_count(S, 3, ROOTS[1][2], 12) == 7^9
    end

    # Completeness at depth 9 is out of brute force's reach, but soundness is
    # not: every cell the automaton emits must really touch the outside.
    @testset "deep rim cells really touch the outside" begin
        root = ROOTS[1][2]
        rim = IGeo7.border_descendants(root, 12)
        rng = MersenneTwister(0x1660)
        bad = 0
        for c in rand(rng, rim, 500)
            any(n -> !IGeo7.z7_is_descendant(n, root), DGG.cell_neighbors(S, 12, c)) ||
                (bad += 1)
        end
        @test bad == 0
    end

    @testset "kernel wiring" begin
        bad = 0
        for (_, root) in ROOTS, depth in 0:4
            level = IGeo7.get_resolution(root)
            DGG.subtree_border(S, level, root, level + depth) ==
            IGeo7.border_descendants(root, level + depth) || (bad += 1)
        end
        @test bad == 0

        root = ROOTS[1][2]
        # `leaf_level < level` is the kernel's own argument mistake, so it is an
        # `ArgumentError` here as in `cell_descendants` — not the native error.
        @test_throws ArgumentError DGG.subtree_border(S, 3, root, 2)
        # Past that, Z7 validity keeps the native error type.
        @test_throws IGeo7.InvalidZ7Error IGeo7.border_descendants(root, 2)
        @test_throws IGeo7.InvalidZ7Error IGeo7.border_descendants(root, 20)
        @test_throws IGeo7.InvalidZ7Error DGG.subtree_border(S, 3, root, 20)
        @test_throws IGeo7.InvalidZ7Error IGeo7.border_descendants(UInt64(12) << 60, 1)
    end
end

end # module IGeo7BorderTests
