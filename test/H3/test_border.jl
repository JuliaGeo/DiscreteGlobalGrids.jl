# test/H3/test_border.jl — the subtree rim for H3: the digit automaton and its
# kernel wiring (`DGG.subtree_border(::H3DGGS, ...)`, src/H3/H3Kernel.jl).
#
# Ground truth is the definition itself, computed the slow way: enumerate the
# *whole* subtree and keep the cells with an edge neighbor outside it
# (`brute_border` below). That takes its truth from libh3 twice over —
# `cellToChildren` for the subtree and `gridDisk` for the neighbors, neither of
# which knows anything about the arc automaton — so the automaton is never
# checked against itself. It costs `7^depth` neighbor sweeps, which caps it
# around depth 6.
#
# Past that depth the reference is the *generic* kernel fallback
# (src/core/kernel.jl), reached with `invoke` since H3 overrides it. It is an
# independent implementation: it decides membership by asking libh3 for actual
# neighbors at every level, where the native path reads the index digits and
# never touches the grid. The two share only the pruning premise — that a cell
# whose neighbors all lie inside the subtree has no descendant touching the
# outside — and the brute-force depths are what pin that premise down.
#
# This matters more here than it does for IGEO7, whose arc transitions fall out
# of the lattice arithmetic in src/IGeo7/engine.jl. H3's table was fitted from
# observation and has no proof behind it (the STATUS paragraph in
# src/H3/H3Kernel.jl says so), so these comparisons are the whole of the
# evidence, not a regression net over a derivation.
#
# Fixtures deliberately cover what a fitted table could plausibly get wrong:
# roots at even *and* odd resolutions (the transitions are level-parity
# dependent, so a suite that only ever tested res-3 roots would pass on half a
# rule), pentagons at several resolutions plus every one of the 12 pentagon base
# cells, and a root whose rim provably reaches into another base cell.
#
# Big sweeps accumulate a failure counter and assert once (the idiom of
# test_neighbors.jl) so the suite's test count stays readable.

module H3BorderTests

using Test
using Random

using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
const H3N = DGG.H3.H3Native

const S = H3DGGS()

"Descendants of `root` at `res` with an edge neighbor outside the subtree, by
exhaustive enumeration of the subtree — `O(7^depth)`, the definition."
function brute_border(root::UInt64, res::Integer)
    level = H3N.get_resolution(root)
    return UInt64[c for c in H3N.cell_to_children(root, res)
                  if any(n -> H3N.cell_to_parent(n, level) != root,
        DGG.cell_neighbors(S, res, c))]
end

"The generic kernel fallback, which H3 overrides — reached past the override."
generic_border(root::UInt64, level::Integer, res::Integer) =
    invoke(DGG.subtree_border, Tuple{DGG.AbstractDGGS,Integer,Any,Integer},
        S, level, root, res)

border(root::UInt64, res::Integer) =
    DGG.subtree_border(S, H3N.get_resolution(root), root, res)

# A res-1 child of the base cell containing (0, 0): small enough to brute-force
# and adjacent to a base-cell boundary, which `crossing_root` below asserts
# rather than assumes.
const CROSSING = H3N.cell_to_children(H3N.lonlat_to_cell(0.0, 0.0, 0), 1)[4]

const ROOTS = (
    ("hexagon, res 1 (odd)", H3N.lonlat_to_cell(10.0, 45.0, 1)),
    ("hexagon, res 2 (even)", H3N.lonlat_to_cell(10.0, 45.0, 2)),
    ("hexagon, res 3 (odd)", H3N.lonlat_to_cell(10.0, 45.0, 3)),
    ("hexagon, res 4 (even)", H3N.lonlat_to_cell(77.0, 28.0, 4)),
    ("hexagon, res 5 (odd), southern", H3N.lonlat_to_cell(-120.0, -30.0, 5)),
    ("hexagon, res 0 base cell", H3N.lonlat_to_cell(-58.0, -15.0, 0)),
    ("hexagon, rim crosses a base cell", CROSSING),
    ("pentagon, res 0 base cell", H3N.get_pentagons(0)[1]),
    ("pentagon, res 1 (odd)", H3N.get_pentagons(1)[1]),
    ("pentagon, res 2 (even)", H3N.get_pentagons(2)[1]),
    ("pentagon, res 3 (odd)", H3N.get_pentagons(3)[1]),
)

@testset "H3 subtree border" begin

    @testset "matches the definition (brute force)" begin
        bad = 0
        for (_, root) in ROOTS, depth in 0:4
            level = H3N.get_resolution(root)
            border(root, level + depth) == sort(brute_border(root, level + depth)) ||
                (bad += 1)
        end
        @test bad == 0

        # Every pentagon base cell, not just the first: the automaton assumes
        # all 12 delete the same child (digit 1, the K axis), where IGEO7 needs
        # a two-valued table. If that were wrong for even one of them, the rim
        # would contain cells libh3 does not have.
        bad = 0
        for pentagon in H3N.get_pentagons(0), depth in 1:3
            border(pentagon, depth) == sort(brute_border(pentagon, depth)) || (bad += 1)
        end
        @test bad == 0

        # Two depths deeper on a smaller selection: 117_649 descendants swept
        # for the depth-6 case.
        bad = 0
        for root in (ROOTS[3][2], ROOTS[4][2], ROOTS[10][2]), depth in 5:6
            level = H3N.get_resolution(root)
            border(root, level + depth) == sort(brute_border(root, level + depth)) ||
                (bad += 1)
        end
        @test bad == 0
    end

    @testset "matches the generic kernel fallback past brute-force range" begin
        for (name, root) in (ROOTS[3], ROOTS[4])
            level = H3N.get_resolution(root)
            for depth in 7:9
                @test border(root, level + depth) ==
                      generic_border(root, level, level + depth)
            end
        end
        pentagon = ROOTS[10][2]
        @test border(pentagon, 2 + 8) == generic_border(pentagon, 2, 2 + 8)
    end

    @testset "result shape" begin
        bad = 0
        for (_, root) in ROOTS, depth in 0:4
            level = H3N.get_resolution(root)
            rim = border(root, level + depth)
            issorted(rim) || (bad += 1)                       # [contract], not sorted
            allunique(rim) || (bad += 1)
            rim isa Vector{UInt64} || (bad += 1)
            all(H3N.is_valid_cell, rim) || (bad += 1)
            all(c -> H3N.get_resolution(c) == level + depth, rim) || (bad += 1)
            all(c -> H3N.cell_to_parent(c, level) == root, rim) || (bad += 1)
            # Only the root can be a pentagon on a rim path, which is what lets
            # `_h3_fill_border!` stop retesting one level down.
            depth >= 1 && any(H3N.is_pentagon, rim) && (bad += 1)
        end
        @test bad == 0

        # The subtree of depth 0 is the cell itself, and its whole neighborhood
        # lies outside it.
        for (_, root) in ROOTS
            @test border(root, H3N.get_resolution(root)) == [root]
        end
    end

    # The fixture named "rim crosses a base cell" has to actually do that, or
    # the case it is there for is untested.
    @testset "crossing root really crosses" begin
        rim = border(CROSSING, 4)
        home = H3N.get_base_cell(CROSSING)
        outside = [n for c in rim for n in DGG.cell_neighbors(S, 4, c)
                   if H3N.get_base_cell(n) != home]
        @test !isempty(outside)
        @test border(CROSSING, 4) == sort(brute_border(CROSSING, 4))
    end

    # The census `B(d) = 3 B(d-1) + (6 hexagon | 5 pentagon)` of H3Kernel.jl's
    # block comment, in closed form. Independent of the enumeration above:
    # `3^d`, not `7^d`, is the whole point of the operation.
    @testset "rim size" begin
        bad = 0
        for (_, root) in ROOTS, depth in 1:6
            level = H3N.get_resolution(root)
            level + depth <= H3N.MAX_RESOLUTION || continue
            expected = H3N.is_pentagon(root) ? (5 * (3^depth - 1)) ÷ 2 : 3^(depth + 1) - 3
            length(border(root, level + depth)) == expected || (bad += 1)
        end
        @test bad == 0

        # res 3 -> res 15 (the deepest H3 goes): 1_594_320 rim cells out of
        # 7^12 = 13_841_287_201 descendants, which is the size the operation
        # exists to avoid touching. It is also where the last digit slot sits at
        # bit 0, so `_h3_digit_shift` is exercised at both ends.
        root = ROOTS[3][2]
        @test length(border(root, 15)) == 3^13 - 3
        @test DGG.subtree_leaf_count(S, 3, root, 15) == 7^12
    end

    # Completeness at depth 12 is out of brute force's reach, but soundness is
    # not: every cell the automaton emits must really touch the outside.
    @testset "deep rim cells really touch the outside" begin
        root = ROOTS[3][2]
        rim = border(root, 15)
        rng = MersenneTwister(0xf3)
        bad = 0
        for c in rand(rng, rim, 500)
            any(n -> H3N.cell_to_parent(n, 3) != root, DGG.cell_neighbors(S, 15, c)) ||
                (bad += 1)
        end
        @test bad == 0
    end

    @testset "kernel wiring and errors" begin
        root = ROOTS[3][2]                                   # a res-3 hexagon
        @test DGG.subtree_border(S, 3, root, 6) == border(root, 6)
        # An `Integer` id and a non-`Int` level both reach the same method.
        @test DGG.subtree_border(S, Int32(3), Int128(root), UInt8(6)) == border(root, 6)

        # `leaf_level < level` is the kernel's own argument mistake, so it is an
        # `ArgumentError` here as in `cell_descendants`.
        @test_throws ArgumentError DGG.subtree_border(S, 3, root, 2)
        # H3 has no error type of its own; `descendant_range` reports a bad
        # resolution and a level/index mismatch as `ArgumentError` too.
        @test_throws ArgumentError DGG.subtree_border(S, 3, root, 16)
        @test_throws ArgumentError DGG.subtree_border(S, 2, root, 6)
        @test_throws ArgumentError DGG.subtree_border(S, 4, root, 6)

        # Id validation is this method's own work, and has to be: libh3
        # validates neither `cellToChildren` nor `cellToChildrenSize`, so
        # `cell_descendants` — the route the generic fallback takes for its
        # degenerate case, to inherit each system's id checking — hands back a
        # confident subtree of nonexistent cells for a malformed index. The
        # assertions on `cell_descendants` below are not aspirational; they pin
        # the behavior that makes the explicit `isValidCell` necessary.
        bogus = root & ~(UInt64(7) << 30)                    # a padding slot cleared
        @test !H3N.is_valid_cell(bogus)
        @test H3N.get_resolution(bogus) == 3
        @test DGG.cell_descendants(S, 3, bogus, 3) == [bogus]
        @test_throws ArgumentError DGG.subtree_border(S, 3, bogus, 5)
        @test_throws ArgumentError DGG.subtree_border(S, 3, bogus, 3)

        # The other shapes `isValidCell` rejects, all of which read as a
        # plausible res-3 index if you only look at the resolution field.
        digit7 = root | (UInt64(7) << 36)                    # digit 7 in an active slot
        nobase = (root & ~(UInt64(0x7f) << 45)) | (UInt64(125) << 45)
        @test_throws ArgumentError DGG.subtree_border(S, 3, digit7, 5)
        @test_throws ArgumentError DGG.subtree_border(S, 3, nobase, 5)

        # A leading K digit under a pentagon base cell: the child H3 deletes.
        pentagon = H3N.get_pentagons(0)[1]
        res1 = (pentagon & ~(UInt64(0xf) << 52)) | (UInt64(1) << 52)
        deleted = (res1 & ~(UInt64(7) << 42)) | (UInt64(1) << 42)
        @test !H3N.is_valid_cell(deleted)
        @test_throws ArgumentError DGG.subtree_border(S, 1, deleted, 4)
    end
end

end # module H3BorderTests
