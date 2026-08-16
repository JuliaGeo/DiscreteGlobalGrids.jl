# Cross-system laws for the one-arg `neighbors` iterator and the positioned
# handles it mints. The oracle is the two-arg form: the iterator promises the
# SAME cells in the SAME order with positions attached, so every law here is
# an equality against per-cell calls that never see the window cursor.

module NeighborhoodTests

using Test
import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import Extents

using DiscreteGlobalGrids: levelgrid, ncells, cellindex, cellposition,
    neighbors, level, has_sorted_subtrees, cellindextype, PartialGrid,
    CellVector, CellLookup, MultiOrderCoverage, AuthalicSystem, Vertex, Edge,
    query, system, SubsetPositionedCell, cellid, Cells

const FB = DGG.Fallbacks

sysname(sys) = sys isa AuthalicSystem ?
               "Authalic($(nameof(typeof(parent(sys)))))" : string(nameof(typeof(sys)))

# ---------------------------------------------------------------------------
# Systems, and the two shapes each is swept over
#
# The subtree compresses to ONE window, so the cursor never moves; the
# coverage of a 1°×1° tile is DOZENS of disjoint windows, so the walk crosses
# a boundary every few cells — which is where a cursor that resolves against
# a stale window would pair a cell with a wrong position, and where a walk
# that advances wrongly would skip or duplicate a cell. The coverage level is
# per-system so every set stays in the hundreds-to-thousands.
# ---------------------------------------------------------------------------

const SWEEP = [
    (DGG.IGeo7System(), 1, 3, 8),
    (DGG.H3System(), 1, 3, 7),
    (DGG.HEALPixSystem(), 1, 4, 11),
    (DGG.A5System(), 1, 3, 11),
    (DGG.S2System(), 1, 4, 11),
    (DGG.ISEA4RSystem(), 1, 4, 11),
    (AuthalicSystem(DGG.IGeo7System()), 1, 3, 8),
]

@testset "the sweep covers every registered system" begin
    swept = Set(typeof(s) for (s, _, _, _) in SWEEP)
    for s in DGG.systems()
        @test typeof(s) in swept
    end
end

const TILE = Extents.Extent(X=(10.0, 11.0), Y=(46.0, 47.0))

rooted(sys, base, depth) =
    CellVector(PartialGrid(sys, cellindex(levelgrid(sys, base), 3), base + depth))

nwindows(cv) = FB.nwindows(FB.windows(cv))

@testset "$(sysname(sys))" for (sys, base, depth, covlvl) in SWEEP
    subtree = rooted(sys, base, depth)
    coverage = CellVector(query(sys, MultiOrderCoverage(TILE); level=covlvl))
    # The second shape must genuinely be the many-window case, or the cursor
    # laws below degenerate into a rerun of the first.
    @test nwindows(coverage) > 1

    @testset "$label" for (label, cv) in
                          ("one rooted subtree" => subtree,
                           "multi-window coverage" => coverage)
        for conn in (Vertex(), Edge())
            # The iterator is the two-arg calls, cell for cell and ring for
            # ring — a cursor that skips or duplicates at a window boundary,
            # or clips against the wrong window, breaks this at the first
            # affected cell.
            res = collect(neighbors(cv; connectivity = conn))
            @test length(res) == length(cv)
            @test all(eachindex(res)) do k
                c = cv[k]
                cellid(res[k][1]) == c &&
                    [cellid(h) for h in res[k][2]] ==
                    collect(neighbors(cv, c; connectivity = conn))
            end
            # Every minted handle re-resolves to the position it carries — a
            # minting site pairing a cell with any other position breaks this.
            @test all(neighbors(cv; connectivity = conn)) do (c, nbrs)
                cellposition(cv, cellid(c)) == cellposition(c) &&
                    all(h -> cellposition(cv, cellid(h)) == cellposition(h), nbrs)
            end
        end
    end

    # Handle `==` delegates to the cell, so agreement is checked position by
    # position — the part plain equality would wave through.
    @testset "the grid and lookup forms are the vector's" begin
        same(a, b) = cellid(a) == cellid(b) && cellposition(a) == cellposition(b)
        agree(x, y) = length(x) == length(y) &&
                      all(same(cx, cy) && length(nx) == length(ny) &&
                          all(splat(same), zip(nx, ny))
                          for ((cx, nx), (cy, ny)) in zip(x, y))
        pg = PartialGrid(sys, cellindex(levelgrid(sys, base), 3), base + depth)
        @test agree(collect(neighbors(pg)), collect(neighbors(subtree)))
        @test agree(collect(neighbors(CellLookup(subtree))), collect(neighbors(subtree)))
    end
end

@testset "an empty subset iterates to nothing" begin
    sys = DGG.IGeo7System()
    empty = CellVector(PartialGrid(sys, 3, cellindextype(sys)[]))
    @test isempty(collect(neighbors(empty)))
end

# ---------------------------------------------------------------------------
# The handle surface: one system, because it is system-generic code
# ---------------------------------------------------------------------------

@testset "a handle is its cell everywhere but at the fast path" begin
    sys = DGG.IGeo7System()
    cv = rooted(sys, 1, 3)
    h, nbrs = first(neighbors(cv))
    c = cellid(h)

    @test h == c && c == h && h == first(neighbors(cv))[1]
    @test hash(h) == hash(c)
    @test convert(typeof(c), h) === c
    @test sprint(show, h) == sprint(show, c)
    @test level(h) == level(c)
    @test isbitstype(typeof(h)) && sizeof(typeof(h)) == 16

    A = DD.DimArray(collect(1.0:length(cv)), (Cells(CellLookup(cv)),))
    # The fast path reads storage at the minted position...
    @test A[h] == parent(A)[cellposition(h)]
    A[h] = -1.0
    @test parent(A)[cellposition(h)] == -1.0
    # ...trusting it unconditionally: a wrong-but-in-range position reads that
    # slot, not the cell's. This is the documented contract, not a defect — an
    # indexing path that silently re-resolved the cell would hide the cost the
    # handle exists to remove.
    wrong = SubsetPositionedCell(c, 7)
    @test A[wrong] == parent(A)[7]
    # A bare cell keeps the resolved path.
    @test A[DD.At(c)] == A[h]
end

# The whole sweep stays off the heap: the yielded pair is a bits handle and a
# fixed-capacity ring of them, and the cursor is three integers of state. A
# `Vector` anywhere in the loop is invisible to every equality above — the
# cells would still be right, just not for free. Multi-window on purpose: the
# fallback search path must also be free.
@testset "the sweep allocates nothing" begin
    sys = DGG.IGeo7System()
    coverage = CellVector(query(sys, MultiOrderCoverage(TILE); level=8))
    @test nwindows(coverage) > 1
    sweeploop(it) = (n = 0; for (c, nbrs) in it
        n += cellposition(c)
        for h in nbrs
            n += cellposition(h)
        end
    end; n)
    it = neighbors(coverage)
    sweeploop(it)
    @test @allocated(sweeploop(it)) == 0
end

end # module NeighborhoodTests
