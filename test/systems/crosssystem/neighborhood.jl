# Compare the positioned `neighbors` iterator with per-cell neighbour calls.

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
# Systems and subset shapes covered by the sweep.
# ---------------------------------------------------------------------------

const SWEEP = [
    (DGG.IGeo7System(), 1, 3, 8),
    (DGG.H3System(), 1, 3, 7),
    (DGG.HEALPixSystem(), 1, 4, 11),
    (DGG.A5System(), 1, 3, 11),
    (DGG.S2System(), 1, 4, 11),
    (DGG.ISEA4RSystem(), 1, 4, 11),
    (AuthalicSystem(DGG.IGeo7System()), 1, 3, 8),
    # The literature families. Levels are chosen per aperture rather than
    # shared, so that each row's rooted subtree and tile coverage come out the
    # same size as the rows above rather than the same level number: the
    # coverages here run 213-873 cells over 19-49 windows, against 779-1035
    # over 48-80 for the six above. A shared level number would make the
    # aperture-9 rows two orders of magnitude larger than the aperture-3 one.
    (DGG.ISEA3HSystem(), 1, 5, 13),
    (DGG.ISEA4HSystem(), 1, 4, 10),
    (DGG.ISEA4TSystem(), 1, 4, 10),
    (DGG.RHEALPixSystem(), 1, 3, 7),
    (DGG.AusPIXSystem(), 1, 3, 7),
    (DGG.IVEA4RSystem(), 1, 4, 10),
    (DGG.IVEA9RSystem(), 1, 3, 7),
    (DGG.RTEA4RSystem(), 1, 4, 10),
    (DGG.RTEA9RSystem(), 1, 3, 7),
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
    # The coverage must exercise window transitions.
    @test nwindows(coverage) > 1

    @testset "$label" for (label, cv) in
                          ("one rooted subtree" => subtree,
                           "multi-window coverage" => coverage)
        for conn in (Vertex(), Edge())
            # Cells and ring order match the per-cell form.
            res = collect(neighbors(cv; connectivity = conn))
            @test length(res) == length(cv)
            @test all(eachindex(res)) do k
                c = cv[k]
                cellid(res[k][1]) == c &&
                    [cellid(h) for h in res[k][2]] ==
                    collect(neighbors(cv, c; connectivity = conn))
            end
            # Every handle carries the position resolved from its cell.
            @test all(neighbors(cv; connectivity = conn)) do (c, nbrs)
                cellposition(cv, cellid(c)) == cellposition(c) &&
                    all(h -> cellposition(cv, cellid(h)) == cellposition(h), nbrs)
            end
        end
    end

    # Compare stored positions explicitly because handle equality uses cell ids.
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
# System-generic handle behavior.
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
    # Handles read and write their stored position.
    @test A[h] == parent(A)[cellposition(h)]
    A[h] = -1.0
    @test parent(A)[cellposition(h)] == -1.0
    # An in-range position is used without resolving the cell id.
    wrong = SubsetPositionedCell(c, 7)
    @test A[wrong] == parent(A)[7]
    # Bare cells use the selector path.
    @test A[DD.At(c)] == A[h]
end

# The multi-window sweep allocates nothing after warmup.
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
