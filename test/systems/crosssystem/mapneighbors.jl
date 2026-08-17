# Cross-system laws for the closure form — `mapneighbors`, `foreachneighbors`
# — and the materialized sweep, `HaloTable`. The oracle is again the two-arg
# `neighbors`: the frame promises the same rings, delivered to `f` instead of
# yielded, with outputs landing by position no matter how the loop ran.

module MapNeighborsTests

using Test
import DiscreteGlobalGrids as DGG
import Extents
using Random: shuffle, Xoshiro

using DiscreteGlobalGrids: levelgrid, cellindex, cellposition, neighbors,
    mapneighbors, foreachneighbors, StorageOrder, HaloTable, halo_table,
    PartialGrid, CellVector, CellLookup, MultiOrderCoverage, AuthalicSystem,
    Vertex, Edge, query, system, cellid, level

const FB = DGG.Fallbacks

sysname(sys) = sys isa AuthalicSystem ?
               "Authalic($(nameof(typeof(parent(sys)))))" : string(nameof(typeof(sys)))

# The same sweep as `neighborhood.jl`: one shape whose cursor never moves and
# one that crosses a window boundary every few cells.
const SWEEP = [
    (DGG.IGeo7System(), 1, 3, 8),
    (DGG.H3System(), 1, 3, 7),
    (DGG.HEALPixSystem(), 1, 4, 11),
    (DGG.A5System(), 1, 3, 11),
    (DGG.S2System(), 1, 4, 11),
    (DGG.ISEA4RSystem(), 1, 4, 11),
    (AuthalicSystem(DGG.IGeo7System()), 1, 3, 8),
]

const TILE = Extents.Extent(X=(10.0, 11.0), Y=(46.0, 47.0))

rooted_pg(sys, base, depth) =
    PartialGrid(sys, cellindex(levelgrid(sys, base), 3), base + depth)

# A probe that reads everything the frame hands over — the cell's position,
# the ring's positions, and the ring's ORDER — so a frame that stores at the
# wrong slot, clips against the wrong window, or reorders a ring cannot match
# the oracle built from per-cell calls.
probe(c, nbrs) =
    cellposition(c) * 31 + sum(i * cellposition(h) for (i, h) in enumerate(nbrs);
        init = 0)

naive(cv; connectivity = Vertex()) =
    [k * 31 + sum(i * cellposition(cv, x)
                  for (i, x) in enumerate(neighbors(cv, cv[k]; connectivity));
          init = 0)
     for k in eachindex(cv)]

@testset "$(sysname(sys))" for (sys, base, depth, covlvl) in SWEEP
    subtree = CellVector(rooted_pg(sys, base, depth))
    coverage = CellVector(query(sys, MultiOrderCoverage(TILE); level=covlvl))

    @testset "$label" for (label, cv) in
                          ("one rooted subtree" => subtree,
                           "multi-window coverage" => coverage)
        n = length(cv)
        for conn in (Vertex(), Edge())
            want = naive(cv; connectivity = conn)
            # The frame is the per-cell calls: sequential, threaded, and under
            # a permutation the outputs land by position, bit for bit.
            @test mapneighbors(probe, cv; threaded = false,
                connectivity = conn) == want
            @test mapneighbors(probe, cv; threaded = true,
                connectivity = conn) == want
            perm = shuffle(Xoshiro(42), 1:n)
            @test mapneighbors(probe, cv; order = perm, threaded = false,
                connectivity = conn) == want
            @test mapneighbors(probe, cv; order = perm, threaded = true,
                connectivity = conn) == want
        end

        # The data form: `f` gets its own sample and the ring's, gathered by
        # the frame, aligned with the ring — the weight by slot index is what
        # catches a gather that reorders or drops one.
        data = collect(1.0:n)
        metric = (c, v, vals) ->
            3.0v + sum(i * vals[i] for i in eachindex(vals); init = 0.0)
        want_data = [3.0data[k] +
                     sum(i * data[cellposition(cv, x)]
                         for (i, x) in enumerate(neighbors(cv, cv[k])); init = 0.0)
                     for k in eachindex(cv)]
        @test mapneighbors(metric, cv, data; threaded = false) == want_data
        @test mapneighbors(metric, cv, data; threaded = true) == want_data

        # A concrete tuple return is the multi-output frame: one vector per
        # component, each in position order.
        pair(c, nbrs) = (cellposition(c), Float64(length(nbrs)))
        a, b = mapneighbors(pair, cv)
        @test a isa Vector{Int} && b isa Vector{Float64}
        @test a == collect(1:n)
        @test b == [Float64(length(neighbors(cv, c))) for c in cv]

        # The materialized sweep is the iterator, row for row and slot for
        # slot — ring order, not ascending — and its sorted rows are
        # `halo_table`'s.
        t = HaloTable(cv)
        rows = halo_table(cv)
        @test length(t) == n
        @test t.offsets[1] == 1 && t.offsets[end] == length(t.nbrs) + 1
        @test all(zip(1:n, neighbors(cv))) do (p, (c, nbrs))
            collect(t[p]) == [cellposition(h) for h in nbrs]
        end
        @test all(p -> sort(collect(t[p])) == rows[p], 1:n)
        @test rows == [neighbors(cv, p, 1) for p in 1:n]

        # `t` and `rows` above are the default THREADED builds, so the laws
        # already ran against the chunked sweep. What remains is that the
        # stitch adds nothing and loses nothing: the sequential build's
        # arrays, field for field — a seam that shifted an offset, or a chunk
        # landing out of range order, breaks equality at the boundary row.
        s = HaloTable(cv; threaded = false)
        @test t.offsets == s.offsets && t.nbrs == s.nbrs
        @test halo_table(cv; threaded = false) == rows
    end

    @testset "the grid and lookup forms are the vector's" begin
        pg = rooted_pg(sys, base, depth)
        want = naive(subtree)
        @test mapneighbors(probe, pg; threaded = false) == want
        @test mapneighbors(probe, CellLookup(subtree); threaded = false) == want
        @test HaloTable(pg) == HaloTable(subtree)
        @test HaloTable(CellLookup(subtree)) == HaloTable(subtree)
    end

    # The root survives the vector round trip, so `halo_table` on a
    # subtree-built vector reaches the interior/rim fast path again; a
    # derived subset still reports an unrooted grid, because its windows
    # remember no ancestor.
    @testset "the root survives the round trip" begin
        pg = rooted_pg(sys, base, depth)
        cv = CellVector(pg)
        back = PartialGrid(cv)
        @test back === pg
        @test FB._is_rooted(back)
        @test !FB._is_rooted(PartialGrid(cv[1:(length(cv)-1)]))
        @test halo_table(cv) == halo_table(pg)
    end
end

@testset "a visit order is followed, once per cell, or refused" begin
    sys = DGG.IGeo7System()
    cv = CellVector(rooted_pg(sys, 1, 3))
    n = length(cv)
    perm = shuffle(Xoshiro(7), 1:n)
    visits = Int[]
    foreachneighbors((c, nbrs) -> push!(visits, cellposition(c)), cv; order = perm)
    @test visits == perm
    # A duplicate would overwrite one cell's output and leave another's
    # undefined, so every non-permutation is refused up front.
    bad = copy(perm); bad[2] = bad[1]
    @test_throws ArgumentError mapneighbors(probe, cv; order = bad)
    @test_throws ArgumentError mapneighbors(probe, cv; order = perm[1:(n-1)])
    @test_throws ArgumentError mapneighbors(probe, cv;
        order = vcat(perm[1:(n-1)], n + 1))
    @test_throws ArgumentError mapneighbors(probe, cv; order = :storage)
    # ...and data that is not laid out against the collection likewise.
    @test_throws ArgumentError mapneighbors((c, v, vals) -> v, cv,
        collect(1.0:(n-1)))
end

# The frame's own cost is three integers of cursor state and the closure it
# was handed: a sequential sweep touches the heap not at all. The functor
# keeps the accumulator out of the closure so the measurement is the frame's.
struct PositionSum <: Function
    acc::Base.RefValue{Int}
end
(s::PositionSum)(c, nbrs) = (s.acc[] += cellposition(c) + length(nbrs); nothing)

@testset "the sequential sweep allocates nothing" begin
    sys = DGG.IGeo7System()
    coverage = CellVector(query(sys, MultiOrderCoverage(TILE); level=8))
    s = PositionSum(Ref(0))
    foreachneighbors(s, coverage)
    @test @allocated(foreachneighbors(s, coverage)) == 0
end

end # module MapNeighborsTests
