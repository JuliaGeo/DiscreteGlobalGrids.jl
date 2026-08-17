# Compare closure-based and materialized sweeps with per-cell `neighbors`.

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

# Cover one-window subtrees and multi-window coverages.
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

# Include cell position, neighbour positions, and ring order in one result.
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
            # Traversal mode does not change position-ordered results.
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

        # Weight by ring slot to detect reordered or missing gathered values.
        data = collect(1.0:n)
        metric = (c, v, vals) ->
            3.0v + sum(i * vals[i] for i in eachindex(vals); init = 0.0)
        want_data = [3.0data[k] +
                     sum(i * data[cellposition(cv, x)]
                         for (i, x) in enumerate(neighbors(cv, cv[k])); init = 0.0)
                     for k in eachindex(cv)]
        @test mapneighbors(metric, cv, data; threaded = false) == want_data
        @test mapneighbors(metric, cv, data; threaded = true) == want_data

        # Concrete tuple results split into position-ordered vectors.
        pair(c, nbrs) = (cellposition(c), Float64(length(nbrs)))
        a, b = mapneighbors(pair, cv)
        @test a isa Vector{Int} && b isa Vector{Float64}
        @test a == collect(1:n)
        @test b == [Float64(length(neighbors(cv, c))) for c in cv]

        # `HaloTable` preserves ring order; sorted rows match `halo_table`.
        t = HaloTable(cv)
        rows = halo_table(cv)
        @test length(t) == n
        @test t.offsets[1] == 1 && t.offsets[end] == length(t.nbrs) + 1
        @test all(zip(1:n, neighbors(cv))) do (p, (c, nbrs))
            collect(t[p]) == [cellposition(h) for h in nbrs]
        end
        @test all(p -> sort(collect(t[p])) == rows[p], 1:n)
        @test rows == [neighbors(cv, p, 1) for p in 1:n]

        # Threaded and sequential builds have identical storage.
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

    # Rooted vectors preserve their grid; derived subsets do not retain a root.
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
    # Orders must visit every cell exactly once.
    bad = copy(perm); bad[2] = bad[1]
    @test_throws ArgumentError mapneighbors(probe, cv; order = bad)
    @test_throws ArgumentError mapneighbors(probe, cv; order = perm[1:(n-1)])
    @test_throws ArgumentError mapneighbors(probe, cv;
        order = vcat(perm[1:(n-1)], n + 1))
    @test_throws ArgumentError mapneighbors(probe, cv; order = :storage)
    # Data must share the collection's axis.
    @test_throws ArgumentError mapneighbors((c, v, vals) -> v, cv,
        collect(1.0:(n-1)))
end

# Keep mutable state in the functor while measuring sweep allocations.
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
