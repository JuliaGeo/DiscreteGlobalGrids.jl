# Compare closure-based and materialized sweeps with per-cell `neighbors`.

module MapNeighborsTests

using Test
import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import Extents
using Random: shuffle, Xoshiro

using DiscreteGlobalGrids: levelgrid, cellindex, localindex, neighbors, ring,
    mapneighbors, foreachneighbors, StorageOrder, adjacency, subtree,
    PartialGrid, CellVector, CellLookup, Cells, MultiOrderCoverage,
    AuthalicSystem, Vertex, Edge, query, system, cellid, level, neighborcount,
    Neighbors, Values, NeighborSlices, Disc, Ring, Cell

include(joinpath(@__DIR__, "..", "..", "helpers.jl"))
using .DGGTestHelpers: syslabel, sweepcovers, basesystem, ishexwalk,
    FORCED_BOUNDS_CHECKS

const FB = DGG.Fallbacks
const EN = DGG.Engine

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

@testset "the sweep covers every registered system" begin
    sweepcovers(SWEEP)
end

const TILE = Extents.Extent(X=(10.0, 11.0), Y=(46.0, 47.0))

rooted_pg(sys, base, depth) =
    subtree(sys, cellindex(levelgrid(sys, base), 3), base + depth)

# Include cell index, neighbour indices, and ring order in one result.
probe(c, nbrs) =
    localindex(c) * 31 + sum(i * localindex(h) for (i, h) in enumerate(nbrs);
        init = 0)

naive(cv; connectivity = Vertex()) =
    [k * 31 + sum(i * localindex(cv, x)
                  for (i, x) in enumerate(neighbors(cv, cv[k]; connectivity));
          init = 0)
     for k in eachindex(cv)]

@testset "$(syslabel(sys))" for (sys, base, depth, covlvl) in SWEEP
    sub = CellVector(rooted_pg(sys, base, depth))
    coverage = CellVector(query(sys, MultiOrderCoverage(TILE); level=covlvl))

    @testset "$label" for (label, cv) in
                          ("one rooted subtree" => sub,
                           "multi-window coverage" => coverage)
        n = length(cv)
        for conn in (Vertex(), Edge())
            want = naive(cv; connectivity = conn)
            # Traversal mode does not change index-ordered results.
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
                     sum(i * data[localindex(cv, x)]
                         for (i, x) in enumerate(neighbors(cv, cv[k])); init = 0.0)
                     for k in eachindex(cv)]
        @test mapneighbors(metric, cv, data; threaded = false) == want_data
        @test mapneighbors(metric, cv, data; threaded = true) == want_data

        # Concrete tuple results split into index-ordered vectors.
        pair(c, nbrs) = (localindex(c), Float64(length(nbrs)))
        a, b = mapneighbors(pair, cv)
        @test a isa Vector{Int} && b isa Vector{Float64}
        @test a == collect(1:n)
        @test b == [Float64(length(neighbors(cv, c))) for c in cv]

        # `adjacency` is the same rows in CSR form, in the same ring order.
        t = adjacency(cv)
        @test length(t) == n
        @test t.offsets[1] == 1 && t.offsets[end] == length(t.indices) + 1
        @test all(zip(1:n, neighbors(cv))) do (p, (c, nbrs))
            collect(t[p]) == [localindex(h) for h in nbrs]
        end
        @test all(p -> collect(t[p]) == neighbors(cv, p, 1), 1:n)

        # Threaded and sequential builds have identical storage.
        @test adjacency(cv; threaded = false) == t
    end

    @testset "the grid and lookup forms are the vector's" begin
        pg = rooted_pg(sys, base, depth)
        want = naive(sub)
        @test mapneighbors(probe, pg; threaded = false) == want
        @test mapneighbors(probe, CellLookup(sub); threaded = false) == want
        @test adjacency(pg) == adjacency(sub)
        @test adjacency(CellLookup(sub)) == adjacency(sub)
    end

    # Rooted vectors preserve their grid; derived subsets do not retain a root.
    @testset "the root survives the round trip" begin
        pg = rooted_pg(sys, base, depth)
        cv = CellVector(pg)
        back = PartialGrid(cv)
        @test back === pg
        @test EN._is_rooted(back)
        @test !EN._is_rooted(PartialGrid(cv[1:(length(cv)-1)]))
        @test adjacency(cv) == adjacency(pg)
    end
end

@testset "a visit order is followed, once per cell, or refused" begin
    sys = DGG.IGeo7System()
    cv = CellVector(rooted_pg(sys, 1, 3))
    n = length(cv)
    perm = shuffle(Xoshiro(7), 1:n)
    visits = Int[]
    foreachneighbors((c, nbrs) -> push!(visits, localindex(c)), cv; order = perm)
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

# A threaded sweep fans the callback out over tasks; one failure must not come
# back as one exception per task.
@testset "a threaded callback failure is reported once, with its cell" begin
    sys = DGG.IGeo7System()
    cv = CellVector(rooted_pg(sys, 1, 3))
    bad = 5
    boom(c, nbrs) = localindex(c) == bad ? error("callback said no") : 1.0

    err = try
        mapneighbors(boom, cv; threaded = true)
        nothing
    catch e
        e
    end
    @test err isa DGG.NeighborCallbackError
    msg = sprint(showerror, err)
    @test occursin(sprint(show, cv[bad]), msg)
    @test occursin("subset index $bad", msg)
    @test occursin("threaded = false", msg)
    # The callback's own exception is the cause, not something swallowed.
    @test occursin("callback said no", msg)
    # And the sequential path still raises exactly what the callback threw.
    @test_throws "callback said no" mapneighbors(boom, cv; threaded = false)
end

# Keep mutable state in the functor while measuring sweep allocations.
struct IndexSum <: Function
    acc::Base.RefValue{Int}
end
(s::IndexSum)(c, nbrs) = (s.acc[] += localindex(c) + length(nbrs); nothing)

@testset "the sequential sweep allocates nothing" begin
    sys = DGG.IGeo7System()
    coverage = CellVector(query(sys, MultiOrderCoverage(TILE); level=8))
    s = IndexSum(Ref(0))
    foreachneighbors(s, coverage)
    @test @allocated(foreachneighbors(s, coverage)) == 0
    # Spelling the default out reaches the same path: `Disc(1)` is the one-ring
    # call the sweep made before the selector existed, keyword and all.
    foreachneighbors(s, coverage; neighborhood = Disc(1))
    @test @allocated(foreachneighbors(s, coverage; neighborhood = Disc(1))) == 0
end

# A wider sweep is free wherever the system carries its winding: the static
# walk builds one stack buffer per shell, and H3 writes libh3's shells into
# stack buffers of its own. Selected by the `winding` trait, never by name.
@testset "a wider sweep allocates nothing where the winding is carried" begin
    carried(sys) = DGG.winding(levelgrid(sys, 1), Vertex()) isa
                   Union{DGG.CounterClockwise,DGG.Clockwise}
    for sys in filter(carried, DGG.systems())
        cv = CellVector(rooted_pg(sys, 1, 3))
        s = IndexSum(Ref(0))
        @testset "$(syslabel(sys)) $conn $nb" for conn in (Vertex(), Edge()),
                nb in (Disc(2), Ring(2), Disc(3), Ring(3))
            foreachneighbors(s, cv; neighborhood = nb, connectivity = conn)
            @test @allocated(foreachneighbors(s, cv; neighborhood = nb,
                connectivity = conn)) == 0 skip = VERSION < v"1.12" ||
                                                  FORCED_BOUNDS_CHECKS
        end
    end
    @test !isempty(filter(carried, DGG.systems()))
end

# The whole-array entry points resolve the cell dimension and hand back the
# caller's own wrapper. One system suffices: the code is system-generic and
# the numbers are pinned to the CellVector layer's.
@testset "any dimarray in, the same lookups out" begin
    sys = DGG.IGeo7System()
    cv = CellVector(rooted_pg(sys, 1, 3))
    n = length(cv)
    data = collect(1.0:n)
    A = DD.DimArray(copy(data), (Cells(CellLookup(cv)),))
    # Put the cell dimension second to exercise automatic detection; each
    # time slice carries its own numbers.
    cubedata = [j * data[k] + 0.1j for j in 1:3, k in 1:n]
    cube = DD.DimArray(copy(cubedata), (DD.Dim{:time}(1:3), Cells(CellLookup(cv))))

    # The default passes indexed handles: `probe` has no three-argument
    # method, so a default flipped to Values() errors instead of matching
    # the CellVector layer's numbers.
    out = mapneighbors(probe, A; threaded = false)
    @test out isa DD.AbstractDimArray
    @test parent(out) == mapneighbors(probe, cv; threaded = false)
    @test DD.dims(out) === DD.dims(A)
    @test parent(DD.lookup(out, 1)) === cv

    # An N-D handle selects a view, while the output retains the cell lookup.
    colsum = mapneighbors((c, nbrs) -> sum(cube[c]), cube; threaded = false)
    @test size(colsum) == (n,)
    @test parent(DD.lookup(colsum, 1)) === cv
    @test parent(colsum) == vec(sum(cubedata; dims = 1))
    h = DGG.SubsetIndexedCell(cv[5], 5)
    @test cube[h] == cubedata[:, 5]
    @test map(DD.name, DD.dims(cube[h])) == (:time,)
    @test view(cube, h) == cubedata[:, 5]

    # Values() matches the bare data call, and each time slice equals its
    # own 1-D call — a slicer along the wrong dim or with interacting slices
    # fails here.
    metric = (c, v, vals) ->
        3.0v + sum(i * vals[i] for i in eachindex(vals); init = 0.0)
    outV = mapneighbors(metric, A; pass = Values(), threaded = false)
    @test parent(outV) == mapneighbors(metric, cv, data; threaded = false)
    @test DD.dims(outV) === DD.dims(A)
    coutV = mapneighbors(metric, cube; pass = Values(), threaded = false)
    @test DD.dims(coutV) === DD.dims(cube)
    @test all(1:3) do j
        parent(coutV)[j, :] ==
            mapneighbors(metric, cv, cubedata[j, :]; threaded = false)
    end

    # NeighborSlices() passes one-dimension-smaller views — the ndims guard
    # poisons the value on any other shape — and refuses a 1-D array by
    # naming Values().
    sliced = (c, s, ns) ->
        ndims(s) == 1 && all(x -> ndims(x) == 1, ns) ?
        sum(s) + sum(sum, ns; init = 0.0) : NaN
    outS = mapneighbors(sliced, cube; pass = NeighborSlices(), threaded = false)
    @test size(outS) == (n,)
    @test parent(DD.lookup(outS, 1)) === cv
    rings = collect(neighbors(cv))
    @test parent(outS) ≈ [sum(cubedata[:, k]) +
        sum(sum(cubedata[:, localindex(h)]) for h in rings[k][2]; init = 0.0)
        for k in 1:n]
    @test_throws "use Values()" mapneighbors(sliced, A; pass = NeighborSlices())

    # Explicit dimension selectors agree with automatic detection.
    @test parent(mapneighbors(metric, cube; spatialdim = Cells, pass = Values(),
        threaded = false)) == parent(coutV)
    @test parent(mapneighbors(metric, cube; spatialdim = :Cells, pass = Values(),
        threaded = false)) == parent(coutV)

    # Invalid dimensions and pass modes report the failing condition.
    @test_throws "carries a cell lookup" mapneighbors(probe,
        DD.DimArray(collect(1.0:4), (DD.X(1:4),)))
    @test_throws "no dimension matching" mapneighbors(probe, cube;
        spatialdim = DD.Ti)
    @test_throws "not a cell lookup" mapneighbors(probe, cube;
        spatialdim = :time)
    @test_throws "pass must be" mapneighbors(probe, cube; pass = :values)

    # Concrete tuple results produce one wrapped array per component.
    two = (c, v, vals) -> (v, Float64(length(vals)))
    a, b = mapneighbors(two, cube; pass = Values(), threaded = false)
    @test a isa DD.AbstractDimArray && b isa DD.AbstractDimArray
    @test parent(a) == cubedata
    @test DD.dims(b) === DD.dims(cube)

    # foreachneighbors and neighbors use the same cell-dimension resolution.
    # Sums are ≈ where the visit order differs from `sum`'s.
    accH = Ref(0)
    foreachneighbors((c, nbrs) -> (accH[] += length(nbrs)), cube)
    @test accH[] == sum(length(r[2]) for r in rings)
    accV = Ref(0.0)
    foreachneighbors((c, v, vals) -> (accV[] += v), cube; pass = Values())
    @test accV[] ≈ sum(cubedata)
    res = collect(neighbors(cube))
    @test length(res) == n
    @test cellid(res[1][1]) == cv[1]

    # The positional `dims` form names the cell dimension the way
    # `DimensionalData.dims(A, dims)` reads it, and agrees with discovery.
    for dims in (Cells, :Cells, (Cells,))
        @test collect(neighbors(cube, dims)) == res
        @test parent(mapneighbors(metric, cube, dims; pass = Values(),
            threaded = false)) == parent(coutV)
        @test adjacency(cube, dims) == adjacency(cv)
    end
    @test adjacency(cube) == adjacency(cv)
    @test adjacency(A) == adjacency(cv)
    @test adjacency(cube; threaded = false) == adjacency(cv)
    acc = Ref(0)
    foreachneighbors((c, nbrs) -> (acc[] += length(nbrs)), cube, Cells)
    @test acc[] == accH[]
    @test_throws "no dimension matching" neighbors(cube, DD.Ti)
    @test_throws "not a cell lookup" adjacency(cube, :time)
    @test_throws "one dimension" neighbors(cube, (Cells, :time))
    @test_throws "carries a cell lookup" adjacency(
        DD.DimArray(collect(1.0:4), (DD.X(1:4),)))
end

# ===========================================================================
# Radius-k and exact-ring sweeps: `neighborhood = Disc(k) | Ring(k)`.
#
# One law covers every case below — the callback's ring argument is the direct
# query the selector names, clipped to the subset — so the oracle is that verb
# and never a second implementation of the shell walk. Everything the sweep
# could get wrong (order, membership, length, the centre, duplicates) is a
# consequence of that one equality.
# ===========================================================================

# The verb each selector defers to.
selected(cv, c, ::Disc{K}; connectivity = Vertex()) where {K} =
    neighbors(cv, c, K; connectivity)
selected(cv, c, ::Ring{K}; connectivity = Vertex()) where {K} =
    ring(cv, c, K; connectivity)

# Local indices hold order, membership and length in one comparable object.
naive_rings(cv, nb; connectivity = Vertex()) =
    [[localindex(cv, x) for x in selected(cv, cv[k], nb; connectivity)]
     for k in eachindex(cv)]

sweep_rings(cv, nb; kw...) =
    mapneighbors((c, nbrs) -> [localindex(h) for h in nbrs], cv;
        neighborhood = nb, kw...)

# The two laws that hold at every visit, whatever the selector names.
function visits_are_clean(cv, nb; connectivity = Vertex())
    clean = true
    foreachneighbors(cv; neighborhood = nb, connectivity) do c, nbrs
        ids = [localindex(h) for h in nbrs]
        clean &= !(localindex(c) in ids) && allunique(ids)
    end
    return clean
end

# `Vertex()` and `Edge()` coincide where exactly three cells meet at every
# vertex — the icosahedral hex family the directed-walk hook names — and are
# two different neighbourhoods everywhere else, so only those systems are
# swept twice. The main testset checks the trait really does say this.
conns(sys) = ishexwalk(sys) ? (Vertex(),) : (Vertex(), Edge())

const SELECTORS = (Disc(0), Disc(1), Disc(2), Disc(3), Ring(1), Ring(2), Ring(3))

@testset "$(syslabel(sys)): the ring argument is the selector's own query" for
        (sys, base, depth, covlvl) in SWEEP

    sub = CellVector(rooted_pg(sys, base, depth))
    coverage = CellVector(query(sys, MultiOrderCoverage(TILE); level=covlvl))

    # What `conns` skips, and why it is safe to skip it.
    @test ishexwalk(sys) ==
          all(c -> neighbors(sub, c) == neighbors(sub, c; connectivity = Edge()), sub)

    @testset "$label" for (label, cv) in
                          ("one rooted subtree" => sub,
                           "multi-window coverage" => coverage)
        @testset "$nb under $conn" for conn in conns(sys), nb in SELECTORS
            @test sweep_rings(cv, nb; connectivity = conn, threaded = false) ==
                  naive_rings(cv, nb; connectivity = conn)
            @test visits_are_clean(cv, nb; connectivity = conn)
        end
    end
end

@testset "$(syslabel(sys)): Disc(1) written out is the default" for
        (sys, base, depth, covlvl) in SWEEP

    cv = CellVector(rooted_pg(sys, base, depth))
    for conn in conns(sys)
        @test mapneighbors(probe, cv; connectivity = conn,
            neighborhood = Disc(1)) == mapneighbors(probe, cv; connectivity = conn)
    end
end

# A subtree with its middle third removed. Distance stays the system's, so the
# hole shortens rings without lengthening any path around itself: from `k == 2`
# the clip is the only thing that separates this fixture from the whole
# subtree, and it is the only thing the sweep may do.
holed(cv) = cv[vcat(1:(length(cv) ÷ 3), (2 * (length(cv) ÷ 3) + 1):length(cv))]

@testset "$(syslabel(sys)): a hole clips, and only clips" for
        (sys, base, depth, covlvl) in SWEEP

    cv = CellVector(rooted_pg(sys, base, depth))
    hv = holed(cv)
    @test length(hv) < length(cv)

    @testset "$nb under $conn" for conn in conns(sys),
                                   nb in (Disc(2), Disc(3), Ring(2), Ring(3))
        @test sweep_rings(hv, nb; connectivity = conn, threaded = false) ==
              naive_rings(hv, nb; connectivity = conn)
        @test visits_are_clean(hv, nb; connectivity = conn)
    end

    # The fixture earns its place: the removed middle really does shorten rings
    # the whole subtree answers in full. Indices below the hole are the same
    # cell in both vectors.
    whole = mapneighbors((c, r) -> length(r), cv; neighborhood = Disc(2),
        threaded = false)
    cut = mapneighbors((c, r) -> length(r), hv; neighborhood = Disc(2),
        threaded = false)
    @test any(k -> cut[k] < whole[k], 1:(length(cv) ÷ 3))
end

# Pentagons are a per-system expectation. H3 carries twelve at every level, and
# a complete level clips nothing, so a short shell there is the pentagon's own
# deficiency rather than an edge of the subset.
@testset "H3 pentagons answer shorter shells" begin
    sys = DGG.H3System()
    grid = levelgrid(sys, 1)
    cv = CellVector(grid)
    pent = [k for k in eachindex(cv) if neighborcount(grid, cv[k]) == 5]
    @test length(pent) == 12

    @testset "$nb" for nb in (Disc(1), Disc(2), Disc(3), Ring(1), Ring(2), Ring(3))
        rings = sweep_rings(cv, nb; threaded = false)
        @test rings == naive_rings(cv, nb)
        @test visits_are_clean(cv, nb)
        hexmax = maximum(length(rings[k]) for k in eachindex(cv) if !(k in pent))
        @test all(k -> length(rings[k]) < hexmax, pent)
    end
end

# A complete level holds every face seam. A `CustomOrder` system measures one
# centroid azimuth per shell member and sorts, which is the `k ≥ 2` path most
# likely to disagree with the verb the sweep claims to be.
customorder(sys) = DGG.winding(basesystem(sys), Vertex()) isa DGG.CustomOrder

@testset "shells cross seams: $(syslabel(sys))" for
        sys in filter(customorder, DGG.systems())

    cv = CellVector(levelgrid(sys, 2))
    @testset "$nb under $conn" for conn in conns(sys),
                                   nb in (Disc(2), Disc(3), Ring(2), Ring(3))
        @test sweep_rings(cv, nb; connectivity = conn, threaded = false) ==
              naive_rings(cv, nb; connectivity = conn)
        @test visits_are_clean(cv, nb; connectivity = conn)
    end
end

@testset "the seam systems are still selected by winding" begin
    # An empty filter above would pass every law in it vacuously.
    @test !isempty(filter(customorder, DGG.systems()))
end

# Every calling form reads the same neighbourhood; the forms differ only in
# what they hand the callback about it. One system suffices — the selector is
# threaded through system-generic code and the numbers are pinned to the
# handle form's.
@testset "every calling form carries the selector" begin
    sys = DGG.IGeo7System()
    pg = rooted_pg(sys, 1, 3)
    cv = CellVector(pg)
    n = length(cv)
    data = collect(1.0:n)
    cubedata = [j * data[k] + 0.1j for j in 1:3, k in 1:n]
    A = DD.DimArray(copy(data), (Cells(CellLookup(cv)),))
    cube = DD.DimArray(copy(cubedata),
        (DD.Dim{:time}(1:3), Cells(CellLookup(cv))))
    metric = (c, v, vals) ->
        3.0v + sum(i * vals[i] for i in eachindex(vals); init = 0.0)
    sliced = (c, s, ns) -> sum(s) + sum(sum, ns; init = 0.0)

    @testset "$nb" for nb in (Disc(0), Disc(2), Ring(2), Ring(3))
        want = naive_rings(cv, nb)
        @test sweep_rings(cv, nb; threaded = false) == want

        # The iterator yields the same rings in the same order, from every
        # collection that offers it: the vector, the subset, the lookup and
        # the dimarray over the lookup.
        for src in (cv, pg, CellLookup(cv), A)
            @test [[localindex(h) for h in nbrs]
                   for (_, nbrs) in neighbors(src; neighborhood = nb)] == want
        end

        # Positional data gathers the ring, slot for slot.
        want_data = [3.0data[k] +
                     sum(i * data[want[k][i]] for i in eachindex(want[k]);
                         init = 0.0) for k in eachindex(cv)]
        @test mapneighbors(metric, cv, data; neighborhood = nb,
            threaded = false) == want_data

        # A field request answers over the same cells, field-major.
        got = mapneighbors(cv; needs = (Cell(),), neighborhood = nb,
            threaded = false) do center, rings
            [localindex(cv, x) for x in rings[1]]
        end
        @test got == want

        # foreachneighbors fills a store the caller owns.
        store = [Int[] for _ in 1:n]
        foreachneighbors(cv; neighborhood = nb) do c, nbrs
            store[localindex(c)] = [localindex(h) for h in nbrs]
        end
        @test store == want

        # The DimArray forms: values, and one view per selected neighbour.
        @test parent(mapneighbors(metric, A; pass = Values(),
            neighborhood = nb, threaded = false)) == want_data
        outS = mapneighbors(sliced, cube; pass = NeighborSlices(),
            neighborhood = nb, threaded = false)
        @test parent(outS) ≈ [sum(cubedata[:, k]) +
                              sum(sum(cubedata[:, j]) for j in want[k]; init = 0.0)
                              for k in eachindex(cv)]
    end
end

@testset "traversal does not change a Disc(2) sweep" begin
    sys = DGG.IGeo7System()
    cv = CellVector(rooted_pg(sys, 1, 3))
    want = naive_rings(cv, Disc(2))
    perm = shuffle(Xoshiro(13), 1:length(cv))
    @test sweep_rings(cv, Disc(2); threaded = false) == want
    @test sweep_rings(cv, Disc(2); threaded = true) == want
    @test sweep_rings(cv, Disc(2); order = perm, threaded = false) == want
    @test sweep_rings(cv, Disc(2); order = perm, threaded = true) == want
end

# The message, not just the type: a caller who wrote the wrong thing is told
# which spellings are the right ones.
errmsg(f) = try
    (f(); "")
catch e
    sprint(showerror, e)
end

@testset "a neighbourhood is Disc(k) or Ring(k), and says so" begin
    sys = DGG.IGeo7System()
    cv = CellVector(rooted_pg(sys, 1, 3))
    data = collect(1.0:length(cv))
    A = DD.DimArray(copy(data), (Cells(CellLookup(cv)),))

    # The constructors refuse before any sweep is built.
    @test_throws ArgumentError Ring(0)
    @test_throws ArgumentError Disc(-1)
    @test_throws ArgumentError Ring(-1)
    zero_ring = errmsg(() -> Ring(0))
    @test contains(zero_ring, "Disc(0)") && contains(zero_ring, "Ring(1)")
    @test contains(errmsg(() -> Disc(-1)), "non-negative")

    # A value that is neither selector is refused by every entry point.
    for f in (() -> mapneighbors(probe, cv; neighborhood = 3),
              () -> mapneighbors(probe, cv, data; neighborhood = 3),
              () -> mapneighbors(probe, cv; needs = (Cell(),), neighborhood = 3),
              () -> foreachneighbors(probe, cv; neighborhood = 3),
              () -> neighbors(cv; neighborhood = 3),
              () -> mapneighbors(probe, A; neighborhood = 3))
        @test_throws ArgumentError f()
        msg = errmsg(f)
        @test contains(msg, "Disc(k)") && contains(msg, "Ring(k)")
    end
end

@testset "the callback's ring type follows the selected capacity" begin
    sys = DGG.IGeo7System()
    cv = CellVector(rooted_pg(sys, 1, 3))
    H = DGG.SubsetIndexedCell{eltype(cv)}
    @testset "$nb under $conn" for conn in (Vertex(), Edge()), nb in SELECTORS
        @test eltype(mapneighbors((c, r) -> r, cv; neighborhood = nb,
            connectivity = conn)) ==
              EN._ringtype(EN._capacity(cv.grid, nb, conn), H)
    end
end

end # module MapNeighborsTests
