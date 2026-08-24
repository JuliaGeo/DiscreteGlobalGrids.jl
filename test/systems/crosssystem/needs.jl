# What a requested field answers, slot by slot, against the per-cell verbs.

module NeedsTests

using Test
import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import Extents
using Random: shuffle, Xoshiro

using DiscreteGlobalGrids: levelgrid, cellindex, localindex, globalindex,
    neighbors, mapneighbors, foreachneighbors, subtree,
    CellVector, CellLookup, Cells, MultiOrderCoverage, AuthalicSystem, Vertex,
    Edge, query, cellid, cell_centroid, rawid, reindex, cellindextypes,
    Cell, Index, Local, Global, Value, Centroid,
    Neighbors, Values, NeighborSlices

include(joinpath(@__DIR__, "..", "..", "helpers.jl"))
using .DGGTestHelpers: syslabel, sweepcovers

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

# The record the sweep is asked for: the whole center tuple and the whole
# ring tuple, so every slot of every need is compared.
record(center, rings) = (center, rings)

# The same answers assembled one cell at a time from the per-cell verbs.
naive_center(cv, data, k) =
    (cv[k], k, globalindex(cv, cv[k]), data[k], cell_centroid(cv.grid, cv[k]))

function naive_rings(cv, data, k; connectivity = Vertex())
    nbrs = collect(neighbors(cv, cv[k]; connectivity))
    return ([cellid(x) for x in nbrs],
        [localindex(cv, x) for x in nbrs],
        [globalindex(cv, x) for x in nbrs],
        [data[localindex(cv, x)] for x in nbrs],
        [cell_centroid(cv.grid, x) for x in nbrs])
end

# One number that reads the visited cell and every slot of the `Value` ring,
# weighted by slot, so a dropped, refetched or reordered ring fails.
weighted(center, rings) = Float64(rawid(center[1])) + 3.0 * center[2] +
    sum(i * rings[2][i] for i in eachindex(rings[2]); init = 0.0)

@testset "$(syslabel(sys))" for (sys, base, depth, covlvl) in SWEEP
    sub = CellVector(rooted_pg(sys, base, depth))
    coverage = CellVector(query(sys, MultiOrderCoverage(TILE); level=covlvl))

    @testset "$label" for (label, cv) in
                          ("one rooted subtree" => sub,
                           "multi-window coverage" => coverage)
        n = length(cv)
        data = collect(1.0:n)
        needs = (Cell(), Index(Local()), Index(Global()), Value(data),
            Centroid())

        # The subtree fixture has local ≠ global, so `Index(Global())` served
        # from the local index dies here.
        @testset "every slot is the per-cell verb's answer" begin
            for conn in (Vertex(), Edge())
                centers, rings = mapneighbors(record, cv; needs,
                    threaded = false, connectivity = conn)
                @test length(centers) == n && length(rings) == n
                @test all(k -> centers[k] == naive_center(cv, data, k), 1:n)
                @test all(1:n) do k
                    map(collect, rings[k]) ==
                        naive_rings(cv, data, k; connectivity = conn)
                end
                # Slot `i` of every ring names the same neighbour: the rings
                # are field-major, all of one length.
                @test all(k -> allequal(map(length, rings[k])), 1:n)

                # Traversal mode does not change index-ordered results.
                @test mapneighbors(record, cv; needs, threaded = true,
                    connectivity = conn) == (centers, rings)
                perm = shuffle(Xoshiro(42), 1:n)
                @test mapneighbors(record, cv; needs, order = perm,
                    threaded = false, connectivity = conn) == (centers, rings)
                @test mapneighbors(record, cv; needs, order = perm,
                    threaded = true, connectivity = conn) == (centers, rings)
            end
        end

        # Every declared id scheme is reachable, and answers `reindex`.
        @testset "Index(T) is reindex" begin
            for T in cellindextypes(sys)
                ids, ringids = mapneighbors((c, r) -> (c[1], collect(r[1])), cv;
                    needs = (Index(T),), threaded = false)
                @test eltype(ids) === T
                @test ids == [reindex(T, sys, cv[k]) for k in 1:n]
                @test all(1:n) do k
                    ringids[k] == [reindex(T, sys, x)
                                   for x in neighbors(cv, cv[k])]
                end
            end
        end

        # Two arrays of different eltypes give two rings of those eltypes:
        # neither ring's eltype comes from the first array.
        @testset "two Values keep their own eltypes" begin
            f32 = collect(Float32, 1:n)
            i16 = collect(Int16, 1:n)
            types = Ref{Any}(nothing)
            vals, ringvals = mapneighbors(cv; needs = (Value(f32), Value(i16)),
                threaded = false) do center, rings
                types[] = (typeof(center), typeof(rings))
                (center, map(collect, rings))
            end
            CT, RT = types[]
            @test CT === Tuple{Float32,Int16}
            @test map(eltype, fieldtypes(RT)) === (Float32, Int16)
            @test eltype(ringvals) === Tuple{Vector{Float32},Vector{Int16}}
            @test all(k -> vals[k] === (f32[k], i16[k]), 1:n)
            @test all(1:n) do k
                nbrs = neighbors(cv, cv[k])
                ringvals[k] == ([f32[localindex(cv, x)] for x in nbrs],
                    [i16[localindex(cv, x)] for x in nbrs])
            end
        end

        # A concrete tuple result is one vector per component, index ordered.
        @testset "a tuple result splits per component" begin
            a, b = mapneighbors((c, r) -> (c[2], Float64(length(r[1]))), cv;
                needs = (Cell(), Index(Local())), threaded = false)
            @test a isa Vector{Int} && b isa Vector{Float64}
            @test a == collect(1:n)
            @test b == [Float64(length(neighbors(cv, c))) for c in cv]
        end
    end

    # A field request reaches the other two collections through their existing
    # keyword forwards: a `needs` dropped on the way changes the callback's
    # arity, so the answer cannot merely differ, it cannot be produced.
    @testset "the grid and lookup forms are the vector's" begin
        pg = rooted_pg(sys, base, depth)
        data = collect(1.0:length(sub))
        needs = (Cell(), Value(data))
        want = mapneighbors(weighted, sub; needs, threaded = false)
        @test mapneighbors(weighted, pg; needs, threaded = false) == want
        @test mapneighbors(weighted, CellLookup(sub); needs,
            threaded = false) == want
    end
end

@testset "a malformed request is refused" begin
    sys = DGG.IGeo7System()
    cv = CellVector(rooted_pg(sys, 1, 3))
    n = length(cv)
    data = collect(1.0:n)

    @test_throws ArgumentError mapneighbors(record, cv; needs = ())
    @test_throws ArgumentError foreachneighbors(record, cv; needs = ())
    @test_throws ArgumentError mapneighbors(record, cv;
        needs = (Cell(), :centroid))
    @test_throws ArgumentError mapneighbors(record, cv;
        needs = (Value(collect(1.0:(n - 1))),))
    # One request without its trailing comma is not a tuple; the message says
    # so with the value and the form that would have worked.
    @test_throws "did you mean" mapneighbors(record, cv; needs = Cell())

    # A malformed space is refused where it is written, not by a MethodError
    # from inside the sweep. `Index(Local)` names the type where the space was
    # meant — the typo the documented `Index(T)` form invites.
    @test_throws ArgumentError Index(:local)
    @test_throws ArgumentError Index(Float64)
    @test_throws "`Local` is the type" Index(Local)
    @test_throws "`Global` is the type" Index(Global)
    # An id type the system does not declare is refused up front, not by
    # `reindex` at the first cell: IGeo7 answers in `Z7Cell` alone, and an
    # empty collection — which has no first cell — must refuse it too.
    @test !(DGG.LevelIndex in cellindextypes(sys))
    empty = cv[1:0]
    @test length(empty) == 0
    @test_throws ArgumentError mapneighbors(record, cv;
        needs = (Index(DGG.LevelIndex),))
    @test_throws ArgumentError mapneighbors(record, empty;
        needs = (Index(DGG.LevelIndex),))
    @test_throws ArgumentError foreachneighbors(record, empty;
        needs = (Index(DGG.LevelIndex),))
    # The scheme the system does declare sweeps that same empty collection.
    @test isempty(mapneighbors((c, r) -> 1.0, empty;
        needs = (Index(first(cellindextypes(sys))),)))
    # A `Value` is one vector on the cell axis, so a cube is refused at the
    # request, not at the first read.
    @test_throws ArgumentError Value(zeros(2, n))
    @test_throws ArgumentError Value(1.0)
    # A field request and a positional data vector are two ways to say the
    # same thing; asking for both leaves the callback's arity undefined.
    @test_throws ArgumentError mapneighbors(record, cv, data;
        needs = (Cell(),))
    @test_throws ArgumentError foreachneighbors(record, cv, data;
        needs = (Cell(),))
end

# Keep mutable state in the functor while measuring sweep allocations.
struct NeedsSum <: Function
    acc::Base.RefValue{Float64}
end

function (s::NeedsSum)(center, rings)
    c, k, g, v, p = center
    cs, ls, gs, vs, ps = rings
    t = Float64(rawid(c)) + Float64(k) + Float64(g) + v + p[1]
    for i in eachindex(ls)
        t += Float64(rawid(cs[i])) + Float64(ls[i]) + Float64(gs[i]) +
             vs[i] + ps[i][1]
    end
    s.acc[] += t
    return nothing
end

@testset "the sequential needs sweep allocates nothing" begin
    sys = DGG.IGeo7System()
    coverage = CellVector(query(sys, MultiOrderCoverage(TILE); level=8))
    data = collect(1.0:length(coverage))
    needs = (Cell(), Index(Local()), Index(Global()), Value(data), Centroid())
    s = NeedsSum(Ref(0.0))
    foreachneighbors(s, coverage; needs)
    @test @allocated(foreachneighbors(s, coverage; needs)) == 0
end

# A field request on a cube runs the same sweep on the cell axis and hands
# back the caller's own lookups. One system suffices: the forward is
# system-generic and the numbers are pinned to the CellVector layer's.
@testset "a field request sweeps a dimarray's cell axis" begin
    sys = DGG.IGeo7System()
    cv = CellVector(rooted_pg(sys, 1, 3))
    n = length(cv)
    data = collect(1.0:n)
    A = DD.DimArray(copy(data), (Cells(CellLookup(cv)),))
    # Put the cell dimension second to exercise automatic detection; each
    # time slice carries its own numbers, none of which the sweep may read.
    cubedata = [j * data[k] + 0.1j for j in 1:3, k in 1:n]
    cube = DD.DimArray(copy(cubedata), (DD.Dim{:time}(1:3), Cells(CellLookup(cv))))

    want = mapneighbors(weighted, cv; needs = (Cell(), Value(parent(A))),
        threaded = false)

    # A one-dimensional dim array is already a vector laid out against the
    # collection, so `Value(A)` needs no unwrapping.
    out = mapneighbors(weighted, A; needs = (Cell(), Value(A)), threaded = false)
    @test out isa DD.AbstractDimArray
    @test parent(out) == want
    @test DD.dims(out) === DD.dims(A)
    @test parent(DD.lookup(out, 1)) === cv

    # The sweep never reads the array: the same request on a cube whose cell
    # dimension is second gives one result per cell, on the cell dimension.
    outc = mapneighbors(weighted, cube; needs = (Cell(), Value(data)),
        threaded = false)
    @test outc isa DD.AbstractDimArray
    @test size(outc) == (n,)
    @test map(DD.name, DD.dims(outc)) == (:Cells,)
    @test parent(DD.lookup(outc, 1)) === cv
    @test parent(outc) == want

    # The side-effecting form visits every cell of the cell axis once.
    visits = zeros(Int, n)
    foreachneighbors(cube; needs = (Index(Local()),)) do center, rings
        visits[center[1]] += 1
    end
    @test visits == ones(Int, n)

    # `pass` and `needs` both name what the callback receives, so only the
    # default `pass` may accompany a field request.
    @test_throws ArgumentError mapneighbors(weighted, A; pass = Values(),
        needs = (Cell(), Value(data)))
    @test_throws ArgumentError mapneighbors(weighted, cube;
        pass = NeighborSlices(), needs = (Cell(), Value(data)))
    @test_throws ArgumentError foreachneighbors(weighted, cube; pass = Values(),
        needs = (Cell(), Value(data)))
    # The default pass is untouched by the keyword's arrival.
    @test parent(mapneighbors(weighted, A; pass = Neighbors(),
        needs = (Cell(), Value(A)), threaded = false)) == want
end

end # module NeedsTests
