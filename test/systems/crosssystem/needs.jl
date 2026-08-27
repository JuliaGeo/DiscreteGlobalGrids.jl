# What a requested field answers, slot by slot, against the per-cell verbs.

module NeedsTests

using Test
import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import Extents
using Random: shuffle, Xoshiro

using DiscreteGlobalGrids: levelgrid, cellindex, localindex, globalindex,
    neighbors, mapneighbors, mapneighbors!, foreachneighbors, subtree,
    CellVector, CellLookup, Cells, MultiOrderCoverage, AuthalicSystem, Vertex,
    Edge, query, cellid, cell_centroid, rawid, reindex, cellindextypes,
    Cell, Index, Local, Global, Value, Centroid, cellfield, border, ncells,
    chunkplan, nchunks, Neighbors, Values, NeighborSlices

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

# Keep mutable state in the functor while measuring sweep allocations: a
# top-level `@allocated` through a non-const global boxes and reports bytes
# that are not there.
struct FieldSum <: Function
    acc::Base.RefValue{Float64}
end

function (s::FieldSum)(center, rings)
    c, k, g, v = center
    cs, ls, gs, vs = rings
    t = Float64(rawid(c)) + Float64(k) + Float64(g) + v
    for i in eachindex(ls)
        t += Float64(rawid(cs[i])) + Float64(ls[i]) + Float64(gs[i]) + vs[i]
    end
    s.acc[] += t
    return nothing
end

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

# Cells, indices and values are read out of state the sweep already has, so a
# request for them alone still allocates nothing at all.
@testset "the sequential needs sweep allocates nothing" begin
    sys = DGG.IGeo7System()
    coverage = CellVector(query(sys, MultiOrderCoverage(TILE); level=8))
    data = collect(1.0:length(coverage))
    needs = (Cell(), Index(Local()), Index(Global()), Value(data))
    s = FieldSum(Ref(0.0))
    foreachneighbors(s, coverage; needs)
    @test @allocated(foreachneighbors(s, coverage; needs)) == 0
end

# One sweep's centroid bytes, measured inside a function that receives the
# collection.
function centroidbytes(cv)
    data = collect(1.0:length(cv))
    needs = (Cell(), Index(Local()), Index(Global()), Value(data), Centroid())
    s = NeedsSum(Ref(0.0))
    foreachneighbors(s, cv; needs)
    return @allocated foreachneighbors(s, cv; needs)
end

# One `mapneighbors` sweep's bytes, again measured inside a function that
# receives what it sweeps. Both kernels below answer a `Float64`, so the two
# requests allocate the same output vector and their difference is the
# centroid machinery alone.
function sweepbytes(f, cv, needs; threaded)
    mapneighbors(f, cv; needs, threaded)
    return @allocated mapneighbors(f, cv; needs, threaded)
end

valuesum(center, rings) = center[1] + sum(rings[1]; init = 0.0)

# Reads the centre value, the centre centroid and every slot of both rings, so
# nothing the sweep hands it can be elided.
function slopesum(center, rings)
    v0, p0 = center
    vs, ps = rings
    t = 0.0
    for i in eachindex(vs)
        t += (v0 - vs[i]) * (p0[1] - ps[i][1])
    end
    return t
end

# `Centroid()` buys a working set, so the law is per cell rather than per
# sweep: one window per task, whatever the sweep's length. The two sequential
# fixtures are 47,659 and 117,649 cells — 2.47x apart, and both past the
# 16,384-slot cap — so both pay for one window of exactly the same size, and a
# sweep that grew its cache with the collection, or rebuilt it per cell,
# differs here. The threaded sweep then pays one such window per task and no
# more: a window shared between tasks would be a race, and one per cell would
# be the allocation this whole design exists to avoid.
@testset "the centroid working set is per task, not per cell" begin
    sys = DGG.IGeo7System()
    coverage = CellVector(query(sys, MultiOrderCoverage(TILE); level=10))
    sub = CellVector(rooted_pg(sys, 1, 6))
    @test length(coverage) == 47659
    @test length(sub) == 7^6
    perslot = sizeof(Int) + sizeof(typeof(cell_centroid(sub.grid, sub[1])))
    window = 16384 * perslot
    small = centroidbytes(coverage)
    large = centroidbytes(sub)
    @test small == large
    @test window <= small < 1 << 20

    n = length(sub)
    data = collect(1.0:n)
    ranges = DGG.Engine._chunk_ranges(n)
    @test length(ranges) > 1
    # What one window costs: the sequential sweep's centroid bytes over the
    # same sweep without one. That is `window` plus the two array headers,
    # which is what the threaded bound must be stated in.
    seqbase = sweepbytes(valuesum, sub, (Value(data),); threaded = false)
    onewindow = sweepbytes(slopesum, sub, (Value(data), Centroid());
        threaded = false) - seqbase
    @test window <= onewindow < window + 1024
    # The output vector and the sweep's own overhead, priced by the same call
    # with no centroid in it.
    output = sweepbytes(valuesum, sub, (Value(data),); threaded = true)
    bytes = sweepbytes(slopesum, sub, (Value(data), Centroid()); threaded = true)
    # A task's window is sized from its own range, so price it from that range
    # rather than from the sequential sweep's, and allow the 512 bytes a task's
    # callback box costs to be spawned: the field `Centroid()` resolves to is
    # named inside that box. Both are per task, not per cell.
    taskwindow = min(nextpow(2, length(first(ranges))), 16384) * perslot +
                 (onewindow - window)
    @test bytes <= length(ranges) * (taskwindow + 512) + output
    # And at least one window per task: a task's window holds its whole range
    # up to the cap, so nothing smaller than this can be one per task.
    @test bytes >= length(ranges) * min(length(first(ranges)), 16384) * perslot +
                   output
end

# One read of a field, measured inside a function that receives it, and
# summing a component so no read is dead code.
readfield!(acc, a, ks) = (for k in ks
    acc[] += a[k][1]
end; nothing)

function fieldbytes(a, ks)
    acc = Ref(0.0)
    readfield!(acc, a, ks)
    return @allocated readfield!(acc, a, ks)
end

# `Centroid()` is the collection's centroid field asked for by name, so the
# spelled-out field must answer the same numbers however much of itself it
# already knows. Each `known` below kills a different reader: an ignored
# `known` (the table row), a subset probed with the swept collection's index
# instead of the subset's own (the border row — those two index spaces differ
# on every cell but the first), and an entry the subset does not carry read
# rather than computed (the empty row, where every read must fall through).
@testset "a centroid field answers Centroid(), whatever it knows" begin
    sys = DGG.IGeo7System()
    cv = CellVector(rooted_pg(sys, 1, 3))
    n = length(cv)
    data = collect(1.0:n)
    table = [cell_centroid(cv.grid, c) for c in cv]

    bd = cv[collect(DGG.border(cv))]
    @test 0 < length(bd) < n
    bcube = DD.DimArray([cell_centroid(cv.grid, c) for c in bd],
        (Cells(CellLookup(bd)),))
    ecube = DD.DimArray(eltype(table)[], (Cells(CellLookup(cv[1:0])),))

    want = mapneighbors(record, cv; needs = (Value(data), Centroid()),
        threaded = false)
    @testset "known = $label" for (label, known) in
                                  ("nothing" => nothing,
                                   "the whole table" => table,
                                   "the border" => bcube,
                                   "an empty subset" => ecube)
        a = DGG.cellfield(cell_centroid, cv; known)
        @test eltype(a) === eltype(table)
        # The field is a vector in its own right: readable with no sweep
        # around it, and answering the per-cell verb at every index.
        @test collect(a) == table
        @test mapneighbors(record, cv; needs = (Value(data), Value(a)),
            threaded = false) == want
        # The field is shared by every task and the readers are not, so a
        # threaded sweep answers the sequential one's numbers.
        @test mapneighbors(record, cv; needs = (Value(data), Value(a)),
            threaded = true) == want
        # Reading is pure: nothing is remembered, so the second read of the
        # same index answers the same value and costs no allocation.
        @test fieldbytes(a, 1:n) == 0
    end

    # A field is bounds-checked like the vector it is.
    @test_throws BoundsError DGG.cellfield(cell_centroid, cv)[n + 1]

    # A field is read by local index, so one built over OTHER cells of the same
    # length answers a different cell's value in every slot without a single
    # index leaving the axis: the layout check cannot see it and the
    # collection check must. Sibling subtrees at the same level are exactly
    # that pair.
    other = CellVector(subtree(sys, cellindex(levelgrid(sys, 1), 4), 4))
    @test length(other) == n
    @test other != cv
    @test_throws ArgumentError mapneighbors(record, cv;
        needs = (Value(data), Value(DGG.cellfield(cell_centroid, other))),
        threaded = false)

    # The check is over cells, not over object identity: the dim-array route
    # sweeps the lookup's own vector, and a field the caller built by naming
    # the cube's cells is over the same cells in a different object. Rejecting
    # it would make the cube route unusable with a field.
    A = DD.DimArray(copy(data), (Cells(CellLookup(cv)),))
    rebuilt = CellVector(DGG.cellset(DD.lookup(A, 1)))
    @test rebuilt !== parent(DD.lookup(A, 1))
    @test rebuilt == cv
    @test map(parent, mapneighbors(record, A;
        needs = (Value(data), Value(DGG.cellfield(cell_centroid, rebuilt))),
        threaded = false)) == want

    # A complete field is read straight through: the sweep builds no window
    # for it, so a sweep over one allocates exactly what the same sweep with
    # no field in it does. A window built anyway would show up here as half a
    # megabyte.
    big = CellVector(rooted_pg(sys, 1, 5))
    bigdata = collect(1.0:length(big))
    full = DGG.cellfield(cell_centroid, big;
        known = [cell_centroid(big.grid, c) for c in big])
    @test sweepbytes(slopesum, big, (Value(bigdata), Value(full));
        threaded = false) ==
          sweepbytes(valuesum, big, (Value(bigdata),); threaded = false)
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

# A request stated about one collection, answered chunk by chunk. A chunk is a
# partial grid whose `1:ncells` is its own, so every number the sweep reports
# has to be translated back — and the law is that the chunked answer IS the
# whole-axis answer, cell for cell. The fixtures are subsets, where the local
# index is not the global one, so a translation that confused the two dies
# here; the read side of the same route is asserted in `test/io/chunk_sweep.jl`.
@testset "a chunked field request answers in the caller's own indices" begin
    sys = DGG.IGeo7System()
    lvl = 4
    grid = levelgrid(sys, lvl)
    # One rooted subtree, and a scattered subset of the same level: the first
    # is one window, the second is thousands, so the halo of a chunk of it is
    # a very different shape.
    scattered = CellVector(sys, lvl,
        [cellindex(grid, p) for p in
         sort(shuffle(Xoshiro(7), 1:ncells(grid))[1:3000])])
    fixtures = ("one rooted subtree" => CellVector(rooted_pg(sys, 1, 3)),
        "a scattered subset" => scattered)

    # Every translated quantity, in one record: the visited cell's local index
    # and the sum of its ring's, so a chunk-local number fails on every chunk
    # but the first; the same in global indices, which no chunk may change; a
    # value, a re-encoded id, a centroid and the ring's length.
    function record(center, rings)
        _, k, g, t, v, p = center
        cs, ks, gs, ts, vs, ps = rings
        return (k, sum(ks; init = 0), g, sum(gs; init = 0),
            v + sum(vs; init = 0.0),
            rawid(t) + sum(rawid(x) for x in ts; init = zero(rawid(t))),
            p[1] + sum(q[1] for q in ps; init = 0.0), length(cs))
    end

    # Three needs: a value, a centroid field, and the caller's own index.
    function record3(center, rings)
        v, p, k = center
        vs, ps, ks = rings
        return (v + sum(vs; init = 0.0),
            p[1] + sum(q[1] for q in ps; init = 0.0), k, sum(ks; init = 0))
    end

    @testset "$label" for (label, cv) in fixtures
        n = length(cv)
        # The fixtures are subsets: if they were not, an untranslated chunk
        # index could pass by accident.
        @test globalindex(cv, cv[1]) != 1
        data = collect(1.0:n)
        A = DD.DimArray(copy(data), (Cells(CellLookup(cv)),))
        T = last(cellindextypes(sys))
        needs = (Cell(), Index(Local()), Index(Global()), Index(T),
            Value(data), Centroid())
        want = mapneighbors(record, cv; needs, threaded = false)

        plan = chunkplan(A; halo = 1, chunks = cld(n, 4))
        @test nchunks(plan) >= 3

        dest = map(similar, want)
        mapneighbors!(dest, record, A, plan; needs, threaded = false)
        @test all(map(==, dest, want))

        # Threading inside a chunk, and cutting the plan into pieces that each
        # write their own part, change nothing.
        threadeddest = map(similar, want)
        mapneighbors!(threadeddest, record, A, plan; needs, threaded = true)
        @test all(map(==, threadeddest, want))

        splitdest = map(similar, want)
        for p in Base.split(plan, 3)
            mapneighbors!(splitdest, record, A, p; needs, threaded = false)
        end
        @test all(map(==, splitdest, want))

        # The plan-free form plans for itself.
        autodest = map(similar, want)
        mapneighbors!(autodest, record, A; needs, chunks = cld(n, 4),
            threaded = false)
        @test all(map(==, autodest, want))

        # A scalar result lands in a plain destination.
        onesum(center, rings) = center[5] + sum(rings[5]; init = 0.0)
        scalarwant = mapneighbors(onesum, cv; needs, threaded = false)
        scalardest = similar(scalarwant)
        mapneighbors!(scalardest, onesum, A, plan; needs, threaded = false)
        @test scalardest == scalarwant

        # A request never reads the cube it is called on, so a cube with other
        # dimensions is swept the same way and still gives one result per cell,
        # on the cell axis alone — `dest` is a vector however wide `A` is.
        cube = DD.DimArray(repeat(data, 1, 3),
            (Cells(CellLookup(cv)), DD.Dim{:time}(1:3)))
        cubedest = similar(scalarwant)
        mapneighbors!(cubedest, onesum, cube,
            chunkplan(cube; halo = 1, chunks = cld(n, 4)); needs,
            threaded = false)
        @test cubedest == scalarwant

        # A cell field carrying a subset of what the caller already knows says
        # which cells it holds by CELL ID, so no chunk may renumber it: the
        # border's centroids read back the same on the chunk route as off it.
        edge = cv[collect(border(cv))]
        part = cellfield(cell_centroid, cv;
            known = DD.DimArray([cell_centroid(cv.grid, c) for c in edge],
                (Cells(CellLookup(edge)),)))
        table = cellfield(cell_centroid, cv;
            known = [cell_centroid(cv.grid, c) for c in cv])
        fields = ("no" => cellfield(cell_centroid, cv), "a subset" => part,
            "a complete" => table)
        @testset "$kind known" for (kind, field) in fields
            fneeds = (Value(data), Value(field), Index(Local()))
            fwant = mapneighbors(record3, cv; needs = fneeds, threaded = false)
            fdest = map(similar, fwant)
            mapneighbors!(fdest, record3, A, plan; needs = fneeds,
                threaded = false)
            @test all(map(==, fdest, fwant))
        end
    end

    # The request is checked against the collection the caller wrote it about,
    # before a chunk is read: the message names the whole axis's length, not
    # some chunk's.
    cv = CellVector(rooted_pg(sys, 1, 3))
    n = length(cv)
    A = DD.DimArray(collect(1.0:n), (Cells(CellLookup(cv)),))
    plan = chunkplan(A; halo = 1, chunks = cld(n, 4))
    dest = zeros(Float64, n)
    err = try
        mapneighbors!(dest, record3, A, plan;
            needs = (Value(collect(1.0:(n-1))), Index(Local())))
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("1:$n", err.msg)
    # A field over other cells is refused by the same up-front check.
    other = CellVector(subtree(sys, cellindex(levelgrid(sys, 1), 4), 4))
    @test_throws ArgumentError mapneighbors!(dest, record3, A, plan;
        needs = (Value(parent(A)), Value(cellfield(cell_centroid, other))))
end

end # module NeedsTests
