# Following a stored cube's chunk lines.
#
# Two things have to hold at once and they pull in opposite directions: the
# result must equal the whole-axis sweep's, cell for cell, and the reads must
# not. The read count is asserted against a floor computed from the plan itself
# — one read for each chunk's own cells, plus one per DISTINCT storage chunk its
# halo reaches — so any traversal that decodes a chunk twice fails, and so does
# any that gets the answer right by reading everything.

module DGGSIOChunkSweepTests

using Test
import DiscreteGlobalGrids as DGG

const HAS_ZARR = try
    @eval using Zarr
    true
catch err
    @warn "Zarr.jl is not loadable: the chunk-sweep suite is skipped." exception = err
    false
end

if HAS_ZARR

import DimensionalData as DD
using DiscreteGlobalGrids: IGeo7System, levelgrid, cellindex, ncells, Cells,
    CellVector, CellLookup, dggread, dggwrite, mapneighbors, mapneighbors!,
    chunkplan, foreachchunk, chunkcube, localindices, ownedindices,
    axisindices, globalindices, chunkhalo, halowidth, nchunks, region, Values,
    foreachneighbors, localindex, Index, Local, Value, Disc, Ring, neighbors,
    ring

# Counts the chunks read from a store. Metadata keys are not chunks; every other
# key is a chunk of the array whose name prefixes it.
const METADATA_KEYS = (".zarray", ".zattrs", ".zgroup", ".zmetadata")

struct CountingStore{S<:Zarr.AbstractStore} <: Zarr.AbstractStore
    parent::S
    reads::Dict{String,Int}
end
CountingStore(p::Zarr.AbstractStore) = CountingStore(p, Dict{String,Int}())

function Base.getindex(s::CountingStore, k::String)
    if !any(m -> endswith(k, m), METADATA_KEYS)
        name = first(split(k, '/'))
        s.reads[name] = get(s.reads, name, 0) + 1
    end
    return s.parent[k]
end
Base.setindex!(s::CountingStore, v, k::String) = (s.parent[k] = v)
Base.delete!(s::CountingStore, k::String) = delete!(s.parent, k)
Zarr.subdirs(s::CountingStore, p) = Zarr.subdirs(s.parent, p)
Zarr.subkeys(s::CountingStore, p) = Zarr.subkeys(s.parent, p)
Zarr.isinitialized(s::CountingStore, k::AbstractString) =
    Zarr.isinitialized(s.parent, k)
Zarr.storagesize(s::CountingStore, p) = Zarr.storagesize(s.parent, p)

counting(path) = Zarr.zopen(CountingStore(Zarr.DirectoryStore(path)), "r")

# A whole level, so every chunk but the first and last has neighbours on both
# sides and no chunk's halo is empty.
const SYS = IGeo7System()
const LEVEL = 4
const GRID = levelgrid(SYS, LEVEL)
const N = ncells(GRID)
const CHUNK = 512

# The stencil, and what the whole-axis sweep makes of it.
stencil(c, v, vs) = v + sum(vs; init=0.0)

function fixture(dir)
    cells = [cellindex(GRID, p) for p in 1:N]
    A = DD.DimArray(Float64.(1:N),
        Cells(CellLookup(CellVector(SYS, LEVEL, cells))); name=:e)
    path = joinpath(dir, "s.zarr")
    dggwrite(path, DD.DimStack((; e=A)); chunks=CHUNK)
    return A, path
end

# One read per chunk's own cells, plus one per distinct storage chunk its halo
# lands in. This is what "follow the chunk lines" means, counted.
readfloor(plan) =
    sum(1 + length(Set(cld(p, CHUNK) for p in chunkhalo(mc))) for mc in plan)

# The selector's own definition, read back one cell at a time: `Disc(k)` is
# `neighbors(cv, cell, k)` and `Ring(k)` is `ring(cv, cell, k)`. On the complete
# level the local index a direct query answers with IS the axis index, so these
# rings index `parent(A)` directly.
oraclerings(cv, ::Disc{K}) where {K} = [neighbors(cv, p, K) for p in 1:length(cv)]
oraclerings(cv, ::Ring{K}) where {K} = [ring(cv, p, K) for p in 1:length(cv)]

# `stencil` evaluated on those rings.
oracle(v, rs) = [v[p] + sum(q -> v[q], rs[p]; init=0.0) for p in eachindex(v)]

# The halo a selector needs, spelled where the test can read it.
_k(::Disc{K}) where {K} = K
_k(::Ring{K}) where {K} = K

@testset "chunk sweep" begin
    mktempdir() do dir
        A, path = fixture(dir)
        want = parent(mapneighbors(stencil, A; pass=Values(), threaded=false))

        @testset "the plan is the store's own chunk grid, before any read" begin
            B = dggread(counting(path))[:e]
            reads = parent(B).storage.reads
            empty!(reads)
            plan = chunkplan(B; halo=1)
            # Planning walks boundaries, not data.
            @test get(reads, "e", 0) == 0

            @test nchunks(plan) == cld(N, CHUNK)
            @test halowidth(plan) == 1
            # The chunks partition the axis, in order.
            @test reduce(vcat, collect(ownedindices(mc)) for mc in plan) == collect(1:N)
            for mc in plan
                h = chunkhalo(mc)
                @test !isempty(h)                      # a global axis has no isolated chunk
                @test issorted(h) && allunique(h)
                @test all(p -> 1 <= p <= N, h)
                @test all(p -> !(p in ownedindices(mc)), h)
            end
        end

        @testset "a chunk arrives as an ordinary cube" begin
            B = dggread(path)[:e]
            plan = chunkplan(B; halo=1)
            owned = 0
            foreachchunk(B, plan) do cc
                cube = chunkcube(cc)
                # Ordinary: the package's own lookup, not the store's.
                @test DD.lookup(cube, Cells) isa CellLookup
                @test length(cube) ==
                      length(ownedindices(cc)) + length(chunkhalo(cc))
                # Owned and halo partition the block, and owned is contiguous
                # because every halo cell is outside the chunk's own run.
                @test sort(vcat(collect(localindices(cc)), chunkhalo(cc))) ==
                      collect(1:length(cube))
                # The owned rows are the chunk's cells, in the axis's order.
                @test parent(cube)[localindices(cc)] == parent(A)[ownedindices(cc)]
                # Every cube cell knows where it came from in the axis, and the
                # owned pair is one set of cells in the two numberings. A
                # translation that dropped or reordered the halo dies here.
                @test length(axisindices(cc)) == length(cube)
                @test issorted(axisindices(cc)) && allunique(axisindices(cc))
                @test axisindices(cc)[localindices(cc)] ==
                      collect(ownedindices(cc))
                @test sort(vcat(axisindices(cc)[chunkhalo(cc)],
                    axisindices(cc)[localindices(cc)])) == axisindices(cc)
                @test parent(cube) == parent(A)[axisindices(cc)]
                owned += length(localindices(cc))
            end
            @test owned == N
        end

        @testset "globalindices still answers, deprecated" begin
            B = dggread(path)[:e]
            plan = chunkplan(B; halo=1)
            @test globalindices(plan[1]) == ownedindices(plan[1])
            foreachchunk(B, plan) do cc
                @test globalindices(cc) == ownedindices(cc)
            end
        end

        @testset "the chunked result is the whole-axis result" begin
            B = dggread(path)[:e]

            dest = zeros(Float64, N)
            mapneighbors!(dest, stencil, B; threaded=false)
            @test dest == want

            # `mapneighbors` takes the same route by itself, and keeps the axis.
            out = mapneighbors(stencil, B; pass=Values(), threaded=false)
            @test parent(out) == want
            @test DD.lookup(out, Cells) === DD.lookup(B, Cells)
            # Threading inside a chunk changes nothing.
            @test parent(mapneighbors(stencil, B; pass=Values())) == want

            # A concrete tuple result is one array per component here too.
            pair = (c, v, vs) -> (v, Float64(length(vs)))
            a, b = mapneighbors(pair, B; pass=Values(), threaded=false)
            @test parent(a) == parent(A)
            @test parent(b) ==
                  Float64.(parent(mapneighbors((c, nbrs) -> length(nbrs), A)))

            # An N-D cube is swept once per index of the other dimensions.
            C = cat(B, B; dims=DD.Ti(1:2))
            outN = mapneighbors(stencil, C; pass=Values(), threaded=false)
            @test parent(outN)[:, 1] == want
            @test parent(outN)[:, 2] == want
        end

        @testset "results stream back into a store" begin
            B = dggread(path)[:e]
            g = Zarr.zgroup(joinpath(dir, "out.zarr"))
            plan = chunkplan(B; halo=1)

            # Nothing but one chunk of the result is in memory at a time, which
            # is the difference between this and `mapneighbors`.
            Z = Zarr.zcreate(Float64, g, "r", N; chunks=(CHUNK,))
            mapneighbors!(Z, stencil, B, plan; threaded=false)
            @test Z[:] == want

            # Pieces of one plan write disjoint ranges of the same store, so
            # they can run at once without coordinating.
            Z2 = Zarr.zcreate(Float64, g, "r2", N; chunks=(CHUNK,))
            @sync for p in Base.split(plan, 4)
                Threads.@spawn mapneighbors!(Z2, stencil, B, p; threaded=false)
            end
            @test Z2[:] == want
        end

        @testset "splitting the plan is how it parallelises" begin
            B = dggread(path)[:e]
            plan = chunkplan(B; halo=1)
            pieces = Base.split(plan, 4)
            @test sum(nchunks, pieces) == nchunks(plan)
            @test reduce(vcat, [collect(ownedindices(mc)) for p in pieces for mc in p]) ==
                  collect(1:N)
            dest = zeros(Float64, N)
            @sync for p in pieces
                Threads.@spawn mapneighbors!(dest, stencil, B, p; threaded=false)
            end
            @test dest == want
        end

        @testset "each chunk is read once, its halo once per foreign chunk" begin
            B = dggread(counting(path))[:e]
            plan = chunkplan(B; halo=1)
            reads = parent(B).storage.reads
            empty!(reads)
            dest = zeros(Float64, N)
            mapneighbors!(dest, stencil, B, plan; threaded=false)
            @test dest == want
            @test reads["e"] == readfloor(plan)

            # The same sweep cell at a time, which is what the plan replaces.
            C = dggread(counting(path))[:e]
            data, cv = parent(C), region(DD.lookup(C, Cells))
            empty!(data.storage.reads)
            foreachneighbors(cv; threaded=false) do c, nbrs
                data[localindex(c)]
                for h in nbrs
                    data[localindex(h)]
                end
            end
            # Two orders of magnitude is the claim; the exact ratio is the
            # chunk length's and not a law.
            @test reads["e"] * 100 < data.storage.reads["e"]
        end

        # A field request answered chunk by chunk. What is asserted here is the
        # READ side and the route's own entry points; that a translated request
        # equals the whole-axis one on a collection where the local index is
        # not the global one is asserted next to the request itself, in
        # `test/systems/crosssystem/needs.jl`.
        @testset "a field request follows the chunk lines" begin
            # `Index(Local())` in the centre and summed over the ring: a chunk
            # reporting its own numbering fails on every chunk but the first.
            kern(center, rings) = (center[1] + sum(rings[1]; init=0.0),
                Float64(center[2]) + sum(rings[2]; init=0))
            axiscv = region(DD.lookup(A, Cells))
            wa, wb = mapneighbors(kern, axiscv;
                needs=(Value(parent(A)), Index(Local())), threaded=false)

            B = dggread(counting(path))[:e]
            reads = parent(B).storage.reads
            plan = chunkplan(B; halo=1)
            needs = (Value(B), Index(Local()))

            # The automatic route: a chunked parent and storage order, the same
            # rule `Values()` follows.
            empty!(reads)
            ga, gb = mapneighbors(kern, B; needs=needs, threaded=false)
            @test parent(ga) == wa
            @test parent(gb) == wb
            @test DD.lookup(ga, Cells) === DD.lookup(B, Cells)
            # Each storage chunk is decoded a bounded number of times: once for
            # the cube's own cells, once for the request's `Value` over the
            # same store. Per-scalar reading is two orders of magnitude above.
            @test readfloor(plan) <= reads["e"] <= 2 * readfloor(plan)

            # A permutation names a visit order over the whole axis, which a
            # chunked sweep cannot honour, so it keeps the whole-axis path and
            # still answers the same.
            pa, pb = mapneighbors(kern, B; needs=needs, order=collect(N:-1:1),
                threaded=false)
            @test parent(pa) == wa && parent(pb) == wb

            # A `Value` over a store, swept on a cube already in memory: one
            # read per storage chunk it touches, never one per scalar. That is
            # the whole claim of reading a request along the chunk lines.
            C = dggread(counting(path))[:e]
            creads = parent(C).storage.reads
            memplan = chunkplan(A; halo=1, chunks=CHUNK)
            empty!(creads)
            dest = (zeros(Float64, N), zeros(Float64, N))
            mapneighbors!(dest, kern, A, memplan;
                needs=(Value(C), Index(Local())), threaded=false)
            @test dest[1] == wa && dest[2] == wb
            @test creads["e"] == readfloor(memplan)

            # Threading inside a chunk, and cutting the plan into pieces that
            # each write their own part, change nothing.
            d2 = (zeros(Float64, N), zeros(Float64, N))
            mapneighbors!(d2, kern, B, plan; needs=needs)
            @test d2[1] == wa && d2[2] == wb
            d3 = (zeros(Float64, N), zeros(Float64, N))
            @sync for p in Base.split(plan, 4)
                Threads.@spawn mapneighbors!(d3, kern, B, p; needs=needs,
                    threaded=false)
            end
            @test d3[1] == wa && d3[2] == wb

            # The side-effecting form takes the route too, and must not fire
            # for a chunk's halo — those cells are read to complete the owned
            # cells' rings and have incomplete rings of their own, and there is
            # no result to throw away afterwards.
            hits = zeros(Int, N)
            acc = zeros(Float64, N)
            foreachneighbors(B; needs=needs, threaded=false) do center, rings
                hits[center[2]] += 1
                acc[center[2]] = center[1] + sum(rings[1]; init=0.0)
            end
            @test all(==(1), hits)
            @test acc == wa
        end

        # ===================================================================
        # Wider neighbourhoods
        #
        # `neighborhood = Disc(k)` or `Ring(k)` widens what each visit sees,
        # and the plan's halo has to widen with it or the cells at the chunk's
        # edge answer from a truncated ring. Everything below is measured
        # against the same two references: the whole-axis sweep, and a per-cell
        # oracle built from `neighbors`/`ring` themselves — the verbs the
        # selector is DEFINED by, so the oracle is the definition read back.
        # ===================================================================

        axiscv = region(DD.lookup(A, Cells))
        v = parent(A)
        nbrings = Dict(nb => oraclerings(axiscv, nb)
                       for nb in (Ring(2), Disc(2), Disc(3)))
        oracles = Dict(nb => oracle(v, r) for (nb, r) in nbrings)
        wholeaxis = Dict(nb => parent(mapneighbors(stencil, A; pass=Values(),
            neighborhood=nb, threaded=false)) for nb in keys(nbrings))

        @testset "the whole-axis sweep is the two direct queries" begin
            # If this fails nothing below means anything: the oracle and the
            # in-memory sweep would be measuring different things.
            for nb in (Ring(2), Disc(2), Disc(3))
                @test wholeaxis[nb] == oracles[nb]
            end
            # A disc is its rings, so the two selectors are not the same sweep.
            @test wholeaxis[Disc(2)] != wholeaxis[Ring(2)]
            @test wholeaxis[Disc(2)] == want .+ wholeaxis[Ring(2)] .- v
        end

        @testset "a wider neighbourhood is still the whole-axis result" begin
            B = dggread(path)[:e]
            for nb in (Ring(2), Disc(3))
                dest = zeros(Float64, N)
                mapneighbors!(dest, stencil, B; neighborhood=nb, threaded=false)
                @test dest == wholeaxis[nb]
                @test dest == oracles[nb]

                # The plan form, threading inside a chunk, and pieces of one
                # plan each writing their own part: same cells, same answer.
                plan = chunkplan(B; halo=_k(nb))
                @test halowidth(plan) == _k(nb)
                d2 = zeros(Float64, N)
                mapneighbors!(d2, stencil, B, plan; neighborhood=nb)
                @test d2 == wholeaxis[nb]
                d3 = zeros(Float64, N)
                @sync for p in Base.split(plan, 4)
                    Threads.@spawn mapneighbors!(d3, stencil, B, p;
                        neighborhood=nb, threaded=false)
                end
                @test d3 == wholeaxis[nb]
            end
        end

        @testset "a plan narrower than the selector is refused" begin
            B = dggread(path)[:e]
            dest = zeros(Float64, N)
            narrow = chunkplan(B; halo=1)

            # The message names both widths and the plan that would work, so a
            # caller can act on it without reading the source.
            @test_throws ArgumentError mapneighbors!(dest, stencil, B, narrow;
                neighborhood=Ring(2), threaded=false)
            @test_throws "chunkplan(A; halo = 2)" mapneighbors!(dest, stencil, B,
                narrow; neighborhood=Ring(2), threaded=false)
            @test_throws "1 ring" mapneighbors!(dest, stencil, B, narrow;
                neighborhood=Ring(2), threaded=false)
            @test_throws "Ring(2)" mapneighbors!(dest, stencil, B, narrow;
                neighborhood=Ring(2), threaded=false)
            @test_throws "chunkplan(A; halo = 3)" mapneighbors!(dest, stencil, B,
                chunkplan(B; halo=2); neighborhood=Disc(3), threaded=false)
            # Nothing was swept before the refusal.
            @test all(iszero, dest)

            # The needs form checks the same invariant.
            @test_throws ArgumentError mapneighbors!(
                (zeros(Float64, N), zeros(Float64, N)),
                (c, r) -> (0.0, 0.0), B, narrow;
                needs=(Value(B), Index(Local())), neighborhood=Ring(2),
                threaded=false)

            # Wider is fine: the halo carries more context than the sweep asks
            # for and the extra rings are simply not visited.
            wide = chunkplan(B; halo=3)
            mapneighbors!(dest, stencil, B, wide; neighborhood=Ring(2),
                threaded=false)
            @test dest == wholeaxis[Ring(2)]
            # And `Disc(1)`, the default, runs under any plan at all.
            d1 = zeros(Float64, N)
            mapneighbors!(d1, stencil, B, wide; threaded=false)
            @test d1 == want
        end

        @testset "the automatic route carries the selector" begin
            B = dggread(counting(path))[:e]
            reads = parent(B).storage.reads
            # What the route must plan for each selector, if it plans at all.
            twoplan = readfloor(chunkplan(B; halo=2))
            oneplan = readfloor(chunkplan(B; halo=1))

            empty!(reads)
            out = mapneighbors(stencil, B; pass=Values(), neighborhood=Disc(2),
                threaded=false)
            # The route ran a plan, and the plan carried two rings: one ring
            # would read less and the whole-axis fallback would read far more.
            @test reads["e"] == twoplan
            @test twoplan > oneplan
            @test parent(out) == wholeaxis[Disc(2)]
            @test parent(out) == oracles[Disc(2)]
            @test DD.lookup(out, Cells) === DD.lookup(B, Cells)
            # Threading inside a chunk changes nothing here either.
            @test parent(mapneighbors(stencil, B; pass=Values(),
                neighborhood=Disc(2))) == wholeaxis[Disc(2)]

            # `Disc(0)` visits every cell with an empty ring, so the callback
            # sees the cell's own value and nothing else. The plan still
            # carries one ring — it may not carry none — and the answer must
            # not show it.
            empty!(reads)
            z0 = mapneighbors(stencil, B; pass=Values(), neighborhood=Disc(0),
                threaded=false)
            @test parent(z0) == v
            @test reads["e"] == oneplan
            counts = mapneighbors((c, val, vs) -> length(vs), B; pass=Values(),
                neighborhood=Disc(0), threaded=false)
            @test all(iszero, parent(counts))

            # An N-D cube is swept once per index of the other dimensions.
            C = cat(B, B; dims=DD.Ti(1:2))
            outN = mapneighbors(stencil, C; pass=Values(), neighborhood=Disc(2),
                threaded=false)
            @test parent(outN)[:, 1] == wholeaxis[Disc(2)]
            @test parent(outN)[:, 2] == wholeaxis[Disc(2)]
        end

        @testset "a wider field request follows the chunk lines" begin
            # `Index(Local())` summed over a two-ring: a chunk reporting its own
            # numbering, or one whose halo is a single ring, fails here.
            kern(center, rings) = (center[1] + sum(rings[1]; init=0.0),
                Float64(center[2]) + sum(rings[2]; init=0))
            ob = [Float64(p) + sum(r; init=0) for (p, r) in enumerate(nbrings[Ring(2)])]

            B = dggread(counting(path))[:e]
            reads = parent(B).storage.reads
            floor2 = readfloor(chunkplan(B; halo=2))
            needs = (Value(B), Index(Local()))

            empty!(reads)
            ga, gb = mapneighbors(kern, B; needs=needs, neighborhood=Ring(2),
                threaded=false)
            # Each storage chunk is decoded a bounded number of times: once for
            # the cube's own cells, once for the request's `Value` over the same
            # store, both over a two-ring halo.
            @test floor2 <= reads["e"] <= 2 * floor2
            @test parent(ga) == oracles[Ring(2)]
            @test parent(gb) == ob
            @test DD.lookup(ga, Cells) === DD.lookup(B, Cells)

            # The side-effecting form fires once per OWNED cell: a halo cell has
            # an incomplete ring of its own and there is no result to discard.
            hits = zeros(Int, N)
            acc = zeros(Float64, N)
            accb = zeros(Float64, N)
            foreachneighbors(B; needs=needs, neighborhood=Ring(2),
                    threaded=false) do center, rings
                hits[center[2]] += 1
                acc[center[2]] = center[1] + sum(rings[1]; init=0.0)
                accb[center[2]] = Float64(center[2]) + sum(rings[2]; init=0)
            end
            @test all(==(1), hits)
            @test acc == oracles[Ring(2)]
            @test accb == ob
        end

        @testset "a wider halo is still read once per chunk it lands in" begin
            B = dggread(counting(path))[:e]
            plan = chunkplan(B; halo=2)
            onering = chunkplan(B; halo=1)
            @test halowidth(plan) == 2
            # Two rings reach further along the axis than one, so the floor is
            # a different number and not a restatement of the one-ring case.
            @test readfloor(plan) > readfloor(onering)

            reads = parent(B).storage.reads
            empty!(reads)
            dest = zeros(Float64, N)
            mapneighbors!(dest, stencil, B, plan; neighborhood=Ring(2),
                threaded=false)
            @test dest == wholeaxis[Ring(2)]
            @test reads["e"] == readfloor(plan)

            # `Disc(2)` reads the same blocks: the halo is the plan's, not the
            # selector's, so widening the sweep within the plan costs nothing.
            empty!(reads)
            d2 = zeros(Float64, N)
            mapneighbors!(d2, stencil, B, plan; neighborhood=Disc(2),
                threaded=false)
            @test d2 == wholeaxis[Disc(2)]
            @test reads["e"] == readfloor(plan)
        end

        @testset "the halo defaults to the selector's reach" begin
            B = dggread(path)[:e]
            d0 = zeros(Float64, N)
            mapneighbors!(d0, stencil, B; neighborhood=Disc(3), threaded=false)
            d3 = zeros(Float64, N)
            mapneighbors!(d3, stencil, B; neighborhood=Disc(3), halo=3,
                threaded=false)
            d5 = zeros(Float64, N)
            mapneighbors!(d5, stencil, B; neighborhood=Disc(3), halo=5,
                threaded=false)
            @test d0 == d3 == d5 == wholeaxis[Disc(3)]

            # The default is the selector's reach and not some fixed maximum:
            # asking for less is still refused.
            @test_throws ArgumentError mapneighbors!(d0, stencil, B;
                neighborhood=Disc(3), halo=2, threaded=false)

            # `Disc(0)` reaches nothing, and a plan may not carry a haloless
            # chunk, so the default is one ring rather than an error.
            z = zeros(Float64, N)
            mapneighbors!(z, stencil, B; neighborhood=Disc(0), threaded=false)
            @test z == v
        end
    end
end

end # if HAS_ZARR

end # module DGGSIOChunkSweepTests
