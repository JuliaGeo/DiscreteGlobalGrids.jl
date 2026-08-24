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
    chunkplan, foreachchunk, chunkcube, localindices, globalindices, chunkhalo,
    halowidth, nchunks, region, Values, foreachneighbors, localindex

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
            @test reduce(vcat, collect(globalindices(mc)) for mc in plan) == collect(1:N)
            for mc in plan
                h = chunkhalo(mc)
                @test !isempty(h)                      # a global axis has no isolated chunk
                @test issorted(h) && allunique(h)
                @test all(p -> 1 <= p <= N, h)
                @test all(p -> !(p in globalindices(mc)), h)
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
                      length(globalindices(cc)) + length(chunkhalo(cc))
                # Owned and halo partition the block, and owned is contiguous
                # because every halo cell is outside the chunk's own run.
                @test sort(vcat(collect(localindices(cc)), chunkhalo(cc))) ==
                      collect(1:length(cube))
                # The owned rows are the chunk's cells, in the axis's order.
                @test parent(cube)[localindices(cc)] == parent(A)[globalindices(cc)]
                owned += length(localindices(cc))
            end
            @test owned == N
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
            @test reduce(vcat, [collect(globalindices(mc)) for p in pieces for mc in p]) ==
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
    end
end

end # if HAS_ZARR

end # module DGGSIOChunkSweepTests
