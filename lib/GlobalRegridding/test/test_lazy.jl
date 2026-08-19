# Chunk discovery, lazy execution, streaming, and spill storage.
#
# The five laws the `L1`–`L5` testsets below stand for:
# L1 construction is free: a lazy array reads no values and builds no weights.
# L2 locality: a destination chunk reads its connected source chunks, once each.
# L3 residency: the budget bounds what is held at once and changes nothing else.
# L4 plan reuse: each chunk pair is built once and reused by every later read.
# L5 chunking invariance: the answer does not depend on the source's chunking.

import DiskArrays
import DimensionalData as DD

# Source read counter

mutable struct T7Counting{T,N,A<:AbstractArray{T,N}} <: DiskArrays.AbstractDiskArray{T,N}
    data::A
    chunks::DiskArrays.GridChunks{N,NTuple{N,DiskArrays.RegularChunks}}
    reads::Vector{NTuple{N,UnitRange{Int}}}
end

T7Counting(a::AbstractArray{T,N}, cs::Tuple) where {T,N} =
    T7Counting(a, DiskArrays.GridChunks(a, cs), NTuple{N,UnitRange{Int}}[])

Base.size(x::T7Counting) = size(x.data)
DiskArrays.haschunks(::T7Counting) = DiskArrays.Chunked()
DiskArrays.eachchunk(x::T7Counting) = x.chunks
function DiskArrays.readblock!(x::T7Counting, out, r::AbstractUnitRange...)
    push!(x.reads, map(UnitRange{Int}, r))
    out .= view(x.data, r...)
    return out
end

t7_reset!(x::T7Counting) = (empty!(x.reads); x)
t7_spatial(x::T7Counting) = sort([r[1:2] for r in x.reads])

# Weight-build counter

mutable struct T7CountingMethod{M<:AbstractRegriddingMethod} <: AbstractRegriddingMethod
    inner::M
    builds::Int
end

T7CountingMethod(inner) = T7CountingMethod(inner, 0)
GR.support_radius(m::T7CountingMethod, space::RegridSpace) =
    GR.support_radius(m.inner, space)

function build_weights!(coo::WeightCOO, m::T7CountingMethod,
    dst::RegridSpace, dst_inds, src::RegridSpace, src_inds)
    countbuild!(m)
    return build_weights!(coo, m.inner, dst, dst_inds, src, src_inds)
end

# Fixed-radius stencil used to test support discovery.

struct T7RadiusMethod <: AbstractRegriddingMethod
    radius::Float64
end

GR.support_radius(m::T7RadiusMethod, ::RegridSpace) = m.radius

function build_weights!(coo::WeightCOO, m::T7RadiusMethod,
    dst::RegridSpace, dst_inds, src::RegridSpace, src_inds)
    for (j, p) in enumerate(dst_inds)
        centre = cellcentroid(dst, p)
        for (k, q) in enumerate(src_inds)
            US.spherical_distance(centre, cellcentroid(src, q)) <= m.radius || continue
            addweight!(coo, j, k, 1.0)
            adddenom!(coo, j, 1.0)
        end
    end
    return coo
end

t7_plan(method, dst, src; policy = Weighted(0.5), storage = PerChunk(),
    budget = 2^30, chunks = nothing, missingval = nothing) =
    ChunkedPlan(method, policy, dst, src, storage, budget, chunks, missingval)

# Source with a configurable `knownempty` oracle.

mutable struct T8Presence{T,N,A<:AbstractArray{T,N},F} <: DiskArrays.AbstractDiskArray{T,N}
    data::A
    chunks::DiskArrays.GridChunks{N,NTuple{N,DiskArrays.RegularChunks}}
    reads::Vector{NTuple{N,UnitRange{Int}}}
    oracle::F
    enabled::Bool
end

T8Presence(a::AbstractArray{T,N}, cs::Tuple, oracle; enabled::Bool = true) where {T,N} =
    T8Presence(a, DiskArrays.GridChunks(a, cs), NTuple{N,UnitRange{Int}}[], oracle, enabled)

Base.size(x::T8Presence) = size(x.data)
DiskArrays.haschunks(::T8Presence) = DiskArrays.Chunked()
DiskArrays.eachchunk(x::T8Presence) = x.chunks
function DiskArrays.readblock!(x::T8Presence, out, r::AbstractUnitRange...)
    push!(x.reads, map(UnitRange{Int}, r))
    out .= view(x.data, r...)
    return out
end

GR.knownempty(x::T8Presence, nd::Tuple{Vararg{AbstractUnitRange}}) =
    x.enabled && x.oracle(nd)

# Reference connectivity from current pairwise chunk caps.
function t8_pairs(dst, src; radius = 0.0)
    dcaps, scaps = GR.chunkextents(dst), GR.chunkextents(src)
    return [(d, s) for d in eachindex(dcaps) for s in eachindex(scaps)
            if US.spherical_distance(dcaps[d].point, scaps[s].point) <=
               dcaps[d].radius + scaps[s].radius + radius]
end

# The connectivity a flat pairwise cap test gives — the independent reference
# the tree descent is checked against.
function t7_pairwise(dst, dstchunk, src; radius = 0.0)
    dcap = chunktree(dst).caps[dstchunk]
    return [s for (s, scap) in enumerate(chunktree(src).caps)
            if US.spherical_distance(dcap.point, scap.point) <=
               dcap.radius + scap.radius + radius]
end

@testset "Lazy path" begin

    # Source and destination use different chunk layouts.
    srcspace = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
        chunks = (4, 2))
    dstspace = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
        chunks = (8, 2))
    field = collect(reshape(1.0:32.0, 8, 4))

    @testset "discovery" begin
        # Dual-tree discovery agrees with pairwise cap checks.
        for c in 1:nchunks(dstspace)
            @test GR.connectedchunks(dstspace, c, srcspace) == t7_pairwise(dstspace, c, srcspace)
        end

        # Dilation only ever adds pairs.
        wide = GR.connectedchunks(dstspace, 1, srcspace; radius = 0.5)
        @test issubset(GR.connectedchunks(dstspace, 1, srcspace), wide)

        # Chunk extents come back indexed by chunk number, which is what the
        # descent's leaf indices mean.
        @test GR.chunkextents(srcspace) == chunktree(srcspace).caps

        # The batched dual descent finds exactly the union of the per-chunk
        # queries — the form a hierarchical destination tree will use.
        pairs = Tuple{Int,Int}[]
        GR.connectedchunkpairs((d, s) -> push!(pairs, (d, s)), dstspace, srcspace)
        expected = sort([(d, s) for d in 1:nchunks(dstspace)
                         for s in GR.connectedchunks(dstspace, d, srcspace)])
        @test sort(pairs) == expected
    end

    @testset "L1 — construction is free" begin
        source = T7Counting(field, (4, 2))
        plan = t7_plan(ToyDiagonalMethod(), dstspace, srcspace)
        A = LazyRegridArray(source, plan)

        # Construction reads no values and builds no weights.
        @test isempty(source.reads)
        @test GR.nblocks(plan.storage) == 0
        @test size(A) == (32,)
        @test DiskArrays.haschunks(A) isa DiskArrays.Chunked
        # The destination's own chunks are the cell axis's chunks, because this
        # destination's chunks are contiguous runs of cell positions.
        @test collect(DiskArrays.eachchunk(A)) == [(1:16,), (17:32,)]
    end

    @testset "L2 — locality" begin
        source = T7Counting(field, (4, 2))
        plan = t7_plan(ToyDiagonalMethod(), dstspace, srcspace)
        A = LazyRegridArray(source, plan)

        connected = t7_pairwise(dstspace, 1, srcspace)
        values = A[1:16]

        # Read each connected source chunk exactly once.
        @test length(source.reads) == length(connected)
        @test t7_spatial(source) ==
              sort([GR.chunkranges(srcspace, s, (8, 4)) for s in connected])
        @test values == vec(field)[1:16]

        # The second destination chunk reads its own connected set, not the
        # first's.
        t7_reset!(source)
        @test A[17:32] == vec(field)[17:32]
        @test t7_spatial(source) ==
              sort([GR.chunkranges(srcspace, s, (8, 4)) for s in t7_pairwise(dstspace, 2, srcspace)])
    end

    @testset "the wave of concurrent builds is bounded by the weight budget" begin
        # Build waves respect the minimum block-size budget.
        plan = t7_plan(ToyDiagonalMethod(), dstspace, srcspace)
        srcchunks = GR.connectedchunks(dstspace, 1, srcspace)
        srcranges = [GR.chunkranges(srcspace, s, (8, 4)) for s in srcchunks]
        nd = length(cellindices(dstspace, 1))
        @test length(srcchunks) > 1

        tiny = t7_plan(ToyDiagonalMethod(), dstspace, srcspace; budget = 2^10)
        @test GR._wavesize(tiny, nd, srcchunks, srcranges) == 1
        @test GR._wavesize(plan, nd, srcchunks, srcranges) ==
              min(Threads.nthreads(), length(srcchunks))

        # One pair is one build whatever the threads: a wave never spawns for a
        # tile it cannot split.
        @test GR._wavesize(plan, nd, srcchunks[1:1], srcranges[1:1]) == 1
    end

    @testset "L4 — plan reuse" begin
        source = T7Counting(field, (4, 2))
        method = T7CountingMethod(ToyDiagonalMethod())
        plan = t7_plan(method, dstspace, srcspace)
        A = LazyRegridArray(source, plan)

        first_read = A[1:16]
        builds = method.builds
        @test builds == length(t7_pairwise(dstspace, 1, srcspace))
        @test GR.nblocks(plan.storage) == builds

        # Repeated reads reuse weights but reload source data.
        t7_reset!(source)
        @test A[1:16] == first_read
        @test method.builds == builds
        @test length(source.reads) == length(t7_pairwise(dstspace, 1, srcspace))
    end

    @testset "non-spatial slices reuse blocks" begin
        cube = cat(field, 10 .* field, 100 .* field; dims = 3)
        source = T7Counting(cube, (4, 2, 1))
        method = T7CountingMethod(ToyDiagonalMethod())
        plan = t7_plan(method, dstspace, srcspace)
        A = LazyRegridArray(source, plan)
        @test size(A) == (32, 3)

        slice1 = A[1:16, 1:1]
        builds = method.builds
        @test slice1 == reshape(vec(field)[1:16], 16, 1)

        # A second field over the same geometry is a second slice through the
        # same blocks: no build, and nothing of the first slice retained.
        slice2 = A[1:16, 2:2]
        @test method.builds == builds
        @test slice2 == reshape(10 .* vec(field)[1:16], 16, 1)

        # Multiple slices reuse the same spatial blocks.
        three = A[1:16, 1:3]
        @test method.builds == builds
        @test three == hcat(vec(field)[1:16], 10 .* vec(field)[1:16],
            100 .* vec(field)[1:16])

        # New destination tiles build each pair once for all slices.
        every = A[1:32, 1:3]
        @test method.builds == builds + length(t7_pairwise(dstspace, 2, srcspace))
        eager = regrid(cube; to = dstspace, from = srcspace,
            method = ToyDiagonalMethod(), lazy = false)
        @test every == eager
    end

    @testset "dilation" begin
        # A 50° stencil crosses a gap between source chunk extents.
        wide_src = ToyLonLatSpace(8, 1; lon = (-80.0, 80.0), lat = (-5.0, 5.0),
            chunks = (4, 1))
        point_dst = ToyLonLatSpace(1, 1; lon = (-40.0, -30.0), lat = (-5.0, 5.0))
        method = T7RadiusMethod(deg2rad(50))

        @test GR.connectedchunks(point_dst, 1, wide_src) == [1]
        @test GR.connectedchunks(point_dst, 1, wide_src; radius = method.radius) == [1, 2]

        line = collect(reshape(1.0:8.0, 8, 1))
        source = T7Counting(line, (4, 1))
        plan = t7_plan(method, point_dst, wide_src)
        A = LazyRegridArray(source, plan)

        # Dilated discovery includes all five cells within the support radius.
        @test A[1:1] == [3.0]
        @test length(source.reads) == 2
        @test A[1:1] == regrid(line; to = point_dst, from = wide_src, method, lazy = false)
    end

    @testset "chunk equivalence" begin
        source = T7Counting(field, (4, 2))
        # Same field, same spaces, one plan per policy: the lazy answer is the
        # eager answer, bit for bit, for a method whose weights are exact.
        for policy in (Weighted(0.5), Extensive())
            plan = t7_plan(ToyDiagonalMethod(; scale = 2.0), dstspace, srcspace;
                policy)
            A = LazyRegridArray(source, plan)
            eager = regrid(field; to = dstspace, from = srcspace,
                method = ToyDiagonalMethod(; scale = 2.0), missingpolicy = policy,
                lazy = false)
            @test A[1:32] == eager
        end

        # Different source chunking preserves the result.
        coarse = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
            chunks = (8, 1))
        plan = t7_plan(ToyDiagonalMethod(), dstspace, coarse)
        @test LazyRegridArray(T7Counting(field, (8, 1)), plan)[1:32] ==
              regrid(field; to = dstspace, from = srcspace,
            method = ToyDiagonalMethod(), lazy = false)

        # `regrid!` materializes the same values chunk by chunk.
        dest = Vector{Float64}(undef, 32)
        regrid!(dest, field, t7_plan(ToyDiagonalMethod(), dstspace, srcspace))
        @test dest == regrid(field; to = dstspace, from = srcspace,
            method = ToyDiagonalMethod(), lazy = false)
    end

    @testset "LRU" begin
        source = T7Counting(field, (4, 2))
        method = T7CountingMethod(ToyDiagonalMethod())
        storage = PerChunk(1)
        plan = t7_plan(method, dstspace, srcspace; storage)
        A = LazyRegridArray(source, plan)

        reference = regrid(field; to = dstspace, from = srcspace,
            method = ToyDiagonalMethod(), lazy = false)

        # One block resident at a time, whatever the pass costs in rebuilds.
        @test A[1:32] == reference
        @test GR.nblocks(storage) == 1
        builds = method.builds
        @test builds > 1

        # A capacity too small to hold a destination chunk's blocks rebuilds them
        # every pass — the answer is untouched, only the work changes.
        @test A[1:32] == reference
        @test method.builds == 2 * builds
        @test GR.nblocks(storage) == 1

        # Unbounded storage over the same plan builds each pair once and keeps it.
        roomy = T7CountingMethod(ToyDiagonalMethod())
        kept = PerChunk()
        B = LazyRegridArray(source, t7_plan(roomy, dstspace, srcspace; storage = kept))
        @test B[1:32] == reference
        once = roomy.builds
        @test B[1:32] == reference
        @test roomy.builds == once
        @test GR.nblocks(kept) == once
        @test GR.storagebytes(kept) > 0
    end

    @testset "missing data" begin
        # Lazy and eager paths agree on per-chunk missing data.
        holey = copy(field)
        holey[1, 1] = NaN
        source = T7Counting(holey, (4, 2))
        plan = t7_plan(ToyDiagonalMethod(), dstspace, srcspace)
        lazy = LazyRegridArray(source, plan)[1:32]
        eager = regrid(holey; to = dstspace, from = srcspace,
            method = ToyDiagonalMethod(), lazy = false)
        @test isnan(lazy[1])
        @test filter(!isnan, lazy) == filter(!isnan, eager)
    end

    @testset "raster source" begin
        # Raster chunk reads follow array order in both orientations.
        xd = DD.X(-168.75:22.5:168.75)
        yd = DD.Y(-78.75:22.5:78.75)
        dst = ToyLonLatSpace(8, 4; chunks = (8, 2))
        for (shape, dims, chunks) in
            (((16, 8), (xd, yd), (8, 4)), ((8, 16), (yd, xd), (4, 8)))
            raster = DD.DimArray(T7Counting(collect(reshape(1.0:128, shape)), chunks), dims)
            src = RasterGrid(raster)
            plan = t7_plan(Conservative(), dst, src)
            @test LazyRegridArray(raster, plan)[1:32] ==
                  regrid(raster; to = dst, from = src, method = Conservative(),
                lazy = false)
        end
    end

    @testset "budget tiles build one restricted tree each" begin
        # Non-aligned tiles hit the CellCapTree fallback; blocks share it per tile.
        xd = DD.X(-168.75:22.5:168.75)
        yd = DD.Y(-78.75:22.5:78.75)
        dst = ToyLonLatSpace(8, 4; chunks = (8, 2))
        raster = DD.DimArray(T7Counting(collect(reshape(1.0:128, 16, 8)), (8, 4)),
            (xd, yd))
        src = RasterGrid(raster)
        plan = t7_plan(Conservative(), dst, src; chunks = (12,))
        A = LazyRegridArray(raster, plan)
        @test !A.tiling.spacetiled
        ntile = length(A.tiling.runs)
        @test ntile == 3
        b0 = GR.cellcaptree_builds()
        out = A[1:32]
        blocks = GR.residency(A).loads + GR.residency(A).hits
        @test blocks > ntile
        @test GR.cellcaptree_builds() - b0 == ntile
        @test out == regrid(raster; to = dst, from = src, method = Conservative(),
            lazy = false)
    end

    @testset "storage policies" begin
        @test_throws ArgumentError PerChunk(0)
        @test_throws ArgumentError PerChunk(; maxbytes = 0)
    end

    # Budget, spill, nodata, and destination chunking

    reference = regrid(field; to = dstspace, from = srcspace,
        method = ToyDiagonalMethod(), lazy = false)

    @testset "L3 — the budget bounds residency" begin
        # Polar destination bands fan into every source chunk.
        bigsrc = ToyLonLatSpace(64, 32; chunks = (16, 8))
        bigdst = ToyLonLatSpace(64, 32; chunks = (64, 8))
        bigfield = collect(reshape(1.0:2048.0, 64, 32))
        chunkbytes = 8 * 16 * 8
        eager = regrid(bigfield; to = bigdst, from = bigsrc,
            method = ToyDiagonalMethod(), lazy = false)

        function t8_read(budget)
            source = T7Counting(bigfield, (16, 8))
            method = T7CountingMethod(ToyDiagonalMethod())
            plan = t7_plan(method, bigdst, bigsrc; budget,
                storage = PerChunk(; maxbytes = GR.weightbudget(budget)))
            A = LazyRegridArray(source, plan)
            out = Vector{Float64}(undef, ncells(bigdst))
            # One call, so hold-versus-stream is decided once over the whole
            # request and a held chunk has a second tile to be reused by.
            DiskArrays.readblock!(A, out, 1:ncells(bigdst))
            return out, GR.residency(A), plan, method
        end

        held, stats_held, plan_held, _ = t8_read(2^24)
        streamed, stats_streamed, plan_streamed, _ = t8_read(4096)

        # The budget changes what is resident and nothing else.
        @test held == eager
        @test streamed == held

        # Holding reuses chunks; streaming keeps one chunk in flight.
        @test stats_held.hits > 0
        @test stats_streamed.hits == 0
        @test stats_held.loads + stats_held.hits ==
              stats_streamed.loads + stats_streamed.hits
        @test stats_streamed.peakbytes <= 2 * chunkbytes
        @test stats_held.peakbytes >= 8 * chunkbytes

        # The weight cache respects its byte budget.
        @test GR.nblocks(plan_streamed.storage) == 1
        @test GR.nblocks(plan_held.storage) == length(t8_pairs(bigdst, bigsrc))
    end

    @testset "Spilled storage" begin
        dir = mktempdir()
        method = T7CountingMethod(ToyDiagonalMethod())
        # A one-byte memory cache forces disk spill.
        storage = Spilled(dir; maxbytes = 1)
        plan = t7_plan(method, dstspace, srcspace; storage)
        A = LazyRegridArray(T7Counting(field, (4, 2)), plan)

        @test A[1:32] == reference
        built = method.builds
        @test built > 1
        @test length(GR.spilledfiles(storage)) == built
        @test !isempty(readdir(dir))

        # A second pass rebuilds nothing: an evicted block comes back off disk,
        # bit for bit, denominator and all.
        @test A[1:32] == reference
        @test method.builds == built
        @test GR.nblocks(storage) == 1

        # A block without a denominator round-trips as one — `nothing` is a
        # meaningful value, not an empty vector.
        bare = t7_plan(ToyDiagonalMethod(; withdenom = false), dstspace, srcspace;
            storage = Spilled(dir; maxbytes = 1))
        @test LazyRegridArray(T7Counting(field, (4, 2)), bare)[1:32] ==
              regrid(field; to = dstspace, from = srcspace,
            method = ToyDiagonalMethod(; withdenom = false), lazy = false)

        # Separate plans do not share spill files.
        again = T7CountingMethod(ToyDiagonalMethod())
        second = t7_plan(again, dstspace, srcspace; storage = Spilled(dir; maxbytes = 1))
        @test LazyRegridArray(T7Counting(field, (4, 2)), second)[1:32] == reference
        @test again.builds == built
        @test length(readdir(dir)) > built
    end

    @testset "missingval on the lazy path" begin
        # The sentinel is applied per source chunk, so a lazy read must reach the
        # same verdict as the eager one on the same field.
        sentinel = copy(field)
        sentinel[1, 1] = -9999.0
        nanned = copy(field)
        nanned[1, 1] = NaN
        plan = t7_plan(ToyDiagonalMethod(), dstspace, srcspace; missingval = -9999.0)
        lazy = LazyRegridArray(T7Counting(sentinel, (4, 2)), plan)[1:32]
        @test all(isequal.(lazy, regrid(nanned; to = dstspace, from = srcspace,
            method = ToyDiagonalMethod(), lazy = false)))
    end

    @testset "knownempty — per non-spatial chunk" begin
        # Emptiness varies across spatial and time chunks.
        cube = cat(field, 10 .* field, 100 .* field, 1000 .* field; dims = 3)
        cube[1:4, 1:2, 1:2] .= NaN
        cube[5:8, 3:4, :] .= NaN
        oracle(nd) = (nd[1] == 1:4 && nd[2] == 1:2 && nd[3] == 1:2) ||
                     (nd[1] == 5:8 && nd[2] == 3:4)

        function t8_cube(policy, enabled)
            source = T8Presence(cube, (4, 2, 2), oracle; enabled)
            method = T7CountingMethod(ToyDiagonalMethod())
            plan = t7_plan(method, dstspace, srcspace; policy)
            A = LazyRegridArray(source, plan)
            out = Matrix{Float64}(undef, 32, 4)
            DiskArrays.readblock!(A, out, 1:32, 1:4)
            return out, source, method, GR.residency(A)
        end

        for policy in (Weighted(0.5), Extensive())
            skipping, source, method, stats = t8_cube(policy, true)
            plain, plainsource, plainmethod, plainstats = t8_cube(policy, false)

            # Skipping preserves the result.
            @test all(isequal.(skipping, plain))
            # And visible in the reads — an empty combination is never asked for.
            @test !any(oracle, source.reads)
            @test any(oracle, plainsource.reads)
            @test length(source.reads) < length(plainsource.reads)
            @test stats.skipped > 0
        end

        # Reference-using policies retain weights for all-times-empty chunks.
        _, weighted_source, weighted_method, weighted_stats = t8_cube(Weighted(0.5), true)
        _, _, plainmethod, _ = t8_cube(Weighted(0.5), false)
        @test weighted_stats.dropped == 0
        @test weighted_method.builds == plainmethod.builds

        # Policies without references drop all-times-empty chunks.
        _, ext_source, ext_method, ext_stats = t8_cube(Extensive(), true)
        _, _, ext_plain, _ = t8_cube(Extensive(), false)
        @test ext_stats.dropped > 0
        @test ext_method.builds < ext_plain.builds
        @test !any(r -> r[1] == 5:8 && r[2] == 3:4, ext_source.reads)
    end

    @testset "knownempty — dropping is licensed, not assumed" begin
        # Empty-chunk reference weights still affect coverage thresholds.
        wide_src = ToyLonLatSpace(8, 1; lon = (-80.0, 80.0), lat = (-5.0, 5.0),
            chunks = (4, 1))
        point_dst = ToyLonLatSpace(1, 1; lon = (20.0, 30.0), lat = (-5.0, 5.0))
        method = T7RadiusMethod(deg2rad(50))
        line = collect(reshape(1.0:8.0, 8, 1))
        line[5:8, 1] .= NaN
        empty2(nd) = nd[1] == 5:8

        function t8_line(policy, enabled)
            source = T8Presence(line, (4, 1), empty2; enabled)
            plan = t7_plan(method, point_dst, wide_src; policy)
            A = LazyRegridArray(source, plan)
            return A[1:1], GR.residency(A), source
        end

        blanked, wstats, wsource = t8_line(Weighted(0.5), true)
        plain, _, _ = t8_line(Weighted(0.5), false)
        @test isnan(only(blanked))
        @test all(isequal.(blanked, plain))
        @test wstats.dropped == 0
        @test wstats.skipped > 0
        # The read of the empty chunk is skipped even though its weights are kept.
        @test !any(empty2, wsource.reads)

        raw, estats, _ = t8_line(Extensive(), true)
        rawplain, _, _ = t8_line(Extensive(), false)
        @test all(isequal.(raw, rawplain))
        @test estats.dropped > 0
    end

    @testset "destination chunking" begin
        # Explicit `chunks` tiles a destination with one space chunk.
        whole = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0))
        undivided = LazyRegridArray(T7Counting(field, (4, 2)),
            t7_plan(ToyDiagonalMethod(), whole, srcspace))
        @test length(collect(DiskArrays.eachchunk(undivided))) == 1

        divided = LazyRegridArray(T7Counting(field, (4, 2)),
            t7_plan(ToyDiagonalMethod(), whole, srcspace; chunks = (8,)))
        @test collect(DiskArrays.eachchunk(divided)) == [(1:8,), (9:16,), (17:24,), (25:32,)]

        expected = regrid(field; to = whole, from = srcspace,
            method = ToyDiagonalMethod(), lazy = false)
        @test divided[1:32] == expected
        # Each declared destination chunk computes only its cells.
        @test divided[9:16] == expected[9:16]

        # `regrid!` drives the chunk loop, so the same tiling is what gets
        # materialized one piece at a time.
        dest = Vector{Float64}(undef, 32)
        regrid!(dest, field, t7_plan(ToyDiagonalMethod(), whole, srcspace; chunks = (8,)))
        @test dest == expected

        # A GridChunks says the same thing, and a mismatched one is caught.
        gridded = t7_plan(ToyDiagonalMethod(), whole, srcspace;
            chunks = DiskArrays.GridChunks((32,), (8,)))
        @test LazyRegridArray(T7Counting(field, (4, 2)), gridded)[1:32] == expected
        @test_throws ArgumentError LazyRegridArray(T7Counting(field, (4, 2)),
            t7_plan(ToyDiagonalMethod(), whole, srcspace; chunks = (8, 2)))
    end

    @testset "outer parallelism wins" begin
        # Top level threads exactly when the session has threads.
        top = Threads.nthreads() > 1 ? GOCore.True() : GOCore.False()
        @test GR._innerthreaded() === top
        # Inside a declared outer wave the inner build stays serial, spawned
        # tasks included, because they inherit the scope.
        Base.ScopedValues.@with GR.OUTER_PARALLEL => true begin
            @test GR._innerthreaded() === GOCore.False()
            @test fetch(Threads.@spawn GR._innerthreaded()) === GOCore.False()
        end
        @test GR._innerthreaded() === top
    end

    @testset "L5 — chunking invariance" begin
        # Conservative results are invariant to incompatible source chunkings.
        xd = DD.X(-168.75:22.5:168.75)
        yd = DD.Y(-78.75:22.5:78.75)
        dst = ToyLonLatSpace(8, 4; chunks = (8, 2))
        values = collect(reshape(1.0:128, 16, 8))

        results = map(((8, 4), (4, 8), (16, 8))) do cs
            raster = DD.DimArray(T7Counting(copy(values), cs), (xd, yd))
            src = RasterGrid(raster)
            storage = PerChunk()
            A = LazyRegridArray(raster, t7_plan(Conservative(), dst, src; storage))
            out = A[1:32]
            # Connectivity is whatever the two spaces' caps say it is *now*, so
            # this stays honest as the extents are tightened.
            @test GR.nblocks(storage) == length(t8_pairs(dst, src))
            return out
        end

        @test results[1] ≈ results[2] rtol = 1e-12
        @test results[1] ≈ results[3] rtol = 1e-12

        # And all of them are the eager whole-domain answer.
        whole = DD.DimArray(copy(values), (xd, yd))
        @test results[1] ≈ regrid(whole; to = dst, from = RasterGrid(whole),
            method = Conservative(), lazy = false) rtol = 1e-12
    end
end
