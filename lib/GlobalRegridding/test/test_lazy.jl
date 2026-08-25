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

# The shared relation oracles. `test_chunkgraph.jl` runs first and defines the
# module; this file reuses it rather than re-spelling either definition, the
# same way it reuses that file's `G4ProbeSpace`.
using .ChunkGraphOracles: graph_pairs, demanded_pairs

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
GR.supportradius(m::T7CountingMethod, space::RegridSpace) =
    GR.supportradius(m.inner, space)

function buildweights!(coo::WeightCOO, m::T7CountingMethod,
    dst::RegridSpace, dst_inds, src::RegridSpace, src_inds)
    countbuild!(m)
    return buildweights!(coo, m.inner, dst, dst_inds, src, src_inds)
end

# Fixed-radius stencil used to test support discovery.

struct T7RadiusMethod <: AbstractRegriddingMethod
    radius::Float64
end

GR.supportradius(m::T7RadiusMethod, ::RegridSpace) = m.radius

function buildweights!(coo::WeightCOO, m::T7RadiusMethod,
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

# A point method with `BarycentricPoint`'s stencils and the default support
# radius of zero, so its relation is cap overlap and nothing else.

struct T9NoReach <: AbstractRegriddingMethod end

GR.outputsampling(::T9NoReach) = DD.Lookups.Points()
GR.sampler(::T9NoReach, space::RegridSpace) = GR.sampler(BarycentricPoint(), space)

# The source chunks a destination tile's stencils name, taken from an operator
# built over the whole domain by the chunk-pair builder, so nothing a tile build
# does takes part in the answer.
function t9_owners(whole::Matrix{Float64}, space::RegridSpace, dinds)
    out = Set{Int}()
    for i in dinds, j in axes(whole, 2)
        iszero(whole[Int(i), j]) || push!(out, Int(GR.chunkat(space, j)))
    end
    return sort!(collect(out))
end

t7_plan(method, dst, src; policy = Weighted(0.5), storage = PerChunk(),
    budget = 2^30, chunks = nothing, missingval = nothing) =
    ChunkedPlan(method, policy, dst, src, storage, budget, chunks, missingval)

# A relation carrying nothing but the caps a wave costs against. `_wavesize` and
# `_blockcosts!` read per-chunk geometry off the plan's relation rather than off
# private vectors, so a test that wants a contrived pair of cap sets builds the
# relation that carries them. The CSR is deliberately edgeless: the wave
# estimator reads the caps and the caller's own `srcchunks`, never the rows.
t7_capgraph(dstcaps, srccaps) =
    GR.ChunkDependencyGraph(GR.DependencyIdentity(),
        fill(1, length(dstcaps) + 1), Int32[],
        fill(1, length(srccaps) + 1), Int32[];
        dstcaps = collect(dstcaps), srccaps = collect(srccaps))

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
# the tree descent is checked against. It reads `chunkextents`, the one way to
# obtain a space's chunk caps.
function t7_pairwise(dst, dstchunk, src; radius = 0.0)
    dcap = GR.chunkextents(dst)[dstchunk]
    return [s for (s, scap) in enumerate(GR.chunkextents(src))
            if US.spherical_distance(dcap.point, scap.point) <=
               dcap.radius + scap.radius + radius]
end

# One destination chunk's source chunks, read the way the lazy executor reads
# them: off the relation the plan owns. These are the rows the executor takes,
# so an assertion about them is an assertion about the read.
t7_sources(plan::ChunkedPlan, d::Integer) =
    Int.(GR.sourcesof(GR.dependencies(plan), Int(d)))

@testset "Lazy path" begin

    # Source and destination use different chunk layouts.
    srcspace = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
        chunks = (4, 2))
    dstspace = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
        chunks = (8, 2))
    field = collect(reshape(1.0:32.0, 8, 4))

    @testset "discovery" begin
        # Every claim below is about the relation a plan owns, because that
        # is the only relation there is: one `candidatechunks!` query per
        # destination cap defines every edge.
        plan = t7_plan(ToyDiagonalMethod(), dstspace, srcspace)
        graph = GR.dependencies(plan)

        # The generic packed index agrees with pairwise cap checks.
        for c in 1:nchunks(dstspace)
            @test t7_sources(plan, c) == t7_pairwise(dstspace, c, srcspace)
        end

        # Dilation only ever adds pairs. The radius reaches the relation
        # through the plan's method, which is the only way it can.
        wideplan = t7_plan(T7RadiusMethod(0.5), dstspace, srcspace)
        @test GR.dependency_radius(GR.dependencies(wideplan)) == 0.5
        wide = t7_sources(wideplan, 1)
        @test issubset(t7_sources(plan, 1), wide)

        # Chunk extents come back indexed by chunk number, which is what the
        # descent's leaf indices mean — and the relation carries the very
        # vectors it was built from, so there is one copy of each, not two.
        @test GR.sourceextents(graph) == GR.chunkextents(srcspace)
        @test GR.destinationextents(graph) == GR.chunkextents(dstspace)

        # The relation's rows ARE the per-destination-cap queries: one
        # `chunkindex(src)`, one `candidatechunks!` per `chunkextents(dst)`
        # entry. `demanded_pairs` replays exactly that loop, so with no
        # `refine` the two must be equal, not merely nested.
        @test graph_pairs(graph) == demanded_pairs(dstspace, srcspace)

        # The generic index is GeometryOps' packed R-tree. Its leaf boxes are
        # outward-rounded bounds of each cap, including tiny, polar,
        # antimeridian and whole-sphere cases.
        index = GR.chunkindex(srcspace)
        @test index isa GR.FlexibleRTrees.RTree
        srccaps = GR.chunkextents(srcspace)
        srcboxes = map(cap -> convert(GR.Extents.Extent, cap), srccaps)
        for capacity in (2, 3, 16)
            packed = GR.FlexibleRTrees.RTree(GR.FlexibleRTrees.HPR(), srccaps;
                extents = srcboxes, nodecapacity = capacity)
            for c in 1:nchunks(dstspace)
                dcap = GR.chunkextents(dstspace)[c]
                @test GR.candidatechunks!(Int[], packed, dcap) ==
                      t7_pairwise(dstspace, c, srcspace)
            end
        end
        capcases = (
            SphericalCap(toy_point(0, 90), 1e-12),
            SphericalCap(toy_point(180, 0), 0.3),
            SphericalCap(toy_point(25, -40), 1.2),
            SphericalCap(toy_point(0, 0), Float64(pi)),
        )
        for cap in capcases
            box = convert(GR.Extents.Extent, cap)
            @test keys(box) == (:X, :Y, :Z)
            centre = cap.point
            seed = abs(centre[3]) < 0.9 ? USPoint(0.0, 0.0, 1.0) : USPoint(1.0, 0.0, 0.0)
            u = LinearAlgebra.normalize(seed - LinearAlgebra.dot(seed, centre) * centre)
            v = USPoint(centre[2] * u[3] - centre[3] * u[2],
                centre[3] * u[1] - centre[1] * u[3],
                centre[1] * u[2] - centre[2] * u[1])
            r = min(cap.radius, Float64(pi))
            for θ in range(0.0, 2pi; length = 65)
                p = cos(r) * centre + sin(r) * (cos(θ) * u + sin(θ) * v)
                @test box.X[1] <= p[1] <= box.X[2]
                @test box.Y[1] <= p[2] <= box.Y[2]
                @test box.Z[1] <= p[3] <= box.Z[2]
            end
        end
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
        # destination's chunks are contiguous runs of cell indices.
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
        srcchunks = t7_sources(plan, 1)
        srcranges = [GR.chunkranges(srcspace, s, (8, 4)) for s in srcchunks]
        nd = length(ownedindices(dstspace, 1))
        # The tile is destination chunk 1, so its rows are `[1]`; the caps come
        # off the plan's own relation.
        graph = GR.dependencies(plan)
        rows = [1]
        nt = Threads.nthreads()
        @test length(srcchunks) > 1

        # Under a declared outer wave, inner threading is already off, so the
        # wave is the only parallelism left and the budget alone bounds it.
        Base.ScopedValues.@with GR.OUTER_PARALLEL => true begin
            tiny = t7_plan(ToyDiagonalMethod(), dstspace, srcspace; budget = 2^10)
            @test GR._wavesize(tiny, nd, srcchunks, srcranges, graph, rows) == 1
            @test GR._wavesize(plan, nd, srcchunks, srcranges, graph, rows) ==
                  min(nt, length(srcchunks))

            # One pair is one build whatever the threads: a wave never spawns
            # for a tile it cannot split.
            @test GR._wavesize(plan, nd, srcchunks[1:1], srcranges[1:1], graph, rows) == 1
        end
    end

    @testset "a wave that loses a task still waits for the rest" begin
        # The failing method is built from the chunk it must fail on, so the
        # rows are read off a plan over the same pair at the same radius first.
        srcchunks = t7_sources(t7_plan(ToyDiagonalMethod(), dstspace, srcspace), 1)
        j = min(3, length(srcchunks))
        @test j > 1
        dinds = ownedindices(dstspace, 1)
        # The first chunk of the wave throws; the others sleep. A `_fillwave!`
        # that raises on the first `fetch` and abandons the rest would return
        # with those tasks still running against a plan the caller is done with.
        bad = Int(first(ownedindices(srcspace, srcchunks[1])))
        method = WaveFailMethod(bad, 0.25)
        plan = t7_plan(method, dstspace, srcspace)
        # ...and the plan under test really does hold those rows.
        @test t7_sources(plan, 1) == srcchunks
        wave = GR.CachedBlock[]
        @test_throws Exception GR._fillwave!(wave, plan, 1, srcchunks, 1, j,
            dinds, GR.TileCells(plan.dst_space, dinds))
        @test method.finished[] == j - 1
    end

    @testset "at top level the wave is weighed against inner threading" begin
        plan = t7_plan(ToyDiagonalMethod(), dstspace, srcspace)
        srcchunks = t7_sources(plan, 1)
        srcranges = [GR.chunkranges(srcspace, s, (8, 4)) for s in srcchunks]
        nd = length(ownedindices(dstspace, 1))
        nt = Threads.nthreads()

        # Costs weigh the chunk by how much of it the tile can reach, so equal
        # chunk sizes do not make an even wave.
        wide = [SphericalCap(toy_point(0, 0), 1.0)]
        near = SphericalCap(toy_point(0, 0), 0.2)
        far = SphericalCap(toy_point(180, 0), 0.2)
        n = max(nt, 4)
        chunks = collect(1:n)
        rows = [1]
        ranges = [srcranges[1] for _ in 1:n]
        @test GR._blockcosts!(zeros(2), [1, 2], ranges[1:2],
                  t7_capgraph(wide, [near, far]), rows) ==
              [Float64(prod(map(length, ranges[1]))), 0.0]

        if nt > 1
            # Every chunk sits inside the tile: the wave is worth its width.
            @test GR._wavesize(plan, nd, chunks, ranges,
                      t7_capgraph(wide, fill(near, n)), rows) == min(nt, n)
            # One chunk carries all the work and the rest meet the tile only at
            # its extent: the wave cannot beat the threading it would suppress.
            lopsided = [i == 1 ? near : far for i in 1:n]
            @test GR._wavesize(plan, nd, chunks, ranges,
                      t7_capgraph(wide, lopsided), rows) == 1
            # Nothing to weigh at all still falls back to one build at a time.
            @test GR._wavesize(plan, nd, chunks, ranges,
                      t7_capgraph(wide, fill(far, n)), rows) == 1
        else
            @test GR._wavesize(plan, nd, chunks, ranges,
                      t7_capgraph(wide, fill(near, n)), rows) == 1
        end
    end

    @testset "cap overlap areas" begin
        area(r) = 2pi * (1 - cos(r))
        # Disjoint, nested, and identical caps are exact.
        @test GR._capoverlap(0.3, 0.3, 1.0) == 0.0
        @test GR._capoverlap(1.0, 0.2, 0.5) ≈ area(0.2)
        @test GR._capoverlap(0.4, 0.4, 0.0) ≈ area(0.4)
        # A cap against the whole sphere is the cap.
        @test GR._capoverlap(Float64(GR._WHOLE_SPHERE.radius), 0.7, 2.0) ≈ area(0.7)
        # Touching caps meet in nothing, and half-covered ones in half.
        @test GR._capoverlap(0.3, 0.5, 0.8) ≈ 0.0 atol = 1e-12
        @test GR._capoverlap(0.5, 0.5, 0.5) ≈ 0.307515 rtol = 1e-5
        # Symmetric in its two caps, and monotone as they separate.
        @test GR._capoverlap(0.9, 0.4, 0.7) ≈ GR._capoverlap(0.4, 0.9, 0.7)
        @test GR._capoverlap(0.9, 0.4, 0.7) > GR._capoverlap(0.9, 0.4, 0.9)
        # Caps wider than a quarter turn go through their complements.
        @test GR._capoverlap(2.0, 0.5, 1.8) ≈ 0.585781 rtol = 1e-5
        @test GR._capoverlap(1.7, 1.7, 3.0) ≈ 1.619108 rtol = 1e-5
        @test GR._capoverlap(2.5, 2.0, 1.0) ≈ 8.402553 rtol = 1e-5
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

        line = collect(reshape(1.0:8.0, 8, 1))
        source = T7Counting(line, (4, 1))
        plan = t7_plan(method, point_dst, wide_src)
        A = LazyRegridArray(source, plan)

        # The support radius widens the plan's own relation — the rows the read
        # below will take — and it reaches it only through the plan's method.
        @test t7_sources(t7_plan(ToyDiagonalMethod(), point_dst, wide_src), 1) == [1]
        @test t7_sources(plan, 1) == [1, 2]

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
        # Non-aligned tiles hit the fallback tree; blocks share it per tile.
        xd = DD.X(-168.75:22.5:168.75)
        yd = DD.Y(-78.75:22.5:78.75)
        dst = CountingSpace(ToyLonLatSpace(8, 4; chunks = (8, 2)))
        raster = DD.DimArray(T7Counting(collect(reshape(1.0:128, 16, 8)), (8, 4)),
            (xd, yd))
        src = RasterGrid(raster)
        plan = t7_plan(Conservative(), dst, src; chunks = (12,))
        A = LazyRegridArray(raster, plan)
        @test !A.tiling.spacetiled
        ntile = length(A.tiling.runs)
        @test ntile == 3
        b0 = dst.builds[]
        out = A[1:32]
        blocks = GR.residency(A).loads + GR.residency(A).hits
        @test blocks > ntile
        @test dst.builds[] - b0 == ntile
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

    @testset "caches and spills carry no reference state of their own" begin
        method = T7CountingMethod(ToyDiagonalMethod())
        plan = t7_plan(method, dstspace, srcspace)
        dinds = ownedindices(dstspace, 1)

        built = GR.blockfor(plan, (1, 1), dinds)
        @test method.builds == 1
        @test built.block.reference === built.block.denom

        # A hit hands back the entry, so it hands back the same reference object.
        # A cache that kept a reference of its own would answer an equal copy.
        hit = GR.blockfor(plan, (1, 1), dinds)
        @test hit === built
        @test hit.block.reference === built.block.reference
        @test method.builds == 1

        # A tile's cached blocks reference their own blocks too.
        tiled = GR.CachedTile(GR.TileWeights([1], WeightBlock[built.block]))
        @test tiled.entries[1].block.reference === built.block.reference

        # The spill format stores weights and the optional denominator, and
        # nothing else: a denominated block and a block with none, over the same
        # weights, differ on disk by exactly the denominator.
        dir = mktempdir()
        inds = ownedindices(srcspace, 1)
        m = length(inds)
        coo = WeightCOO(m)
        buildweights!(coo, ToyDiagonalMethod(; scale = 2.0), srcspace, inds, srcspace, inds)
        denominated = WeightBlock(coo, m, m)
        bare = WeightCOO(m)
        buildweights!(bare, ToyDiagonalMethod(; scale = 2.0, withdenom = false),
            srcspace, inds, srcspace, inds)
        point = WeightBlock(bare, m, m)
        @test point.weights == denominated.weights

        dpath = GR.writeblockfile(joinpath(dir, "denominated.blk"), denominated)
        ppath = GR.writeblockfile(joinpath(dir, "point.blk"), point)
        @test filesize(dpath) - filesize(ppath) == 8 * m

        # Both come back with the reference reconstructed, aliased to the
        # denominator where there is one and recomputed where there is not.
        dback = GR.readblockfile(dpath)
        @test dback.weights == denominated.weights
        @test dback.denom == denominated.denom
        @test dback.reference === dback.denom
        @test GR._blockbytes(dback) == GR._blockbytes(denominated)
        pback = GR.readblockfile(ppath)
        @test pback.weights == point.weights
        @test pback.denom === nothing
        @test pback.reference == point.reference
        @test GR._blockbytes(pback) == GR._blockbytes(point)

        # A spilled tile's blocks are reconstructed the same way.
        tpath = GR.writetilefile(joinpath(dir, "tile.tile"),
            GR.TileWeights([1, 2], WeightBlock[denominated, point]))
        tback = GR.readtilefile(tpath)
        @test tback.sourcechunks == [1, 2]
        @test tback.blocks[1].reference === tback.blocks[1].denom
        @test tback.blocks[2].reference == point.reference
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

    # ----------------------------------------------------------------------
    # The lazy executor is driven by the plan's dependency rows
    # ----------------------------------------------------------------------

    @testset "a lazy read takes its sources from the plan's relation" begin
        plan = t7_plan(ToyDiagonalMethod(), dstspace, srcspace)
        graph = GR.dependencies(plan)
        source = T7Counting(field, (4, 2))
        A = LazyRegridArray(source, plan)

        # One relation, one object: the array's is its plan's, identically, and
        # two arrays over one plan share it. Whatever schedules against the plan
        # is looking at exactly what the read consults.
        @test GR.dependencies(A) === graph
        @test GR.dependencies(A) === GR.dependencies(plan)
        @test GR.dependencies(LazyRegridArray(T7Counting(field, (4, 2)), plan)) === graph

        # This destination's chunks are its tiles, so a tile's sources are that
        # chunk's row, exactly — not a superset of it, and not a re-derivation.
        @test A.tiling.spacetiled
        out = Int[]
        for t in 1:nchunks(dstspace)
            @test GR._connectedsource!(out, A, t) == Int.(GR.sourcesof(graph, t))
        end
    end

    @testset "a derived tile takes the sorted union of its rows" begin
        # A tile smaller or larger than a destination chunk is a *derived* tile:
        # it spans a set of destination chunks and its sources are the union of
        # their rows, ascending and duplicate-free whatever order the rows are
        # visited in.
        dst = ToyLonLatSpace(12, 6; chunks = (3, 2))
        src = ToyLonLatSpace(8, 4; chunks = (2, 1))
        data = collect(reshape(1.0:32.0, 8, 4))
        plan = t7_plan(ToyDiagonalMethod(), dst, src; chunks = (7, ))
        graph = GR.dependencies(plan)
        A = LazyRegridArray(T7Counting(data, (4, 2)), plan)
        @test !A.tiling.spacetiled

        # At least one tile really does span several destination chunks, or the
        # union branch below is never exercised.
        @test any(t -> length(A.tiling.chunksof[t]) > 1, eachindex(A.tiling.runs))

        out = Int[]
        for t in eachindex(A.tiling.runs)
            rows = A.tiling.chunksof[t]
            union = sort!(unique!(reduce(vcat,
                [Int.(GR.sourcesof(graph, d)) for d in rows])))
            @test GR._connectedsource!(out, A, t) == union
            @test issorted(out) && allunique(out)
            # Order-independence: the same rows in any order give the same set.
            @test GR._unionrows!(Int[], graph, reverse(rows)) == union
        end
    end

    @testset "a lazy read performs no dependency discovery" begin
        # `G4ProbeSpace` counts every `chunkextents` call — the destination
        # caps, the identity stamp and the vector the generic `chunkindex`
        # packs all come through it — so the counter rises exactly when a
        # relation is being derived from the space. It is defined in
        # `test_chunkgraph.jl`, which runs first; reuse rather than a second
        # copy is deliberate.
        probe = G4ProbeSpace(srcspace)
        plan = t7_plan(ToyDiagonalMethod(), dstspace, probe)
        built = probe.queries
        @test built > 0                                  # the plan did derive one

        source = T7Counting(field, (4, 2))
        A = LazyRegridArray(source, plan)
        @test probe.queries == built                     # and the array derived none

        expected = vec(field)
        @test A[1:32] == expected
        @test A[1:16] == expected[1:16]
        @test A[17:32] == expected[17:32]
        @test A[5:20] == expected[5:20]
        @test collect(A) == expected
        # Not one question to the source space through the whole of that.
        @test probe.queries == built

        # Structurally: the array holds the relation and nothing else that could
        # answer a spatial query. No source index, no cap vectors of its own.
        @test :graph in fieldnames(typeof(A))
        @test fieldtype(typeof(A), :graph) <: GR.ChunkDependencyGraph
        @test !any(in(fieldnames(typeof(A))), (:srcindex, :dstcaps, :srccaps, :radius))
        @test !any(T -> T <: AbstractVector{<:SphericalCap}, fieldtypes(typeof(A)))
    end

    @testset "wave costing reads the relation's extents, and copies none" begin
        # The per-chunk caps live on the relation. The point is not only
        # that they are reachable, but that there is ONE of each: a plan, its
        # array and every task in a wave read the same vectors.
        plan = t7_plan(ToyDiagonalMethod(), dstspace, srcspace)
        graph = GR.dependencies(plan)
        A = LazyRegridArray(T7Counting(field, (4, 2)), plan)

        @test GR.hasextents(graph)
        @test GR.destinationextents(graph) == GR.chunkextents(dstspace)
        @test GR.sourceextents(graph) == GR.chunkextents(srcspace)
        @test GR.destinationextent(graph, 2) === GR.destinationextents(graph)[2]
        @test GR.sourceextent(graph, 3) === GR.sourceextents(graph)[3]
        # Shared by reference, not copied, through every derived object.
        @test GR.sourceextents(GR.dependencies(A)) === GR.sourceextents(graph)
        @test GR.destinationextents(GR.restrict(graph, [1])) ===
              GR.destinationextents(graph)

        # A relation assembled from bare CSR arrays carries none, says so, and
        # is refused by the lazy array rather than silently degrading.
        bare = GR.ChunkDependencyGraph(GR.DependencyIdentity(), [1, 1], Int32[],
            [1, 1], Int32[])
        @test !GR.hasextents(bare)
        @test_throws ArgumentError GR.sourceextents(bare)
        @test_throws ArgumentError GR.destinationextent(bare, 1)
    end

    @testset "a plan with no relation cannot back a lazy read" begin
        plan = ChunkedPlan(ToyDiagonalMethod(), Weighted(0.5), dstspace, srcspace;
            dependencies = false)
        @test GR.dependencies(plan) === nothing
        err = try
            LazyRegridArray(T7Counting(field, (4, 2)), plan)
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("dependencies = false", err.msg)
    end

    @testset "a point tile reads exactly the chunks its stencils name" begin
        # The counting point fixture of `test_interpolation.jl`: stencils that
        # bracket each destination across source-chunk seams, on a destination
        # whose two chunks are its tiles, against an 18-chunk source.
        src = ToyLonLatSpace(36, 18; chunks = (6, 6))
        dst = ToyLonLatSpace(20, 10; lon = (-19.0, 21.0), lat = (-20.0, 20.0),
            chunks = (20, 5))
        data = collect(reshape(1.0:648.0, 36, 18))
        source = T7Counting(data, (6, 6))
        plan = t7_plan(T5PlaceCount(), dst, src)
        graph = GR.dependencies(plan)
        A = LazyRegridArray(source, plan)
        @test A.tiling.spacetiled && nchunks(dst) == 2

        eager = regrid(data; to = dst, from = src, method = T5PlaceCount(),
            lazy = false)

        for t in 1:nchunks(dst)
            dinds = ownedindices(dst, t)
            t7_reset!(source)
            @test A[first(dinds):last(dinds)] ≈ eager[first(dinds):last(dinds)]

            manifest = plan.storage.tiles[t].sourcechunks
            row = Int.(GR.sourcesof(graph, t))

            # The manifest is the tile's own stencils, and the read is the
            # manifest: no row chunk with no block was loaded, and no manifest
            # chunk was skipped.
            @test manifest == t5_owners(dst, dinds, src)
            @test length(source.reads) == length(manifest)
            @test t7_spatial(source) ==
                  sort([GR.chunkranges(src, s, (36, 18)) for s in manifest])

            # The relation's row is a superset, and a strict one, so reading it
            # rather than the manifest would show up in the count above.
            @test issubset(manifest, row)
            @test length(manifest) < length(row)
        end

        # Reading the whole destination in one call shares one source hold
        # across the tiles, and still touches nothing outside their manifests.
        t7_reset!(source)
        @test A[1:Int(ncells(dst))] ≈ eager
        held = Set(GR.chunkranges(src, s, (36, 18))
                   for t in 1:nchunks(dst) for s in plan.storage.tiles[t].sourcechunks)
        rowheld = Set(GR.chunkranges(src, s, (36, 18))
                      for t in 1:nchunks(dst) for s in Int.(GR.sourcesof(graph, t)))
        @test !isempty(source.reads)
        @test all(in(held), t7_spatial(source))
        @test length(held) < length(rowheld)
    end

    @testset "a point stencil reaching past cap overlap needs a declared radius" begin
        # A destination lattice twenty times finer than the source, wholly
        # inside one source cell and just past its sample site, so its stencils
        # bracket sites in three further source chunks. Each source cell is its
        # own chunk, so the caps are tight and the destination meets one.
        src = ToyLonLatSpace(36, 18; chunks = (1, 1))
        dst = ToyLonLatSpace(2, 2; lon = (5.0, 6.0), lat = (5.0, 6.0),
            chunks = (2, 2))
        data = collect(reshape(1.0:648.0, 36, 18))
        ndst = Int(ncells(dst))
        @test nchunks(dst) == 1

        owner = Int(GR.chunkat(src, cellat(src, cellcentroid(dst, 1))))
        @test all(Int(GR.chunkat(src, cellat(src, cellcentroid(dst, i)))) == owner
                  for i in 1:ndst)

        weights = GR.tileweights(BarycentricPoint(), GR.TileCells(dst, 1:ndst),
            1:ndst, src, GR.sampler(BarycentricPoint(), src))
        @test length(weights.sourcechunks) == 4
        @test owner in weights.sourcechunks

        # Cap overlap alone names only the chunk the destination sits in: the
        # other three own sites whose own cells the destination never touches,
        # so an intersection with a radius-free row would drop their weights.
        caps, dcap = GR.chunkextents(src), GR.chunkextents(dst)[1]
        @test [s for s in weights.sourcechunks
               if US.spherical_distance(dcap.point, caps[s].point) <=
                  dcap.radius + caps[s].radius] == [owner]

        # The declared reach is what puts all four in the relation, and the read
        # is all four.
        plan = t7_plan(BarycentricPoint(), dst, src)
        @test supportradius(BarycentricPoint(), src) ≈ deg2rad(10.0)
        @test issubset(weights.sourcechunks, Int.(GR.sourcesof(GR.dependencies(plan), 1)))
        source = T7Counting(data, (1, 1))
        A = LazyRegridArray(source, plan)
        @test A[1:ndst] ≈ regrid(data; to = dst, from = src,
            method = BarycentricPoint(), lazy = false)
        @test t7_spatial(source) ==
              sort([GR.chunkranges(src, s, (36, 18)) for s in weights.sourcechunks])

        # Declaring no reach leaves three of them out of the relation. The read
        # refuses, naming the hook and the method, rather than answering with
        # three weights silently dropped.
        bare = t7_plan(T9NoReach(), dst, src)
        @test supportradius(T9NoReach(), src) == 0.0
        @test !issubset(weights.sourcechunks,
            Int.(GR.sourcesof(GR.dependencies(bare), 1)))
        err = try
            LazyRegridArray(T7Counting(data, (1, 1)), bare)[1:ndst]
            nothing
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("supportradius", err.msg)
        @test occursin("T9NoReach", err.msg)
    end

    @testset "a point tile's manifest tracks the source's chunking" begin
        # A chunked raster source under two chunkings, and a finer destination
        # lattice inside its hull straddling both chunk seams.
        xd = DD.X(-177.5:5.0:177.5)
        yd = DD.Y(-87.5:5.0:87.5)
        data = collect(reshape(1.0:2592.0, 72, 36))
        dst = ToyLonLatSpace(17, 11; lon = (-102.0, -17.0), lat = (-42.0, 13.0),
            chunks = (17, 4))
        @test nchunks(dst) == 3

        plainspace = RasterGrid(DD.DimArray(data, (xd, yd)))
        eager = regrid(data; to = dst, from = plainspace,
            method = BarycentricPoint(), lazy = false)
        # The operator the chunk-pair builder gives over the whole domain, which
        # no tile build takes part in.
        whole = Matrix(GR.wholeblock(BarycentricPoint(), dst, plainspace).weights)
        @test all(isapprox(1.0), sum(whole; dims = 2))

        manifests = Dict{Tuple{Int,Int},Vector{Vector{Int}}}()
        for chunks in ((18, 18), (12, 9))
            counting = T7Counting(data, chunks)
            raster = DD.DimArray(counting, (xd, yd))
            space = RasterGrid(raster)
            plan = t7_plan(BarycentricPoint(), dst, space)
            A = LazyRegridArray(raster, plan)
            tiles = Vector{Int}[]
            for t in 1:nchunks(dst)
                dinds = ownedindices(dst, t)
                t7_reset!(counting)
                @test A[first(dinds):last(dinds)] ≈ eager[first(dinds):last(dinds)]

                manifest = plan.storage.tiles[t].sourcechunks
                row = Int.(GR.sourcesof(GR.dependencies(plan), t))
                @test manifest == t9_owners(whole, space, dinds)
                @test issubset(manifest, row)
                @test t7_spatial(counting) ==
                      sort([GR.chunkranges(space, s, (72, 36)) for s in manifest])
                push!(tiles, [Int(x) for x in manifest])
            end
            # Somewhere in this destination the row really is wider than the
            # manifest, so the read is not the row by coincidence.
            @test any(t -> length(tiles[t]) <
                           length(GR.sourcesof(GR.dependencies(plan), t)),
                1:nchunks(dst))
            manifests[chunks] = tiles
        end

        # One answer under both chunkings, and a manifest that is not the same
        # list: the stencils do not move with the chunking and the manifest
        # does.
        @test manifests[(18, 18)] != manifests[(12, 9)]
    end
end
