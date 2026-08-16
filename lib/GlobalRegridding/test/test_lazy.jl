# Chunk discovery, `LazyRegridArray`, streaming and spill. Owned by tasks T7, T8.
#
# These are the package's laziness laws (L1, L2, L4 of the task plan), and every
# one of them is an assertion about *what was read* or *what was built*, not
# about a number: the source is wrapped in a disk array that records each
# `readblock!` call, and the method in a wrapper that counts weight builds.
# Correctness is checked against the eager whole-domain answer on the same field.

import DiskArrays
import DimensionalData as DD

# --- A source that records every block read ---------------------------------

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

# --- A method that counts the weight builds it is asked for -----------------

mutable struct T7CountingMethod{M<:AbstractRegriddingMethod} <: AbstractRegriddingMethod
    inner::M
    builds::Int
end

T7CountingMethod(inner) = T7CountingMethod(inner, 0)
GR.support_radius(m::T7CountingMethod, space::RegridSpace) =
    GR.support_radius(m.inner, space)

function build_weights!(coo::WeightCOO, m::T7CountingMethod,
    dst::RegridSpace, dst_inds, src::RegridSpace, src_inds)
    m.builds += 1
    return build_weights!(coo, m.inner, dst, dst_inds, src, src_inds)
end

# --- A method whose stencil reaches a fixed angular distance ----------------
#
# Destination `j` is the mean of every source cell whose centre lies within
# `radius` of its centroid — a stencil that genuinely straddles a source chunk
# boundary, and reports the radius that makes it discoverable.

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
    budget = 2^30, chunks = nothing) =
    ChunkedPlan(method, policy, dst, src, storage, budget, chunks)

# The connectivity a flat pairwise cap test gives — the independent reference
# the tree descent is checked against.
function t7_pairwise(dst, dstchunk, src; radius = 0.0)
    dcap = chunktree(dst).caps[dstchunk]
    return [s for (s, scap) in enumerate(chunktree(src).caps)
            if US.spherical_distance(dcap.point, scap.point) <=
               dcap.radius + scap.radius + radius]
end

@testset "Lazy path" begin

    # A source chunked 4×2 into four chunks, a destination chunked by latitude
    # band into two contiguous ones: the two chunkings deliberately disagree, so
    # a destination chunk needs several source chunks.
    srcspace = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
        chunks = (4, 2))
    dstspace = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
        chunks = (8, 2))
    field = collect(reshape(1.0:32.0, 8, 4))

    @testset "discovery" begin
        # The dual descent and a flat pairwise cap check agree, chunk for chunk.
        # Over-reporting is allowed by the contract but not by this test: the
        # descent must not prune what the caps say meets.
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

        # Nothing was read and nothing was built: a plan and a lazy array are
        # geometry and bookkeeping only.
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

        # Exactly the connected source chunks, each read exactly once, each read
        # exactly one chunk's worth. A missed pair, a repeated read, or a read of
        # the whole array all fail here.
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

    @testset "L4 — plan reuse" begin
        source = T7Counting(field, (4, 2))
        method = T7CountingMethod(ToyDiagonalMethod())
        plan = t7_plan(method, dstspace, srcspace)
        A = LazyRegridArray(source, plan)

        first_read = A[1:16]
        builds = method.builds
        @test builds == length(t7_pairwise(dstspace, 1, srcspace))
        @test GR.nblocks(plan.storage) == builds

        # Re-reading the same destination chunk builds no weights. It does read
        # the source again — the plan caches weights, never data — and that
        # re-read is the documented cost of a second read.
        t7_reset!(source)
        @test A[1:16] == first_read
        @test method.builds == builds
        @test length(source.reads) == length(t7_pairwise(dstspace, 1, srcspace))
    end

    @testset "L4 — non-spatial slices reuse blocks" begin
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

        # All three slices of the same destination chunk at once: still no
        # build, and one spatial plan carried across the slices.
        three = A[1:16, 1:3]
        @test method.builds == builds
        @test three == hcat(vec(field)[1:16], 10 .* vec(field)[1:16],
            100 .* vec(field)[1:16])

        # The rest of the destination is new geometry, so it builds its own
        # pairs — once, for all three slices — and the whole agrees with the
        # eager N-D answer.
        every = A[1:32, 1:3]
        @test method.builds == builds + length(t7_pairwise(dstspace, 2, srcspace))
        eager = regrid(cube; to = dstspace, from = srcspace,
            method = ToyDiagonalMethod(), lazy = false)
        @test every == eager
    end

    @testset "dilation" begin
        # A destination cell sitting well inside the western source chunk, and a
        # stencil reaching 50° — far enough to take in the first cell of the
        # eastern chunk. The two chunks' extents do not meet, so undilated
        # discovery drops the pair outright.
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

        # Source centres are 20° apart from −70°; the destination centre is at
        # −35°, so cells 1–5 lie within 50° and their mean is 3. Dropping the
        # eastern chunk would renormalize over cells 1–4 and give 2.5 —
        # plausible, silent, and wrong.
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

        # A different source chunking is a different set of chunk pairs and the
        # same answer — the invariance the dilated discovery exists to protect.
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
        # `Weighted` blanks a destination whose source went missing, and the
        # lazy path reaches the same verdict as the eager one — the validity
        # mask is derived per source chunk, not per domain.
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
        # A real `RasterGrid` source and real conservative weights, in both
        # array orientations. A chunk is read as a rectangle of the *array*, so
        # its ranges follow the array's dimension order and not the space's
        # X-then-Y one; get that wrong and a `(Y, X)` raster reads out of bounds
        # or transposed.
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

    @testset "storage policies" begin
        @test_throws ErrorException Spilled(mktempdir())
        @test_throws ArgumentError PerChunk(0)
        @test_throws ArgumentError PerChunk(; maxbytes = 0)
    end
end
