# The CopDEM production driver's scheduling policy: the walk order, the pull
# cursor with its taper, the refcount tile cache, and the prefetcher.
#
# `scripts/copdem_policy.jl` is deliberately free of Copernicus DEM, IGeo7, Zarr
# and GlobalRegridding — it takes the tile-to-chunk adjacency as two functions —
# so everything below runs on hand-written adjacencies in microseconds, with no
# regrid and no store. The one test that does touch GlobalRegridding checks the
# other half: that a real `chunk_dependency_graph` plugs into the same seam.

module CopDEMPolicyTests

using Test
import GlobalRegridding as GR
import GeometryOps as GO
const US = GO.UnitSpherical

include(joinpath(@__DIR__, "..", "..", "scripts", "copdem_policy.jl"))

# A toy adjacency: `rows[d]` is the source chunks destination chunk `d` reads.
function adjacency(rows::Vector{Vector{Int}}, nsrc::Int)
    deg = zeros(Int, nsrc)
    for r in rows, s in r
        deg[s] += 1
    end
    return (d -> rows[d]), (s -> deg[s]), length(rows), nsrc
end

"A cache over `rows` whose loader records what it was asked for."
function toycache(rows::Vector{Vector{Int}}, nsrc::Int; permits = 8, load = nothing)
    sourcesof, consumerdegree, ndst, _ = adjacency(rows, nsrc)
    calls = Threads.Atomic{Int}(0)
    loader = load === nothing ?
             (s -> (Threads.atomic_add!(calls, 1); fill(Float32(s), 4))) : load
    cache = RefCountCache{Vector{Float32}}(nsrc, ndst, sourcesof, consumerdegree,
        loader; permits)
    return cache, calls
end

"""
Drain `s` into `into`, appending every index it hands out.

A named function, not a `while` loop inlined into the `Threads.@spawn` below:
the enclosing testset scope already binds `b`, so an inlined
`while (b = claim!(s)) !== nothing` would make every worker assign the *same*
captured `Core.Box` and read back whichever batch another worker wrote last.
Here `b` is a local of one call, so each task owns its own.
"""
function drain!(s::GuidedSchedule, into::Vector{Int})
    while (b = claim!(s)) !== nothing
        append!(into, b)
    end
    return into
end

"How `got` differs from a clean cover of `1:n`, rather than merely that it does."
function coverage(got::Vector{Int}, n::Int)
    counts = zeros(Int, n)
    stray = Int[]
    for v in got
        1 <= v <= n ? (counts[v] += 1) : push!(stray, v)
    end
    return (; count = length(got), duplicated = findall(>(1), counts),
        uncovered = findall(iszero, counts), stray)
end

@testset "morton2" begin
    @test morton2(0, 0) == 0
    @test morton2(1, 0) == 1
    @test morton2(0, 1) == 2
    @test morton2(1, 1) == 3
    @test morton2(3, 0) == 0b000101
    # The interleave is a bijection on the 10-bit square.
    codes = [morton2(x, y) for x in 0:31, y in 0:31]
    @test allunique(codes)
    @test_throws ArgumentError morton2(1024, 0)
    @test_throws ArgumentError morton2(0, -1)
end

@testset "affinity_order" begin
    # Four chunks over four sources swept in the order 4, 3, 2, 1. A chunk is
    # emitted when its LAST source is swept, so the chunk holding source 4 goes
    # first and the one holding source 1 goes last.
    rows = [[1], [2], [3], [4]]
    order = affinity_order(4, d -> rows[d], [4, 3, 2, 1])
    @test order == [4, 3, 2, 1]

    # Chunks closing on the same source run earliest-opening first, so the
    # source that has been held longest is the first one released.
    rows = [[1, 3], [2, 3], [3]]
    order = affinity_order(3, d -> rows[d], [1, 2, 3])
    @test order == [1, 2, 3]

    # It is always a permutation, and a sourceless chunk sorts last.
    rows = [[2], Int[], [1, 2], [1]]
    order = affinity_order(4, d -> rows[d], [10, 20])
    @test sort(order) == 1:4
    @test order[end] == 2
end

@testset "GuidedSchedule taper" begin
    # Full batches while there is plenty, singles at the end.
    s = GuidedSchedule(100, 4, 8)
    @test claim!(s) == 1:8
    s = GuidedSchedule(10, 4, 8)
    @test claim!(s) == 1:3          # cld(10, 4)
    s = GuidedSchedule(4, 4, 8)
    @test claim!(s) == 1:1          # remaining <= W means one at a time
    s = GuidedSchedule(0, 4, 8)
    @test claim!(s) === nothing

    # `taper = false` is spelled `workers = 1`, which restores flat batching.
    s = GuidedSchedule(100, 1, 8)
    @test claim!(s) == 1:8

    # The tail is what the taper is for: with W workers and n items, the last
    # claim is a single item, so no worker is ever handed a batch while the
    # others idle.
    s = GuidedSchedule(37, 6, 8)
    last = nothing
    while (b = claim!(s)) !== nothing
        last = b
    end
    @test length(last) == 1

    # Concurrent claims cover 1:n exactly once, with a varying batch size. The
    # claim loop lives in `drain!` so that each task owns its own `b`; see there.
    n = 5000
    s = GuidedSchedule(n, 8, 8)
    seen = [Int[] for _ in 1:8]
    Threads.@sync for w in 1:8
        Threads.@spawn drain!(s, seen[w])
    end
    # Named parts, so a miscover says whether it dropped, doubled or invented.
    cover = coverage(vcat(seen...), n)
    @test cover.count == n
    @test isempty(cover.duplicated)
    @test isempty(cover.uncovered)
    @test isempty(cover.stray)
    @test claim!(s) === nothing
end

@testset "RefCountCache: eviction on the last consumer" begin
    # Source 1 is read by chunks 1 and 2; source 2 only by chunk 2.
    cache, calls = toycache([[1], [1, 2]], 2)
    @test cachestats(cache).loads == 0
    v = gettile!(cache, 1)
    @test v == fill(1.0f0, 4)
    @test calls[] == 1
    @test gettile!(cache, 1) === v          # cached, not rebuilt
    @test calls[] == 1
    @test cachestats(cache).hits == 1

    retire_column!(cache, 1)                 # one consumer left
    @test gettile!(cache, 1) === v
    retire_column!(cache, 2)                 # the last one
    @test cachestats(cache).live == 0
    @test completed(cache)
    @test isempty(doubleloaded(cache))
    # Source 2 was credited but never demanded: legitimate, and it was never
    # loaded either.
    @test cachestats(cache).loads == 1
end

@testset "RefCountCache: peak residency and at-most-once" begin
    rows = [[1, 2], [2, 3], [3, 4], [4]]
    cache, calls = toycache(rows, 4)
    for d in 1:4
        for s in rows[d]
            gettile!(cache, s)
        end
        retire_column!(cache, d)
    end
    @test calls[] == 4                       # each source built exactly once
    @test cachestats(cache).loads == 4
    @test isempty(doubleloaded(cache))
    @test completed(cache)
    # Only ever two sources alive at once, never all four.
    @test cachestats(cache).peaktiles == 2
    @test cachestats(cache).peakbytes == 2 * 16
end

@testset "RefCountCache: single-flight loading" begin
    started = Threads.Atomic{Int}(0)
    gate = Base.Event()
    cache, _ = toycache([[1] for _ in 1:20], 1;
        load = function (s)
            Threads.atomic_add!(started, 1)
            wait(gate)
            return fill(7.0f0, 4)
        end)
    tasks = [Threads.@spawn gettile!(cache, 1) for _ in 1:20]
    while started[] == 0
        yield()
    end
    notify(gate)
    vals = fetch.(tasks)
    @test all(v -> v == fill(7.0f0, 4), vals)
    @test all(v -> v === vals[1], vals)      # one array, not twenty
    @test started[] == 1
    # Each of the other nineteen either joined the flight or found it finished;
    # which of the two is a race, that none of them loaded is not.
    @test cachestats(cache).waits + cachestats(cache).hits == 19
    @test cachestats(cache).loads == 1
end

@testset "RefCountCache: resume, failure and uncredited demands" begin
    # Resume: retiring every consumer up front makes the source dead, and a
    # demand for a dead source is served transiently without entering the cache.
    cache, calls = toycache([[1], [1]], 1)
    retire_column!(cache, 1)
    retire_column!(cache, 2)
    @test cachestats(cache).loads == 0
    @test gettile!(cache, 1) == fill(1.0f0, 4)   # served, loudly counted
    @test cachestats(cache).uncredited == 1
    @test cachestats(cache).loads == 0           # never cached, never counted
    @test cachestats(cache).live == 0

    # Retiring twice is a bug, not a no-op.
    @test_throws ErrorException retire_column!(cache, 1)

    # A chunk whose regrid threw still retires, and its sources still free.
    cache, _ = toycache([[1, 2]], 2)
    gettile!(cache, 1)
    try
        error("regrid failed")
    catch
    finally
        retire_column!(cache, 1)
    end
    @test completed(cache)
    @test cachestats(cache).live == 0
end

@testset "RefCountCache: a load that throws" begin
    cache, _ = toycache([[1], [1]], 1; load = s -> error("no such tile $s"))
    err = try
        gettile!(cache, 1)
        nothing
    catch e
        e
    end
    @test err isa TileLoadError
    @test err.chunk == 1
    # Remembered, not retried: the second demand throws the same thing.
    @test_throws TileLoadError gettile!(cache, 1)
    @test isempty(doubleloaded(cache))
    retire_column!(cache, 1)
    retire_column!(cache, 2)
    @test completed(cache)
end

@testset "Prefetcher" begin
    rows = [[d] for d in 1:64]
    cache, calls = toycache(rows, 64)
    order = collect(1:64)
    sched = GuidedSchedule(64, 4, 1)
    pf = Prefetcher(cache, order, d -> rows[d], 64, sched; depth = 8, concurrency = 4)
    # The prefetcher primes before the cursor moves at all.
    while pf.issued[] < 8
        yield()
    end
    @test pf.issued[] >= 8
    while (b = claim!(sched)) !== nothing
        advance!(pf)
        for pos in b
            gettile!(cache, order[pos])
            retire_column!(cache, order[pos])
        end
    end
    stop_prefetch!(pf)
    @test prefetchstats(pf).failure === nothing
    @test quiescent(cache)
    @test all(istaskdone, pf.tasks)
    @test isempty(doubleloaded(cache))       # never both fetched and demanded
    @test completed(cache)
end

@testset "Prefetcher: a speculative failure is deferred, not raised" begin
    rows = [[1], [1]]
    cache, _ = toycache(rows, 1; load = s -> error("fetch failed"))
    sched = GuidedSchedule(2, 1, 1)
    pf = Prefetcher(cache, [1, 2], d -> rows[d], 1, sched; depth = 2, concurrency = 1)
    while pf.issued[] < 1
        yield()
    end
    stop_prefetch!(pf)
    # The prefetch task swallowed it...
    @test prefetchstats(pf).failure === nothing
    @test all(istaskdone, pf.tasks)
    # ...and the demand that follows raises it.
    @test_throws TileLoadError gettile!(cache, 1)
end

@testset "Prefetcher prepares a source before speculative decode" begin
    rows = [[1]]
    prepared = Threads.Atomic{Bool}(false)
    cache, _ = toycache(rows, 1;
        load = s -> (prepared[] || error("decoded before preparation"); fill(1.0f0, 4)))
    sched = GuidedSchedule(1, 1, 1)
    pf = Prefetcher(cache, [1], d -> rows[d], 1, sched;
        depth = 1, concurrency = 1, prepare = s -> (prepared[] = true))
    while !prepared[]
        yield()
    end
    stop_prefetch!(pf)
    @test prefetchstats(pf).failure === nothing
    @test gettile!(cache, 1) == fill(1.0f0, 4)
    retire_column!(cache, 1)
    @test completed(cache)
end

@testset "Prefetcher: a result nobody wants any more is discarded" begin
    gate = Base.Event()
    rows = [[1]]
    cache, _ = toycache(rows, 1;
        load = s -> (wait(gate); fill(1.0f0, 4)))
    sched = GuidedSchedule(1, 1, 1)
    pf = Prefetcher(cache, [1], d -> rows[d], 1, sched; depth = 1, concurrency = 1)
    while pf.issued[] < 1
        yield()
    end
    # The only consumer retires while the speculative fetch is still in the air.
    retire_column!(cache, 1)
    notify(gate)
    stop_prefetch!(pf)
    @test cachestats(cache).live == 0        # published into nothing
    @test quiescent(cache)
    @test completed(cache)
end

@testset "Prefetcher: stopping early leaves nothing running" begin
    rows = [[d] for d in 1:200]
    cache, _ = toycache(rows, 200)
    sched = GuidedSchedule(200, 1, 1)
    pf = Prefetcher(cache, collect(1:200), d -> rows[d], 200, sched;
        depth = 200, concurrency = 4)
    stop_prefetch!(pf)
    @test all(istaskdone, pf.tasks)
    @test quiescent(cache)
    @test prefetchstats(pf).failure === nothing
end

@testset "StripedLRUCache is the escape hatch, not the policy" begin
    calls = Threads.Atomic{Int}(0)
    cache = StripedLRUCache{Vector{Float32}}(
        s -> (Threads.atomic_add!(calls, 1); fill(Float32(s), 4));
        slots = 4, stripes = 1)
    for s in 1:4
        gettile!(cache, s)
    end
    @test calls[] == 4
    @test gettile!(cache, 1) == fill(1.0f0, 4)
    @test calls[] == 4                       # still resident
    gettile!(cache, 5)                       # evicts the least recent
    gettile!(cache, 2)                       # ...which was 2
    @test calls[] == 6
    # It answers the same protocol, and retirement is a no-op it can ignore.
    @test retire_column!(cache, 1) === nothing
    @test quiescent(cache)
    @test gettile!(cache, 1; speculative = true) === nothing
    @test cachestats(cache).policy === :lru
end

@testset "the seam holds a real chunk_dependency_graph" begin
    # Two tiny spaces, an actual graph, and the cache/order built from it exactly
    # the way the driver builds them.
    caps(pts, r) = [GO.UnitSpherical.SphericalCap(
        US.UnitSphereFromGeographic()(p), r) for p in pts]
    srcpts = [(0.0, 0.0), (2.0, 0.0), (4.0, 0.0), (6.0, 0.0)]
    dstpts = [(1.0, 0.0), (3.0, 0.0), (5.0, 0.0)]
    srccaps = caps(srcpts, deg2rad(1.2))
    # There are no `RegridSpace`s here to stamp, only caps, so the graph carries
    # the empty identity: it matches no space, and `validate_dependencies` will
    # refuse to certify it for reuse. That is the honest record for a relation
    # built by hand, and nothing in this testset reuses it.
    graph = GR._chunkgraph(GR.DependencyIdentity(), caps(dstpts, deg2rad(1.2)),
        GR._packedchunkindex(srccaps), length(srccaps), 0.0, nothing)
    @test GR.nsourcechunks(graph) == 4
    @test GR.ndestinationchunks(graph) == 3

    keys = [morton2(round(Int, p[1]) + 180, round(Int, p[2]) + 90) for p in srcpts]
    order = affinity_order(3, d -> GR.sourcesof(graph, d), keys)
    @test sort(order) == 1:3

    cache, calls = toycache([[Int(s) for s in GR.sourcesof(graph, d)] for d in 1:3],
        4)
    for d in order
        for s in GR.sourcesof(graph, d)
            gettile!(cache, Int(s))
        end
        retire_column!(cache, d)
    end
    @test completed(cache)
    @test isempty(doubleloaded(cache))
    @test cachestats(cache).uncredited == 0
end

end # module
