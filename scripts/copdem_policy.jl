# Scheduling policy for the CopDEM -> IGeo7 production run: the walk order, the
# tile cache, the pull cursor and the prefetcher. Included by
# `copdem_production.jl`; unit-tested by `test/scripts/copdem_policy.jl`.
#
# Nothing here knows about Copernicus DEM, IGeo7, Zarr, or GlobalRegridding. The
# only thing it is given about the workload is the tile <-> column adjacency,
# handed over as two plain functions:
#
#     sourcesof(d)       the source chunks destination chunk `d` may read
#     consumerdegree(s)  how many destination chunks may read source chunk `s`
#
# `GR.chunk_dependency_graph` supplies both; so does a hand-written adjacency in
# a test. That is the seam: the *graph query* lives in GlobalRegridding, the
# *policy* lives here, and neither depends on the other's types.
#
# The numbers this file is built to hit are measured in
# `regrid-notes/2026-08-21-chunk-dag-sim.md`. In one line: refcount eviction
# walked in tile-affinity order holds 1.43 GiB of tiles where a 3072-slot LRU
# holds 16.48 GiB, reloads nothing, and costs +0.013 % makespan.

# ===========================================================================
# The walk order
# ===========================================================================

"""
    morton2(x, y) -> Int

The 2-D Morton (Z-order) code of the non-negative lattice point `(x, y)`, ten
bits per axis.

Ten bits is enough for the 1-degree lon/lat lattice this run sweeps (360 x 180
after shifting into the non-negative quadrant) and is the width the offline
simulator used, so orders computed here reproduce its measurements exactly.
"""
function morton2(x::Integer, y::Integer)
    (0 <= x < 1024 && 0 <= y < 1024) ||
        throw(ArgumentError("morton2 takes 10-bit coordinates, got ($x, $y)"))
    z = 0
    for i in 0:9
        z |= ((Int(x) >> i) & 1) << (2i) | ((Int(y) >> i) & 1) << (2i + 1)
    end
    return z
end

"""
    affinity_order(ndst, sourcesof, srckeys) -> Vector{Int}

The destination chunks in tile-affinity order: `order[p]` is the destination
chunk to run at position `p`.

Sort the source chunks by `srckeys` — for this run a Morton sweep of the 1-degree
tile lattice — and give each one its rank in that sweep. A destination chunk is
then emitted **the moment its last source is swept**, ties broken by its *first*
source — so among chunks closing on the same tile, the one that has been holding
its earliest tile the longest runs first, and releases it soonest.

The rule optimises the *closing* of a source chunk, not the opening of one, and
that is exactly what refcount eviction rewards: by the time a destination chunk
runs, every source it needs was touched recently and the last of them is free to
go the moment it finishes. Measured against the alternatives on the real
workload: 1.43 GiB peak residency, against 11.91 GiB for a greedy
maximum-overlap walk and 16.48 GiB for the LRU it replaces.

Do **not** substitute an order by cost or degree. Descending degree with
refcount eviction measured 85 GiB — it front-loads the polar chunks, which
between them touch most of the globe, and nothing can be freed until the cheap
chunks sharing their tiles run at the very end.

A destination chunk with no sources at all sorts last; it has no tiles to be
near.
"""
function affinity_order(ndst::Integer, sourcesof, srckeys::AbstractVector)
    rank = invperm(sortperm(srckeys))
    keys = Vector{NTuple{3,Int}}(undef, Int(ndst))
    for d in 1:Int(ndst)
        ss = sourcesof(d)
        if isempty(ss)
            keys[d] = (typemax(Int), typemax(Int), d)
            continue
        end
        lo = hi = rank[Int(first(ss))]
        for s in ss
            r = rank[Int(s)]
            lo = min(lo, r)
            hi = max(hi, r)
        end
        keys[d] = (hi, lo, d)
    end
    return sortperm(keys)
end

# ===========================================================================
# The pull cursor, with a taper
# ===========================================================================

"""
    GuidedSchedule(n, workers, maxbatch)

One atomic cursor over `n` work positions, handed out in batches that shrink as
the queue drains — guided self-scheduling.

Batching keeps a worker on a geographically connected run, which is what makes
its tile cache worth having. It also owns the run's tail: the production
`batch = 8` hands the *last* pull eight chunks to one worker while the other 22
sit idle, and that measured **301 s** of dead time at the end of a 16.8 h run.
Tapering the batch to `ceil(remaining / W)` gives full batches while there is
plenty of work, singles once fewer than `W` chunks remain, and restores the
batch-1 tail of ~38 s at no cost to makespan.

The cursor stays the sole balancing mechanism. An order is only ever the
*sequence* the cursor walks; it is never a static partition, or the ~5x polar
cost skew stops being absorbed.
"""
struct GuidedSchedule
    cursor::Threads.Atomic{Int}
    n::Int
    workers::Int
    maxbatch::Int
end

GuidedSchedule(n::Integer, workers::Integer, maxbatch::Integer) =
    GuidedSchedule(Threads.Atomic{Int}(1), Int(n), max(1, Int(workers)),
        max(1, Int(maxbatch)))

"""
    claim!(s::GuidedSchedule) -> UnitRange{Int} or nothing

Claim the next batch of positions, or `nothing` when the queue is empty.

A compare-and-swap loop rather than `atomic_add!`, because the batch size now
depends on the cursor value: only the task that wins the CAS moves the cursor,
and it moves it by exactly the batch it is about to keep. Every successful CAS
therefore owns a half-open interval exclusively, so a varying batch size can
neither drop a position nor hand one out twice.
"""
function claim!(s::GuidedSchedule)
    while true
        pos = s.cursor[]
        pos > s.n && return nothing
        batch = min(s.maxbatch, cld(s.n - pos + 1, s.workers))
        stop = min(s.n, pos + batch - 1)
        Threads.atomic_cas!(s.cursor, pos, stop + 1) == pos && return pos:stop
    end
end

"The next position the cursor will hand out; `n + 1` once the queue is empty."
cursorposition(s::GuidedSchedule) = s.cursor[]

# ===========================================================================
# Tile caches
# ===========================================================================

"""
    TileLoadError(chunk, err)

A source chunk that could not be loaded. Stored in the cache rather than thrown
from the loading task, so that a *speculative* load never fails a run and a
later *demand* for the same chunk still fails loudly.
"""
struct TileLoadError <: Exception
    chunk::Int
    err::Any
end

Base.showerror(io::IO, e::TileLoadError) =
    print(io, "TileLoadError: source chunk ", e.chunk, ": ",
        sprint(showerror, e.err))

"A load in progress. `done` is level-triggered: latecomers never miss it."
struct Flight
    done::Base.Event
end

Flight() = Flight(Base.Event())

_nbytes(v::AbstractArray) = sizeof(v)
_nbytes(v) = Base.summarysize(v)

abstract type TileCache end

# ---------------------------------------------------------------------------
# Refcount eviction
# ---------------------------------------------------------------------------

"""
    RefCountCache{V}(nsrc, ndst, sourcesof, consumerdegree, load; permits)

A source-chunk cache that frees a chunk when the last destination chunk that
could read it has finished, and never loads a chunk twice.

The invariant, maintained under one lock, is

    remaining[s] == the number of destination chunks that may read `s`
                    and have not yet been retired

from which everything else follows. A destination chunk that is running has not
been retired, so nothing it may read can be freed underneath it. Once
`remaining[s]` reaches zero no correct future demand for `s` exists, so the
value is dropped — and because it can only reach zero once, `s` is never
reloaded. That is why the run's total load count equals the number of source
chunks it touches, structurally rather than by luck.

**Not a bounded cache.** Residency is bounded by the *order*, not by this type:
see [`affinity_order`](@ref) for the one that was measured, and the 85 GiB
counterexample it warns about.

One lock, not stripes. It is held only across array reads and writes — never
across a load, a wait, or a permit acquisition — so the 64 stripes the LRU
needed buy nothing here.

# Loading

`load(s)` builds source chunk `s`. It is called at most once per chunk, under a
`Base.Semaphore(permits)` shared with the prefetcher so that `permits` bounds
concurrent fetches whoever asked for them. Racing demands do not duplicate the
work: the first installs a [`Flight`](@ref) and the rest wait on it, which
matters once a load is an AWS request rather than 15 ms of arithmetic.

A load that throws is captured as a [`TileLoadError`](@ref) and remembered.
Speculative loads swallow it; the next demand throws it. There is no automatic
retry inside one process — a resumed run gets a fresh attempt.

# Demands the graph did not predict

The adjacency is a conservative superset of the true geometry, so every real
read should be a credited edge. Should one arrive anyway — `remaining[s] == 0` —
the chunk is loaded *transiently*, outside the cache, and `uncredited` is
incremented. Correctness is preserved and the run continues; the counter is
reported and checked at the end. Refusing to serve the read would trade a
reportable accounting bug for twenty lost hours.
"""
mutable struct RefCountCache{V,S,F} <: TileCache
    const sourcesof::S
    const load::F
    const lock::ReentrantLock
    const remaining::Vector{Int}
    const retired::BitVector
    const values::Vector{Union{Nothing,V}}
    const flights::Vector{Union{Nothing,Flight}}
    const failures::Vector{Union{Nothing,TileLoadError}}
    const attempts::Vector{Int}
    const demanded::BitVector
    const permits::Base.Semaphore
    bytes::Int
    tiles::Int
    peakbytes::Int
    peaktiles::Int
    uncredited::Int
    hits::Int
    waits::Int
    const uncreditedof::Vector{Int}     # which chunks, up to a cap, for the log
end

function RefCountCache{V}(nsrc::Integer, ndst::Integer, sourcesof, consumerdegree,
        load; permits::Integer = 16) where {V}
    ns, nd = Int(nsrc), Int(ndst)
    remaining = [Int(consumerdegree(s)) for s in 1:ns]
    return RefCountCache{V,typeof(sourcesof),typeof(load)}(
        sourcesof, load, ReentrantLock(), remaining, falses(nd),
        Union{Nothing,V}[nothing for _ in 1:ns],
        Union{Nothing,Flight}[nothing for _ in 1:ns],
        Union{Nothing,TileLoadError}[nothing for _ in 1:ns],
        zeros(Int, ns), falses(ns), Base.Semaphore(max(1, Int(permits))),
        0, 0, 0, 0, 0, 0, 0, Int[])
end

nsourcechunks(c::RefCountCache) = length(c.remaining)
ndestinationchunks(c::RefCountCache) = length(c.retired)

# Called with the lock held.
function _note_peak!(c::RefCountCache)
    c.bytes > c.peakbytes && (c.peakbytes = c.bytes)
    c.tiles > c.peaktiles && (c.peaktiles = c.tiles)
    return nothing
end

"""
    retire_column!(cache, d)

Record that destination chunk `d` will never read anything again, and free every
source chunk whose last consumer it was.

Called exactly once per destination chunk, on **every** terminal outcome:

  - a chunk skipped on resume, before the run starts;
  - a chunk that regridded successfully, as soon as `regrid` returns and before
    the result is written (the write reads no source data, and holding tiles
    across Zarr I/O only inflates residency);
  - a chunk whose regrid threw, in the same `finally`.

Miss any one of those and the count leaks: the chunk is pinned for the rest of
the run and the residency the whole design is for is gone.
"""
function retire_column!(c::RefCountCache, d::Integer)
    di = Int(d)
    lock(c.lock) do
        c.retired[di] && error("destination chunk $di retired twice")
        c.retired[di] = true
        for s in c.sourcesof(di)
            si = Int(s)
            n = c.remaining[si]
            n > 0 || error("refcount underflow on source chunk $si")
            c.remaining[si] = n - 1
            n == 1 || continue
            v = c.values[si]
            if v !== nothing
                c.bytes -= _nbytes(v)
                c.tiles -= 1
                c.values[si] = nothing
            end
            c.failures[si] = nothing
        end
    end
    return nothing
end

"""
    gettile!(cache, s; speculative = false) -> value or nothing

Source chunk `s`, loading it if this is the first ask.

A `speculative` call is the prefetcher's: it starts a load if nobody else has
one going, returns `nothing` either way, and never throws. A demand call blocks
until the value exists and rethrows a stored [`TileLoadError`](@ref).
"""
function gettile!(c::RefCountCache{V}, s::Integer; speculative::Bool = false) where {V}
    si = Int(s)
    waited = false
    while true
        flight = nothing
        lead = false
        ready = nothing
        failure = nothing
        transient = false
        lock(c.lock)
        try
            if speculative
                # Nothing to speculate on if it is dead, done, failed, or already
                # being fetched by somebody.
                if c.remaining[si] > 0 && c.values[si] === nothing &&
                   c.failures[si] === nothing && c.flights[si] === nothing
                    flight = Flight()
                    c.flights[si] = flight
                    c.attempts[si] == 0 ||
                        error("source chunk $si loaded twice")
                    c.attempts[si] += 1
                    lead = true
                end
            elseif c.remaining[si] == 0
                c.uncredited += 1
                length(c.uncreditedof) < 64 && push!(c.uncreditedof, si)
                transient = true
            else
                c.demanded[si] = true
                ready = c.values[si]
                if ready === nothing
                    failure = c.failures[si]
                    if failure === nothing
                        flight = c.flights[si]
                        if flight === nothing
                            flight = Flight()
                            c.flights[si] = flight
                            c.attempts[si] == 0 ||
                                error("source chunk $si loaded twice")
                            c.attempts[si] += 1
                            lead = true
                        elseif !waited
                            c.waits += 1
                        end
                    end
                elseif !waited
                    # A demand is counted once, for how it was answered. A
                    # waiter that comes round again to collect what it waited
                    # for is the same demand, not a second one.
                    c.hits += 1
                end
            end
        finally
            unlock(c.lock)
        end
        # An uncredited demand is served without ever entering the cache, so it
        # cannot disturb the refcounts or the at-most-once load count.
        transient && return _guarded(c, si)
        ready === nothing || return ready
        failure === nothing || throw(failure)
        if lead && speculative
            # A speculative load reports nothing and raises nothing: the failure
            # is stored, and the next demand for the chunk is what raises it.
            try
                _publish!(c, si, flight)
            catch err
                err isa TileLoadError || rethrow()
            end
            return nothing
        elseif lead
            # The leader keeps what it built. It cannot have been evicted from
            # under it: a demand leader is a live consumer, so the count it is
            # holding cannot have reached zero while it was loading.
            value = _publish!(c, si, flight)
            value === nothing || return value
        elseif flight === nothing
            return nothing            # speculative, and somebody else has it
        else
            wait(flight.done)
            waited = true
        end
    end
end

function _guarded(c::RefCountCache, si::Int)
    Base.acquire(c.permits)
    try
        return c.load(si)
    finally
        Base.release(c.permits)
    end
end

# Run the load outside the lock, then take the lock to publish it. Returns the
# value so that a demand leader keeps what it built rather than reading it back.
function _publish!(c::RefCountCache{V}, si::Int, flight::Flight) where {V}
    value = nothing
    failure = nothing
    Base.acquire(c.permits)
    try
        value = c.load(si)::V
    catch err
        failure = TileLoadError(si, err)
    finally
        Base.release(c.permits)
    end
    lock(c.lock) do
        c.flights[si] = nothing
        # Zero here means every consumer retired while this load was in the air:
        # the result is obsolete and is dropped rather than resurrected.
        c.remaining[si] > 0 || return nothing
        if failure === nothing
            c.values[si] = value
            c.bytes += _nbytes(value)
            c.tiles += 1
            _note_peak!(c)
        else
            c.failures[si] = failure
        end
        return nothing
    end
    notify(flight.done)
    failure === nothing || throw(failure)
    return value
end

"""
    cachestats(cache) -> NamedTuple

What the cache did: `loads` (source chunks loaded, each at most once),
`peakbytes` / `peaktiles` (high-water residency), `hits`, `waits` (demands that
joined somebody else's load), `uncredited` (demands the adjacency did not
predict — expected to be zero), and `live` (chunks still held).
"""
function cachestats(c::RefCountCache)
    lock(c.lock) do
        return (policy = :refcount, loads = sum(c.attempts), peakbytes = c.peakbytes,
            peaktiles = c.peaktiles, hits = c.hits, waits = c.waits,
            uncredited = c.uncredited, live = c.tiles, bytes = c.bytes,
            demanded = count(c.demanded), uncreditedof = copy(c.uncreditedof),
            retired = count(c.retired), pinned = count(>(0), c.remaining))
    end
end

"""
    quiescent(cache) -> Bool

No load is in flight. True after a clean shutdown, whether the run finished or
aborted.
"""
quiescent(c::RefCountCache) = lock(c.lock) do
    all(isnothing, c.flights)
end

"""
    completed(cache) -> Bool

Every destination chunk retired and every source chunk freed. Only meaningful
after a run that walked the whole queue.
"""
completed(c::RefCountCache) = lock(c.lock) do
    all(c.retired) && all(iszero, c.remaining) && all(isnothing, c.values)
end

"Source chunks loaded more than once. Empty, by construction; checked anyway."
doubleloaded(c::RefCountCache) = lock(c.lock) do
    findall(>(1), c.attempts)
end

# ---------------------------------------------------------------------------
# The striped LRU it replaces
# ---------------------------------------------------------------------------

"""
    StripedLRUCache{V}(nsrc, load; slots, stripes)

The pre-DAG cache, kept as the escape hatch `cachepolicy = :lru`: `stripes`
independent LRUs of `slots / stripes` entries each, under one lock apiece,
keyed by source chunk number.

Its peak residency is a property of `slots` and not of the workload — 5.49 MiB
times the slot count for mid-latitude tiles — so the 3072 slots the production
run used cost 16.48 GiB to hold a working set that never exceeded 665 tiles.
Shrinking it to match is not the answer either: at 256 slots the same walk
reloads 12.6 % of its tiles. It is here to fall back to, not to use.

Unlike [`RefCountCache`](@ref) a build happens *outside* the stripe lock and two
racers may both build the same chunk, the loser's copy being dropped. That is
cheap when a build is 15 ms of arithmetic and wrong when it is a paid AWS
request, which is why `prefetch_depth` is refused in this mode.
"""
struct StripedLRUCache{V,F} <: TileCache
    load::F
    caches::Vector{Dict{Int,V}}
    orders::Vector{Vector{Int}}
    locks::Vector{ReentrantLock}
    per::Int
    loads::Threads.Atomic{Int}
    hits::Threads.Atomic{Int}
end

function StripedLRUCache{V}(load; slots::Integer = 3072,
        stripes::Integer = 64) where {V}
    ns = max(1, Int(stripes))
    return StripedLRUCache{V,typeof(load)}(load,
        [Dict{Int,V}() for _ in 1:ns], [Int[] for _ in 1:ns],
        [ReentrantLock() for _ in 1:ns], max(1, Int(slots) ÷ ns),
        Threads.Atomic{Int}(0), Threads.Atomic{Int}(0))
end

function gettile!(c::StripedLRUCache, s::Integer; speculative::Bool = false)
    speculative && return nothing
    si = Int(s)
    st = mod(si, length(c.locks)) + 1
    lk, cache, order = c.locks[st], c.caches[st], c.orders[st]
    hit = lock(lk) do
        v = get(cache, si, nothing)
        v === nothing && return nothing
        push!(order, splice!(order, findfirst(==(si), order)))
        return v
    end
    if hit !== nothing
        Threads.atomic_add!(c.hits, 1)
        return hit
    end
    Threads.atomic_add!(c.loads, 1)
    v = c.load(si)
    return lock(lk) do
        existing = get(cache, si, nothing)
        existing === nothing || return existing
        length(cache) >= c.per && delete!(cache, popfirst!(order))
        cache[si] = v
        push!(order, si)
        return v
    end
end

retire_column!(::StripedLRUCache, d::Integer) = nothing
quiescent(::StripedLRUCache) = true

function cachestats(c::StripedLRUCache)
    live = sum(length, c.caches)
    bytes = sum(cache -> sum(_nbytes, values(cache); init = 0), c.caches)
    return (policy = :lru, loads = c.loads[], peakbytes = -1, peaktiles = -1,
        hits = c.hits[], waits = 0, uncredited = 0, live = live, bytes = bytes,
        demanded = -1, uncreditedof = Int[], retired = -1, pinned = -1)
end

# ===========================================================================
# The prefetcher
# ===========================================================================

"""
    Prefetcher(cache, order, sourcesof, schedule; depth, concurrency, prepare)

A pool that keeps the next `depth` positions of the walk order loaded before a
worker asks for them.

It matters exactly when a source chunk costs real latency to obtain — which the
synthetic run does not and the AWS run will. Modelled at `L` seconds per fetch
with eight concurrent fetches, the depth needed to hide it entirely is
`N ~ 32 L`: 8 at 0.2 s, 23 at 1 s, 64 at 2 s. Running with no prefetcher at
`L = 2 s` costs **+21 % makespan**; running with one costs +0.1 GiB of
residency. Raising `concurrency` to 16 halves the depth and the unavoidable
cold-start fill.

# Shape

One driver task, `concurrency` fetch tasks, one request channel, and the
`permits` semaphore *shared with the cache* — so `concurrency` bounds fetches in
flight however they were asked for, speculative and demanded together.

The driver walks a `frontier` forward through the order, never re-scanning, and
puts each not-yet-requested source chunk into the channel. It deduplicates with
a bitmap, so at most one request per source chunk exists and a channel sized to
the source count can never block. When it catches up with the cursor it waits on
an event the workers signal after each successful claim; the event is
autoresetting and the cursor is monotone, so a coalesced or early wakeup is
harmless and a missed one is impossible.

A fetch task first calls optional `prepare(s)`, then calls
`gettile!(cache, s; speculative = true)`. The real-tile driver uses that seam to
put the GeoTIFF on disk before decode; synthetic mode leaves it unset. The cache
call is a no-op when the chunk is resident, in flight, already failed, or already
dead. It never throws: a failed speculative load is remembered by the cache and
re-thrown at the next demand. So a prefetch can waste work, but it cannot fail a
run.

# Against the refcount cache

The prefetcher holds no consumer credit. A flight is a publication pin, not a
reference: if every consumer of a chunk retires while its speculative load is in
the air, the result is discarded rather than cached. This is the reason a
prefetcher cannot induce evict-then-refetch here — a count of zero means no
correct future demand exists, so there is nothing to refetch.
"""
mutable struct Prefetcher
    const requests::Channel{Int}
    const advanced::Base.Event
    const stop::Threads.Atomic{Bool}
    const tasks::Vector{Task}
    const issued::Threads.Atomic{Int}
    const depth::Int
    failure::Any
end

"Signal the prefetcher that the cursor moved. Cheap; called after every claim."
advance!(::Nothing) = nothing
advance!(pf::Prefetcher) = notify(pf.advanced)

function Prefetcher(cache::TileCache, order::Vector{Int}, sourcesof, nsrc::Integer,
        schedule::GuidedSchedule; depth::Integer, concurrency::Integer,
        prepare = nothing)
    d = Int(depth)
    d > 0 || throw(ArgumentError("prefetch depth must be positive, got $depth"))
    requests = Channel{Int}(max(1, Int(nsrc)))
    pf = Prefetcher(requests, Base.Event(true), Threads.Atomic{Bool}(false),
        Task[], Threads.Atomic{Int}(0), d, nothing)
    push!(pf.tasks, Threads.@spawn _prefetch_driver(pf, order, sourcesof,
        Int(nsrc), schedule))
    for _ in 1:max(1, Int(concurrency))
        push!(pf.tasks, Threads.@spawn _prefetch_worker(pf, cache, prepare))
    end
    return pf
end

function _prefetch_driver(pf::Prefetcher, order::Vector{Int}, sourcesof, nsrc::Int,
        schedule::GuidedSchedule)
    queued = falses(nsrc)
    frontier = 1
    n = length(order)
    try
        while !pf.stop[]
            target = min(n, cursorposition(schedule) + pf.depth - 1)
            while frontier <= target
                for s in sourcesof(order[frontier])
                    si = Int(s)
                    queued[si] && continue
                    queued[si] = true
                    put!(pf.requests, si)
                    Threads.atomic_add!(pf.issued, 1)
                end
                frontier += 1
            end
            frontier > n && break
            wait(pf.advanced)
        end
    catch err
        # The driver is speculative too. Record it, stop prefetching, let the
        # run carry on demand-loading.
        pf.failure === nothing && (pf.failure = err)
        pf.stop[] = true
    finally
        close(pf.requests)
    end
    return nothing
end

function _prefetch_worker(pf::Prefetcher, cache::TileCache, prepare)
    for s in pf.requests
        pf.stop[] && continue          # drain the channel, do no more work
        try
            prepare === nothing || prepare(s)
            gettile!(cache, s; speculative = true)
        catch err
            pf.failure === nothing && (pf.failure = err)
            pf.stop[] = true
        end
    end
    return nothing
end

"""
    stop_prefetch!(pf) -> Prefetcher or nothing

Shut the prefetcher down and join its tasks. Idempotent, and safe to call from a
`finally` while the run is unwinding an exception: a prefetch failure is
recorded in `pf.failure`, never rethrown, so it cannot mask the real one.
"""
stop_prefetch!(::Nothing) = nothing
function stop_prefetch!(pf::Prefetcher)
    pf.stop[] = true
    notify(pf.advanced)
    for t in pf.tasks
        try
            wait(t)
        catch err
            pf.failure === nothing && (pf.failure = err)
        end
    end
    return pf
end

"What the prefetcher did: requests issued, and any error it swallowed."
prefetchstats(::Nothing) = (depth = 0, issued = 0, failure = nothing)
prefetchstats(pf::Prefetcher) =
    (depth = pf.depth, issued = pf.issued[], failure = pf.failure)
