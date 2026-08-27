# Copernicus DEM GLO-90 -> IGEO7 level 12, into one global ancestor-subzone Zarr
# store, written by W worker tasks in one process.
#
# Edit `copdem_config` in `dagger_regrid/copdem_helpers.jl`, then run it:
#
#     julia --project=benchmark -t 26 --gcthreads=8 scripts/copdem_production.jl
#
# The GC thread count carries NO second field. `--gcthreads=N,1` turns on Julia's
# concurrent page sweeper, which madvises freed pages from a background thread
# while the workers run; the 2026-08-21 run died of a SIGSEGV that is exactly what
# a page released out from under a live object looks like. The shared `gcguard`
# refuses to start under it. See regrid-notes/2026-08-21-polar-segfault.md.
#
# The run is resumable: it skips every chunk the ledger or the store already has,
# so re-running the same command continues where it stopped.
#
# Words used throughout:
#
#   tile    a 1x1-degree Copernicus DEM tile. The source unit: one source chunk.
#   chunk   one level-`ancestor` IGeo7 cell with all its level-`level`
#           descendants (7^7 = 823 543 cells at 5 -> 12). The destination work
#           unit, one Zarr chunk, one file. The store API calls it a "column".
#   budget  bytes of intersection weights a lazy regrid may hold at once, per
#           worker.
#
# `copdem_config().source` selects real, lazily downloaded GLO-90 GeoTIFFs or the analytic
# field in `copdem_synthetic.jl`. The tile LIST is real in both modes: Copernicus
# ships ~26 450 of the 64 800 tiles in the 1x1-degree lattice, and the source is a
# `PartialGrid` over exactly those. A destination chunk over open ocean therefore
# pairs with no source at all and stays nodata without a network request. Store
# I/O and resume live in `copdem_store.jl`.
const STARTED = Ref(time())

include(joinpath(@__DIR__, "dagger_regrid", "copdem_helpers.jl"))
const CONFIG = copdem_config()

# Destination chunks sampled at start-up to check that a column's own relation
# demands no tile the global one is missing. See the `check` in `run`.
const GRAPHMISS_SAMPLE = 16

# Atomic: workers report their own chunk failures, so `+= 1` would lose some.
const FAILURES = Threads.Atomic{Int}(0)
# The last run's cache statistics, for a harness that calls `main` in-process
# and wants the numbers rather than the log line.
const LASTCACHE = Ref{Any}(nothing)
# Likewise the lazy provider: the cold-network harness reads its counters after
# `main` returns. Production behavior does not depend on this observation seam.
const LASTPROVIDER = Ref{Any}(nothing)
# A one-second libuv wake keeps Julia's signal-profile report listener schedulable
# while every default-pool thread is occupied by the outer worker wave.
const PROFILEPUMP = Ref{Union{Nothing,Task}}(nothing)

"Report a claim and count the failures."
function check(name, ok; detail = "")
    lock(LOGLOCK) do
        ok || Threads.atomic_add!(FAILURES, 1)
        println(stamp(), "  ", ok ? "PASS  " : "FAIL  ", rpad(name, 56), detail)
        flush(stdout)
    end
    return ok
end

"Enable Julia's SIGUSR1 profile peek and write each report to the requested path."
function configureprofile()
    request = get(ENV, "COPDEM_PROFILE_REQUEST", "")
    isempty(request) && return false
    Profile.init(n = 5_000_000, delay = 0.01)
    Profile.set_peek_duration(60.0)
    Profile.peek_report[] = () -> begin
        try
            target = strip(read(request, String))
            isempty(target) && error("profile request file $request is empty")
            open(target, "w") do io
                println(io, "COPDEM_PROFILE_BEGIN duration_s=60 delay_s=0.01")
                Profile.print(io; format = :flat, sortedby = :count,
                    combine = true, mincount = 20)
                println(io, "COPDEM_PROFILE_END")
            end
            say("profile: wrote $target")
        catch err
            println(stderr, "COPDEM_PROFILE_ERROR ", sprint(showerror, err))
            flush(stderr)
        end
    end
    if PROFILEPUMP[] === nothing || istaskdone(PROFILEPUMP[])
        PROFILEPUMP[] = errormonitor(@async while true
            sleep(1.0)
        end)
    end
    say("profile: SIGUSR1 60-second peek armed; request file $request")
    return true
end

"""
    rssgib() -> Float64

CURRENT resident set size in GiB.

Deliberately not `Sys.maxrss()`, which is `ru_maxrss`: a monotone high-water
mark. A heartbeat printing that reports the largest transient the run has ever
had and never comes back down, so it looks like a staircase of plateaus when it
is really a record of peaks. On the GLO-90 run one sub-minute excursion to
81.65 GiB made every heartbeat for the next half hour claim "83.7 GiB" while the
true resident size was 64 GiB. Peaks are worth reporting — see
[`peakrssgib`](@ref) — but they have to be labelled as peaks.
"""
function rssgib()
    try
        return parse(Int, split(read("/proc/self/statm", String))[2]) * 4096 / 2^30
    catch
        return Sys.maxrss() / 2^30   # non-Linux fallback: a peak, but better than nothing
    end
end

"Peak resident set size in GiB since the process started."
peakrssgib() = Sys.maxrss() / 2^30


mutable struct Progress
    const lock::ReentrantLock
    done::Int
    skipped::Int
    cells::Int
    nan::Int
    const total::Int
    const started::Float64
    const every::Float64     # seconds between heartbeats
    const cold::Union{Nothing,Threads.Atomic{Int}}
    last::Float64
end

Progress(total, every, cold = nothing) =
    Progress(ReentrantLock(), 0, 0, 0, 0, total, time(), Float64(every), cold, time())

"""
    graphmisscheck(sys7, layout, srcspace, plan, todochunks, config) -> Bool

Sample `GRAPHMISS_SAMPLE` columns and `check` that none of them demands a source
tile the global relation does not hold.

Since Task E1 a column's lazy read takes its source tiles from that column's
**own** one-row relation, built in `regrid_chunk`, while the refcount cache
retires tiles by consumer count from the **global** relation `dagplan` owns. A
tile a column demanded but the global relation never held would be a silent
graph miss: the cache would have retired it, and the column would either reload
it or — with a bounded cache and no permits left — wait for a producer that is
never coming.

The two relations hold the same rows by construction. Both are one
`GR.candidatechunks!` query of the same source chunk index against the same
destination cap, and the caps are the same because a rooted subtree grid's chunk
cap is the covering's cap for that chunk (asserted on this exact shape in
`test/systems/crosssystem/regrid.jl`). So this samples the claim rather than
proving it. It costs about half a millisecond per column sampled, once, at
start-up.

Uses the column destination `regrid_chunk` itself builds, through the same
`DGG.columncell` mapping, so it checks the path the run takes rather than a
restatement of it.
"""
function graphmisscheck(sys7, layout, srcspace, plan, todochunks, config)
    n = min(GRAPHMISS_SAMPLE, length(todochunks))
    misses, checked = 0, 0
    for d in (n > 0 ? round.(Int, range(1, length(todochunks); length = n)) : Int[])
        a = DGG.columncell(layout, todochunks[d])
        colspace = DGG.DGGSpace(DGG.subtree(sys7, a, config.level);
            chunklevel = config.ancestor)
        col = GR.chunk_dependency_graph(colspace, srcspace; radius = plan.radius)
        row = GR.sourcesof(plan.graph, d)
        misses += count(s -> !insorted(s, row), GR.sourcesof(col, 1))
        checked += 1
    end
    return check("no column demands a tile the graph does not hold",
        misses == 0; detail = "$checked column(s) sampled, $misses miss(es)")
end


"""
    etaseconds(elapsed, computed, remaining) -> Float64

How much longer the remaining chunks take at the rate THIS SESSION has computed
chunks — `NaN` until the first one lands.

`computed` deliberately excludes the chunks skipped at resume. Those cost a set
lookup, not a regrid, so crediting them to the session's elapsed time makes the
rate look arbitrarily fast: a resumed run of the GLO-90 store reported "ETA
0.30 h" with about three hours left, because 54 k skips and 4.5 k real chunks
were divided into one session's clock as if all 58 k had been computed in it.
"""
etaseconds(elapsed, computed, remaining) =
    computed > 0 ? elapsed * remaining / computed : NaN

"""
    heartbeat!(p, force = false)

One progress line, at most every `p.every` seconds.

Every rate in it is SESSION-scoped and labelled so, because that is all a
`Progress` knows: it counts what this process computed. Only the chunk count is
store-wide, and it names both numbers — `d/t chunks (c computed this session)` —
so neither reading can be mistaken for the other. Store-wide cells and NaNs come
from the ledger instead, at the end of the run; see [`storetotals`](@ref).
"""
function heartbeat!(p::Progress, force = false)
    lock(p.lock) do
        now = time()
        (force || now - p.last >= p.every) || return nothing
        p.last = now
        el = now - p.started
        seen = p.done + p.skipped
        eta = etaseconds(el, p.done, p.total - seen)
        cold = p.cold === nothing ? "" : " | demand-cold downloads $(p.cold[])"
        say(@sprintf("HEARTBEAT  %d/%d chunks (%d computed this session, %d skipped) | session %.3e cells, %.0f cells/s | elapsed %s | ETA %s | RSS %.1f GiB (peak %.1f) | session NaN %.2f%%%s",
            seen, p.total, p.done, p.skipped, Float64(p.cells),
            p.cells / max(el, 1e-9), hours(el), hours(eta),
            rssgib(), peakrssgib(), 100 * p.nan / max(p.cells, 1), cold))
        return nothing
    end
end

"""
    reportchecks(donelog, written)

Dry-run self-checks on the progress ARITHMETIC — the part of the run that has no
regrid to verify it, and that a resumed run twice got wrong: an ETA that counted
chunks skipped at resume as work this session, and a final banner whose totals
were session-scoped without saying so. Both were cosmetic; both were believed.

Cheap enough to run on every dry run, and it reads the run's real ledger, so it
also says whether these totals can be trusted for THIS store.
"""
function reportchecks(donelog, written)
    # 10 chunks computed in 10 s with 90 to go is 90 s. Crediting 100 skips to
    # the same 10 s would say 5 s, which is the shape of the bug.
    check("ETA counts only chunks computed this session",
        etaseconds(10.0, 10, 90) ≈ 90.0;
        detail = @sprintf("%.1f s, not the %.1f s that crediting 100 skips gives",
            etaseconds(10.0, 10, 90), etaseconds(10.0, 110, 90)))
    check("ETA is NaN until the first chunk lands", isnan(etaseconds(10.0, 0, 90)))

    # A chunk logged twice — recomputed because a crash lost its file, not its
    # line — is one chunk, at its LAST line's values.
    tmp = tempname()
    try
        open(tmp, "w") do io
            println(io, "{\"col\":7,\"cells\":100,\"nan\":10,\"secs\":1.00,\"w\":1,\"t\":\"x\"}")
            println(io, "{\"col\":7,\"cells\":100,\"nan\":40,\"secs\":1.00,\"w\":2,\"t\":\"x\"}")
            println(io, "{\"col\":9,\"cells\":100,\"nan\":0,\"secs\":1.00,\"w\":1,\"t\":\"x\"}")
        end
        t = storetotals(tmp, Set([7, 9, 11]))
        check("a chunk logged twice counts once, at its last line",
            length(doneledger(tmp)) == 2 && t.cells == 200 && t.nan == 40;
            detail = "$(length(doneledger(tmp))) entries, $(t.cells) cells, $(t.nan) NaN")
        check("a written chunk with no ledger line is unaccounted",
            t.chunks == 3 && t.unaccounted == 1;
            detail = "$(t.chunks) chunks, $(t.unaccounted) unaccounted")
    finally
        rm(tmp; force = true)
    end

    # And against the ledger this run will actually append to.
    led = doneledger(donelog)
    tot = storetotals(donelog, written)
    check("store totals span the whole resume union",
        tot.chunks == length(written) &&
        tot.unaccounted == count(c -> !haskey(led, c) || led[c].cells < 0, written);
        detail = "$(tot.chunks) chunks written, $(tot.unaccounted) with no ledger line")
    check("store cells are the ledger sum over that union",
        tot.cells == sum(Int[led[c].cells for c in written
                             if haskey(led, c) && led[c].cells >= 0]; init = 0);
        detail = @sprintf("%.4e cells, NaN %.3f%%", Float64(tot.cells),
            100 * tot.nan / max(tot.cells, 1)))
    return nothing
end

"""
    workercount(config) -> (W, shape)

The worker count `W`, sized from the `cores` budget rather than from a thread
count, and the parallel `shape` it assumes. `config.workers` overrides `W` when
set to something positive.

A work unit is one chunk, and how many cores one worker can hold depends
entirely on which of the two nested parallelisms is allowed to run:

  - `shape = :outer` sets [`GR.OUTER_PARALLEL`](@ref) in the worker body. The
    weight build stays serial and the lazy plan's source wave is the only
    parallelism inside a unit. Measured: **1.06 cores per worker**, flat in `W`
    while `W` is well below `nthreads`, and **65 k cells per core-second** — the
    most work per core of the two, but it needs `W ≈ cores` workers, and each
    worker holds a chunk's weights, so resident memory scales with `W`.

  - `shape = :inner` leaves `OUTER_PARALLEL` unset, so a narrow wave stands
    aside and ConservativeRegridding threads the weight build itself. Measured
    at 12 threads: one worker reaches 7.0–7.3 cores, the pool saturates at
    **`W ≈ nthreads/4`** holding ~78 % of it, and throughput is **60 k cells per
    core-second** — about 8 % less work per core, for a third of the workers and
    two thirds of the resident memory.

So at a fixed core budget `:outer` converts cores into cells slightly faster,
while `:inner` reaches the same cores with far fewer workers, and is the only
shape that keeps scaling once `cores` exceeds what `W ≈ cores` workers can be
given threads for. `:outer` is the default because the standing budget is 24
cores, where the throughput edge wins; switch to `:inner` for a larger budget or
when resident memory is the binding constraint.
"""
function workercount(config)
    shape = config.shape
    shape in (:outer, :inner) ||
        error("shape must be :outer or :inner, got $(repr(shape))")
    config.workers > 0 && return (config.workers, shape)
    config.cores > 0 || error("cores must be positive when workers is not set")
    # Cores per worker, measured; see above.
    W = shape === :outer ? cld(config.cores, 1.06) : cld(config.cores, 3.1)
    return (max(1, Int(W)), shape)
end

"""
    runchunks(chunks, order, dem, srcspace, sys7, store, layout, skipped, donelog,
              prefetcher, config) -> Progress

`W` tasks — [`workercount`](@ref) — pulling batches of indices off one
[`GuidedSchedule`](@ref), each computing its chunk and writing it straight to its
own Zarr file. Chunks are disjoint, so no two tasks ever touch the same file.

`order[p]` is the graph destination number to run at index `p`, and
`chunks[order[p]]` is the level-`ancestor` chunk itself. The order is a priority
sequence, never a partition: pulling dynamically is what keeps the polar chunks —
which meet 200-360 tiles and cost ~5x a mid-latitude one — off the critical path,
and tapering the batch is what keeps the last pull from handing one worker eight
of them while everybody else is idle.

Every chunk is retired from the cache exactly once, on every terminal outcome and
before its result is written: the Zarr write reads no source data, so holding
tiles across it only inflates residency. Chunks skipped on resume were retired
before the run started and are not in `order` at all.

Under `shape = :outer` each worker body sets `GR.OUTER_PARALLEL`, so the nested
weight builds do not spawn their own tasks and oversubscribe the box; the lazy
plan's source wave is then the only parallelism inside a unit. Under
`shape = :inner` the scoped value is left alone, and `_wavesize` weighs that wave
against the threading inside the build — standing the wave down whenever one
source chunk carries the destination chunk, which is almost always.
"""
function runchunks(chunks, order, sched, dem, srcspace, sys7, store, layout, skipped,
        donelog, prefetcher, config)
    W, shape = workercount(config)
    outer = shape === :outer
    provider = dem.builder.provider
    p = Progress(length(chunks), config.heartbeat,
        provider === nothing ? nothing : provider.ncold)
    p.skipped = skipped
    log = DoneLog(donelog)
    try
        # `false` is the scoped value's own default, so binding it is the same
        # as leaving it unset — the worker is the outermost scope either way.
        Threads.@sync for w in 1:W
            Threads.@spawn @with GR.OUTER_PARALLEL => outer begin
                while true
                    batch = claim!(sched)
                    batch === nothing && break
                    advance!(prefetcher)
                    for pos in batch
                        d = order[pos]
                        ch = chunks[d]
                        t0 = time()
                        vals = try
                            regrid_chunk(dem, srcspace, sys7, layout, ch, config)
                        catch err
                            say("ERROR worker $w chunk $ch: ",
                                sprint(showerror, err, catch_backtrace())[1:min(end, 1500)])
                            Threads.atomic_add!(FAILURES, 1)
                            nothing
                        finally
                            retire_column!(dem.cache, d)
                        end
                        vals === nothing && continue
                        DGG.dggwrite!(store, ch, vals)
                        el = time() - t0
                        nnan = count(isnan, vals)
                        lock(p.lock) do
                            p.done += 1
                            p.cells += length(vals)
                            p.nan += nnan
                        end
                        record!(log, ch, length(vals), nnan, el, w)
                        say(@sprintf("w%02d chunk %d  %d cells  %.1f s  %.0f cells/s  %.1f%% NaN",
                            w, ch, length(vals), el, length(vals) / el,
                            100 * nnan / length(vals)))
                        heartbeat!(p)
                    end
                end
            end
        end
    finally
        close(log)
    end
    heartbeat!(p, true)
    return p
end

"""
    reportcache(cache, prefetcher, provider, nsrc, complete)

What the tile cache did, and whether it did it lawfully.

Three claims, in decreasing strength:

  - **Every source chunk was loaded at most once.** This is the structural
    property of refcount eviction and the reason it never reloads: a chunk can
    only be freed when no destination chunk that may read it is left, so no
    correct demand for it can follow. If this fails the refcounts are wrong.

  - **Every demand was one the graph predicted.** The adjacency is a
    conservative superset of the true geometry, so a read the graph did not
    carry is a hole in the covering. Such a read is still served — the run does
    not die over an accounting bug — but it is counted here.

  - **Every source chunk was loaded at least once**, reported only for a fresh,
    complete, failure-free run. This is a COVERAGE diagnostic, not a law: the
    graph's 1.35x conservatism, a resumed run, and a failed chunk each leave
    tiles legitimately untouched. The offline simulator could equate loads with
    the tile count because it treated adjacency as demand; production reads an
    unknown subset of it.
"""
function reportcache(cache, prefetcher, provider, nsrc::Int, complete::Bool)
    s = cachestats(cache)
    LASTCACHE[] = s
    if s.policy === :refcount
        say(@sprintf("tile cache: %d loads of %d source chunks, peak %.2f GiB / %d tiles, %d hits, %d joined loads, %d live at end",
            s.loads, nsrc, s.peakbytes / 2^30, s.peaktiles, s.hits, s.waits, s.live))
        check("every source chunk was loaded at most once", isempty(doubleloaded(cache));
            detail = "$(length(doubleloaded(cache))) chunk(s) loaded twice")
        check("every demand was a predicted edge", s.uncredited == 0;
            detail = s.uncredited == 0 ? "" :
                     "$(s.uncredited) uncredited demand(s), chunks " *
                     string(s.uncreditedof))
        check("every destination chunk was retired", s.pinned == 0 && s.live == 0;
            detail = "$(s.pinned) source chunk(s) still pinned, $(s.live) still held")
        complete && say(s.loads == nsrc ?
            "coverage: every one of the $nsrc source chunks was read" :
            "coverage: $(nsrc - s.loads) of $nsrc source chunks were never read " *
            "(the graph is a conservative superset; this is not an error)")
    else
        say(@sprintf("tile cache: LRU, %d loads (%d reloads over %d source chunks), %d hits, %d live at end, %.2f GiB held",
            s.loads, max(0, s.loads - nsrc), nsrc, s.hits, s.live, s.bytes / 2^30))
    end
    provider === nothing || say(
        "source downloads: $(provider.ndownloads[]) total, " *
        "$(provider.ncold[]) demand-cold")
    pf = prefetchstats(prefetcher)
    pf.depth == 0 && return nothing
    say("prefetch: depth $(pf.depth), $(pf.issued) tile requests issued")
    check("the prefetcher raised nothing", pf.failure === nothing;
        detail = pf.failure === nothing ? "" : sprint(showerror, pf.failure))
    return nothing
end

# ===========================================================================
# Main
# ===========================================================================


function main(config = CONFIG)
    gcguard(config)
    regridmethod(config)   # fail fast on a bad `method`, before the land mask
    configureprofile()
    realspec = effective_realspec(config.source, config.real)
    println("="^92)
    println(stamp(), "  copdem_production.jl — GLO-$(config.res) -> IGEO7 level " *
                     "$(config.level), level-$(config.ancestor) chunks, " *
                     "$(uppercase(String(config.source))) elevations, " *
                     "$(uppercase(String(config.method))) method")
    println("  julia $(VERSION)  threads=$(Threads.nthreads())  " *
            "gcmark=$(Base.JLOptions().nmarkthreads)  " *
            "gcsweep=$(Base.JLOptions().nsweepthreads)  pid=$(getpid())")
    println("  ", config)
    println("="^92)
    flush(stdout)

    say(config.malloctrim > 0 ?
        "malloc: M_TRIM_THRESHOLD frozen at $(config.malloctrim >> 20) MiB — " *
        (tunemalloc(config.malloctrim) ? "ok" : "MALLOPT FAILED, expect ~3x resident memory") :
        "malloc: left at glibc defaults (malloctrim = 0); expect ~3x resident memory")

    sys = DGG.CopernicusDEMSystem(config.res)
    storesys7 = DGG.IGeo7System()
    sys7 = config.authalic ? DGG.AuthalicSystem(storesys7) : storesys7
    capacity = 7^(config.level - config.ancestor)

    tiledir = joinpath(config.data, "CopernicusDEM", "$(config.res)m")
    tilelist = copdem_tilelist(config)
    landshp = joinpath(config.data, "naturalearth", "ne_10m_land.shp")
    donelog = donelogpath(config.store)
    chunklist = chunklistpath(config.store)

    # Real tiles carry their own ocean nodata; only a synthetic source masks.
    mask = config.source === :synthetic ? landmask(landshp, config.maskarcsec) : NOMASK

    # --- the source ------------------------------------------------------
    t0 = time()
    tiles = listedtiles(sys, tilelist, config.region)
    isempty(tiles) && error("the tile list and region select no tiles")
    real = realtiles(sys, tiledir, realspec)
    filter!(p -> p.first in Set(tiles), real)
    provider = config.source === :real ? LazyCopernicusTiles(sys, tiles;
        cachedir = config.tilecache, baseurl = config.tilebaseurl,
        retries = config.retries, backoff = config.backoff,
        timeout = config.timeout) : nothing
    LASTPROVIDER[] = provider
    ids = TileIds(sys, tiles)
    say("tile list: $(length(tiles)) listed tiles of $(DGG.ncells(sys, 0)) " *
        "($(round(100 * length(tiles) / DGG.ncells(sys, 0); digits = 1))%), " *
        "$(length(real)) local GeoTIFF overrides, " *
        @sprintf("%.3e pixels, %s", Float64(ids.n), secs(time() - t0)))
    provider === nothing || say("lazy GLO-90 source: $(provider.baseurl), cache $(provider.cachedir)")
    for (o, p) in sort!(collect(real))
        say("  real: $(tilestem(sys, DGG.LevelIndex(0, o)))  " *
            "$(round(filesize(p) / 2^10)) KiB")
    end

    builder = TileBuilder(sys, tiles, real, provider, mask; delay = config.fetchdelay)
    config.fetchdelay > 0 && say(@sprintf(
        "fetchdelay: every tile build sleeps %.3f s — a latency shim, not production",
        config.fetchdelay))
    t0 = time()
    srcgrid = DGG.PartialGrid(sys, 1, ids)
    srcspace = DGG.DGGSpace(srcgrid; chunklevel = 0)
    say("source space: $(GR.nchunks(srcspace)) chunks over " *
        @sprintf("%.3e pixels, built in %s", Float64(DGG.ncells(srcgrid)), secs(time() - t0)))
    check("source chunks are the listed tiles", GR.nchunks(srcspace) == length(tiles);
        detail = "$(GR.nchunks(srcspace)) chunks, $(length(tiles)) tiles")

    # --- the destination chunks ------------------------------------------
    todochunks = load_chunklist(chunklist)
    if todochunks === nothing
        t0 = time()
        todochunks = covering_chunks(sys7, sys, tiles, config.ancestor;
            nthreads = max(1, Threads.nthreads() - 1))
        save_chunklist(chunklist, config.ancestor, todochunks)
        say("chunks: $(length(todochunks)) level-$(config.ancestor) chunks cover the " *
            "tiles, computed in $(secs(time() - t0)), cached at $chunklist")
    else
        say("chunks: $(length(todochunks)) read from $chunklist")
    end
    if !isempty(config.chunks)
        todochunks = sort!(unique(config.chunks))
        say("chunks: overridden to $(length(todochunks)) named by `chunks`")
    end
    if 0 < config.maxchunks < length(todochunks)
        todochunks = todochunks[round.(Int, range(1, length(todochunks);
            length = config.maxchunks))]
        say("chunks: limited to $(config.maxchunks) by `maxchunks`")
    end
    ncell = length(todochunks) * capacity
    say(@sprintf("work: %d chunks x %d cells = %.4e level-%d cells, %.1f GiB dense f32",
        length(todochunks), capacity, Float64(ncell), config.level, 4 * ncell / 2^30))

    # --- the store -------------------------------------------------------
    geometry_tag = config.authalic ? sprint(show, sys7) : nothing
    store = openstore(config, storesys7, capacity; geometry_tag)
    layout = store.layout
    check("store layout is level $(config.level) over level-$(config.ancestor) chunks",
        DGG.system(layout) == storesys7 && DGG.level(layout) == config.level &&
        layout.ancestor_level == config.ancestor &&
        layout.capacity == capacity;
        detail = "system $(DGG.system(layout)), level $(DGG.level(layout)), " *
                 "ancestor $(layout.ancestor_level), " *
                 "capacity $(layout.capacity)")

    done = config.resume ? donechunks(donelog, config.store, "elevation") : Set{Int}()
    todo = count(c -> !(c in done), todochunks)
    say("resume: $(length(done)) chunks already written, $todo of " *
        "$(length(todochunks)) to do")

    # --- the dependency graph, the walk order, the cache ------------------
    W, shape = workercount(config)
    plan = dagplan(sys, sys7, tiles, todochunks, srcspace, config)
    say(@sprintf("graph: %d x %d chunks, %d edges (%.2f tiles per chunk), radius %g; dst space %s, graph %s, order %s",
        GR.nsourcechunks(plan.graph), GR.ndestinationchunks(plan.graph), plan.edges,
        plan.edges / max(length(todochunks), 1), plan.radius,
        secs(plan.tspace), secs(plan.tgraph), secs(plan.torder)))
    # One relation, one owner: everything below reads the object the global
    # plan holds, and there is no second one anywhere in this run.
    check("the graph is the global plan's own relation",
        plan.graph === GR.dependencies(plan.globalplan))
    check("graph destination chunks are the work list",
        GR.ndestinationchunks(plan.graph) == length(todochunks))
    check("graph source chunks are the listed tiles",
        GR.nsourcechunks(plan.graph) == length(tiles))
    check("the walk order is a permutation of the work list",
        length(plan.order) == length(todochunks) && isperm(plan.order))
    graphmisscheck(sys7, layout, srcspace, plan, todochunks, config)

    cache = if config.cachepolicy === :refcount
        RefCountCache{Vector{Float32}}(length(tiles), length(todochunks),
            d -> GR.sourcesof(plan.graph, d),
            s -> GR.consumerdegree(plan.graph, s),
            k -> buildtile(builder, k);
            permits = config.fetchconc > 0 ? config.fetchconc : W)
    elseif config.cachepolicy === :lru
        config.prefetch == 0 || error("prefetch needs cachepolicy = :refcount: a " *
            "bounded LRU can evict a prefetched tile before it is used, and the " *
            "churn that causes was never measured")
        StripedLRUCache{Vector{Float32}}(k -> buildtile(builder, k);
            slots = config.cache, stripes = config.stripes)
    else
        error("cachepolicy must be :refcount or :lru, got $(repr(config.cachepolicy))")
    end
    dem = TiledDEM(ids, builder, cache)

    # A chunk skipped on resume is retired before anything starts, so its tiles
    # are never credited to work that will not happen, and the run's hot loop
    # never has to ask again whether a chunk is done.
    order = plan.order
    if !isempty(done)
        keep = Int[]
        for d in order                       # `order` holds destination numbers
            if todochunks[d] in done
                config.cachepolicy === :refcount && retire_column!(cache, d)
            else
                push!(keep, d)
            end
        end
        order = keep
    end
    nskipped = length(plan.order) - length(order)
    say("schedule: $(config.schedule) order over $(length(order)) chunks " *
        "($nskipped retired on resume), cache=$(config.cachepolicy), " *
        "taper=$(config.taper)")

    if config.dryrun
        say(@sprintf("dryrun: max tiles per chunk %d, max chunks per tile %d",
            maximum(d -> GR.sourcedegree(plan.graph, d), 1:length(todochunks); init = 0),
            maximum(s -> GR.consumerdegree(plan.graph, s), 1:length(tiles); init = 0)))
        reportchecks(donelog, done)
        say("dryrun: stopping before any regrid")
        return FAILURES[]
    end

    # --- the run ---------------------------------------------------------
    say("launching W=$W workers, shape=$shape, sized to hold $(config.cores) cores " *
        "(nthreads=$(Threads.nthreads())), batch=$(config.batch), budget=$(config.budget)")
    shape === :inner && Threads.nthreads() < 4W &&
        say("NOTE shape=:inner saturates near W=nthreads/4; " *
            "$(Threads.nthreads()) threads will not fill $W workers")
    sched = GuidedSchedule(length(order), config.taper ? W : 1, config.batch)
    prefetcher = nothing
    t0 = time()
    p = try
        if config.prefetch > 0
            fc = config.fetchconc > 0 ? config.fetchconc : W
            prepare = provider === nothing ? nothing :
                s -> tilepath!(provider, tiles[s]; demand = false)
            prefetcher = Prefetcher(cache, order, d -> GR.sourcesof(plan.graph, d),
                length(tiles), sched; depth = config.prefetch, concurrency = fc,
                prepare)
            say("prefetch: depth $(config.prefetch) chunks, $fc concurrent fetches")
        end
        runchunks(todochunks, order, sched, dem, srcspace, sys7, store, layout,
            nskipped, donelog, prefetcher, config)
    finally
        stop_prefetch!(prefetcher)
    end
    wall = time() - t0
    # Two banners, never one. `p` counts what THIS process computed; the store
    # holds what every session ever wrote. On a resumed run those differ by
    # whatever the earlier sessions did, and the session numbers on their own
    # read as the whole dataset: session 3 of the GLO-90 run signed off with
    # "9.7606e9 cells ... 44.308% NaN" when the store held 5.4500e10 cells at
    # 26.892% NaN. Same ledger the resume set came from, so the two agree.
    written = donechunks(donelog, config.store, "elevation"; label = "store")
    tot = storetotals(donelog, written)
    say(@sprintf("RUN DONE  session: %d chunks computed, %d skipped, %.4e cells in %s = %.0f cells/s aggregate",
        p.done, p.skipped, Float64(p.cells), hours(wall), p.cells / max(wall, 1e-9)))
    say(@sprintf("RUN DONE  session: NaN %.3f%% of the %.4e cells this session wrote",
        100 * p.nan / max(p.cells, 1), Float64(p.cells)))
    say(@sprintf("RUN DONE  store total: %d of %d chunks written, %.4e cells, NaN %.3f%%%s",
        tot.chunks, length(todochunks), Float64(tot.cells),
        100 * tot.nan / max(tot.cells, 1),
        tot.unaccounted == 0 ? "" :
            " — cells and NaN over the $(tot.chunks - tot.unaccounted) chunks the " *
            "ledger accounts for; $(tot.unaccounted) written chunks have no ledger line"))
    say(@sprintf("source tiles decoded: %d real, %d synthetic; %.3e of %.3e pixels valid (%.2f%% land)",
        builder.nreal[], builder.nsynthetic[], Float64(builder.nland[]),
        Float64(builder.npixels[]),
        100 * builder.nland[] / max(builder.npixels[], 1)))
    reportcache(cache, prefetcher, provider, length(tiles),
        FAILURES[] == 0 && nskipped == 0)
    say(@sprintf("RSS %.2f GiB now, %.2f GiB peak", rssgib(), peakrssgib()))

    if config.checks
        sm = SourceMask(mask, DGG.levelgrid(sys, 0), Set(tiles))
        verify(config, sys7, layout, todochunks, sm, real, dem, written)
    end

    say("total wall $(hours(time() - STARTED[])), " *
        (FAILURES[] == 0 ? "NO FAILURES" : "$(FAILURES[]) FAILURE(S)"))
    return FAILURES[]
end

# Run only when this file IS the program. `include`d — by the A/B harness in
# `copdem_dag_validate.jl`, or by anything else that wants to call `main` with a
# config of its own — it defines the run and starts nothing.
if abspath(PROGRAM_FILE) == @__FILE__
    exit(main() == 0 ? 0 : 1)
end
