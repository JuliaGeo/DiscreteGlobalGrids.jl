# Copernicus DEM GLO-90 -> IGEO7 level 12, into one global ancestor-subzone Zarr
# store, written by W worker tasks in one process.
#
#     COPDEM_SOURCE=real COPDEM_STORE=/path/to/store.zarr \
#         julia --project=benchmark -t 26 --gcthreads=8 scripts/copdem_production.jl
#
# `copdem_config` in `copdem_common.jl` lists every knob and its `COPDEM_*`
# variable. `--gcthreads` takes one field; `gcguard` refuses `--gcthreads=N,1`.
#
# The run is resumable: it skips every chunk the ledger or the store already has,
# so re-running the same command continues where it stopped.
#
#   tile    a 1x1-degree Copernicus DEM tile. The source unit: one source chunk.
#   chunk   one level-`ancestor` IGeo7 cell with all its level-`level`
#           descendants (7^7 = 823 543 cells at 5 -> 12). The destination work
#           unit, one Zarr chunk, one file. The store API calls it a "column".
#   budget  bytes of intersection weights a lazy regrid may hold at once, per
#           worker.
#
# The tile list is real for both sources: Copernicus ships ~26 450 of the 64 800
# tiles in the 1x1-degree lattice, and the source is a `PartialGrid` over exactly
# those. A destination chunk over open ocean pairs with no source and stays
# nodata without a network request.
const STARTED = Ref(time())

import Profile
include("copdem_common.jl")

# A one-second libuv wake keeps Julia's signal-profile report listener schedulable
# while every default-pool thread is occupied by the worker wave.
const PROFILEPUMP = Ref{Union{Nothing,Task}}(nothing)

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

Current resident set size in GiB. `Sys.maxrss()` is a monotone high-water mark,
so a heartbeat printing it shows one transient peak for the rest of the run;
[`peakrssgib`](@ref) reports that, labelled as a peak.
"""
function rssgib()
    try
        return parse(Int, split(read("/proc/self/statm", String))[2]) * 4096 / 2^30
    catch
        return Sys.maxrss() / 2^30   # off Linux the peak is all there is
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
    etaseconds(elapsed, computed, remaining) -> Float64

Time left at the rate this session has computed chunks; `NaN` until the first
one lands. `computed` excludes chunks skipped at resume: those cost a set
lookup, and crediting them to the session clock makes the rate arbitrarily fast.
"""
etaseconds(elapsed, computed, remaining) =
    computed > 0 ? elapsed * remaining / computed : NaN

"""
    heartbeat!(p, force = false)

One progress line, at most every `p.every` seconds. Every rate in it is
session-scoped; only the chunk count is store-wide. Store-wide cells and NaNs
come from the ledger at the end of the run; see [`storetotals`](@ref).
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
    runchunks(chunks, order, sched, dem, cache, src, store, p, donelog, prefetcher, config)

`config.workers` tasks pulling batches off one [`GuidedSchedule`](@ref), each
computing its chunk and writing it straight to its own Zarr file. Chunks are
disjoint, so no two tasks touch the same file.

`order[pos]` is the graph destination to run at position `pos`, and
`chunks[order[pos]]` is the level-`ancestor` chunk itself. The order is a
priority sequence, never a partition: pulling dynamically keeps the polar
chunks, which meet 200-360 tiles and cost ~5x a mid-latitude one, off the
critical path.

Every chunk is retired from the cache exactly once, on every terminal outcome
and before its result is written: the Zarr write reads no source data.

Each worker sets `GR.OUTER_PARALLEL`, so nested weight builds stay serial and
the lazy plan's source wave is the only parallelism inside a chunk.
"""
function runchunks(chunks, order, sched, dem, cache, src, store, p, donelog,
        prefetcher, config)
    layout = store.layout
    log = DoneLog(donelog)
    try
        Threads.@sync for w in 1:config.workers
            Threads.@spawn @with GR.OUTER_PARALLEL => true begin
                while true
                    batch = claim!(sched)
                    batch === nothing && break
                    advance!(prefetcher)
                    for pos in batch
                        d = order[pos]
                        ch = chunks[d]
                        t0 = time()
                        vals = try
                            regrid_chunk(dem, src.srcspace, src.sys7, layout, ch, config)
                        catch err
                            say("ERROR worker $w chunk $ch: ",
                                first(sprint(showerror, err, catch_backtrace()), 1500))
                            Threads.atomic_add!(FAILURES, 1)
                            nothing
                        finally
                            retire_column!(cache, d)
                        end
                        vals === nothing && continue
                        config.write && DGG.dggwrite!(store, ch, vals)
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
    reportcache(cache, prefetcher, tilesource)

What the tile cache did, and the laws it must keep:

  - **Every source chunk was loaded at most once.** A chunk is freed only when
    no destination chunk that may read it is left, so no correct demand for it
    can follow. A second load means the refcounts are wrong.
  - **Every demand was one the graph predicted.** The adjacency is a superset of
    what the executor reads, so an unpredicted read is a hole in the covering.
    It is still served, and counted here.
  - **Every destination chunk was retired**, which leaves nothing pinned.
"""
function reportcache(cache, prefetcher, tilesource)
    s = cachestats(cache)
    say(@sprintf("tile cache: %d loads of %d source chunks, peak %.2f GiB / %d tiles, %d hits, %d joined loads, %d live at end",
        s.loads, nsourcechunks(cache), s.peakbytes / 2^30, s.peaktiles, s.hits,
        s.waits, s.live))
    check("every source chunk was loaded at most once", isempty(doubleloaded(cache));
        detail = "$(length(doubleloaded(cache))) chunk(s) loaded twice")
    check("every demand was a predicted edge", s.uncredited == 0;
        detail = s.uncredited == 0 ? "" :
                 "$(s.uncredited) uncredited demand(s), chunks $(s.uncreditedof)")
    check("every destination chunk was retired", s.pinned == 0 && s.live == 0;
        detail = "$(s.pinned) source chunk(s) still pinned, $(s.live) still held")
    tilesource isa CopernicusTiles && say(
        "source downloads: $(tilesource.ndownloads[]) total, " *
        "$(tilesource.ncold[]) demand-cold")
    pf = prefetchstats(prefetcher)
    pf.depth == 0 && return nothing
    say("prefetch: depth $(pf.depth), $(pf.issued) tile requests issued")
    check("the prefetcher raised nothing", pf.failure === nothing;
        detail = pf.failure === nothing ? "" : sprint(showerror, pf.failure))
    return nothing
end

function main(config = copdem_config())
    gcguard()
    configureprofile()
    println("="^92)
    println(stamp(), "  copdem_production.jl — GLO-$(config.res) -> IGEO7 level " *
                     "$(config.level), level-$(config.ancestor) chunks, " *
                     "$(uppercase(String(config.source))) elevations, " *
                     "$(uppercase(String(config.method))) method")
    println("  julia $(VERSION)  threads=$(Threads.nthreads())  " *
            "gcmark=$(Base.JLOptions().nmarkthreads)  pid=$(getpid())")
    println("  ", config)
    println("="^92)
    flush(stdout)
    Sys.islinux() && !tunemalloc() &&
        say("malloc: mallopt failed, expect ~3x resident memory")

    # --- the source and the work list ------------------------------------
    src = opensource(config)
    tiles = src.tiles
    say("source: $(length(tiles)) land tiles of $(DGG.ncells(src.sys, 0)), " *
        @sprintf("%.3e pixels", Float64(length(src.ids))))
    chunks = destination_chunks(config, src)
    capacity = 7^(config.level - config.ancestor)
    say(@sprintf("work: %d chunks x %d cells = %.4e level-%d cells",
        length(chunks), capacity, Float64(length(chunks) * capacity), config.level))

    # --- the store -------------------------------------------------------
    geometry_tag = config.authalic ? sprint(show, src.sys7) : nothing
    store = openstore(config, src.storesys7, capacity; geometry_tag)
    donelog = donelogpath(config.store)
    done = config.resume ? donechunks(donelog, config.store, "elevation") : Set{Int}()

    # --- the dependency graph, the walk order, the cache ------------------
    t0 = time()
    plan = dagplan(src, chunks, config)
    graph = plan.graph
    say("graph: $(GR.nsourcechunks(graph)) x $(GR.ndestinationchunks(graph)) chunks, " *
        secs(time() - t0))
    W = config.workers
    cache = RefCountCache{Vector{Float32}}(length(tiles), length(chunks),
        d -> GR.sourcesof(graph, d), s -> GR.consumerdegree(graph, s),
        k -> loadtile(src.tilesource, tiles[k]); permits = W)
    dem = TiledDEM(src.ids, k -> gettile!(cache, k))

    # A chunk skipped on resume is retired before anything starts, so its tiles
    # are never credited to work that will not happen.
    order = filter(plan.order) do d
        chunks[d] in done || return true
        retire_column!(cache, d)
        return false
    end
    nskipped = length(plan.order) - length(order)
    say("resume: $nskipped chunks already written, $(length(order)) to do")

    if config.dryrun
        say(@sprintf("dryrun: max tiles per chunk %d, max chunks per tile %d",
            maximum(d -> GR.sourcedegree(graph, d), 1:length(chunks); init = 0),
            maximum(s -> GR.consumerdegree(graph, s), 1:length(tiles); init = 0)))
        return FAILURES[]
    end

    # --- the run ---------------------------------------------------------
    say("launching W=$W workers (nthreads=$(Threads.nthreads())), " *
        "batch=$(config.batch), budget=$(config.budget)")
    sched = GuidedSchedule(length(order), W, config.batch)
    p = Progress(length(chunks), config.heartbeat,
        src.tilesource isa CopernicusTiles ? src.tilesource.ncold : nothing)
    p.skipped = nskipped
    prefetcher = nothing
    t0 = time()
    try
        if config.prefetch > 0
            prepare = src.tilesource isa CopernicusTiles ?
                      (s -> tilepath!(src.tilesource, tiles[s]; demand = false)) : nothing
            prefetcher = Prefetcher(cache, order, d -> GR.sourcesof(graph, d),
                length(tiles), sched; depth = config.prefetch, concurrency = W, prepare)
            say("prefetch: depth $(config.prefetch) chunks, $W concurrent fetches")
        end
        runchunks(chunks, order, sched, dem, cache, src, store, p, donelog,
            prefetcher, config)
    finally
        stop_prefetch!(prefetcher)
    end
    wall = time() - t0

    # The session counters describe this process; the store totals span every
    # session that ever wrote to the store. A resumed run reports both.
    written = donechunks(donelog, config.store, "elevation"; label = "store")
    tot = storetotals(donelog, written)
    say(@sprintf("RUN DONE  session: %d chunks computed, %d skipped, %.4e cells in %s = %.0f cells/s, NaN %.3f%%",
        p.done, p.skipped, Float64(p.cells), hours(wall), p.cells / max(wall, 1e-9),
        100 * p.nan / max(p.cells, 1)))
    say(@sprintf("RUN DONE  store total: %d of %d chunks written, %.4e cells, NaN %.3f%%%s",
        tot.chunks, length(chunks), Float64(tot.cells),
        100 * tot.nan / max(tot.cells, 1),
        tot.unaccounted == 0 ? "" :
            "; $(tot.unaccounted) written chunks have no ledger line and are not counted"))
    reportcache(cache, prefetcher, src.tilesource)
    say(@sprintf("RSS %.2f GiB now, %.2f GiB peak", rssgib(), peakrssgib()))

    config.checks && config.source === :synthetic &&
        verify(config, src, store.layout, chunks)

    say("total wall $(hours(time() - STARTED[])), " *
        (FAILURES[] == 0 ? "NO FAILURES" : "$(FAILURES[]) FAILURE(S)"))
    return FAILURES[]
end

# `include`d, this file defines the run and starts nothing.
if abspath(PROGRAM_FILE) == @__FILE__
    exit(main() == 0 ? 0 : 1)
end
