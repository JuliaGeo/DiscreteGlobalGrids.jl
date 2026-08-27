#!/usr/bin/env julia

"""
An isolated Dagger execution backend for the CopDEM production regrid.

The regular production driver is loaded lazily for its source, graph, regrid,
and store seams; it is not modified and its `main` is never called here. The
coordinator retains the dependency graph, affinity order, admission policy, and
completion ledger. A Dagger task computes and writes a complete destination
chunk (or a small affinity-contiguous batch), preserving the existing ordered
source fold.

Typical local-process canary:

    julia --project=scripts/dagger_regrid -e \
        'using Pkg; Pkg.instantiate()'
    julia --project=scripts/dagger_regrid -p 2 -t 1 \
        scripts/dagger_regrid.jl

Set `COPDEM_MAXCHUNKS`, `COPDEM_STORE`, and the other existing CopDEM variables
as usual. `DAGGER_REGRID_MODE=smoke` runs only process-scoped echo tasks;
`DAGGER_REGRID_MODE=canary` additionally requires `COPDEM_MAXCHUNKS <= 32`.
"""
module DaggerRegrid

import Dagger
import Distributed
import Printf: @sprintf
using Base.ScopedValues: @with

export DaggerRegridConfig, DaggerRegridReport, dagger_regrid, dagger_smoke

const SCRIPT_PATH = abspath(@__FILE__)
const SOURCE_VERSION = "2026-08-26-dagger-regrid-d10"
const ALLOWED_FAILPOINTS = (:before_compute, :before_write, :after_write)
const PRODUCTION_LOCK = ReentrantLock()
const PRODUCTION_LOADED = Ref(false)

function _load_production!()
    PRODUCTION_LOADED[] && return nothing
    lock(PRODUCTION_LOCK) do
        PRODUCTION_LOADED[] && return nothing
        Base.include(@__MODULE__, joinpath(@__DIR__, "copdem_production.jl"))
        PRODUCTION_LOADED[] = true
    end
    return nothing
end

"""Configuration owned by the experimental Dagger coordinator."""
struct DaggerRegridConfig{P}
    production::P
    processes::Vector{Int}
    inflight_per_process::Int
    batch::Int
    worker_cache_slots::Int
    worker_cache_stripes::Int
    worker_tilecache_root::Union{Nothing,String}
    failpoints::Dict{Int,Symbol}
end

function DaggerRegridConfig(; production = nothing, processes = nothing,
        inflight_per_process::Integer = 1,
        batch::Union{Nothing,Integer} = nothing,
        worker_cache_slots::Integer = 256, worker_cache_stripes::Integer = 16,
        worker_tilecache_root::Union{Nothing,AbstractString} = nothing,
        failpoints = Dict{Int,Symbol}())
    if production === nothing
        _load_production!()
        production = CONFIG
    end
    batch === nothing && (batch = production.batch)
    pids = processes === nothing ? Distributed.workers() : Int.(processes)
    isempty(pids) && (pids = [Distributed.myid()])
    return DaggerRegridConfig(production, unique!(collect(pids)),
        Int(inflight_per_process), Int(batch), Int(worker_cache_slots),
        Int(worker_cache_stripes),
        worker_tilecache_root === nothing ? nothing : abspath(worker_tilecache_root),
        Dict{Int,Symbol}(failpoints))
end

struct WorkItem
    destination::Int
    chunk::Int
    failpoint::Union{Nothing,Symbol}
end

struct ChunkReport
    destination::Int
    chunk::Int
    cells::Int
    nan::Int
    seconds::Float64
    stage::Symbol
    error::Union{Nothing,String}
end

struct BatchReport
    worker::Int
    seconds::Float64
    chunks::Vector{ChunkReport}
end

struct WorkerStats
    worker::Int
    loads::Int
    hits::Int
    live::Int
    bytes::Int
    real_tiles::Int
    synthetic_tiles::Int
    pixels::Int
end

"""A compact, coordinator-side result; destination arrays never appear here."""
struct DaggerRegridReport
    total::Int
    completed::Int
    skipped::Int
    failed::Int
    cells::Int
    nan::Int
    batches::Int
    seconds::Float64
    processes::Vector{Int}
    failures::Vector{ChunkReport}
    workers::Vector{WorkerStats}
end

struct WorkerState{D,S,Y,O,L,C}
    dem::D
    srcspace::S
    sys7::Y
    store::O
    layout::L
    config::C
end

mutable struct BatchFlight
    task::Dagger.DTask
    process::Int
    items::Vector{WorkItem}
end

function _validate(config::DaggerRegridConfig)
    config.inflight_per_process > 0 ||
        throw(ArgumentError("inflight_per_process must be positive"))
    config.batch > 0 || throw(ArgumentError("batch must be positive"))
    config.worker_cache_slots > 0 ||
        throw(ArgumentError("worker_cache_slots must be positive"))
    config.worker_cache_stripes > 0 ||
        throw(ArgumentError("worker_cache_stripes must be positive"))
    isempty(config.processes) && throw(ArgumentError("processes must not be empty"))
    available = Set(Distributed.procs())
    missing = filter(p -> !(p in available), config.processes)
    isempty(missing) || throw(ArgumentError(
        "Dagger processes $(missing) do not exist; available processes are $(sort!(collect(available)))"))
    config.production.prefetch == 0 || throw(ArgumentError(
        "the isolated backend does not yet implement graph prefetch; set prefetch = 0"))
    for (chunk, point) in config.failpoints
        point in ALLOWED_FAILPOINTS || throw(ArgumentError(
            "failpoint for chunk $chunk is $point; expected one of $(ALLOWED_FAILPOINTS)"))
    end
    regridmethod(config.production)
    return config
end

function _load_workers!(pids; production = false)
    for pid in pids
        pid == Distributed.myid() && continue
        loaded = Distributed.remotecall_fetch(isdefined, pid, Main, :DaggerRegrid)
        if loaded
            version = Distributed.remotecall_fetch(Core.eval, pid, Main,
                :(DaggerRegrid.SOURCE_VERSION))
            version == SOURCE_VERSION || error(
                "process $pid has DaggerRegrid $version loaded, expected $SOURCE_VERSION; use fresh workers")
        else
            Distributed.remotecall_wait(Base.include, pid, Main, SCRIPT_PATH)
        end
        if production
            # Load production outside a Dagger task. The Dagger call sites also
            # enter through `invokelatest`, because a smoke task may already
            # have started a reusable loop with an older Julia 1.12 world age.
            Distributed.remotecall_wait(Core.eval, pid, Main,
                :(DaggerRegrid._load_production!()))
        end
    end
    return nothing
end

function _source_geometry(config)
    sys = DGG.CopernicusDEMSystem(config.res)
    storesys7 = DGG.IGeo7System()
    sys7 = config.authalic ? DGG.AuthalicSystem(storesys7) : storesys7
    tilelist = joinpath(config.data, "CopernicusDEM", "tileList-glo$(config.res).txt")
    tiles = listedtiles(sys, tilelist, config.region)
    isempty(tiles) && error("the tile list and region select no tiles")
    ids = TileIds(sys, tiles)
    srcgrid = DGG.PartialGrid(sys, 1, ids)
    srcspace = DGG.DGGSpace(srcgrid; chunklevel = 0)
    return (; sys, storesys7, sys7, tiles, ids, srcspace)
end

function _destination_chunks(config, geometry)
    !isempty(config.chunks) && return sort!(unique(copy(config.chunks)))
    path = chunklistpath(config.store)
    chunks = load_chunklist(path)
    if chunks === nothing
        chunks = covering_chunks(geometry.sys7, geometry.sys, geometry.tiles,
            config.ancestor; nthreads = max(1, Threads.nthreads() - 1))
        save_chunklist(path, config.ancestor, chunks)
    end
    if 0 < config.maxchunks < length(chunks)
        chunks = chunks[round.(Int, range(1, length(chunks); length = config.maxchunks))]
    end
    return chunks
end

function _prepare_coordinator(config)
    prod = config.production
    geometry = _source_geometry(prod)
    chunks = _destination_chunks(prod, geometry)
    capacity = 7^(prod.level - prod.ancestor)
    geometry_tag = prod.authalic ? sprint(show, geometry.sys7) : nothing
    store = openstore(prod, geometry.storesys7, capacity; geometry_tag)
    layout = store.layout
    donepath = donelogpath(prod.store)
    done = prod.resume ? donechunks(donepath, prod.store, "elevation") : Set{Int}()
    plan = dagplan(geometry.sys, geometry.sys7, geometry.tiles, chunks,
        geometry.srcspace, prod)
    GR.ndestinationchunks(plan.graph) == length(chunks) ||
        error("dependency graph destination count does not match the work list")
    GR.nsourcechunks(plan.graph) == length(geometry.tiles) ||
        error("dependency graph source count does not match the tile list")
    length(plan.order) == length(chunks) && isperm(plan.order) ||
        error("dependency order is not a permutation of the work list")
    order = filter(d -> !(chunks[d] in done), plan.order)
    return (; geometry..., chunks, order, graph = plan.graph, store, layout,
        donepath, skipped = length(plan.order) - length(order))
end

function _worker_config(config::DaggerRegridConfig)
    return (; production = config.production,
        cache_slots = config.worker_cache_slots,
        cache_stripes = config.worker_cache_stripes,
        tilecache_root = config.worker_tilecache_root)
end

function _make_worker_state(settings)
    _load_production!()
    config = settings.production
    gcguard(config)
    config.malloctrim > 0 && tunemalloc(config.malloctrim)
    geometry = _source_geometry(config)
    tiledir = joinpath(config.data, "CopernicusDEM", "$(config.res)m")
    landshp = joinpath(config.data, "naturalearth", "ne_10m_land.shp")
    mask = landmask(landshp, config.maskarcsec)
    real = realtiles(geometry.sys, tiledir,
        effective_realspec(config.source, config.real))
    tileset = Set(geometry.tiles)
    filter!(p -> p.first in tileset, real)
    provider = if config.source === :real
        root = something(settings.tilecache_root, config.tilecache * ".dagger")
        cachedir = joinpath(root, "worker-$(Distributed.myid())")
        LazyCopernicusTiles(geometry.sys, geometry.tiles;
            cachedir, baseurl = config.tilebaseurl, retries = config.retries,
            backoff = config.backoff, timeout = config.timeout)
    else
        nothing
    end
    builder = TileBuilder(geometry.sys, geometry.tiles, real, provider, mask;
        delay = config.fetchdelay)
    cache = StripedLRUCache{Vector{Float32}}(k -> buildtile(builder, k);
        slots = settings.cache_slots, stripes = settings.cache_stripes)
    dem = TiledDEM(geometry.ids, builder, cache)
    store = DGG.subzonestore(config.store)
    return WorkerState(dem, geometry.srcspace, geometry.sys7, store, store.layout, config)
end

function _error_string(err, bt)
    text = sprint(showerror, err, bt)
    return first(text, min(length(text), 1500))
end

function _worker_batch(state::WorkerState, items::Vector{WorkItem})
    reports = ChunkReport[]
    batch_started = time()
    outer = state.config.shape === :outer
    @with GR.OUTER_PARALLEL => outer begin
        for item in items
            started = time()
            stage = :compute
            vals = nothing
            try
                item.failpoint === :before_compute && error("injected before-compute failure")
                vals = regrid_chunk(state.dem, state.srcspace, state.sys7,
                    state.layout, item.chunk, state.config)
                stage = :write
                item.failpoint === :before_write && error("injected before-write failure")
                DGG.dggwrite!(state.store, item.chunk, vals)
                stage = :report
                item.failpoint === :after_write && error("injected after-write failure")
                push!(reports, ChunkReport(item.destination, item.chunk, length(vals),
                    count(isnan, vals), time() - started, :done, nothing))
            catch err
                ncells = vals === nothing ? 0 : length(vals)
                nnan = vals === nothing ? 0 : count(isnan, vals)
                push!(reports, ChunkReport(item.destination, item.chunk, ncells,
                    nnan, time() - started, stage,
                    _error_string(err, catch_backtrace())))
            end
        end
    end
    return BatchReport(Distributed.myid(), time() - batch_started, reports)
end

function _worker_stats(state::WorkerState)
    stats = cachestats(state.dem.cache)
    builder = state.dem.builder
    return WorkerStats(Distributed.myid(), stats.loads, stats.hits, stats.live,
        stats.bytes, builder.nreal[], builder.nsynthetic[], builder.npixels[])
end

_latest(f, args...) = Base.invokelatest(f, args...)

function _spawn_batch(state, pid, items)
    firstchunk, lastchunk = first(items).chunk, last(items).chunk
    options = Dagger.Options(; scope = Dagger.scope(; worker = pid),
        name = "regrid-$firstchunk-$lastchunk", return_type = BatchReport)
    return BatchFlight(Dagger.spawn(_latest, options, _worker_batch, state, items),
        pid, items)
end

function _launch!(flights, ready, slot, state, pid, items)
    flight = _spawn_batch(state, pid, items)
    flights[slot] = flight
    errormonitor(@async begin
        wait(flight.task)
        put!(ready, slot)
    end)
    return flight
end

function _claim_items(schedule, order, chunks, failpoints)
    positions = claim!(schedule)
    positions === nothing && return nothing
    return [WorkItem(order[p], chunks[order[p]],
        get(failpoints, chunks[order[p]], nothing)) for p in positions]
end

function _record_batch!(log, report, failures)
    completed = cells = nan = 0
    for chunk in report.chunks
        if chunk.error === nothing
            record!(log, chunk.chunk, chunk.cells, chunk.nan, chunk.seconds, report.worker)
            completed += 1
            cells += chunk.cells
            nan += chunk.nan
        else
            push!(failures, chunk)
            say("ERROR process $(report.worker) chunk $(chunk.chunk) at $(chunk.stage): ",
                chunk.error)
        end
    end
    return (; completed, cells, nan)
end

function _transport_failures!(failures, flight, err, bt)
    message = _error_string(err, bt)
    for item in flight.items
        push!(failures, ChunkReport(item.destination, item.chunk, 0, 0, 0.0,
            :dagger, message))
    end
    say("ERROR Dagger task on process $(flight.process): ", message)
    return nothing
end

function _drain_flights!(flights, log, failures)
    for flight in flights
        flight === nothing && continue
        try
            _record_batch!(log, fetch(flight.task)::BatchReport, failures)
        catch err
            _transport_failures!(failures, flight, err, catch_backtrace())
        end
    end
    return nothing
end

function _final_worker_stats(state, pids)
    tasks = [Dagger.spawn(_latest,
        Dagger.Options(; scope = Dagger.scope(; worker = pid),
            name = "regrid-stats-$pid", return_type = WorkerStats),
        _worker_stats, state) for pid in pids]
    return WorkerStats[fetch(task) for task in tasks]
end

"""
    dagger_regrid(config = DaggerRegridConfig()) -> DaggerRegridReport

Run the graph-aware CopDEM destination order through a bounded set of
process-scoped Dagger tasks. Worker state is created once per process as a
`Dagger.Shard`; output arrays are written on the worker, and only compact
reports return to the coordinator. The coordinator alone appends the done
ledger.
"""
function dagger_regrid(config::DaggerRegridConfig = DaggerRegridConfig())
    _load_production!()
    _validate(config)
    gcguard(config.production)
    prepared = _prepare_coordinator(config)
    total = length(prepared.chunks)
    if config.production.dryrun
        return DaggerRegridReport(total, 0, prepared.skipped, 0, 0, 0, 0, 0.0,
            copy(config.processes), ChunkReport[], WorkerStats[])
    end

    # Dagger's eager context sees every Distributed worker, even though our
    # process scopes use only `config.processes`; all workers need the task
    # module, while only selected workers need the production helpers.
    _load_workers!(Distributed.workers())
    _load_workers!(config.processes; production = true)
    settings = _worker_config(config)
    state = Dagger.shard(
        () -> Base.invokelatest(_make_worker_state, settings);
        workers = config.processes)
    nslots = length(config.processes) * config.inflight_per_process
    schedule = GuidedSchedule(length(prepared.order),
        config.production.taper ? nslots : 1, config.batch)
    slotpids = repeat(config.processes; inner = config.inflight_per_process)
    flights = Union{Nothing,BatchFlight}[nothing for _ in slotpids]
    ready = Channel{Int}(nslots)
    failures = ChunkReport[]
    completed = cells = nan = batches = 0
    started = time()
    log = DoneLog(prepared.donepath)
    finished = false

    try
        for i in eachindex(flights)
            items = _claim_items(schedule, prepared.order, prepared.chunks,
                config.failpoints)
            items === nothing && break
            _launch!(flights, ready, i, state, slotpids[i], items)
        end
        while any(!isnothing, flights)
            i = take!(ready)
            flight = something(flights[i])
            try
                report = fetch(flight.task)::BatchReport
                delta = _record_batch!(log, report, failures)
                completed += delta.completed
                cells += delta.cells
                nan += delta.nan
                batches += 1
            catch err
                _transport_failures!(failures, flight, err, catch_backtrace())
            end
            items = _claim_items(schedule, prepared.order, prepared.chunks,
                config.failpoints)
            flights[i] = nothing
            items === nothing || _launch!(flights, ready, i, state, slotpids[i], items)
        end
        finished = true
    finally
        finished || _drain_flights!(flights, log, failures)
        close(log)
    end

    workerstats = _final_worker_stats(state, config.processes)
    wall = time() - started
    say(@sprintf("DAGGER DONE  %d completed, %d skipped, %d failed, %d batches in %s",
        completed, prepared.skipped, length(failures), batches, hours(wall)))
    return DaggerRegridReport(total, completed, prepared.skipped, length(failures),
        cells, nan, batches, wall, copy(config.processes), failures, workerstats)
end

dagger_regrid(production; kwargs...) =
    dagger_regrid(DaggerRegridConfig(; production, kwargs...))

_smoke_state() = (worker = Distributed.myid(), version = SOURCE_VERSION)
_smoke_echo(state, value) = (worker = Distributed.myid(),
    state_worker = state.worker, version = state.version, value)

"""Run only small process-scoped Dagger tasks; no source, graph, or store is opened."""
function dagger_smoke(; processes = nothing)
    pids = processes === nothing ? Distributed.workers() : Int.(processes)
    isempty(pids) && (pids = [Distributed.myid()])
    pids = unique!(collect(pids))
    available = Set(Distributed.procs())
    all(pid -> pid in available, pids) ||
        throw(ArgumentError("smoke processes must be members of $(sort!(collect(available)))"))
    _load_workers!(Distributed.workers())
    state = Dagger.shard(_smoke_state; workers = pids)
    tasks = [Dagger.spawn(_smoke_echo,
        Dagger.Options(; scope = Dagger.scope(; worker = pid), name = "smoke-$pid"),
        state, pid) for pid in pids]
    reports = fetch.(tasks)
    all(report -> report.worker == report.state_worker == report.value &&
        report.version == SOURCE_VERSION, reports) ||
        error("Dagger process scope was not honored: $reports")
    return reports
end

function run_cli()
    mode = Symbol(lowercase(get(ENV, "DAGGER_REGRID_MODE", "run")))
    if mode === :smoke
        reports = dagger_smoke()
        println("Dagger smoke passed: ", reports)
        return 0
    end
    _load_production!()
    mode in (:run, :canary) || error(
        "DAGGER_REGRID_MODE must be smoke, canary, or run; got $(repr(mode))")
    mode === :canary && !(0 < CONFIG.maxchunks <= 32) && error(
        "canary mode requires COPDEM_MAXCHUNKS between 1 and 32")
    config = DaggerRegridConfig(
        inflight_per_process = parse(Int, get(ENV, "DAGGER_INFLIGHT", "1")),
        batch = parse(Int, get(ENV, "DAGGER_BATCH", string(CONFIG.batch))),
        worker_cache_slots = parse(Int, get(ENV, "DAGGER_CACHE_SLOTS", "256")),
        worker_cache_stripes = parse(Int, get(ENV, "DAGGER_CACHE_STRIPES", "16")))
    report = dagger_regrid(config)
    return report.failed == 0 ? 0 : 1
end

end # module DaggerRegrid

if abspath(PROGRAM_FILE) == @__FILE__
    exit(DaggerRegrid.run_cli())
end
