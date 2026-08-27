#!/usr/bin/env julia

# Run with:
#
#   julia --project=scripts/dagger_regrid -p 2 -t 1 \
#       scripts/dagger_regrid/main.jl
#
# DAGGER_REGRID_MODE may be smoke, canary, or run.
using Pkg
Pkg.activate(@__DIR__)
Base.active_project()

using DaggerRegrid
import Dagger
import Distributed
import Printf: @sprintf
Distributed.@everywhere using DaggerRegrid

using DaggerRegrid: BatchFlight, BatchReport, ChunkReport, DoneLog
using DaggerRegrid: GuidedSchedule, WorkerStateFactory, WorkerStats
using DaggerRegrid: claim_items!, collect_worker_stats, drain_flights!
using DaggerRegrid: launch_batch!, prepare_coordinator
using DaggerRegrid: record_batch!, record_transport_failures!
using DaggerRegrid: copdem_config, gcguard, hours, say, validate_config, worker_config

# Choose the execution mode before opening a source, graph, or store.
mode = Symbol(lowercase(get(ENV, "DAGGER_REGRID_MODE", "run")))

mode in (:smoke, :run, :canary) ||
    error("DAGGER_REGRID_MODE must be smoke, canary, or run; got $(repr(mode))")

batch = haskey(ENV, "DAGGER_BATCH") ? parse(Int, ENV["DAGGER_BATCH"]) : nothing

# Describe the run with plain, serializable values.
production = copdem_config()
config = DaggerRegridConfig(;
    production,
    inflight_per_process = parse(Int, get(ENV, "DAGGER_INFLIGHT", "1")),
    batch,
    worker_cache_slots = parse(Int, get(ENV, "DAGGER_CACHE_SLOTS", "256")),
    worker_cache_stripes = parse(Int, get(ENV, "DAGGER_CACHE_STRIPES", "16")),
)

mode === :canary && !(0 < config.production.maxchunks <= 32) &&
    error("canary mode requires COPDEM_MAXCHUNKS between 1 and 32")

if mode === :smoke
    smoke_reports = dagger_smoke(; processes = config.processes)
    println("Dagger smoke passed: ", smoke_reports)
    exit(0)
end

validate_config(config)
gcguard(config.production)

# The coordinator owns the graph, affinity order, store metadata, and ledger.
prepared = prepare_coordinator(config)
total = length(prepared.chunks)

if config.production.dryrun
    report = DaggerRegridReport(
        total, 0, prepared.skipped, 0, 0, 0, 0, 0.0,
        copy(config.processes), ChunkReport[], WorkerStats[])
    exit(0)
end

settings = worker_config(config)
state = Dagger.shard(
    WorkerStateFactory(settings);
    workers = config.processes,
)

nslots = length(config.processes) * config.inflight_per_process
schedule = GuidedSchedule(
    length(prepared.order),
    config.production.taper ? nslots : 1,
    config.batch,
)
slotpids = repeat(config.processes; inner = config.inflight_per_process)
flights = Union{Nothing,BatchFlight}[nothing for _ in slotpids]
ready = Channel{Int}(nslots)

failures = ChunkReport[]
completed = 0
cells = 0
nan = 0
batches = 0
started = time()
log = DoneLog(prepared.donepath)
finished = false

try
    # Fill the bounded admission window.
    for slot in eachindex(flights)
        items = claim_items!(schedule, prepared.order, prepared.chunks,
            config.failpoints)
        items === nothing && break
        launch_batch!(flights, ready, slot, state, slotpids[slot], items)
    end

    # Every completion frees one slot for the next guided batch.
    while any(!isnothing, flights)
        slot = take!(ready)
        flight = something(flights[slot])

        try
            batch_report = fetch(flight.task)::BatchReport
            delta = record_batch!(log, batch_report, failures)
            global completed += delta.completed
            global cells += delta.cells
            global nan += delta.nan
            global batches += 1
        catch err
            record_transport_failures!(
                failures, flight, err, catch_backtrace())
        end

        items = claim_items!(schedule, prepared.order, prepared.chunks,
            config.failpoints)
        flights[slot] = nothing
        items === nothing ||
            launch_batch!(flights, ready, slot, state, slotpids[slot], items)
    end
    global finished = true
finally
    finished || drain_flights!(flights, log, failures)
    close(log)
end

workerstats = collect_worker_stats(state, config.processes)
wall = time() - started
say(@sprintf(
    "DAGGER DONE  %d completed, %d skipped, %d failed, %d batches in %s",
    completed, prepared.skipped, length(failures), batches, hours(wall)))

report = DaggerRegridReport(
    total,
    completed,
    prepared.skipped,
    length(failures),
    cells,
    nan,
    batches,
    wall,
    copy(config.processes),
    failures,
    workerstats,
)

report.failed == 0 ||
    error("Dagger regrid finished with $(report.failed) failed chunks")
report
