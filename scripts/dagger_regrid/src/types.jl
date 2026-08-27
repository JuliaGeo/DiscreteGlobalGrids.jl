# Configuration and compact values that cross the coordinator/worker boundary.

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
    production === nothing && (production = copdem_config())
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

function validate_config(config::DaggerRegridConfig)
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
