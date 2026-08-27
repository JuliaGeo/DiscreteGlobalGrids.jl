# Graph preparation, bounded admission, completion, and ledger ownership.

function _source_geometry(config)
    sys = DGG.CopernicusDEMSystem(config.res)
    storesys7 = DGG.IGeo7System()
    sys7 = config.authalic ? DGG.AuthalicSystem(storesys7) : storesys7
    tilelist = copdem_tilelist(config)
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

function prepare_coordinator(config)
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


function spawn_batch(state, pid, items)
    firstchunk, lastchunk = first(items).chunk, last(items).chunk
    options = Dagger.Options(; scope = Dagger.scope(; worker = pid),
        name = "regrid-$firstchunk-$lastchunk", return_type = BatchReport)
    return BatchFlight(Dagger.spawn(run_batch, options, state, items), pid, items)
end

function launch_batch!(flights, ready, slot, state, pid, items)
    flight = spawn_batch(state, pid, items)
    flights[slot] = flight
    errormonitor(@async begin
        wait(flight.task)
        put!(ready, slot)
    end)
    return flight
end

function claim_items!(schedule, order, chunks, failpoints)
    positions = claim!(schedule)
    positions === nothing && return nothing
    return [WorkItem(order[p], chunks[order[p]],
        get(failpoints, chunks[order[p]], nothing)) for p in positions]
end

function record_batch!(log, report, failures)
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

function record_transport_failures!(failures, flight, err, bt)
    message = _error_string(err, bt)
    for item in flight.items
        push!(failures, ChunkReport(item.destination, item.chunk, 0, 0, 0.0,
            :dagger, message))
    end
    say("ERROR Dagger task on process $(flight.process): ", message)
    return nothing
end

function drain_flights!(flights, log, failures)
    for flight in flights
        flight === nothing && continue
        try
            record_batch!(log, fetch(flight.task)::BatchReport, failures)
        catch err
            record_transport_failures!(failures, flight, err, catch_backtrace())
        end
    end
    return nothing
end

function collect_worker_stats(state, pids)
    tasks = [Dagger.spawn(worker_stats,
        Dagger.Options(; scope = Dagger.scope(; worker = pid),
            name = "regrid-stats-$pid", return_type = WorkerStats), state) for pid in pids]
    return WorkerStats[fetch(task) for task in tasks]
end
