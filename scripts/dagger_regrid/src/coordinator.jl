# Graph preparation, bounded admission, completion, and ledger ownership.

function prepare_coordinator(config)
    prod = config.production
    src = opensource(prod)
    chunks = destination_chunks(prod, src)
    capacity = 7^(prod.level - prod.ancestor)
    geometry_tag = prod.authalic ? sprint(show, src.sys7) : nothing
    store = openstore(prod, src.storesys7, capacity; geometry_tag)
    donepath = donelogpath(prod.store)
    done = prod.resume ? donechunks(donepath, prod.store, "elevation") : Set{Int}()
    plan = dagplan(src, chunks, prod)
    order = filter(d -> !(chunks[d] in done), plan.order)
    return (; src, chunks, order, graph = plan.graph, store, donepath,
        skipped = length(plan.order) - length(order))
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
