# Small process-placement check; no source, graph, or store is opened.

_smoke_state() = (worker = Distributed.myid(),)
_smoke_echo(state, value) = (worker = Distributed.myid(),
    state_worker = state.worker, value)

"""Run only small process-scoped Dagger tasks; no source, graph, or store is opened."""
function dagger_smoke(; processes = nothing)
    pids = processes === nothing ? Distributed.workers() : Int.(processes)
    isempty(pids) && (pids = [Distributed.myid()])
    pids = unique!(collect(pids))
    available = Set(Distributed.procs())
    all(pid -> pid in available, pids) ||
        throw(ArgumentError("smoke processes must be members of $(sort!(collect(available)))"))
    state = Dagger.shard(_smoke_state; workers = pids)
    tasks = [Dagger.spawn(_smoke_echo,
        Dagger.Options(; scope = Dagger.scope(; worker = pid), name = "smoke-$pid"),
        state, pid) for pid in pids]
    reports = fetch.(tasks)
    all(report -> report.worker == report.state_worker == report.value, reports) ||
        error("Dagger process scope was not honored: $reports")
    return reports
end
