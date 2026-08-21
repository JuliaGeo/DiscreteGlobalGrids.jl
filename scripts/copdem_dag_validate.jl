# Does DAG-driven scheduling change the numbers? It must not.
#
# Ordering and caching decide WHEN a source tile is in memory, never what a
# destination cell is worth. So the same chunks regridded under the old
# schedule — canonical chunk order, striped LRU, no prefetch — and under the new
# one — tile-affinity order, refcount eviction, tapered batches, prefetcher —
# must produce byte-identical Zarr chunk files.
#
#     julia --project=bench -t 8 --gcthreads=4 scripts/copdem_dag_validate.jl [outdir]
#
# Writes two scratch stores under `outdir` (default a fresh temp directory) and
# compares them file by file. It never touches the production store: both stores
# are created here, and the only thing read from the real run's directory is
# nothing at all.

include(joinpath(@__DIR__, "copdem_production.jl"))

const OUT = length(ARGS) >= 1 ? ARGS[1] : mktempdir(; prefix = "dagvalidate-")

"""
    pickchunks(n) -> Vector{Int}

The chunks to compare: the highest-degree chunks in the covering — the polar
ones, which touch 100-360 tiles each and are the only place the multi-tile paths
run — plus an evenly spaced sample of the rest, so a single-tile chunk and an
all-ocean (all-`NaN`) chunk are in the set too.

Reads the cached covering and the dependency graph, not the store.
"""
function pickchunks(config, sys, sys7, tiles, todochunks, srcspace, n)
    plan = dagplan(sys, sys7, tiles, todochunks, srcspace, config)
    deg = [GR.sourcedegree(plan.graph, d) for d in 1:length(todochunks)]
    bydeg = sortperm(deg; rev = true)
    picked = Int[bydeg[1], bydeg[2]]                      # the two worst polar ones
    single = findfirst(==(1), deg)
    single === nothing || push!(picked, single)
    for d in round.(Int, range(1, length(todochunks); length = n))
        push!(picked, d)
    end
    ds = sort!(unique!(picked))
    say("validate: $(length(ds)) chunks, tile degrees " *
        string(sort([deg[d] for d in ds]; rev = true)))
    return [todochunks[d] for d in ds]
end

"Every Zarr chunk file a store holds, as path -> bytes."
function storefiles(path)
    dir = joinpath(path, "elevation")
    isdir(dir) || return Dict{String,Vector{UInt8}}()
    return Dict(f => read(joinpath(dir, f))
                for f in readdir(dir) if !startswith(f, "."))
end

function runside(base, tag, chunks, extra)
    store = joinpath(OUT, tag * ".zarr")
    FAILURES[] = 0
    config = merge(base, (store = store, chunks = chunks, resume = false,
        heartbeat = 1_000_000, checks = false, maxchunks = 0), extra)
    say("="^70)
    say("validate: side `$tag` -> $store")
    fails = main(config)
    return (store = store, failures = fails, files = storefiles(store))
end

function main_validate()
    base = merge(CONFIG, (dryrun = false,))
    sys = DGG.CopernicusDEMSystem(base.res)
    sys7 = DGG.IGeo7System()
    tilelist = joinpath(base.data, "CopernicusDEM", "tileList-glo$(base.res).txt")
    tiles = listedtiles(sys, tilelist, base.region)
    ids = TileIds(sys, tiles)
    srcspace = DGG.DGGSpace(DGG.PartialGrid(sys, 1, ids); chunklevel = 0)
    todochunks = load_chunklist(chunklistpath(CONFIG.store))
    todochunks === nothing && error("no cached covering to pick chunks from")

    chunks = pickchunks(base, sys, sys7, tiles, todochunks, srcspace, 8)

    legacy = runside(base, "legacy",
        chunks, (schedule = :canonical, cachepolicy = :lru, taper = false,
            prefetch = 0, batch = 8))
    dag = runside(base, "dag",
        chunks, (schedule = :affinity, cachepolicy = :refcount, taper = true,
            prefetch = 16, fetchconc = 4, batch = 8))

    say("="^70)
    FAILURES[] = 0
    check("the legacy side ran clean", legacy.failures == 0)
    check("the DAG side ran clean", dag.failures == 0)
    check("both sides wrote the same chunk files",
        sort(collect(keys(legacy.files))) == sort(collect(keys(dag.files)));
        detail = "$(length(legacy.files)) vs $(length(dag.files))")
    differing = [f for f in keys(legacy.files)
                 if get(dag.files, f, nothing) != legacy.files[f]]
    check("every chunk file is byte-identical", isempty(differing);
        detail = isempty(differing) ? "$(length(legacy.files)) files" :
                 "differ: " * join(sort(differing), ", "))
    bytes = sum(length, values(legacy.files); init = 0)
    say(@sprintf("validate: %d files, %.2f MiB compared, output at %s",
        length(legacy.files), bytes / 2^20, OUT))
    say(FAILURES[] == 0 ? "VALIDATION PASSED" : "$(FAILURES[]) VALIDATION FAILURE(S)")
    return FAILURES[]
end

exit(main_validate() == 0 ? 0 : 1)
