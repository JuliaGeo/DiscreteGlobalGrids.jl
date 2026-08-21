# Controlled cold-network acceptance run for the real CopDEM DAG prefetcher.
#
# The three phases deliberately run as separate Julia processes so cold/warm
# wall times carry comparable JIT cost:
#
#   nice -n 10 julia --project=bench -t 8 --gcthreads=4 \
#       scripts/copdem_prefetch_coldtest.jl plan
#   nice -n 10 julia --project=bench -t 8 --gcthreads=4 \
#       scripts/copdem_prefetch_coldtest.jl cold
#   nice -n 10 julia --project=bench -t 8 --gcthreads=4 \
#       scripts/copdem_prefetch_coldtest.jl warm
#
# `cold` refuses an existing cache or output store. `warm` reuses that cache,
# writes a second store, and compares its encoded elevation chunks byte for
# byte with the cold store. Nothing here removes an input or output.

include(joinpath(@__DIR__, "copdem_production.jl"))
import JSON3

const COLDTEST_SCRATCH = "/home/asinghvi17/geo/scratch-stores"
const COLDTEST_DATA = "/home/asinghvi17/geo/DiscreteGlobalGrids.jl/bench/data"
const COLDTEST_CACHE = joinpath(COLDTEST_SCRATCH, "prefetch-coldtest-cache")
const COLDTEST_COLD_STORE = joinpath(COLDTEST_SCRATCH, "prefetch-coldtest-cold.zarr")
const COLDTEST_WARM_STORE = joinpath(COLDTEST_SCRATCH, "prefetch-coldtest-warm.zarr")
const COLDTEST_RESULTS = joinpath(COLDTEST_SCRATCH, "prefetch-coldtest-results.ndjson")

# Sorted because a PartialGrid's concatenated ancestor subtrees must be ordered.
# These are real-land Himalayan columns centred near 88.72E/28.75N,
# 84.39E/28.20N, and 86.45E/28.03N.
const COLDTEST_COLUMNS = [70493, 73037, 73054]
const COLDTEST_PREFETCH = 32
const COLDTEST_FETCHCONC = 16
const COLDTEST_WORKERS = 3

coldcheck(name, ok; detail = "") = ok ?
    println("PASS  ", rpad(name, 58), detail) :
    error("FAILED: $name" * (isempty(detail) ? "" : " — $detail"))

"The exact broad-phase dependency set the production executor will schedule."
function coldtest_plan()
    sys = DGG.CopernicusDEMSystem(90)
    sys7 = DGG.IGeo7System()
    tilelist = joinpath(COLDTEST_DATA, "CopernicusDEM", "tileList-glo90.txt")
    tiles = listedtiles(sys, tilelist, nothing)
    ids = TileIds(sys, tiles)
    srcspace = DGG.DGGSpace(DGG.PartialGrid(sys, 1, ids); chunklevel = 0)
    config = merge(CONFIG, (source = :real, real = :none, data = COLDTEST_DATA,
        schedule = :affinity, refinegraph = false))
    dag = dagplan(sys, sys7, tiles, COLDTEST_COLUMNS, srcspace, config)
    percolumn = Dict(c => Int.(collect(GR.sourcesof(dag.graph, d)))
        for (d, c) in enumerate(COLDTEST_COLUMNS))
    sourceindices = sort!(unique!(Int[s for ss in values(percolumn) for s in ss]))
    ordinals = [tiles[s] for s in sourceindices]
    stems = [tilestem(sys, DGG.LevelIndex(0, ordinal)) for ordinal in ordinals]
    return (; sys, sys7, tiles, dag, percolumn, sourceindices, ordinals, stems)
end

"Content-Length from a HEAD request; never creates a cache object."
function coldtest_headbytes(plan)
    provider = LazyCopernicusTiles(plan.sys, plan.tiles; cachedir = COLDTEST_CACHE)
    out = Dict{String,Int}()
    for (stem, ordinal) in zip(plan.stems, plan.ordinals)
        response = Downloads.request(tileurl(provider, ordinal);
            method = "HEAD", timeout = 60.0)
        response.status == 200 || error("HEAD for $stem returned $(response.status)")
        entry = findfirst(p -> lowercase(first(p)) == "content-length",
            response.headers)
        entry === nothing && error("HEAD for $stem returned no Content-Length")
        out[stem] = parse(Int, last(response.headers[entry]))
    end
    return out
end

function print_coldtest_plan(plan, headbytes)
    println("columns: ", join(COLDTEST_COLUMNS, ", "))
    println("graph: ", plan.dag.edges, " edges, ", length(plan.stems),
        " unique listed tiles")
    for (d, column) in enumerate(COLDTEST_COLUMNS)
        stems = [tilestem(plan.sys, DGG.LevelIndex(0, plan.tiles[s]))
            for s in plan.percolumn[column]]
        println("  column $column: $(length(stems)) tiles — ", join(stems, ", "))
    end
    println("unique predicted set:")
    for stem in plan.stems
        println("  ", stem, "  ", headbytes[stem], " bytes")
    end
    println("predicted transfer: ", sum(values(headbytes)), " bytes (",
        round(sum(values(headbytes)) / 1e6; digits = 3), " MB)")
end

mutable struct ColdDownloadWatch
    const stems::Vector{String}
    const started::Float64
    const stop::Threads.Atomic{Bool}
    const partseen::Dict{String,Float64}
    const finished::Dict{String,Float64}
end

function start_download_watch(stems)
    watch = ColdDownloadWatch(copy(stems), time(), Threads.Atomic{Bool}(false),
        Dict{String,Float64}(), Dict{String,Float64}())
    task = Threads.@spawn begin
        while true
            now = time()
            for stem in watch.stems
                part = joinpath(COLDTEST_CACHE, stem * ".tif.part")
                final = joinpath(COLDTEST_CACHE, stem * ".tif")
                isfile(part) && get!(watch.partseen, stem, now)
                if isfile(final) && !haskey(watch.finished, stem)
                    # A sub-poll download is explicitly marked with a zero
                    # lower-bound rather than inventing a latency.
                    get!(watch.partseen, stem, now)
                    watch.finished[stem] = now
                end
            end
            partsremain = isdir(COLDTEST_CACHE) &&
                any(endswith(".part"), readdir(COLDTEST_CACHE))
            watch.stop[] && !partsremain && break
            sleep(0.005)
        end
    end
    return watch, task
end

function stop_download_watch!(watch, task)
    watch.stop[] = true
    wait(task)
    return watch
end

function coldtest_parts()
    isdir(COLDTEST_CACHE) || return String[]
    return sort!([f for f in readdir(COLDTEST_CACHE) if endswith(f, ".part")])
end

function coldtest_cachefiles()
    isdir(COLDTEST_CACHE) || return Dict{String,String}()
    return Dict(splitext(f)[1] => joinpath(COLDTEST_CACHE, f)
        for f in readdir(COLDTEST_CACHE) if endswith(f, ".tif"))
end

function coldtest_storefiles(path)
    dir = joinpath(path, "elevation")
    isdir(dir) || return Dict{String,Vector{UInt8}}()
    return Dict(f => read(joinpath(dir, f))
        for f in readdir(dir) if !startswith(f, "."))
end

function coldtest_sanity(path)
    array = Zarr.zopen(path, "r")["elevation"]
    rows = NamedTuple[]
    for column in COLDTEST_COLUMNS
        vals = Vector{Float32}(array[:, column])
        finite = filter(isfinite, vals)
        length(finite) + count(isnan, vals) == length(vals) ||
            error("column $column contains a non-finite non-NaN value")
        isempty(finite) && error("column $column has no finite elevations")
        push!(rows, (column = column, cells = length(vals), finite = length(finite),
            nan = count(isnan, vals), minimum = Float64(minimum(finite)),
            maximum = Float64(maximum(finite))))
    end
    return rows
end

function coldtest_config(store)
    return merge(CONFIG, (source = :real, real = :none, data = COLDTEST_DATA,
        tilecache = COLDTEST_CACHE, store = store, region = nothing,
        maskarcsec = 0, workers = COLDTEST_WORKERS, cores = COLDTEST_WORKERS,
        shape = :outer, schedule = :affinity, cachepolicy = :refcount,
        taper = true, prefetch = COLDTEST_PREFETCH,
        fetchconc = COLDTEST_FETCHCONC, fetchdelay = 0.0, resume = false,
        checks = false, heartbeat = 1_000_000, maxchunks = 0,
        chunks = COLDTEST_COLUMNS, dryrun = false))
end

function prepare_store(store)
    for path in (store, donelogpath(store), chunklistpath(store))
        ispath(path) && error("fresh-store requirement: $path already exists")
    end
    save_chunklist(chunklistpath(store), CONFIG.ancestor, COLDTEST_COLUMNS)
end

function appendjson(path, mode, rows)
    open(path, mode) do io
        for row in rows
            println(io, JSON3.write(row))
        end
    end
end

function run_coldtest(tag::Symbol)
    Threads.nthreads() <= 8 || error("cold test is capped at 8 Julia threads")
    Base.JLOptions().nsweepthreads == 0 || error(
        "refusing concurrent GC sweeper; use a single-field --gcthreads=N")
    isdir(COLDTEST_SCRATCH) || error("scratch root is missing: $COLDTEST_SCRATCH")

    plan = coldtest_plan()
    headbytes = coldtest_headbytes(plan)
    predictedbytes = sum(values(headbytes))
    predictedbytes < 2_000_000_000 || error(
        "predicted transfer $predictedbytes exceeds the 2 GB hard cap")
    print_coldtest_plan(plan, headbytes)

    store = tag === :cold ? COLDTEST_COLD_STORE : COLDTEST_WARM_STORE
    if tag === :cold
        ispath(COLDTEST_CACHE) && error(
            "cold cache must not exist before the run: $COLDTEST_CACHE")
    else
        files = coldtest_cachefiles()
        sort!(collect(keys(files))) == sort(plan.stems) || error(
            "warm cache does not exactly match the predicted tile set")
        isempty(coldtest_parts()) || error("warm cache contains .part files")
        isdir(COLDTEST_COLD_STORE) || error("cold output store is missing")
    end
    prepare_store(store)

    FAILURES[] = 0
    LASTCACHE[] = nothing
    LASTPROVIDER[] = nothing
    watch = task = nothing
    tag === :cold && ((watch, task) = start_download_watch(plan.stems))
    started = time()
    failures = main(coldtest_config(store))
    wall = time() - started
    tag === :cold && stop_download_watch!(watch, task)

    provider = LASTPROVIDER[]
    provider isa LazyCopernicusTiles || error("driver did not expose its real provider")
    cache = LASTCACHE[]
    cache === nothing && error("driver did not expose its cache statistics")
    actual = coldtest_cachefiles()
    actualstems = sort!(collect(keys(actual)))
    parts = coldtest_parts()
    sanity = coldtest_sanity(store)

    coldcheck("driver returned no failures", failures == 0; detail = "$failures")
    coldcheck("actual cache tile set equals predicted set",
        actualstems == sort(plan.stems); detail = "$(length(actualstems)) tiles")
    coldcheck("no .part files remain", isempty(parts); detail = join(parts, ","))
    expecteddownloads = tag === :cold ? length(plan.stems) : 0
    coldcheck("successful GET count has no duplicate tile fetches",
        provider.ndownloads[] == expecteddownloads;
        detail = "$(provider.ndownloads[]) downloads for $(length(plan.stems)) tiles")
    coldcheck("every predicted source tile was decoded exactly once",
        cache.loads == length(plan.stems);
        detail = "$(cache.loads) loads for $(length(plan.stems)) tiles")
    coldcheck("every demand was a predicted graph edge", cache.uncredited == 0;
        detail = "$(cache.uncredited)")
    coldcheck("real output range is physically plausible",
        all(row -> -500 <= row.minimum <= row.maximum <= 9_000, sanity))

    rows = Any[Dict(
        "kind" => "run", "phase" => String(tag), "columns" => COLDTEST_COLUMNS,
        "predicted_tiles" => sort(plan.stems), "predicted_bytes" => predictedbytes,
        "actual_tiles" => actualstems,
        "actual_bytes" => sum(filesize, values(actual); init = 0),
        "successful_downloads" => provider.ndownloads[],
        "demand_cold_downloads" => provider.ncold[],
        "wall_s" => wall, "cache_loads" => cache.loads,
        "cache_hits" => cache.hits, "cache_joined_loads" => cache.waits,
        "cache_uncredited" => cache.uncredited, "part_files" => parts,
        "prefetch_depth" => COLDTEST_PREFETCH,
        "fetch_concurrency" => COLDTEST_FETCHCONC,
        "workers" => COLDTEST_WORKERS, "failures" => failures,
        "sanity" => sanity)]

    if tag === :cold
        length(watch.finished) == length(plan.stems) || error(
            "watcher saw $(length(watch.finished)) of $(length(plan.stems)) completions")
        latencies = Float64[]
        for stem in sort(plan.stems)
            latency = watch.finished[stem] - watch.partseen[stem]
            push!(latencies, latency)
            push!(rows, Dict("kind" => "tile", "phase" => "cold", "stem" => stem,
                "bytes" => filesize(actual[stem]), "head_bytes" => headbytes[stem],
                "part_seen_s" => watch.partseen[stem] - started,
                "finished_s" => watch.finished[stem] - started,
                "latency_s" => latency))
        end
        rows[1]["download_latency_sum_s"] = sum(latencies)
        rows[1]["download_latency_min_s"] = minimum(latencies)
        rows[1]["download_latency_median_s"] = Statistics.median(latencies)
        rows[1]["download_latency_p90_s"] = Statistics.quantile(latencies, 0.9)
        rows[1]["download_latency_max_s"] = maximum(latencies)
        appendjson(COLDTEST_RESULTS, "w", rows)
    else
        coldchunks = coldtest_storefiles(COLDTEST_COLD_STORE)
        warmchunks = coldtest_storefiles(COLDTEST_WARM_STORE)
        samekeys = sort!(collect(keys(coldchunks))) == sort!(collect(keys(warmchunks)))
        differing = samekeys ? sort!([f for f in keys(coldchunks)
            if coldchunks[f] != warmchunks[f]]) : String[]
        identical = samekeys && isempty(differing)
        rows[1]["cold_chunk_files"] = length(coldchunks)
        rows[1]["warm_chunk_files"] = length(warmchunks)
        rows[1]["byte_identical"] = identical
        rows[1]["differing_chunk_files"] = differing
        coldcheck("warm run downloaded no tiles", provider.ndownloads[] == 0;
            detail = "$(provider.ndownloads[])")
        coldcheck("warm run had no demand-cold downloads", provider.ncold[] == 0;
            detail = "$(provider.ncold[])")
        coldcheck("cold and warm elevation chunks are byte-identical", identical;
            detail = "$(length(coldchunks)) chunk files")
        appendjson(COLDTEST_RESULTS, "a", rows)
    end

    println("RESULT phase=", tag, " wall_s=", round(wall; digits = 3),
        " downloads=", provider.ndownloads[], " demand_cold=", provider.ncold[],
        " bytes=", sum(filesize, values(actual); init = 0))
    for row in sanity
        println("OUTPUT column=", row.column, " finite=", row.finite,
            " nan=", row.nan, " range=[", row.minimum, ", ", row.maximum, "]")
    end
    println("results: ", COLDTEST_RESULTS)
    return 0
end

function main_coldtest()
    phase = isempty(ARGS) ? "plan" : ARGS[1]
    phase in ("plan", "cold", "warm") || error("phase must be plan, cold, or warm")
    if phase == "plan"
        plan = coldtest_plan()
        headbytes = coldtest_headbytes(plan)
        print_coldtest_plan(plan, headbytes)
        ispath(COLDTEST_CACHE) && error(
            "planning requires the fresh cache path to remain absent")
        println("PLAN ONLY: cache remains absent; no GeoTIFF body was downloaded")
        return 0
    end
    return run_coldtest(Symbol(phase))
end

exit(main_coldtest())
