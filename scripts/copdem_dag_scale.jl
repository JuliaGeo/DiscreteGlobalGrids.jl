# What DAG-driven scheduling costs and saves, measured two ways.
#
#     julia --project=benchmark -t 12 --gcthreads=4 scripts/copdem_dag_scale.jl OUTDIR
#
# **Part 1, the whole globe, statically.** Build the real 66 178 x 26 475
# dependency graph, take the tile-affinity walk order, and replay refcount
# eviction over it with `W` columns in flight — no regrid, no pixels, seconds of
# arithmetic. Tile sizes are the source grid's own, so the peak is in real bytes.
# This is the number to compare against the offline simulator's 1.43 GiB, and
# against the 3072-slot LRU's 16.48 GiB.
#
# **Part 2, a few hundred chunks, live.** Regrid a real region twice — once the
# way the run used to work, once DAG-driven — into scratch stores, and report
# residency, reloads, wall time and the tail. Then once more with fake fetch
# latency, with and without the prefetcher, which is the only way to see the
# prefetcher do anything at all on a synthetic source.
#
# Writes `<outdir>/scale.ndjson`. Never touches the production store.

include(joinpath(@__DIR__, "copdem_production.jl"))
import JSON3

const OUT = length(ARGS) >= 1 ? ARGS[1] : mktempdir(; prefix = "dagscale-")
const RECORDS = Dict{String,Any}[]

emit(; kw...) = (push!(RECORDS, Dict{String,Any}(String(k) => v for (k, v) in kw));
                 RECORDS[end])

# ---------------------------------------------------------------------------
# Part 1: replay the globe
# ---------------------------------------------------------------------------

"""
    replay(graph, order, bytes, W) -> (peakbytes, peaktiles, loads)

Refcount residency over `order` with `W` columns in flight at once.

The window is the point. A column's tiles must be resident from the moment it
starts until it retires, so with `W` workers the live set is the union of `W`
consecutive columns' tiles plus everything not yet released behind them — which
is what makes peak residency a property of the ORDER and not of a cache size.
Costs are ignored: every column is one step. That understates the spread a real
5x polar cost skew produces and overstates nothing, so it is a floor on the
peak, not a ceiling.
"""
function replay(sourcesof, consumerdegree, nsrc::Int, order::Vector{Int},
        bytes::Vector{Int}, W::Int)
    remaining = [consumerdegree(s) for s in 1:nsrc]
    resident = falses(nsrc)
    live = 0
    tiles = 0
    peak = 0
    peaktiles = 0
    loads = 0
    for (i, _) in enumerate(order)
        # Start column i.
        for s in sourcesof(order[i])
            si = Int(s)
            resident[si] && continue
            resident[si] = true
            live += bytes[si]
            tiles += 1
            loads += 1
        end
        peak = max(peak, live)
        peaktiles = max(peaktiles, tiles)
        # Retire the column that leaves the window.
        j = i - W + 1
        j >= 1 || continue
        for s in sourcesof(order[j])
            si = Int(s)
            remaining[si] -= 1
            remaining[si] == 0 || continue
            if resident[si]
                resident[si] = false
                live -= bytes[si]
                tiles -= 1
            end
        end
    end
    return peak, peaktiles, loads
end

"Peak bytes of an `n`-slot LRU over the same walk: the `n` largest tiles it can hold."
lrupeak(bytes::Vector{Int}, slots::Int) =
    sum(sort(bytes; rev = true)[1:min(slots, length(bytes))])

function globalreplay(config, W)
    say("="^70)
    say("part 1: the whole globe, replayed")
    sys = DGG.CopernicusDEMSystem(config.res)
    sys7 = DGG.IGeo7System()
    tiles = listedtiles(sys, joinpath(config.data, "CopernicusDEM",
        "tileList-glo$(config.res).txt"), nothing)
    ids = TileIds(sys, tiles)
    srcspace = DGG.DGGSpace(DGG.PartialGrid(sys, 1, ids); chunklevel = 0)
    chunks = load_chunklist(chunklistpath(config.store))
    chunks === nothing && error("no cached covering at $(chunklistpath(config.store))")
    plan = dagplan(sys, sys7, tiles, chunks, srcspace, config)
    # A tile's real size: its pixel count from the source grid, four bytes each.
    bytes = [4 * ((k < length(tiles) ? ids.offsets[k + 1] : ids.n) - ids.offsets[k])
             for k in 1:length(tiles)]
    say(@sprintf("graph: %d x %d, %d edges; all tiles at once %.2f GiB",
        length(tiles), length(chunks), plan.edges, sum(bytes) / 2^30))
    canonical = collect(1:length(chunks))
    sourcesof = d -> GR.sourcesof(plan.graph, d)
    consumerdegree = s -> GR.consumerdegree(plan.graph, s)
    for (name, ord) in ("affinity" => plan.order, "canonical" => canonical)
        t = @elapsed pk, pt, ld = replay(sourcesof, consumerdegree, length(tiles),
            ord, bytes, W)
        say(@sprintf("%-10s refcount peak %.2f GiB / %d tiles, %d loads (%d tiles), replay %s",
            name, pk / 2^30, pt, ld, length(tiles), secs(t)))
        emit(part = "globalreplay", order = name, workers = W,
            peakbytes = pk, peaktiles = pt, loads = ld, tiles = length(tiles),
            chunks = length(chunks), edges = plan.edges, allbytes = sum(bytes))
    end
    for slots in (3072, 1024)
        say(@sprintf("LRU %d slots would hold %.2f GiB", slots,
            lrupeak(bytes, slots) / 2^30))
        emit(part = "lrubound", slots = slots, peakbytes = lrupeak(bytes, slots))
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Part 2: a live region
# ---------------------------------------------------------------------------

"Worker idle at the tail: makespan minus the first worker to run out of work."
function tailidle(donelog, wall)
    isfile(donelog) || return NaN
    last = Dict{Int,Float64}()
    t0 = Inf
    for line in eachline(donelog)
        j = try JSON3.read(line) catch; continue end
        t = Dates.datetime2unix(Dates.DateTime(j.t))
        t0 = min(t0, t - j.secs)
        last[j.w] = max(get(last, j.w, 0.0), t)
    end
    isempty(last) && return NaN
    return (t0 + wall) - minimum(values(last))
end

function liveside(base, tag, extra)
    store = joinpath(OUT, tag * ".zarr")
    for f in (store, store * ".done.ndjson", store * ".columns.txt")
        rm(f; recursive = true, force = true)
    end
    FAILURES[] = 0
    LASTCACHE[] = nothing
    config = merge(base, (store = store, resume = false, heartbeat = 1_000_000,
        checks = false), extra)
    say("="^70)
    say("part 2: side `$tag`")
    t0 = time()
    fails = main(config)
    wall = time() - t0
    st = LASTCACHE[]
    rec = emit(part = "live", side = tag, wall = wall, failures = fails,
        tail = tailidle(store * ".done.ndjson", wall),
        loads = st === nothing ? -1 : st.loads,
        peakbytes = st === nothing ? -1 : st.peakbytes,
        peaktiles = st === nothing ? -1 : st.peaktiles,
        liveend = st === nothing ? -1 : st.live,
        uncredited = st === nothing ? -1 : st.uncredited,
        schedule = String(config.schedule), cachepolicy = String(config.cachepolicy),
        taper = config.taper, prefetch = config.prefetch,
        fetchdelay = config.fetchdelay, batch = config.batch)
    say(@sprintf("side `%s`: wall %s, tail idle %.1f s, %d failure(s)",
        tag, secs(wall), rec["tail"], fails))
    return rec
end

function main_scale()
    W = 12
    base = merge(CONFIG, (cores = W, workers = W, dryrun = false,
        region = [(0.0, 20.0, 35.0, 50.0)], maxchunks = 400))
    globalreplay(merge(CONFIG, (region = nothing,)), 23)

    liveside(base, "legacy", (schedule = :canonical, cachepolicy = :lru,
        taper = false, prefetch = 0, batch = 8))
    liveside(base, "dag", (schedule = :affinity, cachepolicy = :refcount,
        taper = true, prefetch = 0, batch = 8))
    # The prefetcher only has anything to hide when a tile costs time to obtain,
    # and a fabricated tile costs none, so impose some.
    liveside(base, "delay-nopf", (schedule = :affinity, cachepolicy = :refcount,
        taper = true, prefetch = 0, batch = 8, fetchdelay = 0.5, maxchunks = 150))
    liveside(base, "delay-pf", (schedule = :affinity, cachepolicy = :refcount,
        taper = true, prefetch = 32, fetchconc = 8, batch = 8, fetchdelay = 0.5,
        maxchunks = 150))

    path = joinpath(OUT, "scale.ndjson")
    open(path, "w") do io
        for r in RECORDS
            println(io, JSON3.write(r))
        end
    end
    say("wrote $path")
    return 0
end

exit(main_scale())
