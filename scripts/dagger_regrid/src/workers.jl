# Process-local state and complete-destination execution.

function worker_config(config::DaggerRegridConfig)
    return (; production = config.production,
        cache_slots = config.worker_cache_slots,
        cache_stripes = config.worker_cache_stripes,
        tilecache_root = config.worker_tilecache_root)
end

struct WorkerStateFactory{S}
    settings::S
end

(factory::WorkerStateFactory)() = make_worker_state(factory.settings)

function make_worker_state(settings)
    config = settings.production
    gcguard(config)
    config.malloctrim > 0 && tunemalloc(config.malloctrim)
    geometry = _source_geometry(config)
    tiledir = joinpath(config.data, "CopernicusDEM", "$(config.res)m")
    landshp = joinpath(config.data, "naturalearth", "ne_10m_land.shp")
    # Real tiles carry their own ocean nodata; only a synthetic source masks.
    mask = config.source === :synthetic ? landmask(landshp, config.maskarcsec) : NOMASK
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

function run_batch(state::WorkerState, items::Vector{WorkItem})
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

function worker_stats(state::WorkerState)
    stats = cachestats(state.dem.cache)
    builder = state.dem.builder
    return WorkerStats(Distributed.myid(), stats.loads, stats.hits, stats.live,
        stats.bytes, builder.nreal[], builder.nsynthetic[], builder.npixels[])
end
