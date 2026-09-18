# Shared by both CopDEM drivers: the configuration, the source, the dependency
# graph with its walk order, and the one-chunk regrid.
#
# An include fragment. `copdem_production.jl` (threads) and the `DaggerRegrid`
# module (processes) each load it into their own module. Tile access, the lazy
# DEM array and the synthetic source live in `lib/CopernicusUtils`.

import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import GeometryOps as GO
import DimensionalData as DD
import Zarr
import Statistics
import Dates
import Printf: @sprintf
using Base.ScopedValues: @with
using CopernicusUtils

const CD = DGG.CopernicusDEM
const US = GO.UnitSpherical

_env(key, default::Int) = parse(Int, get(ENV, key, string(default)))
_env(key, default::Bool) = get(ENV, key, default ? "1" : "0") != "0"
_env(key, default::Symbol) = Symbol(get(ENV, key, String(default)))

"`COPDEM_REGION=w,e,s,n` as a one-box region list; `nothing` when unset."
function _envregion()
    haskey(ENV, "COPDEM_REGION") || return nothing
    w, e, s, n = parse.(Float64, split(ENV["COPDEM_REGION"], ','))
    return [(w, e, s, n)]
end

"""
    copdem_config(; overrides...) -> NamedTuple

The run configuration: defaults, then `COPDEM_*` environment variables, then
`overrides`. A driver's `main` takes the result, so a script or REPL session
runs a variant with `main(copdem_config(source = :real, maxchunks = 2))`.
"""
function copdem_config(; overrides...)
    res = _env("COPDEM_RES", 90)                # 90 for GLO-90, 30 for GLO-30
    level = _env("COPDEM_LEVEL", 12)            # IGeo7 output level
    source = _env("COPDEM_SOURCE", :synthetic)  # :real (AWS COGs) or :synthetic
    data = get(ENV, "RASTERDATASOURCES_PATH", joinpath(@__DIR__, "..", "benchmark", "data"))
    config = (;
        res, level, source, data,
        ancestor    = _env("COPDEM_ANCESTOR", 5),  # chunk root level
        authalic    = true,     # compute on the WGS84 authalic geometry
        method      = _env("COPDEM_METHOD", :conservative), # or :point, :nearest, :nearest-direct
        store       = get(ENV, "COPDEM_STORE",
                          joinpath(data, "stores", "glo$res-$source-igeo7-l$level.zarr")),
        region      = _envregion(),  # nothing for the globe, or [(w, e, s, n), ...]
        maskarcsec  = _env("COPDEM_MASKARCSEC", 15), # synthetic land mask; 0 disables it
        tilecache   = get(ENV, "COPDEM_TILE_CACHE", joinpath(data, "CopernicusDEM", "tiles")),
        tilebaseurl = get(ENV, "COPDEM_TILE_BASEURL", CopernicusUtils.bucketurl(res)),
        download    = _env("COPDEM_DOWNLOAD", true), # false: tilecache is complete and read-only
        workers     = _env("COPDEM_WORKERS", 2 * Threads.nthreads()), # concurrent chunk tasks
        batch       = 8,        # chunks handed out per pull, at most
        budget      = 2^30,     # lazy-regrid byte budget, per worker
        prefetch    = _env("COPDEM_PREFETCH", 0), # chunks of tile lookahead; 0 disables
        resume      = true,     # skip chunks already written
        write       = _env("COPDEM_WRITE", true), # false: compute and ledger, never write
        checks      = _env("COPDEM_CHECKS", false), # run the synthetic oracle afterwards
        checkchunks = 6,        # chunks the oracle reads back
        heartbeat   = 300,      # seconds between progress lines
        maxchunks   = _env("COPDEM_MAXCHUNKS", 0), # 0 = no limit
        chunks      = Int[],    # explicit chunk indices, overriding the covering
        dryrun      = _env("COPDEM_DRYRUN", false), # plan and report, compute nothing
    )
    config = merge(config, NamedTuple(overrides))
    config.source in (:real, :synthetic) ||
        error("source must be :real or :synthetic, got $(repr(config.source))")
    regridmethod(config)
    return config
end

"""
    regridmethod(config) -> AbstractRegriddingMethod

| `config.method` | method | a destination cell holds |
|---|---|---|
| `:conservative` | `Conservative()` | the coverage-normalised mean of the pixels it overlaps |
| `:point` | `BarycentricPoint()` | a sample at its centroid, interpolated between posts |
| `:nearest` | `NearestCell()` | the nearest post |
| `:nearest-direct` | `DirectNearest()` | the nearest post, with no weight assembly |

Every plan in a run takes its method from here. The global plan's dependency
relation is built at the method's own `supportradius`, so the per-chunk plans
and the global one must agree.
"""
function regridmethod(config)
    config.method === :conservative && return DGG.Conservative()
    config.method === :point && return DGG.BarycentricPoint()
    config.method === :nearest && return DGG.NearestCell()
    config.method in (Symbol("nearest-direct"), :nearest_direct) && return DGG.DirectNearest()
    return error("method must be :conservative, :point, :nearest or :nearest-direct, " *
                 "got $(repr(config.method))")
end

"""
    regridpolicy(config) -> AbstractMissingPolicy

`Weighted(1)` for `:point`: a stencil naming a post with no elevation blanks the
cell. `Weighted(0.5)` otherwise: an area row blanks below half coverage, and a
nearest row is one weight of 1.0 and passes through.
"""
regridpolicy(config) = config.method === :point ? DGG.Weighted(1) : DGG.Weighted(0.5)

secs(t) = @sprintf("%.1f s", t)
hours(t) = @sprintf("%.2f h", t / 3600)

"""
    tunemalloc() -> Bool

Freeze glibc's malloc trim threshold at 32 MiB; `false` off glibc.

glibc grows its mmap and trim thresholds as large blocks are freed, after which
decoded tiles and weight blocks (both ~5 MiB) come from arena heap and are never
returned: resident memory settles near three times the live heap. Setting any
threshold explicitly stops the growth.
"""
function tunemalloc()
    Sys.islinux() || return false
    return try
        # M_TRIM_THRESHOLD is -1 in glibc's malloc.h; mallopt returns 1 on success.
        ccall(:mallopt, Cint, (Cint, Cint), Cint(-1), Cint(32 * 2^20)) == 1
    catch
        false
    end
end

"""
    gcguard()

Refuse to run under `--gcthreads=N,1`. The concurrent page sweeper `madvise`s
freed pages from a background thread while workers run, and this workload has
segfaulted under it on a live, task-local array. Pass `--gcthreads=N`.
"""
function gcguard()
    nsweep = Base.JLOptions().nsweepthreads
    nsweep == 0 || error("Julia's concurrent page sweeper is on " *
        "(--gcthreads=…,$nsweep) and segfaults this run; pass " *
        "--gcthreads=$(Base.JLOptions().nmarkthreads)")
    return nothing
end

const LOGLOCK = ReentrantLock()
# Atomic: workers report their own chunk failures.
const FAILURES = Threads.Atomic{Int}(0)

stamp() = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")

function say(parts...)
    lock(LOGLOCK) do
        println(stamp(), "  ", parts...)
        flush(stdout)
    end
end

"Report a claim and count the failures."
function check(name, ok; detail = "")
    ok || Threads.atomic_add!(FAILURES, 1)
    say(ok ? "PASS  " : "FAIL  ", rpad(name, 56), detail)
    return ok
end

include("copdem_store.jl")
include("copdem_policy.jl")

# ===========================================================================
# The source
# ===========================================================================

"""
    opensource(config) -> NamedTuple

Everything about the source that reads no pixel: both grid systems, the listed
tiles, the tile source, and the source space whose chunk `k` is `tiles[k]`.
`mask` is the synthetic land mask, [`NOMASK`](@ref) for real tiles.
"""
function opensource(config)
    sys = DGG.CopernicusDEMSystem(config.res)
    storesys7 = DGG.IGeo7System()
    sys7 = config.authalic ? DGG.AuthalicSystem(storesys7) : storesys7
    list = tilelist(config.data, config.res; baseurl = config.tilebaseurl)
    tiles = listedtiles(sys, list, config.region)
    isempty(tiles) && error("the tile list and region select no tiles")
    mask = NOMASK
    tilesource = if config.source === :real
        CopernicusTiles(sys, tiles; cachedir = config.tilecache,
            download = config.download, baseurl = config.tilebaseurl)
    else
        shp = joinpath(config.data, "naturalearth", "ne_10m_land.shp")
        mask = landmask(shp, config.maskarcsec)
        SyntheticTiles(sys, mask)
    end
    ids = TileIds(sys, tiles)
    srcspace = DGG.DGGSpace(DGG.PartialGrid(sys, 1, ids); chunklevel = 0)
    return (; sys, storesys7, sys7, tiles, tilesource, mask, ids, srcspace)
end

"""
    destination_chunks(config, src) -> Vector{Int}

The level-`ancestor` chunks this run covers: `config.chunks` when given,
otherwise the covering of the listed tiles, cached beside the store, thinned to
`config.maxchunks` evenly spaced chunks when that is set.
"""
function destination_chunks(config, src)
    isempty(config.chunks) || return sort!(unique(config.chunks))
    path = chunklistpath(config.store)
    chunks = load_chunklist(path)
    if chunks === nothing
        chunks = covering_chunks(src.sys7, src.sys, src.tiles, config.ancestor;
            nthreads = max(1, Threads.nthreads() - 1))
        save_chunklist(path, config.ancestor, chunks)
    end
    if 0 < config.maxchunks < length(chunks)
        chunks = chunks[round.(Int, range(1, length(chunks); length = config.maxchunks))]
    end
    return chunks
end

# ===========================================================================
# The dependency graph, and the walk order it implies
# ===========================================================================

"""
    dagplan(src, chunks, config) -> (; graph, order)

The tile-to-chunk dependency graph over exactly the work this run does, and the
order to walk it in.

Graph destination `d` is `chunks[d]` and graph source `s` is `src.tiles[s]`, the
numbering the source space and the tile caches use. `order[p]` is the
destination to run at position `p`; see [`affinity_order`](@ref).

A `GR.ChunkedPlan` over the whole covering owns the relation, so the schedule
and the refcount cache read the same pairing the executor uses. Building it
reads no DEM data.
"""
function dagplan(src, chunks::Vector{Int}, config)
    g = DGG.levelgrid(src.sys7, config.ancestor)
    dstids = SubtreeIds(src.sys7, [DGG.cellindex(g, c) for c in chunks], config.level)
    dstspace = DGG.DGGSpace(DGG.PartialGrid(src.sys7, config.level, dstids);
        chunklevel = config.ancestor)
    plan = GR.ChunkedPlan(regridmethod(config), regridpolicy(config), dstspace,
        src.srcspace; budget = config.budget, dependencies = true)
    graph = GR.dependencies(plan)
    # Morton sweep of the 1-degree tile lattice.
    keys = map(src.tiles) do t
        lat, lon = CD.tilecorner(src.sys, DGG.LevelIndex(0, t))
        morton2(Int(lon) + 180, Int(lat) + 90)
    end
    order = affinity_order(length(chunks), d -> GR.sourcesof(graph, d), keys)
    return (; graph, order)
end

"""
    regrid_chunk(dem, srcspace, sys7, layout, chunk, config) -> Vector{Float32}

One work unit: the level-`config.level` values of one chunk.

`DGG.subtree` hands the destination over as a rooted grid, so the regridder
knows it is one chunk without scanning the level-`ancestor` grid. The plan
builds its own one-row dependency relation, which costs a fraction of a
millisecond against a chunk that takes seconds.

Pass no `chunks` to `GR.regrid`: supplying it defeats the plan's pairing, which
is what prunes the source tiles this chunk does not meet.
"""
function regrid_chunk(dem, srcspace, sys7, layout, chunk::Int, config)
    dstgrid = DGG.subtree(sys7, DGG.columncell(layout, chunk), config.level)
    dstspace = DGG.DGGSpace(dstgrid; chunklevel = config.ancestor)
    out = GR.regrid(dem; to = dstspace, from = srcspace,
        method = regridmethod(config), missingpolicy = regridpolicy(config),
        lazy = true, budget = config.budget)
    return Float32.(vec(collect(out)))
end

include("copdem_verify.jl")
