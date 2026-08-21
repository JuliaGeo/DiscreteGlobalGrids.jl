# Copernicus DEM GLO-90 -> IGEO7 level 12, into one global ancestor-subzone Zarr
# store, written by W worker tasks in one process.
#
# Edit `CONFIG` below, then run it:
#
#     julia --project=bench -t 26 --gcthreads=8,1 scripts/copdem_production.jl
#
# The run is resumable: it skips every chunk the ledger or the store already has,
# so re-running the same command continues where it stopped.
#
# Words used throughout:
#
#   tile    a 1x1-degree Copernicus DEM tile. The source unit: one source chunk.
#   chunk   one level-`ancestor` IGeo7 cell with all its level-`level`
#           descendants (7^7 = 823 543 cells at 5 -> 12). The destination work
#           unit, one Zarr chunk, one file. The store API calls it a "column".
#   budget  bytes of intersection weights a lazy regrid may hold at once, per
#           worker.
#
# Elevations are FABRICATED — see `copdem_synthetic.jl`. What is real is the tile
# LIST: Copernicus ships ~26 450 of the 64 800 tiles the 1x1-degree lattice has,
# and the source is a `PartialGrid` over exactly those, so a destination chunk
# over open ocean pairs with no source at all. Store I/O and resume live in
# `copdem_store.jl`.

import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import GeometryOps as GO
import DimensionalData as DD
import DiskArrays
import ArchGDAL
import GeoInterface as GI
import Extents
import Zarr
import Statistics
import Dates
import Printf: @sprintf
using Base.ScopedValues: @with

const CD = DGG.CopernicusDEM
const US = GO.UnitSpherical

# ===========================================================================
# Configuration — edit this, there is no command line
# ===========================================================================

# These are the values the production run is using.
const CONFIG = (
    res         = 90,       # 90 for GLO-90, 30 for GLO-30
    level       = 12,       # IGeo7 output level
    ancestor    = 5,        # chunk root level
    store       = "/home/asinghvi17/geo/dggstores/copdem90-igeo7-l12-synthetic.zarr",
    region      = nothing,  # nothing for the globe, or [(w, e, s, n), ...] boxes
    maskarcsec  = 15,       # land-mask lattice, arcseconds; 0 disables the mask
    real        = :auto,    # :auto (every GeoTIFF found), :none, or ["stem", ...]
    workers     = 21,       # concurrent worker tasks
    batch       = 8,        # chunks handed out per pull
    budget      = 2^30,     # lazy-regrid byte budget, per worker
    cache       = 3072,     # decoded source tiles held across all workers
    stripes     = 64,       # independent locks over that cache
    resume      = true,     # skip chunks already written
    checks      = false,    # run the synthetic oracle after the run
    checkchunks = 6,        # chunks to verify when `checks`
    heartbeat   = 300,      # seconds between summary lines
    maxchunks   = 0,        # 0 = no limit; a smoke-test knob
    chunks      = Int[],    # explicit chunk indices, overriding the covering
    dryrun      = false,    # plan and report, compute nothing
    data        = get(ENV, "RASTERDATASOURCES_PATH",
                      joinpath(@__DIR__, "..", "bench", "data")),
)

# ===========================================================================
# Logging
# ===========================================================================

const LOGLOCK = ReentrantLock()
const STARTED = Ref(time())
const FAILURES = Ref(0)

stamp() = Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS")

function say(parts...)
    lock(LOGLOCK) do
        println(stamp(), "  ", parts...)
        flush(stdout)
    end
end

"Report a claim and count the failures."
function check(name, ok; detail = "")
    lock(LOGLOCK) do
        ok || (FAILURES[] += 1)
        println(stamp(), "  ", ok ? "PASS  " : "FAIL  ", rpad(name, 56), detail)
        flush(stdout)
    end
    return ok
end

rssgib() = Sys.maxrss() / 2^30
secs(t) = @sprintf("%.1f s", t)
hours(t) = @sprintf("%.2f h", t / 3600)

include("copdem_store.jl")
include("copdem_synthetic.jl")

# ===========================================================================
# The tile list: which tiles exist at all
# ===========================================================================

"The tile an AWS stem names, e.g. `Copernicus_DSM_COG_30_S90_00_E000_00_DEM`."
function stemtile(sys, stem::AbstractString)
    m = match(r"_([NS])(\d{2})_00_([EW])(\d{3})_00_DEM$", stem)
    m === nothing && throw(ArgumentError("$stem is not a Copernicus DEM stem"))
    lat = parse(Int, m[2]) * (m[1] == "S" ? -1 : 1)
    lon = parse(Int, m[4]) * (m[3] == "W" ? -1 : 1)
    return CD.tilecell(sys, lat, lon)
end

"The AWS object stem of `tile`. Inverse of [`stemtile`](@ref)."
function tilestem(sys, tile)
    lat, lon = CD.tilecorner(sys, tile)
    tag = CD.lat_intervals(sys) == 3600 ? "10" : "30"
    return string("Copernicus_DSM_COG_", tag, "_", lat < 0 ? "S" : "N",
        lpad(abs(lat), 2, '0'), "_00_", lon < 0 ? "W" : "E",
        lpad(abs(lon), 3, '0'), "_00_DEM")
end

"True unless `regions` is a list of (w, e, s, n) boxes and (`lon`, `lat`) misses them all."
inregions(::Nothing, lon, lat) = true
inregions(regions, lon, lat) =
    any(r -> r[1] <= lon <= r[2] && r[3] <= lat <= r[4], regions)

"""
    listedtiles(sys, path, regions) -> Vector{Int}

The ordinals of every tile named in the Copernicus tile list, ascending, kept to
`regions` if those are given.

A tile NOT on the list does not exist here: it is absent from the source grid, so
a destination cell over open ocean never pairs with a source chunk and is never
computed.
"""
function listedtiles(sys, path::AbstractString, regions)
    isfile(path) || throw(ArgumentError("no tile list at $path"))
    out = Int[]
    for line in eachline(path)
        stem = strip(line)
        isempty(stem) && continue
        t = stemtile(sys, stem)
        lat, lon = CD.tilecorner(sys, t)
        inregions(regions, lon, lat) || continue
        push!(out, Int(t.index))
    end
    return sort!(unique!(out))
end

"""
    realtiles(sys, dir, spec) -> Dict{Int,String}

The tile ordinals backed by a GeoTIFF in `dir`: `:auto` takes every file found
there, `:none` gives an all-synthetic globe, a vector of stems takes just those.
A named stem that is not on disk is an error.
"""
function realtiles(sys, dir, spec)
    spec === :none && return Dict{Int,String}()
    isdir(dir) || (spec === :auto && return Dict{Int,String}())
    stems = spec === :auto ?
            [splitext(f)[1] for f in readdir(dir) if endswith(f, ".tif")] :
            String.(spec)
    out = Dict{Int,String}()
    for stem in stems
        path = joinpath(dir, stem * ".tif")
        isfile(path) || throw(ArgumentError("no GeoTIFF for $stem at $path"))
        out[Int(stemtile(sys, stem).index)] = path
    end
    return out
end

# ===========================================================================
# The source: one lazy vector of pixels, one chunk per listed tile
# ===========================================================================

"""
    TileIds(sys, tiles)

The pixel (level-1) cell ids of the listed tiles as one lazy ascending vector.

Lazy because the source grid stores its ids by reference and never materialises
them, so 3.8e10 pixels cost an offsets table of 26 450 integers. It also declares
itself sorted, which saves the O(n) scan the grid would otherwise run: the tiles
are sorted and their descendant ranges are disjoint, so the concatenation ascends
by construction.
"""
struct TileIds{G} <: AbstractVector{DGG.LevelIndex}
    complete::G
    starts::Vector{Int}      # first level-1 position of each tile
    offsets::Vector{Int}     # cells before each tile; `offsets[1] == 0`
    n::Int
end

function TileIds(sys, tiles::Vector{Int})
    complete = DGG.levelgrid(sys, 1)
    starts = Vector{Int}(undef, length(tiles))
    offsets = Vector{Int}(undef, length(tiles))
    acc = 0
    for (k, t) in enumerate(tiles)
        r = DGG.descendant_range(sys, DGG.LevelIndex(0, t), 1)
        starts[k] = Int(first(r))
        offsets[k] = acc
        acc += length(r)
    end
    return TileIds(complete, starts, offsets, acc)
end

Base.size(v::TileIds) = (v.n,)
Base.IndexStyle(::Type{<:TileIds}) = IndexLinear()

@inline function Base.getindex(v::TileIds, i::Int)
    @boundscheck (1 <= i <= v.n) || throw(BoundsError(v, i))
    k = searchsortedlast(v.offsets, i - 1)
    return DGG.cellindex(v.complete, v.starts[k] + (i - 1 - v.offsets[k]))
end

DGG.Helpers.strictly_increasing(::TileIds) = true

"Which listed tile holds source position `p`, and the offset within it."
@inline function tileat(v::TileIds, p::Int)
    k = searchsortedlast(v.offsets, p - 1)
    return k, p - v.offsets[k]
end

"""
    TiledDEM(sys, ids, tiles; realtiles, mask, cachesize, stripes)

Every listed tile's pixels as one `Float32` vector in the source grid's own
position order, chunk `k` being listed tile `k` — so a read is always tile
aligned and the regridder's source chunks are the DEM's own tiles.

A tile named in `realtiles` decodes from its GeoTIFF and brings its own nodata.
Every other tile is fabricated by `synthetic_tile` (see `copdem_synthetic.jl`).
The two are indistinguishable downstream: `NaN` is the regridder's invalid
sentinel whatever produced it.

The cache is striped — `stripes` independent LRUs under `stripes` locks — and a
tile is built OUTSIDE its lock, so the workers do not serialise on one mutex for
the ~15 ms a 1200x1200 tile takes.
"""
struct TiledDEM{S<:DGG.CopernicusDEMSystem,I,C,M} <: DiskArrays.AbstractDiskArray{Float32,1}
    sys::S
    ids::I
    tiles::Vector{Int}
    chunks::C
    realtiles::Dict{Int,String}
    mask::M
    caches::Vector{Dict{Int,Vector{Float32}}}
    orders::Vector{Vector{Int}}
    locks::Vector{ReentrantLock}
    per::Int
    nreal::Threads.Atomic{Int}
    nsynthetic::Threads.Atomic{Int}
    nland::Threads.Atomic{Int}
    npixels::Threads.Atomic{Int}
end

function TiledDEM(sys, ids::TileIds, tiles::Vector{Int}; realtiles::Dict{Int,String},
    mask, cachesize::Integer = 3072, stripes::Integer = 64)
    widths = [ids.offsets[k + 1] - ids.offsets[k] for k in 1:(length(tiles) - 1)]
    push!(widths, ids.n - ids.offsets[end])
    chunks = DiskArrays.GridChunks(DiskArrays.IrregularChunks(; chunksizes = widths))
    ns = Int(stripes)
    return TiledDEM(sys, ids, tiles, chunks, realtiles, mask,
        [Dict{Int,Vector{Float32}}() for _ in 1:ns], [Int[] for _ in 1:ns],
        [ReentrantLock() for _ in 1:ns], max(1, Int(cachesize) ÷ ns),
        Threads.Atomic{Int}(0), Threads.Atomic{Int}(0),
        Threads.Atomic{Int}(0), Threads.Atomic{Int}(0))
end

Base.size(A::TiledDEM) = (A.ids.n,)
DiskArrays.eachchunk(A::TiledDEM) = A.chunks
DiskArrays.haschunks(::TiledDEM) = DiskArrays.Chunked()

"The decoded band of the GeoTIFF at `path`, flattened into position order."
readtile(path) = vec(ArchGDAL.read(ds -> ArchGDAL.read(ds, 1), path))

"One tile's pixels, from the stripe cache or freshly built."
function tilevalues!(A::TiledDEM, ordinal::Int)
    s = mod(ordinal, length(A.locks)) + 1
    lk, cache, order = A.locks[s], A.caches[s], A.orders[s]
    hit = lock(lk) do
        v = get(cache, ordinal, nothing)
        v === nothing && return nothing
        push!(order, splice!(order, findfirst(==(ordinal), order)))
        return v
    end
    hit === nothing || return hit
    # Built outside the lock: two tasks racing the same tile both build it and
    # the loser's copy is dropped, which is cheaper than blocking the pool.
    path = get(A.realtiles, ordinal, nothing)
    v = if path === nothing
        Threads.atomic_add!(A.nsynthetic, 1)
        vals, nland = synthetic_tile(A.sys, DGG.LevelIndex(0, ordinal), A.mask)
        Threads.atomic_add!(A.nland, nland)
        Threads.atomic_add!(A.npixels, length(vals))
        vals
    else
        Threads.atomic_add!(A.nreal, 1)
        vals = readtile(path)
        Threads.atomic_add!(A.nland, count(isfinite, vals))
        Threads.atomic_add!(A.npixels, length(vals))
        vals
    end
    return lock(lk) do
        existing = get(cache, ordinal, nothing)
        existing === nothing || return existing
        length(cache) >= A.per && delete!(cache, popfirst!(order))
        cache[ordinal] = v
        push!(order, ordinal)
        return v
    end
end

function DiskArrays.readblock!(A::TiledDEM, out, r::AbstractUnitRange)
    p = first(r)
    while p <= last(r)
        k, off = tileat(A.ids, p)
        width = (k < length(A.tiles) ? A.ids.offsets[k + 1] : A.ids.n) - A.ids.offsets[k]
        stop = min(last(r), A.ids.offsets[k] + width)
        seg = (p - first(r) + 1):(stop - first(r) + 1)
        v = tilevalues!(A, A.tiles[k])
        out[seg] .= view(v, off:(off + length(seg) - 1))
        p = stop + 1
    end
    return out
end

# ===========================================================================
# The destination chunks
# ===========================================================================

"""
    covering_chunks(sys7, sys, tiles, ancestor; nthreads) -> Vector{Int}

The level-`ancestor` chunks that cover the listed tiles, ascending.

Each tile's 1x1-degree extent is queried AT `ancestor`, so the answer is whole
level-`ancestor` cells and a work unit's destination is a complete subtree — what
the subzone store requires and what makes a chunk exactly one Zarr chunk.

Every chunk here meets at least one listed tile, so pruning empty work is
structural: a chunk with no source is never enqueued.
"""
function covering_chunks(sys7, sys, tiles::Vector{Int}, ancestor::Int; nthreads = 1)
    g = DGG.levelgrid(sys7, ancestor)
    parts = [Set{Int}() for _ in 1:nthreads]
    Threads.@sync for w in 1:nthreads
        Threads.@spawn for k in w:nthreads:length(tiles)
            t = DGG.LevelIndex(0, tiles[k])
            lat, lon = CD.tilecorner(sys, t)
            ex = Extents.Extent(X = (Float64(lon), Float64(lon) + 1),
                Y = (Float64(lat), Float64(lat) + 1))
            set = DGG.query(sys7, DGG.MultiOrderCoverage(ex); level = ancestor)
            for c in DGG.CellVector(set; level = ancestor)
                push!(parts[w], DGG.cellposition(g, c))
            end
        end
    end
    return sort!(collect(union(parts...)))
end

# ===========================================================================
# The run
# ===========================================================================

mutable struct Progress
    const lock::ReentrantLock
    done::Int
    skipped::Int
    cells::Int
    nan::Int
    const total::Int
    const started::Float64
    const every::Float64     # seconds between heartbeats
    last::Float64
end

Progress(total, every) =
    Progress(ReentrantLock(), 0, 0, 0, 0, total, time(), Float64(every), time())

"""
    regrid_chunk(dem, srcspace, sys7, layout, chunk, config) -> Vector{Float32}

One work unit: the level-`config.level` values of one chunk.

`DGG.subtree` hands the destination over as a ROOTED grid, so the regridder knows
it is one chunk without scanning the global level-`ancestor` grid for it. The
regrid is lazy and `chunks` is deliberately not passed — supplying it defeats the
plan's pairing, which is what prunes the source tiles this chunk does not meet.
"""
function regrid_chunk(dem, srcspace, sys7, layout, chunk::Int, config)
    a = DGG.columncell(layout, chunk)
    dstgrid = DGG.subtree(sys7, a, config.level)
    dstspace = DGG.DGGSpace(dstgrid; chunklevel = config.ancestor)
    out = GR.regrid(dem; to = dstspace, from = srcspace,
        method = DGG.Conservative(), missingpolicy = DGG.Weighted(0.5),
        lazy = true, budget = config.budget)
    return Float32.(vec(collect(out)))
end

function heartbeat!(p::Progress, force = false)
    lock(p.lock) do
        now = time()
        (force || now - p.last >= p.every) || return nothing
        p.last = now
        el = now - p.started
        done = p.done + p.skipped
        rate = p.cells / max(el, 1e-9)
        left = p.total - done
        eta = done > 0 ? el * left / done : NaN
        say(@sprintf("HEARTBEAT  %d/%d chunks (%d skipped) | %.3e cells | %.0f cells/s | elapsed %s | ETA %s | RSS %.1f GiB | NaN %.2f%%",
            done, p.total, p.skipped, Float64(p.cells), rate, hours(el), hours(eta),
            rssgib(), 100 * p.nan / max(p.cells, 1)))
        return nothing
    end
end

"""
    runchunks(chunks, dem, srcspace, sys7, store, layout, done, donelog, config) -> Progress

`config.workers` tasks pulling contiguous batches of chunks off one atomic
cursor, each computing its chunk and writing it straight to its own Zarr file.
Chunks are disjoint, so no two tasks ever touch the same file.

Batches are contiguous in canonical Z7 order, so a worker walks a geographically
connected run and its tile cache stays hot; pulling dynamically keeps the polar
chunks — which meet 200-360 tiles and cost ~5x a mid-latitude one — off the
critical path.

`GR.OUTER_PARALLEL` is set inside each worker body: the outer loop is already one
task per chunk, so the nested weight builds must not spawn their own and
oversubscribe the box. The lazy plan's own source wave still spawns, which is the
parallelism a polar chunk can actually use.
"""
function runchunks(chunks, dem, srcspace, sys7, store, layout, done, donelog, config)
    W, B = config.workers, config.batch
    W >= 1 || error("workers must be at least 1, got $W")
    p = Progress(length(chunks), config.heartbeat)
    cursor = Threads.Atomic{Int}(1)
    log = DoneLog(donelog)
    try
        Threads.@sync for w in 1:W
            Threads.@spawn @with GR.OUTER_PARALLEL => true begin
                while true
                    i = Threads.atomic_add!(cursor, B)
                    i > length(chunks) && break
                    for ch in chunks[i:min(i + B - 1, length(chunks))]
                        if ch in done
                            lock(p.lock) do; p.skipped += 1; end
                            continue
                        end
                        t0 = time()
                        vals = try
                            regrid_chunk(dem, srcspace, sys7, layout, ch, config)
                        catch err
                            say("ERROR worker $w chunk $ch: ",
                                sprint(showerror, err, catch_backtrace())[1:min(end, 1500)])
                            FAILURES[] += 1
                            continue
                        end
                        DGG.dggwrite!(store, ch, vals)
                        el = time() - t0
                        nnan = count(isnan, vals)
                        lock(p.lock) do
                            p.done += 1
                            p.cells += length(vals)
                            p.nan += nnan
                        end
                        record!(log, ch, length(vals), nnan, el, w)
                        say(@sprintf("w%02d chunk %d  %d cells  %.1f s  %.0f cells/s  %.1f%% NaN",
                            w, ch, length(vals), el, length(vals) / el,
                            100 * nnan / length(vals)))
                        heartbeat!(p)
                    end
                end
            end
        end
    finally
        close(log)
    end
    heartbeat!(p, true)
    return p
end

# ===========================================================================
# Main
# ===========================================================================

function main(config = CONFIG)
    println("="^92)
    println(stamp(), "  copdem_production.jl — GLO-$(config.res) -> IGEO7 level " *
                     "$(config.level), level-$(config.ancestor) chunks, SYNTHETIC elevations")
    println("  julia $(VERSION)  threads=$(Threads.nthreads())  gc=$(Threads.ngcthreads())  pid=$(getpid())")
    println("  ", config)
    println("="^92)
    flush(stdout)

    sys = DGG.CopernicusDEMSystem(config.res)
    sys7 = DGG.IGeo7System()
    capacity = 7^(config.level - config.ancestor)

    tiledir = joinpath(config.data, "CopernicusDEM", "$(config.res)m")
    tilelist = joinpath(config.data, "CopernicusDEM", "tileList-glo$(config.res).txt")
    landshp = joinpath(config.data, "naturalearth", "ne_10m_land.shp")
    donelog = donelogpath(config.store)
    chunklist = chunklistpath(config.store)

    mask = landmask(landshp, config.maskarcsec)

    # --- the source ------------------------------------------------------
    t0 = time()
    tiles = listedtiles(sys, tilelist, config.region)
    isempty(tiles) && error("the tile list and region select no tiles")
    real = realtiles(sys, tiledir, config.real)
    filter!(p -> p.first in Set(tiles), real)
    ids = TileIds(sys, tiles)
    say("tile list: $(length(tiles)) listed tiles of $(DGG.ncells(sys, 0)) " *
        "($(round(100 * length(tiles) / DGG.ncells(sys, 0); digits = 1))%), " *
        "$(length(real)) backed by a real GeoTIFF, " *
        @sprintf("%.3e pixels, %s", Float64(ids.n), secs(time() - t0)))
    for (o, p) in sort!(collect(real))
        say("  real: $(tilestem(sys, DGG.LevelIndex(0, o)))  " *
            "$(round(filesize(p) / 2^10)) KiB")
    end

    dem = TiledDEM(sys, ids, tiles; realtiles = real, mask = mask,
        cachesize = config.cache, stripes = config.stripes)
    t0 = time()
    srcgrid = DGG.PartialGrid(sys, 1, ids)
    srcspace = DGG.DGGSpace(srcgrid; chunklevel = 0)
    say("source space: $(GR.nchunks(srcspace)) chunks over " *
        @sprintf("%.3e pixels, built in %s", Float64(DGG.ncells(srcgrid)), secs(time() - t0)))
    check("source chunks are the listed tiles", GR.nchunks(srcspace) == length(tiles);
        detail = "$(GR.nchunks(srcspace)) chunks, $(length(tiles)) tiles")

    # --- the destination chunks ------------------------------------------
    todochunks = load_chunklist(chunklist)
    if todochunks === nothing
        t0 = time()
        todochunks = covering_chunks(sys7, sys, tiles, config.ancestor;
            nthreads = max(1, Threads.nthreads() - 1))
        save_chunklist(chunklist, config.ancestor, todochunks)
        say("chunks: $(length(todochunks)) level-$(config.ancestor) chunks cover the " *
            "tiles, computed in $(secs(time() - t0)), cached at $chunklist")
    else
        say("chunks: $(length(todochunks)) read from $chunklist")
    end
    if !isempty(config.chunks)
        todochunks = sort!(unique(config.chunks))
        say("chunks: overridden to $(length(todochunks)) named by `chunks`")
    end
    if 0 < config.maxchunks < length(todochunks)
        todochunks = todochunks[round.(Int, range(1, length(todochunks);
            length = config.maxchunks))]
        say("chunks: limited to $(config.maxchunks) by `maxchunks`")
    end
    ncell = length(todochunks) * capacity
    say(@sprintf("work: %d chunks x %d cells = %.4e level-%d cells, %.1f GiB dense f32",
        length(todochunks), capacity, Float64(ncell), config.level, 4 * ncell / 2^30))

    # --- the store -------------------------------------------------------
    store = openstore(config, sys7, capacity)
    layout = store.layout
    check("store layout is level $(config.level) over level-$(config.ancestor) chunks",
        DGG.level(layout) == config.level && layout.ancestor_level == config.ancestor &&
        layout.capacity == capacity;
        detail = "level $(DGG.level(layout)), ancestor $(layout.ancestor_level), " *
                 "capacity $(layout.capacity)")

    done = config.resume ? donechunks(donelog, config.store, "elevation") : Set{Int}()
    todo = count(c -> !(c in done), todochunks)
    say("resume: $(length(done)) chunks already written, $todo of " *
        "$(length(todochunks)) to do")

    if config.dryrun
        say("dryrun: stopping before any regrid")
        return FAILURES[]
    end

    # --- the run ---------------------------------------------------------
    say("launching $(config.workers) workers on $(Threads.nthreads()) threads, " *
        "batch=$(config.batch), budget=$(config.budget)")
    t0 = time()
    p = runchunks(todochunks, dem, srcspace, sys7, store, layout, done, donelog, config)
    wall = time() - t0
    say(@sprintf("RUN DONE  %d chunks computed, %d skipped, %.4e cells in %s = %.0f cells/s aggregate",
        p.done, p.skipped, Float64(p.cells), hours(wall), p.cells / max(wall, 1e-9)))
    say(@sprintf("source tiles decoded: %d real, %d synthetic; %.3e of %.3e pixels valid (%.2f%% land)",
        dem.nreal[], dem.nsynthetic[], Float64(dem.nland[]), Float64(dem.npixels[]),
        100 * dem.nland[] / max(dem.npixels[], 1)))
    say(@sprintf("destination NaN fraction: %.3f%% of %.4e written cells",
        100 * p.nan / max(p.cells, 1), Float64(p.cells)))
    say(@sprintf("peak RSS %.2f GiB", rssgib()))

    if config.checks
        sm = SourceMask(mask, DGG.levelgrid(sys, 0), Set(tiles))
        verify(config, sys7, layout, todochunks, sm, real, dem,
            donechunks(donelog, config.store, "elevation"))
    end

    say("total wall $(hours(time() - STARTED[])), " *
        (FAILURES[] == 0 ? "NO FAILURES" : "$(FAILURES[]) FAILURE(S)"))
    return FAILURES[]
end

exit(main() == 0 ? 0 : 1)
