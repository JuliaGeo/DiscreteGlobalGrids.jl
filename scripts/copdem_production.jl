# Copernicus DEM GLO-90 -> IGEO7 level 12, into one global ancestor-subzone Zarr
# store, written by W worker tasks in one process.
#
# Edit `CONFIG` below, then run it:
#
#     julia --project=bench -t 26 --gcthreads=8 scripts/copdem_production.jl
#
# The GC thread count carries NO second field. `--gcthreads=N,1` turns on Julia's
# concurrent page sweeper, which madvises freed pages from a background thread
# while the workers run; the 2026-08-21 run died of a SIGSEGV that is exactly what
# a page released out from under a live object looks like. `gcguard` below refuses
# to start under it. See regrid-notes/2026-08-21-polar-segfault.md.
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
# `CONFIG.source` selects real, lazily downloaded GLO-90 GeoTIFFs or the analytic
# field in `copdem_synthetic.jl`. The tile LIST is real in both modes: Copernicus
# ships ~26 450 of the 64 800 tiles in the 1x1-degree lattice, and the source is a
# `PartialGrid` over exactly those. A destination chunk over open ocean therefore
# pairs with no source at all and stays nodata without a network request. Store
# I/O and resume live in `copdem_store.jl`.

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
import Downloads
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
    source      = :real,    # :real (lazy AWS tiles) or :synthetic
    store       = "/home/asinghvi17/geo/dggstores/copdem90-igeo7-l12-real.zarr",
    region      = nothing,  # nothing for the globe, or [(w, e, s, n), ...] boxes
    maskarcsec  = 15,       # land-mask lattice, arcseconds; 0 disables the mask
    real        = :auto,    # local GeoTIFF overrides: :auto, :none, or ["stem", ...]
    tilecache   = get(ENV, "COPDEM_TILE_CACHE",
                      joinpath(@__DIR__, "..", "bench", "data", "CopernicusDEM", "tiles")),
    tilebaseurl = "https://copernicus-dem-90m.s3.amazonaws.com",
    retries     = 4,        # total GET attempts for a transient failure
    backoff     = 1.0,      # seconds; doubled between attempts
    timeout     = 600.0,    # seconds per tile GET
    workers     = 0,        # concurrent worker tasks; 0 = size from `cores`
    cores       = 24,       # the core budget `workers` is sized to hold
    shape       = :outer,   # :outer or :inner; see `workercount`
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
    allowsweeper = false,   # run anyway under `--gcthreads=N,1`; see `gcguard`
    malloctrim  = 32 * 2^20, # glibc M_TRIM_THRESHOLD; 0 leaves glibc alone
    data        = get(ENV, "RASTERDATASOURCES_PATH",
                      joinpath(@__DIR__, "..", "bench", "data")),
)

# ===========================================================================
# Logging
# ===========================================================================

const LOGLOCK = ReentrantLock()
const STARTED = Ref(time())
# Atomic: workers report their own chunk failures, so `+= 1` would lose some.
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
    lock(LOGLOCK) do
        ok || Threads.atomic_add!(FAILURES, 1)
        println(stamp(), "  ", ok ? "PASS  " : "FAIL  ", rpad(name, 56), detail)
        flush(stdout)
    end
    return ok
end

"""
    rssgib() -> Float64

CURRENT resident set size in GiB.

Deliberately not `Sys.maxrss()`, which is `ru_maxrss`: a monotone high-water
mark. A heartbeat printing that reports the largest transient the run has ever
had and never comes back down, so it looks like a staircase of plateaus when it
is really a record of peaks. On the GLO-90 run one sub-minute excursion to
81.65 GiB made every heartbeat for the next half hour claim "83.7 GiB" while the
true resident size was 64 GiB. Peaks are worth reporting — see
[`peakrssgib`](@ref) — but they have to be labelled as peaks.
"""
function rssgib()
    try
        return parse(Int, split(read("/proc/self/statm", String))[2]) * 4096 / 2^30
    catch
        return Sys.maxrss() / 2^30   # non-Linux fallback: a peak, but better than nothing
    end
end

"Peak resident set size in GiB since the process started."
peakrssgib() = Sys.maxrss() / 2^30

secs(t) = @sprintf("%.1f s", t)
hours(t) = @sprintf("%.2f h", t / 3600)

"""
    tunemalloc(trim) -> Bool

Freeze glibc's malloc thresholds at `trim` bytes; `false` if that did not take.

glibc raises its `mmap` and `trim` thresholds every time a large mmapped block
is freed — up to a 32 MiB cap — so that later allocations of that size come from
arena heap rather than a fresh mapping. This pipeline's hot allocations sit
squarely inside that range: a mid-latitude decoded GLO-90 tile is
`1200 x 1200 x 4 B` = 5.49 MiB, and the lazy plan's weight blocks are the same
order. Within seconds of starting, every one of them is served from arena heap,
and freeing it returns nothing to the operating system.

Measured (`regrid-notes/2026-08-21-memory-attribution.md`): resident memory
settles at about three times the live Julia heap and stays there through two
full `GC.gc(true)` passes — 15.9 GiB resident against 3.0 GiB live. It is not
the garbage collector; `--heap-size-hint` buys 14 % of it for half the
throughput.

Setting any one of the thresholds explicitly sets glibc's `no_dyn_threshold`,
which stops the growth, and that flag — not the number — is the whole effect:
8 MiB and 32 MiB measure identically. So pass a large one, which trims just as
well as the 128 KiB default with far fewer syscalls. `malloctrim = 0` opts out.
"""
function tunemalloc(trim::Integer)
    (trim > 0 && Sys.islinux()) || return false
    return try
        # M_TRIM_THRESHOLD is -1 in glibc's malloc.h; mallopt returns 1 on success.
        ccall(:mallopt, Cint, (Cint, Cint), Cint(-1), Cint(trim)) == 1
    catch
        false
    end
end

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
# Lazy real GLO-90 acquisition
# ===========================================================================

"The production AWS Open Data root for Copernicus DEM GLO-90 COGs."
const COPDEM_GLO90_BASEURL = "https://copernicus-dem-90m.s3.amazonaws.com"

"""
    LazyCopernicusTiles(sys, listed; cachedir, baseurl, retries, backoff, timeout)

A single-flight, persistent source for real Copernicus GLO-90 GeoTIFFs.

No request is made by construction. On the first access to a listed tile,
[`tilepath!`](@ref) downloads its COG to `cachedir/<stem>.tif.part` and renames
it atomically to `cachedir/<stem>.tif`. A later access trusts only the final
name. One lock per listed tile makes concurrent worker requests single-flight
within the production process.

An unlisted tile is ocean in this source. [`loadtile`](@ref) returns an all-NaN
tile for it without constructing a URL or touching the network. The production
`PartialGrid` omits these chunks altogether, which has the same downstream
nodata semantics without allocating them.
"""
struct LazyCopernicusTiles{S<:DGG.CopernicusDEMSystem}
    sys::S
    listed::Set{Int}
    cachedir::String
    baseurl::String
    retries::Int
    backoff::Float64
    timeout::Float64
    locks::Dict{Int,ReentrantLock}
    downloader::Downloads.Downloader
    ndownloads::Threads.Atomic{Int}
end

function LazyCopernicusTiles(sys::DGG.CopernicusDEMSystem, listed;
        cachedir::AbstractString,
        baseurl::AbstractString = COPDEM_GLO90_BASEURL,
        retries::Integer = 4, backoff::Real = 1.0, timeout::Real = 600.0)
    CD.lat_intervals(sys) == 1200 || throw(ArgumentError(
        "the lazy AWS provider is for GLO-90; got $sys"))
    retries >= 1 || throw(ArgumentError("retries must be at least 1"))
    backoff >= 0 || throw(ArgumentError("backoff must be non-negative"))
    timeout > 0 || throw(ArgumentError("timeout must be positive"))
    tiles = Set(Int.(listed))
    locks = Dict(t => ReentrantLock() for t in tiles)
    return LazyCopernicusTiles(sys, tiles, abspath(String(cachedir)),
        String(rstrip(String(baseurl), '/')), Int(retries), Float64(backoff),
        Float64(timeout), locks, Downloads.Downloader(), Threads.Atomic{Int}(0))
end

"The final local cache path of a tile; no filesystem or network access."
function tilecachepath(provider::LazyCopernicusTiles, ordinal::Int)
    stem = tilestem(provider.sys, DGG.LevelIndex(0, ordinal))
    return joinpath(provider.cachedir, stem * ".tif")
end

"The verified public S3 URL of a tile."
function tileurl(provider::LazyCopernicusTiles, ordinal::Int)
    stem = tilestem(provider.sys, DGG.LevelIndex(0, ordinal))
    return string(provider.baseurl, "/", stem, "/", stem, ".tif")
end

"Whether an HTTP status is worth retrying."
_transientstatus(status::Integer) = status == 0 || status == 408 || status == 429 ||
    500 <= status < 600

"""
    tilepath!(provider, ordinal) -> Union{String,Nothing}

Return the cached path for a listed tile, downloading it atomically if needed.
Return `nothing` immediately for an unlisted (ocean) tile.
"""
function tilepath!(provider::LazyCopernicusTiles, ordinal::Int)
    ordinal in provider.listed || return nothing
    path = tilecachepath(provider, ordinal)
    isfile(path) && return path
    return lock(provider.locks[ordinal]) do
        isfile(path) && return path
        mkpath(provider.cachedir)
        part = path * ".part"
        url = tileurl(provider, ordinal)
        stem = splitext(basename(path))[1]
        last_error = nothing
        for attempt in 1:provider.retries
            try
                # `Downloads.download` truncates an existing output path, so a
                # stale .part left by a killed run is safe to reuse.
                Downloads.download(url, part; timeout = provider.timeout,
                    downloader = provider.downloader)
                filesize(part) > 0 || error("downloaded an empty object from $url")
                Base.Filesystem.rename(part, path)
                Threads.atomic_add!(provider.ndownloads, 1)
                return path
            catch err
                last_error = err
                status = err isa Downloads.RequestError ? err.response.status : 0
                rm(part; force = true)
                if status == 403 || status == 404
                    error("listed Copernicus GLO-90 tile $stem returned HTTP $status from $url")
                end
                if !_transientstatus(status) || attempt == provider.retries
                    break
                end
                delay = provider.backoff * 2.0^(attempt - 1)
                @warn "Copernicus tile download failed; retrying" stem attempt attempts=provider.retries delay status exception=(err, catch_backtrace())
                sleep(delay)
            end
        end
        detail = sprint(showerror, last_error)
        error("failed to download listed Copernicus GLO-90 tile $stem after " *
              "$(provider.retries) attempt(s) from $url: $detail")
    end
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
    TiledDEM(sys, ids, tiles; realtiles, provider, mask, cachesize, stripes)

Every listed tile's pixels as one `Float32` vector in the source grid's own
position order, chunk `k` being listed tile `k` — so a read is always tile
aligned and the regridder's source chunks are the DEM's own tiles.

A tile named in `realtiles` decodes from its local GeoTIFF. Otherwise a non-null
`provider` lazily fetches and decodes the public GLO-90 COG; with no provider it
is fabricated by `synthetic_tile` (see `copdem_synthetic.jl`). The sources are
indistinguishable downstream: `NaN` is the regridder's invalid sentinel whatever
produced it. The Natural Earth mask remains exclusively in the synthetic path,
exactly as before; real COGs bring their own values/nodata.

The cache is striped — `stripes` independent LRUs under `stripes` locks — and a
tile is built OUTSIDE its lock, so the workers do not serialise on one mutex for
the ~15 ms a 1200x1200 tile takes.
"""
struct TiledDEM{S<:DGG.CopernicusDEMSystem,I,C,P,M} <: DiskArrays.AbstractDiskArray{Float32,1}
    sys::S
    ids::I
    tiles::Vector{Int}
    chunks::C
    realtiles::Dict{Int,String}
    provider::P
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
    provider = nothing, mask, cachesize::Integer = 3072, stripes::Integer = 64)
    widths = [ids.offsets[k + 1] - ids.offsets[k] for k in 1:(length(tiles) - 1)]
    push!(widths, ids.n - ids.offsets[end])
    chunks = DiskArrays.GridChunks(DiskArrays.IrregularChunks(; chunksizes = widths))
    ns = Int(stripes)
    return TiledDEM(sys, ids, tiles, chunks, realtiles, provider, mask,
        [Dict{Int,Vector{Float32}}() for _ in 1:ns], [Int[] for _ in 1:ns],
        [ReentrantLock() for _ in 1:ns], max(1, Int(cachesize) ÷ ns),
        Threads.Atomic{Int}(0), Threads.Atomic{Int}(0),
        Threads.Atomic{Int}(0), Threads.Atomic{Int}(0))
end

Base.size(A::TiledDEM) = (A.ids.n,)
DiskArrays.eachchunk(A::TiledDEM) = A.chunks
DiskArrays.haschunks(::TiledDEM) = DiskArrays.Chunked()

# GDAL is entered from whichever worker first wants a real tile, and
# `tilevalues!` deliberately builds outside the stripe lock, so several workers
# can be in `readtile` at once. ArchGDAL states no thread-safety guarantee, and
# the whole holding is four 1x1-degree tiles decoded once each, so serialising
# them costs nothing measurable and removes the only native-library race in the
# run. (Not the 2026-08-21 crash: the cursor never reached the chunks that
# carry those tiles. It would have reached them ~1300 chunks later.)
const GDALLOCK = ReentrantLock()

"The decoded band of the GeoTIFF at `path`, validated and flattened into position order."
function readtile(path, sys, tile::DGG.LevelIndex)
    band = lock(GDALLOCK) do
        ArchGDAL.read(ds -> ArchGDAL.read(ds, 1), path)
    end
    lat, _ = CD.tilecorner(sys, tile)
    expected = (Int(CD.ncols_at(sys, lat)), Int(CD.lat_intervals(sys)))
    size(band) == expected || error(
        "Copernicus tile $(tilestem(sys, tile)) at $path has raster size " *
        "$(size(band)); expected $expected")
    return Float32.(vec(band))
end

"Decode a listed real tile, or return ocean nodata without a request."
function loadtile(provider::LazyCopernicusTiles, ordinal::Int)
    path = tilepath!(provider, ordinal)
    tile = DGG.LevelIndex(0, ordinal)
    if path === nothing
        lat, _ = CD.tilecorner(provider.sys, tile)
        return fill(NaN32, Int(CD.ncols_at(provider.sys, lat) *
                               CD.lat_intervals(provider.sys)))
    end
    return readtile(path, provider.sys, tile)
end

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
    v = if path !== nothing
        Threads.atomic_add!(A.nreal, 1)
        vals = readtile(path, A.sys, DGG.LevelIndex(0, ordinal))
        Threads.atomic_add!(A.nland, count(isfinite, vals))
        Threads.atomic_add!(A.npixels, length(vals))
        vals
    elseif A.provider === nothing
        Threads.atomic_add!(A.nsynthetic, 1)
        vals, nland = synthetic_tile(A.sys, DGG.LevelIndex(0, ordinal), A.mask)
        Threads.atomic_add!(A.nland, nland)
        Threads.atomic_add!(A.npixels, length(vals))
        vals
    else
        Threads.atomic_add!(A.nreal, 1)
        vals = loadtile(A.provider, ordinal)
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
        say(@sprintf("HEARTBEAT  %d/%d chunks (%d skipped) | %.3e cells | %.0f cells/s | elapsed %s | ETA %s | RSS %.1f GiB (peak %.1f) | NaN %.2f%%",
            done, p.total, p.skipped, Float64(p.cells), rate, hours(el), hours(eta),
            rssgib(), peakrssgib(), 100 * p.nan / max(p.cells, 1)))
        return nothing
    end
end

"""
    workercount(config) -> (W, shape)

The worker count `W`, sized from the `cores` budget rather than from a thread
count, and the parallel `shape` it assumes. `config.workers` overrides `W` when
set to something positive.

A work unit is one chunk, and how many cores one worker can hold depends
entirely on which of the two nested parallelisms is allowed to run:

  - `shape = :outer` sets [`GR.OUTER_PARALLEL`](@ref) in the worker body. The
    weight build stays serial and the lazy plan's source wave is the only
    parallelism inside a unit. Measured: **1.06 cores per worker**, flat in `W`
    while `W` is well below `nthreads`, and **65 k cells per core-second** — the
    most work per core of the two, but it needs `W ≈ cores` workers, and each
    worker holds a chunk's weights, so resident memory scales with `W`.

  - `shape = :inner` leaves `OUTER_PARALLEL` unset, so a narrow wave stands
    aside and ConservativeRegridding threads the weight build itself. Measured
    at 12 threads: one worker reaches 7.0–7.3 cores, the pool saturates at
    **`W ≈ nthreads/4`** holding ~78 % of it, and throughput is **60 k cells per
    core-second** — about 8 % less work per core, for a third of the workers and
    two thirds of the resident memory.

So at a fixed core budget `:outer` converts cores into cells slightly faster,
while `:inner` reaches the same cores with far fewer workers, and is the only
shape that keeps scaling once `cores` exceeds what `W ≈ cores` workers can be
given threads for. `:outer` is the default because the standing budget is 24
cores, where the throughput edge wins; switch to `:inner` for a larger budget or
when resident memory is the binding constraint.
"""
function workercount(config)
    shape = config.shape
    shape in (:outer, :inner) ||
        error("shape must be :outer or :inner, got $(repr(shape))")
    config.workers > 0 && return (config.workers, shape)
    config.cores > 0 || error("cores must be positive when workers is not set")
    # Cores per worker, measured; see above.
    W = shape === :outer ? cld(config.cores, 1.06) : cld(config.cores, 3.1)
    return (max(1, Int(W)), shape)
end

"""
    runchunks(chunks, dem, srcspace, sys7, store, layout, done, donelog, config) -> Progress

`W` tasks — [`workercount`](@ref) — pulling contiguous batches of chunks off one
atomic cursor, each computing its chunk and writing it straight to its own Zarr
file. Chunks are disjoint, so no two tasks ever touch the same file.

Batches are contiguous in canonical Z7 order, so a worker walks a geographically
connected run and its tile cache stays hot; pulling dynamically keeps the polar
chunks — which meet 200-360 tiles and cost ~5x a mid-latitude one — off the
critical path.

Under `shape = :outer` each worker body sets `GR.OUTER_PARALLEL`, so the nested
weight builds do not spawn their own tasks and oversubscribe the box; the lazy
plan's source wave is then the only parallelism inside a unit. Under
`shape = :inner` the scoped value is left alone, and `_wavesize` weighs that wave
against the threading inside the build — standing the wave down whenever one
source chunk carries the destination chunk, which is almost always.
"""
function runchunks(chunks, dem, srcspace, sys7, store, layout, done, donelog, config)
    W, shape = workercount(config)
    outer = shape === :outer
    B = config.batch
    p = Progress(length(chunks), config.heartbeat)
    cursor = Threads.Atomic{Int}(1)
    log = DoneLog(donelog)
    try
        # `false` is the scoped value's own default, so binding it is the same
        # as leaving it unset — the worker is the outermost scope either way.
        Threads.@sync for w in 1:W
            Threads.@spawn @with GR.OUTER_PARALLEL => outer begin
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
                            Threads.atomic_add!(FAILURES, 1)
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

"""
    gcguard(config)

Refuse to run under Julia's concurrent page sweeper.

`--gcthreads=N,M` asks for `N` mark threads and `M` (0 or 1) threads for the
*concurrent sweeping phase*. `M` defaults to 0; `M = 1` is an opt-in that hands a
background thread the job of releasing swept pages — `jl_concurrent_gc_threadfun`
→ `gc_free_pages` → `jl_gc_free_page` → `madvise` — while the worker threads keep
running. `madvise(MADV_DONTNEED)` on anonymous memory does not unmap the page; it
drops its contents, so the next read of anything still living there returns zeros.

On 2026-08-21 this run died at chunk ~121700 with

    [3718607] signal 11 (128): Segmentation fault
    getindex at ./essentials.jl:920 [inlined]
    _child_extent at GeometryOps .../dual_depth_first_search.jl:37 [inlined]
    dual_depth_first_search at .../dual_depth_first_search.jl:78

which is the `memoryrefget` *after* a `@boundscheck` that passed, on a freshly
built seven-element `Vector{SphericalCap}` that no other task can reach — an
`Array` whose length still read correctly while its data pointer did not. A page
released under a live object produces precisely that, and the previous run's
SIGTERM backtrace (same log, line 28153) caught the sweeper thread inside
`madvise` in this very workload. Two independent audits and ~1000 Antarctic
chunks re-run under `--check-bounds=yes` found no out-of-bounds access, no
unsafe/`ccall` code, and no unsynchronized shared state on the path.

So the configuration is the defect, and this is where the configuration is
written down. `allowsweeper = true` overrides, for someone deliberately testing
it.
"""
function gcguard(config)
    nsweep = Base.JLOptions().nsweepthreads
    nsweep == 0 && return nothing
    config.allowsweeper && return say(
        "WARNING: running with $nsweep concurrent GC sweep thread(s) by " *
        "`allowsweeper = true`; this is the 2026-08-21 SIGSEGV configuration")
    error("""
        this run was launched with Julia's concurrent page sweeper enabled \
        (--gcthreads=…,$nsweep). It released a page under a live object on \
        2026-08-21 and killed the global run at chunk ~121700. Drop the second \
        field: `--gcthreads=$(Base.JLOptions().nmarkthreads)`. Set \
        `allowsweeper = true` to run under it anyway.""")
end

function main(config = CONFIG)
    gcguard(config)
    println("="^92)
    config.source in (:real, :synthetic) || error(
        "source must be :real or :synthetic, got $(repr(config.source))")
    println(stamp(), "  copdem_production.jl — GLO-$(config.res) -> IGEO7 level " *
                     "$(config.level), level-$(config.ancestor) chunks, " *
                     "$(uppercase(String(config.source))) elevations")
    println("  julia $(VERSION)  threads=$(Threads.nthreads())  " *
            "gcmark=$(Base.JLOptions().nmarkthreads)  " *
            "gcsweep=$(Base.JLOptions().nsweepthreads)  pid=$(getpid())")
    println("  ", config)
    println("="^92)
    flush(stdout)

    say(config.malloctrim > 0 ?
        "malloc: M_TRIM_THRESHOLD frozen at $(config.malloctrim >> 20) MiB — " *
        (tunemalloc(config.malloctrim) ? "ok" : "MALLOPT FAILED, expect ~3x resident memory") :
        "malloc: left at glibc defaults (malloctrim = 0); expect ~3x resident memory")

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
    provider = config.source === :real ? LazyCopernicusTiles(sys, tiles;
        cachedir = config.tilecache, baseurl = config.tilebaseurl,
        retries = config.retries, backoff = config.backoff,
        timeout = config.timeout) : nothing
    ids = TileIds(sys, tiles)
    say("tile list: $(length(tiles)) listed tiles of $(DGG.ncells(sys, 0)) " *
        "($(round(100 * length(tiles) / DGG.ncells(sys, 0); digits = 1))%), " *
        "$(length(real)) local GeoTIFF overrides, " *
        @sprintf("%.3e pixels, %s", Float64(ids.n), secs(time() - t0)))
    provider === nothing || say("lazy GLO-90 source: $(provider.baseurl), cache $(provider.cachedir)")
    for (o, p) in sort!(collect(real))
        say("  real: $(tilestem(sys, DGG.LevelIndex(0, o)))  " *
            "$(round(filesize(p) / 2^10)) KiB")
    end

    dem = TiledDEM(sys, ids, tiles; realtiles = real, provider = provider, mask = mask,
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
    W, shape = workercount(config)
    say("launching W=$W workers, shape=$shape, sized to hold $(config.cores) cores " *
        "(nthreads=$(Threads.nthreads())), batch=$(config.batch), budget=$(config.budget)")
    shape === :inner && Threads.nthreads() < 4W &&
        say("NOTE shape=:inner saturates near W=nthreads/4; " *
            "$(Threads.nthreads()) threads will not fill $W workers")
    t0 = time()
    p = runchunks(todochunks, dem, srcspace, sys7, store, layout, done, donelog, config)
    wall = time() - t0
    say(@sprintf("RUN DONE  %d chunks computed, %d skipped, %.4e cells in %s = %.0f cells/s aggregate",
        p.done, p.skipped, Float64(p.cells), hours(wall), p.cells / max(wall, 1e-9)))
    say(@sprintf("source tiles decoded: %d real, %d synthetic; %.3e of %.3e pixels valid (%.2f%% land)",
        dem.nreal[], dem.nsynthetic[], Float64(dem.nland[]), Float64(dem.npixels[]),
        100 * dem.nland[] / max(dem.npixels[], 1)))
    provider === nothing || say("source tile downloads this process: $(provider.ndownloads[])")
    say(@sprintf("destination NaN fraction: %.3f%% of %.4e written cells",
        100 * p.nan / max(p.cells, 1), Float64(p.cells)))
    say(@sprintf("RSS %.2f GiB now, %.2f GiB peak", rssgib(), peakrssgib()))

    if config.checks
        sm = SourceMask(mask, DGG.levelgrid(sys, 0), Set(tiles))
        verify(config, sys7, layout, todochunks, sm, real, dem,
            donechunks(donelog, config.store, "elevation"))
    end

    say("total wall $(hours(time() - STARTED[])), " *
        (FAILURES[] == 0 ? "NO FAILURES" : "$(FAILURES[]) FAILURE(S)"))
    return FAILURES[]
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main() == 0 ? 0 : 1)
end
