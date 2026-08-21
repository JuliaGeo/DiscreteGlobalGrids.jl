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
    workers     = 0,        # concurrent worker tasks; 0 = size from `cores`
    cores       = 24,       # the core budget `workers` is sized to hold
    shape       = :outer,   # :outer or :inner; see `workercount`
    batch       = 8,        # chunks handed out per pull, at most; see `taper`
    taper       = true,     # shrink the batch as the queue drains
    budget      = 2^30,     # lazy-regrid byte budget, per worker
    schedule    = :affinity, # :affinity (tile-affinity walk) or :canonical
    cachepolicy = :refcount, # :refcount (graph-driven) or :lru
    cache       = 3072,     # tiles held across all workers, `cachepolicy = :lru` only
    stripes     = 64,       # independent locks over that cache, ditto
    refinegraph = false,    # narrow the graph past the regridder's own pairing
    prefetch    = 0,        # columns of lookahead; 0 disables the prefetcher
    fetchconc   = 0,        # concurrent source loads; 0 means "as many as workers"
    fetchdelay  = 0.0,      # seconds of fake latency per tile build; a test knob
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
# The last run's cache statistics, for a harness that calls `main` in-process
# and wants the numbers rather than the log line.
const LASTCACHE = Ref{Any}(nothing)

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
include("copdem_policy.jl")

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
    SubtreeIds(sys, parents, level)

The level-`level` descendants of `parents` as one lazy ascending vector of cell
ids, subtree `k` being `parents[k]`.

Lazy because a grid stores its ids by reference and never materialises them, so
the source's 3.8e10 pixels cost an offsets table of 26 450 integers and the
destination's 5.4e10 cells cost 66 178. It also declares itself sorted, which
saves the O(n) scan the grid would otherwise run: `parents` is ascending and
descendant ranges are disjoint, so the concatenation ascends by construction.

Both sides of this run are built from it. The source is the listed 1x1-degree
tiles resolved to pixels ([`TileIds`](@ref)); the destination is the covering
level-`ancestor` chunks resolved to level-`level` cells, which is the space the
chunk dependency graph is built against.
"""
struct SubtreeIds{G,ID} <: AbstractVector{ID}
    complete::G
    starts::Vector{Int}      # first level-`level` position of each subtree
    offsets::Vector{Int}     # cells before each subtree; `offsets[1] == 0`
    n::Int
end

function SubtreeIds(sys, parents::AbstractVector, level::Int)
    complete = DGG.levelgrid(sys, level)
    starts = Vector{Int}(undef, length(parents))
    offsets = Vector{Int}(undef, length(parents))
    acc = 0
    for (k, p) in enumerate(parents)
        r = DGG.descendant_range(sys, p, level)
        starts[k] = Int(first(r))
        offsets[k] = acc
        acc += length(r)
    end
    # A grid only accepts its own system's cell type, and that differs by system
    # — `LevelIndex` for Copernicus, `Z7Cell` for IGeo7 — so take it from the
    # grid rather than naming one.
    ID = typeof(DGG.cellindex(complete, first(starts)))
    return SubtreeIds{typeof(complete),ID}(complete, starts, offsets, acc)
end

"The pixel (level-1) cell ids of the listed tiles."
TileIds(sys, tiles::Vector{Int}) =
    SubtreeIds(sys, [DGG.LevelIndex(0, t) for t in tiles], 1)

Base.size(v::SubtreeIds) = (v.n,)
Base.IndexStyle(::Type{<:SubtreeIds}) = IndexLinear()

@inline function Base.getindex(v::SubtreeIds, i::Int)
    @boundscheck (1 <= i <= v.n) || throw(BoundsError(v, i))
    k = searchsortedlast(v.offsets, i - 1)
    return DGG.cellindex(v.complete, v.starts[k] + (i - 1 - v.offsets[k]))
end

DGG.Helpers.strictly_increasing(::SubtreeIds) = true

"Which subtree holds position `p`, and the offset within it."
@inline function tileat(v::SubtreeIds, p::Int)
    k = searchsortedlast(v.offsets, p - 1)
    return k, p - v.offsets[k]
end

"""
    TileBuilder(sys, tiles, realtiles, mask, delay)

Where a source tile's pixels come from, and nothing about when.

A tile named in `realtiles` decodes from its GeoTIFF and brings its own nodata.
Every other tile is fabricated by `synthetic_tile` (see `copdem_synthetic.jl`).
The two are indistinguishable downstream: `NaN` is the regridder's invalid
sentinel whatever produced it.

Split out from [`TiledDEM`](@ref) so the cache can hold the *loader* rather than
the array that holds the cache. `delay` sleeps before every build, which is how
the AWS run's fetch latency is imposed on the synthetic one for a scheduling
measurement; it is a test knob and defaults to zero.
"""
struct TileBuilder{S<:DGG.CopernicusDEMSystem,M}
    sys::S
    tiles::Vector{Int}
    realtiles::Dict{Int,String}
    mask::M
    delay::Float64
    nreal::Threads.Atomic{Int}
    nsynthetic::Threads.Atomic{Int}
    nland::Threads.Atomic{Int}
    npixels::Threads.Atomic{Int}
end

TileBuilder(sys, tiles::Vector{Int}, realtiles::Dict{Int,String}, mask;
    delay::Real = 0.0) =
    TileBuilder(sys, tiles, realtiles, mask, Float64(delay),
        Threads.Atomic{Int}(0), Threads.Atomic{Int}(0),
        Threads.Atomic{Int}(0), Threads.Atomic{Int}(0))

# GDAL is entered from whichever worker first wants a real tile, and the caches
# deliberately build outside their locks, so several workers can be in `readtile`
# at once. ArchGDAL states no thread-safety guarantee, and the whole holding is
# four 1x1-degree tiles decoded once each, so serialising them costs nothing
# measurable and removes the only native-library race in the run. (Not the
# 2026-08-21 crash: the cursor never reached the chunks that carry those tiles.
# It would have reached them ~1300 chunks later.)
const GDALLOCK = ReentrantLock()

"The decoded band of the GeoTIFF at `path`, flattened into position order."
readtile(path) = lock(GDALLOCK) do
    vec(ArchGDAL.read(ds -> ArchGDAL.read(ds, 1), path))
end

"Build source chunk `k` — listed tile `k` — from its GeoTIFF or from thin air."
function buildtile(b::TileBuilder, k::Int)
    b.delay > 0 && sleep(b.delay)
    ordinal = b.tiles[k]
    path = get(b.realtiles, ordinal, nothing)
    vals = if path === nothing
        Threads.atomic_add!(b.nsynthetic, 1)
        v, nland = synthetic_tile(b.sys, DGG.LevelIndex(0, ordinal), b.mask)
        Threads.atomic_add!(b.nland, nland)
        v
    else
        Threads.atomic_add!(b.nreal, 1)
        v = readtile(path)
        Threads.atomic_add!(b.nland, count(isfinite, v))
        v
    end
    Threads.atomic_add!(b.npixels, length(vals))
    return vals
end

"""
    TiledDEM(ids, builder, cache)

Every listed tile's pixels as one `Float32` vector in the source grid's own
position order, chunk `k` being listed tile `k` — so a read is always tile
aligned and the regridder's source chunks are the DEM's own tiles, which is what
makes the chunk dependency graph's source side the DEM's own tile list.

`cache` decides which tiles are resident and for how long: a
[`RefCountCache`](@ref) driven by that graph, or the [`StripedLRUCache`](@ref)
this run used before it had one.
"""
struct TiledDEM{I,C,K} <: DiskArrays.AbstractDiskArray{Float32,1}
    ids::I
    builder::TileBuilder
    cache::K
    chunks::C
end

function TiledDEM(ids::SubtreeIds, builder::TileBuilder, cache)
    n = length(builder.tiles)
    widths = [ids.offsets[k + 1] - ids.offsets[k] for k in 1:(n - 1)]
    push!(widths, ids.n - ids.offsets[end])
    chunks = DiskArrays.GridChunks(DiskArrays.IrregularChunks(; chunksizes = widths))
    return TiledDEM(ids, builder, cache, chunks)
end

Base.size(A::TiledDEM) = (A.ids.n,)
DiskArrays.eachchunk(A::TiledDEM) = A.chunks
DiskArrays.haschunks(::TiledDEM) = DiskArrays.Chunked()

function DiskArrays.readblock!(A::TiledDEM, out, r::AbstractUnitRange)
    ntiles = length(A.builder.tiles)
    p = first(r)
    while p <= last(r)
        k, off = tileat(A.ids, p)
        width = (k < ntiles ? A.ids.offsets[k + 1] : A.ids.n) - A.ids.offsets[k]
        stop = min(last(r), A.ids.offsets[k] + width)
        seg = (p - first(r) + 1):(stop - first(r) + 1)
        # `k` is the source CHUNK number, which is what the dependency graph and
        # both caches are keyed by; only the builder knows it is tile `tiles[k]`.
        v = gettile!(A.cache, k)
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
# The dependency graph, and the walk order it implies
# ===========================================================================

# A conservative lon/lat narrow phase for this specific pair, available to
# `chunk_dependency_graph` as its `refine` hook — and OFF by default. Read on.
#
# The broad phase is cap-versus-cap, and a circle circumscribing a square already
# inflates by its half-diagonal, so on this workload the cap relation carries
# 1.75x the edges the exact geometry has. A Copernicus tile is exactly a
# 1-degree lon/lat box — always 1x1 in EXTENT; it is the pixel count that shrinks
# poleward, not the footprint — and a spherical cap has an exact lon/lat bounding
# box, so two boxes that do not overlap cannot intersect. That takes the
# inflation to 1.35x for no measurable build cost, and it drops none of the
# exact edges: the narrow phase was checked against an independently built exact
# adjacency and kept all 186 069 of them.
#
# **And it is still the wrong graph to schedule on.** Refcount eviction does not
# need a superset of the *true geometry*; it needs a superset of what the
# *executor actually reads*, and the executor pairs its chunks by exactly the
# cap test the broad phase uses. Every pair the narrow phase removes is a pair
# the regridder will still read the tile for. Measured on eleven chunks: the
# refined graph credited chunk 100 with one tile and the regrid then demanded
# three, eight such demands over the eleven. No data was wrong — an uncredited
# demand is served — but the tile had already been freed, so the reload the
# design exists to prevent came back.
#
# So the narrow phase stays here, behind `refinegraph`, for the day the
# regridder's own pairing gets tighter. Turning it on today trades the whole
# at-most-once property for 25 % fewer edges.
#
# `PAD` absorbs the half-pixel outer frame `node_extent` adds to a tile box and
# any rounding, so the test only ever errs toward keeping an edge.
const PAD = 0.01

"Half-width in longitude degrees of a cap's bounding box; 180 at a pole."
function lonhalfwidth(latdeg::Float64, rdeg::Float64)
    abs(latdeg) + rdeg >= 90 - PAD && return 180.0
    s = sind(rdeg) / cosd(latdeg)
    s >= 1 && return 180.0
    return rad2deg(asin(s))
end

"Circular longitude separation, in degrees."
circulardlon(a::Float64, b::Float64) = abs(mod(a - b + 180.0, 360.0) - 180.0)

function boxesoverlap(dlat::Float64, dlon::Float64, dlathalf::Float64,
        dlonhalf::Float64, tlat::Float64, tlon::Float64)
    abs(dlat - tlat) <= dlathalf + 0.5 + PAD || return false
    dlonhalf >= 180 - PAD && return true
    return circulardlon(dlon, tlon) <= dlonhalf + 0.5 + PAD
end

"""
    dagplan(sys, sys7, tiles, chunks, srcspace, config) -> NamedTuple

The tile-to-chunk dependency graph over exactly the work this run will do, and
the order to walk it in.

The destination space is built over `chunks` in the order they are given, so
graph destination `d` is `chunks[d]` and graph source `s` is `tiles[s]` — the
same numbering the source space and both caches use. The order is a *separate*
permutation, `order[p] -> d`, applied by the cursor; the graph itself is never
built in it.

Returns `(graph, order, seconds, edges, dstspace)`. Building it costs about a
second on the real pair, almost all of it the destination space; the graph
itself is 0.12 s for 66 178 x 26 475 chunks.
"""
function dagplan(sys, sys7, tiles::Vector{Int}, chunks::Vector{Int}, srcspace, config)
    t0 = time()
    g5 = DGG.levelgrid(sys7, config.ancestor)
    dstids = SubtreeIds(sys7, [DGG.cellindex(g5, c) for c in chunks], config.level)
    dstgrid = DGG.PartialGrid(sys7, config.level, dstids)
    dstspace = DGG.DGGSpace(dstgrid; chunklevel = config.ancestor)
    tspace = time() - t0

    # Tile centres, and each destination cap as a lon/lat box, for `refine`.
    corner = [CD.tilecorner(sys, DGG.LevelIndex(0, t)) for t in tiles]
    tlat = [Float64(c[1]) + 0.5 for c in corner]
    tlon = [Float64(c[2]) + 0.5 for c in corner]
    caps = GR.chunkextents(dstspace)
    dlat = [rad2deg(atan(c.point[3], hypot(c.point[1], c.point[2]))) for c in caps]
    dlon = [rad2deg(atan(c.point[2], c.point[1])) for c in caps]
    drad = [rad2deg(c.radius) for c in caps]
    dhalf = [lonhalfwidth(dlat[i], drad[i]) for i in eachindex(caps)]
    refine = config.refinegraph ?
             ((d, s) -> boxesoverlap(dlat[d], dlon[d], drad[d], dhalf[d],
                 tlat[s], tlon[s])) : nothing

    t1 = time()
    radius = Float64(GR.support_radius(DGG.Conservative(), srcspace))
    graph = GR.chunk_dependency_graph(dstspace, srcspace; radius, refine)
    tgraph = time() - t1

    # Sweep the tiles in Morton order over the 1-degree lattice and emit a chunk
    # when its last tile is swept. `:canonical` keeps the ascending chunk order
    # the run used before, as the control for an A/B.
    t2 = time()
    order = if config.schedule === :affinity
        keys = [morton2(Int(c[2]) + 180, Int(c[1]) + 90) for c in corner]
        affinity_order(length(chunks), d -> GR.sourcesof(graph, d), keys)
    elseif config.schedule === :canonical
        collect(1:length(chunks))
    else
        error("schedule must be :affinity or :canonical, got $(repr(config.schedule))")
    end
    torder = time() - t2

    return (graph = graph, order = order, dstspace = dstspace, radius = radius,
        tspace = tspace, tgraph = tgraph, torder = torder,
        edges = length(order) == 0 ? 0 : sum(d -> GR.sourcedegree(graph, d),
            1:length(chunks)))
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
    runchunks(chunks, order, dem, srcspace, sys7, store, layout, skipped, donelog,
              prefetcher, config) -> Progress

`W` tasks — [`workercount`](@ref) — pulling batches of positions off one
[`GuidedSchedule`](@ref), each computing its chunk and writing it straight to its
own Zarr file. Chunks are disjoint, so no two tasks ever touch the same file.

`order[p]` is the graph destination number to run at position `p`, and
`chunks[order[p]]` is the level-`ancestor` chunk itself. The order is a priority
sequence, never a partition: pulling dynamically is what keeps the polar chunks —
which meet 200-360 tiles and cost ~5x a mid-latitude one — off the critical path,
and tapering the batch is what keeps the last pull from handing one worker eight
of them while everybody else is idle.

Every chunk is retired from the cache exactly once, on every terminal outcome and
before its result is written: the Zarr write reads no source data, so holding
tiles across it only inflates residency. Chunks skipped on resume were retired
before the run started and are not in `order` at all.

Under `shape = :outer` each worker body sets `GR.OUTER_PARALLEL`, so the nested
weight builds do not spawn their own tasks and oversubscribe the box; the lazy
plan's source wave is then the only parallelism inside a unit. Under
`shape = :inner` the scoped value is left alone, and `_wavesize` weighs that wave
against the threading inside the build — standing the wave down whenever one
source chunk carries the destination chunk, which is almost always.
"""
function runchunks(chunks, order, sched, dem, srcspace, sys7, store, layout, skipped,
        donelog, prefetcher, config)
    W, shape = workercount(config)
    outer = shape === :outer
    p = Progress(length(chunks), config.heartbeat)
    p.skipped = skipped
    log = DoneLog(donelog)
    try
        # `false` is the scoped value's own default, so binding it is the same
        # as leaving it unset — the worker is the outermost scope either way.
        Threads.@sync for w in 1:W
            Threads.@spawn @with GR.OUTER_PARALLEL => outer begin
                while true
                    batch = claim!(sched)
                    batch === nothing && break
                    advance!(prefetcher)
                    for pos in batch
                        d = order[pos]
                        ch = chunks[d]
                        t0 = time()
                        vals = try
                            regrid_chunk(dem, srcspace, sys7, layout, ch, config)
                        catch err
                            say("ERROR worker $w chunk $ch: ",
                                sprint(showerror, err, catch_backtrace())[1:min(end, 1500)])
                            Threads.atomic_add!(FAILURES, 1)
                            nothing
                        finally
                            retire_column!(dem.cache, d)
                        end
                        vals === nothing && continue
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

"""
    reportcache(cache, prefetcher, nsrc, complete)

What the tile cache did, and whether it did it lawfully.

Three claims, in decreasing strength:

  - **Every source chunk was loaded at most once.** This is the structural
    property of refcount eviction and the reason it never reloads: a chunk can
    only be freed when no destination chunk that may read it is left, so no
    correct demand for it can follow. If this fails the refcounts are wrong.

  - **Every demand was one the graph predicted.** The adjacency is a
    conservative superset of the true geometry, so a read the graph did not
    carry is a hole in the covering. Such a read is still served — the run does
    not die over an accounting bug — but it is counted here.

  - **Every source chunk was loaded at least once**, reported only for a fresh,
    complete, failure-free run. This is a COVERAGE diagnostic, not a law: the
    graph's 1.35x conservatism, a resumed run, and a failed chunk each leave
    tiles legitimately untouched. The offline simulator could equate loads with
    the tile count because it treated adjacency as demand; production reads an
    unknown subset of it.
"""
function reportcache(cache, prefetcher, nsrc::Int, complete::Bool)
    s = cachestats(cache)
    LASTCACHE[] = s
    if s.policy === :refcount
        say(@sprintf("tile cache: %d loads of %d source chunks, peak %.2f GiB / %d tiles, %d hits, %d joined loads, %d live at end",
            s.loads, nsrc, s.peakbytes / 2^30, s.peaktiles, s.hits, s.waits, s.live))
        check("every source chunk was loaded at most once", isempty(doubleloaded(cache));
            detail = "$(length(doubleloaded(cache))) chunk(s) loaded twice")
        check("every demand was a predicted edge", s.uncredited == 0;
            detail = s.uncredited == 0 ? "" :
                     "$(s.uncredited) uncredited demand(s), chunks " *
                     string(s.uncreditedof))
        check("every destination chunk was retired", s.pinned == 0 && s.live == 0;
            detail = "$(s.pinned) source chunk(s) still pinned, $(s.live) still held")
        complete && say(s.loads == nsrc ?
            "coverage: every one of the $nsrc source chunks was read" :
            "coverage: $(nsrc - s.loads) of $nsrc source chunks were never read " *
            "(the graph is a conservative superset; this is not an error)")
    else
        say(@sprintf("tile cache: LRU, %d loads (%d reloads over %d source chunks), %d hits, %d live at end, %.2f GiB held",
            s.loads, max(0, s.loads - nsrc), nsrc, s.hits, s.live, s.bytes / 2^30))
    end
    pf = prefetchstats(prefetcher)
    pf.depth == 0 && return nothing
    say("prefetch: depth $(pf.depth), $(pf.issued) tile requests issued")
    check("the prefetcher raised nothing", pf.failure === nothing;
        detail = pf.failure === nothing ? "" : sprint(showerror, pf.failure))
    return nothing
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
    println(stamp(), "  copdem_production.jl — GLO-$(config.res) -> IGEO7 level " *
                     "$(config.level), level-$(config.ancestor) chunks, SYNTHETIC elevations")
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
    ids = TileIds(sys, tiles)
    say("tile list: $(length(tiles)) listed tiles of $(DGG.ncells(sys, 0)) " *
        "($(round(100 * length(tiles) / DGG.ncells(sys, 0); digits = 1))%), " *
        "$(length(real)) backed by a real GeoTIFF, " *
        @sprintf("%.3e pixels, %s", Float64(ids.n), secs(time() - t0)))
    for (o, p) in sort!(collect(real))
        say("  real: $(tilestem(sys, DGG.LevelIndex(0, o)))  " *
            "$(round(filesize(p) / 2^10)) KiB")
    end

    builder = TileBuilder(sys, tiles, real, mask; delay = config.fetchdelay)
    config.fetchdelay > 0 && say(@sprintf(
        "fetchdelay: every tile build sleeps %.3f s — a latency shim, not production",
        config.fetchdelay))
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

    # --- the dependency graph, the walk order, the cache ------------------
    W, shape = workercount(config)
    plan = dagplan(sys, sys7, tiles, todochunks, srcspace, config)
    say(@sprintf("graph: %d x %d chunks, %d edges (%.2f tiles per chunk), radius %g; dst space %s, graph %s, order %s",
        GR.nsourcechunks(plan.graph), GR.ndestinationchunks(plan.graph), plan.edges,
        plan.edges / max(length(todochunks), 1), plan.radius,
        secs(plan.tspace), secs(plan.tgraph), secs(plan.torder)))
    check("graph destination chunks are the work list",
        GR.ndestinationchunks(plan.graph) == length(todochunks))
    check("graph source chunks are the listed tiles",
        GR.nsourcechunks(plan.graph) == length(tiles))
    check("the walk order is a permutation of the work list",
        length(plan.order) == length(todochunks) && isperm(plan.order))

    cache = if config.cachepolicy === :refcount
        RefCountCache{Vector{Float32}}(length(tiles), length(todochunks),
            d -> GR.sourcesof(plan.graph, d),
            s -> GR.consumerdegree(plan.graph, s),
            k -> buildtile(builder, k);
            permits = config.fetchconc > 0 ? config.fetchconc : W)
    elseif config.cachepolicy === :lru
        config.prefetch == 0 || error("prefetch needs cachepolicy = :refcount: a " *
            "bounded LRU can evict a prefetched tile before it is used, and the " *
            "churn that causes was never measured")
        StripedLRUCache{Vector{Float32}}(k -> buildtile(builder, k);
            slots = config.cache, stripes = config.stripes)
    else
        error("cachepolicy must be :refcount or :lru, got $(repr(config.cachepolicy))")
    end
    dem = TiledDEM(ids, builder, cache)

    # A chunk skipped on resume is retired before anything starts, so its tiles
    # are never credited to work that will not happen, and the run's hot loop
    # never has to ask again whether a chunk is done.
    order = plan.order
    if !isempty(done)
        keep = Int[]
        for d in order                       # `order` holds destination numbers
            if todochunks[d] in done
                config.cachepolicy === :refcount && retire_column!(cache, d)
            else
                push!(keep, d)
            end
        end
        order = keep
    end
    nskipped = length(plan.order) - length(order)
    say("schedule: $(config.schedule) order over $(length(order)) chunks " *
        "($nskipped retired on resume), cache=$(config.cachepolicy), " *
        "taper=$(config.taper)")

    if config.dryrun
        say(@sprintf("dryrun: max tiles per chunk %d, max chunks per tile %d",
            maximum(d -> GR.sourcedegree(plan.graph, d), 1:length(todochunks); init = 0),
            maximum(s -> GR.consumerdegree(plan.graph, s), 1:length(tiles); init = 0)))
        say("dryrun: stopping before any regrid")
        return FAILURES[]
    end

    # --- the run ---------------------------------------------------------
    say("launching W=$W workers, shape=$shape, sized to hold $(config.cores) cores " *
        "(nthreads=$(Threads.nthreads())), batch=$(config.batch), budget=$(config.budget)")
    shape === :inner && Threads.nthreads() < 4W &&
        say("NOTE shape=:inner saturates near W=nthreads/4; " *
            "$(Threads.nthreads()) threads will not fill $W workers")
    sched = GuidedSchedule(length(order), config.taper ? W : 1, config.batch)
    prefetcher = nothing
    t0 = time()
    p = try
        if config.prefetch > 0
            fc = config.fetchconc > 0 ? config.fetchconc : W
            prefetcher = Prefetcher(cache, order, d -> GR.sourcesof(plan.graph, d),
                length(tiles), sched; depth = config.prefetch, concurrency = fc)
            say("prefetch: depth $(config.prefetch) chunks, $fc concurrent fetches")
        end
        runchunks(todochunks, order, sched, dem, srcspace, sys7, store, layout,
            nskipped, donelog, prefetcher, config)
    finally
        stop_prefetch!(prefetcher)
    end
    wall = time() - t0
    say(@sprintf("RUN DONE  %d chunks computed, %d skipped, %.4e cells in %s = %.0f cells/s aggregate",
        p.done, p.skipped, Float64(p.cells), hours(wall), p.cells / max(wall, 1e-9)))
    say(@sprintf("source tiles decoded: %d real, %d synthetic; %.3e of %.3e pixels valid (%.2f%% land)",
        builder.nreal[], builder.nsynthetic[], Float64(builder.nland[]),
        Float64(builder.npixels[]),
        100 * builder.nland[] / max(builder.npixels[], 1)))
    say(@sprintf("destination NaN fraction: %.3f%% of %.4e written cells",
        100 * p.nan / max(p.cells, 1), Float64(p.cells)))
    reportcache(cache, prefetcher, length(tiles), FAILURES[] == 0 && nskipped == 0)
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

# Run only when this file IS the program. `include`d — by the A/B harness in
# `copdem_dag_validate.jl`, or by anything else that wants to call `main` with a
# config of its own — it defines the run and starts nothing.
if abspath(PROGRAM_FILE) == @__FILE__
    exit(main() == 0 ? 0 : 1)
end
