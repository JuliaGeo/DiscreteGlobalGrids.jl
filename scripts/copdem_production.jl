# Copernicus DEM GLO-90 -> IGEO7 level 12, into ONE global ancestor-subzone Zarr
# store, written by W worker tasks in one process.
#
#     julia --project=bench -t 63 --gcthreads=8,1 scripts/copdem_production.jl \
#         workers=21 region=all store=/home/asinghvi17/geo/dggstores/copdem90.zarr
#
# The DEM itself is never downloaded. What IS real is the tile LIST: Copernicus
# ships ~26 450 land tiles of the 64 800 the 1x1-degree lattice has, and the
# holding is exactly that list. Every listed tile's pixels are synthesised from
# the analytic field `SYNTHETIC`, except the handful cached on disk, which decode
# for real. Tiles off the list do not exist — the source is a `PartialGrid` over
# the listed tiles only, so an over-covered open-ocean destination pairs with no
# source chunk at all.
#
# Within a listed tile, ocean pixels are NODATA. A land mask rasterised once from
# Natural Earth coastlines decides, so the run exercises the missing-data
# machinery at scale: conservative weights renormalise over the land fraction of
# every coastal cell, and a destination cell that is all ocean comes back NaN
# inside a column that is otherwise written.
#
# Shape, per `regrid-notes/2026-08-20-copdem-zarr-spike.md` sections 10, 13, 14
# and `2026-08-20-la-choice.md`:
#
#   * SOURCE chunks are DEM tiles: `DGGSpace(PartialGrid(land pixels);
#     chunklevel = 0)`, one chunk per listed tile, taking the windowed
#     `BlockCursor` fast path (section 13.1).
#   * DESTINATION work unit is ONE level-5 ancestor column — 7^7 = 823 543
#     level-12 cells, the measured optimum (la-choice).
#   * WRITE unit is the same column: one Zarr chunk, one file, written by
#     `dggwrite!` from whichever task computed it. Columns are disjoint, so no
#     two tasks ever touch the same file and nothing shared is rewritten.
#   * Columns are queued in contiguous batches in canonical Z7 order, so a
#     worker walks a geographically connected run and its DEM tile cache stays
#     hot; workers pull batches dynamically so the polar columns, which cost ~5x
#     a mid-latitude one, do not set the wall clock.
#
# Configuration is `key=value` in ARGS, or `config=<file>` of the same lines
# (ARGS win). See `DEFAULTS`.

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
# Configuration
# ===========================================================================

const DEFAULTS = Dict{String,String}(
    "res" => "90",                 # 90 for GLO-90, 30 for GLO-30
    "level" => "12",               # IGEO7 output level
    "ancestor" => "5",             # La: one work unit and one Zarr chunk
    "store" => "/home/asinghvi17/geo/dggstores/copdem90-igeo7-l12-synthetic.zarr",
    "tilelist" => "",              # "" -> <data>/CopernicusDEM/tileList-glo90.txt
    "tiles" => "",                 # real-tile directory; "" -> <data>/CopernicusDEM/<res>m
    "real" => "auto",              # "auto" = every GeoTIFF there, "none", or stems
    "land" => "",                  # land shapefile; "" -> <data>/naturalearth/ne_10m_land.shp
    "maskarcsec" => "15",          # land-mask resolution; 0 disables the mask
    "region" => "all",             # "all", or "box:w,e,s,n" joined by ";"
    "workers" => "0",              # W: concurrent worker tasks; 0 = size from `cores`
    "cores" => "24",               # the core budget W is sized to hold
    "shape" => "outer",            # "outer" or "inner"; see `workercount`
    "batch" => "8",                # columns handed out per pull
    "budget" => string(2^30),      # lazy regrid byte budget, per worker
    "cache" => "3072",             # decoded source tiles held across all workers
    "stripes" => "64",             # independent locks over that cache
    "resume" => "true",            # skip columns the done log or the store already has
    "checks" => "false",           # run the synthetic oracle after the columns
    "checkcolumns" => "6",         # columns to verify when `checks`
    "heartbeat" => "300",          # seconds between summary lines
    "columncache" => "",           # "" -> <store>.columns.txt
    "donelog" => "",               # "" -> <store>.done.ndjson
    "maxcolumns" => "0",           # 0 = no limit; a smoke-test knob
    "columns" => "",               # explicit column indices, overriding the covering
    "dryrun" => "false",           # plan and report, compute nothing
)

function parseargs(args)
    cfg = copy(DEFAULTS)
    pairs = Dict{String,String}()
    for a in args
        k, v = split(a, "=", limit = 2)
        pairs[String(k)] = String(v)
    end
    if haskey(pairs, "config")
        for line in eachline(pairs["config"])
            s = strip(first(split(line, "#", limit = 2)))
            isempty(s) && continue
            k, v = split(s, "=", limit = 2)
            cfg[strip(String(k))] = strip(String(v))
        end
    end
    merge!(cfg, pairs)
    return cfg
end

const CONFIG = parseargs(ARGS)

cfg(k) = CONFIG[k]
cfgint(k) = parse(Int, CONFIG[k])
cfgbool(k) = parse(Bool, CONFIG[k])

datadir() = get(ENV, "RASTERDATASOURCES_PATH",
    joinpath(@__DIR__, "..", "bench", "data"))

const RES = cfgint("res")
const LEVEL = cfgint("level")
const ANCESTOR = cfgint("ancestor")
const STORE = cfg("store")
const TILEDIR = isempty(cfg("tiles")) ?
                joinpath(datadir(), "CopernicusDEM", "$(RES)m") : cfg("tiles")
const TILELIST = isempty(cfg("tilelist")) ?
                 joinpath(datadir(), "CopernicusDEM", "tileList-glo$(RES).txt") :
                 cfg("tilelist")
const LANDSHP = isempty(cfg("land")) ?
                joinpath(datadir(), "naturalearth", "ne_10m_land.shp") : cfg("land")
const DONELOG = isempty(cfg("donelog")) ? STORE * ".done.ndjson" : cfg("donelog")
const COLUMNCACHE = isempty(cfg("columncache")) ? STORE * ".columns.txt" :
                    cfg("columncache")

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

# ===========================================================================
# The synthetic field
# ===========================================================================

"""
    SYNTHETIC(lon, lat) -> Float64

The analytic stand-in for a DEM tile, in degrees in and metres out:
`1000 sin(3λ) cos(2φ) + 500 cos(7λ) sin(5φ) + 100`.

Smooth on the scale of a level-$(LEVEL) cell by four orders of magnitude, so a
conservative cell mean reduces to the field's value at the cell centre — that is
the oracle the `checks` pass uses. It stays the oracle over a partly-ocean cell,
because a weighted mean of a nearly-constant field is that constant however the
weights renormalise.
"""
SYNTHETIC(lon, lat) = let λ = deg2rad(lon), φ = deg2rad(lat)
    1000 * sin(3λ) * cos(2φ) + 500 * cos(7λ) * sin(5φ) + 100
end

SYNTHETIC(p::GO.UnitSphericalPoint) = SYNTHETIC(US.GeographicFromUnitSphere()(p)...)

# ===========================================================================
# The land mask
# ===========================================================================

"""
    LandMask(ncols, nrows, bits)

A global land/ocean raster on the plate-carrée lattice, `bits[col, row]` true
over land. Column 1 starts at 180°W and row 1 at 90°N, so a lookup is two
divisions and a bit read.

It exists to make ocean pixels INSIDE a listed tile nodata. Copernicus ships
land tiles, so most listed tiles are mostly land and the mask matters at the
coast and on island tiles — which is exactly where the conservative weights have
to renormalise over a partial cell.
"""
struct LandMask
    ncols::Int
    nrows::Int
    dlon::Float64
    dlat::Float64
    bits::BitMatrix
end

const NOMASK = nothing

@inline function island(m::LandMask, lon::Real, lat::Real)
    c = clamp(floor(Int, (lon + 180) / m.dlon) + 1, 1, m.ncols)
    r = clamp(floor(Int, (90 - lat) / m.dlat) + 1, 1, m.nrows)
    return @inbounds m.bits[c, r]
end

island(::Nothing, lon, lat) = true

"""
    landrings(path) -> Vector{Vector{Tuple{Float64,Float64}}}

Every ring of every polygon in a land shapefile, exterior and interior alike.
Holes need no special treatment: the scanline fill below is even-odd, which
reads an interior ring as the hole it is.
"""
function landrings(path::AbstractString)
    out = Vector{Vector{Tuple{Float64,Float64}}}()
    ArchGDAL.read(path) do ds
        layer = ArchGDAL.getlayer(ds, 0)
        for feat in layer
            _collectrings!(out, ArchGDAL.getgeom(feat))
        end
    end
    return out
end

function _collectrings!(out, g)
    t = GI.geomtrait(g)
    if t isa GI.PolygonTrait
        for r in GI.getring(g)
            push!(out, [(Float64(GI.x(p)), Float64(GI.y(p))) for p in GI.getpoint(r)])
        end
    elseif t isa GI.MultiPolygonTrait
        for p in GI.getgeom(g)
            _collectrings!(out, p)
        end
    end
    return out
end

"""
    rasterize_land(rings, arcsec) -> LandMask

Scanline (even-odd) rasterisation of `rings` onto the global `arcsec` lattice.

Edges are bucketed by the first scanline they reach and swept with an active
list, so the cost is proportional to the number of (edge, row) incidences rather
than to rings x rows. Natural Earth is clipped to [-180, 180], so no edge wraps
the antimeridian and the sweep needs no seam handling.
"""
function rasterize_land(rings, arcsec::Int)
    nrows = (180 * 3600) ÷ arcsec
    ncols = (360 * 3600) ÷ arcsec
    dlat = 180 / nrows
    dlon = 360 / ncols
    X1 = Float64[]; Y1 = Float64[]; X2 = Float64[]; Y2 = Float64[]
    for r in rings, k in 1:(length(r) - 1)
        (xa, ya), (xb, yb) = r[k], r[k + 1]
        ya == yb && continue
        if ya < yb
            push!(X1, xa); push!(Y1, ya); push!(X2, xb); push!(Y2, yb)
        else
            push!(X1, xb); push!(Y1, yb); push!(X2, xa); push!(Y2, ya)
        end
    end
    ne = length(X1)
    rowof(y) = clamp(floor(Int, (90 - y) / dlat) + 1, 1, nrows)
    starts = [rowof(Y2[e]) for e in 1:ne]     # `Y2` is the northern endpoint
    order = sortperm(starts)
    bits = falses(ncols, nrows)
    active = Int[]
    xs = Float64[]
    p = 1
    for j in 1:nrows
        ylat = 90 - (j - 0.5) * dlat
        while p <= ne && starts[order[p]] <= j
            push!(active, order[p]); p += 1
        end
        isempty(active) && continue
        empty!(xs)
        i = 1
        while i <= length(active)
            e = active[i]
            if Y1[e] > ylat                    # wholly north of the scanline now
                active[i] = active[end]; pop!(active); continue
            end
            Y1[e] <= ylat < Y2[e] &&
                push!(xs, X1[e] + (ylat - Y1[e]) * (X2[e] - X1[e]) / (Y2[e] - Y1[e]))
            i += 1
        end
        length(xs) < 2 && continue
        sort!(xs)
        for k in 1:2:(length(xs) - 1)
            c1 = clamp(floor(Int, (xs[k] + 180) / dlon) + 1, 1, ncols)
            c2 = clamp(ceil(Int, (xs[k + 1] + 180) / dlon), 1, ncols)
            c1 <= c2 && (@inbounds bits[c1:c2, j] .= true)
        end
    end
    return LandMask(ncols, nrows, dlon, dlat, bits)
end

# ===========================================================================
# The tile list: which tiles exist at all
# ===========================================================================

"The level-0 tile an AWS stem names, e.g. `Copernicus_DSM_COG_30_S90_00_E000_00_DEM`."
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

"""
    listedtiles(sys, path, regions) -> Vector{Int}

The ordinals of every tile named in the Copernicus tile list, ascending, kept to
`regions` where those are given.

This is the one thing about the run that is real data. Tiles NOT on the list do
not exist: they are absent from the source `PartialGrid`, so a destination cell
over open ocean never pairs with a source chunk and is never computed.
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

inregions(::Nothing, lon, lat) = true
inregions(regions, lon, lat) =
    any(r -> r[1] <= lon <= r[2] && r[3] <= lat <= r[4], regions)

"`box:w,e,s,n` boxes joined by `;`, or `nothing` for the whole globe."
function parseregions(spec::AbstractString)
    spec in ("all", "") && return nothing
    out = NTuple{4,Float64}[]
    for part in split(spec, ";")
        kind, rest = split(strip(part), ":", limit = 2)
        kind == "box" || throw(ArgumentError("region is `all` or `box:w,e,s,n`, not $part"))
        w, e, s, n = parse.(Float64, split(rest, ","))
        push!(out, (w, e, s, n))
    end
    return out
end

# ===========================================================================
# The source: a lazy id vector, a lazy value vector, one chunk per listed tile
# ===========================================================================

"""
    TileIds(complete, starts, widths)

The level-1 (pixel) ids of the listed tiles, as one lazy ascending vector.

`PartialGrid` stores its ids by reference and never materialises them, so the
holding's 3.8e10 pixels cost an offsets table of 26 450 integers and nothing
else. `strictly_increasing` is answered by construction rather than by the O(n)
scan `PartialGrid` would otherwise run — the tiles are sorted and their
descendant ranges are disjoint, so the concatenation ascends.
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

"Which listed tile holds source position `p`, and its offset within it."
@inline function tileat(v::TileIds, p::Int)
    k = searchsortedlast(v.offsets, p - 1)
    return k, p - v.offsets[k]
end

"""
    TiledDEM(sys, ids; realtiles, mask, cachesize, stripes)

Every listed tile's pixels as one `Float32` vector in the SOURCE space's own
position order, with chunk `k` equal to listed tile `k` — so a read is always
tile aligned and the regridder's source chunks are the DEM's own tiles.

Tiles named in `realtiles` decode from their GeoTIFF; every other tile is
synthesised from [`SYNTHETIC`](@ref) at the pixel posts, with ocean posts set to
`NaN32` by `mask`. The two kinds are indistinguishable downstream: `NaN` is the
regridder's own invalid sentinel whatever produced it.

The cache is striped: `stripes` independent LRUs under `stripes` locks, and a
tile is generated OUTSIDE its lock, so twenty-one workers do not serialise on
one mutex for the ~15 ms a 1200x1200 synthesis takes.
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

"""
    synthetic_tile(sys, tile, mask) -> Vector{Float32}

[`SYNTHETIC`](@ref) at every post of `tile`, in the tile's own position order —
raster row `j` north first, column `i` west to east, `i` fastest — with ocean
posts set to `NaN32`.

The posts are the pixel-is-point lattice: column `i` sits at `lon_w + i/ncols`
and row `j` at `lat_s + 1 - j/N`.
"""
function synthetic_tile(sys, tile::DGG.LevelIndex, mask)
    lat_s, lon_w = CD.tilecorner(sys, tile)
    nc = Int(CD.ncols_at(sys, lat_s))
    nrows = Int(CD.lat_intervals(sys))
    out = Vector{Float32}(undef, nc * nrows)
    lons = [lon_w + i / nc for i in 0:(nc - 1)]
    s3 = [sin(3 * deg2rad(l)) for l in lons]
    c7 = [cos(7 * deg2rad(l)) for l in lons]
    nland = 0
    @inbounds for j in 0:(nrows - 1)
        lat = (lat_s + 1) - j / nrows
        φ = deg2rad(lat)
        c2, s5 = cos(2φ), sin(5φ)
        base = j * nc
        for i in 1:nc
            if island(mask, lons[i], lat)
                out[base + i] = 1000 * s3[i] * c2 + 500 * c7[i] * s5 + 100
                nland += 1
            else
                out[base + i] = NaN32
            end
        end
    end
    return out, nland
end

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
    # Built outside the lock: a synthesis is milliseconds and a decode is longer,
    # and neither is worth blocking twenty other workers on. Two tasks racing the
    # same tile both build it; the loser's copy is dropped.
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

"""
    realtiles(sys, dir, spec) -> Dict{Int,String}

The tile ordinals backed by a GeoTIFF in `dir`: `"auto"` for every file found
there, `"none"` for an all-synthetic globe, or a comma-separated list of stems.
A named stem that is not on disk is an error.
"""
function realtiles(sys, dir, spec::AbstractString)
    spec == "none" && return Dict{Int,String}()
    isdir(dir) || (spec == "auto" && return Dict{Int,String}())
    stems = spec == "auto" ?
            [splitext(f)[1] for f in readdir(dir) if endswith(f, ".tif")] :
            String.(split(spec, ","))
    out = Dict{Int,String}()
    for stem in stems
        path = joinpath(dir, stem * ".tif")
        isfile(path) || throw(ArgumentError("no GeoTIFF for $stem at $path"))
        out[Int(stemtile(sys, stem).index)] = path
    end
    return out
end

# ===========================================================================
# The destination columns
# ===========================================================================

"""
    covering_columns(sys7, sys, tiles, la) -> Vector{Int}

The level-`la` columns that cover the listed tiles, ascending.

Each tile's 1x1-degree extent is queried as a `MultiOrderCoverage` AT `la`, so
the answer is whole level-`la` cells and the destination of a work unit is one
complete subtree by construction — which is what the subzone store demands
(`DGGSFormatError(check = :incomplete_subtree)` otherwise) and what makes a
column exactly one Zarr chunk.

Every column here meets at least one listed tile, so the skip-pruning the brief
asks for is structural: a column with no source is never enqueued at all.
"""
function covering_columns(sys7, sys, tiles::Vector{Int}, la::Int; nthreads = 1)
    g = DGG.levelgrid(sys7, la)
    parts = [Set{Int}() for _ in 1:nthreads]
    chunks = [k:nthreads:length(tiles) for k in 1:nthreads]
    Threads.@sync for (w, ks) in enumerate(chunks)
        Threads.@spawn for k in ks
            t = DGG.LevelIndex(0, tiles[k])
            lat, lon = CD.tilecorner(sys, t)
            ex = Extents.Extent(X = (Float64(lon), Float64(lon) + 1),
                Y = (Float64(lat), Float64(lat) + 1))
            set = DGG.query(sys7, DGG.MultiOrderCoverage(ex); level = la)
            for c in DGG.CellVector(set; level = la)
                push!(parts[w], DGG.cellposition(g, c))
            end
        end
    end
    return sort!(collect(union(parts...)))
end

"""
    load_columns(path) / save_columns(path, cols)

The column list cached beside the store, because computing it is ~26 000
coverage queries and it is the same list on every resume.
"""
function load_columns(path)
    isfile(path) || return nothing
    cols = Int[]
    for line in eachline(path)
        s = strip(line)
        (isempty(s) || startswith(s, "#")) && continue
        push!(cols, parse(Int, s))
    end
    return cols
end

function save_columns(path, cols)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "# level-$(ANCESTOR) column indices covering the listed tiles")
        for c in cols
            println(io, c)
        end
    end
    return path
end

# ===========================================================================
# Resume
# ===========================================================================

"""
    donecolumns(logpath, storepath, layer) -> Set{Int}

Which columns are already written, from the append-only done log UNION whatever
chunk files the store itself holds.

Two sources because neither is complete on its own. The log is the record of
what this script did and can be deleted or truncated; the chunk listing is the
record of what Zarr HAS, and Zarr does not store a chunk whose every value is
the fill value — so a column that came out all-`NaN`, which over the ocean-side
of the covering is a normal outcome, leaves no file behind and is
indistinguishable from one nobody computed.

The union is therefore what is skipped: recomputing a written column would be
correct but wasteful, and skipping an unwritten one would leave a hole.
"""
function donecolumns(logpath, storepath, layer)
    fromlog = Set{Int}()
    if isfile(logpath)
        for line in eachline(logpath)
            m = match(r"\"col\":(\d+)", line)
            m === nothing || push!(fromlog, parse(Int, m[1]))
        end
    end
    fromdisk = storecolumns(storepath, layer)
    if fromdisk !== nothing
        only_log = length(setdiff(fromlog, fromdisk))
        only_disk = length(setdiff(fromdisk, fromlog))
        (only_log == 0 && only_disk == 0) ||
            say("resume: $only_log logged columns have no chunk file (all-NaN " *
                "columns are stored as nothing at all), $only_disk chunk files " *
                "have no log line; taking the union")
        return union(fromlog, fromdisk)
    end
    say("resume: no chunk listing available at $storepath, using the done log alone")
    return fromlog
end

"""
    storecolumns(path, layer) -> Set{Int} or nothing

The columns a directory store already has a chunk file for.

Zarr v2 names a chunk by its indices in ZARR order, which is the reverse of
Julia's, and this layout is `(capacity, ncolumns)` in Julia — so a column's
chunk is `"<column-1>.0"` and the leading field is the column. The format is
CHECKED rather than assumed: a listing that does not parse returns `nothing` and
the caller falls back to the log.
"""
function storecolumns(path, layer)
    dir = joinpath(path, layer)
    isdir(dir) || return nothing
    out = Set{Int}()
    for f in readdir(dir)
        startswith(f, ".") && continue
        m = match(r"^(\d+)\.(\d+)$", f)
        m === nothing && return nothing
        parse(Int, m[2]) == 0 || return nothing
        push!(out, parse(Int, m[1]) + 1)
    end
    return out
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
    last::Float64
end

Progress(total) = Progress(ReentrantLock(), 0, 0, 0, 0, total, time(), time())

"""
    regrid_column(dem, srcspace, sys7, layout, col; budget) -> Vector{Float32}

One work unit: the level-12 values of one level-`ANCESTOR` column.

`DGG.subtree` gives the column's complete subtree as a ROOTED `PartialGrid`, so
`DGGSpace` knows the destination is one chunk without scanning the global
level-`ANCESTOR` grid for it. `chunks` is never passed — supplying it defeats
the plan's pairing, which is what prunes the source tiles that do not meet this
column.
"""
function regrid_column(dem, srcspace, sys7, layout, col::Int; budget::Int)
    a = DGG.columncell(layout, col)
    dstgrid = DGG.subtree(sys7, a, LEVEL)
    dstspace = DGG.DGGSpace(dstgrid; chunklevel = ANCESTOR)
    out = GR.regrid(dem; to = dstspace, from = srcspace,
        method = DGG.Conservative(), missingpolicy = DGG.Weighted(0.5),
        lazy = true, budget = budget)
    return Float32.(vec(collect(out)))
end

function heartbeat!(p::Progress, force = false)
    hb = cfgint("heartbeat")
    lock(p.lock) do
        now = time()
        (force || now - p.last >= hb) || return nothing
        p.last = now
        el = now - p.started
        done = p.done + p.skipped
        rate = p.cells / max(el, 1e-9)
        left = p.total - done
        eta = done > 0 ? el * left / done : NaN
        say(@sprintf("HEARTBEAT  %d/%d columns (%d skipped) | %.3e cells | %.0f cells/s | elapsed %s | ETA %s | RSS %.1f GiB | NaN %.2f%%",
            done, p.total, p.skipped, Float64(p.cells), rate, hours(el), hours(eta),
            rssgib(), 100 * p.nan / max(p.cells, 1)))
        return nothing
    end
end

"""
    workercount() -> (W, shape)

The worker count `W`, sized from the `cores` budget rather than from a thread
count, and the parallel `shape` it assumes. `workers` overrides `W` when set.

A work unit is one column, and how many cores one worker can hold depends
entirely on which of the two nested parallelisms is allowed to run:

  - `shape = "outer"` sets [`GR.OUTER_PARALLEL`](@ref) in the worker body. The
    weight build stays serial and the lazy plan's source-chunk wave is the only
    parallelism inside a unit. Measured: **1.06 cores per worker**, flat in `W`
    while `W` is well below `nthreads`, and **65 k cells per core-second** — the
    most work per core of the two, but it needs `W ≈ cores` workers, and each
    worker holds a column's weights, so resident memory scales with `W`.

  - `shape = "inner"` leaves `OUTER_PARALLEL` unset, so a narrow wave stands
    aside and ConservativeRegridding threads the weight build itself. Measured
    at 12 threads: one worker reaches 7.0–7.3 cores, the pool saturates at
    **`W ≈ nthreads/4`** holding ~78 % of it, and throughput is **60 k cells per
    core-second** — about 8 % less work per core, for a third of the workers and
    two thirds of the resident memory.

So at a fixed core budget `outer` converts cores into cells slightly faster,
while `inner` reaches the same cores with far fewer workers, and is the only
shape that keeps scaling once `cores` exceeds what `W ≈ cores` workers can be
given threads for. `outer` is the default because the standing budget is 24
cores, where the throughput edge wins; switch to `inner` for a larger budget or
when resident memory is the binding constraint.
"""
function workercount()
    named = cfgint("workers")
    shape = cfg("shape")
    shape in ("outer", "inner") ||
        error("shape must be \"outer\" or \"inner\", got $(repr(shape))")
    named > 0 && return (named, shape)
    c = cfgint("cores")
    c > 0 || error("cores must be positive when workers is not set")
    # Cores per worker, measured; see the docstring.
    W = shape == "outer" ? cld(c, 1.06) : cld(c, 3.1)
    return (max(1, Int(W)), shape)
end

"""
    runcolumns(cols, dem, srcspace, sys7, store, layout, done) -> Progress

`W` worker tasks pulling contiguous batches of columns off one atomic cursor.

Under `shape = "outer"` each worker body sets `GR.OUTER_PARALLEL`, so the nested
weight builds do not spawn their own tasks and oversubscribe the box; the lazy
plan's source-chunk wave is then the only parallelism inside a unit. Under
`shape = "inner"` the scoped value is left alone, and `_wavesize` weighs that
wave against the threading inside the build — standing the wave down whenever
one source chunk carries the column, which is almost always.
"""
function runcolumns(cols, dem, srcspace, sys7, store, layout, done)
    W, shape = workercount()
    outer = shape == "outer"
    B = cfgint("batch")
    budget = cfgint("budget")
    p = Progress(length(cols))
    cursor = Threads.Atomic{Int}(1)
    donelock = ReentrantLock()
    donio = open(DONELOG, "a")
    try
        # `false` is the scoped value's own default, so binding it is the same
        # as leaving it unset — the worker is the outermost scope either way.
        Threads.@sync for w in 1:W
            Threads.@spawn @with GR.OUTER_PARALLEL => outer begin
                while true
                    i = Threads.atomic_add!(cursor, B)
                    i > length(cols) && break
                    for col in cols[i:min(i + B - 1, length(cols))]
                        if col in done
                            lock(p.lock) do; p.skipped += 1; end
                            continue
                        end
                        t0 = time()
                        vals = try
                            regrid_column(dem, srcspace, sys7, layout, col; budget)
                        catch err
                            say("ERROR worker $w column $col: ",
                                sprint(showerror, err, catch_backtrace())[1:min(end, 1500)])
                            FAILURES[] += 1
                            continue
                        end
                        DGG.dggwrite!(store, col, vals)
                        el = time() - t0
                        nnan = count(isnan, vals)
                        lock(p.lock) do
                            p.done += 1
                            p.cells += length(vals)
                            p.nan += nnan
                        end
                        lock(donelock) do
                            println(donio, @sprintf("{\"col\":%d,\"cells\":%d,\"nan\":%d,\"secs\":%.2f,\"w\":%d,\"t\":\"%s\"}",
                                col, length(vals), nnan, el, w, stamp()))
                            flush(donio)
                        end
                        say(@sprintf("w%02d col %d  %d cells  %.1f s  %.0f cells/s  %.1f%% NaN",
                            w, col, length(vals), el, length(vals) / el,
                            100 * nnan / length(vals)))
                        heartbeat!(p)
                    end
                end
            end
        end
    finally
        close(donio)
    end
    heartbeat!(p, true)
    return p
end

# ===========================================================================
# The synthetic oracle
# ===========================================================================

const BOUNDARY_SAMPLES = 8

"""
    SourceMask(mask, g0, listed)

The predicate the oracle actually needs: whether a point has a VALID SOURCE
PIXEL under it, which takes both of the run's two independent kinds of absence.

A pixel is valid when its tile is on the Copernicus list AND the land mask calls
the post land. Either alone is the wrong test. Natural Earth calls a point land
that Copernicus never listed a tile for — a lake island, a coastline the two
datasets disagree about, or, in a regional run, a tile the region filter cut —
and Copernicus lists tiles whose interior is mostly sea. A cell over either is
`NaN` and correctly so.
"""
struct SourceMask{M,G}
    mask::M
    g0::G
    listed::Set{Int}
end

@inline function hassource(s::SourceMask, p::GO.UnitSphericalPoint)
    lon, lat = US.GeographicFromUnitSphere()(p)
    island(s.mask, lon, lat) || return false
    t = DGG.cellat(s.g0, p)
    return t !== nothing && Int(t.index) in s.listed
end

"""
    cellsource(g12, c, sm) -> (:full, :none, :mixed)

Whether a destination cell is fed everywhere, nowhere, or partly, sampled at its
centroid and around its densified boundary.

`:mixed` is the case the run exists to exercise, and the case the analytic
oracle is only approximately right about: those cells are checked against the
field's own range across the cell rather than against its value at the centre,
because a renormalised mean over the valid fraction is a mean over SOME subset
of the cell and the field is not quite constant on one.
"""
function cellsource(g12, c, sm::SourceMask)
    has = hassource(sm, DGG.cell_centroid(g12, c))
    ring = DGG.cell_boundary(g12, c)
    for j in eachindex(ring)
        a, b = ring[j], ring[mod1(j + 1, length(ring))]
        for s in 0:(BOUNDARY_SAMPLES - 1)
            hassource(sm, US.slerp(a, b, s / BOUNDARY_SAMPLES)) == has || return :mixed
        end
    end
    return has ? :full : :none
end

"""
    verify(store, sys7, layout, cols, sm, real, dem, written) -> nothing

Read written columns back and hold them against the analytic field.

Five claims, in the order they can fail:

  * every value is either finite or `NaN` — nothing else exists;
  * a cell with no valid source pixel anywhere under it is `NaN`, and a cell
    fed everywhere by synthetic tiles is `SYNTHETIC` at its centre to within the
    field's own curvature across the cell;
  * a cell straddling the coast — fed over part of itself — still lies inside
    the field's range across that cell, which is what says the conservative
    weights renormalised over the valid fraction rather than diluting toward
    zero or blanking the cell;
  * a REAL tile's own cells are finite, the S90 pole tiles included;
  * a column nobody wrote reads back all `NaN` at zero bytes on disk.
"""
function verify(store, sys7, layout, cols, sm, real, dem, written)
    g12 = DGG.levelgrid(sys7, LEVEL)
    g0 = DGG.levelgrid(dem.sys, 0)
    n = min(cfgint("checkcolumns"), length(cols))
    picked = cols[round.(Int, range(1, length(cols); length = n))]
    say("verify: reading back $n of $(length(cols)) columns")
    anymixed = false
    for col in picked
        a = DGG.columncell(layout, col)
        stack = DGG.dggread(STORE; ancestors = [a])
        vals = collect(stack[:elevation])
        cells = collect(DD.lookup(stack[:elevation], DGG.Cells))
        h = DGG.columnlength(layout, col)
        check("column $col: read back $h values",
            length(vals) == h && length(cells) == h;
            detail = "$(length(vals)) values, $(length(cells)) cells")

        nfin = count(isfinite, vals)
        nnan = count(isnan, vals)
        check("column $col: every value is finite or NaN", nfin + nnan == length(vals);
            detail = "$nfin finite, $nnan NaN")

        # A sparse sample: the full boundary walk is 8 slerps per edge and a
        # column is 823 543 cells, so the oracle takes a stride rather than
        # every cell. The stride is prime to 7 so it does not land on one
        # subtree's worth of siblings.
        stride = max(1, length(vals) ÷ 400)
        nfull = nnone = nmixed = 0
        errs = Float64[]
        bracketfail = 0
        nanfail = Int[]
        finfail = Int[]
        for k in 1:stride:length(vals)
            c = cells[k]
            # Only cells fed exclusively by synthetic tiles have an oracle.
            Int(DGG.cellat(g0, DGG.cell_centroid(g12, c)).index) in keys(real) && continue
            kind = cellsource(g12, c, sm)
            if kind === :none
                nnone += 1
                isnan(vals[k]) || push!(nanfail, k)
            elseif kind === :full
                nfull += 1
                isfinite(vals[k]) || (push!(finfail, k); continue)
                push!(errs, abs(vals[k] - SYNTHETIC(DGG.cell_centroid(g12, c))))
            else
                nmixed += 1
                isfinite(vals[k]) || continue
                samples = SYNTHETIC.(DGG.cell_boundary(g12, c))
                push!(samples, SYNTHETIC(DGG.cell_centroid(g12, c)))
                lo, hi = extrema(samples)
                tol = 1e-3 + 0.05 * (hi - lo)
                (lo - tol <= vals[k] <= hi + tol) || (bracketfail += 1)
            end
        end
        nmixed > 0 && (anymixed = true)
        say("column $col sampled: $nfull fully fed / $nnone unfed / $nmixed coastal")
        check("column $col: unfed cells are NaN", isempty(nanfail);
            detail = "$(length(nanfail)) of $nnone unfed cells are not NaN")
        check("column $col: fully-fed cells are finite", isempty(finfail);
            detail = "$(length(finfail)) of $nfull fully-fed cells are NaN")
        if !isempty(errs)
            mx = maximum(errs)
            check("column $col: fully-fed cells match the analytic field", mx < 1.0;
                detail = @sprintf("max %.3e m, RMS %.3e m over %d cells", mx,
                    sqrt(Statistics.mean(abs2, errs)), length(errs)))
        end
        check("column $col: coastal cells bracket the field", bracketfail == 0;
            detail = "$bracketfail of $nmixed coastal cells outside their own field range")
        if nnan > 0 && nfin > 0
            say("column $col holds BOTH land and nodata: $nfin finite, $nnan NaN " *
                @sprintf("(%.1f%%)", 100 * nnan / length(vals)))
        end
    end
    # A global run must meet the coast; a region need not — an all-Antarctica
    # box has 100 % valid pixels and nothing to renormalise, and saying so is
    # not the same as failing to exercise the path.
    if cfg("region") == "all"
        check("some sampled cell straddles the coast", anymixed;
            detail = anymixed ? "the renormalising path ran" : "no coastal cell was sampled")
    else
        say("coastal cells sampled in this region: " *
            (anymixed ? "yes, the renormalising path ran" :
             "none — a region with no coastline has nothing to renormalise"))
    end

    # The real tiles, and the pole in particular. A GLO tile keeps its full
    # degree of longitude but that degree pinches to nothing at 90 degrees, so
    # the S90 row is 360 slivers meeting at a point — the one place the source
    # geometry is degenerate, and the one worth naming a check after.
    for (o, _) in sort!(collect(real))
        t = DGG.LevelIndex(0, o)
        lat, lon = CD.tilecorner(dem.sys, t)
        c = DGG.cellat(g12, US.UnitSphereFromGeographic()((lon + 0.5, lat + 0.5)))
        c === nothing && continue
        col = DGG.columnindex(layout, DGG.ancestor(sys7, c, ANCESTOR))
        col in written || continue
        stack = DGG.dggread(STORE; ancestors = [DGG.columncell(layout, col)])
        vals = collect(stack[:elevation])
        k = findfirst(==(c), collect(DD.lookup(stack[:elevation], DGG.Cells)))
        check("real tile $(tilestem(dem.sys, t)) is finite at its centre",
            k !== nothing && isfinite(vals[k]);
            detail = k === nothing ? "its centre cell is not in column $col" :
                     @sprintf("%.2f m in column %d", vals[k], col))
    end

    # A column nobody wrote: the first level-`ANCESTOR` column that is not in the
    # covering, which by construction meets no listed tile.
    inset = Set(cols)
    empty = findfirst(i -> !(i in inset), 1:DGG.ncells(sys7, ANCESTOR))
    if empty !== nothing
        a = DGG.columncell(layout, empty)
        vals = collect(DGG.dggread(STORE; ancestors = [a])[:elevation])
        check("unwritten ocean column $empty reads back NaN", all(isnan, vals);
            detail = "$(count(isnan, vals)) of $(length(vals)) NaN")
    end
    return nothing
end

# ===========================================================================
# Main
# ===========================================================================

function main()
    println("="^92)
    println(stamp(), "  copdem_production.jl — GLO-$RES -> IGEO7 level $LEVEL, " *
                     "level-$ANCESTOR ancestor columns, SYNTHETIC data")
    println("  julia $(VERSION)  threads=$(Threads.nthreads())  gc=$(Threads.ngcthreads())  pid=$(getpid())")
    for k in sort!(collect(keys(CONFIG)))
        print("  $k=$(CONFIG[k])")
    end
    println()
    println("="^92)
    flush(stdout)

    sys = DGG.CopernicusDEMSystem(RES)
    sys7 = DGG.IGeo7System()
    capacity = 7^(LEVEL - ANCESTOR)

    # --- the land mask ---------------------------------------------------
    arcsec = cfgint("maskarcsec")
    mask = if arcsec > 0
        t0 = time()
        rings = landrings(LANDSHP)
        m = rasterize_land(rings, arcsec)
        frac = count(m.bits) / length(m.bits)
        say(@sprintf("land mask: %d x %d at %d arcsec from %d rings / %d vertices, %.2f%% land by lattice cell, %.0f MiB, %s",
            m.ncols, m.nrows, arcsec, length(rings), sum(length, rings),
            100 * frac, sizeof(m.bits.chunks) / 2^20, secs(time() - t0)))
        m
    else
        say("land mask: DISABLED (maskarcsec=0), every pixel of a listed tile is valid")
        NOMASK
    end

    # --- the source ------------------------------------------------------
    regions = parseregions(cfg("region"))
    t0 = time()
    tiles = listedtiles(sys, TILELIST, regions)
    isempty(tiles) && error("the tile list and region select no tiles")
    real = realtiles(sys, TILEDIR, cfg("real"))
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
        cachesize = cfgint("cache"), stripes = cfgint("stripes"))
    t0 = time()
    srcgrid = DGG.PartialGrid(sys, 1, ids)
    srcspace = DGG.DGGSpace(srcgrid; chunklevel = 0)
    say("source space: $(GR.nchunks(srcspace)) chunks over " *
        @sprintf("%.3e pixels, built in %s", Float64(DGG.ncells(srcgrid)), secs(time() - t0)))
    check("source chunks are the listed tiles", GR.nchunks(srcspace) == length(tiles);
        detail = "$(GR.nchunks(srcspace)) chunks, $(length(tiles)) tiles")

    # --- the destination columns ----------------------------------------
    cols = load_columns(COLUMNCACHE)
    if cols === nothing
        t0 = time()
        cols = covering_columns(sys7, sys, tiles, ANCESTOR;
            nthreads = max(1, Threads.nthreads() - 1))
        save_columns(COLUMNCACHE, cols)
        say("columns: $(length(cols)) level-$ANCESTOR columns cover the tiles, " *
            "computed in $(secs(time() - t0)), cached at $COLUMNCACHE")
    else
        say("columns: $(length(cols)) read from $COLUMNCACHE")
    end
    if !isempty(cfg("columns"))
        cols = sort!(unique!(parse.(Int, split(cfg("columns"), ","))))
        say("columns: overridden to $(length(cols)) named by `columns=`")
    end
    lim = cfgint("maxcolumns")
    0 < lim < length(cols) &&
        (cols = cols[round.(Int, range(1, length(cols); length = lim))];
         say("columns: limited to $lim by maxcolumns"))
    ncell = length(cols) * capacity
    say(@sprintf("work: %d columns x %d cells = %.4e level-%d cells, %.1f GiB dense f32",
        length(cols), capacity, Float64(ncell), LEVEL, 4 * ncell / 2^30))

    # --- the store -------------------------------------------------------
    fresh = !isdir(STORE)
    if fresh
        mkpath(dirname(STORE))
        t0 = time()
        store = DGG.subzonestore(STORE, sys7, LEVEL; ancestor_level = ANCESTOR,
            layers = ("elevation" => Float32,), capacity = capacity,
            fill_value = NaN, ancestor_coordinate = true,
            attrs = Dict{String,Any}("title" => "Copernicus DEM GLO-$RES (SYNTHETIC) on IGEO7 level $LEVEL",
                "source" => "synthetic analytic field over the real Copernicus GLO-$RES tile list",
                "created" => stamp()))
        say("store: created $STORE, $(DGG.ncells(sys7, ANCESTOR)) columns of " *
            "$capacity, $(secs(time() - t0))")
    else
        store = DGG.subzonestore(STORE)
        say("store: reopened $STORE")
    end
    layout = store.layout
    check("store layout is level $LEVEL over level-$ANCESTOR columns",
        DGG.level(layout) == LEVEL && layout.ancestor_level == ANCESTOR &&
        layout.capacity == capacity;
        detail = "level $(DGG.level(layout)), ancestor $(layout.ancestor_level), " *
                 "capacity $(layout.capacity)")

    done = cfgbool("resume") ? donecolumns(DONELOG, STORE, "elevation") : Set{Int}()
    todo = count(c -> !(c in done), cols)
    say("resume: $(length(done)) columns already written, $todo of $(length(cols)) to do")

    if cfgbool("dryrun")
        say("dryrun: stopping before any regrid")
        return FAILURES[]
    end

    # --- the run ---------------------------------------------------------
    W, shape = workercount()
    say("launching W=$W workers, shape=$shape, sized to hold $(cfg("cores")) cores " *
        "(nthreads=$(Threads.nthreads())), batch=$(cfg("batch")), budget=$(cfg("budget"))")
    shape == "inner" && Threads.nthreads() < 4W &&
        say("NOTE shape=inner saturates near W=nthreads/4; " *
            "$(Threads.nthreads()) threads will not fill $W workers")
    t0 = time()
    p = runcolumns(cols, dem, srcspace, sys7, store, layout, done)
    wall = time() - t0
    say(@sprintf("RUN DONE  %d columns computed, %d skipped, %.4e cells in %s = %.0f cells/s aggregate",
        p.done, p.skipped, Float64(p.cells), hours(wall), p.cells / max(wall, 1e-9)))
    say(@sprintf("source tiles decoded: %d real, %d synthetic; %.3e of %.3e pixels valid (%.2f%% land)",
        dem.nreal[], dem.nsynthetic[], Float64(dem.nland[]), Float64(dem.npixels[]),
        100 * dem.nland[] / max(dem.npixels[], 1)))
    say(@sprintf("destination NaN fraction: %.3f%% of %.4e written cells",
        100 * p.nan / max(p.cells, 1), Float64(p.cells)))
    say(@sprintf("peak RSS %.2f GiB", rssgib()))

    if cfgbool("checks")
        sm = SourceMask(mask, DGG.levelgrid(sys, 0), Set(tiles))
        verify(store, sys7, layout, cols, sm, real, dem,
            donecolumns(DONELOG, STORE, "elevation"))
    end

    say("total wall $(hours(time() - STARTED[])), " *
        (FAILURES[] == 0 ? "NO FAILURES" : "$(FAILURES[]) FAILURE(S)"))
    return FAILURES[]
end

exit(main() == 0 ? 0 : 1)
