# Copernicus DEM -> IGEO7 level-12 -> chunked Zarr, tree-aligned on both sides.
#
#     julia --project=bench -t 8 scripts/copdem_to_igeo7_zarr.jl \
#         region=cap:-89.9 ancestor=7 method=conservative out=/tmp/spike
#
# The chunking is the point of this script. Nothing here mosaics rasters:
#
#   * SOURCE chunks are Copernicus DEM tiles — the level-0 cells of
#     `CopernicusDEMSystem`. The DEM enters as ONE lazy `DiskArrays` vector over
#     the complete level-1 (pixel) grid whose chunks are exactly the tiles'
#     `descendant_range`s, so a "read" is always a whole 1°x1° tile.
#   * DESTINATION chunks are IGEO7 ancestor subtrees at level `ancestor`. The
#     region is queried as a `MultiOrderCoverage` AT THAT LEVEL and only then
#     expanded to `level`, so every cell of the destination belongs to a
#     COMPLETE ancestor subtree and every ancestor run is exactly
#     `7^(level - ancestor)` cells long.
#   * WRITE chunks are the same runs: `dggwrite(...; chunk_target = 7^(level -
#     ancestor))` makes `:auto` pick that length, and the store's manifest
#     marker records `ancestor_level` and `ancestor_aligned` so the alignment is
#     checked rather than assumed.
#
# Only the tiles named by `real=` are read from disk; every other tile is filled
# with the analytic field `SYNTHETIC`, which makes the whole globe available
# without downloading it and gives the output an exact oracle.
#
# Configuration is `key=value` in ARGS; see `CONFIG` below for the keys.

import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import GeometryOps as GO
import DimensionalData as DD
import DiskArrays
import ArchGDAL
import Extents
import Zarr
import Statistics
import Printf: @sprintf

const CD = DGG.CopernicusDEM
const US = GO.UnitSpherical

# ===========================================================================
# Configuration
# ===========================================================================

const DEFAULTS = Dict{String,String}(
    "res" => "90",                 # 90 for GLO-90, 30 for GLO-30
    "level" => "12",               # IGEO7 output level
    "ancestor" => "7",             # IGEO7 ancestor level: chunk = one subtree
    "method" => "conservative",    # conservative | bilinear | both
    "region" => "cap:-89.9",       # cap:<max latitude> | box:<w>,<e>,<s>,<n>
    "tiles" => "",                 # real-tile directory; "" uses `datadir()`
    "real" => "auto",              # "auto" = every tile found there, or stems
    "out" => "",                   # output directory; "" uses `datadir()/spike-out`
    "budget" => string(2^30),      # lazy regrid byte budget
    "cache" => "128",              # decoded source tiles held in the LRU
    "checks" => "true",            # run the validation pass
    "pairs" => "4",                # destination chunks to log source pairings for
    "merge" => "step",             # id-range merge rule: step | rank
)

parseargs(args) = merge(DEFAULTS, Dict(String(first(p)) => String(last(p))
                                       for p in split.(args, "=", limit = 2)))

const CONFIG = parseargs(ARGS)

cfg(k) = CONFIG[k]
cfgint(k) = parse(Int, CONFIG[k])
cfgbool(k) = parse(Bool, CONFIG[k])

datadir() = get(ENV, "RASTERDATASOURCES_PATH",
    joinpath(@__DIR__, "..", "bench", "data"))

const RES = cfgint("res")
const LEVEL = cfgint("level")
const ANCESTOR = cfgint("ancestor")
const SUBTREE = 7^(LEVEL - ANCESTOR)          # level-`LEVEL` cells per chunk
const TILEDIR = isempty(cfg("tiles")) ?
                joinpath(datadir(), "CopernicusDEM", "$(RES)m") : cfg("tiles")
const OUTDIR = isempty(cfg("out")) ? joinpath(datadir(), "spike-out") : cfg("out")

# ===========================================================================
# Reporting
# ===========================================================================

const FAILURES = Ref(0)
const STARTED = Ref(time())

function check(name, ok; detail = "")
    ok || (FAILURES[] += 1)
    println(ok ? "PASS  " : "FAIL  ", rpad(name, 54), detail)
    flush(stdout)
    return ok
end

note(text) = (println("      ", text); flush(stdout))
section(text) = (println(); println("--- ", text, " ", "-"^max(0, 62 - length(text)));
                 flush(stdout))

rss() = @sprintf("%.2f GiB", Sys.maxrss() / 2^30)
secs(t) = @sprintf("%.1f s", t)

# ===========================================================================
# The synthetic field
# ===========================================================================

"""
    SYNTHETIC(lon, lat) -> Float64

The analytic stand-in for a DEM tile, in degrees in and metres out:
`1000 sin(3λ) cos(2φ) + 500 cos(7λ) sin(5φ) + 100`, with λ and φ in radians.

Smooth on the scale of a level-$(LEVEL) cell by four orders of magnitude, so a
conservative cell mean and a bilinear sample both reduce to the field's value at
the cell centre — that is the oracle every synthetic check below uses.
"""
SYNTHETIC(lon, lat) = let λ = deg2rad(lon), φ = deg2rad(lat)
    1000 * sin(3λ) * cos(2φ) + 500 * cos(7λ) * sin(5φ) + 100
end

SYNTHETIC(p::GO.UnitSphericalPoint) = SYNTHETIC(US.GeographicFromUnitSphere()(p)...)

# ===========================================================================
# The source: one lazy vector over every GLO pixel, one chunk per tile
# ===========================================================================

"""
    TiledDEM(sys; realtiles, cachesize)

Every Copernicus DEM pixel as one `Float32` vector in level-1 cell order, with
chunk `k` equal to tile `k`'s `descendant_range` — so a read is always tile
aligned and the regridder's source chunks ARE the DEM's own tiles.

`realtiles` maps a tile ordinal to the GeoTIFF that holds it; every other tile
is generated from [`SYNTHETIC`](@ref) at the pixel posts. Decoded tiles are held
in a `cachesize`-entry LRU behind a lock, because the lazy regridder reads source
chunks from several tasks at once.
"""
struct TiledDEM{S<:DGG.CopernicusDEMSystem,C} <: DiskArrays.AbstractDiskArray{Float32,1}
    sys::S
    chunks::C
    realtiles::Dict{Int,String}
    cache::Dict{Int,Vector{Float32}}
    order::Vector{Int}                  # LRU order, oldest first
    lock::ReentrantLock
    loaded::Vector{Int}                 # every decode, in decode order
    nreal::Base.RefValue{Int}
    nsynthetic::Base.RefValue{Int}
    cachesize::Int
end

function TiledDEM(sys::DGG.CopernicusDEMSystem; realtiles::Dict{Int,String},
    cachesize::Integer = 128)
    widths = [length(DGG.descendant_range(sys, t, 1)) for t in DGG.rootcells(sys)]
    chunks = DiskArrays.GridChunks(DiskArrays.IrregularChunks(; chunksizes = widths))
    return TiledDEM(sys, chunks, realtiles, Dict{Int,Vector{Float32}}(), Int[],
        ReentrantLock(), Int[], Ref(0), Ref(0), Int(cachesize))
end

Base.size(A::TiledDEM) = (DGG.ncells(A.sys, 1),)
DiskArrays.eachchunk(A::TiledDEM) = A.chunks
DiskArrays.haschunks(::TiledDEM) = DiskArrays.Chunked()

"The AWS object stem of `tile`, e.g. `Copernicus_DSM_COG_30_S90_00_E000_00_DEM`."
function tilestem(sys, tile)
    lat, lon = CD.tilecorner(sys, tile)
    tag = lat_intervals_tag(sys)
    return string("Copernicus_DSM_COG_", tag, "_", lat < 0 ? "S" : "N",
        lpad(abs(lat), 2, '0'), "_00_", lon < 0 ? "W" : "E",
        lpad(abs(lon), 3, '0'), "_00_DEM")
end

lat_intervals_tag(sys) = CD.lat_intervals(sys) == 3600 ? "10" : "30"

"""
    synthetic_tile(sys, tile) -> Vector{Float32}

[`SYNTHETIC`](@ref) at every post of `tile`, in the tile's own position order —
raster row `j` north first, column `i` west to east, `i` fastest.

The posts are the pixel-is-point lattice: column `i` sits at `lon_w + i/ncols`
and row `j` at `lat_s + 1 - j/N`, which is `cell_box`'s centre by construction.
"""
function synthetic_tile(sys::DGG.CopernicusDEMSystem, tile::DGG.LevelIndex)
    lat_s, lon_w = CD.tilecorner(sys, tile)
    nc = Int(CD.ncols_at(sys, lat_s))
    nrows = Int(CD.lat_intervals(sys))
    out = Vector{Float32}(undef, nc * nrows)
    lons = [deg2rad(lon_w + i / nc) for i in 0:(nc - 1)]
    s3, c7 = sin.(3 .* lons), cos.(7 .* lons)
    @inbounds for j in 0:(nrows - 1)
        φ = deg2rad((lat_s + 1) - j / nrows)
        c2, s5 = cos(2φ), sin(5φ)
        base = j * nc
        for i in 1:nc
            out[base + i] = 1000 * s3[i] * c2 + 500 * c7[i] * s5 + 100
        end
    end
    return out
end

"The decoded band of the GeoTIFF at `path`, flattened into position order."
readtile(path) = vec(ArchGDAL.read(ds -> ArchGDAL.read(ds, 1), path))

# One tile's values, from the LRU or freshly built. The lock covers the whole
# body: two tasks asking for the same tile must not both decode it.
function tilevalues!(A::TiledDEM, ordinal::Int)
    lock(A.lock) do
        cached = get(A.cache, ordinal, nothing)
        if cached !== nothing
            push!(A.order, splice!(A.order, findfirst(==(ordinal), A.order)))
            return cached
        end
        path = get(A.realtiles, ordinal, nothing)
        v = if path === nothing
            A.nsynthetic[] += 1
            synthetic_tile(A.sys, DGG.LevelIndex(0, ordinal))
        else
            A.nreal[] += 1
            readtile(path)
        end
        push!(A.loaded, ordinal)
        length(A.cache) >= A.cachesize && delete!(A.cache, popfirst!(A.order))
        A.cache[ordinal] = v
        push!(A.order, ordinal)
        return v
    end
end

function DiskArrays.readblock!(A::TiledDEM, out, r::AbstractUnitRange)
    p = first(r)
    while p <= last(r)
        tile = Base.parent(A.sys, DGG.LevelIndex(1, p - 1))
        window = DGG.descendant_range(A.sys, tile, 1)
        stop = min(last(r), last(window))
        seg = (p - first(r) + 1):(stop - first(r) + 1)
        v = tilevalues!(A, Int(tile.index))
        out[seg] .= view(v, (p - first(window) + 1):(stop - first(window) + 1))
        p = stop + 1
    end
    return out
end

"""
    realtiles(sys, dir, spec) -> Dict{Int,String}

The tile ordinals backed by a GeoTIFF in `dir`. `spec` is `"auto"` for every
file found there, `"none"` for an all-synthetic globe, or a comma-separated list
of stems; a named stem that is not on disk is an error, because a spike that
silently synthesises the tile it meant to test against is a spike that tests
nothing.
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

"The level-0 tile an AWS stem names. Inverse of [`tilestem`](@ref)."
function stemtile(sys, stem::AbstractString)
    m = match(r"_([NS])(\d{2})_00_([EW])(\d{3})_00_DEM$", stem)
    m === nothing && throw(ArgumentError("$stem is not a Copernicus DEM stem"))
    lat = parse(Int, m[2]) * (m[1] == "S" ? -1 : 1)
    lon = parse(Int, m[4]) * (m[3] == "W" ? -1 : 1)
    return CD.tilecell(sys, lat, lon)
end

# ===========================================================================
# The region and the destination
# ===========================================================================

"""
    region(spec) -> SphericalCap or Extents.Extent

`cap:<lat>` is the spherical cap around the SOUTH pole reaching up to latitude
`lat` — a cap, not a lat/lon box, so its rim is a small circle rather than four
great-circle edges that bow poleward between their corners.

`box:<w>,<e>,<s>,<n>` is the lon/lat rectangle, for sector and band-edge runs.
"""
function region(spec::AbstractString)
    kind, rest = split(spec, ":", limit = 2)
    if kind == "cap"
        lat = parse(Float64, rest)
        return US.SphericalCap(US.UnitSphericalPoint(0.0, 0.0, -1.0),
            deg2rad(90 + lat))
    elseif kind == "box"
        w, e, s, n = parse.(Float64, split(rest, ","))
        return Extents.Extent(X = (w, e), Y = (s, n))
    end
    throw(ArgumentError("region is cap:<lat> or box:<w>,<e>,<s>,<n>, not $spec"))
end

"""
    destination(sys, target, level, ancestor) -> (grid, set)

The level-`level` cells of the ancestor subtrees that cover `target`.

The query runs AT `ancestor`, so every cell the coverage emits sits at that
level or coarser and expanding it to `level` yields only COMPLETE subtrees. That
is what makes the ancestor runs uniform — `7^(level - ancestor)` cells each —
and therefore what lets one uniform Zarr chunk length land on every one of them.
The cost is over-coverage: the destination is the smallest union of whole
subtrees containing the region, not the region.
"""
function destination(sys, target, level::Int, ancestor::Int)
    set = DGG.query(sys, DGG.MultiOrderCoverage(target); level = ancestor)
    cv = DGG.CellVector(set; level = level)
    return DGG.PartialGrid(cv), set
end

# ===========================================================================
# The conservative path
# ===========================================================================

"""
    conservative_cube(dem, srcspace, dstspace; budget) -> lazy DimArray

The whole regrid as one lazy cube. The plan pairs each destination ancestor
subtree with each source DEM tile whose cap meets it and drops the rest, which
is the pairing this script exists to exercise; `chunks` is deliberately NOT
passed, since supplying it defeats that pruning.
"""
function conservative_cube(dem, srcspace, dstspace; budget)
    cube = DD.DimArray(dem, DGG.Cells(DGG.CellLookup(DGG.levelgrid(dem.sys, 1)));
        name = :elevation)
    out = DGG.regrid(cube; to = dstspace, from = srcspace,
        method = DGG.Conservative(), missingpolicy = DGG.Weighted(0.5),
        lazy = true, budget = budget)
    return DD.rebuild(out; name = :elevation)
end

# ===========================================================================
# The bilinear path
# ===========================================================================
#
# `BilinearPoint` needs `hascellchart(source)`, which only `RasterGrid` answers
# `true` to, so the DGGSpace source of the conservative path cannot be reused.
# The chunk shape is kept instead of the source: walk the SAME destination
# ancestor subtrees, split each one by the DEM tile its cell centres fall in,
# and run one small raster-to-cells regrid per (subtree, tile) pair.

"""
    tileraster(dem, tile; halo = true) -> DimArray

`tile` as a lon/lat raster of its posts, optionally ringed by one pixel of its
neighbours so that a bilinear stencil at the tile border still straddles real
data.

East and west neighbours always share the tile's band and so its column count;
north and south neighbours only sometimes do. Where they do not, the edge row is
replicated and the interpolation degrades to clamping across that border —
`BAND_EDGE_CLAMPS` counts it.
"""
function tileraster(dem::TiledDEM, tile::DGG.LevelIndex; halo::Bool = true)
    sys = dem.sys
    lat_s, lon_w = CD.tilecorner(sys, tile)
    nc = Int(CD.ncols_at(sys, lat_s))
    nrows = Int(CD.lat_intervals(sys))
    dlon, dlat = 1 / nc, 1 / nrows
    core = reshape(tilevalues!(dem, Int(tile.index)), nc, nrows)

    pad = halo ? 1 : 0
    A = Matrix{Float32}(undef, nc + 2pad, nrows + 2pad)
    A[(1 + pad):(nc + pad), (1 + pad):(nrows + pad)] .= core

    if halo
        west = _neighbourtile(sys, lat_s, lon_w - 1)
        east = _neighbourtile(sys, lat_s, lon_w + 1)
        A[1, (1 + pad):(nrows + pad)] .= _column(dem, west, nc, nrows, nc)
        A[nc + 2, (1 + pad):(nrows + pad)] .= _column(dem, east, nc, nrows, 1)
        # Rows: `j = 0` is the NORTH row, so the north halo comes from the tile
        # one degree poleward and the south halo from one degree equatorward.
        A[:, 1] .= _row(dem, sys, lat_s + 1, lon_w, nc, A, pad, nrows, :north)
        A[:, nrows + 2] .= _row(dem, sys, lat_s - 1, lon_w, nc, A, pad, nrows, :south)
    end

    lons = [lon_w + (i - 1 - pad) * dlon for i in 1:(nc + 2pad)]
    lats = [(lat_s + 1) - (j - 1 - pad) * dlat for j in 1:(nrows + 2pad)]
    return DD.DimArray(A, (DD.X(lons), DD.Y(lats)); name = :elevation)
end

_neighbourtile(sys, lat_s, lon_w) =
    CD.tilecell(sys, lat_s, mod(lon_w + 180, 360) - 180)

function _column(dem, tile, nc, nrows, i)
    v = tilevalues!(dem, Int(tile.index))
    lat_s, _ = CD.tilecorner(dem.sys, tile)
    ncn = Int(CD.ncols_at(dem.sys, lat_s))
    ncn == nc || return fill(NaN32, nrows)   # unreachable: E/W share the band
    return [v[(j - 1) * nc + i] for j in 1:nrows]
end

const BAND_EDGE_CLAMPS = Ref(0)

# A north/south halo row, or the replicated edge row where the neighbour tile
# has a different column count and its posts do not line up.
function _row(dem, sys, lat_s, lon_w, nc, A, pad, nrows, side)
    inner = side === :north ? 2 : nrows + 1
    if -90 <= lat_s <= 89 && Int(CD.ncols_at(sys, lat_s)) == nc
        tile = CD.tilecell(sys, lat_s, lon_w)
        v = tilevalues!(dem, Int(tile.index))
        j = side === :north ? nrows - 1 : 0     # neighbour's adjacent post row
        row = Vector{Float32}(undef, nc + 2pad)
        row[(1 + pad):(nc + pad)] .= view(v, (j * nc + 1):(j * nc + nc))
        row[1], row[end] = row[1 + pad], row[nc + pad]
        return row
    end
    BAND_EDGE_CLAMPS[] += 1
    return A[:, inner]
end

"""
    bilinear_chunk(dem, sys7, ids) -> Vector{Float32}

One destination ancestor subtree, assembled from one bilinear regrid per DEM
tile. The cells are partitioned by the tile their CENTRE falls in, so no cell is
ever sampled from a raster it lies outside — which is what keeps `BilinearPoint`
from silently clamp-extrapolating.
"""
function bilinear_chunk(dem::TiledDEM, sys7, ids::AbstractVector)
    g12 = DGG.levelgrid(sys7, LEVEL)
    g0 = DGG.levelgrid(dem.sys, 0)
    groups = Dict{Int,Vector{Int}}()
    for (k, c) in enumerate(ids)
        t = DGG.cellat(g0, DGG.cell_centroid(g12, c))
        push!(get!(() -> Int[], groups, Int(t.index)), k)
    end
    out = Vector{Float32}(undef, length(ids))
    for (ordinal, ks) in groups
        raster = tileraster(dem, DGG.LevelIndex(0, ordinal))
        sub = DGG.PartialGrid(sys7, LEVEL, ids[ks])
        vals = DGG.regrid(raster; to = sub, method = DGG.BilinearPoint(),
            missingpolicy = DGG.Weighted(0.5), lazy = false)
        out[ks] .= Float32.(vec(parent(vals)))
    end
    return out
end

"""
    BilinearCells(dem, sys7, grid, ranges)

The bilinear result as a lazy vector over the destination cell axis whose chunks
are the ancestor runs — the same write unit the conservative path streams, so
`dggwrite` sees one interface for both methods and never materialises the whole
level.
"""
struct BilinearCells{D,S,G,C} <: DiskArrays.AbstractDiskArray{Float32,1}
    dem::D
    sys7::S
    grid::G
    ranges::Vector{UnitRange{Int}}
    chunks::C
end

function BilinearCells(dem, sys7, grid, ranges)
    widths = [length(r) for r in ranges]
    return BilinearCells(dem, sys7, grid, ranges,
        DiskArrays.GridChunks(DiskArrays.IrregularChunks(; chunksizes = widths)))
end

Base.size(A::BilinearCells) = (DGG.ncells(A.grid),)
DiskArrays.eachchunk(A::BilinearCells) = A.chunks
DiskArrays.haschunks(::BilinearCells) = DiskArrays.Chunked()

function DiskArrays.readblock!(A::BilinearCells, out, r::AbstractUnitRange)
    ids = [DGG.cellindex(A.grid, p) for p in r]
    out .= bilinear_chunk(A.dem, A.sys7, ids)
    return out
end

# ===========================================================================
# Writing and reading back
# ===========================================================================

"""
    worker_groups(dstspace, sys, level) -> Vector{Vector{Int}}

The destination chunks grouped by their level-`level` ancestor: the partition an
outer, worker-parallel run would hand out, one group per worker task.

Not wired up here — this spike writes one store from one `dggwrite` call — but
the partition is the mechanical half of that design and belongs beside the
chunking it depends on. Two properties make it the right unit:

  - **Spatially clustered.** IGEO7 ids sort into contiguous subtrees, so the
    chunks under one coarse ancestor are neighbours on the sphere. A worker's
    DEM-tile LRU therefore stays hot: the mid-latitude runs here decoded 2 tiles
    for their first chunk and 0 for the next three.
  - **Disjoint on the write side.** Each group is a contiguous run of the cell
    axis and an ancestor run in its own right, so a group is exactly one store
    shard — no two workers ever touch the same Zarr chunk.

`level` is the knob: coarser gives fewer, larger groups. See the spike record for
the sizing rule (`workers x threads <= CPU threads`, and `workers x per-worker
footprint <= RAM`).
"""
function worker_groups(dstspace, sys, level::Integer)
    groups = Dict{Any,Vector{Int}}()
    order = Any[]
    for k in eachindex(dstspace.chunkids)
        a = DGG.ancestor(sys, dstspace.chunkids[k], Int(level))
        haskey(groups, a) || push!(order, a)
        push!(get!(() -> Int[], groups, a), k)
    end
    return [groups[a] for a in order]
end

"""
    writestore(path, cube) -> path

`chunk_target` is set to one ancestor subtree, so `chunks = :auto` picks exactly
that length: it is the largest whole number of level-`ANCESTOR` runs at or under
the target, and every run here is that long.

`encoding` is left at `:auto`, which resolves to ranges here — the axis is
sorted and unique, which is all `RangesEncoding` asks. `merge = :rank` shrinks
that coordinate by two orders of magnitude on an ancestor-snapped destination
(see the spike record) but is read back only by a rank-aware reader, so `:step`
stays the default. Chunk bytes are Zarr.jl's own Blosc default.
"""
function writestore(path, cube)
    ispath(path) && rm(path; recursive = true)
    mkpath(dirname(path))
    DGG.dggwrite(path, cube; chunks = :auto, chunk_target = SUBTREE,
        merge = Symbol(cfg("merge")))
    return path
end

"The chunk grid and manifest marker a store ended up with."
function storeplan(path, layer = "elevation")
    z = Zarr.zopen(path)
    arr = z[layer]
    marker = nothing
    for name in keys(z.arrays)
        a = z[name]
        haskey(a.attrs, "dggs_chunk_manifest") &&
            (marker = a.attrs["dggs_chunk_manifest"])
    end
    names = sort!(collect(keys(z.arrays)))
    coord = "cell_id_ranges" in names ? "cell_id_ranges" :
            ("cell_ids" in names ? "cell_ids" : "(implicit)")
    rows = coord == "cell_id_ranges" ? size(z[coord][:, :], 2) : 0
    return (chunks = arr.metadata.chunks, shape = arr.metadata.shape[],
        marker = marker, arrays = names, coord = coord, rows = rows,
        encoding = get(get(z.attrs, "dggs", Dict()), "compression", "?"),
        coordbytes = coord == "(implicit)" ? 0 : dirsize(joinpath(path, coord)),
        layerbytes = dirsize(joinpath(path, layer)))
end

dirsize(path) = sum(filesize(joinpath(root, f))
                    for (root, _, files) in walkdir(path) for f in files; init = 0)

"""
    axisnote(plan, ncells) -> nothing

What the id axis cost, which is the point of the `encoding` keyword: the store
names every cell once, and `:auto` picks ranges over a dense id array wherever
the axis is sorted and unique.

The row count is worth reading. `merge = :step` merges only ids ADJACENT AS
INTEGERS, and an IGEO7 `Z7Cell` packs each base-7 digit into four bits, so only
seven of every sixteen integer values name a cell: seven siblings is the longest
integer-contiguous run that exists, and the coordinate is `ncells / 7` rows
however tidy the cell set is. `merge = :rank` merges runs of consecutive CELLS
and collapses an ancestor-snapped destination to a handful of rows, at the price
of needing a rank-aware reader.
"""
function axisnote(plan, ncells)
    dense = 8 * ncells
    note("id axis: encoding $(repr(plan.encoding)), coordinate $(plan.coord), " *
         "$(plan.rows) rows, $(plan.coordbytes) B on disk against $dense B " *
         @sprintf("for a dense UInt64 axis (%.2f%%); layer %d B", 100 * plan.coordbytes / max(dense, 1),
        plan.layerbytes))
    return nothing
end

# ===========================================================================
# Validation
# ===========================================================================

"""
    synthetic_only(g0, g12, cells, real) -> BitVector

Which destination cells draw only on synthetic tiles, where `g0` is the DEM's
level-0 (tile) grid.

Sampling the centroid and the boundary VERTICES is not enough near a pole. A
DEM tile keeps its full degree of longitude but that degree pinches to nothing
at the pole, so a tile is about 190 m wide at 89.9° and 19 m at 89.99° — narrower
than the 61 m level-$(LEVEL) cell it crosses. A real tile can therefore slice
through a cell without containing any of its vertices. Each edge is densified to
`$(BOUNDARY_SAMPLES)` samples, and a cell whose cap reaches the pole — where all
360 tiles of the row meet — is treated as touching every one of them.
"""
const BOUNDARY_SAMPLES = 16

function synthetic_only(g0, g12, cells, real::Dict{Int,String})
    out = trues(length(cells))
    isempty(real) && return out
    poles = (US.UnitSphericalPoint(0.0, 0.0, 1.0), US.UnitSphericalPoint(0.0, 0.0, -1.0))
    sys7 = g12.system
    for (k, c) in enumerate(cells)
        cap = DGG.node_extent(sys7, c)
        if any(p -> US.spherical_distance(cap.point, p) <= cap.radius, poles)
            out[k] = false
            continue
        end
        ring = DGG.cell_boundary(g12, c)
        hit = Int(DGG.cellat(g0, DGG.cell_centroid(g12, c)).index) in keys(real)
        for j in eachindex(ring)
            hit && break
            a, b = ring[j], ring[mod1(j + 1, length(ring))]
            for s in 0:(BOUNDARY_SAMPLES - 1)
                p = US.slerp(a, b, s / BOUNDARY_SAMPLES)
                if Int(DGG.cellat(g0, p).index) in keys(real)
                    hit = true
                    break
                end
            end
        end
        out[k] = !hit
    end
    return out
end

"""
    polar_outliers(name, errs, ks, cells, g12) -> nothing

Where the worst analytic errors sit, measured in cell radii from the nearer
pole. If the tail is a pole artefact this prints single-digit distances; if it
prints large ones the regrid, not the oracle, is wrong.
"""
function polar_outliers(name, errs, ks, cells, g12)
    order = sortperm(errs; rev = true)
    worst = order[1:min(5, length(order))]
    ncap = 0
    for j in eachindex(errs)
        c = cells[ks[j]]
        cap = DGG.node_extent(g12.system, c)
        d = min(US.spherical_distance(cap.point, US.UnitSphericalPoint(0.0, 0.0, 1.0)),
            US.spherical_distance(cap.point, US.UnitSphericalPoint(0.0, 0.0, -1.0)))
        errs[j] > 1.0 && d < 10 * cap.radius && (ncap += 1)
    end
    over = count(>(1.0), errs)
    note("$name: $over cells over 1 m, of which $ncap sit within 10 cell radii " *
         "of a pole")
    for j in worst
        c = cells[ks[j]]
        cap = DGG.node_extent(g12.system, c)
        d = min(US.spherical_distance(cap.point, US.UnitSphericalPoint(0.0, 0.0, 1.0)),
            US.spherical_distance(cap.point, US.UnitSphericalPoint(0.0, 0.0, -1.0)))
        lon, lat = US.GeographicFromUnitSphere()(DGG.cell_centroid(g12, c))
        note(@sprintf("  err %8.2f m at lon %10.5f lat %10.5f, %.1f cell radii from the pole",
            errs[j], lon, lat, d / cap.radius))
    end
    return nothing
end

"""
    validate(name, values, cells, g0, g12, real) -> nothing

Finiteness, range, and — over the cells fed only by synthetic tiles — the
absolute error against [`SYNTHETIC`](@ref) at the cell centre.
"""
function validate(name, values, cells, g0, g12, real)
    nan = count(!isfinite, values)
    check("$name: every cell is finite", nan == 0;
        detail = nan == 0 ? "$(length(values)) cells" : "$nan non-finite")
    if nan > 0
        bad = findall(!isfinite, values)[1:min(5, nan)]
        for k in bad
            lon, lat = US.GeographicFromUnitSphere()(DGG.cell_centroid(g12, cells[k]))
            note("non-finite at position $k, cell $(cells[k]), " *
                 @sprintf("lon %.5f lat %.5f", lon, lat))
        end
    end

    syn = synthetic_only(g0, g12, cells, real)
    nsyn = count(syn)
    check("$name: some cells are synthetic-only", nsyn > 0;
        detail = "$nsyn of $(length(cells)) ($(round(100nsyn / length(cells); digits = 1))%)")
    if nsyn > 0
        ks = [k for k in findall(syn) if isfinite(values[k])]
        # TWO oracles, because one of them is only valid away from a pole.
        #
        # `errs` is the field at the cell CENTRE, which equals the cell's mean
        # only where the field is flat across the cell. It is not, near a pole:
        # a level-$(LEVEL) cell is about 61 m across but a degree of longitude
        # pinches to nothing there, so a cell at 89.9° spans a third of a DEGREE
        # of longitude, over which `SYNTHETIC` swings tens of metres. That is the
        # oracle failing, not the regrid.
        #
        # `bracket` is valid everywhere: a conservative cell mean is an average
        # of the field over the cell, so it must lie between the field's smallest
        # and largest values there. `spread` is how much room that leaves.
        errs, rel, outside = Float64[], Float64[], Int[]
        for k in ks
            c = cells[k]
            samples = push!(SYNTHETIC.(DGG.cell_boundary(g12, c)),
                SYNTHETIC(DGG.cell_centroid(g12, c)))
            lo, hi = extrema(samples)
            spread = hi - lo
            e = abs(values[k] - samples[end])
            push!(errs, e)
            push!(rel, e / max(spread, 1e-12))
            tol = 1e-3 + 0.05 * spread
            (lo - tol <= values[k] <= hi + tol) || push!(outside, k)
        end
        mx, rms = maximum(errs), sqrt(Statistics.mean(abs2, errs))
        spreadmax = maximum(rel)
        check("$name: the cell mean brackets the analytic field", isempty(outside);
            detail = @sprintf("%d of %d cells outside their own cell's field range",
                length(outside), length(ks)))
        note(@sprintf("%s: |value - field(centre)| max %.3e m, RMS %.3e m; worst is %.2f x the field's own spread across that cell",
            name, mx, rms, spreadmax))
        polar_outliers(name, errs, ks, cells, g12)
    end

    finite = filter(isfinite, values)
    isempty(finite) || note("$name: range " *
                            @sprintf("[%.2f, %.2f] m, mean %.2f", minimum(finite),
        maximum(finite), Statistics.mean(finite)))
    return nothing
end

"""
    readback(name, path, cells) -> Vector{Float64}

The store as `dggread` sees it, with its cell axis checked against the axis the
regrid was asked for. Validation runs on THIS vector rather than on the lazy
cube, so the numbers that are checked are the numbers on disk and the regrid is
not recomputed to check itself.
"""
function readback(name, path, cells)
    back = DGG.dggread(path; lazy = true)
    layer = back[:elevation]
    ids = collect(DD.lookup(layer, DGG.Cells))
    check("$name round trip: cell count and ids",
        length(ids) == length(cells) && ids == collect(cells);
        detail = "$(length(ids)) cells")
    return Float64.(collect(layer))
end

"""
    recompute_chunk(name, cube, r, stored) -> nothing

One destination chunk computed a second time and compared against the store.
This is what says the store holds the regrid's output; `readback` only says the
axis survived.
"""
function recompute_chunk(name, cube, r, stored)
    fresh = Float32.(collect(cube[r]))
    same = count(k -> isequal(fresh[k], Float32(stored[first(r) + k - 1])), eachindex(fresh))
    check("$name: a recomputed chunk is bit-identical to the store",
        same == length(fresh);
        detail = "$same of $(length(fresh)) cells in $(first(r)):$(last(r))")
    return nothing
end

# ===========================================================================
# Main
# ===========================================================================

function run_conservative(dem, srcspace, dstgrid, dstspace, cells, g0, g12, real)
    section("conservative")
    t0 = time()
    cube = conservative_cube(dem, srcspace, dstspace; budget = cfgint("budget"))
    note("plan built in $(secs(time() - t0)), result $(size(cube))")

    # Log the (destination subtree x source tile) pairing for a few chunks.
    npairs = cfgint("pairs")
    if npairs > 0
        for k in 1:min(npairs, length(dstspace.ranges))
            r = dstspace.ranges[k]
            before = length(dem.loaded)
            v = cube[r]
            tiles = unique(dem.loaded[(before + 1):end])
            stems = [tilestem(dem.sys, DGG.LevelIndex(0, t)) for t in tiles]
            note("chunk $k (cells $(first(r)):$(last(r)), ancestor $(dstspace.chunkids[k])) " *
                 "newly decoded $(length(tiles)) tiles" *
                 (isempty(stems) ? " (all cached)" : ": " * join(first(stems, 3), ", ") *
                                                     (length(stems) > 3 ? " ..." : "")))
            length(v) == length(r) || error("chunk read returned $(length(v)) values")
        end
    end

    t0 = time()
    path = joinpath(OUTDIR, "out_conservative.zarr")
    writestore(path, cube)
    wall = time() - t0
    plan = storeplan(path)
    note("wrote $path in $(secs(wall)), $(round(dirsize(path) / 2^20; digits = 1)) MiB, " *
         "peak RSS $(rss())")
    check("write chunks are one ancestor subtree",
        plan.chunks == (SUBTREE,);
        detail = "chunks $(plan.chunks), target $SUBTREE, shape $(plan.shape)")
    if plan.marker !== nothing
        # `:auto` reports the COARSEST level whose runs still fit the target, not
        # the level asked for: where each level-`ANCESTOR` subtree of the
        # destination happens to be the only one under its parent, the parent's
        # runs are the same length and the plan reports the parent. That is a
        # stronger statement, not a wrong one, so `<=` is the check.
        al = get(plan.marker, "ancestor_level", nothing)
        check("the manifest records ancestor alignment",
            al isa Integer && al <= ANCESTOR &&
            get(plan.marker, "ancestor_aligned", false) === true &&
            get(plan.marker, "chunk_length", 0) == SUBTREE;
            detail = "ancestor_level=$al aligned=$(get(plan.marker, "ancestor_aligned", nothing)) " *
                     "chunk_length=$(get(plan.marker, "chunk_length", nothing))")
    end
    note("store arrays: " * join(plan.arrays, ", "))
    axisnote(plan, length(cells))

    values = readback("conservative", path, cells)
    validate("conservative", values, cells, g0, g12, real)
    recompute_chunk("conservative", cube, dstspace.ranges[1], values)
    note("source tiles decoded: $(dem.nreal[]) real, $(dem.nsynthetic[]) synthetic")
    note("wall $(secs(wall)) for $(length(values)) cells = " *
         @sprintf("%.1f cells/s", length(values) / wall))
    return values
end

function run_bilinear(dem, sys7, dstgrid, dstspace, cells, g0, g12, real)
    section("bilinear")
    t0 = time()
    lazy = BilinearCells(dem, sys7, dstgrid, dstspace.ranges)
    cube = DD.DimArray(lazy, DGG.Cells(DGG.CellLookup(dstgrid)); name = :elevation)
    path = joinpath(OUTDIR, "out_bilinear.zarr")
    writestore(path, cube)
    wall = time() - t0
    plan = storeplan(path)
    note("wrote $path in $(secs(wall)), $(round(dirsize(path) / 2^20; digits = 1)) MiB, " *
         "peak RSS $(rss()), band-edge clamps $(BAND_EDGE_CLAMPS[])")
    check("write chunks are one ancestor subtree", plan.chunks == (SUBTREE,);
        detail = "chunks $(plan.chunks), target $SUBTREE")
    axisnote(plan, length(cells))

    values = readback("bilinear", path, cells)
    validate("bilinear", values, cells, g0, g12, real)
    return values
end

function main()
    println("="^78)
    println("copdem_to_igeo7_zarr.jl — GLO-$RES -> IGEO7 level $LEVEL, " *
            "chunked on level-$ANCESTOR subtrees")
    println("julia $(VERSION)  threads=$(Threads.nthreads())")
    for k in sort!(collect(keys(CONFIG)))
        print("  $k=$(CONFIG[k])")
    end
    println()
    println("="^78)

    sys = DGG.CopernicusDEMSystem(RES)
    sys7 = DGG.IGeo7System()
    g12 = DGG.levelgrid(sys7, LEVEL)

    section("source")
    real = realtiles(sys, TILEDIR, cfg("real"))
    note("real tiles: $(length(real)) from $TILEDIR")
    for (o, p) in sort!(collect(real))
        note("  $(tilestem(sys, DGG.LevelIndex(0, o)))  " *
             "$(round(filesize(p) / 2^10; digits = 0)) KiB")
    end
    dem = TiledDEM(sys; realtiles = real, cachesize = cfgint("cache"))
    t0 = time()
    srcspace = DGG.DGGSpace(DGG.levelgrid(sys, 1); chunklevel = 0)
    note("source space: $(GR.nchunks(srcspace)) chunks (DEM tiles) over " *
         "$(DGG.ncells(sys, 1)) pixels, built in $(secs(time() - t0))")
    check("source chunks are the DEM's own tiles",
        GR.nchunks(srcspace) == 64_800 &&
        srcspace.ranges[1] == 1:length(DGG.descendant_range(sys, DGG.LevelIndex(0, 0), 1)))

    section("destination")
    target = region(cfg("region"))
    t0 = time()
    dstgrid, set = destination(sys7, target, LEVEL, ANCESTOR)
    cells = DGG.CellVector(dstgrid)
    note("region $(cfg("region")) -> $(length(set)) coverage cells at level " *
         "<= $ANCESTOR, $(DGG.ncells(dstgrid)) level-$LEVEL cells, " *
         "$(secs(time() - t0))")
    t0 = time()
    dstspace = DGG.DGGSpace(dstgrid; chunklevel = ANCESTOR)
    note("destination space: $(GR.nchunks(dstspace)) chunks at level $ANCESTOR, " *
         "built in $(secs(time() - t0))")
    widths = unique(length.(dstspace.ranges))
    check("every ancestor run is a complete subtree", widths == [SUBTREE];
        detail = "run lengths $(widths), 7^($LEVEL-$ANCESTOR) = $SUBTREE")
    note(@sprintf("output %.2f MiB dense f32", 4 * DGG.ncells(dstgrid) / 2^20))
    if ANCESTOR >= 2
        g = worker_groups(dstspace, sys7, ANCESTOR - 2)
        note("outer-parallel partition at level $(ANCESTOR - 2): $(length(g)) groups " *
             "of $(extrema(length.(g))) chunks (not run here; see the spike record)")
    end

    g0 = DGG.levelgrid(sys, 0)
    method = cfg("method")
    cons = method in ("conservative", "both") ?
           run_conservative(dem, srcspace, dstgrid, dstspace, cells, g0, g12, real) :
           nothing
    bilin = method in ("bilinear", "both") ?
            run_bilinear(dem, sys7, dstgrid, dstspace, cells, g0, g12, real) : nothing

    if cons !== nothing && bilin !== nothing
        section("agreement")
        ok = findall(k -> isfinite(cons[k]) && isfinite(bilin[k]), eachindex(cons))
        d = [bilin[k] - cons[k] for k in ok]
        check("bilinear and conservative agree", !isempty(ok);
            detail = @sprintf("RMS %.4e m, max %.4e m, corr %.8f over %d cells",
                sqrt(Statistics.mean(abs2, d)), maximum(abs, d),
                Statistics.cor(cons[ok], bilin[ok]), length(ok)))
    end

    println()
    println("  total wall $(secs(time() - STARTED[])), peak RSS $(rss())")
    println(FAILURES[] == 0 ? "ALL CHECKS PASSED" : "$(FAILURES[]) CHECK(S) FAILED")
    return FAILURES[]
end

exit(main() == 0 ? 0 : 1)
