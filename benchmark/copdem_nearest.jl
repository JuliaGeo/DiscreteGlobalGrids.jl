# What `NearestCell` costs on real Copernicus DEM tiles regridded to IGeo7.
#
#     julia -t auto --project=benchmark benchmark/copdem_nearest.jl \
#         res=30 box=10,12,46,48 method=nearest
#
# `benchmark/point_tile_baseline.jl` prices the two nearest routes on a toy
# raster: one 750x500 destination tile, a source whose `cellat` is two binary
# searches, four candidate chunks, and everything resident. This file prices the
# same method on the workload the toy stands in for — real 1x1-degree Copernicus
# DEM GeoTIFFs into an IGeo7 grid — where `cellat` is a hierarchy descent, the
# destination is thousands of tiles rather than one, the source chunks are files
# on disk, and a destination cell is finer than a source pixel.
#
# The source is built exactly as `scripts/copdem_production.jl` builds it, by
# including that file: one `TiledDEM` over a `PartialGrid` of the listed tiles'
# pixels, chunk `k` being tile `k`, decoded through `TileBuilder` and held in a
# `StripedLRUCache`. Nothing here reimplements the source.
#
# ARMS. Each is a phase of one lazy regrid, and each phase is timed on a plan in
# a known state rather than by subtraction:
#
#   1. `destination`   building the IGeo7 cells the box covers at `level`, the
#                      `PartialGrid` over them and the `DGGSpace` that tiles it.
#                      Not part of a regrid, but a caller pays it.
#   2. `source`        the same for the Copernicus side: tile ids, `PartialGrid`,
#                      `DGGSpace` at chunk level 0 (one chunk per tile).
#   3. `plan`          `plan_regrid`, which builds the chunk dependency relation
#                      and reads no data.
#   4. `first tile`    on a FRESH plan, one read of destination tile 1: the
#                      first `TileWeights` build, the source loads its manifest
#                      names, and the application. This is the per-tile cold cost
#                      with no cache of any kind warm.
#   5. `cold sweep`    on a second fresh plan, every destination tile read once
#                      in order. Weights are budgeted, so a tile evicted here is
#                      rebuilt if a later tile needs it; the reported cache
#                      counters say whether that happened.
#   6. `warm sweep`    the same sweep again against the plan arm 5 left warm.
#
# `method=pairnearest` swaps in `PairNearest`, which is `NearestCell`'s own
# `buildweights!` with no `sampler`, so a plan around it builds one
# `(destination tile, source chunk)` pair at a time. Nothing but the routing
# differs, so the two methods in one session price the routing alone on real
# data. `method=conservative` runs `Conservative()` over the same arms for scale.
#
# WHAT IS COUNTED, not asserted:
#
#   - source chunks read and how many times each, from a counting `DiskArrays`
#     wrapper that records every `readblock!` range;
#   - GeoTIFFs decoded and pixels decoded, from `TileBuilder`'s own atomics, so
#     a tile decoded twice is visible as a cache miss;
#   - source residency and loads/hits/skips from `residency(::LazyRegridArray)`;
#   - resident weight bytes from the plan's own `PerChunk` accounting, with the
#     manifest lengths of the cached tiles beside them;
#   - peak RSS from `Sys.maxrss()`, which is a high-water mark for the whole
#     process and therefore includes the warm-up.
#
# CORRECTNESS. `check=N` verifies N random destination cells against the source
# directly: the sample site is the destination cell centroid, `cellat` on the
# source space names the pixel, and the regridded value must be BIT-IDENTICAL to
# that pixel (`isequal`, so a NaN pixel must give a NaN cell). A destination
# cell whose site the source covers nowhere must be missing. A nearest stencil
# is one entry of exactly 1.0 and `Weighted` divides by exactly one, so bit
# identity is the right assertion here and approximate equality would be a
# weaker test. `check` is skipped for `conservative`, which averages.
#
# WARM-UP AND STATISTIC. Every arm is a one-shot measurement of a phase that
# takes seconds to minutes, so there is no min-of-n. Compilation is removed
# instead: the process first runs one complete miniature regrid over a
# 0.02-degree box at the same level with the same method, which instantiates
# every type the measured arms use. That warm-up's cost is reported so it can be
# seen not to be in the arms. `provenance()` stamps the run with the Julia
# version, thread and CPU counts and the machine's power mode; the same code
# clocks differently under low power, so compare runs only through a
# same-session ratio.
#
# CONFIGURATION. All `key=value`, all optional, order-free:
#
#   res=90|30          Copernicus product. GLO-90 is 1200 rows a degree, GLO-30
#                      3600. Default 90.
#   box=w,e,s,n        whole degrees, east and north EXCLUSIVE, so `10,11,46,47`
#                      is the single tile N46 E010. Default 10,11,46,47.
#   level=L            IGeo7 destination level. Default: the bare-system rule,
#                      `levelfor(IGeo7System(), source space)` — the level
#                      `to = IGeo7System()` would pick.
#   chunk=A            destination chunk (tile) level. Default `level - 7`, so a
#                      full tile is 7^7 = 823,543 cells, the production chunk.
#   method=nearest|pairnearest|conservative
#   budget=B           lazy byte budget. Default 2^31, the API default.
#   slots=N            source tiles held resident. Default: all of them.
#   data=DIR           directory of `Copernicus_DSM_COG_*_DEM.tif` files.
#                      Default `$COPDEM_TILE_CACHE`, else
#                      `$RASTERDATASOURCES_PATH/CopernicusDEM/<res>m`, else
#                      `bench/data/CopernicusDEM/tiles`. Tiles are NOT
#                      downloaded here; a missing one is an error naming its URL.
#   check=N            destination cells to verify. Default 1000; 0 disables.
#   profile=K          profile a cold read and a warm read of the first K
#                      destination tiles and print the top 25 frames by SELF
#                      time. 0 (the default) disables. The two listings are on
#                      a plan of their own, after the arms.
#   sweeps=0           run arms 1-4 only, skipping both full sweeps.
#
# No data is committed and no path here is absolute.

include(joinpath(@__DIR__, "..", "scripts", "copdem_production.jl"))

using Printf
import Random
import DiskArrays
import Profile

import GlobalRegridding: AbstractRegriddingMethod, WeightCOO, NearestCell

# ---------------------------------------------------------------------------
# Nearest weights on the chunk-pair route
# ---------------------------------------------------------------------------

"""
    PairNearest()

`NearestCell`'s weights on the chunk-pair route: it reports `Points()` and
forwards `buildweights!`, but supplies no `sampler`, so a plan around it builds
one `(destination tile, source chunk)` pair at a time. `placed` counts the
destination cells its builds were asked to place, one `cellat` query each.

The same wrapper as `benchmark/point_tile_baseline.jl`'s, so the two files
measure the same two routes.
"""
struct PairNearest <: AbstractRegriddingMethod
    placed::Threads.Atomic{Int}
end

PairNearest() = PairNearest(Threads.Atomic{Int}(0))

GR.outputsampling(::PairNearest) = DD.Lookups.Points()

GR.supportradius(::PairNearest, space::GR.RegridSpace) =
    GR.supportradius(NearestCell(), space)

function GR.buildweights!(coo::WeightCOO, method::PairNearest, dst::GR.RegridSpace,
    dst_inds, src::GR.RegridSpace, src_inds)
    Threads.atomic_add!(method.placed, length(dst_inds))
    return GR.buildweights!(coo, NearestCell(), dst, dst_inds, src, src_inds)
end

# ---------------------------------------------------------------------------
# A source that records what is read from it
# ---------------------------------------------------------------------------

"""
    CountingTiles(parent)

The production `TiledDEM`, wrapped so that every `readblock!` is recorded
against the source chunk — the Copernicus tile — it lands in. `reads` therefore
gives both how many tiles a phase touched and how many times each was fetched,
independently of how many times the GeoTIFF behind it was decoded.
"""
struct CountingTiles{P} <: DiskArrays.AbstractDiskArray{Float32,1}
    parent::P
    reads::Dict{Int,Int}
    pixels::Threads.Atomic{Int}
    lock::ReentrantLock
end

CountingTiles(parent) =
    CountingTiles(parent, Dict{Int,Int}(), Threads.Atomic{Int}(0), ReentrantLock())

Base.size(A::CountingTiles) = size(A.parent)
DiskArrays.eachchunk(A::CountingTiles) = DiskArrays.eachchunk(A.parent)
DiskArrays.haschunks(::CountingTiles) = DiskArrays.Chunked()

function DiskArrays.readblock!(A::CountingTiles, out, r::AbstractUnitRange)
    k, _ = tileat(A.parent.ids, first(r))
    lock(A.lock) do
        A.reads[k] = get(A.reads, k, 0) + 1
    end
    Threads.atomic_add!(A.pixels, length(r))
    return DiskArrays.readblock!(A.parent, out, r)
end

resetreads!(A::CountingTiles) = lock(A.lock) do
    empty!(A.reads)
    A.pixels[] = 0
    return A
end

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

"Parse `key=value` arguments into a NamedTuple over the defaults."
function options(args = ARGS)
    kv = Dict{String,String}()
    for a in args
        i = findfirst('=', a)
        i === nothing && error("arguments are key=value; got $a")
        kv[a[1:(i - 1)]] = a[(i + 1):end]
    end
    known = ("res", "box", "level", "chunk", "method", "budget", "slots", "data",
        "check", "profile", "sweeps")
    for k in keys(kv)
        k in known || error("unknown option $k; known: $(join(known, ", "))")
    end
    get1(k, d) = parse(Int, get(kv, k, string(d)))
    box = Tuple(parse.(Int, split(get(kv, "box", "10,11,46,47"), ',')))
    length(box) == 4 || error("box is w,e,s,n in whole degrees")
    res = get1("res", 90)
    res in (30, 90) || error("res is 90 or 30")
    method = get(kv, "method", "nearest")
    method in ("nearest", "pairnearest", "conservative") ||
        error("method is nearest, pairnearest or conservative")
    return (; res, box, method,
        level = haskey(kv, "level") ? parse(Int, kv["level"]) : nothing,
        chunk = haskey(kv, "chunk") ? parse(Int, kv["chunk"]) : nothing,
        budget = get1("budget", 2^31),
        slots = haskey(kv, "slots") ? parse(Int, kv["slots"]) : nothing,
        data = get(kv, "data", ""),
        check = get1("check", 1000),
        profile = get1("profile", 0),
        sweeps = get1("sweeps", 1) != 0)
end

"The directory the tiles are read from, by the same precedence the header states."
function tiledir(opt)
    isempty(opt.data) || return abspath(opt.data)
    env = get(ENV, "COPDEM_TILE_CACHE", "")
    isempty(env) || return abspath(env)
    root = get(ENV, "RASTERDATASOURCES_PATH", "")
    isempty(root) ||
        return abspath(joinpath(root, "CopernicusDEM", string(opt.res) * "m"))
    return abspath(joinpath(@__DIR__, "..", "bench", "data", "CopernicusDEM", "tiles"))
end

buildmethod(name) = name == "nearest" ? NearestCell() :
                    name == "pairnearest" ? PairNearest() : DGG.Conservative()

# `LazyCopernicusTiles`' own URL shape, for the error message on a missing tile.
tileurl(res::Int, stem::AbstractString) =
    string(res == 90 ? "https://copernicus-dem-90m.s3.amazonaws.com" :
           "https://copernicus-dem-30m.s3.amazonaws.com", "/", stem, "/", stem, ".tif")

# ---------------------------------------------------------------------------
# The two spaces
# ---------------------------------------------------------------------------

"""
    sourceside(opt, dir) -> NamedTuple

The Copernicus side, built the way `scripts/copdem_production.jl` builds it: one
lazy `TiledDEM` whose chunk `k` is tile `k`, over a `PartialGrid` of those tiles'
pixels, chunked at level 0 so a source chunk is exactly a tile.
"""
function sourceside(opt, dir)
    w, e, s, n = opt.box
    sys = DGG.CopernicusDEMSystem(opt.res)
    tiles = sort!([Int(CD.tilecell(sys, lat, lon).index)
                   for lat in s:(n - 1) for lon in w:(e - 1)])
    paths = Dict{Int,String}()
    for t in tiles
        stem = tilestem(sys, DGG.LevelIndex(0, t))
        path = joinpath(dir, stem * ".tif")
        isfile(path) || error("no GeoTIFF for $stem at $path; fetch it from " *
                              tileurl(opt.res, stem))
        paths[t] = path
    end
    builder = TileBuilder(sys, tiles, paths, nothing, NOMASK)
    slots = opt.slots === nothing ? length(tiles) : opt.slots
    cache = StripedLRUCache{Vector{Float32}}(k -> buildtile(builder, k);
        slots = slots, stripes = 1)
    ids = TileIds(sys, tiles)
    dem = CountingTiles(TiledDEM(ids, builder, cache))
    grid = DGG.PartialGrid(sys, 1, ids)
    space = DGG.DGGSpace(grid; chunklevel = 0)
    bytes = sum(filesize(paths[t]) for t in tiles)
    return (; sys, tiles, paths, builder, cache, dem, grid, space, slots, bytes)
end

"""
    destinationside(sys7, box, level, chunklevel) -> NamedTuple

The IGeo7 cells the box covers at `level`, as a `PartialGrid` tiled at
`chunklevel`. The cover is the cells `MultiOrderCoverage` names, so the
destination is the box and nothing more; taking whole ancestor subtrees the way
a subzone store must would add cells outside the box that no source can reach,
and those would dilute every per-cell number here.
"""
function destinationside(sys7, box, level::Int, chunklevel::Int)
    w, e, s, n = box
    ex = Extents.Extent(X = (Float64(w), Float64(e)), Y = (Float64(s), Float64(n)))
    set = DGG.query(sys7, DGG.MultiOrderCoverage(ex); level = level)
    cells = DGG.CellVector(set; level = level)
    grid = DGG.PartialGrid(cells)
    space = DGG.DGGSpace(grid; chunklevel = chunklevel)
    return (; grid, space)
end

# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------

"""
    phase(label, f) -> (; time, bytes, gctime, value)

One shot of `f`, reported and returned. There is no min-of-n: every arm here
runs for seconds at least, and compilation has been removed by the warm-up
rather than by repetition.
"""
function phase(label, f)
    GC.gc()
    r = @timed f()
    @printf("  %-26s %10.3f s  %14s B alloc  %6.1f%% gc\n",
        label, r.time, commas(r.bytes), 100 * r.gctime / max(r.time, eps()))
    flush(stdout)
    return (; time = r.time, bytes = r.bytes, gctime = r.gctime, value = r.value)
end

commas(n::Integer) = replace(string(n), r"(?<=[0-9])(?=(?:[0-9]{3})+$)" => ",")
commas(x) = string(x)

"""
    provenance() -> NamedTuple

The machine state the timings belong to. Power mode is part of it: the same code
on the same machine runs at different clocks under low power, so a number here is
comparable only against another carrying the same stamp.
"""
function provenance()
    shell(cmd) = try
        strip(read(cmd, String))
    catch
        ""
    end
    power = Sys.isapple() ?
            let s = shell(`pmset -g`)
                m = match(r"\b(low)?powermode\s+(\S+)", s)
                m === nothing ? "unreported" :
                string(m.captures[1] === nothing ? "powermode" : "lowpowermode",
                    " ", m.captures[2])
            end : "unreported"
    return (; julia = string(VERSION), threads = Threads.nthreads(),
        gcthreads = Threads.ngcthreads(),
        ncpu = Sys.isapple() ? shell(`sysctl -n hw.ncpu`) : string(Sys.CPU_THREADS),
        power)
end

"Read every destination tile once, in order, into `out`."
function sweep!(out, A, space)
    for c in 1:Int(GR.nchunks(space))
        r = GR.ownedindices(space, c)
        isempty(r) && continue
        out[r] .= A[r]
    end
    return out
end

"Resident weight accounting, on whichever route the plan took."
function weightstate(plan)
    st = plan.storage
    tiles = hasproperty(st, :tiles) ? collect(values(st.tiles)) : []
    blocks = hasproperty(st, :blocks) ? collect(values(st.blocks)) : []
    manifests = [length(t.sourcechunks) for t in tiles]
    return (; accounted = GR.storagebytes(st), budget = GR.weightbudget(plan.budget),
        ntiles = length(tiles), nblocks = length(blocks),
        manifest_min = isempty(manifests) ? 0 : minimum(manifests),
        manifest_max = isempty(manifests) ? 0 : maximum(manifests),
        manifest_sum = sum(manifests; init = 0))
end

# ---------------------------------------------------------------------------
# Profiling: a flat listing by SELF time
# ---------------------------------------------------------------------------

"""
    flatself(io, label; top = 25)

The top `top` frames by SELF time in the current profile buffer: the frame a
sample was taken IN, not the frames it was called from. `Profile.print`'s flat
format counts a frame in every sample it appears anywhere in, which for a deep
regrid call stack puts the entry points on top and says nothing about where the
time went.

Only samples whose backtrace reaches this package are counted. A profile taken
under `-t auto` samples every thread, and the threads with no work to do park in
`__psynch_cvwait`, which on a single-threaded read is more than nine samples in
ten and says nothing about the read. Both totals are printed, so the fraction
discarded is visible rather than assumed.
"""
function flatself(io, label; top = 25)
    data = Profile.fetch(include_meta = false)
    lidict = Profile.getdict(data)
    # The package root as it is on disk, so a checkout under any directory name
    # is recognised; matching the string "DiscreteGlobalGrids" would miss a
    # worktree that is not called that.
    root = normpath(joinpath(@__DIR__, ".."))
    ours(f) = let path = string(f.file)
        startswith(path, root) || occursin("GlobalRegridding", path)
    end
    counts = Dict{String,Int}()
    total = working = 0
    block = UInt64[]
    function flush!()
        isempty(block) && return nothing
        total += 1
        mine = false
        for ip in block
            for f in get(lidict, ip, ())
                if ours(f)
                    mine = true
                    break
                end
            end
            mine && break
        end
        if mine
            working += 1
            frames = get(lidict, block[1], nothing)
            if frames !== nothing && !isempty(frames)
                f = frames[1]
                key = string(f.func, " at ", basename(string(f.file)), ":", f.line)
                counts[key] = get(counts, key, 0) + 1
            end
        end
        empty!(block)
        return nothing
    end
    for ip in data
        ip == 0 ? flush!() : push!(block, ip)
    end
    flush!()
    @printf(io, "\n%s\n", label)
    @printf(io, "  %d samples reached this package, of %d over all threads (%.1f%%)\n",
        working, total, 100working / max(total, 1))
    working == 0 && return nothing
    for (k, v) in first(sort!(collect(counts); by = last, rev = true), top)
        @printf(io, "  %6.2f%%  %7d  %s\n", 100v / working, v, k)
    end
    # Self time names the leaf, which on a regrid is usually arithmetic in Base.
    # This second table says which of OUR functions those leaves were reached
    # through: a frame is credited once per sample it appears anywhere in, so
    # the numbers are inclusive and nest rather than sum to 100%.
    incl = Dict{String,Int}()
    seen = Set{String}()
    block = UInt64[]
    function tally!()
        isempty(block) && return nothing
        empty!(seen)
        for ip in block, f in get(lidict, ip, ())
            ours(f) || continue
            push!(seen, string(f.func, " at ", basename(string(f.file)), ":", f.line))
        end
        isempty(seen) || for k in seen
            incl[k] = get(incl, k, 0) + 1
        end
        empty!(block)
        return nothing
    end
    for ip in data
        ip == 0 ? tally!() : push!(block, ip)
    end
    tally!()
    println(io, "  -- inclusive, this package's frames only --")
    for (k, v) in first(sort!(collect(incl); by = last, rev = true), top)
        @printf(io, "  %6.2f%%  %7d  %s\n", 100v / working, v, k)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Correctness on real data
# ---------------------------------------------------------------------------

"""
    verify(out, dst, src, dem, n; seed) -> NamedTuple

Check `n` random destination cells against the source read directly. The sample
site is the destination cell's centroid, the source cell is `cellat` of that
site, and the regridded value must equal the source pixel BIT for bit. Where
`cellat` finds nothing the cell must be missing.
"""
function verify(out, dst, src, dem, n::Int; seed = 20260825)
    n <= 0 && return (; checked = 0, mapped = 0, outside = 0, bad = 0, firstbad = nothing)
    rng = Random.MersenneTwister(seed)
    idx = rand(rng, 1:length(out), n)
    mapped = outside = bad = 0
    firstbad = nothing
    for i in idx
        p = GR.cellcentroid(dst, i)
        j = GR.cellat(src, p)
        if j === nothing
            outside += 1
            ok = !isfinite(out[i])
        else
            mapped += 1
            ok = isequal(out[i], dem[j])
        end
        if !ok
            bad += 1
            firstbad === nothing && (firstbad = (; i, j, got = out[i],
                want = j === nothing ? missing : dem[j]))
        end
    end
    return (; checked = n, mapped, outside, bad, firstbad)
end

# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------

function run(opt = options())
    dir = tiledir(opt)
    sys7 = DGG.IGeo7System()
    println(provenance())
    @printf("\nres=GLO-%d  box=(%d,%d,%d,%d)  method=%s  budget=%s B  tiles from %s\n",
        opt.res, opt.box..., opt.method, commas(opt.budget), dir)

    src = sourceside(opt, dir)
    level = opt.level === nothing ? DGG.levelfor(sys7, src.space) : opt.level
    chunklevel = opt.chunk === nothing ? max(first(DGG.levels(sys7)), level - 7) : opt.chunk
    method = buildmethod(opt.method)

    @printf("source   %d tiles, %s pixels, %s B of GeoTIFF, %d cache slots\n",
        length(src.tiles), commas(DGG.ncells(src.grid)), commas(src.bytes), src.slots)
    @printf("         pixel size %.2f m (median over the box)\n",
        DGG.cellsize(src.space))
    @printf("dest     IGeo7 level %d%s, chunk level %d, %.2f m cells\n",
        level, opt.level === nothing ? " (bare-system rule)" : " (override)",
        chunklevel, DGG.cellsize(sys7, level))
    @printf("         IGeo7 level %d holds %s cells globally\n",
        level, commas(DGG.ncells(sys7, level)))
    flush(stdout)

    # --- warm-up: one complete miniature regrid, to take compilation out ----
    # On a source of its own, so no tile it decodes is resident for any arm and
    # the arms below start with an empty tile cache.
    warm = @timed begin
        w = sourceside(opt, dir)
        ex = Extents.Extent(X = (Float64(opt.box[1]) + 0.40, Float64(opt.box[1]) + 0.42),
            Y = (Float64(opt.box[3]) + 0.40, Float64(opt.box[3]) + 0.42))
        set = DGG.query(sys7, DGG.MultiOrderCoverage(ex); level = level)
        g = DGG.PartialGrid(DGG.CellVector(set; level = level))
        sp = DGG.DGGSpace(g; chunklevel = chunklevel)
        p = GR.plan_regrid(w.dem; to = sp, from = w.space, method = method,
            missingpolicy = DGG.Weighted(0.5), lazy = true, budget = opt.budget)
        collect(GR.regrid(w.dem, p))
    end
    @printf("\nwarm-up  %.3f s over %s cells (compilation; not in any arm below)\n",
        warm.time, commas(length(warm.value)))
    flush(stdout)
    @printf("         the GeoTIFFs it touched stay in the OS page cache; every arm\n")
    @printf("         below starts with an EMPTY decoded-tile cache of its own\n")
    method isa PairNearest && (method.placed[] = 0)

    println("\narms")
    a1 = phase("destination", () -> destinationside(sys7, opt.box, level, chunklevel))
    dst = a1.value
    ndst = DGG.ncells(dst.grid)
    a2 = phase("source", () -> DGG.DGGSpace(src.grid; chunklevel = 0))
    a3 = phase("plan", () -> GR.plan_regrid(src.dem; to = dst.space, from = src.space,
        method = method, missingpolicy = DGG.Weighted(0.5), lazy = true,
        budget = opt.budget))
    plan = a3.value
    graph = GR.dependencies(plan)
    nch = Int(GR.nchunks(dst.space))
    cand = [length(GR.sourcesof(graph, c)) for c in 1:nch]
    @printf("\ndestination %s cells in %s tiles; candidates per tile min %d, median %.1f, max %d\n",
        commas(ndst), commas(nch), minimum(cand), Statistics.median(cand), maximum(cand))
    flush(stdout)

    # arm 4: one tile on a fresh plan over a source whose tile cache is empty
    src1 = sourceside(opt, dir)
    fresh = GR.plan_regrid(src1.dem; to = dst.space, from = src1.space, method = method,
        missingpolicy = DGG.Weighted(0.5), lazy = true, budget = opt.budget)
    A1 = GR.regrid(src1.dem, fresh)
    tile1 = GR.ownedindices(dst.space, 1)
    a4 = phase("first tile", () -> A1[tile1])
    @printf("  tile 1 is %s cells; %d source chunk reads, %s pixels, %d GeoTIFF decodes\n",
        commas(length(tile1)), sum(values(src1.dem.reads); init = 0),
        commas(src1.dem.pixels[]), src1.builder.nreal[])
    method isa PairNearest && (method.placed[] = 0)
    A1 = nothing
    fresh = nothing
    src1 = nothing

    out = nothing
    a5 = a6 = nothing
    if opt.sweeps
        out = Vector{Float32}(undef, ndst)
        A = GR.regrid(src.dem, plan)
        a5 = phase("cold sweep", () -> sweep!(out, A, dst.space))
        coldreads = copy(src.dem.reads)
        coldpixels = src.dem.pixels[]
        coldstats = deepcopy(GR.residency(A))
        colddecodes = src.builder.nreal[]
        coldweights = weightstate(plan)
        coldplaced = method isa PairNearest ? method.placed[] : 0
        resetreads!(src.dem)
        a6 = phase("warm sweep", () -> sweep!(out, A, dst.space))

        nfinite = count(isfinite, out)
        @printf("\nvalues     %s finite of %s (%.1f%%)",
            commas(nfinite), commas(ndst), 100nfinite / ndst)
        if nfinite > 0
            lo, hi = extrema(Iterators.filter(isfinite, out))
            @printf("  range [%.2f, %.2f] m", lo, hi)
        end
        println()
        @printf("throughput cold %s cells/s, warm %s cells/s\n",
            commas(round(Int, ndst / a5.time)), commas(round(Int, ndst / a6.time)))
        @printf("           cold %s finite cells/s\n",
            commas(round(Int, nfinite / a5.time)))
        @printf("source     %d of %d chunks read cold, %s times, %s pixels, %d GeoTIFF decodes\n",
            length(coldreads), length(src.tiles),
            commas(sum(values(coldreads); init = 0)), commas(coldpixels), colddecodes)
        @printf("           warm sweep read %d chunks, %s times\n",
            length(src.dem.reads), commas(sum(values(src.dem.reads); init = 0)))
        println("           residency cold ", coldstats)
        println("           residency warm ", GR.residency(A))
        @printf("weights    %s B of a %s B budget in %d cached tiles / %d blocks; manifest %d..%d, %s total\n",
            commas(coldweights.accounted), commas(coldweights.budget),
            coldweights.ntiles, coldweights.nblocks, coldweights.manifest_min,
            coldweights.manifest_max, commas(coldweights.manifest_sum))
        if method isa PairNearest
            @printf("locations  %s destination cells placed by the pair route on the cold sweep, %s on the warm\n",
                commas(coldplaced), commas(method.placed[] - coldplaced))
        end
    end

    if out !== nothing && opt.check > 0 && !(method isa DGG.Conservative)
        v = verify(out, dst.space, src.space, src.dem, opt.check)
        @printf("\ncheck      %s cells: %s mapped, %s outside coverage, %d wrong\n",
            commas(v.checked), commas(v.mapped), commas(v.outside), v.bad)
        v.firstbad === nothing || println("           first mismatch ", v.firstbad)
    end

    if opt.profile > 0
        k = min(opt.profile, nch)
        rs = [GR.ownedindices(dst.space, c) for c in 1:k]
        # Every source tile is decoded and resident before either listing, so
        # neither is a profile of GDAL. What the arms above cost in decoding is
        # reported there, in seconds and in decode counts.
        for j in 1:length(src.tiles)
            gettile!(src.cache, j)
        end
        resetreads!(src.dem)
        pplan = GR.plan_regrid(src.dem; to = dst.space, from = src.space,
            method = method, missingpolicy = DGG.Weighted(0.5), lazy = true,
            budget = opt.budget)
        AP = GR.regrid(src.dem, pplan)
        Profile.init(n = 40_000_000, delay = 0.002)
        Profile.clear()
        Profile.@profile for r in rs
            AP[r]
        end
        flatself(stdout,
            "cold read of $k destination tiles ($(commas(sum(length, rs))) cells): weights built here")
        Profile.clear()
        Profile.@profile for r in rs
            AP[r]
        end
        flatself(stdout,
            "warm read of the same $k tiles: weights cached, source loads and application only")
        Profile.clear()
    end

    @printf("\npeak RSS   %.2f GiB\n", Sys.maxrss() / 2^30)
    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    run()
end
