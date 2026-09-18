# Real Copernicus DEM tiles regridded to IGeo7: what each method costs, and how
# far apart the methods' answers are.
#
#     julia -t auto --project=benchmark benchmark/copdem.jl cost \
#         res=90 box=10,11,46,47 methods=nearest,point,conservative
#     julia -t auto --project=benchmark benchmark/copdem.jl semantics \
#         res=90 box=10,11,46,47 arms=real,linear,quadratic,rim,polar
#
# The source is `scripts/copdem_production.jl`'s: one `TiledDEM` over a
# `PartialGrid` of the listed tiles' pixels, chunk `k` being tile `k`, decoded by
# `CopernicusTiles` and held in a `StripedLRUCache`. Methods and missing policies
# come from that script's `regridmethod` and `regridpolicy`, so what runs here is
# what a production run configured the same way computes. `pairnearest` is the
# one extra name: `NearestCell`'s weights with no `sampler`, which prices the
# chunk-pair route against the sampled one.
#
# COST. One destination, and per method a source, plan and tile cache of its
# own, all in one process so the ratios belong to one machine state:
#
#   destination   the IGeo7 cells the box covers at `level`, as a `PartialGrid`
#                 tiled at `chunk`. Paid once, shared.
#   source        the Copernicus `DGGSpace` at chunk level 0 (one chunk per tile).
#   plan          `plan_regrid`: builds the chunk dependency relation, reads no data.
#   first tile    one read of destination tile 1 on a fresh plan and an empty
#                 tile cache: the per-tile cold cost.
#   stencils      `point` only: one `weightsat!` per destination site with no
#                 source read, against the holding and against the complete level.
#   cold sweep    every destination tile read once, in order.
#   warm sweep    the same sweep against the plan the cold sweep left warm.
#
# Counted, not asserted: source chunk reads and pixels from a counting
# `DiskArrays` wrapper; GeoTIFF decodes from the tile cache's load counter;
# `residency(::LazyRegridArray)`; resident weight bytes from the plan's own
# accounting; peak RSS from `Sys.maxrss()`, a process high-water mark that
# includes the warm-up.
#
# `check=N` verifies N random destination cells per method. A nearest value must
# be BIT-IDENTICAL to the pixel `cellat` names at the cell centroid. A point
# value must be the weighted sum of the posts `weightsat!` names. Where the
# source has nothing, the cell must be missing. Conservative has no pointwise
# oracle and is skipped.
#
# Every phase is a one-shot measurement of seconds of work; compilation is
# removed by a complete miniature regrid per method before any arm.
#
# SEMANTICS. Each arm regrids one source under every method, the first method
# being the baseline the others are compared against:
#
#   real       the listed tiles of `box`, on real elevations.
#   linear     a field affine in lon/lat on the scaled lattice twin, destination
#              finer than the source. Point stencils reproduce it to roundoff;
#              the area mean is off by about half a source pixel of slope.
#   quadratic  a field with curvature on the same twin: no method reproduces it.
#   rim        the twin with a destination reaching past its edge. A point
#              stencil needing an absent post stays unmapped; an area weight exists.
#   polar      the strip between the twin's polemost post row and the pole, run
#              under both settings of the point method's `poles` policy.
#
# Reported per arm and method: cold sweep time, chunk reads, residency, error
# against the analytic field where there is one, and against the baseline the
# fraction both placed, mean and max absolute difference, and cells only one
# placed.
#
# OPTIONS, all `key=value`, order-free:
#
#   res=90|30        Copernicus product. Default 90.
#   box=w,e,s,n      whole degrees, east and north exclusive. Default 10,11,46,47.
#   level=L          IGeo7 destination level. Default `levelfor(IGeo7System(), source)`.
#   chunk=A          destination chunk level. Default `level - 7`.
#   methods=a,b      cost default nearest,point,conservative; semantics default
#                    conservative,point (first is the baseline).
#   budget=B         lazy byte budget. Default 2^31.
#   slots=N          source tiles held resident. Default all of them.
#   data=DIR         directory of `Copernicus_DSM_COG_*_DEM.tif`. Default
#                    `$COPDEM_TILE_CACHE`, else
#                    `$RASTERDATASOURCES_PATH/CopernicusDEM/<res>m`, else
#                    `benchmark/data/CopernicusDEM/tiles`. A missing tile is an
#                    error naming its URL; nothing is downloaded.
#   check=N          cost: cells verified per method. Default 1000; 0 disables.
#   profile=K        cost: flat self-time profile of a cold and a warm read of the
#                    first K destination tiles, per method. Default 0.
#   sweeps=0         cost: skip the sweeps.
#   arms=a,b         semantics: which arms. Default all five.
#   twin=W,E,S,N     semantics: the twin's box. Default 10,12,49,51, which
#                    straddles the 50-degree band edge.
#   fine=K           semantics: destination levels finer than the rule's on the
#                    twin arms. Default 1.

include(joinpath(@__DIR__, "..", "scripts", "copdem_common.jl"))

using Printf
import Random
import DiskArrays
import Extents
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

const SEMANTIC_ARMS = (:real, :linear, :quadratic, :rim, :polar)

"""
    methodpair(name) -> (method, missingpolicy)

The pair `scripts/copdem_production.jl` runs under `method = name`, plus
`pairnearest` for [`PairNearest`](@ref).
"""
function methodpair(name)
    name == "pairnearest" && return (PairNearest(), DGG.Weighted(0.5))
    config = (; method = Symbol(name))
    return (regridmethod(config), regridpolicy(config))
end

"Parse `key=value` arguments into a NamedTuple over `command`'s defaults."
function options(command::Symbol, args)
    kv = Dict{String,String}()
    for a in args
        i = findfirst('=', a)
        i === nothing && error("arguments are key=value; got $a")
        kv[a[1:(i - 1)]] = a[(i + 1):end]
    end
    known = ("res", "box", "level", "chunk", "methods", "budget", "slots", "data",
        "check", "profile", "sweeps", "arms", "twin", "fine")
    for k in keys(kv)
        k in known || error("unknown option $k; known: $(join(known, ", "))")
    end
    int(k, d) = parse(Int, get(kv, k, string(d)))
    maybeint(k) = haskey(kv, k) ? parse(Int, kv[k]) : nothing
    degreebox(k, d) = let b = Tuple(parse.(Int, split(get(kv, k, d), ',')))
        length(b) == 4 || error("$k is w,e,s,n in whole degrees")
        b
    end
    res = int("res", 90)
    res in (30, 90) || error("res is 90 or 30")
    methods = Tuple(String.(split(get(kv, "methods", command === :cost ?
        "nearest,point,conservative" : "conservative,point"), ',')))
    foreach(methodpair, methods)          # reject an unknown name before any work
    command === :semantics && length(methods) < 2 &&
        error("semantics needs a baseline method and at least one other")
    arms = Tuple(Symbol.(split(get(kv, "arms", join(SEMANTIC_ARMS, ',')), ',')))
    for a in arms
        a in SEMANTIC_ARMS || error("unknown arm $a; known: $(join(SEMANTIC_ARMS, ", "))")
    end
    return (; res, methods, arms,
        box = degreebox("box", "10,11,46,47"),
        twin = degreebox("twin", "10,12,49,51"),
        level = maybeint("level"), chunk = maybeint("chunk"), slots = maybeint("slots"),
        budget = int("budget", 2^31), data = get(kv, "data", ""),
        check = int("check", 1000), profile = int("profile", 0),
        sweeps = int("sweeps", 1) != 0, fine = int("fine", 1))
end

"The directory the tiles are read from, by the precedence the header states."
function tiledir(opt)
    isempty(opt.data) || return abspath(opt.data)
    env = get(ENV, "COPDEM_TILE_CACHE", "")
    isempty(env) || return abspath(env)
    root = get(ENV, "RASTERDATASOURCES_PATH", "")
    isempty(root) ||
        return abspath(joinpath(root, "CopernicusDEM", string(opt.res) * "m"))
    return abspath(joinpath(@__DIR__, "data", "CopernicusDEM", "tiles"))
end

"The destination level and chunk level `opt` asks for against `srcspace`."
function destinationlevels(opt, sys7, srcspace)
    level = opt.level === nothing ? DGG.levelfor(sys7, srcspace) : opt.level
    chunk = opt.chunk === nothing ? max(first(DGG.levels(sys7)), level - 7) : opt.chunk
    return level, chunk
end

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
    # Local files only: a missing tile is an error naming its URL.
    tilesource = CopernicusTiles(sys, tiles; cachedir = dir, download = false)
    paths = Dict(t => tilepath!(tilesource, t) for t in tiles)
    slots = opt.slots === nothing ? length(tiles) : opt.slots
    cache = StripedLRUCache{Vector{Float32}}(k -> loadtile(tilesource, tiles[k]);
        slots = slots, stripes = 1)
    ids = TileIds(sys, tiles)
    dem = CountingTiles(TiledDEM(ids, cache))
    grid = DGG.PartialGrid(sys, 1, ids)
    space = DGG.DGGSpace(grid; chunklevel = 0)
    bytes = sum(filesize(paths[t]) for t in tiles)
    return (; sys, tiles, paths, cache, dem, grid, space, slots, bytes)
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
    verifynearest(out, dst, src, dem, n; seed) -> NamedTuple

Check `n` random destination cells against the source read directly. The sample
site is the destination cell's centroid, the source cell is `cellat` of that
site, and the regridded value must equal the source pixel BIT for bit. Where
`cellat` finds nothing the cell must be missing.
"""
function verifynearest(out, dst, src, dem, n::Int; seed = 20260825)
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

"""
    stencilpass(smp, sites, n) -> (mapped, entries)

One `weightsat!` per destination sample site, into a row reused across all of
them and nothing else: the stencil's own cost, with no source read and no weight
stored. `mapped` counts the sites the source could place and `entries` the
nonzeros they emitted.
"""
function stencilpass(smp, sites, n::Int)
    row = GR.WeightRow()
    mapped = 0
    entries = 0
    for i in 1:n
        GR.ismapped(GR.weightsat!(row, smp, sites[i])) || continue
        mapped += 1
        entries += length(row)
    end
    return (mapped, entries)
end

"""
    verifypoint(out, dst, src, dem, n) -> (checked, mapped, missing, wrong)

`n` random destination cells of the barycentric arm, against the stencil
recomputed here and the source pixels read directly.
"""
function verifypoint(out, dst, src, dem, n::Int; seed = 20260826)
    n <= 0 && return (0, 0, 0, 0)
    rng = Random.MersenneTwister(seed)
    smp = GR.sampler(GR.BarycentricPoint(), src.space)
    row = GR.WeightRow()
    sites = GR.samplesites(dst.space)
    checked = mapped = absent = wrong = 0
    for i in Random.rand(rng, 1:DGG.ncells(dst.grid), n)
        checked += 1
        if !GR.ismapped(GR.weightsat!(row, smp, sites[i]))
            absent += 1
            isfinite(out[i]) && (wrong += 1)
            continue
        end
        mapped += 1
        want = 0.0
        good = true
        for k in 1:length(row)
            v = dem[row.indices[k]]
            isfinite(v) || (good = false; break)
            want += row.weights[k] * v
        end
        good || continue
        isapprox(out[i], Float32(want); rtol = 1e-6) || (wrong += 1)
    end
    return (checked, mapped, absent, wrong)
end

# ---------------------------------------------------------------------------
# cost
# ---------------------------------------------------------------------------

plan_for(name, src, dstspace, opt) = let (method, policy) = methodpair(name)
    GR.plan_regrid(src.dem; to = dstspace, from = src.space, method,
        missingpolicy = policy, lazy = true, budget = opt.budget)
end

"A cold and a warm profiled read of the first `opt.profile` destination tiles."
function profilearm(name, opt, src, dst)
    k = min(opt.profile, Int(GR.nchunks(dst.space)))
    rs = [GR.ownedindices(dst.space, c) for c in 1:k]
    # Every source tile is resident first, so neither listing profiles GDAL.
    foreach(src.cache, 1:length(src.tiles))
    AP = GR.regrid(src.dem, plan_for(name, src, dst.space, opt))
    Profile.init(n = 40_000_000, delay = 0.002)
    for label in ("cold read of $k destination tiles ($(commas(sum(length, rs))) cells): weights built here",
                  "warm read of the same $k tiles: weights cached, source loads and application only")
        Profile.clear()
        Profile.@profile for r in rs
            AP[r]
        end
        flatself(stdout, label)
    end
    Profile.clear()
    return nothing
end

"One method's plan, first tile and two sweeps, on a source whose tile cache is empty."
function arm(name, opt, dir, dst)
    method, _ = methodpair(name)
    ndst = Int(DGG.ncells(dst.grid))
    src = sourceside(opt, dir)
    println("\n$name — $(repr(method))")
    p = phase("plan", () -> plan_for(name, src, dst.space, opt))
    plan = p.value
    graph = GR.dependencies(plan)
    cand = [length(GR.sourcesof(graph, c)) for c in 1:Int(GR.nchunks(dst.space))]
    @printf("  %-26s min %d, median %.1f, max %d\n", "candidates per tile",
        minimum(cand), Statistics.median(cand), maximum(cand))

    src1 = sourceside(opt, dir)
    A1 = GR.regrid(src1.dem, plan_for(name, src1, dst.space, opt))
    tile1 = GR.ownedindices(dst.space, 1)
    phase("first tile", () -> A1[tile1])
    @printf("  %-26s %s cells; %d chunk reads, %s pixels, %d GeoTIFF decodes\n",
        "tile 1", commas(length(tile1)), sum(values(src1.dem.reads); init = 0),
        commas(src1.dem.pixels[]), src1.cache.loads[])
    src1 = A1 = nothing
    opt.sweeps || return (; name, cold = p, warm = p)

    if method isa GR.BarycentricPoint
        sites = GR.samplesites(dst.space)
        # The complete level makes a post's index arithmetic; against the holding
        # it is a search of the holding's ids, four times a query.
        whole = DGG.DGGSpace(DGG.levelgrid(src.sys, 1); chunklevel = 0)
        for (label, space) in (("stencils", src.space), ("stencils, complete level", whole))
            smp = GR.sampler(method, space)
            stencilpass(smp, sites, min(ndst, 1024))          # warm the row's buffers
            st = phase(label, () -> stencilpass(smp, sites, ndst))
            @printf("  %-26s %.1f ns and %.1f B each, %s entries over %s placed\n",
                "one stencil", 1e9 * st.time / ndst, st.bytes / ndst,
                commas(st.value[2]), commas(st.value[1]))
        end
    end

    method isa PairNearest && (method.placed[] = 0)
    out = Vector{Float32}(undef, ndst)
    A = GR.regrid(src.dem, plan)
    cold = phase("cold sweep", () -> sweep!(out, A, dst.space))
    coldreads = copy(src.dem.reads)
    coldpixels = src.dem.pixels[]
    colddecodes = src.cache.loads[]
    coldstats = deepcopy(GR.residency(A))
    w = weightstate(plan)
    coldplaced = method isa PairNearest ? method.placed[] : 0
    resetreads!(src.dem)
    warm = phase("warm sweep", () -> sweep!(out, A, dst.space))

    nfinite = count(isfinite, out)
    lo, hi = nfinite > 0 ? extrema(Iterators.filter(isfinite, out)) : (NaN, NaN)
    @printf("  %-26s %s of %s (%.2f%%), range [%.2f, %.2f] m\n", "finite values",
        commas(nfinite), commas(ndst), 100nfinite / ndst, lo, hi)
    @printf("  %-26s %d of %d chunks, %s reads, %s pixels, %d GeoTIFF decodes\n",
        "cold source", length(coldreads), length(src.tiles),
        commas(sum(values(coldreads); init = 0)), commas(coldpixels), colddecodes)
    @printf("  %-26s %d chunks, %s reads\n", "warm source", length(src.dem.reads),
        commas(sum(values(src.dem.reads); init = 0)))
    println("  residency cold ", coldstats)
    println("  residency warm ", GR.residency(A))
    @printf("  %-26s %s B of a %s B budget in %d cached tiles / %d blocks; manifest %d..%d, %s total\n",
        "weights", commas(w.accounted), commas(w.budget), w.ntiles, w.nblocks,
        w.manifest_min, w.manifest_max, commas(w.manifest_sum))
    method isa PairNearest && @printf("  %-26s %s cells placed cold, %s warm\n",
        "pair-route locations", commas(coldplaced), commas(method.placed[] - coldplaced))

    # Verified while this arm's source and values are alive, so no arm holds
    # another's plan or weights and the peak below is one arm's.
    if opt.check > 0 && method isa GR.BarycentricPoint
        c = verifypoint(out, dst, src, src.dem, opt.check)
        @printf("  %-26s %s cells, %s mapped, %s outside, %s wrong\n", "check",
            commas(c[1]), commas(c[2]), commas(c[3]), commas(c[4]))
    elseif opt.check > 0 && !(method isa DGG.Conservative)
        v = verifynearest(out, dst.space, src.space, src.dem, opt.check)
        @printf("  %-26s %s cells, %s mapped, %s outside, %d wrong\n", "check",
            commas(v.checked), commas(v.mapped), commas(v.outside), v.bad)
        v.firstbad === nothing || println("  first mismatch ", v.firstbad)
    end
    opt.profile > 0 && profilearm(name, opt, src, dst)
    @printf("  %-26s %.2f GiB\n", "peak RSS after this arm", Sys.maxrss() / 2^30)
    flush(stdout)
    return (; name, cold, warm)
end

function cost(opt)
    dir = tiledir(opt)
    sys7 = DGG.IGeo7System()
    println(provenance())
    @printf("\ncost  res=GLO-%d  box=(%d,%d,%d,%d)  methods=%s  budget=%s B  tiles from %s\n",
        opt.res, opt.box..., join(opt.methods, ","), commas(opt.budget), dir)

    probe = sourceside(opt, dir)
    level, chunklevel = destinationlevels(opt, sys7, probe.space)
    @printf("source   %d tiles, %s pixels, %s B of GeoTIFF, %d cache slots, %.2f m pixels\n",
        length(probe.tiles), commas(DGG.ncells(probe.grid)), commas(probe.bytes),
        probe.slots, DGG.cellsize(probe.space))
    @printf("dest     IGeo7 level %d%s, chunk level %d, %.2f m cells\n", level,
        opt.level === nothing ? " (bare-system rule)" : " (override)", chunklevel,
        DGG.cellsize(sys7, level))
    flush(stdout)

    # Compilation, on sources of their own, so every arm starts with an empty
    # decoded-tile cache. The GeoTIFFs stay in the OS page cache either way.
    warm = @timed for name in opt.methods
        w = sourceside(opt, dir)
        ex = Extents.Extent(X = (Float64(opt.box[1]) + 0.40, Float64(opt.box[1]) + 0.42),
            Y = (Float64(opt.box[3]) + 0.40, Float64(opt.box[3]) + 0.42))
        set = DGG.query(sys7, DGG.MultiOrderCoverage(ex); level)
        sp = DGG.DGGSpace(DGG.PartialGrid(DGG.CellVector(set; level)); chunklevel)
        collect(GR.regrid(w.dem, plan_for(name, w, sp, opt)))
    end
    @printf("\nwarm-up  %.3f s for %d method(s) (compilation; in no arm below)\n",
        warm.time, length(opt.methods))

    println("\narms")
    dst = phase("destination", () -> destinationside(sys7, opt.box, level, chunklevel)).value
    ndst = Int(DGG.ncells(dst.grid))
    phase("source", () -> DGG.DGGSpace(probe.grid; chunklevel = 0))
    @printf("  %-26s %s cells in %s tiles\n", "destination",
        commas(ndst), commas(Int(GR.nchunks(dst.space))))
    flush(stdout)

    probe = nothing                       # the sizing source is not an arm's
    results = [arm(name, opt, dir, dst) for name in opt.methods]

    if opt.sweeps
        println("\n$(Threads.nthreads()) threads, $(commas(ndst)) destination cells")
        @printf("%-14s %10s %12s %12s %10s %12s %12s %8s\n",
            "method", "cold s", "cold ns/cell", "cold cells/s",
            "warm s", "warm ns/cell", "warm cells/s", "cold x")
        for r in results
            @printf("%-14s %10.3f %12.1f %12s %10.3f %12.1f %12s %8.2f\n", r.name,
                r.cold.time, 1e9 * r.cold.time / ndst,
                commas(round(Int, ndst / r.cold.time)),
                r.warm.time, 1e9 * r.warm.time / ndst,
                commas(round(Int, ndst / r.warm.time)),
                r.cold.time / results[1].cold.time)
        end
    end
    @printf("\npeak RSS   %.2f GiB\n", Sys.maxrss() / 2^30)
    return results
end

# ---------------------------------------------------------------------------
# semantics
# ---------------------------------------------------------------------------

# The scaled conformance lattice: 30 latitude intervals to the degree instead of
# 1200 or 3600, and every band ratio of a real Copernicus product. A tile of it
# is 900 posts rather than 1.4 million, so a whole arm fits in a second.
const TWIN = CD.CopernicusDEMSystem{30}()

"""
    haspolarpolicy() -> Bool

Whether the point method in this build carries a polar policy, which it states
by taking one as a field.
"""
haspolarpolicy() = :poles in fieldnames(DGG.BarycentricPoint)

# ---------------------------------------------------------------------------
# Comparing two fields over the same cells
# ---------------------------------------------------------------------------

"""
    difference(base, other) -> NamedTuple

What separates two regridded fields over the same destination cells: what
fraction each placed, what fraction both placed, the spread between them where
both did, and how many cells one placed and the other did not — which for a
point method against an area one is the rim.
"""
function difference(base::AbstractVector, other::AbstractVector)
    n = length(base)
    n == length(other) || throw(DimensionMismatch(
        "the two fields cover $(n) and $(length(other)) cells"))
    both = 0
    total = 0.0
    worst = 0.0
    bmiss = 0
    omiss = 0
    bonly = 0
    oonly = 0
    @inbounds for i in 1:n
        a, p = Float64(base[i]), Float64(other[i])
        ga, gp = isfinite(a), isfinite(p)
        ga || (bmiss += 1)
        gp || (omiss += 1)
        ga && !gp && (bonly += 1)
        gp && !ga && (oonly += 1)
        (ga && gp) || continue
        both += 1
        d = abs(a - p)
        total += d
        worst = max(worst, d)
    end
    return (; n, both, bothfrac = both / n,
        baseplaced = n - bmiss, otherplaced = n - omiss,
        baseonly = bonly, otheronly = oonly,
        basemissing = bmiss / n, othermissing = omiss / n,
        meanabs = both == 0 ? NaN : total / both,
        maxabs = both == 0 ? NaN : worst)
end

"""
    fielderror(values, sites, f) -> (mean, max, n)

How far a regridded field is from the analytic field it was built from, at the
destination sample sites, over the cells it placed.
"""
function fielderror(values::AbstractVector, sites, f)
    n = 0
    total = 0.0
    worst = 0.0
    for i in eachindex(values)
        isfinite(values[i]) || continue
        n += 1
        d = abs(Float64(values[i]) - f(sites[i]))
        total += d
        worst = max(worst, d)
    end
    return (n == 0 ? NaN : total / n, n == 0 ? NaN : worst, n)
end

"Every method's field against the first one's, one block per comparison."
function reportdifferences(results, unit)
    base = results[1]
    println("\nagainst $(base.name), in $unit")
    out = Any[]
    for r in results[2:end]
        d = difference(base.out, r.out)
        @printf("  %-26s %s of %s cells (%.2f%%)\n", "$(r.name): both placed",
            commas(d.both), commas(d.n), 100 * d.bothfrac)
        @printf("  %-26s %s %.4f%%, %s %.4f%%\n", "  missing", base.name,
            100 * d.basemissing, r.name, 100 * d.othermissing)
        @printf("  %-26s %s only %s, %s only %s (%.4f%% of the %s %s placed)\n",
            "  placed by one alone", base.name, commas(d.baseonly), r.name,
            commas(d.otheronly),
            d.baseplaced == 0 ? NaN : 100 * d.baseonly / d.baseplaced,
            commas(d.baseplaced), base.name)
        @printf("  %-26s mean %.6g, max %.6g\n", "  |difference|", d.meanabs, d.maxabs)
        push!(out, (; name = r.name, d))
    end
    flush(stdout)
    return out
end

# ---------------------------------------------------------------------------
# One method over one pair
# ---------------------------------------------------------------------------

"""
    sweepmethod(name, dem, srcspace, dstspace, ndst, budget) -> NamedTuple

One method's cold sweep over every destination tile, and what it read.

The source is handed in already built, so its cache state is the caller's to
control: each method here gets one of its own, decoded from nothing.
"""
function sweepmethod(name, dem, srcspace, dstspace, ndst::Int, budget::Int)
    method, policy = methodpair(name)
    println("\n$name — $(repr(method)) under $(repr(policy))")
    plan = GR.plan_regrid(dem; to = dstspace, from = srcspace, method = method,
        missingpolicy = policy, lazy = true, budget = budget)
    A = GR.regrid(dem, plan)
    # The source's own precision: rounding an affine field to Float32 costs more
    # than the interpolation that reproduces it does.
    out = Vector{eltype(A)}(undef, ndst)
    cold = phase("cold sweep", () -> sweep!(out, A, dstspace))
    st = GR.residency(A)
    reads = dem isa CountingTiles ? sum(values(dem.reads); init = 0) : nreads(dem)
    @printf("  %-26s %s chunk reads, %s loads, %s cache hits, %s peak source B\n",
        "source", commas(reads), commas(st.loads), commas(st.hits),
        commas(st.peakbytes))
    @printf("  %-26s %s of %s placed, %.2f GiB peak RSS\n", "values",
        commas(count(isfinite, out)), commas(ndst), Sys.maxrss() / 2^30)
    flush(stdout)
    return (; name, out, time = cold.time, reads, loads = st.loads, hits = st.hits,
        peakbytes = st.peakbytes)
end

# ---------------------------------------------------------------------------
# The scaled lattice twin, and analytic fields over it
# ---------------------------------------------------------------------------

"A source over a holding of twin posts that counts the reads made of it."
struct CountingPosts{C} <: DiskArrays.AbstractDiskArray{Float64,1}
    values::Vector{Float64}
    chunks::C
    reads::Threads.Atomic{Int}
end

Base.size(A::CountingPosts) = size(A.values)
DiskArrays.eachchunk(A::CountingPosts) = A.chunks
DiskArrays.haschunks(::CountingPosts) = DiskArrays.Chunked()

function DiskArrays.readblock!(A::CountingPosts, out, r::AbstractUnitRange)
    Threads.atomic_add!(A.reads, 1)
    out .= view(A.values, r)
    return out
end

nreads(A::CountingPosts) = A.reads[]

"""
    twinside(box, f) -> NamedTuple

The twin lattice's posts over `box`, held as a partial grid chunked one 1-degree
tile per chunk, carrying `f` evaluated at each post's own sample site.

The holding is a strict subset of the twin's complete level, so a stencil
reaching a post outside it is a rim — which is what the `rim` arm looks at, and
what every arm's interior is measured away from.
"""
function twinside(box, f)
    w, e, s, n = box
    tiles = sort!([Int(CD.tilecell(TWIN, lat, lon).index)
                   for lat in s:(n - 1) for lon in w:(e - 1)])
    ids = TileIds(TWIN, tiles)
    grid = DGG.PartialGrid(TWIN, 1, ids)
    space = DGG.DGGSpace(grid; chunklevel = 0)
    sites = GR.samplesites(space)
    values = [f(sites[i]) for i in 1:Int(GR.ncells(space))]
    widths = [length(GR.ownedindices(space, k)) for k in 1:Int(GR.nchunks(space))]
    chunks = DiskArrays.GridChunks(DiskArrays.IrregularChunks(; chunksizes = widths))
    dem = CountingPosts(values, chunks, Threads.Atomic{Int}(0))
    return (; tiles, grid, space, dem, sites)
end

"""
    capdestination(sys7, south, level, chunklevel) -> NamedTuple

The IGeo7 cells of the cap north of latitude `south`, as a `PartialGrid` tiled
at `chunklevel`.

`destinationside` takes whole degrees, and the region this arm needs is a strip
one part in a few hundred of a degree deep — the gap between a lattice's
polemost post row and the pole itself.
"""
function capdestination(sys7, south::Float64, level::Int, chunklevel::Int)
    ex = Extents.Extent(X = (-180.0, 180.0), Y = (south, 90.0))
    set = DGG.query(sys7, DGG.MultiOrderCoverage(ex); level = level)
    grid = DGG.PartialGrid(DGG.CellVector(set; level = level))
    return (; grid, space = DGG.DGGSpace(grid; chunklevel = chunklevel))
end

"Longitude and latitude in degrees of a unit-sphere point."
lonlat(p) = let g = US.GeographicFromUnitSphere()(p)
    (Float64(g[1]), Float64(g[2]))
end

"An affine field in longitude and latitude, in metres."
affinefield(p) = let (lon, lat) = lonlat(p)
    1000.0 + 37.0 * lon - 61.0 * lat
end

"A field with curvature, in metres, over the same box."
curvedfield(p) = let (lon, lat) = lonlat(p)
    1000.0 + 90.0 * (lon - 11.0)^2 + 140.0 * (lat - 50.0)^2
end

"""
    polarfield(p) -> Float64

A field for the polar cap, in metres: a ramp in the plane tangent at the pole.

Longitude is not a coordinate at a pole — every meridian meets there — so a
field affine in longitude is discontinuous exactly where this arm looks, and
would report its own jump rather than the methods'. This one is affine in the
sphere's own `x` and `y`, which is affine across the cap and single-valued on
the whole sphere.
"""
polarfield(p) = 1000.0 + 7.0e5 * p[1] + 4.0e5 * p[2]

# ---------------------------------------------------------------------------
# The arms
# ---------------------------------------------------------------------------

"""
    realarm(opt, dir, sys7) -> NamedTuple

The listed GLO-90 tiles of `opt.box` into IGeo7, under every named method.
"""
function realarm(opt, dir, sys7)
    println("\n", "="^78, "\nreal — GLO-$(opt.res) $(opt.box) elevations")
    probe = sourceside(opt, dir)
    level, chunklevel = destinationlevels(opt, sys7, probe.space)
    @printf("source   %d tiles, %s posts, %.2f m spacing\n", length(probe.tiles),
        commas(DGG.ncells(probe.grid)), DGG.cellsize(probe.space))
    probe = nothing
    dst = destinationside(sys7, opt.box, level, chunklevel)
    ndst = Int(DGG.ncells(dst.grid))
    @printf("dest     IGeo7 level %d, chunk level %d, %.2f m cells, %s cells in %s tiles\n",
        level, chunklevel, DGG.cellsize(sys7, level), commas(ndst),
        commas(Int(GR.nchunks(dst.space))))
    flush(stdout)

    results = map(opt.methods) do name
        src = sourceside(opt, dir)
        sweepmethod(name, src.dem, src.space, dst.space, ndst, opt.budget)
    end
    return (; arm = :real, level, ndst, results,
        diffs = reportdifferences(results, "metres"))
end

"""
    twinarm(label, opt, sys7, f; box = opt.twin, dest = opt.twin) -> NamedTuple

An analytic field on the twin holding `box`, regridded onto the IGeo7 cells of
`dest` under every named method, and each method's own error against the field.
"""
function twinarm(label, opt, sys7, f; box = opt.twin, dest = opt.twin)
    println("\n", "="^78, "\n$label — the scaled twin over $(box), destination $(dest)")
    src = twinside(box, f)
    level = (opt.level === nothing ? DGG.levelfor(sys7, src.space) : opt.level) + opt.fine
    chunklevel = opt.chunk === nothing ? max(first(DGG.levels(sys7)), level - 3) : opt.chunk
    dst = destinationside(sys7, dest, level, chunklevel)
    ndst = Int(DGG.ncells(dst.grid))
    @printf("source   %d tiles, %s posts, %.2f m spacing in %s chunks\n",
        length(src.tiles), commas(DGG.ncells(src.grid)), DGG.cellsize(src.space),
        commas(Int(GR.nchunks(src.space))))
    @printf("dest     IGeo7 level %d, %.2f m cells, %s cells in %s tiles\n",
        level, DGG.cellsize(sys7, level), commas(ndst),
        commas(Int(GR.nchunks(dst.space))))
    flush(stdout)

    results = map(opt.methods) do name
        one = twinside(box, f)
        sweepmethod(name, one.dem, one.space, dst.space, ndst, opt.budget)
    end

    # The strict threshold against the loose one on a source with no missing
    # values: a complete point row's valid weight is its whole row, so
    # `Weighted(1)` blanks nothing `Weighted(0)` keeps.
    strict = findfirst(r -> r.name == "point", results)
    if strict !== nothing
        loose = twinside(box, f)
        looseplan = GR.plan_regrid(loose.dem; to = dst.space, from = loose.space,
            method = DGG.BarycentricPoint(), missingpolicy = DGG.Weighted(0),
            lazy = true, budget = opt.budget)
        looseA = GR.regrid(loose.dem, looseplan)
        loosevals = Vector{eltype(looseA)}(undef, ndst)
        sweep!(loosevals, looseA, dst.space)
        @printf("\n  %-26s %s under Weighted(1), %s under Weighted(0)\n",
            "point cells placed", commas(count(isfinite, results[strict].out)),
            commas(count(isfinite, loosevals)))
    end

    sites = GR.samplesites(dst.space)
    println("\nagainst the analytic field, at the destination sample sites")
    errs = map(results) do r
        m, x, n = fielderror(r.out, sites, f)
        @printf("  %-26s mean %.6g m, max %.6g m over %s placed cells\n",
            r.name, m, x, commas(n))
        (; name = r.name, mean = m, max = x, placed = n)
    end
    return (; arm = Symbol(label), level, ndst, results, errs,
        diffs = reportdifferences(results, "metres"))
end

"""
    polararm(opt, sys7) -> NamedTuple

The strip between the twin lattice's polemost post row and the pole.

No dual cell of source sample sites covers it: the natural one has that whole
row as its corners. The point method's `poles` policy decides what happens
there, and this arm runs it both ways against the area method.
"""
function polararm(opt, sys7)
    println("\n", "="^78, "\npolar — between the polemost post row and the pole")
    if !haspolarpolicy()
        println("  the point method in this build carries no polar policy, so a\n" *
                "  query poleward of the polemost post row is unmapped by\n" *
                "  construction and there is nothing to compare. Rerun once the\n" *
                "  policy is in the build.")
        return (; arm = :polar, ran = false)
    end
    # The polemost row of tiles, all the way round: the policy takes the nearest
    # site of that row, so a holding short in longitude would make its own rim
    # the thing measured instead of the policy.
    box = (-180, 180, 89, 90)
    src = twinside(box, polarfield)
    st = CD._pointstate(TWIN)
    edge = CD._sitelat(st, Int64(0))
    # The strip is a quarter of a post's latitude spacing deep, so the rule's
    # level resolves it with a cell or two. Three levels finer is 343 times the
    # cells, which is what makes a fraction a number.
    level = (opt.level === nothing ? DGG.levelfor(sys7, src.space) : opt.level) +
            opt.fine + 2
    chunklevel = opt.chunk === nothing ? max(first(DGG.levels(sys7)), level - 3) : opt.chunk
    dst = capdestination(sys7, edge, level, chunklevel)
    ndst = Int(DGG.ncells(dst.grid))
    @printf("source   %d tiles, %s posts in %s chunks; polemost row at %.6f N\n",
        length(src.tiles), commas(DGG.ncells(src.grid)),
        commas(Int(GR.nchunks(src.space))), edge)
    @printf("dest     IGeo7 level %d, %.2f m cells, %s cells north of it in %s tiles\n",
        level, DGG.cellsize(sys7, level), commas(ndst),
        commas(Int(GR.nchunks(dst.space))))
    flush(stdout)

    results = map(opt.methods) do name
        one = twinside(box, polarfield)
        sweepmethod(name, one.dem, one.space, dst.space, ndst, opt.budget)
    end

    # The same point method with the policy switched off: the difference between
    # the two point columns is the policy's whole effect.
    off = twinside(box, polarfield)
    plan = GR.plan_regrid(off.dem; to = dst.space, from = off.space,
        method = DGG.BarycentricPoint(; poles = nothing),
        missingpolicy = DGG.Weighted(1), lazy = true, budget = opt.budget)
    bareA = GR.regrid(off.dem, plan)
    bare = Vector{eltype(bareA)}(undef, ndst)
    sweep!(bare, bareA, dst.space)
    point = findfirst(r -> r.name == "point", results)
    if point !== nothing
        @printf("\n  %-26s %s of %s with poles = NearestCell(), %s with poles = nothing\n",
            "cells placed", commas(count(isfinite, results[point].out)),
            commas(ndst), commas(count(isfinite, bare)))
    end
    sites = GR.samplesites(dst.space)
    println("\nagainst the analytic field, at the destination sample sites")
    for r in results
        m, x, n = fielderror(r.out, sites, polarfield)
        @printf("  %-26s mean %.6g m, max %.6g m over %s placed cells\n",
            r.name, m, x, commas(n))
    end
    let (m, x, n) = fielderror(bare, sites, polarfield)
        @printf("  %-26s mean %.6g m, max %.6g m over %s placed cells\n",
            "point, poles = nothing", m, x, commas(n))
    end
    return (; arm = :polar, ran = true, level, ndst, results, bare,
        diffs = reportdifferences(results, "metres"))
end

# ---------------------------------------------------------------------------
# The run
# ---------------------------------------------------------------------------

function semantics(opt)
    dir = tiledir(opt)
    sys7 = DGG.IGeo7System()
    println(provenance())
    @printf("\narms=%s  methods=%s  res=GLO-%d  box=%s  twin=%s  fine=+%d  budget=%s B\n",
        join(opt.arms, ","), join(opt.methods, ","), opt.res, string(opt.box),
        string(opt.twin), opt.fine, commas(opt.budget))
    :real in opt.arms && println("tiles from $dir")

    # Compilation, on a pair of its own, so no arm below starts warm.
    warm = @timed for name in opt.methods
        w = twinside((0, 1, 0, 1), affinefield)
        lvl = DGG.levelfor(sys7, w.space)
        sp = destinationside(sys7, (0, 1, 0, 1), lvl,
            max(first(DGG.levels(sys7)), lvl - 3))
        method, policy = methodpair(name)
        collect(GR.regrid(w.dem; to = sp.space, from = w.space, method = method,
            missingpolicy = policy, lazy = true, budget = opt.budget))
    end
    @printf("warm-up  %.3f s for %d method(s) (compilation; in no arm below)\n",
        warm.time, length(opt.methods))
    flush(stdout)

    out = Any[]
    for a in opt.arms
        a === :real && push!(out, realarm(opt, dir, sys7))
        a === :linear && push!(out, twinarm("linear", opt, sys7, affinefield))
        a === :quadratic && push!(out, twinarm("quadratic", opt, sys7, curvedfield))
        # The destination reaches a whole degree past the holding on every side,
        # so its edge cells need posts the holding does not have.
        a === :rim && push!(out, twinarm("rim", opt, sys7, affinefield;
            dest = (opt.twin[1] - 1, opt.twin[2] + 1, opt.twin[3] - 1, opt.twin[4] + 1)))
        a === :polar && push!(out, polararm(opt, sys7))
    end

    println("\n", "="^78, "\nsummary — each method against $(opt.methods[1])")
    @printf("%-11s %-14s %9s %10s %12s %12s %12s\n", "arm", "method", "cold s",
        "missing", "base only", "mean |diff|", "max |diff|")
    for r in out
        get(r, :ran, true) || continue
        @printf("%-11s %-14s %9.3f %9.3f%% %12s %12s %12s\n", r.arm,
            r.results[1].name, r.results[1].time,
            100 * (1 - r.diffs[1].d.baseplaced / r.diffs[1].d.n), "-", "-", "-")
        for (k, e) in enumerate(r.diffs)
            @printf("%-11s %-14s %9.3f %9.3f%% %12s %12.6g %12.6g\n", "", e.name,
                r.results[k + 1].time, 100 * e.d.othermissing,
                commas(e.d.baseonly), e.d.meanabs, e.d.maxabs)
        end
    end
    @printf("\n%d threads, peak RSS %.2f GiB\n", Threads.nthreads(),
        Sys.maxrss() / 2^30)
    return out
end

function main(args = ARGS)
    command = !isempty(args) && !occursin('=', args[1]) ? Symbol(args[1]) : :cost
    rest = command === :cost && (isempty(args) || occursin('=', args[1])) ? args : args[2:end]
    command in (:cost, :semantics) || error("usage: copdem.jl [cost|semantics] key=value...")
    opt = options(command, rest)
    return command === :cost ? cost(opt) : semantics(opt)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
