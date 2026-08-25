# What one large destination tile costs a point method, on each of the two
# build routes it can take.
#
#     julia -t auto --project=benchmark benchmark/point_tile_baseline.jl
#
# A `ChunkedPlan` builds weights one `(destination tile, source chunk)` pair at
# a time, and each pair's `buildweights!` walks every destination cell of the
# tile. A point method therefore locates each destination point once per
# CANDIDATE source chunk, not once, while the stencils it finally keeps are four
# entries per destination however many chunks were searched. This file prices
# that gap on a raster large enough for it to matter, and records what a tile's
# weights and reads cost beside it.
#
# Four arms over the same tile, all with `BilinearPoint`:
#
#   1. `one pass`        one `buildweights!` over the whole source index range.
#                        Every destination is located once and every entry is
#                        kept, so this is the point-location floor: the work a
#                        build that partitioned by source chunk after the fact
#                        would do.
#   2. `pair loop`       one `buildblock` per candidate chunk, which is what the
#                        plan does. Same stencils, `k` passes over the tile.
#   3. `cold tile`       a fresh plan and one lazy read of the whole tile:
#                        arm 2 plus source loads and the weight application.
#   4. `warm tile`       the same read against the plan arm 3 left warm, so the
#                        weights are cached and only application and source
#                        loads remain.
#
# Four more arms over that same tile with `BarycentricPoint`, which supplies a
# `sampler` and is therefore built one destination tile at a time:
#
#   5. `tile weights`    one `tileweights` build of the whole tile, sampler
#                        preparation included. Every destination is located once
#                        and every entry is filed under the source chunk owning
#                        it, so this replaces arm 2, not arm 1.
#   6. `plan build`      the same spaces, chunk index and relation as arm 3's,
#                        prepared for the other method.
#   7. `cold read`       a fresh plan and one lazy read of the whole tile: arm 5
#                        plus the loads its manifest names and the application.
#   8. `warm read`       the same read against the plan arm 7 left warm.
#
# The tile-route accounting rows beside them come from one further instrumented
# read of a clean plan, and describe what that route does rather than what the
# pair route did: the point locations one build performs, counted through a
# sampler that records every `weightsat!` query and answers with
# `BarycentricPoint`'s own stencil, rather than asserted; the nonzeros and bytes
# the single cached `TileWeights` holds, manifest included; the length of that
# manifest against the number of candidates the relation names; and the source
# chunks the read touched with how many times each was fetched. The counted
# build's manifest is checked against the plan's, so the counted numbers belong
# to the tile the plan built.
#
# SIZES. The source is a global 0.125-degree raster, 2880x1440 = 4,147,200 cells
# in 8x4 = 32 chunks of 360x360, which is the shape of an ordinary global
# climate or DEM mosaic and large enough that a chunk does not fit in cache. The
# destination is a single 750x500 = 375,000-cell tile over 60 degrees of
# longitude and 40 of latitude: one tile, so the tiling contributes nothing to
# the count, and a region whose caps reach several source chunks rather than one
# or all of them. The number of candidates is read from the plan's relation, not
# assumed.
#
# Every arm is warmed once and reported as the minimum of three samples;
# allocations come from `@allocated` on one further call. `provenance()` stamps
# each run with the Julia version, thread and CPU counts and the machine's power
# mode, because the same code clocks differently under low power. Compare a
# later run with this one only through a same-session ratio, never by carrying
# an absolute across sessions or power modes. Weight residency is the plan's own
# block accounting, `storagebytes` over `PerChunk` (`_blockbytes`: 16 bytes a
# nonzero, 8 a column pointer, 8 a reference entry), with `Base.summarysize` of
# the same blocks beside it. Source residency is
# `residency(::LazyRegridArray)`, and the per-chunk read counts come from a
# counting `DiskArrays` source that records every `readblock!` range.
#
# Storage is unbounded `PerChunk()` so nothing is evicted: the reported weight
# bytes are what the tile's blocks actually occupy, not what a budget allows.
# `weightbudget(plan.budget)` is printed beside them for comparison.
#
# RECORDED 2026-08-25, M-series macOS, Julia 1.12.6, 8 threads of 12 CPUs,
# `powermode 2`, warm-up plus minimum of three samples:
#
#   arm                                  time        allocated
#   one pass over the tile          67.801 ms      101,422,496 B
#   pair loop (4 candidates)       324.564 ms      221,472,320 B
#   plan construction                0.066 ms            4,320 B
#   cold tile read                 331.528 ms      254,579,792 B
#   warm tile read                   4.478 ms       19,275,648 B
#   tile weights, one pass         105.406 ms      212,468,720 B
#   plan construction, tile route    0.062 ms           78,304 B
#   cold read, tile route          105.509 ms      245,848,432 B
#   warm read, tile route            3.021 ms       19,349,536 B
#
#   375,000 destination cells, 4 candidate source chunks [12, 13, 20, 21], all
#   four contributing, 1,500,000 stencil entries on both routes, 4.00 per
#   destination cell
#
#   chunk-pair route  1,500,000 point locations where 375,000 would do
#                     4 cached blocks, weight bytes 40,147,488 (plan
#                     accounting) / 28,147,912 (summarysize)
#                     4 source chunks read, each once; peak source residency
#                     4,147,200 B
#   tile route          375,000 point locations, one per destination
#                     1 cached tile of 4 blocks over a 4-chunk manifest, weight
#                     bytes 40,147,584 / 28,148,112
#                     4 source chunks read, each once; peak source residency
#                     4,147,200 B
#
# The pair loop is 4.79x one pass at k = 4: the destination pass is repeated per
# candidate and nothing else about the build changes. Weight construction is
# 98.6% of a cold tile read, so this repetition, not source I/O or the weight
# application, is what a large point tile costs on that route.
#
# The fused build locates each of the 375,000 destinations once, and the manifest
# it emits equals the relation's candidate list here, so the reads are the same
# four: 3.08x the pair loop for the same 1,500,000 entries, and a cold read 3.14x
# the pair route's. Its 1.55x over the one-pass floor is not overhead over the
# same work — the floor arm is a different kernel emitting one unpartitioned COO,
# where the fused build also files every entry under its owning chunk and
# assembles four blocks. Weight construction is 97.1% of its own cold read, and
# the 96 extra bytes it holds are the manifest.
#
# These absolutes belong to that machine state. A later comparison must be a
# ratio inside one session: re-run this file on the unchanged tree, then again
# on the changed one, and compare those two.

using Printf
import DimensionalData as DD
import DiskArrays
using DimensionalData: X, Y
using GlobalRegridding
import GlobalRegridding as GR

const NSAMPLES = parse(Int, get(ENV, "DGG_POINT_TILE_SAMPLES", "3"))

# Source shape: a global raster in tens of chunks.
const NX, NY = 2880, 1440
const CX, CY = 360, 360

# Destination shape: one tile over a region reaching several source chunks. The
# 0.08-degree step shares no cell centre with the source's 0.125-degree one, so
# every stencil has four entries and the mean is not diluted by destinations
# sitting exactly on a source sample site.
const DNX, DNY = 750, 500
const DLON, DLAT = (-30.0, 30.0), (-20.0, 20.0)

# --- a source that records what is read from it ----------------------------

"""
    CountingSource(parent, chunks)

A chunked `DiskArrays` array over `parent` that counts `readblock!` calls per
requested range. `reads` therefore gives both how many source chunks a read
touched and how many times each was fetched.
"""
struct CountingSource{T,N,A<:AbstractArray{T,N}} <: DiskArrays.AbstractDiskArray{T,N}
    parent::A
    chunks::DiskArrays.GridChunks{N}
    reads::Dict{NTuple{N,UnitRange{Int}},Int}
end

CountingSource(parent::AbstractArray{T,N}, chunks) where {T,N} =
    CountingSource{T,N,typeof(parent)}(parent, chunks,
        Dict{NTuple{N,UnitRange{Int}},Int}())

Base.size(A::CountingSource) = size(A.parent)

DiskArrays.haschunks(::CountingSource) = DiskArrays.Chunked()
DiskArrays.eachchunk(A::CountingSource) = A.chunks

function DiskArrays.readblock!(A::CountingSource{T,N}, out,
    r::AbstractUnitRange...) where {T,N}
    key = ntuple(i -> UnitRange{Int}(r[i]), N)
    A.reads[key] = get(A.reads, key, 0) + 1
    out .= view(A.parent, r...)
    return out
end

DiskArrays.writeblock!(A::CountingSource, v, r::AbstractUnitRange...) =
    (view(A.parent, r...) .= v; v)

# --- a sampler that records how often it is asked -------------------------

"""
    CountingSites(inner, placed)

Sampler state that counts point locations and answers with `inner`'s stencil.
`placed` is therefore the number of `weightsat!` queries a build performed,
however many blocks came out of them, and the stencils are the ones the wrapped
sampler produces.
"""
struct CountingSites{S}
    inner::S
    placed::Threads.Atomic{Int}
end

function GR.weightsat!(row::GR.WeightRow,
    s::GR.Sampler{<:GR.RegridSpace,<:AbstractVector,<:CountingSites}, p)
    Threads.atomic_add!(s.state.placed, 1)
    return GR.weightsat!(row, s.state.inner, p)
end

countingsampler(space, placed) =
    GR.Sampler(space, GR.samplesites(space),
        CountingSites(GR.sampler(BarycentricPoint(), space), placed))

# --- fixture ---------------------------------------------------------------

centres(lo, hi, n) = range(lo + (hi - lo) / 2n, hi - (hi - lo) / 2n; length = n)

# A smooth field the regrid reproduces exactly, so no arm can be measuring a
# computation the others skip.
surface(lon, lat) = 2.0 + 0.01lon + 0.03lat + 0.0007lon * lat

function fixture()
    lons, lats = centres(-180.0, 180.0, NX), centres(-90.0, 90.0, NY)
    raw = Array{Float64}(undef, NX, NY)
    @inbounds for j in 1:NY, i in 1:NX
        raw[i, j] = surface(lons[i], lats[j])
    end
    source = CountingSource(raw, DiskArrays.GridChunks(raw, (CX, CY)))
    data = DD.DimArray(source, (X(lons), Y(lats)))
    src = RasterGrid(data)

    dlons, dlats = centres(DLON..., DNX), centres(DLAT..., DNY)
    dst = RasterGrid(DD.DimArray(zeros(DNX, DNY), (X(dlons), Y(dlats))))
    want = [surface(dlons[i], dlats[j]) for j in 1:DNY for i in 1:DNX]
    return (; src, dst, data, source, want)
end

newplan(f) = ChunkedPlan(BilinearPoint(), Weighted(0.5), f.dst, f.src;
    storage = PerChunk())

# The same fixture on the fused route: `BarycentricPoint` supplies a `sampler`,
# so its build unit is the destination tile and its reads are its manifest.
newtileplan(f) = ChunkedPlan(BarycentricPoint(), Weighted(0.5), f.dst, f.src;
    storage = PerChunk())

# --- measurement -----------------------------------------------------------

"""
    measure(label, f) -> (time, alloc, value)

Warm `f` once, take the minimum wall time of `NSAMPLES` calls, and allocate one
further call to price it. The returned value is the warm-up's, so nothing
measured is dead code.
"""
function measure(label, f)
    value = f()
    best = Inf
    for _ in 1:NSAMPLES
        GC.gc(false)
        best = min(best, @elapsed f())
    end
    alloc = @allocated f()
    @printf("%-34s %10.3f ms %16s B\n", label, best * 1e3, alloc)
    return (time = best, alloc = alloc, value = value)
end

# Bytes the tile's blocks occupy, by the plan's accounting and by walking them.
function weightbytes(plan)
    blocks = [e.block for e in values(plan.storage.blocks)]
    return (accounted = GR.storagebytes(plan.storage),
        summarysize = Base.summarysize(blocks),
        nblocks = GR.nblocks(plan.storage),
        nonzeros = sum(b -> count(!iszero, b.weights), blocks; init = 0),
        contributing = count(b -> any(!iszero, b.weights), blocks))
end

# The same accounting for a plan whose cache entry is one whole tile: its
# blocks are the tile's, and its manifest is the exact source-chunk list.
function tilebytes(plan)
    tiles = collect(values(plan.storage.tiles))
    blocks = [e.block for t in tiles for e in t.entries]
    return (accounted = GR.storagebytes(plan.storage),
        summarysize = Base.summarysize(blocks),
        ntiles = length(tiles),
        nblocks = length(blocks),
        nonzeros = sum(b -> count(!iszero, b.weights), blocks; init = 0),
        manifest = isempty(tiles) ? Int[] : tiles[1].sourcechunks)
end

"""
    provenance() -> NamedTuple

The machine state the timings belong to. Power mode is part of it: the same
code on the same machine runs at different clocks under low power, so a number
recorded here is comparable only against another carrying the same stamp.
"""
function provenance()
    shell(cmd) = try
        strip(read(cmd, String))
    catch
        ""
    end
    # macOS reports `powermode` (0 normal, 1 low, 2 high) or, on machines that
    # only offer the toggle, `lowpowermode`.
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
        power, samples = NSAMPLES)
end

function main()
    println(provenance())
    f = fixture()
    tile = GR.ownedindices(f.dst, 1)
    plan = newplan(f)
    candidates = collect(Int, GR.sourcesof(GR.dependencies(plan), 1))
    k = length(candidates)

    @printf("\nsource       %d x %d = %d cells, %d chunks of %d x %d\n",
        NX, NY, GR.ncells(f.src), GR.nchunks(f.src), CX, CY)
    @printf("destination  %d x %d = %d cells, %d tile\n",
        DNX, DNY, GR.ncells(f.dst), GR.nchunks(f.dst))
    @printf("candidates   %d of %d source chunks: %s\n\n",
        k, GR.nchunks(f.src), string(candidates))

    println("BilinearPoint, chunk-pair route")
    whole = 1:Int(GR.ncells(f.src))
    onepass = measure("one pass over the tile", function ()
        coo = WeightCOO(length(tile))
        buildweights!(coo, BilinearPoint(), f.dst, tile, f.src, whole)
        return length(coo.vals)
    end)

    # `buildblock` consults no storage, so each call is a fresh build of one
    # pair against the plan the executor would use.
    pairs = measure("pair loop ($k candidates)", function ()
        n = 0
        for s in candidates
            n += count(!iszero,
                GR.buildblock(plan, tile, GR.ownedindices(f.src, s)).weights)
        end
        return n
    end)

    # Plan construction — spaces, chunk index and dependency relation — reads no
    # data and builds no weights, but a cold read pays it, so it is priced apart
    # from the two arms that contain it.
    planbuild = measure("plan construction", () -> newplan(f))

    # A fresh plan per sample, so every sample builds the weights it applies.
    cold = measure("cold tile read", function ()
        p = newplan(f)
        A = GR.LazyRegridArray(f.data, p)
        return A[1:Int(GR.ncells(f.dst))]
    end)

    warm = measure("warm tile read", function ()
        A = GR.LazyRegridArray(f.data, plan)
        return A[1:Int(GR.ncells(f.dst))]
    end)

    cold.value ≈ f.want || error("the cold read did not reproduce the source surface")
    warm.value ≈ f.want || error("the warm read did not reproduce the source surface")

    # The fused route over the same tile. The sampler is prepared inside the
    # arm, so nothing the build needs is excluded from it; a plan prepares one
    # per read, not one per tile.
    println("\nBarycentricPoint, fused tile route")
    tilebuild = measure("tile weights, one pass", function ()
        tw = GR.tileweights(BarycentricPoint(), GR.TileCells(f.dst, tile), tile,
            f.src, GR.sampler(BarycentricPoint(), f.src))
        return sum(b -> count(!iszero, b.weights), tw.blocks; init = 0)
    end)

    tileplanbuild = measure("plan construction, tile route",
        () -> newtileplan(f))

    coldtile = measure("cold read, tile route", function ()
        p = newtileplan(f)
        A = GR.LazyRegridArray(f.data, p)
        return A[1:Int(GR.ncells(f.dst))]
    end)

    tileplan = newtileplan(f)
    warmtile = measure("warm read, tile route", function ()
        A = GR.LazyRegridArray(f.data, tileplan)
        return A[1:Int(GR.ncells(f.dst))]
    end)

    coldtile.value ≈ f.want ||
        error("the cold tile-route read did not reproduce the source surface")
    warmtile.value ≈ f.want ||
        error("the warm tile-route read did not reproduce the source surface")

    # One instrumented read, on a clean plan and a clean read counter, so the
    # residency and per-chunk numbers describe exactly one tile.
    empty!(f.source.reads)
    fresh = newplan(f)
    A = GR.LazyRegridArray(f.data, fresh)
    A[1:Int(GR.ncells(f.dst))]
    stats = GR.residency(A)
    bytes = weightbytes(fresh)
    counts = sort(collect(values(f.source.reads)))

    # The same instrumented read on the fused route, again on a clean plan and a
    # clean read counter.
    empty!(f.source.reads)
    freshtile = newtileplan(f)
    At = GR.LazyRegridArray(f.data, freshtile)
    At[1:Int(GR.ncells(f.dst))]
    tilestats = GR.residency(At)
    tbytes = tilebytes(freshtile)
    tilecounts = sort(collect(values(f.source.reads)))

    # Point locations, counted rather than asserted: one build of the same tile
    # with the same stencils, through a sampler that records every query.
    placed = Threads.Atomic{Int}(0)
    counted = GR.tileweights(BarycentricPoint(), GR.TileCells(f.dst, tile), tile,
        f.src, countingsampler(f.src, placed))
    counted.sourcechunks == tbytes.manifest ||
        error("the counted build and the plan's tile disagree on the manifest")

    println()
    @printf("%-46s %14d\n", "destination cells in the tile", length(tile))
    @printf("%-46s %14d\n", "candidate source chunks (k)", k)
    @printf("%-46s %14d\n", "point locations, one pass", length(tile))
    @printf("%-46s %14d\n", "point locations, pair loop", length(tile) * k)
    @printf("%-46s %14d\n", "stencil entries kept", onepass.value)
    @printf("%-46s %14.2f\n", "entries per destination cell",
        onepass.value / length(tile))
    @printf("%-46s %14.2fx\n", "pair loop / one pass", pairs.time / onepass.time)
    @printf("%-46s %14.2fx\n", "cold / warm tile read", cold.time / warm.time)
    @printf("%-46s %14.3f\n", "plan construction, ms", planbuild.time * 1e3)

    println()
    @printf("%-46s %14d\n", "weight blocks held", bytes.nblocks)
    @printf("%-46s %14d\n", "blocks carrying a nonzero", bytes.contributing)
    @printf("%-46s %14d\n", "weight nonzeros held", bytes.nonzeros)
    @printf("%-46s %14d\n", "weight bytes, plan accounting", bytes.accounted)
    @printf("%-46s %14d\n", "weight bytes, summarysize", bytes.summarysize)
    @printf("%-46s %14d\n", "weightbudget(plan.budget)", GR.weightbudget(fresh.budget))

    println()
    @printf("%-46s %14d\n", "source chunks read", length(counts))
    @printf("%-46s %14s\n", "reads per chunk (min, max)",
        string((minimum(counts), maximum(counts))))
    @printf("%-46s %14d\n", "readblock! calls", sum(counts))
    @printf("%-46s %14d\n", "source loads (residency)", stats.loads)
    @printf("%-46s %14d\n", "source cache hits (residency)", stats.hits)
    @printf("%-46s %14d\n", "peak source bytes (residency)", stats.peakbytes)

    println()
    @printf("%-46s %14d\n", "point locations, tile route", placed[])
    @printf("%-46s %14d\n", "stencil entries kept, tile route", tbytes.nonzeros)
    @printf("%-46s %14.2f\n", "entries per destination cell, tile route",
        tbytes.nonzeros / length(tile))
    @printf("%-46s %14.2fx\n", "pair loop / tile build", pairs.time / tilebuild.time)
    @printf("%-46s %14.2fx\n", "cold pair read / cold tile-route read",
        cold.time / coldtile.time)
    @printf("%-46s %14.2fx\n", "warm pair read / warm tile-route read",
        warm.time / warmtile.time)
    @printf("%-46s %14.3f\n", "plan construction, tile route, ms",
        tileplanbuild.time * 1e3)

    println()
    @printf("%-46s %14d\n", "tiles held", tbytes.ntiles)
    @printf("%-46s %14d\n", "weight blocks in the tile", tbytes.nblocks)
    @printf("%-46s %14d\n", "manifest length", length(tbytes.manifest))
    @printf("%-46s %14s\n", "manifest", string(tbytes.manifest))
    @printf("%-46s %14d\n", "candidate source chunks (k)", k)
    @printf("%-46s %14d\n", "weight nonzeros held, tile route", tbytes.nonzeros)
    @printf("%-46s %14d\n", "weight bytes, plan accounting, tile route",
        tbytes.accounted)
    @printf("%-46s %14d\n", "weight bytes, summarysize, tile route",
        tbytes.summarysize)

    println()
    @printf("%-46s %14d\n", "source chunks read, tile route", length(tilecounts))
    @printf("%-46s %14s\n", "reads per chunk, tile route (min, max)",
        string((minimum(tilecounts), maximum(tilecounts))))
    @printf("%-46s %14d\n", "readblock! calls, tile route", sum(tilecounts))
    @printf("%-46s %14d\n", "source loads, tile route", tilestats.loads)
    @printf("%-46s %14d\n", "source cache hits, tile route", tilestats.hits)
    @printf("%-46s %14d\n", "peak source bytes, tile route", tilestats.peakbytes)

    println()
    println("=== reading it ===")
    @printf("%d destination points, located %d times to keep %d stencil entries.\n",
        length(tile), length(tile) * k, onepass.value)
    @printf("One pass costs %.3f ms; the %d pair builds cost %.3f ms (%.2fx).\n",
        onepass.time * 1e3, k, pairs.time * 1e3, pairs.time / onepass.time)
    @printf("Weight construction is %.1f%% of a cold tile read.\n",
        100 * (cold.time - warm.time) / cold.time)
    @printf("Source chunks read: %d, of which %d carry a weight.\n",
        length(counts), bytes.contributing)
    @printf("The fused build locates each of those %d points once for the same %d\n",
        placed[], tbytes.nonzeros)
    @printf("entries: %.3f ms against the pair loop's %.3f ms (%.2fx), and a cold read\n",
        tilebuild.time * 1e3, pairs.time * 1e3, pairs.time / tilebuild.time)
    @printf("%.2fx the pair route's. Its manifest of %d chunks took %d reads.\n",
        cold.time / coldtile.time, length(tbytes.manifest), sum(tilecounts))
end

main()
