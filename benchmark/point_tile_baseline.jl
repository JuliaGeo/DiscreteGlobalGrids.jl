# What one large destination tile costs a point method today.
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
#   arm                              time        allocated
#   one pass over the tile      72.950 ms       96,081,312 B
#   pair loop (4 candidates)   320.903 ms      244,901,440 B
#   plan construction            0.058 ms            4,096 B
#   cold tile read             308.418 ms      252,154,736 B
#   warm tile read               2.834 ms       19,275,648 B
#
#   375,000 destination cells, 4 candidate source chunks, all four contributing
#   1,500,000 stencil entries, 4.00 per destination cell
#   1,500,000 point locations performed where 375,000 would do
#   weight bytes 40,147,488 (plan accounting) / 28,147,912 (summarysize)
#   4 source chunks read, each once; peak source residency 4,147,200 B
#
# The pair loop is 4.40x one pass at k = 4: the destination pass is repeated per
# candidate and nothing else about the build changes. Weight construction is
# 99.1% of a cold tile read, so this repetition, not source I/O or the weight
# application, is what a large point tile costs. Source reads are already exact
# here — four chunks, one read each — so what the fused build has to remove is
# the repeated point location, not a read.
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

    # One instrumented read, on a clean plan and a clean read counter, so the
    # residency and per-chunk numbers describe exactly one tile.
    empty!(f.source.reads)
    fresh = newplan(f)
    A = GR.LazyRegridArray(f.data, fresh)
    A[1:Int(GR.ncells(f.dst))]
    stats = GR.residency(A)
    bytes = weightbytes(fresh)
    counts = sort(collect(values(f.source.reads)))

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
    @printf("%-46s %14d\n", "source chunks read", length(f.source.reads))
    @printf("%-46s %14s\n", "reads per chunk (min, max)",
        string((minimum(counts), maximum(counts))))
    @printf("%-46s %14d\n", "readblock! calls", sum(counts))
    @printf("%-46s %14d\n", "source loads (residency)", stats.loads)
    @printf("%-46s %14d\n", "source cache hits (residency)", stats.hits)
    @printf("%-46s %14d\n", "peak source bytes (residency)", stats.peakbytes)

    println()
    println("=== reading it ===")
    @printf("%d destination points, located %d times to keep %d stencil entries.\n",
        length(tile), length(tile) * k, onepass.value)
    @printf("One pass costs %.3f ms; the %d pair builds cost %.3f ms (%.2fx).\n",
        onepass.time * 1e3, k, pairs.time * 1e3, pairs.time / onepass.time)
    @printf("Weight construction is %.1f%% of a cold tile read.\n",
        100 * (cold.time - warm.time) / cold.time)
    @printf("Source chunks read: %d, of which %d carry a weight.\n",
        length(f.source.reads), bytes.contributing)
end

main()
