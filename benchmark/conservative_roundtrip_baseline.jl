# Exact reproduction of the P1/P4 raster workload used to expose the
# Conservative CSC -> COO -> CSC peak-memory round trip. Run modes in separate
# processes because `Sys.maxrss()` cannot be reset.
#
#     julia -t 8 --project=benchmark benchmark/conservative_roundtrip_baseline.jl direct
#     julia -t 8 --project=benchmark benchmark/conservative_roundtrip_baseline.jl chunked
#
# Baseline on 2026-08-22, Julia 1.12.6, 8 threads: both modes produce 7,120,800
# nonzeros, a 166,291,416-byte block, and sum to 4pi. Direct: 3.14 s,
# 12,333,000,896 bytes allocated, 2,640.2 MiB max RSS. Chunked: 3.16 s,
# 13,143,046,416 bytes allocated, 2,784.4 MiB max RSS. The current round trip
# therefore adds 810,045,520 allocated bytes and 144.2 MiB peak RSS.

import DimensionalData as DD
import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import SparseArrays

function rasterspace(nx, ny; T = Float64, chunks = nothing)
    dx, dy = 360 / nx, 180 / ny
    data = DD.DimArray(zeros(T, nx, ny),
        (DD.X(range(-180 + dx / 2, 180 - dx / 2; length = nx)),
         DD.Y(range(-90 + dy / 2, 90 - dy / 2; length = ny))))
    return chunks === nothing ? GR.RasterGrid(data) : GR.RasterGrid(data; chunks)
end

chunkranges(n, width) = [i:min(i + width - 1, n) for i in 1:width:n]

function build(mode, dst, src)
    if mode === :direct
        return GR.wholeblock(DGG.Conservative(), dst, src)
    end
    plan = GR.ChunkedPlan(DGG.Conservative(), GR.Weighted(0.5), dst, src)
    return GR.buildblock(plan, 1:GR.ncells(dst), 1:GR.ncells(src))
end

function main(mode::Symbol)
    mode in (:direct, :chunked) || error("mode must be `direct` or `chunked`")

    # Compile the selected construction path before measuring the full pair.
    build(mode, rasterspace(12, 6), rasterspace(36, 18))
    GC.gc()

    dst = rasterspace(360, 180;
        chunks = (chunkranges(360, 360), chunkranges(180, 20)))
    src = rasterspace(3600, 1800;
        chunks = (chunkranges(3600, 512), chunkranges(1800, 512)))
    GC.gc()
    timed = @timed build(mode, dst, src)
    block = timed.value
    maxrss = Int(Sys.maxrss())
    println((
        mode,
        destination_cells = GR.ncells(dst),
        source_cells = GR.ncells(src),
        nonzeros = SparseArrays.nnz(block.weights),
        weight_sum = sum(block.weights),
        seconds = timed.time,
        allocated = timed.bytes,
        block_summarysize = Base.summarysize(block),
        process_maxrss_bytes = maxrss,
        process_maxrss_mib = maxrss / 2^20,
    ))
end

main(Symbol(only(ARGS)))
