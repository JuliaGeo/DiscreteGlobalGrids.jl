# Production-shaped Conservative weight benchmark and CSC -> COO -> CSC memory
# baseline. Run each mode in a fresh process because `Sys.maxrss()` is a
# process-lifetime high-water mark.
#
#     julia -t 8 --project=benchmark benchmark/conservative_block_baseline.jl direct
#     julia -t 8 --project=benchmark benchmark/conservative_block_baseline.jl chunked
#
# Baseline on 2026-08-22, Julia 1.12.6, 8 threads: both modes produce 487,174
# nonzeros and sum to 4pi. Direct: 3.36 s, 754.5 MiB max RSS. Chunked:
# 2.31 s, 755.4 MiB max RSS. This workload gates the production tree path;
# `conservative_roundtrip_baseline.jl` is the larger sparse-copy stress case.

import ConservativeRegridding
import DimensionalData as DD
import DiscreteGlobalGrids as DGG
import GlobalRegridding as GR
import SparseArrays

function spaces()
    dgg = DGG.levelgrid(DGG.IGeo7System(), 5)
    dst = DGG.DGGSpace(dgg)
    raster = DD.DimArray(zeros(Float32, 360, 180),
        (DD.X(-179.5:1.0:179.5), DD.Y(-89.5:1.0:89.5)))
    src = GR.RasterGrid(raster)
    return dst, src
end

function main(mode::Symbol)
    mode in (:direct, :chunked) || error("mode must be `direct` or `chunked`")
    dst, src = spaces()
    GC.gc()
    timed = if mode === :direct
        @timed GR.wholeblock(DGG.Conservative(), dst, src)
    else
        plan = GR.ChunkedPlan(DGG.Conservative(), GR.Weighted(0.5), dst, src)
        @timed GR.buildblock(plan, 1:GR.ncells(dst), 1:GR.ncells(src))
    end
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
