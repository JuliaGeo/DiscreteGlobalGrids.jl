# Raster range-extent and small-spatial-chunk traversal benchmark.
#
#     nice -n 10 julia --project=benchmark --threads=8 \
#         benchmark/raster_extent_baseline.jl
#
# The first comparison is the A4 retention gate for the precisely geographic
# analytic specialization: 2,000 deterministic rectangles on a 720x360 global
# grid, against ConservativeRegridding's general perimeter walk over the same
# RasterGridView. The second measurement exercises candidate discovery with
# 2x2-cell DiskArrays ownership chunks.
#
# A4 checkpoint on 2026-08-22, Julia 1.12.6, 8 threads: geographic analytic
# 0.861 ms, zero allocations; CR perimeter 18.884 ms, zero allocations (21.9x
# slower). The 64,800-chunk / 128-query traversal took 5.732 ms and allocated
# 58,288 bytes in 107 allocations.

import BenchmarkTools
import ConservativeRegridding as CR
import DimensionalData as DD
import GeometryOps as GO
import GlobalRegridding as GR

function global_raster(nx, ny; chunkwidth = nothing)
    dx, dy = 360 / nx, 180 / ny
    data = DD.DimArray(zeros(Float32, nx, ny),
        (DD.X(range(-180 + dx / 2, 180 - dx / 2; length = nx)),
         DD.Y(range(-90 + dy / 2, 90 - dy / 2; length = ny))))
    chunkwidth === nothing && return GR.RasterGrid(data)
    chunks(n, width) = [i:min(i + width - 1, n) for i in 1:width:n]
    return GR.RasterGrid(data;
        chunks = (chunks(nx, chunkwidth), chunks(ny, chunkwidth)))
end

function rectangles(nx, ny, n)
    return map(1:n) do k
        ilo = mod(37k - 1, nx) + 1
        jlo = mod(53k - 1, ny) + 1
        ihi = min(nx, ilo + mod(29k, 181))
        jhi = min(ny, jlo + mod(31k, 91))
        (ilo:ihi, jlo:jhi)
    end
end

function extent_checksum(f, grid, rects)
    total = 0.0
    for (irange, jrange) in rects
        total += f(grid, irange, jrange).radius
    end
    return total
end

geographic_extent(grid, irange, jrange) =
    CR.Trees.cell_range_extent(grid, irange, jrange)
perimeter_extent(grid, irange, jrange) =
    GR._curvilinear_range_extent(grid, irange, jrange)

function query_batch(index, caps)
    out = Int[]
    checksum = 0
    for cap in caps
        GR.candidatechunks!(out, index, cap)
        checksum += length(out)
    end
    return checksum
end

function estimate(f; seconds = 5)
    trial = BenchmarkTools.run(BenchmarkTools.@benchmarkable($f());
        samples = 10, evals = 1, seconds)
    e = BenchmarkTools.median(trial)
    return (; milliseconds = e.time / 1e6, allocated = e.memory, allocations = e.allocs)
end

function main()
    geographic = global_raster(720, 360)
    grid = GR.RasterGridView(geographic)
    rects = rectangles(720, 360, 2_000)
    analytic = () -> extent_checksum(geographic_extent, grid, rects)
    perimeter = () -> extent_checksum(perimeter_extent, grid, rects)
    analytic() == perimeter() && error(
        "benchmark paths unexpectedly produced identical aggregate cap radii")

    smallchunks = global_raster(720, 360; chunkwidth = 2)
    index = GR.chunkindex(smallchunks)
    to_sphere = GO.UnitSpherical.UnitSphereFromGeographic()
    caps = [GO.UnitSpherical.SphericalCap(
        to_sphere(((mod(41k, 720) + 0.5) / 2 - 180,
                   (mod(43k, 360) + 0.5) / 2 - 90)), deg2rad(0.75))
        for k in 1:128]
    traversal = () -> query_batch(index, caps)
    traversal()

    GC.gc()
    println((
        rectangles = length(rects),
        geographic_analytic = estimate(analytic),
        cr_perimeter = estimate(perimeter),
        small_chunk_shape = (2, 2),
        small_chunk_count = GR.nchunks(smallchunks),
        traversal_queries = length(caps),
        small_chunk_traversal = estimate(traversal),
    ))
end

main()
