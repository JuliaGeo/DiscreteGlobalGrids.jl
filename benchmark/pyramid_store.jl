# Benchmark the pyramid layout: writing a partially covered earth, opening it,
# and the two reads a viewer does per frame.
#
#   julia --project=benchmark benchmark/pyramid_store.jl
#
# The case is the one the layout exists for: half an earth, so half the chunks
# are never written and the coarse levels are half fill. The synthetic field is
# random, which is the worst case for the compressor and makes the on-disk size
# comparable with the raw bytes.
#
# Output on `feat/pyramid-layout`, Julia 1.13.0, M-series macOS, 2026-09-17,
# level 7 (4.1 M cells of a 8.2 M-cell level), default chunking (7^6):
#
#   write 4.1 M cells              0.22 s
#   on disk                       17.0 MB   (16.5 MB of Float32 leaf data)
#   open the store                 0.0010 s
#   collect the whole leaf level   0.106 s
#   cellvalues, 2 000 cells @ L4   0.0001 s
#   cellvalues, 400 000 cells      0.016 s
#   holdsdata, one cell            38 us
#
# The last two are the frame path, and the gap between them is the point: a
# descent that asks `holdsdata` per candidate pays a chunk decompression per
# call, and one that asks `cellvalues` for the whole level's candidates pays one
# per chunk. Batch the frame.
#
# The real case, from `docs/src/tutorials/hydrology.jl` — one 1x1 degree
# Copernicus DEM tile at 30 m regridded onto IGEO7 level 13, same machine:
#
#   16 172 725 cells written        1.07 s
#   on disk                        55.7 MB  in 227 chunk files
#   read the whole region back      2.07 s   (bit-identical to what went in)
#   open the store                  0.0014 s
#
# Level 13 addresses 968 890 104 072 cells in 9 882 516 chunks. The store has
# 170 of those chunks. Sparsity is not a special case here; it is the layout.

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
using Zarr
using Printf
using Statistics

const SYS = DGG.IGeo7System()
const LEVEL = parse(Int, get(ENV, "DGG_PYRAMID_BENCH_LEVEL", "7"))
const COVERAGE = 6                      # of the twelve base subtrees

"Half an earth at `LEVEL`, as a cell vector."
function halfearth(level)
    roots = collect(DGG.rootcells(SYS))[1:COVERAGE]
    cells = sort!(reduce(vcat,
        [collect(DGG.descendants(SYS, r, level)) for r in roots]))
    return DGG.CellVector(SYS, level, cells)
end

report(label, t; unit="s") = @printf("  %-32s %8.4f %s\n", label, t, unit)

function main()
    cv = halfearth(LEVEL)
    n = length(cv)
    @printf("half an earth at IGEO7 level %d: %d of %d cells\n",
        LEVEL, n, DGG.ncells(SYS, LEVEL))
    cube = DD.DimArray(rand(Float32, n),
        (DGG.Cells(DGG.CellLookup(cv)),); name=:elevation)

    # Warm every path once; the first call through any of them is compilation.
    warm = mktempdir()
    DGG.dggwrite(warm, cube; layout=:pyramid)
    warmed = DGG.dggread(warm)[:elevation]
    collect(warmed[0][:elevation])
    DGG.holdsdata(warmed, first(DGG.rootcells(SYS)))
    DGG.cellvalues(warmed, LEVEL, [DGG.cellindex(DGG.levelgrid(SYS, LEVEL), 1)])

    dest = mktempdir()
    report("write $(round(n / 1e6; digits=1)) M cells",
        (@timed DGG.dggwrite(dest, cube; layout=:pyramid)).time)
    bytes = sum(filesize(joinpath(r, f)) for (r, _, fs) in walkdir(dest) for f in fs)
    @printf("  %-32s %8.1f MB  (%.1f MB of leaf data)\n", "on disk",
        bytes / 1e6, 4n / 1e6)

    report("open the store", (@timed DGG.dggread(dest)).time)
    pyr = DGG.dggread(dest)[:elevation]

    report("collect the whole leaf level",
        (@timed collect(pyr[LEVEL][:elevation])).time)

    # A frame's worth: the cells a zoomed-out view draws, and the cells a
    # zoomed-in one does.
    coarse = [DGG.cellindex(DGG.levelgrid(SYS, 4), i)
              for i in 1:min(2000, DGG.ncells(SYS, 4))]
    report("cellvalues, $(length(coarse)) cells @ L4",
        (@timed DGG.cellvalues(pyr, 4, coarse)).time)

    fine = [DGG.cellindex(DGG.levelgrid(SYS, LEVEL), i)
            for i in 1:min(400_000, DGG.ncells(SYS, LEVEL))]
    report("cellvalues, $(length(fine)) cells @ L$LEVEL",
        (@timed DGG.cellvalues(pyr, LEVEL, fine)).time)

    t = (@timed [DGG.holdsdata(pyr, c) for c in coarse]).time
    report("holdsdata, one cell", t / length(coarse) * 1e6; unit="us")
    return nothing
end

main()
