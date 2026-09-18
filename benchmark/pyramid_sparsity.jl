# Where each layout degenerates: the same level-13 cells, arranged from a solid
# rectangle to scattered points.
#
#   julia --project=benchmark benchmark/pyramid_sparsity.jl
#
# No network: the cell sets are geometry, and the values are random `Float32`
# so they do not compress and every difference below is the layout's own.
#
# # What the arithmetic says
#
# With `n` cells stored, `r` maximal runs of consecutive indices, `C` chunks
# touched, `W = 7^k` slots per chunk and a value type of `v` bytes, the four
# layouts lay down, before the compressor:
#
#   :cells dense     v·n + 8n          index is Θ(n) always -- 2x the data at Float32
#   :cells ranges    v·n + 16r         beats dense iff r < n/2, i.e. mean run > 2
#   :subzones        v·C'·capacity     and refuses unless every column is whole
#   :pyramid         v·C·W·7/6         INDEPENDENT of n; the ratio that matters
#                                      is the fill n/(C·W)
#
# So each fails differently. `dense` never blows up, it just never gets cheap.
# `ranges` is unbeatable on runs and becomes WORSE than dense once the mean run
# drops below two. `subzones` cannot express a partly covered subtree at all.
# The pyramid pays for chunks it touches whether or not they are full, so its
# cost is flat in `n` and its enemy is scatter, not size.
#
# # What actually happens
#
# Measured on `feat/pyramid-layout`, Julia 1.13.0, M-series macOS, 2026-09-18,
# IGEO7 level 13, default chunking (7^6). Megabytes on disk / files:
#
#   case                            cells    runs  chunks   fill    dense    ranges     pyramid
#   tile, 1x1 degree             16181892    9909     170  80.9%  58.1/35   57.5/19    67.1/227
#   straddling a base cell       16177029   10502     167  82.3%      -     57.5/19    67.1/233
#   1x1 box at the north pole      210558    5239      20   8.9%      -      0.8/3      1.1/42
#   thin strip, 1 x 0.01 deg       165737    5250      14  10.1%   0.6/3     0.6/3      0.9/38
#   every 7th cell of the tile    2311699 2311699     170  11.6%   8.3/7     8.5/5     35.2/227
#   every 49th cell                330243  330243     169   1.7%   1.2/3     1.2/3     10.2/226
#   1000 scattered in the tile       1000    1000     157   0.0%   0.0/3     0.0/3      0.5/212
#   1000 scattered on the earth      1000    1000    1000   0.0%   0.0/3     0.0/3    10.6/5411
#
# Five things that arithmetic alone gets wrong:
#
#  1. A RECTANGLE IS FINE WHEREVER YOU PUT IT. Straddling a base-cell boundary
#     costs nothing measurable -- 82.3% fill against 80.9% mid-face. At the
#     chunk root level a one-degree box already spans about 170 cells, so it
#     crosses subtree boundaries wherever it sits; the level-0 ones are not
#     special. What DOES hurt is a shape with a large perimeter for its area:
#     the polar box and the thin strip fall to 9-10% fill.
#  2. COMPRESSION RESCUES MOST OF THE AMPLIFICATION. The earth-scatter case is
#     549 MB of logical slots and 10.6 MB on disk, because solid fill is a
#     constant and Blosc eats it.
#  3. WHAT IT DOES NOT RESCUE IS INTERLEAVED SPARSITY. Every 7th cell is 4.2x
#     dense on disk, every 49th is 8.5x, because a chunk of alternating NaN and
#     value has no runs to squeeze -- unlike a chunk that is all fill.
#  4. THE REAL WALL IS FILE COUNT, not bytes: 1000 cells scattered over the
#     earth cost 5411 files, five per stored value, one object-store request
#     each. That is the number to watch, and the coarse levels contribute most
#     of it because scattered cells stay scattered all the way up.
#  5. RANGES REALLY DOES CROSS OVER. Every 7th cell makes every run length one,
#     and `ranges` lands at 8.5 MB against dense's 8.3 -- the predicted flip at
#     mean run length two, muted on disk because both halves compress.
#
# The chunk exponent trades bytes against files, and which way depends on
# whether the sparsity is finer or coarser than a chunk:
#
#   case                    7^2            7^3            7^4            7^6
#   the whole tile     82.0MB/387441  73.1MB/55914   69.8MB/8213    67.1MB/227
#   every 49th cell    32.3MB/386000  19.9MB/55758   14.3MB/8187    10.2MB/226
#   1000 over earth     0.7MB/9363     1.0MB/8375     1.1MB/7387    10.6MB/5411
#
# Thinning touches every chunk whatever the chunk size, so a bigger chunk is
# strictly better there. Isolated points waste one chunk each, so a smaller
# chunk is better on bytes -- and worse on files. Tune it to the sparsity you
# have, and if that sparsity is scatter rather than shape, use `:cells`.

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import Extents
using Zarr
using Printf
using Random

const SYS = DGG.IGeo7System()
const L = 13
const K = 6
const W = 7^K
const GRID = DGG.levelgrid(SYS, L)

bytes(d) = sum(filesize(joinpath(r, f)) for (r, _, fs) in walkdir(d) for f in fs; init=0)
files(d) = sum(count(f -> !startswith(f, "."), fs) for (r, _, fs) in walkdir(d); init=0)

box(x0, x1, y0, y1) = Extents.Extent(X=(x0, x1), Y=(y0, y1))

"The level-`L` cells covering a lon/lat box."
rect(e) = DGG.CellVector(DGG.query(SYS, DGG.MultiOrderCoverage(e); level=L); level=L)

indices(cv) = reduce(vcat,
    [collect(lo:hi) for (lo, hi) in DGG.Engine.intervals(DGG.region(cv).windows)])
fromindices(ix) = DGG.CellVector(SYS, L, [DGG.cellindex(GRID, i) for i in ix])

"""
    shape(cv) -> (n, runs, chunks, fill)

What a cell set looks like to a layout: how many cells, how many maximal runs of
consecutive indices (what `:cells ranges` stores), how many level-`L` chunks it
touches, and how full it leaves them (what `:pyramid` pays for).
"""
function shape(cv)
    ivs = DGG.Engine.intervals(DGG.region(cv).windows)
    n = sum(hi - lo + 1 for (lo, hi) in ivs)
    chunks = Set{Int}()
    for (lo, hi) in ivs
        p = lo
        while p <= hi
            a = DGG.ancestor(SYS, DGG.cellindex(GRID, p), L - K)
            push!(chunks, DGG.slotindex(SYS, a))
            p = min(hi, last(DGG.descendant_range(SYS, a, L))) + 1
        end
    end
    return n, length(ivs), length(chunks), n / (length(chunks) * W)
end

function main()
    Random.seed!(1)
    tile = rect(box(10.0, 11.0, 46.0, 47.0))
    base = indices(tile)
    # The Z7 base cell is the top four bits of an id, and its boundaries are not
    # where the twelve level-0 cell CENTRES would put them; 54.8613E at this
    # latitude is a real one, found by scanning those bits.
    cases = [
        "tile, 1x1 degree" => tile,
        "straddling a base cell" => rect(box(54.3613, 55.3613, 46.0, 47.0)),
        "1x1 box at the north pole" => rect(box(-0.5, 0.5, 89.0, 90.0)),
        "thin strip, 1 x 0.01 deg" => rect(box(10.0, 11.0, 46.50, 46.51)),
        "every 7th cell of the tile" => fromindices(base[1:7:end]),
        "every 49th cell" => fromindices(base[1:49:end]),
        "1000 scattered in the tile" => fromindices(sort!(unique(rand(base, 1000)))),
        "1000 scattered on the earth" =>
            fromindices(sort!(unique(rand(1:DGG.ncells(SYS, L), 1000)))),
    ]

    @printf("%-28s %9s %9s %7s %6s | %11s %11s %13s\n", "case", "cells", "runs",
        "chunks", "fill", "dense", "ranges", "pyramid")
    println("-"^104)
    for (label, cv) in cases
        n, runs, chunks, fill = shape(cv)
        A = DD.DimArray(rand(Float32, n), (DGG.Cells(DGG.CellLookup(cv)),);
            name=:elevation)
        cols = String[]
        for kw in ((; encoding=:dense), (; encoding=:ranges), (; layout=:pyramid))
            d = joinpath(mktempdir(), "s.zarr")
            try
                DGG.dggwrite(d, A; kw...)
                push!(cols, @sprintf("%7.1fMB/%d", bytes(d) / 1e6, files(d)))
            catch err
                push!(cols, "    refused")
            end
        end
        @printf("%-28s %9d %9d %7d %5.1f%% | %11s %11s %13s\n",
            label, n, runs, chunks, 100fill, cols[1], cols[2], cols[3])
    end
    return nothing
end

main()
