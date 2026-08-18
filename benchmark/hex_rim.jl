# Benchmark the hexagonal subtree automata on H3 and IGeo7: rim, interior, halo.
#
#   julia --project=test benchmark/hex_rim.jl
#
# Each cell reports the minimum wall time over nine runs and the bytes one run
# allocates, both after a warmup. The engines are consumed by a counting loop,
# so the only allocations reported are the walk's own.
#
# Roots are one hexagon and one pentagon per system at level 1, walked to a mid
# depth and a deep one. The interior walk is `O(7^d)` where the rim is `O(3^d)`,
# so it uses shallower depths to keep the run under a minute.
#
# This script exists because the two hexagonal automata — H3's `_h3_border_step`
# stack walk and IGeo7's `_border_step` stack walk — are the same automaton with
# different digit tables, and unifying them behind one `HexRimEngine{C}` was
# tried and measured. Minimum wall time in microseconds, Julia 1.12.6, M-series
# macOS, 2026-08-17, eight interleaved A/B rounds against `f4bce1d`, per-cell
# minimum. Allocations were 0 B in every cell of both columns.
#
#                            separate   unified
#   H3       hex rim d=5        34.75     15.12   -56.5%
#   H3       hex rim d=9      3109.12   1211.96   -61.0%
#   H3       hex interior d=4   66.79     25.96   -61.1%
#   H3       hex interior d=6 3626.96   1292.83   -64.4%
#   H3       hex halo d=4       79.71     70.71   -11.3%
#   H3       hex halo d=6      467.29    392.71   -16.0%
#   H3       pent rim d=5       28.83     11.25   -61.0%
#   H3       pent rim d=9     2622.83   1041.71   -60.3%
#   H3       pent interior d=4  55.00     21.71   -60.5%
#   H3       pent interior d=6 2994.17  1078.67   -64.0%
#   H3       pent halo d=4      65.42     60.29    -7.8%
#   H3       pent halo d=6     340.12    295.96   -13.0%
#   IGeo7    hex rim d=5        16.92     14.88   -12.1%
#   IGeo7    hex rim d=9      1392.88   1389.83    -0.2%
#   IGeo7    hex interior d=4   27.42     30.08    +9.7%
#   IGeo7    hex interior d=6 1497.17   1550.08    +3.5%
#   IGeo7    hex halo d=4       58.75     60.33    +2.7%
#   IGeo7    hex halo d=6      351.12    354.17    +0.9%
#   IGeo7    pent rim d=5       14.21     12.42   -12.6%
#   IGeo7    pent rim d=9     1108.38   1154.08    +4.1%
#   IGeo7    pent interior d=4  23.25     24.79    +6.6%
#   IGeo7    pent interior d=6 1186.50  1309.54   +10.4%
#   IGeo7    pent halo d=4      46.25     45.58    -1.4%
#   IGeo7    pent halo d=6     249.12    244.04    -2.0%
#
# The unification was reverted on the IGeo7 interior regressions. It was not
# algorithmic: the shared and separate walks compile to the same typed IR with
# no dynamic dispatch, and the gap closes under `--pkgimages=no` (1372 vs 1390)
# or `JULIA_CPU_TARGET=native` (1372 vs 1340), appearing only in the default
# multiversioned pkgimage at frame-stack capacity 20.
#
# The H3 column points at a standing cost that outlives the revert: H3's frame
# carries a fifth `full::Bool` field, which makes `H3Frame` five bytes where
# `Z7Frame` is four, and `small_push`/`small_setlast` rebuild the whole tuple on
# every descent. `full` is redundant with `L == 0`, and the unified engine that
# won those rows dropped it; whether dropping it from `H3Frame` alone recovers
# the same margin is unmeasured.

import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: Vertex, levelgrid, cellindex, level

const RUNS = 9

const H3 = DGG.H3
const IGeo7 = DGG.IGeo7

function drain(e)
    n = 0
    for _ in e
        n += 1
    end
    return n
end

rim(sys, c, target) = drain(DGG.rim_engine(sys, c, target, Vertex()))
interior(sys, c, target) = drain(DGG.interior_engine(sys, c, target, Vertex()))
halo(sys, c, target) = drain(DGG.halo_engine(sys, c, target, Vertex()))

function measure(f, sys, c, target)
    n = f(sys, c, target)
    best = Inf
    for _ in 1:RUNS
        t = @elapsed f(sys, c, target)
        best = min(best, t)
    end
    return (n, best, @allocated f(sys, c, target))
end

function report(label, f, sys, c, target)
    n, t, b = measure(f, sys, c, target)
    println(rpad(label, 34), lpad(n, 9), " cells  ",
        lpad(round(t * 1e6; digits = 2), 11), " us  ", lpad(b, 8), " B")
end

# A hexagon and a pentagon root, so the deleted-digit path is measured too.
function roots(sys, pentagon)
    grid = levelgrid(sys, 1)
    hex = cellindex(grid, DGG.ncells(grid) ÷ 2)
    return (("hex", hex), ("pentagon", pentagon))
end

function sweep(name, sys, pentagon, rimdepths, interiordepths, halodepths)
    for (kind, c) in roots(sys, pentagon)
        base = level(c)
        for d in rimdepths
            report("$name $kind rim d=$d", rim, sys, c, base + d)
        end
        for d in interiordepths
            report("$name $kind interior d=$d", interior, sys, c, base + d)
        end
        for d in halodepths
            report("$name $kind halo d=$d", halo, sys, c, base + d)
        end
    end
end

const H3_PENTAGON = H3.H3Cell(first(H3.H3Native.get_pentagons(1)))
const Z7_PENTAGON = IGeo7.Z7Cell(IGeo7.z7_from_string("000"))

sweep("H3", H3.H3System(), H3_PENTAGON, (5, 9), (4, 6), (4, 6))
sweep("IGeo7", IGeo7.IGeo7System(), Z7_PENTAGON, (5, 9), (4, 6), (4, 6))
