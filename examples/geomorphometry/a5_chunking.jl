#!/usr/bin/env julia
# Measure subtree chunk addressing and halo construction on A5.
#
#   julia --project=test examples/geomorphometry/a5_chunking.jl
#
# The report separates two operations:
#
#   1. `descendant_range` does not exist for A5 (`has_sorted_subtrees` is
#      false), so there is no O(1) way to name the chunk's position block. The
#      fallback materializes `descendants`, takes their position hull, and
#      verifies that the hull is contiguous, using O(subtree) time and memory.
#   2. A5 has no specialised halo engine, so `halo` on a subtree runs the
#      generic outside-first geometry walk.

include("SphericalTerrain.jl")
using .SphericalTerrain: chunk_range
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: Vertex
using Printf

function contiguity(sys, rootlevel, target, nroots)
    g = DGG.levelgrid(sys, target)
    groot = DGG.levelgrid(sys, rootlevel)
    n = DGG.ncells(groot)
    bad = 0
    for i in 1:nroots
        p = 1 + (i - 1) * max(1, n ÷ nroots)
        p > n && break
        ps = sort!([DGG.cellposition(g, c)
                    for c in DGG.descendants(sys, DGG.cellindex(groot, p), target)])
        ps == first(ps):last(ps) || (bad += 1)
    end
    return bad
end

function main()
    println("### is a subtree a contiguous position block on A5?")
    for (rl, tl, nr) in ((0, 6, 12), (1, 7, 12), (2, 8, 10), (3, 9, 6), (4, 10, 4))
        b = contiguity(DGG.A5System(), rl, tl, nr)
        @printf("  L%d roots -> L%-2d : %d of %d NOT contiguous\n", rl, tl, b, nr)
    end

    println()
    println("### cost of naming the chunk's position block")
    for (sys, rl, tl) in ((DGG.HEALPixSystem(), 4, 12), (DGG.S2System(), 4, 12),
                          (DGG.ISEA4RSystem(), 4, 12), (DGG.A5System(), 4, 12))
        root = DGG.cellindex(DGG.levelgrid(sys, rl), 3)
        nm = string(nameof(typeof(sys)))
        if DGG.has_sorted_subtrees(sys)
            DGG.descendant_range(sys, root, tl)
            t = @elapsed r = DGG.descendant_range(sys, root, tl)
            a = @allocated DGG.descendant_range(sys, root, tl)
            @printf("  %-16s descendant_range L%d->L%d: %d cells, %.6f s, %d B\n",
                nm, rl, tl, length(r), t, a)
        else
            @printf("  %-16s descendant_range: NO METHOD (has_sorted_subtrees = false)\n", nm)
            for t2 in (8, 10, 12)
                t2 > tl && continue
                el = @elapsed r = chunk_range(sys, root, t2)
                al = @allocated chunk_range(sys, root, t2)
                @printf("      fallback hull-of-descendants L%d->L%-2d: %d cells, %.3f s, %.1f MB\n",
                    rl, t2, length(r), el, al / 2^20)
            end
        end
    end

    println()
    println("### cost of the margin: halo(subtree), one root per system")
    for (sys, rl) in ((DGG.HEALPixSystem(), 6), (DGG.S2System(), 6),
                      (DGG.ISEA4RSystem(), 6), (DGG.IGeo7System(), 6),
                      (DGG.H3System(), 6), (DGG.A5System(), 6))
        root = DGG.cellat(DGG.levelgrid(sys, rl), 12.3, 41.7)
        nm = string(nameof(typeof(sys)))
        print(rpad("  " * nm, 18))
        for d in 1:6
            tl = rl + d
            tl > DGG.maxlevel(sys) && continue
            t = @elapsed h = collect(DGG.halo(DGG.subtree(sys, root, tl)))
            @printf(" d%d:%d/%.3fs", d, length(h), t)
            t > 60 && (print(" [stopped: over a minute]"); break)
        end
        println()
    end
end

main()
