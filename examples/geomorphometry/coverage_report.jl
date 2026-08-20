#!/usr/bin/env julia
# Count the geometry labels reached by the deterministic sweep.
#
#   julia --project=test examples/geomorphometry/coverage_report.jl
#
# The case set includes every level-0 and level-1 root, depths 0–3, and both
# connectivity modes.

include("Harness.jl")
using .Harness
using .Harness: modal_degree, pole_positions, nrange
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: Vertex, Edge
using Printf

const MAXCELLS = 30_000

function main()
    println("### per-system geometry inventory")
    for sys in DGG.systems(), conn in (Vertex(),)
        k = gridctx(sys, 3, conn)
        degs = [length(nrange(k, p)) for p in 1:length(k)]
        m = modal_degree(k)
        np, sp = pole_positions(k)
        @printf("  %-16s L3 %5d cells  degrees %s (modal %d)  low=%d high=%d  poles at %d,%d\n",
            nameof(typeof(sys)), length(k), string(sort(unique(degs))), m,
            count(<(m), degs), count(>(m), degs), np, sp)
    end

    println()
    println("### sweep case labels (all roots at L0/L1, depths 0..3, both connectivities)")
    total = Dict{Symbol,Int}()
    for sys in DGG.systems()
        counts = Dict{Symbol,Int}()
        n = 0
        for rootlevel in 0:1
            DGG.ncells(sys, rootlevel) > 200 && continue
            groot = DGG.levelgrid(sys, rootlevel)
            for depth in 0:3
                target = rootlevel + depth
                (target > DGG.maxlevel(sys) || DGG.ncells(sys, target) > MAXCELLS) && continue
                for conn in (Vertex(), Edge()), p in 1:DGG.ncells(groot)
                    root = DGG.cellindex(groot, p)
                    n += 1
                    for t in classify_root(sys, root, target, conn)
                        counts[t] = get(counts, t, 0) + 1
                        total[t] = get(total, t, 0) + 1
                    end
                end
            end
        end
        @printf("  %-16s %5d cases: %s\n", nameof(typeof(sys)), n,
            join(["$k=$v" for (k, v) in sort(collect(counts); by = first)], " "))
    end
    println("  TOTAL: ", join(["$k=$v" for (k, v) in sort(collect(total); by = first)], " "))
end

main()
