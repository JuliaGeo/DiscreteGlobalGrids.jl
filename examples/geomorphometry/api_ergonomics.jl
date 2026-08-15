#!/usr/bin/env julia
# Report iterator sizing, ID-to-position conversion, and wider-halo costs.
#
#   julia --project=test examples/geomorphometry/api_ergonomics.jl

import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: Vertex, Edge
using Printf

function iterator_size()
    println("### 1. IteratorSize / length: HasLength only at depth 0")
    for sys in DGG.systems()
        r = DGG.cellindex(DGG.levelgrid(sys, 2), 1)
        line = String[]
        for l in (2, 3, 6)
            it = DGG.SubtreeHaloIterator(sys, r, l)
            len = try string(length(it)) catch e; string(nameof(typeof(e))) end
            push!(line, "L$l:$(Base.IteratorSize(typeof(it)))/$len")
        end
        @printf("  %-16s %s\n", nameof(typeof(sys)), join(line, "  "))
    end
    println("  => a caller cannot pre-size a read buffer for any real chunk.")
end

function id_to_position()
    println()
    println("### 2. the halo yields cell ids; arrays want positions")
    for (sys, rl, tl) in ((DGG.IGeo7System(), 4, 10), (DGG.H3System(), 4, 10),
                          (DGG.HEALPixSystem(), 4, 12), (DGG.S2System(), 4, 12),
                          (DGG.ISEA4RSystem(), 4, 12))
        g = DGG.levelgrid(sys, tl)
        root = DGG.cellindex(DGG.levelgrid(sys, rl), 3)
        DGG.subtree_halo(sys, root, tl)
        tw = @elapsed h = DGG.subtree_halo(sys, root, tl)
        DGG.cellposition.(Ref(g), h)
        tp = @elapsed DGG.cellposition.(Ref(g), h)
        ap = @allocated DGG.cellposition.(Ref(g), h)
        @printf("  %-16s halo=%4d  walk %8.1f us   id->pos %6.1f us (%4.1f%%, %d B)\n",
            nameof(typeof(sys)), length(h), tw * 1e6, tp * 1e6,
            100 * tp / (tw + tp), ap)
    end
    println("  => cheap, so ids are defensible on cost; the friction is that every")
    println("     caller writes the same cellposition broadcast and its own sorted map.")
end

function wider_k()
    println()
    println("### 3. there is no k>1 halo; the workaround leaves the fast engine")
    for sys in (DGG.IGeo7System(), DGG.H3System(), DGG.HEALPixSystem(),
                DGG.S2System(), DGG.ISEA4RSystem())
        tl = 6
        g = DGG.levelgrid(sys, tl)
        root = DGG.cellindex(DGG.levelgrid(sys, 2), 1)
        r = DGG.descendant_range(sys, root, tl)
        build() = begin
            h1 = DGG.subtree_halo(sys, root, tl)
            ids = vcat([DGG.cellindex(g, p) for p in r], h1)
            sort!(ids; by = c -> DGG.cellposition(g, c))
            (h1, ids)
        end
        h1, ids = build(); collect(DGG.halo(DGG.CellVector(sys, tl, ids)))   # warm
        t1 = @elapsed h1 = DGG.subtree_halo(sys, root, tl)
        _, ids = build()
        t2 = @elapsed h2 = collect(DGG.halo(DGG.CellVector(sys, tl, ids)))
        # The first and second rings together equal the chunk's two-ring halo.
        want = Set{Int}()
        for p in r, m in DGG.neighbors(g, DGG.cellindex(g, p), 2; connectivity = Vertex())
            q = DGG.cellposition(g, m)
            (first(r) <= q <= last(r)) || push!(want, q)
        end
        got = Set(DGG.cellposition(g, c) for c in vcat(h1, h2))
        @printf("  %-16s L6  k=1: %4d cells %8.2f ms | k=2 extra: %4d cells %9.2f ms (%.0fx)  correct=%s\n",
            nameof(typeof(sys)), length(h1), t1 * 1e3, length(h2), t2 * 1e3,
            t2 / t1, got == want)
    end
    println("  => halo(CellVector) is the generic OutsideWalkEngine, not the")
    println("     hexagonal/square subtree walk, so asking for a two-ring margin")
    println("     costs three orders of magnitude more on IGeo7/H3.")
end

iterator_size()
id_to_position()
wider_k()
