#!/usr/bin/env julia
# Report iterator sizing, ID-to-position conversion, and wider-halo costs.
#
#   julia --project=test examples/geomorphometry/api_ergonomics.jl

import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: Vertex, Edge
using Printf

function iterator_size()
    println("### 1. IteratorSize / length / sizehint: HasLength only at depth 0")
    for sys in DGG.systems()
        r = DGG.cellindex(DGG.levelgrid(sys, 2), 1)
        line = String[]
        for l in (2, 3, 6)
            it = DGG.halo(DGG.subtree(sys, r, l); cells = true)
            len = try string(length(it)) catch e; string(nameof(typeof(e))) end
            push!(line, "L$l:$(Base.IteratorSize(typeof(it)))/$len/$(DGG.sizehint(it))")
        end
        @printf("  %-16s %s\n", nameof(typeof(sys)), join(line, "  "))
    end
    println("  => `length` is exact only at depth 0, so a read buffer is reserved")
    println("     against `sizehint`, which may over- or undershoot.")
end

function id_to_position()
    println()
    println("### 2. the halo answers in positions; ids cost a conversion")
    for (sys, rl, tl) in ((DGG.IGeo7System(), 4, 10), (DGG.H3System(), 4, 10),
                          (DGG.HEALPixSystem(), 4, 12), (DGG.S2System(), 4, 12),
                          (DGG.ISEA4RSystem(), 4, 12))
        g = DGG.levelgrid(sys, tl)
        root = DGG.cellindex(DGG.levelgrid(sys, rl), 3)
        pg = DGG.subtree(sys, root, tl)
        collect(DGG.halo(pg; cells = true))
        tw = @elapsed h = collect(DGG.halo(pg; cells = true))
        DGG.cellposition.(Ref(g), h)
        tp = @elapsed DGG.cellposition.(Ref(g), h)
        ap = @allocated DGG.cellposition.(Ref(g), h)
        collect(DGG.halo(pg))
        tn = @elapsed collect(DGG.halo(pg))
        @printf("  %-16s halo=%4d  ids %8.1f us   id->pos %6.1f us (%4.1f%%, %d B)   positions %8.1f us\n",
            nameof(typeof(sys)), length(h), tw * 1e6, tp * 1e6,
            100 * tp / (tw + tp), ap, tn * 1e6)
    end
    println("  => the conversion is cheap, and the position walk folds it in, so a")
    println("     caller indexing an array writes no broadcast and no sorted map.")
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
        pg = DGG.subtree(sys, root, tl)
        build() = begin
            ids = vcat([DGG.cellindex(g, p) for p in r],
                collect(DGG.halo(pg; cells = true)))
            sort!(ids; by = c -> DGG.cellposition(g, c))
            ids
        end
        ids = build(); collect(DGG.halo(DGG.CellVector(sys, tl, ids)))   # warm
        t1 = @elapsed h1 = collect(DGG.halo(pg))
        ids = build()
        t2 = @elapsed h2 = collect(DGG.halo(DGG.CellVector(sys, tl, ids)))
        # The first and second rings together equal the chunk's two-ring halo.
        want = Set{Int}()
        for p in r, m in DGG.neighbors(g, DGG.cellindex(g, p), 2; connectivity = Vertex())
            q = DGG.cellposition(g, m)
            (first(r) <= q <= last(r)) || push!(want, q)
        end
        got = Set(vcat(h1, h2))
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
