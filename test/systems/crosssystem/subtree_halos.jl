# ---------------------------------------------------------------------------
# The OUTSIDE face of a subtree boundary.
#
# `subtree_iterators.jl` walks the inside of that boundary — the descendants
# with a neighbour that is not one. This file walks the outside: the level-`l`
# cells that are NOT descendants but have a neighbour that is. Same boundary,
# opposite side, so the two files share a fixture vocabulary and nothing else.
#
# Written against the generic interface only: no system module is imported, and
# every law runs against every system in `systems()`.
# ---------------------------------------------------------------------------

module SubtreeHaloTests

using Test
using DiscreteGlobalGrids
using DiscreteGlobalGrids: systems, levelgrid, level, max_level, ncells,
    cellindex, neighbors, ancestor, Vertex, Edge,
    SubtreeHaloIterator, subtree_halo
import DiscreteGlobalGrids as DGG

# ---------------------------------------------------------------------------
# Depth zero: a cell's own one-ring
# ---------------------------------------------------------------------------

@testset "depth zero is the cell's own one-ring" begin
    for sys in systems()
        grid = levelgrid(sys, 1)
        c = cellindex(grid, 1)
        for conn in (Vertex(), Edge())
            expected = sort!(collect(neighbors(grid, c, 1; connectivity = conn)))
            it = SubtreeHaloIterator(sys, c, 1; connectivity = conn)
            @test collect(it) == expected
            @test subtree_halo(sys, c, 1; connectivity = conn) == expected
            @test eltype(it) == DGG.cellindextype(sys)
        end
    end
end

@testset "level validation" begin
    for sys in systems()
        grid = levelgrid(sys, 1)
        c = cellindex(grid, 1)
        @test_throws ArgumentError SubtreeHaloIterator(sys, c, 0)
        @test_throws ArgumentError SubtreeHaloIterator(sys, c, max_level(sys) + 1)
    end
end

# ---------------------------------------------------------------------------
# The defining law — the only oracle that tests the ENUMERATION
# ---------------------------------------------------------------------------

# The law itself, computed the slow honest way: every level-l cell that is not a
# descendant and has a descendant neighbour. O(ncells) and unusable in anger,
# which is exactly why it is the oracle.
#
# It is worth being precise about what this buys that the forced-geometry
# testset below does not. This law shares NO CODE with the halo walk: it decides
# which cells to consider by scanning positions 1:ncells in order, where the
# walk decides by a pruned depth-first descent of the hierarchy. So it pins
# three separate things at once — which cells the walk emits, in what order, and
# that the pruning threw nothing away. Comparing the walk against ITSELF under a
# different adjacency provider, which is what `forced_geometry_halo` does for
# today's generic engine, pins only the last question in the adjacency test and
# is blind to all three of these. That is why this testset sweeps wide and why
# widening it is not gold-plating: a bug in `_admit`'s enumeration, or in the
# cap prune's soundness margin, is invisible everywhere else in this file.
function law_halo(sys, c, l; connectivity = Vertex())
    grid = levelgrid(sys, l)
    lc = level(c)
    out = DGG.cellindextype(sys)[]
    for p in 1:ncells(grid)
        x = cellindex(grid, p)
        ancestor(sys, x, lc) == c && continue
        any(nb -> ancestor(sys, nb, lc) == c,
            neighbors(grid, x, 1; connectivity)) && push!(out, x)
    end
    return out
end

check_law(sys, c, l, conn) =
    @test collect(SubtreeHaloIterator(sys, c, l; connectivity = conn)) ==
          law_halo(sys, c, l; connectivity = conn)

# A deterministic spread: no RNG, so a failure names the same cell every run.
function sample_cells(grid, n::Int)
    total = ncells(grid)
    step = max(1, total ÷ n)
    return [cellindex(grid, i) for i in 1:step:total]
end

# The cells whose one-ring is not the modal size — pentagons, face corners,
# poles — found by DEGREE, so a sweep needs no system knowledge and a system
# added later is covered without anyone remembering to list its oddities.
function irregular_cells(grid, limit::Int)
    cells = [cellindex(grid, i) for i in 1:ncells(grid)]
    degrees = [length(neighbors(grid, c, 1)) for c in cells]
    counts = Dict{Int,Int}()
    for d in degrees
        counts[d] = get(counts, d, 0) + 1
    end
    modal = argmax(k -> counts[k], keys(counts))
    odd = [c for (c, d) in zip(cells, degrees) if d != modal]
    return odd[1:min(limit, length(odd))]
end

# The budget. `law_halo` is O(ncells) per call and the deep levels of the
# aperture-7 systems are where it bites, so the sweep SAMPLES ROOTS as it goes
# deeper rather than dropping the depth: depth is what exposes the cap prune
# (an under-covering root cap only starts dropping cells once the nodes it
# prunes are smaller than the subtree's own overhang), so a shallow-only sweep
# would be the cheap half of the coverage and the useless half.
#
# One caveat this file should not pretend away: on A5 the shipped engine is
# `ScanHaloEngine`, which enumerates by the same ascending position scan the law
# does. There the law is a check on the adjacency test and the descendant skip,
# not on the enumeration — nothing here can be independent of an engine that is
# already the naive one. A5's independent check is the forced-geometry testset.

@testset "$(nameof(typeof(sys))): the defining law" for sys in systems()
    grid0 = levelgrid(sys, 0)
    n0 = ncells(grid0)
    mx = max_level(sys)

    # EVERY level-0 root at depth 1, both connectivities. The full generation
    # rather than a sample because it is cheap and it contains every awkward
    # cell at once: all twelve IGeo7 and H3 pentagons, HEALPix's polar faces,
    # S2's six cube faces, ISEA4R's diamonds, A5's dodecahedral roots.
    for i in 1:n0, conn in (Vertex(), Edge())
        check_law(sys, cellindex(grid0, i), 1, conn)
    end

    # A spread of level-0 roots deeper. Depth 2 everywhere; depth 3 on the five
    # systems with sorted subtrees, where the walk is the pruned descent whose
    # prune this is testing (A5 has no `descendant_range`, so it scans and there
    # is no prune to break).
    if mx >= 2
        for c in sample_cells(grid0, 6), conn in (Vertex(), Edge())
            check_law(sys, c, 2, conn)
        end
    end
    if DGG.has_sorted_subtrees(sys) && mx >= 3
        for c in sample_cells(grid0, 3), conn in (Vertex(), Edge())
            check_law(sys, c, 3, conn)
        end
    end
    # Depth 4 from a level-0 root, on the two aperture-7 systems only. This is
    # the configuration where descendants overhang their parent's drawn polygon
    # by the largest margin, so it is the one that notices a root cap that has
    # stopped covering them. Not hypothetical: swapping the walk's `rootcap`
    # from `node_extent` to the under-covering `cell_cap` — which changes no
    # arithmetic, only the covering margin the prune's soundness rests on —
    # is caught HERE and nowhere else in this file, by the two H3 arms. Two
    # roots is all the runtime affords: H3's level-4 grid is 288k cells and the
    # law visits every one.
    if (sys isa DGG.IGeo7System || sys isa DGG.H3System) && mx >= 4
        for c in sample_cells(grid0, 2), conn in (Vertex(), Edge())
            check_law(sys, c, 4, conn)
        end
    end

    # Roots that are no longer whole faces or whole pentagon fans: sampled and
    # irregular-degree cells one and two levels down, each to depth 2.
    for base in 1:min(2, mx)
        gridb = levelgrid(sys, base)
        roots = unique(vcat(sample_cells(gridb, 4), irregular_cells(gridb, 2)))
        for c in roots, l in (base + 1):min(base + 2, mx), conn in (Vertex(), Edge())
            check_law(sys, c, l, conn)
        end
    end
end

# ---------------------------------------------------------------------------
# The independent oracle
# ---------------------------------------------------------------------------

# The generic walk forced onto unit-sphere boundary comparison, reached through
# the POSITIONAL constructor so the keyword one keeps choosing whatever engine
# the system ships.
forced_geometry_halo(sys, c, l, conn) = DGG.collect_subtree(
    DGG.SubtreeHaloIterator(sys, c, Int(l), conn,
        DGG.Fallbacks.geometry_halo_engine(sys, c, Int(l), conn)))

# WHAT THIS TESTSET DOES AND DOES NOT TEST, while the generic walk is the only
# engine there is.
#
# Both sides of every comparison below are the SAME `OutsideWalkEngine`: the
# same `_admit`, the same descendant-range skip, the same cap prune, the same
# depth-first descent. Only the adjacency PROVIDER differs — the system's native
# one-ring on one side, unit-sphere boundary comparison on the other. So what is
# under test here is the adjacency PREDICATE, and only that: whether the two
# definitions of "touches" agree, at pentagons, at poles, at cube corners, and
# under both connectivities. It is a real question and this is the right way to
# ask it.
#
# It is NOT an end-to-end oracle for today's engine. It cannot see a bug in the
# enumeration — a candidate the walk never considers is a candidate neither side
# considers, and a cap prune that has stopped being sound prunes both sides
# identically. `law_halo` above is what covers that, by enumerating from an
# ascending position scan that shares nothing with the walk at all.
#
# This becomes a genuine end-to-end oracle with the SPECIALIZED engines of the
# later tasks — the square band walk, the seam-aware stream merge, the
# calibrated hexagonal walks. Those enumerate differently from the generic
# descent, so comparing one against the other tests the enumeration and the
# predicate together. That is the reason this helper exists now rather than
# later: it is being stood up and pinned before it has to carry that weight.
#
# The sweep is EVERY level-0 root at depth 1, and a spread of them below. The
# full root generation rather than a sample, deliberately: it is the only cheap
# sweep that contains every structurally awkward cell at once — all twelve IGeo7
# and H3 pentagons, HEALPix's polar faces, S2's cube corners, ISEA4R's
# icosahedral-vertex diamonds and A5's twelve dodecahedral roots. Those are
# exactly the configurations where a *drawn* boundary and the *hierarchy's*
# adjacency could legitimately part company, and where A5's `Vertex()`/`Edge()`
# split is widest.
#
# NO EXCLUSION IS NEEDED. Two were anticipated — A5's connectivity split, and
# the aperture-7 systems at pentagons — and neither materialised: the two
# providers agree element for element everywhere this testset looks. If a future
# system or refinement does disagree here, the native indexed walk is
# authoritative and the exclusion belongs in this comment, named — never as a
# silent `skip`.

# Depth zero is the one case where the geometry provider has no subtree to
# descend: `root` is its own only target-level descendant. The cursor cannot
# express that — seeded at the target level it would descend past it to
# `max_level` and throw on a cell with no children — so the provider answers it
# directly against the root's own boundary. This pins that it does.
@testset "forced geometry at depth zero" begin
    for sys in systems()
        for base in 0:1
            grid = levelgrid(sys, base)
            for c in (cellindex(grid, 1), cellindex(grid, ncells(grid))),
                conn in (Vertex(), Edge())
                @test forced_geometry_halo(sys, c, base, conn) ==
                      collect(SubtreeHaloIterator(sys, c, base; connectivity = conn))
            end
        end
    end
end

@testset "$(nameof(typeof(sys))): geometry agrees with topology" for sys in systems()
    grid0 = levelgrid(sys, 0)
    n0 = ncells(grid0)
    for i in 1:n0, conn in (Vertex(), Edge())
        c = cellindex(grid0, i)
        @test forced_geometry_halo(sys, c, 1, conn) ==
              collect(SubtreeHaloIterator(sys, c, 1; connectivity = conn))
    end
    if max_level(sys) >= 2
        for i in 1:max(1, n0 ÷ 4):n0, conn in (Vertex(), Edge())
            c = cellindex(grid0, i)
            @test forced_geometry_halo(sys, c, 2, conn) ==
                  collect(SubtreeHaloIterator(sys, c, 2; connectivity = conn))
        end
        # And once more one level down, where a root is no longer a whole face
        # or a whole pentagon fan: sampled level-1 roots plus the ones whose
        # degree marks them as a seam, pole or pentagon child.
        grid1 = levelgrid(sys, 1)
        n1 = ncells(grid1)
        roots1 = unique(vcat([cellindex(grid1, i) for i in 1:max(1, n1 ÷ 4):n1],
            irregular_cells(grid1, 2)))
        for c in roots1, conn in (Vertex(), Edge())
            @test forced_geometry_halo(sys, c, 2, conn) ==
                  collect(SubtreeHaloIterator(sys, c, 2; connectivity = conn))
        end
    end
end

end # module
