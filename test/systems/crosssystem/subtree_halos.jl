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
using DiscreteGlobalGrids: systems, levelgrid, level, levels, max_level, ncells,
    cellindex, cellposition, neighbors, ancestor, descendants, descendant_range,
    subtree_border, has_sorted_subtrees, Vertex, Edge, Connectivity,
    SubtreeHaloIterator, subtree_halo, halo
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
# The defining law
# ---------------------------------------------------------------------------

# The law itself, computed the slow honest way: every level-l cell that is not a
# descendant and has a descendant neighbour. O(ncells) and unusable in anger,
# which is exactly why it is the oracle.
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

@testset "$(nameof(typeof(sys))): the defining law" for sys in systems()
    grid0 = levelgrid(sys, 0)
    for c in (cellindex(grid0, 1), cellindex(grid0, ncells(grid0)))
        for l in 1:min(2, max_level(sys)), conn in (Vertex(), Edge())
            @test collect(SubtreeHaloIterator(sys, c, l; connectivity = conn)) ==
                  law_halo(sys, c, l; connectivity = conn)
        end
    end
end

# ---------------------------------------------------------------------------
# The independent oracle
# ---------------------------------------------------------------------------

# The generic walk forced onto unit-sphere boundary comparison, reached through
# the POSITIONAL constructor so the keyword one keeps choosing whatever engine
# the system ships. This shares no index arithmetic with `neighbors`, which is
# what makes it an oracle: a specialization tested against `neighbors` or
# `subtree_border` is tested against its own topology.
forced_geometry_halo(sys, c, l, conn) = DGG.collect_subtree(
    DGG.SubtreeHaloIterator(sys, c, Int(l), conn,
        DGG.Fallbacks.geometry_halo_engine(sys, c, Int(l), conn)))

# EVERY level-0 root at depth 1, and a spread of them at depth 2.
#
# The full root generation rather than a sample, and deliberately: it is the
# only cheap sweep that contains every structurally awkward cell at once — all
# twelve IGeo7 and H3 pentagons, HEALPix's polar faces, S2's cube corners,
# ISEA4R's icosahedral-vertex diamonds and A5's twelve dodecahedral roots. Those
# are exactly the configurations where a *drawn* boundary and the *hierarchy's*
# adjacency could legitimately part company, and where A5's `Vertex()`/`Edge()`
# split is widest.
#
# NO EXCLUSION IS NEEDED. Two were anticipated — A5's connectivity split, and
# the aperture-7 systems at pentagons — and neither materialised: the two
# providers agree element for element everywhere this testset looks. If a future
# system or refinement does disagree here, the native indexed walk is
# authoritative and the exclusion belongs in this comment, named — never as a
# silent `skip`.

# The cells whose one-ring is not the modal size — pentagons, face corners,
# poles — found by DEGREE, so the sweep needs no system knowledge and a system
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
