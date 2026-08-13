# ---------------------------------------------------------------------------
# T7 — cross-system tests.
#
# Everything here is a law stated ONCE and run against EVERY system in
# `systems()`. That is the point: each system's own suite checks it against its
# own oracle, in its own vocabulary, and a law that is restated three times in
# three vocabularies is a law that can drift in one of them. These tests are
# written against the generic interface only — no system module is imported,
# nothing here knows what a Z7 digit or an nside is — so adding a system to
# `systems()` automatically subjects it to all of them.
#
# Two families:
#
#   * the ROTATIONAL NEIGHBOUR CONTRACT (the T7 owner decision): rings
#     concatenated outward, the tail-block law, and the counter-clockwise
#     winding, measured geometrically rather than restated from the lattice the
#     system used to produce it;
#   * the SUBTREE RIM HOOK: `subtree_border` / `subtree_interior` against the
#     brute-force definition. All three systems override the border with an
#     automaton, so the fallback in `src/fallbacks/subtree.jl` is never reached
#     for them — the definition is spelled out again here, independently, and
#     the automata are checked against it.
# ---------------------------------------------------------------------------

module CrossSystemTests

using Test
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: systems, levels, levelgrid, ncells, cellindex,
    cell_centroid, neighbors, ring, level, children, descendants,
    subtree_border, subtree_interior, Vertex, Edge

# A deterministic spread of cells: no RNG, so a failure names the same cell on
# every run and on every machine.
function sample_cells(grid, n::Int)
    total = ncells(grid)
    step = max(1, total ÷ n)
    return [cellindex(grid, i) for i in 1:step:total]
end

# ---------------------------------------------------------------------------
# Winding, measured
#
# The number of counter-clockwise turns a closed sequence of directions makes
# about the subject cell's centroid, seen from OUTSIDE the sphere. A single
# rotational cycle gives exactly 1.0; a sequence sorted by id gives something
# else, and a cycle wound the wrong way gives -1.0.
#
# Written from the raw coordinates rather than by reusing any system's tangent
# frame, so that it tests the order rather than restating how the order was
# produced.
# ---------------------------------------------------------------------------

function tangent_frame(p)
    # Any reference direction not parallel to `p`; the turn count does not
    # depend on which, only on the handedness.
    a = abs(p[3]) < 0.9 ? (0.0, 0.0, 1.0) : (1.0, 0.0, 0.0)
    radial = a[1] * p[1] + a[2] * p[2] + a[3] * p[3]
    t = (a[1] - radial * p[1], a[2] - radial * p[2], a[3] - radial * p[3])
    n = sqrt(t[1]^2 + t[2]^2 + t[3]^2)
    e1 = (t[1] / n, t[2] / n, t[3] / n)
    # `e2 = p x e1` makes (e1, e2) right-handed as seen from outside, so
    # increasing atan(u.e2, u.e1) is counter-clockwise from outside.
    e2 = (p[2] * e1[3] - p[3] * e1[2],
          p[3] * e1[1] - p[1] * e1[3],
          p[1] * e1[2] - p[2] * e1[1])
    return e1, e2
end

function azimuth(p, e1, e2, q)
    d = (q[1] - p[1], q[2] - p[2], q[3] - p[3])
    return atan(d[1] * e2[1] + d[2] * e2[2] + d[3] * e2[3],
                d[1] * e1[1] + d[2] * e1[2] + d[3] * e1[3])
end

"Net counter-clockwise turns the shell makes about `c`'s centroid, as a Float64."
function winding_turns(grid, c, shell)
    length(shell) < 3 && return 1.0
    p = cell_centroid(grid, c)
    e1, e2 = tangent_frame(p)
    as = [azimuth(p, e1, e2, cell_centroid(grid, x)) for x in shell]
    total = 0.0
    for i in eachindex(as)
        d = as[mod1(i + 1, length(as))] - as[i]
        # Fold each step into (-pi, pi]: a single traversal never turns more
        # than half a circle between adjacent members of a shell.
        while d <= -π
            d += 2π
        end
        while d > π
            d -= 2π
        end
        total += d
    end
    return total / (2π)
end

# The rim of `c`'s subtree at `l`, straight from the definition: a descendant
# with a neighbour that is not a descendant. Slow and obviously correct — the
# whole point is that it shares no code with the automata it checks.
function brute_force_border(sys, c, l; connectivity = Vertex())
    lc = level(c)
    l == lc && return [c]
    grid = levelgrid(sys, l)
    inside = Set(descendants(sys, c, l))
    return [d for d in descendants(sys, c, l)
            if any(nb -> !(nb in inside), neighbors(grid, d, 1; connectivity))]
end

# ---------------------------------------------------------------------------

@testset "cross-system" begin

    @testset "systems() registry" begin
        ss = systems()
        @test ss isa Tuple
        @test !isempty(ss)
        @test allunique(typeof.(ss))
        for sys in ss
            @test sys isa DGG.AbstractHierarchicalGridSystem
            @test !isempty(levels(sys))
            g = levelgrid(sys, first(levels(sys)))
            @test ncells(g) > 0
            @test DGG.system(g) == sys
        end
    end

    for sys in systems()
        name = string(nameof(typeof(sys)))

        @testset "$name: rotational neighbour contract" begin
            # A level deep enough that a k=3 disc is not most of the sphere,
            # but shallow enough to stay quick.
            l = min(3, last(levels(sys)))
            grid = levelgrid(sys, l)
            for c in sample_cells(grid, 24), conn in (Vertex(), Edge())
                shells = [collect(ring(grid, c, j; connectivity = conn)) for j in 1:3]
                for k in 1:3
                    disc = collect(neighbors(grid, c, k; connectivity = conn))
                    shell = shells[k]

                    # Rings concatenated outward, element for element.
                    @test disc == reduce(vcat, shells[1:k]; init = eltype(disc)[])

                    # ...and therefore the ring is the disc's tail block. This
                    # is the law that a separate walk for `ring` breaks while
                    # still agreeing as a set.
                    if !isempty(shell)
                        @test disc[(end - length(shell) + 1):end] == shell
                    end

                    # No duplicates, and the subject is never its own neighbour.
                    @test allunique(disc)
                    @test !(c in disc)
                end

                # Each ring is ONE counter-clockwise cycle seen from outside.
                # Sorting by id would not be, and a clockwise cycle would give
                # -1. The tolerance absorbs the shells' geometric irregularity
                # near pentagons and face seams.
                for shell in shells
                    length(shell) >= 3 || continue
                    @test winding_turns(grid, c, shell) ≈ 1.0 atol = 0.05
                end
            end
        end

        @testset "$name: subtree rim hook" begin
            root_level = first(levels(sys))
            grid0 = levelgrid(sys, root_level)
            probes = sample_cells(grid0, 6)
            for c in probes
                lc = level(c)

                # A depth-0 subtree is the cell itself, and it is all rim.
                @test collect(subtree_border(sys, c, lc)) == [c]
                @test isempty(collect(subtree_interior(sys, c, lc)))

                for depth in 1:2
                    l = lc + depth
                    l <= last(levels(sys)) || continue

                    border = collect(subtree_border(sys, c, l))
                    interior = collect(subtree_interior(sys, c, l))
                    kids = collect(descendants(sys, c, l))

                    # The automaton against the definition. Every system here
                    # overrides `subtree_border`, so this is the differential
                    # test that keeps the fast path honest.
                    @test Set(border) == Set(brute_force_border(sys, c, l))

                    # Border and interior partition the subtree.
                    @test isempty(intersect(Set(border), Set(interior)))
                    @test union(Set(border), Set(interior)) == Set(kids)
                    @test length(border) + length(interior) == length(kids)

                    # The rim is a small minority once there is any depth to
                    # speak of — the property that makes the hook worth having.
                    depth >= 2 && @test length(border) < length(kids)

                    @test allunique(border)
                    @test eltype(border) === DGG.cellindextype(sys)

                    # The interface documents the border's order as ascending
                    # canonical order unless a system says otherwise, and none
                    # of the three does. Verified across all three automatons
                    # before pinning it here: the Z7 digit automaton, H3's
                    # digit-arc automaton and HEALPix's Morton rim walk all
                    # emit ascending by construction (the Morton walk visits
                    # quadrants in id order precisely so that it can). Left
                    # unpinned, an automaton could start emitting a rim in walk
                    # order and only the docs would be wrong.
                    @test issorted(border)
                    @test issorted(interior)
                end

                # Asking for a level above the cell's own is an error, not an
                # empty answer, and uniformly so across systems.
                lc > first(levels(sys)) &&
                    @test_throws ArgumentError subtree_border(sys, c, lc - 1)
            end
        end
    end
end

end # module CrossSystemTests
