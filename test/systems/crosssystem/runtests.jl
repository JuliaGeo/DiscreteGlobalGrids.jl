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
#     brute-force definition, spelled out here independently of both.
#
#     Five registered systems override the border with an `O(rim)` automaton —
#     IGeo7's Z7 digit predicate, H3's digit-arc walk, and the shared
#     aperture-4 square walk under HEALPix's and ISEA4R's Morton curve and S2's
#     Hilbert one — and for those this is the differential test that keeps the
#     fast path honest.
#
#     The remaining systems deliberately retain the fallback. A5 has no
#     `descendant_range`; the new systems have not yet added a family-specific
#     rim walker. For them the same assertions still bite, just one rung lower:
#     they check
#     `src/fallbacks/subtree.jl` against the definition. That is not circular —
#     the fallback decides membership by walking each neighbour up to the root
#     with `ancestor`, while `brute_force_border` below materialises the
#     descendant set and asks it — so the two agree only if both are right.
#
#     The lazy iterators those walks are made of get their own file,
#     `subtree_iterators.jl`.
# ---------------------------------------------------------------------------

module CrossSystemTests

using Test
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: systems, levels, levelgrid, ncells, cellindex,
    cell_boundary, cell_centroid, cellat, neighbors, ring, level, children,
    descendants, subtree_border, subtree_interior, Vertex, Edge, PartialGrid

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

    @testset "who ships a subtree-rim automaton" begin
        # The header comment names four systems with an `O(rim)` automaton and
        # two that keep the generic fallback. Pinned here so the claim cannot
        # rot the way its predecessor did ("all three systems override the
        # border" survived two systems being added that do not): a system that
        # gains or loses an automaton must come and edit this list, which is
        # exactly the moment to re-read the comment.
        # T20 moved S2 across: its subtree is an aligned square block on one
        # face and its rim is that block's perimeter, so it now shares the
        # aperture-4 walk with HEALPix and ISEA4R under its own Hilbert curve.
        # A5 is the only system left on the fallback, and structurally so — no
        # `descendant_range`, and no closed-form child adjacency to build an
        # automaton from.
        automaton = Set([:IGeo7System, :H3System, :HEALPixSystem, :ISEA4RSystem,
            :S2System])
        fallback = Set([:A5System, :ISEA3HSystem, :ISEA4HSystem, :ISEA4TSystem,
            :RHEALPixSystem, :AuthalicSystem, :IVEA4RSystem, :IVEA9RSystem,
            :RTEA4RSystem, :RTEA9RSystem])
        for sys in systems()
            n = nameof(typeof(sys))
            c = cellindex(levelgrid(sys, first(levels(sys))), 1)
            m = which(subtree_border, Base.typesof(sys, c, level(c)))
            overrides = m.module !== DGG.Fallbacks
            # Module provenance, not specificity, and deliberately: the
            # question here is "did this system write a rim walker", and
            # `Fallbacks` owning the method is precisely the answer "no".
            @test (n in automaton) ⊻ (n in fallback)   # nobody unaccounted for
            @test overrides == (n in automaton)
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

        # -------------------------------------------------------------------
        # The generic point-in-cell test, against the one point that is
        # unarguably inside a cell.
        #
        # This is a REGRESSION LAW for a substrate bug, not a system contract,
        # which is why it lives here and not in the conformance package: what
        # it guards is `Fallbacks.point_in_cell`, and the way to guard a
        # generic is to sweep it over every system in the registry, which is
        # precisely what this file is for. Putting it in the conformance
        # package would have restated a defect in the shared substrate as an
        # obligation each system owes separately.
        #
        # The bug (T2 vintage, found by T9, fixed in T13): `point_in_cell`
        # asked `spherical_ring_encloses` first, whose test arc runs from the
        # query point to the antipode of the ring's vertex mass. For a point
        # INSIDE a small cell that arc is a near-half-turn, and its
        # between-ness test then admits every point on the sphere — so the
        # cell's own centroid was reported OUTSIDE its own boundary for 282 of
        # 3072 HEALPix level-4 cells, and for cells of H3, S2 and ISEA4R
        # besides, while IGeo7's twelve exactly-symmetric pentagons returned
        # `nothing`. See `anchor_arc_is_conditioned`.
        #
        # Both halves matter. The predicate is the bug; `cellat` on a
        # `PartialGrid` is the symptom a user would actually hit, because a
        # partial grid has no native point location to hide the fallback.
        # -------------------------------------------------------------------
        @testset "$name: fallback point-in-cell accepts a cell's own centroid" begin
            l = min(2, last(levels(sys)))
            grid = levelgrid(sys, l)
            probes = sample_cells(grid, 200)

            offenders = [c for c in probes
                         if DGG.Fallbacks.point_in_cell(cell_boundary(grid, c),
                                                        cell_centroid(grid, c)) !== true]
            # Named, not merely counted: a bare count makes a regression here
            # a number to argue with rather than a cell to go and look at.
            @test isempty(offenders)
            isempty(offenders) || @info "$name: centroid rejected" first(offenders, 5)

            # The other side of the same predicate: the ANTIPODE of a centroid
            # is not inside the cell.
            #
            # `point_in_cell`'s last resort is a winding number, which measures
            # how often the ring encircles the probe — and a ring encircles the
            # antipode of a point exactly as often as it encircles the point
            # itself, with the opposite sign. A bare magnitude test therefore
            # calls both contained. That is why `ring_winding_verdict` carries a
            # sub-hemispheric guard, and why it is pinned here directly as well
            # as through the cascade: the mode is LATENT, because for these
            # probes one of the two parity algorithms answers first and the
            # winding verdict is never reached, so a test that only went through
            # `point_in_cell` would pass with the guard removed and leave the
            # safety resting on an unstated cascade-ordering invariant.
            wrong_side = [c for c in probes
                          if DGG.Fallbacks.point_in_cell(cell_boundary(grid, c),
                                                         -cell_centroid(grid, c)) !== false]
            @test isempty(wrong_side)

            @test all(probes) do c
                r, m = DGG.Fallbacks.open_ring(cell_boundary(grid, c))
                # inside: decided; antipode of inside: abstains outright
                DGG.Fallbacks.ring_winding_verdict(r, m, cell_centroid(grid, c)) === true &&
                DGG.Fallbacks.ring_winding_verdict(r, m, -cell_centroid(grid, c)) === nothing
            end

            # ...and the consequence. A `PartialGrid` has no native `cellat`,
            # so this is the generic descend-and-test path end to end: the
            # centroid of a cell in the subset must locate that same cell.
            subset = sort(probes)
            pg = PartialGrid(sys, l, subset)
            missed = [c for c in subset if cellat(pg, cell_centroid(grid, c)) !== c]
            @test isempty(missed)
            isempty(missed) || @info "$name: cellat missed its own cell" first(missed, 5)
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

                    # Against the definition — the automaton's differential
                    # test for the four systems that ship one, and the generic
                    # fallback's for A5 and S2. See the header for which is
                    # which and why.
                    @test Set(border) == Set(brute_force_border(sys, c, l))

                    # Border and interior partition the subtree.
                    @test isempty(intersect(Set(border), Set(interior)))
                    @test union(Set(border), Set(interior)) == Set(kids)
                    @test length(border) + length(interior) == length(kids)

                    # The rim is a small minority once there is any depth and
                    # more than one descendant. ISEA3H/4H's two polar roots
                    # are intentionally singleton all-zero chains, so their
                    # subtrees have no possible interior at any depth.
                    depth >= 2 && length(kids) > 1 &&
                        @test length(border) < length(kids)

                    @test allunique(border)
                    @test eltype(border) === DGG.cellindextype(sys)

                    # The interface documents the border's order as ascending
                    # canonical order unless a system says otherwise, and none
                    # of the registered systems does. Verified across the fast paths
                    # before pinning it here: the Z7 digit automaton, H3's
                    # digit-arc automaton, HEALPix's Morton rim walk and
                    # ISEA4R's edge walk all emit ascending by construction
                    # (the Morton walk visits quadrants in id order precisely
                    # so that it can). A5 and S2 inherit it from the fallback,
                    # which preserves `descendants` order. Left unpinned, an
                    # automaton could start emitting a rim in walk order and
                    # only the docs would be wrong.
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
