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
    cellindex, cellposition, neighbors, ancestor, subtree_border, Vertex, Edge,
    SubtreeHaloIterator, subtree_halo
import DiscreteGlobalGrids as DGG

# ---------------------------------------------------------------------------
# The one fixture that cannot live inside the outer testset
#
# Everything below runs inside ONE outer `@testset`, so a failure anywhere is
# recorded and every later section still runs: a TOP-LEVEL testset throws when
# it finishes with failures, which would abort the rest of the file, while a
# nested one only reports upwards. A `struct` cannot be declared in that local
# scope, so the lying engine the guard testset needs is declared out here —
# next to the reason it is not declared where it is used.
# ---------------------------------------------------------------------------

# Claims three, yields one — exactly the shape `collect_subtree` exists to
# catch. Without it `collect`'s own `HasLength` route sizes the vector from the
# claim and hands back two `undef` slots as cell ids; the square band walk is an
# engine that claims a closed-form count, so this guard is load-bearing for it
# and not a curiosity.
struct MiscountingEngine end
Base.iterate(::MiscountingEngine) = (DGG.LevelIndex(0, 0), 1)
Base.iterate(::MiscountingEngine, ::Int) = nothing
Base.eltype(::Type{MiscountingEngine}) = DGG.LevelIndex
Base.IteratorSize(::Type{MiscountingEngine}) = Base.HasLength()
Base.length(::MiscountingEngine) = 3

# Everything from here down is ONE testset, so a failure in an early section is
# recorded and the rest of the file still runs. See the fixture note above.
@testset "subtree halos" begin

    # -----------------------------------------------------------------------
    # Depth zero: a cell's own one-ring
    # -----------------------------------------------------------------------

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

    # -----------------------------------------------------------------------
    # The defining law — the only oracle that tests the ENUMERATION
    # -----------------------------------------------------------------------

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
        # is caught here, by the two H3 arms — and by the H3 arms only; the IGeo7
        # half of this same loop passes, which is why both systems are swept
        # rather than one. Two roots is all the runtime affords: H3's level-4 grid
        # is 288k cells and the law visits every one.
        #
        # SINCE TASK 6 THIS IS NO LONGER THE ONLY ARM THAT CATCHES IT, and the
        # strengthening is worth recording because it was a side effect rather
        # than a design. `forced_geometry_halo` IS the generic walk — the same
        # `_admit`, the same `rootcap` — so once H3 and IGeo7 grew a directed walk
        # that shares none of that, every comparison of the two became a comparison
        # of a sound enumeration against a broken one. The same mutation now also
        # fails "the directed walk against forced geometry" (118 + 88 assertions)
        # and "the directed walk at depth four" (2 + 2). What is still true is the
        # narrow claim: the arm below is where an under-covering cap is caught with
        # the generic walk on BOTH sides of the comparison, which is the only shape
        # that would survive the specializations being removed.
        #
        # THE GENERIC WALK IS BUILT EXPLICITLY, not reached. Every system now
        # ships a specialization, so the keyword constructor no longer returns
        # `OutsideWalkEngine` anywhere and this arm would otherwise have stopped
        # testing the cap prune the moment Task 6 landed — silently, since both
        # engines answer correctly. The law is computed once and both walks are
        # held to it.
        if (sys isa DGG.IGeo7System || sys isa DGG.H3System) && mx >= 4
            for c in sample_cells(grid0, 2), conn in (Vertex(), Edge())
                want = law_halo(sys, c, 4; connectivity = conn)
                @test collect(SubtreeHaloIterator(sys, c, 4; connectivity = conn)) ==
                      want
                @test collect(SubtreeHaloIterator(sys, c, 4, conn,
                    DGG.Fallbacks.generic_halo_engine(sys, c, 4, conn))) == want
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

    # -----------------------------------------------------------------------
    # The independent oracle
    # -----------------------------------------------------------------------

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

    # -----------------------------------------------------------------------
    # The sweep harness
    #
    # `law_halo` above answers "is this the right SET, in the right order?" by a
    # scan that shares nothing with the walk, and it is the strongest oracle here —
    # but it is `O(ncells)` per call, so it can only be afforded shallow. The bundle
    # below is the other half: a fixed list of the iterator contract's laws, cheap
    # enough to run at every root, level and connectivity the sweep reaches, and the
    # thing every specialization added later is put through unchanged. It is not an
    # oracle — it cannot tell a walk that drops a cell from one that never should
    # have emitted it — but it pins uniqueness, sortedness, outside ancestry, both
    # adjacency directions, `subtree_halo`/`collect` agreement, and any `length` an
    # engine claims.
    # -----------------------------------------------------------------------

    # Bases 0, 1 and 2, so a root is a whole face, then a quarter of one, then a
    # sixteenth — the last is the first that can be nowhere flush with its face
    # edge, which is the configuration the square band walk needs.
    #
    # A5 STOPS AT BASE 1, on purpose and not to hide anything: it is the one system
    # with no `descendant_range`, so it takes `ScanHaloEngine`, which is `O(ncells)`
    # per halo — and its `neighbors` allocates a `Set` per call, so the constant is
    # large too. Sweeping it at base 2 (targets up to level 4, 3840 cells scanned
    # twice per case) costs more than the other five systems together and pins
    # nothing they do not. Bases 0 and 1 still run every law at every level, and
    # `law_halo` and the forced-geometry testset above both cover A5 at full width.
    sweep_bases(sys) = filter(b -> b <= max_level(sys),
        sys isa DGG.A5System ? (0, 1) : (0, 1, 2))

    sweep_roots(sys, base::Int) = (grid = levelgrid(sys, base);
        unique(vcat(sample_cells(grid, 4), irregular_cells(grid, 2))))

    # The deepest target whose grid is within `budget` times the root generation's,
    # so "deep" means the same amount of WORK on an aperture-4 and an aperture-7
    # system rather than the same number of levels.
    function deep_depth(sys, base::Int, budget::Int = 70_000)
        d = 0
        while base + d + 1 <= max_level(sys) &&
            ncells(levelgrid(sys, base + d + 1)) ÷ ncells(levelgrid(sys, base)) <= budget
            d += 1
        end
        return d
    end

    # The per-case law bundle. Every engine, every specialization, goes through it.
    function check_halo_case(sys, c, l, conn)
        it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
        h = collect(it)
        @test h == subtree_halo(sys, c, l; connectivity = conn)
        @test allunique(h)
        lc = level(c)
        @test all(x -> ancestor(sys, x, lc) != c, h)          # outside ancestry
        grid = levelgrid(sys, l)
        @test issorted([cellposition(grid, x) for x in h])    # canonical order
        # Both adjacency directions, under the same connectivity: every border cell
        # reaches the halo, and every halo cell reaches the border.
        border = subtree_border(sys, c, l; connectivity = conn)
        hs, bs = Set(h), Set(border)
        @test all(r -> any(in(hs), neighbors(grid, r, 1; connectivity = conn)), border)
        @test all(x -> any(in(bs), neighbors(grid, x, 1; connectivity = conn)), h)
        Base.IteratorSize(typeof(it)) isa Base.HasLength && @test length(it) == length(h)
        return h
    end

    @testset "$(nameof(typeof(sys))) at level $base" for sys in systems(),
            base in sweep_bases(sys)
        for c in sweep_roots(sys, base), l in base:min(base + 2, max_level(sys))
            hv = check_halo_case(sys, c, l, Vertex())
            he = check_halo_case(sys, c, l, Edge())
            # `Edge()` is `Vertex()` minus the cells that touch at a point only.
            @test issubset(Set(he), Set(hv))
        end
    end

    # One deep case per system, where `deep_depth` puts the target grid within a
    # fixed factor of the root generation. One root and one connectivity: the point
    # is that the laws still hold when the halo is hundreds of cells and the descent
    # is long, not to re-sweep at depth. A5 is excluded for the reason `sweep_bases`
    # gives — its `deep_depth` from a level-0 root is level 7, a million-cell scan.
    @testset "$(nameof(typeof(sys))) at depth" for sys in
            filter(s -> !(s isa DGG.A5System), systems())
        d = deep_depth(sys, 0)
        d >= 1 || continue
        check_halo_case(sys, cellindex(levelgrid(sys, 0), 1), d, Vertex())
    end

    # -----------------------------------------------------------------------
    # The square band walk — HEALPix, S2 and ISEA4R away from a face edge
    # -----------------------------------------------------------------------

    # WHY THE ORACLE HERE IS `forced_geometry_halo` AND NOT `law_halo` ALONE.
    #
    # For the generic engine the forced-geometry comparison was only a test of the
    # adjacency predicate: both sides ran the same `OutsideWalkEngine`, the same
    # `_admit`, the same descendant-range skip, the same depth-first descent. That
    # stops being true here. `SquareBandEngine` enumerates by descending the FACE's
    # quadtree in curve order and pruning by lattice overlap — it never calls
    # `_admit`, never calls `descendant_range`, never looks at a cap, and never asks
    # the system for a neighbour at all. So comparing it against the forced-geometry
    # walk now pins the enumeration and the predicate together, end to end, which is
    # exactly what the design asks specializations to be tested by. `law_halo` is
    # kept as one arm below because it is cheaper to be sure than to argue: it
    # enumerates from an ascending position scan and shares nothing with either.
    #
    # What must NOT be the only oracle is `neighbors` or `subtree_border` — the band
    # walk is index arithmetic on the same lattice those are built from, so they
    # would agree with a wrong band for the same reason it was wrong.

    # Which walk a root gets is a fact about the lattice, not something a test
    # should hard-code, so roots are CLASSIFIED by the engine the constructor
    # actually chose. Three classes, because there are three walks:
    #
    #   `inface`   — `SquareBandEngine` under `NoCheck`: the block is nowhere
    #                flush with its face edge, the band IS the halo, and the
    #                count is closed form.
    #   `seam`     — `SquareBandEngine` under `NativeCheck`: the block is flush
    #                somewhere, the rectangles are a conservative superset, and
    #                every candidate goes through the native one-ring.
    #   `fallback` — anything else, i.e. the generic outside-first walk.
    #
    # The two CLAIMED classes are checked against the oracle below, cell for cell.
    # `fallback` is only asserted EMPTY (`check_root_classes`) — on these three
    # systems it is supposed to have no members at all, so there is nothing here
    # to oracle. The generic walk it names is still oracled, on the systems that
    # genuinely take it: see "the generic fallback still agrees with the oracle".
    function classify_roots(sys, base::Int, l::Int, conn)
        grid = levelgrid(sys, base)
        C = DGG.cellindextype(sys)
        inface, seam, fallback = C[], C[], C[]
        for i in 1:ncells(grid)
            c = cellindex(grid, i)
            e = SubtreeHaloIterator(sys, c, l; connectivity = conn).engine
            if e isa DGG.Fallbacks.SquareBandEngine
                push!(e.check isa DGG.Fallbacks.NoCheck ? inface : seam, c)
            else
                push!(fallback, c)
            end
        end
        return inface, seam, fallback
    end

    # Evenly spaced picks, so a sample of a face-ordered list crosses faces
    # instead of staying on the first one.
    function spread(v, n::Int)
        isempty(v) && return v
        length(v) <= n && return v
        step = length(v) ÷ n
        return [v[1 + (i - 1) * step] for i in 1:n]
    end

    # HOW MANY ROOTS EACH WALK MUST CLAIM, in closed form — and why that has to be
    # pinned separately from everything else in this section.
    #
    # `classify_roots` reads the classes OFF THE CODE UNDER TEST: whatever the
    # three `halo_engine` methods do is what it reports. So a guard that grows
    # STRICTER is invisible to every other assertion here — the blocks it stops
    # claiming fall to the seam walk, which probes nothing on a non-flush block
    # and so reduces to the same band box filtered by the native one-ring. Right
    # answer, slower walk, and every oracle comparison still passes element for
    # element. Only a count notices. Narrowing HEALPix's interval test to
    # `side <= x0 && x0 + 2side <= n - 1`, say, misroutes 16 blocks per base and
    # fails exactly the two counts below.
    #
    # A depth-`d` block at base `b` sits at lattice origin `(ix, iy) · 2^d` with
    # `ix, iy ∈ [0, 2^b)`, and the width-one band fits inside `[0, n-1]²` exactly
    # when `1 <= ix` and `ix <= 2^b - 2`, likewise `iy` — independent of `d`. So
    # each face contributes `max(0, 2^b - 2)²` in-face blocks and the count is that
    # times the number of faces, which is the level-0 generation: 12 on HEALPix, 6
    # on S2, 10 on ISEA4R. At bases 0 and 1 that is ZERO — a level-0 block is the
    # whole face and a level-1 block is flush on one side per axis — so those two
    # bases are entirely the seam walk, which is the point of them being swept.
    #
    # Every remaining root is the seam walk and NOTHING falls back: that is the
    # third assertion, and it is the one that would fail if a seam configuration
    # were quietly handed to `generic_halo_engine`.
    #
    # (One mutation this CANNOT catch, because it is not one: tightening the guard
    # to `2 <= x0 && x0 + side <= n - 2` admits exactly the same blocks. `x0` is
    # `ix << d` with `d >= 1`, so `x0` and `x0 + side` are even and `n` is a power
    # of two — `x0 >= 1` and `x0 >= 2` are the same predicate here, as are
    # `<= n - 1` and `<= n - 2`. The guard has a spare parity of slack in it.)
    inface_root_count(sys, base::Int) =
        ncells(levelgrid(sys, 0)) * max(0, (1 << base) - 2)^2

    function check_root_classes(sys, base::Int, inface, seam, fallback)
        total = ncells(levelgrid(sys, base))
        @test length(inface) == inface_root_count(sys, base)
        @test length(seam) == total - inface_root_count(sys, base)
        @test isempty(fallback)
    end

    SQUARE_SYSTEMS = (HEALPixSystem(), S2System(), ISEA4RSystem())

    # BASE 2 IS NOT ENOUGH, and this is measured, not defensive. At base 2 exactly
    # four of the sixteen cells per face are non-flush, and they are the four at
    # lattice `(1,1)`, `(1,2)`, `(2,1)`, `(2,2)` — every one of which is mapped to
    # itself by the square symmetry that a wrong S2 orientation seed induces
    # (`SWAP` fixes the two diagonal blocks, `SWAP|INVERT` the two anti-diagonal
    # ones). Seeding the descent with `_hilbert_orientation(c.index, ...)` instead
    # of the face root's `isodd(face) ? SWAP_MASK : 0x0` — the exact error PR #19
    # recorded — therefore passes every base-2 arm of this file element for
    # element, on all 144 band cases, while being wrong.
    #
    # At base 3 the non-flush set is 36 of 64 per face and most of those blocks are
    # not fixed by either symmetry: the same mutation fails 576 of 864 cases. So
    # base 3 is in the sweep, and a comment is not a substitute for it.
    BAND_BASES = (2, 3)

    @testset "$(nameof(typeof(sys))): the band walk against forced geometry" for sys in
            SQUARE_SYSTEMS
        for base in BAND_BASES, d in 1:2, conn in (Vertex(), Edge())
            l = base + d
            l <= max_level(sys) || continue
            inface, seam, fallback = classify_roots(sys, base, l, conn)
            # The specialization was reached, and reached on exactly the blocks the
            # lattice says it should be. See `inface_root_count`.
            check_root_classes(sys, base, inface, seam, fallback)
            # Six claimed roots: enough to cross faces (four per face at base 2, so
            # this reaches at least two of them) without paying for the oracle 64
            # times.
            for c in inface[1:min(6, length(inface))]
                it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
                @test it.engine isa DGG.Fallbacks.SquareBandEngine
                @test it.engine.check isa DGG.Fallbacks.NoCheck
                @test collect(it) == forced_geometry_halo(sys, c, l, conn)
                check_halo_case(sys, c, l, conn)
            end
        end
    end

    # One arm against the O(ncells) brute force as well. `law_halo` shares nothing
    # with either the band walk or the geometry walk — it decides which cells to
    # consider by scanning positions 1:ncells — so it is the one check that a wrong
    # band and a wrong oracle cannot pass together. Depth 2 only: `law_halo` visits
    # every cell of the target level, and levels 4 and 5 are where that is still
    # cheap on all three systems (at most 12288 cells).
    @testset "$(nameof(typeof(sys))): the band walk against the law" for sys in
            SQUARE_SYSTEMS
        for base in BAND_BASES, conn in (Vertex(), Edge())
            l = base + 2
            l <= max_level(sys) || continue
            inface, _, _ = classify_roots(sys, base, l, conn)
            for c in inface[1:min(4, length(inface))]
                @test collect(SubtreeHaloIterator(sys, c, l; connectivity = conn)) ==
                      law_halo(sys, c, l; connectivity = conn)
            end
        end
    end

    # The closed-form count, VERIFIED rather than declared: `4·side + 4` band cells
    # under `Vertex()` and `4·side` under `Edge()`, on every block size the sweep
    # can reach. `IteratorSize() == HasLength()` rests on this, and `collect_subtree`
    # turns a miscount into an `error` rather than an `undef` tail — but an engine
    # that counted wrong AND walked wrong by the same amount would slip past that,
    # so the formula is pinned against a real collect here.
    @testset "$(nameof(typeof(sys))): the band count is closed form" for sys in
            SQUARE_SYSTEMS
        for base in BAND_BASES, d in 1:4
            l = base + d
            l <= max_level(sys) || continue
            side = 1 << d
            for conn in (Vertex(), Edge())
                inface, _, _ = classify_roots(sys, base, l, conn)
                isempty(inface) && continue
                for c in inface[1:min(3, length(inface))]
                    it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
                    @test Base.IteratorSize(typeof(it)) isa Base.HasLength
                    @test length(it) == (conn isa Vertex ? 4side + 4 : 4side)
                    @test length(collect(it)) == length(it)
                end
            end
        end
    end

    # The depth check PR #19's review found necessary: a 64x64 block, whose band is
    # 260 cells, element for element against the oracle on both connectivities. The
    # shallow cases above all have `side <= 8`, where a wrong `_restore_code` on the
    # way back up the face descent can still land on the right cell by accident; at
    # nine levels of descent it cannot. Three roots from base 3 rather than one from
    # base 2, for `BAND_BASES`' reason: a base-2 block is symmetric under the very
    # transforms a wrong descent applies.
    @testset "the band walk at 64x64" begin
        sys = S2System()
        l = 3 + 6
        inface, _, _ = classify_roots(sys, 3, l, Vertex())
        @test !isempty(inface)
        for c in inface[1:min(3, length(inface))], conn in (Vertex(), Edge())
            it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
            @test it.engine isa DGG.Fallbacks.SquareBandEngine
            @test length(it) == (conn isa Vertex ? 260 : 256)
            @test collect(it) == forced_geometry_halo(sys, c, l, conn)
            # The contract bundle at nine levels of descent too: sortedness and both
            # adjacency directions are where a deep walk can go wrong differently
            # from a shallow one, and the oracle comparison alone does not say so.
            check_halo_case(sys, c, l, conn)
        end
    end

    # The corner law. `Edge()` drops exactly the four cells that touch the block at
    # a vertex only, so the two halos differ by four and nothing else. Pinned to a
    # NON-FLUSH block deliberately: a flush block can lose a corner across a cube
    # corner, where three faces meet and the diagonal neighbour does not exist, so
    # the count would be three there and the law would read as broken.
    @testset "Edge drops exactly the four diagonal corners" begin
        sys = S2System()
        inface, _, _ = classify_roots(sys, 2, 4, Vertex())
        c = first(inface)
        hv = collect(SubtreeHaloIterator(sys, c, 4; connectivity = Vertex()))
        he = collect(SubtreeHaloIterator(sys, c, 4; connectivity = Edge()))
        @test length(setdiff(Set(hv), Set(he))) == 4
        @test issubset(Set(he), Set(hv))
    end

    # -----------------------------------------------------------------------
    # The seam walk — the same engine where the block touches a face edge
    # -----------------------------------------------------------------------

    # WHAT IS DIFFERENT HERE, AND WHY THE ORACLE MATTERS MORE.
    #
    # The in-face band is exact by construction: the band and the halo are the
    # same set, so a comparison against geometry is checking arithmetic. The seam
    # band is a conservative SUPERSET filtered by the native one-ring, so there
    # are two ways to be wrong and only one of them is loud. Yielding a cell that
    # is not a halo cell fails the filter and would fail everything below.
    # MISSING one is silent: the rectangles simply never propose it, the filter
    # never sees it, and `neighbors` and `subtree_border` would agree with the
    # gap, because the missing cell is missing from the same index arithmetic
    # they are built from. `forced_geometry_halo` is the only oracle in this file
    # that can see a candidate the walk never considered, which is why every arm
    # of this section goes through it, and why one arm goes through `law_halo` as
    # well.
    #
    # BASES 0 AND 1 ARE THE POINT. A level-0 block is the whole face and is flush
    # on all four sides; a level-1 block is flush on one side per axis and its
    # corner cell is a face corner. Those are precisely the configurations the
    # in-face guard used to send to the generic walk, so they are the surface
    # this section exists to retire — and `check_root_classes` asserts that all
    # of them are claimed, so a quiet re-routing back to the generic walk fails
    # here rather than passing everywhere.
    #
    # DEPTH 3 rather than 2 because the derivation's monotonicity argument is
    # about the INTERIOR rim cells of a flush side, and only two cells of each
    # side are ever probed. At depth 1 a side is two cells and both are probed,
    # so the argument is not exercised at all; at depth 2 it is two of four; at
    # depth 3 it is two of eight, and six cells of every flush side reach faces
    # no probe ever asked about.
    seam_roots(sys, base::Int, seam) = unique(vcat(spread(seam, 6),
        filter(in(Set(seam)), irregular_cells(levelgrid(sys, base), 4))))

    @testset "$(nameof(typeof(sys))): the seam walk against forced geometry" for sys in
            SQUARE_SYSTEMS
        for base in (0, 1, 2, 3), d in 1:3, conn in (Vertex(), Edge())
            l = base + d
            l <= max_level(sys) || continue
            inface, seam, fallback = classify_roots(sys, base, l, conn)
            check_root_classes(sys, base, inface, seam, fallback)
            @test !isempty(seam)                 # the seam path was reached
            for c in seam_roots(sys, base, seam)
                it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
                @test it.engine isa DGG.Fallbacks.SquareBandEngine
                @test it.engine.check isa DGG.Fallbacks.NativeCheck
                @test collect(it) == forced_geometry_halo(sys, c, l, conn)
                check_halo_case(sys, c, l, conn)
            end
        end
    end

    # And one arm against the O(ncells) brute force, which shares nothing with
    # either the seam walk or the geometry walk. Bases 0 and 1 at depth 2, where
    # the target grids are at most 3072 cells on all three systems.
    @testset "$(nameof(typeof(sys))): the seam walk against the law" for sys in
            SQUARE_SYSTEMS
        for base in (0, 1), conn in (Vertex(), Edge())
            l = base + 2
            l <= max_level(sys) || continue
            _, seam, _ = classify_roots(sys, base, l, conn)
            for c in spread(seam, 4)
                @test collect(SubtreeHaloIterator(sys, c, l; connectivity = conn)) ==
                      law_halo(sys, c, l; connectivity = conn)
            end
        end
    end

    # The count contract, in the negative. No perimeter formula survives a seam —
    # a cube corner is three cells where the in-face rule wants four, an ISEA4R
    # icosahedral vertex is five — so the seam engine declares `SizeUnknown()`
    # and defines NO `length`. The `MethodError` is the contract; a `length` that
    # silently walked the halo to answer would be the thing the design forbids.
    @testset "the seam walk declares no length" begin
        for sys in SQUARE_SYSTEMS, base in (0, 1)
            l = base + 2
            l <= max_level(sys) || continue
            _, seam, _ = classify_roots(sys, base, l, Vertex())
            isempty(seam) && continue
            it = SubtreeHaloIterator(sys, first(seam), l)
            @test Base.IteratorSize(typeof(it)) isa Base.SizeUnknown
            @test_throws MethodError length(it)
        end
    end

    # Deeper than the sweep can afford at every root, on one root per system: a
    # 32x32 whole-face block from a level-0 root, where each of the four flush
    # sides is thirty-two cells long and only two of those were ever probed. If
    # the monotonicity argument were wrong anywhere, this is where the gap would
    # be widest.
    @testset "the seam walk at depth five" begin
        for sys in SQUARE_SYSTEMS, conn in (Vertex(), Edge())
            c = cellindex(levelgrid(sys, 0), 1)
            it = SubtreeHaloIterator(sys, c, 5; connectivity = conn)
            @test it.engine isa DGG.Fallbacks.SquareBandEngine
            @test collect(it) == forced_geometry_halo(sys, c, 5, conn)
            check_halo_case(sys, c, 5, conn)
        end
    end

    # THE DEEP REGIME, where the shallow sweep's arithmetic stops being the same
    # arithmetic.
    #
    # Everything above tops out at level 9 (`the band walk at 64x64`), and the
    # seam derivation's claim that the maps are affine at every `n` is an
    # argument plus an offline sweep, neither of which leaves a runnable artifact
    # in this repo. Two things can only go wrong deep:
    #
    #   * `FaceRect` stores its bounds as `Int32`, which holds a lattice
    #     coordinate through LEVEL 31 and overflows at 32. S2's `max_level` is
    #     30, so the headroom is exactly one level, and the failure mode is an
    #     `InexactError` thrown by `FaceRect`'s constructor from inside iterator
    #     construction — loud, but from a place that names neither the level nor
    #     the field. A `max_level` bump must be evaluated against 31, not against
    #     30 and not against `_SQUARE_CAP`; see `FaceRect`'s docstring.
    #   * The face-quadtree descent is 30 levels long here rather than nine, so a
    #     `code`/`x`/`y` restore that is off by a level has thirty chances to
    #     show rather than nine.
    #
    # A MAX-LEVEL CORNER BLOCK is the sharpest cheap case: root at
    # `max_level - 1`, target `max_level`, so the block is 2x2, flush on two
    # sides, and its corner is a face corner — three faces meet there on S2 and
    # the diagonal candidate does not exist, which is why the counts differ by
    # system. It costs a few hundred microseconds because the halo is a dozen
    # cells however deep the descent; the level-20 arm is the same shape one
    # decade shallower, so a failure that is about the DEPTH and not about the
    # corner separates the two.
    @testset "the band walk at max_level and at level 20" begin
        for sys in SQUARE_SYSTEMS
            mx = max_level(sys)
            for base in (mx - 1, 19), conn in (Vertex(), Edge())
                l = base + 1
                # Position 1 is lattice (0, 0) of face 0 under both curves: the
                # Morton systems because min-code is min-corner, S2 because the
                # Hilbert curve enters face 0 at its origin.
                c = cellindex(levelgrid(sys, base), 1)
                it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
                @test it.engine isa DGG.Fallbacks.SquareBandEngine
                @test it.engine.check isa DGG.Fallbacks.NativeCheck
                @test collect(it) == forced_geometry_halo(sys, c, l, conn)
                check_halo_case(sys, c, l, conn)
            end
        end
        # The counts, pinned so a walk that agreed with a wrong oracle would
        # still have to explain itself. A 2x2 corner block has five in-face band
        # cells; the rest come across the two seams, and the diagonal one exists
        # on HEALPix and not at an S2 cube corner or an ISEA4R icosahedral
        # vertex. `Edge()` drops all four diagonal contacts on every system.
        for (sys, nv) in ((HEALPixSystem(), 12), (S2System(), 11),
                          (ISEA4RSystem(), 11))
            mx = max_level(sys)
            c = cellindex(levelgrid(sys, mx - 1), 1)
            @test length(collect(SubtreeHaloIterator(sys, c, mx))) == nv
            @test length(collect(SubtreeHaloIterator(sys, c, mx;
                connectivity = Edge()))) == 8
        end
    end

    # -----------------------------------------------------------------------
    # The calibrated directed walk — H3 and IGeo7
    # -----------------------------------------------------------------------

    # WHY THE ORACLE HERE IS `forced_geometry_halo`, AND WHY IT IS A REAL ONE.
    #
    # `HexArcHaloEngine` enumerates by seeding each NEIGHBOUR's rim automaton
    # with an arc calibrated from observation and walking it down. It shares no
    # step with the generic descent — no `_admit`, no cap, no `rootcells` — and
    # it never asks about the root's own subtree at all. So comparing it against
    # the forced-geometry walk pins the enumeration and the adjacency predicate
    # together, end to end. `law_halo` runs on one arm below for the same reason
    # it does on the square walks: it decides which cells to consider by scanning
    # positions `1:ncells`, so a wrong walk and a wrong oracle cannot pass it
    # together.
    #
    # What must NOT be the only oracle is `subtree_border` — the seeded automaton
    # IS the border automaton, so the two would agree with a wrong arc for
    # precisely the reason it was wrong.
    #
    # THE TWO STRUCTURAL RISKS THIS SECTION EXISTS FOR.
    #
    #   * Pentagons. Twelve per level per system, and they are where the arc-3
    #     calibration lives: the only configuration in 52,182 measured pairs
    #     where the minimal covering arc is three directions wide is a pentagon
    #     neighbour whose deleted direction falls between the two touching
    #     children. They are pinned BY NAME (`ispentagon`, `z7_is_pentagon`)
    #     rather than found by degree, so a change to `irregular_cells` cannot
    #     quietly stop covering them, and the count is asserted to be twelve at
    #     every base.
    #   * Parity. The two systems' automata have exchanged parity branches and
    #     test their `L < 6` guards in the opposite order, so a seed that is
    #     right on H3 at an even level can be wrong on IGeo7 at the same level.
    #     Every arm therefore runs BOTH systems at consecutive root levels, which
    #     is what makes an even and an odd seeding of the same shape appear.

    HEX_SYSTEMS = (H3System(), IGeo7System())

    hex_ispentagon(sys, c) = sys isa DGG.H3System ?
        DGG.H3.ispentagon(c) : DGG.IGeo7.z7_is_pentagon(c.id)

    # The twelve pentagons of a level, by NAME. Deep levels have billions of
    # cells and must never be enumerated: the centre child of a pentagon is a
    # pentagon, so descending the first child from each level-0 pentagon names
    # the whole level-`base` pentagon set in twelve short walks.
    function hex_pentagons(sys, base::Int)
        grid0 = levelgrid(sys, 0)
        out = DGG.cellindextype(sys)[]
        for i in 1:ncells(grid0)
            c = cellindex(grid0, i)
            hex_ispentagon(sys, c) || continue
            for _ in 1:base
                c = first(DGG.children(sys, c))
            end
            hex_ispentagon(sys, c) && push!(out, c)
        end
        return out
    end

    # All twelve pentagons, the ring around two of them (the arc-3 neighbours),
    # and a spread of ordinary cells by position.
    function hex_roots(sys, base::Int, nhex::Int)
        grid = levelgrid(sys, base)
        pents = hex_pentagons(sys, base)
        @test length(pents) == 12
        n = ncells(grid)
        step = max(1, n ÷ nhex)
        rest = [cellindex(grid, i) for i in 1:step:n][1:min(nhex, end)]
        around = unique(vcat([collect(neighbors(grid, p, 1)) for p in pents[1:2]]...))
        return unique(vcat(pents, around, rest))
    end

    # Which walk a root gets is read OFF THE CODE UNDER TEST, exactly as
    # `classify_roots` does for the square systems. Three classes:
    #
    #   `child`    — `HexChildHaloEngine`: `target == level(root) + 1`, where the
    #                calibration is already the answer and no automaton runs.
    #   `arc`      — `HexArcHaloEngine`: one seeded rim automaton per neighbour.
    #   `fallback` — anything else, i.e. the generic outside-first walk.
    function classify_hex_roots(sys, roots, l::Int, conn)
        C = DGG.cellindextype(sys)
        child, arc, fallback = C[], C[], C[]
        for c in roots
            e = SubtreeHaloIterator(sys, c, l; connectivity = conn).engine
            if e isa DGG.Fallbacks.HexChildHaloEngine
                push!(child, c)
            elseif e isa DGG.Fallbacks.HexArcHaloEngine
                push!(arc, c)
            else
                push!(fallback, c)
            end
        end
        return child, arc, fallback
    end

    # HOW MANY ROOTS EACH WALK MUST CLAIM, and why a count is needed on top of
    # every oracle comparison below.
    #
    # `hex_halo_engine` has four guards that send a case back to the generic
    # walk, and none of them was observed to fire anywhere in the spike that
    # measured this design. A guard that grew stricter — a calibration that
    # started rejecting arc-3 pentagons, say — would therefore be INVISIBLE to
    # every oracle arm in this file: the rejected roots would fall to the generic
    # walk, which answers correctly, and every comparison would still pass
    # element for element. Only a count notices.
    #
    # The rule is not statistical: depth one is `HexChildHaloEngine` and every
    # deeper target is `HexArcHaloEngine`, for EVERY root of both systems at
    # every level, so the two claimed classes partition the generation and the
    # fallback class is empty. That is asserted exhaustively over the level-0 and
    # level-1 generations below, which is 122 + 842 H3 roots and 12 + 72 IGeo7
    # ones.
    function check_hex_classes(sys, roots, base::Int, l::Int, conn)
        child, arc, fallback = classify_hex_roots(sys, roots, l, conn)
        if l == base + 1
            @test length(child) == length(roots)
            @test isempty(arc)
        else
            @test isempty(child)
            @test length(arc) == length(roots)
        end
        @test isempty(fallback)
    end

    @testset "$(nameof(typeof(sys))): every root takes the directed walk" for sys in
            HEX_SYSTEMS
        for base in (0, 1), conn in (Vertex(), Edge())
            grid = levelgrid(sys, base)
            roots = [cellindex(grid, i) for i in 1:ncells(grid)]
            for l in (base + 1):min(base + 4, max_level(sys))
                check_hex_classes(sys, roots, base, l, conn)
            end
        end
        # AND ONE DEEP BASE, because the purpose of this testset is to notice a
        # guard that grew stricter and the guards are not all base-independent:
        # `_hex_validate` runs only from depth three, `_hex_calibrate` reads a
        # ring whose shape at base 0 is a whole base cell's and at base 8 an
        # ordinary hexagon's, and the seeded frames sit at the other parity. The
        # generation cannot be enumerated — H3's level-8 grid is 7e8 cells — so
        # the roots are the twelve pentagons BY NAME, the ring around two of them
        # (the arc-3 neighbours), and a positional spread. Depth 4 is included at
        # every base for the same reason: it is the first depth at which a seeded
        # arc has been through three transitions, and the shallow arms would not
        # notice a guard that only fires there.
        base = 8
        if base + 1 <= max_level(sys)
            roots = hex_roots(sys, base, 4)
            for conn in (Vertex(), Edge()),
                    l in (base + 1):min(base + 4, max_level(sys))
                check_hex_classes(sys, roots, base, l, conn)
            end
        end
    end

    # The differential sweep. Bases 0 and 1 are whole base cells and their
    # children; bases 5 and 8 are ordinary cells deep in the hierarchy, where the
    # ring is six neighbours of an ordinary hexagon and the seeded frames sit at
    # both parities. Depths 1-3 everywhere: depth 1 is the automaton-free path,
    # depth 2 is the first seeded walk and runs WITHOUT the depth-two validation
    # (`hex_halo_engine` skips it there), and depth 3 is the first depth that
    # validates.
    @testset "$(nameof(typeof(sys))): the directed walk against forced geometry" for
            sys in HEX_SYSTEMS
        for base in (0, 1, 5, 8)
            base + 1 <= max_level(sys) || continue
            roots = hex_roots(sys, base, 4)
            for c in roots, d in 1:3, conn in (Vertex(), Edge())
                l = base + d
                l <= max_level(sys) || continue
                it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
                @test collect(it) == forced_geometry_halo(sys, c, l, conn)
            end
        end
    end

    # THE EXACTNESS CONTRACT, WHICH NOTHING ELSE IN THIS FILE PINS.
    #
    # `HexArcHaloEngine` runs `_touches_subtree(IndexedNeighbors(), e, x)` on
    # every candidate before yielding it, and `halo.jl` says in as many words that
    # the check is what makes the engine EXACT rather than TRUSTED. But the
    # calibrated walk is already exact — candidate-to-halo ratio 1.0000 at every
    # depth — so the check never rejects anything, and deleting both call sites
    # leaves this file at its full pass count with nothing red. A future
    # simplification pass would find a green suite saying the check may go.
    #
    # So the check is pinned by making the band conservative on purpose. Widening
    # every calibrated arc `(L, s)` to `(L + 1, s - 1)` keeps the original arc
    # inside the new one, so the widened walk is a SUPERSET of the calibrated
    # one — asserted below by counting the raw automaton output, which is a
    # strict inequality and not a pinned number — and the engine's answer must be
    # unchanged, because the check filters the surplus back out. Delete either
    # `_touches_subtree` call in `HexArcHaloEngine` and this testset fails; every
    # other arm in this file stays green.
    #
    # Deliberately NOT a count of the surplus: the point is the invariant "the
    # emitted set is the halo whatever the band proposes", and a number would pin
    # the widening rather than the check.
    function widen_hex_arcs(e)
        ring = e.ring
        for i in 1:length(ring)
            h = ring[i]
            ring = DGG.Helpers.small_setindex(ring,
                DGG.Fallbacks.HexNeighbour(h.cell, h.lo, h.arclen + Int8(1),
                    Int8(mod(Int(h.start) - 1, 6))), i)
        end
        return DGG.Fallbacks.HexArcHaloEngine(e.system, e.grid, e.root,
            e.rootlevel, e.target, e.connectivity, ring)
    end

    # The engine's candidate stream BEFORE the check: the seeded automata alone,
    # which is what "conservative band" names.
    function hex_candidate_count(e)
        n = 0
        for i in 1:length(e.ring)
            nb = e.ring[i]
            for _ in DGG.seeded_rim_engine(e.system, nb.cell, e.target,
                    Int(nb.arclen), Int(nb.start))
                n += 1
            end
        end
        return n
    end

    @testset "$(nameof(typeof(sys))): the check filters a widened arc" for sys in
            HEX_SYSTEMS
        for base in (0, 1, 2), d in 2:3, conn in (Vertex(), Edge())
            l = base + d
            l <= max_level(sys) || continue
            for c in spread(hex_roots(sys, base, 3), 5)
                it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
                @test it.engine isa DGG.Fallbacks.HexArcHaloEngine
                wide = widen_hex_arcs(it.engine)
                # The band really did widen, so the equality below is the check
                # doing work and not the mutation being a no-op.
                @test hex_candidate_count(wide) > hex_candidate_count(it.engine)
                @test collect(SubtreeHaloIterator(sys, c, l, conn, wide)) ==
                      collect(it)
            end
        end
    end

    # The contract bundle — uniqueness, sortedness, outside ancestry, both
    # adjacency directions — on a smaller spread, because it builds two `Set`s
    # and a `subtree_border` per case. Pentagons first, so the expensive laws run
    # where the walk is least ordinary.
    @testset "$(nameof(typeof(sys))): the directed walk keeps the contract" for sys in
            HEX_SYSTEMS
        for base in (0, 2), conn in (Vertex(), Edge())
            base + 1 <= max_level(sys) || continue
            for c in spread(hex_roots(sys, base, 2), 6),
                    l in (base + 1):min(base + 3, max_level(sys))
                check_halo_case(sys, c, l, conn)
            end
        end
    end

    # Depth four, where the halo is 246 cells around a hexagon and 205 around a
    # pentagon, and where a seeded arc has been through three transitions rather
    # than one. The oracle costs about 50 ms a call here, so this is two roots —
    # a pentagon and a hexagon — per system rather than a sweep.
    @testset "$(nameof(typeof(sys))): the directed walk at depth four" for sys in
            HEX_SYSTEMS
        grid = levelgrid(sys, 2)
        pent = first(hex_pentagons(sys, 2))
        hex = first(filter(c -> !hex_ispentagon(sys, c),
            collect(neighbors(grid, pent, 1))))
        for c in (pent, hex), conn in (Vertex(), Edge())
            it = SubtreeHaloIterator(sys, c, 6; connectivity = conn)
            @test it.engine isa DGG.Fallbacks.HexArcHaloEngine
            @test collect(it) == forced_geometry_halo(sys, c, 6, conn)
        end
    end

    # One arm against the O(ncells) brute force, which shares nothing with either
    # the directed walk or the geometry walk. Level-0 roots at depths 1 and 2,
    # where the target grids are at most 5882 cells on H3 and 492 on IGeo7.
    @testset "$(nameof(typeof(sys))): the directed walk against the law" for sys in
            HEX_SYSTEMS
        for c in spread(hex_roots(sys, 0, 4), 6), d in 1:2, conn in (Vertex(), Edge())
            @test collect(SubtreeHaloIterator(sys, c, d; connectivity = conn)) ==
                  law_halo(sys, c, d; connectivity = conn)
        end
    end

    # THE COUNTS, PINNED BUT NOT PROMISED. `3^(d+1) + 3` around a hexagon and
    # `5(3^d + 1)/2` around a pentagon — verified in 176/176 configurations by
    # the spike and derived from a per-neighbour census of `(3^d + 1)/2` that
    # holds for both arc lengths, but derived by ENUMERATION rather than from the
    # transition recurrence. So the numbers are pinned here against a real
    # `collect`, and `IteratorSize` still says `SizeUnknown()` with no `length`
    # method at all — the next testset is what says so. Pinning them is not the
    # same as promising them: this arm turns "the census changed" into a failure
    # instead of a silent drift, which is what the design asks for while a count
    # is evidence.
    @testset "$(nameof(typeof(sys))): the directed walk's census" for sys in HEX_SYSTEMS
        for base in (0, 3)
            grid = levelgrid(sys, base)
            pents = hex_pentagons(sys, base)
            hexes = filter(c -> !hex_ispentagon(sys, c),
                collect(neighbors(grid, first(pents), 1)))
            # IGeo7's base tessellation is twelve pentagons and nothing else, so
            # at base 0 there is no hexagon to take; base 3 supplies both.
            cells = isempty(hexes) ? [first(pents)] : [first(pents), first(hexes)]
            for d in 1:4
                base + d <= max_level(sys) || continue
                for c in cells
                    n = length(collect(SubtreeHaloIterator(sys, c, base + d)))
                    @test n == (hex_ispentagon(sys, c) ? (5 * (3^d + 1)) ÷ 2 :
                                3^(d + 1) + 3)
                end
            end
        end
    end

    # The count contract in the negative, as for the seam walk. The census above
    # is evidence, not an API promise, so neither hex engine declares a length
    # and the `MethodError` is the contract being kept.
    @testset "the directed walk declares no length" begin
        for sys in HEX_SYSTEMS, d in 1:2
            c = cellindex(levelgrid(sys, 1), 1)
            it = SubtreeHaloIterator(sys, c, 1 + d)
            @test Base.IteratorSize(typeof(it)) isa Base.SizeUnknown
            @test_throws MethodError length(it)
        end
    end

    # THE LAZINESS LAW, which is the other thing Task 6 delivers. The generic
    # walk descends from `rootcells` and prunes by cap, so reaching the first
    # halo cell of a deep target costs a traversal that grows with the target
    # level: taking ten cells of an IGeo7 L=5, d=7 halo cost 42 ms and 779 KB of
    # allocation — allocation in proportion to the whole halo, which the design's
    # verification section forbids in as many words. The directed walk holds one
    # seeded automaton and its frame stack, both isbits, so the same prefix is a
    # flat 256 bytes at every depth.
    #
    # The assertion that matters is INDEPENDENCE, not a threshold: the same
    # number of bytes at depth three, five and seven is the `O(depth)` claim, and
    # a walk that started materialising would break it at once.
    @testset "a prefix of a deep halo costs O(depth), not O(halo)" begin
        prefix10 = it -> begin
            n = 0
            for _ in it
                n += 1
                n >= 10 && break
            end
            n
        end
        for sys in HEX_SYSTEMS
            base = 5
            grid = levelgrid(sys, base)
            root = cellindex(grid, ncells(grid) ÷ 2 + 1)
            depths = filter(d -> base + d <= max_level(sys), [3, 5, 7])
            allocs = map(depths) do d
                prefix10(SubtreeHaloIterator(sys, root, base + d))     # compile
                @allocated prefix10(SubtreeHaloIterator(sys, root, base + d))
            end
            @test all(a -> a <= 4096, allocs)
            @test length(unique(allocs)) == 1
            @test all(d -> prefix10(SubtreeHaloIterator(sys, root, base + d)) == 10,
                depths)
        end
        # And the walk it replaced, on the system where the violation was
        # measured: the generic prefix allocates hundreds of kilobytes where the
        # directed one allocates hundreds of bytes.
        sys = IGeo7System()
        grid = levelgrid(sys, 5)
        root = cellindex(grid, ncells(grid) ÷ 2 + 1)
        gen = () -> prefix10(SubtreeHaloIterator(sys, root, 12, Vertex(),
            DGG.Fallbacks.generic_halo_engine(sys, root, 12, Vertex())))
        gen()
        dir = () -> prefix10(SubtreeHaloIterator(sys, root, 12))
        dir()
        @test @allocated(gen()) > 100 * @allocated(dir())
    end

    # -----------------------------------------------------------------------
    # The generic walk, still oracled where it is still the walk
    # -----------------------------------------------------------------------

    # `check_root_classes` and `check_hex_classes` assert that no root of any
    # system falls back, which says the specializations were reached but says
    # nothing about the walk they would have fallen back TO. That walk is not
    # dead code: it is what every one of `hex_halo_engine`'s and
    # `_seam_band_engine`'s guards returns, and what a system added later
    # inherits until it writes an engine of its own. Since Task 6 no system
    # reaches it through the keyword constructor at all, so it is BUILT here and
    # held to both oracles — the geometry walk, and the `O(ncells)` law, which
    # enumerates from an ascending position scan and shares nothing with either.
    #
    # Without this arm the generic walk would be exercised only as the oracle's
    # own carrier: `forced_geometry_halo` runs the same `OutsideWalkEngine` under
    # a different adjacency provider, so a bug in the enumeration or the cap
    # prune would move both sides together and show up nowhere. (The one other
    # place it is still pinned is the depth-4 aperture-7 arm of "the defining
    # law", which is where an under-covering root cap is caught.)
    @testset "the generic fallback still agrees with the oracle" begin
        for sys in HEX_SYSTEMS, base in (0, 1), conn in (Vertex(), Edge())
            grid = levelgrid(sys, base)
            roots = unique(vcat(sample_cells(grid, 3), irregular_cells(grid, 2)))
            l = base + 1
            for c in roots
                it = SubtreeHaloIterator(sys, c, l, conn,
                    DGG.Fallbacks.generic_halo_engine(sys, c, l, conn))
                @test it.engine isa DGG.Fallbacks.OutsideWalkEngine
                h = collect(it)
                @test h == forced_geometry_halo(sys, c, l, conn)
                @test h == law_halo(sys, c, l; connectivity = conn)
            end
        end
    end

    # -----------------------------------------------------------------------
    # Resumability, on every engine this file can reach
    # -----------------------------------------------------------------------

    # THE PREFIX-EQUALITY LAW, mirroring `subtree_iterators.jl`'s: "the walk is
    # resumable, not restarted". Four cells taken off the front must be the first
    # four of the collected form, and a second `collect` must reproduce the whole
    # walk — which is the observable half of the house rule that an engine is
    # immutable and ALL walk state travels in the value `iterate` threads. An
    # engine that cached a cursor in a mutable field would pass every oracle in
    # this file and fail exactly here, on the second pass.
    #
    # Run against EVERY engine type rather than one per system: the seven walks
    # thread seven different state types — a selection emit, a frame stack, a
    # position counter, a quadtree descent under two emit rules, a child cursor,
    # and a seeded automaton plus its own stack — and resumability is a property
    # of the state, not of the system. The set assertion at the end is what says
    # the list did not quietly stop covering one.
    @testset "the walk is resumable, not restarted" begin
        seen = Set{Symbol}()
        # `SquareBandEngine` is two walks wearing one name, so the tag reads the
        # emit rule as well: the exact band and the seam band differ in what they
        # do between yields, which is exactly what resumability is about.
        engine_tag(e) = e isa DGG.Fallbacks.SquareBandEngine ?
            (e.check isa DGG.Fallbacks.NoCheck ? :SquareBandNoCheck :
             :SquareBandNativeCheck) : nameof(typeof(e))
        function check_prefix(it)
            push!(seen, engine_tag(it.engine))
            full = collect(it)
            @test length(full) >= 4
            prefix = eltype(it)[]
            for x in it
                push!(prefix, x)
                length(prefix) >= 4 && break
            end
            @test prefix == full[1:4]
            @test collect(it) == full        # a second pass gives the same walk
        end
        for sys in systems()
            c = cellindex(levelgrid(sys, 1), 1)
            l = min(level(c) + 2, max_level(sys))
            check_prefix(SubtreeHaloIterator(sys, c, l))              # shipped
            check_prefix(SubtreeHaloIterator(sys, c, level(c)))       # one-ring
            # The generic walk is no longer reachable through the keyword
            # constructor on any system, so it is built explicitly — and on A5
            # that same call is the scan, which is how both fallbacks are covered
            # without naming either system.
            check_prefix(SubtreeHaloIterator(sys, c, l, Vertex(),
                DGG.Fallbacks.generic_halo_engine(sys, c, l, Vertex())))
        end
        # Depth one on the aperture-7 systems is the automaton-free child walk,
        # which the `l = level(c) + 2` cases above never reach.
        for sys in HEX_SYSTEMS
            check_prefix(SubtreeHaloIterator(sys, cellindex(levelgrid(sys, 1), 1), 2))
        end
        # Both square emit rules, on the blocks the constructor actually claims.
        for sys in SQUARE_SYSTEMS
            inface, seam, _ = classify_roots(sys, 2, 4, Vertex())
            check_prefix(SubtreeHaloIterator(sys, first(inface), 4))
            check_prefix(SubtreeHaloIterator(sys, first(seam), 4))
        end
        @test seen == Set((:RingHaloEngine, :OutsideWalkEngine, :ScanHaloEngine,
            :SquareBandNoCheck, :SquareBandNativeCheck, :HexChildHaloEngine,
            :HexArcHaloEngine))
    end

    # -----------------------------------------------------------------------
    # The wrapper, and the guard on a lying count
    #
    # Last rather than first, because the wrapper test needs `classify_roots`: a
    # level-0 root is flush on all four sides on every square system, so a forwarding
    # test that only used one would never reach a specialized engine at all.
    # -----------------------------------------------------------------------

    # The authalic transform moves where a cell is DRAWN, not which cells are
    # adjacent, so the halo through the wrapper must be the halo without it — the
    # same cell ids, in the same order. `halo_engine(::AuthalicSystem, ...)` is one
    # forwarding line, and this is what says the line is there.
    @testset "AuthalicSystem forwards the halo walk" begin
        for sys in systems()
            wrapped = DGG.AuthalicSystem(sys)
            grid0 = levelgrid(sys, 0)
            c = cellindex(grid0, 1)
            for l in level(c):min(level(c) + 2, max_level(sys))
                @test collect(SubtreeHaloIterator(wrapped, c, l)) ==
                      collect(SubtreeHaloIterator(sys, c, l))
            end
        end
        # And on a root the SPECIALIZATION claims. Forwarding that only ever ran on
        # a level-0 root would be forwarding that only ever reached the generic
        # walk — the wrapper would be free to lose the fast path, and every
        # assertion above would still pass.
        for sys in SQUARE_SYSTEMS
            wrapped = DGG.AuthalicSystem(sys)
            for base in BAND_BASES
                l = base + 2
                l <= max_level(sys) || continue
                inface, _, _ = classify_roots(sys, base, l, Vertex())
                isempty(inface) && continue
                c = first(inface)
                for conn in (Vertex(), Edge())
                    it = SubtreeHaloIterator(wrapped, c, l; connectivity = conn)
                    @test it.engine isa DGG.Fallbacks.SquareBandEngine
                    @test collect(it) ==
                          collect(SubtreeHaloIterator(sys, c, l; connectivity = conn))
                end
            end
        end
        # And on the aperture-7 specialization, for the same reason.
        for sys in HEX_SYSTEMS
            wrapped = DGG.AuthalicSystem(sys)
            c = cellindex(levelgrid(sys, 1), 1)
            for d in 1:2, conn in (Vertex(), Edge())
                it = SubtreeHaloIterator(wrapped, c, 1 + d; connectivity = conn)
                @test it.engine isa (d == 1 ? DGG.Fallbacks.HexChildHaloEngine :
                                     DGG.Fallbacks.HexArcHaloEngine)
                @test collect(it) ==
                      collect(SubtreeHaloIterator(sys, c, 1 + d; connectivity = conn))
            end
        end
    end

    @testset "collect is the guarded path" begin
        sys = HEALPixSystem()
        c = cellindex(levelgrid(sys, 1), 1)
        lying = SubtreeHaloIterator(sys, c, 1, Vertex(), MiscountingEngine())
        @test_throws ErrorException collect(lying)
    end

end  # @testset "subtree halos"

end # module
