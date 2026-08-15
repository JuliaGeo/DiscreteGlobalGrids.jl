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

# ---------------------------------------------------------------------------
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
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# The wrapper, and the guard on a lying count
# ---------------------------------------------------------------------------

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
end

# Claims three, yields one — exactly the shape `collect_subtree` exists to
# catch. Without it `collect`'s own `HasLength` route sizes the vector from the
# claim and hands back two `undef` slots as cell ids; the specializations below
# are the engines that claim a closed-form count, so this guard is load-bearing
# for them and not a curiosity.
struct MiscountingEngine end
Base.iterate(::MiscountingEngine) = (DGG.LevelIndex(0, 0), 1)
Base.iterate(::MiscountingEngine, ::Int) = nothing
Base.eltype(::Type{MiscountingEngine}) = DGG.LevelIndex
Base.IteratorSize(::Type{MiscountingEngine}) = Base.HasLength()
Base.length(::MiscountingEngine) = 3

@testset "collect is the guarded path" begin
    sys = HEALPixSystem()
    c = cellindex(levelgrid(sys, 1), 1)
    lying = SubtreeHaloIterator(sys, c, 1, Vertex(), MiscountingEngine())
    @test_throws ErrorException collect(lying)
end

# ---------------------------------------------------------------------------
# The square band walk — HEALPix, S2 and ISEA4R away from a face edge
# ---------------------------------------------------------------------------

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

# Which roots the band walk claims is a fact about the lattice, not something a
# test should hard-code: a block is claimed exactly when it is nowhere flush
# with its face edge, which at base 2 is four of the sixteen cells per face on
# all three systems. So roots are CLASSIFIED by the engine the constructor
# actually chose, and both classes are checked — the claimed ones because that
# is the walk under test, the unclaimed ones because a fallback that quietly
# stopped agreeing would otherwise be invisible.
#
# A level-0 block is flush on all four sides and a level-1 block on two per
# axis, so the engine is reachable only for roots at level >= 2. Base 2 is
# therefore the shallowest base this section can use — and, on its own, not a
# sufficient one. See `BAND_BASES`.
function classify_roots(sys, base::Int, l::Int, conn)
    grid = levelgrid(sys, base)
    C = DGG.cellindextype(sys)
    band, fallback = C[], C[]
    for i in 1:ncells(grid)
        c = cellindex(grid, i)
        it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
        push!(it.engine isa DGG.Fallbacks.SquareBandEngine ? band : fallback, c)
    end
    return band, fallback
end

const SQUARE_SYSTEMS = (HEALPixSystem(), S2System(), ISEA4RSystem())

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
const BAND_BASES = (2, 3)

@testset "$(nameof(typeof(sys))): the band walk against forced geometry" for sys in
        SQUARE_SYSTEMS
    lazies = 0
    for base in BAND_BASES, d in 1:2, conn in (Vertex(), Edge())
        l = base + d
        l <= max_level(sys) || continue
        band, fallback = classify_roots(sys, base, l, conn)
        @test !isempty(band)                    # the specialization was reached
        lazies += length(band)
        # Six claimed roots: enough to cross faces (four per face, so this
        # reaches at least two of them) without paying for the oracle 64 times.
        for c in band[1:min(6, length(band))]
            it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
            @test it.engine isa DGG.Fallbacks.SquareBandEngine
            @test collect(it) == forced_geometry_halo(sys, c, l, conn)
            check_halo_case(sys, c, l, conn)
        end
        # And three flush roots, which must still take the generic walk and must
        # still agree with the oracle: the guard is a runtime interval test, and
        # a guard that admitted a flush block would be caught here.
        for c in fallback[1:min(3, length(fallback))]
            it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
            @test !(it.engine isa DGG.Fallbacks.SquareBandEngine)
            @test collect(it) == forced_geometry_halo(sys, c, l, conn)
        end
    end
    @test lazies > 0
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
        band, _ = classify_roots(sys, base, l, conn)
        for c in band[1:min(4, length(band))]
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
            band, _ = classify_roots(sys, base, l, conn)
            isempty(band) && continue
            for c in band[1:min(3, length(band))]
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
    band, _ = classify_roots(sys, 3, l, Vertex())
    @test !isempty(band)
    for c in band[1:min(3, length(band))], conn in (Vertex(), Edge())
        it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
        @test it.engine isa DGG.Fallbacks.SquareBandEngine
        @test length(it) == (conn isa Vertex ? 260 : 256)
        @test collect(it) == forced_geometry_halo(sys, c, l, conn)
    end
end

# The corner law. `Edge()` drops exactly the four cells that touch the block at
# a vertex only, so the two halos differ by four and nothing else. Pinned to a
# NON-FLUSH block deliberately: a flush block can lose a corner across a cube
# corner, where three faces meet and the diagonal neighbour does not exist, so
# the count would be three there and the law would read as broken.
@testset "Edge drops exactly the four diagonal corners" begin
    sys = S2System()
    band, _ = classify_roots(sys, 2, 4, Vertex())
    c = first(band)
    hv = collect(SubtreeHaloIterator(sys, c, 4; connectivity = Vertex()))
    he = collect(SubtreeHaloIterator(sys, c, 4; connectivity = Edge()))
    @test length(setdiff(Set(hv), Set(he))) == 4
    @test issubset(Set(he), Set(hv))
end

end # module
