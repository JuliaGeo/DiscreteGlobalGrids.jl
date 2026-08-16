# Cross-system laws for the outside boundary of a subtree: target-level cells
# outside the subtree that touch one of its descendants.
#
# The suite compares halo iterators with an ascending full-grid scan and a
# forced-geometry oracle. It also checks ordering, uniqueness, adjacency,
# container delegation, iterator protocol, and allocation behavior. Generic
# laws sweep `systems()`; specialization-specific laws use explicit system sets.

module SubtreeHaloTests

using Test
using DiscreteGlobalGrids
using DiscreteGlobalGrids: systems, levelgrid, level, max_level, ncells,
    cellindex, cellposition, neighbors, ancestor, subtree_border, Vertex, Edge,
    SubtreeHaloIterator, subtree_halo, halo_positions, halo_sizehint
import DiscreteGlobalGrids as DGG

# Fixture types must be declared outside the local scope of the outer testset.

struct MiscountingEngine end
Base.iterate(::MiscountingEngine) = (DGG.LevelIndex(0, 0), 1)
Base.iterate(::MiscountingEngine, ::Int) = nothing
Base.eltype(::Type{MiscountingEngine}) = DGG.LevelIndex
Base.IteratorSize(::Type{MiscountingEngine}) = Base.HasLength()
Base.length(::MiscountingEngine) = 3

# Eager reference with the same values and iterator protocol as its wrapped
# engine. It materializes during construction and again for each fresh walk, so
# allocation laws can distinguish eager and incremental implementations without
# changing the enumeration.
struct EagerHaloEngine{E,C}
    inner::E
    cells::Vector{C}
end

EagerHaloEngine(inner) = EagerHaloEngine(inner, DGG.collect_subtree(inner))

Base.iterate(e::EagerHaloEngine{E,C}) where {E,C} =
    iterate(e, (DGG.collect_subtree(e.inner), 1))
Base.iterate(::EagerHaloEngine{E,C}, s::Tuple{Vector{C},Int}) where {E,C} =
    s[2] > length(s[1]) ? nothing : (@inbounds(s[1][s[2]]), (s[1], s[2] + 1))
Base.eltype(::Type{<:EagerHaloEngine{E,C}}) where {E,C} = C
Base.IteratorSize(::Type{<:EagerHaloEngine}) = Base.HasLength()
Base.length(e::EagerHaloEngine) = length(e.cells)

# Subset wrapper that counts membership-position and span queries. This measures
# traversal work without relying on wall-clock time.
mutable struct CountingSubset{S}
    inner::S
    calls::Int
end

CountingSubset(inner) = CountingSubset{typeof(inner)}(inner, 0)

DGG.cellposition(cs::CountingSubset, c::DGG.AbstractCellIndex) =
    (cs.calls += 1; DGG.cellposition(cs.inner, c))

DGG.Fallbacks.subset_span(cs::CountingSubset, lo::Int, hi::Int) =
    (cs.calls += 1; DGG.Fallbacks.subset_span(cs.inner, lo, hi))

# `wrap` optionally replaces the shipped walk with its eager equivalent.
counting_iterator(sys, cs::CountingSubset, complete, l, conn, wrap = identity) =
    DGG.Fallbacks.SubsetHaloIterator(cs, conn, wrap(
        DGG.Fallbacks.subset_halo_engine(sys, cs, complete, Int(l), conn)))

counting_halo(sys, cs::CountingSubset, complete, l, conn) =
    DGG.collect_subtree(counting_iterator(sys, cs, complete, l, conn))

# Return query counts for an `n`-cell prefix and for complete collection. Reset
# the counter after construction so this helper measures traversal only.
function counted_walk(it, cs::CountingSubset, n::Int)
    cs.calls = 0
    take_n(it, n)
    prefix = cs.calls
    cs.calls = 0
    h = DGG.collect_subtree(it)
    return (prefix, cs.calls, h)
end


take_n(it, n::Int) = (seen = 0; for _ in it
    seen += 1
    seen >= n && break
end; seen)

build_and_take(sys, c, l, n::Int) = take_n(SubtreeHaloIterator(sys, c, l), n)

lazy_bytes(sys, c, l, n::Int) =
    (build_and_take(sys, c, l, n); @allocated build_and_take(sys, c, l, n))
eager_bytes(sys, c, l) = (subtree_halo(sys, c, l); @allocated subtree_halo(sys, c, l))

const SINK = Ref{Any}(nothing)

construct!(sys, c, l) = (SINK[] = SubtreeHaloIterator(sys, c, l); nothing)
construct_bytes(sys, c, l) = (construct!(sys, c, l); @allocated construct!(sys, c, l))

subset_construct!(sub) = (SINK[] = halo(sub); nothing)
subset_construct_bytes(sub) =
    (subset_construct!(sub); @allocated subset_construct!(sub))

# Build the generic outside-first walk explicitly; specialized constructors do
# not select it.
generic_iterator(sys, c, l) = SubtreeHaloIterator(sys, c, Int(l), Vertex(),
    DGG.Fallbacks.generic_halo_engine(sys, c, Int(l), Vertex()))

generic_take(sys, c, l, n::Int) = take_n(generic_iterator(sys, c, l), n)
generic_collect(sys, c, l) = DGG.collect_subtree(generic_iterator(sys, c, l))
generic_construct!(sys, c, l) = (SINK[] = generic_iterator(sys, c, l); nothing)
generic_construct_bytes(sys, c, l) =
    (generic_construct!(sys, c, l); @allocated generic_construct!(sys, c, l))

# Apply the same measurements to the eager fixture.
fixture_iterator(sys, c, l) = SubtreeHaloIterator(sys, c, Int(l), Vertex(),
    EagerHaloEngine(DGG.Fallbacks.halo_engine(sys, c, Int(l), Vertex())))

fixture_collect(sys, c, l) = DGG.collect_subtree(fixture_iterator(sys, c, l))
fixture_construct!(sys, c, l) = (SINK[] = fixture_iterator(sys, c, l); nothing)
fixture_take(sys, c, l, n::Int) = take_n(fixture_iterator(sys, c, l), n)
fixture_construct_bytes(sys, c, l) =
    (fixture_construct!(sys, c, l); @allocated fixture_construct!(sys, c, l))
fixture_prefix_bytes(sys, c, l, n::Int) =
    (fixture_take(sys, c, l, n); @allocated fixture_take(sys, c, l, n))
fixture_collect_bytes(sys, c, l) =
    (fixture_collect(sys, c, l); @allocated fixture_collect(sys, c, l))

# The level-0 cell these laws hang a subtree from.
#
# It is root 1 on thirteen of the fifteen systems, and has to be asked for
# rather than written as `law_root(sys)` because of the other
# two: ISEA3H and ISEA4H seal their two polar roots as single-child chains, so
# root 1 has a ONE-cell subtree at every depth. Every law below that reads a
# growth ratio, an allocation trend, or a halo big enough to slice reads it as
# 5 — the pentagon's one-ring — at every level, and fails on a degeneracy
# rather than on the property it means to test.
#
# Picking the first root whose subtree actually branches keeps all thirteen
# unchanged (so the constants pinned against them still mean what they meant)
# and gives the two hexagon systems a representative subtree.
law_root(sys) = law_root(sys, 4)
function law_root(sys, probe::Int)
    g0 = levelgrid(sys, 0)
    l = min(probe, max_level(sys))
    for i in 1:ncells(g0)
        c = cellindex(g0, i)
        count(_ -> true, DGG.descendants(sys, c, l)) > 1 && return c
    end
    return cellindex(g0, 1)
end

# The depth ladders below are written as level NUMBERS, which silently means
# wildly different amounts of WORK once the registry spans apertures 3 to 9:
# level 7 is a 16k-cell subtree at aperture 4 and a 4.8M-cell one at aperture 9.
# That matters here because the four aperture-9 systems have no specialized halo
# engine — their walk is the generic one, which costs what the full-grid
# reference scan costs (measured at depth 7 on rHEALPix: 4.46s for the scan,
# 4.46s for the shipped walk, against 0.006s for S2, which has one) — so the
# cost lands on every law in this file, and the ladders are walked repeatedly.
#
# Cap the ladder at the largest subtree these laws already walked before the
# literature systems were registered: H3 at level 7. Every system main swept is
# under it and keeps its exact depths; only the aperture-9 newcomers are pulled
# in, and their ladders still grow by a factor of 9 per rung, which is more than
# any growth ratio asserted here needs.
const LAW_SUBTREE_BUDGET = 7^7

function law_depths(sys, ladder)
    b = max(2, length(children(sys, law_root(sys))))
    top = floor(Int, log(LAW_SUBTREE_BUDGET) / log(b))
    # The other end of the same problem. These ladders are written for aperture
    # 4, and a level moves the halo by roughly `sqrt(aperture)`, so at aperture 3
    # the same span of levels grows the halo 4.5x where the law asks for 5x —
    # a fail on arithmetic rather than on behaviour. Stretch the top rung by the
    # levels the coarser aperture needs to cover the same ground; systems at
    # aperture 4 and above are already there and keep their exact ladders.
    stretch = max(0, ceil(Int, (last(ladder) - first(ladder)) * (log(4) / log(b) - 1)))
    stretched = (ladder[1:(end - 1)]..., last(ladder) + stretch)
    return unique(min(d, top) for d in stretched)
end

# One outer testset records failures while allowing later nested testsets to run.
@testset "subtree halos" begin

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

    # All engines expose the same validation messages through the public iterator.
    @testset "level validation" begin
        msg(f) = try
            f()
            nothing
        catch e
            e isa ArgumentError ? e.msg : rethrow()
        end
        for sys in systems()
            grid = levelgrid(sys, 1)
            c = cellindex(grid, 1)
            mx = max_level(sys)
            @test msg(() -> SubtreeHaloIterator(sys, c, 0)) ==
                "subtree_halo: level 0 is above the cell's own level 1"
            @test msg(() -> SubtreeHaloIterator(sys, c, mx + 1)) ==
                "subtree_halo: level $(mx + 1) is past max_level $mx"
        end
    end

    # Full-grid oracle: scan target-level positions and retain outside cells with
    # a descendant neighbour. This O(ncells) definition shares no traversal code
    # with the halo engines and preserves canonical position order.
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

    # Deterministic, evenly spaced root sample.
    sample_cells(grid, n::Int) =
        [cellindex(grid, i) for i in 1:max(1, ncells(grid) ÷ n):ncells(grid)]

    # Select cells whose one-ring size differs from the grid's modal degree.
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


    @testset "$(nameof(typeof(sys))): the defining law" for sys in systems()
        grid0 = levelgrid(sys, 0)
        n0 = ncells(grid0)
        mx = max_level(sys)

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

    # The geometry engine decides adjacency by counting shared boundary points,
    # so it states a law about the DRAWN rings and not about the lattice. The
    # two central-place hexagon systems draw their Snyder edges as a documented
    # finite approximation (see `src/systems/IVEARTEA/system.jl` for the same
    # tradeoff stated for the rhombic family, and the ISEA3H/4H note in
    # `regridding_conservation.jl`), so on them the two relations are not the
    # same relation and geometry comes out a strict subset of topology:
    # measured over this file's own roots, ISEA4H disagrees in 6 of the level-1
    # cases and 5 of the level-2 ones, ISEA3H in 1 level-2 case, and every
    # disagreement is a cell topology holds and geometry drops. The other
    # thirteen registered systems — the seven other literature ones included —
    # agree exactly, which is what keeps this an exclusion rather than a
    # weakened law. `stencils.jl` excludes the same two from its geometric
    # member-adjacency law for the related reason that they refine
    # non-congruently.
    exact_rings(sys) = !(sys isa DGG.ISEA3HSystem || sys isa DGG.ISEA4HSystem)

    @testset "forced geometry at depth zero" begin
        for sys in filter(exact_rings, systems())
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

    @testset "$(nameof(typeof(sys))): geometry agrees with topology" for sys in
                                                                        filter(exact_rings, systems())
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

    function check_halo_case(sys, c, l, conn)
        it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
        h = try
            collect(it)
        catch err
            @test err === nothing     # fails naming the throw, rather than aborting
            return eltype(it)[]
        end
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
        return h
    end

    @testset "$(nameof(typeof(sys))) at level $base" for sys in systems(),
            base in sweep_bases(sys)
        for c in sweep_roots(sys, base), l in base:min(base + 2, max_level(sys))
            hv = check_halo_case(sys, c, l, Vertex())
            he = check_halo_case(sys, c, l, Edge())
            @test issubset(Set(he), Set(hv))
        end
    end

    @testset "$(nameof(typeof(sys))) at depth" for sys in
            filter(s -> !(s isa DGG.A5System), systems())
        d = deep_depth(sys, 0)
        d >= 1 || continue
        check_halo_case(sys, law_root(sys), d, Vertex())
    end

    # -----------------------------------------------------------------------
    # The square band walk — HEALPix, S2 and ISEA4R away from a face edge
    # -----------------------------------------------------------------------


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

    function spread(v, n::Int)
        isempty(v) && return v
        length(v) <= n && return v
        step = length(v) ÷ n
        return [v[1 + (i - 1) * step] for i in 1:n]
    end

    inface_root_count(sys, base::Int) =
        ncells(levelgrid(sys, 0)) * max(0, (1 << base) - 2)^2

    function check_root_classes(sys, base::Int, inface, seam, fallback)
        total = ncells(levelgrid(sys, base))
        @test length(inface) == inface_root_count(sys, base)
        @test length(seam) == total - inface_root_count(sys, base)
        @test isempty(fallback)
    end

    SQUARE_SYSTEMS = (HEALPixSystem(), S2System(), ISEA4RSystem())

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
            for c in spread(inface, 6)
                it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
                @test it.engine isa DGG.Fallbacks.SquareBandEngine
                @test it.engine.check isa DGG.Fallbacks.NoCheck
                @test collect(it) == forced_geometry_halo(sys, c, l, conn)
                check_halo_case(sys, c, l, conn)
            end
        end
    end

    # One arm against the `O(ncells)` brute force as well: the one check a wrong
    # band and a wrong oracle cannot pass together. Depth 2 only, where a visit
    # to every cell of the target level is still cheap.
    @testset "$(nameof(typeof(sys))): the band walk against the law" for sys in
            SQUARE_SYSTEMS
        for base in BAND_BASES, conn in (Vertex(), Edge())
            l = base + 2
            l <= max_level(sys) || continue
            inface, _, _ = classify_roots(sys, base, l, conn)
            for c in spread(inface, 4)
                @test collect(SubtreeHaloIterator(sys, c, l; connectivity = conn)) ==
                      law_halo(sys, c, l; connectivity = conn)
            end
        end
    end

    @testset "$(nameof(typeof(sys))): the band count is closed form" for sys in
            SQUARE_SYSTEMS
        for base in BAND_BASES, d in 1:4
            l = base + d
            l <= max_level(sys) || continue
            side = 1 << d
            for conn in (Vertex(), Edge())
                inface, _, _ = classify_roots(sys, base, l, conn)
                isempty(inface) && continue
                for c in spread(inface, 3)
                    it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
                    @test Base.IteratorSize(typeof(it)) isa Base.HasLength
                    @test length(it) == (conn isa Vertex ? 4side + 4 : 4side)
                end
            end
        end
    end

    # A 64x64 block, whose band is 260 cells, element for element against the
    # oracle on both connectivities. The shallow cases above all have
    # `side <= 8`, where a wrong `_restore_code` on the way back up the face
    # descent can still land on the right cell by accident; at nine levels of
    # descent it cannot. Three roots from base 3 rather than one from base 2,
    # for `BAND_BASES`' reason: a base-2 block is symmetric under the very
    # transforms a wrong descent applies.
    @testset "the band walk at 64x64" begin
        sys = S2System()
        l = 3 + 6
        inface, _, _ = classify_roots(sys, 3, l, Vertex())
        @test !isempty(inface)
        for c in spread(inface, 3), conn in (Vertex(), Edge())
            it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
            @test it.engine isa DGG.Fallbacks.SquareBandEngine
            @test length(it) == (conn isa Vertex ? 260 : 256)
            @test collect(it) == forced_geometry_halo(sys, c, l, conn)
            # The contract bundle at nine levels of descent too: sortedness and
            # both adjacency directions are where a deep walk can go wrong
            # differently from a shallow one.
            check_halo_case(sys, c, l, conn)
        end
    end

    @testset "Edge drops exactly the four diagonal corners" begin
        sys = S2System()
        inface, _, _ = classify_roots(sys, 2, 4, Vertex())
        # Four blocks on four different faces, not the first four of face 0:
        # `_band_corner` works in lattice coordinates and the four in-face
        # blocks of one S2 face are exactly the set a wrong orientation seed
        # maps to itself.
        for c in spread(inface, 4)
            hv = collect(SubtreeHaloIterator(sys, c, 4; connectivity = Vertex()))
            he = collect(SubtreeHaloIterator(sys, c, 4; connectivity = Edge()))
            @test length(setdiff(Set(hv), Set(he))) == 4
            @test issubset(Set(he), Set(hv))
        end
    end

    # -----------------------------------------------------------------------
    # The seam walk — the same engine where the block touches a face edge
    # -----------------------------------------------------------------------

    seam_roots(sys, base::Int, seam) = unique(vcat(spread(seam, 6),
        filter(in(Set(seam)), irregular_cells(levelgrid(sys, base), 4))))

    function band_candidate_count(e)
        n = 0
        for _ in DGG.Fallbacks.SquareBandEngine(e.curve, DGG.Fallbacks.NoCheck(),
                e.level, e.faceside, e.homeface, e.x0, e.y0, e.side, true, e.rects)
            n += 1
        end
        return n
    end

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
                h = collect(it)
                @test h == forced_geometry_halo(sys, c, l, conn)
                @test band_candidate_count(it.engine) ==
                      (conn isa Vertex ? length(h) : length(subtree_halo(sys, c, l)))
                check_halo_case(sys, c, l, conn)
            end
        end
    end

    # And one arm against the `O(ncells)` brute force. Bases 0 and 1 at depth 2,
    # where the target grids are at most 3072 cells on all three systems.
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
    # sides is thirty-two cells long and only two were ever probed. If the
    # monotonicity argument were wrong anywhere, the gap would be widest here.
    @testset "the seam walk at depth five" begin
        for sys in SQUARE_SYSTEMS, conn in (Vertex(), Edge())
            c = law_root(sys)
            it = SubtreeHaloIterator(sys, c, 5; connectivity = conn)
            @test it.engine isa DGG.Fallbacks.SquareBandEngine
            @test collect(it) == forced_geometry_halo(sys, c, 5, conn)
            # And still tight thirty-two cells along a flush side, which is
            # where a bound taken lazily would have the most room to be slack.
            @test band_candidate_count(it.engine) == length(subtree_halo(sys, c, 5))
            check_halo_case(sys, c, 5, conn)
        end
    end

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
                # And tight at level 30, where the rectangles' `Int32` bounds
                # are one level from overflowing: a bound derived a level too
                # late would be wide here and nowhere else.
                @test band_candidate_count(it.engine) == length(subtree_halo(sys, c, l))
                check_halo_case(sys, c, l, conn)
            end
        end
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
        # AND ONE DEEP BASE, because the guards are not all base-independent:
        # `_hex_validate` runs only from depth three, `_hex_calibrate` reads a
        # ring that is a whole base cell's at base 0 and an ordinary hexagon's
        # at base 8, and the seeded frames sit at the other parity. The
        # generation cannot be enumerated — H3's level-8 grid is 7e8 cells — so
        # the roots are the twelve pentagons BY NAME, the ring around two of
        # them, and a positional spread. Depth 4 is in every base because it is
        # the first at which a seeded arc has been through three transitions.
        base = 8
        if base + 1 <= max_level(sys)
            roots = hex_roots(sys, base, 4)
            for conn in (Vertex(), Edge()),
                    l in (base + 1):min(base + 4, max_level(sys))
                check_hex_classes(sys, roots, base, l, conn)
            end
        end
    end

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
                it.engine isa DGG.Fallbacks.HexArcHaloEngine || continue
                wide = widen_hex_arcs(it.engine)
                @test hex_candidate_count(wide) > hex_candidate_count(it.engine)
                @test collect(SubtreeHaloIterator(sys, c, l, conn, wide)) ==
                      collect(it)
            end
        end
    end

    # The contract bundle on a smaller spread, because it builds two `Set`s and
    # a `subtree_border` per case. Pentagons first.
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
    # pentagon, and a seeded arc has been through three transitions rather than
    # one. The oracle costs about 50 ms a call, so this is two roots — a
    # pentagon and a hexagon — per system rather than a sweep.
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

    # One arm against the `O(ncells)` brute force: level-0 roots at depths 1 and
    # 2, where the target grids are at most 5882 cells on H3 and 492 on IGeo7.
    @testset "$(nameof(typeof(sys))): the directed walk against the law" for sys in
            HEX_SYSTEMS
        for c in spread(hex_roots(sys, 0, 4), 6), d in 1:2, conn in (Vertex(), Edge())
            @test collect(SubtreeHaloIterator(sys, c, d; connectivity = conn)) ==
                  law_halo(sys, c, d; connectivity = conn)
        end
    end

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

    # The count contract in the negative, as for the seam walk: the census above
    # is evidence, not an API promise, so neither hex engine declares a length.
    @testset "the directed walk declares no length" begin
        for sys in HEX_SYSTEMS, d in 1:2
            c = cellindex(levelgrid(sys, 1), 1)
            it = SubtreeHaloIterator(sys, c, 1 + d)
            @test Base.IteratorSize(typeof(it)) isa Base.SizeUnknown
            @test_throws MethodError length(it)
        end
    end

    @testset "a prefix of a deep halo costs O(depth), not O(halo)" begin
        prefix10(it) = take_n(it, 10)
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
    # A5 — the one system with no specialization, and the assertion that says so
    # -----------------------------------------------------------------------

    # A5 lacks sorted subtree ranges and a boundary automaton, so it uses the
    # scan engine. Targets stop at level 2 because the geometry oracle performs
    # a pruned boundary walk for every scanned cell.
    @testset "A5 stays on the linear scan" begin
        sys = A5System()
        for base in (0, 1), conn in (Vertex(), Edge())
            grid = levelgrid(sys, base)
            for c in sample_cells(grid, 2), l in base:2
                it = SubtreeHaloIterator(sys, c, l; connectivity = conn)
                # Depth zero is the native one-ring on every system, A5 included;
                # everything deeper is the scan.
                @test it.engine isa (l == level(c) ? DGG.Fallbacks.RingHaloEngine :
                                     DGG.Fallbacks.ScanHaloEngine)
                @test collect(it) == forced_geometry_halo(sys, c, l, conn)
                check_halo_case(sys, c, l, conn)
            end
        end
    end

    # -----------------------------------------------------------------------
    # The generic walk, still oracled where it is still the walk
    # -----------------------------------------------------------------------

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
    # Subset halos — `halo` on the three containers
    # -----------------------------------------------------------------------

    function holed_halo_oracle(sys, c, l, removed, conn, whole)
        grid = levelgrid(sys, l)
        lc = level(c)
        gone = Set(removed)
        held(x) = ancestor(sys, x, lc) == c && !(x in gone)
        out = eltype(whole)[]
        for x in vcat(whole, collect(removed))
            any(held, neighbors(grid, x, 1; connectivity = conn)) && push!(out, x)
        end
        sort!(out; by = x -> cellposition(grid, x))
        return out
    end

    @testset "$(nameof(typeof(sys))): halo on subsets" for sys in systems()
        l = min(2, max_level(sys))
        c = law_root(sys)
        pg = PartialGrid(sys, c, l)
        cv = CellVector(pg)
        expected = subtree_halo(sys, c, l)
        @test collect(halo(pg)) == expected
        @test collect(halo(cv)) == expected
        @test collect(halo(CellLookup(cv))) == expected
        # And against the geometry oracle, on both connectivities. `cv` is a
        # `CellVector`, so this reaches `SubsetHaloIterator` on every system
        # including the five whose rooted grid delegates.
        geom = Dict(conn => forced_geometry_halo(sys, c, l, conn)
                    for conn in (Vertex(), Edge()))
        for conn in (Vertex(), Edge())
            @test collect(halo(cv; connectivity = conn)) == geom[conn]
        end

        @test halo(pg) isa (DGG.has_sorted_subtrees(sys) ? SubtreeHaloIterator :
                            DGG.Fallbacks.SubsetHaloIterator)
        @test halo(cv) isa DGG.Fallbacks.SubsetHaloIterator

        # The same cells with the root forgotten must give the same answer.
        loose = PartialGrid(sys, l, collect(pg.ids))
        @test halo(loose) isa DGG.Fallbacks.SubsetHaloIterator
        @test collect(halo(loose)) == expected

        @test collect(halo(pg; connectivity = Edge())) == geom[Edge()]
        @test all(x -> cellposition(pg, x) === nothing, collect(halo(pg)))
        @test all(x -> cellposition(loose, x) === nothing, collect(halo(loose)))

        interior = collect(DGG.subtree_interior(sys, c, l))
        # one interior cell, the whole interior, a border patch. A shallow
        # subtree can be all rim — ISEA3H at this level is — so the single-cell
        # slice is bounded the same way the border patch already is; the empty
        # entries are skipped below.
        holes = (interior[1:min(1, end)], interior,
            subtree_border(sys, c, l)[1:min(3, end)])
        for removed in holes
            isempty(removed) && continue
            ids = filter(!in(Set(removed)), collect(pg.ids))
            holed = PartialGrid(sys, l, ids)
            for conn in (Vertex(), Edge())
                want = holed_halo_oracle(sys, c, l, removed, conn, geom[conn])
                hh = collect(halo(holed; connectivity = conn))
                @test hh == want
                @test hh == collect(halo(CellVector(holed); connectivity = conn))
                @test hh == collect(halo(CellLookup(CellVector(holed));
                    connectivity = conn))
                @test allunique(hh)
                @test all(x -> cellposition(holed, x) === nothing, hh)
                rooted = PartialGrid(sys, l, ids; root = c)
                @test halo(rooted) isa DGG.Fallbacks.SubsetHaloIterator
                @test collect(halo(rooted; connectivity = conn)) == hh
            end
        end
    end

    # The prune probes `coarse_probe_rings(sys)` rings, not always one, because
    # on a non-congruent refinement a descendant can leave its ancestor's
    # footprint. This is that trait's law: at the declared ring count nothing
    # escapes. Counting the escapes rather than asserting a bare `false` is what
    # makes the failure legible — and one below the declared count the counts
    # are nonzero, so this kills a trait lowered under the geometry.
    escape_count(sys, l, lc, k) = begin
        grid = levelgrid(sys, l)
        coarse = levelgrid(sys, lc)
        out = 0
        for p in 1:ncells(grid)
            x = cellindex(grid, p)
            a = ancestor(sys, x, lc)
            ring = neighbors(coarse, a, k; connectivity = Vertex())
            for y in neighbors(grid, x, 1; connectivity = Vertex())
                b = ancestor(sys, y, lc)
                (b == a || any(==(b), ring)) || (out += 1)
            end
        end
        out
    end

    @testset "the coarse-containment law the subset prune rests on" begin
        for sys in systems()
            k = DGG.coarse_probe_rings(sys)
            escaped = Tuple{Int,Int}[]
            for l in 1:6
                l <= max_level(sys) || continue
                ncells(levelgrid(sys, l)) <= 300_000 || continue
                out = escape_count(sys, l, l - 1, k)
                out == 0 || push!(escaped, (l, out))
            end
            @test escaped == Tuple{Int,Int}[]
        end
    end

    # The trait is a claim about geometry, so it has to be wrong to raise as
    # well as to lower. ISEA3H is the one registered system that needs more than
    # one ring; every other system is pinned at one so that a gratuitous raise —
    # which costs query time, since the prune then retires less — fails here.
    @testset "two rings on ISEA3H, and one everywhere else" begin
        for sys in systems()
            @test DGG.coarse_probe_rings(sys) == (sys isa DGG.ISEA3HSystem ? 2 : 1)
        end
        # The escapes ISEA3H would suffer at one ring, and the aperture-4
        # counterpart that has none. Only ODD coarse levels escape: those are
        # the even-depth ancestors, the ones the alternating-digit maximiser
        # drives to the full covering factor 2 rather than to `sqrt(3)`.
        s3 = DGG.ISEA3HSystem()
        @test [escape_count(s3, l, l - 1, 1) for l in 2:6] == [20, 0, 180, 0, 1620]
        @test all(escape_count(s3, l, l - 1, 2) == 0 for l in 2:6)
        s4 = DGG.ISEA4HSystem()
        @test all(escape_count(s4, l, lc, 1) == 0 for l in 1:5, lc in 0:4 if lc < l)
    end

    # A rooted subtree with one interior cell removed: the smallest departure
    # from a subtree there is, and the one that forces the subset walk — a
    # complete subtree would delegate and measure the wrong engine.
    function holed_subtree(sys, c, l)
        r = DGG.descendant_range(sys, c, l)
        grid = levelgrid(sys, l)
        ids = [cellindex(grid, p) for p in r]
        deleteat!(ids, length(ids) ÷ 2)
        return PartialGrid(sys, l, ids; root = c)
    end

    @testset "the subset walk's cost follows the halo, not the subset" begin
        for sys in systems()
            DGG.has_sorted_subtrees(sys) || continue
            c = law_root(sys)
            # Choose depths with comparable subtree sizes across apertures.
            depths = any(h -> sys isa typeof(h), HEX_SYSTEMS) ? (3, 6) : (4, 9)
            all(l -> l <= max_level(sys), depths) || continue
            calls = Int[]; members = Int[]; halos = Int[]
            for l in depths
                cv = CellVector(holed_subtree(sys, c, l))
                counted = CountingSubset(cv)
                h = counting_halo(sys, counted, levelgrid(sys, l), l, Vertex())
                @test h == collect(halo(cv))
                push!(calls, counted.calls)
                push!(members, length(cv))
                push!(halos, length(h))
            end
            # Member count grows much faster than halo size.
            @test members[2] / members[1] >= 4 * (halos[2] / halos[1])
            # Query growth follows halo size. This is the law's actual content,
            # and it holds on every system including the two-ring one.
            @test calls[2] <= 1.6 * (halos[2] / halos[1]) * calls[1]
            # The walk inspects only a fraction of members — an absolute bound,
            # and one that a wider coarse probe cannot meet. ISEA3H at these
            # depths: the one-ring probe makes 8987 calls, 0.46 of its members,
            # and returns 596 of the 649 halo cells — it is cheap because it is
            # wrong. Two rings makes it right at 43603 calls, 2.2 of its member
            # count. So the fraction is pinned per probe width rather than
            # dropped, and a regression to a third ring still fails here.
            if DGG.coarse_probe_rings(sys) == 1
                @test calls[2] <= 0.6 * members[2]
            else
                @test calls[2] <= 2.5 * members[2]
            end
        end
    end

    @testset "a mixed-level set has no halo" begin
        for sys in systems()
            c = cellindex(levelgrid(sys, 1), 1)
            set = DGG.MultiOrderCellSet(sys, [c], [1], trues(1), level(c))
            @test_throws MethodError halo(set)
        end
    end


    engine_tag(e) = e isa DGG.Fallbacks.SquareBandEngine ?
        (e.check isa DGG.Fallbacks.NoCheck ? :SquareBandNoCheck :
         :SquareBandNativeCheck) : nameof(typeof(e))

    ALL_ENGINE_TAGS = Set((:RingHaloEngine, :OutsideWalkEngine, :ScanHaloEngine,
        :SquareBandNoCheck, :SquareBandNativeCheck, :HexChildHaloEngine,
        :HexArcHaloEngine))

    inface_root(sys, base::Int, l::Int) =
        last(spread(first(classify_roots(sys, base, l, Vertex())), 5))

    @testset "the walk is resumable, not restarted" begin
        seen = Set{Symbol}()
        wrappers = Set{Symbol}()
        function check_prefix(it)
            push!(seen, engine_tag(it.engine))
            push!(wrappers, nameof(typeof(it)))
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
            # that same call is the scan, which covers both fallbacks without
            # naming either system.
            check_prefix(SubtreeHaloIterator(sys, c, l, Vertex(),
                DGG.Fallbacks.generic_halo_engine(sys, c, l, Vertex())))
        end
        # Depth one on the aperture-7 systems is the automaton-free child walk,
        # which the `l = level(c) + 2` cases above never reach.
        for sys in HEX_SYSTEMS
            check_prefix(SubtreeHaloIterator(sys, cellindex(levelgrid(sys, 1), 1), 2))
        end
        for sys in SQUARE_SYSTEMS
            inface, seam, _ = classify_roots(sys, 2, 4, Vertex())
            check_prefix(SubtreeHaloIterator(sys, last(spread(inface, 5)), 4))
            check_prefix(SubtreeHaloIterator(sys, last(spread(seam, 5)), 4))
        end
        # And the subset wrapper on both containers: `SubsetHaloIterator`
        # forwards the whole protocol itself, so resumability is a property of
        # the wrapper as much as of the engine inside it.
        for sys in systems()
            l = min(2, max_level(sys))
            c = law_root(sys)
            loose = PartialGrid(sys, l, collect(PartialGrid(sys, c, l).ids))
            check_prefix(halo(loose))
            check_prefix(halo(CellVector(loose)))
        end
        @test seen == ALL_ENGINE_TAGS
        @test wrappers == Set((:SubtreeHaloIterator, :SubsetHaloIterator))
    end

    # -----------------------------------------------------------------------
    # Type stability in `eltype`, on every engine and every system
    # -----------------------------------------------------------------------

    @testset "eltype is the system's cell index type, on every engine" begin
        seen = Set{Symbol}()
        wrappers = Set{Symbol}()
        function check_eltype(sys, it)
            push!(seen, engine_tag(it.engine))
            push!(wrappers, nameof(typeof(it)))
            C = DGG.cellindextype(sys)
            @test eltype(typeof(it)) === C
            @test isconcretetype(eltype(typeof(it)))
            @test collect(it) isa Vector{C}
        end
        for sys in systems()
            mx = max_level(sys)
            c0 = law_root(sys)
            check_eltype(sys, SubtreeHaloIterator(sys, c0, 0))       # the one-ring
            for l in 1:min(2, mx)
                check_eltype(sys, SubtreeHaloIterator(sys, c0, l))   # what it ships
            end
            l = min(2, mx)
            check_eltype(sys, generic_iterator(sys, c0, l))          # walk, or scan
            pg = PartialGrid(sys, c0, l)
            loose = PartialGrid(sys, l, collect(pg.ids))
            check_eltype(sys, halo(pg))
            check_eltype(sys, halo(loose))
            check_eltype(sys, halo(CellVector(loose)))
            check_eltype(sys, halo(CellLookup(CellVector(loose))))
        end
        # The counted square emit rule, which no level-0 or level-1 root can
        # reach: those blocks are flush with their face edge on all three
        # systems.
        for sys in SQUARE_SYSTEMS
            check_eltype(sys, SubtreeHaloIterator(sys, inface_root(sys, 2, 4), 4))
        end
        @test seen == ALL_ENGINE_TAGS
        @test wrappers == Set((:SubtreeHaloIterator, :SubsetHaloIterator))
    end

    @testset "subtree_halo's return type is inferred, not just correct" begin
        for sys in systems()
            c = law_root(sys)
            T = Tuple{typeof(sys),typeof(c),Int}
            @test only(Base.return_types(subtree_halo, T)) ===
                Vector{DGG.cellindextype(sys)}
        end
    end

    # -----------------------------------------------------------------------
    # The count contract, for ALL SEVEN ENGINES AT ONCE
    # -----------------------------------------------------------------------

    @testset "length is truthful where it exists and absent where it does not" begin
        counted = Set{Symbol}()
        refusing = Set{Symbol}()
        function check_count(it)
            tag = engine_tag(it.engine)
            h = collect(it)
            @test !isempty(h)
            if Base.IteratorSize(typeof(it)) isa Base.HasLength
                push!(counted, tag)
            else
                @test Base.IteratorSize(typeof(it)) isa Base.SizeUnknown
                push!(refusing, tag)
                @test_throws MethodError length(it)
            end
        end
        for sys in systems()
            mx = max_level(sys)
            c0 = law_root(sys)
            check_count(SubtreeHaloIterator(sys, c0, 0))
            for l in 1:min(2, mx)
                check_count(SubtreeHaloIterator(sys, c0, l))
            end
            l = min(2, mx)
            check_count(generic_iterator(sys, c0, l))
            loose = PartialGrid(sys, l, collect(PartialGrid(sys, c0, l).ids))
            check_count(halo(loose))
            check_count(halo(CellVector(loose)))
        end
        for sys in SQUARE_SYSTEMS
            for d in 1:3
                l = 2 + d
                l <= max_level(sys) || continue
                check_count(SubtreeHaloIterator(sys, inface_root(sys, 2, l), l))
            end
        end
        @test counted == Set((:RingHaloEngine, :SquareBandNoCheck))
        @test refusing == Set((:OutsideWalkEngine, :ScanHaloEngine,
            :SquareBandNativeCheck, :HexChildHaloEngine, :HexArcHaloEngine))
    end


    @testset "halo_positions is the ascending position stream, on every engine" begin
        seen = Set{Symbol}()
        wrappers = Set{Symbol}()
        function check_positions(grid, it)
            push!(seen, engine_tag(it.engine))
            push!(wrappers, nameof(typeof(it)))
            ids = collect(it)
            hp = halo_positions(it)
            ps = collect(hp)
            @test !isempty(ps)
            # The stream IS the conversion of the id stream, element for element
            # — including on the subset walks, where the grid a position means
            # is the COMPLETE level and not the subset.
            @test ps == [cellposition(grid, x) for x in ids]
            # The contract.
            @test issorted(ps)
            @test allunique(ps)
            # Answered in the type domain, for `collect`'s sake: an iterator
            # that only knew `Int` at run time would hand back a `Vector{Any}`
            # with nothing red anywhere.
            @test eltype(typeof(hp)) === Int
            @test ps isa Vector{Int}
            # Reading positions instead of ids changes nothing about counting.
            @test Base.IteratorSize(typeof(hp)) === Base.IteratorSize(typeof(it))
            if Base.IteratorSize(typeof(hp)) isa Base.HasLength
                @test length(hp) == length(ps)
            else
                @test_throws MethodError length(hp)
            end
            # Resumable, not restarted — the wrapper threads the walk's own
            # state, so a second pass must reproduce the first.
            @test collect(Iterators.take(hp, 4)) == ps[1:min(4, length(ps))]
            @test collect(hp) == ps
        end
        for sys in systems()
            mx = max_level(sys)
            c0 = law_root(sys)
            for l in 0:min(2, mx)
                check_positions(levelgrid(sys, l), SubtreeHaloIterator(sys, c0, l))
            end
            l = min(2, mx)
            check_positions(levelgrid(sys, l), generic_iterator(sys, c0, l))
            loose = PartialGrid(sys, l, collect(PartialGrid(sys, c0, l).ids))
            check_positions(levelgrid(sys, l), halo(loose))
            check_positions(levelgrid(sys, l), halo(CellVector(loose)))
        end
        for sys in SQUARE_SYSTEMS
            check_positions(levelgrid(sys, 4),
                SubtreeHaloIterator(sys, inface_root(sys, 2, 4), 4))
        end
        @test seen == ALL_ENGINE_TAGS
        @test wrappers == Set((:SubtreeHaloIterator, :SubsetHaloIterator))
        # The three-argument form is the two-argument one, and reaches the same
        # walk the id verb does.
        for sys in systems()
            l = min(2, max_level(sys))
            c0 = law_root(sys)
            @test collect(halo_positions(sys, c0, l)) ==
                  [cellposition(levelgrid(sys, l), x)
                   for x in subtree_halo(sys, c0, l)]
            @test collect(halo_positions(sys, c0, l; connectivity = Edge())) ==
                  [cellposition(levelgrid(sys, l), x)
                   for x in subtree_halo(sys, c0, l; connectivity = Edge())]
        end
    end

    # -----------------------------------------------------------------------
    # `halo_sizehint`: bounding, approximate, and never a count
    # -----------------------------------------------------------------------

    @testset "halo_sizehint bounds the walk, and stays out of the count" begin
        # --- the aperture-4 band, over whole generations -------------------
        excess = Dict{Tuple{Symbol,Any},Set{Int}}()
        unbounded = 0
        for sys in SQUARE_SYSTEMS, conn in (Vertex(), Edge())
            s = Set{Int}()
            for base in 0:2, d in 1:3
                l = base + d
                l <= max_level(sys) || continue
                g = levelgrid(sys, base)
                for i in 1:ncells(g)
                    it = SubtreeHaloIterator(sys, cellindex(g, i), l;
                        connectivity = conn)
                    n = length(collect(it))
                    h = halo_sizehint(it)
                    unbounded += (h === nothing || n > h)
                    push!(s, n - 4 * 2^d)
                end
            end
            excess[(nameof(typeof(sys)), conn)] = s
        end
        @test unbounded == 0
        @test excess[(:HEALPixSystem, Vertex())] == Set((2, 3, 4))
        @test excess[(:S2System, Vertex())] == Set((0, 3, 4))
        @test excess[(:ISEA4RSystem, Vertex())] == Set((2, 3, 4, 5))
        for sys in SQUARE_SYSTEMS
            @test excess[(nameof(typeof(sys)), Edge())] == Set((0,))
        end
        # The fifth cell is a fact about the icosahedral vertex and not about
        # the depth: it is still exactly five at side 4096, where a fitted
        # constant that drifted with the perimeter would be thousands off.
        let sys = ISEA4RSystem(), g = levelgrid(sys, 2)
            for d in (6, 12)
                it = SubtreeHaloIterator(sys, cellindex(g, 11), 2 + d)
                @test length(collect(it)) == 4 * 2^d + 5
                @test halo_sizehint(it) == 4 * 2^d + 8
            end
        end

        # --- the aperture-7 census -----------------------------------------
        for sys in HEX_SYSTEMS, base in 0:1, d in 1:4
            l = base + d
            l <= max_level(sys) || continue
            g = levelgrid(sys, base)
            step = max(1, ncells(g) ÷ 8)
            for i in 1:step:ncells(g)
                it = SubtreeHaloIterator(sys, cellindex(g, i), l)
                n = length(collect(it))
                @test n in (3^(d + 1) + 3, (5 * (3^d + 1)) ÷ 2)
                @test halo_sizehint(it) == 3^(d + 1) + 3
                @test n <= halo_sizehint(it)
            end
        end

        # --- where nothing is known ----------------------------------------
        for sys in systems()
            l = min(2, max_level(sys))
            c0 = law_root(sys)
            @test halo_sizehint(generic_iterator(sys, c0, l)) === nothing
            loose = PartialGrid(sys, l, collect(PartialGrid(sys, c0, l).ids))
            @test halo_sizehint(halo(loose)) === nothing
            @test halo_sizehint(halo(CellVector(loose))) === nothing
        end
        let sys = A5System(), g = levelgrid(sys, 1)
            c = cellindex(g, 1)
            @test halo_sizehint(SubtreeHaloIterator(sys, c, 2)) === nothing
            @test [length(subtree_halo(sys, c, 1 + d)) for d in 1:3] == [12, 22, 42]
        end

        # --- and it is not a count -----------------------------------------
        # Each of the four bounded engines, asked both questions. A hint that
        # had become a `length` would answer the second here instead of
        # throwing, and `SizeUnknown` would have quietly turned into
        # `HasLength` on the two that refuse.
        hinted = Set{Symbol}()
        function check_not_a_count(it)
            push!(hinted, engine_tag(it.engine))
            @test halo_sizehint(it) isa Int
            @test halo_sizehint(it) >= length(collect(it))
            @test halo_sizehint(halo_positions(it)) == halo_sizehint(it)
            if Base.IteratorSize(typeof(it)) isa Base.SizeUnknown
                @test_throws MethodError length(it)
            else
                @test length(it) == length(collect(it))
            end
        end
        for sys in SQUARE_SYSTEMS
            c0 = law_root(sys)
            check_not_a_count(SubtreeHaloIterator(sys, c0, 0))     # the one-ring
            check_not_a_count(SubtreeHaloIterator(sys, c0, 3))     # the seam band
            check_not_a_count(SubtreeHaloIterator(sys, inface_root(sys, 2, 4), 4))
        end
        for sys in HEX_SYSTEMS
            c1 = cellindex(levelgrid(sys, 1), 1)
            check_not_a_count(SubtreeHaloIterator(sys, c1, 2))
            check_not_a_count(SubtreeHaloIterator(sys, c1, 4))
        end
        @test hinted == Set((:RingHaloEngine, :SquareBandNoCheck,
            :SquareBandNativeCheck, :HexChildHaloEngine, :HexArcHaloEngine))
    end

    # -----------------------------------------------------------------------
    # "Consumable incrementally or in caller-selected batches"
    # -----------------------------------------------------------------------

    # The prefix law above says a partial walk is a prefix. This says the two
    # standard ways a caller cuts a lazy stream into pieces — `Iterators.take`
    # and `Iterators.partition` — put the pieces back together into exactly the
    # walk, at every chunk size including ones that do not divide it. Not free
    # given `SizeUnknown()`: `partition` takes a different route for a sized
    # iterator, and an engine whose `iterate` mutated shared state would
    # reassemble into something shorter than the collect.
    @testset "consumable incrementally and in caller-chosen batches" begin
        function check_batches(it)
            full = collect(it)
            @test length(full) >= 6
            @test collect(Iterators.take(it, 3)) == full[1:3]
            @test collect(Iterators.take(it, length(full) + 5)) == full
            @test isempty(collect(Iterators.take(it, 0)))
            for k in (1, 2, 5)
                parts = collect.(Iterators.partition(it, k))
                @test reduce(vcat, parts) == full
                @test all(p -> 1 <= length(p) <= k, parts)
            end
        end
        for sys in systems()
            l = min(2, max_level(sys))
            c0 = law_root(sys)
            check_batches(SubtreeHaloIterator(sys, c0, l))
            check_batches(generic_iterator(sys, c0, l))
            loose = PartialGrid(sys, l, collect(PartialGrid(sys, c0, l).ids))
            check_batches(halo(loose))
        end
        for sys in SQUARE_SYSTEMS
            check_batches(SubtreeHaloIterator(sys, inface_root(sys, 2, 4), 4))
        end
    end

    # -----------------------------------------------------------------------
    # The one awkward cell that DEGREE does not find
    # -----------------------------------------------------------------------

    @testset "the two HEALPix poles, pinned by location rather than by degree" begin
        sys = HEALPixSystem()
        for base in (1, 2, 3), lat in (89.999, -89.999)
            grid = levelgrid(sys, base)
            p = DGG.cellat(grid, 0.0, lat)
            @test p !== nothing
            for d in 1:2, conn in (Vertex(), Edge())
                l = base + d
                l <= max_level(sys) || continue
                @test collect(SubtreeHaloIterator(sys, p, l; connectivity = conn)) ==
                      forced_geometry_halo(sys, p, l, conn)
                check_halo_case(sys, p, l, conn)
            end
        end
        grid = levelgrid(sys, 1)
        for lat in (89.999, -89.999), conn in (Vertex(), Edge())
            p = DGG.cellat(grid, 0.0, lat)
            @test collect(SubtreeHaloIterator(sys, p, 3; connectivity = conn)) ==
                  law_halo(sys, p, 3; connectivity = conn)
        end
    end



    DEPTH_FLAT_SYSTEMS = (HEALPixSystem(), S2System(), ISEA4RSystem(),
        H3System(), IGeo7System())

    @testset "construction does not allocate in proportion to the halo" begin
        for sys in systems()
            c = law_root(sys)
            # A5's targets stop at 3: its `subtree_halo` at level 4 is a
            # 3840-cell scan whose per-cell cost is a `Set`-allocating
            # `neighbors`, and the law here needs only two comparable points.
            depths = filter(l -> l <= max_level(sys),
                sys isa DGG.A5System ? (1, 2, 3) : law_depths(sys, (3, 5, 7)))
            ship = [construct_bytes(sys, c, l) for l in depths]
            gen = [generic_construct_bytes(sys, c, l) for l in depths]
            sizes = [length(subtree_halo(sys, c, l)) for l in depths]
            @test last(sizes) >= 2 * first(sizes)
            @test maximum(ship) - minimum(ship) <= 64
            @test maximum(gen) - minimum(gen) <= 64
        end
    end

    @testset "the specialized prefix costs the same at every depth" begin
        for sys in DEPTH_FLAT_SYSTEMS
            c = law_root(sys)
            depths = filter(l -> l <= max_level(sys), (3, 5, 7))
            allocs = [lazy_bytes(sys, c, l, 4) for l in depths]
            sizes = [length(subtree_halo(sys, c, l)) for l in depths]
            # Halo cardinality increases substantially across these levels.
            @test last(sizes) >= 8 * first(sizes)
            # Four-cell prefix allocation remains level-independent.
            @test maximum(allocs) - minimum(allocs) <= 64
            @test all(>(0), allocs)
        end
        for sys in SQUARE_SYSTEMS
            c = inface_root(sys, 3, 4)
            depths = filter(l -> l <= max_level(sys), (4, 6, 9))
            allocs = [lazy_bytes(sys, c, l, 4) for l in depths]
            sizes = [length(subtree_halo(sys, c, l)) for l in depths]
            @test last(sizes) >= 8 * first(sizes)
            @test maximum(allocs) - minimum(allocs) <= 64
        end
    end

    @testset "a short prefix is a small, non-growing fraction of the collect" begin
        for sys in systems()
            c = law_root(sys)
            depths = sys isa DGG.A5System ? (1, 4) : law_depths(sys, (3, 7))
            all(l -> l <= max_level(sys), depths) || continue
            fracs = [lazy_bytes(sys, c, l, 4) / eager_bytes(sys, c, l)
                     for l in depths]
            sizes = [length(subtree_halo(sys, c, l)) for l in depths]
            # The complete halo grows across these target levels.
            @test last(sizes) >= 5 * first(sizes)
            # Prefix allocation remains bounded.
            @test last(fracs) <= first(fracs)
            @test last(fracs) < 0.15
        end
    end

    @testset "the generic walk obeys the same fraction law" begin
        for sys in systems()
            sys isa DGG.A5System && continue
            c = law_root(sys)
            depths = law_depths(sys, (3, 6))
            all(l -> l <= max_level(sys), depths) || continue
            fracs = Float64[]
            sizes = Int[]
            for l in depths
                h = generic_collect(sys, c, l)             # warm up, and count
                eb = @allocated generic_collect(sys, c, l)
                generic_take(sys, c, l, 4)                 # warm up
                lb = @allocated generic_take(sys, c, l, 4)
                push!(fracs, lb / eb)
                push!(sizes, length(h))
            end
            # Collection size grows substantially between target levels.
            @test last(sizes) >= 5 * first(sizes)
            @test last(fracs) <= first(fracs)
            @test last(fracs) < 0.15
        end
    end

    # The negative control for the three allocation laws above: wrap the real
    # engine in one that materializes, and check the laws REFUSE it. Without it
    # they could be passing because nothing measurable is allocated either way.
    #
    # Three systems and not `systems()`. Eagerness is a property of the wrapper,
    # not of the system under it — an eager engine over IVEA4R fails for the
    # same reason it fails over S2 — so the only axis worth spanning is the
    # engine being wrapped: a square specialization, a hex one, and a system on
    # the generic walk. Sweeping all fifteen re-proved that at depth 7 on
    # 4.8M-cell aperture-9 subtrees, which was 170 of the 230 minutes the suite
    # took when this file first met the literature registry.
    EAGER_CONTROL_SYSTEMS = (S2System(), IGeo7System(), IVEA4RSystem())

    @testset "an eager engine with the same surface fails every allocation law" begin
        for sys in EAGER_CONTROL_SYSTEMS
            c = law_root(sys)
            depths = filter(l -> l <= max_level(sys), law_depths(sys, (3, 5, 7)))
            @test collect(fixture_iterator(sys, c, last(depths))) ==
                  subtree_halo(sys, c, last(depths))
            ctor = [fixture_construct_bytes(sys, c, l) for l in depths]
            pref = [fixture_prefix_bytes(sys, c, l, 4) for l in depths]
            @test maximum(ctor) - minimum(ctor) > 64      # arm 1 refuses it
            @test maximum(pref) - minimum(pref) > 64      # arm 2 refuses it
            @test fixture_prefix_bytes(sys, c, last(depths), 4) >=
                  0.15 * fixture_collect_bytes(sys, c, last(depths))   # arm 3
        end
    end

    # The absolute query count of one complete subset walk, per system. A ratio
    # cannot see a walk that gets uniformly more expensive, which is what these
    # pin.
    #
    # The spread across the second group is the interesting part, and it is not
    # noise: IVEA9R and RTEA9R ask over half a million questions where rHEALPix,
    # at the same aperture and the same subset size, asks 9612. They are the
    # systems with no specialized halo engine AND no sorted subtrees, so the
    # walk cannot prune by descendant range and falls back to asking about
    # cells one at a time. That is a standing performance gap, recorded here so
    # that closing it shows up as a failure rather than passing unnoticed.
    SUBSET_QUESTION_TOTALS = Dict(
        :IGeo7System => 4622, :H3System => 6005, :HEALPixSystem => 2143,
        :A5System => 28010, :S2System => 1572, :ISEA4RSystem => 1865,
        :ISEA3HSystem => 2641, :ISEA4HSystem => 1162, :ISEA4TSystem => 3296,
        :RHEALPixSystem => 9612, :AuthalicSystem => 9612, :IVEA4RSystem => 20548,
        :IVEA9RSystem => 535818, :RTEA4RSystem => 20538, :RTEA9RSystem => 535821)

    @testset "the subset walk is lazy, and its construction is O(1)" begin
        for sys in systems()
            mx = max_level(sys)
            ctor = (grid = Int[], vector = Int[])
            sizes = Int[]
            for l in unique((1, min(4, mx)))
                c = law_root(sys)
                loose = PartialGrid(sys, l, collect(PartialGrid(sys, c, l).ids))
                push!(sizes, ncells(loose))
                for (built, sub) in ((ctor.grid, loose),
                        (ctor.vector, CellVector(loose)))
                    push!(built, subset_construct_bytes(sub))
                end
            end
            # The larger subset has over ten times as many members.
            @test last(sizes) >= 10 * first(sizes)
            # Construction allocation remains independent of member count.
            @test maximum(ctor.grid) - minimum(ctor.grid) <= 64
            @test maximum(ctor.vector) - minimum(ctor.vector) <= 64
            # Check traversal work on the larger subset.
            deep = min(4, mx)
            root = law_root(sys)
            cv = CellVector(PartialGrid(sys, deep,
                collect(PartialGrid(sys, root, deep).ids)))
            complete = levelgrid(sys, deep)
            lazy = CountingSubset(cv)
            lp, lc, lh = counted_walk(
                counting_iterator(sys, lazy, complete, deep, Vertex()), lazy, 4)
            eager = CountingSubset(cv)
            ep, ec, eh = counted_walk(
                counting_iterator(sys, eager, complete, deep, Vertex(),
                    EagerHaloEngine), eager, 4)
            @test lh == eh == collect(halo(cv))
            @test lp <= 0.2 * lc         # the law
            @test ep > 0.2 * ec          # and the engine it refuses
            # Pin the absolute query count as well as its ratio.
            @test lc == SUBSET_QUESTION_TOTALS[nameof(typeof(sys))]
        end
    end

    @testset "the size hint is what preallocates the collect" begin
        function hint_bytes(sys, c, l)
            it = SubtreeHaloIterator(sys, c, l)
            h = halo_sizehint(it)
            @test h !== nothing
            DGG.collect_subtree(it, nothing)
            DGG.collect_subtree(it, h)
            grown = @allocated DGG.collect_subtree(it, nothing)
            hinted = @allocated DGG.collect_subtree(it, h)
            out = DGG.collect_subtree(it, h)
            return (grown, hinted, length(out) * sizeof(eltype(it)))
        end
        # A seam block on two of the three aperture-4 systems, and a depth-8
        # halo on both aperture-7 ones: 16387, 16389, 19686 and 16405 cells,
        # which is large enough that the growth doubling has happened many
        # times over.
        for (sys, base, i, l) in ((HEALPixSystem(), 2, 6, 14),
                                  (ISEA4RSystem(), 2, 11, 14),
                                  (H3System(), 0, 1, 8),
                                  (IGeo7System(), 0, 1, 8))
            c = cellindex(levelgrid(sys, base), i)
            grown, hinted, answer = hint_bytes(sys, c, l)
            @test 2 * hinted < grown
            @test hinted < 1.4 * answer
            @test grown > 2 * answer
        end
        let sys = S2System()
            grown, hinted, _ = hint_bytes(sys, cellindex(levelgrid(sys, 2), 6), 14)
            @test hinted < grown
        end
    end

    # -----------------------------------------------------------------------
    # The wrapper, and the guard on a lying count
    #
    # Last rather than first, because the wrapper test needs `classify_roots`.
    # -----------------------------------------------------------------------

    # AusPIX is a registered `AuthalicSystem`, so the sweeps below have to skip
    # the systems that are already wrapped: re-wrapping one throws by contract,
    # and the law here is about what the wrapper forwards, which its own parent
    # already covers.
    rewrappable(sys) = !(sys isa DGG.AuthalicSystem)

    @testset "AuthalicSystem forwards the halo walk" begin
        seen = Set{Symbol}()
        for sys in filter(rewrappable, systems())
            wrapped = DGG.AuthalicSystem(sys)
            c = law_root(sys)
            for l in level(c):min(level(c) + 2, max_level(sys))
                it = SubtreeHaloIterator(wrapped, c, l)
                push!(seen, engine_tag(it.engine))
                @test collect(it) == collect(SubtreeHaloIterator(sys, c, l))
            end
        end
        for sys in filter(rewrappable, systems())
            wrapped = DGG.AuthalicSystem(sys)
            c = law_root(sys)
            l = min(2, max_level(sys))
            it = SubtreeHaloIterator(wrapped, c, l, Vertex(),
                DGG.Fallbacks.generic_halo_engine(wrapped, c, l, Vertex()))
            push!(seen, engine_tag(it.engine))
            @test collect(it) == collect(SubtreeHaloIterator(sys, c, l, Vertex(),
                DGG.Fallbacks.generic_halo_engine(sys, c, l, Vertex())))
        end
        for sys in SQUARE_SYSTEMS
            wrapped = DGG.AuthalicSystem(sys)
            for base in BAND_BASES
                l = base + 2
                l <= max_level(sys) || continue
                inface, _, _ = classify_roots(sys, base, l, Vertex())
                isempty(inface) && continue
                c = last(spread(inface, 5))
                for conn in (Vertex(), Edge())
                    it = SubtreeHaloIterator(wrapped, c, l; connectivity = conn)
                    push!(seen, engine_tag(it.engine))
                    @test engine_tag(it.engine) === :SquareBandNoCheck
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
                push!(seen, engine_tag(it.engine))
                @test it.engine isa (d == 1 ? DGG.Fallbacks.HexChildHaloEngine :
                                     DGG.Fallbacks.HexArcHaloEngine)
                @test collect(it) ==
                      collect(SubtreeHaloIterator(sys, c, 1 + d; connectivity = conn))
            end
        end
        @test seen == ALL_ENGINE_TAGS
        for sys in filter(rewrappable, systems())
            wrapped = DGG.AuthalicSystem(sys)
            l = min(2, max_level(sys))
            c = law_root(sys)
            pg = PartialGrid(wrapped, c, l)
            loose = PartialGrid(wrapped, l, collect(pg.ids))
            expected = subtree_halo(sys, c, l)
            @test halo(pg) isa (DGG.has_sorted_subtrees(sys) ? SubtreeHaloIterator :
                                DGG.Fallbacks.SubsetHaloIterator)
            @test collect(halo(pg)) == expected
            @test collect(halo(loose)) == expected
            @test collect(halo(CellVector(pg))) == expected
            @test collect(halo(CellLookup(CellVector(pg)))) == expected
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
