# ---------------------------------------------------------------------------
# T22 — the lazy subtree halo.
#
# `NeighborCellIterator` is the third member of the T20 family, and the one that
# walks OUTSIDE the subtree: the level-`l` cells that are not descendants of `c`
# but touch one. `subtree_halo` is `collect` of it, and `halo(sub)` is the same
# question asked of a subset container. The defining law is executable:
#
#     halo == sort!(unique!([nb for r in subtree_border(sys, c, l; κ)
#                               for nb in neighbors(levelgrid(sys, l), r, 1; κ)
#                               if ancestor(sys, nb, level(c)) != c]))
#
# Written against the generic interface only, like its siblings here: no system
# module is imported, every law runs against every system in `systems()`, plus
# an `AuthalicSystem`-wrapped run — with one internal touch, the engine-type
# checks against `DGG.Fallbacks.EagerEngine`, which are how the sweeps PROVE
# the automaton was reached rather than the generic engine tested twice.
# Helpers are deliberate copies of the ones in `subtree_iterators.jl` rather
# than imports — the two files must be able to disagree.
#
# The laws, in order:
#
#   * DEFINING LAW — the iterator, the eager verb, and the brute-force oracle
#     agree, all systems × both connectivities × depths 0-2 and one deep.
#   * CONSEQUENCES — halo is disjoint from the subtree, deduplicated, in the
#     complete level's position order; every border cell has a neighbour in the
#     halo and every halo cell one in the border; `Edge()`'s halo is inside
#     `Vertex()`'s, strictly so somewhere concrete.
#   * THE AUTOMATON — the exterior-perimeter walk on the three aligned-block
#     systems, which the bases-0-and-1 sweep cannot reach: a level-0 block is
#     flush with its face on every side and a level-1 block on two of its four
#     (one per axis), and any flush block takes the generic engine (its halo
#     crosses a seam). The level-2 sweep is where the lazy walk actually runs,
#     so it gets its own oracle pass, a deep element-for-element case, and the
#     T20-style allocation ratio.
#   * SUBSETS — `halo(sub)` agrees between a rooted `PartialGrid`, the same ids
#     with the root forgotten, and a `CellVector`; every element is out of set;
#     a hole's cells join the halo.
#   * GUARD — `collect` routes through `collect_subtree`, so a counted engine
#     that under-delivers errors instead of handing back `undef` as cell ids.
#
# IGeo7 and H3 take the materialising generic engine everywhere — an aperture-7
# subtree boundary is a fractal, not a block, so there is no perimeter to walk —
# and A5 does too, for T20's reason: no `descendant_range`, no automaton. Their
# iterators are honest `O(rim · degree)`-memory iterators, excluded from the
# allocation law exactly as A5 is excluded in `subtree_iterators.jl`: they would
# measure nothing, not fail.
# ---------------------------------------------------------------------------

module SubtreeHaloTests

using Test
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: systems, levels, max_level, levelgrid, ncells,
    cellindex, cellposition, neighbors, level, ancestor, descendants, rootcells,
    subtree_border, subtree_interior, subtree_halo, halo,
    NeighborCellIterator, PartialGrid, CellVector,
    AuthalicSystem, Vertex, Edge, Connectivity

# The defining law, verbatim. Independent of the AUTOMATON wherever an
# automaton runs — the square systems' band descent shares none of this
# arithmetic — but on the generic-engine paths (IGeo7, H3, A5, flush blocks)
# it is a re-derivation from the same T20-proven primitives the engine itself
# uses, not a second opinion. The genuinely independent checks there are the
# two coverage directions, sortedness by the level's own positions, and the
# subset laws below. `sort!` before `unique!` so the dedup is the cheap
# adjacent one and the result is already in canonical order.
function brute_force_halo(sys, c, l; connectivity = Vertex())
    grid = levelgrid(sys, l)
    lc = level(c)
    out = [nb for r in subtree_border(sys, c, l; connectivity)
              for nb in neighbors(grid, r, 1; connectivity)
              if ancestor(sys, nb, lc) != c]
    return unique!(sort!(out))
end

# A deterministic spread: no RNG, so a failure names the same cell every run.
function sample_cells(grid, n::Int)
    total = ncells(grid)
    step = max(1, total ÷ n)
    return [cellindex(grid, i) for i in 1:step:total]
end

# The cells whose one-ring is not the modal size — pentagons, face corners —
# found by degree so the sweep needs no system knowledge.
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

sweep_bases(sys) = filter(b -> b <= max_level(sys), (0, 1))

function sweep_roots(sys, base::Int)
    grid = levelgrid(sys, base)
    return unique(vcat(sample_cells(grid, 4), irregular_cells(grid, 2)))
end

function deep_depth(sys, base::Int, budget::Int = 70_000)
    d = 0
    while base + d + 1 <= max_level(sys) &&
        ncells(levelgrid(sys, base + d + 1)) ÷ ncells(levelgrid(sys, base)) <= budget
        d += 1
    end
    return d
end

# The allocation harness from `subtree_iterators.jl`, verbatim: construction is
# inside the measurement, or an engine that materialises in its constructor
# measures zero.
function take_n(it, n::Int)
    seen = 0
    for _ in it
        seen += 1
        seen >= n && break
    end
    return seen
end

build_and_take(T, sys, c, l, n::Int) = take_n(T(sys, c, l), n)
build_and_collect(T, sys, c, l) = length(collect(T(sys, c, l)))

function lazy_bytes(T, sys, c, l, n::Int)
    build_and_take(T, sys, c, l, n)
    return @allocated build_and_take(T, sys, c, l, n)
end

function eager_bytes(T, sys, c, l)
    build_and_collect(T, sys, c, l)
    return @allocated build_and_collect(T, sys, c, l)
end

# Claims three, yields one — exactly the shape `collect_subtree` exists to
# catch, wired in as an engine through the positional constructor.
struct MiscountingEngine end
Base.iterate(::MiscountingEngine) = (DGG.LevelIndex(0, 0), 1)
Base.iterate(::MiscountingEngine, ::Int) = nothing
Base.eltype(::Type{MiscountingEngine}) = DGG.LevelIndex
Base.IteratorSize(::Type{MiscountingEngine}) = Base.HasLength()
Base.length(::MiscountingEngine) = 3

# One halo, all of its per-case laws. Factored because the base sweep and the
# level-2 automaton sweep assert the same things about different roots.
function check_halo_case(sys, c, l, conn)
    it = NeighborCellIterator(sys, c, l; connectivity = conn)
    h = collect(it)
    @test h == subtree_halo(sys, c, l; connectivity = conn)
    @test h == brute_force_halo(sys, c, l; connectivity = conn)
    @test allunique(h)
    lc = level(c)
    @test all(x -> ancestor(sys, x, lc) != c, h)
    grid = levelgrid(sys, l)
    @test issorted([cellposition(grid, x) for x in h])
    # Both coverage directions: border cells reach the halo and halo cells
    # reach the border, under the same connectivity.
    border = subtree_border(sys, c, l; connectivity = conn)
    hs, bs = Set(h), Set(border)
    @test all(r -> any(in(hs), neighbors(grid, r, 1; connectivity = conn)), border)
    @test all(x -> any(in(bs), neighbors(grid, x, 1; connectivity = conn)), h)
    # Wherever the engine claims a count, the walk delivers it (the guarded
    # `collect` above would have errored otherwise; this pins `length` too).
    Base.IteratorSize(typeof(it)) isa Base.HasLength && @test length(it) == length(h)
    return h
end

@testset "subtree halos" begin

    for sys in systems(), base in sweep_bases(sys)
        sysname = string(nameof(typeof(sys)))
        name = "$sysname at level $base"
        roots = sweep_roots(sys, base)

        @testset "$name: the defining law and its consequences" begin
            for c in roots
                deep = base + deep_depth(sys, base)
                ls = c == first(roots) ? (base:base+2..., deep) : (base:base+2...,)
                for l in unique(ls)
                    l <= max_level(sys) || continue
                    hv = check_halo_case(sys, c, l, Vertex())
                    he = check_halo_case(sys, c, l, Edge())
                    # An edge-neighbour is a vertex-neighbour, so the halos nest.
                    @test issubset(Set(he), Set(hv))
                end
            end
        end

        @testset "$name: depth zero, and validation" begin
            c = first(roots)
            lc = level(c)
            grid = levelgrid(sys, lc)
            for conn in (Vertex(), Edge())
                # A cell's halo at its own level is exactly its one-ring, sorted.
                expected = sort!([nb for nb in neighbors(grid, c, 1; connectivity = conn)])
                @test collect(NeighborCellIterator(sys, c, lc; connectivity = conn)) ==
                      expected
                @test subtree_halo(sys, c, lc; connectivity = conn) == expected
            end
            if lc > 0
                @test_throws ArgumentError NeighborCellIterator(sys, c, lc - 1)
            end
            @test_throws ArgumentError NeighborCellIterator(sys, c, max_level(sys) + 1)
        end
    end

    # The exterior-perimeter automaton. Only the three aligned-block systems
    # have one, and only for blocks not flush with a face edge — which first
    # exist at level 2 — so this sweep is the automaton's oracle pass, and the
    # `lazies` counter proves the sweep reached it rather than testing the
    # generic engine twice.
    for sys in systems()
        sysname = string(nameof(typeof(sys)))
        sysname in ("HEALPixSystem", "S2System", "ISEA4RSystem") || continue
        grid2 = levelgrid(sys, 2)
        cells2 = [cellindex(grid2, i) for i in 1:ncells(grid2)]

        @testset "$sysname: the exterior-perimeter automaton, against the oracle" begin
            lazies = 0
            for c in cells2, d in 1:2, conn in (Vertex(), Edge())
                it = NeighborCellIterator(sys, c, 2 + d; connectivity = conn)
                it.engine isa DGG.Fallbacks.EagerEngine || (lazies += 1)
                check_halo_case(sys, c, 2 + d, conn)
            end
            @test lazies > 0
        end

        @testset "$sysname: lazy" begin
            c = first(x for x in cells2 if !(NeighborCellIterator(sys, x, 3).engine
                isa DGG.Fallbacks.EagerEngine))
            deep = deep_depth(sys, 2)
            @test deep >= 4                    # or the measurement is not one
            deep_l = 2 + deep
            shallow_l = 2 + deep - 3

            lazy = lazy_bytes(NeighborCellIterator, sys, c, deep_l, 4)
            eager = eager_bytes(NeighborCellIterator, sys, c, deep_l)

            # Building the walk and taking four cells off the front costs a
            # fraction of the walk they are the front of...
            @test lazy * 8 < eager
            # ...and does not grow with the perimeter, only with the depth.
            @test lazy <= lazy_bytes(NeighborCellIterator, sys, c, shallow_l, 4) + 64

            # Resumable, not restarted: the prefix is the collected prefix.
            it = NeighborCellIterator(sys, c, deep_l)
            prefix = eltype(it)[]
            for x in it
                push!(prefix, x)
                length(prefix) >= 4 && break
            end
            @test prefix == collect(it)[1:4]
        end
    end

    # The level-2 sweep above element-checks the automaton only up to side 4;
    # the deep runs in "lazy" check the count and a prefix, not the elements.
    # One deep block, element for element against the oracle, closes that gap —
    # the engine type is the proof the block is non-flush and the automaton is
    # the code under test.
    @testset "the automaton at depth, element for element" begin
        sys = only(s for s in systems() if string(nameof(typeof(s))) == "S2System")
        grid2 = levelgrid(sys, 2)
        c = first(x for x in (cellindex(grid2, i) for i in 1:ncells(grid2))
            if !(NeighborCellIterator(sys, x, 3).engine isa DGG.Fallbacks.EagerEngine))
        l = 2 + 6                              # a 64 x 64 block, halo 260 cells
        for conn in (Vertex(), Edge())
            it = NeighborCellIterator(sys, c, l; connectivity = conn)
            @test !(it.engine isa DGG.Fallbacks.EagerEngine)
            @test collect(it) == brute_force_halo(sys, c, l; connectivity = conn)
        end
    end

    # `Edge() ⊆ Vertex()` is strict somewhere concrete: a non-flush square
    # block's vertex halo has exactly the four diagonal-contact corners its
    # edge halo lacks. Pinned on the first S2 level-2 cell whose block is not
    # flush with its face — the engine type says which — so a regression names
    # a specific cell. (A flush block can lose a corner across a cube corner,
    # which is why the count is pinned here and not in the sweep above.)
    @testset "Edge halo is strictly inside Vertex halo on a square block" begin
        sys = only(s for s in systems() if string(nameof(typeof(s))) == "S2System")
        grid2 = levelgrid(sys, 2)
        c = first(x for x in (cellindex(grid2, i) for i in 1:ncells(grid2))
            if !(NeighborCellIterator(sys, x, 3).engine isa DGG.Fallbacks.EagerEngine))
        hv = subtree_halo(sys, c, 3; connectivity = Vertex())
        he = subtree_halo(sys, c, 3; connectivity = Edge())
        @test issubset(Set(he), Set(hv))
        @test length(setdiff(Set(hv), Set(he))) == 4
    end

    # `halo(sub)`: the subset verb. On a rooted whole-subtree `PartialGrid` it
    # IS the lazy iterator; on anything else it is the eager clipped scan — and
    # the two must agree exactly on the same cells, which is the law that keeps
    # the fast path honest.
    @testset "halo on subsets" begin
        for sys in systems()
            base = min(1, max_level(sys))
            l = min(base + 2, max_level(sys))
            c = first(sweep_roots(sys, base))
            want = subtree_halo(sys, c, l)

            pg = PartialGrid(sys, c, l)
            ids = descendants(sys, c, l)
            loose = PartialGrid(sys, l, ids)   # same cells, root forgotten
            cv = CellVector(sys, l, ids)

            DGG.has_sorted_subtrees(sys) && @test halo(pg) isa NeighborCellIterator
            @test collect(halo(pg)) == want
            @test collect(halo(loose)) == want
            @test collect(halo(cv)) == want
            @test collect(halo(pg; connectivity = Edge())) ==
                  subtree_halo(sys, c, l; connectivity = Edge())
            @test all(x -> cellposition(pg, x) === nothing, collect(halo(pg)))

            # A hole's cells touch the subset from outside, so they are halo.
            interior = subtree_interior(sys, c, l)
            if !isempty(interior)
                hole = interior[max(1, length(interior) ÷ 2)]
                kept = [x for x in ids if x != hole]
                holed_pg = PartialGrid(sys, l, kept)
                holed_cv = CellVector(sys, l, kept)
                @test hole in collect(halo(holed_pg))
                @test hole in collect(halo(holed_cv))
                @test collect(halo(holed_pg)) == collect(halo(holed_cv))
            end
        end
    end

    # The wrapper is a system too, and the halo of a subtree is a question
    # about the discrete hierarchy that the ellipsoid does not touch.
    @testset "AuthalicSystem forwards the halo walk" begin
        for sys in systems()
            wrapped = AuthalicSystem(sys)
            base = min(1, max_level(sys))
            c = first(sweep_roots(sys, base))
            for d in 0:2
                l = base + d
                l <= max_level(sys) || continue
                @test collect(NeighborCellIterator(wrapped, c, l)) ==
                      collect(NeighborCellIterator(sys, c, l))
                @test subtree_halo(wrapped, c, l) == subtree_halo(sys, c, l)
            end
            if base + 2 <= max_level(sys)
                @test collect(halo(PartialGrid(wrapped, c, base + 2))) ==
                      collect(halo(PartialGrid(sys, c, base + 2)))
            end
        end
    end

    # `collect` must BE the guarded path, not merely parallel to it: given
    # `HasLength`, Base's own `collect` sizes a vector from `length` and an
    # under-delivering walk would return an `undef` tail as cell ids. The lying
    # engine below is the regression test for the routing itself.
    @testset "collect is the guarded path" begin
        sys = first(systems())
        c = first(rootcells(sys))
        lying = NeighborCellIterator(sys, c, Int(level(c)), Vertex(),
            MiscountingEngine())
        @test_throws ErrorException collect(lying)
    end
end

end # module SubtreeHaloTests
