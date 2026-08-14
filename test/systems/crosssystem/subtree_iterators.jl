# ---------------------------------------------------------------------------
# T20 — the lazy subtree walkers.
#
# `EdgeCellIterator` and `InnerCellIterator` are the resumable form of
# `subtree_border` / `subtree_interior`, and since T20 they are also the only
# form: the eager verbs are `collect` of them, so a system writes ONE walk and
# both faces of it are that walk. That is what makes this file's first law worth
# stating even though it looks tautological — it pins that the two faces cannot
# drift apart, because there is nothing left for them to drift between.
#
# Written against the generic interface only, like its neighbours here: no
# system module is imported, and every law runs against every system in
# `systems()` plus an `AuthalicSystem`-wrapped one.
#
# The laws, in order:
#
#   * AGREEMENT — `collect` of each iterator is the eager verb, and both are
#     the brute-force definition. Hexagon AND pentagon / face-corner roots,
#     depths 0-3, both connectivities.
#   * GENERIC AGREEMENT — the same, against the generic engine reached by
#     `invoke`. Every system but A5 overrides the walk, so without this the
#     fallback in `src/fallbacks/subtree_iterators.jl` would be unreachable
#     code on all six; with it, each automaton is differentially tested against
#     a walk that shares none of its arithmetic. This is the oracle for S2's
#     new Hilbert fast path in particular.
#   * PARTITION — rim and interior are disjoint, ascending, and together the
#     descendants.
#   * COUNTED — `length` agrees with the walk wherever `IteratorSize` promises
#     one, and the `HasLength`/`SizeUnknown` split is where it claims to be.
#   * LAZY — the measurement the whole task is for: taking four cells off a
#     deep rim allocates a bounded amount that does NOT grow with the rim,
#     while `collect` of it grows linearly. Asserted as a ratio, not as bytes.
#     Excluded on A5 — see below.
#
# A5 is the one system with no rim automaton, and structurally so: its four
# Hilbert children cover the parent's area but not its footprint, so there is no
# child-adjacency predicate to prune on, and `has_sorted_subtrees` is false so
# there is no position range to walk instead. Its iterators are honest
# iterators over an internally materialised vector — correct, `HasLength`, and
# `O(subtree)` in memory. Every law here applies to it except LAZY, and LAZY is
# excluded because on A5 it would MEASURE NOTHING rather than fail: the subtree
# is built in the constructor, so iterating the result allocates zero and the
# ratio comes out as favourable as a real automaton's while the memory has
# already been spent. An exclusion, not a known-failure.
# ---------------------------------------------------------------------------

module SubtreeIteratorTests

using Test
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: systems, levels, max_level, levelgrid, ncells,
    cellindex, neighbors, level, descendants, rootcells, cellindextype,
    subtree_border, subtree_interior, EdgeCellIterator, InnerCellIterator,
    AuthalicSystem, Vertex, Edge, Connectivity

# The rim of `c`'s subtree at `l`, straight from the definition: a descendant
# with a neighbour that is not a descendant. Slow and obviously correct — the
# whole point is that it shares no code with the walks it checks. Deliberately
# a copy of the one in `runtests.jl` rather than an import: these two files must
# be able to disagree.
function brute_force_border(sys, c, l; connectivity = Vertex())
    lc = level(c)
    l == lc && return [c]
    grid = levelgrid(sys, l)
    inside = Set(descendants(sys, c, l))
    return [d for d in descendants(sys, c, l)
            if any(nb -> !(nb in inside), neighbors(grid, d, 1; connectivity))]
end

# A deterministic spread: no RNG, so a failure names the same cell every run.
function sample_cells(grid, n::Int)
    total = ncells(grid)
    step = max(1, total ÷ n)
    return [cellindex(grid, i) for i in 1:step:total]
end

# The cells whose one-ring is not the modal size: IGeo7's and H3's twelve
# pentagons, S2's 24 face corners, ISEA4R's icosahedral vertices. Stated by
# degree rather than by name so the sweep finds them without knowing what a
# pentagon is, and finds nothing (harmlessly) on a system that has none.
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

# The generic engine, reached past every system's override. `invoke` on the
# least-specific signature is the only way to ask a system "what would the
# fallback have said?", which is exactly the question a differential test on an
# automaton needs answered.
const GENERIC = Tuple{DGG.AbstractHierarchicalGridSystem,DGG.AbstractCellIndex,
    Int,Connectivity}

generic_rim(sys, c, l, conn) =
    collect(invoke(DGG.rim_engine, GENERIC, sys, c, Int(l), conn))
generic_interior(sys, c, l, conn) =
    collect(invoke(DGG.interior_engine, GENERIC, sys, c, Int(l), conn))

# Iterating without collecting: the shape the LAZY law measures.
function take_n(it, n::Int)
    seen = 0
    for _ in it
        seen += 1
        seen >= n && break
    end
    return seen
end

# `@allocated` behind a function barrier, and always after a warm-up call, so
# what is measured is the walk rather than its compilation.
function lazy_bytes(it, n::Int)
    take_n(it, n)
    return @allocated take_n(it, n)
end

function eager_bytes(it)
    collect(it)
    return @allocated collect(it)
end

# The base level each system is swept from, and the roots swept at it. Level 1
# everywhere: deep enough that a subtree is not the whole sphere, shallow enough
# that depth 3 stays cheap on an aperture-7 system.
function sweep_roots(sys)
    base = min(1, max_level(sys))
    grid = levelgrid(sys, base)
    roots = vcat(sample_cells(grid, 4), irregular_cells(grid, 2))
    return base, unique(roots)
end

# The deepest subtree under `base` that stays under a cell budget — derived from
# the system's own aperture rather than hardcoded, so aperture 4 gets depth 8 and
# aperture 7 gets depth 5 without this file knowing either number.
function deep_depth(sys, base::Int, budget::Int = 70_000)
    d = 0
    while base + d + 1 <= max_level(sys) &&
        ncells(levelgrid(sys, base + d + 1)) ÷ ncells(levelgrid(sys, base)) <= budget
        d += 1
    end
    return d
end

@testset "subtree iterators" begin

    for sys in systems()
        name = string(nameof(typeof(sys)))
        base, roots = sweep_roots(sys)

        @testset "$name: agreement with the eager verbs and the definition" begin
            for c in roots, d in 0:3
                l = base + d
                l <= max_level(sys) || continue
                for conn in (Vertex(), Edge())
                    rim = collect(EdgeCellIterator(sys, c, l; connectivity = conn))
                    inner = collect(InnerCellIterator(sys, c, l; connectivity = conn))

                    # The eager verbs ARE these, so this pins the wrapper.
                    @test rim == subtree_border(sys, c, l; connectivity = conn)
                    @test inner == subtree_interior(sys, c, l; connectivity = conn)

                    # ...and both are the definition.
                    oracle = brute_force_border(sys, c, l; connectivity = conn)
                    @test rim == oracle
                    @test inner == [x for x in descendants(sys, c, l)
                                    if !(x in Set(oracle))]
                end
            end
        end

        @testset "$name: agreement with the generic engine" begin
            for c in roots, d in 0:3
                l = base + d
                l <= max_level(sys) || continue
                for conn in (Vertex(), Edge())
                    @test collect(EdgeCellIterator(sys, c, l; connectivity = conn)) ==
                          generic_rim(sys, c, l, conn)
                    @test collect(InnerCellIterator(sys, c, l; connectivity = conn)) ==
                          generic_interior(sys, c, l, conn)
                end
            end
        end

        @testset "$name: partition" begin
            for c in roots, d in 0:3
                l = base + d
                l <= max_level(sys) || continue
                for conn in (Vertex(), Edge())
                    rim = collect(EdgeCellIterator(sys, c, l; connectivity = conn))
                    inner = collect(InnerCellIterator(sys, c, l; connectivity = conn))
                    all_of_them = descendants(sys, c, l)

                    @test issorted(rim)
                    @test issorted(inner)
                    @test isempty(intersect(Set(rim), Set(inner)))
                    @test sort(vcat(rim, inner)) == sort(all_of_them)
                    @test length(rim) + length(inner) == length(all_of_them)
                end
            end
        end

        @testset "$name: counted, and inferred" begin
            for c in roots, d in 0:3
                l = base + d
                l <= max_level(sys) || continue
                for it in (EdgeCellIterator(sys, c, l), InnerCellIterator(sys, c, l))
                    @test eltype(it) == cellindextype(sys)
                    if Base.IteratorSize(typeof(it)) isa Base.HasLength
                        @test length(it) == length(collect(it))
                    else
                        # The contract is that no `length` walks the subtree to
                        # answer, so an uncounted iterator must not have one.
                        @test_throws MethodError length(it)
                    end
                    # Type stability where there is anything to yield.
                    isempty(collect(it)) || @test (@inferred first(it)) isa
                                                  cellindextype(sys)
                end
            end
        end

        @testset "$name: validation" begin
            c = first(roots)
            lc = level(c)
            if lc > 0
                @test_throws ArgumentError EdgeCellIterator(sys, c, lc - 1)
                @test_throws ArgumentError InnerCellIterator(sys, c, lc - 1)
            end
            @test_throws ArgumentError EdgeCellIterator(sys, c, max_level(sys) + 1)
            @test_throws ArgumentError InnerCellIterator(sys, c, max_level(sys) + 1)

            # Depth zero: a cell is its own rim, and has no interior at all.
            @test collect(EdgeCellIterator(sys, c, lc)) == [c]
            @test isempty(collect(InnerCellIterator(sys, c, lc)))
        end

        # A5 is excluded by construction, not by accident: with no automaton its
        # engine is a vector built in the constructor, so these assertions would
        # pass while measuring nothing. See this file's header.
        if name != "A5System"
            @testset "$name: lazy" begin
                c = first(roots)
                deep = deep_depth(sys, base)
                @test deep >= 4                    # or the measurement is not one
                shallow_l = base + deep - 3
                deep_l = base + deep

                for T in (EdgeCellIterator, InnerCellIterator)
                    shallow_it = T(sys, c, shallow_l)
                    deep_it = T(sys, c, deep_l)
                    length(collect(deep_it)) >= 64 || continue

                    lazy = lazy_bytes(deep_it, 4)
                    eager = eager_bytes(deep_it)

                    # The point of the whole task, as a ratio rather than a byte
                    # count: four cells off the front cost a fraction of the
                    # walk they are the front of.
                    @test lazy * 8 < eager

                    # ...and the sharper form of it. Between these two depths the
                    # rim grows by the aperture cubed, and a walk whose state is
                    # its depth must not notice. This is the assertion that
                    # would fail if the iterator ever went back to building the
                    # rim and handing out its elements.
                    @test lazy <= lazy_bytes(shallow_it, 4) + 64

                    # The walk is resumable, not restarted: the prefix agrees
                    # with the collected form element for element.
                    prefix = eltype(deep_it)[]
                    for x in deep_it
                        push!(prefix, x)
                        length(prefix) >= 4 && break
                    end
                    @test prefix == collect(deep_it)[1:4]
                end
            end
        end
    end

    # The wrapper is a system too, and the rim of a subtree is a question about
    # the hierarchy that the ellipsoid does not touch — so both walks must
    # forward to the wrapped system and answer identically.
    @testset "AuthalicSystem forwards both walks" begin
        for sys in systems()
            wrapped = AuthalicSystem(sys)
            base, roots = sweep_roots(sys)
            for c in roots, d in 0:2
                l = base + d
                l <= max_level(sys) || continue
                @test collect(EdgeCellIterator(wrapped, c, l)) ==
                      collect(EdgeCellIterator(sys, c, l))
                @test collect(InnerCellIterator(wrapped, c, l)) ==
                      collect(InnerCellIterator(sys, c, l))
                @test subtree_border(wrapped, c, l) == subtree_border(sys, c, l)
                @test subtree_interior(wrapped, c, l) == subtree_interior(sys, c, l)
            end
        end
    end
end

end # module SubtreeIteratorTests
