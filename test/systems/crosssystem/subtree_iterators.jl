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

# The generic engine, wrapped back into the public iterator through its
# positional constructor, so the fallback can be asked every question the fast
# path is asked — including the ones about `IteratorSize` and `length`.
generic_edge(sys, c, l, conn = Vertex()) = EdgeCellIterator(sys, c, Int(l), conn,
    invoke(DGG.rim_engine, GENERIC, sys, c, Int(l), conn))
generic_inner(sys, c, l, conn = Vertex()) = InnerCellIterator(sys, c, Int(l), conn,
    invoke(DGG.interior_engine, GENERIC, sys, c, Int(l), conn))

generic_rim(sys, c, l, conn) = collect(generic_edge(sys, c, l, conn))
generic_interior(sys, c, l, conn) = collect(generic_inner(sys, c, l, conn))

# Iterating without collecting: the shape the LAZY law measures.
function take_n(it, n::Int)
    seen = 0
    for _ in it
        seen += 1
        seen >= n && break
    end
    return seen
end

# CONSTRUCTION IS INSIDE THE MEASUREMENT, deliberately. Measuring only the
# iteration of an already-built iterator cannot fail: an engine that materialises
# its subtree in the constructor — which is exactly the regression this law
# exists to catch, and exactly what A5 does — then walks a ready vector and
# measures zero. Building and taking four is the whole claim.
function build_and_take(T, sys, c, l, n::Int)
    return take_n(T(sys, c, l), n)
end

function build_and_collect(T, sys, c, l)
    return length(collect(T(sys, c, l)))
end

# `@allocated` behind a function barrier, and always after a warm-up call, so
# what is measured is the walk rather than its compilation.
function lazy_bytes(T, sys, c, l, n::Int)
    build_and_take(T, sys, c, l, n)
    return @allocated build_and_take(T, sys, c, l, n)
end

function eager_bytes(T, sys, c, l)
    build_and_collect(T, sys, c, l)
    return @allocated build_and_collect(T, sys, c, l)
end

# BOTH base levels, because they are different shapes of problem. A level-0 root
# is a whole face or base cell: its block is flush with every seam the system
# has, every mask bit stays set the whole way down, and on the square systems the
# rim is the face boundary whose neighbours live on another face. A level-1 root
# is an interior block with seams on some sides and siblings on others. Level 1
# alone would leave the flush-everywhere arithmetic unexercised.
sweep_bases(sys) = filter(b -> b <= max_level(sys), (0, 1))

function sweep_roots(sys, base::Int)
    grid = levelgrid(sys, base)
    return unique(vcat(sample_cells(grid, 4), irregular_cells(grid, 2)))
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

    for sys in systems(), base in sweep_bases(sys)
        sysname = string(nameof(typeof(sys)))
        name = "$sysname at level $base"
        roots = sweep_roots(sys, base)

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
                    # Every system's own walk is counted — the five automata by
                    # closed form, A5 by the vector it materialised.
                    @test Base.IteratorSize(typeof(it)) isa Base.HasLength
                    @test length(it) == length(collect(it))
                    # Type stability where there is anything to yield.
                    isempty(collect(it)) || @test (@inferred first(it)) isa
                                                  cellindextype(sys)
                end
            end
        end

        # The other half of COUNTED. The generic walk is a lazy scan with no
        # closed-form count, and the contract is that such a walk carries NO
        # `length` rather than one that would traverse to answer — so the
        # `MethodError` is the assertion, not an accident. A5 is excluded
        # because its fallback materialises instead of scanning, which makes it
        # legitimately counted.
        if DGG.has_sorted_subtrees(sys) && base + 2 <= max_level(sys)
            @testset "$name: the generic walk is uncounted" begin
                l = base + 2
                for c in roots
                    for git in (generic_edge(sys, c, l), generic_inner(sys, c, l))
                        @test Base.IteratorSize(typeof(git)) isa Base.SizeUnknown
                        @test_throws MethodError length(git)
                    end
                    # ...and it still answers what the fast path answers.
                    @test collect(generic_edge(sys, c, l)) ==
                          collect(EdgeCellIterator(sys, c, l))
                    @test collect(generic_inner(sys, c, l)) ==
                          collect(InnerCellIterator(sys, c, l))
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
        if sysname != "A5System"
            @testset "$name: lazy" begin
                c = first(roots)
                deep = deep_depth(sys, base)
                @test deep >= 4                    # or the measurement is not one
                shallow_l = base + deep - 3
                deep_l = base + deep

                for T in (EdgeCellIterator, InnerCellIterator)
                    deep_it = T(sys, c, deep_l)
                    length(collect(deep_it)) >= 64 || continue

                    lazy = lazy_bytes(T, sys, c, deep_l, 4)
                    eager = eager_bytes(T, sys, c, deep_l)

                    # The point of the whole task, as a ratio rather than a byte
                    # count: building the walk and taking four cells off the
                    # front costs a fraction of the walk they are the front of.
                    @test lazy * 8 < eager

                    # ...and the sharper form of it. Between these two depths the
                    # rim grows by the aperture cubed, and a walk whose state is
                    # its depth must not notice. This is the assertion that would
                    # fail if the iterator ever went back to building the rim and
                    # handing out its elements.
                    @test lazy <= lazy_bytes(T, sys, c, shallow_l, 4) + 64

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
        for sys in systems(), base in sweep_bases(sys)
            wrapped = AuthalicSystem(sys)
            roots = sweep_roots(sys, base)
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
