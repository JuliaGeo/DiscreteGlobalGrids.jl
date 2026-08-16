# Cross-system laws for lazy subtree-border and subtree-interior iterators.
#
# The tests compare each public iterator and eager collector with a brute-force
# definition, exercise the generic fallback directly, verify that border and
# interior partition the descendants, and check the iterator-size contract.
#
# Systems without sorted subtrees — A5 and the IVEA/RTEA rhombic family —
# materialize the subtree in the constructor, having neither sorted subtree
# ranges nor a boundary automaton. Their value and iterator-size laws are
# tested, but they are excluded from prefix-allocation checks because
# construction has already paid the output-sized allocation. The exclusion is
# keyed on `has_sorted_subtrees`, not on a system name, so registering another
# materializing system does not silently make the law measure nothing.

module SubtreeIteratorTests

using Test
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: systems, levels, max_level, levelgrid, ncells,
    cellindex, neighbors, level, descendants, rootcells, cellindextype,
    subtree_border, subtree_interior, EdgeCellIterator, InnerCellIterator,
    AuthalicSystem, Vertex, Edge, Connectivity

# Brute-force subtree rim: descendants with at least one outside neighbour.
# This implementation does not share traversal code with the iterators.
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

# Call the least-specific method to bypass system-specific engines.
const GENERIC = Tuple{DGG.AbstractHierarchicalGridSystem,DGG.AbstractCellIndex,
    Int,Connectivity}

# Wrap the generic engine in the public iterator to test the same protocol.
generic_edge(sys, c, l, conn = Vertex()) = EdgeCellIterator(sys, c, Int(l), conn,
    invoke(DGG.rim_engine, GENERIC, sys, c, Int(l), conn))
generic_inner(sys, c, l, conn = Vertex()) = InnerCellIterator(sys, c, Int(l), conn,
    invoke(DGG.interior_engine, GENERIC, sys, c, Int(l), conn))

generic_rim(sys, c, l, conn) = collect(generic_edge(sys, c, l, conn))
generic_interior(sys, c, l, conn) = collect(generic_inner(sys, c, l, conn))

# Consume the first `n` elements without materializing the full iterator.
function take_n(it, n::Int)
    seen = 0
    for _ in it
        seen += 1
        seen >= n && break
    end
    return seen
end

# Include iterator construction so eager materialization is measured.
function build_and_take(T, sys, c, l, n::Int)
    return take_n(T(sys, c, l), n)
end

function build_and_collect(T, sys, c, l)
    return length(collect(T(sys, c, l)))
end

# Warm each function before measuring allocations.
function lazy_bytes(T, sys, c, l, n::Int)
    build_and_take(T, sys, c, l, n)
    return @allocated build_and_take(T, sys, c, l, n)
end

function eager_bytes(T, sys, c, l)
    build_and_collect(T, sys, c, l)
    return @allocated build_and_collect(T, sys, c, l)
end

# Level-0 roots exercise whole-face seams; level-1 roots exercise blocks with a
# mixture of seam and sibling boundaries.
sweep_bases(sys) = filter(b -> b <= max_level(sys), (0, 1))

function sweep_roots(sys, base::Int)
    grid = levelgrid(sys, base)
    return unique(vcat(sample_cells(grid, 4), irregular_cells(grid, 2)))
end

# Select the deepest subtree that remains within the cell budget.
function deep_depth(sys, base::Int, budget::Int = 70_000)
    # A scan fallback pays for every descendant's topology while collecting.
    # Keep that oracle large enough to demonstrate depth/resumability without
    # turning the cross-system law into a multi-minute geometry benchmark.
    root = first(sweep_roots(sys, base))
    probe_level = min(base + 1, max_level(sys))
    probe = EdgeCellIterator(sys, root, probe_level)
    if Base.IteratorSize(typeof(probe)) isa Base.SizeUnknown
        budget = min(budget, 7_000)
    end
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
                    # Closed-form automata and materialized fallbacks are counted;
                    # a lazy generic scan advertises SizeUnknown and carries no
                    # traversal-cost `length`.
                    size_trait = Base.IteratorSize(typeof(it))
                    @test size_trait isa Union{Base.HasLength,Base.SizeUnknown}
                    if size_trait isa Base.HasLength
                        @test length(it) == length(collect(it))
                    else
                        @test_throws MethodError length(it)
                    end
                    # `first` is defined only for nonempty iterators.
                    isempty(collect(it)) || @test (@inferred first(it)) isa
                                                  cellindextype(sys)
                end
            end
        end

        # The generic lazy scan has no closed-form count and therefore no
        # `length`; the `MethodError` is the assertion, not an accident. Systems
        # without sorted subtrees materialize the fallback instead of scanning,
        # which makes it legitimately counted, so they are excluded.
        if DGG.has_sorted_subtrees(sys) && base + 2 <= max_level(sys)
            @testset "$name: the generic walk is uncounted" begin
                l = base + 2
                for c in roots
                    for git in (generic_edge(sys, c, l), generic_inner(sys, c, l))
                        @test Base.IteratorSize(typeof(git)) isa Base.SizeUnknown
                        @test_throws MethodError length(git)
                    end
                    # Generic and specialized walks return the same cells.
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

        # Systems without sorted subtrees materialize during construction, so
        # iteration-only allocation is not a valid laziness measurement there.
        if DGG.has_sorted_subtrees(sys)
            @testset "$name: lazy" begin
                c = first(roots)
                deep = deep_depth(sys, base)
                @test deep >= 4                    # provides a three-level comparison
                shallow_l = base + deep - 3
                deep_l = base + deep

                for T in (EdgeCellIterator, InnerCellIterator)
                    deep_it = T(sys, c, deep_l)
                    length(collect(deep_it)) >= 64 || continue

                    lazy = lazy_bytes(T, sys, c, deep_l, 4)
                    eager = eager_bytes(T, sys, c, deep_l)

                    # A four-cell prefix must allocate much less than collection.
                    @test lazy * 8 < eager

                    # Prefix allocation stays flat as the rim grows by the
                    # aperture cubed over three levels. Closed-form automata are
                    # effectively constant-sized; a generic DFS scan may add a
                    # small O(depth) stack, but never O(subtree) storage.
                    slack = Base.IteratorSize(typeof(deep_it)) isa Base.SizeUnknown ?
                            2_048 : 64
                    @test lazy <= lazy_bytes(T, sys, c, shallow_l, 4) + slack

                    # Incremental iteration returns the same prefix as collection.
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

    # Authalic wrapping changes geometry but not subtree hierarchy.
    @testset "AuthalicSystem forwards both walks" begin
        for sys in systems()
            sys isa DGG.AuthalicSystem && continue
            for base in sweep_bases(sys)
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
end

end # module SubtreeIteratorTests
