# Interface-layer tests.
#
# T1 ships declarations, contracts and trait defaults — no algorithms — so this
# suite checks exactly that: the module loads with the new surface and none of
# the old one, `LevelIndex` behaves, the trait defaults are the documented ones,
# and every generic that is only *declared* raises a `MethodError` on a type
# that has not implemented it (rather than silently answering `nothing`, which
# is the failure mode a bodiless method definition would have produced).

module InterfaceTests

using Test
using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG

# Types that implement nothing at all. Every generic must bounce off them.
struct UnimplementedGrid <: AbstractGrid end
struct UnimplementedSystem <: AbstractHierarchicalGridSystem end
struct UnimplementedIndex <: AbstractCellIndex end

"""
    DetachedDocProbe

The positive control for the docstring-coverage test below, and a live
specimen of a trap that cost the A5 port a review cycle: **a comment between a
docstring and the thing it documents silently detaches it.** No warning, no
error — the docstring is simply parsed as a free-standing string and thrown
away.

This one is detached on purpose. `docstring` must report it as undocumented,
which is what makes the coverage assertions mean anything.
"""
# THIS COMMENT IS THE POINT — do not remove it, it is what detaches the
# docstring above.
struct DetachedDocProbe end

const EXPORTED = filter(!=(:DiscreteGlobalGrids), names(DiscreteGlobalGrids))

# The public surface each shipped system contributes: its singleton, the grid
# type `levelgrid` returns, and its canonical id type. Derived from the
# registry rather than listed, so a newly registered system is covered without
# anyone remembering to add it here.
const SYSTEM_TYPE_NAMES = sort!(unique!(reduce(vcat,
    [[nameof(typeof(s)),
      nameof(typeof(DGG.levelgrid(s, first(DGG.levels(s))))),
      nameof(DGG.cellindextype(s))] for s in DGG.systems()])))

"""
    docstring(mod, name) -> String

Every docstring attached to `mod.name`, concatenated, or `""` if it has none.

Reads the doc metadata directly rather than calling `Base.Docs.doc`, whose
`Binding` method only exists once the `REPL` stdlib is loaded — true in an
interactive session, false in the process `Pkg.test` spawns. `Binding` resolves
through imports, so a re-exported name (`ncells`, `Intersects`) is found in the
module that defines it.
"""
function docstring(mod::Module, name::Symbol)
    binding = Base.Docs.Binding(mod, name)
    io = IOBuffer()
    for m in Base.Docs.modules
        meta = Base.Docs.meta(m; autoinit=false)
        meta === nothing && continue
        multidoc = get(meta, binding, nothing)
        multidoc === nothing && continue
        for sig in multidoc.order
            for chunk in multidoc.docs[sig].text
                print(io, chunk)
            end
        end
    end
    return String(take!(io))
end

@testset "module surface" begin
    @test AbstractGrid isa Type
    @test AbstractHierarchicalGridSystem isa Type
    @test AbstractCellIndex isa Type
    @test UnimplementedGrid() isa AbstractGrid
    @test UnimplementedSystem() isa AbstractHierarchicalGridSystem
    @test UnimplementedIndex() isa AbstractCellIndex

    # The interface names are all there...
    for n in (:ncells, :cellindex, :cell_boundary, :cell_centroid, :cellposition,
              :rawid, :reindex, :cellindextypes, :cell_polygon, :cell_area,
              :cell_extent, :getcell, :cellat, :neighbors, :ring, :treeify,
              :query, :system, :level, :cellindextype, :levels, :max_level,
              :levelgrid, :rootcells, :children, :node_extent, :cap_inflation,
              :max_neighbors, :has_sorted_subtrees, :ancestor, :descendants,
              :descendant_range, :LevelIndex, :Connectivity, :Vertex, :Edge,
              :authalic_sphere)
        @test n in EXPORTED
    end

    # ... and the old architecture is gone, not merely shadowed.
    for n in (:AbstractDGGS, :all_systems, :DGGSGrid, :DGGSCursor, :cell_neighbors,
              :cell_children, :cell_parent, :has_exact_subtree_cap,
              :supports_prefix_ranges, :subtree_grid, :HEALPixDGGS)
        @test !isdefined(DiscreteGlobalGrids, n)
    end

    # `treeify`/`ncells`/`getcell` are `Trees`' own bindings, not lookalikes:
    # a method added here is a method `ConservativeRegridding` will dispatch to.
    Trees = DGG.Trees
    @test DGG.treeify === Trees.treeify
    @test DGG.ncells === Trees.ncells
    @test DGG.getcell === Trees.getcell

    # DE9IM predicates are re-exported, and `parent` unwraps the target.
    @test Intersects <: DE9IMPredicate
    @test parent(Covers(:target)) === :target
    @test Intersects(:t) isa DE9IMPredicate

    # Every exported name carries a docstring...
    @test docstring(@__MODULE__, :UnimplementedGrid) == ""  # the probe can fail
    for n in EXPORTED
        @test !isempty(docstring(DiscreteGlobalGrids, n))
    end

    # ...and the coverage check can actually see the failure mode that matters.
    #
    # `UnimplementedGrid` above only shows that a name which never had a
    # docstring reads as empty. The shape that bit the A5 port is different and
    # much quieter: a docstring that WAS written, and was detached from its
    # definition by an interposed comment. Julia reports nothing at all; for a
    # method on a shared generic, `@doc` then shows the *interface* docstring,
    # so the name still looks documented. This pins that a detached docstring
    # reads as absent, which is what gives the loop above its teeth.
    @test docstring(@__MODULE__, :DetachedDocProbe) == ""

    # For the system types the loop above is enough — there is no generic for
    # their text to fall through from, so detachment shows up as empty. Pinning
    # that each one's docstring NAMES it closes the remaining gap: a name whose
    # only documentation is inherited or interface-level text would pass an
    # emptiness test and fail this one.
    #
    # This does not reach a detached docstring on a *method* of a shared
    # generic — `cellat`, `neighbors` — where the interface text keeps showing
    # through and no package-level check can tell the difference. Those are
    # pinned in the systems' own suites, which is where the specific method is
    # in scope; A5's suite does exactly that for all of its public names.
    for n in SYSTEM_TYPE_NAMES
        doc = docstring(DiscreteGlobalGrids, n)
        @test !isempty(doc)
        @test occursin(string(n), doc)
    end

    # Two contracts are load-bearing enough to pin in the docs themselves: the
    # covering law (without it, generic tree pruning is silently unsound) and
    # the position-vs-identity rule (without it, `Int` arguments are a coin
    # flip). Both are what a third-party implementor reads instead of the
    # design document.
    node_extent_doc = docstring(DiscreteGlobalGrids, :node_extent)
    @test occursin("covering law", lowercase(node_extent_doc))
    @test occursin("every descendant", node_extent_doc)
    @test occursin("every depth", node_extent_doc)

    grid_doc = docstring(DiscreteGlobalGrids, :AbstractGrid)
    @test occursin("position", lowercase(grid_doc))
    @test occursin("1:ncells(grid)", grid_doc)
end

@testset "LevelIndex" begin
    c = LevelIndex(3, 17)

    @test c isa AbstractCellIndex
    @test isbits(c)
    @test sizeof(LevelIndex) == 16

    # `level` is total on ids and returns an `Int` regardless of field width.
    @test level(c) === 3
    @test rawid(c) == 17

    # Integer arguments convert, and the fields are wide enough that a level
    # past an `Int8` (or a count past an `Int32`) is representable.
    @test LevelIndex(Int8(3), Int32(17)) === c
    @test level(LevelIndex(200, 1)) === 200
    @test rawid(LevelIndex(0, 5_000_000_000)) == 5_000_000_000

    # Value equality and hashing agree, as required of a cell index.
    @test c == LevelIndex(3, 17)
    @test hash(c) == hash(LevelIndex(3, 17))
    @test c != LevelIndex(3, 18)
    @test c != LevelIndex(4, 17)

    # Total order, lexicographic in (level, index).
    @test c < LevelIndex(3, 18)
    @test !(c < c)
    @test LevelIndex(3, typemax(Int32)) < LevelIndex(4, 0)
    @test sort([LevelIndex(4, 0), LevelIndex(3, 18), LevelIndex(3, 17)]) ==
          [LevelIndex(3, 17), LevelIndex(3, 18), LevelIndex(4, 0)]

    @test repr(c) == "LevelIndex(3, 17)"
end

@testset "connectivity singletons" begin
    @test Vertex() isa Connectivity
    @test Edge() isa Connectivity
    @test Vertex() === Vertex()
    @test Edge() === Edge()
    @test isbits(Vertex()) && isbits(Edge())
    @test Vertex() != Edge()
end

@testset "trait defaults" begin
    s = UnimplementedSystem()
    g = UnimplementedGrid()

    # Documented defaults, live.
    @test has_sorted_subtrees(s) === false
    @test cap_inflation(s) === 1.2

    # A standalone grid has no hierarchy, and says so rather than erroring.
    @test system(g) === nothing
    @test level(g) === nothing
end

@testset "unimplemented base generics throw MethodError" begin
    g = UnimplementedGrid()
    c = LevelIndex(0, 0)

    @test_throws MethodError ncells(g)
    @test_throws MethodError cellindex(g, 1)
    @test_throws MethodError cellindex(g, 1, LevelIndex)
    @test_throws MethodError cell_boundary(g, c)
    @test_throws MethodError cell_centroid(g, c)
    @test_throws MethodError cellposition(g, c)
    @test_throws MethodError cell_polygon(g, c)
    @test_throws MethodError cell_area(g, c)
    @test_throws MethodError cell_extent(g, c)
    @test_throws MethodError getcell(g, 1)
    @test_throws MethodError cellat(g, 0.0, 0.0)
    @test_throws MethodError neighbors(g, c)
    @test_throws MethodError neighbors(g, c, 2; connectivity=Edge())
    @test_throws MethodError ring(g, c, 1)
    @test_throws MethodError treeify(g)
    @test_throws MethodError query(g, Intersects(nothing))

    # Ids that implement nothing: `level`/`rawid` are required, not defaulted.
    @test_throws MethodError level(UnimplementedIndex())
    @test_throws MethodError rawid(UnimplementedIndex())
end

@testset "unimplemented system generics throw MethodError" begin
    s = UnimplementedSystem()
    c = LevelIndex(1, 0)

    @test_throws MethodError cellindextype(s)
    @test_throws MethodError levels(s)
    @test_throws MethodError levelgrid(s, 0)
    @test_throws MethodError rootcells(s)
    @test_throws MethodError parent(s, c)
    @test_throws MethodError children(s, c)
    @test_throws MethodError node_extent(s, c)
    @test_throws MethodError ancestor(s, c, 0)
    @test_throws MethodError descendants(s, c, 2)
    @test_throws MethodError descendant_range(s, c, 2)
    @test_throws MethodError reindex(LevelIndex, s, c)

    # The derived trait defaults are total only once their primitive is wired.
    @test_throws MethodError max_level(s)          # -> levels(s)
    @test_throws MethodError cellindextypes(s)     # -> cellindextype(s)
    @test_throws MethodError max_neighbors(s)      # -> max_neighbors(s, Vertex())
    @test_throws MethodError max_neighbors(s, Vertex())
end

end # module InterfaceTests
