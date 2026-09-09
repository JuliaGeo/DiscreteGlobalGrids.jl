# Interface-layer tests.
#
# The suite checks the exported surface, `LevelIndex`, documented trait
# defaults, and the `MethodError` behavior of unimplemented interface methods.

module InterfaceTests

using Test
using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG

# Types that implement nothing at all. Every generic must bounce off them.
struct UnimplementedGrid <: AbstractGrid end
struct UnimplementedSystem <: AbstractHierarchicalGridSystem end
struct UnimplementedIndex <: AbstractCellIndex end

# A system that declares its id scheme and its levels, and stops there. Its
# methods are below, beside the testset that reads them.
struct IdentifiedSystem <: AbstractHierarchicalGridSystem end

"""
    DetachedDocProbe

A positive control for the docstring-coverage test. The intervening comment
deliberately detaches this docstring from the type, so `docstring` must report
the type as undocumented.
"""
# This comment must remain between the docstring and the type definition.
struct DetachedDocProbe end

# `names` returns both tiers on Julia >= 1.11, so the split is made explicitly:
# `EXPORTED` is what `using DiscreteGlobalGrids` brings in, `PUBLIC` is that
# plus the names reachable only through the module path.
const PUBLIC = filter(!=(:DiscreteGlobalGrids), names(DiscreteGlobalGrids))
const EXPORTED = filter(n -> Base.isexported(DiscreteGlobalGrids, n), PUBLIC)

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

    # Required interface names are exported.
    for n in (:ncells, :cellindex, :cell_boundary, :cell_centroid, :localindex,
              :globalindex, :rawid, :reindex, :cellindextypes,
              :cell_area, :cell_extent, :getcell, :cellat, :neighbors, :ring,
              :treeify, :query, :system, :level, :cellindextype, :levels,
              :maxlevel, :levelgrid, :rootcells, :children, :node_extent,
              :maxneighbors, :has_sorted_subtrees, :has_congruent_refinement,
              :has_direct_location,
              :ancestor, :descendants,
              :descendant_range, :LevelIndex, :Connectivity, :Vertex, :Edge,
              :cellsize, :levelfor, :subtree, :halo, :border, :interior,
              :adjacency, :AdjacencyTable, :halocells, :haloindices)
        @test n in EXPORTED
    end

    # The `public` tier: documented and reachable by module path, deliberately
    # absent from `using`. Machinery a caller names only to talk ABOUT it.
    for n in (:EdgeCellIterator, :InnerCellIterator, :SubtreeHaloIterator,
              :SubsetHaloIterator, :HaloIndexIterator, :RegionSide,
              :sizehint, :halo_indices, :SubsetIndexedCell,
              :HierarchicalGridCursor, :StorageOrder, :NeighborCallbackError,
              :cap_inflation, :directioncode, :authalic_sphere,
              :StoreSnapshot, :StoreDescription, :ArrayEntry, :ChunkManifest,
              :GridReference, :CONVENTION_REGISTRY, :DEFAULT_WRITE_CONVENTIONS,
              :ENCODING_REGISTRY, :GRID_REFERENCE)
        @test n in PUBLIC
        @test !Base.isexported(DiscreteGlobalGrids, n)
    end

    # The type every boundary and centroid method returns, and the
    # index-space covering verb the zonal recipe calls, are reachable from
    # the top namespace — `using DiscreteGlobalGrids` and no module path.
    @test :UnitSphericalPoint in EXPORTED
    @test UnitSphericalPoint === DGG.GO.UnitSpherical.UnitSphericalPoint
    @test :covering_indices in EXPORTED

    # Retired interface names are not defined.
    for n in (:AbstractDGGS, :all_systems, :DGGSGrid, :DGGSCursor, :cell_neighbors,
              :cell_children, :cell_parent, :has_exact_subtree_cap,
              :supports_prefix_ranges, :subtree_grid, :HEALPixDGGS,
              :subtree_halo, :subtree_border, :subtree_interior, :halo_table,
              :stencil_table, :HaloTable, :StencilTable, :halo_sizehint,
              :is_contained, :max_level, :max_neighbors)
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

    # Every name of either tier carries a docstring.
    @test docstring(@__MODULE__, :UnimplementedGrid) == ""  # the probe can fail
    for n in PUBLIC
        @test !isempty(docstring(DiscreteGlobalGrids, n))
    end

    # An interposed comment detaches a docstring, which must read as absent.
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
    # the index-vs-identity rule (without it, `Int` arguments are a coin
    # flip). Both must be explicit in the public documentation.
    node_extent_doc = docstring(DiscreteGlobalGrids, :node_extent)
    @test occursin("covering law", lowercase(node_extent_doc))
    @test occursin("every descendant", node_extent_doc)
    @test occursin("every depth", node_extent_doc)

    grid_doc = docstring(DiscreteGlobalGrids, :AbstractGrid)
    @test occursin("index", lowercase(grid_doc))
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

@testset "the radix-4 quad-face family" begin
    # S2, HEALPix and ISEA4R get their hierarchy, level-grid arithmetic and
    # subtree engines from `AbstractQuadFaceGridSystem` alone. A system that
    # silently detached from the supertype would keep compiling and start
    # answering with the generic fallbacks instead.
    @test AbstractQuadFaceGridSystem <: AbstractHierarchicalGridSystem
    for s in (S2System(), HEALPixSystem(), ISEA4RSystem())
        @test s isa AbstractQuadFaceGridSystem
    end

    # No other registered system reaches the family's arithmetic: a non-member
    # selects its own `rootcells`, never the supertype's. Asked of the method
    # table rather than of the type names, so a system renamed or added is
    # covered without anyone editing a list.
    for s in DGG.systems()
        s isa AbstractQuadFaceGridSystem && continue
        m = which(DGG.rootcells, Base.typesof(s))
        @test m.sig.parameters[2] !== AbstractQuadFaceGridSystem
    end
end

@testset "trait defaults" begin
    s = UnimplementedSystem()
    g = UnimplementedGrid()

    # Documented defaults, live.
    @test has_sorted_subtrees(s) === false
    @test has_congruent_refinement(s) === false
    @test DGG.cap_inflation(s) === 1.2

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
    @test_throws MethodError globalindex(g, c)
    @test_throws MethodError DGG.cell_polygon(g, c)
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

# `AbstractHierarchicalGridSystem`'s docstring splits the implementor surface
# into a required half and a defaulted half, and `IdentifiedSystem` is the split
# made executable: it declares identity and levels and nothing else, so every
# required method must refuse to answer and every defaulted one must answer
# anyway. A default handed to a required method, or taken from a defaulted one,
# rewrites the documented contract without touching the documentation.
DGG.cellindextype(::IdentifiedSystem) = LevelIndex
DGG.levels(::IdentifiedSystem) = 0:2

@testset "system required/defaulted split" begin
    s = IdentifiedSystem()
    c = LevelIndex(1, 0)

    # Declared, so they answer.
    @test cellindextype(s) === LevelIndex
    @test levels(s) == 0:2

    # Required: identity and hierarchy.
    @test_throws Exception rootcells(s)
    @test_throws Exception parent(s, c)
    @test_throws Exception children(s, c)

    # Required: the five level-grid primitives, in their system-level arity.
    @test_throws Exception ncells(s, 1)
    @test_throws Exception cellindex(s, 1, 1)
    @test_throws Exception globalindex(s, c)
    @test_throws Exception cell_boundary(s, c)
    @test_throws Exception cell_centroid(s, c)

    # Defaulted, and answering for a system that declared none of them.
    @test maxneighbors(s, Vertex()) === nothing
    @test levelgrid(s, 1) === HierarchicalLevelGrid(s, 1)
    @test DGG.cap_inflation(s) === 1.2
    @test maxlevel(s) == 2
    @test has_sorted_subtrees(s) === false
    # `node_extent`'s default needs the geometry this system does not have, so
    # only its existence is pinned here; the covering law it must satisfy is in
    # the fallbacks suite.
    @test hasmethod(node_extent, Tuple{AbstractHierarchicalGridSystem,AbstractCellIndex})
end

@testset "unimplemented system generics throw MethodError" begin
    s = UnimplementedSystem()
    c = LevelIndex(1, 0)

    @test_throws MethodError cellindextype(s)
    @test_throws MethodError levels(s)
    # The default `levelgrid` validates against the unimplemented `levels(s)`.
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
    @test_throws MethodError maxlevel(s)          # -> levels(s)
    @test_throws MethodError cellindextypes(s)     # -> cellindextype(s)

    # `maxneighbors` is total: a system that declares no static degree bound
    # answers `nothing`, and the subset neighbour machinery buffers its
    # one-rings in a `Vector` instead of a `SmallVector`.
    @test maxneighbors(s) === nothing             # -> maxneighbors(s, Vertex())
    @test maxneighbors(s, Vertex()) === nothing
    @test maxneighbors(s, Edge()) === nothing

    # Standalone grids have no system declaration to forward, so the grid
    # form preserves the same explicit "no bound" answer.
    g = UnimplementedGrid()
    @test maxneighbors(g) === nothing
    @test maxneighbors(g, Vertex()) === nothing
    @test maxneighbors(g, Edge()) === nothing
end

end # module InterfaceTests

include("sizing.jl")
