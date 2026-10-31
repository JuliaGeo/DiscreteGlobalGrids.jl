# Correctness certification: the conformance laws run against the generic
# fallbacks, which are what answer for a system that implements no fast path.
# Not part of `Pkg.test()` — the geometric adjacency is a spatial-tree query and
# a boundary comparison per cell, against a system's own O(1) arithmetic — so it
# runs as its own CI job:
#
#   julia --project=test test/fallbacks/generic_fallbacks.jl
#
# Two systems, one per adjacency family: HEALPix is the quad-face lattice (its
# engines are shared with S2 and ISEA4R) and IGeo7 the aperture-7 hex with
# pentagons (shared with H3). Both declare `has_sorted_subtrees`, so the cursor
# behind the fallbacks prunes by `descendant_range`.
#
# The wrapper hides a system's fast paths by TYPE rather than by `invoke`: a
# system's specializations hang on `HierarchicalLevelGrid{TheSystem}` or on
# `TheSystem`, and the wrapped system appears as
# `HierarchicalLevelGrid{GenericFallbackSystem{TheSystem}}` — outside every one
# of those signatures. Only the required surface (identity, boundary, centroid,
# hierarchy, the covering methods and the traits) is forwarded; `cellat`,
# `neighbors`, `ring`, `ancestor`, `descendants` and the traversal engines are
# left to the interface-wide implementations.

module GenericFallbackCertification

using Test
using DiscreteGlobalGrids
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGrids: AbstractHierarchicalGridSystem, AbstractCellIndex, Connectivity
import GridInterface
include("../conformance_helpers.jl")
using .ConformanceHelpers: conformance_grid, conformance_system

struct GenericFallbackSystem{S<:AbstractHierarchicalGridSystem} <: AbstractHierarchicalGridSystem
    system::S
end

# The wrapped system's own complete-level grid: the entry point every consumer
# uses, so both implementor styles forward identically.
_inner(w::GenericFallbackSystem, l::Integer) = DGG.levelgrid(w.system, Int(l))

DGG.cellindextype(w::GenericFallbackSystem) = DGG.cellindextype(w.system)
DGG.levels(w::GenericFallbackSystem) = DGG.levels(w.system)
DGG.rootcells(w::GenericFallbackSystem) = DGG.rootcells(w.system)
Base.parent(w::GenericFallbackSystem, c::AbstractCellIndex) = Base.parent(w.system, c)
DGG.children(w::GenericFallbackSystem, c::AbstractCellIndex) = DGG.children(w.system, c)

DGG.node_extent(w::GenericFallbackSystem, c::AbstractCellIndex) =
    DGG.node_extent(w.system, c)
DGG.cap_inflation(w::GenericFallbackSystem) = DGG.cap_inflation(w.system)
DGG.has_sorted_subtrees(w::GenericFallbackSystem) = DGG.has_sorted_subtrees(w.system)
DGG.has_congruent_refinement(w::GenericFallbackSystem) =
    DGG.has_congruent_refinement(w.system)
DGG.descendant_range(w::GenericFallbackSystem, c::AbstractCellIndex, l::Integer) =
    DGG.descendant_range(w.system, c, l)
DGG.maxneighbors(w::GenericFallbackSystem, conn::Connectivity) =
    DGG.maxneighbors(w.system, conn)

# The shape traits are declarations about the cells, not fast paths: the
# wrapped level has its system's, so the ring and dual laws run over the
# fallbacks rather than being skipped as `GeneralCells`.
const _WrappedLevel = DGG.HierarchicalLevelGrid{<:GenericFallbackSystem}
GridInterface.cellshape(g::_WrappedLevel) = GridInterface.cellshape(_inner(DGG.system(g), DGG.level(g)))
GridInterface.edgekind(g::_WrappedLevel) = GridInterface.edgekind(_inner(DGG.system(g), DGG.level(g)))
GridInterface.tessellation(g::_WrappedLevel) = GridInterface.tessellation(_inner(DGG.system(g), DGG.level(g)))

DGG.ncells(w::GenericFallbackSystem, l::Integer) = DGG.ncells(_inner(w, l))
DGG.cellindex(w::GenericFallbackSystem, l::Integer, i::Int) = DGG.cellindex(_inner(w, l), i)
DGG.globalindex(w::GenericFallbackSystem, c::AbstractCellIndex) =
    DGG.globalindex(_inner(w, DGG.level(c)), c)
DGG.cell_boundary(w::GenericFallbackSystem, c::AbstractCellIndex) =
    DGG.cell_boundary(_inner(w, DGG.level(c)), c)
DGG.cell_centroid(w::GenericFallbackSystem, c::AbstractCellIndex) =
    DGG.cell_centroid(_inner(w, DGG.level(c)), c)

Base.show(io::IO, w::GenericFallbackSystem) =
    print(io, "GenericFallbackSystem(", w.system, ")")

# The coarsest levels: enough for the geometric fallbacks, and above the cell
# size at which the harness's own area tolerances give out.
function certify(sys; levels=(1, 2, 3))
    wrapped = GenericFallbackSystem(sys)
    grid = DGG.levelgrid(wrapped, first(levels))
    probe = DGG.cellindex(grid, 1)
    @testset "generic fallbacks: $(nameof(typeof(sys)))" begin
        # The wrapper reaches no specialization: every answer below is the
        # interface-wide one, while the contract the fallbacks consume is still
        # the system's.
        @test DGG.cell_boundary(grid, probe) ==
              DGG.cell_boundary(DGG.levelgrid(sys, first(levels)), probe)
        @test DGG.node_extent(wrapped, probe) == DGG.node_extent(sys, probe)
        @test !DGG.has_direct_location(wrapped)
        # The fallbacks' known costs: the extent of a cell is the cap of its
        # boundary (allocating), and the area of a cell is its straight-edged
        # polygon, whose chords differ between a parent and its children, so an
        # antichain mixing levels sums to 4πr² (1 + 8e-5).
        known = ["forwarding sweep over extent allocates nothing", "total area at most 4πr²"]
        for l in levels
            conformance_grid(DGG.levelgrid(wrapped, l); nsamples=8, known,
                name="$(nameof(typeof(sys))) level $l (generic fallbacks)")
        end
        conformance_system(wrapped; levels=levels, nsamples=4, known,
            name="$(nameof(typeof(sys))) (generic fallbacks)")
    end
end

@testset "generic fallbacks" begin
    certify(HEALPixSystem())
    certify(IGeo7System())
end

end # module GenericFallbackCertification
