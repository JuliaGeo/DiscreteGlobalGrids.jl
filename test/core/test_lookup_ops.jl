# test/core/test_lookup_ops.jl — the generic lookup operations layer
# (src/core/lookup_ops.jl) and the neighbor kernel generics
# (src/core/kernel.jl): export/binding identities, the accessor pair, the
# unported-system error paths, and the allocation posture of the stencil
# sweep. Per-system neighbor *correctness* lives in the system suites
# (test/H3/test_neighbors.jl, test/IGeo7/test_neighbors.jl, the HEALPix
# suite); the spherical tree query behind `zonal` lives in
# test_tree_queries.jl.

module LookupOpsTests

using Test
using Random
import DimensionalData as DD
import GeoInterface as GI
using SmallCollections: SmallVector

using DiscreteGlobalGrids
const DGG = DiscreteGlobalGrids
using DiscreteGlobalGrids.HEALPix.HealpixLookups: HealpixLookup, Cells, HealpixLookups
using DiscreteGlobalGrids.H3.H3Lookups: H3Lookup
using DiscreteGlobalGrids.A5.A5Lookups: A5Lookup
using DiscreteGlobalGrids.IGeo7.IGeo7Lookups: IGeo7Lookup

@testset "lookup operations layer" begin
    @testset "export surface and shared bindings" begin
        exported = names(DiscreteGlobalGrids)
        for name in (:max_neighbors, :cell_neighbors, :dggs_system, :dggs_level,
                     :neighbor_indices, :stencil, :zonal)
            @test name in exported
        end
        # `zonal` and `stencil` are exported from `HealpixLookups` too — as
        # the SAME bindings, so a `using` of both namespaces (which every
        # HEALPix test file does) can never make them ambiguous.
        for name in (:zonal, :stencil)
            @test getproperty(DiscreteGlobalGrids, name) ===
                  getproperty(HealpixLookups, name)
        end
    end

    @testset "dggs_system / dggs_level over every lookup type" begin
        hl = HealpixLookup(Int64[3, 5]; level=2)
        @test dggs_system(hl) === HEALPixDGGS() && dggs_level(hl) == 2
        h3root = sort(DGG.root_ids(H3DGGS()))
        h3 = H3Lookup(h3root; resolution=0)
        @test dggs_system(h3) === H3DGGS() && dggs_level(h3) == 0
        ig = IGeo7Lookup(collect(DGG.root_ids(IGEO7DGGS())); resolution=0)
        @test dggs_system(ig) === IGEO7DGGS() && dggs_level(ig) == 0
        a5 = A5Lookup(sort(DGG.root_ids(A5DGGS())); resolution=0)
        @test dggs_system(a5) === A5DGGS() && dggs_level(a5) == 0
    end

    @testset "unported systems fail at the operation, not the accessor" begin
        # S2 is the remaining instance: no `cell_neighbors`, and no lookup type
        # either, so the gap can only be reached at the kernel.
        @test_throws NotPortedError max_neighbors(S2DGGS())
        # A5 was this testset's lookup-level example until `cell_neighbors` was
        # wired for it (`src/A5/A5Kernel.jl`). It now answers the whole chain,
        # which is what the assertions below pin — the accessor/operation split
        # this testset is named for currently has no lookup-reachable instance,
        # and re-acquiring one means a *new* system, not a regression in this.
        a5 = A5Lookup(sort(DGG.root_ids(A5DGGS())); resolution=0)
        @test max_neighbors(A5DGGS()) == 5
        halo = neighbor_indices(a5)
        @test length(halo) == length(a5)
        # The 12 res-0 cells are the dodecahedron's faces, so the lookup is
        # closed under adjacency: every neighbor resolves to a position in it
        # and none falls through to the `0` sentinel.
        @test all(nbs -> length(nbs) == 5 && all(!iszero, nbs), halo)
        A = DD.DimArray(fill(1.0, length(a5)), DD.Dim{:cells}(a5))
        @test stencil((c, nbs) -> c + sum(nbs), A) == fill(6.0, length(a5))
        # `zonal` used to be a NotPortedError here too (the old kernel
        # descent required descendant_range); the spherical tree query runs
        # on the selection-cursor fallback instead, so A5 answers.
        box = GI.Polygon([GI.LinearRing([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)])])
        @test zonal(sum, A; of=[box], boundary=:touches)[1] >= 1.0
    end

    @testset "the cell-dimension finder" begin
        hl = HealpixLookup(Int64[3, 5, 9]; level=2)
        box = GI.Polygon([GI.LinearRing([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0), (0.0, 0.0)])])
        # no DGGS dimension at all
        plain = DD.DimArray(rand(3), DD.X([1.0, 2.0, 3.0]))
        @test_throws ArgumentError stencil((c, nbs) -> c, plain)
        @test_throws ArgumentError zonal(sum, plain; of=[box])
        # stencil is 1-D over the cell dimension by contract
        A2 = DD.DimArray(rand(3, 2), (Cells(hl), DD.Ti([1.0, 2.0])))
        @test_throws ArgumentError stencil((c, nbs) -> c, A2)
        # ...while zonal selects along it and keeps the rest (covered per
        # system; here just the dim discovery on a 2-D array)
        @test length(zonal(sum, A2; of=[box])) == 1
    end

    # The geometry cap itself (now spherical, prepared-query backed) is
    # covered next door in test_tree_queries.jl, beside the descent it prunes.

    @testset "sorted insertion helper" begin
        v = SmallVector{6,Int}()
        for x in (5, 1, 3, 2, 4)
            v = DGG._insert_sorted(v, x)
        end
        @test collect(v) == [1, 2, 3, 4, 5]
        @test issorted(DGG._insert_sorted(v, 0)) && issorted(DGG._insert_sorted(v, 6))
    end

    @testset "stencil sweep does not allocate per cell" begin
        level = 4
        ids = collect(Int64, 0:(12 * 4^level - 1))
        l = HealpixLookup(ids; level)
        A = DD.DimArray(rand(length(ids)), Cells(l))
        nbi = neighbor_indices(l)
        smooth(c, nbs) = c + sum(nbs; init=zero(c))
        stencil(smooth, A; nbidx=nbi)               # compile
        allocated = @allocated stencil(smooth, A; nbidx=nbi)
        # the output vector plus rebuild overhead; a per-cell heap container
        # (the old Vector-based sweep) costs an order of magnitude more
        @test allocated < 4 * sizeof(Float64) * length(ids)
    end
end

end # module LookupOpsTests
