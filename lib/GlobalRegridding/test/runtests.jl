using Test
using GlobalRegridding
import GlobalRegridding as GR
import LinearAlgebra

# Shared analytic test doubles.
include("toyspaces.jl")

# Method used to test the fallback error.
struct UnimplementedMethod <: AbstractRegriddingMethod end

@testset "GlobalRegridding" begin

    @testset "qualified space extension contract" begin
        hooks = (
            :subtree,
            :chunkextents, :chunkextent, :chunkindex, :candidatechunks!,
            :chunkranges,
            :chartaxes, :chartcoords, :chartlocalindex, :chartperiod, :chartspacing,
            :destinationdims, :dimsource, :_asspace,
        )
        integration_hooks = (
            :resolvespatialdims,
            :_prepare_raster_transform_pair, :_task_prepared_raster_transform,
        )
        @test all(name -> Base.ispublic(GR, name), (hooks..., integration_hooks...))
        docs = Base.Docs.meta(GR)
        @test all(name -> haskey(docs, Base.Docs.Binding(GR, name)), hooks)
    end

    @testset "toy space contract" begin
        space = ToyLonLatSpace(8, 4; chunks = (4, 2))

        # Chunks partition cell indices.
        covered = reduce(vcat, [collect(ownedindices(space, c)) for c in 1:nchunks(space)])
        @test sort(covered) == collect(1:ncells(space))

        # Full-width chunks are contiguous.
        rows = ToyLonLatSpace(8, 4; chunks = (8, 2))
        @test ownedindices(rows, 2) isa AbstractUnitRange
        @test ownedindices(rows, 2) == 17:32

        # Point location inverts cell centroids.
        @test all(cellat(space, cellcentroid(space, i)) == i for i in 1:ncells(space))
        @test cellat(space, toy_point(0, 0)) isa Int

        # Partial spaces return `nothing` outside coverage.
        patch = ToyLonLatSpace(4, 2; lat = (0.0, 40.0))
        @test cellat(patch, toy_point(0, -10)) === nothing

        # Generic `chunkat` inverts `ownedindices` for indices and points.
        @test all(GR.chunkat(space, i) == c
                  for c in 1:nchunks(space) for i in ownedindices(space, c))
        @test GR.chunkat(space, cellcentroid(space, 5)) == GR.chunkat(space, 5)
        @test GR.chunkat(patch, toy_point(0, -10)) === nothing

        # Analytic and polygon areas both cover the sphere.
        global_space = ToyLonLatSpace(8, 4)
        @test sum(graticule_area(global_space, i) for i in 1:ncells(global_space)) ≈ 4pi
        @test sum(GO.area(manifold(global_space), getcell(global_space, i))
                  for i in 1:ncells(global_space)) ≈ 4pi

        # Cell rings are counter-clockwise from outside.
        oriented = GOCore.Spherical(; radius = 1.0, oriented = true)
        @test GO.area(oriented, getcell(global_space, 1)) < 2pi

        # Regional chunk extents cover their cell corners.
        region = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
            chunks = (4, 2))
        caps = chunktree(region).caps
        @test all(GR.US._contains(caps[c], p)
                  for c in 1:nchunks(region)
                  for i in ownedindices(region, c)
                  for p in cellcorners(region, i))
        @test maximum(cap.radius for cap in caps) < pi / 2
    end

    @testset "weight blocks" begin
        space = ToyLonLatSpace(4, 2)
        inds = ownedindices(space, 1)
        coo = WeightCOO(length(inds))
        buildweights!(coo, ToyDiagonalMethod(; scale = 2.0), space, inds, space, inds)
        block = WeightBlock(coo, length(inds), length(inds))

        # Weights use chunk-local index order.
        @test size(block) == (8, 8)
        @test Matrix(block.weights) == 2.0 * Matrix(LinearAlgebra.I, 8, 8)
        @test block.denom == fill(2.0, 8)

        # Denominators remain optional.
        bare = WeightCOO(length(inds))
        buildweights!(bare, ToyDiagonalMethod(; withdenom = false), space, inds, space, inds)
        @test WeightBlock(bare, length(inds), length(inds)).denom === nothing
    end

    @testset "unimplemented hooks" begin
        space = ToyLonLatSpace(4, 2)
        inds = ownedindices(space, 1)
        # Missing method implementations produce a useful error.
        @test_throws "UnimplementedMethod" buildweights!(
            WeightCOO(length(inds)), UnimplementedMethod(), space, inds, space, inds)
        @test supportradius(UnimplementedMethod(), space) == 0.0
        @test_throws ArgumentError Weighted(1.5)
    end

    include("test_rastergrid.jl")
    include("test_proj.jl")
    include("test_conservative.jl")
    include("test_interpolation.jl")
    include("test_executor.jl")
    include("test_chunkgraph.jl")
    include("test_lazy.jl")
    include("test_integration.jl")
end
