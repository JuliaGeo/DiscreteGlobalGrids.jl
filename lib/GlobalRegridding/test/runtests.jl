using Test
using GlobalRegridding
import GlobalRegridding as GR
import LinearAlgebra

# Shared analytic test doubles.
include("toyspaces.jl")

# Method used to test the fallback error.
struct UnimplementedMethod <: AbstractRegriddingMethod end

@testset "GlobalRegridding" begin

    @testset "execution hook contract" begin
        # A regrid side is any GridInterface grid; its cells, trees, chunks,
        # charts and dual cells are read through GridInterface's verbs. What a
        # package may still extend here are the execution hooks: plain
        # functions with `Any` defaults, public and documented.
        hooks = (
            :chunkranges, :knownempty, :samplerstate, :supportradius,
            :destinationdims, :dimsource, :rebuildoutput,
            :sourcemissingval, :outputmissingval, :destinationmissingval,
            :overlappairs, :locatecell, :chunkindex, :candidatechunks!,
            :_asgrid,
        )
        integration_hooks = (
            :resolvespatialdims,
            :_prepare_raster_transform_pair, :_task_prepared_raster_transform,
            :_crs_native_to_unit_sphere, :_projected_crs_native_to_unit_sphere,
        )
        # What a package's own point sampler (`samplerstate` + `weightsat!`)
        # writes its rows and statuses with.
        sampler_hooks = (
            :weightsat!, :Sampler, :WeightRow, :addentry!, :snapunit,
            :meanvalueweights!, :WeightStatus, :ismapped,
        )
        @test all(name -> Base.ispublic(GR, name),
            (hooks..., integration_hooks..., sampler_hooks...,
             :WeightsMapped, :WeightsRim, :WeightsOutside, :WeightsDegenerate))
        docs = Base.Docs.meta(GR)
        @test all(name -> haskey(docs, Base.Docs.Binding(GR, name)),
            (hooks..., integration_hooks..., sampler_hooks...))
        # The `Any` defaults, callable on anything that is a grid.
        space = ToyLonLatSpace(8, 4; chunks = (4, 2))
        @test GR.knownempty(zeros(2, 2), (1:2, 1:2)) === false
        @test GR.supportradius(ToyDiagonalMethod(), space) == 0.0
        @test GR.destinationdims(space, DD.Lookups.Points()) === nothing
        @test GR.dimsource(1:4) === nothing
        # `chunkranges` reads a chunk's positions off `chunks(grid)`.
        @test GR.chunkranges(space, 2, (8, 4)) == (5:8, 1:2)
        @test_throws ArgumentError GR.chunkranges(space, 1, (ncells(space),))
        bands = ToyLonLatSpace(8, 4; chunks = (8, 2))
        @test GR.chunkranges(bands, 1, (ncells(bands),)) == (1:16,)
        # The grid-shaped hooks are gone, not merely unexported.
        for name in (:RegridSpace, :celltree, :getcell, :cellcentroid, :cellneighbors,
                     :cellcap, :cellcorners, :hasanalyticlocation, :hascellchart,
                     :hasdualcells, :chartaxes, :chartcoords, :chartlocalindex,
                     :chartperiod, :chartspacing, :nchunks, :ownedindices,
                     :cellindices, :chunkextents, :chunkextent, :manifold,
                     :expensivecellgeometry, :_asspace, :RasterGridView,
                     :CellSpaceRTree)
            @test !isdefined(GR, name)
        end
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
        @test all(chunkat(space, i) == c
                  for c in 1:nchunks(space) for i in ownedindices(space, c))
        @test chunkat(space, cellcentroid(space, 5)) == chunkat(space, 5)
        @test chunkat(patch, toy_point(0, -10)) === nothing

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
        caps = chunkextents(region)
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
    include("test_rastercrs.jl")
    include("test_conservative.jl")
    include("test_interpolation.jl")
    include("test_barycentric.jl")
    include("test_executor.jl")
    include("test_chunkgraph.jl")
    include("test_lazy.jl")
    include("test_directnearest.jl")
    include("test_integration.jl")
end
