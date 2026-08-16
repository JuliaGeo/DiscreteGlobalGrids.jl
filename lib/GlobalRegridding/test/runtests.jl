using Test
using GlobalRegridding
import GlobalRegridding as GR
import LinearAlgebra

# The analytic test doubles every other file here builds on.
include("toyspaces.jl")

# A method that never implements `build_weights!`, for the fallback error path.
struct UnimplementedMethod <: AbstractRegriddingMethod end

@testset "GlobalRegridding" begin

    @testset "toy space contract" begin
        space = ToyLonLatSpace(8, 4; chunks = (4, 2))

        # Chunks partition the cell positions exactly. Everything downstream —
        # discovery, block assembly, accumulation — double-counts or drops cells
        # if this slips.
        covered = reduce(vcat, [collect(cellindices(space, c)) for c in 1:nchunks(space)])
        @test sort(covered) == collect(1:ncells(space))

        # A chunk spanning the full longitude row is contiguous in position
        # space, and says so.
        rows = ToyLonLatSpace(8, 4; chunks = (8, 2))
        @test cellindices(rows, 2) isa AbstractUnitRange
        @test cellindices(rows, 2) == 17:32

        # `cellat` inverts `cellcentroid`: the lattice mapping is not transposed
        # and the cell boxes tile without gaps.
        @test all(cellat(space, cellcentroid(space, i)) == i for i in 1:ncells(space))
        @test cellat(space, toy_point(0, 0)) isa Int

        # A partial space does not answer outside its coverage.
        patch = ToyLonLatSpace(4, 2; lat = (0.0, 40.0))
        @test cellat(patch, toy_point(0, -10)) === nothing

        # `chunkat` inverts `cellindices` for a space that defines no method of
        # its own — the fallback scans the chunks, and this is the space that
        # exercises it. The point form is `cellat` composed with it and is
        # `nothing` exactly where `cellat` is.
        @test all(GR.chunkat(space, i) == c
                  for c in 1:nchunks(space) for i in cellindices(space, c))
        @test GR.chunkat(space, cellcentroid(space, 5)) == GR.chunkat(space, 5)
        @test GR.chunkat(patch, toy_point(0, -10)) === nothing

        # A global space is exactly the sphere, by its analytic areas and by its
        # polygons, which tile without gap or overlap.
        global_space = ToyLonLatSpace(8, 4)
        @test sum(graticule_area(global_space, i) for i in 1:ncells(global_space)) ≈ 4pi
        @test sum(GO.area(manifold(global_space), getcell(global_space, i))
                  for i in 1:ncells(global_space)) ≈ 4pi

        # Rings are counter-clockwise seen from outside: read with orientation
        # honoured, a cell is the region it bounds and not its complement.
        oriented = GOCore.Spherical(; radius = 1.0, oriented = true)
        @test GO.area(oriented, getcell(global_space, 1)) < 2pi

        # Every chunk extent covers the geometry of the cells the chunk owns.
        # Discovery prunes on these, so a cap that does not cover silently drops
        # pairs. Checked on a regional space, where the caps are tighter than
        # the whole sphere and the law can actually fail.
        region = ToyLonLatSpace(8, 4; lon = (-40.0, 40.0), lat = (-20.0, 20.0),
            chunks = (4, 2))
        caps = chunktree(region).caps
        @test all(GR.US._contains(caps[c], p)
                  for c in 1:nchunks(region)
                  for i in cellindices(region, c)
                  for p in cellcorners(region, i))
        @test maximum(cap.radius for cap in caps) < pi / 2
    end

    @testset "weight blocks" begin
        space = ToyLonLatSpace(4, 2)
        inds = cellindices(space, 1)
        coo = WeightCOO(length(inds))
        build_weights!(coo, ToyDiagonalMethod(; scale = 2.0), space, inds, space, inds)
        block = WeightBlock(coo, length(inds), length(inds))

        # Chunk-local indices, in the order the builder was handed them.
        @test size(block) == (8, 8)
        @test Matrix(block.weights) == 2.0 * Matrix(LinearAlgebra.I, 8, 8)
        @test block.denom == fill(2.0, 8)

        # A builder that reports no denominator produces a block without one.
        bare = WeightCOO(length(inds))
        build_weights!(bare, ToyDiagonalMethod(; withdenom = false), space, inds, space, inds)
        @test WeightBlock(bare, length(inds), length(inds)).denom === nothing
    end

    @testset "unimplemented hooks" begin
        space = ToyLonLatSpace(4, 2)
        inds = cellindices(space, 1)
        # A method with no `build_weights!` says which method and which spaces.
        @test_throws "UnimplementedMethod" build_weights!(
            WeightCOO(length(inds)), UnimplementedMethod(), space, inds, space, inds)
        @test support_radius(UnimplementedMethod(), space) == 0.0
        @test_throws ArgumentError Weighted(1.5)
    end

    include("test_rastergrid.jl")
    include("test_conservative.jl")
    include("test_interpolation.jl")
    include("test_executor.jl")
    include("test_lazy.jl")
    include("test_integration.jl")
end
