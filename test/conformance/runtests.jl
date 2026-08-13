# T3 replaces this placeholder with the conformance-harness self-test.
module TestConformance

using Test
using DiscreteGlobalGrids

@testset "conformance (placeholder)" begin
    @test test_grid_interface isa Function
    @test test_hierarchical_system isa Function
end

end # module TestConformance
