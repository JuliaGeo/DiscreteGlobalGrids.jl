# T2 replaces this placeholder with the fallback-substrate suite.
module TestFallbacks

using Test
using DiscreteGlobalGrids

@testset "fallbacks (placeholder)" begin
    @test isdefined(DiscreteGlobalGrids, :Fallbacks)
end

end # module TestFallbacks
