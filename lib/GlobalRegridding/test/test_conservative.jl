# `Conservative()` weight construction. Owned by task T3.

@testset "Conservative weights" begin
    @test Conservative() isa AbstractRegriddingMethod
end
