# `NearestCell()` and `BilinearPoint()` weight construction. Owned by task T4.

@testset "Interpolation weights" begin
    @test NearestCell() isa AbstractRegriddingMethod
    @test BilinearPoint() isa AbstractRegriddingMethod
end
