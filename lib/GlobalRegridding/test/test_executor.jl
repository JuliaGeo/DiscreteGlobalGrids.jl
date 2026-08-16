# Plans, the eager executor, and the user API. Owned by task T5.

@testset "Executor" begin
    @test_throws ErrorException regrid(zeros(2, 2); to = ToyLonLatSpace(2, 2))
end
