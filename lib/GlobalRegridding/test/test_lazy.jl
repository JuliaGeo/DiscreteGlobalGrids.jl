# Chunk discovery, `LazyRegridArray`, streaming and spill. Owned by tasks T7, T8.

@testset "Lazy path" begin
    @test_throws ErrorException PerChunk(4)
end
