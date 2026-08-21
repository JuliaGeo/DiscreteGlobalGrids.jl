module CopDEMSourceModeTests

using Test

include(joinpath(@__DIR__, "..", "..", "scripts", "copdem_source_mode.jl"))

@testset "synthetic CopDEM source is absolute" begin
    @test effective_realspec(:synthetic, :none) === :none
    @test_throws ArgumentError effective_realspec(:synthetic, :auto)
    @test_throws ArgumentError effective_realspec(:synthetic,
        ["Copernicus_DSM_COG_30_S90_00_E000_00_DEM"])
end

@testset "real CopDEM source selection is unchanged" begin
    stems = ["Copernicus_DSM_COG_30_N00_00_E006_00_DEM"]
    @test effective_realspec(:real, :auto) === :auto
    @test effective_realspec(:real, :none) === :none
    @test effective_realspec(:real, stems) === stems
    @test_throws ArgumentError effective_realspec(:hybrid, :none)
end

end # module
