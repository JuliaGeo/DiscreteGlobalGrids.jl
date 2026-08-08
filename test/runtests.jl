using Test
using DiscreteGlobalGrids

# Each unit suite is wrapped in its own module (see e.g. test/A5/runtests.jl)
# so that same-named exports (cell_center, lonlat_to_cell, ...) from different
# systems never collide in a shared namespace. `include` evaluates at the top
# level of Main, so the module definitions inside these files are legal even
# from within the @testset.
@testset "DiscreteGlobalGrids.jl" begin
    include("test_helpers.jl")
    include("core/runtests.jl")
    include("A5/runtests.jl")
    include("H3/runtests.jl")
    include("HEALPix/runtests.jl")
    include("IGeo7/runtests.jl")
    include("ISEA4R/runtests.jl")
    include("S2/runtests.jl")
end
