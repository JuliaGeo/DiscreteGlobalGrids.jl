using Test
using DiscreteGlobalGrids

# Each unit suite is wrapped in its own module so that same-named exports from
# different systems never collide in a shared namespace. `include` evaluates at
# the top level of Main, so the module definitions inside these files are legal
# even from within the @testset.
#
# Only the redesigned suites run here. The old-architecture suites (test/core,
# test/A5, test/H3, test/HEALPix, test/IGeo7, test/ISEA4R, test/S2,
# test/test_helpers.jl) still exist on disk as reference for the ports and are
# expected to fail against this branch until T4-T7 land them; T8 deletes what
# the kill list retires.
@testset "DiscreteGlobalGrids.jl" begin
    include("interface/runtests.jl")
    # TODO(T2): include("fallbacks/runtests.jl")
    # TODO(T3): include("conformance/runtests.jl")
    # TODO(T4-T6): include("systems/{IGeo7,H3,HEALPix}/runtests.jl")
end
