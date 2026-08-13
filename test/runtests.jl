using Test
using DiscreteGlobalGrids

# Each unit suite is wrapped in its own module so that same-named exports from
# different systems never collide in a shared namespace. `include` evaluates at
# the top level of Main, so the module definitions inside these files are legal
# even from within the @testset.
#
# Only the redesigned suites run here. The old-architecture suites (test/core,
# test/A5, test/H3, test/HEALPix, test/IGeo7, test/ISEA4R, test/S2,
# test/test_helpers.jl) still exist on disk as reference for the remaining
# ports and are expected to fail against this branch; T8 deletes what the kill
# list retires.
#
# The conformance harness self-tests live in their own workspace package under
# lib/. Each system suite imports that package and runs its two public suites
# on its own system, alongside that system's own oracle vectors.
@testset "DiscreteGlobalGrids.jl" begin
    include("interface/runtests.jl")
    include("fallbacks/runtests.jl")
    include("fallbacks/authalic.jl")
    include("systems/IGeo7/runtests.jl")
    include("systems/H3/runtests.jl")
    include("systems/HEALPix/runtests.jl")
    include("systems/crosssystem/runtests.jl")
    # TODO(T10): include("systems/A5/runtests.jl")
    # TODO(T11): include("systems/S2/runtests.jl")
    # TODO(T12): include("systems/ISEA4R/runtests.jl")
end
