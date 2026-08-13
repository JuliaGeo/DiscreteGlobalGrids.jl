using Test
using DiscreteGlobalGrids

# Each unit suite is wrapped in its own module so that same-named exports from
# different systems never collide in a shared namespace. `include` evaluates at
# the top level of Main, so the module definitions inside these files are legal
# even from within the @testset.
#
# These are all the suites there are: T8 deleted the old-architecture ones
# (test/core, test/A5, test/H3, test/HEALPix, test/IGeo7, test/ISEA4R,
# test/S2, test/test_helpers.jl) once their assertions had been ported. The
# DGGRID oracle vectors they carried now live in test/systems/IGeo7/vectors/.
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
    include("systems/A5/runtests.jl")
    include("systems/S2/runtests.jl")
    include("systems/ISEA4R/runtests.jl")
    # Last on purpose: the cross-system laws sweep `systems()`, so they are the
    # suite that grows when a system is registered rather than when one is
    # written. Running them after every per-system suite means a failure here is
    # unambiguously a *contract* failure, not a port that never worked.
    include("systems/crosssystem/runtests.jl")
end
