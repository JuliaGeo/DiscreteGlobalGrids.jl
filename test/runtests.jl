using Test
using DiscreteGlobalGrids

# Each unit suite is wrapped in its own module so that same-named exports from
# different systems never collide in a shared namespace. `include` evaluates at
# the top level of Main, so the module definitions inside these files are legal
# even from within the @testset.
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
    # Run interface-wide laws after each system's implementation tests.
    include("systems/crosssystem/runtests.jl")
    # Multi-order suites share the committed California outline fixture.
    include("systems/crosssystem/multiorder_polygons.jl")
    include("systems/crosssystem/multiorder_budget.jl")
    # Test the compressed cell container before the lookup that wraps it.
    include("systems/crosssystem/cell_vector.jl")
    # MOC storage, in dependency order: the mixed-level container, the
    # aggregation verbs that build one, then its cube face.
    include("systems/crosssystem/multiorder_vector.jl")
    include("systems/crosssystem/aggregate.jl")
    include("systems/crosssystem/multiorder_data.jl")
    include("systems/crosssystem/subtree_iterators.jl")
    include("systems/crosssystem/stencils.jl")
    include("systems/crosssystem/subtree_halos.jl")
    include("systems/crosssystem/dimensionaldata.jl")
    include("systems/crosssystem/regridding_conservation.jl")
    include("plotting/runtests.jl")
end
