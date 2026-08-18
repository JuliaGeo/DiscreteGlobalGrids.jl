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
    # CopernicusDEM is not in `systems()`, so the cross-system sweeps below
    # never see it; its laws are stated only here.
    include("systems/CopernicusDEM/runtests.jl")
    # Run interface-wide laws after each system's implementation tests.
    include("systems/crosssystem/runtests.jl")
    # Multi-order suites share the committed California outline fixture.
    include("systems/crosssystem/multiorder_polygons.jl")
    include("systems/crosssystem/multiorder_budget.jl")
    # Test the compressed cell container before the lookup that wraps it.
    include("systems/crosssystem/cell_vector.jl")
    include("systems/crosssystem/subtree_iterators.jl")
    include("systems/crosssystem/stencils.jl")
    include("systems/crosssystem/neighborhood.jl")
    include("systems/crosssystem/mapneighbors.jl")
    include("systems/crosssystem/subtree_halos.jl")
    # The algebra over regions reads the halo walk and both cell containers.
    include("systems/crosssystem/region_algebra.jl")
    include("systems/crosssystem/dimensionaldata.jl")
    include("systems/crosssystem/regridding_conservation.jl")
    # The GlobalRegridding face reads the grids, the cell containers, and the
    # cube axis, so it runs after all three.
    include("systems/crosssystem/regrid.jl")
    # Acceptance: the tiled-DEM, south-pole, streaming-and-spill case, on the
    # face the file above unit-tests.
    include("systems/crosssystem/regrid_acceptance.jl")
    include("io/runtests.jl")
    include("plotting/runtests.jl")
end
