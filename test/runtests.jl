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
    # T15: the same treatment for multi-order coverage, against a committed
    # real-world outline rather than a mock's box. Last because it is the
    # slowest, and separate from the file above because its laws are about one
    # operation rather than about the interface at large.
    include("systems/crosssystem/multiorder_polygons.jl")
    # T18: the same outline, the same seven systems, and the other mode of the
    # same query — a cell BUDGET instead of a depth. Straight after the file
    # above because it shares its fixture and its samplers, and because the two
    # laws it cannot state everywhere are the ones that file already measures.
    include("systems/crosssystem/multiorder_budget.jl")
    # T19: the compressed cell collection the cube axis is made of, exercised
    # with no DimensionalData in sight. Before the axis file because it is the
    # axis's substrate: a failure here is the compression, and the file below
    # would then fail for a reason that is not its own.
    include("systems/crosssystem/cell_vector.jl")
    # T20: the lazy form of the subtree rim and interior, which since T20 is the
    # only form — the eager verbs are `collect` of these. After the file above
    # for the same reason it is after the per-system suites: those suites check
    # each automaton against its own oracle, so a failure here is the iterator
    # protocol or the sharing, not the walk.
    include("systems/crosssystem/subtree_iterators.jl")
    # T21: adjacency on a SUBSET, which is the complete level's answer clipped
    # to membership. Straight after the file above because the rooted halo table
    # is built from those two iterators, so a failure there is the walk and a
    # failure here is the clipping.
    include("systems/crosssystem/stencils.jl")
    # T22: the OUTSIDE face of the boundary `subtree_iterators.jl` walks the
    # inside of. After both files above on purpose: the halo's laws are stated
    # against `subtree_border` and against the subset adjacency the stencil file
    # pins, so a failure upstream of here is the rim walk or the clipping, and a
    # failure here is the halo walk itself.
    include("systems/crosssystem/subtree_halos.jl")
    # T16: the DimensionalData cell axis, which is the multi-order set above
    # read as a cube dimension. After that file for the same reason it is after
    # the per-system suites — a failure here is the lookup, not the coverage.
    include("systems/crosssystem/dimensionaldata.jl")
    # T17: and the same treatment for regridding, in BOTH directions. The
    # source direction was the only one anything checked; the destination
    # direction is silently non-conservative on every system whose cell rings
    # are non-convex, for a reason that lives in GeometryOps' clipper. Those
    # arms are `@test_broken` and the file says why.
    include("systems/crosssystem/regridding_conservation.jl")
    include("plotting/runtests.jl")
end
