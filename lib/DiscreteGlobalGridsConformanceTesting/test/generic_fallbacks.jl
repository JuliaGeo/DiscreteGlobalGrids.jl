# Correctness certification: the conformance laws run against the generic
# fallbacks, which are what answer for a system that implements no fast path.
# Not part of `Pkg.test()` — the geometric adjacency is a spatial-tree query and
# a boundary comparison per cell, against a system's own O(1) arithmetic — so it
# runs as its own CI job.
#
# Two systems, one per adjacency family: HEALPix is the quad-face lattice (its
# engines are shared with S2 and ISEA4R) and IGeo7 the aperture-7 hex with
# pentagons (shared with H3). Both declare `has_sorted_subtrees`, so the cursor
# behind the fallbacks prunes by `descendant_range`; A5, the one system that
# does not, sends it down the selection-mode path and is the next system to add
# here.

using Test
using DiscreteGlobalGrids
using DiscreteGlobalGridsConformanceTesting

@testset "generic fallbacks" begin
    test_generic_fallbacks(HEALPixSystem())
    test_generic_fallbacks(IGeo7System())
end
