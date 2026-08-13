# ---------------------------------------------------------------------------
# T3 — Conformance harness. STUB.
#
# The executable form of the interface contracts: `test_grid_interface(grid)`
# and `test_hierarchical_system(sys)`, property suites a third-party
# implementor runs with two calls. Cursor-free by design — it exercises
# interface primitives only, so it does not depend on the T2 fallbacks.
#
# The function names are declared (and re-exported from the package) here so
# that T3 only ever edits this file and its test file — never the shared
# module or test entry points.
# ---------------------------------------------------------------------------

module Conformance

# TODO(T3): implement the conformance property suites as methods of these.

"""
    test_grid_interface(grid; n_samples, rng)

Property-test `grid` against the base-interface contracts: the
`cellindex`/`cellposition` bijection over sampled positions, boundary rings of
unit-norm points (implicitly closed, CCW from outside), centroids inside their
cell's covering cap, `cellat(cell_centroid(grid, c)) == c` where `cellat` is
implemented, and determinism of repeated calls.

Runs as a `Test.@testset`; a third-party grid implementor calls this and
[`test_hierarchical_system`](@ref) to validate an implementation.
"""
function test_grid_interface end

"""
    test_hierarchical_system(sys; levels, n_samples, rng)

Property-test `sys` against the hierarchical contracts: `parent`/`children`
inverses, child-count sanity, `levelgrid` round-trips, the **covering law**
(sampled descendants' boundaries contained in every ancestor's
[`node_extent`](@ref), several levels down), neighbor determinism and symmetry
under both connectivities, and `descendant_range` completeness when
`has_sorted_subtrees(sys)`.

Runs as a `Test.@testset`; sampled, seeded, reproducible.
"""
function test_hierarchical_system end

export test_grid_interface, test_hierarchical_system

end # module Conformance
