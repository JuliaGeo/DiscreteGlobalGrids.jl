# ---------------------------------------------------------------------------
# T12 — ISEA4R. STUB.
#
# Canonical id is `LevelIndex` over `diamond * 4^level + morton` (0-based,
# dense), `ISEA4RSystem <: AbstractHierarchicalGridSystem`. Port from
# `src/ISEA4R/chart.jl` (xyd <-> morton, xyd <-> row-major) and
# `src/ISEA4R/Isea4rKernel.jl` (`ordinal_to_cell`, boundary/center/cap/
# polygon); the old `test/ISEA4R/` suite is the oracle. ISEA9R stays parked.
#
# The shared ISEA basis is ALREADY PORTED at `src/systems/ISEA/` (T4 brought it
# over with IGeo7). Do not re-copy it — include it the way `IGeo7.jl` does, via
# the `isdefined(@__MODULE__, :ISEA) || include(...)` guard, so the include
# order of `src/systems/` stays irrelevant.
#
#   - Morton nesting makes positions subtree-sorted:
#     `has_sorted_subtrees = true`, parent/children by shift.
#   - Neighbours: the rhombus lattice plus the 10-diamond icosahedral seam
#     topology (derive the adjacency tables, verify geometrically via shared
#     boundary points). Moore 8 / Edge 4, CCW cycle from a documented start,
#     per the rotational contract in `src/interface/grid.jl`.
#   - Children tile the parent rhombus in chart space and Snyder is a per-face
#     homeomorphism, so nesting is exact: same `node_extent` choice and
#     verification as T11.
# ---------------------------------------------------------------------------

module ISEA4R

# TODO(T12): port the ISEA4R system onto the new interface.

end # module ISEA4R
