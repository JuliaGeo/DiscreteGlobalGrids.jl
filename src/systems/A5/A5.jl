# ---------------------------------------------------------------------------
# T10 — A5. STUB.
#
# `A5Cell <: AbstractCellIndex` wrapping the native `UInt64` serialization
# (origin + quintant + Hilbert bits from bit 58, `WORLD_CELL = 0`, the res-30
# special case), `A5System <: AbstractHierarchicalGridSystem`. Port from
# `src/A5/A5Native.jl` (the pure-Julia serialization) and `src/A5/A5Kernel.jl`
# (children/parent/descendants, pentagon neighbours, boundary/center,
# ordinals); the old `test/A5/` suite is the oracle.
#
# Two things make A5 the odd one out, deliberately:
#
#   - `has_sorted_subtrees = false`. A5 is the first real system to run on the
#     cursor's selection mode and on `MultiOrderCellSet`'s `(level, ordinal)`
#     fallback sort, so its suite must drive treeify/query/coverage through
#     those paths rather than leaving the substrate's mocks as the only
#     coverage.
#   - `cap_inflation` overrides to 1.75, and pentagons mean `max_neighbors` is
#     5 with `Edge()` == `Vertex()`.
#
# Neighbour order follows the rotational contract in `src/interface/grid.jl`:
# derive the CCW cycle from geometry (azimuth about the cell centre) and VERIFY
# the winding — the old kernel sorted by id, and both hex systems' "natural"
# orders turned out to need fixing.
# ---------------------------------------------------------------------------

module A5

# TODO(T10): port the A5 system onto the new interface.

end # module A5
