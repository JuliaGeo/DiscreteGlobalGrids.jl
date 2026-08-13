# ---------------------------------------------------------------------------
# T11 — S2. STUB.
#
# Canonical id is `LevelIndex` over the scaffold ordinal
# `face * 4^level + hilbert` (0-based, dense), `S2System <:
# AbstractHierarchicalGridSystem`. Port from `src/S2/chart.jl` (xyf <-> hilbert
# and xyf <-> row-major codecs) and `src/S2/S2Kernel.jl` (`ordinal_to_cell`,
# boundary/center/cap/polygon); the old `test/S2/` suite is the oracle.
#
# This port COMPLETES the surface rather than merely rewiring it: the old
# wiring had no `cell_to_ordinal`, no neighbours, and no hierarchy.
#
#   - Hilbert nesting should make positions subtree-sorted
#     (`has_sorted_subtrees = true`, parent/children = `div(p, 4)` / `4p + k`) —
#     VERIFY the child-order claim against the chart codec before relying on it.
#   - Neighbours: within-face lattice plus cube-edge adjacency for seam
#     crossings. Moore 8 / Edge 4, CCW cycle from a documented start, verified
#     geometrically (edge neighbours share >= 2 boundary points, corner
#     neighbours exactly 1).
#   - S2 quads subdivide exactly, so the cell's own cap may serve as
#     `node_extent`: verify the nesting, then document whether this ships the
#     exact cap or the inflated default.
#   - The native S2 cellid (lsb-sentinel encoding) as an alternate scheme via
#     `reindex` is OPTIONAL; the scaffold ordinal stays canonical either way.
# ---------------------------------------------------------------------------

module S2

# TODO(T11): port the S2 system onto the new interface.

end # module S2
