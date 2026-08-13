# ---------------------------------------------------------------------------
# T2 — Fallback substrate. STUB.
#
# Every generic implementation of the interface declared in `src/interface/`
# lives under this directory, split into files `include`d from here:
# `HierarchicalGridCursor` and its SpatialTreeInterface / CR.Trees wiring, the
# fallback position-space quadtree, `PartialGrid`, `MultiOrderCellSet`, the
# geometry generics (`cell_polygon`, `cell_area`, `cell_extent`, the default
# `node_extent`), generic `cellat`, geometric `neighbors`/`ring`, the derived
# `ancestor`/`descendants`, and the `query` engine.
#
# Extend the parent's generics from in here, e.g.
#
#     import ..DiscreteGlobalGrids as DGG
#     import ..DiscreteGlobalGrids: AbstractGrid, cell_polygon, cell_boundary
#
# and define methods on those imported bindings; they are the same functions
# the package exports, so a method added here is the package's method.
# ---------------------------------------------------------------------------

module Fallbacks

# TODO(T2): implement the fallback substrate here.

end # module Fallbacks
