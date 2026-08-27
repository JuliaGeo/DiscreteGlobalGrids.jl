# Regridding simplification Phase 1

Date: 2026-08-22

## Outcome

Lazy chunk discovery now queries the spatial hierarchy already owned by each
space instead of requiring every space to construct a flat cap-tree adapter.
The common private seam is:

```julia
index = chunkindex(source_space)
candidatechunks!(out, index, destination_cap; radius)
```

The index is intentionally private and query-oriented. Its internal node
extents do not have to share one type.

## Implemented paths

- Raster storage ownership remains in the `DiskArrays.eachchunk` ranges already
  captured by `RasterGrid`. A zero-copy `RasterGridView` presents the real
  task-local coordinate transform and storage-order numbering through
  ConservativeRegridding's existing `AbstractCurvilinearGrid` interface.
  `TopDownQuadtreeCursor` supplies the implicit spatial hierarchy; traversal
  stops when a node belongs to one DiskArrays chunk.
- Generic spaces pack their chunk caps into GeometryOps'
  `FlexibleRTrees.RTree`. Leaves retain caps for the exact angular predicate;
  outward-rounded XYZ boxes provide the hierarchical broad phase. Empty input
  has a private zero-node index.
- DGG spaces are their own chunk index. Ordinary complete and partial grids use
  `HierarchicalGridCursor` over the original grid and stop at `chunklevel`, so
  ancestor extents continue to cover all stored leaves. CopernicusDEM reuses
  its existing block cursor over the complete tile lattice and filters frontier
  IDs through the selected `PartialGrid` IDs, avoiding a 26,475-root scan per
  query.

`CapQuery` and `DGGChunkTree` were removed. `LazyRegridArray.srctree` became
`srcindex`; no eager weight, cell-tree, or source-data behavior changed.
`RasterFlatTree` remains for scattered-cell trees and the legacy `chunktree`
compatibility method; production chunk discovery does not construct it.

## Correctness gates

Tests cover:

- R-tree result identity with brute-force cap intersection under node
  capacities 2, 3, and 16;
- outward XYZ bounds for tiny, polar, antimeridian, and whole-sphere caps;
- raster whole-sphere coverage, both storage numberings, arbitrary callable
  transforms, DiskArrays-derived chunk ranges, and zero source reads;
- raster boundary-point no-false-negative checks;
- DGG whole-sphere coverage and self-reachability;
- exact DGG frontier ID/range correspondence for complete and sparse grids; and
- CopernicusDEM's separate block-cursor path.

Verification through the Julia MCP:

- GlobalRegridding: 3,118 passed, 1 expected broken;
- DGG cross-system regridding: passed;
- regridding conservation: 76 passed, 12 expected broken;
- regridding acceptance: 22 passed; and
- CopernicusDEM system suite: 16,258 passed, 3 expected broken.

## Graph performance gate retained for Phase 2

Phase 1 does not replace the materialized dependency graph's latitude-sorted
builder. The exact four-thread production baseline after Phase 1 is:

- 26,475 source chunks;
- 66,175 destination chunks;
- 326,386 edges;
- 3,352,520-byte graph;
- 0.0581 s minimum and 0.0594 s median construction over five samples; and
- 10,654,960 minimum allocated bytes.

A provisional graph builder using one Phase 1 hierarchy query per destination
was slower: 0.1500 s median and 13,973,280 minimum allocated bytes. It produced
326,064 edges and a 3,349,944-byte graph. The 322-edge difference is consistent
with hierarchy pruning of unused cap over-coverage, but it is not itself proof
that no geometrically real pair was dropped.

That prototype was backed out. Phase 2 must add actual-cell no-miss
instrumentation, retain deterministic CSR output, and at least recover the
current graph-build performance before the dependency graph adopts the native
indexes. Until then, lazy discovery uses the native hierarchy while scheduling
keeps the established conservative cap join.
