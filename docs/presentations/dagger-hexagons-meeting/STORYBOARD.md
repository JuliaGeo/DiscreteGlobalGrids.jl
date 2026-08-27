# Dagger for Global Grids — storyboard

Ten-minute meeting primer. By the end, the room should understand how DGGS
storage creates chunk-aware work, see the two workloads that expose it, and be
ready to choose a small Dagger integration experiment.

## Story

| Beat | Point | Graphic |
| --- | --- | --- |
| 1. Many global grids | A DGGS covers the globe with cells that have stable addresses. The cell shape and topology vary by system. | **WebGL Makie:** comparable globes for Tripolar, IGeo7, HEALPix and ISEA. |
| 2. The hierarchy | Every IGeo7 address belongs to a nested hierarchy. Coarse cells contain successively finer cells. | **WebGL Makie:** three IGeo7 levels overlaid on one globe, with finer outlines drawn at decreasing opacity. |
| 3. Hierarchy becomes storage | One level-5 cell can be a storage chunk whose contents are its level-12 descendants. The logical chunk is a polygon, while its cells occupy a contiguous address interval. | **WebGL Makie:** the level-5 parent beside the detailed level-12 footprint, with the parent outline overlaid and the corresponding address interval below. |
| 4. Address order is geometric | The contiguous Z7 interval is also a space-filling walk through the descendants. | **WebGL Makie:** successive level-12 centroids inside one level-5 chunk, beside a zoomed global Z7 curve. |
| 5. Chunks are geometric | Across the globe, storage is partitioned into polygonal level-5 chunks, each resolved at level 12. These are not rectangular raster blocks. | **WebGL Makie:** an autoplaying camera moves from the globe to one chunk and the level-12 outlines of all its neighboring chunks. |
| 6a. One-map regridding | Construct a regridder once, then apply it to a single two-dimensional map. | **WebGL Makie:** source raster globe → sparse overlap weights → destination DGGS globe. |
| 6b. Long-cube regridding | Reuse the spatial plan across many variables and timestamps. NASA's 120 m ITS_LIVE data is the motivating scale. | **WebGL Makie:** a fixed spatial plan while a field rotates once around both the source and regridded destination globes. |
| 7. The graph we already have | Regridding already maps each destination chunk to the source chunks that may contribute. Today that graph drives ordering, caching, prefetch and validation inside one process. | **WebGL Makie:** source-chunk graph beside the IGeo7 destination, with contributing source footprints outlined in red on the globe. |
| 8. A stencil really runs | For an out-of-memory grid, process owned cells in storage order, read the one-ring halo, and write the computed value. Spatial neighbours may live in distant chunks. | **Video:** a real one-ring mean advances monotonically through one owned interval; the globe highlights each neighbourhood and the storage strip shows its halo reads. |
| 9. The grid supplies the halo | For any chunk, DiscreteGlobalGrids can identify the chunks—and cells—needed for its halo. | **WebGL Makie:** owned cells, halo cells and supplying storage chunks linked in one view. |
| 10. The runtime question | Both workloads already describe work and dependencies. Dagger could own placement, movement and multiprocessing. | **Video:** regrid dependencies and stencil halos collapse into the same task shape, then fan out to Dagger workers. |
| 11. The meeting decision | Choose one thin vertical slice and decide what dependency information should cross the boundary. | **WebGL Makie:** three candidate contracts: chunk graph, chunk-plus-cell map, or higher-level task plan. |

## Graphic rule

Every beat has one explanatory visual. WebGL Makie is for spatial inspection
and browser-driven camera movement; video is for streaming and execution order.
FlyThroughPaths defines camera paths for both forms.

## Scope

- Introductory, not a Dagger tutorial.
- Regridding and stencils are the two driving examples.
- No proposed public API yet.
- Close with a concrete experiment for the meeting to choose.
