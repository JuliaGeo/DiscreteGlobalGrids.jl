# The tiled raster cursor

- Date: 2026-08-25
- Scope: `src/engine/tiled_raster.jl` (new), `src/interface/grid.jl`,
  `src/systems/CopernicusDEM/cursor.jl`, and the Copernicus DEM tests.
- Commit: `Descend a quadtree inside each raster tile`

## 1. A grid of raster tiles

- Neither a lattice nor a hierarchy: a collection of raster tiles, each a
  rectangle of pixels, whose pixels together are the grid's cells. A Copernicus
  DEM holding is one (1x1-degree tiles of 1200x1200 or 3600x3600 pixels); so is
  a STAC catalogue, whose items may overlap and carry their own projections.
- It has exactly two structures — a flat set of tiles with footprints, and a
  rectangle of pixels inside each — and the tree is those two, in that order.

## 2. The interface

A grid opts in by answering `raster_tiles`; the other three are then required.
All four are declared in `src/interface/grid.jl`.

| hook | contract |
| --- | --- |
| `raster_tiles(grid, inds) -> tiles or nothing` | The tiles covering grid indices `inds`, as an indexable collection of opaque **tile handles**; `nothing` (the default) when this grid's cells are not raster pixels. The tiles must hold every cell of `inds` exactly once, and no cell outside it. |
| `raster_shape(grid, tile) -> (nrows, ncols)` | The tile's rectangle in pixels. Both positive. Rows and columns are numbered `0:nrows-1` and `0:ncols-1` everywhere else. |
| `raster_localindex(grid, tile, j, i) -> Int` | The index **in `grid`** of the pixel at row `j`, column `i`. It is the tile's offset in the grid plus the pixel's row-major index within the tile, so a tile's pixels are one contiguous block of grid indices in row-major order. |
| `raster_cap(grid, tile, j0, j1, i0, i1) -> SphericalCap` | A cap containing the geometry — boundaries, not merely centres — of every pixel in the closed sub-rectangle `j0:j1` x `i0:i1`. |

- A handle is opaque to the engine, which reads nothing out of it and names tiles
  by caps and shapes alone: **no lattice coordinate reaches the engine**.
- `raster_cap` is asked for whole tiles at build time, interior rectangles during
  descent, and single pixels at the leaves. A closed-form box (Copernicus DEM has
  one, from its latitude bands) answers in constant time; without one a grid
  merges per-pixel caps, correct but linear, and wants a memo of its own.

## 3. The tree

- **Tile layer.** `RasterTileTree` sorts tiles by the Morton key of their cap
  centres and splits them `RASTER_TILE_ARITY`-ways per level down to one tile per
  node, storing each node's cap (the merge of its children's) and pixel count.
  It is `IndexTree`'s packing reused, not re-derived. Packed caps rather than
  `FlexibleRTrees.RTree`, whose Cartesian `Extent`s would need a conversion at
  every node, as `CellSpaceRTree` pays. It assumes nothing about arrangement.
- **Pixel layer.** A node holding one tile *is* that tile's whole raster, so no
  level is spent on the tile itself; descending bisects the rectangle along its
  longer axis until a node holds `LEAF_CELLS` pixels or fewer. That is
  `BlockCursor`'s descent and the same code: `bisect_parts`, `rect_part`,
  `LEAF_CELLS` and `LeafCells` now live in the engine.
- **Reuse, not duplication.** `BlockCursor`'s law `index = id - origin` holds
  only where a rectangle's cells are one consecutive id run, which a tile set
  crossing a latitude row is not. The law was generalised instead:
  `raster_localindex` subsumes it as the one-tile case. `BlockCursor` survives
  for the complete lattice and its windows, where a tiled tree would need a table
  of 64 800 tile caps against an `O(1)` box from two divisions; its `PartialGrid`
  constructor is gone, so no shape has two trees.
- **Extents.** Tile node extents are stored and nest, each an `Extents.union` of
  its children's. Raster node extents derive from `raster_cap` and are memoized
  per task by `RasterExtentMemo` — direct-mapped, 1024 slots, whole-key compared,
  keyed on `(tile slot, j0, j1, i0, i1)`, cleared when a task turns to another
  tree — so `node_extent_is_expensive` is `false`. Node caps bound cell
  **geometry**: they cover the pixels' boundaries, not those pixels' own caps.

## 4. `treeify` and `subcursor`

- `treeify(::PartialGrid{<:CopernicusDEMSystem})` builds the tiled tree for any
  holding — one tile, a row, a square, an L, four tiles on four continents. The
  `HierarchicalGridCursor` fallback stays reachable only for an id the lattice
  does not have, which `raster_tiles` answers with `nothing`.
- `subcursor(…, inds)` builds the tree over `inds` alone, with leaf indices still
  the grid's own at every node, which is `subcursor`'s contract.
- `raster_tiles` walks the ids by runs and cuts each into at most three
  rectangles — a part row, the whole rows under it, a part row. A run is found by
  binary search: ids ascend strictly, so `id(p+k) - id(p) - k` never decreases.
  A level-0 holding is a set of 1x1 rectangles.

## 5. What a STAC-backed system implements

Four methods and one handle type; nothing else changes, because overlapping
footprints, per-item projections and heterogeneous sizes are already what a layer
that packs caps packs.

- handle: the item's identifier, raster size and offset in the grid — the grid's
  cells are the catalogue's pixels in item order;
- `raster_tiles`: the items whose pixel blocks meet `inds`, in that order, from
  the catalogue rather than any lattice;
- `raster_shape`: the item's raster dimensions;
- `raster_localindex`: the item's offset plus `j * ncols + i + 1`;
- `raster_cap`: the sub-rectangle's corners through the item's own geotransform
  and CRS, padded for edge bow as the Copernicus lattice's is — a cap fitted to
  four corners alone does not cover the pixels' boundaries.

## 6. The root's shape

Packed root against a flat one (a single node with every tile as a child), on a
scaled twin lattice of 30x30-pixel tiles, over 2000 random cell centroids per
point query and a dual-tree join against a complete HEALPix level, best of two
runs:

| tiles | cells | root | ns per point query | dual-tree join |
| --- | --- | --- | --- | --- |
| 2x2 (4 tiles) | 3 600 | packed | 3 999 | 0.0766 s |
| 2x2 (4 tiles) | 3 600 | flat | 4 068 | 0.0717 s |
| 4x4 (16 tiles) | 14 400 | packed | 5 661 | 0.0825 s |
| 4x4 (16 tiles) | 14 400 | flat | 5 594 | 0.0825 s |
| 32x32 (1024 tiles) | 643 200 | packed | 9 028 | - |
| 32x32 (1024 tiles) | 643 200 | flat | 30 691 | - |

- They cost the same at the tile counts a holding of tens has — a four-tile
  packed root is one level over its tiles and tests about as many caps as a flat
  root — and their joins hold the same stored entries at both sizes.
- They diverge with the tile count, a flat root testing every tile on every
  query. At 1024 tiles the packed root is 3.4 times faster on a point query, and
  the ratio grows; a global Copernicus holding is 26 450 tiles, so the root shape
  is the difference between 26 450 cap tests per query and a few dozen.

## 7. The gate

Source: Copernicus DEM GLO-90 tiles. Destination: the IGeo7 level-12 cells
covering the same box, chunked at their level-5 ancestors, taking the chunk with
the most candidate source chunks. Two routes to the source tree: **chunked**, the
pair being one source chunk's indices, and **whole**, the whole source space with
the plan reaching `treeify(grid)` itself. Nanoseconds per destination cell, best
of two runs, before and after interleaved on an otherwise idle machine:

| shape | method | route | cells | before | after | after / before |
| --- | --- | --- | --- | ---: | ---: | ---: |
| one tile | Conservative | chunked | 512 | 6 578 | 6 691 | 1.02 |
| one tile | Conservative | whole | 512 | 6 406 | 6 587 | 1.03 |
| one tile | NearestCell | chunked | 512 | 349 | 335 | 0.96 |
| one tile | NearestCell | whole | 512 | 344 | 333 | 0.97 |
| two-tile strip | Conservative | chunked | 16 384 | 3 147 | 3 082 | 0.98 |
| two-tile strip | Conservative | whole | 16 384 | 2 033 | 1 999 | 0.98 |
| two-tile strip | NearestCell | chunked | 16 384 | 701 | 690 | 0.98 |
| two-tile strip | NearestCell | whole | 16 384 | 375 | 350 | 0.93 |
| 2x2 | Conservative | chunked | 16 384 | 5 800 | 5 567 | 0.96 |
| **2x2** | **Conservative** | **whole** | **16 384** | **268 539** | **2 313** | **0.0086** |
| 2x2 | NearestCell | chunked | 16 384 | 2 332 | 2 335 | 1.00 |
| 2x2 | NearestCell | whole | 16 384 | 597 | 588 | 0.98 |

- The 2x2 whole-source conservative block is 116 times cheaper and allocates
  4 435 822 576 bytes before against 130 397 632 after; its 2 313 ns per cell sit
  between the strip's 1 999 and the single tile's 6 587, so the shape no longer
  decides the price.
- Every other arm is within five per cent in both directions, this benchmark's
  run-to-run spread. The chunked route never had the cliff — a chunk of a holding
  is one whole tile, which the block cursor could address — so those columns
  compare two trees over one rectangle and stay under the noise.

## 8. Locating a point

`cellat` on the 2x2 holding, over 32 sample sites through the grid's tree:

| | before | after |
| --- | ---: | ---: |
| per point | 943.1 ms | 0.0192 ms |
| allocated, 32 points | 23 316 538 000 B | 166 576 B |

- The generic cursor's tile node holds all 1 440 000 of a tile's pixels as
  children and derives a cap for each; the tiled tree bisects to nine. Building
  it costs 4 microseconds and 5 088 bytes for four tiles, against a constant 160
  for the wrapper it replaces.

## 9. Suites

| suite | before | after |
| --- | --- | --- |
| `test/runtests.jl` | 1 046 865 pass, 0 fail, 0 error, 17 broken | 1 047 316 pass, 0 fail, 0 error, 17 broken |
| the regridding bridge cross-check | 383 pass, 0 fail, 0 error, 12 broken | 383 pass, 0 fail, 0 error, 12 broken |

- The root suite runs in 8m12s with exactly the two pre-existing `globalindices`
  deprecation warnings and no others.
- Its 451 extra passes are all new Copernicus DEM assertions; that file alone is
  16 779 pass and 3 broken, net of the block-cursor rows that moved from holdings
  onto the complete lattice and its windows.
- They cover the index law and its round trip over five tile shapes; the explicit
  offset-plus-row-major law over 400 random pixels of a square and an L; leaf
  coverage with no duplicate and no outside index; every leaf cap against its
  pixel's boundary; cap nesting through the tile layer; point-query parity
  against the per-tile block cursors over 300+ probes per shape including cell
  edges and corners; the STI and `Trees` surface the regridder reads; type
  inference at every node touched and hook called; an empty holding; a level-0
  holding; and intersection matrices identical to the generic cursor's.
