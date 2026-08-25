# The tiled raster cursor

- Date: 2026-08-25
- Scope: `src/engine/tiled_raster.jl` (new), `src/interface/grid.jl`,
  `src/systems/CopernicusDEM/cursor.jl`, and the Copernicus DEM tests.
- Commit: `Descend a quadtree inside each raster tile`

## 1. What a grid of raster tiles is

Some grids are not a lattice and are not a hierarchy of parents and children.
They are **a collection of raster tiles**: a set of rectangles, each a grid of
pixels, whose pixels together are the grid's cells. A Copernicus DEM holding is
one — a set of 1x1-degree tiles, each 1200x1200 or 3600x3600 pixels. So is a
STAC catalogue: a list of items, each a raster with its own footprint, size and
projection, overlapping its neighbours or not.

Such a grid has no single lattice to bisect and no hierarchy to descend, so
neither of the trees this package had suited it. It has instead exactly two
structures: a flat set of tiles, each with a footprint, and a rectangle of
pixels inside each tile. The tree built here is those two structures, in that
order — a packed tree over the tiles' caps, and a bisection quadtree inside each
tile — and it is built from an interface of four hooks that name nothing else.

## 2. The interface

A grid opts in by answering `raster_tiles`; the other three hooks are then
required. All four are declared in `src/interface/grid.jl`.

| hook | contract |
| --- | --- |
| `raster_tiles(grid, inds) -> tiles or nothing` | The tiles covering grid indices `inds`, as an indexable collection of opaque **tile handles**; `nothing` (the default) when this grid's cells are not raster pixels. The tiles must hold every cell of `inds` exactly once, and no cell outside it. |
| `raster_shape(grid, tile) -> (nrows, ncols)` | The tile's rectangle in pixels. Both positive. Rows and columns are numbered `0:nrows-1` and `0:ncols-1` everywhere else. |
| `raster_localindex(grid, tile, j, i) -> Int` | The index **in `grid`** of the pixel at row `j`, column `i`. It is the tile's offset in the grid plus the pixel's row-major index within the tile, so a tile's pixels are one contiguous block of grid indices in row-major order. |
| `raster_cap(grid, tile, j0, j1, i0, i1) -> SphericalCap` | A cap containing the geometry — boundaries, not merely centres — of every pixel in the closed sub-rectangle `j0:j1` x `i0:i1`. |

A tile handle is opaque to the engine: it is produced by `raster_tiles` and
passed back to the other three, and the engine reads nothing out of it. A grid
therefore carries in the handle whatever its geometry hooks want — an
identifier, the rectangle's origin, its offset in the grid — without decoding
anything twice. The engine names tiles by nothing but their caps and their
shapes, so **no lattice coordinate reaches it**.

`raster_cap` is asked for whole tiles once at build time, for interior
rectangles during descent, and for single pixels at the leaves. A grid with a
closed-form box for any sub-rectangle — Copernicus DEM has one, from its
latitude bands — answers in constant time; a grid without one merges the
sub-rectangle's per-pixel caps, which is correct but linear in the rectangle,
and wants a memo of its own above the engine's.

## 3. The tree

### The tile layer

`RasterTileTree` sorts the tiles by the Morton key of their cap centres and
splits them `RASTER_TILE_ARITY`-ways per level down to one tile per node,
storing each node's cap (the merge of its children's) and pixel count. This is
the packing `IndexTree` uses for cells, reused rather than re-derived: the same
`_morton_key` on the cap centre's longitude and latitude, the same 4-ary split,
the same stored nesting extents.

It is `IndexTree`'s packing and not `FlexibleRTrees.RTree` because every node
extent on this side of the boundary is a `SphericalCap`, and the R-tree owns
Cartesian `Extent`s — which is why `GlobalRegridding`'s `CellSpaceRTree` has to
convert one to the other at every node. A packed spherical-cap tree needs no
conversion, and its merge, `Extents.union` on two caps, is the operation the
DGG side already uses for node extents.

The tile layer assumes nothing about the tiles' arrangement. Copernicus DEM's
tiles happen to form a lattice; the root neither knows nor exploits that. A
STAC catalogue's arbitrary footprints pack the same way.

### The pixel layer

A node holding one tile *is* that tile's whole raster, so the tree spends no
level on the tile itself: descending into it bisects the rectangle — the longer
axis in two, near-equal parts — until a node holds `LEAF_CELLS` pixels or
fewer. That is exactly `BlockCursor`'s descent, and it is the same code:
`bisect_parts`, `rect_part`, `LEAF_CELLS` and the `LeafCells` leaf container now
live in the engine, and `BlockCursor` calls them.

### Reuse or generalise

`BlockCursor` bisects a rectangle of a lattice and names a cell by
`index = id - origin`. That law is what stopped it being the per-tile subtree:
it holds only where the rectangle's cells are one consecutive id run, which a
tile set crossing a latitude row is not.

So the **law was generalised, not the quadtree duplicated**. The engine's raster
node splits by the same rule, partitions by the same helper, and fills the same
leaf container; only the index law differs, and `raster_localindex` — a tile's
offset plus the pixel's row-major index within it — subsumes `id - origin`,
which is the special case of one tile whose offset is `origin` and whose
row-major order is the id order.

`BlockCursor` survives for the complete lattice and for the windows `subcursor`
cuts out of it, where it is strictly better: a lattice-wide tiled tree would
need a table of 64 800 tile caps, where the lattice derives any node's box in
O(1) from two divisions. Its `PartialGrid` constructor is gone, so no shape has
two trees: a complete level and its windows get `BlockCursor`, a holding gets
`TiledRasterCursor`, and nothing gets both.

### The extents

Tile node extents are stored in the tree and nest, because each is an
`Extents.union` of its children's. Raster node extents are derived from
`raster_cap` and memoized per task by `RasterExtentMemo` — the direct-mapped,
1024-slot, whole-key-compared table `MemoBlockCursor` uses, keyed on
`(tile slot, j0, j1, i0, i1)` and cleared when a task turns to another tree.
`node_extent_is_expensive` is therefore `false`: a hit is a key compare and a
load.

Node caps bound cell **geometry**, as everywhere else here: a node's cap covers
the boundaries of the pixels beneath it, not necessarily those pixels' own caps.

## 4. `treeify` and `subcursor`

`treeify(::PartialGrid{<:CopernicusDEMSystem})` builds the tiled raster tree for
any holding — one tile, a row, a square, an L, four tiles on four continents.
The `HierarchicalGridCursor` fallback remains reachable only when the grid names
an id the lattice does not have, which `raster_tiles` answers with `nothing`.

`subcursor(::PartialGrid{<:CopernicusDEMSystem}, inds)` builds the tree over
`inds` alone. Leaf indices are the grid's own at every node, so a window's tree
names its cells exactly as the whole grid's does, which is `subcursor`'s
contract.

Copernicus DEM's `raster_tiles` walks the holding's ids by runs and cuts each
run into at most three rectangles — a part row, the whole rows under it, a part
row — so each rectangle is one contiguous block of grid indices and the
row-major law holds on it. The run itself is found by binary search: ids ascend
strictly, so `id(p+k) - id(p) - k` never decreases. A level-0 holding, where the
cell is the tile, is a set of 1x1 rectangles.

## 5. What a STAC-backed system implements

Four methods and one handle type:

- a handle carrying the item's identifier, its raster size and its offset in the
  grid — the grid's cells are the catalogue's pixels in item order;
- `raster_tiles`: the items whose pixel blocks meet `inds`, in that order,
  taken from the catalogue rather than from any lattice;
- `raster_shape`: the item's raster dimensions;
- `raster_localindex`: the item's offset plus `j * ncols + i + 1`;
- `raster_cap`: the sub-rectangle's corners through the item's own geotransform
  and CRS, as a cap padded for the bow of its edges the way the Copernicus
  lattice's is — the contract is that the cap covers the pixels' boundaries,
  and a cap fitted to four corners alone does not.

Nothing else changes. Overlapping footprints, per-item projections and
heterogeneous tile sizes are already what the tile layer packs, because it packs
caps.

## 6. The root's shape

The alternative to a packed root is a flat one: a single node with every tile as
a child. Measured on a scaled twin lattice whose tiles are 30x30 pixels, over
2000 random cell centroids per point query and a dual-tree join against a
complete HEALPix level, best of two runs:

| tiles | cells | root | ns per point query | dual-tree join |
| --- | --- | --- | --- | --- |
| 2x2 (4 tiles) | 3 600 | packed | 3 999 | 0.0766 s |
| 2x2 (4 tiles) | 3 600 | flat | 4 068 | 0.0717 s |
| 4x4 (16 tiles) | 14 400 | packed | 5 661 | 0.0825 s |
| 4x4 (16 tiles) | 14 400 | flat | 5 594 | 0.0825 s |
| 32x32 (1024 tiles) | 643 200 | packed | 9 028 | - |
| 32x32 (1024 tiles) | 643 200 | flat | 30 691 | - |

The two roots' joins hold the same number of stored entries at both sizes, and
they cost the same at the tile counts a holding of tens has: a four-tile packed root
is one level over its tiles and tests about as many caps as a flat root does.

They diverge with the tile count, because a flat root tests every tile on every
query while a packed one discards blocks of them. At 1024 tiles the packed root
answers a point query 3.4 times faster, and the ratio grows from there. A global
Copernicus holding is 26 450 tiles, so the root shape is the difference between
26 450 cap tests per query and a few dozen.

## 7. The gate

Source: Copernicus DEM GLO-90 tiles. Destination: the IGeo7 level-12 cells
covering the same box, chunked at their level-5 ancestors, taking the chunk with
the most candidate source chunks. Two routes to the source tree: **chunked**,
where the pair is one source chunk's indices, and **whole**, where it is the
whole source space and the plan reaches `treeify(grid)` itself. Nanoseconds per
destination cell, best of two runs, before and after interleaved on an otherwise
idle machine:

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

The 2x2 whole-source conservative block is 116 times cheaper and allocates
4 435 822 576 bytes before against 130 397 632 after, and its 2 313 nanoseconds
per cell sit between the strip's 1 999 and the single tile's 6 587 — the shape
no longer decides the price. Every other arm is within five per cent of what it
cost before, in both directions, which is this benchmark's run-to-run spread.

The chunked route never had the cliff: a chunk of a holding is one whole tile,
which the block cursor could address, so those columns compare two trees over
the same rectangle and measure the difference between arithmetic on a lattice
and two loads from a tile table. It does not show up above the noise.

## 8. Locating a point

`cellat` on the 2x2 holding, driven through the grid's tree over 32 sample
sites:

| | before | after |
| --- | ---: | ---: |
| per point | 943.1 ms | 0.0192 ms |
| allocated, 32 points | 23 316 538 000 B | 166 576 B |

The generic cursor's tile node holds all 1 440 000 of a tile's pixels as
children and derives a cap for each; the tiled tree bisects to nine. Building
the tree itself costs 4 microseconds and 5 088 bytes for four tiles, against a
constant 160 for the wrapper it replaces.

## 9. Suites

| suite | before | after |
| --- | --- | --- |
| `test/runtests.jl` | 1 046 865 pass, 0 fail, 0 error, 17 broken | 1 047 316 pass, 0 fail, 0 error, 17 broken |
| the regridding bridge cross-check | 383 pass, 0 fail, 0 error, 12 broken | 383 pass, 0 fail, 0 error, 12 broken |

The root suite runs in 8m12s with exactly the two pre-existing `globalindices`
deprecation warnings and no others. Its 451 extra passes are all new assertions
on the Copernicus DEM side, whose file alone — the root suite includes it — is
16 779 pass and 3 broken, net of the block-cursor rows that moved from holdings
onto the complete lattice and its windows: the index law and its round trip
over five tile shapes, the explicit offset-plus-row-major law over 400 random
pixels of a square and an L, leaf coverage with no duplicate and no outside
index, every leaf cap against its pixel's boundary, cap nesting through the
tile layer, point-query parity against the per-tile block cursors over more
than 300 probes per shape including cell edges and corners, the STI and `Trees`
surface the regridder reads, type inference at every node the search touches
and on every hook it calls, an empty holding, a level-0 holding, and
intersection matrices identical to the generic cursor's for a square and an L.
