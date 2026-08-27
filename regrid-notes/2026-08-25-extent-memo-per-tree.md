# One extent table per tree

- Date: 2026-08-25
- Scope: `src/engine/extent_memo.jl` (new), `src/engine/engine.jl`,
  `src/engine/tiled_raster.jl`, `src/systems/CopernicusDEM/cursor_memo.jl`,
  and the Copernicus DEM tests.
- Commit: `Keep one extent table per tree`

## 1. What the memo is for

Two cursors derive a node's extent instead of storing one:

- `MemoBlockCursor` over `BlockCursor` — the lazy tree over the complete
  Copernicus DEM lattice and the windows `subcursor` cuts out of it. A cap comes
  from `_node_box` through `_box_cap`: five `sincosd` pairs, four spherical
  distances.
- `TiledRasterCursor` — a bisection quadtree inside each raster tile. A raster
  node's cap comes from the grid's `raster_cap` hook. Its tile-layer caps are
  stored in the tree, not derived.

Both are asked for the same node many times — a dual-tree join asks each node
once per opposing node, and a tile's tree is rewalked for every source block,
destination column and worker — so both memoize in 1024 direct-mapped slots in
`task_local_storage`. A key hashes to one slot; a hit is a compare and a load, a
miss derives and overwrites. The slot holds the whole key, so a collision costs a
re-derive and never a wrong extent, and memory per task is bounded whatever the
lattice or tile size.

## 2. Why one table is the wrong number

- A node's extent is fixed by its rectangle **and its tree**: the same tile
  rectangle has one cap where cells are tiles and another where they are pixels
  (the leaf pad differs), and the same `(tile slot, rows, columns)` names
  different pixels in two holdings. One table must therefore clear all 1024 slots
  whenever the task turns to another tree.
- That turn is the common case: a **dual-tree join** alternates between two trees
  on one task (a CopDEM-to-CopDEM regrid puts two cursors of a kind on one task);
  a **window against its holding** joins the tree `subcursor` cut for one chunk
  against the whole holding's; a **self-join** reads one tree through two cursors.
- Under one table each alternation clears on every ask, so the memo never serves
  a hit: it pays a hash, a 1024-slot `fill!` and a box for the new identity, then
  derives anyway — a pessimisation exactly where it was meant to pay, worse the
  deeper the join. Measured: alternating asks cost 681 ns each on the raster side
  and 1069 ns on the lattice side, against 84 ns and 55 ns per tree, and 160 and
  16 bytes per ask against none.

## 3. The facility

`src/engine/extent_memo.jl` holds the slot, hash and compare logic once; both
cursors call it.

- **`ExtentTable{K,I}`** — one tree's `EXTENT_MEMO_SLOTS` (1024) slots. `K` is
  the node-key type: `NTuple{9,Int}` for a block node (tile and pixel rectangles
  plus the pixel flag), `NTuple{5,Int}` for a raster node (tile slot and
  rectangle). `I` is the tree-identity type. `misses` counts derivations since
  the table was keyed, written on the miss path only.
- **`ExtentMemo{K,I}`** — one task's tables for that type pair, at most
  `EXTENT_MEMO_TABLES` (4). Four covers a join's two trees and a per-chunk window
  tree read beside the holding it came from.
- **`extent_table(K, id)`** — the calling task's table for `id`, found by an
  identity (`===`) scan, so a turn costs at most four compares and leaves every
  keyed table intact. A miss past the count re-keys the **least recently used**
  table: the pattern that motivates eviction is one long-lived tree against a
  stream of short-lived ones, and a rotating victim would evict the long-lived
  tree every fourth window.
- **`memoized_extent(derive, table, key)`** — the slot arithmetic and compare,
  returning a `Cap`. The hit path writes only the recency stamp.

The tree identity is whatever pins the geometry, compared with `===` and nothing
else:

| cursor | identity | key |
| --- | --- | --- |
| `MemoBlockCursor` | `(grid, sys, level)` | `NTuple{9,Int}` |
| `TiledRasterCursor` | the `RasterTileTree` | `NTuple{5,Int}` |

Both identities are concrete, so the compare is inline and nothing is boxed for
it. `STI.node_extent` on both cursors returns `Cap`, is `@inferred`, and
allocates nothing on the hit path. Memory per task is four tables per `(K, I)`
pair — a block table is 112 KiB and a raster table 80 KiB, so at most 448 KiB and
320 KiB, and nothing until a task asks for an extent.

## 4. Behaviour

- Unchanged: every answer is the memo-free cursor's, bit for bit, whoever asks
  and in whatever order; a collision is still a miss that overwrites; tables are
  still task-local and lock-free.
- Changed: a task alternating between trees derives each distinct node once
  instead of once per ask, and a task touching more trees than it holds tables
  for re-keys one rather than answering wrongly.

## 5. Gate

Same machine, same session, arms interleaved, best of two. `before` is the
single-table design; `after` is this one.

### Copernicus DEM to Copernicus DEM

`Conservative` weights for a 16,384-cell destination window — a prefix of the
GLO-30 tile N46 E010 — against all four chunks of a GLO-90 2x2 holding (N46-N47,
E010-E011), on the `buildweights!` path the chunked route takes. Both sides are
`TiledRasterCursor`s. Real tiles from disk; 21,968 nonzero weights before and
after, 16 builds per timed run.

| arm | before | after | ratio |
| --- | --- | --- | --- |
| ns per destination cell, one thread | 5178.5 | 3598.9 | 0.695 |
| bytes per destination cell, one thread | 4712.6 | 4447.3 | 0.944 |
| ns per destination cell, `-t auto` | 1495.7 | 1369.1 | 0.915 |
| bytes per destination cell, `-t auto` | 6612.9 | 7743.8 | 1.171 |

The one-thread arm is the memo's own effect: one task asks every extent and sees
every repeat. Under `-t auto` the search spreads over worker tasks that see fewer
repeats, and the extra bytes are the tables — up to four per key type per task
instead of one, 80 KiB each allocated once per task, so a per-task constant
landing here on a fresh set of tasks sixteen times over.

Micro-arms: `node_extent` over two trees turn and turn about, 102,400 asks,
`-t auto`.

| arm | before | after | ratio |
| --- | --- | --- | --- |
| two `TiledRasterCursor`s, ns per ask | 681.0 | 83.8 | 0.123 |
| two `TiledRasterCursor`s, bytes per ask | 160 | 0 | — |
| two `MemoBlockCursor`s, ns per ask | 1069.2 | 55.3 | 0.052 |
| two `MemoBlockCursor`s, bytes per ask | 16 | 0 | — |

The bytes are the box the single table's `Any` identity field needed for the tree
it was being re-keyed to, once per ask.

### Copernicus DEM to IGeo7, unchanged

The previous gate's arms: `Conservative` and `NearestCell` on one tile, a
two-tile strip and a 2x2 holding, on the whole and chunked source routes, plus
the tree's own point query. Nanoseconds per destination cell, `-t auto`.

| arm | before | after | ratio |
| --- | --- | --- | --- |
| 1tile conservative chunked | 6335.8 | 6370.2 | 1.005 |
| 1tile conservative whole | 6341.5 | 6301.4 | 0.994 |
| 1tile nearest chunked | 355.0 | 332.4 | 0.936 |
| 1tile nearest whole | 352.9 | 339.8 | 0.963 |
| strip conservative chunked | 3105.4 | 3012.5 | 0.970 |
| strip conservative whole | 2036.1 | 1995.3 | 0.980 |
| strip nearest chunked | 737.4 | 690.1 | 0.936 |
| strip nearest whole | 360.2 | 345.8 | 0.960 |
| 2x2 conservative chunked | 5627.2 | 5526.7 | 0.982 |
| 2x2 conservative whole | 810052.1 | 810416.8 | 1.000 |
| 2x2 nearest chunked | 2341.0 | 2304.4 | 0.984 |
| 2x2 nearest whole | 600.2 | 580.3 | 0.967 |
| 2x2 point query, ns per point | 19600.0 | 18200.0 | 0.929 |

Every arm is within the benchmark's spread, none slower by more than half a per
cent. Allocated bytes are unchanged to three figures except `strip conservative
chunked`, which moves by five per cent.

## 6. Suites

- `test/systems/CopernicusDEM/runtests.jl` standalone: 16,800 pass, 0 fail,
  0 error, 3 broken (16,779 before; the difference is the tests added here).
- Root `test/runtests.jl`: 1,047,337 pass, 0 fail, 0 error, 17 broken, with the
  two pre-existing `globalindices` deprecation warnings and no others.
- Regrid cross-check: `XCHECK TOTAL passed=383 failed=0 errored=0 broken=12`.

## 7. Tests

In `test/systems/CopernicusDEM/runtests.jl`, testset `one extent table per tree`:

- **alternating derives each node once.** Two trees of one kind — two lattice
  cursors differing in level, two `TiledRasterCursor`s over different holdings —
  walked turn and turn about three times on a task with no tables of its own.
  Derivations after each pass must equal the distinct node count and must not
  grow; under one table it would be the ask count, six times larger. The lattice
  arm also asserts the two trees give some tile rectangle the same key, so the
  alternation is the one a shared table would clear on.
- **the same weights.** `intersection_areas` for a holding against a holding, and
  for a level-0 lattice window against a level-1 one, must equal the memo-free
  cursors' matrix — `HierarchicalGridCursor` for the first, bare `BlockCursor`s
  for the second.
- **eviction.** Five trees round-robined on one task, three passes: every answer
  is still the `raster_cap` of the tree the node came from, the store holds
  exactly `EXTENT_MEMO_TABLES` tables, and a re-keyed table's second ask hits.

The existing memo testsets stay: every answer bit-identical to the bare cursor's,
a real slot collision resolved as a miss, and eight tasks walking one grid's
nodes with no shared slot.
