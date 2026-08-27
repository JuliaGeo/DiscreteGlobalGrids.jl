# D1 — one prepared destination per tile

- Date: 2026-08-25
- Card: `regrid-notes/2026-08-21-regridding-simplification-plan.md`,
  "Phase 7 — prepared destination geometry", Task D1
- Commit: `Consolidate destination geometry preparation`
- Benchmarks: `benchmark/conservative_block_baseline.jl`,
  `benchmark/conservative_roundtrip_baseline.jl`

## 1. What a destination tile is now

- `DestinationCache` (`conservative.jl`) holds one tile's index set, chunk-local
  index map, restricted tree, one polygon slot per tile row, and a fill flag per
  row under one lock. Not a `RegridSpace`: `plan.dst_space` still answers.
- `preparedestination(method, dst_space, dst_inds, budget)` answers the index set
  or a cache, and that answer *is* the destination a build addresses:

```julia
pairblock(::Conservative, dst_space, dst_inds,   src_space, src_inds)
pairblock(::Conservative, dst_space, dst_cache,  src_space, src_inds)
```

- The first prepares nothing: `CellMemo` on both sides, its own restricted tree.
  The second reuses the tile's tree, reads polygons by the block row already
  computed, and memoizes the source alone. Same values, CSC layout and
  denominator order, asserted entry for entry.
- Slots fill from the candidate pairs, never from a prediction:
  `work_items(op, pairs)` runs once per `intersection_areas`, after the dual
  search and before partitioning or spawning, filling unfilled rows under the
  lock and returning the list unchanged.
- Hence: a cell no source overlaps is never synthesized, none twice, serial and
  threaded assembly leave the same cache, and each block extends the same array.
- The clipping loop does one array index and holds no lock. A block reads after
  releasing the fill lock, so another block's slot is visible; rows are distinct
  memory and the array never resizes.
- The fill flag is a `Bool` per row, not `isassigned`: an isbits polygon is
  stored inline and leaves no empty slot. An IGeo7 cell is 184 bytes, a
  `SmallList` ring inside the polygon.
- The tile's tree is opaque to inference (untyped field, `_destinationtree`'s
  barrier, one at the executor's choice), so a block's assembly compiles for the
  tree type the run reaches rather than every type the space can produce.
- Deleted: `TileCells <: RegridSpace` and its fourteen forwards, its `subtree`,
  `_tiletree!`, `_synthesizecells`, `_cachedtree`, `CachedCellTree` with its nine
  spatial-tree forwards and `_task_prepared_raster_tree`, `_TILE_CELL_CACHE_MAX`.
- `blockfor`/`buildblock` take the prepared destination where a wrapper space
  went; `tilefor` drops the argument, the point-method tile route reading sample
  sites and no polygon.
- Untouched: `chunkextents` and `CapCachedTree` cache extents, not polygons.

## 2. The spike, and the default it justifies

Same session, three repetitions per arm, Julia 1.12.6, 8 threads, idle machine.

### Cost of one destination cell

A `@belapsed` sweep of every cell divided by the cell count, so a synthesized
polygon's allocation is included; the prepared column sweeps a slot array.

| destination | one `getcell` | through the tree | one prepared slot | polygon bytes |
|---|---:|---:|---:|---:|
| IGeo7 level 5 (168,072 cells) | 421.5 ns | 427.3 ns | 0.5 ns | 184 |
| `RasterGrid` 360×180 (64,800 cells) | 27.9 ns | 27.0 ns | 0.6 ns | 192 |

A DGGS cell costs fifteen times a raster cell (four chart evaluations off two
axes). A prepared slot is free next to either.

### How many syntheses each way

Counted exactly by replaying the 64-slot memo's policy over the candidate list.
The memo beats one-synthesis-per-pair because the candidate order repeats a
destination within a source-leaf pairing.

| workload | destination cells | candidate pairs | memo syntheses | prepared syntheses |
|---|---:|---:|---:|---:|
| IGeo7 level 5 ← 360×180 | 168,072 | 732,772 | 192,996 | 168,072 |
| 360×180 ← 3600×1800 | 64,800 | 33,285,272 | 102,720 | 64,800 |

Within one block preparing saves 13 % of the syntheses on the first and 37 % on
the second, and pays a fill pass over the whole candidate list. The saving that
matters is across blocks: every block starts with a cold memo.

### Time and peak memory, prepared against not

`whole` is one block over the whole domain; `tile` is one destination tile
against every source chunk, the shape the lazy executor builds. Builds are warmed
on a small workload of the same types, so these are not first-call latency.

| destination | shape | prepared | seconds, median (min–max) | max RSS, mean |
|---|---|---|---:|---:|
| IGeo7 | whole (1 block) | no | 0.195 (0.188–0.196) | 941 MiB |
| IGeo7 | whole (1 block) | yes | 0.245 (0.236–0.248) | 981 MiB |
| IGeo7 | tile (16 blocks) | no | 0.374 (0.369–0.378) | 904 MiB |
| IGeo7 | tile (16 blocks) | yes | 0.244 (0.242–0.249) | 913 MiB |
| raster | whole (1 block) | no | 10.957 (10.555–11.640) | 2,880 MiB |
| raster | whole (1 block) | yes | 11.527 (11.410–11.834) | 2,868 MiB |
| raster | tile (32 blocks) | no | 11.312 (11.008–11.355) | 1,076 MiB |
| raster | tile (32 blocks) | yes | 11.583 (11.366–12.062) | 1,073 MiB |

- **Preparing pays only where more than one block reads it**: 35 % faster over a
  DGGS tile's sixteen blocks (0.374 s → 0.244 s) and 26 % slower on one block
  (0.195 s → 0.245 s), ranges apart in both.
- One block's memo working set is 64 polygons, 12 KiB and cache-resident; a
  prepared tile is a 30 MiB array read in candidate order, costing more traffic
  than the 25,000 syntheses it saves.
- So a one-block build prepares nothing — `wholeblock`, a bare `buildblock`, and
  a tile with one source chunk. The executor prepares once per tile where it has
  several, before the tile's blocks start.
- **A cheap destination never pays**: on a raster destination preparing is 2–5 %
  slower with overlapping ranges in both shapes, and adds 12 MiB of resident
  polygons and a fill pass over 33 million pairs. Raster destinations prepare
  nothing.

## 3. The trait: measured, not inherited

| space | index set | tree | `node_extent_is_expensive` |
|---|---|---|---|
| `DGGSpace` (IGeo7) | whole space | `CapCachedTree` | `true` |
| `DGGSpace` (IGeo7) | a range | `CellSpaceRTree` | `false` |
| `RasterGrid` | whole space | `TopDownQuadtreeCursor` | `true` |
| `RasterGrid` | a lattice rectangle | `TopDownQuadtreeCursor` | `true` |
| `RasterGrid` | scattered indices | `CellSpaceRTree` | `false` |

- `STI.node_extent_is_expensive` was the obvious candidate and classifies the
  wrong thing: `true` for both production destination trees and `false` only for
  the packed R-tree fallback, which both spaces reach depending on the index set.
- It is a property of the tree an index set produces, and "extents derived or
  stored?" is orthogonal to "is a cell polygon expensive to synthesize?".
- So the space property is stated directly: `expensivecellgeometry(space)`,
  `true` by default, `false` on `RasterGrid`; documented, public, qualified.
- `preparesdestination(method, dst_space)` is the method side: `false` by
  default, `expensivecellgeometry(dst_space)` for `Conservative`.

## 4. The budget rule, and why the tiling did not change

- `DEST_BUDGET_SHARE` still reserves one eighth of the budget for a destination
  tile's per-cell state, now named `destcellbudget(budget)`.
- `destcellbytes(method, dst_space, ndst)` is what one cell holds: the executor's
  `DEST_BYTES_PER_CELL` accumulators plus one prepared polygon and its flag where
  kept, from one `Base.summarysize` probe.
- `destcellsfit(method, dst_space, ndst, budget)` is the one rule and replaces
  `_TILE_CELL_CACHE_MAX`: a tile prepares when its whole per-cell state fits
  `destcellbudget`. A single-block build prepares nothing whatever the budget.
- The tiling charges only the executor's own per-cell cost; `_defaulttilesizes`
  is arithmetically unchanged (`fld(budget, 8 * 40)` became
  `fld(destcellbudget(budget), 40)`, the same integer).
- Charging prepared polygons in the tiling was measured and rejected: at 225
  bytes per IGeo7 cell it divides the tile size by 5.6, costing one more
  restricted tree and one more read of every shared source chunk, and it moves
  the cliff `2026-08-23-dggs-subcursor-scope.md` documents — 823,543 cells at
  `budget = 2^30` is one run today, two if polygons are charged, dropping it from
  an `O(1)` hierarchy cursor onto a packed R-tree over 412 k cells.

| workload | destination | chunks | tiles before | tiles after | tile cells | per-cell bytes |
|---|---|---:|---:|---:|---|---:|
| block | IGeo7 level 5 | 72 | 72 | 72 | 2,001 / 2,401 | 225 |
| roundtrip | `RasterGrid` 360×180 | 9 | 9 | 9 | 7,200 | 40 |

Tile counts and every other tiling are unchanged: the same function of the same
inputs, both destinations space-tiled before the budget is reached. The largest
block tile holds 2,401 × 225 = 528 KiB against a 256 MiB share, so a lazy tile of
that shape prepares; the round-trip's raster destination prepares nothing and is
charged the executor's 40 bytes alone.

## 5. What it verifies

`test_conservative.jl` and `test_proj.jl`, 56 assertions net (4,657 → 4,713
passing, one broken unchanged):

- **prepared and unprepared are the same weights**, raster and cell-space: equal
  matrices, `colptr`, `rowvals`, `nzrange` per column, `===` on every value and
  denominator, and the same against the generic route. An empty source keeps the
  index form's degenerate contract.
- **one restricted tree per tile, one synthesis per prepared row**: a tile built
  against two source chunks builds one tree and advances the synthesis counter by
  exactly the filled-row count over both blocks; rows neither chunk reaches are
  never synthesized. Preparing per block builds two trees with identical weights.
- **concurrent blocks share one cache**: four blocks under
  `OUTER_PARALLEL => true` build one tree, fill each row once, and each equals
  its serial reference.
- **the budget rule**: the charge per destination kind and method; that tile size
  uses the executor's cost alone; that a tile at the cap is refused and a smaller
  one granted; and that a cheap destination, a point method, too small a budget
  and an empty tile each answer with the index set.
- **a method that reads no prepared geometry still builds**: it names the tile's
  cells through the space, fills no slot, and produces the index-form block.
- **nothing added to the clipping loop**: one operator call against a prepared
  destination allocates zero bytes, and no more than against a memo.
- `test_lazy.jl`'s "budget tiles build one restricted tree each" now proves the
  executor prepares: unprepared, each tile's blocks would build their own trees.
- `test_proj.jl`: a prepared destination's stored tree retains the safe chart
  across block tasks, and a serial build still derives a private prepared cursor.

Mutants killed: differing values, layout or denominators; a fill that predicts
rows (untouched-row count); a fill per block or per task (tree and synthesis
counts); a cache shared without a lock (concurrent test); a prepared destination
handed to a method that cannot read one; a clipping loop synthesizing behind the
slot (allocation floor); a budget rule replaced by a constant; a tiling that
charges polygons.

Dropped with the wrappers they tested: "one tile's cell geometry, synthesized
once", "the destination tile cell cache is bounded", and "one tile's restricted
tree, built once", whose claim the counting test makes more precisely.

## 6. Gate

Two states interleaved rep by rep out of two checkouts, three repetitions per
arm, one session, Julia 1.12.6, 8 threads, machine idle of other benchmarks
before each arm (an earlier attempt was discarded for an overlap). Both states
produce identical nonzero counts and weight sums. Seconds are the median, peak
RSS the mean.

| workload / mode | seconds before | seconds after | ratio | RSS before | RSS after | ratio |
|---|---:|---:|---:|---:|---:|---:|
| block / direct | 4.881 | 2.625 | 0.538 | 866.0 MiB | 809.8 MiB | 0.935 |
| block / chunked | 2.527 | 2.548 | 1.008 | 844.6 MiB | 819.0 MiB | 0.970 |
| roundtrip / direct | 11.289 | 11.776 | 1.043 | 2,716.0 MiB | 2,794.2 MiB | 1.029 |
| roundtrip / chunked | 11.856 | 12.104 | 1.021 | 2,594.8 MiB | 2,583.8 MiB | 0.996 |

- Neither benchmark warms its build, so each arm times a cold call: the build
  plus compiling it.
- The round-trip's arms move by 2 % and 4 % against within-state ranges spanning
  11.2–12.4 s, and peak RSS by at most 3 % against this workload's usual
  200–320 MiB spread.
- The block workload's chunked arm is level and its eager arm 46 % faster: that
  route used to compile its assembly for every tree type the destination space
  can produce, and §1's barrier leaves it compiling the one it reaches.

Timing the same build three times in one process separates build from compile:

| state / mode | first call | second | third | first allocated | second allocated |
|---|---:|---:|---:|---:|---:|
| before / direct | 5.453 s | 0.123 s | 0.128 s | 1,424,871,584 | 255,315,040 |
| before / chunked | 3.240 s | 0.126 s | 0.126 s | 1,329,831,584 | 255,315,536 |
| after / direct | 3.449 s | 0.136 s | 0.139 s | 1,392,899,440 | 255,315,040 |
| after / chunked | 3.393 s | 0.128 s | 0.124 s | 1,451,921,888 | 255,315,136 |

The steady-state build is 0.123–0.139 s in all four states, allocating the same
bytes exactly on the eager route and within 400 bytes of 255 MB on the chunked
one. First calls hold no surprise: the eager route drops 2.0 s, the chunked route
is level.

## 7. Suites

- `lib/GlobalRegridding/test/runtests.jl`: 4,713 passed, 0 failed, 0 errored, 1
  broken, zero deprecation warnings, against a baseline of 4,657 / 0 / 0 / 1
  measured this session on the untouched checkout. The 56 new assertions are §5's
  six testsets (81) less the three that tested the deleted wrappers (25).
- Bridge cross-check (`test/systems/crosssystem/regridding_conservation.jl`,
  `regrid.jl`, `regrid_acceptance.jl`): `XCHECK TOTAL passed=380 failed=0
  errored=0 broken=12`, unchanged, zero deprecation warnings. This is where
  `regrid.jl` runs; it is not a standalone entry point.
- Eager and chunked weights and values for `Conservative`, `NearestCell` and
  `BilinearPoint` over the lib fixtures: 54 serialized entries, all `==` and all
  `isequal` to the untouched checkout's reference.
