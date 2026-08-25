# D1 — one prepared destination per tile

- Date: 2026-08-25
- Card: `regrid-notes/2026-08-21-regridding-simplification-plan.md`,
  "Phase 7 — prepared destination geometry", Task D1
- Commit: `Consolidate destination geometry preparation`
- Benchmarks: `benchmark/conservative_block_baseline.jl`,
  `benchmark/conservative_roundtrip_baseline.jl`

## 1. What a destination tile is now

`DestinationCache` (`lib/GlobalRegridding/src/conservative.jl`) holds one
tile's index set, its chunk-local index map, its restricted tree, one polygon
slot per tile row, and a fill flag per row under one lock. It is not a
`RegridSpace`, forwards nothing, and answers no question about the
destination: `plan.dst_space` still does.

`preparedestination(method, dst_space, dst_inds, budget)` answers either the
index set or a cache, and that answer *is* the destination a build addresses:

```julia
pairblock(::Conservative, dst_space, dst_inds,   src_space, src_inds)
pairblock(::Conservative, dst_space, dst_cache,  src_space, src_inds)
```

The first is the builder that prepares nothing. It keeps the task-local
`CellMemo` on both sides and builds its own restricted tree, exactly as
before. The second reuses the tile's tree, takes its destination polygons from
the cache by the block row it has already computed for the matrix, and keeps
the memo for the source alone. Both produce the same weights — values, CSC
layout, and denominator accumulation order — and the tests assert that entry
for entry.

Slots fill from the candidate pairs, never from a prediction.
`ConservativeRegridding.work_items(op, candidate_pairs)` runs once per
`intersection_areas` call, after the dual search has produced the whole pair
list and before the list is partitioned or any assembly task is spawned. The
operator's `work_items` walks that list under the cache lock, synthesizes the
polygon of every destination row not yet filled, and returns the list
unchanged. So a destination cell no source overlaps is never synthesized, no
cell is synthesized twice, serial and threaded assembly leave the same cache,
and each block of the tile extends the same array. The clipping loop then does
one array index and holds no lock.

A block's own reads happen after it releases the fill lock, so a slot another
block filled is visible to it; distinct rows are distinct memory, and the
array never resizes.

The fill flag is a `Bool` per row rather than `isassigned` on an `undef`
vector, because an isbits polygon is stored inline and leaves no empty slot to
test. An IGeo7 cell is exactly that: 184 bytes, a `SmallList` ring stored in
the polygon itself.

The tile's tree is opaque to inference on purpose. `DestinationCache` holds
it in an untyped field, `_destinationtree` puts an unprepared destination's
tree behind the same barrier, and the executor's choice between prepared and
unprepared crosses one too. `subtree` answers one of several tree types for a
space, and a block's assembly is compiled for the tree it is handed, so
inferring the answer compiles that assembly for every type the space can
produce — latency every fresh process pays, whatever it goes on to build. One
dynamic dispatch per block compiles the type the run reaches instead. The
clipping loop sits inside that specialization and pays nothing for it.

Deleted: `TileCells <: RegridSpace` and its fourteen forwarding methods, its
`subtree`, `_tiletree!`, `_synthesizecells`, `_cachedtree`, `CachedCellTree`
with its nine spatial-tree forwards and its `_task_prepared_raster_tree`
method, and `_TILE_CELL_CACHE_MAX`. `plans.jl`'s `blockfor` and `buildblock`
take the prepared destination in the argument that used to take a wrapper
space; `tilefor` drops that argument entirely, because the point-method tile
route reads destination *sample sites* and no cell polygon. Nothing on that
route touches the cache.

The raster and DGG chunk-extent caches (`chunkextents`, `CapCachedTree`) are
untouched. They cache extents, not polygons.

## 2. The spike, and the default it justifies

Same session, three repetitions of each arm, Julia 1.12.6, 8 threads, machine
otherwise idle.

### Cost of one destination cell

Each figure is a `@belapsed` sweep of every cell of the space divided by the
cell count, so it includes the allocation a synthesized polygon makes; the
prepared column sweeps the same polygons out of a slot array instead.

| destination | one `getcell` | through the tree | one prepared slot | polygon bytes |
|---|---:|---:|---:|---:|
| IGeo7 level 5 (168,072 cells) | 421.5 ns | 427.3 ns | 0.5 ns | 184 |
| `RasterGrid` 360×180 (64,800 cells) | 27.9 ns | 27.0 ns | 0.6 ns | 192 |

A DGGS cell is fifteen times the cost of a raster cell, which is four chart
evaluations off two coordinate axes. A prepared slot is free next to either.

### How many syntheses each way

The 64-slot direct-mapped memo is far better than "one synthesis per candidate
pair": the candidate order repeats a destination within a source-leaf pairing,
which is what the memo is shaped for. Counted exactly, by replaying the memo's
policy over the candidate list:

| workload | destination cells | candidate pairs | memo syntheses | prepared syntheses |
|---|---:|---:|---:|---:|
| IGeo7 level 5 ← 360×180 | 168,072 | 732,772 | 192,996 | 168,072 |
| 360×180 ← 3600×1800 | 64,800 | 33,285,272 | 102,720 | 64,800 |

Within *one* block, therefore, preparing saves 13 % of the destination
syntheses on the first and 37 % on the second — and pays a fill pass over the
whole candidate list to get them. The saving that matters is across blocks:
every block of a tile starts with a cold memo, so an unprepared tile
re-synthesizes its destinations once per source chunk.

### Time and peak memory, prepared against not

`whole` is one block over the whole domain; `tile` is one destination tile
against every source chunk, which is the shape the lazy executor builds.

Each build is warmed on a small workload of the same types first, so these are
build times and not first-call latency.

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

Two facts decide the default.

**Preparing pays only where more than one block reads it.** On a DGGS
destination it is 35 % faster across a tile's sixteen blocks (0.374 s → 0.244 s,
ranges apart) and 26 % *slower* on a single block (0.195 s → 0.245 s, ranges
apart). A single block's memo working set is 64 polygons — 12 KiB, resident in
cache — while a prepared tile is a 30 MiB array read in candidate order, and
that traffic costs more than the 25,000 syntheses it saves. So a build of one
block prepares nothing: the eager whole domain (`wholeblock`), a bare
`buildblock`, and a tile the executor finds has exactly one source chunk. The
executor prepares where it has several, once per tile, before the tile's blocks
start.

**A cheap destination never pays.** On a raster destination preparing is 2–5 %
slower with overlapping ranges in both shapes, and it adds 12 MiB of resident
polygons and a fill pass over 33 million candidate pairs. Raster destinations
prepare nothing.

## 3. The trait: measured, not inherited

The spatial-tree trait `STI.node_extent_is_expensive` was the obvious
candidate. It does not classify what the measurements classify:

| space | index set | tree | `node_extent_is_expensive` |
|---|---|---|---|
| `DGGSpace` (IGeo7) | whole space | `CapCachedTree` | `true` |
| `DGGSpace` (IGeo7) | a range | `CellSpaceRTree` | `false` |
| `RasterGrid` | whole space | `TopDownQuadtreeCursor` | `true` |
| `RasterGrid` | a lattice rectangle | `TopDownQuadtreeCursor` | `true` |
| `RasterGrid` | scattered indices | `CellSpaceRTree` | `false` |

It answers `true` for both production destination trees and `false` only for
the packed R-tree fallback — which both spaces reach, depending on the index
set. It is a property of the tree an index set produces, not of the space, and
what it describes (are extents derived or stored?) is orthogonal to what
decides this (is a cell polygon expensive to synthesize?): the R-tree stores
extents *because* it derived them from cells it had to build.

So GlobalRegridding states the space property directly.
`expensivecellgeometry(space)` defaults to `true` — a space that derives cell
boundaries from an index pays for every call — and `RasterGrid` overrides it
to `false`. It is a documented, public, qualified space hook, so a space whose
cells are cheap can say so.

`preparesdestination(method, dst_space)` is the method side: `false` by
default, and `expensivecellgeometry(dst_space)` for `Conservative`. An area
method reads a destination polygon once per overlapping source leaf; a point
method evaluates at a sample site and reads none.

## 4. The budget rule, and why the tiling did not change

`DEST_BUDGET_SHARE` still reserves one eighth of the budget for one
destination tile's per-cell state, now named: `destcellbudget(budget)`.
`destcellbytes(method, dst_space, ndst)` is what one destination cell holds
inside it — the executor's `DEST_BYTES_PER_CELL` accumulators, plus one
prepared polygon and its flag where the build keeps them, from one
`Base.summarysize` probe.

`destcellsfit(method, dst_space, ndst, budget)` is the one rule, and it
replaces `_TILE_CELL_CACHE_MAX`: a tile prepares its geometry when the tile's
whole per-cell state fits `destcellbudget`. It is what refuses a lazy tile
whose polygons the budget cannot hold; a single-block build prepares nothing
whatever the budget says. The constant is gone, and nothing is capped by a
number that does not name a budget.

The tiling itself charges only the executor's own per-cell cost, exactly as
before, and `_defaulttilesizes` is arithmetically unchanged
(`fld(budget, 8 * 40)` became `fld(destcellbudget(budget), 40)`, which is the
same integer). Charging prepared polygons there instead was measured and
rejected: at 225 bytes per IGeo7 cell it divides the budget-derived tile size
by 5.6, and a smaller tile costs one more restricted tree and one more read of
every source chunk it shares — both larger than the syntheses preparing saves.
It also moves the cliff `regrid-notes/2026-08-23-dggs-subcursor-scope.md`
documents: a whole-space DGGS destination of 823,543 cells at `budget = 2^30`
is one run today, and charging polygons would split it into two, which drops
it from an `O(1)` hierarchy cursor onto a packed R-tree over 412 k cells.
Preparing takes what is left of the share instead, and is refused where there
is nothing left.

Tile counts on both benchmark workloads are therefore unchanged, and so is
every other tiling: the destination tiling is the same function of the same
inputs. Both workloads' destinations are space-tiled in any case, which the
budget never reaches:

| workload | destination | chunks | tiles before | tiles after | tile cells | per-cell bytes |
|---|---|---:|---:|---:|---|---:|
| block | IGeo7 level 5 | 72 | 72 | 72 | 2,001 / 2,401 | 225 |
| roundtrip | `RasterGrid` 360×180 | 9 | 9 | 9 | 7,200 | 40 |

The largest block tile holds 2,401 × 225 = 528 KiB of per-cell state against
a 256 MiB share, so a lazy tile of that shape prepares. The round-trip's
destination is a raster and prepares nothing, so its cells are charged the
executor's 40 bytes alone.

## 5. What it verifies

`lib/GlobalRegridding/test/test_conservative.jl` and `test_proj.jl`, 56
assertions net (4,657 → 4,713 passing, one broken unchanged):

- **prepared and unprepared are the same weights**, for a raster destination
  and a cell-space destination: equal sparse matrices, equal `colptr`, equal
  `rowvals`, equal `nzrange` per column, `===` on every stored value and
  denominator, and the same again against the generic coordinate-list route.
  An empty source keeps the degenerate contract of the index form.
- **one restricted tree per tile, one synthesis per prepared row**: a tile
  built against two source chunks builds one tree, and the wrapper space's
  synthesis counter advances by exactly the number of filled rows over both
  blocks together. Every row a block weighs is filled; the rows neither source
  chunk reaches are never synthesized. Preparing per block instead builds two
  trees, and its weights are identical bit for bit.
- **concurrent blocks share one cache**: four blocks of one tile spawned under
  `OUTER_PARALLEL => true` build one tree, fill each row exactly once between
  them, and each equals its serially built reference.
- **the budget rule**: what a cell is charged, per destination kind and per
  method; that a tile is sized by the executor's cost alone; that a tile at
  the size cap is refused prepared geometry and a smaller one is granted it;
  and that a cheap destination, a point method, too small a budget and an
  empty tile each answer with the index set.
- **a method that reads no prepared geometry still builds**: handed a
  prepared destination, a method on the generic coordinate-list route names
  the tile's cells through the space, fills no slot, and produces the
  index-form block.
- **nothing added to the clipping loop**: one operator call against a prepared
  destination allocates zero bytes, and no more than the same call against a
  memo.
- `test_lazy.jl`'s "budget tiles build one restricted tree each" is unchanged
  and now proves the executor prepares: without preparation each of the tiles'
  several blocks would build its own tree, and the count would exceed the tile
  count.
- `test_proj.jl` keeps the Proj obligation on the new object: a prepared
  destination's stored tree retains the safe chart across block tasks, and a
  serial build still derives a private prepared cursor from it.

The mutants these kill: a prepared build whose values, layout or denominators
differ from the unprepared one; a fill that predicts rows instead of reading
the candidate list, which the untouched-row count catches; a fill that runs
per block or per task, which the tree and synthesis counts catch; a cache
shared without a lock, which the concurrent test exercises; a prepared
destination handed to a method that cannot read one; a clipping loop
that synthesizes behind the prepared slot, which the allocation floor catches;
a budget rule replaced by a constant, and a tiling that quietly charges
polygons.

Dropped with the wrappers they tested: "one tile's cell geometry, synthesized
once" (`TileCells` forwarding and `CachedCellTree` fall-through), "the
destination tile cell cache is bounded" (`_TILE_CELL_CACHE_MAX`), and "one
tile's restricted tree, built once", whose claim the new counting test makes
more precisely.

## 6. Gate

One session, the two states interleaved rep by rep out of two checkouts, three
repetitions of each arm, Julia 1.12.6, 8 threads. Every arm waited for the
machine to be idle of other benchmark processes before it started, and an
earlier attempt at this table was discarded because another benchmark
overlapped it. Both states produce identical nonzero counts and weight sums on
both workloads. Seconds are the median of the three repetitions, peak RSS their
mean.

| workload / mode | seconds before | seconds after | ratio | RSS before | RSS after | ratio |
|---|---:|---:|---:|---:|---:|---:|
| block / direct | 4.881 | 2.625 | 0.538 | 866.0 MiB | 809.8 MiB | 0.935 |
| block / chunked | 2.527 | 2.548 | 1.008 | 844.6 MiB | 819.0 MiB | 0.970 |
| roundtrip / direct | 11.289 | 11.776 | 1.043 | 2,716.0 MiB | 2,794.2 MiB | 1.029 |
| roundtrip / chunked | 11.856 | 12.104 | 1.021 | 2,594.8 MiB | 2,583.8 MiB | 0.996 |

Neither benchmark warms its build, so each arm times one cold call: the build
plus what compiling it costs. The round-trip's two arms move by 2 % and 4 %
against within-state ranges that together span 11.2–12.4 s, and peak RSS by at
most 3 % against the 200–320 MiB spread this workload has always had. The block
workload's chunked arm is level, and its eager arm is 46 % faster: that route
used to compile its assembly for every tree type the destination space can
produce, and the barrier in §1 leaves it compiling the one it reaches.

Timing the same build three times in one process separates the build from
compiling it:

| state / mode | first call | second | third | first allocated | second allocated |
|---|---:|---:|---:|---:|---:|
| before / direct | 5.453 s | 0.123 s | 0.128 s | 1,424,871,584 | 255,315,040 |
| before / chunked | 3.240 s | 0.126 s | 0.126 s | 1,329,831,584 | 255,315,536 |
| after / direct | 3.449 s | 0.136 s | 0.139 s | 1,392,899,440 | 255,315,040 |
| after / chunked | 3.393 s | 0.128 s | 0.124 s | 1,451,921,888 | 255,315,136 |

The steady-state build is 0.123–0.139 s in all four states, allocating the same
bytes exactly on the eager route and within 400 bytes of 255 MB on the chunked
one: the weights cost what they cost before. First calls hold no surprise
either — the eager route drops 2.0 s and the chunked route is level, which is
what the benchmark arms report.

## 7. Suites

- `lib/GlobalRegridding/test/runtests.jl`: 4,713 passed, 0 failed, 0 errored,
  1 broken, zero deprecation warnings, against a baseline of 4,657 / 0 / 0 / 1
  measured in this session on the untouched checkout. The 56 new assertions are
  the six testsets in §5 (81 assertions) less the three testsets that tested
  the deleted wrappers (25).
- Bridge cross-check (`test/systems/crosssystem/`
  `regridding_conservation.jl`, `regrid.jl`, `regrid_acceptance.jl`):
  `XCHECK TOTAL passed=380 failed=0 errored=0 broken=12`, unchanged, zero
  deprecation warnings. This is where `test/systems/crosssystem/regrid.jl`
  runs; it is not a standalone entry point.
- Eager and chunked weights and values for `Conservative`, `NearestCell` and
  `BilinearPoint` over the lib fixtures: 54 serialized entries, all `==` and
  all `isequal` to the reference taken on the untouched checkout.
