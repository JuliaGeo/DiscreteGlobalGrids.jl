# Nearest regridding of Copernicus DEM tiles to IGeo7

> Superseded for the locator cost by `2026-08-25-partialgrid-cellat.md`.

- Date: 2026-08-25
- Benchmark: `benchmark/copdem_nearest.jl`
- Nothing in the library changed. This is what `NearestCell` costs on real
  Copernicus DEM data, and where that cost is.

`benchmark/point_tile_baseline.jl` priced the two nearest routes on a toy raster.
This note prices the same method on real 1x1-degree GeoTIFFs into IGeo7, where
source and destination are both `PartialGrid`s.

Headline: `cellat` on the source `PartialGrid` costs 6,367 ns and 347 bytes per
call where the same answer is available in 102 ns and no bytes, and about 0.7 s
whenever the source tiles do not form one contiguous run. §3–4 establish both,
§5 shows the effect on a regrid, §14 says what the second makes unreachable.

## 1. The workload

- Source: one `TiledDEM` over a `PartialGrid` of the listed tiles' pixels, built
  by `scripts/copdem_production.jl`, which the benchmark includes. Chunk `k` is
  tile `k`, decoded through `TileBuilder` into a `StripedLRUCache` of one slot
  per tile; `DGGSpace(grid; chunklevel = 0)`, so a chunk is one Copernicus tile.
- Destination: `query(sys, MultiOrderCoverage(box); level)` as a `CellVector` and
  a `PartialGrid`, tiled at `level - 7`, so a full tile is 7^7 = 823,543 cells —
  the chunk the production run writes. Whole ancestor subtrees would add cells no
  source can reach and dilute every per-cell number.
- `NearestCell()`, `Weighted(0.5)`, `lazy = true`, the default 2 GiB budget,
  `-t auto`.

## 2. The level the bare-system rule picks

`to = IGeo7System()` names no cells until a level is chosen; as a destination it
takes `levelfor(sys, src_space)`, the level whose median cell area is closest to
the source's.

| product | rows per degree | median pixel over the box | IGeo7 level | IGeo7 cell | cells globally |
|---|---|---|---|---|---|
| GLO-90 | 1200 | 76.88 m | **12** | 60.71 m | 138,412,872,012 |
| GLO-30 | 3600 | 25.63 m | **13** | 22.94 m | 968,890,104,072 |

- `scripts/copdem_production.jl` fixes level 12 for GLO-90, which the rule picks
  anyway, so GLO-90 has one arm. GLO-30 differs and both are measured: level 13
  is the rule's, level 12 the production setting applied to 30 m data.
- A bare system as destination names the whole globe at that level — 1.4e11 cells
  at level 12, 9.7e11 at level 13. No arm sweeps that; each names its box's cells
  and §11 extrapolates.

## 3. Where the time goes: `cellat` on the source

One destination tile of 823,543 cells, GLO-90 N46 E010 into IGeo7 level 12, each
phase in its own compiled function so a top-level loop over globals cannot price
boxing instead of work. One shot per phase after a warm-up on a 1,024-cell slice.

| phase | s | ns per destination cell | bytes |
|---|---|---|---|
| `CellVector` decode alone | 0.050 | 60.3 | 0 |
| `sites[i]` — decode + IGeo7 centroid | 0.106 | 128.3 | 0 |
| `weightsat!` — site read + `cellat` | 5.249 | 6,373.5 | 285,488,416 |
| `cellat` alone, sites precomputed | 5.011 | **6,085.2** | 285,488,416 |
| `chunkat` over the located indices | 0.001 | 1.1 | 0 |
| **whole `tileweights`** | **5.180** | **6,290.1** | 401,288,656 |
| residual (`tileweights` less the three above) | 0.062 | 75.6 | |
| cold read of the same tile | 5.453 | 6,621.0 | 511,855,280 |
| warm read of the same tile | 0.014 | 16.8 | 40,408,096 |

- The phases sum: 128.3 + 6,085.2 + 1.1 + 75.6 = 6,290.2 ns, and `tileweights` is
  95.0 % of the cold read. The remaining 331 ns per cell is the GeoTIFF decode,
  the source load and the weight application together.
- A warm read is 16.8 ns per cell, so applying weights is 0.25 % of a cold read:
  every question about this workload is a question about `cellat`.
- Allocation follows: 622 bytes per destination cell over the cold read, 487 in
  `tileweights`, 347 in `cellat` alone. A nearest stencil is one entry of one;
  none of those bytes is the stencil.
- 46.55 s wall against 52.55 s user, a ratio of 1.13, and `tileweights` is a
  serial loop, so per-core and aggregate throughput are the same number here.

## 4. The same answer costs 102 ns

`CopernicusDEMSystem` defines a closed-form `cellat(::LevelGrid, p)` — a
latitude-band lookup and two floors — and `DGGSpace`'s `cellat` reduces whatever
the grid answers through `localindex` anyway. Composing the two directly, over
100,000 sample sites of the same tile:

| route | ns per call | bytes per call |
|---|---|---|
| `cellat(levelgrid(sys, 1), p)` then `localindex(grid, c)` | **102.2** | **0** |
| `cellat(space, p)` on the `DGGSpace` over that grid, as shipped | **6,366.6** | **346.7** |
| `treeify(grid)` on its own | 47.1 | 352 |

- **62x for the same answer.** `treeify` is not the cost: memoized, 47 ns.
- The cost is what `cellat(grid::AbstractGrid, p)` (`src/fallbacks/locate.jl`)
  does with the tree — an `STI.query` over the tile's 1,440,000 pixel cells, a
  `sort!` of survivors into a fresh vector, then an exact `point_in_cell` against
  a `cell_boundary` polygon built per candidate.
- That is the generic locator for a grid with no closed-form method, and
  `PartialGrid` defines none, so a subset over a system that has one still takes
  it.

## 5. `NearestCell` is slower than `Conservative` on this data

Same destination tile against the same source chunk, one session, 20,000
destination cells each:

| method, `buildweights!` on one (tile, chunk) pair | ns per destination cell | bytes |
|---|---|---|
| `Conservative()` | 3,186.5 | 84,690,640 |
| `NearestCell()` | 6,758.9 | 8,307,376 |

The lookup is **2.1x the clip**. `Conservative` never calls `cellat`: it descends
the source tree once per destination cell, where the nearest locator pays a tree
query, a sort and a polygon test for a point it could compute arithmetically.

## 6. The cliff: a source that is not one contiguous run

`treeify(::PartialGrid{<:CopernicusDEMSystem})` gives a `BlockCursor` only when
the cells are one lattice rectangle held as one contiguous id run, and a
`HierarchicalGridCursor` otherwise. Tiles in one latitude row are such a run;
tiles spanning two rows are not, and `cellat` inherits the split. 200 calls per
fast shape, 3 per slow one:

| source tiles | chunks | tree `treeify` gives | ms per `cellat` |
|---|---|---|---|
| N46 E010 | 1 | `MemoBlockCursor` | 0.007 |
| N46 E010–E011 | 2 | `MemoBlockCursor` | 0.008 |
| N46 E010–E012 | 3 | `MemoBlockCursor` | 0.007 |
| N46 E010–E013 | 4 | `MemoBlockCursor` | 0.007 |
| N46–N47 E010 | 2 | `HierarchicalGridCursor` | **662.3** |
| N46–N47 E010–E011 | 4 | `HierarchicalGridCursor` | **731.7** |

- Four orders of magnitude, decided by whether the tile set crosses a latitude
  row. A 996-cell destination against the 2x2 block did not finish one read in
  30 seconds of profiling, 95.6 % of package samples under
  `tileweights` → `weightsat!`.
- Not a corner case: every production destination is a chunk of a global run, and
  a 55 km IGeo7 level-5 cell straddles a latitude row wherever it crosses a whole
  degree. A 2x2-degree box is the smallest interesting one and is already on the
  wrong side.

## 7. Throughput, GLO-90 and GLO-30

Only east-west strips within one latitude row are reachable (§6), so boxes grow
along longitude from latitude row 46. Each arm is a fresh process; `cold sweep`
reads every destination tile once on a plan never read, `warm sweep` repeats it
against the plan that sweep left behind.

### GLO-90 into IGeo7 level 12

| source tiles | destination cells | destination tiles | candidates per tile (med/max) | destination arm | cold sweep | cold cells/s | warm sweep | warm cells/s | peak RSS |
|---|---|---|---|---|---|---|---|---|---|
| 1 (E010) | 2,313,802 | 8 | 1.0 / 1 | 1.163 s | 15.809 s | 146,364 | 0.039 s | 59,300,001 | 1.27 GiB |
| 2 (E010-E011) | 4,624,703 | 13 | 2.0 / 2 | 1.543 s | 31.216 s | 148,153 | 0.048 s | 96,922,022 | 1.48 GiB |
| 3 (E010-E012) | 6,935,841 | 18 | 2.0 / 3 | 1.717 s | 48.295 s | 143,613 | 0.079 s | 87,629,025 | 1.65 GiB |
| 4 (E010-E013) | 9,246,816 | 21 | 2.0 / 3 | 2.070 s | 64.063 s | 144,338 | **63.974 s** | **144,540** | 1.67 GiB |

| source tiles | source chunk reads | pixels copied | GeoTIFF decodes | resident weight bytes | budget | destination tiles cached |
|---|---|---|---|---|---|---|
| 1 | 8 | 11,520,000 | 1 | 147,619,072 | 536,870,912 | 8 of 8 |
| 2 | 16 | 23,040,000 | 2 | 309,841,656 | 536,870,912 | 13 of 13 |
| 3 | 25 | 36,000,000 | 3 | 486,996,576 | 536,870,912 | 18 of 18 |
| 4 | 33 | 47,520,000 | 4 | 525,094,312 | 536,870,912 | **16 of 21** |

`Conservative()` over the identical one-tile destination: 6.284 s and 368,199
cells/s against nearest's 15.809 s and 146,364 — **2.52x faster end to end**, the
whole-pipeline form of §5.

### GLO-30 into IGeo7

| destination | source tiles | destination cells | destination tiles | candidates (med/max) | destination arm | cold sweep | cold cells/s | warm sweep | warm cells/s | peak RSS |
|---|---|---|---|---|---|---|---|---|---|---|
| level 13 (the rule's) | 1 | 16,181,892 | 33 | 1.0 / 1 | 5.678 s | 119.084 s | 135,886 | 118.608 s | 136,432 | 2.05 GiB |
| level 13 | 2 | 32,351,955 | 59 | 1.0 / 2 | 7.058 s | 251.166 s | 128,807 | 246.933 s | 131,015 | 2.15 GiB |
| level 13 | 3 | 48,523,116 | 81 | 2.0 / 3 | 7.652 s | 371.226 s | 130,710 | 369.284 s | 131,398 | 2.37 GiB |
| level 12 (production's) | 1 | 2,313,802 | 8 | 1.0 / 1 | 1.147 s | 18.218 s | 127,009 | 17.876 s | 129,438 | 1.97 GiB |

- The level-12 arm is every GLO-90 arm's destination fed by a nine times finer
  source, at 127,009 cells/s against 146,364: nine times finer costs 13 %, so
  `cellat`'s cost is the tree query and the polygon test, not the pixels.
- **Every GLO-30 arm sweeps warm at cold speed.** A weight block's reference
  vector is sized by the source chunk it names, so the level-12 destination holds
  451,737,104 B against GLO-90's 147,619,072 for identical cells — three times,
  for a 12,960,000-pixel tile against a 1,440,000-pixel one. Only 4 of 8 tiles
  fit the 512 MiB budget, and an in-order sweep against an under-sized LRU reuses
  nothing: the pair-route arm places 16,181,892 cells cold and as many warm.
- Same boundary in GLO-90's four-tile warm sweep: 21 tiles needing 525,094,312 B
  against a 536,870,912 B budget keeps 16, and warm costs what cold did (63.974 s
  against 64.063 s), where three tiles fit and sweep warm in 0.079 s, eight
  hundred times faster. The cliff is at the budget, not the box.
- Source reads amplify with the destination tiling, not the box: the level-13
  three-tile arm read its 3 chunks 98 times and copied 1,270,080,000 pixels,
  5.1 GB of memcpy for 128 MB of GeoTIFF decoded three times.
- Throughput is flat across box size — 146,364, 148,153, 143,613, 144,338 cells/s
  over four-fold growth on GLO-90, 135,886, 128,807, 130,710 over three-fold on
  GLO-30 — so the cost is per destination cell. At 145,000 cells/s a cell costs
  6.9 us, of which §3 accounts 6.6 us in `cellat`.
- Destination-space construction grows sublinearly, the coverage set compacting
  whole subtrees while only the boundary refines: 1.163, 1.543, 1.717, 2.070 s
  for 2.3 to 9.2 million cells at level 12, allocating 760 MB to 1.70 GB, and
  5.678 to 7.652 s for 16.2 to 48.5 million at level 13. No regrid pays it, but a
  caller does.
- `plan_regrid` never exceeded 8.5 MB or a measurable fraction of a second.

## 8. The two routes on real data

`PairNearest` is `NearestCell`'s own `buildweights!` with no `sampler`, so a plan
around it builds one `(destination tile, source chunk)` pair at a time. Nothing
but the routing differs, and each pair of arms ran on the same data.

| product | source tiles | candidates (med/max) | pair locations | tile locations | location ratio | pair cold sweep | tile cold sweep | time ratio |
|---|---|---|---|---|---|---|---|---|
| GLO-90 | 1 | 1.0 / 1 | 2,313,802 | 2,313,802 | 1.00x | 15.111 s | 15.809 s | 0.96x |
| GLO-90 | 2 | 2.0 / 2 | 8,561,833 | 4,624,703 | **1.85x** | 56.876 s | 31.216 s | **1.82x** |
| GLO-30 | 1 | 1.0 / 1 | 16,181,892 | 16,181,892 | 1.00x | 117.745 s | 119.084 s | 0.99x |
| GLO-30 | 2 | 1.0 / 2 | 50,372,164 | 32,351,955 | **1.56x** | 377.992 s | 251.166 s | **1.51x** |

- At one source tile the routes do identical work and the pair route is 1–4 %
  ahead: the tile route's `chunkat` per entry and its manifest, nothing else.
- At two source tiles the pair route performs 1.85 and 1.56 locations per
  destination cell against the tile route's one, and costs 1.82x and 1.51x the
  time. **Location ratio and time ratio agree to within 3 %.** The toy raster said
  4.00x locations, 3.13x time; here, where location is 92 % of the read, the two
  coincide instead of diverging.
- The multiplicity is small because a destination tile is smaller than a source
  tile — a 55 km IGeo7 level-5 cell against a 111 km by 77 km tile — so most
  destination tiles sit inside one source tile. A coarser destination tiling, or
  a finer-chunked source, would raise it.

## 9. The profile

Two `Profile` listings over the same four destination tiles (526,893 cells) of
the two-tile GLO-30 strip, on a plan of its own after every source tile is
decoded and resident, so neither is a profile of GDAL. Only samples whose
backtrace reaches this package are counted; under `-t auto` the idle threads park
in `__psynch_cvwait` and are 94 % of all samples. The east-west strip stands in
for the 2x2 box, which §6 makes unprofilable.

**Cold read, weights built here** (1,273 samples reaching the package, of 22,032):

| inclusive | frame |
|---|---|
| 98.90% | `readblock!` at `lazy.jl:422` |
| 97.96% | `_readdestination!` → `tilefor` → `gettile!` |
| 94.27% | `tileweights` at `plans.jl:1000` |
| 92.54% | `weightsat!` at `interpolation.jl:76` |
| 89.95% | `cellat` at `regridding.jl:158` (the `DGGSpace` method) |
| **84.68%** | **`cellat` at `locate.jl:14` (the generic fallback)** |
| 26.16% | `node_extent` at `cursor.jl:247` |
| 22.94% | `child_indices_extents` at `cursor_memo.jl:101` |
| 21.45% | `_box_cap` at `cursor.jl:200` |
| 20.42% | `_leafcell` at `cursor.jl:289` |
| 6.91% | `cap_contains` at `caps.jl:25` |

- Top self time: `sind` 12.25 %, `cosd` 8.32 % across two entries, `_taskmemo`
  5.89 % across three, `node_extent` 3.46 %, `getchild` 1.65 % — trigonometry and
  node-extent bookkeeping, with no frame belonging to a nearest stencil.
- The chain reads plainly: 92.5 % of the cold read is `weightsat!`, 90 % the
  `cellat` it calls, and 84.7 % inside the generic locator alone.
- **Warm read, weights cached** (11 samples of 208): `_applygroup!`/`anyinvalid`
  54.55 %, `_readsource!` and the `TiledDEM` `readblock!` 27.27 %, `applyblock!`
  9.09 %. Too cheap to profile at this delay — the same statement as §3's 16.8 ns.

## 10. Correctness on real data

`check=1000` verifies random destination cells against the source read directly:
the sample site is the cell's centroid, `cellat` on the source space names the
pixel, and the regridded value must equal it bit for bit (`isequal`, so a NaN
pixel must give a NaN cell); where `cellat` finds nothing the cell must be
missing. A nearest stencil is one entry of exactly 1.0 and `Weighted` divides by
exactly one, so bit identity is the right assertion.

| arm | cells checked | mapped | outside coverage | wrong |
|---|---|---|---|---|
| GLO-90, 1 tile | 1,000 | 996 | 4 | **0** |
| GLO-90, 2 tiles | 1,000 | 997 | 3 | **0** |
| GLO-90, 3 tiles | 1,000 | 998 | 2 | **0** |
| GLO-90, 4 tiles | 1,000 | 996 | 4 | **0** |
| GLO-90, 1 tile, pair route | 1,000 | 996 | 4 | **0** |
| GLO-90, 2 tiles, pair route | 1,000 | 997 | 3 | **0** |
| GLO-30, level 13, 1 tile | 1,000 | 999 | 1 | **0** |
| GLO-30, level 13, 2 tiles | 1,000 | 997 | 3 | **0** |
| GLO-30, level 13, 3 tiles | 1,000 | 997 | 3 | **0** |
| GLO-30, level 13, 1 tile, pair route | 1,000 | 999 | 1 | **0** |
| GLO-30, level 13, 2 tiles, pair route | 1,000 | 997 | 3 | **0** |
| GLO-30, level 12 | 1,000 | 996 | 4 | **0** |

Twelve thousand cells, none wrong. The cells outside coverage are the box's own
boundary — the coverage query keeps every IGeo7 cell the box touches, so a cell
whose centroid falls just outside a tile's 1-degree extent has no source, and
every one came back missing as `Weighted(0.5)` requires. Finite fractions run
99.7 % to 99.9 %, elevations 22.5 m to 3,874.1 m, which is the Alps.

## 11. What a full globe would take

The production GLO-90 destination is 66,175 IGeo7 level-5 chunks of 823,543
cells, 5.45e10 cells; at level 13 the same land is seven times that, 3.81e11.
Extrapolation assumes the rate stays flat, which §7 establishes over four-fold
growth in box size and no further, and no scheduling loss.

| run | destination cells | at the measured rate, one core | at 40 workers, perfect scaling |
|---|---|---|---|
| GLO-90 to level 12 | 5.45e10 | 145,000 cells/s -> **104 hours** | 2.6 hours |
| GLO-30 to level 13 | 3.81e11 | 130,000 cells/s -> **815 hours** | 20 hours |

- **Neither is reachable today**: both rates belong to a source that is one
  contiguous run, and a 55 km level-5 cell crosses a whole degree of latitude
  nearly everywhere, so nearly every production chunk takes the
  `HierarchicalGridCursor` path at 0.7 s per destination cell — 3.8e10 seconds of
  core time for GLO-90 alone, about 1,200 years.
- With `cellat` composed from §4's closed form a destination cell would cost
  128 ns for its site, 102 ns for the lookup, 1 ns for `chunkat` and 76 ns for
  row and accumulator work — 307 ns of weights — plus 331 ns of source read and
  application: about 640 ns, or 1.57 million cells/s, roughly 9.6 hours for
  GLO-90 and 67 for GLO-30 on one core, cliff included, the closed form
  consulting no tree. Arithmetic on measured parts, not a measurement.
- Data is the other wall: 26,475 tiles of about 5 MB is roughly 132 GB for
  GLO-90, the same tiles at about 43 MB roughly 1.1 TB for GLO-30. Neither fits
  the 23 GiB free here, which is why every arm is a box.

## 12. Data footprint

Downloaded into the directory `data=` named, outside the repository, and safe to
delete: `copdem-data/`, 280 MB total, 288,835,243 bytes of GeoTIFF.

- `tiles-glo90/` — 77 MB, 15 tiles (N46–N49 x E010–E013) plus one symlink to a
  tile already in `$RASTERDATASOURCES_PATH`.
- `tiles-glo30/` — 203 MB, 5 tiles (E012 for N46–N47, and N48 x E010–E012) plus
  four symlinks to tiles already there.
- The symlinked tiles cost nothing and are untouched; nothing was written into
  `$RASTERDATASOURCES_PATH`. The 6 GiB download cap was never approached: time,
  not disk, bounded the arms, and §6 bounded them harder than either.
- `~/Downloads/Copernicus_DSM_10_N46_00_E010_00_DEM.tif` was read but not used:
  3601x3601, the ESA distribution shape, where the AWS COG is 3600x3600 and is
  what `CopernicusDEMSystem(30)` and `readtile` require. The COG was fetched.

## 13. Metadata

- Julia 1.12.6, macOS 26.5, M-series, `-t auto` = 8 threads and 8 GC threads of
  12 CPUs, `pmset -g` reports `powermode 2`, 24 GiB of RAM.
- Every timed arm ran with the machine otherwise idle: the driver waits until no
  other Julia process matching a benchmark or scratch path is running, and each
  arm's log records the count before and after it.
- One shot per arm after a compilation warm-up. Ratios inside one process are the
  portable part; absolutes belong to this machine state and power mode.
- Source reads come from a `DiskArrays` wrapper recording every `readblock!`
  range against its tile; decodes and decoded pixels from `TileBuilder`'s
  atomics; residency from `residency(::LazyRegridArray)`; resident weight bytes
  from the plan's `PerChunk` accounting with the cached tiles' manifest lengths;
  peak RSS from `Sys.maxrss()`, warm-up included.
- Profile listings are self time over samples whose backtrace reaches this
  package; both totals are printed so the discarded fraction is visible.

## 14. What this does not show

- **No box that crosses a latitude row** (§6). The 2x2-degree box the plan called
  for, and every larger block, is unreachable: at 0.7 s per destination cell a
  9.3-million-cell destination is about 75 days.
- **No full-globe or production-shaped run.** A production destination chunk is a
  level-5 subtree, not a box, and its source is a multi-row tile set for most
  chunks.
- **One machine state, one process per arm.** Every arm is a one-shot measurement
  of a phase running for seconds; there is no min-of-n. Compilation is removed
  instead, by one complete miniature regrid before any arm runs.
- **Nothing about accuracy beyond §10's bit-identity check.** No arm compares
  against an independent implementation.
- **Nothing about the destination being IGeo7 specifically.** The destination side
  costs 128 ns per cell (§3) and does not dominate; another system would move
  that number and leave the conclusion.
