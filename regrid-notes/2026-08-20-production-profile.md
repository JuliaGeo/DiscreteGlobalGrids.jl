# Production profile: where the CopDEM→IGeo7 L12 run's time actually goes (2026-08-20)

Companion data: `2026-08-20-production-profile.ndjson` (73 measurements).
Methodology: `regrid-notes/profile-performance.md` / the `profile-performance`
skill — warm up, one-shot `@timed`/`withcpu` baselines, `Profile.@profile` +
`FlameGraphs` traversal and aggregation, `Profile.Allocs`. No GUI, no
ProfileView.

The live run (PID 2915199) was **not touched**. Every experiment ran in a
separate process at `-t 12` or below against its own scratch store
(`/home/asinghvi17/geo/dggstores/profile-scratch/prof.zarr`), on the same
worktree code (`/home/asinghvi17/geo/DGG-subzone-store`, branch
`claude/copdem-production`, git clean), the same `bench` environment, CR
`6a4b997` / GO `2825c47`, Julia 1.12.6. Load was checked before each run and
stayed between 47 and 68.

---

## 1. Headline

Three claims, each measured:

1. **The 63-thread process was never going to hold more than ~22 cores.** The
   parallel shape gives each of the 21 workers a *hard* ceiling of ~1.3 cores,
   because the only parallelism inside a work unit is the lazy plan's
   source-chunk wave, and the wave's blocks are 70–100 % one block. Measured
   ideal wave speedup at mid latitude: **1.05–1.44**. Measured actual:
   **1.23–1.27 cores** for one worker on 12 threads. The run holds 22.2 cores
   with 21 workers — i.e. it is running at **~95 % of its own design ceiling**.
2. **There is no per-core regression.** Production converts 823 543 cells in
   ~28.0 core-seconds; a controlled single-threaded run of the identical work
   unit on this box today takes 28.3 core-seconds. The "25.5 k cells/s per busy
   core" in §15.9.1 is aggregate-cells/s ÷ instantaneous-process-cores, which
   charges the 8 GC threads and the write/log/heartbeat path to the compute; the
   per-worker figure from the run's own done log is **29 688 cells per
   worker-second**, and per core-second **28.0 k**.
3. **Two-thirds of every CPU second is spent *finding* candidate cell pairs, and
   half of that is recomputing spherical caps for source-tree nodes that never
   change.** `node_extent` on the CopernicusDEM `BlockCursor` is **47.7 %** of
   all CPU; the actual intersection-area arithmetic is only **24.7 %**. The
   destination tree got a cap cache in PR #44; the *source* tree did not.

So: at a 24-core budget the parallel shape is a dead end (≤ 8 % left), and the
whole remaining win is in the per-cell cost of the hot path.

---

## 2. The live run as measured

`ps`, `/proc/2915199/task/*/stat` deltas, and the run's own heartbeats:

| | 19:15Z | 20:12Z |
|---|---:|---:|
| columns done | ~1 930 | 4 420 / 66 178 |
| aggregate cells/s | 585 000 | 636 379 |
| process cores (10 s Σ utime+stime delta) | 22.30 | 22.24 |
| threads above 0.5 core | 23 | 23 |
| threads above 0.9 core | **0** | — |
| busiest single thread | 0.79 core | — |
| OS threads in process | 136 | — |
| RSS | 29.3 GiB | 35.3 GiB |
| 1-min load (box, shared) | 65.4 | 47.7 |
| ETA | 25.5 h | 22.2 h |

Cells per core-second: 636 379 / 22.24 = **28 614**.

Done-log distribution over the first 1 930 columns: mean 27.74 s, median 27.73,
p10 17.5, p90 38.6, max 74.2. Rate is flat in the NaN fraction (30–31 k cells/s
for every band from 0 % to 80 % NaN) and jumps only for the **fully empty**
columns (100 % NaN → 68 214 cells/s, 311 of 1 930). A 100 %-NaN column still
costs ~12 s, but it is 2.2× cheaper than a fed one, not free — §15.9.1's
candidate 2 is real but small.

The box is 2× AMD EPYC 9354, 64 physical cores, **no SMT**, 2 NUMA nodes, shared
with other tenants at load 47–68 throughout.

---

## 3. Q1 — Where the worker threads idle

### 3.1 The per-column time budget: everything except weight building is free

`timed_column` phase timing, 3 mid-latitude columns, `-t 1` (`exp:"phase"`,
tag `serial_t1`):

| phase | seconds for 3 columns | share |
|---|---:|---:|
| `subtree` + `DGGSpace(...; chunklevel=5)` (destination space) | 0.000047 | 0.00005 % |
| `GR.regrid(...; lazy=true)` (plan, tiling, caps, source tree) | 0.0071 | 0.008 % |
| `collect` (the whole read: covering, weights, source, assembly) | **87.98** | **99.98 %** |
| `dggwrite!` (Zarr chunk, 823 543 f32, Blosc) | 0.011 | 0.013 % |

§15.5's rooted-`PartialGrid` fix worked completely: the destination space costs
**4 µs**, not the 0.30 s la-choice measured. Plan build is 2 ms. The write is
4 ms. **There is no serial per-column section worth attacking outside the
regrid itself.**

Inside `collect`, the flame graph's topmost-inclusive aggregation
(`exp:"profile"`, `serial_t1`, 17 988 samples at 5 ms):

| frame | inclusive |
|---|---:|
| `_readdestination!` | 99.7 % |
| `blockfor` (weight-block construction) | **99.6 %** |
| `_applygroup!` (sparse apply) | 0.06 % |
| `_sourcefor!` → `readblock!` → `tilevalues!` → `synthetic_tile` | 0.02 % |
| `_writechunk!` (missing policy, renormalise, fill out) | 0.01 % |
| `_connectedsource!` / `connectedchunks!` | 0.00 % |
| `TileCells` / cap-tree | 0.00 % |

Synthetic tile decode, mask sampling and TiledDEM churn — the things §15.9.1
worried about — are **0.02 % of the run**. The striped cache is doing its job;
this is not an I/O or a source-materialisation problem.

### 3.2 The wave is bounded, and worse, it is imbalanced

Destination tiling is **one tile per column**: `_defaulttilesizes` gives
`min(budget/320, 823 543) = 823 543`, so `runs = [1:823543]`, one
`TileCells`, one `_connectedsource!`, one wave loop, one `_writechunk!`.
The wave is therefore the *only* parallelism in a work unit, and
`w = _wavesize = min(nthreads, n_srcchunks, weightbudget ÷ floorbytes)`.

**Source-chunk fan-out** over 3 000 randomly sampled production columns
(`exp:"fanout_dist"`, plan build is 2 ms so this is cheap and exact):

| n source chunks | 1 | 2 | 3 | 4 | 5 | 6 | 7–10 | 11–20 | >20 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| columns | 138 | 182 | 389 | 1277 | 361 | 214 | 291 | 90 | 58 |

mean 5.03, median 4, 66 % of columns at ≤ 4. By band: −90..−60 mean 12.7
(max 224), 60..90 mean 7.1, everything between 3.4–4.3. **Harmonic-mean
speedup cap if `blockfor` were perfectly parallel: 3.61×.**

But the blocks in a wave are not equal. Serial per-block timing
(`exp:"blocktimes"`, `-t 1`), seconds and nonzeros per source chunk:

| column | lat | n src | per-block seconds | nnz | max block | **ideal wave speedup** |
|---|---:|---:|---|---|---:|---:|
| 728 | 47.3 | 5 | 0.00 **20.17** 7.19 0.01 0.01 | 0, 2.16 M, 0.76 M, 0, 0 | 74 % | **1.36** |
| 729 | 47.3 | 4 | 0.00 0.00 **25.89** 1.41 | 0, 0, 2.81 M, 0.13 M | 95 % | **1.05** |
| 730 | 47.3 | 5 | 2.80 0.01 0.01 **23.38** 0.82 | 0.31 M, 0, 0, 2.51 M, 0.09 M | 87 % | **1.16** |
| 98241 | −33.8 | 1 | **26.36** | 2.55 M | 100 % | **1.00** |
| 34818 | 9.6 | 4 | 0.16 0.28 **18.57** 7.67 | 0, 0, 1.75 M, 0.70 M | 70 % | **1.44** |
| 34819 | 9.6 | 4 | 0.00 1.14 0.00 **25.24** | 0, 0.11 M, 0, 2.35 M | 96 % | **1.05** |
| 115426 | −76.7 | 10 | 0.01×5, 8.82, **13.54**, 5.43, 0, 0 | … 3.6 K, 0.95 M, 1.41 M, 0.60 M … | 49 % | **2.06** |
| 115427 | −76.7 | 8 | 0.00 0.00 1.79 12.20 **12.35** 1.62 0 0 | … | 44 % | **2.26** |

Two things at once:

* **2 to 6 of the "connected" source chunks produce zero pairs** and cost
  ~0.01 s. The cap/tree connectivity query over-reports; `_wavesize` counts
  those chunks as parallel work that does not exist.
* Of the chunks that *do* produce pairs, **one carries 44–100 %** of the block
  time, because block cost tracks nnz and a level-5 column's overlap with the
  1° tile lattice is lopsided.

Mean ideal wave speedup over these twelve columns: **1.30**.

### 3.3 The accounting closes

Directly measured, one worker (W = 1) on a 12-thread process with
`GR.OUTER_PARALLEL => true` exactly as the worker body sets it
(`exp:"phase"`, tag `W1_t12`):

| batch | wall, 3 cols | **cores held** | cells/s | cells/s/core |
|---|---:|---:|---:|---:|
| midlat_N (nsrc 4–5) | 74.5 s | **1.27** | 33 141 | 26 144 |
| midlat_S (nsrc 1) | 39.6 s | **1.17** | 62 439 | 53 593 |
| polar_S (nsrc 3–10) | 38.5 s | **1.96** | 64 136 | 32 741 |

against 88.0 s / 0.99 cores for the same mid-latitude batch at `-t 1`: the
12-thread wave buys **1.18× wall for 1.27 cores**. In the flame graph of that
run, `poptask` (idle worker threads) is **90.2 % of samples** — the other 11
threads are asleep, not contending.

W-sweep through the real `runcolumns` loop, 33 contiguous mid-latitude columns,
batch = 1, `-t 12` (`exp:"wsweep"`):

| W | wall | cores | cells/s | **cells/s/core** | cores/worker |
|---:|---:|---:|---:|---:|---:|
| 4 | 198.8 s | 4.92 | 136 733 | 27 793 | 1.23 |
| 8 | 122.2 s | 7.88 | 222 357 | 28 233 | 0.98 |
| 11 | 109.7 s | 8.85 | 247 793 | 27 997 | 0.80 |

Per-core-second throughput is **flat at 27.8–28.2 k** across W: no
memory-bandwidth or lock cliff, and the tile cache, `donelock` and `p.lock` are
not contended. `cores/worker` falls only because W approaches `nthreads` and the
wave has nowhere to spawn — at W = 21 on 63 threads there is plenty of room,
which is why the live run gets 1.06.

**63 threads → 22.2 busy, explained:**
21 workers × (1 base thread + ~0.06–0.3 of a wave that is one dominant block)
≈ 22–25 cores. The remaining 40 threads have nothing to be given.

### 3.4 The one lever that *does* scale, and why the code cannot reach it

CR's own per-weight-build threading (`_intersectionareas` → `_innerthreaded`,
plus `_leafcaps` in `src/cap_cached_tree.jl`) is what `OUTER_PARALLEL` turns
off. Measured on the identical work, all source chunks of a column built
through `buildblock`, `-t 12` (`exp:"nowave"`):

| column | lat | nnz | inner OFF: wall / cores | inner ON: wall / cores | speedup |
|---|---:|---:|---|---|---:|
| 728 | 47.3 | 2.92 M | 28.88 s / 1.24 | **3.32 s / 9.43** | **8.7×** |
| 729 | 47.3 | 2.93 M | 28.78 s / 1.21 | 3.74 s / 8.35 | 7.7× |
| 34818 | 9.6 | 2.46 M | 27.60 s / 1.18 | 3.24 s / 9.48 | 8.5× |
| 34819 | 9.6 | 2.46 M | 27.81 s / 1.19 | 3.22 s / 9.18 | 8.6× |
| 115424 | −76.7 | 0.85 M | 9.61 s / 1.14 | 1.37 s / 8.33 | 7.0× |
| 115426 | −76.7 | 2.95 M | 28.82 s / 1.23 | 3.53 s / 9.06 | 8.2× |
| 98241 | −33.8 | 2.55 M | 27.24 s / 1.23 | 3.03 s / 9.44 | 9.0× |
| 98242 | −33.8 | 0.28 M | 4.12 s / 1.32 | 0.43 s / 9.55 | 9.6× |

**8.7× on 12 threads, uniform across latitude and across fan-out, ~92 %
efficiency per held core** (823 543 / 3.32 / 9.43 = 26.3 k cells/s/core vs 28.3 k
serial). This is the parallelism that works, and it is unreachable from
configuration: `_fillwave!` wraps every spawned block build in
`@with OUTER_PARALLEL => true` (lazy.jl:477), so **whenever a column has more
than one source chunk the wave itself suppresses inner threading**, regardless
of what the caller sets. Removing `GR.OUTER_PARALLEL => true` from the worker
body changes nothing for 95 % of columns. Verified: `W1_t12` with `OUTER=0`
reproduces `OUTER=1` to within noise (1.257 vs 1.268 cores).

---

## 4. Q2 — The hot path per core

Flame graph of the block builds for two mid-latitude columns, `-t 1`, 2 ms
sampling, 29 628 samples (`exp:"hotpath"`). Topmost-inclusive by function:

| frame | midlat_N | polar_S |
|---|---:|---:|
| `intersection_areas` (CR) | **93.0 %** | 90.0 % |
| ├ `get_all_candidate_pairs` / `dual_depth_first_search` | **67.8 %** | 64.9 % |
| │  ├ `node_extent` (source `BlockCursor`) | **47.7 %** | 46.0 % |
| │  │   └ `_box_cap` | 41.9 % | 40.0 % |
| │  │   └ `UnitSphereFromGeographic` / `sincosd` | 19.5 / 18.9 % | 19.5 / 18.9 % |
| │  ├ `child_indices_extents` (leaf caps) | 21.3 % | 20.2 % |
| │  └ `spherical_distance` (cap overlap test) | 24.7 % | 23.4 % |
| └ `assemble_sparse_matrix_coo` / `_assemble_chunk` | **24.7 %** | 24.0 % |
|    └ `BlockAreaOperator` / `IntersectionAreaOperator` | 24.6 / 17.2 % | 23.9 / 16.3 % |
| `TileCells` restricted-tree build (once per column) | ~5 % (1.5 s of 29) | ~5 % |

Top self-cost frames are all scalar float work inside those two phases —
`*` float.jl:497 16.5–17.1 %, `+` 8.5–9.1 %, `GenericMemory` boot.jl:588
**7.0–7.5 %** (allocation), `<` 6.7–7.1 %, `-` 6.0–6.6 %, `muladd` 5.1–5.6 %,
then `node_extent` cursor.jl:246/247 at 1.0–1.2 % self each,
`memmove` 1.0–1.2 %. **Runtime dispatch is 0.006 % — zero.** GC is 4.5–4.8 %.

### What this says

* **The pair search costs 2.7× the useful work.** Discovering which
  (destination cell, source pixel) pairs intersect is 67.8 % of CPU; computing
  their intersection areas is 24.7 %.
* **The source tree recomputes its geometry on every visit.**
  `STI.node_extent(::BlockCursor)` (`src/systems/CopernicusDEM/cursor.jl:245–248`)
  calls `_node_box` then `_box_cap`, and `_box_cap` runs `sincosd` /
  `UnitSphereFromGeographic` / `spherical_distance` per node. It is invoked for
  every interior node the dual DFS touches, on every block build, for every
  column, in every worker — for a tile whose tree is **identical every time**.
  Interior-node extents alone are 47.7 − 21.3 ≈ **26 % of all CPU**; the leaf
  caps built by `child_indices_extents` are the other ~21 %.
  The destination (IGeo7) side is already cached — `CapCachedTree`, PR #44. The
  source side is the half that was never done.
* **`child_indices_extents` allocates a fresh `Vector{Tuple{Int,Cap}}` per
  leaf.** Allocation profile of one column (`exp:"alloc_profile"`,
  `sample_rate=0.01`): **3.22 GiB allocated per column, 4 198 B per destination
  cell, GC 4.5 %**, of which
  **71.6 % is `Memory{Tuple{Int64, SphericalCap{Float64}}}`** — that vector,
  and nothing else. `Memory{UnitSphericalPoint}` (the cell polygons handed to
  the clipper) is 11.5 %, `Memory{SphericalCap}` 5.3 %.

### Against the 2026-08-19 post-stack plateau structure

That campaign's plateaus were *cap-tree build + assembly* and a *serial tail*.
Both are gone from this workload:

* the destination cap tree is now built once per column under `TileCells`'
  lock and costs 1.5 s of 29 (~5 %), not a plateau;
* the "serial tail" (assembly, `_fillcoo!`, `_writechunk!`) is **0.07 %** here,
  because the destination is one chunk written by one task;
* what replaced them is a single new plateau — the **source-side tree walk** —
  which that campaign could not see, since its source was a `RasterGrid` with a
  memoised `MemoRasterTree` (`raster_tree_memo.jl`), not a `DGGSpace` over a
  `PartialGrid` of CopernicusDEM tiles. The §13 source-path work made a
  `DGGSpace` source *usable*; it did not give it the tree memo the raster source
  already has.

New relative to that campaign, and measured to be irrelevant: subzone write
path (0.013 %), synthetic tile decode (0.02 %), mask sampling (inside that
0.02 %), TiledDEM churn (nil — 14 `loads`, 0 evictions per 3 columns).

---

## 5. Q3 — The 33 k → 25.5 k per-busy-core question

The drop is largely an artefact of how the two numbers were formed, and the
residual is the shared box.

| figure | value | what it divides by |
|---|---:|---|
| la-choice La = 5, `box:10,10.05,45,45.05`, `-t 8` | 33 026 cells/s/core | warm wall × 2.36 "util of 8" |
| §15.9.1 production | 25 400 | 585 000 ÷ 23 instantaneous process cores |
| production, done log | **29 688 cells / worker-second** | Σ per-column wall over 1 930 columns |
| production, process CPU (20:12Z) | **28 614 cells / core-second** | 636 379 ÷ 22.24 |
| this campaign, controlled `-t 1`, same code, same box | **28 277 cells / core-second** | one thread, no workers |

* **The 25.5 k figure charges the GC threads and the write/log path to compute.**
  The like-for-like number is 28.6 k, and a single thread doing nothing else on
  this box today gets 28.3 k. So production is at **99 % of the achievable
  single-thread rate** — there is no production-specific loss.
* **Polar mix: not a cause.** polar_S runs at 32 741 cells/s/core, *better* than
  mid-latitude, because more of its columns are partly unfed.
* **Coastal NaN handling: not a cause.** The done log shows 30–31 k cells/s at
  every NaN fraction from 0 % to 80 %; only the fully empty columns differ
  (68 k), and they are 16 % of columns and pull the average *up*.
* **Write contention: not a cause.** 4 ms per column, disjoint files.
* **Tile-cache misses: not a cause.** 0.02 % of CPU in the whole source path.
* **GC at 63 threads: minor.** 4.5–4.8 % measured at `-t 1`; `--gcthreads=8,1`
  is a reasonable setting for 3.2 GiB of churn per column × 21 workers, and the
  GC threads are part of why "process cores" (22.2) exceeds "worker cores"
  (21 × 1.06 = 22.3 — they coincide).
* **What is left (28.3 k vs 33.0 k, −14 %)** is the box: la-choice ran at 1-min
  load 19.9, this campaign at 47–68, on a 2-socket NUMA machine where 40+
  foreign threads share the LLC and memory controllers. The live run's busiest
  thread never exceeds **0.79 core-seconds per wall second** — it is descheduled
  21 % of the time — and its allocation rate (3.2 GiB/column) makes it
  bandwidth-sensitive. Secondarily, la-choice's box selected a source holding of
  one or two tiles against production's 26 475, so its column paired against a
  single source chunk (production averages 5, of which 2–6 are empty).

**Verdict: no regression. The hot path costs what it costs; the difference is
the shared machine and the accounting.**

---

## 6. Q4 — Ranked fixes, with projected wall clock at a 24-core budget

Baseline for every projection: **28.6 k cells per core-second** measured, and
**61 758 columns × 823 543 = 5.086 × 10¹⁰ cells** left at 20:12Z (4 420 of
66 178 done). Current: 22.24 cores → 636 k cells/s → **22.2 h**.

### (a) Config-only, applicable to a resume of the live run

| # | change | measured basis | effect | new ETA |
|---|---|---|---|---:|
| a1 | `workers=24` (from 21) | cores/worker 1.06 is stable while W ≪ nthreads; W-sweep is linear until W→nthreads | 22.24 → ~25.4 cores, +14 % | **19.5 h** |
| a2 | `workers=23` | as above, exactly the 24-core budget | 22.24 → ~24.4 cores, +10 % | **20.2 h** |
| a3 | `-t 26 --gcthreads=4,1` with `workers=23` | 37 of 63 threads are never runnable; `_wavesize` caps at `min(nt, nsrc, fits)` and nsrc ≤ 24 for 99.9 % of columns | neutral throughput, less scheduler/GC-thread overhead, lower RSS | 20.2 h |
| a4 | `batch`, `cache`, `budget`, `stripes` | tile path is 0.02 % of CPU; `budget` < 2.63e8 would *split* the destination into several tiles and cost more | **neutral or harmful — do not change** | — |

A resume costs ~20 s of startup plus up to one partial column per worker
(≈ 21 × 28 s ≈ 10 core-minutes). a1/a2 buy 2.0–2.7 h. **This is the whole
config-only budget; nothing else in the conf file moves the number.**

### (b) Script changes (`scripts/copdem_production.jl`)

| # | change | effect | new ETA (24 cores) |
|---|---|---|---:|
| b1 | drop `@with GR.OUTER_PARALLEL => true` from the worker body | **≈ nil.** `_fillwave!` re-imposes it for every column with > 1 source chunk (95 % of them). Measured: 1.257 vs 1.268 cores. Only worth doing *together with* c3 | 20.6 h |
| b2 | set `workers` from a core budget rather than a thread count (`W = ceil(budget / 1.06)`) | makes a1 the default and stops the `threadsper` fiction | 19.5 h |

### (c) Library changes — where the real win is

| # | change | evidence | CPU cut | new ETA (24 cores) |
|---|---|---|---:|---:|
| **c1** | **Memoise source-tree node extents per DEM tile.** Give `DGGSpace`-over-`PartialGrid` sources the equivalent of `raster_tree_memo.jl`'s `MemoRasterTree`: cache `STI.node_extent` for the interior nodes of a source chunk's `BlockCursor` tree, keyed by chunk, shared across columns and workers. A tile's interior tree is a few thousand caps (~100 KB); 26 475 tiles is ~2.6 GB, or an LRU of the hot few hundred | `node_extent` 47.7 % inclusive, `_box_cap` 41.9 %; interior-only ≈ 26 % | **−26 %** | **15.3 h** |
| **c2** | **Stop allocating a fresh `Vector{Tuple{Int,Cap}}` per leaf in `child_indices_extents`.** Fill a per-task reusable buffer (`STI` allows returning any indexable) | 71.6 % of 3.22 GiB/column is exactly this type; `GenericMemory` 7.5 % self, `memmove` 1.2 %, GC 4.5 % | **−9 %** | **14.0 h** (with c1) |
| c3 | **Make `_wavesize` estimate real work, and don't set `OUTER_PARALLEL` inside `_fillwave!` when the wave is narrow.** Fall back to CR's inner threading, which measures 8.7× on 12 threads at 92 % per-core efficiency, versus the wave's 1.2× | §3.4 | no CPU cut; lets W drop from 21 to ~8 for the same 24 cores | 14.0 h, **RSS ~3× lower**, and the only path past 24 cores |
| c4 | Cache leaf caps per tile too (1.44 M caps ≈ 46 MB/tile, small LRU, hot within a batch of 8 contiguous columns) | `child_indices_extents` 21.3 % | −10 to −15 % more | ~12 h |
| c5 | Prune zero-pair source chunks before the wave (a cheap cap-vs-cap reject) | 2–6 of 4–10 chunks yield nnz = 0 but cost only ~0.01 s each | **< 0.5 % — not worth it** | — |

Ordering rationale: c1 and c2 are the only changes that cut *work*, they are
independent, and both sit in code the destination side already has a precedent
for (`CapCachedTree`, `MemoRasterTree`). c3 changes the *shape* — at a 24-core
budget it is worth nothing in wall clock but it is what makes the design honest
and it drops the 35 GiB RSS, so it belongs before c4.

### Projected ETA summary for the remaining 5.086 × 10¹⁰ cells

| configuration | cells/s | ETA |
|---|---:|---:|
| as running (22.24 cores) | 636 k | **22.2 h** |
| + a1 `workers=24` (25.4 cores) | 727 k | 19.4 h |
| a2 at exactly 24 cores | 687 k | 20.6 h |
| 24 cores + c1 | 928 k | 15.2 h |
| 24 cores + c1 + c2 | 1 020 k | **13.9 h** |
| 24 cores + c1 + c2 + c4 | ~1 180 k | ~12.0 h |

---

## 7. Reproducing

Scratch driver (not committed; regenerate from this note):

```bash
S=<scratch>
head -n -1 /home/asinghvi17/geo/DGG-subzone-store/scripts/copdem_production.jl > $S/prod.jl
# setup.jl sets ARGS to the production conf with store=<scratch>, includes prod.jl,
# and exports lazycolumn / timed_column / fanout / withcpu / emit.
cd /home/asinghvi17/geo/DGG-subzone-store
RASTERDATASOURCES_PATH=/home/asinghvi17/geo/DiscreteGlobalGrids.jl/bench/data \
  julia --project=bench -t 12 $S/e7_sweep.jl
```

| script | what it produced | `exp` |
|---|---|---|
| `e1_fanout.jl` | latitude bands, representative batches, per-column tiles/chunks/wave | `fanout` |
| `e2_dist.jl` | fan-out over 3 000 sampled columns | `fanout_dist`, `fanout_band` |
| `e3_phase.jl` | phase timing + flame graph, `-t 1` and `-t 12`, `OUTER` on/off | `phase`, `profile` |
| `e4_block.jl` | per-source-chunk block build times and nnz | `blocktimes` |
| `e5_inner.jl` | CR inner threading on/off, cores and wall | `nowave` |
| `e6_hot.jl` | topmost-inclusive hot-path breakdown | `hotpath` |
| `e7_sweep.jl` | W-sweep through the real `runcolumns` | `wsweep` |
| `e8_alloc.jl` | `@timed` baseline + `Profile.Allocs` by type | `alloc_baseline`, `alloc_profile` |

Live-run observations came from `/proc/2915199/task/*/stat` deltas and the run
log; nothing wrote to the live store, log, done log or column cache.
