# 2026-08-19 — Thread utilization during `regrid`: where the parallelism goes

Measures how much work each thread actually does during `regrid`, phase by
phase, and explains the maintainer's observation of **~135 % CPU at
`--threads 8`**. Measurement only; no code changes shipped.

## Trees, machine, method

- **DGG worktree**: `a8fec729fea445497e8a7857ea9035178adad514` (origin/main,
  "Merge pull request #40 … cr-main-repin"), checked out fresh under the
  session scratchpad and pinned for the whole campaign; the user's checkout
  was not modified except to add the three files this campaign delivers.
- **CR arm `main`**: `89ec5ca` (v0.2.9, the rev the repo pins). **CR arm
  `frontier`**: `508a637` (branch `claude/budget-frontier`, PR #137 HEAD).
  GeometryOps main (git-tree `c30ec9107`, v0.1.44), GeometryOpsCore v0.1.11.
  Julia 1.12.6; GC threads = nthreads at every thread count (default, on).
- **Machine**: the shared 64-core box, load average 19–35 throughout
  (recorded per stage row in the NDJSON). No `taskset`: with 64 cores and
  ~25 load there were always free cores, and the metrics below are
  contention-robust (process/thread CPU clocks and deterministic work
  counters, not wall-clock exclusivity). At most two Julia suites ran
  concurrently, staggered.
- **Instrumentation**: a probe module compiled into private clones of both CR
  arms and into the worktree's GlobalRegridding (full diffs in the logs
  tarball; the shared `CRpr` clone and the user's checkout were untouched).
  Phases record wall, `CLOCK_PROCESS_CPUTIME_ID`, `CLOCK_THREAD_CPUTIME_ID`
  and GC-time deltas; every spawned traversal/clip task records pairs
  handled, weights stored, thread-CPU time and start/end thread ids. **Zero
  task migrations** were observed across all ~200k instrumented tasks, so
  per-task thread-CPU readings are valid. Overhead: the mean instrumented
  task costs 769 µs of CPU against ~1 µs of probe work; instrumented e512
  build wall (2.6–3.0 s) matches the uninstrumented post-fix campaign
  (3.37 s under heavier load). No MCP Julia server was available, so each
  (arm × thread count) ran as one script process; first-call vs steady-state
  is separated by repetition number inside each process.
- **Workloads** (bench/harness.jl conventions): GLO-30 tile N45_00_E010_00
  windows → IGEO7 level-13 covering, Conservative + Weighted(0.5).
  `eN` = in-memory raster, eager (whole-domain block). `lN` = the same
  window behind a 128×128-chunked CountingDisk, lazy read of the full
  destination, budget 2^30, default destination tiling (which is ONE tile at
  ≤2048² — the known `chunkcells` defaults issue). `dN` = DGG-native
  CopernicusDEM source, eager. `moc` = the whole 3600² tile, lazy, 512×512
  source chunks, defaults otherwise.
- **Raw data**: `2026-08-19-thread-utilization.ndjson` beside this file
  (26,950 rows: stage rows, phase rows, per-task summaries, cold-process
  rows, profile frames) and `2026-08-19-thread-utilization-logs.tar.gz`
  (full per-task NDJSON, run logs, `/usr/bin/time -v` outputs,
  instrumentation diffs, runner and analysis scripts).

## 1. Utilization curve — was 135 % expected? (Q3)

Warm steady-state, per stage; `util` = process CPU time / wall.

| workload | arm | t1 | t2 | t4 | t8 |
|:--|:--|--:|--:|--:|--:|
| e512 (plan) | main | 16.80 s · 1.00 | 8.85 s · 1.96 | 4.69 s · 3.70 | **2.60 s · 7.12** |
| e512 (plan) | frontier | 16.70 s · 1.00 | 10.99 s · 1.58 | 7.66 s · 2.33 | **6.58 s · 2.77** |
| l512 (read) | main | 18.77 s · 0.99 | 9.83 s · 1.99 | 5.20 s · 3.83 | 2.84 s · 7.29 |
| l512 (read) | frontier | 18.66 s · 0.99 | 9.96 s · 1.96 | 5.35 s · 3.75 | 2.94 s · 7.23 |
| l1024 (read) | main | 75.81 s · 0.99 | 40.76 s · 1.96 | 22.14 s · 3.77 | 12.17 s · 7.16 |
| l1024 (read) | frontier | 75.64 s · 1.00 | 42.47 s · 1.93 | 22.97 s · 3.63 | 13.09 s · 6.77 |
| e1024 (plan) | main | 61.23 s · 1.00 | 31.39 s · 1.96 | 17.03 s · 3.70 | 9.24 s · 6.92 |
| e1024 (plan) | frontier | 60.75 s · 1.00 | 35.37 s · 1.81 | 26.14 s · 2.45 | **23.32 s · 2.84** |

DGG-native source at t8: d512 plan 18.2 s · util 7.72 (main), 19.2 s · 7.03
(frontier); d1024 73.9 s · 7.70 (main), 77.1 s · 7.12 (frontier). MOC scale
is §6 (util ~3.9–4.2, a different bottleneck).

**Verdict on the curve.** On the pinned tree the hot paths are healthy: the
main arm runs at 690–730 % CPU at 8 threads on every workload from 512² to
1024², speedup 6.2–6.6× over t1. **135 % is not what this tree does in
steady state.** The regimes that do produce it are in §5.

**The frontier arm is a real regression on eager raster→DGG builds**: eager
utilization collapses to 230–280 % and e1024 takes 2.5× the main arm's wall
(23.3 s vs 9.2 s). Its lazy and DGG→DGG cases are unaffected (§3). The
user's 135 % predates the frontier, so this is a finding about the intended
future default, not the explanation of the observation.

## 2. Phase decomposition at 8 threads (Q1)

Eager e512, warm (rep 3): one `plan_regrid` + `regrid`, phases sequential in
wall time; util = phase pcpu / phase wall.

| phase | main wall | main util | frontier wall | frontier util | serial by construction? |
|:--|--:|--:|--:|--:|:--|
| dst/src `subtree` | ~0.000 s | – | ~0.000 s | – | yes (O(1) cursors) |
| candidate query (`cr_query`) | 1.362 s | 7.9 | 5.259 s | 2.0 | no — see §3 |
| — query spawn walk / frontier descent | 0.278 s | overlaps tasks | 0.001 s | 1.0 | yes |
| — query fetch | 1.044 s | 7.8 | 5.236 s | 2.0 | no |
| — query vcat merge | 0.040 s | 1.0 | (in fetch) | – | yes |
| clip + COO assembly (`cr_clip`) | 1.077 s | 6.5 | 1.117 s | 6.3 | no (4×nt chunk tasks) |
| — COO vcat merge | 0.007 s | 1.0 | 0.007 s | 1.0 | yes |
| CR `sparse()` | 0.026 s | 1.0 | 0.067 s | 1.0 | yes |
| `_fillcoo!` | 0.052 s | 1.0 | 0.108 s | 1.0 | yes |
| WeightBlock `sparse()` | 0.082 s | 1.0 | 0.030 s | 1.0 | yes |
| **plan total** | **2.60 s** | **7.12** | **6.58 s** | **2.77** | |
| apply (`applyplan!`) | 0.004 s warm, 0.17 s first | 1.0 | same | 1.0 | yes |

e1024 (main, t8) has the same shape: query 4.78 s · 7.7, clip 3.94 s · 6.6,
serial tail (`_fillcoo!` 0.27 + two `sparse()` 0.25 + merges 0.15) ≈ 0.67 s
of 9.24 s. At t1 the same e512 build splits query 10.3 s / clip 6.1 s — the
candidate query is ~60 % of the build even serially.

On the campaign question about `assemble_sparse_matrix_coo`: **it does
parallelize**, over candidate pairs — ChunkSplitters partitions the pair
list into 4×nthreads chunks (32 at t8), one StableTask each, measured util
6.3–6.6. Only its final vcat merges and the `sparse()` calls are serial.

Lazy l1024, t8, main (read wall 12.17 s):

| piece | wall | util | note |
|:--|--:|--:|:--|
| `wave_fill` (8 waves × 8 builds) | 12.06 s | 7.05–7.30 per wave | 99.1 % of the read |
| serial residue (`apply_group` 0.055 s + `add_ref` 0.041 s + `src_read` 0.010 s + `write_chunk` 0.008 s + `connected_src` ≈0) | 0.11 s | 1.0 | 0.9 % |

Inside the waves each of the 64 block builds nests the threaded CR path, so
the wave keeps all 8 threads busy even where wave width < 8. The candidate
query is 74 % of build task-seconds (62.2 of 84.4 s) because **each source
chunk re-queries the same whole-tile destination tree** — at 1024² that is
~1.7× the eager path's total query CPU: a throughput tax, not a utilization
one. The `TileCells` geometry cache never engaged in any measured case (the
single default tile always exceeds the 2^16-cell cap), so destination
polygons are re-synthesized per candidate pair inside the threaded clip —
which is why t1's 37 % "dstgeom" attribution does NOT become a serial
bottleneck at t8: it lives inside the parallel region.

**Serial-by-construction inventory**: plan side — spawn recursion / frontier
descent, both vcat merges, both `sparse()` calls, `_fillcoo!`; apply side —
all of `applyplan!`; lazy driver — source reads, `applyblock!` accumulation,
`_writechunk!`, connected-chunk discovery; and outside `regrid` proper, the
destination covering query and `DGGSpace` construction. At t8 these sum to
≈8 % of e512 plan wall and ≈1 % of l1024 read wall.

## 3. Per-task and per-thread accounting (Q2)

Clip tasks (identical on both arms; ChunkSplitters even split):

| case (t8, warm) | n tasks | pairs/task min/med/max | pairs CV | tcpu med/max | per-thread tcpu |
|:--|--:|:--|--:|--:|:--|
| e512 | 32 | 58,981 / 58,981 / 58,982 | 0.00 | 194 / 335 ms | balanced |
| e1024 | 32 | 236,090 / … / 236,091 | 0.00 | 754 / 915 ms | balanced |
| l1024 (64 builds) | 2,048 | 3,683 / 3,690 / 3,698 | 0.00 | 12.7 / 133 ms | balanced |

Pair counts are exactly even; per-pair clip cost varies ~1.7× max/median,
which 4-chunks-per-thread granularity absorbs.

Query tasks — the arms diverge completely:

| case (t8) | arm | n tasks | tcpu sum | tcpu max | tcpu CV | max task's pair share | per-thread tcpu min/max |
|:--|:--|--:|--:|--:|--:|--:|--:|
| e512 | main | 13,288 | 10.22 s | 0.036 s | 1.45 | 0.1 % | 1.27 / 1.31 s |
| e512 | frontier | 117 | 10.44 s | **5.20 s** | 6.15 | **50 %** | 0.14 / 5.20 s |
| e1024 | main | 5,939 | 35.85 s | 0.114 s | 1.47 | 0.4 % | 4.43 / 4.55 s |
| e1024 | frontier | 141 | 36.27 s | **18.13 s** | 6.78 | **50 %** | 0.30 / 18.26 s |

The main arm's recursive spawn produces thousands of small tasks whose
per-thread CPU balances to within 3 %. The frontier planner stops at 117–141
pairs believing them balanced: **every pair's estimated weight is nearly
identical (estw ≈ 2.6–2.8×10⁵, CV 0.02) while actual work spans five orders
of magnitude** (30 → 943,951 pairs; corr(estw, tcpu) ≈ 0.20). The realized
pair sizes halve down a chain (50 %, 25 %, 12.5 %, 6.3 % of all pairs): the
planner split the wrong pairs and left one branch holding half the search.
`split_weight` itself is answered O(1) by both trees (ncells fallback); the
failure is `pair_weight`'s cap-overlap × cell-density model, which saturates
on fine RasterGrid × DGGSpace whole-space pairs whose SphericalCaps nearly
coincide. DGG→DGG (d1024 util 7.12) and the lazy path's many small
per-chunk queries are unaffected — the failure needs one huge foreign-tree
pair, exactly the eager raster→DGG case. The 2026-08-18 threading-policy
campaign's coarse `dgg_n144` case never exercised this.

## 4. Attribution of the gap (Q4)

Amdahl account of e512 plan at t8, main arm (wall 2.60 s, CPU 18.5 s, util
7.12): ideal wall at 8 threads is 2.31 s; the 0.29 s (11 %) gap is ≈0.18 s
serial-phase excess (the 0.21 s serial tail minus its 1/8 share — merges,
two `sparse()` calls, `_fillcoo!`, spawn walk = 8 % of wall), the remainder
GC stop-the-world and scheduler idle. In-phase GC time totals 0.46 s (18 %
of wall; GC threads on, = nthreads). The CPU profile agrees:
`GenericMemory` allocation is 16–17 % of CPU samples, `poptask`/
`try_yieldto` (idle scheduler) 17–20 %, and no runtime-dispatch hotspots in
the clip loop. l1024 lazy at t8: 99.1 % of wall inside threaded waves at
util 7.05–7.30, 0.9 % serial driver, GC 2.35 s of 12.17 s (19 %).

**Ranked levers:**

1. **Frontier arm: fix `pair_weight` or re-split reactively** — worth 2.5×
   wall on eager raster→DGG builds (e1024 23.3 s → 9.2 s; util 2.8 → 6.9).
   As measured, the frontier would ship a regression exactly where DGG
   lives. (Largest lever overall, but only on the frontier arm.)
2. **Allocation/GC pressure** (both arms, every workload): 16–17 % of CPU
   samples allocate and 18–20 % of wall is GC. The per-pair `push!` COO
   vectors and per-task result vectors are the sites. Worth roughly one
   thread of the eight at t8.
3. **MOC-scale defaults** (§6): budget-derived tiles rebuild a
   `CellCapTree` over ~3.3M destination cells 320 times (2,856 s of task
   wall) and the `_wavesize` weight-budget bound caps the wave at 4
   concurrent builds — together they hold the whole-tile default path at
   util ~3.9 and ~900 s against 140–175 s measured previously with
   `chunklevel = 7`. Chunk-aligned tiles (or a cached/cheaper restricted
   subtree, or a bigger weight budget) recover both.
4. **Serial tail of the eager build** (~8 % of plan wall at t8): bounded win
   (≤0.7 s at 1024²) but it caps utilization at ~7.3/8 as builds shrink.
5. **Lazy repeated-query overhead** (throughput, not utilization): one
   destination-tree query per source chunk costs ~1.7× eager query CPU at
   1024²; a shared candidate index or coarser source chunks recovers it.
6. **Not levers** (measured, killed): destination geometry synthesis (lives
   inside the threaded region at t8); `assemble_sparse_matrix_coo` (already
   parallel, exactly even splits); the DGG-native source path (util 7.7);
   GC threads off (they are on); the single-tile lazy default as a
   *serializer* at ≤2048² (waves + nested threading cover it).

## 5. Reproducing the 135 % (Q5)

Whole-process `%CPU` from `/usr/bin/time -v`, one cold `julia -t 8` process
per cell — package load + destination covering + first-call compile + one
regrid, i.e. exactly "run regrid locally once":

| workload | main | frontier |
|:--|--:|--:|
| e128 | **160 %** | **158 %** |
| e256 | 171 % | 162 % |
| l256 | 170 % | – |
| e512 | 217 % | 179 % |
| l512 | 217 % | 215 % |
| e1024 | 309 % | 213 % |
| l1024 | 337 % | 332 % |

Utilization of the **regrid call alone** inside those processes (space
setup + plan + apply, excluding Julia startup and package load):

| workload | main | frontier |
|:--|--:|--:|
| e128 | **1.30** | **1.29** |
| e256 | 1.79 | 1.47 |
| l256 | 1.69 | – |
| e512 | 2.99 | 1.97 |

The structure: the destination covering query + `DGGSpace` construction +
its first-call compile is a fixed ~5 s at util ~1.07 (serial), while the
parallel weight build shrinks with the workload (0.7 s at 128², 3.2 s at
512²). At ≤256² the serial setup dominates and the whole call measures
**129–179 % CPU at 8 threads**; the 135 % observation lands squarely in
this band. Two further depressants on a maintainer's laptop: a first-ever
run also pays package precompilation at ~100 % for minutes (this box had
warm caches), and the pre-fix tree (eager 512² was 51 s before the post-fix
campaign) had proportionally more serial work.

**Verdict**: ~135 % at 8 threads is the signature of a **first-call or
small-workload invocation, dominated by the serial destination-space setup
and compilation** — not of the steady-state regrid kernel, which measures
690–730 % on this tree. Closest single reproductions: cold e128 (call-only
1.29–1.30, whole-process 158–160 %). Nothing in the warm suite comes near:
the minimum warm t8 utilization measured is 2.77 (frontier eager, its own
regression) and 6.77 for any main-arm case at ≤1024². (At MOC scale the
warm default path does drop to 3.85–4.21 — §6 — so "watched `top` during a
whole-tile lazy run with default chunking" is the one warm scenario in the
same neighborhood, though it reads ~390 %, not 135 %.)

## 6. MOC-scale workload (3600², lazy, t8, defaults)

Read wall **901 s (main) / 840 s (frontier)** at util **3.85 / 4.21**, GC
176/180 s (~20 %). Structure: 16.49M destination cells → **5 budget-derived
tiles** (the single space chunk is unusable, so `_desttiling` falls to the
budget path); 64 source chunks; 320 block builds in 77 waves → wave width
≈4, exactly the `_wavesize` weight-budget bound (`fits` =
weightbudget/floorbytes ≈ 4 at 3.3M-cell tiles). Dominant cost:
`subtree(dst_space, inds)` falls back to `CellCapTree` for every build,
because budget tiles are not chunk-aligned — **320 rebuilds of a
~3.3M-cell cap tree, 2,856 s of summed task wall** (11,055 CPU-s), i.e.
per-cell polygon+cap synthesis repeated 64× per tile. The prior harness
sweep's 139–175 s full-tile figure used `chunklevel = 7` destination
chunking, whose chunk-aligned tiles reuse the grid hierarchy in O(1); the
default path measured here is ~5–6× slower and holds utilization near 4.
This is the same defaults issue the postfix campaign flagged, now with the
utilization mechanism attached: half the loss is the wave-width budget
bound, half is that the redundant `CellCapTree` builds do not fill all
eight threads.

## Caveats

- Shared box: wall times carry load noise (loadavg per row in the NDJSON);
  utilization ratios were stable across repetitions (±0.15 at t8).
- The instrumented trees differ from the pinned trees only by the probe
  diffs in the tarball (timing wrappers and counters; identical algorithms;
  identical outputs across arms and thread counts by digest).
- The 135 % reproduction is a band, not a point: the exact figure depends
  on core count, precompile-cache state, and workload size; this box
  brackets it at 129–179 % for cold ≤512² eager invocations.
