# 2026-08-19 — Post-stack profile of the eager MOC regrid: where the 50.8 s go

Four-worker profiling campaign on the post-stack tree, answering the
maintainer's observation of a CPU curve that reads **~750 % → ~600 % → ~100 %**
in `top` during a single warm `regrid`, and ranking what to fix next.
Measurement only; no code changes shipped.

## Trees, machine, method

- **DGG**: `1bb575d1852d37b84079370bd97147337b577b6b` (main), checked out in a
  private worktree under the session scratchpad; the user's checkout was not
  modified except to add the two files this note delivers.
- **ConservativeRegridding**: pinned to `claude/budget-frontier`,
  git-tree-sha1 `873cc732e77638eb44d7d92e06880585a31200a9`, resolved commit
  `6a4b997ab2e66dea45b5b62b5b2f32f0a3b279b0` (v0.2.9). **This pin does not
  carry the pair-gen threading regression** measured at `508a637` in the
  2026-08-19 thread-utilization campaign (§3 of that note): the dual-tree
  descent spreads across all 8 workers here.
- **Julia 1.12.6**, `-t 8` with default `--gcthreads=8,1` unless a row says
  otherwise; `-t 64` probe as marked.
- **Machine**: the shared 64-core box, 1-minute load **23–46** throughout the
  t8 work (51–62 during the t64 probe), recorded per run in the NDJSON. No
  `taskset`; metrics are process/thread CPU clocks, `GC_Diff` counters and
  `/proc` sampling, not wall-clock exclusivity.
- **Workload** (identical for all four workers): a 3600×3600 CopernicusDEM
  GLO-30 tile as an in-memory raster, **eager** conservative regrid onto a
  level-13 IGEO7 MOC coverage of **16,181,892 cells**. **93,079,031 candidate
  pairs**, **66,129,312 nnz**. Warm (rep 2/3 in-process); the profiled call is
  always the *second or later* regrid against the same grid object.
- **Baseline**: **50.8–51.8 s** warm at t8 (pre-stack: 102–110 s).
- **Method per worker**: (A) `/proc/<pid>/stat` utilization timeline at 0.5 s
  with phase marks from an instrumented decomposed replay; (B) `Profile`
  sampling of the whole call plus one profile per phase, `delay=0.01`,
  1,136,682 slots, idle/sleeping samples dropped; (C) deterministic
  `GC_Diff`/`@timed` per phase plus `Profile.Allocs` for site attribution, and
  a `gcthreads` × `nthreads` sweep.
- **Raw data** (session scratchpad, not in the repo):
  `scratchpad/utilrun/` (A: `t8-run1/`, `t64-run1/` — `samples.txt`,
  `marks.ndjson`, `plateau.py`), `scratchpad/profB/` (B: `whole.jlprof`,
  `ph_*.jlprof`, `an_*.txt`), `scratchpad/allocC/` (C: `c1`–`c5` run dirs,
  `summary.txt`, `alloc-tables.txt`). Machine-readable summary of this note:
  `2026-08-19-post-stack-profile.ndjson` beside it.

## 1. Consolidated phase table (t8, warm)

Spine is C's deterministic per-phase `@timed`/`GC_Diff` split; A's columns are
independent `/proc` phase-mark readings from a different process; B's wall is
the profiled replay (under `Profile` overhead). Util = phase CPU-s / phase wall.

| phase | wall (C) | CPU-s (C) | util | share of CPU | GiB | GC s | wall (A) | wall (B) |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|
| dst `CapCachedTree` build | 3.73 | 24.82 | **6.65** | 8.1 % | 0.72 | 0.00 | 3.68 | 3.90 |
| candidate-pair gen (dual-tree) | 19.55 | 127.14 | **6.50** | 41.5 % | 15.34 | 2.54 | 19.98 | 20.16 |
| clip + COO assembly | 20.20 | 144.04 | **7.13** | 47.0 % | 11.69 | 1.62 | 20.60 | 20.58 |
| CR `sparse()` | 1.62 | 2.47 | 1.52 | 0.8 % | 2.29 | 0.15 | 1.62 | 2.09 |
| `_fillcoo!` | 3.83 | 5.27 | 1.38 | 1.7 % | 5.59 | 0.59 | 3.96 | 4.11 |
| WeightBlock `sparse()` | 1.94 | 2.41 | 1.24 | 0.8 % | 2.41 | 0.25 | 1.95 | 1.86 |
| `applyplan!` | 0.38 | 0.43 | 1.13 | 0.1 % | 0.42 | 0.02 | 0.37 | 0.38 |
| **sum of phases** | **51.25** | **306.58** | 5.98 | 100 % | **38.46** | **5.15** | | |
| monolithic warm call | 50.82 | 308.39 | 6.07 | | 38.58 | 5.43 | 50.76–51.80 | 53.21 |

Also measured but too small to table: `subtree(src)` 0.018 s, operator
construction 0.021 s.

**Reconciliation.** Decomposed/monolithic = **100.9 % of wall, 99.4 % of CPU**
— the split is exact and the instrumentation costs ~0.4 s. Phase-by-phase A vs
C: **−3.7 % … +3.4 % on wall**, all within tolerance. B vs C: within ±7.3 % on
every phase except `CR sparse()` (+29 %, see flags). B's total allocated bytes
41.43 GB and C's 38.581 GiB are **the same number** (38.581 GiB = 41.42 GB),
measured by two independent mechanisms; B's GC 5.30 s vs C's 5.43 s agree to
2.4 %. Sum of phase CPU-s **306.6 vs 308.4 monolithic (99.4 %)**.

## 2. The plateau narrative: 750 → 600 → 100

A's `/proc` utilization timeline (0.5 s samples, warm rep 2; rep 3 identical to
within a sample):

| segment | window | wall | CPU-s | mean %CPU | sample range |
|:--|:--|--:|--:|--:|:--|
| 1 dst cap-tree build | 0.0–4.0 s | 4.0 | 27.4 | **686 %** | 98 → **802 %** |
| 2 candidate-pair gen (core) | 4.0–19.8 s | 15.8 | 113.3 | **717 %** | 469–781 % |
| 3 pair-gen wind-down + `reduce(vcat)` | 19.8–23.4 s | 3.6 | 14.9 | **415 %** | **106 %** trough |
| 4 clip + COO assembly | 23.4–41.9 s | 18.5 | 138.2 | **747 %** | 298–802 % |
| 5 serial tail | 41.9–51.8 s | 9.9 | 14.2 | **143 %** | 98–298 % |
| whole call | | 51.8 | 308.0 | **595 %** | |

**The three numbers the maintainer saw.** *~750 %* is the two genuinely
parallel plateaus — the destination cap-tree build touches **802 %** at its
peak and the clip/assembly plateau holds **747 %** for 18.5 s, so any glance at
`top` during the first 4 s or the middle 18 s reads 750–800 %. *~600 %* is
`top`'s smoothed whole-call mean, **595–602 %** measured over both warm reps,
which is what a several-second `top` refresh reports once it has averaged the
plateaus with the trough and the tail (the pair-gen phase's own mean, 622–668 %
including its wind-down, sits in the same band). *~100 %* is the last **9–10 s
of every call**, a strictly serial tail — `sparse()` → `_fillcoo!` →
WeightBlock `sparse()` → `apply` — running on exactly one mutator thread; the
readings that oscillate to 143 % mean and 200–298 % peaks are **not** a second
worker but the GC threads: C confirmed this two ways (counters — tail non-GC
CPU 6.75 s over 7.78 s wall = one thread, with 3.81 GC CPU-s over the tail =
3.8 effective GC threads at t8; and direct `/proc/<pid>/task` sampling at t64,
where one thread shows 6.91 CPU-s over 6.94 s and 63 others show 0.24–0.28 s
each).

**The trough.** Between the two big plateaus utilization collapses to **106 %
for ~1.4 s** (segment 3, mean 415 % over its 3.6 s window). B identified it as
the serial `reduce(vcat, ...)` at `ConservativeRegridding/intersection_areas.jl:181`
that merges the per-task candidate-pair vectors: **85 % `memmove`, 1.59 GB
moved at ~1.3 GB/s**.

**What the plateaus are made of** (B's flat profile, whole call, 25,336 active
samples): **zero** `jl_apply_generic`/`jl_invoke` frames anywhere — no runtime
dispatch, no type-instability lever. The top self-time frames are `*`
(12.24 %), `+` (6.65 %), `<` (5.92 %), `muladd` (3.11 %), `-` (2.57 %), `/`
(1.91 %) — **~35 % of all samples are scalar `Float64` arithmetic**. This is a
spherical-geometry FLOP workload, not abstraction overhead. Region rollup:
dual-tree search 42.4 %, clip/assembly 41.2 %, cap build 8.5 %, serial thread-1
3.5 %, GC threads 4.4 %.

**Assembly is co-dominant with pair-gen, not subordinate.** At 144.0 vs 127.1
CPU-s (47.0 % vs 41.5 %) the clip/assembly phase is the single largest CPU
consumer, correcting the earlier campaign's "the candidate query dominates"
claim — which held at 512²–1024², where the query is 60–74 % of build CPU, but
not at MOC scale.

**Where the missing 2 threads went.** Whole call: 50.82 s wall against
306.58 CPU-s, i.e. an ideal 8-thread wall of 38.32 s and a **deficit of
98.2 idle-thread-s**:

| bucket | idle-thread-s | share |
|:--|--:|--:|
| serial tail (7.77 s × 7 threads) | 54.4 | 55 % |
| serial `vcat` trough (excess over its own CPU/8) | 13.9 | 14 % |
| stop-the-world GC inside the parallel regions | 17.4 | 18 % |
| scheduler idle, load imbalance, box contention | 12.5 | 13 % |

**Serial regions are 69 % of the gap to ideal utilization**; GC is real but
second-order — excluding every GC window, pair-gen still only reaches ~6.3 of
8 threads.

## 3. Scaling: t8 vs t64, and the `gcthreads` sweep

| phase | t8 wall | t8 util | t64 wall | t64 util | speedup |
|:--|--:|--:|--:|--:|--:|
| dst cap-tree build | 3.73 | 6.65 | 1.41 | 17.7 | 2.6× |
| candidate-pair gen | 19.55 | 6.50 | 7.53 | 21.0 | 2.6× |
| clip + COO assembly | 20.20 | 7.13 | 6.16 | 26.9 | 3.3× |
| CR `sparse()` | 1.62 | 1.52 | 1.59 | 3.2 | 1.0× |
| `_fillcoo!` | 3.83 | 1.38 | 4.00 | 3.4 | 1.0× |
| WeightBlock `sparse()` | 1.94 | 1.24 | 2.08 | 3.4 | 0.9× |
| `applyplan!` | 0.38 | 1.13 | 0.70 | 4.1 | 0.5× |
| **parallel region** | **43.48** | | **15.10** | | **2.88×** |
| **serial tail** | **7.77** | | **8.37** | | **0.93×** |
| whole warm call | 50.8–51.8 | 5.95–6.07 | 22.3–23.6 | 15.6–16.9 | **2.3×** |

**8× the threads buys 2.3× the wall.** The parallel region scales 2.9× (not 8×;
the box is shared and the phases are memory-bandwidth-heavy), and the serial
tail is **exactly invariant** — 7.77 s at t8, 7.78 s at t64 — so it grows from
**15 % of the call at t8 to 37 % at t64**. The tail's apparent 301 % at t64 is
GC threads: non-GC tail CPU 6.75 s (one thread) plus 16.19 GC CPU-s over 1.03 s
of GC wall = **15.7 effective GC threads**; the same arithmetic at t8 gives 3.8,
which is exactly A's 131–146 % reading.

`gcthreads` sweep (2 warm reps each unless noted):

| config | wall s | CPU-s | util | GC s | GC % | pauses | ms/pause |
|:--|--:|--:|--:|--:|--:|--:|--:|
| t8 default (8,1) | 50.98 | 304.5 | 5.97 | 5.05 | 9.9 | 131 | 38.5 |
| t8 `4,1` | 51.67 | 310.9 | 6.02 | 4.29 | 8.3 | 116 | 36.8 |
| t8 `1,0` (serial GC) | **55.54** | 298.4 | 5.37 | 9.54 | 17.2 | 80 | 119.2 |
| t8 `64,1` | 50.02 | **388.4** | 7.76 | 3.91 | 7.8 | 148 | 26.4 |
| t64 default (64,1), 1 rep | 22.33 | 377.1 | 16.89 | 5.59 | 25.0 | 102 | 54.8 |
| **t64 `--gcthreads=8,1`** | **20.83** | **328.4** | 15.76 | 3.55 | 17.1 | 78 | 45.9 |

At t8 **the default is already right**: serial GC is a real 4.5 s regression,
and `64,1` buys 1 s of wall for +84 CPU-s of pure GC-thread burn (its 7.76
"utilization" is fake). At t64 `--gcthreads=8,1` is a **genuine win: −7 % wall,
−36 % GC time, −49 CPU-s**. Note the GC% doubling from t8 to t64 is arithmetic,
not a regression: **absolute GC wall is identical** (5.59 s vs 5.43 s) while
non-GC wall shrinks 45.4 → 16.7 s. GC does not scale with mutator threads
(mark/sweep is bandwidth-bound) and per-pause cost *rises* with mutator count
(26.4 ms/pause at t8/64gc vs 54.8 at t64/64gc — safepoint synchronization of 64
mutators costs ~2×).

## 4. Allocation and GC state

**41.43 GB churn for 93.1 M candidate pairs = 445 B/pair**, against 2440 B/pair
and 221.17 GB in the pre-stack campaign — **5.34× lower for identical work**.
185,076,141 allocations, GC 5.43 s (10.7 % of wall), 131 pauses (109
incremental, 22 full), 33.9 GiB freed, max pause 386 ms.

**Destination geometry allocation is now exactly zero.** `getcell(dst_space, i)`
allocates **0.0 B/cell**: `cell_boundary` returns an isbits `SmallList`, the
polygon stays stack-resident (`fallbacks/geometry.jl:186-193`). The pre-stack
campaign attributed **52.9 % of all bytes (~109 GiB)** to destination geometry;
**none of it survives**. What is left:

| site | bytes | share of sampled | note |
|:--|--:|--:|:--|
| GO `dual_depth_first_search.jl:33` `_child_extents` | ~7 GB | 42.0 % of pair-gen | `Memory{Tuple{Int,SphericalCap}}` per descent node — **upstream GeometryOps** |
| DGG `cap_cached_tree.jl:43` | — | 29.0 % of pair-gen | cap bookkeeping |
| GR `raster_tree_memo.jl:82` | — | 14.6 % of pair-gen | |
| `rastergrid.jl:444` `_cellring` + `:413` `getcell` | ~6 GB | 73.9 % of assembly | **source pixel** geometry, 240 B/cell heap |
| CR `conservative.jl:315-317` COO `push!` growth | ~3–5 GB | ~26 % of assembly | no `sizehint!` |

Pair-gen's 16.47 GB is **100 % cap/extent bookkeeping with zero polygons**;
assembly's 12.55 GB is ~100 % source-pixel geometry plus COO growth, with
**zero samples from Sutherland-Hodgman or IGeo7 boundary code**. The serial tail
holds **27.8 % of all bytes in 528 allocations** — large-buffer traffic, not
object churn.

Per-phase GC sample share (B): `_fillcoo!` 6.4 % (worst), pair-gen 4.3 %,
assembly 2.4 %, both `sparse()` calls 0.0 %.

## 5. Anatomy of the two fixable regions

**The serial tail (9.4 s), and the decisive redundancy.** B established that
`nnz_coo == nnz(block) == 66,129,312` exactly: **the dual-tree emits each pair
once**, `sparse()`'s dedup merges zero entries, and with the identity indexmap
of the eager whole-block path `block` **already is** the final weights matrix.
The pipeline therefore runs COO → CSC → COO → CSC and **the middle two steps are
provably redundant on this path** (they are *not* redundant on the chunked
path, where the indexmap is non-trivial — any fix must be path-aware).
Breakdown: `_fillcoo!` (4.11 s) is 91.3 % `addweight!`
(`conservative.jl:370`), of which **53.7 % is `push!` → `_growend!` → `memmove`
(2.21 s of Vector growth with no `sizehint!`)**, 22.5 % the actual store,
15.7 % the CSC read loop, 8.1 % GC. CR `sparse()` (2.09 s) is 100 %
`SparseArrays.sparse!`: three passes over 66.1 M entries, dedup/repack ~33 %,
CSR scatter ~24 %, final CSC scatter ~23 %, GC zero. WeightBlock `sparse()`
(1.86 s) is 94.6 % `sparse!` plus 5 % `copy(coo.denom)`.

**Clip + assembly (20.6 s, 144 CPU-s).** 74.8 % is the real clip kernel
(`SphericalClipAreaOperator`): Sutherland-Hodgman clip-to-edge 62.6 %, of which
`spherical_orient` predicates 30.8 % and `robust_cross_product` ~19 %; ring
dedup via `isapprox` 4.4 %; **the area integral itself is only 1.7 %** — the
cost is clipping, not integrating. The other 25 % is geometry materialization
(`_memocell` 19.2 %, ~27 CPU-s), dominated by **destination** cells (15.1 %:
`getcell` → `cell_polygon` → `cell_boundary_cartesian` (`z7grid.jl:205`) →
`dev_to_xyz` (`snyder.jl:359`) = 12.6 %) with source raster cells at 4.1 %. The
64-slot direct-mapped `CellMemo` works: only 19 % of lookups reach `getcell`.

**Cap-tree build (3.90 s @ 767 % = 29 CPU-s) — real work, doubly redundant.**
The split is a perfectly balanced 8-task static partition producing exactly
`ncells` caps; 96.7 % is `_fillcaps!`, of which the `cell_boundary` chain is
79 % → `dev_to_xyz` (Snyder inverse) 68.6 %, whose self time is raw
transcendentals. But: **(a) across calls it is 100 % redundant** — `_leafcaps`
is not memoized on the grid, so every regrid rebuilds it, and the profiled call
was the *second* regrid on the same grid object and still paid the full build;
**(b) within the call it is ~half duplicated with assembly** — the build
computes each cell's 6–7 XYZ boundary vertices, keeps a 4-number cap, discards
the vertices, and assembly then re-runs the identical chain for the same
16.18 M cells (~23 CPU-s). Whole-call, `cell_boundary` is **16.3 % of all
regrid CPU** and `dev_to_xyz` alone is **14.0 %**.

## 6. Ranked levers

Ceilings are measured, at t8, against the 50.8 s warm baseline.

| # | lever | t8 ceiling | t64 relevance | where | decision? |
|--:|:--|--:|:--|:--|:--|
| 1 | **Fuse `_fillcoo!` + WeightBlock `sparse()`** on the eager whole-block path: keep the CR `sparse()` result as the weights matrix, take one row-sum pass for `denom` | **−5.5 s (10.8 %)** | −5.5 s = **26 %** of the 20.8 s t64 wall | CR `conservative.jl` + DGG WeightBlock | **yes** — must be path-aware (invalid on the chunked/non-identity-indexmap path); lives partly in CR, not DGG |
| 2 | **Destination geometry reuse.** (a) memoize `_leafcaps`/`CapCachedTree` on the grid object: **−3.7 s on every repeat regrid** (0 on a first call); (b) keep the boundary rings the cap build already computes and let assembly read them: **−2.7 s** (21.7 CPU-s at util 7.1) | **−3.7 s + −2.7 s = −6.4 s (12.6 %)** on a repeat call | (a) shrinks to −1.4 s at t64 (the build parallelizes, util 17.7); (b) scales with it | DGG `z7grid.jl` / `cap_cached_tree.jl` | **yes** — (b) costs **~2.7 GB f64 / ~0.9 GB f32** of resident cache; (a) needs a cache-lifetime/invalidation policy on the grid object |
| 3 | **Rest of the serial tail**: parallel COO→CSC (bandwidth-bound, 3–4×) −1.15 s; parallel concat replacing `reduce(vcat)` −1.0 s; uniqueness-assuming `sparse` dropping the dead dedup pass −0.5 s | **−2.65 s (5.2 %)** | worth ~13 % of the t64 wall — the tail is 37 % of the call there | CR `intersection_areas.jl` | no, but each needs a correctness argument for the uniqueness assumption |
| 3′ | `sizehint!` the COO vectors in `addweight!` — **subsumed by #1**, but the cheap standalone version if #1 is not taken | −2.5 s | | CR `conservative.jl:370` | no |
| 4 | **`--gcthreads=8,1` at t64** (documentation/default, not code) | 0 at t8 (default is already optimal) | **−1.5 s (−7 %), −36 % GC, −49 CPU-s** | launcher docs | no |
| 5 | **GO `_child_extents` allocation pattern** — the single largest allocation site (~7 GB, 42 % of pair-gen bytes). Ceiling is bounded by pair-gen's whole GC deficit: **−1.3 s** | −1.3 s | proportional | **upstream GeometryOps** | **yes** — upstream API/ownership question |

**Combined ceiling: 50.8 s → ~36.3 s (−28.6 %) on a repeat call**, of which
levers 1+2 are 78 %.

**Measured non-levers.** Type instability (zero dispatch frames in the entire
call). The clip kernel itself (75 % of assembly is genuine Sutherland-Hodgman
FLOPs; the area integral is 1.7 %). Destination-geometry *allocation* (already
0 bytes). `assemble_sparse_matrix_coo`'s threading (util 7.13, the best in the
call). GC tuning at t8 (default optimal). Pair-gen's task balance at this CR pin
(the `508a637` estimator regression is absent here).

## Caveats

- Shared box, load 23–46 (51–62 at t64): wall times carry load noise;
  utilization ratios and CPU-second splits are contention-robust and reproduced
  across three warm reps to within ±0.15 threads.
- B's profiled call is 53.21 s against the 50.8 s baseline — **+4.7 % `Profile`
  overhead**, unevenly distributed (worst on short memory-bound serial phases,
  see the `sparse()` flag below). B's *shares* are used, never its absolute walls.
- `Profile.Allocs` byte coverage is **56–65 % of deterministic bytes** (object
  counts match to 0.2 %), so allocation byte shares are shares-of-sampled; every
  conclusion above rests on deterministic `GC_Diff` plus `@allocated` micros.
- The t64 figures are 1–2 reps each, on a box whose load rose to 51–62 during
  them; treat them as a scaling probe, not a benchmark.
- **Reconciliation flags** (all adjudicated, none affecting a conclusion):
  1. **CR `sparse()`: B 2.09 s vs A/C 1.62 s (+29 %)**. Adjudicated to
     A/C — the phase is a 1.6 s memory-bound `sparse!` where `Profile`'s
     per-sample cost lands hardest; absolute discrepancy 0.47 s. All other
     phases agree to ±7.3 %.
  2. **Region ordering.** B's flat sample rollup puts dual-tree search
     (42.4 %) above clip/assembly (41.2 %); A's and C's independent CPU clocks
     both put assembly above pair-gen (144.0 vs 127.1 CPU-s, 47.0 % vs
     41.5 %). Adjudicated to A/C: deterministic CPU-second windows beat
     sample-share bucketing, which carries several points of ambiguity from
     frames shared between the two regions (`raster_tree_memo`, `_memocell`).
     Both methods agree the two regions together are ~89 % of call CPU.
  3. **Plateau percentages.** A's headline 767 % / 745 % / 688 % are the
     *tightest inner segments* of each phase; the whole-phase means including
     ramp-up, GC dips and wind-down are 665 % / 713 % / 650 %. Both are in this
     note; the inner-segment figures are what `top` shows mid-plateau, the
     phase means are what the deterministic table reports.
  4. **"Serial regions are 71 % of the gap"** (A) counts the serial tail *plus*
     the `vcat` trough and uses raw serial wall; recomputed as excess over each
     region's own CPU/8 it is **69 %** (tail alone 61 %, trough 14 %, less the
     rounding overlap), consistent with C's 55 % for the tail alone in
     idle-thread-seconds. No conflict, only a definitional difference.
