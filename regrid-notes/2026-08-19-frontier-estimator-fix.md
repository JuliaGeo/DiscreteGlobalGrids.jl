# 2026-08-19 — Frontier estimator fix + nesting-parallelism audit

Follow-up to `2026-08-19-thread-utilization.md` §3's finding that the budget
frontier's `pair_weight` was flat on eager RasterGrid→DGGSpace builds. Shared
64-core box, loadavg 20–36 throughout (per-row in the NDJSON where recorded);
all comparisons back-to-back in the same session. DGG base: origin/main
9fc8514 (per-arm worktrees); CR arms: main = 89ec5ca (v0.2.9, with DGG
a8fec729 since main CR predates `split_weight`), frontier-before = 508a637,
frontier-fixed = 6a4b997.

## Job 1 — nesting audit: verdict (a), one site

DGG nests parallelism in exactly one place: the lazy wave.
`_fillwave!` (lib/GlobalRegridding/src/lazy.jl:476 at 9fc8514) spawns one
task per chunk-pair block; each block's conservative build reaches
`_intersectionareas` (conservative.jl:294), which passed threaded=True
whenever `Threads.nthreads() > 1`. Every other CR invocation — eager
`wholeblock` (api.jl:215), interpolation methods, discovery — is top-level
or serial; DGG-core threaded regions (neighborhood.jl:374, adjacency.jl:165)
never reach CR.

Policy implemented per the maintainer's directive ("outer parallelism
wins"): a `ScopedValue` (`GR.OUTER_PARALLEL`) declares the wave and
`_innerthreaded()` resolves the inner CR flag to False inside it; spawned
tasks inherit the scope. PR #46, merged (0bc2b02), lib suite 486 pass/1
broken = 475 baseline + 7 (PR #45 tests) + 4 (its own). Measured on the
affected workloads (t8, fixed CR, two rounds, identical digests):
l512 2.78–2.86 s (nested) → 2.86–3.02 s (outer); l1024 11.94–12.06 →
12.14–12.79 s; round 2 within noise. Cost 0–6% wall at t8; the win is that a
wave no longer multiplies into nthreads × (frontier + clip chunks) tasks.

## Job 2 — the estimator

**Root cause.** `Trees.ncells(::RasterCellTree)` (and `RasterFlatTree`,
`CellCapTree`) answers the whole space — the matrix-size contract — so CR's
`split_weight` fallback read N_total at every node. In `pair_weight` that
density's 1/cap_area(node) scaling cancels against the overlap area exactly:
any raster window facing a containing DGG subtree weighs ≈ N_total (measured
estw 2.6–2.8e5 ≈ 512², CV 0.02, at every depth 1–9). The DGG side's honest
count is drowned (its system caps dilute density ~6000×). Splits never
reduced the heap max, every split inflated the estimated total, and the
share test stalled: 117 pairs, one depth-1 pair holding the full 334k-cell
covering and 50% of all 1.89M candidate pairs.

**Fix, CR side** (6a4b997 on `claude/budget-frontier`, PR #137 — commented,
not merged): a split may only redistribute its parent's weight — children
whose raw `pair_weight`s sum past the parent's are scaled back, never up.
Honest subadditive estimates pass through unchanged (sibling ratios
preserved, d4ef556's tied-pair behavior intact); a saturating estimate
degrades to even breadth-first splitting. Regression test: a lying-weight
wrapper + deterministic pair-count bound (≤10% max share; 23–25% before,
5–6% after). CR suite 10053 pass / 0 fail at t8 = 10049 baseline + 4 new.

**Fix, DGG side** (PR #45, merged 7f801d4): exact `Trees.split_weight` for
`RasterCellTree` (index rectangle), `RasterFlatTree` (stored entries),
`CellCapTree` (cap window) — same shape as `PositionTreeNode`/`BlockCursor`.

**Balance and estimates** (deterministic, nchunks=64, from
`frontier`+per-pair serial DFS):

| case | arm | ntasks | max task share | corr(est, npairs) |
|:--|:--|--:|--:|--:|
| e512 | before | 117 | 50.0 % | 0.20 |
| e512 | CR fix | 190 | 1.6 % | 0.20 |
| e512 | CR+DGG fix | 190 | 1.6 % | 0.42 |
| e1024 | before | 141 | 50.0 % | 0.21 |
| e1024 | CR fix | 154 | 1.6 % | 0.28 |
| e1024 | CR+DGG fix | 154 | 1.6 % | 0.36 |

Residual corr ceiling: ~half the frontier pairs are geometric false
positives (caps intersect, zero pairs) which no O(1) model sees; they cost
~0 CPU.

**Wall/util** (eager plan, t8, warm, same session, load 29–36):

| size | frontier@508a637 | fixed | CR main |
|:--|--:|--:|--:|
| 512² | 6.54 s · 2.78 | **2.63 s · 6.78** | 2.67 s · 6.98 |
| 1024² | 22.83 s · 2.84 | **9.24 s · 6.78** | 9.74 s · 6.78 |

DGG-native source unchanged (d512 19.6→20.0 s, d1024 75.8→75.6 s, util
6.9–7.2). Identity: pair lists == serial at chunks_per_thread ∈ {1,2,8,32}
and weight matrices bit-identical (colptr/rowval/nzval) on e512 and d512 at
t8; regrid digests 424312876.1135559 (512²) / 1311105293.6084747 (1024²)
identical across all arms and both jobs' trees.

Raw data: session scratchpad `work/` (eager2-t8.ndjson, dgg-t8.ndjson,
lazy-t8.ndjson, diag-*.ndjson, cr-suite-fix.log, dgg*-tests logs).
