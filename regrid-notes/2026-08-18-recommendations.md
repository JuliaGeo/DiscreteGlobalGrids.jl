# 2026-08-18 — Consolidated recommendations

The single read for "what next and why", now that all measurement campaigns are
done. Evidence lives in `regrid-notes/perf-P3-brief.md` (15 ranked defects),
`regrid-notes/alloc-attribution.md` (churn decomposition),
`regrid-notes/2026-08-18-threading-policy.md` (budget-frontier verdict),
`regrid-notes/2026-08-18-integration-plan.md` (PR mechanics),
`bench/results/postfix-REPORT.md` (mode/size/thread campaign),
`bench/results/mem-postfix-REPORT.md` (memory envelope), and
`bench/results/cache-summary.txt` (CR #131 isolation). This document cites
their numbers; it does not re-derive them.

## 1. Verdict

As of 2026-08-18 the stack regrids one Copernicus GLO-30 tile onto IGEO7 L13
correctly, conservatively, and on small hardware: Conservative output is
bit-identical across 80/82 mode configurations (the two exceptions are
`DGGSpace`-source runs off by half a Float32 ulp) and across 1–64 threads,
conservation closes to 8.3e-10 source-side at the full tile, and the lazy path
fits in 2 GB (floor 1.42 GB sampled / 1.60 GB kernel) at ~160 s on 8 threads.
What blocks calling it ready for global Copernicus scale is not throughput but
three specific items: `BilinearPoint` returns missing or wrong values whenever
the destination is tiled (root cause found in `_onbranch`, verified 2-line fix,
PR being authored now); the library's default destination chunking produces one
destination chunk at every size and does not finish the full tile within 2400 s
that `chunklevel = 6` completes in 156 s; and the open release gates —
GeometryOps registration, a ConservativeRegridding release carrying #131, then
the DGG unpin. Beyond those, truly global work has two absolute walls no
per-tile fix moves — the ~51 TB global weight matrix (no single-matrix
formulation exists) and the measured 27.7 GB / 813 s complete-globe L13
`DGGSpace` — but neither binds tile-at-a-time processing, which is the workable
formulation today.

## 2. Correctness, ranked

**C1 — `BilinearPoint` regional-raster branch fold (top item; fix PR in
flight).** `_onbranch` (`lib/GlobalRegridding/src/rastergrid.jl:961-965`) folds
longitudes to the low end of the chart branch, throwing points just west of a
regional raster a full period east; the non-periodic `_locate` then clamps to
the last column. Untiled runs fill the west pad ring with far-side garbage;
tiled runs lose the cells outright — 812–7,751 across the measured sizes, and
one tiling returns finite values wrong by up to 3.848 in the data's own units
(`postfix-REPORT.md` §5.6). Until it lands there is no configuration where
bilinear is both correct and tractable at 3600² (untiled DNF-adjacent at
1671 s vs 28 s tiled). The 2-line nearest-representative fix drives every
tiling delta to zero and changes untiled output only for out-of-coverage cells.
Orthogonal design question to settle in review, not in this PR:
clamp-extrapolation vs NaN outside coverage.

**C2 — the default destination chunking is a defaults bug.**
`chunkcells = DEFAULT_CHUNK_CELLS = 4096` (`src/regridding.jl:9`) yields **one**
destination chunk at every size 256²–3600²: 19× at 2048² with (128,128) source
chunks (1104 s vs 58 s, bit-identical output), DNF at the full tile vs 156 s at
`chunklevel = 6`. Root cause is P3-brief D3: `_chunklevel` compares the
subset's cell count against global level counts (a 4000× miss). Fix D3 and D4's
`_chunkwindows` global-level scan **together or not at all** — a corrected
chunk level that `_chunkwindows` cannot enumerate is not buildable.

**C3 — `chunks = (n,)` defeats empty-pair pruning.** Pruning goes 94% → exactly
0% (`builds = ntiles × nsrcchunks` to the unit); 3.2× at 1024², 14.6× at 2048²,
non-terminating at 3600² at all three values tried. Same weights, same bytes
read. Either make flat runs participate in pruning or warn/refuse; today an
ordinary keyword silently turns off the mechanism that makes lazy viable.

**C4 — `Spilled` on tmpfs is an OOM trap.** Spill files on tmpfs are charged to
the cgroup: 2 GB cap OOM-kills on `/tmp`, runs at 169.6 s / 2.03 GB peak on
btrfs (P3-brief D15). The obvious idiom `Spilled(mktempdir())` lands on tmpfs
on most Linux systems. Cheap fix: check the filesystem at construction, warn or
refuse on tmpfs/ramfs.

**C5 — `NearestCell` with a DGG source effectively hangs** (P3-brief D6):
196.5 s of `plan_regrid` for a 16×16 window, 0.42 s per destination cell.
Throw or fix; do not leave the hang.

Watch items, not defects: S2's extent margin is one ULP by design (P3-brief
D1b) — any change to corner/centre arithmetic there needs the conformance
covering law run first; conformance reach should grow per D1b's four
recommendations (depth ≥ 6 / branch ≥ 3, pentagon one-rings, a
`tree_covering_problems` for the cursor layer, and **no** cap-nesting
assertion); the CR depwarn traces to CR's own test file (4 manifold-less
`GO.area` calls, pre-existing) — a test-only upstream PR.

## 3. Performance, ranked by measured value

All P3 figures are the 1024² lazy read (61.17 s baseline, 1 thread) unless
noted; churn shares from `alloc-attribution.md` (eager 1024² plan, 18.44 GB,
full-tile scale ×11.99).

| rank | item | measured value | where |
|---|---|---|---|
| P1 | Budget-frontier traversal, `chunk_factor = 8`, balance refinement off | real regrid run t8 14.9 → 11.6 s (1.29×), t64 9.9 → 7.2 s; 1.4–7.6× on isolated t64 traversals; digests identical | threading-policy doc; own CR PR |
| P2 | DGG #37 inline cell polygon (open) | erases the 5.1 GiB clip-loop boundary chain ≈ 27% of plan churn (~61 GiB at full tile) | alloc-attribution; merges clean |
| P3 | Subtree/tile boundary cache (D2/R1) | −17.5 s standalone ≈ −28% of the read; with P5, −18.1 s → 43.1 s (1.42×) | P3-brief D2 |
| P4 | `cell_cap` overridable hook + `PartialGrid`/`AuthalicGrid` forwards (D2a) | prerequisite: 50.2% of 14.2 M `cell_boundary` calls enter there and bypass any `node_extent` override | P3-brief D2a |
| P5 | Analytic IGEO7 extent, `K_cell = 1.18`, `K_node = 1.25` (D2b) | 5.44× per call, 0 B; marginal −0.6 s once P3 exists, but covers the 10/22 tiles where the cache disengages and unlocks `node_extent_is_expensive = false`. `K = 1` under-bounds 13.6% — never ship it | P3-brief D2b |
| P6 | `spherical_distance` once per cap (D12) | ~−12% (12.6 s) for ~10 lines in `caps.jl:67-70` + `rastergrid.jl:610-613` | P3-brief D12 |
| P7 | Materialise source-leaf caps (D13) | 14.42% = 8.8 s → ~0.79 s, ≈ −13%; 16 caps × 32 B per leaf | P3-brief D13 |
| P8 | Materialise `PartialGrid.ids` (D11) | 8.34% = 5.1 s of 96× redundant decode; 10.6 MB at 1024². Counting overlap with P3: R1+R4 ≈ −20 s combined, not additive | P3-brief D11 |
| P9 | D3+D4 `_chunklevel`/`_chunkwindows` fix | same change as C2; also removes the 29–30 GB global-`DGGSpace` chunk blowup path | P3-brief D3/D4 |
| P10 | Descent-side cap memoization (`src/fallbacks/cursor.jl:269`) | 4.0 GiB = 22% of plan churn (leaf caps rebuilt once per opposing leaf, 3232 B/cell); **not** covered by #37 | alloc-attribution |
| P11 | Upstream area-only intersection operator + streaming `area.jl:275` | area 14.9% (~31 GiB full tile, one `collect`); clip result polygons 13.0% (~27 GiB); shared degenerate polygon for 2.17 M disjoint pairs, 2.5% | alloc-attribution; GO-owned |
| P12 | Lazy build-wave scaling | summed `buildtime` flat at ≥625 s full tile at any thread count; lazy stalls at 8.3× vs eager's 9.9× at 64 threads | postfix-REPORT gap #3 |

Two campaign facts that shape defaults rather than patches: **eager vs lazy is
a memory choice, not a time choice** (end-to-end 161.9 s eager vs 158.7 s lazy
`chunklevel = 6`; 5.4 GiB vs 248 MiB RSS growth — recommend lazy, (128,128)
source chunks, `chunklevel = 6`, eager only when the operator is reused), and
**`budget = 2^26` is the small-machine sweet spot** (1.82 GB / 163.9 s, within
2% of unconstrained time; 2 GB cap holds, 1 GB does not; eager OOMs at 6 GB
regardless of thread count).

## 4. Sequencing

Merge order per `2026-08-18-integration-plan.md`, unchanged:
**#34 → (fix conformance-CI `[sources]` GeometryOpsCore subdir on b-series) →
#35 → #37 → (merge b-series into b4-adopt) → #36 → (fix `toys.jl` `max_level`;
sign off `halo` retype + `subtree` name clash) → #38 → perf work written
against `src/engine/` and the new names.** One integrated suite run on the
assembled tree before #38 merges. The bilinear fix PR (C1) slots in
independently — `rastergrid.jl` is outside every stack conflict site — land it
as soon as review clears; likewise the D3+D4 defaults fix (C2/P9).

Perf implementation order (dependencies, then value): D2 cache → `cell_cap`
hook → analytic IGEO7 extent → D12 → D13 → D11 → D3+D4 → D5 only with its A/B.
Re-measure the P3 baseline once on the post-stack tree before writing patches.

External gates: GeometryOps release is cut and awaiting registration; once it
registers, cut the ConservativeRegridding release **with #131 merged**
(isolated: 7–10% wall, −33% churn, ~103 GiB at full tile, digests identical —
`cache-summary.txt`) — that also moots the pin branch's missing #133 vcat fix.
**Do not merge CR #134** (see §5). The budget-frontier traversal (P1) ships as
its own subsequent CR PR: `chunk_factor = 8`, balance refinement off,
`should_parallelize` deprecated, optional `split_weight`. Then unpin CR + GO in
the three pin files, reconcile the GO 0.1.43 compat line, and cut the
registrable DGG release. CompatHelper PRs get batch-regenerated at unpin time,
not before.

## 5. Explicitly not worth doing

- **Merging CR #134.** Measured 27–51% *slower* than shipping `&&` on the real
  DGG regrid (512² t8: 5.85 s vs 4.54 s; 1024² t64: 14.91 s vs 9.85 s); its
  22–48-task ceiling costs 8–12% max-share imbalance on foreign trees. The
  frontier subsumes its one virtue.
- **The CHEAP alternative** (#134 + hand-written per-tree policies): also a
  regression on the regrid path (6.19 s vs 4.54 s at 512²/t8), and it converts
  "extensions need nothing" into "every extension tunes a policy".
- **Frontier balance refinement / `chunk_factor = 32`.** `CR_BALANCE=0` is
  equal or faster on 6/6 cells; 8 beats 32 at t64 on 6/6 traversals and both
  regrid windows.
- **Further clip-scratch hunting.** The sampler finds 0 B of cache-buffer
  allocation in the shipped arm — the remaining 13% of clip churn is the result
  polygon itself; only an upstream area-only operator (P11) can take it.
- **`cells_cap` tightening** — 0.14% of wall clock.
- **D7 IGEO7 `cell_area` closed form for this pipeline** — zero profile samples
  (`Conservative` never calls it); worth doing only for `Extensive`/`cellarea`,
  and it is not a drop-in (−11.76% at pentagons).
- **D5 bucket-size change without the A/B** at `bucket_size ∈ {0, 7, 16, 49}` —
  after D2, extents are near-free and bucketing may be a net loss.
- **A cap-nesting conformance assertion** — fails 100% of
  `BlockCursor`/`RasterCellTree` nodes, correctly, for a benign reason
  (P3-brief D1b).
- **A `has_analytic_extent` trait** — duplicates `node_extent_is_expensive`,
  which already has a consumer and folds at compile time.
- **Rescuing eager with fewer threads or heap hints** — eager churns the same
  ~206 GB at any thread count, OOMs under 6 GB at 1, 2, and 8 threads, and
  `--heap-size-hint` moves GC fraction, not peak.
- **One-sided descent as a standalone traversal change** — the 55 dropped pairs
  are all zero-area, but it over-tasks (2,801 vs a 256 target) and the
  frontier's one-sided *splitting* already covers its degenerate case.
- **Chasing eager-vs-lazy on time** — 2% apart end-to-end at the full tile,
  inside the ±20% noise band; only the 22× memory gap is real.
