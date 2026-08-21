# 2026-08-18 — Threading policy for CR's dual-tree traversal: the wave-3 verdict

Decides the granularity policy for `ConservativeRegridding`'s multithreaded
dual depth-first search, to ship as its own CR PR (the "later CR PR" of
`2026-08-18-integration-plan.md` §5). Base tree for every arm:
`260f9199` (slug `08cPv`, v0.2.8) — the same tree the post-fix DGG campaign
(`bench/results/postfix-REPORT.md`) measured. Machine: the shared 64-core box,
load average ~25–45 throughout; every NDJSON row records `loadavg` at emit.
Raw record: `regrid-notes/2026-08-18-threading-policy.ndjson` (540 rows, one
per repetition) and `regrid-notes/2026-08-18-threading-policy-logs.tar.gz`
(per-run logs + the arms' unified diffs against the base tree); the working
copies lived in the session scratchpad (`wave3/`).

## The arms

- **A** — the shipped policy: recurse and spawn a task whenever both nodes'
  `should_parallelize` predicates fire; foreign node types hit the `::Any`
  cap-area fallback, which fires at the top of the tree. 18,944–59,888 tasks,
  up to 49 % of them empty (dgg: 9,219 of 18,944).
- **E** — same recursive spawn, leaf-count default policy in `src` only;
  foreign trees under-split (22 tasks on `clima_cs>ocea_latlon`, one task
  holding 12 % of all pairs at any thread count).
- **CHEAP** — E plus hand-written per-extension `should_parallelize` methods
  (leaf-count rule against the *global* cell total) for the cubed-sphere,
  HEALPix and ring-grid node types.
- **S1 / S2** — **budget frontier**: phase 1 walks both trees serially,
  keeping a max-heap of node *pairs* keyed by an overlap-area × cell-density
  weight, splitting the heaviest pair until `nchunks =
  CR_CHUNKS_PER_THREAD × nthreads` pairs exist; phase 2 spawns one task per
  pair; results concatenate in DFS pre-order, reproducing the serial emission
  order exactly. S1 = 32 chunks/thread, S2 = 8. `subtree_cells(node)` (O(1),
  answered by `ncells`) orders the heap; a wrong answer costs balance, never
  correctness.

## Correctness: identical output, identical order, everywhere

- Every arm × workload × thread count in the timing pass returns the same
  `npairs` and the same `pairhash` (e.g. `0xe70b92ad7559b63c` for the DGG
  case, all five arms, t1/t8/t64).
- The explicit gate (`CHECK=1`: serial re-traversal, order and set compared
  element-wise) ran at t1, t8 and t64 for all five arms × 6 workloads:
  **55 checks, `order_identical=true` and `set_identical=true` in all 55,
  zero violations** (`res/gate_*.log`, `res/gate.out`).
- End-to-end weight matrices (`E2E=1`) agree per workload across arms on
  nnz, sum and sighash (e.g. dgg: nnz 106,309, sighash `0x439cccea9d5dc993`
  for A, E, CHEAP, S1, S2 alike; `res/e2e_*.log`).
- The real DGG regrid (pass B, `bench/harness.jl` workload) returns the
  identical digest per size across all arms and thread counts
  (`n1326354_s1311105293.6084747` at 1024², 30 runs).
- `subtree_cells` fallback probe: **0 fallback hits** over the five foreign
  sweat workloads — every node type in the suite answers `ncells` O(1)
  (`res/probe_fallback.log`).

## Traversal wall time (min over 2 rounds × 3 reps, seconds)

| case | nt | A | E | CHEAP | S1 | S2 |
|:--|--:|--:|--:|--:|--:|--:|
| dgg_n144 | 1 | 4.265 | 4.270 | 4.216 | 4.236 | 4.266 |
| dgg_n144 | 8 | 0.625 | 0.620 | 0.625 | 0.628 | 0.645 |
| dgg_n144 | 64 | 0.580 | 0.254 | 0.249 | 0.273 | 0.253 |
| clima_cs>ocea_latlon | 1 | 1.338 | 1.330 | 1.332 | 1.326 | 1.328 |
| clima_cs>ocea_latlon | 8 | 0.206 | 0.233 | 0.200 | 0.207 | 0.224 |
| clima_cs>ocea_latlon | 64 | **0.845** | 0.194 | **0.820** | 0.124 | **0.083** |
| ocea_latlon>hpx_nested | 1 | 1.891 | 1.889 | 1.892 | 1.889 | 1.892 |
| ocea_latlon>hpx_nested | 8 | 0.338 | 0.395 | 0.337 | 0.319 | 0.382 |
| ocea_latlon>hpx_nested | 64 | **1.178** | 0.268 | 0.430 | 0.216 | **0.186** |
| ocea_latlon>ocea_tripolar | 1 | 2.094 | 2.102 | 2.097 | 2.086 | 2.078 |
| ocea_latlon>ocea_tripolar | 8 | 0.337 | 0.337 | 0.328 | 0.327 | 0.344 |
| ocea_latlon>ocea_tripolar | 64 | 0.772 | 0.760 | 0.725 | 0.194 | **0.167** |
| sw_clenshaw>sw_gaussian | 1 | 0.894 | 0.894 | 0.896 | 0.888 | 0.892 |
| sw_clenshaw>sw_gaussian | 8 | 0.147 | 0.147 | 0.145 | 0.147 | 0.166 |
| sw_clenshaw>sw_gaussian | 64 | 0.665 | 0.809 | 0.816 | 0.107 | **0.078** |
| hpx_nested>hpx_ring | 1 | 2.353 | 2.355 | 2.364 | 2.341 | 2.335 |
| hpx_nested>hpx_ring | 8 | 0.545 | 0.536 | 0.484 | 0.478 | 0.512 |
| hpx_nested>hpx_ring | 64 | 0.385 | 0.392 | 0.723 | 0.341 | **0.314** |

Reading it:

- **t1: all five arms identical to <1 %.** No policy costs anything serial.
- **t8: all within ~10 % of each other**; A is competitive here (its
  oversubscription has not started to hurt), S2 trails A by 5–13 % on three
  cases, S1 stays within 5 %.
- **t64 is the verdict.** A loses monotonicity on four of six workloads
  (clima t64 is 4.1× its own t8 time; hpx_nested 3.5×), burning 12–28 CPU-s
  where S2 burns 1.1–5.3. E and CHEAP each blow up on the workloads their
  hand policies don't fit (E: clenshaw 0.81 s; CHEAP: clima 0.82 s,
  hpx_ring 0.72 s). **S2 is fastest or tied-fastest in all six cases at t64**
  (up to 10.1× over A on clima, 8.6× on clenshaw), S1 second in five.
  Frontier phase-1 cost is 4–22 ms — noise.
- The recursive arms' failure is structural, not tuning: A's task count is
  fixed by tree shape (18,944–59,888 tasks with up to half empty, max-task
  share 8–21 % of all pairs on the E side), so no `should_parallelize`
  spelling fixes all six workloads at once — E and CHEAP are exactly such
  spellings and each fails somewhere. The frontier sizes tasks by measured
  weight, so its worst task share falls from 12 % (E, clima) to 0.1–0.8 %.

## End-to-end: the real DGG regrid (pass B, min over 3 reps, seconds)

512² and 1024² GLO-30 windows through `bench/harness.jl`'s lazy conservative
path (srcchunk (128,128), chunklevel 7); identical digests everywhere.

| size | nt | A | E | CHEAP | S1 | S2 |
|--:|--:|--:|--:|--:|--:|--:|
| 512² | 1 | 29.25 | 29.25 | 28.98 | 29.16 | 29.20 |
| 512² | 8 | 19.01 | 20.17 | 19.89 | 17.30 | **17.31** |
| 512² | 64 | 18.42 | 19.50 | 19.06 | 16.73 | **16.31** |
| 1024² | 1 | 76.58 | 76.41 | 76.16 | 76.75 | 76.58 |
| 1024² | 8 | 30.20 | 35.83 | 35.07 | 27.92 | **27.45** |
| 1024² | 64 | 25.80 | 30.46 | 29.82 | 24.21 | **22.20** |

The traversal is a minority of regrid wall, so the 10× microbenchmark gaps
compress — but the ordering survives: **S2 is fastest at every thread count
measured, 9 % over A at t8/1024² and 14 % at t64**, and never loses at t1.
(The A t8/t64 rows and part of the E t8 run were first measured concurrently
with another wave3 pass by mistake; those logs are quarantined as
`res/B_*.log.contended` / `designspace/wave3-regrid-contended.ndjson` and every
number above is from the serialized rerun.)

## The `balance` refinement

`CR_BALANCE=0` (stop at exactly `nchunks` frontier pairs, no post-split of
outliers) at t8/t64: within noise of S1 everywhere, occasionally faster
(clima t64 0.103 s vs S1's 0.124 s), with worse worst-task share (0.6 % vs
0.1 %). The refinement costs milliseconds and buys balance headroom; keep it
on by default, keep the switch.

## Verdict

1. **Adopt the budget-frontier policy** (arm S) as CR's threading strategy:
   serial heap-split to `chunks_per_thread × nthreads` pairs, one task per
   pair, DFS pre-order concatenation. It is the only arm that is never
   pathological on any of the six workloads at any thread count, and it
   preserves the serial emission order bit-for-bit (55/55 checks).
2. **Default `chunks_per_thread = 8`** (S2). End-to-end it beats 32 (S1) at
   both t8 and t64; in the traversal microbenchmarks it wins t64 outright and
   its worst t8 deficit (13 %, clenshaw) is on a 0.15 s traversal. Keep the
   env override.
3. **Keep `balance = true`**; the refinement is milliseconds.
4. **Ship `subtree_cells` with the generic `ncells` answer** plus the four
   one-line methods for wrapper/multi-tree roots; the probe measured zero
   fallback hits across every suite grid, and a fallback costs balance only.
5. **Drop the per-extension `should_parallelize` tuning route** (arms E and
   CHEAP): each spelling measured here fixes some workloads and regresses
   others, and the frontier makes the predicate irrelevant to task granularity.

