# 2026-08-19 — The serial sparse-assembly tail: racing assembly strategies on the real triplets

The eager MOC regrid ends in ~9 s of single-threaded sparse assembly. The
2026-08-19 post-stack profile ranked *fusing* that tail (drop the COO → CSC →
COO → CSC round trip) as lever #1, worth −5.5 s. But an earlier campaign had
already tried exactly that and **measured it slower**, so this campaign captured
the real per-task COO arrays, raced assembly strategies against each other on
them, and then went back to find out why the two results disagreed.

**Verdict: both measurements were right, and the earlier one was reading a bug.**
The fused route is 5.1 s faster in isolation; the earlier implementation put the
denominator loop *inline* in `wholeblock`, where the assembled block's type is
not inferred, so the loop ran boxed and cost 13.9 s instead of 0.17 s. One
function barrier turns the earlier campaign's +7 s regression into **−5.7 s
(−11.2 %) end to end**, with 8.0 GiB less churn and a bit-identical result.

## Machine, trees, method

- **DGG**: main `7b3c727`, in a private scratchpad worktree; branch
  `claude/tail-assembly`. **ConservativeRegridding** pinned
  `claude/budget-frontier`, git-tree-sha1 `873cc732…` (v0.2.9, resolved
  `6a4b997`); **GeometryOps**/**GeometryOpsCore** `main` (v0.1.44 / v0.1.11).
  Julia 1.12.6, `-t 8` (and a `-t 1` arm), default `--gcthreads=8,1`.
- **Machine**: the shared 64-core box, 1-minute load recorded with every number
  (20–33 for the decisive runs, 38–58 for one diagnostic pass, noted per table).
- **Workload**: `bench/moc_probe.jl` conventions — a 3600×3600 Copernicus GLO-30
  tile N45/E010 as an in-memory raster, eager conservative regrid onto a
  level-13 IGEO7 MOC covering of 16,181,892 cells. 93,079,031 candidate pairs,
  **66,129,312 nnz**, digest `5b696475a3665634`.
- **Method**: a capture driver replays `CR.assemble_sparse_matrix_coo`'s
  partition loop through CR's documented developer entry points — no source
  edits — and serialises the 32 per-task `(rows, cols, vals)` triples
  (1.48 GiB) plus a reference CSC and denominator. The strategy race then runs
  entirely on that file: warm-up call per strategy, then reps interleaved
  strategy-by-strategy so drift is shared, `GC.gc(true)` twice between calls,
  `@timed` one-shots (every piece is ≥0.1 s, so BenchmarkTools' sampling buys
  nothing over repeated one-shots here) and `/proc/loadavg` on both sides of
  every measurement.
- **Raw data**: `2026-08-19-tail-assembly.ndjson` beside this note (167 rows:
  isolated race t8/t1, piece breakdown, both three-arm diagnostics, the
  end-to-end target passes, and the structural facts).

## 1. History audit — what was tried before, and what won

| when | where | what | verdict then |
|:--|:--|:--|:--|
| `regrid-notes/perf-P1.md` §P4 | 512²-scale | "remove the double sparse assembly" — `_fillcoo!` + second `sparse()` | **deferred**: 0.24 % of a one-shot regrid, filed as a memory win, not a time win |
| `regrid-notes/perf-P3.md` R6 | eager 512² | CR's `reduce(vcat)` COO merge | filed upstream, not taken |
| **`regrid-notes/2026-08-19-eager-moc-perf.md` item 5**, branch `claude/moc-wholeblock-adopt`, commit **`b9bb01d`** | **the target workload** | `wholeblock(::Conservative, …)` adopts the assembled CR matrix and reads denominators off it — the exact fusion | **DROPPED.** Same-session alternating A/B: fast 52.5 / 53.1 s vs generic 46.0 / 46.6 s, +0.47 GiB churn. "The fast path is reproducibly slower. Mechanism unresolved (post-stage should be strictly less work; suspect allocator/locality interaction inside the otherwise-identical assemble). Branch deleted." No PR was opened. |

So the user's recollection is exact: the fusion was tried on this workload and
the current route won. `b9bb01d` survives in the object store and its variant is
raced below as **C**, so the comparison is apples-to-apples.

## 2. The structural facts, measured

Capture run at t8, load 17.8–22.6 (`kind: "structure"` in the NDJSON):

| claim | evidence |
|:--|:--|
| 32 partitions, balanced | per-partition nnz 2,060,669 … 2,072,742; total **66,129,312** |
| every value is stored | `minimum(vals) = 3.523e-29`, **0 non-positive** — nothing is dropped by `_fillcoo!`'s `w > 0` guard |
| **no duplicate coordinates** | linear indices `(col-1)*nrows + row` sorted: **0 adjacent repeats**, 66,129,312 unique |
| `sparse()` merges nothing | `nnz(block) = 66,129,312 = nnz_coo` |
| the round trip is the **identity** | `_fillcoo!` emits 66,129,312 entries; re-`sparse()`ing them gives `colptr ==`, `rowval ==`, `nzval ===` against the first CSC, all `true` |
| and it is already sorted | the COO `_fillcoo!` emits is `issorted` by `cols` and by `(cols, rows)` — it hands `sparse()` back its own output |

The middle two steps of COO → CSC → COO → CSC are therefore provably redundant
**on the eager whole-block path**. They are *not* redundant on the chunked path,
where `indexmap` is non-trivial and blocks from several source chunks are summed
into one COO — any fix has to be path-aware.

## 3. The strategy race, on the captured triplets

Nine strategies, 3 interleaved reps each. **Correctness gate, passed by all
nine**: `colptr ==`, `rowval ==`, `nzval ===`, `denom ===` against the reference.

- **A** `current` — `reduce(vcat)` → `sparse` → `_fillcoo!` → `sparse` (baseline)
- **B** `sizehint` — A with the `WeightCOO` vectors pre-sized (lever 3′)
- **C** `adopt` — **`b9bb01d`'s variant**: keep the assembled CSC, one denominator pass
- **D** `adopt + parallel concat` — threaded `copyto!` per task segment
- **E** `current + parallel concat`
- **F** `adopt + no-dedup CSC` — serial counting sort by column, rows ordered in place (uniqueness verified in §2)
- **G** `adopt + parallel CSC` — coarse column histogram picks entry-balanced windows, one bucketing pass straight off the per-task arrays (**no concat at all**), then windows sorted independently
- **H** `G + parallel denominators` — each task scans the whole CSC in CSC order but accumulates only the destination rows it owns, so the summation order per row, and hence the rounding, is unchanged
- **I** `current + parallel CSC`

| strategy | t8 reps (s) | t8 GiB | t1 reps (s) | t1 GiB |
|:--|:--|--:|:--|--:|
| **A current** | 8.32 / 8.18 / 7.90 | 11.88 | 8.58 / 8.62 / 8.52 | 11.88 |
| B sizehint | 5.84 / 5.73 / 5.92 | 7.77 | 6.15 / 5.97 / 5.74 | 7.77 |
| **C adopt (`b9bb01d`)** | 2.63 / 2.45 / 2.50 | 3.88 | 2.39 / 2.44 / 2.38 | 3.88 |
| D adopt + par concat | 1.98 / 1.92 / 2.26 | 3.88 | 2.39 / 2.36 / 2.34 | 3.88 |
| E current + par concat | 7.86 / 7.96 / 7.77 | 11.88 | 8.60 / 8.48 / 8.38 | 11.88 |
| F adopt + no-dedup CSC | 1.84 / 1.85 / 1.84 | 2.87 | 2.05 / 1.78 / 1.76 | 2.87 |
| G adopt + par CSC | 0.51 / 0.46 / 0.52 | 2.87 | 1.98 / 1.98 / 1.98 | 2.87 |
| **H adopt + par CSC + par denom** | **0.50 / 0.47 / 0.49** | 2.87 | 2.00 / 2.00 / 2.00 | 2.87 |
| I current + par CSC | 6.03 / 6.23 / 6.03 | 10.87 | 7.73 / 7.80 / 7.76 | 10.87 |

t8 load 21.6–30.0, t1 load 19.7–31.8. Rep spread is ≤5 % everywhere except
D-rep3 (a load spike to 30). **The ordering is unambiguous at both thread
counts**: the gap between A and any adopt variant (≥5.3 s) is two orders of
magnitude larger than the noise.

Piece breakdown (3 reps each; t8 load 24–27, t1 load 24–27):

| piece | t8 (s) | t1 (s) | GiB |
|:--|--:|--:|--:|
| `reduce(vcat)` ×3 | 0.704–0.739 | 0.716–0.758 | 1.48 |
| threaded concat ×3 | **0.162–0.167** | 0.713–0.751 | 1.48 |
| `sparse()` (dedup) | 1.535–1.544 | 1.520–1.538 | 2.28 |
| uniqueness-assuming CSC, serial | 0.903–0.909 | 0.896–0.905 | 1.28 |
| threaded CSC off the per-task arrays | **0.296–0.333** | 1.868–2.008 | 2.75 |
| `_fillcoo!` | 3.485–3.648 | 3.702–4.038 | 5.71 |
| `_fillcoo!` pre-sized | 1.376–1.438 | 1.359–1.424 | 1.60 |
| denominator pass, serial | 0.166–0.168 | 0.163–0.206 | 0.12 |
| denominator pass, threaded | 0.115–0.116 | 0.157–0.158 | 0.12 |
| `WeightBlock` `sparse()` | 1.739–1.762 | 1.694–1.706 | 2.41 |

Reading: `_fillcoo!` is 3.5 s, **60 % of it `push!` growth** (pre-sizing alone
takes it to 1.4 s); the dedup pass in `sparse()` is worth 0.63 s on a matrix with
zero duplicates; the concat parallelises 4.4×, the COO→CSC 3.0×. The tail as
shipped is **7.6 s of the 7.9 s** strategy-A total, matching the post-stack
profile's in-situ split (1.62 + 3.83 + 1.94 s plus the `vcat` trough) to within
4 %, so the harness is faithful.

## 4. Why the earlier campaign measured the opposite

Strategy C is `b9bb01d` verbatim and it is 5.4 s *faster* than A in isolation.
Re-running the earlier same-session A/B on the current tree reproduced the
earlier verdict exactly (t8, load 29.6→49.4, ABBA order, both arms digesting
`5b696475a3665634` and agreeing bit for bit on `colptr`/`rowval`/`nzval`/`denom`):

| arm | round 1 | round 2 | churn |
|:--|--:|--:|--:|
| fast (`b9bb01d`) | 58.01 s | 60.00 s | 38.62 GiB |
| generic | 50.58 s | 54.04 s | 38.16 GiB |

The tell is the churn: the fast path does 8 GiB *less* work in the tail and yet
allocates **more** over the whole call. Adding a **control arm** — the parallel
region alone (operator, both subtrees, `_intersectionareas`, no tail at all) —
localises it (t8, load 38–58, ABBA):

| arm | wall | CPU-s | util | churn | tail vs core |
|:--|--:|--:|--:|--:|--:|
| core (no tail) | 46.26 / 43.81 s | 298.9 / 296.3 | 6.46 / 6.76 | 30.04 GiB | — |
| generic | 53.75 / 51.36 s | 306.6 / 306.1 | 5.70 / 5.96 | 38.16 GiB | +7.5 s, +8.12 GiB |
| fast (`b9bb01d`) | 59.87 / 57.71 s | 312.3 / 311.2 | 5.22 / 5.39 | 38.62 GiB | **+13.7 s, +8.58 GiB** |

A denominator pass that costs 0.17 s and 0.12 GiB standalone was costing 13.7 s
and 8.58 GiB in situ. The cause is inference:

```
_intersectionareas      -> SparseArrays.SparseMatrixCSC      # no type parameters
CR.intersection_areas   -> SparseArrays.SparseMatrixCSC
```

(`Base.infer_return_type`, real space types.) CR's assembly deliberately takes
one dynamic dispatch per chunk (`task_local_operator`'s return type is not
inferrable), so `assemble_sparse_matrix_coo` returns a `vals` vector of unknown
eltype and `sparse()` hands back an **unparameterised** `SparseMatrixCSC`. With
the denominator loop written *inline* in `wholeblock`, `rowvals(block)` and
`nonzeros(block)` are `Any`, every `vals[t]` boxes a `Float64` and every
`denom[rows[t]] += w` is a dynamic dispatch — 66.1 M times.

The generic route never had this problem: `_fillcoo!(coo, block)` is a
**function barrier**, so its loop specialises on the runtime type. The earlier
campaign's "suspect allocator/locality interaction" was a missing barrier.

**Fix**: move the loop into `_blockdenom(block::AbstractSparseMatrixCSC, ndst)`.
Same arithmetic, same CSC order, one dispatch instead of 66.1 million.

| arm (post-fix, t8, load 30–43, ABBA) | round 1 | round 2 | churn | util |
|:--|--:|--:|--:|--:|
| core (no tail) | 44.59 s | 44.27 s | 30.04 GiB | 6.63 / 6.65 |
| **fast (adopt + barrier)** | **43.78 s** | **44.50 s** | **30.16 GiB** | 6.77 / 6.61 |
| generic | 50.17 s | 49.80 s | 38.16 GiB | 6.06 / 6.11 |

The adopt tail now costs **+0.12 GiB and nothing measurable in wall** over the
core, exactly as the isolated race predicted.

## 5. End to end, against main, one window

`bench/moc/target.jl`, 3 warm reps per session, sessions interleaved
main → branch → branch → main, t8, load 27–33 throughout.

| pass | arm | rep 1 (cold) | rep 2 | rep 3 | churn (warm) | GC | maxrss |
|--:|:--|--:|--:|--:|--:|--:|--:|
| 1 | main `7b3c727` | 61.73 s | 51.16 s | 50.35 s | 38.58 GiB | 5.1–5.4 s | 10.59 GiB |
| 2 | branch | 53.51 s | 44.89 s | 45.30 s | 30.58 GiB | 4.1–4.4 s | 9.91 GiB |
| 3 | branch | 53.96 s | 46.10 s | 44.80 s | 30.58 GiB | 4.1–4.4 s | 9.68 GiB |
| 4 | main `7b3c727` | 59.80 s | 51.14 s | 51.24 s | 38.58 GiB | 5.2–5.3 s | 10.73 GiB |

**Warm t8: main 50.35–51.24 s → branch 44.80–46.10 s, −5.7 s (−11.2 %).** The
two ranges do not overlap across four independent sessions. Churn −8.00 GiB
(−20.7 %), GC −1.0 s, maxrss −0.9 GiB. Cold first call −6.9 s. Digest
`sum = 3.368822268175383e10  sha = 5b696475a3665634  n = 16181892` on **every**
pass, both arms.

## 6. What ships, and what does not

Shipped as **PR #56** (`claude/tail-assembly`, 45 lines in `lib/GlobalRegridding/src/conservative.jl`):
`wholeblock(::Conservative, dst_space, src_space)` adopts the assembled matrix
and takes the denominators off it through `_blockdenom`. Degenerate sides
`invoke` the generic method, which is the only place the two disagree
(`ndst == 0` leaves the generic block without denominators). The chunked path,
which needs the round trip, is untouched.

**Not shipped, and why.** The race's winner is **H** at 0.47 s — another 2.0 s
below the adopt path. But everything below C lives in ConservativeRegridding,
not here: the concat and the COO→CSC are inside
`assemble_sparse_matrix_coo`/`intersection_areas`, and DGG only ever sees the
finished `SparseMatrixCSC`. The pins stay on `claude/budget-frontier` for this
campaign, so those are filed upstream, with the isolated numbers to justify them:

| upstream lever | isolated t8 win | note |
|:--|--:|:--|
| threaded concat replacing `reduce(vcat)` (CR `intersection_areas.jl:204-206`) | −0.55 s | trivially safe, no assumption |
| uniqueness-assuming COO→CSC (skip the dedup pass) | −0.63 s | needs a per-operator uniqueness contract; **true** for `BlockAreaOperator`, not for operators that emit several entries per work item |
| threaded COO→CSC straight off the per-task triplets (subsumes the concat) | −1.2 s | column-window bucketing; bit-identical, verified on 8 randomised shapes and on the real matrix |
| threaded denominator pass | −0.05 s | row-range partition keeps the summation order |

A fifth non-lever worth recording: `sizehint!`ing the `WeightCOO` (lever 3′,
strategy B) is worth 2.4 s on the *chunked* path, which still runs `_fillcoo!`
and is not helped by the fusion at all. It is the cheap standalone follow-up.

## Caveats

- Shared box, load 20–33 for the decisive runs and 38–58 for the pre-fix
  three-arm pass; that pass is used only for its *ratios* (fast vs generic vs
  core in the same session), never for absolute wall.
- The isolated race replays the tail on serialised triplets, so it excludes the
  heap state a real call arrives with. The end-to-end numbers in §5 are the ones
  that decide; the isolated numbers explain them.
- The uniqueness result in §2 is measured on this workload, not proved for every
  operator; the shipped change does not depend on it (it keeps CR's `sparse()`).
  The upstream levers that do depend on it are flagged above.
- `_intersectionareas`'s abstract return type is the real defect here; the
  barrier is a fix at the call site, not at the source. Making
  `assemble_sparse_matrix_coo` infer its eltype upstream would remove the trap
  for every future caller.
