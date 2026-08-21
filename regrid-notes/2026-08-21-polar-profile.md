# Polar regime profile: CopDEM GLO-90 → IGEO7 L12

**The premise this campaign was commissioned on is false. There is no polar
penalty.** A deep-Antarctic column costs the same as a mid-latitude one, to
within 5 %, and the phase split of its CPU is the same to within 2 points. The
5× seen in the production done log is an artefact of *when* polar columns were
scheduled, not of *where* they are.

Worktree `/home/asinghvi17/geo/DGG-polarprof`, detached at `e2f90d1` — the exact
code the production run's final epoch executed — with the run's known-good
`Manifest.toml` copied from `/home/asinghvi17/geo/DGG-perf-ladder`
(ConservativeRegridding 0.2.9 @ `6a4b997`, as in the run). Measurements in
`2026-08-21-polar-profile.ndjson`. Harness in the session scratchpad, built on
the perf-ladder campaign's `ladder-scratch/setup.jl`; see §10.

**Box conditions.** The production run exited at 13:53Z; nothing held its
done log. Other users and other agents held 15–25 cores of 64 throughout. All
single-column numbers are **core-seconds** (`utime+stime` of the process across
the call), which is insensitive to that background load; spread across 3 reps
was under 5 %. Between-process spread is larger, ~15 %, and is the noise floor
for §9.

---

## 1. Headline: intrinsic cost is flat in latitude

On a quiet box at `-t 1`, every full-land column between 26 °N and 83 °S costs
**9.6–10.6 core-seconds**, flat in latitude and flat in source-tile count from 4
to 24:

| column | tiles | lat | production `secs` | wall s | **core-s** | alloc MiB | GC % |
|---|---:|---:|---:|---:|---:|---:|---:|
| 12322 midlat | 4 | 26.0 | 12.51 | 8.79 | **9.55** | 1704 | 4.4 |
| 76268 midlat | 4 | 34.8 | 12.51 | 9.17 | **9.89** | 1656 | 3.4 |
| 728 midlat | 5 | 47.1 | 33.25 | 9.96 | **11.05** | 1919 | 4.5 |
| 26176 arctic | 8 | 67.1 | 49.19 | 9.22 | **10.20** | 1903 | 4.6 |
| 122487 antarctic | 9 | −71.1 | 47.57 | 8.34 | **9.09** | 1731 | 4.1 |
| 121611 antarctic | 9 | −75.3 | 47.54 | 9.21 | **9.89** | 1774 | 3.3 |
| 121618 antarctic | 10 | −74.9 | 47.52 | 9.44 | **10.43** | 1809 | 4.3 |
| 123160 deep | 21 | −82.1 | 47.66 | 9.41 | **10.56** | 2229 | 4.7 |
| 122718 deep | 24 | −83.1 | 342.11 | 9.83 | **10.61** | 2452 | 3.6 |
| 123204 pole | 267 | −89.5 | 198.39 | 45.97 | **52.13** | 16679 | 5.9 |
| 123203 pole | 360 | −89.5 | 246.57 | 91.27 | **101.67** | 30597 | 5.1 |

Median of the three mid-latitude columns 9.89 core-s; of the four |lat| > 65
columns at 8–24 tiles, 10.32 core-s. **Intrinsic polar / mid-latitude = 1.04.**
Not 5×, not the 2–2.5× the brief expected to survive. Only the two dozen columns
that sit on the pole itself cost more (§4).

The production `secs` field is not a measure of column cost: the two 360-tile
columns with identical geometry recorded 246.57 s and 350.18 s, and column
122718 (24 tiles) recorded 342 s against an intrinsic 10.61.

## 2. Where the 5× came from: epoch mixing

The done log spans three separately-launched epochs with different code and GC
flags. Bucketing all 66 178 rows by epoch **and** by latitude band inside each
epoch:

| epoch | code / flags | columns | wall h | cols/s | median mid-lat s | median polar s | **polar / mid** |
|---|---|---:|---:|---:|---:|---:|---:|
| A 18:32–04:00 | subzone-store code, W=21, `-t 63`, `--gcthreads=8,1` | 27 404 | 9.47 | 0.804 | 28.63 | 23.28 | **0.81** |
| B 04:15–07:05 | perf-ladder `7627ee3`, cores=24 outer, `-t 26`, `--gcthreads=4,1` | 26 315 | 2.83 | 2.58 | 9.34 | 8.57 | **0.92** |
| C 08:13–13:55 | perf-ladder `e2f90d1`, cores=24 outer, `-t 26`, `--gcthreads=4` | 11 852 | 5.70 | 0.578 | 40.71 | 42.32 | **1.04** |

**Within every epoch, polar costs the same as mid-latitude or less.** The
47.6-vs-9.5 headline compares polar columns — which canonical Z7 order pushed
into the slow epoch C — against mid-latitude columns mostly retired in the fast
epoch B. Epoch B's 9.34 s median matches this campaign's quiet-box 9.9 core-s to
within 6 %, so epoch B *is* the healthy regime.

## 3. What actually regressed between B and C (not polar, and not solved)

B and C ran identical configuration (`cores=24 shape=outer batch=8
budget=1073741824 cache=3072 stripes=64`, `-t 26`) and differ in two things: the
commit and the GC flag.

**The commit is not the cause.** The entire `7627ee3 → e2f90d1` library delta is
defensive: `_fillwave!` now drains every spawned task instead of abandoning them
after the first failure, and three `PerChunk` observers (`show`, `nblocks`,
`storagebytes`) take the lock they were missing. No hot-path arithmetic changed.

**The workers were blocked, not slow.** The done log's worker ids show both
epochs ran 23 workers essentially never idle, and `ps` `pcpu` — a *lifetime*
average, differentiated here into instantaneous cores — shows what they drew:

| epoch | workers | occupancy | **instantaneous cores** (p10/p90) | core-s per column | external demand |
|---|---:|---:|---:|---:|---:|
| A | 21 | 100 % | 25.19 (22.1 / 28.2) | — | 24.2 |
| B | 23 | 99 % | **21.50** (20.4 / 22.6) | **8.3** | 18.8 |
| C | 23 | 99 % | **9.58** (8.3 / 11.3) | **16.5** | 31.4 (p90 67.7) |

Epoch C's 4.5× decomposes into **2.25× from drawing 9.6 cores instead of 21.5**
and **2.0× from burning 16.5 core-s per column instead of 8.3**. Queue
starvation is excluded — every worker held a column 99 % of the time. A
memory-bandwidth wall is excluded — that would show ~24 cores busy and slow.

Three candidate causes were tested and two were **refuted**:

- **Wave over-spawn** (polar columns spawn wider waves and poison the scheduler)
  — refuted in §8: deep-polar columns cost only 1.18× mid-latitude at fixed W on
  a quiet box.
- **Large heap forcing stop-the-world sweeps** — refuted in §8: pre-warming the
  tile cache to a 22.9 GiB live heap left cost flat at 10.69 core-s/column and
  *lowered* GC share to 12.9 %.
- **The GC flag itself** (`--gcthreads=4,1 → --gcthreads=4`, the concurrent page
  sweeper disabled after the SIGSEGV) — **not testable here**: this campaign is
  forbidden from running a nonzero second `--gcthreads` field, so the A/B could
  not be performed. It remains the leading unexplained candidate, alongside the
  1.7× higher external load epoch C ran under.

**This is the run's real cost driver and it is not a polar problem.** It is the
highest-value follow-up, and it needs a 23-worker / 24-core box to reproduce.

## 4. The degree tail

Cost is flat to 24 tiles then turns up, superlinearly:

| column | tiles | core-s | marginal above the 10.14 core-s flat baseline |
|---|---:|---:|---:|
| 123204 | 267 | 52.13 | 164 ms/tile |
| 123203 | 360 | 101.67 | 262 ms/tile |

**It is not discovery and not cache thrash.** The lazy plan's own accounting
(`GR.residency`) plus a standalone timing of `connectedchunks!`:

| column | tiles | discovery | plan build | loads | hold hits | dropped | source bytes |
|---|---:|---:|---:|---:|---:|---:|---:|
| 12322 midlat | 4 | 0.5 ms | 2.0 ms | 4 | 0 | 0 | 22.0 MiB |
| 123160 deep | 21 | 0.5 ms | 2.0 ms | 21 | 0 | 0 | 23.1 MiB |
| 123204 pole | 267 | 0.5 ms | 2.0 ms | 267 | 0 | 0 | 146.7 MiB |
| 123203 pole | 360 | 0.5 ms | 2.0 ms | 360 | 0 | 0 | 197.8 MiB |

Discovery costs **0.5 ms at 360 tiles** — 0.0005 % of the column. Every tile is
loaded exactly once; the 0.75 GiB data budget is never stressed. The pre-ladder
note's "discovery cost 2.7× useful work" is dead.

### How much of the run is in that tail

| tile degree | columns | share |
|---|---:|---:|
| 1–9 | 61 987 | 93.67 % |
| 10–24 | 3 755 | 5.67 % |
| 25–49 | 324 | 0.49 % |
| 50–99 | 88 | 0.13 % |
| 100–199 | 12 | 0.02 % |
| 200+ | 12 | 0.02 % |

Costing every column at the measured flat 10.14 core-s up to 24 tiles and
10.14 + 0.20 × (degree − 24) above — a *linear* model, which understates the two
pole columns and so is conservative **against** the tail — the whole run is
**187 core-hours** and the entire excess above the flat baseline is **0.62
core-hours, 0.33 %**. The 24 columns of degree ≥ 100 account for 0.24 core-h.

## 5. What sets a column's cost: loaded source cells

Joining `-t 1` cost to `residency` byte counts gives a single-variable model
holding across a 9× range of source volume and a 90× range of tile degree:

| column | tiles | MiB/tile | source Mcells | core-s | **µs / source cell** |
|---|---:|---:|---:|---:|---:|
| 12322 midlat 26 °N | 4 | 5.49 | 5.76 | 9.55 | 1.66 |
| 76268 midlat 35 °N | 4 | 5.49 | 5.76 | 9.89 | 1.72 |
| 728 midlat 47 °N | 5 | 5.49 | 7.20 | 11.05 | 1.53 |
| 26176 arctic 67 °N | 8 | 2.75 | 5.76 | 10.20 | 1.77 |
| 122487 71 °S | 9 | 1.83 | 4.32 | 9.09 | 2.10 |
| 121611 75 °S | 9 | 1.83 | 4.32 | 9.89 | 2.29 |
| 121618 75 °S | 10 | 1.83 | 4.80 | 10.43 | 2.17 |
| 123160 82 °S | 21 | 1.10 | 6.05 | 10.56 | 1.75 |
| 122718 83 °S | 24 | 1.10 | 6.91 | 10.61 | 1.54 |
| 123204 90 °S | 267 | 0.55 | 38.45 | 52.13 | 1.36 |
| 123203 90 °S | 360 | 0.55 | 51.84 | 101.67 | 1.96 |

**1.36–2.29 µs per loaded source cell, with no trend in latitude or degree.**

This kills the per-tile-overhead reading of §4. The Copernicus band table
(`src/systems/CopernicusDEM/bands.jl:35-51`) reduces tile width poleward — 1200
columns below 50°, then 800 / 600 / 400 / 240 / **120** — so a tile at 89.5 °S
holds 1200 × 120 = 144 000 cells, 0.55 MiB. At 1.7 µs/cell that tile is worth
**245 ms**, precisely the 164–262 ms "marginal cost per tile" of §4.
**Per-tile fixed overhead is not resolvable above the cost of the tile's own
pixels; it is consistent with zero.**

A pole column is expensive because it reads 51.8 M source cells against a
mid-latitude column's 5.76 M — nine times the data for the same 823 543
destination cells — not because it was handed 360 chunks instead of 4.

**Why the pole reads so much.** The band table's coarsest step is 10×, applied
across the whole 85–90° band, but at 89.5° meridian convergence is 115×. The
last band's pixels are 8 m × 92 m: **11× oversampled in longitude**. That
surplus is the pole columns' entire excess cost.

## 6. Polar CPU split vs mid-latitude: identical

`Profile.@profile` at 2 ms, restricted to samples inside `regrid_column`
(Julia was launched with `--gcthreads=4` and Profile samples every thread, so
~75 % of raw samples are the four idle GC threads parked in `pthread_cond_wait`
and must be excluded). Partition is leaf-first, first bucket matched wins, so
the columns sum to 100 %:

| phase | mid-lat (4 tiles) | polar 75–82 °S (9–21 tiles) | pole (360 tiles) |
|---|---:|---:|---:|
| clipping + sparse assembly | 42.09 % | **42.58 %** | 47.11 % |
| boundary generation | 30.64 % | **28.36 %** | 11.42 % |
| tree search | 26.22 % | **27.30 %** | 38.93 % |
| destination geometry | 0.65 % | 1.04 % | 0.10 % |
| tile decode / synthesis | 0.40 % | 0.71 % | 2.44 % |
| samples in `regrid_column` | 24 875 | 26 048 | 42 708 |

**The polar split is the mid-latitude split** — every phase within 2.3 points.
It shifts only at the 360-tile pole column, where tree search rises to 38.9 %
and boundary generation collapses to 11.4 %: boundary work is destination-side
over a constant 823 543 cells, so its share falls as the total grows 10×, while
tree search is per (destination, source-chunk) pair.

Cross-check against the brief's mid-latitude split: `snyder_inv*` inclusive is
**22.97 %** here against the brief's "~22 % snyder_inv_xyz boundary
generation" — an independent confirmation of the methodology. The brief's
"~10 % tree search" is narrower than this bucket, which also carries
`dual_depth_first_search` (28.1 % inclusive).

## 7. Candidate verdicts

### Candidate 1 — latitude-adaptive E-W source chunk merging: **REJECT**

Two independent bounds, and they agree.

*By tail share:* the whole excess of the degree tail above the flat baseline is
**0.62 of 187 core-hours, 0.33 %** (§4). 436 of 66 178 columns (0.66 %) are above
24 tiles at all.

*By mechanism:* §5 shows the tail's cost is its **pixels, not its chunk count** —
per-tile fixed overhead is consistent with zero. Merging 360 chunks into 45
leaves the same 51.8 M source cells to be paired, so it would not remove even
that 0.33 %.

No prototype was built. The second bound makes one pointless: the change does
not act on the quantity that sets the cost. At L13 with `La = 6` a column covers
one seventh the area and meets *fewer* 1° tiles, so the tail shrinks further.

**It does have a non-CPU justification, recorded separately:** the 360-tile
column churns **30.6 GiB of allocation** to produce 3.3 MB of output, 18× a
mid-latitude column. That is a GC-pressure and peak-memory hazard for the 64 GB
laptop target. Rank it as a memory fix, not a speed fix — and note that capping
wave width would achieve the same thing more cheaply.

### Candidate 2 — conservative E-W pixel pre-aggregation: **bounded at 13.7 % globally, 2.8 % of it polar**

Scope only, as commissioned, but the bound can be made exact: the only quantity
that matters is loaded source cells (§5) and the band table is a closed form.
Taking every column's tile set from the chunk-DAG adjacency, counting each
tile's cells from the band table, and comparing against the same tiles at
**square pixels for the tile's mid-latitude** (`min(1200 cos lat, current)`, so
it never adds cells):

| band | tile loads | cells (T) | square (T) | saving | core-h saved at 1.75 µs/cell |
|---|---:|---:|---:|---:|---:|
| < 50° | 168 209 | 0.242 | 0.215 | 11.4 % | 13.44 |
| 50–60 | 29 581 | 0.028 | 0.024 | 14.0 % | 1.93 |
| 60–70 | 38 094 | 0.027 | 0.023 | 16.4 % | 2.19 |
| 70–80 | 46 333 | 0.023 | 0.018 | 21.5 % | 2.38 |
| 80–85 | 22 159 | 0.007 | 0.005 | 33.4 % | 1.15 |
| 85–90 | 22 016 | 0.004 | 0.002 | **55.2 %** | 0.99 |
| **total** | 326 392 | **0.332** | **0.286** | **13.7 %** | **22.1** |

Ceiling **22.1 of 161 core-hours of source-side work (~12 % of a GLO-90/L12
run)** — but only **4.5 core-hours, 2.8 %, from |lat| > 60**. The polar bands
have the largest percentage saving and the smallest absolute one, because they
hold 5 % of the source cells. Copernicus already does most of this aggregation:
a polar column reads 25 % *fewer* source cells than the mid-latitude control.

Caveats, both against the estimate: `residency.bytes` counts whole loaded tiles,
so the 1.75 µs/cell constant absorbs work aggregation would not remove; and the
saving assumes aggregation all the way to square, which below 50° means
resampling the product's native grid.

**Not a polar measure.** Its real content is a global source-resolution decision
that breaks bit-identity, and it belongs in the next-store-generation
conversation on its own merits.

### Candidate 3 — c4 leaf-cap cache and wider c1 memo in the polar regime: **still null**

Answered from the profile, which has far better statistics than the timings
(between-process noise here is ~15 %, larger than the effect). Inclusive share
of the target frames, restricted to samples inside `regrid_column`:

| target | mid-lat | polar 75–82 °S | pole (360 tiles) |
|---|---:|---:|---:|
| **c4** source `child_indices_extents` → `LeafCells` | 6.31 % | **6.37 %** | 11.05 % |
| — its `_leafcell` cap derivation | 6.12 % | **6.17 %** | 10.75 % |
| **c1** `node_extent` (interior) | 8.76 % | **8.90 %** | 14.92 % |
| `_box_cap` (the un-memoised fallback) | 6.53 % | 6.63 % | 11.77 % |
| `_taskmemo` (memo lookup itself) | 0.26 % | 0.25 % | 0.25 % |

The `_leafcell` share is **6.17 % in the polar regime against 6.12 %
mid-latitude** — and against the perf-ladder's recorded 6.17 % for the same
frame, an exact reproduction. **Polar tile multiplicity does not revive c4**: its
ceiling is the same −6.4 % it was when the ladder dropped it, and it only
doubles on the 24 pole columns, 0.03 % of the run.

Throwaway patches were built and measured anyway (c1 `_MEMO_EXTENT_SLOTS`
1024 → 16384; c4 ported from the adopted `MemoRasterTree` template into
`MemoBlockCursor`, by value since `LeafCells` is `isbits`):

| column | tiles | base | c1 wide | c4 | c1 Δ | c4 Δ |
|---|---:|---:|---:|---:|---:|---:|
| 12322 midlat | 4 | 10.90 | 9.63 | 9.33 | −11.7 % | −14.4 % |
| 121611 | 9 | 10.08 | 11.37 | 9.57 | **+12.8 %** | −5.1 % |
| 123160 | 21 | 10.43 | 10.91 | 10.17 | +4.6 % | −2.5 % |
| 123204 | 267 | 54.28 | 52.63 | 50.98 | −3.0 % | −6.1 % |

c1's sign flips across columns — null. c4 is consistently negative but its
mid-latitude −14.4 % exceeds its own theoretical ceiling of −6.2 %, so the
timings are noise-dominated; the profile's −6.4 % ceiling is the real number.
Both patches were reverted; the worktree is clean at `e2f90d1`.

### Candidate 4 — cost-ordering of columns: dead, not revisited (per brief).

## 8. Concurrency, utilization, and the two refutations

### W workers on one cursor, quiet box, `shape = outer`, `-t 12`

| set | W | cores held | core-s/col | GC % | pauses |
|---|---:|---:|---:|---:|---:|
| mid | 1 | 1.40 | 10.12 | 4.6 | 363 |
| mid | 8 | 7.61 | **9.46** | 13.2 | 822 |
| ant | 1 | 2.05 | 10.94 | 6.8 | 404 |
| ant | 8 | 9.68 | **10.25** | 17.1 | 936 |
| deep | 1 | 2.42 | 10.86 | 9.5 | 387 |
| deep | 8 | 9.48 | **11.12** | 22.0 | 1056 |

Deep-polar costs **1.18×** mid-latitude at W=8. Wave over-spawn cannot produce
epoch C's 4.5× — **refuted**.

### Heap size, same 160 deep columns, same W=8

| pre-warmed tiles | live heap after GC | core-s/col | GC % | end RSS |
|---:|---:|---:|---:|---:|
| 0 | 2.13 GiB | 10.85 | 17.1 % | 11.6 GiB |
| 1 500 | 4.02 GiB | 10.71 | 13.0 % | 17.5 GiB |
| 3 000 | 12.10 GiB | 10.65 | 11.6 % | 25.3 GiB |
| 9 000 | 22.89 GiB | **10.69** | **12.9 %** | 20.0 GiB |

Cost is flat and GC share *falls* as the heap grows (more headroom, fewer
collections). A large heap alone does not stall the multi-worker regime —
**refuted**. Production's 40–58 GiB at 23 workers on 26 threads is still
untested; that needs a 24-core box.

### One column's utilization at `-t 12`

| column | tiles | production `outer` | default shape | wall speedup |
|---|---:|---|---|---:|
| 12322 midlat | 4 | 2.47 cores, 10.11 core-s | 7.37 cores, 9.61 core-s | 3.11× |
| 76268 midlat | 4 | 1.16 cores, 10.00 core-s | **7.54 cores**, 9.52 core-s | 6.82× |
| 26176 arctic | 8 | 1.54 cores, 9.94 core-s | 7.66 cores, 9.61 core-s | 5.13× |
| 121611 | 9 | 1.81 cores, 10.21 core-s | **7.73 cores**, 9.89 core-s | 4.41× |
| 122487 | 9 | 2.01 cores, 9.50 core-s | 7.14 cores, 9.16 core-s | 3.68× |
| 123160 | 21 | 2.63 cores, 10.45 core-s | **7.16 cores**, 10.19 core-s | 2.77× |
| 122718 | 24 | 2.78 cores, 11.45 core-s | **6.60 cores**, 10.88 core-s | 2.49× |
| 123204 | 267 | 6.24 cores, 51.49 core-s | 6.23 cores, 49.44 core-s | 1.04× |

**Correcting a briefed assumption: c3 *does* engage on multi-chunk polar
columns.** The 21- and 24-tile columns reach 7.16 and 6.60 cores under the
default shape. It fails to engage only on the 267-tile pole column, where the
wave is genuinely wide enough that the estimator prefers it.

Under production's `shape = outer` a column is deliberately capped at 1.2–2.8
cores and the box is filled with W=23 workers instead. Inner threading reaches
7.5 cores at the **same or slightly fewer core-seconds** (9.52 vs 10.00), so it
is a free 2.5–6.8× latency win. Polar columns hold more cores under `outer`
(2.0–2.8 vs 1.2–1.5) because more source chunks make a wider wave.

## 9. Recommendation, ranked

The 53.4 h GLO-30/L13 projection is computed *from* `-t 1` intrinsic rates
(1 281 core-h at 24 cores) — that is, it already assumes epoch-B-class
efficiency. **There are no hours to take off it from polar work. There are ~187
hours to avoid adding by not repeating epoch C** (53.4 h × 4.5 = 240 h).

1. **Find and fix the epoch-C collapse (≈4.5×, i.e. 53 h → ~240 h if repeated).**
   Not a polar problem; wave over-spawn and heap size are refuted (§8). Leading
   remaining candidate is the concurrent-page-sweeper flag, which this campaign
   is forbidden from A/B-testing, compounded by 1.7× external load. **Do not
   launch GLO-30 until a 2 000-column canary at W=23 on 24 cores reproduces
   epoch B's 8.3 core-s/column.** This is the only item on the list worth more
   than a few percent.
2. **Switch `shape` from `outer` to `inner` (W ≈ 8 rather than 23).** Free on
   throughput (same core-seconds), 2.5–6.8× lower per-column latency, and it cuts
   in-flight columns ~3×, which is the dominant term in RSS — 23 × 1 GiB of
   per-worker weight budget becomes 8 × 1 GiB. Directly serves the 64 GB laptop
   target. Verify at 24 cores before adopting; §8's numbers are at `-t 12`.
3. **Cap wave width on high-degree columns.** Removes the 30.6 GiB allocation
   spike on the 24 pole columns for a few lines. Memory fix, ~0 % speed.
4. **Candidate 2, if and only if a new store generation is being cut anyway** —
   ~12 % globally, 2.8 % from the polar bands, breaks bit-identity.
5. **Candidates 1, 3 and 4: closed.** ≤0.33 %, null, and dead respectively.

## 10. Reproducing

```bash
git worktree add --detach /home/asinghvi17/geo/DGG-polarprof e2f90d1
cp /home/asinghvi17/geo/DGG-perf-ladder/Manifest.toml /home/asinghvi17/geo/DGG-polarprof/
# bench/ is a workspace member, so the manifest lives at the worktree ROOT.
# A fresh instantiate re-resolves ConservativeRegridding to a tree without
# `Trees.split_weight` (2026-08-21-chunk-dag-sim.md §9). CR must read 0.2.9 @ 6a4b997.

S=<scratchpad>/pp
head -n -1 /home/asinghvi17/geo/DGG-polarprof/scripts/copdem_production.jl > $S/prod.jl
cp /home/asinghvi17/geo/dggstores/copdem90-igeo7-l12-synthetic.zarr.columns.txt $S/cols.txt
cd /home/asinghvi17/geo/DGG-polarprof
RASTERDATASOURCES_PATH=/home/asinghvi17/geo/DiscreteGlobalGrids.jl/bench/data \
  nice -n 10 julia --project=bench -t 1 --gcthreads=4 $S/s1_intrinsic.jl
```

**Never a nonzero second `--gcthreads` field.** `--gcthreads=4` only; the
concurrent page sweeper is what the 2026-08-21 SIGSEGV was
(`2026-08-21-polar-segfault.md`), and this harness bypasses `main`, so the
script's own `gcguard` never runs to protect it.

| script | what it produced |
|---|---|
| `setup.jl` | the production world without `main()`; `world`, `lazycolumn`, `timecolumn`, `withcpu`, `emit`, the reference column table |
| `s1_intrinsic.jl` | §1 — 11 reference columns × 3 reps at `-t 1` |
| `s2_phase.jl` + `s2b_offline.jl` + `s2c_watch.jl` | §6, §7 c3 — profile, then offline partition and watched frames restricted to `regrid_column` |
| `s3_concurrent.jl` | §8 — W workers on one cursor, `shape = outer`, per band |
| `s4_tiles.jl` | §4–5 — `GR.residency` and standalone `connectedchunks!` timing |
| `s5_cand3.jl` | §7 — c1 memo width and the c4 leaf-cap memo |
| `s6_heap.jl` | §8 — the same W=8 sweep with the tile cache pre-warmed |
| `s7_util.jl` | §8 — one column at `-t 12`, `OUTER_PARALLEL` true vs false |

Column selection joined the production done log to the chunk-DAG adjacency
(`2026-08-21-chunk-dag-adjacency.ndjson`), taking each tile's latitude from its
Copernicus stem for a per-column mean latitude and tile degree.

The concurrency runs wrote to a fresh
`/home/asinghvi17/geo/dggstores/polarprof-scratch/pp.zarr`. The production store,
its done log, its column cache, and the worktrees `DGG-perf-ladder`,
`DGG-chunk-dag`, `DGG-snyder-impl` and `DGG-verify` were not written to.
`DGG-polarprof` is left clean at `e2f90d1` — the c1/c4 patches of §7 are applied
and `git checkout`-reverted inside `run_s5.sh`.
