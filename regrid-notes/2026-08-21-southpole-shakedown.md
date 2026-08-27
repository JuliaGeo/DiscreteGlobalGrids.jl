# South-pole DAG-executor shakedown

## Outcome

The throughput-collapse regime was **not reproduced**.  All three safe runs are
much nearer the healthy epoch-B reference (8.3 core-s/column) than the collapsed
epoch-C reference (16.5 core-s/column): 9.179, 9.422, and 10.018
core-s/column.  Changing only the first `--gcthreads` field from 4 to 8 moved
core cost by +2.65%, far below both the campaign's 15% trigger and the old 2x
regime split.  Optional config D (`--gcthreads=16`) was therefore not run.

The DAG residency result is strong.  Peak decoded-tile residency was 0.70-0.95
GiB, every demanded tile loaded at most once, all tiles retired, and there were
zero uncredited demands and zero downloads.  Kernel high-water RSS was 19.38
GiB (A), 23.76 GiB (B), and 10.67 GiB (C), all comfortably below 64 GiB and
far below the old 83.7 GiB S2/W=21 LRU-era peak.

The production-store cross-check formally **fails** its 1e-3 m gate, but the
failure is localized to an invalid production pole artifact, not nondeterminism
in the new executor.  A, B, and C are bit-identical in both affected pole
columns.  Production alone contains values outside the synthetic source's
exact global range; details are in the cross-check section.

Harness: `scripts/southpole_shakedown.jl`, commit
`c8287c563467795416ef83fafb7739fda790cc40`, pushed append-only to
`origin/claude/perf-ladder`.

## Campaign definition and safety

- Worktree/branch: `/home/asinghvi17/geo/DGG-perf-ladder`,
  `claude/perf-ladder`, starting at `22cf2cf`.
- Julia 1.12.6, `-t 26`, `cores=24`, `workers=0`, `nice -n 10`.
- A/B resolve to W=23 with `shape=:outer`; C resolves to W=8 with
  `shape=:inner`.
- Synthetic source, `real=:none`, no network provider, `prefetch=0`.
- New-stack defaults retained: affinity order, `RefCountCache`, guided taper,
  batch 8, 1 GiB per-worker lazy budget, 32 MiB malloc trim threshold.
- Every launch carried only a first-field GC setting.  `gcsweep=0` was checked
  before each run and independently enforced by the production driver's
  `gcguard`.  The concurrent sweeper was never tested.
- Each run had a 60-minute `SIGINT` timeout plus a five-minute forced-kill
  backstop.  None approached it.
- Stores and logs were written only below
  `/home/asinghvi17/geo/scratch-stores/`.  The production store was opened only
  with Zarr mode `"r"` during the cross-check.
- Runs were sequential: A, then B, then C.  A preflight stopped before creating
  a store because the worktree-local data root lacked the Natural Earth
  shapefile.  The harness was pinned to the already-resolved read-only data tree
  in the main checkout, the unused scratch column sidecar was removed, and A
  was relaunched.  The preflight error remains in the NDJSON and in
  `southpole-A-gc4-outer-preflight-error.log`; final results are the
  `run_final` records.

The box was shared.  Observed one-minute load averages during the runs were
roughly 25-42 on 64 cores, including long-lived unrelated jobs.  Wall rates are
therefore secondary to process CPU (`utime+stime`) and its core-s/column ratio.

## Exact column set

The production covering contains 66,178 level-5 columns.  Of those, 5,380 have
centroid latitude below -60 degrees (range -89.735117 to -60.009448 degrees), so
the harness selected exactly 2,000:

1. retain all 1,087 production-ledger `nan==0` columns below -80 degrees;
2. retain the two southernmost topological IGeo7 pentagons, columns 112049 and
   154067;
3. sort the remaining eligible columns by `(latitude, column)` and take 911
   deterministic evenly spaced ranks.

The two pentagons are an explicit exception to the strict latitude band: their
centroids are both -58.282525588539 degrees.  No topological IGeo7 pentagon lies
below -60 degrees, so this exception is necessary to satisfy the requested
pentagon coverage.  Neither pentagon belongs to the production land covering;
both stores correctly read them as fill.

The exact 2,000-element integer array and the complete 1,087-element mandatory
deep-land array are in the NDJSON `column_set` record.  The canonical comma-join
of the selected integers has SHA-256:

```text
62cef2b46fcac83f308f4a6314e8acd80ebae71f9b96b72fe8d3380184dd1635
```

Every config's `run_start` and `run_final` repeats that digest.

## Results

`wall` is harness wall including setup; `compute` is inferred from the earliest
ledger start through the final completion.  `cores` is aggregate process CPU
divided by harness wall.  The p10/median/p90 values are the differentiated
five-second `/proc/<pid>/stat` samples.  Rates use the 1,646,811,486 cells
actually written (the two pentagon subtrees are shorter than a hexagon column).
Peak RSS is kernel `VmHWM`, not a sampled current-RSS approximation.

The sampler requested a five-second cadence and every record carries its actual
interval.  Under W=23, the in-process sampler task was itself scheduler-delayed:
A's interval median/p90/max was 6.32/9.26/19.68 s and B's was
6.33/9.80/25.69 s.  C, with W=8, achieved 5.07/5.23/6.86 s.  Instantaneous
cores always use the actual interval's `/proc/self/stat` tick delta, so the
values remain differentiated process CPU rather than `ps` lifetime averages,
but A/B's percentiles are more smoothed than a strict five-second external
sampler would be.  Total CPU, mean cores, and core-s/column use run endpoint
deltas and are unaffected by this cadence slippage.

| config | shape / W | GC mark | wall / compute | cells/s | core-s/col | cores mean; p10/med/p90 | **peak RSS** | cache peak | loads / uncredited / cold | last-10 span | regime |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| A | outer / 23 | 4 | 18.14 / 17.98 min | 1,526,819 | **9.179** | 16.87; 15.98/17.08/17.95 | **19.38 GiB** | 0.932 GiB / 1,120 tiles | 6,775 / 0 / 0 | 21 s | epoch-B-like |
| B | outer / 23 | 8 | 17.81 / 17.65 min | 1,555,357 | **9.422** | 17.64; 16.96/17.85/18.66 | **23.76 GiB** | 0.951 GiB / 1,117 tiles | 6,775 / 0 / 0 | 24 s | epoch-B-like |
| C | inner / 8 | 4 | 23.79 / 23.63 min | 1,161,608 | **10.018** | 14.04; 13.49/14.34/15.20 | **10.67 GiB** | 0.699 GiB / 920 tiles | 6,775 / 0 / 0 | 69 s | epoch-B-like |

All configs computed 2,000/2,000 columns with zero driver failures.  All loaded
the same 6,775 source tiles exactly once, ended with zero live/pinned tiles, and
reported zero uncredited demands.  A/B had zero joined loads; B/C each had two,
which is ordinary single-flight sharing and not a reload.

### Regime verdict and first-field GC effect

- A is 10.6% above the 8.3 epoch-B reference and 44.4% below the 16.5 epoch-C
  reference.
- B is 13.5% above epoch B and 42.9% below epoch C.
- C is 20.7% above epoch B and 39.3% below epoch C.  It is still nearer B by a
  wide margin and never showed the old 9.6-core/16.5-core-s combination.
- B versus A: +2.65% core-s/column, +4.54% mean cores, -1.83% compute wall,
  +1.87% cells/s, and +22.6% peak RSS.  This one-replicate difference is small
  and internally mixed; it does not support the mark-thread first field as a
  regime switch.
- D's >15% condition was not met.

The collapse is therefore ruled out for this Antarctic canary under safe
single-field GC configurations.  This does **not** test whether the unsafe
second-field sweeper caused the old epoch split; that configuration remains
forbidden.

### Inner shape

C confirms the memory side of the earlier polar profile but not its
throughput-neutral claim at this launch shape.  Versus A, inner used 44.9% less
peak RSS and 25.0% less cache residency, but 9.1% more core-s/column and 23.9%
fewer cells/s.  The driver warned at launch that `shape=:inner` saturates near
four Julia threads per worker: W=8 would like roughly 32 threads, while this
campaign retained the production `-t26` process shape.  It consequently drew
only 14.0 mean cores.  The 69-second last-10 span is the same underfill made
visible at the expensive pole tail.

This is a measurement result, not a request to change the driver.  A separate
`-t32`/inner or W=6 sizing run is needed before calling inner
throughput-neutral at a 24-core target.

## RSS and the 64 GiB question

All three runs fit 64 GiB with large headroom: 44.62 GiB for A, 40.24 GiB for
B, and 53.33 GiB for C.  Relative to 83.7 GiB, peaks fell 76.8%, 71.6%, and
87.2%, respectively.

The tile cache is no longer the RSS floor.  Its peak is below 1 GiB and it is
empty at shutdown.  Nevertheless shutdown RSS is 8.67, 7.93, and 8.37 GiB for
A/B/C, against final Julia GC-live sizes of only 1.79, 1.92, and 1.25 GiB.
That leaves a roughly 6-7 GiB non-live residual with no cache tiles: the fixed
445 MiB land mask and runtime/Zarr state plus allocator-retained pages from the
regrid working buffers.  During compute, per-worker lazy weight/scratch working
sets and their allocator high-water mark dominate the transient peak.  Reducing
in-flight workers from 23 to 8 cuts that component sharply, which is why C peaks
at 10.67 GiB even though all configs load the same tiles.

The five-second current-RSS samples missed short peaks (sampled maxima 17.08,
18.66, and 9.89 GiB), so the table correctly uses exact kernel `VmHWM` values.

## Production-store cross-check

Twenty-four exact Zarr columns were read from config A and the production store:
two pentagons, six below -80 degrees, and sixteen spread through the remaining
band.  NaN masks match in every cell.

| column | latitude | class | max abs residual (m) | non-bit-equal Float32 | NaN-mask mismatch |
|---:|---:|---|---:|---:|---:|
| 112049 | -58.282526 | pentagon | 0 | 0 | 0 |
| 114833 | -60.009448 | band | 0 | 0 | 0 |
| 115619 | -79.999087 | band | 0 | 0 | 0 |
| 121338 | -80.010794 | south of -80 | 0 | 0 | 0 |
| 121404 | -75.375343 | band | 0 | 0 | 0 |
| 121503 | -81.214597 | south of -80 | 3.051758e-5 | 1 | 0 |
| 121641 | -76.253030 | band | 0 | 0 | 0 |
| 121654 | -73.639439 | band | 0 | 0 | 0 |
| 123003 | -85.805880 | south of -80 | 2.384186e-7 | 2 | 0 |
| **123203** | **-89.735117** | **south of -80 / pole** | **3,301.082123** | **8,825** | **0** |
| 123481 | -69.934350 | band | 0 | 0 | 0 |
| 123546 | -71.881576 | band | 0 | 0 | 0 |
| 154067 | -58.282526 | pentagon | 0 | 0 | 0 |
| 156654 | -66.256770 | band | 0 | 0 | 0 |
| 157532 | -77.151687 | band | 0 | 0 | 0 |
| 157634 | -74.507485 | band | 6.103516e-5 | 1 | 0 |
| 157767 | -78.037876 | band | 0 | 0 | 0 |
| 157870 | -68.938059 | band | 0 | 0 | 0 |
| 157897 | -67.767404 | band | 0 | 0 | 0 |
| 158043 | -72.761212 | band | 0 | 0 | 0 |
| 158205 | -70.934102 | band | 0 | 0 | 0 |
| 164719 | -83.976298 | south of -80 | 0 | 0 | 0 |
| 164831 | -82.536904 | south of -80 | 0 | 0 | 0 |
| 164945 | -78.986181 | band | 0 | 0 | 0 |

Aggregate: 19,765,032 stored Float32 positions, 8,829 non-bit-equal values,
zero NaN-mask mismatches, and maximum residual 3,301.082123 m.  Thus the literal
production-store gate fails.  If column 123203 is excluded, only four cells are
non-bit-equal and the maximum residual is 6.103516e-5 m, comfortably below the
1e-3 m threshold.

### Pole diagnosis

The cross-check failure is reproducible and isolated:

| column | A vs B bit differences | A vs C | production differences | >1e-3 m | max residual | production values outside [-1400,1600] m | fresh outside range |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 123203 | 0 | 0 | 8,825 | 8,800 | 3,301.082 m | 7,754 | 0 |
| 123204 | 0 | 0 | 1,322 | 1,308 | 3,270.526 m | 1,233 | 0 |

The analytic field is
`1000 sin(3lon) cos(2lat) + 500 cos(7lon) sin(5lat) + 100`, so its exact global
range is [-1400, 1600] m.  A valid positive-weight conservative mean cannot
leave that range.  At the worst row of 123203 (1.8931 E, 89.7812 S), A/B/C all
hold -485.54794 m and the analytic centroid is -485.55489 m, while production
holds 2815.5342 m.  At the worst row of 123204, fresh is -486.47617 m versus an
analytic -486.47829 m, while production is 2784.05 m.

This rules out scheduler order, GC first-field count, and inner/outer shape as
the source of the discrepancy.  It does not yet distinguish corruption in the
old store from an invalid old pole-weight path.  The production pole chunks
must be treated as suspect and regenerated before they are used as a numerical
oracle.

## Carry into the real GLO-90 run

1. Keep single-field GC settings only.  `--gcthreads=4` is the conservative
   choice from this campaign: it had the lowest core cost and lower peak RSS
   than 8, while the throughput difference was noise-sized.
2. Carry the affinity/refcount/taper DAG stack.  It held decoded tiles below 1
   GiB, loaded each demanded tile once, retired fully, and removed the old RSS
   ramp.
3. Use outer/W=23 at `-t26` when throughput is the priority.  Inner/W=8 is the
   clear memory choice but needs enough Julia threads or different worker
   sizing before being adopted as throughput-neutral.
4. The 64 GiB laptop-fit question is answered yes for this synthetic Antarctic
   canary, including outer.  The memory floor is now working buffers and
   allocator/runtime retention, not decoded-tile residency.
5. Do not use production columns 123203/123204 as numerical references until
   regenerated.  Retain the rest of the 24-column comparison as evidence: 23
   sampled columns passed the 1e-3 m gate and all NaN masks matched.
6. This synthetic campaign says nothing new about network prefetch.  Real
   GLO-90 should use the separately measured fetch policy; cold-download
   behavior was deliberately absent here.

## Open questions

- Was the old pole artifact persisted corruption, or an invalid historical
  pole-weight calculation?  Regenerating just 123203/123204 with an independent
  oracle is the next diagnostic.
- Does inner recover outer throughput at `-t32`, or should a 26-thread launch
  resolve to W=6 instead of W=8?  This campaign measured the current driver
  behavior and did not tune it.
- A/B were one replicate each on a shared box.  The 2.65% core-cost difference
  is not large enough to establish a mark-thread optimum, only enough to rule
  out the old 2x collapse.

Raw record: `regrid-notes/2026-08-21-southpole-shakedown.ndjson`.  It contains
the exact column arrays, every requested five-second CPU/RSS/GC-live sample with
its actual elapsed interval, run summaries, all 24 cross-check rows, and both
pole diagnostic rows.
