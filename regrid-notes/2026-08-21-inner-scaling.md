# CopDEM DAG driver inner-shape thread scaling

## Outcome

`shape=:inner` is **not throughput-neutral at 32 Julia threads**.  With the
driver's default fixed 24-core sizing budget, outer resolved to W=23 and inner
to W=8 at both thread counts.  At `-t32`, default inner used 10.380
core-s/column versus outer's 9.381 (+10.7%), drew 14.61 versus 17.72 mean cores
(-17.5%), and took 23.68 versus 17.65 minutes (+34.2%).  The extra threads did
not recover the t26 gap.

The requested W boost was available through the existing `config.workers`
override.  Inner/W=15 improved over inner/W=8 by 5.5% in wall and 1.5% in core
cost, but still took 26.8% longer than outer and drew only 15.23 cores.  It is
not evidence that a larger W can make inner neutral; in fact, the very small
utilization gain from seven more workers argues against W being the principal
limit.

At the deployment point, **use outer at `-t21`** when wall time is the primary
criterion.  It finished in 22.88 minutes versus inner's 27.72 minutes, a 4.85
minute / 17.5% wall saving.  Its exact kernel peak was 18.45 GiB and its
shutdown RSS floor was 7.93 GiB.  Inner is the memory fallback at 9.90 GiB
peak and 6.60 GiB shutdown, but it is 21.2% slower relative to outer.

All five configurations completed before the 75-minute limit, with zero driver
failures.  All required DAG cache checks passed.  The five stores were
bit-identical on the 12-column correctness sample, including both pentagons.

## Campaign definition and guards

- Worktree/branch: `/home/asinghvi17/geo/DGG-perf-ladder`,
  `claude/perf-ladder`, starting at `fe16005`.
- Harness: `scripts/southpole_shakedown.jl`; scaling support commit
  `0d0f4d9f748720d3be4cdeda176cd93fd1f7a959`, pushed append-only to
  `origin/claude/perf-ladder`.
- Julia 1.12.6; `nice -n 10`; `--gcthreads=4` with no second field.  Every
  process reported `gcmark=4`, `gcsweep=0`.
- The prior campaign's fixed `cores=24` budget was retained so the only primary
  changes were Julia thread count and driver shape.  The W-boost run alone set
  an explicit worker count.
- Synthetic source only, `real=:none`, `prefetch=0`, no network provider and no
  downloads.  The source inventory reported zero local GeoTIFF overrides.
- Fresh stores were created as
  `/home/asinghvi17/geo/scratch-stores/innerscale-<config>.zarr`.  Unrelated
  pre-existing `prefetch-coldtest-*` artifacts in the scratch root were left
  untouched.
- Runs were strictly sequential in this order: outer-t32, inner-t32,
  inner-t32-Wboost, outer-t21, inner-t21.
- Each launch used an external 75-minute `SIGINT` timeout and five-minute
  forced-kill backstop.  No timeout fired.
- Telemetry differentiated `/proc/<pid>/stat` utime+stime ticks at a requested
  five-second cadence and read current RSS plus kernel `VmHWM`.  No `ps pcpu`
  value was used.

The original production column and completion sidecars had been removed since
the shakedown.  The harness therefore recovered the complete canonical column
and mandatory-deep-land arrays from the prior shakedown NDJSON, then freshly
checked count, uniqueness, containment, pentagon topology, and the comma-join
SHA-256 before any run.  The verified set is exactly 2,000 columns with digest:

```text
62cef2b46fcac83f308f4a6314e8acd80ebae71f9b96b72fe8d3380184dd1635
```

The pentagons are columns 112049 and 154067.  A first selector preflight stopped
before creating a store when it found the missing source sidecars; no column
set was substituted without the digest check.

## Results

`wall` includes setup.  `compute` runs from the earliest ledger start through
the final completion.  `cores` is endpoint process CPU divided by harness wall.
The p10/median/p90 values are differentiated process-CPU samples.  Tail is the
completion span of the final ten columns.  Every run wrote the same
1,646,811,486 cells.

| config | shape / threads / W | wall / compute | cells/s | core-s/col | cores mean; p10/med/p90 | peak RSS | shutdown RSS | tail |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| outer-t32 | outer / 32 / 23 | 17.65 / 17.33 min | 1,583,412 | **9.381** | 17.72; 16.57/18.12/19.23 | **24.60 GiB** | 9.08 GiB | 13 s |
| inner-t32 | inner / 32 / 8 | 23.68 / 23.39 min | 1,173,486 | **10.380** | 14.61; 13.78/14.95/15.93 | **10.90 GiB** | 6.87 GiB | 73 s |
| inner-t32-Wboost | inner / 32 / 15 | 22.37 / 22.07 min | 1,243,552 | **10.223** | 15.23; 14.45/15.54/16.53 | **13.54 GiB** | 7.70 GiB | 59 s |
| outer-t21 | outer / 21 / 23 | 22.88 / 22.59 min | 1,215,108 | **9.411** | 13.71; 12.96/13.88/14.70 | **18.45 GiB** | 7.93 GiB | 9 s |
| inner-t21 | inner / 21 / 8 | 27.72 / 27.41 min | 1,001,473 | **10.077** | 12.12; 11.60/12.32/13.00 | **9.90 GiB** | 6.60 GiB | 39 s |

The in-process sampler was delayed most under the outer worker pool.  Sample
interval median/p90/max was 5.42/6.92/14.32 s for outer-t32 and
5.70/7.20/18.76 s for outer-t21.  The three inner runs were near the requested
cadence: medians 5.07-5.10 s and p90 5.27-5.40 s.  Instantaneous cores use each
sample's actual interval.  Endpoint CPU, mean cores, core-s/column, and kernel
peak RSS are unaffected.  Short RSS excursions were material: outer-t32's
sampled maximum was only 17.84 GiB against exact `VmHWM` 24.60 GiB, and
outer-t21 sampled 15.61 GiB against exact 18.45 GiB.

## DAG cache invariants

Every configuration demanded and loaded the same 6,775 source tiles exactly
once.  Reloads, uncredited demands, cold downloads, live tiles at end, and
pinned tiles at end were all zero.

| config | cache peak | loads / demanded / reloads | uncredited / cold / live / pinned | joined loads |
|---|---:|---:|---:|---:|
| outer-t32 | 0.922 GiB / 1,010 tiles | 6,775 / 6,775 / 0 | 0 / 0 / 0 / 0 | 2 |
| inner-t32 | 0.669 GiB / 900 tiles | 6,775 / 6,775 / 0 | 0 / 0 / 0 / 0 | 0 |
| inner-t32-Wboost | 0.839 GiB / 1,044 tiles | 6,775 / 6,775 / 0 | 0 / 0 / 0 / 0 | 2 |
| outer-t21 | 0.893 GiB / 1,019 tiles | 6,775 / 6,775 / 0 | 0 / 0 / 0 / 0 | 1 |
| inner-t21 | 0.687 GiB / 896 tiles | 6,775 / 6,775 / 0 | 0 / 0 / 0 / 0 | 1 |

## t32 throughput-neutrality verdict

The answer is **no** for both default and boosted inner shapes.

- Default inner versus outer: +10.65% core-s/column, -17.54% mean cores,
  -25.89% cells/s, and +34.20% wall.  Its memory benefit is real: -55.67% peak
  RSS and -24.30% shutdown RSS.
- W=15 versus W=8: -1.51% core-s/column, +4.26% mean cores, +5.97% cells/s,
  and -5.54% wall, at +24.14% peak RSS.
- W=15 versus outer: +8.98% core-s/column, -14.03% mean cores, +26.77% wall,
  and -44.97% peak RSS.

The t26 shakedown measured outer/W=23 at 9.179 core-s/column and 16.87 mean
cores, and inner/W=8 at 10.018 and 14.04.  At t32 those became 9.381/17.72 and
10.380/14.61 respectively.  Outer gained 0.85 mean cores and 2.7% wall; inner
gained only 0.57 mean cores and 0.5% wall.  Core cost moved +2.2% for outer and
+3.6% for inner, small run-to-run changes in the wrong direction for a claimed
recovery.  Six more Julia threads therefore did not change the regime or close
the shape gap.

## W-resolution story

The current formula is a function of `config.cores`, not of
`Threads.nthreads()`:

```text
outer: W = ceil(cores / 1.06)
inner: W = ceil(cores / 3.1)
```

With this campaign's retained `cores=24`, W is therefore 23 outer and 8 inner
at both `-t32` and `-t21`.  `config.workers=15` supplied the supported W-boost
override.  If the budget were instead explicitly set equal to the Julia thread
count, the same formula would resolve to outer/inner W=31/11 at t32 and W=20/7
at t21; the driver does not make that substitution automatically.

There are idle Julia threads at both operating points, but the measurements do
not support one universal W correction:

- At t32, outer/W=23 drew 17.72 cores against the configured 24-core target.
  Pure proportional extrapolation gives W about 31 to draw 24 cores
  (`23 * 24 / 17.72`), conveniently the same order as a 32-core formula
  budget.  W=30-31 is the defensible next outer sizing measurement, not a
  validated fix; memory and near-`nthreads` contention can break linearity.
- At t32, inner/W=8 already equals the driver's `nthreads/4` saturation
  heuristic.  Raising W by 87.5% to 15 increased mean utilization by only 0.62
  core.  Linear extrapolation of that observed slope would require an absurd
  W about 43 merely to match outer's 17.72 cores, beyond the 32-thread pool.
  No practical W-only fix is supported.  W=8 remains the efficient memory
  point; W=15 buys only a small wall improvement.
- At t21, outer/W=23 already exceeds the Julia thread count and inner/W=8
  exceeds the heuristic saturation point of about W=5.  Neither shape is
  worker-starved there, so increasing W is not a plausible utilization fix.
  A deployment-specific 21-core budget would produce W=20 outer and W=7 inner;
  W around 5 may be enough for inner by the heuristic, but those smaller-W
  memory/contention points were not measured and should not be assigned the
  measured peak-RSS numbers.

The fixed-budget coefficients therefore deserve a later sizing sweep, but this
campaign supplies no basis for a production code change.  In particular, the
inner shortfall should not be described as merely a low default-W cap.

## t21 deployment recommendation and 64 GB headroom

Outer-t21 wins wall by 290.8 seconds.  Inner spends 7.08% more CPU per column,
draws 11.64% fewer cores, and delivers 17.58% fewer cells/s.  It saves 8.55 GiB
of peak RSS (46.32%) and 1.34 GiB at shutdown (16.86%).

Recommend **`shape=:outer`, Julia `-t21`, `--gcthreads=4`, current W=23** for
the 21-thread user machine, with an expected synthetic-canary peak of about
**18.5 GiB** and a post-compute floor of about **7.9 GiB**.  Against the
campaign's 64 GiB convention that leaves 45.55 GiB at peak.  A marketed 64 GB
machine exposes about 59.6 GiB before firmware reservations, still leaving
roughly 41 GiB for the OS, browser, and other processes at this measured peak.
Even a 20 GiB co-tenant working set would retain about 21 GiB of margin.

This is enough headroom to prefer outer for throughput, but the exact kernel
peak—not the sampled RSS—must be used for capacity planning.  If the user's
browser/OS workload is exceptionally memory-heavy or swap avoidance dominates
latency, inner-t21 is the explicit fallback: about 9.9 GiB peak and 6.6 GiB
shutdown, in exchange for 4.85 minutes on this canary.

## Bit identity

Outer-t32 was compared against every other store on these 12 deterministic
columns:

```text
112049, 114833, 115619, 115717, 121338, 121993,
123033, 123203, 154067, 157714, 163382, 164693
```

The sample contains both pentagons 112049 and 154067.  Across four candidate
comparisons and 39,530,064 stored Float32 positions there were **zero bit
differences**.  The STOP condition did not fire.

## Carry-forward

1. For the real GLO-90 run and the eventual 21-thread machine, start from
   outer/`-t21`/single-field GC4 when throughput matters.  Budget about 18.5
   GiB peak and 7.9 GiB shutdown for this synthetic Antarctic canary, then
   confirm with a real-source canary before treating those as hard limits.
2. Keep inner-t21 as the memory-pressure fallback, with the measured 9.9 GiB
   peak and the measured 21.2% wall penalty relative to outer.
3. Do not expect `-t32` alone to make inner throughput-neutral.  W=15 also does
   not close the gap.
4. Keep the affinity/refcount/taper stack and the no-concurrent-sweeper rule.
   Every run retained one-time tile loads and fully retired the DAG cache.
5. Real GLO-90 still needs the separately validated fetch/prefetch policy.  This
   campaign deliberately had zero downloads and does not measure network time
   or network-related memory.
6. If worker sizing is revisited later, measure outer W=30-31 at t32 and smaller
   deployment-budget counts around outer W=20 / inner W=5-7.  These are
   measurement candidates, not changes justified by this campaign.

Raw telemetry and checks are in
`regrid-notes/2026-08-21-inner-scaling.ndjson`.  Per-run stdout logs and stores
are under `/home/asinghvi17/geo/scratch-stores/`.
