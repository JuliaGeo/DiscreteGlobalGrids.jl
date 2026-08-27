# CopDEM DAG driver t21 utilization attribution and tuning

## Outcome

The old `outer/-t21/W=23/gc4` result is reproducible.  On the new 800-column
expensive Antarctic set it drew **13.64 cores** versus the earlier 2,000-column
result's 13.71, took 10.25 minutes wall, and peaked at 20.18 GiB.

The missing cores are now attributed:

- Julia scheduler thread 1 is parked in the driver's `Threads.@sync` and is
  essentially unused (0.20% busy in the steady window).
- Threads 2--21 are not split into busy and idle populations.  Their aggregate
  busy fractions are tightly clustered at 63.77--64.96%, while GC occupies
  36.0% of baseline wall.  They work outside GC and stop together for GC.
- The 20 worker threads lose about 7.14 average cores to those stops.  The four
  GC TIDs do about 1.22 cores of GC work during the same window, making the net
  GC loss about **5.92 cores**.  The parked main thread loses another **1.00**.
  That accounts for 6.92 of the 6.91-core steady-state gap; setup and tail lower
  the whole-run mean by a further 0.44 core, giving the observed 7.36-core gap.

This is primarily **stop-the-world allocation/GC time plus one parked scheduler
thread**, not seven worker threads starved of columns and not a Zarr, done-log,
cache, or pull-cursor lock.

Increasing W helps, but no allowed safe configuration gets close to all 21
cores.  The completed configuration with the highest whole-run draw was
W32/GC8 at 14.80 cores; it was slower and more CPU-expensive than W40/GC4.
W40/GC8 reached the 40 GiB guard after 358 columns and was cleanly aborted.

For throughput, recommend **outer, `-t21`, `workers=40`, single-field
`--gcthreads=4`**.  It drew 14.74 cores, finished the canary in 9.32 minutes
(9.0% / 55.34 seconds faster than W23), used 10.305 core-s/column, and peaked at
34.40 GiB.  On a marketed 64 GB machine (about 59.6 GiB) that leaves about 25.2
GiB for the OS, browser, and other processes.  Use W32/GC4 as the headroom
fallback: 26.45 GiB peak, only 13.17 seconds / 2.3% slower than W40.

There is therefore no truthful knob-only answer that “makes it draw close to
all 21.”  W40/GC4 is the best measured wall point; getting substantially closer
requires reducing allocation/GC in weight assembly and making the `runchunks`
caller thread participate in work.

## Campaign definition and guards

- Worktree/branch: `/home/asinghvi17/geo/DGG-perf-ladder`,
  `claude/perf-ladder`, starting at `71de327`.
- Harness extension: `scripts/southpole_shakedown.jl`, commit `aa917ca`.
  Julia 1.12 exposes no `Base.SIGUSR1` constant; the Linux profile-signal fix is
  commit `366f6f1`.  Both commits were pushed append-only before accepted runs.
- Julia 1.12.6; `nice -n 10`; outer shape; exactly 21 Julia threads.  Every
  child reported the requested first-field GC count and `gcsweep=0`.  No launch
  contained a second `--gcthreads` field.
- Synthetic source only, `real=:none`, `prefetch=0`, no network provider and no
  downloads.
- Runs were sequential: baseline, W32/GC4, W40/GC4, W40/GC8, W32/GC8.  The
  supervisor sampled the child externally every three seconds and imposed a
  40-minute wall limit plus a 40 GiB kernel-HWM limit.
- Fresh stores were `/home/asinghvi17/geo/scratch-stores/t21util-<config>.zarr`.
  Logs have the same prefix and a `.log` suffix.  The W40/GC8 partial store is
  retained; it was not used as a completed measurement or correctness oracle.
- Raw records are local-only at
  `/home/asinghvi17/geo/DiscreteGlobalGrids.jl/regrid-notes/2026-08-22-t21-utilization.ndjson`.
  The NDJSON is not committed.

The first baseline attempt stopped at the requested profile signal because of
the harness-only `Base.SIGUSR1` name error.  Its 193-column store and log are
recoverably preserved under
`t21util-baseline-failed-sigusr1-20260822*`.  It has no `run_final`, is excluded
from every result below, and motivated `366f6f1`.  The fresh rerun completed.

## Exact 800-column set

The selector recovers and revalidates the canonical 2,000-column set and its
1,087-column `< -80 degrees`, production-ledger `nan==0` full-land block.  It
then keeps:

1. both canonical southernmost IGeo7 pentagons, 112049 and 154067; and
2. 798 evenly spaced ranks, including both endpoints, after sorting the 1,087
   deep full-land columns by `(latitude, column)`.

The result has exactly 800 unique columns and SHA-256:

```text
f0dd12b7f7fb35be7fce9a7d9718821a9aeb2bd536706abeb63f6a11d2fe0313
```

All completed runs wrote the same 658,559,886 Float32 positions.  The two
pentagons are intentional all-fill topology checks; the other 798 columns keep
the set expensive and land-heavy.

## W resolution and overrides

The production policy does not use `Threads.nthreads()` when resolving W:

```text
outer: W = ceil(config.cores / 1.06)
inner: W = ceil(config.cores / 3.1)
```

The retained fixed `cores=24` budget therefore resolves outer to W=23 at
`-t21`.  `config.workers` is the supported positive override used for W32 and
W40.  Batch 8, guided taper, affinity order, refcount cache, 1 GiB lazy budget
per worker, and the other validated driver settings were unchanged.

Simply replacing `config.cores` with `Threads.nthreads()` in the current
formula would resolve W=20, which is opposite the measured direction.  A later
default-policy change needs a retuned oversubscription factor and a memory cap,
not a blind variable substitution.

## Baseline attribution

### Per-thread busy distribution

The supervisor, running outside the t21 child, read
`/proc/<pid>/task/<tid>/stat` every three seconds.  Before compute the child
mapped Julia default-pool thread IDs 1--21 to unique Linux TIDs.  This avoids
confusing the Julia pool with the additional GC, runtime, I/O, and BLAS TIDs in
`/proc`.

Over the steady 61--605 second window:

| population | aggregate busy fraction |
|---|---:|
| Julia thread 1 (driver caller) | **0.20%** |
| threads 2--21 minimum | 63.77% |
| threads 2--21 p10 / median / p90 | **63.89 / 64.28 / 64.83%** |
| threads 2--21 maximum | 64.96% |

Across individual three-second samples, the median p10/median/p90 thread busy
fractions were 61.92/62.58/67.22%.  The low-busy count was one, not seven: it
was always thread 1.  All 20 worker threads rose and fell together.  This rules
out a stable `~14 busy + ~7 idle` runnable-task shortage.

The same externally sampled steady window averaged 14.09 process cores:

| accounting component | mean cores |
|---|---:|
| threads 2--21 | 12.86 |
| Julia thread 1 | 0.002 |
| four post-pool GC TIDs plus other runtime work | 1.22 |
| process total | **14.09** |

Relative to 21 cores, the accounting is:

```text
parked caller thread                         1.00 core lost
20 workers: 20.00 - 12.86                   7.14 cores lost during GC
GC/runtime TIDs doing the collection       -1.22 cores recovered as GC work
                                            ----
steady-state gap                            6.92 cores (measured 6.91)
```

The whole-run 13.64-core mean includes graph/store setup and queue drain, adding
about 0.44 core to the gap.  That is where the headline 7.36 missing cores go.

### GC share

`Base.gc_num()` was differentiated on every in-process heartbeat and over the
run endpoints.  The endpoint values are not affected by heartbeat delay.

| baseline GC metric | result |
|---|---:|
| GC time / wall | **221.29 / 614.76 s = 36.0%** |
| total mark time | 79.97 s (13.0% of wall) |
| total sweep time | 141.31 s (23.0% of wall) |
| time reaching safepoints | 4.19 s |
| pauses / full sweeps | 2,807 / 451 |
| allocated bytes during measured run | 2,089.76 GiB |

Internal GC intervals had p10/median/p90 wall fractions of
32.25/37.10/41.38%.  The internal sampler was scheduler-delayed (interval
median/p90/max 6.43/13.20/26.11 seconds rather than the requested five), but
every sample uses its actual interval and the independent external sampler
stayed at three seconds.

The match is direct: worker-thread median busy is 64.28%, while
`1 - whole-run GC fraction` is 64.00%.  The four GC TIDs consumed about 1.22
cores while the 20 worker TIDs were stopped, so the net loss is about 5.92
cores rather than the full `20 * 0.36 = 7.20`.

### Profile and lock suspects

At 182 seconds the supervisor requested one 60-second Julia SIGUSR1 profile at
10 ms sampling delay.  The combined flat report contained 95,425 snapshots.
Its largest inclusive stacks were:

| frame family | inclusive snapshots |
|---|---:|
| lazy `blockfor` / `_fillwave!` | 57--58 k |
| `build_weights!` | 54.5 k |
| conservative `_intersectionareas` / `intersection_areas` | 49.9 k |
| candidate-pair / depth-first intersection walk | 25.5 k |
| COO chunk assembly | 23.8 k |

`poptask` had 3,817/95,425 snapshots (4.0%), all on the caller waiting in the
`runchunks` `Threads.@sync`.  The profile's printed 96% utilization is
conditional on stacks the CPU sampler could collect; stopped worker TIDs are
accounted by `/proc` and GC wall time instead.  The signal stack also caught GC
helper threads in `jl_parallel_gc_threadfun` and worker threads in
`jl_safepoint_wait_gc`.

No `DGG.dggwrite!`, `DoneLog.record!`, or `claim!` frame appeared in the flat
profile.  The refcount-cache retire lock had only 53/95,425 inclusive samples
(0.056%); `gettile!` had 24.  No lock is remotely large enough to explain a
seven-core gap.  Zarr serialization, done-log append, the pull cursor, and the
cache lock are therefore ruled out as the cap.

## Sweep results

Wall includes setup; compute is earliest ledger start through final completion.
Core-s/column and mean cores use process CPU endpoint deltas.  The p10/median/p90
cores are differentiated in-process samples using their real intervals.  Peak
RSS is kernel `VmHWM`, not the lower sampled current RSS.

| config | W / GC mark | wall / compute | cells/s | core-s/col | mean cores; p10/med/p90 | GC wall | peak / end RSS | last-10 tail | outcome |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| baseline | 23 resolved / 4 | 10.25 / 9.97 min | 1,100,571 | 10.484 | 13.64; 12.95/13.94/14.69 | 36.0% | 20.18 / 7.86 GiB | 20 s | complete |
| W32-gc4 | 32 override / 4 | 9.54 / 9.27 min | 1,184,056 | **10.177** | 14.22; 13.07/14.54/15.25 | 32.0% | 26.45 / 8.78 GiB | 27 s | complete |
| W40-gc4 | 40 override / 4 | **9.32 / 9.05 min** | **1,212,170** | 10.305 | 14.74; 13.36/15.08/16.46 | **29.5%** | 34.40 / 12.68 GiB | **16 s** | complete |
| W40-gc8 | 40 override / 8 | >4.66 min | n/a | n/a | n/a | n/a | **40.19 GiB HWM** | n/a | guard abort at 358/800 |
| W32-gc8 | 32 override / 8 | 9.39 / 9.12 min | 1,203,113 | 10.424 | **14.80**; 13.90/15.16/16.03 | 31.0% | 33.40 / 8.25 GiB | 22 s | complete |

All four completed runs demanded and loaded the same 3,585 source tiles exactly
once.  Reloads, uncredited demands, live tiles at end, and pinned tiles at end
were zero.  Cache peak was 0.777/1,046 tiles at W23, 0.913/1,233 at W32/GC4,
1.087/1,414 at W40/GC4, and 0.956/1,294 at W32/GC8.  Joined loads were zero,
zero, one, and one respectively.

### W effect and memory price

Relative to baseline:

- W32/GC4: 6.86% less wall, 2.93% less CPU per column, +0.58 mean cores,
  -4.0 percentage points GC wall, and +6.27 GiB peak RSS.
- W40/GC4: 9.00% less wall, 1.70% less CPU per column, +1.09 mean cores,
  -6.46 percentage points GC wall, and +14.21 GiB peak RSS.

The observed RSS price was 0.696 GiB per added worker from W23 to W32 and
0.993 GiB per added worker from W32 to W40 (0.836 GiB per worker over W23 to
W40).  The marginal W32-to-W40 gain is only 2.30% / 13.17 seconds wall for
7.95 GiB, so W40 is the measured wall knee rather than evidence to keep raising
W.  Extrapolating the latest slope puts only about five more workers below the
40 GiB guard, with no basis to expect a large utilization jump.

### First-field GC effect

At W32, GC8 versus GC4 changed:

- wall by -1.56% (8.94 seconds);
- mean cores by +4.06% (14.22 to 14.80);
- CPU per column by **+2.43%**;
- GC wall only from 32.01% to 30.97%; and
- peak RSS by **+6.95 GiB** (26.45 to 33.40 GiB).

This is more CPU draw, not a compelling throughput win.  W40/GC8 exceeded the
40 GiB limit before halfway, so it cannot be the deployment answer.  A GC16
variant was not run: the baseline plus four allowed variants were exhausted,
GC8 already had the wrong CPU/RSS trade, and W40's remaining GC was dominated
by 120.28 seconds of sweep versus 44.96 seconds of mark.  More first-field mark
threads cannot safely promise a near-21-core result.

## Deployment recommendation and 64 GB frame

Use:

```text
nice -n 10 julia --project=benchmark -t 21 --gcthreads=4 ...
shape=:outer
workers=40                 # explicit override; do not leave workers=0 / W=23
batch=8, taper=true, schedule=:affinity, cachepolicy=:refcount
```

Never add a nonzero second GC field.

Expected from this exact synthetic canary: **14.74 whole-run cores** (15.25 in
the externally sampled steady window), **9.32 minutes wall**, **10.305
core-s/column**, and **34.40 GiB peak / 12.68 GiB shutdown RSS**.  Against the
matched W23 baseline that is 55.34 seconds / 9.0% faster.  Applying only that
measured ratio to the earlier 22.88-minute 2,000-column result gives a rough
20.82-minute estimate, not a measured 2,000-column W40 time.

Memory headroom at the measured peak:

| frame | headroom after W40/GC4 peak |
|---|---:|
| 64 GiB convention | 29.60 GiB |
| marketed 64 GB = about 59.60 GiB | **25.20 GiB** |

That is adequate for an ordinary OS and browser, but not unlimited.  If their
combined working set can approach 20 GiB or a real-source canary adds more than
about 5 GiB, choose W32/GC4.  It leaves about 33.15 GiB on a 59.6 GiB machine,
saves 7.95 GiB versus W40, and costs only 13.17 seconds / 2.3% on this canary.
Do not choose GC8 merely to make CPU monitoring show a larger number.

## What needs code change

Config alone recovered only 1.15 of the 7.36 missing whole-run cores.  Two code
seams remain:

1. **GC/allocation:** the dominant seam is transient allocation in
   `GlobalRegridding.build_weights!` and the
   `ConservativeRegridding._assemble_chunk` / `_run_and_store!` COO and
   intersection-area path named by the profile.  Reusing per-worker COO,
   candidate-pair, and weight buffers or otherwise reducing their allocation is
   preferable to increasing GC threads.  Halving W40's measured 29.5% GC wall
   is worth roughly 15% wall and should move draw into the 17--18 core range.
   Removing nearly all of it is an upper bound near 30% wall and about 20 cores.
2. **Parked caller:** `scripts/copdem_production.jl:runchunks` spawns all W
   workers and waits in `Threads.@sync`; Linux TID data show that caller never
   participates.  Running one worker loop inline while spawning the others can
   recover at most about 0.7 mean core at the current GC fraction (roughly
   4--5% wall after GC is improved).

Together those are the path toward 21.  Parallelizing a per-column assembly
tail is not the first seam here: 20 worker TIDs are uniformly busy outside GC,
and W40 already hides serial column phases.  Zarr/done-log/cursor locking is
also not a candidate.  The forbidden concurrent GC sweeper might appear to
target the 120-second sweep component, but its known page-corruption failure
makes it categorically unavailable.

The `workercount` default should eventually be retuned so the deployment does
not silently resolve from a stale fixed 24-core budget.  The new policy must
combine thread count, an empirically measured outer oversubscription factor,
and a memory ceiling; this campaign supports W40 specifically at t21, not a
universal `2*nthreads` rule.

## Implication for this box's real GLO-90 run

At report time this host exposed 64 CPUs, 751 GiB total RAM, and 204 GiB
available RAM (swap was full).  The measured 34.40 GiB W40/GC4 peak fits its
current frame comfortably.  For the imminent real run on this box, use
outer/`-t21`/W40/single-field GC4 and keep the affinity/refcount/taper stack.
Do not raise W beyond 40 without another measurement; W40 is already the wall
knee and the deployment-memory slope is steep.

This campaign had no network and does not supersede the accepted real-source
policy.  Retain the separately tested prefetch depth 32 / fetch concurrency 16,
watch `demand-cold downloads` (expected zero) and kernel HWM early, and fall
back to W32/GC4 if the real provider adds unexpected memory.  The real run must
still omit the unsafe GC sweep field.

## Bit identity

Baseline was compared with all three completed variants on six deterministic
columns:

```text
112049, 121338, 122678, 123072, 123203, 154067
```

The set includes both pentagons and four evenly spaced deep full-land columns.
Across three candidates and 14,823,774 compared Float32 positions there were
**zero bit differences**.  Configuration tuning did not change values.

