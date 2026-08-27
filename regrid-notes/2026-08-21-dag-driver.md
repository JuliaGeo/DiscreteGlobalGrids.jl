# DAG-driven scheduling in the CopDEM production driver

Branch `claude/dag-driver` (worktree `/home/asinghvi17/geo/DGG-dag-driver`), off
`claude/perf-ladder` @ `0d0f069`, with `claude/chunk-dag` and
`claude/driver-reporting` merged in. Measurements in
`2026-08-21-dag-driver.ndjson`.

Companions: `2026-08-21-chunk-dag-api.md` (the graph this consumes) and
`2026-08-21-chunk-dag-sim.md` (the offline simulator whose policy table this
implements). Read the second one first; everything below is either its
prediction realised or its prediction corrected.

## 1. What landed, and where

| piece | file | why there |
|---|---|---|
| the graph query | `lib/GlobalRegridding/src/chunkgraph.jl` | already merged; untouched |
| Morton key, affinity order | `scripts/copdem_policy.jl` | policy |
| `GuidedSchedule` + batch taper | `scripts/copdem_policy.jl` | policy |
| `RefCountCache`, `StripedLRUCache` | `scripts/copdem_policy.jl` | policy |
| `Prefetcher` | `scripts/copdem_policy.jl` | policy |
| `dagplan`, the `refine` narrow phase, `reportcache` | `scripts/copdem_production.jl` | workload-specific |
| policy tests (74) | `test/scripts/copdem_policy.jl` | wired into the root suite |
| byte-identity A/B | `scripts/copdem_dag_validate.jl` | standing proof |
| residency / tail / prefetch measurement | `scripts/copdem_dag_scale.jl` | reproducible |

Nothing new went into `GlobalRegridding`. The settled principle — graph query in
the library, policy in the driver — survived contact: the policy module takes the
adjacency as two functions, `sourcesof(d)` and `consumerdegree(s)`, so it never
mentions `ChunkDependencyGraph`, and 73 of the 74 tests run on hand-written
`Vector{Vector{Int}}` adjacencies in milliseconds.

## 2. Codex advice, and where I deviated

One consultation, `gpt-5.6-sol` at `model_reasoning_effort=xhigh`, on ownership,
the cache's locking discipline, the prefetcher's task structure, the test seam,
and the lock-free taper.

**Adopted, essentially verbatim:** all five features are driver policy and none
is a regridding concept; build the destination space in canonical chunk order and
keep the Morton walk as a separate permutation `order[p] -> d`; one lock over the
whole cache rather than stripes, held only across array reads and writes;
single-flight loading through a level-triggered `Base.Event`; a shared
`Base.Semaphore` bounding demand and speculative fetches together; retire a chunk
on every terminal outcome, in a `finally`, *before* the Zarr write; retire
resume-skipped chunks up front and drop them from the order rather than testing
`ch in done` in the hot loop; discard a speculative result whose refcount reached
zero mid-flight; a fixed prefetch pool plus a deduplicated request channel rather
than a task per tile, with a `frontier` the driver never rewinds; prefetch errors
captured into the cache and re-raised only at the next demand; `stop_prefetch!`
that never rethrows and so cannot mask the real exception; the CAS taper loop;
`scripts/copdem_policy.jl` as a side-effect-free module tested from
`test/scripts/`, with the driver made includable so a harness can call
`main(config)`.

**Deviated:**

1. **The "loads == tiles" assertion the brief asked for is not a law, and Codex
   was right to say so.** Refcount guarantees *at most one* load per source
   chunk. It cannot guarantee at least one: the graph is a conservative superset,
   a resumed run retires chunks before they load anything, and a failed chunk
   never reads what it would have. So the hard check is
   `all(attempts .<= 1)` — enforced at the load site, not just at the end — and
   `loads == nsrc` is reported as a coverage diagnostic on a fresh, complete,
   failure-free run only. The simulator could equate the two because it treated
   adjacency as demand.
2. **An uncredited demand is served, not refused.** Codex specified
   `@assert remaining[s] > 0` on the demand path. A wrong assertion there
   destroys twenty hours of production work to report an accounting bug, so the
   cache loads the chunk transiently — outside the cache, so refcounts and the
   at-most-once count stay clean — increments `uncredited`, and the run's closing
   `check` fails on it. Same signal, no lost run. This turned out to matter
   within the hour (§3).
3. **The order is the simulator's tile-sweep affinity rule, not "Morton over
   destination-column centroids"** as the brief put it. The sweep is over the
   *tiles*; a chunk is emitted when its last tile is swept, ties broken by its
   first. That is the rule that measured 1.43 GiB, and ordering columns directly
   by their own Morton code was never measured.
4. **No `refinegraph` by default** — §3, the one real finding here.

## 3. The `refine` hook is the wrong tool for eviction

The plan said to build the driver's graph with the lon/lat box narrow phase,
which takes the cap relation's 1.75x edge inflation down to 1.35x. I did, and the
first A/B run reported **8 uncredited demands over 11 chunks**.

Isolated per chunk, with a one-chunk destination space each:

| chunk | refined edges | tiles the regrid demanded | uncredited |
|---|---|---|---|
| 100 | 1 | 3 | 5901, 6130 |
| 72197 | 2 | 4 | 13357, 13359 |
| 1 | 4 | 5 | 5439 |
| 32646 | 4 | 5 | 7057 |
| 57594 | 2 | 3 | 16855 |
| 136095 | 2 | 3 | 18694 |
| 123203 (360 tiles) | 360 | 360 | — |

Every one of those tiles **is** in the cap-only graph. So the narrow phase is not
buggy: it is correct about the geometry and wrong about the executor.

The reason is that refcount eviction does not need a superset of the true
geometry. It needs a superset of what the regridder actually **reads**, and the
regridder pairs its own chunks with exactly the cap test the broad phase uses.
Every pair the narrow phase removes is a pair the regridder still reads the tile
for. The data was never at risk — byte-identity passed even with 8 uncredited
demands, because the transient path serves them — but the tile had already been
freed, so the reload the whole design exists to prevent came straight back.

`refinegraph` therefore defaults to `false`. The narrow phase stays in the driver
behind the knob, for the day the regridder's own pairing gets tighter, with the
finding written where the knob is. The cost of not using it is §4's 2.12 GiB
instead of the simulator's 1.43 GiB.

Generalised: **a `refine` that is tighter than the consumer's own pairing rule is
unsafe for scheduling, however correct it is about geometry.** The
`chunk_dependency_graph` docstring warns that a wrong `refine` silently corrupts
results; this is the scheduling-side version of the same hazard, and it is
quieter, because it costs performance rather than correctness.

## 4. The globe, replayed

The real 66 178 x 26 475 graph, the affinity order, refcount eviction, 23 chunks
in flight, tile sizes taken from the source grid's own pixel counts. Pure
arithmetic — no regrid, no pixels, 0.03 s.

| order | eviction | peak | peak tiles | loads |
|---|---|---|---|---|
| affinity | refcount | **2.12 GiB** | 1 162 | 26 475 = every tile, once |
| canonical | refcount | 4.29 GiB | 1 621 | 26 475 |
| any | LRU 3072 (the run's setting) | 16.48 GiB | 3 072 | — |
| any | LRU 1024 | 5.49 GiB | 1 024 | — |

All 26 475 tiles resident at once would be 93.62 GiB.

Against the simulator: it measured 1.43 GiB for affinity+refcount and 3.29 GiB
for baseline+refcount using the *exact* `MultiOrderCoverage` adjacency. This graph
carries 326 392 edges against that adjacency's 186 069, a 1.75x superset, and the
peaks scale almost exactly with it — 1.48x and 1.30x. The prediction held; the
extra is the conservatism §3 says we have to keep.

**The 16.48 GiB LRU bound is reproduced to the digit**, which is the check that
the two calculations are measuring the same thing.

So: **7.8x less tile residency than the run's current LRU, and zero reloads by
construction.** The tile cache is now 2 GiB of a 64 GiB budget.

## 5. Live, 400 chunks

Region 0-20E / 35-50N, 279 tiles, 400 chunks, W = 12, `nice -n 10`,
`--gcthreads=4`, on a box shared with another agent.

| side | schedule | cache | taper | wall | tail idle | loads | reloads | peak |
|---|---|---|---|---|---|---|---|---|
| legacy | canonical | LRU 3072 | no | 384.3 s | **67.6 s** | 279 | 0 | 1.46 GiB held at end |
| dag | affinity | refcount | yes | **354.6 s** | **10.8 s** | 279 | 0 | **0.38 GiB / 73 tiles** |

- **Residency 3.8x lower**, and the working set is 73 tiles of the region's 279 —
  the number the LRU was paying 3072 slots to not know.
- **Zero reloads on both sides**, and zero uncredited demands. At this scale the
  LRU's 3072 slots exceed the 279 tiles that exist, so it cannot reload either;
  the reload comparison only bites at globe scale, where §4 has it.
- **The tail: 67.6 s -> 10.8 s.** This is the taper doing exactly what the
  simulator said it would. 400 chunks is 1/165 of the real run, and the effect
  there was 301 s -> 38 s.
- Makespan −7.7 %, which on a shared box is noise-adjacent; the simulator's
  prediction is ±0.02 %, and nothing here contradicts that.

## 6. The prefetcher

It cannot help a source that costs nothing to obtain, and a fabricated tile costs
nothing, so latency was imposed with `fetchdelay` — a `sleep` in the tile
builder, config-gated, defaulting to zero.

**At 0.5 s per tile** (150 chunks, ~3.7 tiles per chunk) the prefetcher made no
difference at all: 136.5 s without, 136.4 s with. ~2 s of fetch per chunk against
~10 s of compute, spread over 12 workers, is already hidden by demand loading.
Residency, however, rose 0.32 -> 1.02 GiB: a depth of 32 over a 150-chunk queue
is a fifth of the whole run in flight. **Prefetch depth is not free when the
queue is short.**

**At 3.0 s per tile**, `fetchconc = 8`, 100 chunks, both sides on the same
permit count:

| depth | wall | peak | peak tiles |
|---|---|---|---|
| 0 | 166.0 s | 0.23 GiB | 43 |
| 16 | **105.5 s** | 0.34 GiB | 64 |
| 48 | 104.1 s | 0.36 GiB | 69 |

**−36 % makespan for +0.11 GiB**, and depth 48 buys nothing over 16 — the knee is
where the simulator put it (`N ~ 32L` at `fc = 8`; here fetch is only ~11 s per
chunk against ~10 s of compute, so a shallower depth suffices). The mechanism
works and the depth knob behaves. `prefetch` stays 0 for the synthetic run and
should be 32 with `fetchconc = 16` for the AWS one, per the simulator's §6.

## 7. Byte identity

`scripts/copdem_dag_validate.jl`, run twice into scratch stores: 11 chunks — the
two worst polar ones at 360 tiles each, a single-tile one, three all-ocean ones,
and an even sample between — under `(canonical, LRU, no taper, no prefetch)` and
under `(affinity, refcount, taper, prefetch 16, fetchconc 4)`.

**8 chunk files, 11.98 MiB, byte-identical.** The three missing files are the
all-`NaN` chunks, which Zarr stores as nothing at all — themselves a case worth
having in the set, since they exercise the path where a chunk consumes tiles and
writes no file.

Both sides ran with zero failures, zero uncredited demands, and the DAG side
loaded each of its 389 tiles exactly once.

This is not a comparison against the production store: the closed-form Snyder
inverse landed on this base, so the numbers moved for unrelated and intended
reasons. Only the within-branch A/B is meaningful.

## 8. Tests

- `test/scripts/copdem_policy.jl`, **74 tests**, wired into `test/runtests.jl`.
  Concurrent claims covering `1:n` exactly once under a varying batch; the taper
  reaching size 1; twenty concurrent demands invoking one loader; a value
  surviving the first of two retirements and not the second; resume retirement
  making a source dead; a failed regrid still retiring; a load that throws being
  remembered rather than retried; a speculative failure deferred to the next
  demand; a refcount reaching zero mid-flight discarding the publication; an
  early `stop_prefetch!` leaving no live task; the LRU escape hatch; and one test
  that plugs a real `chunk_dependency_graph` into the same seam.
- `lib/GlobalRegridding`: **2 295 pass, 1 broken**, unchanged by this branch.
- Driver dry checks extended: the graph's bipartition matches the work list, the
  order is a permutation of it, and the dry run reports both degree extremes.
  The closing report checks at-most-once, no uncredited demand, and full
  retirement.

## 9. Merges

- `claude/chunk-dag` @ `0419988` — clean, no conflicts.
- `claude/driver-reporting` @ `aeceb55` — two conflicts in `main`, both resolved
  by keeping both sides: the dry-run block now prints the graph's degree extremes
  *and* runs their `reportchecks`; the closing banner keeps their session/store
  split, drops the single `destination NaN fraction` line it supersedes, and
  gains the tile-cache report. Their loader counters moved from `TiledDEM` to
  `TileBuilder` on this branch and were repointed.

The bench environment was resolved from `DGG-perf-ladder`'s known-good
`Manifest.toml`; only `Graphs` and its three dependencies were added, and
ConservativeRegridding stayed pinned at `6a4b997`.

## 10. Not done

- No run at globe scale. §4 is a static replay, not a regrid; the residency
  number is arithmetic over the real graph and the real tile sizes, and the tail
  and reload behaviour are measured only at 400 chunks.
- The prefetcher has never met a real network. `fetchdelay` is a `sleep`, so it
  models latency and nothing else — no bandwidth term, no failures, no retries,
  no S3 client. The retry policy in particular is deliberately absent: a failed
  load is remembered for the process and retried only by resuming.
- `cachepolicy = :lru` refuses `prefetch > 0` rather than supporting it. A
  bounded LRU can evict a prefetched tile before it is used and that churn was
  never measured.
- The 1.75x conservatism of §3 is real cost — 2.12 GiB where 1.43 would do.
  Closing it needs the *regridder's* pairing to get tighter, not the graph's.
