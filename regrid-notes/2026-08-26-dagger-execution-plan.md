# Dagger as an isolated regridding execution backend

Status: **isolated D0/D1 scaffold implemented on 2026-08-26 and reorganized on
2026-08-27**. The regular regridding entry points and behavior remain unchanged;
their reusable CopDEM definitions now live in one shared include. A one-chunk,
one-process synthetic vertical slice passes; no production-data canary,
byte-identity comparison, or speed claim has been made yet.

Companions: `2026-08-21-chunk-dag-api.md`, `2026-08-21-chunk-dag-sim.md`,
`2026-08-21-dag-driver.md`, `2026-08-22-t21-utilization.md`, and
`2026-08-23-full-run-attribution.md`.

## 1. Decision and boundary

Test Dagger as a **separate, optional execution backend**, reached through
`scripts/dagger_regrid/main.jl`. Do not add Dagger calls, conditionals, options,
or dependencies to `GlobalRegridding.regrid!`, `LazyRegridArray`,
`_readdestination!`, or the regular CopDEM production loop while the idea is
being evaluated.

The domain scheduler remains ours. It owns the source/destination graph,
destination priority, source affinity, batching, tapering, cache policy, retry
policy, and completion ledger. Dagger receives only bounded batches that our
scheduler has declared ready and provides process placement, execution, failure
delivery, and instrumentation.

This is only approximately "Dagger as an execution engine": an eager Dagger
`DTask` still passes through Dagger's own scheduler. The separation is that the
Dagger scheduler never decides the domain-level order or receives the complete
regridding graph as a task DAG.

The execution unit is a **complete destination chunk or a small,
affinity-contiguous batch of destination chunks**. It is not a
`(destination, source)` graph edge and not a partial floating-point reduction.
The existing executor may construct weights concurrently, but applies source
contributions in ascending source-chunk order. Keeping the complete destination
fold inside one worker preserves the current byte-identity contract.

## 2. Isolation in the tree

The implementation is a small, explicit script directory:

```text
scripts/dagger_regrid/main.jl            readable top-level smoke/canary/run flow
scripts/dagger_regrid/README.md           short purpose and usage guide
scripts/dagger_regrid/src/DaggerRegrid.jl package module and include order
scripts/dagger_regrid/src/types.jl        configuration and transport values
scripts/dagger_regrid/src/workers.jl      process-local source/cache/write work
scripts/dagger_regrid/src/coordinator.jl  graph preparation, admission, and ledger
scripts/dagger_regrid/src/smoke.jl        process-placement smoke
scripts/dagger_regrid/copdem_helpers.jl  shared CopDEM source/graph/chunk helpers
scripts/dagger_regrid/Project.toml        isolated Dagger dependency environment
```

The project file is dependency metadata rather than another implementation
layer. It avoids adding Dagger to the root environment or coupling the
experiment to the presentation environment that happens to contain Dagger
today. Its package metadata lets the entry point and distributed workers load
the implementation with `using DaggerRegrid`; no absolute source path or
worker-side `include` protocol is needed. The executable `main.jl` contains
top-level control flow, not a CLI function or `let` block, and leaves its
intermediate values bound for interactive inspection.

Both drivers eagerly include `copdem_helpers.jl`. The file holds the common
configuration, source construction, graph planning, one-chunk regrid, store,
and scheduler/cache definitions. `copdem_production.jl` retains only its
threaded run orchestration, profiling, progress, checks, and reporting. This
keeps one implementation of the domain work while ensuring every Dagger worker
has the complete method table before its first task. `smoke` opens no source,
graph, or store, but it pays the geospatial package-load cost. The Dagger
backend reuses three boundaries:

1. `dagplan` and the graph queries for destination order and source affinity;
2. `regrid_chunk` as the complete, deterministically ordered work unit;
3. the disjoint chunk writer, with ledger commits kept on the coordinator.

Worker processes load the `DaggerRegrid` package and reconstruct their own
provider, source space, cache state, and output handle from plain configuration.
The top-level script shows the coordinator lifecycle directly: preparation,
Shard construction, bounded admission, completion, teardown, and final
statistics. Only the mechanical operations live in helper functions.

The top-level script exposes `smoke`, `canary`, and `run`, so the first
serialization checks and fixed-core comparison use the exact implementation
being tested rather than parallel test/harness scaffolding. Dedicated test or
benchmark files come later only if the backend survives the promotion gates.

These ownership boundaries remain:

- `scripts/dagger_regrid/main.jl` does not become another mode inside
  `scripts/copdem_production.jl`.
- Dagger is not added to the root `Project.toml` during the experiment. The
  separate script environment tracks Dagger's upstream `master` branch while
  the experiment is active.
- Production scheduler policy may be called or included from
  `scripts/copdem_policy.jl`; it is not copied into GlobalRegridding.
- Smoke checks and canaries invoke `scripts/dagger_regrid/main.jl` explicitly.
  Nothing selects it automatically from `regrid!`.
- No extension module is added until the experiment passes its promotion gates.

`DaggerRegridConfig` is the narrow boundary between the readable script and its
workers. It contains serializable descriptions of the input, output,
destination range, method, missing policy, scheduling limits, and worker
topology. It must not contain an open store, downloader, `LazyRegridArray`, live
`ChunkedPlan`, lock, condition, GDAL dataset, or PROJ pointer. This is an
experimental script contract, not a new public package API.

### Implemented slice

The files under `scripts/dagger_regrid/` now provide:

- `DaggerRegridConfig` and a compact `DaggerRegridReport` around the top-level
  execution flow;
- a coordinator-owned affinity order and guided rolling admission window;
- process-scoped DTasks, with one completion-channel watcher per admitted batch;
- one process-local provider, source space, bounded LRU, and output handle held
  by a `Dagger.Shard` rather than serialized per task;
- complete-destination computation and disjoint worker-side Zarr writes;
- coordinator-only done-ledger appends after a successful compact report;
- process-private real-source download cache directories, avoiding cross-process
  `.part` races at the deliberate cost of measuring duplicate downloads;
- `:before_compute`, `:before_write`, and `:after_write` failure injection;
- safe exceptional teardown that drains side-effecting write tasks instead of
  treating Dagger cancellation as process preemption; and
- `smoke`, `canary`, and `run` command modes.

The implementation assigns a **process**, not an exact thread. Dagger still
chooses a thread within that process. With the first intended topology of one
Julia thread per worker process, process and exact-thread placement coincide.
This slice reports failed chunks but does not retry them in the same session;
bounded retry/admission policy remains D4 work, with the existing store/ledger
union providing safe resume behavior meanwhile.

The isolated project selects Dagger directly from its upstream `master` branch.
The validation below resolved `master` at commit
`c6ff43da9c02cad11e9940950f5bc6ffdbbb41f8` (`Dagger v0.22.1`). The repository
ignores generated `Manifest.toml` files, so `master` is deliberately a moving
target: record the resolved commit with every performance result. The project
also follows the current ConservativeRegridding `main` selection and resolves
DGG and GlobalRegridding from this checkout.

### Validation completed

- Julia parsing and `git diff --check` pass.
- The isolated project resolves and loads Dagger 0.22.1 from upstream `master`.
- Executor-only smoke passes on the coordinator process.
- Executor-only smoke passes on two local single-thread worker processes and
  verifies both that each process-scoped task ran on its requested process and
  that it received the `Dagger.Shard` state piece owned by that process.
- A one-worker synthetic nearest-direct canary constructs worker-local state,
  builds one CopDEM tile, computes all seven finite cells of one level-5 to
  level-6 destination chunk, writes the disjoint temporary Zarr chunk, commits
  its coordinator ledger entry, and returns compact worker/cache statistics.
- The first canary exposed a Julia 1.12 world-age boundary when production
  helpers were loaded from inside a Dagger task. The shared helper is now loaded
  while the coordinator and worker modules are defined, removing the lazy-load
  lock, worker branch, and all `invokelatest` calls.

The next action is still the D1 8--32-chunk byte-identity canary below, not a
production run.

## 3. Ownership and data flow

```text
coordinator process
  graph + affinity order + guided cursor + reservations + completion ledger
                       |
              bounded ready batches
                       v
              process-scoped DTasks
                       |
                       v
worker-local source / plan / cache / output handle
                       |
        ordered complete destination computation
                       |
             disjoint Zarr chunk write
                       v
          small BatchReport to coordinator
```

### Coordinator state

The coordinator owns one graph and the canonical mapping between graph rows and
destination chunk IDs. It keeps at most a configured number of batches in flight
per worker process, selects the target process, and spawns an eager DTask with a
process scope. It submits a rolling window rather than all 66,175 destination
chunks at once.

The coordinator is also the single owner of store metadata and the completion
ledger. A worker may write its disjoint physical Zarr chunk directly, but only a
successful `BatchReport` advances the ledger. Reports contain destination IDs,
timings, counters, output checksums or byte counts, and structured failures; they
do not contain destination arrays.

### Worker state

Each process reconstructs its state once from plain configuration, preferably as
a Dagger `Shard`. It reopens the source provider and output store, reconstructs
spaces and plans, and owns process-local source/weight caches. A job sends only a
small vector of destination IDs and immutable run options.

The worker explicitly establishes the nested-parallelism mode. It must not rely
on a process-local scoped value such as `OUTER_PARALLEL` propagating over a
distributed task boundary.

Source download coordination also becomes process-aware. The current CopDEM
single-flight locks protect threads in one process; they do not prevent two
processes from fetching or renaming the same `.part` file.

## 4. What the chunk graph does, and does not do

`ChunkDependencyGraph` is a compact scheduling oracle. `sourcesof(d)` supplies
the affinity signature of a destination; `consumerdegree(s)` supplies the global
number of possible consumers. The graph relation is deliberately conservative:
an edge means that execution may read a source chunk, not that the pair has a
positive final weight.

Consequences for this backend:

1. Do not create one DTask per edge. The production graph has about 326,000
   edges, while the useful independent outputs number about 66,000.
2. Do not ask Dagger to reduce source-edge results. A different floating-point
   reduction tree would leave the existing byte-identity contract.
3. Do not submit the whole graph and expect Dagger locality to rediscover the
   existing Morton/affinity policy. Our scheduler must make placement explicit.
4. A destination can be retried as a complete unit because no destination writes
   another destination's output.

## 5. Cache policy across processes

The current process-wide refcount cache can guarantee at most one load of a
source chunk in a run. That guarantee does not automatically survive multiple
worker processes. A global consumer count says when all destinations are done;
it does not say when a particular process has finished with its local copy.

The first slice therefore uses a bounded worker-local LRU, or no persistent
source cache, and measures duplicate source reads. It must not naively clone the
existing `RefCountCache` once per worker.

The next slice gives each worker coarse affinity reservations, called
superblocks here. For the destinations reserved to one worker, derive local
consumer counts and retire that worker's source copies when those counts reach
zero. Dynamic pull remains inside a reservation. Work stealing happens only at
unstarted superblock boundaries, where ownership and local counts can move
together.

Only after measuring this design should the experiment consider source-load
DTasks whose results are Dagger chunks. That version exposes source arrays to
Dagger's data movement and may replace duplicated storage reads with duplicated
network transfers, so it is not the starting point.

## 6. Staged experiment

### D0 — serialization boundary, no regridding

- Start two local worker processes from the isolated environment.
- Initialize worker-local state from a plain synthetic-source configuration.
- Round-trip the complete small graph once per worker, not once per job.
- Submit no-op and counter-only batches to measure DTask creation, dispatch,
  report, and process-scope overhead.
- Prove that plans, providers, locks, and open stores are reconstructed rather
  than captured.

### D1 — tiny vertical slice

Run 8--32 synthetic destination chunks through one complete-destination task or
small batch at a time. Include:

- `Weighted` and `Extensive` reducers;
- ordinary, high-degree, missing, and all-NaN destinations;
- deliberately shuffled task completion;
- one injected failure before output commit and one after physical write but
  before ledger commit;
- one and two in-flight batches per process.

Compare every physical output chunk byte-for-byte with the regular path. Confirm
that destination arrays are not gathered into the coordinator and that graph
and worker state are not serialized with each task.

### D2 — fixed-core local process canary

Use the standing representative chunk set and compare, at a fixed total CPU
budget:

1. the current threaded scheduler;
2. Dagger with one process;
3. Dagger with 2 and 4 processes;
4. many single-threaded processes versus fewer multithreaded processes.

Record cold and warm wall time, compilation, task-dispatch overhead, CPU use, GC
wall per process, aggregate and peak RSS, source loads and duplicate loads,
weight builds, bytes read and transferred, cache hit rate, output throughput,
tail time, and byte identity.

The hypothesis is specifically that separate Julia heaps recover time currently
lost to stop-the-world GC. Thread-backed Dagger alone is not expected to improve
the mature threaded driver.

### D3 — graph-aware process ownership

- Add affinity superblocks and worker-local consumer counts.
- Preserve guided batch taper inside each worker reservation.
- Add stealing at superblock boundaries and record every ownership transfer.
- Compare no cache, bounded LRU, and reservation-local refcount caches.
- Measure the trade between load balance and duplicated source reads, including
  polar destinations with much higher graph degree.

### D4 — I/O, failure, and resume

- Keep store creation and metadata mutation on the coordinator.
- Permit workers to write only disjoint destination chunks.
- Make ledger commits centralized and idempotent.
- Make source download temporary paths and locking process-safe.
- Kill a worker during source acquisition, computation, physical write, and
  report delivery; verify that retry neither loses nor double-commits a chunk.
- Cover all-NaN chunks, whose successful output may have no physical chunk file
  and therefore depends on the ledger.

### D5 — multi-node canary

Proceed only if local multiprocess execution wins. Reopen source and output
stores independently on every node and measure shared-store bandwidth, source
read amplification, Dagger data movement, aggregate memory, and scaling
efficiency. The scheduler continues to assign affinity superblocks; it does not
submit the source/destination graph as a Dagger DAG.

## 7. Promotion gates

The separate backend remains experimental unless all of these hold:

- **Correctness:** every canary output is byte-identical to the regular path;
  graph demand remains a superset of actual reads; source accumulation order is
  unchanged.
- **Isolation:** regular regridding behavior remains unchanged; Dagger stays out
  of the root dependency graph; the backend runs only through an explicit
  `scripts/dagger_regrid/main.jl` invocation.
- **Overhead:** one-process Dagger is within approximately 5% of the current
  driver on representative batches.
- **Speed:** multiple processes produce a material repeatable improvement,
  proposed as at least 20% on the fixed-core canary before multi-node work.
- **Memory:** both peak-per-process and aggregate RSS fit the declared budget;
  no cache grows without a bound.
- **I/O:** output arrays do not travel through the coordinator and source-load
  amplification is reported and bounded.
- **Reliability:** failures retire scheduler ownership exactly once, retries are
  idempotent, and resume works for both physical and all-NaN chunks.

Passing those gates justifies a separate decision about an optional package
extension and an executor-neutral internal seam. It does **not** automatically
justify putting Dagger branches into the regular `regrid!` implementation.

## 8. Expected result

The likely win is conservative, compute-heavy regridding across multiple Julia
processes or machines. The present full-run attribution says clipping and weight
construction dominate while store time is small, and the utilization runs show
material single-process GC loss. Separate heaps can plausibly turn more of the
available cores into useful work.

The likely non-win is Dagger on the same threads, point/direct methods with little
computation per destination, or a cluster whose shared source store saturates
before its CPUs. This plan is structured to reject those cases cheaply without
touching the established regridding path.
