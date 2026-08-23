# Regridding simplification execution plan

- Original date: 2026-08-21
- Rewritten: 2026-08-22

This is the authoritative execution plan. The superseded evidence-rich version
is preserved in
`regrid-notes/2026-08-21-regridding-simplification-plan-detailed-archive.md`.
Use that archive for measurements, code citations, rejected alternatives, and
dependency-version history; use this file for task order and scope.

## Outcome

Global regridding should have one representation at each layer:

```text
space geometry
|- raster: DiskArrays ranges over CR's curvilinear quadtree
|- DGG: native grid cursor stopped at the chunk frontier
`- arbitrary subset: GeometryOps FlexibleRTree
                         |
                         v
              one chunk-query interface
                         |
                         v
             plan-owned dependency CSR
                         |
              +----------+----------+
              |          |          |
              v          v          v
          executor   scheduler   cache policy
              |
              v
      on-demand final WeightBlocks
              |
              v
       on-demand source data
```

Planning is eager about geometry and chunk dependencies. Weight construction
and source reads remain lazy.

## Current checkpoint

- Phase 0, environment and benchmark baseline: complete.
- Phase 1, native hierarchical chunk indexes: complete.
- Phases 1A, 1B, and 1C: **complete**. They finished the spatial layer before
  graph integration.
- Current implementation phase: **2**.
- The nearest/bilinear method redesign is deferred to a separate plan.

Landed cards, with the commit each one shipped as:

| card | commit | subject |
|---|---|---|
| A1 | `47a7cf7` | Simplify raster coordinate transformations |
| A2 | `41fb204` | Add task-local Proj raster transformations |
| A3 | `a07de52` | Converge raster cells on the CR quadtree |
| A4 | `e821c00` | Simplify raster range extent construction |
| B1/B2/B3 | `714f101` + `93e836d` | Adopt public spherical extent operations; Pin public spherical extent dependencies |
| B4 | `d194e2d` | Reuse packed trees for fallback cell geometry |
| C1 | `5e91b0c` | Consolidate the regridding space interface |
| G1 | PR #70 | Add dependency graph correctness and performance gates |
| G2 | `da9e737` | landed early via PR #69; see the card |

B1 and B2 are upstream commits in GeometryOps and ConservativeRegridding;
`93e836d` is the commit that pinned them and records their SHAs.

Evidence:

- `regrid-notes/2026-08-22-regridding-phase-0-baseline.md`
- `regrid-notes/2026-08-22-regridding-phase-1.md`
- `regrid-notes/2026-08-23-chunk-dag-coverage.md` — why G2 landed early
- `regrid-notes/2026-08-23-g1-graph-oracles.md` — G1's oracles, and the
  retroactive verdict on G2's gates
- `regrid-notes/2026-08-21-simplification-plan-review.md`
- `regrid-notes/2026-08-21-regridding-simplification-plan-detailed-archive.md`

## Decisions that tasks must preserve

1. GlobalRegridding is unit-spherical internally. A raster's native CRS may be
   geographic or projected, but its explicit transform maps to the unit sphere.
2. Proj integration uses Proj.jl's exposed API only. No direct `ccall`,
   `PROJ_jll` call, or reproduced binding is permitted.
3. Proj contexts and transformation clones are task-local, because Julia tasks
   may migrate between OS threads.
4. Raster storage chunk ownership comes from `DiskArrays.eachchunk`, never from
   coordinate lookups.
5. DGG chunk extents cover all stored leaves. Traverse the original grid cursor
   to `chunklevel`; do not rebuild chunk ancestors as literal leaves.
6. A `ChunkedPlan` owns or validates one immutable dependency relation.
   Scheduler, executor, prefetch, and cache policy consume that same object.
7. Planning reads no source values, builds no weights, and performs no provider
   network metadata requests.
8. `DirectPlan` and `ChunkedPlan` remain separate. `ShapedRegridArray` remains.
9. `cellindices` and `chunkranges` remain separate contracts: cell ownership is
   not the same as a rectangular storage read.
10. Do not create a general cache hierarchy. Simplify individual cache state
    only where a representation is genuinely duplicated.

## Standing acceptance laws

Every task must preserve all applicable laws:

1. No spatial index or dependency graph has geometric false negatives.
2. Source chunks are applied in deterministic ascending order.
3. Eager, lazy, and differently chunked plans produce equivalent values.
4. Reusing a plan reuses its dependency graph and built weight blocks.
5. Plan and lazy-array construction read no source values.
6. Declared budgets bound mutable residency, apart from documented immutable
   planning metadata and the one streamed scratch buffer.
7. Graph and executor agree exactly on permitted reads by sharing one graph.
8. A final `WeightBlock` carries one reusable reference vector; caches do not
   copy it.
9. Every cross-package `RegridSpace` specialization is part of the documented
   qualified extension interface.
10. Point-method behavior is unchanged by this plan.

## Sub-agent execution protocol

Each task below is a complete handoff unit for one sub-agent.

- Assign exactly one task ID and its listed evidence to the agent.
- Start from the commit produced by all prerequisites.
- Use one write agent at a time in the shared worktree. Tasks may run in
  parallel only in isolated worktrees and only when their file scopes do not
  overlap.
- The agent may make only the decisions explicitly delegated by its task.
  Escalate a failed gate or new architectural choice instead of expanding scope.
- Preserve unrelated user changes and begin with a clean status check.
- Run the narrow verification listed on the card. The phase-ending task also
  runs the phase gate.
- One task produces one focused commit. The final task commit in a phase is the
  required post-phase checkpoint.
- The handoff report states files changed, deletions achieved, tests and
  benchmarks run, results, and any remaining risks.

A suitable assignment prompt is:

```text
Execute task <ID> from
regrid-notes/2026-08-21-regridding-simplification-plan.md only.
Read its prerequisites and evidence, stay inside its scope, run its verification,
and stop if an exit condition cannot be demonstrated. Do not begin the next task.
```

## Dependency order

```text
A1 -> A2 -> A3 -> A4 ----+
                          +-> B3 -> B4 -> C1              Phases 1A-1C
B1 -> B2 -----------------+
C1 -> G1 -> G2 -> G3 -> G4 -> E1 -> E2                   Phases 2-4
E2 -> W1 -> W2                                             Phase 5
W2 -> O1 -> O2                                             Phase 6
O2 -> P1                                                   Phase 7
P1 -> L1 -> L2                                             Phase 8
```

B1 may be prepared in an isolated GeometryOps worktree while A tasks run, but
B2 must use the committed B1 API and B3 waits for both the upstream work and
the completed raster geometry. Do not integrate or pin provisional upstream
work.

The `G1 -> G2` edge is historical only: G2 landed *before* G1, as a correctness
fix (PR #69), and G1/PR #70 ran its gates retroactively. Read the order as
`C1 -> G2 -> G1 -> G3 -> G4` for what actually happened. Both cards are closed.

## Phase 1A — raster geometry and CRS

### Task A1 — normalize the raster transform contract

**Prerequisite:** completed Phase 1.
**Owns:** `rastergrid.jl` transform names, constructors, chart metadata, and
raster tests.

Actions:

- Replace local `LonLatToSphere` and `SphereToLonLat` with GeometryOps'
  `UnitSphereFromGeographic` and `GeographicFromUnitSphere`.
- Rename internal fields/helpers to `native_to_unit_sphere` and
  `unit_sphere_to_native`; name geographic-only and radian-returning helpers
  precisely.
- Normalize internal calls to the one-coordinate transformation convention.
- Adapt the documented legacy two-argument closure once at construction.
- Keep geographic degrees as the backwards-compatible default. Require an
  explicit transform or CRS for projected native coordinates.

Do not add Proj support or change point-method algorithms in this task.

Verify default cell rings, inverse coordinates, output lookups, arbitrary
callables, eager/lazy parity, and zero source reads during construction.

**Done when:** local geographic transformations are deleted and one explicit
native-to-unit-sphere contract serves every raster geometry caller.
**Commit:** `Simplify raster coordinate transformations`.

### Task A2 — add task-local Proj integration

**Prerequisite:** A1.
**Owns:** the Proj extension, dependency metadata, cleanup code, and extension
tests.

Actions:

- Add a core hook that leaves ordinary Julia/GeometryOps transformations
  unchanged.
- In a Proj extension, cache one cloned context and forward/reverse
  transformation pair per Julia task and template.
- Use only `Proj.proj_context_clone`, `Proj.proj_clone`,
  `Proj.Transformation`, and `Proj.proj_context_destroy`. Use the exposed
  `Proj.proj_destroy` only to release a non-null raw clone if construction fails
  before ownership transfers to a `Proj.Transformation`.
- Retain the mutable template in one shared chart state and serialize only the
  brief clone-pair construction from it. Operational contexts and clones are
  task-local and require no lock.
- Release transformations before their context; make cleanup idempotent and
  correct under partial construction failure.
- Construct projected raster pipelines with `always_xy = true` and compose
  native-to-geographic Proj conversion with GeometryOps geographic-to-sphere
  conversion.

Verify task isolation, reuse within one task, concurrent numerical identity,
round trips, failure cleanup, and absence of a direct PROJ `ccall` in this
package.

**Done when:** no operational Proj context or transformation clone is shared
between tasks and every task owns the complete lifetime of its clones. The
shared grid may retain the template through its locked chart state.
**Commit:** `Add task-local Proj raster transformations`.

### Task A3 — converge the structured raster tree

**Prerequisites:** A1 and A2.
**Owns:** raster cell-tree structure, `RasterGridView`,
`raster_tree_memo.jl`, and structural tests/benchmarks.

Actions:

- Make `RasterGridView` the sole structured raster lattice adapter.
- Make traversal create a task-prepared cursor view that hoists the current
  task's chart pair once. Do not make `Trees.getvertex` perform task-local
  storage lookup for every vertex, and do not store the prepared clone on a
  cursor or index that may be shared between tasks.
- Return restricted CR `TopDownQuadtreeCursor`s for whole and rectangular cell
  ranges; preserve global storage-order numbering.
- Delete `RasterCellTree` and remove `MemoRasterTree` if the committed benchmark
  shows no material regression. If caching is earned, adapt one small cache to
  the CR cursor instead of retaining another tree.
- Preserve the existing node-extent behavior during this structural cutover;
  A4 owns its semantic simplification.
- Leave scattered subsets on the temporary fallback for B4.

Verify whole and restricted cursor cell numbering, both storage orders,
single-build behavior, eager/lazy identity, allocation, and the existing
traversal benchmark.

**Done when:** one structured raster adapter replaces `RasterCellTree`, with no
tree or memoization regression.
**Commit:** `Converge raster cells on the CR quadtree`.

### Task A4 — simplify and prove raster range extents

**Prerequisite:** A3.
**Owns:** general raster range extents, geographic cap specializations, edge
tables, coverage laws, and extent benchmarks.

Actions:

- Delegate the general curvilinear path to CR's perimeter-based
  `cell_range_extent`.
- Remove finite `_boxcap` sampling as a purportedly generic conservative bound
  for arbitrary projected transforms.
- Retain geographic wide-longitude, antimeridian, polar, or sine/cosine edge
  specializations only when a coverage law or committed benchmark earns them;
  name them as geographic optimizations.
- Fix generally useful extent gaps in GeometryOps or CR instead of adding
  another general cap builder here.

Verify full-longitude, polar, antimeridian, projected, large-index rectangles,
boundary-point coverage, small-spatial-chunk traversal, allocation, eager/lazy
identity, and zero source reads.

**Phase 1A gate:** one structured raster adapter; GeometryOps owns geographic
conversion; Proj.jl owns all PROJ calls; every extent covers its leaves.
**Commit:** `Simplify raster range extent construction`.

## Phase 1B — spherical extents and fallback cell trees

### Task B1 — publish GeometryOps spherical extent primitives

**Prerequisite:** none beyond the pinned GeometryOps baseline.
**Repository:** GeometryOps.
**Owns:** unit-spherical cap utilities and tests.

Actions:

- Add public conservative cap-cap `Extents.intersects`, delegating to the
  existing robust predicate.
- Add one public outward-rounded Cartesian XYZ bound for a `SphericalCap`.
  Choose `Extents.extent` only if its semantics are honest; otherwise use an
  explicitly named unit-spherical helper.
- Publish the existing containment and cap-merge operations used across DGG and
  regridding. Add angular dilation only if it clarifies multiple callers.
- Do not publish exact cap-lens area solely for this refactor.

Verify ULP-near contact, coincident centres, tiny/polar/antimeridian/
hemisphere/whole-sphere bounds, containment, merge coverage, and cap-box
symmetry.

**Done when:** downstream packages need no private GeometryOps cap call or
local cap-to-XYZ implementation.
**Commit:** in GeometryOps, with its SHA recorded for B3.

### Task B2 — make ConservativeRegridding extent-generic

**Prerequisite:** committed B1.
**Repository:** ConservativeRegridding.
**Owns:** dual-search predicate selection, mixed-extent tests, and frontier cost
estimation.

Actions:

- Use `Extents.intersects` uniformly instead of selecting a private predicate
  for spherical manifolds.
- Test cap-cap, cap-box, box-cap, and box-box tree pairings.
- Measure the multithreaded frontier when a box node takes the generic
  `pair_weight`. Add a cheap mixed-extent estimate only if the production-like
  balance gate demonstrates it is needed.
- Do not change clipping geometry or sparse assembly.

**Done when:** mixed trees produce the same candidate/weight result without a
material task-balance regression.
**Commit:** in ConservativeRegridding, with its SHA recorded for B3.

### Task B3 — integrate the public spherical extent protocol

**Prerequisites:** A4 and B2; B2 transitively requires B1.
**Owns:** root workspace pins, GlobalRegridding/DGG cap helpers, discovery
predicates, and protocol integration tests.

Actions:

- Pin the committed B1/B2 revisions and re-resolve the environment.
- Replace local cap-to-XYZ, robust cap intersection, and cap merging with the
  public GeometryOps operations.
- Make `DilatedIntersects` delegate to a dilated cap and the public robust
  cap-cap predicate.
- Remove private GeometryOps cap calls from DiscreteGlobalGrids wrappers; delete
  wrappers that no longer add a DGG-specific contract.
- Leave fallback tree types in place for B4 so protocol and representation
  failures remain bisectable.

Verify generic chunk candidate identity, DGG cap coverage, tiny/contact/polar/
antimeridian support dilation, environment resolution, and all existing
fallback cell-tree results.

**Done when:** this workspace calls one public spherical extent protocol and
records the exact upstream SHAs.
**Commit:** `Adopt public spherical extent operations`.

### Task B4 — replace custom fallback cell trees

**Prerequisite:** B3.
**Owns:** generic/scattered cell fallbacks and their tests/benchmarks.

Actions:

- Add one thin cell-space adapter over `FlexibleRTrees.RTree`, retaining the
  owning space and original global cell positions.
- Use it only when a native restricted cursor is unavailable.
- Delete `CellCapTree` and `RasterFlatTree`'s scattered-cell role. Keep the
  latter only if the legacy `chunktree` bridge still requires it until E2.
- **Make "when a native restricted cursor is unavailable" true on the DGGS
  side** (folded in 2026-08-23; scope
  `regrid-notes/2026-08-23-dggs-subcursor-scope.md`, decided by the user). It is
  vacuous today: `subcursor` has one concrete method, CopernicusDEM
  (`src/systems/CopernicusDEM/cursor.jl:412`); every hierarchical DGGS takes the
  `= nothing` default at `src/interface/grid.jl:570`, so a non-chunk-aligned
  window packs an R-tree instead of windowing the hierarchy. Implement
  `subcursor` for hierarchical DGGS grids over the existing `SubtreeIds` O(1)
  contiguous-run view, routed through a `PartialGrid` for the descent clamp.
  This **subsumes `_chunkcursor` (`src/regridding.jl:250`) entirely** — a net
  deletion. Two hard requirements from the scope: the position shift must be
  **bidirectional** (a one-directional shift passes leaf-set and split-weight
  assertions while producing silently wrong weights), and the window must carry
  **cached leaf caps** — a bare cursor loses the join by 3-4x, worse than the
  fallback it replaces.
- Guard the destination-tiling cliff independently of the rest of B4, and first.
  Nothing in production reaches the fallback today (all four `subtree` call
  sites take case 1 or case 3, with 4.07x headroom in `_defaulttilesizes`), but
  at `budget = 2^27` every destination unit packs a ~419k-cell R-tree on 40
  workers: measured 1.213 s against 0.501 s for a cached window, at 4.3x the
  memory. That cliff, not a steady-state win, is the reason this is on the list.
- Not folded in: `_cachedcelltree`'s missing size guard (unlike
  `_cachedchunktree`). User reviewed 2026-08-23 and accepted it as-is.

Verify packing/node-capacity invariance, original leaf numbering, native versus
fallback Conservative block identity, threaded determinism, production
candidate count, task balance, time, and peak memory.

**Phase 1B gate:** GeometryOps owns reusable cap primitives, CR accepts mixed
extents, and GlobalRegridding has no custom general-purpose fallback cell tree.

**Commit:** `Reuse packed trees for fallback cell geometry`, recording both
upstream SHAs in the commit message.

## Phase 1C — the `RegridSpace` contract

### Task C1 — consolidate the qualified extension interface

**Prerequisite:** B4.
**Owns:** interface declarations/docs, module public/export lists, DGG extension
methods, and conformance tests.

Actions:

- Gather and document the contracts by responsibility:
  `ncells/getcell/celltree/subtree`,
  `nchunks/cellindices/chunkextent/chunkindex/candidatechunks!`,
  `chunkranges`, point/chart hooks, and
  `destinationdims/dimsource/_asspace`.
- Mark load-bearing qualified hooks public even when unexported.
- Correct the module comment: DGG already extends `chunkindex`,
  `candidatechunks!`, and `chunkranges`.
- Keep `cellindices` distinct from rectangular `chunkranges`.
- Do not add an abstract chunk-index type, trait bundle, plan core, or forwarding
  space wrapper.
- Leave `chunktree` documented as a temporary compatibility bridge for E2.

**Phase 1C gate:** every external specialization is named by one coherent
contract and no new interface abstraction exists.
**Commit:** `Consolidate the regridding space interface`.

## Phase 2 — one authoritative dependency graph

### Task G1 — establish graph correctness and performance oracles

**Prerequisite:** C1.
**Owns:** no-miss instrumentation, graph benchmark harness, and oracle tests.
Production graph construction remains unchanged.

Actions:

- Add an actual-cell small-space oracle that proves every geometrically
  contributing pair is retained.
- Add brute-force cap identity checks for the generic R-tree path. **Do not
  treat the flat cap relation as an upper bound for raster/DGG native
  hierarchies** — this card originally said to, and #70 measured and pinned the
  opposite. The generic packed-R-tree index equals the cap join exactly; the DGG
  hierarchical descent is a strict subset of it; the `RasterGrid` quadtree and
  the CopernicusDEM level-0 frontier both *cross* it, holding pairs the cap join
  rejects as well as missing pairs it holds. The cap join bounds the graph in
  neither direction on either shipped native hierarchy. What both relations do
  satisfy is `truth ⊆ ·`, and that is the invariant to assert.
- Commit a repeatable comparison harness for the latitude join and native-index
  row builder, including complete plan time, graph time/bytes, allocations, and
  peak memory.
- Record current results on small spatial raster chunks, complete/rooted/sparse
  DGGs, polar/antimeridian cases, nonzero support, nonuniform coverage, and the
  production Copernicus-to-IGeo7 pair.

**Done when:** the next task can change builders and obtain an objective
correctness/performance verdict without adding instrumentation.
**Commit:** `Add dependency graph correctness and performance gates`.
**Landed:** PR #70. Records: `regrid-notes/2026-08-23-g1-graph-oracles.md` (the
oracles, the harness, and the retroactive verdict on G2) and
`regrid-notes/2026-08-23-chunk-dag-coverage.md` (the coverage defect that made
#69 land G2 early). Post-#69 the contract the gates assert is an **equality**
between the graph and the demanded relation, not a containment; G3 and G4 use
`benchmark/chunk_graph_gates.jl` and the shared oracles in
`lib/GlobalRegridding/test/graphoracles.jl` to prove they did not change it.

### Task G2 — cut over to indexed graph construction — CLOSED

**Status: LANDED EARLY, OUT OF ORDER, AND CLOSED.** Not an implementation card.
Do not reopen it.

It shipped as PR #69 (`da9e737`, builder commit `ba2bbfa`, *Build the chunk
dependency graph from the source space's own chunk index*) — before G1 existed,
and therefore without either of the gates this card said to pass first. It landed
as a **correctness fix**, not as the planned performance-neutral cutover: a graph
built from `chunkextents` directly *crossed* the executor's own candidate
relation instead of dominating it, so a refcount derived from it retired source
tiles a subsequent read still demanded. Record:
`regrid-notes/2026-08-23-chunk-dag-coverage.md`.

What shipped: destination-major rows from one source `chunkindex` and independent
`candidatechunks!` queries; deterministic sorted rows and the bidirectional Int32
CSR preserved; the latitude-sorted join deleted.

Both gates were then run retroactively by G1/PR #70, on the case matrix this card
asked for (small spatial raster chunks, complete/rooted/sparse DGGs, polar and
antimeridian cases, nonzero support, nonuniform coverage, and the production
Copernicus-to-IGeo7 pair). Verdict in
`regrid-notes/2026-08-23-g1-graph-oracles.md` §5:

- **Correctness gate: PASSED**, retroactively. The indexed relation holds every
  geometrically contributing pair in every oracle-checked case at every radius,
  and every demanded pair on the production pair — where the deleted cap join
  misses 72.
- **Performance gate: FAILED, AND WAIVED.** This card required recovering the
  archived 0.0594 s latitude-join median before removing that builder. It was not
  recovered and cannot be; the two builders do different work. Measured on the
  production pair: **5.7× slower at t8** (0.1219 s vs 0.0213 s), **6.1× at t4**
  (0.2289 s vs 0.0376 s), and **1.6× the allocations** (17.1 MB vs 10.6 MB). Per
  destination it is worse on a shallow raster source — up to 60× on
  4320×2160/162 chunks.

The waiver, and its reasoning: 0.12 s against a run whose recorded wall time is
8.81 h — about 4 × 10⁻⁶ of it, and 1/90th of the space construction that precedes
it — bought in exchange for a relation a refcount can actually be derived from
(72 → 0 demanded-but-unheld pairs) and 322 fewer edges. The gate is recorded as
*failed and waived*, never as passed. If a later task wants the factor back, the
measurement points at the `RasterGrid` arm's per-query cursor copy and quadtree
descent; `benchmark/chunk_graph_gates.jl` with
`DGG_GRAPH_GATE_CASES=raster-4320-162chunks` is the reproducer, and its
`:latjoin` arm keeps the deleted builder runnable for as long as the waiver
stands.

**Commit:** `Build the chunk dependency graph from the source space's own chunk
index` (`ba2bbfa`, merged as `da9e737`).

### Task G3 — give dependency objects identity and row views

**Prerequisite:** G2.
**Owns:** `ChunkDependencyGraph` identity, immutable restriction/view support,
and tests.

Actions:

- Stamp chunk counts, support radius, source/destination geometry identity or
  fingerprints, and a serializable narrow-phase tag.
- Add validated `restrict(graph, destination)` row views for per-column plans,
  sharing the parent relation rather than rebuilding it.
- Retain `sourcesof`, `consumersof`, degrees, counts, and bidirectional CSR.
- Keep Graphs.jl compatibility over the same CSR; do not create a second graph.

Row views are worth **strictly more** than when this card was written. Post-#69,
`chunkindex(src_space)` is rebuilt on *every* `chunk_dependency_graph` call
(`lib/GlobalRegridding/src/discovery.jl:40` — the generic path is a fresh packed
R-tree over `chunkextents`), so a per-column rebuild now costs a full index build
per column, not a `sortperm` over caps that were already sorted. Sharing the
parent relation avoids that entirely.

**Done when:** an invalid graph reuse fails at construction and row views retain
global destination identity for refinement.
**Commit:** `Add reusable dependency graph identities`.

### Task G4 — make `ChunkedPlan` the sole graph owner

**Prerequisite:** G3.
**Owns:** `plan_regrid`, `ChunkedPlan`, graph constructors/accessors, production
plan construction, and tests. Executor discovery remains unchanged.

Actions:

- Make `refine` a keyword of lazy `plan_regrid` and nowhere else.
- Construct or validate dependencies once in the plan; expose
  `dependencies(plan)` as a non-building accessor.
- Remove `refine` from post-plan graph builders.
- Make production global/per-column plans share the global graph or validated
  row views; never rebuild a one-destination graph per column.
- Delete `connectedchunkpairs` (`lib/GlobalRegridding/src/discovery.jl:203-217`).
  The archive scheduled it for Phase 4 alongside `connectedchunks`, but #69 made
  it a free deletion here: it is now line-for-line the same loop as
  `_chunkgraph`'s rows — one `chunkindex(src)`, one `candidatechunks!` per
  `chunkextents(dst)` entry — and its only caller is
  `lib/GlobalRegridding/test/test_lazy.jl:145`. Move that assertion onto the
  graph's rows and drop the function.
- Preserve zero source reads, zero weight builds, deterministic construction,
  and no network metadata work.

**Phase 2 gate:** one logical plan exposes exactly one validated relation and a
narrow phase cannot be supplied after the plan exists.
**Commit:** `Make chunked plans own dependency graphs`.

## Phase 3 — graph-backed lazy execution

### Task E1 — consume graph rows in the lazy executor

**Prerequisite:** G4.
**Owns:** lazy source selection, tile-row unions, wave metadata ownership,
executor tests, and production graph-miss validation.

Actions:

- For a space-aligned destination tile, use `sourcesof(dependencies(plan), d)`.
- For a derived tile spanning several destination chunks, take a deterministic
  sorted union of its graph rows.
- Keep data-dependent `knownempty` filtering after adjacency selection.
- Remove per-read spatial index queries and `LazyRegridArray.srcindex`.
- Move cap metadata needed only for wave costing onto the immutable dependency
  object or an on-demand accessor; do not replicate vectors per worker/column.
- Prove `nonzero eager edges ⊆ graph rows ⊆ broad cap relation` on small spaces.

Verify repeated reads/slices, explicit and budget tiling, no silent graph
misses, source order, eager/lazy values, residency, and production policy.

**Phase 3 gate:** no lazy read performs geometric dependency discovery; executor
and scheduler consume the same object.
**Commit:** `Drive lazy regridding from dependency rows`.

## Phase 4 — delete legacy discovery

### Task E2 — retire duplicate discovery state and `chunktree`

**Prerequisite:** E1.
**Owns:** discovery compatibility functions/types, module exports, tests, and
scripts.

Actions:

- Remove `connectedchunks` and `connectedchunks!`. (`connectedchunkpairs` moved
  forward to G4 — after #69 it duplicates `_chunkgraph`'s rows exactly, so it is
  a free deletion there.)
- Remove any latitude join left after G2, post-plan graph builder, redundant cap
  vectors, and compatibility source-index state.
- Keep a cheap `chunkextent(space, chunk)` and `chunkextents` materialization
  only for real diagnostics/planning consumers.
- Audit repository and known downstream use of `chunktree`. If none remains,
  remove its export, declaration, Raster compatibility tree, and tests.
- Rewrite oracle tests against `sourcesof(dependencies(plan), d)`.

**Phase 4 gate:** one query implementation defines graph edges, and neither the
executor nor interface translates that relation back through a compatibility
tree.
**Commit:** `Remove duplicate chunk discovery paths`.

## Phase 5 — one final weight block representation

### Task W1 — normalize denominator and reference storage

**Prerequisite:** E2.
**Owns:** `WeightCOO`, `WeightBlock`, `CachedBlock`, executor reference helpers,
spill reconstruction, and tests.

Actions:

- Change `WeightCOO` to `denom::Union{Nothing,Vector}` and remove `hasdenom`.
  Allocate on `markdenominated!` or first `adddenom!`.
- Make each `WeightBlock` own its final reusable reference vector. Alias it to
  `denom` for denominated blocks; otherwise compute row sums once.
- Remove `CachedBlock.ref`; use `entry.block.reference`.
- Reuse the stored reference across eager plan applications.
- Reconstruct, rather than serialize, a redundant reference in `Spilled` and
  count aliased arrays once in `_blockbytes`.

Verify empty sides, zero-coverage rows, point blocks, sparse/dense blocks,
repeated application allocations, cache hits, spill round trips, and byte
accounting.

**Done when:** weights, optional denominator semantics, and one reference vector
live in the final block; no cache copies numerical reference state.
**Commit:** `Store reference weights in WeightBlock`.

### Task W2 — unify eager and chunked block construction

**Prerequisite:** W1.
**Owns:** method construction seam, Conservative specialization, plan builders,
and weight tests/benchmark.

Actions:

- Add `weightblock(method, dst_space, dst_inds, src_space, src_inds)`.
- Keep `build_weights!` plus `WeightCOO` as the compatible generic fallback.
- Let Conservative adopt its assembled CSC directly with denominators computed
  once; remove CSC-to-COO-to-CSC conversion.
- Make eager and chunked builders call the same seam; remove the special
  `wholeblock(::Conservative)` path.
- Preserve wrapper/instrumentation dispatch through an explicit forwarding
  rule or trait.

Verify sparse structure/value identity, partition invariance, denominators,
empty sides, third-party emitters, and the committed peak-RSS benchmark.

**Phase 5 gate:** one block-builder path and one final block representation,
with no measured eager regression and the chunked round-trip memory removed.

**Commit:** `Unify final weight block construction`.

## Phase 6 — generic DGG output and API forwarding

### Task O1 — remove DGG-specific output application

**Prerequisite:** W2.
**Owns:** `DGGSpace.destinationdims`, eager/lazy wrappers, DGG regrid methods,
and cross-system tests.

Actions:

- Define the generic DGG destination dimension as
  `(Cells(CellLookup(space.grid)),)`.
- Remove `_DirectToDGG`, `_ChunkedToDGG`, specialized `GR.regrid` methods, and
  `_ascube`.
- Add the one-axis lazy short circuit so DGG output does not create a pointless
  `ShapedRegridArray`.
- Preserve non-spatial dimensions, lookups, parent type, values, and eager/lazy
  parity.

**Done when:** DGG destinations use only generic output wrapping.
**Commit:** `Use generic output dimensions for DGG regridding`.

### Task O2 — remove target and keyword duplication

**Prerequisite:** O1.
**Owns:** DGG target resolution, `regrid`/`regrid!` convenience methods,
`plan_regrid` validation, docs, and tests.

Actions:

- Replace `regridgrid` plus `RegridTarget` with direct `DGGSpace` constructors
  or `_asspace` methods for each target spelling.
- Retain the meaningful bare-system error and source-dependent level choice.
- Make `plan_regrid` the single executable owner of keyword defaults and
  validation; forward `kwargs...` from `regrid` and `regrid!`.
- Preserve explicit user-facing keyword documentation and error quality.

**Phase 6 gate:** no DGG-specific plan application or intermediate conversion
generic remains, and API defaults occur once.
**Commit:** `Simplify regridding target and keyword resolution`.

## Phase 7 — prepared destination geometry

### Task P1 — replace the destination cache wrapper stack

**Prerequisite:** O2.
**Owns:** `TileCells`, `CachedCellTree`, destination `CellMemo`, block
construction, cache tests, and Conservative benchmark.

Actions:

- Introduce one private prepared destination tile holding its index set/local
  map, one restricted tree, and optional prebuilt polygons.
- Share it across all source blocks for the executor tile.
- Let `BlockAreaOperator` obtain destination polygons directly from it.
- Remove `TileCells <: RegridSpace` and its forwarding methods,
  `CachedCellTree`, and destination `CellMemo` when polygons are prepared.
- Retain a task-local memo for unprepared source geometry and oversized
  destinations.
- Remove stale raster/DGG extent caches only when their individual benchmark
  gate earns deletion; do not merge polygon and extent caches by name.

Verify one tree build per tile, cached/uncached identity, concurrent safety,
cache threshold, allocations, time, and peak memory.

**Phase 7 gate:** destination geometry follows one preparation/cache path with
no production regression.
**Commit:** `Consolidate destination geometry preparation`.

## Phase 8 — private lazy state

### Task L1 — give `SourceHold` one entry dictionary

**Prerequisite:** P1.
**Owns:** source-hold entries, eviction, residency statistics, and lazy tests.

Actions:

- Replace parallel `held` and `used` dictionaries with entries containing
  `value`, `bytes`, and `used`.
- Touch and evict through the entry; subtract its stored byte count.
- Do not unify `SourceHold` and `PerChunk`. Share a tiny oldest-key helper only
  if both call sites become simpler without parameterizing their semantics.

Verify load/hit/skip/drop counts, current/peak bytes, oversize streaming,
eviction order, budget limits, and acceptance residency.

**Done when:** one key has one source-residency record and no generic cache
abstraction was introduced.
**Commit:** `Simplify lazy source residency entries`.

### Task L2 — normalize DiskArrays metadata once

**Prerequisite:** L1.
**Owns:** lazy constructor chunk metadata, pass-through/source group helpers,
chunked-source naming, and tests.

Actions:

- Normalize a valid source `GridChunks` once per lazy-array application.
- Derive output pass-through chunks and non-spatial read groups from that one
  value.
- Reuse it for inferred raster spatial chunks only if this does not couple
  `RasterGrid` to `LazyRegridArray` or put data-specific state on the plan.
- Rename `isdiskbacked` to reflect that it tests declared chunking. Preserve the
  separate true-disk predicate used for DiskArrays reads.

Use a counting fixture to verify one `eachchunk` interpretation and test
regular, irregular, absent, malformed/non-grid chunk descriptions, output
tiling, read groups, defaults, and zero source reads.

**Phase 8 gate:** lazy source metadata and residency each have one private
representation, with unchanged behavior and no new framework.
**Commit:** `Normalize lazy DiskArrays chunk metadata`.

## Phase gates

At the final task of every phase, run the narrow tests above and the complete
applicable gate:

| Gate | Required verification |
|---|---|
| GlobalRegridding | complete package tests |
| DGG integration | cross-system regrid and acceptance tests |
| Conservation | conservation suite with only documented broken cases |
| CopernicusDEM | system suite and zero-source-read construction |
| Spatial changes | no-false-negative property tests and small-chunk benchmark |
| Graph changes | deterministic CSR, graph bytes/time, actual-cell no-miss proof |
| Weight/cache changes | production Conservative time, allocations, and peak RSS |
| Proj changes | concurrent task isolation and lifetime/failure tests |

Performance reports must state thread count, shape/parallelism mode, input
sizes, warmup, statistic, allocations, and peak-memory method. The detailed
archive lists the current reference measurements and exclusions, including the
hybrid-source production columns that are invalid byte-identity references.

## Definition of done

- Raster geometry uses GeometryOps transformations, task-local Proj.jl clones,
  and one CR curvilinear-grid adapter over DiskArrays-owned ranges.
- GeometryOps owns public spherical extent primitives; CR accepts mixed extent
  trees; arbitrary cell subsets use one `FlexibleRTree` adapter.
- The qualified `RegridSpace` contract names every external specialization.
- A reusable `ChunkedPlan` owns one dependency CSR consumed by execution,
  scheduling, prefetch, and cache policy.
- Lazy reads perform no geometric dependency discovery.
- Eager and chunked methods build the same final `WeightBlock`; Conservative
  performs no CSC-to-COO-to-CSC round trip; reference vectors are not copied by
  caches.
- DGG output and target resolution use generic paths, and `plan_regrid` alone
  owns keyword defaults.
- Destination polygon geometry has one prepared-tile path.
- Lazy source residency and DiskArrays chunk metadata each have one private
  representation.
- The point-method redesign remains independently actionable.
- All correctness suites and representative production/small-chunk performance
  gates pass, and every phase-ending commit records its evidence.
