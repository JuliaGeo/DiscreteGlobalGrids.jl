# Regridding simplification execution plan

- Original date: 2026-08-21
- Rewritten: 2026-08-22

This is the authoritative execution plan for conservative regridding and for
the shared spatial, planning, and execution infrastructure every method uses.
The superseded evidence-rich version is preserved in
`regrid-notes/2026-08-21-regridding-simplification-plan-detailed-archive.md`.
Use that archive for measurements, code citations, rejected alternatives, and
dependency-version history; use this file for task order and scope.

The point-method redesign is owned by
`regrid-notes/2026-08-23-barycentric-regridding-plan.md`. That plan names the
point API and its kernels; this plan owns the seams in the shared code that
admit them, as Phase 9 below. Neither plan changes conservative behaviour.

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
      +-------+--------+
      |                |
      v                v
 area methods     point methods
 one WeightBlock  one TileWeights per destination
 per chunk pair   tile, with an exact source-chunk list
      |                |
      +-------+--------+
              v
       on-demand source data
```

Planning is eager about geometry and chunk dependencies. Weight construction
and source reads remain lazy. A point method is the one exception to the
assumption that an exact chunk relation can always be materialized before
weights: its exact dependency list and its weights fall out of the same
destination pass, so they are built together, per tile, on first request.
Planning still reads no source values.

## Vocabulary

The index words are `DiscreteGlobalGrids`' own and are used identically here
and in the barycentric plan's `## Vocabulary`, which is the reference text for
the point-method terms. "Position" names a geometric location and never an
index.

| term | meaning |
|---|---|
| local index | a collection's `1:ncells(space)`. What `cellat(space, p)` answers, what `cellcentroid` and `chunkat` take, what `ownedindices` lists, and the space every stencil is written in. A `DGGSpace` may wrap a partial grid, so this is not the complete level's numbering. |
| chunk-local index | the local index of one chunk read as a partial grid: `j` with `ownedindices(space, k)[j] == i`, which is `i - first(ownedindices(space, k)) + 1` for a contiguous chunk. `WeightCOO` rows and columns and a `WeightBlock`'s axes are chunk-local. |
| global index | the complete level's numbering (`globalindex`, `Index(Global())`). It equals the local index only on a complete grid. |
| axis index | the output cube's full cell axis read as a local index. |
| chunk | a space's own unit of ownership, `ownedindices(space, k)` for `k` in `1:nchunks(space)`. |
| tile | the lazy array's write block: one destination chunk, or a run of them (`DestTiling`, `lib/GlobalRegridding/src/lazy.jl:56`). |
| partial grid | a collection holding fewer cells than its level. |
| sample site | where a cell's value is taken to sit, on either side: `cellcentroid(space, i)` unless the space says otherwise. A point method interpolates between source sample sites and evaluates at destination sample sites; `samplesites(space)` answers both. |
| rim | an *area* — the part of border cells lying outside the dual complex, where no dual cell exists. Distinct from `border`/`halo`, which are sets of cells. |

Four hooks were renamed after the cards below were written; each old name is
kept as a deprecation shim and appears in this plan only as history:
`cellindices` → `ownedindices` (`spaces.jl:79`, shim at `:94`),
`build_weights!` → `buildweights!` (`methods.jl:131`, shim at `:159`),
`support_radius` → `supportradius` (`methods.jl:145`, shim at `:171`), and
`chartposition` → `chartlocalindex` (`spaces.jl:258`, shim at `:287`).

## Current checkpoint

- Phase 0, environment and benchmark baseline: complete.
- Phase 1, native hierarchical chunk indexes: complete.
- Phases 1A, 1B, and 1C: **complete**. They finished the spatial layer before
  graph integration. One action was folded into B4 after B4 had already
  shipped and has not landed; it is carried below as B5, the only open item in
  Phases 1A-1C.
- Phase 2, one authoritative dependency graph: **complete**. G1, G2, G3 and G4
  are all closed. `ChunkedPlan` owns the relation; `dependencies(plan)` is the
  only way to read it and builds nothing; `refine` is a keyword of lazy
  `plan_regrid` alone.
- Phase 3, graph-backed lazy execution: **complete**. E1 is closed. The lazy
  executor takes a tile's source chunks from `sourcesof(dependencies(plan), d)`
  — on the chunk-pair route; a tile with `TileWeights` reads its manifest, see
  Phase 9 — and performs no dependency discovery; `LazyRegridArray.srcindex`
  and its cap vectors are gone, per-chunk caps now live on the relation, and a
  chunked plan therefore owns a relation **by default** — `dependencies = false`
  is the opt-out, and a plan that takes it cannot back a `LazyRegridArray`.
  G4's one open action closed with it, and production declined the result; see
  E1.
- Phase 4, delete legacy discovery: **complete**. E2 is closed. `chunktree`,
  `RasterFlatTree`, `connectedchunks`, `connectedchunks!` and the
  `chunktree`-collecting `chunkextents` fallback are all gone; `chunkextents` is
  a required hook with no fallback, and one `candidatechunks!` query per
  destination cap defines every edge.
- Phase 5, one final weight block representation: **complete**. W1 and W2 are
  closed. A final `WeightBlock` holds weights, the optional denominator and the
  one reference vector its values are normalized against — aliased to the
  denominator where there is one, row sums where there is not, stored once and
  copied by no cache. `weightblock` is the one builder path, for the eager whole
  domain and a chunk pair alike, and `Conservative` reaches it by specializing
  `pairblock`: it adopts the matrix it assembles and reads each denominator off
  that matrix once, so no conservative block is copied through a coordinate list
  on any route. `buildweights!` with `WeightCOO` remains the generic assembly
  and the only hook a method must supply, and a method wrapping another takes
  the inner method's build by forwarding `pairblock`.
- Phase 6, generic DGG output and API forwarding: **complete**. O1 and O2 are
  closed. A DGG destination supplies one output hook — `destinationdims`,
  answering the cells' own axis — and one `_asspace` method per target
  spelling; this package applies no plan, converts no target through an
  intermediate generic, and labels no result of its own on either route, and a
  destination whose one axis is the whole of its shape is labelled rather than
  reshaped or wrapped. `plan_regrid` states every keyword, every default and
  every check; `regrid` and `regrid!` forward `kwargs...` to it and declare
  nothing, refusing only the three keywords that describe a relation a kept
  plan owns.
- Phase 9, point-method admission: **complete**. S1, S2 and S3 are closed, and
  the barycentric plan's P0, P1 and P2 landed inside it. `outputsampling` names
  the build path by trait alone; a point method supplying a `sampler` builds one
  `TileWeights` per destination tile, cached and locked by tile number; and a
  tile reads exactly the source chunks its stencils name — its own manifest,
  neither the relation's row nor an intersection of the two. A point plan still
  owns a relation, for tile order, wave costing, refcounts and prefetch, and its
  rows are a documented superset that decides no read; `supportradius` is the
  declared bound that keeps them a superset, and a manifest naming a chunk no
  row of its tile holds is refused rather than silently trimmed. Conservative
  weights, keys, read counts and values are unchanged.
- **Next in this plan's own order: Phase 7** (D1), behind the user's gate.
  Phase 9 branched from E2 beside Phase 5 and is closed, so nothing runs
  alongside Phase 7; Phase 8 follows it in order.

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
| G3 | PR #71 | Add reusable dependency graph identities |
| G4 | PR #72 | Make chunked plans own dependency graphs |
| E1 | PR #74 | Drive lazy regridding from dependency rows |
| E2 | PR #75 | Remove duplicate chunk discovery paths |
| S1 | `8057af2` | Split weight construction on output sampling |
| P0 | `df1d8b8` | Record the point-method baseline |
| P1 | `e57bd82` | Add the barycentric point contracts and kernels |
| P2 | `874bf74` | Specialise barycentric points on chart axes |
| S2 | `8c05a2f` | Build point weights per destination tile |
| S3 | `2800ba3` | Read exact source chunks for point tiles |
| W1 | `d01da5c` | Store reference weights in WeightBlock |
| P3 | `d227de2` | Record the fused tile-weight execution |
| W2 | `6a4314f` | Unify final weight block construction |
| O1 | `619db4e` | Use generic output dimensions for DGG regridding |
| O2 | `fad29f6` | Simplify regridding target and keyword resolution |

B1 and B2 are upstream commits in GeometryOps and ConservativeRegridding;
`93e836d` is the commit that pinned them and records their SHAs. P0, P1, P2 and P3
are the barycentric plan's own cards, landed here because Phase 9 admits them.

Evidence:

- `regrid-notes/2026-08-22-regridding-phase-0-baseline.md`
- `regrid-notes/2026-08-22-regridding-phase-1.md`
- `regrid-notes/2026-08-23-chunk-dag-coverage.md` — why G2 landed early
- `regrid-notes/2026-08-23-g1-graph-oracles.md` — G1's oracles, and the
  retroactive verdict on G2's gates
- `regrid-notes/2026-08-23-g3-graph-identity.md` — G3's identity design, what it
  catches and what it cannot, and the row-view measurements
- `regrid-notes/2026-08-23-g4-plan-owns-graph.md` — G4's ownership design, the
  Phase 2 gate as tests rather than prose, and what a per-column relation would
  have cost production
- `regrid-notes/2026-08-23-e1-graph-backed-lazy.md` — E1's row-driven executor,
  where cap metadata moved, `subspace_dependencies` and why production declined
  it, and the Phase 3 gate as tests
- `regrid-notes/2026-08-23-e2-delete-legacy-discovery.md` — E2's deletion audit,
  the `chunktree` verdict, what the card listed that was not dead, and the
  Phase 4 gate as tests
- `regrid-notes/2026-08-25-s1-split-on-sampling.md` — S1's trait-selected build
  path, and what stayed on `buildweights!`
- `regrid-notes/2026-08-25-p0-point-baseline.md` — the pinned `BilinearPoint`
  behaviour and the counted repeat of point location per candidate chunk
- `regrid-notes/2026-08-25-p1-point-contracts.md` — the point seam: `WeightRow`,
  `sampler`, `weightsat!`, dual cells and the coordinate kernels
- `regrid-notes/2026-08-25-p2-raster-q1.md` — the chart specialization, every
  intentional difference from `BilinearPoint`, and why the support radius is
  declared on prepared state
- `regrid-notes/2026-08-25-s2-tile-weights.md` — S2's `TileWeights`, the tile as
  cache and locking unit, and the budget and spill round trips
- `regrid-notes/2026-08-25-s3-exact-reads.md` — S3's selection: what decides a
  read on each route, what the relation is for, and why `supportradius` is a
  bound rather than an approximation
- `regrid-notes/2026-08-25-w1-block-reference.md` — W1's final block: the one
  reference vector, when it aliases a denominator, and what a builder hands it
- `regrid-notes/2026-08-25-w2-one-builder.md` — W2's one build path,
  Conservative's adopted assembly, and the wrapper forwarding rule
- `regrid-notes/2026-08-25-o1-generic-output.md` — O1's one output hook, the
  two wrappers' one-axis short circuits, and the comparator behind them
- `regrid-notes/2026-08-25-o2-target-keywords.md` — O2's one `_asspace` method
  per target spelling, the bare-system level rule and its error, and where
  every keyword default now lives
- `regrid-notes/2026-08-23-barycentric-regridding-plan.md` — the point-method
  plan Phase 9 admits, and the authoritative naming for its parts
- `regrid-notes/generic-barycentric-patch-regridding.md` and
  `regrid-notes/esmf-regridding-methods.md` — the point-method use cases
- `regrid-notes/2026-08-24-go-fh-clipping-eval.md` — PR #77's GeometryOps pin:
  GO-WITH-CAVEATS as an off-driver improvement, NO-GO as a clipping
  optimisation, because the conservative weight kernel is
  `GO.ConvexConvexSutherlandHodgman`
  (`lib/GlobalRegridding/src/intersection_area.jl:31-32`), which that branch
  does not touch. The pin informs the Phase 5 gate and does not open it.
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
   For a point method that relation orders, costs and prefetches; what is
   actually read is the tile's exact manifest, which must be a subset of it.
7. Planning reads no source values, builds no weights, and performs no provider
   network metadata requests.
8. `DirectPlan` and `ChunkedPlan` remain separate. `ShapedRegridArray` remains.
9. `ownedindices` and `chunkranges` remain separate contracts: cell ownership is
   not the same as a rectangular storage read.
10. Do not create a general cache hierarchy. Simplify individual cache state
    only where a representation is genuinely duplicated.
11. Conservative regridding's behaviour is fixed. Every place the point path
    diverges from it is named at the seam where it diverges, and the
    conservative branch of that seam is the one that exists today.

## Standing acceptance laws

Every task must preserve all applicable laws:

1. No spatial index or dependency graph has geometric false negatives.
2. Source chunks are applied in deterministic ascending order.
3. Eager, lazy, and differently chunked plans produce equivalent values.
4. Reusing a plan reuses its dependency graph and built weight blocks.
5. Plan and lazy-array construction read no source values.
6. Declared budgets bound mutable residency, apart from documented immutable
   planning metadata and the one streamed scratch buffer.
7. Graph and executor agree exactly on permitted reads. An area method shares
   one graph; a point tile's exact manifest is a subset of its graph rows and
   bounds its reads.
8. A final `WeightBlock` carries one reusable reference vector; caches do not
   copy it.
9. Every cross-package `RegridSpace` specialization is part of the documented
   qualified extension interface.
10. Conservative values, read counts and weights are unchanged by every card
    here. Point-method behavior is unchanged by Phases 1A-8; Phase 9 changes
    where a point method's weights are built and which chunks its tiles read,
    never a weight.

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
B1 -> B2 -----------------+          `-> B5
C1 -> G1 -> G2 -> G3 -> G4 -> E1 -> E2                   Phases 2-4
E2 -> W1 -> W2                                             Phase 5
W2 -> O1 -> O2                                             Phase 6
O2 -> D1                                                   Phase 7
D1 -> L1 -> L2                                             Phase 8
E2 -> S1 -> S2 -> S3                                       Phase 9
```

B1 may be prepared in an isolated GeometryOps worktree while A tasks run, but
B2 must use the committed B1 API and B3 waits for both the upstream work and
the completed raster geometry. Do not integrate or pin provisional upstream
work.

The `G1 -> G2` edge is historical only: G2 landed *before* G1, as a correctness
fix (PR #69), and G1/PR #70 ran its gates retroactively. Read the order as
`C1 -> G2 -> G1 -> G3 -> G4` for what actually happened. Both cards are closed.

Phase 9 branches from E2 beside Phase 5 and depends on nothing in Phases 5-8,
but S2 and S3 share `plans.jl` and `lazy.jl` with W1, W2 and L1, so they take
the shared worktree in sequence rather than concurrently. S2 also requires the
barycentric plan's P1 (`weightsat!`, `WeightRow`, `sampler`) to exist.

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
**Landed:** `47a7cf7`.

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
**Landed:** `41fb204`. The two preparation hooks it added are public:
`_prepare_raster_transform_pair` and `_task_prepared_raster_transform`
(`lib/GlobalRegridding/src/GlobalRegridding.jl:120`).

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
  ranges; preserve storage-order numbering of the space's local indices.
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
**Landed:** `a07de52`. `RasterGridView` is
`lib/GlobalRegridding/src/rastergrid.jl:783`; `RasterCellTree` and
`MemoRasterTree` are gone.

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
**Landed:** `e821c00`. `Trees.cell_range_extent(::RasterGridView, ...)` is
`rastergrid.jl:829`; `_boxcap` is gone.

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
**Landed:** upstream; the SHA is in `93e836d`.

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
**Landed:** upstream; the SHA is in `93e836d`.

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
**Landed:** `714f101` + `93e836d`. `DilatedIntersects` is
`lib/GlobalRegridding/src/discovery.jl:6`.

### Task B4 — replace custom fallback cell trees

**Prerequisite:** B3.
**Owns:** generic/scattered cell fallbacks and their tests/benchmarks.

Actions:

- Add one thin cell-space adapter over `FlexibleRTrees.RTree`, retaining the
  owning space and the cells' original local indices.
- Use it only when a native restricted cursor is unavailable.
- Delete `CellCapTree` and `RasterFlatTree`'s scattered-cell role. Keep the
  latter only if the legacy `chunktree` bridge still requires it until E2.
- Not folded in: `_cachedcelltree`'s missing size guard (unlike
  `_cachedchunktree`). User reviewed 2026-08-23 and accepted it as-is.

Verify packing/node-capacity invariance, original leaf numbering, native versus
fallback Conservative block identity, threaded determinism, production
candidate count, task balance, time, and peak memory.

**Commit:** `Reuse packed trees for fallback cell geometry`, recording both
upstream SHAs in the commit message.
**Landed:** `d194e2d`. The adapter is `CellSpaceRTree`
(`lib/GlobalRegridding/src/conservative.jl:25`); `CellCapTree` is gone, and E2
took `RasterFlatTree` with the `chunktree` bridge. One folded-in action did not
ship with it and is carried as B5.

### Task B5 — window the DGGS hierarchy for a non-chunk-aligned subtree

**Prerequisite:** B4.
**Owns:** `subcursor` for hierarchical DGGS grids, the `DGGSpace` subtree
route, and its cache/benchmark. Scope:
`regrid-notes/2026-08-23-dggs-subcursor-scope.md`, decided by the user.

B4's "use the R-tree adapter only when a native restricted cursor is
unavailable" is vacuous today: `subcursor` has one concrete method, CopernicusDEM
(`src/systems/CopernicusDEM/cursor.jl:412`), and every hierarchical DGGS takes
the `= nothing` default at `src/interface/grid.jl:665`, so a non-chunk-aligned
window packs an R-tree instead of windowing the hierarchy.

Actions:

- Implement `subcursor` for hierarchical DGGS grids over the existing
  `SubtreeIds` O(1) contiguous-run view (`src/engine/partial_grid.jl:10`),
  routed through a `PartialGrid` for the descent clamp.
- Subsume `_chunkcursor` (`src/regridding.jl:254`) entirely — a net deletion.
  It is the third branch of `GR.subtree(::DGGSpace, inds)`
  (`src/regridding.jl:243-251`); the window replaces it.
- Two hard requirements from the scope: the index shift must be
  **bidirectional** — a one-directional shift passes leaf-set and split-weight
  assertions while producing silently wrong weights — and the window must carry
  **cached leaf caps**, because a bare cursor loses the join by 3-4x, worse than
  the fallback it replaces.
- Guard the destination-tiling cliff independently and first. Nothing in
  production reaches the fallback today (all four `subtree` call sites take case
  1 or case 3, with 4.07x headroom in `_defaulttilesizes`,
  `lib/GlobalRegridding/src/lazy.jl:85`), but at `budget = 2^27` every
  destination unit packs a ~419k-cell R-tree on 40 workers: measured 1.213 s
  against 0.501 s for a cached window, at 4.3x the memory. That cliff, not a
  steady-state win, is the reason this is on the list.

Verify leaf-set identity against the R-tree fallback, weight identity in both
shift directions, cached-cap join time, the destination-tiling cliff at
`budget = 2^27`, and threaded determinism.

**Phase 1B gate:** GeometryOps owns reusable cap primitives, CR accepts mixed
extents, and GlobalRegridding has no custom general-purpose fallback cell tree.
**Commit:** `Window the DGGS hierarchy for restricted cell trees`.

## Phase 1C — the `RegridSpace` contract

### Task C1 — consolidate the qualified extension interface

**Prerequisite:** B4.
**Owns:** interface declarations/docs, module public/export lists, DGG extension
methods, and conformance tests.

Actions:

- Gather and document the contracts by responsibility:
  `ncells/getcell/celltree/subtree`,
  `nchunks/ownedindices/chunkextent/chunkindex/candidatechunks!`,
  `chunkranges`, point/chart hooks, and
  `destinationdims/dimsource/_asspace`.
- Mark load-bearing qualified hooks public even when unexported.
- Correct the module comment: DGG already extends `chunkindex`,
  `candidatechunks!`, and `chunkranges`.
- Keep `ownedindices` distinct from rectangular `chunkranges`.
- Do not add an abstract chunk-index type, trait bundle, plan core, or forwarding
  space wrapper.
- Leave `chunktree` documented as a temporary compatibility bridge for E2.

**Phase 1C gate:** every external specialization is named by one coherent
contract and no new interface abstraction exists.
**Commit:** `Consolidate the regridding space interface`.
**Landed:** `5e91b0c`. The grouped contract is
`lib/GlobalRegridding/src/spaces.jl`; the `public` lists are
`lib/GlobalRegridding/src/GlobalRegridding.jl:101-137`. `ownedindices` was then
spelled `cellindices`; `7a92dee` renamed it and left the shim at `spaces.jl:94`.

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

**Status: LANDED EARLY, OUT OF ORDER, AND CLOSED.** Do not reopen it.

It shipped as PR #69 (`ba2bbfa`, merged as `da9e737`) before G1 existed, and so
without either gate this card said to pass first. It landed as a **correctness
fix**, not as the planned performance-neutral cutover: a graph built from
`chunkextents` directly *crossed* the executor's own candidate relation instead
of dominating it, so a refcount derived from it retired source tiles a
subsequent read still demanded (`regrid-notes/2026-08-23-chunk-dag-coverage.md`).
What shipped: destination-major rows from one source `chunkindex` and
independent `candidatechunks!` queries; deterministic sorted rows and the
bidirectional Int32 CSR preserved; the latitude-sorted join deleted.

G1/PR #70 then ran both gates retroactively, on this card's own case matrix.
Verdict in `regrid-notes/2026-08-23-g1-graph-oracles.md` §5:

- **Correctness gate: PASSED.** The indexed relation holds every geometrically
  contributing pair in every oracle-checked case at every radius, and every
  demanded pair on the production pair — where the deleted cap join misses 72.
- **Performance gate: FAILED, AND WAIVED.** The archived 0.0594 s latitude-join
  median this card required was not recovered and cannot be; the two builders do
  different work. On the production pair: **5.7× slower at t8** (0.1219 s vs
  0.0213 s), **6.1× at t4** (0.2289 s vs 0.0376 s), **1.6× the allocations**
  (17.1 MB vs 10.6 MB), and up to 60× per destination on a shallow raster source
  (4320×2160/162 chunks).

The waiver: 0.12 s against a run whose recorded wall time is 8.81 h — about
4 × 10⁻⁶ of it, and 1/90th of the space construction that precedes it — for a
relation a refcount can be derived from (72 → 0 demanded-but-unheld pairs) and
322 fewer edges. The gate is recorded as *failed and waived*, never as passed.
The factor, if a later task wants it back, is in the `RasterGrid` arm's
per-query cursor copy and quadtree descent; `benchmark/chunk_graph_gates.jl`
with `DGG_GRAPH_GATE_CASES=raster-4320-162chunks` reproduces it, and its
`:latjoin` arm keeps the deleted builder runnable while the waiver stands.

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

**Done when:** an invalid graph reuse fails at construction and row views retain
the destination's own chunk numbering for refinement.
**Commit:** `Add reusable dependency graph identities`.
**Landed:** PR #71, stacked on #70. Record:
`regrid-notes/2026-08-23-g3-graph-identity.md`. `DependencyIdentity` is
`lib/GlobalRegridding/src/chunkgraph.jl:129`, `restrict` is `:979`. Two
corrections to this card's premise:

- `chunkindex` is **not** a fresh packed R-tree per column on either shipped
  native space. `chunkindex(::DGGSpace) = space` is a field read and the
  `RasterGrid` one is a cursor; both measure below 0.5 µs in every harness case.
  The R-tree claim holds only on the generic path, which neither shipped native
  space takes. The row view's saving is real ([measured] 26× for one destination
  and 80× for a 4 136-chunk column on the production pair, with 16× less
  allocation), but it comes from the destination caps and the per-row queries.
- `restrict` pays `O(nsourcechunks)` for its source-major transpose, because a
  row view's refcounts are its own. Restricting to a *single* row on a
  26 475-source problem is [measured] slower than re-querying it. Restrict a
  column, not a row.

The identity is a fingerprint, not a proof: equal source chunk caps do not imply
equal source relations, because the relation comes from `chunkindex` rather than
from the caps — the divergence #69 exists to fix. That hole is documented on
`spacestamp` and left open; closing it needs an index-identity hook on the
qualified space interface, which no card owns.

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
- Delete `connectedchunkpairs`. The archive scheduled it for Phase 4 alongside
  `connectedchunks`, but #69 made it a free deletion here: it is now
  line-for-line the same loop as `_chunkgraph`'s rows — one `chunkindex(src)`,
  one `candidatechunks!` per `chunkextents(dst)` entry — and its only caller is
  one lazy test assertion. Move that assertion onto the graph's rows and drop
  the function.
- Preserve zero source reads, zero weight builds, deterministic construction,
  and no network metadata work.

**Phase 2 gate:** one logical plan exposes exactly one validated relation and a
narrow phase cannot be supplied after the plan exists.
**Commit:** `Make chunked plans own dependency graphs`.
**Landed:** PR #72, stacked on #71. Record:
`regrid-notes/2026-08-23-g4-plan-owns-graph.md`. Five of the six actions landed
as written. One correction to the fourth, which G3 predicted and G4 measured:

- **Production's per-column plans cannot take a validated row view, and do not
  need one.** "Never rebuild a one-destination graph per column" is met by
  construction: `dependencies` then defaulted to `nothing`, so `regrid_chunk`'s
  plans built no relation and the default cost [measured] 1.00 µs and 784 B
  against 512 µs / 542 KB for `dependencies = true` and 193 µs / 424 KB for a
  one-row `restrict` — 511× and 193× on 66 175 columns. "Share the global graph
  or validated row views" is **not possible at this destination shape**:
  `regrid_chunk`'s destination is a rooted one-chunk subtree grid, a different
  space from `dagplan`'s 66 175-chunk `PartialGrid`, and a row view still stamps
  the *whole* destination space, so `validate_dependencies` correctly refuses
  it. G4 proved the refusal and left the shape open; closing it needs
  destination-**sub**space re-stamping, which E1 added and production declined.

The global plan does own the relation: `dagplan`
(`scripts/copdem_production.jl:762`) constructs a `GR.ChunkedPlan` with
`dependencies = true` and the narrow phase, and the schedule, refcount cache,
prefetcher and closing validator all read `GR.dependencies(globalplan)` —
asserted at run time by a new `check(...)` in the driver. The relation is
[measured] unchanged: 326 064 edges, equal `DependencyIdentity`, 0.118 s warm.

E1 later reversed this card's `dependencies` default. The current behaviour is
in `_plandependencies` (`lib/GlobalRegridding/src/plans.jl:486`): `nothing` and
`true` both build, `false` is the opt-out.

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
**Landed:** PR #74, stacked on #72. Record:
`regrid-notes/2026-08-23-e1-graph-backed-lazy.md`. All six actions landed;
`_connectedsource!` is `lib/GlobalRegridding/src/lazy.jl:641` and the caps it
costs waves from are read off the relation in `_wavesize` (`lazy.jl:508`).
Three things beyond the card's actions:

- **A chunked plan now owns a relation by default.** It had to: a lazy read *is*
  a read of the rows, so `dependencies = nothing` — meaning "own none" — would
  have made the default plan unable to back a `LazyRegridArray`. `nothing` and
  `true` both build; `false` is the opt-out and is refused by the lazy array, by
  name, in `_lazygraph` (`lazy.jl:275`). That reverses G4's default, and G4's
  measured reason for it (nothing consumed a plan's relation) is exactly what
  this card removed.
- **G4's open action is closed, and production declined it.**
  `subspace_dependencies(g, subspace, destinations)` (`chunkgraph.jl:1047`)
  re-stamps a row view onto a sub-space whose chunk caps it reproduces cap for
  cap — sound because the destination half of the relation is a function of
  those caps alone — and a per-column plan can adopt it. But [measured] on the
  production pair it is only **0.82×** the time and **0.79×** the bytes of
  simply building the one-row relation (408 µs / 425 KB against 500 µs /
  542 KB), because both pay the same `O(nsourcechunks)` transpose and adoption
  still stamps the source space. `regrid_chunk` therefore keeps building its
  own, and E1 added a sampled `graphmisscheck`
  (`scripts/copdem_production.jl:841`) to the driver instead.
- **The card's proof obligation is true only where the source index is the
  generic one.** `graph rows` are not inside the cap join on a native
  hierarchy — the two relations cross, which is why #69 exists — so the raster
  arm of the proof test asserts the eager containment and says so, rather than
  asserting something false.

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
**Landed:** PR #75, stacked on #74. Record:
`regrid-notes/2026-08-23-e2-delete-legacy-discovery.md`. All five actions
landed; what the three spellings of the relation are now is recorded in place at
`lib/GlobalRegridding/src/discovery.jl:133-157`. Three corrections to this card:

- **The `chunktree` verdict is remove, and the audit is in the record.** No
  in-repo or `DiscreteGlobalGrids` space implemented it except the bridge
  itself and three test toys; `DGGSpace` never had one. `chunkextents` is now a
  **required hook with no fallback** (`spaces.jl:125`) — asserted as
  `!hasmethod(chunkextents, Tuple{RegridSpace})` — which is the structural half
  of the gate. `RasterFlatTree` went with it.
- **Two things this card listed were not dead, and were kept.** There is no
  latitude join left in `src/` at all — G2 deleted it — and what remains is the
  gates harness's deliberately archived `:latjoin` arm, whose documented sunset
  condition (G2's waiver retired) has **not** fired, so deleting it would have
  destroyed both the waiver's audit trail and half of E2's own before/after
  evidence. There is likewise no post-plan graph builder left and no redundant
  cap vector: E1 moved the graph's two cap vectors onto the relation on purpose
  and they are aliases to the two `DGGSpace`s' own arrays.
- **`chunkextent` was kept and given the *cheap* method the card asks for.**
  After `connectedchunks` went it had zero callers and zero specializations, but
  it is `public`, documented, and asserted by the qualified-contract testset.
  `RasterGrid` now specializes it at `O(1)` in the chunk count
  (`rastergrid.jl:934`) instead of materializing the whole vector, which is the
  difference between keeping the card's letter and its point.

## Phase 5 — one final weight block representation

### Task W1 — normalize denominator and reference storage

**Prerequisite:** E2.
**Owns:** `WeightCOO`, `WeightBlock`, `CachedBlock`, executor reference helpers,
spill reconstruction, and tests.

Actions:

- Change `WeightCOO` (`lib/GlobalRegridding/src/methods.jl:69`) to
  `denom::Union{Nothing,Vector}` and remove `hasdenom`. Allocate on
  `markdenominated!` or first `adddenom!`.
- Make each `WeightBlock` own its final reusable reference vector. Alias it to
  `denom` for denominated blocks; otherwise compute row sums once.
- Remove `CachedBlock.ref` (`plans.jl:108`); use `entry.block.reference`.
- Reuse the stored reference across eager plan applications.
- Reconstruct, rather than serialize, a redundant reference in `Spilled` and
  count aliased arrays once in `_blockbytes` (`plans.jl:116`).

The `TileWeights` Phase 9 adds will hold ordinary `WeightBlock`s, so this card's
representation is the one it caches per tile. Point blocks carry no denominator;
the empty-`denom` path is theirs and must stay allocation-free.

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

- Add `weightblock(method, dst_space, dst_inds, src_space, src_inds)` (landed
  with S1, `8057af2`).
- Keep `buildweights!` plus `WeightCOO` as the compatible generic fallback.
- Let Conservative adopt its assembled CSC directly with denominators computed
  once; remove CSC-to-COO-to-CSC conversion.
- Make eager and chunked builders call the same seam; remove the special
  `wholeblock(::Conservative)` path
  (`lib/GlobalRegridding/src/conservative.jl:371`), leaving the generic
  `wholeblock` (`api.jl:236`).
- Preserve wrapper/instrumentation dispatch through an explicit forwarding
  rule or trait.

The new seam is per chunk pair, which is the area-method unit. It does not
serve a point method, whose unit is a destination tile (S2); the two seams
coexist and are selected by `outputsampling` (S1).

Verify sparse structure/value identity, partition invariance, denominators,
empty sides, third-party emitters, and the committed peak-RSS benchmark.

**Phase 5 gate:** one block-builder path and one final block representation,
with no measured eager regression and the chunked round-trip memory removed.
Met: interleaved same-session runs put eager production time within 0.9 % on
medians against a 9.5 % within-state spread, and the chunked route's allocations
fall by 812,885,051 bytes on `conservative_roundtrip_baseline.jl` — Phase 0's
810,045,520-byte round trip — and by 6.6 % on `conservative_block_baseline.jl`,
after which the chunked route allocates what the eager route allocates.

**Commit:** `Unify final weight block construction`.

## Phase 6 — generic DGG output and API forwarding

### Task O1 — remove DGG-specific output application

**Prerequisite:** W2.
**Owns:** `DGGSpace.destinationdims`, eager/lazy wrappers, DGG regrid methods,
and cross-system tests.

Actions:

- Define the generic DGG destination dimension as
  `(Cells(CellLookup(space.grid)),)`.
- Remove `_DirectToDGG` and `_ChunkedToDGG` (`src/regridding.jl:306-308`), the
  specialized `GR.regrid` methods, and `_ascube` (`src/regridding.jl:319`).
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

- Replace `regridgrid` plus `RegridTarget` (`src/regridding.jl:274-284`) with
  direct `DGGSpace` constructors or `_asspace` methods for each target spelling.
- Retain the meaningful bare-system error and source-dependent level choice.
- Make `plan_regrid` the single executable owner of keyword defaults and
  validation; forward `kwargs...` from `regrid` and `regrid!`.
- Preserve explicit user-facing keyword documentation and error quality.

**Phase 6 gate:** no DGG-specific plan application or intermediate conversion
generic remains, and API defaults occur once.
Met: this package defines no `regrid`, `regrid!` or `plan_regrid` method and no
`_ascube`, `regridgrid` or `RegridTarget`; each target spelling is one
`_asspace` method; and every keyword's default is declared on `plan_regrid`'s
one method, with `regrid` and `regrid!` declaring only `kwargs...` and the lazy
budget's value named once as `DEFAULT_BUDGET`.

**Commit:** `Simplify regridding target and keyword resolution`.

## Phase 7 — prepared destination geometry

This phase's card is `D1`; `P1` in this plan always names the barycentric
plan's card.

### Task D1 — replace the destination cache wrapper stack

**Prerequisite:** O2.
**Owns:** `TileCells`, `CachedCellTree`, destination `CellMemo`, block
construction, cache tests, and Conservative benchmark.

Actions:

- Introduce one private prepared destination tile holding its index set and its
  chunk-local map, one restricted tree, and optional prebuilt polygons.
- Share it across all source blocks for the executor tile.
- Let `BlockAreaOperator` (`lib/GlobalRegridding/src/conservative.jl:294`)
  obtain destination polygons directly from it.
- Remove `TileCells <: RegridSpace` (`conservative.jl:126`) and its forwarding
  methods, `CachedCellTree` (`:209`), and destination `CellMemo` (`:258`) when
  polygons are prepared.
- Retain a task-local memo for unprepared source geometry and oversized
  destinations.
- Remove stale raster/DGG extent caches only when their individual benchmark
  gate earns deletion; do not merge polygon and extent caches by name.

The prepared tile is destination *polygon* geometry, which a point method does
not use: a point tile needs its destination sample sites, not its cells. Keep
the sample-site path off this wrapper.

Verify one tree build per tile, cached/uncached identity, concurrent safety,
cache threshold, allocations, time, and peak memory.

**Phase 7 gate:** destination geometry follows one preparation/cache path with
no production regression.
**Commit:** `Consolidate destination geometry preparation`.

## Phase 8 — private lazy state

### Task L1 — give `SourceHold` one entry dictionary

**Prerequisite:** D1.
**Owns:** source-hold entries, eviction, residency statistics, and lazy tests.

Actions:

- Replace `SourceHold`'s parallel `held` and `used` dictionaries
  (`lib/GlobalRegridding/src/lazy.jl:144`) with entries containing `value`,
  `bytes`, and `used`.
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
- Rename `isdiskbacked` (`lib/GlobalRegridding/src/api.jl:4`) to reflect that it
  tests declared chunking. Preserve the separate true-disk predicate used for
  DiskArrays reads (`lazy.jl:827`).

Use a counting fixture to verify one `eachchunk` interpretation and test
regular, irregular, absent, malformed/non-grid chunk descriptions, output
tiling, read groups, defaults, and zero source reads.

**Phase 8 gate:** lazy source metadata and residency each have one private
representation, with unchanged behavior and no new framework.
**Commit:** `Normalize lazy DiskArrays chunk metadata`.

## Phase 9 — point-method admission

The barycentric plan's own cards are P0-P6. This phase owns only the shared code
they must enter through: the method seam in `methods.jl`/`plans.jl`, the build
and cache unit in `plans.jl`, and the tile loop in `lazy.jl`. Nothing here
changes a conservative weight, key, read or value; every card names where the
conservative branch stays.

Use the barycentric plan's names exactly: `weightsat!(row, sampler, p)`,
`WeightRow`, `sampler(method, space) -> Sampler` holding a `cellfield` of sample
sites (`src/engine/cellfield.jl:61`), `dualcellat(sampler, p) -> DualCell`,
basis kinds `Linear`/`Bilinear`/`MeanValue`, and `TileWeights(sourcechunks,
blocks)`.

### Task S1 — split weight construction on output sampling

**Prerequisite:** E2.
**Owns:** the point/area split in `plans.jl` and `methods.jl`, and the tests
that pin it. No kernel and no weight changes.

Actions:

- Select the build path with `outputsampling(method) isa DD.Lookups.Points`, the
  trait every method already carries (`lib/GlobalRegridding/src/methods.jl:55-57`:
  `NearestCell` and `BilinearPoint` report `Points()`, everything else
  `Intervals(Center())`). Never dispatch on a concrete method type.
- `Conservative` and every method reporting `Intervals` keep `buildweights!`
  (`methods.jl:131`) with `WeightCOO` and one `WeightBlock` per chunk pair —
  `blockfor`/`buildblock` (`plans.jl:521-548`) unchanged.
- Keep `buildweights!` as the generic fallback for a `Points` method that
  supplies no `sampler`, so a third-party point method and the current
  `NearestCell`/`BilinearPoint` builders (`interpolation.jl`) keep working
  byte-identically. Migrating either of them is the barycentric plan's decision,
  not this card's.
- Add no new trait and no method-type union.

Verify eager and chunked values for `Conservative`, `NearestCell` and
`BilinearPoint` unchanged; a `Points` method with no sampler still routes to
`buildweights!`; dispatch is on the trait in a test that defines its own method.

**Done when:** the shared builders name a point path by trait alone, and no
conservative call site moved.
**Commit:** `Split weight construction on output sampling`.

### Task S2 — build point weights per destination tile

**Prerequisites:** S1, and the barycentric plan's P1 (`WeightRow`, `sampler`,
`weightsat!`).
**Owns:** `TileWeights`, the point branch of `blockfor`/`buildblock`
(`plans.jl:521-548`), block storage keys in `PerChunk` (`plans.jl:133`) and
`Spilled` (`plans.jl:240`), and the tile loop's build step in
`_readdestination!` (`lib/GlobalRegridding/src/lazy.jl:424`).

Actions:

- Add `TileWeights(sourcechunks::Vector{Int}, blocks::Vector{WeightBlock})`:
  the exact sorted source chunk numbers and their blocks in the same order.
  Rows are the tile's chunk-local destination indices and columns are the source
  chunk's chunk-local indices — the convention `WeightCOO` already documents
  (`methods.jl:61-68`).
- Build one tile in one pass over its destination sample sites: `weightsat!` per
  point, then partition each nonzero entry by `chunkat(src_space, i)`
  (`spaces.jl:154`) and convert its local index to the chunk-local one by
  `ownedindices` arithmetic when the chunk is contiguous, and through
  `indexmap`/`localindex` (`lib/GlobalRegridding/src/shared.jl:18-37`)
  otherwise. Combine duplicate source cells deterministically and finalize
  nonempty blocks in ascending source-chunk order.
- Make the tile the cache and locking unit for a point method: concurrent
  requests for one tile build its stencils once. Its whole footprint counts
  against `weightbudget` (`plans.jl:94`), and a spilled tile stores its manifest
  beside its blocks so an evicted dependency list is never recovered by rerunning
  the destination pass.
- Keep temporary memory bounded by the tile's nonzeros plus named worker
  scratch, never by destination cells times candidate chunks.
- Leave `Conservative`'s `(destination tile, source chunk)` key, its `PerChunk`
  eviction order, and its `Spilled` files exactly as they are.

Verify one stencil query per requested destination cell; identical values eager,
lazy and under a different source chunking; ascending source-chunk application;
budget and spill round trips for a point tile; conservative cache statistics
unchanged.

**Done when:** a point method's atomic build unit is one destination tile with an
exact sorted source-chunk manifest, and the conservative unit is still one chunk
pair. This is half of the barycentric plan's P3 gate.
**Commit:** `Build point weights per destination tile`.

### Task S3 — read exactly a point tile's source chunks

**Prerequisite:** S2.
**Owns:** `_connectedsource!` (`lazy.jl:641`), `_lazygraph` (`lazy.jl:275`),
`_plandependencies` (`plans.jl:486`), `_builddependencies`
(`chunkgraph.jl:740`), and the read-count and residency tests.

Actions:

- For a tile that has a `TileWeights`, take its source chunks from
  `TileWeights.sourcechunks`. The graph row is a superset; the manifest is
  exact. Every other tile — conservative, or a `Points` method still on
  `buildweights!` — keeps taking `sourcesof(dependencies(plan), d)`, or the
  ascending union of its rows for a derived tile.
- Apply `knownempty` filtering after selection on both paths, as
  `_connectedsource!` already does: data-dependent filtering may drop a chunk
  the manifest holds and may never add one it does not.
- Buffer nothing: what a tile with `TileWeights` reads is `sourcechunks`
  exactly, neither an intersection with the relation's row nor the row. A
  sampler-bearing method still declares `supportradius` (`methods.jl:145`) — a
  true bound on its stencils, not an approximation of them, because cap overlap
  alone is no superset of a stencil reaching sample sites whose own cells the
  destination never touches — and a manifest naming a chunk outside the tile's
  rows is refused with an `ArgumentError` naming `supportradius` and the method.
  No card may dilate a destination cap or join neighbouring chunks to
  approximate a point stencil's reach.
- State what a point plan's relation is for. `LazyRegridArray` requires one
  (`_lazygraph`) and `_wavesize` (`lazy.jl:508`) costs waves from the caps it
  carries, so a point plan still owns one — for ordering, wave costing,
  refcounts and prefetch only. Document it as a superset that never decides a
  read.
- Assert standing law 7 for point tiles: the chunks read equal the sorted union
  owning the tile's nonzero stencil entries, and that union is inside the tile's
  graph rows.

Verify a per-tile read-count assertion; no chunk outside the manifest reaches
`_sourcefor!` (`lazy.jl:720`); conservative read counts, residency and values
unchanged; a rechunked source gives the same values and a manifest that tracks
the new chunking.

**Phase 9 gate:** a point tile reads exactly the chunks its stencils name, with
no cap superset and no support dilation, and every conservative read count,
weight and value is unchanged.
**Commit:** `Read exact source chunks for point tiles`.

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
| Point-path changes | one stencil query per destination, per-tile read-count assertion, eager/lazy/rechunked equivalence, unchanged conservative results |

Performance reports must state thread count, shape/parallelism mode, input
sizes, warmup, statistic, allocations, and peak-memory method. The detailed
archive lists the current reference measurements and exclusions, including the
hybrid-source production columns that are invalid byte-identity references.

## Definition of done

- Raster geometry uses GeometryOps transformations, task-local Proj.jl clones,
  and one CR curvilinear-grid adapter over DiskArrays-owned ranges.
- GeometryOps owns public spherical extent primitives; CR accepts mixed extent
  trees; arbitrary cell subsets use one `FlexibleRTree` adapter, taken only when
  the source grid cannot window its own hierarchy.
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
- The shared path admits a point method on its own terms: dispatch by
  `outputsampling`, one `TileWeights` per destination tile, and reads that are
  that tile's exact source-chunk manifest rather than a cap superset — the
  relation's rows remain a superset that decides no read of such a tile, kept
  one by the method's declared `supportradius` and never by dilating a cap or
  buffering an adjacent chunk — with conservative behaviour unchanged, and the
  kernels themselves owned by the barycentric plan.
- All correctness suites and representative production/small-chunk performance
  gates pass, and every phase-ending commit records its evidence.
