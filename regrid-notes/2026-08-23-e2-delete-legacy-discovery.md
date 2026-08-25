# E2 — the duplicate discovery paths and the `chunktree` bridge are gone

- Date: 2026-08-23
- Card: `regrid-notes/2026-08-21-regridding-simplification-plan.md`,
  "Phase 4 — delete legacy discovery", Task E2
- Commit: `Remove duplicate chunk discovery paths` (`526c057`), on top of
  `0e3fb70`.

Phase 4's gate: **one query implementation defines graph edges, and neither the
executor nor the interface translates that relation back through a compatibility
tree.** Both halves are tests, in
`test_chunkgraph.jl`, testset `one query implementation defines every edge`.

This card is almost entirely removal, so §1 states what was deleted and the
audit that licensed each deletion, and §2 states what was **kept** and why —
including the two things the card lists that turned out not to be dead.

## 1. Deleted, and the audit for each

### 1.1 `connectedchunks` / `connectedchunks!`

`lib/GlobalRegridding/src/discovery.jl`. Four methods: the space form, the two
`(cap, src_space)` / `(cap, srcindex)` in-place forms, and the
`(caps::Vector, srcindex)` union form.

**Audit before deleting.** Neither name was exported, `public`, or documented
outside its own docstring. Repo-wide, after E1 there were **no callers in
`src/`, `ext/`, `scripts/` or `benchmark/`** — E1's record says so and a
repo-wide grep confirmed it. The nine remaining call sites were all in
`lib/GlobalRegridding/test/test_lazy.jl`.

The `(cap, srcindex)` and `(caps, srcindex)` forms were the card's
"compatibility source-index state": they let a caller carry a *prebuilt source
index* around outside any relation and re-query it, which is exactly the state
E1 removed from `LazyRegridArray` and which a plan's relation now owns. They
went with the function.

### 1.2 The `chunktree` bridge — the verdict is **remove**

The card says remove it *only if no use remains*. The audit, in full:

| site | what it was | verdict |
|---|---|---|
| `GlobalRegridding.jl:65` | `export chunktree` | removed |
| `spaces.jl:143-153` | the declaration + docstring, marked "temporary compatibility bridge … E2 decides" | removed |
| `discovery.jl:131-155` | `chunkextents(::RegridSpace)` fallback + `_collectextents!`, the **only** consumer in `src/` | removed |
| `rastergrid.jl:889-926, 968-971` | `RasterFlatTree` + 15 trait methods, and `chunktree(::RasterGrid)` that built it | removed |
| `conservative.jl:156` | `chunktree(tc::TileCells)` forwarder | replaced by `chunkextents`/`chunkextent`/`chunkindex` forwarders |
| `toyspaces.jl:308, 354` | `ToyLonLatSpace` and `CountingSpace` methods — the toys' *only* way of reporting chunk caps | became `chunkextents` methods, same caps |
| `test_chunkgraph.jl:46` | `G4ProbeSpace`'s counter seam | moved up one level to `chunkextents` |
| 8 test sites | `chunktree(x).caps` | `GR.chunkextents(x)` |

**Downstream, checked, not assumed.** `DiscreteGlobalGrids` — the one in-repo
downstream — defines **no** `GR.chunktree` method at all: `DGGSpace` specializes
`chunkextents`, `chunkindex` and `candidatechunks!` directly
(`src/regridding.jl:159-221`). `src/cap_cached_tree.jl`'s `_cachedchunktree` and
`test/systems/crosssystem/regrid.jl:520`'s use of it are a **private DGG name**
for caching a `HierarchicalGridCursor`'s caps, unrelated to `GR.chunktree`; both
are untouched. `regrid.jl:84-86` binds a **local variable** called `chunktree`
to a `GR.subtree` result; also untouched. `ext/` has no reference. `docs/` does
not build GlobalRegridding at all, so no docs page names it.

So the only thing that still *needed* `chunktree` was `chunktree` itself: a
space packed its caps into a flat tree so that the generic `chunkextents` could
walk the tree and collect the same caps straight back out. `chunkextents` is now
a **required hook with no fallback**, which is asserted structurally:
`!hasmethod(GR.chunkextents, Tuple{RegridSpace})`.

`RasterFlatTree` went with it, closing the note at `rastergrid.jl:969-971` that
the bridge was its only remaining role. `_WHOLE_SPHERE`, which its constructor
used, has three other callers and stays.

**Behaviour change, stated:** an out-of-repo space that implemented only
`chunktree` no longer works. That is the deletion the card asked for, it is a
`0.x` unregistered workspace sub-package, and nothing in the repo or in
`DiscreteGlobalGrids` was such a space.

### 1.3 What the deletion did *not* touch, having looked

The caps a `RasterGrid` reports are computed by the identical
`Trees.cell_range_extent` call as before; `chunktree(::RasterGrid)` was built
*from* `chunkextents(::RasterGrid)`, so collecting them back out was a no-op and
removing the round trip cannot move a value. The toys' caps are the identical
expression, now returned as a `Vector{Cap}` instead of wrapped in a
`ToyCapTree`. `ToyCapTree` itself survives — it is still `celltree`.

## 2. Kept, and why — the card lists things that are not dead

### 2.1 `chunkextents` — kept, with many real consumers

Named, not assumed: `spacestamp` (`chunkgraph.jl:96`), `_builddependencies`
(`chunkgraph.jl:747-748`, both sides, which is where the relation's own cap
vectors come from), `subspace_dependencies` (`chunkgraph.jl:1050`), the generic
`chunkindex` (`discovery.jl:40`), `chunkextent`, and in `scripts/`:
`copdem_production.jl:774` (the lon/lat boxes the `refine` narrow phase is built
from), `copdem_download_count.jl:169`, `chunk_dag_poc.jl:229-230`. Plus the
shared oracles and the gates harness. This is the card's "real
diagnostics/planning consumers", and there are many.

### 2.2 `chunkextent` — kept, and given the cheap method the card asks for

**Audit:** after `connectedchunks` went, `chunkextent` had **zero callers in the
repo** and **zero specializations anywhere**. It is `public`, documented, and
asserted by the qualified-contract testset in
`lib/GlobalRegridding/test/runtests.jl:17`.

The card says to keep *a cheap* `chunkextent`. The generic fallback
`chunkextents(space)[chunk]` was not cheap — it materialized every chunk's cap
to take one — so keeping the name unchanged would have kept the letter and not
the point. `RasterGrid` now specializes it (`rastergrid.jl`): one `chunkbox` and
one `cell_range_extent`, `O(1)` in the number of chunks instead of `O(nchunks)`.
On a `DGGSpace` the fallback is already `O(1)` — `chunkextents` is the field read
`space.caps`. Both `TileCells` forwarders are new.

The docstring now also points a caller at `destinationextent`/`sourceextent` on
the relation, which is the right answer whenever a relation exists and is the
reason the singular form has no hot caller.

### 2.3 The `:latjoin` arm in `benchmark/chunk_graph_gates.jl` — kept

The card says "remove any latitude join left after G2". **In `src/` there is
none**: G2/PR #69 deleted the production builder and a grep of
`lib/GlobalRegridding/src`, `src` and `ext` finds no latitude join. What remains
is the harness's deliberately archived *reimplementation*, which carries an
explicit SUNSET CONDITION (`chunk_graph_gates.jl:38-49`): delete it when G2's
waived performance gate is retired, or when the production builder stops being
comparable to a flat cap join. **Neither has happened** — the waiver in
`2026-08-23-g1-graph-oracles.md` §5 still stands and is the thing the arm exists
to keep auditable. Deleting it would also have destroyed half of this card's own
before/after evidence. Kept, unchanged, and its 13 rows are part of §4's
comparison.

### 2.4 The generic `candidatechunks!(out, index, dstcap)` STI fallback — kept

`discovery.jl:120-129`, "compatibility fallback for existing extension trees".
After 1.2 nothing in the repo reaches it: every `chunkindex` in the repo returns
an `RTree`, a `TopDownQuadtreeCursor`, an `EmptyChunkIndex`, a `DGGSpace`, or a
`TileCells` forward of one of those. It is nonetheless kept, because it is not a
second relation — it is one *method* of the one query function, dispatching on
the index type, and it is what makes the documented `chunkindex` promise
("structured spaces may return any native hierarchy") true. The Phase 4 gate is
about one query *implementation* defining the edges, and it does: `_chunkgraph`
calls `candidatechunks!` once per destination cap and nothing else, which §3
asserts by counting.

### 2.5 Nothing to do: post-plan graph builder, redundant cap vectors

- **Post-plan graph builder.** G4 already removed `refine` from post-plan
  builders, and `api.jl:183` documents that `chunk_dependency_graph` has no
  `plan` method. `chunk_dependency_graph(dst, src; radius, refine)` is the only
  builder and `dependencies(plan)` is a non-building accessor. Nothing left.
- **Redundant cap vectors.** E1 removed `LazyRegridArray`'s `dstcaps`/`srccaps`.
  The two vectors on `ChunkDependencyGraph` are the relation's *own inputs*, are
  aliases to the two `DGGSpace`s' `caps` fields in production, and are what
  E1 deliberately moved there; `lazy.jl:543-544` binds local names to them by
  reference and copies nothing. Not redundant, not removed.

## 3. The Phase 4 gate, as tests

`lib/GlobalRegridding/test/test_chunkgraph.jl`, testset **`one query
implementation defines every edge`** — 19 assertions.

**Structurally**, that the duplicate paths are gone from the module rather than
merely unexported:

```julia
for name in (:connectedchunks, :connectedchunks!, :connectedchunkpairs,
             :chunktree, :RasterFlatTree, :_collectextents!)
    @test !isdefined(GR, name)
end
@test :chunktree ∉ names(GlobalRegridding)
@test !Base.ispublic(GR, :chunktree)
@test !hasmethod(GR.chunkextents, Tuple{RegridSpace})   # no tree to fall back to
```

The `hasmethod` line is the "no compatibility tree" half: there is no method on
the abstract space type, so no space can be answered by collecting a tree.

**Behaviourally**, that one query defines every edge. `E2QuerySpace` wraps a
space with an `E2QueryIndex` that counts `candidatechunks!` calls:

```julia
probe = E2QuerySpace(src)
g = planned_dependencies(dst, probe)
@test probe.index.queries == nchunks(dst)          # one query per destination cap
@test graph_pairs(g) == demanded_pairs(dst, src)   # and the answers ARE the relation
```

then, that a lazy read over that plan adds **no** query of its own — the E1 gate
re-asserted at the query seam rather than at `chunkextents`.

Plus that the cheap `chunkextent` agrees with the vector it is cheaper than, on
the misaligned raster.

The oracles are `graph_pairs`/`demanded_pairs` from
`lib/GlobalRegridding/test/graphoracles.jl`, **reused, not copied**.

### The rewritten oracle tests

`test_lazy.jl`'s `discovery` testset now reads every claim off
`sourcesof(dependencies(plan), d)`, through one helper:

```julia
t7_sources(plan::ChunkedPlan, d::Integer) =
    Int.(GR.sourcesof(GR.dependencies(plan), Int(d)))
```

- the per-chunk agreement with the pairwise cap reference is `t7_sources(plan, c)`;
- dilation now reaches the relation **through the plan's method**, which since G4
  is the only way it can, and the testset asserts `dependency_radius == 0.5`
  before asserting the containment;
- `chunkextents(src) == chunktree(src).caps` became
  `sourceextents(graph) == chunkextents(src)` and the destination equivalent —
  a stronger claim, that the relation carries the very vectors it was built from;
- `sort(pairs) == expected` became `graph_pairs(graph) == demanded_pairs(...)`,
  the shared oracle, which is the same equality without a third spelling of the
  loop.

The four wave/dilation testsets that used `connectedchunks` to obtain a tile's
sources now take them off the plan they then execute — so those assertions are
about the rows the executor will actually use. The wave-failure testset builds
its failing method *from* a chunk number, so it reads the rows off a radius-0
plan first and then asserts the plan under test holds the same row.

`t7_pairwise` keeps its independent pairwise-cap definition; only its input
moved from `chunktree(space).caps` to `chunkextents(space)`.

## 4. `benchmark/chunk_graph_gates.jl` — the relation is identical

Run at the branch point (`0e3fb70`, before any edit) and on this branch. 13
cases including the production `copdem90-igeo7-l12` pair from a **local** tile
list, both arms, **26 ndjson rows each**, `-t 8 --gcthreads=4`, 5 samples,
Julia 1.12.6, GeometryOps `697c7cc`, ConservativeRegridding `e8fc67f`.

Both runs print:

```
cases run: 13.  Oracle-checked: 9.  Oracle skipped: 4.
verdict: PASS on 9 oracle-checked case(s); 4 case(s) unchecked
```

A real geometric verdict on nine cases, not a `NOT CHECKED` run.

**Every relation field is identical on all 26 rows**, compared field by field:
`edges`, `demanded_pairs`, `demand_missing`, `oracle_pairs`, `oracle_missing`,
`only_here`, `missing_here`, `destination_chunks`, `source_chunks`, `radius`,
`identity_bytes`, `graph_allocated_bytes`, `graph_summarysize_bytes`. That
includes the production pair's 326 064 indexed edges, 326 386 latjoin edges and
the 72 / 394 crossing.

`graph_allocated_bytes` is **byte-identical on every one of the 26 rows** —
17 123 800 B indexed / 10 614 384 B latjoin on the production pair — which is
the direct evidence that removing the bridge removed no allocation and added
none: `chunkextents(::RasterGrid)` and `chunkextents(::DGGSpace)` were already
what the builder called, and `chunktree` was never on the build path.

**Timings.** Production `:indexed` +1.8 % (0.11581 → 0.11788 s), within the
±2 % band G3, G4 and E1 all recorded on unchanged code. Sub-millisecond cases
swing between −28.5 % and +19.6 % in both directions on both arms, the same
spread as before, on a shared 64-core box with other work on it. Nothing on the
build path changed at all.

## 5. Values and residency

`scratchpad/residency.jl` — E1's script, reused verbatim, not committed (it is a
measurement, not a suite). One 720×360 tiled raster onto IGeo7 L4 with level-2
chunks, 24 012 destination cells, 492 chunks ← 72 chunks, a 4 MiB budget,
`-t 4 --gcthreads=4`. Measured at `0e3fb70` and again after `526c057`.

| | before (`0e3fb70`) | after |
|---|---|---|
| `LazyStats` | loads 328, hits 3654, skipped 0, dropped 0, **peak 2 073 600 B** | **identical** |
| source `readblock!` calls | 400 | **400** |
| repeated whole read identical | true | **true** |
| all seven slices identical | true | **true** |
| differently tiled plan identical | true | **true** |
| max abs(eager − lazy) | **2.44249e-15** | **2.44249e-15** |
| live at construction: plan | 82 072 B | **82 072 B** |
| live at construction: array (minus source and plan) | 43 376 B | **43 376 B** |

Every number is identical, including the two `summarysize` figures — the plan
and the array hold exactly the same bytes as before, because nothing removed was
ever retained by either. The 2.44e-15 eager/lazy difference is the pre-existing
accumulation-order difference E1 measured and is **bit-identical**, which is the
evidence that E2 moved no value.

Peak process RSS was 919.4 MiB before and 913.2 MiB after; the whole-process
deltas within a run were 0.0 and 55.7 MiB in the two runs, which is compilation
noise at this size, not a signal, and is not reported as one. The number the
budget actually bounds — `LazyStats.peak` — is byte-identical.

## 6. Suites

`-t 8 --gcthreads=4`, Julia 1.12.6, except the GlobalRegridding suite at `-t 4`.
The "before" column was re-measured at `0e3fb70` and reproduces E1's reported
numbers exactly.

| suite | before (`0e3fb70`) | after | delta |
|---|---|---|---|
| `lib/GlobalRegridding/test` | 3 999 / 1 broken / 0 fail | **4 021 / 1 / 0** | **+22 pass** |
| `test/systems/crosssystem/regrid.jl` | 236 / 0 | **236 / 0** | none |
| `test/systems/crosssystem/regrid_acceptance.jl` | 22 / 0 | **22 / 0** | none |
| `test/scripts/copdem_policy.jl` | 89 / 0 | **89 / 0** | none |
| `test/scripts/copdem_source_mode.jl` | 7 / 0 | **7 / 0** | none |
| `test/systems/CopernicusDEM/runtests.jl` | 16 258 / 3 broken / 0 fail | **16 258 / 3 / 0** | none |

**No count dropped.** Every assertion that named a deleted function was
*rewritten* onto the relation rather than removed, so nothing was deleted along
with the code it covered — which is why the four suites outside GlobalRegridding
are unchanged to the assertion.

**+22 in GlobalRegridding**, exactly:

| where | net | what |
|---|---:|---|
| `one query implementation defines every edge` (new) | +19 | 6 `!isdefined`, 2 export/public, 3 `hasmethod`, 2 `chunkextent` specialization, 3 query-count + oracle, 2 lazy-read-adds-none, 1 cheap-`chunkextent` agreement |
| `discovery` | +2 | `dependency_radius == 0.5` for the dilated plan (+1); `chunkextents == chunktree(...).caps` (−1) became `sourceextents` and `destinationextents` (+2) |
| `a wave that loses a task still waits for the rest` | +1 | the plan under test really does hold the rows the failing method was built from |
| `discovery` loop, `dilation`, `test_conservative`, `test_rastergrid`, `runtests` | 0 | same assertions, `chunkextents`/`t7_sources` instead of `chunktree(...).caps`/`connectedchunks` |

The one pre-existing broken GlobalRegridding test is still broken; the three
pre-existing broken CopernicusDEM conformance skips are unchanged. No failures
or errors in any suite.

## 7. Residual uncertainty

- **`chunkextent` has no in-repo caller.** It is kept because the card names it,
  because it is `public`, documented and asserted by the qualified-contract
  testset, and because it now has a genuinely cheap `RasterGrid` method — but
  its only exercise is that testset and the new gate's agreement check. If a
  future card wants the public surface smaller, this is the name to reconsider,
  and the honest reason to keep it today is API stability rather than use.
- **The generic STI `candidatechunks!` fallback is now unreached in-repo.** §2.4
  argues it is one method of the one query rather than a second relation. That
  argument is a judgement about what the gate means, not a measurement; a
  reviewer who reads the gate as "one *method*" would delete it.
- **The `:latjoin` arm outlives this card.** It is archived code kept for a
  standing waiver. Every task after this one that runs the gates harness pays
  its build time and has to read past it.
- **No out-of-repo consumer was audited, because none is known.** The
  `chunktree` verdict rests on the repo plus `DiscreteGlobalGrids`, which is the
  only downstream in this workspace. GlobalRegridding is an unregistered `0.x`
  path dependency, so this is a complete audit of everything that can resolve
  it, but it is not a search of a registry.
- **The production driver was not run end to end**, the same as under G4 and E1.
  `scripts/copdem_production.jl` is untouched by this card; the production pair
  is covered by the gates harness's `copdem90-igeo7-l12` case, whose relation is
  identical field for field.
- **Sub-millisecond gate timings are noisy on this box.** The ±28 % swings in
  §4 are shared-machine noise, in both directions on both arms, on a build path
  that did not change and whose allocations are byte-identical. The production
  case's +1.8 % is the only timing worth quoting and it is inside the band the
  three previous cards recorded on unchanged code.
