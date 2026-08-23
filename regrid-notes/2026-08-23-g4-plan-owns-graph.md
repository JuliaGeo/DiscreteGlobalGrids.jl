# G4 — `ChunkedPlan` is the sole owner of the dependency graph

- Date: 2026-08-23
- Card: `regrid-notes/2026-08-21-regridding-simplification-plan.md`, Phase 2,
  Task G4. Last card of Phase 2.
- Branch: `claude/g4-plan-owns-graph`, cut from `claude/g3-graph-identity`
  @ `f9f268a` (PR #71's head). PR **#72**, stacked into
  `claude/g3-graph-identity`; #70 and #71 stay open.
- Commit: `Make chunked plans own dependency graphs` (`71fefab`).
- Reproducers: `benchmark/chunk_graph_gates.jl` (the G1 gate, unmodified) and
  `benchmark/plan_dependency_ownership.jl` (new, this card's measurement). Both
  take `DGG_COPDEM_TILELIST` pointing at a **local** tile list and neither ever
  downloads.

The card asked for six things: make `refine` a keyword of lazy `plan_regrid` and
nowhere else; construct or validate dependencies once in the plan and expose
`dependencies(plan)` as a non-building accessor; remove `refine` from post-plan
graph builders; make production's global and per-column plans share the global
graph or validated row views and never rebuild a one-destination graph per
column; delete `connectedchunkpairs`; and preserve zero source reads, zero
weight builds, deterministic construction and no network metadata work.

Five landed as written. The fourth landed as **half the requirement met by
construction and the other half proved impossible at production's shape**; §4
says exactly what and why, with the measurement.

---

## 1. The ownership design

`ChunkedPlan` gained one field and the module gained one accessor.

```julia
struct ChunkedPlan{...,G<:Union{Nothing,ChunkDependencyGraph}} <: AbstractRegriddingPlan
    ...
    dependencies::G
end

dependencies(plan::ChunkedPlan) = plan.dependencies
dependencies(::DirectPlan)      = nothing
```

`chunkgraph.jl` now precedes `plans.jl` in the include order, so the plan's
field type can name the graph's type rather than be an unconstrained parameter.

The relation is chosen **once**, at construction, by one keyword on lazy
`plan_regrid` (and on the `ChunkedPlan` keyword constructor it forwards to):

| `dependencies` | what the plan holds | what construction costs |
|---|---|---|
| `nothing` (default) | none | **nothing at all** — a dedicated `_plandependencies(::Nothing,::Nothing,::Nothing,…)` method returns `nothing` without even calling `support_radius` |
| `true` | a relation built here, at `support_radius(method, src_space)` | one `chunk_dependency_graph` |
| a `ChunkDependencyGraph` | that object, by reference, after `validate_dependencies` certifies it | two `spacestamp` calls |
| `false` | none, explicitly | nothing; and it *rejects* `refine`/`narrow` rather than dropping them |

`refine` (the narrow phase) and `narrow` (its Symbol tag) are accepted only
alongside a construction; passing either with `dependencies === nothing` implies
`true`. `refine` with a supplied graph is an error, and the message says why: a
narrow phase applies while a relation is built, so reapplying it afterwards
would prove nothing. Name the phase the supplied graph already carries with
`narrow` instead.

The default is `nothing`, not `true`, and that is a deliberate design decision
rather than caution. Nothing in tree *consumes* a plan's relation yet — the lazy
executor still discovers a tile's sources by querying the source index, which is
E1's card to change — so building one by default would add cost to every lazy
plan for a benefit that has not arrived, and would build 66 175 of them in the
production run. §4 measures what that would have cost.

### What moved out of the graph builder

- `chunk_dependency_graph(plan::ChunkedPlan; refine, narrow)` is **deleted**.
  `dependencies(plan)` replaces it, and there is now no method of
  `chunk_dependency_graph` that accepts a plan at all.
- `chunk_dependency_graph(dst_space, src_space; radius)` **lost `refine` and
  `narrow`**. The narrow phase moved to a private `_builddependencies(dst, src,
  radius, refine, narrow)` that only `plans.jl` calls, once, during
  construction. The public builder is now purely "the relation these two spaces
  have at this radius", with no way to thin it.

## 2. What `dependencies(plan)` guarantees

1. **It builds nothing.** It is a field read. It queries neither space, touches
   no index, and reads no data. Asserted with a `G4ProbeSpace` that counts every
   `chunktree` call — the single funnel through which both `chunkextents` (the
   destination caps *and* the identity stamp) and the generic `chunkindex` pass:
   the counter is `0` after constructing a default plan, and is unchanged by
   five subsequent `dependencies` calls on a plan that owns a relation.
2. **It is the same object every time.** `dependencies(p) === dependencies(p)`
   and `=== p.dependencies`. "The plan's relation" is therefore a phrase that
   denotes one object, which is what lets a scheduler, a refcount cache, a
   prefetcher and a validator agree without any of them re-deriving it.
3. **At most one, ever.** A plan holds zero or one relation and there is no
   second route to a different one: no plan method on the builder, and the
   struct is immutable so the field cannot be swapped afterwards.
4. **If it is non-`nothing`, it was validated or built here.** Built: from the
   plan's own spaces at the plan's own radius. Adopted: only after
   `validate_dependencies` agreed on both space stamps in their own roles, the
   radius (`plan radius <= graph radius`), the narrow-phase tag exactly, and
   that it is a whole-destination-space relation rather than a row view.
5. **It never changed the plan's other guarantees.** Construction still reads no
   source data, builds no weights, is deterministic in thread count (the
   `_chunkgraph` body is untouched: rows are written by index), and makes no
   network metadata request.
6. **`nothing` is honest, not an error.** A default plan owns no relation and
   says so, rather than lazily materializing one on first ask — which would make
   the accessor building after all.

## 3. The Phase 2 gate, enforced and tested

> One logical plan exposes exactly one validated relation, and a narrow phase
> cannot be supplied after the plan exists.

Both halves are asserted, in `lib/GlobalRegridding/test/test_chunkgraph.jl`, in
four new testsets totalling **62 assertions**.

**"Exactly one"** — `a plan owns exactly one relation` (21):

- the default owns none and asks the probe space nothing (`queries == 0`), for
  both the keyword and the positional constructor;
- `dependencies = true` builds one at the plan's *own* radius, and it equals
  `chunk_dependency_graph(dst, src; radius)` pair for pair;
- the accessor is idempotent and adds no query;
- `isempty(methods(chunk_dependency_graph, Tuple{ChunkedPlan}))`, and both
  `chunk_dependency_graph(plan)` and `chunk_dependency_graph(plan; refine)`
  throw `MethodError`;
- `dependencies(::DirectPlan) === nothing`.

**"Not after the fact"** — `a narrow phase cannot be supplied after the plan
exists` (17): the gate is a **method error**, not a runtime throw. For each of
`:refine` and `:narrow`, `Base.kwarg_decl` over every method of
`chunk_dependency_graph`, `restrict`, `regrid` and `regrid!` shows the keyword is
absent, and over `plan_regrid` shows it is present — so `plan_regrid` (and the
`ChunkedPlan` constructor it forwards to) is the only spelling in the package.
`regrid(data; …, refine = f)` throws `MethodError`. The plan's field is on an
immutable struct: `setfield!` throws.

**Validated, not merely stored** — `a plan validates a relation it did not
build` (18): adoption is by reference (`=== g`); a wider-support method, either
space moving, the pair swapped, an unclaimed narrow phase, a claimed one that is
not there, a row view, and a non-graph value are each an `ArgumentError` at
plan construction.

**The API boundary** — `` `refine` reaches the plan through `plan_regrid` only``
(6): a lazy plan carries the tag; an eager plan names each of `dependencies`,
`refine`, `narrow` in the error it throws rather than ignoring them.

## 4. Production: who owns what, and the measurement

### The global plan owns the relation

`scripts/copdem_production.jl:764` — the call site this card exists to remove —
is gone. `dagplan` now constructs the plan and reads the relation off it:

```julia
globalplan = GR.ChunkedPlan(DGG.Conservative(), DGG.Weighted(0.5),
    dstspace, srcspace; budget = config.budget, dependencies = true,
    refine, narrow)
graph  = GR.dependencies(globalplan)
radius = GR.dependency_radius(graph)
```

The radius is now the plan's method's own `support_radius` instead of a
hand-computed `Float64(GR.support_radius(DGG.Conservative(), srcspace))` — the
same number, from the same call, one layer down. `dagplan` returns `globalplan`
beside `graph`, and the driver gained one structural assertion,
`plan.graph === GR.dependencies(plan.globalplan)`, so "one relation, one owner"
is checked at run time and not only asserted here.

Everything else is untouched: the schedule (`affinity_order` over
`sourcesof(graph, d)`), the refcount cache, the prefetcher and the closing
validator all read the same object they read before, by the same names.
`dagplan` takes no `dem` argument, so plan-time still cannot read a source value
or issue a network request even in principle.

[measured] the whole `dagplan` path on the real GLO-90 × IGeo7-L12 pair
(66 175 destination chunks × 26 475 tiles, `-t 8 --gcthreads=4`, Julia 1.12.6,
local tile list, no downloads):

- `graph === GR.dependencies(globalplan)`: **true**
- **326 064 edges**, radius 0, narrow `:none` — and the relation equals
  `chunk_dependency_graph(dstspace, srcspace; radius)` **pair for pair**, with
  **equal `DependencyIdentity` records**
- the graph build is **0.118 s warm** (0.572 s on the first call, compilation
  included) — the same number `dagplan`'s docstring has always claimed and
  within noise of the gates harness's 0.1143 s median for the same pair
- with `refinegraph = true`, 0.159 s
- `refinegraph = true` still narrows (**250 769 edges**), still tags
  `:copdem_tile_lonlat_box`, and is still a subset of the cap relation

### The per-column plans own none — and cannot own a view

The card asked for per-column plans to "share the global graph or validated row
views". **The first half — never rebuild a one-destination graph per column — is
met by construction: `regrid_chunk`'s plans take the default and build nothing.
The second half is not possible at this destination shape, and the identity G3
built is what proves it.**

`regrid_chunk` builds its destination as `DGG.subtree(sys7, a, level)` — a
*rooted one-chunk grid over one ancestor cell's cells*. That is a different
space from `dagplan`'s 66 175-chunk `PartialGrid`: different `ncells`, different
`nchunks`, therefore a different `SpaceStamp`. A row view of the global graph
still stamps the **whole** destination space (by design — a view knows what it
is a view of), so `validate_dependencies` refuses it, correctly. Offering the
per-column plan either the row view or the whole global graph raises
`ArgumentError` naming the mismatch; both refusals are measured below and both
are asserted in the suite on toy spaces of the same shape.

Making it possible would need a genuinely new mechanism — re-stamping a row
view's destination half against a *sub-space* whose chunk caps equal the
selected parent caps. That is sound in principle (G3 §1.3: the destination half
of the relation is a function of the caps alone) but it is not a thing G3 built,
it needs the parent's caps retained to check, and it would buy nothing today
because no per-column consumer reads a relation yet. It is E1's decision, not
G4's, and it is recorded in §7 rather than smuggled in here.

### What a per-column relation would have cost

`benchmark/plan_dependency_ownership.jl`, on the production pair, over 25
columns sampled across the covering, 5 samples each, medians, `-t 8
--gcthreads=4`. The destination in every arm is the space `regrid_chunk` itself
builds: `DGG.subtree(sys7, a, 12)` at `chunklevel = 5`, one chunk. Arm A is what
production does today; B and C are the two things the card told G4 not to do.

| per-column arm | median | bytes | × 66 175 | × 66 175 alloc |
|---|---:|---:|---:|---:|
| **A** `dependencies = nothing` — today | **1.00 µs** | **784 B** | **0.07 s** | **0.05 GB** |
| B `dependencies = true` — a one-destination graph | 512 µs | 541 664 B | 33.8 s | 35.8 GB |
| C `restrict(graph, [d])` — a one-row view | 193 µs | 424 304 B | 12.7 s | 28.1 GB |

**B costs 511× the time and 691× the allocation of A; C costs 193× and 541×.**
Over the run's 66 175 columns that is 33.8 s and 35.8 GB of churn for B, 12.7 s
and 28.1 GB for C, against A's 0.07 s and 0.05 GB. Neither buys anything today,
because nothing in the per-column path reads a relation.

The 784 bytes and 1 µs arm A does pay are the `PerChunk` storage object and the
plan struct — the plan itself, not a relation. `dependencies(plan) === nothing`
is asserted inside the timed loop, so arm A provably builds none.

C's cost is dominated by the `O(nsourcechunks)` source-major transpose G3
measured and warned about: a 26 476-entry offset array to describe a row holding
5 edges. This is that warning, on the exact shape it warned about.

And neither is legal anyway. Offering the per-column plan either relation is
refused at construction, by name:

```
per-column plan offered a one-row view of the global graph:
  refused: ChunkDependencyGraph(26475 source × 1 destination rows, 5 edges,
  radius 0.0, row view of 66175) was not built against this destination space:
  1 chunks against 66175

per-column plan offered the whole global graph:
  refused: ChunkDependencyGraph(26475 source × 66175 destination chunks,
  326064 edges, radius 0.0) was not built against this destination space:
  1 chunks against 66175
```

## 5. `connectedchunkpairs` is gone

Deleted from `lib/GlobalRegridding/src/discovery.jl`, replaced by a comment
saying what took its place. Since PR #69 it ran line for line the loop
`_chunkgraph` runs to fill its destination-major rows — one
`chunkindex(src_space)`, one `candidatechunks!` per `chunkextents(dst_space)`
entry — so it was a second spelling of the same relation with none of the
identity, the CSR, or the reverse direction.

Its only caller was `lib/GlobalRegridding/test/test_lazy.jl:145`, asserting that
batched queries find exactly the union of the per-chunk ones. That assertion now
reads the graph's own rows against the same `connectedchunks` expectation, so
the claim is unchanged and the assertion count is unchanged.

## 6. Gates and suites

### `benchmark/chunk_graph_gates.jl` — the relation is identical

13 cases including the production GLO-90 × IGeo7-L12 pair from a local tile
list, both arms, 26 ndjson rows, `-t 8 --gcthreads=4`, 5 samples.

**Across all 26 rows, every relation field is unchanged**: `edges`,
`demanded_pairs`, `demand_missing`, `oracle_pairs`, `oracle_missing`,
`only_here`, `missing_here`, `destination_chunks`, `source_chunks`, `radius`,
`identity_bytes` and `graph_summarysize_bytes`. That includes the production
pair's 326 064 indexed edges, 326 386 latjoin edges and the 72 / 394 crossing.

Both runs print:

```
cases run: 13.  Oracle-checked: 9.  Oracle skipped: 4.
verdict: PASS on 9 oracle-checked case(s); 4 case(s) unchecked
```

That is a real geometric verdict on nine cases, not a `NOT CHECKED` run.

Timings move within the run-to-run band of a shared box: the production case is
**+0.6 %** (0.1136 → 0.1143 s) on the indexed arm and +1.0 % on latjoin; the
sub-millisecond cases swing ±25 % in both directions, which is the same spread
G3 recorded on unchanged code. Nothing on the build path changed — the same
`_builddependencies` body runs, one call frame deeper.

### Suites

| suite | before (`f9f268a`) | after | delta |
|---|---|---|---|
| `lib/GlobalRegridding/test` | 3 828 pass / 1 broken / 0 fail | **3 891 / 1 / 0** | **+63 pass** |
| `test/systems/crosssystem/regrid.jl` | 202 / 0 fail | **202 / 0** | none |
| `test/scripts/copdem_policy.jl` | 87 / 0 fail | **87 / 0** | none |
| `test/scripts/copdem_source_mode.jl` | 7 / 0 fail | **7 / 0** | none |
| `test/systems/crosssystem/regrid_acceptance.jl` | not baselined | **22 / 0** | after-state |
| `test/systems/CopernicusDEM/runtests.jl` | not baselined | **16 258 / 3 broken / 0 fail** | after-state |

Every delta explained. The +63 is exactly 62 + 1:

- **62** — the four new testsets in §3: 21 + 17 + 18 + 6.
- **1** — one assertion added to the existing `the narrow-phase tag is a claim
  about the relation` testset, checking that a plan-built relation with no
  `refine` is tagged `:none` just as the space-form builder's is.

Nothing else changed count. The existing `refine` testsets kept their assertion
counts exactly; only the spelling moved from `chunk_dependency_graph(dst, src;
refine)` to a plan (`planned_dependencies`, a two-line helper over a
`G4RadiusMethod` that gives the plan the support radius a test wants, since a
plan takes its radius from its method rather than from a keyword). The one
pre-existing broken test is still broken.

`regrid_acceptance.jl` was not in G3's baseline set; it is the other root-suite
file that touches `plan_regrid`/`ChunkedPlan`, so it is reported here as an
after-state. The CopernicusDEM system suite is run because Phase 2's gate table
asks for it at the end of a phase ("system suite and zero-source-read
construction"); it does not touch the plan API at all and its three broken cases
are the pre-existing conformance skips.

Phase 2 gate table, since G4 is the phase's final task:

| gate row | verification | result |
|---|---|---|
| GlobalRegridding | complete package tests | 3 891 / 1 broken / 0 fail |
| DGG integration | cross-system regrid and acceptance tests | 202 / 0 and 22 / 0 |
| CopernicusDEM | system suite; zero-source-read construction | 16 258 / 3 broken / 0 fail; `dagplan` takes no `dem` argument and the measured run opened no raster |
| Graph changes | deterministic CSR, graph bytes/time, actual-cell no-miss proof | relation identical on all 26 harness rows including `graph_summarysize_bytes`; PASS on 9 oracle-checked cases |

## 7. Residual uncertainty

- **The card's "share a row view" is unmet, and §4 says why.** Production's
  per-column destination is a different space, so the identity refuses the view.
  Closing it needs destination-subspace re-stamping, which no card owns and
  which nothing would consume until E1.
- **The source-side fingerprint hole G3 documented is now reachable from a plan.**
  `dependencies = g` runs `validate_dependencies`, so two source spaces with
  equal chunk caps over different hierarchies can be certified for each other's
  relations. G4 did **not** make it materially more dangerous — no code in tree
  adopts a relation across spaces, and the only new adoption path is one a caller
  must ask for explicitly — but it did move the hole from "a function nobody
  calls" to "a public plan keyword". Closing it still needs a `chunkindexstamp`
  hook on the qualified space interface. Recommend it be given a card before a
  caller starts relying on cross-plan reuse.
- **Nothing in tree consumes a plan's relation yet.** The lazy executor still
  runs its own `candidatechunks!` per tile (E1's card). Until then
  `dependencies = true` is only useful to a caller like production's driver that
  reads the relation itself, and the default stays `nothing` for that reason.
- **The production edits are not covered by any test suite.** No suite includes
  `scripts/copdem_production.jl`. `dagplan` was exercised on the real pair by the
  measurement in §4 (both `refinegraph` settings) and the driver's new
  `check(...)` asserts single ownership at run time, but the driver itself was
  not run end to end.
- `benchmark/regridding_plan_baseline.jl` was not re-measured; it times
  `chunk_dependency_graph` on the production pair and will read within the ±1 %
  band §6 records.
- `benchmark/chunk_graph_gates.jl`'s archived `latjoin_graph` still has its own
  `refine` parameter. It is a private reimplementation of the builder #69
  deleted, kept verbatim so the waiver stays auditable, and it is always called
  with `refine = nothing`. It is not a library keyword and was left alone.
- `scripts/copdem_dag_validate.jl` builds three graphs. They are three
  *different* relations over three different destination sets, not a duplicated
  one, so there was nothing for this card to dedupe.
