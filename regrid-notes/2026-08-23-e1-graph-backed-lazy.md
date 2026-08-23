# E1 — the lazy executor is driven by dependency rows

- Date: 2026-08-23
- Card: `regrid-notes/2026-08-21-regridding-simplification-plan.md`,
  "Phase 3 — graph-backed lazy execution", Task E1
- Branch: `claude/e1-graph-backed-lazy`, cut from `claude/g4-plan-owns-graph`
  @ `1228317` (PR #72's head). PR opens into `claude/g4-plan-owns-graph`.
- Commit: `Drive lazy regridding from dependency rows`

Phase 3's gate: **no lazy read performs geometric dependency discovery;
executor and scheduler consume the same object.** Both halves are tests.

## 1. Where a tile's rows come from

`LazyRegridArray` gained one field, `graph`, and lost three: `srcindex`,
`dstcaps`, `srccaps`, plus the cached `radius`. `graph === dependencies(plan)`,
by reference, and `dependencies(A)` is a new accessor that says so.

```julia
function _connectedsource!(out::Vector{Int}, A::LazyRegridArray, t::Int)
    rows = A.tiling.chunksof[t]
    if length(rows) == 1
        _copyrow!(out, sourcesof(A.graph, @inbounds rows[1]))
    else
        _unionrows!(out, A.graph, rows)
    end
    if A.dropempty
        before = length(out)
        filter!(s -> !_allempty(A, s), out)
        A.stats.dropped += before - length(out)
    end
    return out
end
```

**Space-aligned tile.** When the destination's chunks partition the cell axis
into contiguous runs the tiling is `spacetiled`, tile `t` *is* destination chunk
`t`, and `chunksof[t] == [t]`. That tile takes row `t` verbatim, widened from
the CSR's `Int32` to the `Int` chunk numbers the rest of the executor uses.

**Derived tile.** An explicitly chunked or budget-derived tile spans a set of
destination chunks — including the `_overlappingchunks` fallback, where a run
that overlaps no span is given *every* chunk. Its sources are the ascending
union of those rows: concatenate, `sort!`, `unique!`. One pass, one buffer, and
the same answer for any order the rows are visited in, which a test asserts by
feeding the rows reversed.

**Data-dependent filtering stays after adjacency selection.** `knownempty` may
only *drop* a chunk the relation holds; it can never add one. That ordering is
what keeps "what a read may load" equal to "what the relation permits", which is
the property a refcount is derived from.

`DestTiling.capsof` was renamed `chunksof`. It held chunk indices, never caps;
now that the caps live elsewhere the old name was actively misleading.

## 2. What happened to `srcindex` and the per-read queries

Deleted. `lazy.jl` no longer names `chunkindex`, `candidatechunks!` or
`connectedchunks!` anywhere, and `LazyRegridArray` construction no longer calls
`chunkextents` on either space.

The gate is asserted three ways in
`lib/GlobalRegridding/test/test_lazy.jl`, testset `a lazy read performs no
dependency discovery`:

1. **Behaviourally.** The source is a `G4ProbeSpace` (G3/G4's counter, reused
   rather than copied), which counts every `chunktree` call — the single funnel
   both `chunkextents` and the generic `chunkindex` pass through. After plan
   construction the counter is `built > 0`; after building the array, five whole
   reads, three slices and a `collect`, it is still exactly `built`.
2. **Structurally.** `:graph ∈ fieldnames(typeof(A))` with
   `fieldtype(..., :graph) <: ChunkDependencyGraph`; none of `:srcindex`,
   `:dstcaps`, `:srccaps`, `:radius` is a field; and no field of the array has a
   type `<: AbstractVector{<:SphericalCap}`. The array cannot answer a spatial
   query because it holds nothing that could.
3. **By identity.** `dependencies(A) === dependencies(A.plan)`, and two arrays
   over one plan share the identical object — so a scheduler reading the plan is
   reading what the read consults.

`connectedchunks`/`connectedchunks!` still exist in `discovery.jl`. They are
E2's deletion, not E1's; nothing in `src/` calls them any more.

### The default had to move

A lazy read is now a read of the plan's rows, so a plan that owns none cannot
back one. `dependencies` therefore **defaults to building**: `nothing` (the
default) and `true` both build, a `ChunkDependencyGraph` is validated and
adopted, and `false` is the explicit opt-out. Every `ChunkedPlan` constructor
spelling builds one — keyword and both positional forms — leaving the
nine-argument field constructor as the one way to assemble a plan around an
existing relation or none.

`LazyRegridArray` refuses a plan with no relation at construction, naming the
keyword that caused it, rather than letting a `nothing` reach a tile query.

## 3. Where wave-costing cap metadata lives

On the relation. `ChunkDependencyGraph` gained `extents::Bool`,
`dstcaps::Vector{Cap}` and `srccaps::Vector{Cap}`, read through
`hasextents`, `destinationextents`/`sourceextents` and
`destinationextent`/`sourceextent`. `_wavesize` and `_blockcosts!` take
`(graph, rows)` instead of two cap vectors.

**They cost nothing to keep, because they are the relation's own inputs.**
`_builddependencies` already took `chunkextents` of both sides — the destination
caps to query with, the source caps inside `spacestamp(src)` — and threw the
source vector away. It now takes each once and keeps it. Measured on the
production pair (`benchmark/chunk_graph_gates.jl`), `graph_allocated_bytes` is
**identical on every one of the 13 indexed rows**.

**They are not replicated per worker or per column.**

- On a `DGGSpace`, `chunkextents(space) = space.caps` is a field read and
  `convert(Vector{Cap}, ::Vector{Cap})` aliases, so the graph holds the space's
  own arrays. [measured] on an IGeo7 L4←L3 pair: `sourceextents(g) ===
  chunkextents(src)` and `destinationextents(g) === chunkextents(dst)`, both
  `true`. Production's source and destination are both `DGGSpace`s.
- A `restrict` row view and a `subspace_dependencies` view share the parent's
  two vectors by reference; asserted with `===` in both suites.
- A wave's spawned tasks read the graph's vectors; nothing per task is built.
- The lazy array, which used to hold two cap vectors of its own, holds none.

`destinationextents` is indexed by the *destination space's* chunk number, like
`dstoff`, so a row view indexes it through `globaldestination`;
`destinationextent(g, d)` does that for you.

A graph assembled from the bare-CSR constructor carries no extents (the two cap
keywords default to empty). It answers `hasextents(g) === false`, the extent
accessors throw an `ArgumentError` saying why, and `LazyRegridArray` refuses it
rather than silently degrading its wave costing. `benchmark/chunk_graph_gates.jl`'s
archived `latjoin` arm is the one caller of that constructor and is unchanged.

## 4. The subspace-stamping blocker

G4 left this open: `regrid_chunk`'s destination is a rooted one-chunk subtree
grid, a different space from `dagplan`'s 66 175-chunk `PartialGrid`, so a
`restrict` row view — which still stamps the whole destination — is refused by
`validate_dependencies`, correctly. E1 closed it, and then measured that
production does not want it.

### `subspace_dependencies(g, subspace, destinations)`

Returns `restrict(g, destinations)` re-stamped as the relation of `subspace`.
The soundness argument is one line and it is the reason the check is the check
it is:

> The destination half of the relation is a function of the destination chunk
> caps alone — a row is one `candidatechunks!` query of the source index against
> one cap. So if `subspace`'s chunk `k` has *the same cap* as `g`'s destination
> chunk `destinations[k]`, the row `subspace` would produce for `k` is the row
> `g` already holds.

So cap equality is checked, **exactly**, cap for cap, and a mismatch is an
`ArgumentError` rather than a re-stamp: two nearly-equal caps are two different
relations. The count must match. The source half is untouched and its stamp is
carried over verbatim, because it is literally the parent's columns. The
destination stamp is *computed from the sub-space*, not asserted.

`validate_dependencies` needed one change, and it is a simplification: the
`destinations === nothing` branch used to reject any restricted graph and then
check `g.ndst == g.id.dst.nchunks`. The second check subsumes the first — a row
view has fewer rows than the space it stamps, a sub-space view has exactly as
many as the space *it* stamps — so the `isrestricted` rejection is gone and the
count check now carries the whole meaning, with an error message that names
`subspace_dependencies` as the way out.

**What it does not establish, stated plainly.** Equal caps do not imply equal
cell sets. A caller that hands over an unrelated space with coincidentally equal
chunk caps gets a certified relation for a destination it did not mean. That is
the same hole `spacestamp` documents, in the place it matters most, and it is
why the docstring says passing a genuine sub-space is the caller's obligation.

Proved on toy spaces (`test_chunkgraph.jl`, `a row view re-stamped onto a
sub-space is that space's relation`, 25 assertions) and on the real DGG shape
(`test/systems/crosssystem/regrid.jl`, `a column adopts the whole covering's
relation and reads the same`, 34 assertions over three columns of an IGeo7
L4/level-2-chunk covering):

- row for row it equals the parent's rows for those chunks;
- `graph_pairs` and the whole `DependencyIdentity` equal what the sub-space
  builds for itself with `chunk_dependency_graph`;
- a plan over the sub-space adopts it **by reference**, with no `destinations`
  argument, while a plain `restrict` view and the whole covering's relation are
  both still refused;
- **the values do not move**: a lazy regrid through the adopted view is `==` to
  one through the column's own relation, and `≈` to the eager whole-domain
  answer for that column at `rtol = 1e-12`;
- refusals: wrong chunk count, a destination outside the parent's rows,
  descending destinations, a space whose caps are not the parent's caps for
  those chunks, and a graph carrying no extents.

### Production does not adopt it, and here is the number

`benchmark/plan_dependency_ownership.jl`, real GLO-90 × IGeo7-L12 pair from a
local tile list, 66 175 destination chunks × 26 475 tiles, 25 columns sampled
across the covering, 5 samples each, medians, `-t 8 --gcthreads=4`:

| per-column arm | median | bytes | × 66 175 | × 66 175 alloc |
|---|---:|---:|---:|---:|
| A `dependencies = false` — the floor | 1.14 µs | 784 B | 0.1 s | 0.1 GB |
| **B default: builds one — today** | **500 µs** | **541 568 B** | **33.1 s** | **35.8 GB** |
| C `restrict(graph, [d])` alone | 190 µs | 424 368 B | 12.6 s | 28.1 GB |
| D subspace view, adopted by a plan | 408 µs | 425 520 B | 27.0 s | 28.2 GB |

**D/B = 0.82× time, 0.79× bytes.** Re-stamping saves 92 µs and 116 KB per
column — 6.1 s and 7.7 GB over the whole run, against a recorded run wall time
of 8.81 h. Both arms pay the same `O(nsourcechunks)` transpose over 26 475
tiles, and adoption still pays `spacestamp(srcspace)` to certify the source
half; what a rebuild adds on top is one `candidatechunks!` query.

So `regrid_chunk` keeps the default and builds its own one-row relation.
Threading the global relation and each column's row number through the worker
loop would buy 0.02 % of the run. The benchmark also asserts, inside the timed
loop, that B's and D's relations hold the same sources for every sampled column,
and prints the three adoption outcomes: the plain row view **refused**, the
whole global graph **refused**, the re-stamped view **ADOPTED**.

The cost E1 does impose on production is arm B itself: **33.1 s and 35.8 GB of
churn over 66 175 columns**, 0.1 % of the run's wall time and 1.1 MB/s of
allocation. It is not free and it is not hidden.

## 5. The proof obligation

`nonzero eager edges ⊆ graph rows ⊆ broad cap relation`, as a test, on small
spaces: `test_chunkgraph.jl`, `eager weights ⊆ graph rows ⊆ the broad cap
relation`.

The left-hand term is new and is in the shared oracle file
(`lib/GlobalRegridding/test/graphoracles.jl`, reused not copied):

```julia
eager_pairs(dst, src; method = GR.Conservative())
```

It builds the **eager** whole-domain `WeightBlock` and reads its nonzero entries
out as `(dstchunk, srcchunk)`. It consults no cap, no chunk extent and no chunk
index, so it is an oracle for the relation rather than a second opinion, and it
differs from the existing `contributing_pairs` in what it is evidence about:
that one is independent geometric truth, this one is the weight builder's own
answer. A relation missing a pair from it is a relation under which a lazy read
never loads a source the eager path weighted — the precise way the two could
return different numbers.

The middle term is read the way the executor reads it: `graph_pairs` of
`dependencies(ChunkedPlan(Conservative(), …, dst, src))`.

Cases: four generic-index pairs (a plain toy pair, a polar source band, an
antimeridian regional source, and a coarser destination than source) assert both
containments; a misaligned-raster source asserts the left-hand one and
`eager ⊆ capjoin` instead of `rows ⊆ capjoin`, **because the right-hand
containment is false on a native hierarchy** and the suite already pins that it
is (`the cap join is an identity only on the generic index`: the two relations
cross in both directions, which is why #69 exists). Stating the obligation
without that carve-out would have been stating something untrue.

## 6. Gates and suites

### `benchmark/chunk_graph_gates.jl` — the relation is identical

Run on the branch point (`1228317`, in a detached worktree) and on this branch:
13 cases including the production pair from a local tile list, both arms,
26 ndjson rows each, `-t 8 --gcthreads=4`, 5 samples.

**Every relation field is identical on all 26 rows**: `edges`,
`demanded_pairs`, `demand_missing`, `oracle_pairs`, `oracle_missing`,
`only_here`, `missing_here`, `destination_chunks`, `source_chunks`, `radius`,
and `identity_bytes` (still 1040 B for the record plus its stamps). That
includes the production pair's 326 064 indexed edges, 326 386 latjoin edges and
the 72 / 394 crossing.

Both runs print:

```
cases run: 13.  Oracle-checked: 9.  Oracle skipped: 4.
verdict: PASS on 9 oracle-checked case(s); 4 case(s) unchecked
```

A real geometric verdict on nine cases, not a `NOT CHECKED` run.

What did move, and why:

- **`graph_summarysize_bytes`, indexed arm: 1.8×–3.8×** (production 3 350 064 →
  7 056 168). `Base.summarysize` counts everything *reachable*, and the graph now
  reaches the two cap vectors it always derived itself from. On `DGGSpace`s
  those are the spaces' own arrays — [measured] `===` — so nothing new is
  allocated and nothing extra is held: `graph_allocated_bytes` is unchanged on
  every indexed row. On a `RasterGrid` the vectors are genuinely built, and the
  graph now retains what the lazy array used to build and retain twice over.
- **`graph_summarysize_bytes`, latjoin arm: +88 B on every row.** The two empty
  `Cap[]` headers and the `Bool`. That arm assembles from bare CSR and carries no
  extents, which is the point of the flag.
- **Timings.** Production `:indexed` **+0.7 %** (0.11643 → 0.11722 s), within
  the ±1 % band G3 and G4 both recorded on unchanged code. Sub-millisecond cases
  swing ±25 % in both directions, the same spread as before. Nothing on the
  build path changed except keeping two references.
- `restrict_*`/`rebuild_*` moved because `rebuild_graph` in the harness was
  updated to mirror `_builddependencies`: it takes each side's caps **once** and
  uses them for both the stamp and the relation, instead of paying
  `dependency_identity` and `chunkextents(dst)` separately. That is what a
  rebuild costs now; the previous spelling would have measured a rebuild nobody
  performs.

### Suites

Everything below at `-t 8 --gcthreads=4`, Julia 1.12.6, except the
GlobalRegridding suite at `-t 4`.

| suite | before (`1228317`) | after | delta |
|---|---|---|---|
| `lib/GlobalRegridding/test` | 3 891 / 1 broken / 0 fail | **3 999 / 1 / 0** | **+108 pass** |
| `test/systems/crosssystem/regrid.jl` | 202 / 0 | **236 / 0** | **+34** |
| `test/systems/crosssystem/regrid_acceptance.jl` | 22 / 0 | **22 / 0** | none |
| `test/scripts/copdem_policy.jl` | 87 / 0 | **89 / 0** | **+2** |
| `test/scripts/copdem_source_mode.jl` | 7 / 0 | **7 / 0** | none |
| `test/systems/CopernicusDEM/runtests.jl` | 16 258 / 3 broken / 0 fail | **16 258 / 3 / 0** | none |

Every delta accounted for, assertion by assertion.

**+108 in GlobalRegridding**, exactly:

| where | net | what |
|---|---:|---|
| `a plan owns exactly one relation` | +1 | −5 asserts about the old `nothing` default, +6 about the new build-by-default in all three constructor spellings and the `false` opt-out |
| `` `refine` reaches the plan through `plan_regrid` only`` | +1 | the default-plan assertion became `narrowphase == :none`, and one was added for `dependencies = false` |
| `eager weights ⊆ graph rows ⊆ the broad cap relation` (new) | +15 | 4 generic cases × 3, plus the raster case × 3 |
| `a row view re-stamped onto a sub-space…` (new) | +25 | §4 |
| `a lazy read takes its sources from the plan's relation` (new) | +6 | identity, `spacetiled`, and both destination chunks' rows |
| `a derived tile takes the sorted union of its rows` (new) | +35 | 2 setup, then 3 per tile over 11 derived tiles |
| `a lazy read performs no dependency discovery` (new) | +12 | §2 |
| `wave costing reads the relation's extents, and copies none` (new) | +10 | §3 |
| `a plan with no relation cannot back a lazy read` (new) | +3 | the refusal and its message |
| the wave testsets | 0 | same assertions, `(graph, rows)` instead of two cap vectors |

**+34 in `regrid.jl`**: one `hasextents`, then 11 assertions on each of three
columns (§4).

**+2 in `copdem_policy.jl`**: the hand-built relation in `the seam holds a real
chunk_dependency_graph` now also asserts `hasextents` and that its source
extents are the very vector it was built from. Its `GR._chunkgraph` call was
updated for the private builder's new argument (source caps rather than a source
chunk count) — the one place outside `chunkgraph.jl` that calls it.

The one pre-existing broken GlobalRegridding test is still broken; the three
pre-existing broken CopernicusDEM conformance skips are unchanged.

## 7. Values, repeated reads, and residency

`scratchpad/residency.jl` (not committed — it is a measurement, not a suite):
one 720×360 tiled raster onto IGeo7 L4 with level-2 chunks, 24 012 destination
cells, 492 chunks ← 72 chunks, a 4 MiB budget, `-t 4 --gcthreads=4`. Run at the
branch point and on this branch.

| | before | after |
|---|---|---|
| `LazyStats` | loads 328, hits 3654, skipped 0, dropped 0, **peak 2 073 600 B** | **identical** |
| source `readblock!` calls | 400 | **400** |
| repeated whole read identical | true | **true** |
| all seven slices identical | true | **true** |
| differently tiled plan identical | true | **true** |
| max abs(eager − lazy) | **2.44249e-15** | **2.44249e-15** |
| live at construction: plan | 63 688 B | 82 072 B |
| live at construction: array (minus source and plan) | 46 424 B | 43 376 B |

The executor's own residency counter is **byte-identical**, which is the number
the budget bounds. Whole-process `Sys.maxrss` deltas were 0.0–61.6 MiB in both
directions across runs at this size — compilation noise, not a signal, and not
reported as one.

What changed about what is held live: the plan now holds a relation
(+18 384 B here) and the array no longer holds a source index or two cap vectors
(−3 048 B). Net **+15 336 B** on a 492 × 72-chunk pair.

At production scale the graph's genuinely-new bytes are its four CSR arrays —
`dstoff` 529 KB, `srcof` 1.30 MB, `srcoff` 212 KB, `dstof` 1.30 MB ≈ 3.35 MB,
unchanged from before, because both cap vectors are aliases to the two
`DGGSpace`s' own `caps` fields. For a per-column plan the retained cost is
dominated by the `O(nsourcechunks)` transpose: **≈ 212 KB live per column while
that column's plan is alive**, so ≈ 8.5 MB across 40 concurrent workers. On the
64 GB target that is not a constraint; it is the one new live cost E1 adds and it
is stated rather than assumed.

The 2.44e-15 eager/lazy difference is pre-existing accumulation order — the lazy
path sums chunk pair by chunk pair — and it is **bit-identical before and
after**, which is the evidence that E1 moved no value.

## 8. Production

`scripts/copdem_production.jl`:

- `regrid_chunk`'s docstring now records that its plan owns a one-row relation,
  that this is forced rather than chosen, and both routes to adopting the global
  one with the numbers from §4.
- `dagplan`'s docstring records the same, and that its graph is still the only
  relation over the whole covering.
- **New: `graphmisscheck`**, called from `run` beside the existing graph
  `check`s. It samples `GRAPHMISS_SAMPLE = 16` columns, builds each column's own
  relation through the same `DGG.columncell` mapping `regrid_chunk` uses, and
  checks that no column demands a tile the global relation does not hold. This
  is the card's "production graph-miss validation": the refcount cache retires
  tiles by the global relation's consumer counts while a column now reads its
  own relation's row, and a miss would be a retired tile the column reloads —
  or, with no permits left, waits forever for.

**What was exercised, and how.** No suite includes this script, so the
measurements are the exercise:

- `benchmark/plan_dependency_ownership.jl` ran the whole `dagplan` path twice
  on the real pair (both `refinegraph` settings) with a local tile list and no
  downloads: `graph === dependencies(globalplan)` true, 326 064 edges, relation
  equal to `chunk_dependency_graph(dst, src; radius)` pair for pair with equal
  `DependencyIdentity`, 0.130 s warm, and `refinegraph = true` still narrowing
  to 250 769 edges tagged `:copdem_tile_lonlat_box` and still a subset.
- `graphmisscheck` was run on the real pair from a scratch harness that includes
  the production script and builds the layout with `DGG.SubzoneLayout(sys7, 12,
  5)` — the same object the store holds — after asserting `columncell(layout, c)
  == cellindex(levelgrid(sys7, 5), c)` on a sample, so the mapping under test is
  the mapping the run uses. Result: **PASS, 16 columns sampled, 0 misses, in
  0.430 s**, with the `FAILURES` counter left at 0. A **negative control** — the
  same call against a deliberately edgeless relation — correctly **FAILED** with
  60 misses and incremented `FAILURES`, so the check can fail.
- The driver itself was still not run end to end, and no store was opened.

## 9. Residual uncertainty

- **The per-column transpose is the cost, and nothing here removes it.** Arm B's
  500 µs is dominated by `spacestamp(srcspace)` over 26 475 caps and the
  `O(nsourcechunks)` source-major transpose, for a relation with one row and a
  handful of edges. The lazy executor never calls `consumersof`. A relation that
  built its source-major direction on demand would cut both B and D, and would
  make a per-column plan nearly free — but it means memoized state on an object
  documented as immutable, and it is not E1's card.
- **The source-side fingerprint hole is unchanged and is now on one more path.**
  `subspace_dependencies` carries the parent's source stamp over verbatim, which
  is right — the columns are literally the parent's — but the plan that adopts
  the result still certifies the source half by `spacestamp`, so two source
  spaces with equal chunk caps over different hierarchies remain
  interchangeable. E1 did not widen it: the destination half, which is what
  re-stamping touches, is the half `spacestamp` documents as *sound*. The fix is
  still a `chunkindexstamp` hook on the qualified space interface, and it still
  has no card.
- **Equal caps are not equal cells.** `subspace_dependencies` cannot check that
  the sub-space's chunk `k` holds the parent's chunk `destinations[k]`'s cells,
  only that it reproduces its cap. Documented on the function; the obligation is
  the caller's.
- **`hasextents` is a flag, not a derivation.** A caller constructing a graph
  through the bare-CSR constructor with one side's caps and not the other gets an
  `ArgumentError`; with neither, a graph no lazy read can use. That is deliberate
  — degrading wave costing silently would be worse — but it is a new way for a
  hand-built graph to be rejected late.
- **`GRAPHMISS_SAMPLE = 16` is a sample.** It is not a proof over 66 175
  columns; the proof is the test in `regrid.jl` on the same shape plus the
  construction argument in `graphmisscheck`'s docstring. Raising it costs
  ~0.5 ms per column.
- **`benchmark/regridding_plan_baseline.jl` was not re-measured.** It times
  `chunk_dependency_graph` on the production pair and will read within the
  +0.7 % the gates harness recorded.
- The `latjoin` arm's `graph_summarysize_bytes` is now 88 B larger on every row.
  It is an archived builder kept for the G2 waiver; the number is noise against
  its purpose, but it is not literally the number G4 recorded.
