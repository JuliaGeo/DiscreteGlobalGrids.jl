# G1 — dependency graph correctness and performance gates

Task G1 of `2026-08-21-regridding-simplification-plan.md` (Phase 2). Production
graph construction is **unchanged**: this task adds gates and instrumentation
only, so that a later task can swap builders and get an objective verdict without
writing new instrumentation first.

It also settles the question **G2 skipped**. PR #69 (`da9e737`) cut
`chunk_dependency_graph` over to indexed construction as a correctness fix,
before G1 existed and therefore without G2's performance gate — "recover the
recorded latitude-join performance before removing it". The old builder is gone,
so the gate cannot be run as written. This note re-creates the deleted builder
inside the harness and answers the question retroactively. **The gate fails on
its literal wording** — see §5 — and the cutover is still right.

Landed as `48b9815`, with the review findings applied in `07d24ea`, on top of
`9adef54`.

## 1. The oracles

Three relations over the same `(dstchunk, srcchunk)` pair space:

| relation | how it is built | what it is for |
|---|---|---|
| **truth** | real cell rings, clipped through the same spherical kernel the conservative weight builder uses; positive area, or (at nonzero support) cells within `radius` of each other, measured **cell to cell** — every vertex of one against every edge of the other, in both directions | the geometric fact. No cap, no chunk extent, no chunk index takes part |
| **demand** | one `candidatechunks!` per destination cap on `chunkindex(src)` | what a lazy read can ask for. A relation missing any of this cannot be a refcount |
| **cap join** | brute-force `DilatedIntersects` over `chunkextents(dst) × chunkextents(src)` | the definition of the deleted builder's relation |

The invariant every builder must keep is **truth ⊆ graph**. The invariant that
makes the graph a *refcount* is **demand ⊆ graph** — and post-#69, with no
`refine`, **demand = graph** exactly, because the builder issues precisely the
queries `demanded_pairs` replays. The tests assert the equality, not the
containment: `truth ⊆ graph` on its own would pass on a builder that returned
the complete bipartite relation, and G3/G4 need a gate that catches a *changed*
relation, not only a lossy one.

### Where they live

All four relations are defined **once**, in
`lib/GlobalRegridding/test/graphoracles.jl` (module `ChunkGraphOracles`), and
`include`d by path from the three sites that check them. That file is the
innermost of the three — the root suite and `benchmark/` both already depend on
GlobalRegridding — so it adds no dependency edge in the wrong direction, and it
follows the `test/helpers.jl` convention: a module, included by path, defined
inside the including module. It exports `contributing_pairs` (the truth oracle,
`O(ncells²)`), `demanded_pairs`, `capjoin_pairs` and `graph_pairs`.

- `lib/GlobalRegridding/test/test_chunkgraph.jl`
  - testset **"no geometrically contributing pair is dropped"** — seven cases:
    generic (packed R-tree) source at radius 0 and 0.2; the shipped `RasterGrid`
    quadtree at radius 0 and 0.1; nonuniform raster chunks; a polar source band;
    an antimeridian regional source. Each asserts `truth ⊆ graph`,
    `graph == demand`, **and** `truth ⊆ capjoin` — the last catches a chunk cap
    that is too tight, which is a `chunkextents` bug the graph would inherit.
  - testset **"the cap join is an identity only on the generic index"** — the
    brute-force cap identity gate, plus the measured refutation of the "upper
    bound" reading (below).
- `test/systems/crosssystem/regrid.jl`
  - the same truth oracle over real DGG geometry, testset **"no geometrically
    contributing pair is dropped"**: complete IGeo7 L2/chunk-1 from IGeo7
    L3/chunk-2 (all twelve pentagons, both poles, aperture-7 children that reach
    outside their parent's boundary), IGeo7 from S2 (two systems sharing no cell
    edges), a rooted subtree destination, a scattered non-rooted subset, and a
    chunked raster source into a DGG destination at radius 0 and 0.05. Each
    asserts `truth ⊆ graph`, `graph == demand`, `truth ⊆ capjoin`, and — on the
    four DGG-source cases — `graph ⊆ capjoin`.
  - the same file's existing **"the dependency graph holds every pair the chunk
    index answers"** covers demand-domination on the CopDEM level-0 frontier.
    Post-#69 that sweep compares `candidatechunks!` against a graph *built from*
    `candidatechunks!`, so it still catches a mis-assembled CSR but can no
    longer catch a wrong *choice* of index; its own comment now says so. The
    same testset gained the cap-crossing assertion described below.
- `benchmark/chunk_graph_gates.jl` — all four relations again, over the cases
  too large to assert on every CI run, including the production pair.

Suites after the additions and the review fixes: GlobalRegridding **[measured]
3737 pass / 1 broken / 0 fail** (baseline on this tip was 3702/1/0; the first
G1 landing took it to 3730, and the `graph == demand` equality added the other
seven. The broken one is pre-existing). `test/systems/crosssystem/regrid.jl` runs
clean in **[measured] 47 s**, the new testset contributing **[measured] 30**
assertions in **[measured] 4.4 s** and the cap-crossing window taking "the
dependency graph holds every pair the chunk index answers" from 7 to
**[measured] 9**; `crosssystem/runtests.jl` **[measured] 5222 pass / 0 fail**
(untouched — `regrid.jl` is included by the root `runtests.jl`, not by it).

### The brute-force cap identity gate, and its limits

The card asks for "brute-force cap identity checks for the generic R-tree path"
and to "treat the flat cap relation as an upper bound for raster/DGG native
hierarchies". The first holds exactly. **The second is false as stated, and the
tests now pin that.**

| source index | relation vs the cap join |
|---|---|
| generic `_packedchunkindex` (any space with no `chunkindex` of its own) | **equal**, at every radius tried — it is a packed R-tree over the very caps `chunkextents` reports |
| `DGGSpace` via `_dggcandidatechunks!` | **subset**. The descent can prune a chunk whose own cap intersects, never add one whose cap does not |
| `RasterGrid` quadtree | **crosses in both directions**. A straddling leaf pushes every chunk in its rectangle with no per-chunk cap test, so the index holds pairs the cap join rejects; and the cap join holds pairs the descent never reaches |
| `CopernicusDEM` level-0 frontier | **crosses in both directions**. Pinned in CI on the whole level-0 GLO-90 frontier into IGeo7 L3/chunk-2 — no tile list, no download, **[measured] 0.36 s** for the brute-force join over 492 × 64 800 chunk pairs: **[measured] 6** pairs the index holds and the cap join rejects, **[measured] 10** the other way. At production scale (GLO-90 × IGeo7-L12) the first number is **[measured] 72**; *that* figure is harness-only, reproducible only with a local tile list, and is not asserted anywhere |

So the cap join bounds the graph in *neither* direction on the two shipped
native hierarchies. What both relations do satisfy is `truth ⊆ ·`, and that is
the invariant the tests assert everywhere.

## 2. The harness

```
julia -t 8 --project=benchmark benchmark/chunk_graph_gates.jl
```

Environment:

| variable | meaning |
|---|---|
| `DGG_GRAPH_GATE_SAMPLES` | timed samples per arm (default 5) |
| `DGG_GRAPH_GATE_CASES` | comma-separated case-name filter |
| `DGG_GRAPH_GATE_NDJSON` | append one JSON line per `(case, arm)` |
| `DGG_GRAPH_GATE_RAW` | set to `1` to add the unprefiltered `:latjoin_raw` arm (default off) |
| `DGG_COPDEM_TILELIST` | a **local** Copernicus tile list. The production case is skipped without it; the harness never downloads |

Two arms by default:

- `:indexed` — production `chunk_dependency_graph`.
- `:latjoin` — the deleted latitude-sorted cap join, reimplemented in the
  harness from `ba2bbfa^:lib/GlobalRegridding/src/chunkgraph.jl`, with its
  private helpers copied in so it cannot drift as `chunkgraph.jl` changes.
- `:latjoin_raw` — the same with the latitude prefilter off, which is what
  isolates the prefilter's contribution. Opt-in behind `DGG_GRAPH_GATE_RAW=1`:
  it is an attribution aid, not a candidate builder, and it is an order of
  magnitude the slowest arm in the matrix.

`:latjoin` has a **sunset condition**, stated in the file header: it goes when
G2's waiver is retired, or when the production builder stops being comparable to
a flat cap join, whichever comes first. It reconstructs a builder that no longer
exists and is not maintenance the repo owes indefinitely. The oracles have no
sunset — G3 and G4 use them to prove they did not change the relation.

### What the gate line means

The summary distinguishes **oracle-checked** cases from **skipped** ones. It has
to: `oracle_missing` is `-1` when the `O(ncells²)` sweep was skipped as too
large, and `demand_missing` is 0 for `:indexed` *by construction* — post-#69 the
builder issues exactly the queries `demanded_pairs` replays, so that column can
only read 0 on this arm and is not evidence. A run over nothing but skipped
cases (`DGG_GRAPH_GATE_CASES=copdem90-igeo7-l12`, for instance) now prints
`verdict: NOT CHECKED`, never `PASS`.

Per row: destination/source chunk and cell counts, radius, edges, space-build
seconds, `ChunkedPlan` seconds, graph seconds (min and median), complete plan
seconds, allocated bytes, GC seconds, `Base.summarysize` of the graph, peak RSS
growth, microseconds per destination, and the four correctness columns
(`demanded_pairs`, `demand_missing`, `oracle_pairs`, `oracle_missing`) plus the
two relation-difference columns (`only_here`, `missing_here`, both against
`:indexed`).

### Provenance stamping

Every ndjson row carries the Julia version, the thread count, the GC thread
count, the repo HEAD, and the **git revision and tree hash** of GeometryOps,
GeometryOpsCore and ConservativeRegridding, read from the active manifest at
runtime rather than hardcoded. Those three are pinned to branches rather than
releases, and spherical-predicate/Foster-Hormann work in flight on the
GeometryOps side will move clipping cost materially. Every number below is a
**pre-change baseline**; a re-run after a pin bump relabels itself, and rows with
different stamps are not comparable.

This run:

```
julia 1.12.6   threads 8   gcthreads 8
GeometryOps            697c7cc81b3d61666cf42bf697bf186cd7f6b2e2  tree 9944edf3893eb267769150fb563b2ccf01df7cf1
GeometryOpsCore        main                                      tree 0896ad28d34338853bd9d335fee15053d748dc58
ConservativeRegridding e8fc67f420da0362d59b354e2db6fcde8107fec0  tree c571dff5cb9bc1fb2a8c3112c6870baf13b91946
repo 9adef543b7283383726378a7fbb4874697dbcbc0   samples 5
```

## 3. Recorded results

Julia 1.12.6, 8 threads, 5 samples, medians, `DGG_GRAPH_GATE_RAW=1` (the `raw s`
column is opt-in). `idx-only`/`lat-only` are the relation difference in each
direction. `demand-miss` is the number of pairs a lazy read would ask for that
the arm does not hold — for `:indexed` it is 0 in every case by construction, so
only the `:latjoin` column is shown.

| case | dst chunks | src chunks | indexed s | latjoin s | raw s | idx/lat | idx alloc B | lat alloc B | idx graph B | lat graph B | idx edges | lat edges | idx-only | lat-only | truth pairs | idx miss | lat demand-miss |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| raster-small | 72 | 24 | 0.0003 | 0.0001 | 0.0001 | 2.3× | 79 112 | 52 216 | 4 240 | 3 560 | 405 | 320 | 92 | 7 | 193 | 0 | **92** |
| raster-support (r=0.05) | 72 | 24 | 0.0003 | 0.0001 | 0.0001 | 2.1× | 83 656 | 53 208 | 4 608 | 3 760 | 451 | 345 | 114 | 8 | 222 | 0 | **114** |
| raster-nonuniform | 72 | 49 | 0.0005 | 0.0001 | 0.0001 | 3.1× | 90 040 | 58 024 | 4 680 | 4 664 | 435 | 433 | 33 | 31 | 264 | 0 | **33** |
| dgg-complete | 72 | 492 | 0.0012 | 0.0002 | 0.0002 | 7.0× | 96 584 | 104 856 | 18 280 | 18 280 | 1 692 | 1 692 | 0 | 0 | 1 332 | 0 | 0 |
| dgg-crosssystem | 72 | 24 | 0.0002 | 0.0001 | 0.0001 | 1.6× | 92 520 | 51 496 | 3 256 | 3 368 | 282 | 296 | 0 | 14 | 194 | 0 | 0 |
| dgg-rooted | 7 | 72 | 0.0002 | 0.0001 | 0.0001 | 2.3× | 6 832 | 9 104 | 1 048 | 1 048 | 23 | 23 | 0 | 0 | 19 | 0 | 0 |
| dgg-sparse | 175 | 24 | 0.0002 | 0.0001 | 0.0001 | 2.2× | 143 880 | 58 168 | 5 168 | 5 336 | 418 | 439 | 0 | 21 | 267 | 0 | 0 |
| polar-source | 72 | 15 | 0.0001 | 0.0001 | 0.0002 | 1.4× | 53 832 | 41 704 | 1 704 | 1 536 | 97 | 76 | 25 | 4 | 41 | 0 | **25** |
| antimeridian-source | 72 | 12 | 0.0001 | 0.0001 | 0.0001 | 1.5× | 52 360 | 40 984 | 1 528 | 1 360 | 78 | 57 | 22 | 1 | 34 | 0 | **22** |
| raster-4320-162chunks | 3 432 | 162 | 0.0303 | 0.0005 | 0.0008 | **60.3×** | 871 232 | 437 848 | 82 104 | 106 424 | 6 640 | 9 680 | 8 | 3 048 | — | — | **8** |
| raster-4320-1800chunks | 3 432 | 1 800 | 0.0626 | 0.0017 | 0.0060 | **35.9×** | 4 150 744 | 700 688 | 170 808 | 198 424 | 16 090 | 19 542 | 120 | 3 572 | — | — | **120** |
| dgg-large | 3 432 | 492 | 0.0226 | 0.0006 | 0.0016 | **36.0×** | 485 560 | 503 048 | 135 400 | 135 400 | 12 972 | 12 972 | 0 | 0 | — | — | 0 |
| **copdem90-igeo7-l12** | 66 175 | 26 475 | **0.1219** | **0.0213** | 1.0011 | **5.7×** | 17 123 768 | 10 614 224 | 3 349 944 | 3 352 520 | 326 064 | 326 386 | 72 | 394 | — | — | **72** |

Cases with `—` in the truth columns are too large for the `O(ncells²)` sweep;
their correctness is covered by the demand column and by the small cases of the
same shape.

The production pair at 4 threads, for comparison against the archived figure:

| arm | seconds (median) | allocated | graph bytes | edges |
|---|---:|---:|---:|---:|
| indexed | 0.2289 | 13 410 144 | 3 349 944 | 326 064 |
| latjoin | 0.0376 | 10 591 960 | 3 352 520 | 326 386 |
| latjoin_raw | 1.9822 | 10 591 960 | 3 352 520 | 326 386 |

The `:latjoin` arm reproduces the archived baseline's relation exactly —
**[measured] 66 175** destination chunks, **[measured] 26 475** source chunks,
**[measured] 326 386** edges, **[measured] 3 352 520**-byte graph — the numbers
`benchmark/regridding_plan_baseline.jl`'s header used to carry as its own. That
header now records them as superseded and states what the script measures today
(326 064 edges, 3 349 944 bytes, 0.229 s at t4), since #69 changed the builder
underneath it; its behaviour is unchanged. The archived 0.0594 s median was taken
on a differently loaded machine; **[measured] 0.0376 s** here at the same thread
count.

Space construction dwarfs both arms on the production pair: **[measured] 11.1 s**
to build the two spaces, against 0.12 s for the graph. `ChunkedPlan`
construction is **[measured] ~2 µs** — lazy plans hold no weights — so "complete
plan time" is graph time plus noise.

## 4. What the oracles caught

Nothing wrong with the current builder: `oracle_missing` and `demand_missing` are
**[measured] 0** for `:indexed` in every case, at every radius, on every space
shape. That is the result G1 wanted and it is not a null result — it is the first
time the relation has been checked against real cell geometry rather than against
a second cap construction.

What the oracles caught in *themselves* is worth recording, because it is the
kind of error an oracle can make silently. The nonzero-support branch originally
tested vertex-to-vertex distance, which under-approximates: two cells can come
within `radius` edge-to-edge with no vertex pair that close. Replacing it with a
cell-to-cell distance (every vertex of one against every edge of the other, both
directions) grew `raster-support`'s truth set from 209 to **[measured] 222**
pairs — 13 genuinely contributing pairs the `r > 0` oracle had not been demanding.
All three arms hold all 222, so nothing was wrong with any builder; the *gate*
was 6% weaker than it read.

What the oracles *do* catch is the deleted builder, on cases nobody had run it
against:

- **The cap join misses demanded pairs on the shipped raster path at every
  size**: 92 of 405 demanded pairs on a 24-chunk raster, 8 of 6 640 on the
  162-chunk 4320×2160 raster, 120 of 16 090 on the 1800-chunk one, and
  **[measured] 72** of 326 064 on the production pair — the
  exact defect `2026-08-23-chunk-dag-coverage.md` diagnosed, now reproduced by a
  committed harness rather than a one-off script.
- **It misses them on polar and antimeridian regional sources too** (25 and 22),
  which no previous measurement covered.
- **The latitude prefilter was worth 47× on the production pair** (0.0213 s vs
  1.0011 s at t8) and only ~1.5–3.5× at small sizes. Read that as a diagnosis of
  the deleted builder, not as a requirement on a future one: it is the deleted
  builder measured against itself, in a configuration — an unprefiltered flat
  cap join — that ships nowhere. What the 47× says is that a flat
  `O(ndst × nsrc)` scan over every source cap *is* the defect, and a latitude
  band is a bandaid over it that recovers most but not all of the loss. The
  shipped builder needs no band, because the pruning happens inside the source
  space's own index instead of being bolted onto a scan: the indexed arm **beats
  the unprefiltered cap join by 8.2×** (0.1219 s vs 1.0011 s) while losing to the
  prefiltered one by 5.7×. A future builder should be judged on whether it prunes
  structurally, not on whether it kept the band.

## 5. Retroactive verdict on G2's performance gate

G2's action list says: *"Recover the recorded latitude-join performance before
removing it. The archive records the 0.0594 s median... Remove the
latitude-sorted join only after both gates pass."*

**Correctness gate: PASS.** The indexed relation holds every geometrically
contributing pair on every case in the matrix, and every demanded pair on the
production pair, where the cap join misses 72.

**Performance gate: FAIL, on the literal wording.** The indexed builder is
**[measured] 5.7× slower** than the latitude join on the production pair at 8
threads (0.1219 s vs 0.0213 s) and **[measured] 6.1×** at 4 threads (0.2289 s vs
0.0376 s). It also allocates **[measured] 1.6×** as much (17.1 MB vs 10.6 MB).
It did not recover the recorded performance and cannot be made to; the two
builders do different work.

**The raster path is worse than the production pair, and that is the finding
#69 flagged and nobody had measured.** Per destination cap:

| case | indexed µs/dst | latjoin µs/dst | ratio |
|---|---:|---:|---:|
| raster-4320-162chunks | 8.8 | 0.15 | **60×** |
| raster-4320-1800chunks | 18.3 | 0.51 | **36×** |
| dgg-large | 6.6 | 0.18 | **36×** |
| copdem90-igeo7-l12 | 1.84 | 0.32 | 5.7× |

#69's residual-risk note estimated ~58× per query on a 4320×2160 raster with 162
chunks; measured here at **[measured] 60.3×**, so that estimate was right. But
the ratio *falls* as the source side grows, within each index family:

- raster quadtree: 60.3× at 162 source chunks, **[measured] 35.9×** at 1 800.
- DGG hierarchy: 36.0× at 492 source chunks (`_dggcandidatechunks!`),
  **[measured] 5.7×** at 26 475 (the CopernicusDEM level-0 frontier — a
  different descent, so read this pair as a trend, not a controlled comparison).

That is the expected shape: the cap join is `O(ndst × nsrc)` before its
latitude prefilter, while an indexed query is logarithmic in the source side.
The indexed builder is the one that scales; it is losing a large constant factor
at every size measured so far. The 2.3 s serial / 0.3 s at t8 extrapolation
#69 worried about does **not** materialise on any case in the matrix.

Two things make the failed gate acceptable rather than blocking:

1. **Absolute cost.** The worst absolute number in the matrix is 0.12 s (t8) on
   a run whose recorded wall time is 8.81 h — about 4 × 10⁻⁶ of it. Space
   construction alone is 90× the whole graph build.
2. **What the extra time buys.** A relation that a refcount can be derived
   from — 72 → 0 demanded-but-unheld pairs — plus 322 fewer edges (326 064 vs
   326 386), so refcounts are tighter as well as sound.

**Recommendation.** Record the gate as *failed and waived*, not as passed. If a
later task wants the factor back, the measurements point at one place: the
`RasterGrid` arm's `_task_prepared_raster_tree(index)` cursor copy per query and
the quadtree descent it drives, which is where the 60× lives. `benchmark/
chunk_graph_gates.jl` with `DGG_GRAPH_GATE_CASES=raster-4320-162chunks` is the
one-line reproducer.

## 6. Residual uncertainty

1. **`peak_rss_growth_bytes` is uninformative once a case has already touched
   its high-water mark.** `Sys.maxrss` never falls, so the production rows read
   0 — space construction set the mark before the arms ran. The column is only
   meaningful for a case run in a fresh process. Use `graph_summarysize_bytes`
   and `graph_allocated_bytes` for graph memory; they are exact.
2. **The truth oracle is `O(ncells²)`** and therefore only runs on small spaces.
   The production pair's correctness rests on demand-domination, which is a
   weaker statement — it says the graph holds what the executor asks for, not
   what geometry requires. The two coincide only because the executor's own
   candidate relation is itself conservative; that has not been proven, only
   relied upon.
3. **The truth oracle uses `ConvexConvexSutherlandHodgman`** through
   `_intersectionoperator`, the same kernel the conservative weight builder uses.
   If that kernel is wrong on some cell shape, the oracle inherits the error —
   though in the direction that matters it inherits it *conservatively*: a
   spurious positive area only makes the oracle demand more edges.
4. **All timings predate the GeometryOps spherical-predicate/Foster-Hormann
   work.** Clipping is 24.1% of the run's time mix, so the truth oracle's own
   cost and, through `IntersectionAreaOperator`, much of the surrounding
   workload will move. The graph builders themselves touch no clipping, so the
   §3 graph timings should be stable across that change; the oracle sweep times
   will not be.
5. **The `raster-4320-*` cases use a `zeros` in-memory array,** so the quadtree
   descent is measured without the disk-array chunk machinery a real store would
   put behind it. A far larger raster with a deeper cursor was still not
   measured — this narrows #69's residual risk 1 but does not close it.
6. **The two relation-difference columns are against `:indexed`, not against
   truth.** A pair in `lat-only` is not necessarily spurious; it may be a
   genuinely contributing pair that both relations hold. Read `only_here` as
   "how the arms differ", never as "how wrong the other arm is".
7. **The nonzero-support branch measures a cell-to-cell distance, and its
   exactness rests on one geometric argument rather than on a test.** The branch
   is reached only when the clip already returned zero intersection area, and for
   two boundary rings with no positive intersection the minimum separation is
   attained at a vertex of one and a point on an edge of the other — so testing
   every vertex against every edge in both directions is exact there. It was
   previously vertex-to-vertex, which *under*-approximated the truth exactly in
   the `r > 0` cases the oracle exists to cover: two cells can lie within
   `radius` edge-to-edge with no vertex pair that close. The direction of the old
   error was safe (a weaker assertion, never a wrong one), and the fix makes the
   `r > 0` assertions strictly stronger; they still pass. What remains unproven
   is the argument itself — a ring pathological enough to break it (self-touching
   or degenerate) would make the oracle under-demand again, silently.
8. **`graph == demand` is asserted on toy, raster and DGG sources but not at
   production scale.** The CopDEM-scale equality is a construction property of
   `_fillrow!`, not a measurement; the harness's `demand_missing` column reads 0
   on `:indexed` for the same reason and is not independent evidence.

STATUS: gates landed; production graph construction untouched.
