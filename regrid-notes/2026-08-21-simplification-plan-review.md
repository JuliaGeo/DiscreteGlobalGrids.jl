# Review: regridding simplification plan (2026-08-21)

Reviewed: `regrid-notes/2026-08-21-regridding-simplification-plan.md` at branch
`claude/perf-ladder` `ce67738`. Every claim below carries a file:line; paths
are repo-relative, `GO/` is the GeometryOps checkout the root Manifest resolves
(`~/.julia/packages/GeometryOps/zXf0f`), `CR/` is ConservativeRegridding
`~/.julia/packages/ConservativeRegridding/slNrk` (tree `873cc732…` = rev
`66ed54c`, confirmed via `Base.version_slug`).

## Verdict

The direction is right and the description of the code is mostly accurate:
one-leaf chunk trees, the latitude cap join, per-read rediscovery, the
CSC→COO→CSC round trip, the destination-side wrapper stack and the two
evictors all exist as described. The plan is not executable as written,
for four reasons that are each cheap to fix in the document and expensive
to discover mid-execution:

1. Two premises are stale (the baseline story, the scheduler/executor
   disagreement) and one is true only at the correctly resolved dependency
   tree, which is not what is installed here (cached-dual-DFS redundancy).
2. Phase 1's DGG section was revised during review and now states the
   covering law correctly (`node_extent` bounds descendant geometry, not
   descendant caps — see P4, retracted); one CopernicusDEM cursor gap remains.
3. Two motivations are mis-calibrated: Phase 3 is an architecture change with
   no throughput payoff, and Phase 5 is a 1.3 GiB memory fix presented as tidiness.
4. Nothing is cited. Every number and incident is unattributed, ten test paths
   are wrong, and three gates name artifacts that do not exist.

Status 2026-08-22: the rulings on the four premises were applied to the plan
directly; see "Follow-up edits" at the end.

## Premises that do not hold

### P1. "The recorded test baselines disagree" (plan:117-123)

The 496-pass / 22-error record exists nowhere in the repo — not in any
`regrid-notes/*.md`, not in any blob reachable from any ref. The green record
is real: GlobalRegridding 2,295 / 0 / 1 broken at `merge-train-2.md:110-115`,
corroborated at `copdem-lazy-source.md:192` and `dag-driver.md:227`.

The failure the plan describes is nevertheless reproducible on this machine,
for a reason the plan does not name:

- `lib/GlobalRegridding/src/intersection_area.jl:31` calls
  `GO.intersection_area` on the live Conservative path
  (`conservative.jl:339`, `:370`).
- Project.toml pins GeometryOps at `36c853e0` (all five `[sources]` blocks,
  commit `1183e6f`), whose tree is `a7de7c62` (recorded in `7627ee3`'s
  message) and which defines `intersection_area` — `1183e6f`'s message
  records a throwaway-depot instantiate proving `isdefined(GeometryOps,
  :intersection_area)` and 2,202 / 1 broken.
- The local, gitignored `Manifest.toml` (`.gitignore:2`) carries
  `repo-rev = 36c853e0` but `git-tree-sha1 = c30ec910`, and `c30ec910` is the
  tree of `02750768` ("Bump patch version to v0.1.44"), the commit that
  *predates* `intersection_area` — verified with `git rev-parse
  02750768^{tree}` in `~/.julia/clones/18278735739139359408`. No installed
  GeometryOps tree on this machine contains the symbol; slug `AhIww` (tree
  `a7de7c62`) is absent. Runtime `isdefined(GeometryOps, :intersection_area)
  == false` from the repo root.

So: the pin is right, the local Manifest is stale. Two further facts change
what Phase 0 should do:

- `36c853e0` is now an orphan. GitHub compare `36c853e0...32c60581` reports
  the branch `claude/perf-ladder-predicates` diverged (4 ahead / 2 behind):
  it was rebased. The commit is still fetchable today; orphans are not
  guaranteed to stay so.
- `intersection_area` merged to GeometryOps main on 2026-08-21 as `6bdb983e`
  ("Add an `intersection_area` prototype (#473)"). The two perf predicates
  (`32c60581`, `f7194a35`) are the remaining unmerged delta.

The same staleness holds for ConservativeRegridding, found on the
2026-08-22 re-check: the Manifest pairs `repo-rev = 66ed54c` with tree
`873cc732`, which is `git rev-parse 6a4b997^{tree}`; the true tree of
`66ed54c` is `9469fdc5`. So `~/.julia/packages/ConservativeRegridding/slNrk`
is `6a4b997`'s code under a `66ed54c` label, which is why it has neither
`CachedDualDepthFirstSearch.jl` nor `_sparse_from_chunks` — and why every
suite count recorded "at the new pin" since `5c3960b` ran against `6a4b997`.
The three successive CR pin commits (`1183e6f` → `6a4b997`; `ea0df21` →
`3b50a2f`; `5c3960b` → `66ed54c`) reconcile the two notes: each recorded a
different moment. `66ed54c` is reachable (tip of `claude/cached-dual-dfs`)
and is the only CR revision carrying `split_weight`, sparse-CSC
`intersection_areas`, `_sparse_from_chunks` and the cached dual DFS; CR main
`89ec5ca5` lacks `split_weight` and cannot precompile GlobalRegridding.

**Change (applied to the plan, Phase 0):** repin GeometryOps to `32c60581`
(main `6bdb983e` + exactly the two perf-predicate commits; #476 open, green;
`36c853e0` confirmed orphan by `branches-where-head = []`), keep CR at
`66ed54c`, `rm Manifest.toml`, re-instantiate, assert the three `isdefined`s,
record tree SHAs in the pin commit. Verified in a clean scratchpad depot:
both resolve (trees `63274824…`, `9469fdc5…`), precompile, and all three
symbols load. Upstream `src/` is unchanged between `36c853e0` and either
candidate, so no API on GlobalRegridding's call surface moves. Not yet
executed against the working repo.

### P2. "Scheduling and execution can disagree" — ruled: `refine` is set once

With `refine === nothing` the two relations are the *same predicate over the
same caps*:

- executor: `DilatedIntersects(radius)(dstcap, srccap)` over
  `chunktree` leaf caps, `lazy.jl:585-591` → `discovery.jl:109-141`,
  predicate at `discovery.jl:47-48`;
- graph: `DilatedIntersects(radius)(dstcaps[d], srccaps[s]) && (refine ===
  nothing || refine(d, s))`, `chunkgraph.jl:444, 460-461`;
- `chunktree` leaves == `chunkextents` asserted at
  `lib/GlobalRegridding/test/test_lazy.jl:141`.

The divergence the plan describes happened only with the driver's lon/lat-box
`refine` supplied (`copdem_production.jl:722-724`), is documented as such at
`dag-driver.md:81-119` and `copdem_production.jl:640-668`, and `refinegraph =
false` is already the default (`copdem_production.jl:81`). The executor-side
safety net the plan asks for in Phase 3 — serve transiently, count, fail the
closing check — already exists at `scripts/copdem_policy.jl:227, 350-351, 381`,
asserted at `copdem_production.jl:1047`, tested at
`test/scripts/copdem_policy.jl:175-183, 359`. The adjacency-injection seam
Phase 2 asks for also exists there: the policy takes `sourcesof` and
`consumerdegree` as functions, and `test/scripts/copdem_policy.jl:333-366`
plugs a real `chunk_dependency_graph` into it.

**Ruling (applied):** `refine` is settable exactly once. The plan now has a
Phase 2 subsection "`refine` is set once": the narrow phase is a keyword of
`plan_regrid` only; `ChunkedPlan` owns the graph; `ChunkDependencyGraph`
stamps chunk counts, radius, a space fingerprint and a narrow-phase tag (the
tag, not the callable — `chunkgraph.jl:39-41`, `Base.zero` at `:209-211`);
`chunk_dependency_graph(plan; refine)` becomes `dependencies(plan)`; the
space-form builder loses `refine`. The contract after Phase 3 is
`E* ⊆ rows(graph) ⊆ capjoin(r)` with `E*` = cell geometry dilated by
`support_radius`, tested against the nonzero columns of the eager whole-domain
block, because once the executor reads graph rows a dropped edge is a silently
missing contribution that the `uncredited` counter cannot see. Production
keeps one global plan and constructs per-column plans from
`restrict(dependencies(globalplan), d)` — a rebuild would call `refine` with
the local `d`. `connectedchunks`/`connectedchunks!`/`connectedchunkpairs`
(`discovery.jl:88-154`) are removed as a third relation; `copdem_policy.jl`'s
two-function seam is unchanged; `copdem_dag_validate.jl:31, 61` stops
building the graph three times; the `uncredited` counter stays as the
structural assertion that cache and executor hold the same object.

### P3. "CR's cached dual DFS makes `MemoRasterTree`'s node-extent memo redundant" — true at the pin, false at what is installed

The first version of this review called the premise false at the pinned CR.
Re-checked 2026-08-22: it is false at the *installed* tree
(`~/.julia/packages/ConservativeRegridding/slNrk`, which has no
`CachedDualDepthFirstSearch.jl`) because that tree is `6a4b997`'s, not
`66ed54c`'s — see P1. At the true `66ed54c` tree (`9469fdc5`) the file
exists and `cached_dual_depth_first_search` loads (verified in the clean
scratchpad depot). The plan's Phase 7 paragraph therefore stands, conditional
on Phase 0's re-resolve, and was amended to say so rather than dropped.

Two adjacent facts:

- The two notes were both right about different moments: `cr-pr138-repin.md:11`
  records the pin after `ea0df21` (`3b50a2f`), `cr-cached-dfs.md:74` after
  `5c3960b` (`66ed54c`). No contradiction once the pin history is read in
  order.
- "Remeasure the DGG leaf-cap cache" was already done and is null:
  `polar-profile.md:281-315` widened `_MEMO_EXTENT_SLOTS` (sign flips across
  columns) and ported the leaf-cap memo into `MemoBlockCursor` (−6.4 %
  profile ceiling, timings inside noise); both reverted at `e2f90d1`.
  `cr-cached-dfs.md:215-218` flags `_CACHED_BUCKET_SIZE = 49` for
  re-measurement once side-2 extents are stack-cached; that is what remains.

**Change (applied):** Phase 7 now cites the file, the re-resolve precondition,
the null result, and the one measurement still open.

### P4. `node_extent` covers descendant polygons, not descendant caps — retracted as a defect (plan:180-239)

The first version of this review reported a law violation: `node_extent(parent)`
does not contain `node_extent(child)` or `cell_cap(child)` (measured
`(d + r_child) / r_parent`: S2 1.115 at L0→L1, ISEA4R 1.020, IGeo7 1.012, H3
0.99974), and recommended materialising internal extents bottom-up with
`merge_caps`. That recommendation is withdrawn. The covering law
(`src/interface/system.jl:147-174`) is a bound over every descendant cell's
*geometry* through `maxlevel`; caps are inflated bounds, not geometry, so a
hierarchy of `node_extent`s is not nested and is not meant to be. Pruning a
node by `node_extent` against a geometry-derived query is sound: it can only
drop pairs whose overlap existed between two caps' unused over-coverage, never
a pair with a geometrically possible contribution.

The plan as revised already states this correctly (plan:191-200, laws at
225-234, Phase 2 gate at 270-274, definition of done at 600-602) and traverses
`treeify(space.grid)` to the chunk frontier instead of building a `PartialGrid`
at `chunklevel`, which is the right construction: frontier nodes stay internal
nodes whose extent is `node_extent`, and the tight `cell_cap` appears only at
literal cell leaves, where there is nothing below to cover. The only gate that
follows is the one the plan now has — no geometric false negatives, with the
flat cap join as an upper bound — not identity with the cap join.

Because three reviewers in a row misread this, the distinction is now stated
as a contract in the `node_extent`, `cell_cap`, `cells_cap`, `merge_caps` and
cursor docstrings (see "Follow-up edits" below), so that future readers cannot
take `node_extent` for a cap over caps.

Adapter facts that still apply to the revised Phase 1 (verified empirically
on IGeo7/S2/CopDEM spaces); the first two are already handled by traversing
the original tree rather than a new `PartialGrid`, and are kept here as the
reason that choice is necessary:

- the cursor's leaf `child_indices_extents` returns tight `cell_cap`s
  (`src/engine/cursor.jl:265-275`); IGeo7 L0 `cell_cap` radius 0.652 vs
  `node_extent` 0.783, a 17 % shortfall that does not cover aperture-7
  overhang. The adapter must override leaf extents, not only `node_extent`.
  Because `chunkextents(::DGGSpace) = space.caps` (`src/regridding.jl:161`),
  the shortfall is invisible to the existing `capscover` test and would only
  surface in live queries — keep that override after cutover.
- sparse tightening is real and would prune valid chunks: 320 internal nodes
  tighten on a chunklevel-3 IGeo7 subset (`cursor.jl:232-257`); the plan's
  bypass is necessary, not optional.
- `treeify` returns a `BlockCursor`/`MemoBlockCursor` for CopernicusDEM
  complete grids and lattice-rectangle `PartialGrid`s
  (`src/systems/CopernicusDEM/cursor.jl:348-373`); CopDEM has sorted subtrees
  (`system.jl:12`) so it *is* chunked and the plan's fallback list ("systems
  without sorted subtrees" = A5 only, `src/systems/A5/system.jl:115`) does not
  cover it. Existing code guards this at `src/regridding.jl:191-193` and
  `src/cap_cached_tree.jl:123-125`.
- root handling is inverted in practice: the auto `chunklevel` is chosen from
  cell counts independent of the root, so `DGGSpace(subtree(IGeo7, L2, 4))`
  gives `chunklevel = 0`, `nchunks = 1`, and passing `root` to `PartialGrid`
  throws (`partial_grid.jl:87-96`). Guard on `grid.root_level <=
  space.chunklevel` (mirror `src/regridding.jl:106-115`); build the trivial
  index off `nchunks == 1`, not `chunklevel < 0`; guard on `_ischunked`, not
  `isempty(chunkids)` (`_wholechunk` stores `ID[]` with `ID = Any` when there
  is no system, `src/regridding.jl:50-56`, which fails `PartialGrid`'s eltype
  check at `partial_grid.jl:75-76`).

## Per-phase findings

### Phase 0

Covered by P1 and P3. Add to the phase: create the artifacts the later gates
assume (see Hygiene H3) — they do not exist and cannot be "compared against".

### Phase 1 — Raster / generic (`FlexibleRTrees`)

Every API premise holds at the pin, in the resolved tree:

- `GO/src/utils/FlexibleRTrees/types.jl:78-101` `RTree(HPR(), items;
  extents, indices)`, no `GI.extent` call when `extents` is passed;
  `leafindices = indices[perm]` at `:99` preserves original indices.
- `GO/src/utils/UnitSpherical/cap.jl:131` and `:141-142` —
  `Extents.intersects(cap, ::Extent{(:X,:Y,:Z)})` in both orders, "may have
  false positives but never false negatives" (`:129`).
- No `Extents.extent(::SphericalCap)`; no cap→XYZ-box helper anywhere in
  GO/GOCore/UnitSpherical. `cap_xyz_extent` is genuinely new.
- FlexibleRTrees has been in GO since 0.1.42 (`f57277b08`, #431); compat
  `"0.1.43"` already admits it. No pin move needed for this phase.

Corrections to the write-up:

- `extents` must be a concrete `Vector{<:Extents.Extent}` and the tree takes
  ownership (`types.jl:79, 91-93`).
- A bare cap is not a predicate: `STI.query(tree, cap)` throws
  (`sanitize_predicate`, `SpatialTreeInterface.jl:67-70`). Use
  `depth_first_search(Base.Fix1(Extents.intersects, cap), tree)`.
- Unqualified `GeometryOps.query` is an `UndefVarError`: `SpatialTreeInterface`
  and `FlexibleRTrees` export *different* `query`s.
- `Extents.extent(cap)` silently returns `nothing` via Extents' universal
  fallback (`Extents.jl:124`) — worth an explicit negative test.
- `Extents.intersects` dispatches on the exact `(:X, :Y, :Z)` name order.
- Singleton collections work; only the empty case needs a fallback
  (`types.jl:86`). Delete "and singleton" at plan:174.
- Radius ≥ π is already a whole-sphere query: `chord = 2 sin(½ min(r, π))`,
  `cap.jl:137`. plan:173 can cite rather than require.
- CR's spherical broad phase is hard-wired to cap–cap
  (`CR/src/regridder/intersection_areas.jl:259-263`), so the box R-tree can
  only ever be the private chunk index — the plan says this; it is worth
  saying *why*.

### Phase 2 — graph

Confirmed: `_caplatitude` `chunkgraph.jl:376`, sorted lats `:386-388`, global
band `:389`, `_fillrow!` `:440-467`, parallel over destination blocks
`:394-408`, counting-sort transpose `:470-492`. The band is widened globally by
the single largest source cap (`:445`). Already done, which the plan asks for
as future work: Int32 edge indices with a range check
(`chunkgraph.jl:380-383`). Already existing, which the plan could reuse:
`prefilter = false` is the documented brute-force reference
(`chunkgraph.jl:305-306`, exercised at `test_chunkgraph.jl:28`).

Wrong: "the public `connectedchunks` query" (plan:260) — it is neither
exported nor `public` (`GlobalRegridding.jl:61-97`). `chunktree` *is*
exported (`:62`).

Gate scope: the 0.122 s / 326,392 edges / 26,475 × 66,178 figures are real
(`copdem-lazy-source.md:110-113`, full graph, 4 threads) but are 0.12 s of a
~1 s plan build "almost all of it the destination space"
(`copdem_production.jl:701-702`). A gate on 0.122 s measures 12 % of planning.
No note records graph bytes — plan:130 and plan:312 are new measurements, not
comparisons; `chunk_dag_poc.jl:351`'s `adjacency_bytes` is NDJSON file size.

Refine deferral costs nothing the notes can price: refined build is slower
(0.172 vs 0.122 s, `copdem-lazy-source.md:112-113`), saves zero downloads and
zero GB (`:112-118, 171-173`), and broke at-most-once (`dag-driver.md:85-97`).
The 1.43 GiB residency figure at `dag-driver.md:136` is for the *exact*
186,069-edge adjacency, not the refined graph; residency under `refine =
true` was never measured.

Per-column waste the plan names but does not price (plan:264-267): each of
66,178 per-column lazy regrids (`copdem_production.jl:781-789, 981`) runs
`foldl(merge_caps, space.caps)` over 26,475 caps (`src/regridding.jl:210-212`)
and allocates a 26,475-byte `emptymemo` (`lazy.jl:268`); `srccaps` is shared
by reference, not copied.

### Phase 3 — executor

Confirmed: `DestTiling.capsof` `lazy.jl:66-70, 86, 91`; `knownempty`
`executor.jl:30-37` applied per (chunk, non-spatial group) `lazy.jl:434` and
memoised `:601-612`; `othergroups` `:236, 361-369`; ascending source order
`discovery.jl:112-113`, `lazy.jl:426`. Everything named "executor" lives in
`lazy.jl`; `executor.jl` is accumulation/finalisation/shaping.

`_connectedsource!` (`lazy.jl:585-598`) is called per tile per `readblock!`
(`:405, 415`) with nothing cached across reads except emptiness. But the
"dual-tree query" is a 1 × nchunks linear scan: `CapQuery` is a one-entry
leaf (`discovery.jl:16-26`), every current `chunktree` is a one-entry leaf,
so `dual_depth_first_search` takes the leaf×leaf branch
(`GO/src/utils/SpatialTreeInterface/dual_depth_first_search.jl:46-56`).

The payoff is architectural, not throughput. `polar-profile.md:124-136`
measured discovery at 0.5 ms per column at 360 tiles ("0.0005 % of the
column"), plan build at 2.0 ms, against ~10 core-s columns: ~33 s over 187
core-hours. The exit criterion "no chunk-tree search after plan build" reads
as a performance claim the notes contradict. State the real win: one relation
owned by the plan, consumed by scheduler and executor.

`_wavesize`/`_blockcosts!` (`lazy.jl:469-485, 502-518`) read both cap vectors,
so `srccaps` is not retained "solely for per-read discovery"; Phase 4's
caveat is correct and Phase 3's removal list should be only `srctree`.

### Phase 4

Generic `chunkextent(space, chunk)` is `chunkextents(space)[chunk]`
(`discovery.jl:86`): O(nchunks) per call on `RasterGrid`
(`rastergrid.jl:912-919`), and `connectedchunks` pays it every call
(`discovery.jl:94-96`). "Keep cheap leaf extent access" describes an API that
exists but is not cheap; making it cheap is new work.

`RasterFlatTree`'s scattered-cell role is real (`rastergrid.jl:907-910,
928-931`; tests `test_rastergrid.jl:357`, `test_integration.jl:294`) and
separable from the chunk-lattice constructor (`:912-919`).

`DGGChunkTree`'s `Trees.ncells/getcell` (`src/regridding.jl:226-229`) are
dead — `chunktree` is consumed only by `discovery.jl` and forwarded by
`TileCells` (`conservative.jl:157-158`); nothing asks it for cell geometry.

### Phase 5 — `WeightBlock`

Confirmed round trip: `buildblock` (`plans.jl:419-425`) → `WeightCOO` →
`build_weights!(::Conservative)` (`conservative.jl:328-346`) → CSC at `:342`
→ `_fillcoo!` (`:406-416`) → `WeightBlock(coo, …)` (`plans.jl:26-29`) calls
`sparse` a second time. Fast path `wholeblock(::Conservative)`
(`conservative.jl:357-377`) adopts the CSC and computes denominators via
`_blockdenom` (`:381-391`). Bit-identity between the two is already tested
(`test_conservative.jl:385-397`); partition invariance at `:89-114` and
`test_lazy.jl:691-717`.

Stronger justification than the plan gives: `docs/design/regrid-notes/perf-P1.md:354-368`
("P4") measures the round trip at 0.17 % of time but **1335 MiB `maxrss`
delta, 8.4× the result**, from CR matrix + COO + final matrix being live
together — on the chunked path, which is the one with no fast route. That
file is gitignored (`.gitignore:9` `docs/*`): the number has no committed
artifact and should be re-measured in Phase 0.

Facts to carry into the design:

- Phase 5 step 3 already exists: `WeightBlock(coo::WeightCOO, ndst, nsrc)`
  (`plans.jl:26-29`). The seam is a dispatch wrapper.
- Required method for a new regridder is `build_weights!` only
  (`methods.jl:130`); a `weightblock` with generic fallback is compatible.
- Type-dispatch overrides are lost by wrapper methods:
  `T6CountingMethod(Conservative())` at `test_integration.jl:133` and the
  `P1CountingMethod` harness already miss `wholeblock`; Phase 5 extends that
  to the chunked path, so counter-instrumented benchmarks measure the slow
  path. Use a trait or a forwarding requirement.
- `wholeblock(::Conservative)` `invoke`s the generic path when either side is
  empty (`conservative.jl:360-363`) so `hasdenom` matches; a `weightblock`
  override must too. There is an untested asymmetry: empty `dst_inds` returns
  before `markdenominated!` (`:330`), empty `src_inds` after (`:331`).
- An alternative the plan should adjudicate: `BlockAreaOperator` already
  holds both index maps (`conservative.jl:280-286`); pushing straight into the
  `WeightCOO` (`perf-P1.md:364-366`) removes the CSC instead of adopting it and
  is a smaller change than a new seam.
- `wholeblock` is private (`GlobalRegridding.jl:60-99`); its removal is
  internal.
- Denominators are accumulated positive intersected area per block
  (`_fillcoo!` `:411-413`, `_blockdenom` `:385-389`), summed by the executor
  (`lazy.jl:429`, `executor.jl:432`) and thresholded in `finalize!`
  (`executor.jl:305-320`). Per-block denominators are the only correct
  design; the plan's law is well-posed.

### Phase 6 — DGG output

Confirmed: `_DirectToDGG` `src/regridding.jl:270-271`, `_ChunkedToDGG`
`:272-273`, the two `GR.regrid` methods `:275-280`, `_ascube` `:283-291`.
The eager method is pure duplication (it `invoke`s the generic path, which
already labelled axis 1 via `wrapoutput`, then relabels).
`destinationdims(space, sampling)` exists with that signature
(`spaces.jl:132-138`, `public` at `GlobalRegridding.jl:89`) and is consumed by
`wrapoutput` (`executor.jl:397-410`), `wraplazy` (`lazy.jl:835-861`) and both
`regrid!` shape checks.

Measured with the proposed method (NearestCell, since Conservative is
un-runnable locally — P1): eager types/dims/values identical; lazy dims/values
identical, but `parent` changes `LazyRegridArray → ShapedRegridArray`
(`lazy.jl:843`). So the phase is not deletion-only: `wraplazy` needs a
1-tuple short-circuit inside GlobalRegridding, and no existing test would
catch the type change (`regrid_acceptance.jl:96`, `regrid.jl:362-365` check
only `Cells` and `parent` behaviour). Also fix the now-false message at
`api.jl:184-186` ("a lazy regrid returns an unlabelled disk array").

### Phase 7 — caches

The destination-side overlap is real: with prebuilt polygons,
`_memocell(op.dstmemo, dst_tree, i2)` (`conservative.jl:310`) memoises what
`CachedCellTree.getcell` (`:228-232`) already reads from an array; on DGG
destinations the stack is `CachedCellTree(CapCachedTree(cursor))` + `CellMemo`
(`src/regridding.jl:176-183`). The source-side `CellMemo` earns its keep, as
the plan says.

Required laws, current status: restricted tree built once — tested
(`test_conservative.jl:312-343`, `test_lazy.jl:447-468`); cached/uncached
bit-identical — tested (`test_conservative.jl:265-271, 337-342`); no shared
mutation — by construction (`ReentrantLock` `conservative.jl:176-186`),
untested; tile threshold `_TILE_CELL_CACHE_MAX` (`:119`) — appears nowhere in
`test/`; production Conservative benchmark — does not exist in-repo
(`benchmark/` holds toys/geomorphometry/hex_border/maxneighbors; `scripts/bench/`
holds `halo_subset_scaling.jl`).

"Measured performance justification" is half true: `CapCachedTree` carries
its sweep in-tree (`src/cap_cached_tree.jl:78-108`); `MemoRasterTree`'s
docstring (`raster_tree_memo.jl:3-13`) is an argument, and its only
measurements are in gitignored `docs/design/regrid-notes/`. See P3 for the
cached-DFS paragraph.

### Phase 8 — drop it

Genuine duplication between `_evict!` (`plans.jl:212-229`) and
`_evictoldest!` (`lazy.jl:198-211`) is ~11 lines (scan for min stamp,
subtract bytes, delete). Everything else differs: lock vs none, count+bytes
vs bytes-only, protected `keep` key vs none, insert-then-retain-newest vs
refuse-oversized (`lazy.jl:187`), struct field vs parallel `Dict`. A shared
helper needs more parameters than the lines it removes. Additionally
`Spilled` wraps a `PerChunk` (`plans.jl:239-251`) and the open "Spilled
fingerprint redesign" appears in no note; touching `PerChunk` without
naming it is a scope leak.

## Calibration summary

| plan says | evidence says |
|---|---|
| lazy/graph "can disagree" (P2) | identical predicate unless `refine` supplied; already off; counter + closing check already implemented |
| Phase 3 removes a search cost | 0.5 ms/column of ~10 core-s (`polar-profile.md:124-136`); win is one owned relation |
| Phase 5 avoids a round trip | 1335 MiB peak-RSS on the chunked path (`perf-P1.md:354-368`, uncommitted) |
| Phase 2 gate on 0.122 s | that is 12 % of a ~1 s plan build (`copdem_production.jl:701-702`) |
| graph bytes "baseline" | no note records it; new measurement |
| deferring `refine` has a cost | never measured; refined graph is slower to build and saves nothing |
| compact edge indices "where counts permit" | already Int32 (`chunkgraph.jl:380-383`) |
| `connectedchunks` is public | it is not; `chunktree` is |
| Phase 7 "remeasure leaf-cap cache" | already measured null (`polar-profile.md:281-315`) |
| Phase 8 shares "substantially duplicated" eviction | ~11 lines |

## Hygiene

- **H1 Provenance.** The plan cites no note, script, or line. Every number
  above has one; carry them into the plan so a reviewer audits instead of
  re-deriving.
- **H2 Paths.** All ten `test/test_*.jl` entries in the verification matrix
  (plan:524-532) live under `lib/GlobalRegridding/test/`. The root-level
  paths do not exist.
- **H3 Gates without artifacts.** Graph bytes (Phase 2), tile-threshold bound
  and a production Conservative benchmark (Phase 7), and the P4 memory figure
  (Phase 5) have no committed artifact. Phase 0 must create them before any
  phase can claim "no regression".
- **H4 The "1 broken"** is `test_conservative.jl:147`, the GeometryOps
  non-convex destination clip, with the same arm in
  `regridding_conservation.jl:86-99` (three `@test_broken` per non-convex
  system). Phase 5's "denominators identical including zero-coverage rows" is
  asserted over a path with a known non-conserving arm; say so.
- **H5 Pole-artifact exclusions.** Eleven production columns are hybrid-source
  (`pole-artifact.md:196-208`: 122975, 122977, 122979, 122981, 123147, 123172,
  123176, 123199, 123202, 123203, 123204) and are not valid byte-identity
  references for any production-subset gate.
- **H6 Memory gates must name their shape.** `shape = outer → inner` is the
  ranked-#2 recommendation and the dominant RSS term
  (`polar-profile.md:387-391`, 23 × 1 GiB → 8 × 1 GiB per-worker budget); the
  unexplained epoch-C 4.5× collapse (`polar-profile.md:373-386`) bounds any
  "no material regression" claim measured on the shared box.
- **H7 Upstream decisions in flight.** GO PR #474 / CR cached-DFS ownership
  (`cr-cached-dfs.md:171-188, 225-226`) decides whether Phase 7's memo
  paragraph ever becomes true; #473 is now merged (P1).
- **H8 Scripts the matrix omits** but which cover its concerns:
  `test/scripts/copdem_source_mode.jl` (synthetic-source absoluteness) and
  `scripts/copdem_prefetch_coldtest.jl:247-248` (the only check that every
  demand was a graph edge against a real network).

## Recommended edits, in order

1. Rewrite plan:117-123 and Phase 0 per P1; add "re-instantiate, assert
   `intersection_area`, repin to a reachable commit, record tree SHAs".
2. Rewrite plan:66-69 and 284-288 per P2; cite `copdem_policy.jl` as the
   existing seam and counter.
3. Fix Phase 1's DGG law/prescription conflict per P4; add the four adapter
   facts (leaf override, tightening bypass, BlockCursor, root guard).
4. Apply the Phase 1 Raster API corrections (predicate wrapping, qualified
   `query`, vector ownership, drop singleton fallback, cite `cap.jl:137`).
5. Restate Phase 3's payoff as architectural; keep `srccaps`, remove only
   `srctree`; keep the uncredited-demand counter as the acceptance mechanism.
6. Widen the Phase 2 gate to whole plan build; mark graph bytes as a new
   measurement; add the per-column fold cost.
7. Promote Phase 5's motivation to memory; adjudicate push-into-COO vs adopt
   CSC; add the wrapper-method dispatch rule and the empty-side contract.
8. Phase 6: add the `wraplazy` 1-tuple short-circuit and a type assertion;
   fix `api.jl:184-186`.
9. Phase 7: delete the cached-DFS paragraph; cite the null c1/c4 result;
   create the tile-threshold test and the production benchmark in Phase 0.
10. Delete Phase 8; reference the Spilled redesign as out of scope.
11. Fix the ten matrix paths; add provenance throughout; add H4–H8.

## Method

Five parallel read-only investigations (discovery/graph/executor;
GeometryOps/CR API at the pin; DGG hierarchy and output wrappers; weight
construction and caches; session notes, scripts and git history), each
returning file:line verdicts, followed by direct checks of the Manifest tree
hash against the GeometryOps clone, the pin-moving commits `acf6909`,
`7627ee3`, `1183e6f`, and the GitHub compare for `36c853e0`. Empirical
DGG measurements (leaf-index identity, tightening counts, cap ratios,
Phase 6 dims/types) were run with `julia --project` on this checkout using
`NearestCell`, because Conservative cannot load here (P1).

## Follow-up edits (2026-08-22)

Rulings received on the four premises: repin as things have merged; `refine`
settable once; re-check the cached DFS; `node_extent` covers descendant
polygons, not caps, and that must be unmissable in the docs. Applied:

**Plan** (`2026-08-21-regridding-simplification-plan.md`, untracked):

- Current state: the disagreement paragraph restated as a property of a
  narrow phase with two owners; the baselines section replaced by the
  Manifest-staleness account for both sources and the 2026-08-22 upstream
  state.
- Phase 0 rewritten: repin GeometryOps → `32c60581`, keep CR `66ed54c`,
  `rm Manifest.toml`, re-instantiate, three `isdefined`s, tree SHAs in the
  commit message, reference counts from `1183e6f`, four baseline artifacts.
- Phase 1: Raster API details (vector ownership, predicate wrapping,
  qualified `query`, silent `Extents.extent`, singleton needs no fallback,
  radius ≥ π already handled, CR's cap–cap broad phase); DGG: CopernicusDEM
  `BlockCursor` caveat, `nchunks == 1`/`_ischunked` dispatch, why tightening
  is sound on the original tree and unsound on a `PartialGrid` at
  `chunklevel`.
- Phase 2: gate widened to whole plan build; bytes marked as a new
  measurement; per-column fold cost priced; `Int32` edges noted as existing;
  `connectedchunks` demoted to a third relation; row-view `restrict` chosen
  with the local-`d` trap; new subsection "`refine` is set once".
- Phase 3: payoff restated as architectural (0.5 ms/column); caps retained
  until Phase 4; the post-Phase-3 `refine` contract, its weight-block test,
  the 186,069-edge certificate as a runnable test, the `uncredited` counter's
  new meaning, and the driver consequences.
- Phase 4: `chunkextent` is O(nchunks) today; removals of
  `connectedchunks*`, the plan-form graph entry, `DGGChunkTree`'s dead
  `getcell`; positional `_chunkgraph` call in `test/scripts/copdem_policy.jl`.
- Phase 5: memory motivation (1335 MiB, uncommitted source), existing
  finalizer, wrapper-dispatch rule, empty-side contract, push-into-COO
  alternative, existing tests cited per law, RSS exit criterion.
- Phase 6: `wraplazy` one-axis short-circuit, `parent`-type assertion,
  `api.jl:184-186` message.
- Phase 7: cached-DFS paragraph kept and conditioned on the re-resolve; null
  c1/c4 result cited; `_CACHED_BUCKET_SIZE` as the open measurement; law
  status and missing artifacts listed.
- Phase 8 dropped; scope bullet and commit item removed.
- Verification matrix: ten paths fixed to `lib/GlobalRegridding/test/`; two
  scripts added; gate hygiene (the broken test, pole exclusions, `shape`,
  provenance rule); standing law 7 and definition of done updated.

**Docstrings** (uncommitted, docstrings and comments only; package loads):
`node_extent` and `cap_inflation` in `src/interface/system.jl`; `cell_cap`,
`cells_cap`, `merge_caps` in `src/fallbacks/caps.jl`; the default
`node_extent` in `src/fallbacks/geometry.jl`; `HierarchicalGridCursor`,
`STI.node_extent(cursor)` and `child_indices_extents` in
`src/engine/cursor.jl`; the halo soundness comment in `src/engine/halo.jl`;
`PositionTree` in `src/engine/position_tree.jl`; the trait table row in
`src/interface/types.jl`; and the per-system comments in S2, ISEA4R, HEALPix,
CopernicusDEM, H3 and IGeo7. Each states the contract — contains every
descendant cell's geometry; not required to contain descendants' caps or
`node_extent`s, and for real systems does not — and what pruning may assume.
No numeric ratios were added, since no committed artifact reproduces them.

**Not executed:** the repin itself (five `Project.toml` edits, `rm
Manifest.toml`, instantiate, suites). The recipe and its clean-depot
verification are in Phase 0.
