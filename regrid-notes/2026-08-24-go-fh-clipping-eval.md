# Evaluating the GeometryOps spherical Foster-Hormann branch as a DGG pin

- Date: 2026-08-24
- Branch: `claude/go-fh-clipping-eval`, cut from `claude/e2-delete-legacy-discovery`
  @ `04ef77b` (PR #75's head). `claude/perf-ladder` was not touched.
- Machine: 64-core shared box, Julia 1.12.6, `--gcthreads=4` everywhere.
- Verdict: **GO-WITH-CAVEATS as a correctness-neutral, off-driver improvement.
  NO-GO as a clipping optimisation** — the hypothesis this evaluation was set up
  to test is disproven by the code, and then by measurement.

## 0. Headline

**DiscreteGlobalGrids does not use Foster-Hormann clipping.** Its conservative
weight kernel is
`GO.intersection_area(GO.ConvexConvexSutherlandHodgman(manifold), …)`
(`lib/GlobalRegridding/src/intersection_area.jl:31`). The branch does not touch
`src/methods/clipping/sutherland_hodgman.jl` — it is byte-identical to the base
pin — and `spherical_orient`, the only predicate that kernel calls, is
byte-identical too. So the 24.1% `IntersectionAreaOperator` bucket in
`2026-08-23-full-run-attribution.md` **cannot move**, and it did not: both block
builders are within +1.3 %, allocating byte-identically, and every weight and
every regridded value is bit-identical. (One micro arm is a reproducible +8 % on
byte-identical source; §5 records it and does not explain it away.)

What the branch *does* buy DGG is on two other surfaces, both real and both
large: the RelateNG region-query path (`DGG.query`, −57% to −63% on `Within`)
and the generic `point_in_cell` / `cellat` fallback (−74% to −79% wall,
640 B/call → **0 B/call**). Neither is in the production driver's hot loop.

Separately, the evaluation found a **bug in the branch** (§7): spherical
`FosterHormannClipping` `intersection_area` returns a whole polygon's area for
two rings that only touch along a shared edge. It is reported with a
self-contained reproducer; nothing was changed in the GeometryOps tree.

## 1. The pin — the requested SHA no longer exists

The card asked for `as/spherical-fh-clipping` tip `8af92af81`, described as 13
commits ahead of the current pin `697c7cc8` with `697c7cc8` as an ancestor.

That is no longer the state of the world. `git ls-remote` shows the branch has
been **rebased and force-pushed**; `8af92af81` is not on the remote and cannot
be fetched (`Pkg.update` fails with `Did not find rev 8af92af81… in repository`,
and a stale remote-tracking ref in a local clone is what makes it *look*
present). The card's ancestry claim no longer holds either.

**Pinned SHA: `35996798bcc34d325452cfa73926b37918ece65c`**
(manifest tree `1e8708fd530e62235dc711da64449e9c2df82ba0`), the live tip as of
this session. Same 13 commit subjects, rebased onto GeometryOps `main`
(`3b1878db0`).

Consequences of the rebase, stated because they widen the blast radius the card
assumed:

- `697c7cc8` is **not** an ancestor of the new tip. This is not a fast-forward.
- The effective diff is **26 commits**, not 13. Beyond the FH/predicate work it
  now carries `main`'s resolution of the spherical-cap primitives
  (`UnitSpherical/cap.jl`, `methods/extent.jl`) — the very files the card
  predicted would not conflict — plus GO PR #472
  (`relateng/relate_ng.jl`, `relate_geometry.jl`, "Reduce the per-call B-side
  cost of prepared `RelateNG` queries") and `RobustCrossProduct.jl`.
- Some of the observed query win is therefore #472's, not the FH branch's. The
  measurement cannot separate them, and this record does not claim it can.

### 1.1 There are five pin sites, not four

The card names four. `test/Project.toml:29` is a fifth and carries the same
`rev`. Leaving it would have run every non-GlobalRegridding suite against the
old GeometryOps. All five were bumped:

| file | line |
|---|---|
| `Project.toml` | 33 |
| `lib/GlobalRegridding/Project.toml` | 29 |
| `lib/DiscreteGlobalGridsConformanceTesting/Project.toml` | 20 |
| `docs/Project.toml` | 35 |
| `test/Project.toml` | 29 |

`benchmark/` and `lib/DiscreteGlobalGridsVisualization/` carry no GeometryOps
source; they inherit through the workspace manifest.

### 1.2 Resolution

No version conflict. GeometryOps stays `v0.1.44`; the branch changes neither
`Project.toml` nor `GeometryOpsCore/`, so the `GeometryOpsCore = "main"` pin is
**unaffected** — its manifest tree is `0896ad28d34338853bd9d335fee15053d748dc58`
before and after, and ConservativeRegridding stays at `e8fc67f4`. The two
comparison arms therefore differ in exactly one dependency.

One trap worth recording: `Pkg.resolve()` rewrites `repo-rev` but **leaves
`git-tree-sha1` pointing at the old tree**, so the old source stays loaded and
the "after" run silently measures the "before" code. `Pkg.update("GeometryOps")`
is what actually moves the tree. Caught by asserting
`isdefined(GO.UnitSpherical, :exact_spherical_orient)`.

The baseline manifest was also mildly stale against the branch point's
`Project.toml`s (`Pkg` warned; the resolve added `Chairmarks` and a
`DimensionalData` weakdep to the visualization ext). Neither is loaded by
anything measured here, and the GeometryOps / GeometryOpsCore /
ConservativeRegridding entries were correct in the baseline manifest, so no
"before" number is affected.

## 2. Suites — re-measured at the branch point, then after

`--gcthreads=4`; `-t 4` for GlobalRegridding, `-t 8` for the rest.

| suite | before (re-measured at `04ef77b`) | after (`35996798b`) | delta |
|---|---|---|---|
| `lib/GlobalRegridding/test` | 4 021 pass / 1 broken / 0 fail | **4 021 / 1 / 0** | none |
| `test/systems/crosssystem/regrid.jl` | 236 / 0 | **236 / 0** | none |
| `test/systems/crosssystem/regrid_acceptance.jl` | 22 / 0 | **22 / 0** | none |
| `test/scripts/copdem_policy.jl` | 89 / 0 | **89 / 0** | none |
| `test/scripts/copdem_source_mode.jl` | 7 / 0 | **7 / 0** | none |
| `test/systems/CopernicusDEM/runtests.jl` | 16 258 / 3 broken / 0 fail | **16 258 / 3 / 0** | none |

Every "before" column reproduces E2's reported baseline exactly. **No delta in
any suite, so there is nothing to explain.** No failures, no errors; the one
pre-existing broken GlobalRegridding test and the three pre-existing broken
CopernicusDEM conformance skips are unchanged in both columns.

One operational note: the GlobalRegridding suite must be run with
`--project=lib/GlobalRegridding/test`, not `--project=lib/GlobalRegridding` —
the latter cannot see `Proj` and errors out at 583 assertions.

## 3. Do values move? No — bit-identical, everywhere measured

A dedicated script dumps weights and regridded fields as **raw `Float64` bytes**
and compares them with `cmp`, not with a tolerance. Nine dumps, all cases:

| case | quantity | before vs after |
|---|---|---|
| eager whole-domain block, 1° raster → IGeo7 L4 | 184 076 nonzeros | **byte-identical** |
| | CSC `rowval`, `colptr` | **byte-identical** |
| | denominator vector (24 012) | **byte-identical** |
| | Σ weights | `12.566370614359169` both |
| regrid, 0.5° raster → IGeo7 L4, eager | 24 012 values | **byte-identical** |
| regrid, same, lazy (4 MiB budget) | 24 012 values | **byte-identical** |
| **pentagons** (the 12 five-vertex L4 cells), eager and lazy | 12 values + areas | **byte-identical** |
| **polar** — the south-pole chunk of a 0.1° tiled DEM onto IGeo7 L7, lazy | 2 401 values, 406 NaN | **byte-identical** |
| polar, eager reference over the reachable latitude band | 2 401 values, 406 NaN | **byte-identical** |

Magnitudes, stated as the card asks — of the change in the values themselves:

- `max |after − before|` = **0.0 exactly**, on every one of the nine dumps, and
  the SHA-256 of every dump matches.
- `max |eager − lazy|` = **0.0** before and **0.0** after. (This case is
  stricter than E1's 2.44249e-15; that figure came from a *differently tiled*
  720×360 → L4 plan with level-2 chunks. The point stands either way: the
  eager/lazy gap did not move.)
- polar `max |lazy − eager|` = **0.0** absolute and **0.0** relative, before and
  after, over the 1 995 non-NaN cells; the NaN count is 406 in all four runs.

No shift at all, not even at the ulp level. That is the expected result once
§0's code argument is granted: the clipper and its predicate are byte-identical
source, so the arithmetic cannot differ.

## 4. Gates harness — flat, as predicted

`benchmark/chunk_graph_gates.jl`, `-t 8 --gcthreads=4`, 5 samples, 13 cases
including the production `copdem90-igeo7-l12` pair from a **local** tile list
(`bench/data/CopernicusDEM/tileList-glo90.txt`, 26 475 tiles; nothing was
downloaded). 26 ndjson rows each run, self-stamped:
GeometryOps `697c7cc8`/tree `9944edf3` → `35996798b`/tree `1e8708fd`,
GeometryOpsCore tree `0896ad28` unchanged, ConservativeRegridding `e8fc67f4`
unchanged.

**Every relation field is identical on all 26 rows**, compared field by field:
`edges`, `demanded_pairs`, `demand_missing`, `oracle_pairs`, `oracle_missing`,
`only_here`, `missing_here`, `destination_chunks`, `source_chunks`, `radius`,
`identity_bytes`, `graph_allocated_bytes`, `graph_summarysize_bytes`. That
includes the production pair's 326 064 indexed / 326 386 latjoin edges and the
72 / 394 crossing. `graph_allocated_bytes` is byte-identical on every row. Both
runs print the same verdict, character for character:

```
cases run: 13.  Oracle-checked: 9.  Oracle skipped: 4.
verdict: PASS on 9 oracle-checked case(s); 4 case(s) unchecked
```

Timings, the two cases large enough to read:

| case | arm | before | after | delta |
|---|---|---:|---:|---:|
| `copdem90-igeo7-l12` | indexed | 0.113268 s | 0.111151 s | **−1.9 %** |
| `copdem90-igeo7-l12` | latjoin | 0.021335 s | 0.020900 s | **−2.0 %** |
| `raster-4320-1800chunks` | indexed | 0.060275 s | 0.060042 s | −0.4 % |
| `dgg-large` | indexed | 0.032867 s | 0.033111 s | +0.7 % |
| `raster-4320-162chunks` | indexed | 0.034943 s | 0.038972 s | +11.5 % |

The sub-millisecond cases swing between **−23.6 % and +50.8 %** in both
directions on both arms — the same shared-box noise E2, E1, G3 and G4 all
recorded on unchanged code, a little wider today because the box was busier.
The graph builders touch no clipping, so flat is the expected and observed
result; nothing here is a win or a regression.

`space_seconds` for the production case is **+1.4 %** (7.210 → 7.307 s). That
number matters because it contains `covering_chunks`, i.e. 26 475
`DGG.query(sys7, MultiOrderCoverage(extent); level = 5)` calls — an extent
target, not a polygon target, and therefore *not* on the RelateNG path §6
measures. The production planner does not collect the query win.

## 5. The clipping measurement

The gates harness does not exercise clipping, so this is a purpose-built
harness. Method:

1. Build the production `BlockAreaOperator` exactly as `build_weights!` does,
   but with a **recording inner operator** that returns 0 (so the block
   assembles nothing) and pushes the two geometries it is handed. Run under
   `@with GR.OUTER_PARALLEL => true`, which forces `_innerthreaded()` serial so
   the recorder is race-free. This harvests the *real* candidate pairs, in the
   real representation, straight off the production code path.
2. Stride-subsample to ≤ 40 000 pairs (keeps the geographic mix).
3. Time `GR.IntersectionAreaOperator(manifold)` over the pair list, `GC.gc()`
   before each sample, `@timed` for wall / bytes / GC.
4. Separately time the two block builders end to end.
5. Every arm carries a value fingerprint (the summed area, or `nnz` plus
   Σ weights) so a "faster" arm that moved a number would be visible in its own
   row.

`-t 8 --gcthreads=4`. The whole "after" arm was run **twice**, as a noise floor.

| arm | pairs (candidates) | reps | before min | after min | after₂ min | Δ vs before | after vs after₂ | bytes | GC frac |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `micro.global` 1° → L4 | 37 647 (263 524) | 15 | 0.031130 s | 0.033621 s | 0.033220 s | **+8.0 %** | −1.2 % | 80 B both | 0.000 |
| `micro.pentagon` 12 five-vertex cells | 137 (137) | 51 | 85.75 µs | 86.83 µs | 87.13 µs | +1.3 % | +0.3 % | 80 B both | 0.000 |
| `micro.polar` south-pole column, 0.1° src | 38 514 (269 598) | 15 | 0.029190 s | 0.029528 s | 0.028636 s | +1.2 % | −3.0 % | 80 B both | 0.000 |
| `block.whole` eager whole domain | — | 11 | 0.091195 s | 0.091229 s | 0.091663 s | **+0.0 %** | +0.5 % | 115 459 808 B both | 0.0057 → 0.0058 |
| `block.column` one destination column | — | 11 | 0.643202 s | 0.651846 s | 0.663159 s | +1.3 % | +1.7 % | 588 306 744 B both | 0.0385 → 0.0364 |

Per-pair cost is 827 ns → 893 ns (`micro.global`) and 758 ns → 767 ns
(`micro.polar`); the operator's cache keeps the sweep allocation-free (the 80 B
is `@timed`'s own overhead, not per call).

**Reading.** The clipping cost did not change. `block.whole` is +0.0 % and
allocates byte-identically; `block.column` is +1.3 % and allocates
byte-identically; the two large micro arms are +8.0 % and +1.2 %. Every value
fingerprint is identical across all three runs.

The one thing that is *not* noise: `micro.global` is +8.0 % / +6.7 % on two
independent "after" runs that agree with each other to 1.2 %, so it reproduces.
Since `sutherland_hodgman.jl` and `spherical_orient` are byte-identical source
and the summed area is bit-identical, this can only be a code-layout /
inlining / specialisation artefact of the larger module, not a change in the
work done. It is confined to that arm — the same kernel on the polar geometry is
−1.9 % to +1.2 %, and neither block builder sees it. Recorded, not explained.

## 6. What the branch *does* speed up

Measured with the same discipline, because the honest answer to "what is the pin
worth" is not visible in the clipping arms alone.

### 6.1 RelateNG region queries — `DGG.query` (`src/engine/query.jl`)

California multipolygon, `-t 8`, 7 samples per arm, min:

| arm | before | after | delta | bytes before → after |
|---|---:|---:|---:|---|
| `query.intersects` IGeo7 L7 | 0.130487 s | 0.120948 s | −7.3 % | 21.6 MB → 13.7 MB |
| `query.within` IGeo7 L7 | 1.969130 s | 0.737109 s | **−62.6 %** | 1 347 MB → 493 MB |
| `query.intersects` IGeo7 L9 | 5.573261 s | 5.497855 s | −1.4 % | 541 MB → 477 MB |
| `query.within` IGeo7 L9 | 99.515544 s | 42.288359 s | **−57.5 %** | 64.6 GB → 22.2 GB |

Result sets are identical (same cell counts). The repeat "after" run reproduces
to within 1.8 %. Attributable to `rk_orient`'s exact fallback moving from
`ExactPredicates.orient(a, b, c, origin)` to the new
`UnitSpherical.exact_spherical_orient`, to `_spherical_region_extent`'s new
vertex-cap prefilter, and to GO PR #472 — which the rebase folded in and which
this measurement cannot separate.

### 6.2 The generic `point_in_cell` / `cellat` fallback

`-t 1`, 11 samples, 4 000 (ring, point) pairs. IGeo7 has no native `cellat`, so
this is the path it takes.

| arm | before | after | delta | bytes/call |
|---|---:|---:|---:|---|
| `point_in_cell.centroid` L5 | 856 ns/call | 182 ns/call | **−78.7 %** | 640 B → **0 B** |
| `point_in_cell.foreign` L5 | 913 ns/call | 239 ns/call | **−73.9 %** | 640 B → **0 B** |
| `point_in_cell.centroid` L9 | 847 ns/call | 185 ns/call | **−78.1 %** | 640 B → **0 B** |
| `point_in_cell.foreign` L9 | 920 ns/call | 241 ns/call | **−73.8 %** | 640 B → **0 B** |

This is the `_ring_contains` / `_ring_encloses_parity` specialisation fix in
`UnitSpherical/predicates.jl`: Julia declined to specialise on the `Function`-
typed predicates the keyword bodies only forwarded, so every edge dispatched
dynamically and allocated. Binding them to `where` type parameters removes the
allocation entirely.

### 6.3 The new exact orient predicate, directly

3 000 IGeo7 L11 vertex triples:

| | ns/call | bytes/call |
|---|---:|---:|
| `ExactPredicates.orient(a, b, c, origin)` (before and after) | 57.6 / 58.8 | 109.3 |
| `UnitSpherical.exact_spherical_orient` (new) | **4.23** | **0** |

**13.6× faster, allocation-free**, same sign function. Matches the docstring's
own claim in direction.

## 7. A bug in the branch — spherical Foster-Hormann, reproducer included

While answering "would switching DGG's operator to the now-spherical FH help?",
the two clippers were run over the same 37 647 harvested pairs. They disagree:

- Σ area: SH `1.7954127730867251`, FH `1.8009722998004403` — **3.10e-3
  relative**, far outside floating-point noise.
- 24 795 pairs differ at all; **108** differ by more than 1e-9 relative.
- **95 pairs where SH returns exactly 0 and FH returns a positive area.** In
  every one, FH's answer equals *the clip polygon's whole area, to the last
  digit*.
- 0 pairs the other way round.

`OverlayNG` agrees with SH (0.0). The failing geometry is two rings that meet
only along the equator: a 1° raster cell spanning latitude [−1, 0] and the
IGeo7 L4 hexagon immediately north of it. FH's "no crossings at all: either one
polygon contains the other, or they are disjoint" branch appears to classify
edge-touching as containment. The branch's tip commit is
"Rebuild no-crossing results in the caller's representation", which is the
natural place to look.

Self-contained reproducer (no DiscreteGlobalGrids needed; coordinates lifted
verbatim from the failing pair):

```julia
import GeometryOps as GO, GeoInterface as GI
const USP = GO.UnitSpherical.UnitSphericalPoint
m = GO.Spherical(radius = 1.0)

subject = GI.Polygon([GI.LinearRing([          # 1° raster cell, lat [-1, 0]
    USP(0.5591077356830366, -0.828911306117208, -0.01745240643728351),
    USP(0.57348907788161,   -0.8190272834649945, -0.01745240643728351),
    USP(0.573576436351046,  -0.8191520442889918,  0.0),
    USP(0.5591929034707468, -0.8290375725550417,  0.0),
    USP(0.5591077356830366, -0.828911306117208, -0.01745240643728351)])])

clip = GI.Polygon([GI.LinearRing([             # IGeo7 L4 hexagon, lat [0, +1.5]
    USP(0.5650152296594431, -0.825080475016157,  1.0888033474565051e-16),
    USP(0.5720775736396,    -0.820092382092957,  0.013256491682489636),
    USP(0.5672946026105719, -0.8230973309611949, 0.026222463909228798),
    USP(0.5560304420970006, -0.8307544552348473, 0.026022731771791912),
    USP(0.5499227585688582, -0.8351125525975779, 0.013114270928092375),
    USP(0.5541748359488557, -0.8324002950510404, 1.2195901602416452e-16),
    USP(0.5650152296594431, -0.825080475016157,  1.0888033474565051e-16)])])

GO.intersection_area(GO.ConvexConvexSutherlandHodgman(m), subject, clip)  # 0        ✓
GO.intersection_area(GO.OverlayNG(; manifold = m),        subject, clip)  # 0.0      ✓
GO.intersection_area(GO.FosterHormannClipping(m),         subject, clip)  # 5.2374128241001073e-4  ✗
GO.area(m, clip)                                                          # 5.2374128241001073e-4
```

Symmetric in the argument order. **A synthetic pair of lon/lat rectangles
meeting on the equator does *not* reproduce it** — FH returns 0 there — so the
trigger involves the hexagon's geometry (a vertex at z ≈ 1.09e-16, i.e. exactly
on the shared great circle, is the obvious suspect), not merely edge contact.

Nothing in the GeometryOps tree was modified.

**This does not affect DiscreteGlobalGrids today**, which never calls FH; it is
reported because it is the direct answer to whether DGG should switch.

### 7.1 And on cost, FH is not competitive here anyway

Same 37 647 pairs, `-t 8`:

| clipper | ns/pair | bytes/pair |
|---|---:|---:|
| `ConvexConvexSutherlandHodgman`, cached — same run, boxed closure | 1 012 | 336 |
| `ConvexConvexSutherlandHodgman` through `IntersectionAreaOperator` (§5, this pin) | **893** | ~0 |
| `FosterHormannClipping`, no cache | 4 313 | 2 782 |
| `FosterHormannClipping`, `FosterHormannCache` (new on this branch) | 3 822 | 772 |

(The three-clipper arms all go through a boxed closure so they share one harness;
the closure costs SH ~120 ns and 336 B per pair that the production operator does
not pay. FH's numbers carry the same overhead, so the *ratio* is the fair read.)

The new cache is worth **−11 % wall and −72 % allocation on FH**, which is a real
improvement to FH. It still leaves FH **3.8× slower than Sutherland-Hodgman in
the same run** (4.3× against the production operator's own 893 ns/pair), which is
the expected ordering: DGG's cell pairs are convex-convex, exactly SH's sweet
spot.

## 8. Driver-scale implication — **[estimated]**, from a bounded measurement

Extrapolated from §5, using `2026-08-23-full-run-attribution.md`'s buckets
(clipping **[measured] 70 093 / 291 450 samples = 24.050 %**,
**[estimated] 34.826 core-h** of **[estimated] 148.10**, at
**[estimated] 16.81** integrated mean cores, **[measured] 8.81 h** wall).

The measured clipping delta is **0 %** (`block.whole` +0.0 %, `block.column`
+1.3 %, `micro.polar` −1.9 % to +1.2 %, all with byte-identical allocations).
So:

| | value |
|---|---:|
| Clipping core-hours saved | **[estimated] 0.00 of 34.826** |
| Whole-run wall saved | **[estimated] 0.00 h of 8.81 h** |
| Envelope if the true clipping delta were the worst measured block figure (+1.3 %) | **[estimated] −0.45 core-h**, i.e. **+1.6 min** of the 8.81 h (a *regression* of 0.3 %) |
| Envelope at ±3 % (the block arms' run-to-run band) | **[estimated] ±1.04 core-h = ±3.7 min = ±0.7 % of wall** |

**This is an extrapolation from a bounded measurement, not a measured run.** No
production run was performed. The §6 wins do not enter this table because
neither surface appears in the attribution report's bucket list, and the one
place the driver does call `DGG.query` — `covering_chunks` at planning — was
measured flat (§4, +1.4 %).

## 9. Recommendation

**GO-WITH-CAVEATS.**

Adopt the pin if the goal is (a) staying current with the branch the user is
developing, (b) the `DGG.query` and `cellat` wins, which are large and real, or
(c) removing 640 B/call of allocation from a fallback IGeo7 actually uses. It is
correctness-neutral by the strongest test available: six suites unchanged to the
assertion, 26 gate rows identical field for field, and nine value dumps
bit-identical.

**Do not adopt it as a clipping optimisation.** The premise — that clipping is
24.1 % of the run and this branch attacks clipping — does not connect: DGG's
clipping is Sutherland-Hodgman, which the branch does not touch. Expect **zero**
production-driver improvement.

Caveats attaching to a GO decision:

1. **The pin is not what the card described.** The requested SHA is gone; the
   live tip is a rebase that pulls in `main` (including PR #472 and `main`'s
   cap-primitive resolution). Anyone re-running this must re-fetch, and anyone
   comparing to it must use `35996798b`, not a branch name — the branch has been
   force-pushed once already **today**.
2. **The branch carries a live bug** (§7). It is not on DGG's path, but pinning
   a branch means pinning its defects, and this one silently reports a whole
   cell's area as an intersection.
3. `micro.global` is a reproducible +8 % on identical clipping source (§5).
   Small, confined to one arm, unexplained.

## 10. Residual uncertainty

- **The query win cannot be attributed.** The rebase folded GO PR #472 into the
  same diff as the exact-orient and extent-prefilter changes. §6.1's −57 % to
  −63 % is the *pin's* effect, which is what a pin decision needs, but it is not
  evidence about the FH branch's own commits.
- **No production run.** §8 is arithmetic on a bounded measurement, as the card
  required. The attribution report's own core-hour split is itself
  sample-attributed and labelled `[estimated]` there; §8 inherits that.
- **`micro.pentagon` is 137 pairs.** Fifty-one reps make the 1.3 % figure
  stable, but it is a small, geographically concentrated sample and should not
  carry weight on its own. The pentagon *values* are bit-identical, which is the
  claim that matters.
- **`query.within.L9` is one 100 s → 42 s arm at 7 samples.** Reproduced across
  two "after" runs to 0.7 %, but it is a single shape (one multipolygon, one
  predicate, one level).
- **`GeometryOpsCore = "main"` was left alone and did not move**, so the two
  arms differ in one dependency. If `main` moves before this is merged, that
  invariant lapses and the comparison must be re-run.
- **Shared box.** Sub-millisecond gate cases swung ±50 % on unchanged relations
  in both directions. Every conclusion above rests on arms of ≥ 30 ms with a
  repeat run, not on those.
- **The FH bug's trigger is not isolated.** §7 shows a failing case and a
  passing near-neighbour; which geometric feature separates them was not
  determined, and finding that out is GeometryOps' work, not this card's.

## 11. Artifacts

Everything in this record is reproducible from the scripts used; per the
repo policy only this `.md` is committed. The ndjson and dumps stayed local:
`clip-{before,after,after2}.ndjson`, `pred-{before,after}.ndjson`,
`gates-{before,after}.ndjson` (26 rows each), `values-{before,after}/` (nine
`.f64` dumps plus SHA-256 summaries), `fh.ndjson`, and the suite logs.
