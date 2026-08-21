# What is left after c1–c4: the next ranked improvements for the CopDEM→IGeo7 run

Companion data: `2026-08-20-more-improvements.ndjson`.
Predecessors: `2026-08-20-production-profile.md` (the profile whose §6 list this
goes beyond), `2026-08-20-perf-ladder-spec.md` (c1–c4, being implemented
elsewhere), `2026-08-19-post-stack-profile.md`.

**Scope.** c1 (source-tree cap memo), c2 (leaf-cap allocation reuse), c3
(wave/inner-threading shape) and c4 (leaf-cap cache) are out of scope and are
not re-argued here. Everything below is either orthogonal to them or acts on
what survives them.


## Headline

Three changes, none of them on the c1–c4 list, none of them touching an
algorithm, together take a production column from **28.8 s to 13.2 s (2.18×)**
at `-t 1` with **bit-identical output** and **44 % less allocation** — measured
back-to-back in one process, on the real plan, on both a fed column and a
100 %-NaN one:

| | change | what it is | alone |
|---|---|---|---:|
| **N1** | `bucket_size` on the destination cursor | the IGeo7 dual-tree descent currently runs all seven levels down to **823 543 leaves of one cell each**, because `PartialGrid`'s `bucket_size` defaults to 0. Stopping two levels early is a **default**, not a rewrite. | **−57.5 %** |
| **N2** | `spherical_orient` | pays two `isUnitLength` asserts and a `normalize` (`sqrt` + 3 divides) per call, then reads only `sign(n ⋅ c)`. Called **114 212 848 times per column**. | **−5.97 %** |
| **N3** | `cell_boundary` for CopernicusDEM | returns a heap `Vector`, which `closed_ring` copies into a second heap `Vector`. IGeo7 already returns an isbits `SmallList`; the source system never got the same treatment. | **−2.57 %** |

The single most useful thing in this note for the c1–c4 implementer is **not**
any of those. It is §2.1's finding that **leaf sizes tuned for a tree that
rebuilds its caps are badly wrong for a tree that caches them** — the isolated
sweep says a coarse destination leaf is +25 %, the production tree says −57 %,
and the difference is `CapCachedTree`. **c4 turns the *source* tree into a
cached one, so `BlockCursor`'s `LEAF_CELLS = 9` must be re-swept the moment c4
lands.**

---

## 0. Method, and why the numbers here are trustworthy

Same worktree code as the live run (`claude/copdem-production` @ `47961ea`,
private worktree `/home/asinghvi17/geo/DGG-hunter`, branch `claude/perf-hunter`,
root `Manifest.toml` copied from the production worktree so the CR
git-tree-sha1 pin is preserved). CR `6a4b997`, GO `2825c47`, Julia 1.12.6.
Live run PID 2915199 untouched; every experiment in its own process at `-t 1`
against its own scratch, box load 44–68 throughout and recorded per row.

The reference workload is one level-5 IGeo7 column of 823 543 level-12 cells
regridded conservatively from the `DGGSpace`-over-`PartialGrid` of 26 475
CopernicusDEM tiles — the production work unit. Reference columns: **728**
(mid-latitude, fed, 5 source chunks), **34818** (equatorial, fed),
**98241** (single source chunk, 100 % NaN output), plus fully-empty
columns 98/106/117/118/167.

**Harness validation.** `-t 1` on column 728 measured **28.87 s / 28 528 cells
per core-second**, against the production profile's **28.3 k** and the live
run's **28.6 k**. Allocation **3.3015 GiB/column** against the profile's
3.22 GiB. GC 4.02 % against 4.5 %. The harness reproduces the campaign it
builds on.

**Two protocol traps this campaign fell into, and how they were fixed** — both
matter for anyone re-running this:

1. *Cross-process timing on this box is worthless at this resolution.* The same
   column, same code, measured 29.07 s, 31.73 s and 28.87 s in three different
   processes at loads 67, 55 and 51. All patch attribution below is therefore
   **same-process, back-to-back**, patch applied by `@eval` between timed
   stages. The in-process reproducibility, measured as a null stage (a patch
   that turned out not to dispatch — see the method note below), is
   **0.23 %**.
2. *A `LazyRegridArray` shares its plan's weight-block cache.* Collecting the
   same array twice measures the cache, not the work — the second read came
   back **700× faster** and bit-identical. Every variant must be built from a
   **fresh `GR.regrid(...)` call**. (An earlier version of the empty-column
   experiment was invalidated by this and re-run.)

**A method note that matters more than it sounds.** Redefining a GO/DGG method
from `Main` only takes effect if the redefinition's signature is *at least as
specific* as the original. `GO.UnitSpherical.spherical_orient` is declared
`(a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint)`; an
untyped `(a, b, c)` redefinition is silently **less** specific and never
dispatches, with no warning. This produced a false negative that was only
caught by counting calls. Anyone prototyping this way must assert the patch is
live (call count, or `methods(...)`) before believing a null result.

---

## 1. Where the column's 28.87 s actually go

Measured on column 728, `-t 1`, by direct instrumentation rather than by
sampling:

| quantity | measured | how |
|---|---:|---|
| candidate pairs emitted by the dual DFS | **4 217 328** | counted (`e3_dfs`, two fed chunks) |
| nonzeros surviving the clip | 2 920 000 (report §3.2) | — |
| **candidate → nonzero hit rate** | **69 %** | ⇒ 31 % of clips produce nothing |
| cap-overlap predicate calls (`_intersects`) | **51 495 549** | counted (`e9_counts`) |
| **predicate calls per candidate pair** | **12.2** | — |
| `spherical_orient` calls | **114 212 848** | counted (`e9_counts`) |
| **`spherical_orient` calls per candidate pair** | **27.1** | — |
| leaf/leaf visits | **2 315 991** | counted (`e4_prune`) |
| mean source leaf size | **6.79 cells** | counted |
| **mean destination leaf size** | **1.00 cell** | counted |
| dual-DFS wall (two fed chunks, instrumented) | ~20.0 s of 28.9 | `e3_dfs` |
| ⇒ clip + assembly | ~8.9 s of 28.9 (**31 %**) | by difference |

Two facts here were not in the previous campaign and drive most of what
follows.

**The destination tree descends to single cells.** `STI.isleaf` for
`HierarchicalGridCursor` (`src/engine/cursor.jl:194–199`) stops only at
`leaf_level` unless `bucket_size > 0`, and `PartialGrid`'s `bucket_size`
defaults to **0** (`src/engine/partial_grid.jl:69`). So the IGeo7 side is a
7-way tree walked seven levels down to 823 543 one-cell leaves. Every one of
those levels costs interior node visits on **both** trees.

**The clip is called 27 times per candidate pair.** `spherical_orient` is one
call per subject vertex per clip edge (`GO .../sutherland_hodgman.jl:310,316`),
i.e. ~6 hexagon edges × ~4.5 live vertices.

---

## 2. Ranked improvements

Percentages are **of the measured single-thread column** (28.87 s, column 728,
`-t 1`). Where an item acts on work that c1–c4 do not touch, its absolute
saving survives them and its *share* grows.

| # | improvement | owner | measured on today's column | risk |
|---|---|---|---:|---|
| **N1** | **destination-cursor `bucket_size`: stop the IGeo7 descent early** | DGG | **−57.5 % / −61.9 %** at `bucket_size=49` (2.4–2.6×); −50 % at `8` | low |
| **N2** | strip `spherical_orient` of its unit-length asserts and its irrelevant `normalize` | GO | **−5.97 %** | low |
| **N3** | CopernicusDEM `cell_boundary` → inline `SmallList` (isbits polygon) | DGG | **−2.57 %**, allocation −9.0 % | low |
| **N4** | transcendental-free cap-overlap predicate | GO | **−0.32 %** | low |
| **N6** | `_sh_spherical_intersection`'s three redundant normalizations | GO | not isolated, ≈ 2–4 % argued | low |
| **N5** | value-aware source pruning at **sub-tile** granularity | DGG + GR | argued; per-*tile* variant measured **zero** | high |
| **N7** | latent gap: `knownempty`/`dropempty` never exercised on this path | GR | correctness | — |

Every one of N1–N4 was verified **bit-identical** (`ndiff = 0`, exact
`isequal` including NaN positions) against the unpatched baseline on every
reference column it was run on.

### 2.0 Projected share of a post-c1..c4 column

The brief asks for each item as a share of the column *after* c1–c4. That
requires a judgement about overlap, so here it is explicitly rather than buried
in a percentage.

c1 (source cap memo) and c4 (source leaf-cap cache) both make the **source
tree walk** cheaper. **N1 removes most of that walk outright** — it is not
"another 57 % on top", it is largely the *same* work c1/c4 would have made
cheap. So N1 and c1/c4 are closer to `max` than to `sum`, and the ordering
matters: whichever lands first collects most of the credit.

N2, N3 and N6 act on the **clip**, which c1–c4 and N1 all leave alone. Their
absolute savings are invariant, so their share rises as everything else falls.
This is measured, not assumed: N2 is −5.97 % of a 28.9 s column standalone and
−6.5 % of a 14.6 s column after N1 — the same ~0.95 s both times (§2.1a).

| item | absolute saving on column 728 | share of today (28.9 s) | share of a post-N1 column (12.5 s) | overlaps c1/c4? |
|---|---:|---:|---:|---|
| **N1** | 16.9 s | 57.5 % | — | **yes, heavily** |
| **N2** | ~0.95 s | 5.97 % | **≈ 7.6 %** | no |
| **N3** | ~0.40 s | 2.57 % | **≈ 3.2 %** | partly (allocation) |
| **N4** | ~0.09 s | 0.32 % | ≈ 0.7 % | no |
| **N6** | not isolated | ≈ 2–4 % argued | ≈ 3–5 % | no |
| c2-equivalent (measured here) | ~0.96 s | 3.44 % | — | **is c2** |

The practical consequence: **measure c1, c2 and c4 against a `bucket_size = 49`
baseline.** Against today's baseline they will look better than they will
actually deliver.

### 2.1 N1 — the destination tree descends one level too far

This is the headline. `STI.isleaf` for `HierarchicalGridCursor`
(`src/engine/cursor.jl:194–199`) stops only at `leaf_level` unless
`bucket_size > 0`, and `PartialGrid`'s `bucket_size` defaults to **0**
(`src/engine/partial_grid.jl:69`). So the destination side of the dual DFS is a
7-way tree walked all seven levels from L5 to L12, ending in **823 543 leaves
of one cell each** (measured mean destination leaf size: **1.00**). The last
level alone holds 6/7 of all the tree's nodes, and every node visited there
costs an interior predicate on both trees and a `_box_cap` on the source side.

**End-to-end, on the real plan and the real `CapCachedTree`**, fresh plan per
read, 2 reps, min wall (`e12_bucket`):

| `bucket_size` | col 728 wall | GiB | col 34818 wall | GiB | output |
|---:|---:|---:|---:|---:|---|
| **0 (current)** | 29.373 s | 3.3015 | 28.012 s | 2.8474 | — |
| 4 | 28.886 s | 3.3015 | 28.200 s | 2.8474 | identical |
| **8** | **14.627 s** | **2.1582** | **13.536 s** | **1.7817** | **identical** |
| 16 | 14.579 s | 2.1582 | 13.211 s | 1.7817 | identical |
| **49** | **12.473 s** | 2.3388 | **10.677 s** | 1.8676 | **identical** |
| 343 | 18.528 s | 4.8298 | — | — | identical |

* `bucket_size = 4` is a **null control**: an IGeo7 node has 7 children, so
  `7 > 4` and no node ever becomes a leaf early. Allocation is identical to 13
  significant figures and the ±0.7 % wall difference is the noise floor.
* `bucket_size = 8` makes level-11 nodes (7 cells) leaves. **−50.2 % and
  −51.7 % wall, −34.6 % and −37.4 % allocation, `ndiff = 0`.**
* `bucket_size = 16` reproduces `8` to 13 significant figures of allocation —
  the tree is the same, as it must be for any value in `7:48`. The two wall
  readings differ by 0.3 % and 2.4 %, which brackets the noise on this box.

`bucket_size = 49` — **two** IGeo7 levels early — is better still: **−57.5 %
and −61.9 %**, again bit-identical, and `343` (three levels) turns back up to
18.53 s. So the end-to-end curve is 29.4 → 14.6 → **12.5** → 18.5 and the
optimum is **49**, i.e. `7²`, two refinement levels above the destination
resolution.

**And here is the part worth reading twice.** The *isolated* dual-DFS sweep
(`e6_bucket`, same source chunk, counted, but against a **bare
`HierarchicalGridCursor`** rather than the production `CapCachedTree`) says the
exact opposite about coarse leaves:

| `bucket_size` | dst leaf | predicate calls | preds/pair | candidate pairs | DFS wall |
|---:|---|---:|---:|---:|---:|
| **0 (current)** | 1 cell | 38 091 199 | 12.21 | 3 120 315 | **19.77 s** |
| **8** | ≤ 8 cells | **30 128 838** (−20.9 %) | 9.63 | 3 130 187 (+0.3 %) | **10.73 s (−45.7 %)** |
| 50 | ≤ 50 cells | 82 149 366 (+116 %) | 26.24 | 3 130 187 | 24.66 s (+25 %) |
| 350 | ≤ 350 cells | 364 021 699 (+856 %) | 116.29 | 3 130 187 | 106.86 s (+441 %) |

Isolated, `bucket_size = 50` is **+25 %** and `350` is **+441 %**; end-to-end,
`49` is **−57.5 %**. The two disagree because the isolated tree rebuilds a
leaf's cell caps on every visit, so a 49-cell leaf costs 49 `cell_cap` calls
per visit, whereas production's `CapCachedTree`
(`src/cap_cached_tree.jl:47–55`) has them precomputed and a big leaf is nearly
free. **The current leaf sizes are tuned for a tree that does not cache its
leaf caps, and the production tree does.**

That generalises, and it is the real lesson of N1: *once caps are cached,
leaves should be bigger.* The same argument applies to the **source** side the
moment c4 lands — `BlockCursor`'s `LEAF_CELLS = 9`
(`src/systems/CopernicusDEM/cursor.jl:4`) was chosen under the same
uncached assumption, and c4 invalidates it. **Whoever lands c4 should re-sweep
`LEAF_CELLS` immediately afterwards**; on this evidence the win from doing so
could be comparable to c4's own.

The candidate-pair set grows by 0.3 % going from `0` to `8`, and those extra
pairs clip to zero area, so the answer does not move at any setting tested
(`ndiff = 0` throughout).

**Implementation.** `bucket_size` already exists end to end and is already
plumbed: `PartialGrid(...; bucket_size)`
(`src/engine/partial_grid.jl:69`) → `_grid_bucket_size`
(`src/engine/cursor.jl:57`) → `HierarchicalGridCursor` → `_chunkcursor`
(`src/regridding.jl:188–198`, which forwards `root.bucket_size`) →
`CapCachedTree`. Nothing new has to be built; the default is simply wrong for a
dual-tree join against a source of comparable resolution. The principled fix is
to **derive it** rather than hard-code 8: the right stopping level is the one
where a destination node's cell count is comparable to the opposing tree's leaf
size (`LEAF_CELLS = 9` for `BlockCursor`), which is what a `bucket_size` of
7–9 expresses. A blunt `bucket_size = 8` default for `subtree`-derived
destination grids captures the whole measured win.

**Why this is orthogonal to c1–c4.** c1, c2 and c4 make each source-node visit
and each leaf-cap block *cheaper*; N1 removes a whole layer of visits. The
allocation drop (−34.6 %) is itself mostly `Tuple{Int,Cap}` vectors — the same
bytes c2 targets — so the two overlap partially in allocation but not in the
tree-walk arithmetic. N1 composed with N2 and N3 is measured in §2.1a.

**Latitude and regime coverage — the caveat is resolved.** The sweep above is
mid/low-latitude with 4–5 source chunks, so `e14_polar` repeats `0` vs `49` on
two −76.7° polar columns (8–10 source chunks, source leaves spanning far more
longitude) and on the 100 %-NaN column:

| column | latitude | NaN | `bucket_size = 0` | `bucket_size = 49` | speedup | output |
|---:|---:|---:|---:|---:|---:|---|
| 115426 | −76.7° | 4.2 % | 36.883 s | **13.616 s** | **2.71×** | identical |
| 115427 | −76.7° | 56.4 % | 37.390 s | **13.870 s** | **2.71×** | identical |
| 98241 | −33.8° | 100 % | 27.455 s | **11.481 s** | **2.39×** | identical |
| 728 | 47.3° | 0 % | 29.373 s | 12.473 s | 2.35× | identical |
| 34818 | 9.6° | 0 % | 28.012 s | 10.677 s | 2.62× | identical |

**N1 holds at every latitude and every fed/empty regime tested, and it is
*better* at the poles** (2.71× against 2.35–2.62×), where the base column is
most expensive. Every row is bit-identical, NaN positions included.

### 2.1a Do N1 and N2/N3 compose? Yes.

Same process, fresh plan per read, 2 reps, min wall (`e13_combo`). Column 98241
is the 100 %-NaN single-source-chunk column, included to check that N1 is not a
fed-column-only effect:

| stage | col 728 | Δ vs previous | col 98241 (100 % NaN) | Δ vs previous |
|---|---:|---:|---:|---:|
| base (`bucket_size = 0`) | 28.805 s | — | 27.155 s | — |
| `bucket_size = 8` | 14.576 s | **−49.4 %** | 13.379 s | **−50.7 %** |
| + N2 `spherical_orient` | 13.628 s | **−6.5 %** | 12.295 s | **−8.1 %** |
| + N3 `SmallList` boundary | 13.226 s | **−2.9 %** | 11.587 s | **−5.8 %** |

Cumulative: column 728 **28.805 s → 13.226 s (2.18×)**, allocation 3.3015 →
1.8594 GiB (**−43.7 %**); column 98241 **27.155 s → 11.587 s (2.34×)**,
allocation 2.8722 → 1.5840 GiB (**−44.9 %**). Bit-identical output at every
stage on both columns — from three changes, none of which touches an algorithm,
and none of which is on the c1–c4 list. And `bucket_size = 49` beats the `8`
used here by a further ~15 %, so the realistic combined figure is **≈ 2.5×**.

`ndiff = 0` at every stage on both columns. Two things follow:

* **N1 is not fed-column-specific.** The 100 %-NaN column halves too — it pays
  the same tree walk, it just multiplies by NODATA at the end.
* **N2's absolute saving survives N1, so its share grows.** Standalone it is
  −5.97 % of a 28.9 s column; on top of N1 it is −6.5 % of a 14.6 s column,
  i.e. the same ~0.95 s. That is the general shape of every clip-side item on
  this list, and the reason N2/N3/N6 should not be discounted just because c1–c4
  and N1 are bigger: they act on the part that is left.

### 2.2 N2 — `spherical_orient` pays for work it throws away

`GO/src/utils/UnitSpherical/predicates.jl:30–42` calls
`robust_cross_product(a, b)` and then reads only `sign(n ⋅ c)`. But
`robust_cross_product`
(`.../robustcrossproduct/RobustCrossProduct.jl:98–190`)

* opens with **two `@boundscheck @assert isUnitLength(...)`** which, because the
  function is not `@inline`d into an `@inbounds` context, execute on every call
  — each an `isapprox` on `sum(abs2, v)`;
* returns **`normalize(result)`** — a `sqrt` and three divides — whose
  magnitude is then discarded by the sign test;
* recomputes `kMinNorm`, containing two `sqrt(3)`, per invocation.

Stripped version (keep the stable `cross(a−b, a+b)` formula and the
degeneracy branch, drop the asserts, drop the normalize, and scale the
tolerance by `norm_sqr` so the `< 16 eps` degeneracy test is unchanged —
compare `dp² < tol²·nn` and no `sqrt` is needed):

* microbenchmark **34.70 ns → 2.85 ns, 12.2×**, 0 disagreements in 512 cases
  (`e3_kernel`);
* **114 212 848 calls per column** (`e9_counts`);
* end-to-end, same process, 3 reps each: **29.416 s → 27.659 s, −5.97 %**,
  output **bit-identical**, allocation unchanged (`e10_summary`).

The clip is ~31 % of today's column and is untouched by c1–c4, so this is
**≈ 10 % of the post-c1..c4 column** and it is the largest single item that
survives the ladder.

Files: `GeometryOps/src/utils/UnitSpherical/predicates.jl`,
`GeometryOps/src/utils/UnitSpherical/robustcrossproduct/`. **GO-side**, so it
needs a GO branch and a pin bump (precedent: the existing GO PR-473 pin).

### 2.3 N3 — the CopernicusDEM cell boundary is the only heap polygon left

`DGG.cell_boundary(::CopernicusDEMSystem, ::LevelIndex)`
(`src/systems/CopernicusDEM/system.jl:220–228`) returns a **heap
`Vector{UnitSphericalPoint}`**, which `closed_ring`
(`src/fallbacks/geometry.jl:13–26`) then copies into a **second** heap vector.
IGeo7 does the opposite — it returns a `Helpers.SmallList`, `closed_ring` has a
`SmallList` method (`geometry.jl:31–40`), and the destination polygon is
therefore isbits and never touches the heap. The source side simply never got
the same treatment.

The source cell polygon is rebuilt roughly **once per candidate pair** (the
`CellMemo` in `conservative.jl:249–276` is 64 direct-mapped slots, and the
source index is the *outer* loop of the leaf/leaf pass, so it misses once per
source cell per leaf-pair visit). That matches the profile's
`Memory{UnitSphericalPoint}` = 11.5 % of allocation.

Returning a `SmallList{4,USPoint}` (3 for the pole triangles) measured
**28.711 s → 27.973 s, −2.57 %**, allocation **3.3015 → 3.0031 GiB, −9.0 %**,
output bit-identical (`e7_ladder` stage `orient+cap+bnd`).

Files: `src/systems/CopernicusDEM/system.jl:220`. **DGG-side, self-contained.**
Risk: `cell_boundary`'s return type changes for this system; anything asserting
`Vector` in the CopernicusDEM tests needs updating. The interface already
documents both shapes.

### 2.4 N4 — the cap-overlap predicate

`GO.UnitSpherical._intersects` (`cap.jl:99–101`) is
`spherical_distance(x.point, y.point) <= x.radius + y.radius`, and
`spherical_distance` (`point.jl:144`) is `atan(norm(cross(x,y)), x⋅y)` — a
cross, a norm, a dot and an `atan` per call, **51.5 M times per column**.

The test can be answered from the squared chord with no transcendental and no
new fields on `SphericalCap`, because the chord is a lower bound on the arc and
`arc ≤ chord·(1 + chord²/20)` for `chord ≤ 1` (the bound is tight at
`chord = 1`: `2·asin(0.5) = 1.0472 ≤ 1.05`):

```
s = x.radius + y.radius;  s >= π && return true
c2 = |x.point − y.point|²;  ss = s*s
c2 > ss && return false                       # chord ≤ arc ⇒ arc > s
c2 <= 1 && (f = 1 + c2/20; c2*f*f <= ss) && return true
return spherical_distance(x.point, y.point) <= s   # razor-thin band, ~never
```

Microbenchmark **13.09 ns → 2.17 ns, 6.0×**, 0 disagreements (`e3_kernel`).
End-to-end, same process: **28.803 s → 28.711 s, −0.32 %**, bit-identical,
allocation unchanged (`e7_ladder`).

The gap between 6× on the kernel and 0.3 % end-to-end is the honest finding:
**the predicate is not where the DFS's time is** — cap *construction*
(`_box_cap`: five `UnitSphereFromGeographic` and four `spherical_distance` per
cap) is, and that is exactly what c1/c4 remove. Do N4 anyway (it is a handful of
lines and it is the last transcendental in the search), but expect ≈ 0.6 % of
the post-ladder column, not 6×.

### 2.5 N5 — source regions with no valid data (the per-tile version is a dead end)

**The prize is real.** **937 of the first 7 470 completed columns (12.5 %) come
back 100 % NaN, and they account for 8.7 % of the run's wall clock** (done-log
mean 17.91 s against 27.38 s for a fed column). They are not cheap: column
98241 is 100 % NaN and costs the **full 27.7 core-seconds** at `-t 1`, building
2.55 M weights against pixels that are every one of them NODATA.

**The obvious mechanism does not reach it.** `GR.knownempty(data, ndchunk)`
(`lib/GlobalRegridding/src/executor.jl:37`) already exists as the "this storage
chunk holds nothing valid, skip it" hook, and `_connectedsource!` already
filters on it via `_allempty` (`lazy.jl:520–528`). `TiledDEM` does not
implement it, and `dropempty` is off anyway because `LazyRegridArray` sets
`dropempty = !usesreference(policy)` (`lazy.jl:265`) while production runs
`Weighted(0.5)`.

I implemented `knownempty(::TiledDEM, ranges)` exactly — a listed tile is
all-NODATA iff its 1° box contains no land cell of the run's own 15-arcsec
mask, which is exact here because the GLO-90 post lattice (1/1200°) is finer
than the mask (1/240°) — and forced `dropempty = true`. Result
(`e11_dropempty`, fresh plan per variant):

| column | NaN | source chunks | chunks `knownempty` called on | **chunks reported empty** | base | with drop | speedup |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 98241 | 100 % | 1 | 2 | **0** | 27.71 s | 27.63 s | **1.00** |
| 98 | 100 % | 3 | 6 | **0** | 3.64 s | 3.30 s | 1.11 |
| 117 | 100 % | 1 | 2 | **0** | 4.44 s | 4.45 s | 1.00 |

**Zero.** The census says **977 of 26 475 listed tiles (3.7 %) are 100 % ocean**
(`e6_census`, 12.2 s for the whole census), but the fully-NaN *columns do not
meet those tiles*. They meet ordinary **land** tiles and happen to fall on the
ocean part of them. Chunk granularity is the tile, and the tile is not empty.

*(An earlier run of this experiment reported 700–6900× speedups. That was the
plan's weight-block cache — see §0 trap 2 — not the drop. The table above is
the corrected measurement.)*

**What would actually work**, and why it is expensive:

* **Sub-tile validity.** Attach a coarse validity bitmap to the source — say
  64×64 blocks over a 1200×1200 tile, 512 bytes per tile, 13 MB for the whole
  holding — and let `BlockCursor`'s descent (`src/systems/CopernicusDEM/cursor.jl`)
  reject a node whose blocks are all invalid. This prunes the fully-NaN columns
  outright and, inside fed columns, the ~19 % of in-tile pixels that are ocean
  (33 % global land over a 40.9 % tile holding). Ceiling: the 8.7 % plus a part
  of the 19 %.
* **Making the drop exact under `Weighted`.** Whatever the granularity,
  `dropempty` is disabled because `addreference!` needs every connected chunk's
  row sums for the coverage denominator. That denominator is recoverable
  without the block: the row sum of a chunk's weight block is *exactly* the area
  of the destination cell intersected with the chunk's **footprint**, because
  the chunk's cells tile that footprint disjointly. So a skipped region's
  reference contribution is one clip per destination cell against one box
  instead of thousands of per-pixel clips.

Owner: DGG (the source cursor and the `TiledDEM` adapter) plus GR (the
reference path in `lazy.jl`). **Risk: high** — it makes the weight builder
value-aware, which it deliberately is not today, and it must be proved not to
move any coastal cell across the `Weighted(0.5)` threshold. This is the one
item on the list that is a design change rather than a fix; rank it last
despite the size of the prize.

### 2.6 N6 — redundant normalization in the clip (argued, not isolated)

`_sh_spherical_intersection` (`GO .../sutherland_hodgman.jl:253–294`) performs
**two `norm`s and two vector divisions (:266–267, :274–275), a cross, a third
`norm` (:279) and a third division (:287)** per computed intersection vertex —
two of which repeat the normalization `robust_cross_product` already did
internally. `_sh_clip_spherical!`'s prologue additionally runs
`point ≈ clip_points[1]` (an `isapprox`, hence a `sqrt`) per clip vertex per
pair (:368) and unconditionally `append!`s the whole subject ring into
`cache.subject` (:388) for a containment check consulted only in the empty-result
branch (:404–407).

With ~4 intersection vertices per pair and 4.22 M pairs, the redundant
normalizations alone are of the same order as N2's win. Not measured in
isolation here because it needs a real edit rather than a method swap; flagged
as the obvious follow-up **once N2 lands**, since both live in the same file
and the same profile line.

### 2.7 N7 — a latent bug worth fixing while nearby

`_allempty` (`lazy.jl:520–528`) loops over `A.othergroups`. For a source with
**no** non-spatial dimensions — which is exactly this workload, a 1-D
`TiledDEM` — `_groupgrid` returns a single empty tuple, so the loop runs once
and the semantics are fine. But the surrounding contract is never exercised on
the production path at all, because `dropempty` is `false` under `Weighted`.
Any test of `knownempty` on this stack is currently vacuous. Whoever
implements N5 should add a test that the hook is actually consulted.

---

## 3. Measured dead ends — do not re-investigate

| what | result | evidence |
|---|---|---|
| **Pre-pruning the leaf/leaf double loop** (test each source cell against the opposing *node* extent before looping over its cells) | **Worse.** Variant A: **+8.2 %** predicate calls. Variant B (also pre-filtering the destination cells against the source node): **+12.7 %**. Both produce the identical pair set. | `e4_prune`. The reason is §1: the destination leaf holds **exactly one cell**, so the "pre-prune" is one extra predicate per source cell to save at most one. Revisit *only* if N1 lands, which makes destination leaves plural. |
| **Better dual-DFS pruning quality in general** (the deferred F3/D1 row-band bisection) | **Not the problem on this path.** The CopDEM source already uses `BlockCursor`, an index *rectangle* tree that bisects the longer spatial axis, so the linear-range pathology the F3/D1 comment describes (`conservative.jl:43–48`) never arises here. Measured pruning quality: 12.2 predicates per emitted pair, 69 % of emitted candidates have nonzero area. The cap test is already within ~1.44× of the theoretical minimum for circumscribing circles (2π/4 = 1.57). | `e3_dfs`, `e4_prune` |
| **`_fillcoo!` + WeightBlock sparse fusion** (the −5.5 s idea from 2026-08-19) | **Not applicable.** That campaign's serial tail was 9–10 s of a 50.8 s eager MOC regrid. Here the destination is one tile written by one task and the whole assembly/`_fillcoo!`/`_writechunk!` tail is **0.07 %** of the column. | production profile §4 |
| **Source tile-cache locality / batch and column ordering** | **Nothing to win.** The entire source path — `_sourcefor!` → `readblock!` → `tilevalues!` → `synthetic_tile` — is **0.02 %** of CPU, with 14 loads and 0 evictions per 3 columns. | production profile §3.1 |
| **"Compute the cell polygons once per tile instead of per visit"** (the brief's own candidate) | **Already true — nothing to win.** Counted `cell_polygon` calls for one whole column (`e6_polycount`): destination IGeo7 **834 942 builds for 823 543 cells = 1.014 per cell**; source CopernicusDEM **953 935 builds = 1.16 per destination cell, i.e. 0.23 per candidate pair**. The 64-slot direct-mapped `CellMemo` (`conservative.jl:249–276`) is doing its job almost perfectly, because the dual DFS emits pairs in tree order. Precomputing a whole tile's 823 543 IGeo7 polygons would cost ~150 MB per worker to save ~1.4 % of the builds. What *is* worth fixing is not the count but the **cost of each build** — which is N3. | `e6_polycount` |
| **Float32 in the clipper** | **Do not.** Cell separations at IGeo7 L12 are ~3×10⁻⁵ rad on the unit sphere; Float32's 6×10⁻⁸ relative resolution puts the rounding error at ~0.2 % of a cell edge, which is a visible area error on a conservative weight. The narrowing already happens at the right place — `finalize!` (`executor.jl:296–320`) stores into a Float32 output from Float64 accumulators. | argued |
| **`spherical_orient` with an untyped redefinition** | Measures **0.2 %** — because it never dispatches. This is the reproducibility floor, not a result. | §0 |

---

## 4. What this says about the ladder itself

One number the implementer should have: a **c2-equivalent** (task-local
reusable buffer for `BlockCursor.child_indices_extents`, replacing the
`push!`-grown fresh `Vector{Tuple{Int,Cap}}` per leaf visit) measured, in the
same process, on top of N2+N3+N4:

* **27.973 s → 27.010 s, −3.44 %**
* allocation **3.0031 → 1.4846 GiB, −50.6 %**
* GC **4.08 % → 2.49 %**
* output bit-identical

So c2 does everything the spec expects to the *allocation* — it more than halves
it — but converts that into only **−3.4 % wall, not the −9 % the spec
projects**. Allocation is not this workload's bottleneck; the `push!`-growth
(a fresh vector grown 1→2→4→8 per call, ~3 reallocations for 6.79 entries) is
what makes the byte count look worse than the time. Worth recalibrating the
ladder's arithmetic before c2's result is compared against its target.

Full same-process ladder, column 728, 3 reps each, min wall (`e7_ladder`):

| stage | wall | Δ | GiB | GC % |
|---|---:|---:|---:|---:|
| base | 28.870 | — | 3.3015 | 4.02 |
| + (null stage) | 28.803 | −0.23 % | 3.3015 | 4.00 |
| + cap predicate (N4) | 28.711 | −0.32 % | 3.3015 | 5.09 |
| + `SmallList` boundary (N3) | 27.973 | −2.57 % | 3.0031 | 4.08 |
| + c2-equivalent leaf-cap buffer | 27.010 | −3.44 % | 1.4846 | 2.49 |

and separately, N2 on its own base (`e10_summary`): 29.416 → 27.659, **−5.97 %**.

---

## 5. Recommended order

1. **N1 — `bucket_size` on the destination cursor.** A default, not a rewrite,
   and it **more than halves the column** (29.4 s → 12.5 s at `49`,
   bit-identical on both fed reference columns). It removes tree-walk *visits*,
   which is a different axis from c1/c2/c4's cheaper-per-visit, so it should be
   done **first** — every later measurement is then taken against the tree the
   run will actually use.
   *Do before shipping:* confirm on a polar column (`e14_polar` — a −76.7°
   column costs 36.9 s at `bucket_size = 0`, so the absolute prize there is
   larger still), and decide whether to hard-code `49` or derive the stopping
   level from the opposing tree's leaf size.
2. **Re-sweep `LEAF_CELLS` on the source side the moment c4 lands.** N1's
   result is that leaf sizes tuned for an uncached tree are badly wrong for a
   cached one, and c4 turns the source tree into a cached one. This costs one
   sweep and may be worth as much as c4 itself.
3. **N2 — `spherical_orient`.** −5.97 %, bit-identical, and it lives in the
   clip, which c1–c4 do not touch. Needs a GO branch and a pin bump.
4. **N3 — `SmallList` cell boundary for CopernicusDEM.** −2.57 % and −9 %
   allocation, entirely DGG-local, and it brings the source system in line with
   what IGeo7 already does.
5. **N6 — the redundant normalizations in `_sh_spherical_intersection`.**
   Free-riding on N2's file and profile line.
6. **N4 — the cap predicate.** A handful of lines for ~0.3 %; do it while N2 is
   open, not on its own.
7. **N5 — sub-tile value-aware pruning.** Largest remaining prize (the 8.7 % of
   wall spent on 100 %-NaN columns) but the only genuine design change on the
   list. Last.

**Recalibrating the ladder's own arithmetic.** The projections in
`2026-08-20-perf-ladder-spec.md` are anchored on a 28.9 s column. If N1 goes in
first, the column is ~12.5 s and c1/c2/c4's *absolute* savings shrink with it —
they are savings on the tree walk, and N1 removed most of the tree walk. The
honest expectation is that **N1 and c1/c4 substantially overlap**, and that
their combination is much closer to `max` than to `sum`. Measure c1 against a
`bucket_size = 49` baseline, not against today's.

## 6. Reproducing

Scratch driver (not committed; regenerate from this note). Same shape as the
production profile's §7, plus a `world()` helper that builds the production
source/mask/tile list without running `main()`:

```bash
S=/home/asinghvi17/geo/dggstores/hunter-scratch
head -n -1 <worktree>/scripts/copdem_production.jl > $S/prod.jl
cp <live store>.columns.txt $S/hunt.zarr.columns.txt   # skip the covering rebuild
# setup.jl sets ARGS to the production conf with store=$S/hunt.zarr, includes
# prod.jl, and exports world() / lazycolumn() / withcpu() / emit().
cd /home/asinghvi17/geo/DGG-hunter
RASTERDATASOURCES_PATH=/home/asinghvi17/geo/DiscreteGlobalGrids.jl/bench/data \
  julia --project=bench -t 1 $S/e7_ladder.jl
```

| script | what it produced | `exp` | status |
|---|---|---|---|
| `e1_empty.jl` | harness validation against the production profile; fed vs empty columns | `e1_column` | valid |
| `e3_count.jl` | dual-DFS predicate/pair census; kernel microbenchmarks | `e3_dfs`, `e3_kernel` | valid |
| `e4_prune.jl` | leaf/leaf pre-prune variants A/B; baselines; leaf sizes | `e4_prune`, `e4_baseline` | valid (A/B are dead ends) |
| `e5_patch.jl` | cross-process patch runs | `e5_patch` | superseded — cross-process drift |
| `e6_misc.jl` | ocean-tile census; `bucket_size` sweep; cell-polygon build counts | `e6_census`, `e6_bucket`, `e6_polycount` | valid |
| `e7_ladder.jl` | same-process cumulative patch ladder | `e7_ladder`, `e7_summary` | valid |
| `e8_clip.jl` | isolated clip kernel | `e8_*` | **invalid** — degenerate synthetic geometry |
| `e9_counts.jl` | `spherical_orient` / `_intersects` call counts on a real column | `e9_counts` | valid |
| `e10_orient.jl` | same-process orient ladder with the **correct** signature | `e10_summary` | valid |
| `e11_empty.jl` | `knownempty` + `dropempty`, fresh plan per variant | `e11_dropempty` | valid |
| `e12_bucket.jl` | end-to-end `bucket_size` sweep on the real plan | `e12_bucket` | valid |

The worktree `/home/asinghvi17/geo/DGG-hunter` (branch `claude/perf-hunter`,
based on `claude/copdem-production` @ `47961ea`) is **clean** — every patch was
applied at runtime by `@eval`, so there is no source diff to push and no PR to
open. The exact patch bodies are in Appendix A below, ready to be lifted into
GeometryOps and DiscreteGlobalGrids.

Nothing wrote to the live store, log, done log or column cache; the live run's
own scratch, `ladder-scratch` and `profile-scratch` were not touched.


---

## Appendix A — the exact patches, ready to lift

All three were verified bit-identical on every reference column they were run
on. They are transcribed here because the scratch scripts are transient.

**N1 — no code at all.** Give the destination `PartialGrid` a `bucket_size`.
In the production script's `regrid_column`, the destination grid comes from
`DGG.subtree(sys7, a, LEVEL)`; the measurement rebuilt it as

```julia
g0 = DGG.subtree(w.sys7, a, LEVEL)
g  = DGG.PartialGrid(w.sys7, LEVEL, g0.ids; bucket_size = 49, root = a)
dstspace = DGG.DGGSpace(g; chunklevel = ANCESTOR)
```

The principled home for this is a non-zero default in DGG rather than a
caller-side rebuild — `_grid_bucket_size` (`src/engine/cursor.jl:57`) and
`PartialGrid`'s keyword (`src/engine/partial_grid.jl:69`) — ideally derived
from the opposing tree's leaf size rather than hard-coded.

**N2 — `GeometryOps/src/utils/UnitSpherical/predicates.jl:30`.** Keep the
stable `cross(a−b, a+b)` formula and the degeneracy branch; drop the two
`isUnitLength` asserts and the `normalize`, and scale the tolerance by
`norm_sqr` so the `|dot| < 16 eps` degeneracy test is preserved exactly without
a `sqrt`:

```julia
function spherical_orient(a::UnitSphericalPoint, b::UnitSphericalPoint, c::UnitSphericalPoint)
    d1 = a[1] - b[1]; d2 = a[2] - b[2]; d3 = a[3] - b[3]
    s1 = a[1] + b[1]; s2 = a[2] + b[2]; s3 = a[3] + b[3]
    n1 = d2 * s3 - d3 * s2; n2 = d3 * s1 - d1 * s3; n3 = d1 * s2 - d2 * s1
    dp = n1 * c[1] + n2 * c[2] + n3 * c[3]
    nn = n1 * n1 + n2 * n2 + n3 * n3
    tol = eps(Float64) * 16
    dp * dp < tol * tol * nn && return 0
    return dp > 0 ? 1 : -1
end
```

(The real `robust_cross_product` also has an exact-arithmetic fallback for
`norm_sqr(result) < kMinNorm^2`, i.e. `a` and `b` nearly equal or antipodal.
The clip guards `edge_start == edge_end` before calling
(`sutherland_hodgman.jl:397`), so that branch is unreachable on this path — but
a principled upstream patch should keep it rather than delete it.)

**N3 — `DiscreteGlobalGrids/src/systems/CopernicusDEM/system.jl:220`.**

```julia
function DGG.cell_boundary(sys::CopernicusDEMSystem, c::DGG.LevelIndex)
    west, east, south, north = cell_box(sys, c)
    H = DGG.Helpers
    e = H.empty_small_list(Val(4), NORTH_POLE)
    if north == 90.0
        e = H.small_push(e, TO_SPHERE((west, south)))
        e = H.small_push(e, TO_SPHERE((east, south)))
        return H.small_push(e, NORTH_POLE)
    elseif south == -90.0
        e = H.small_push(e, SOUTH_POLE)
        e = H.small_push(e, TO_SPHERE((east, north)))
        return H.small_push(e, TO_SPHERE((west, north)))
    end
    e = H.small_push(e, TO_SPHERE((west, south)))
    e = H.small_push(e, TO_SPHERE((east, south)))
    e = H.small_push(e, TO_SPHERE((east, north)))
    return H.small_push(e, TO_SPHERE((west, north)))
end
```

`closed_ring` already has a `SmallList` method (`src/fallbacks/geometry.jl:31`)
that keeps the result inline, so the polygon becomes isbits end to end, exactly
as IGeo7's already is.

**N4 — `GeometryOps/src/utils/UnitSpherical/cap.jl:99`.**

```julia
function _intersects(x::SphericalCap, y::SphericalCap)
    s = x.radius + y.radius
    s >= pi && return true
    d1 = x.point[1] - y.point[1]
    d2 = x.point[2] - y.point[2]
    d3 = x.point[3] - y.point[3]
    c2 = d1 * d1 + d2 * d2 + d3 * d3      # squared chord
    ss = s * s
    c2 > ss && return false                # chord <= arc, so arc > s
    if c2 <= 1.0
        f = 1.0 + c2 / 20                  # arc <= chord*(1 + chord^2/20), chord <= 1
        c2 * f * f <= ss && return true
    end
    return spherical_distance(x.point, y.point) <= s   # razor-thin band
end
```

The bound is `2·asin(c/2) = c·(1 + c²/24 + 3c⁴/640 + …) ≤ c·(1 + c²/20)` for
`c ≤ 1`, tight at `c = 1` (`1.0472 ≤ 1.05`). NaN-radius (whole-sphere) caps
behave as before: every comparison against NaN is false and the function falls
through to the original test.
