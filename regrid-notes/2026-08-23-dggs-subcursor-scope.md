# Scoping: a `subcursor` for hierarchical DGGS grids

**Date:** 2026-08-23
**Branch:** `claude/dggs-subcursor-scope`, cut from `claude/chunk-dag-coverage` (`7a06ce0`).
**Status:** scoping only. No production implementation, no PR.
**Prototype:** `/home/asinghvi17/geo/scratch-subcursor/` (untracked, throwaway).

## Verdict

**Build it, but not now, and not for the reason the gap was framed.**

The gap is real and a windowed cursor is both correct and measurably better than
the packed fallback — 1.5× to 2.4× faster and 4.3× leaner across the whole size
range, with an exact leaf set over pentagons and both poles. But **no workload we
run reaches the fallback today**: production, the lazy smoke test, and every
example land on case 1 or case 3. The measured win is currently worth zero
seconds of wall clock.

What justifies keeping it on the list is not the win, it is a **latent cliff**:
the production destination clears case 4 by a single comparison inside
`_defaulttilesizes`, with 4.07× headroom in `budget`. Below `budget = 2^28` that
comparison flips and every destination work unit builds a packed R-tree over
~419 k cells — measured at **1.21 s and 68 MiB per unit, on 40 workers**.

Recommendation is **option (b′)**, in §5: fold the windowed cursor into **Task
B4**, which already owns this dispatch decision, and add a cheap guard against
the cliff independently. Do not create a new phase for it.

---

## 1. Reachability — which call paths actually reach case 4

`GR.subtree(space::DGGSpace, inds)` (`src/regridding.jl:239`) tries: whole space →
`subcursor` → exact chunk range → `GR.CellSpaceRTree`.

### `subtree` has exactly four call sites, all in the conservative method

| Site | `inds` | Reaches case 4 for a DGGSpace? |
|---|---|---|
| `lib/GlobalRegridding/src/conservative.jl:381-382` (`wholeblock`, eager) | literally `1:ndst` / `1:nsrc` | **No** — always case 1 |
| `lib/GlobalRegridding/src/conservative.jl:351` (source side of `build_weights!`) | `cellindices(src_space, chunk)` = `space.ranges[k]` verbatim | **No** — always case 2 or 3 |
| `lib/GlobalRegridding/src/conservative.jl:350` (destination side) | `dst_inds` from `_tileindices` | **Yes, possible** |
| `lib/GlobalRegridding/src/conservative.jl:176` (`TileCells` memo) | the same `dst_inds` | same as above |

`NearestCell` / `BilinearPoint` never call `subtree` (`interpolation.jl:12`, `:126`
use `cellat`/`chartcoords`). `chunkgraph.jl` and `discovery.jl` contain no
`subtree` calls. The halo/radius machinery widens *which chunks are selected*,
never the index set handed to `subtree` (`lazy.jl:573-587`).

### The one live route is the destination tiling

`_tileindices` (`lib/GlobalRegridding/src/lazy.jl:571-572`):

```julia
_tileindices(A::LazyRegridArray, t::Int) =
    A.tiling.spacetiled ? cellindices(A.plan.dst_space, t) : A.tiling.runs[t]
```

`spacetiled = true` gives exactly `space.ranges[t]` → case 3. `spacetiled = false`
gives a budget-derived run from `_runs` with no relation to chunk boundaries →
case 4. `_spacetileable` (`lazy.jl:62-81`) makes `spacetiled` false whenever the
user passed `chunks =`, **or** the destination `DGGSpace` has exactly one chunk.

### Two things I expected to be gaps and are not

**`PartialGrid`-backed DGGSpaces are fine.** I checked this directly: a
`PartialGrid` over IGeo7 yields a `WindowCursor` (`selection === nothing`, because
`has_sorted_subtrees(IGeo7System())` is true), `_chunkwindows` gives it real chunk
ranges via `_subsetwindow`, and `_chunkcursor` succeeds. Both a contiguous 2 M-cell
slab and a scattered every-third-cell subset returned a `CapCachedTree`, not a
`CellSpaceRTree`. `_iswholespace` also fires correctly for `1:n` over a
`PartialGrid`, because its positions are `1:length(ids)` over the subset
(`src/engine/partial_grid.jl:142-143`) — asserted at
`test/systems/crosssystem/regrid.jl:341-346`.

**A5 never reaches case 4 either.** `has_sorted_subtrees(A5System())` is false
(`src/systems/A5/system.jl:115`), so `DGGSpace` takes `_wholechunk`
(`src/regridding.jl:41`): one chunk, so `inds` is always `1:n`, so `subtree` always
takes case 1 and `_cachedcelltree` hands back `GR.celltree(space)`.

### The production driver does not reach case 4

`scripts/copdem_production.jl`:

- **Source** (`:1218-1219`): `DGGSpace(PartialGrid(CopernicusDEMSystem, 1, ids); chunklevel = 0)`; `sinds` is always one whole 1° tile, which the CopDEM `subcursor` recognizes as one rectangular id run → **case 2**.
- **Destination** (`:817-824`): `DGGSpace(DGG.subtree(sys7, a, 12); chunklevel = 5)`. The `PartialGrid` is rooted at level 5 and `chunklevel = 5`, so `_ancestorpositions` returns one position → `nchunks == 1`, `ranges == [1:823543]`. `_spacetileable` fails on `length(spans) > 1`, so `spacetiled = false`, and `_defaulttilesizes(823543, 1, 2^30)` returns `min(fld(2^30, 320), 823543) = 823543` → a **single** run → `_iswholespace` → **case 1**.

`scripts/copdem_lazy_smoke.jl:95-101` is the same shape at `budget = 2^29`. Examples
use eager `DirectPlan`s → `wholeblock` → case 1 on both sides.

### The only existing constructions of `CellSpaceRTree` over a DGGSpace are tests

`test/systems/CopernicusDEM/runtests.jl:1037`, `:1039`, `:1071`, and the deliberate
A/B in `benchmark/fallback_cell_tree.jl:45`. **No test drives a DGGSpace into case 4
through the executor** — that is itself a coverage gap.

### The cliff

Production clears case 4 by one comparison in `_defaulttilesizes`
(`lazy.jl:84-88`, constants at `:44-46`: `DEST_BUDGET_SHARE = 8`,
`DEST_BYTES_PER_CELL = 40`):

```
frombudget = budget / 320  must be >=  ndst = 7^(level - ancestor) = 823_543
=> budget >= 263_533_760 bytes (~251 MiB)
```

| `budget` | headroom | outcome |
|---|---|---|
| `2^30` (production) | 4.07× | case 1 |
| `2^29` (lazy smoke) | 2.04× | case 1 |
| `2^28` | 1.02× | case 1, by 1.9% |
| `2^27` | **fails** | 2 runs of ~419 k → **case 4, both** |

The same cliff trips on shape rather than budget if `level - ancestor` reaches 8
(`7^8 = 5_764_801 > 3_355_443` at the default budget), or if anyone passes
`chunks =` to `GR.regrid` in `regrid_chunk`. Note `copdem_production.jl:814-816`
already warns against `chunks =` there — for the unrelated reason that it defeats
source pruning. The guardrail exists but not for this hazard.

**This is the finding that sets the priority.** The gap costs nothing today and
is one configuration change away from costing 1.2 s and 68 MiB per work unit.

---

## 2. Cost of the fallback — measured

IGeo7 level 6 (1,176,492 cells, chunklevel 2, 492 chunks), destination IGeo7
level 4, `julia --threads=4`, Julia 1.12.6. Harness:
`/home/asinghvi17/geo/scratch-subcursor/check2.jl`. "Join" is the real consumer —
`GR._intersectionareas` against a destination tree, i.e. `build_weights!`'s body,
so it includes discovery *and* clipping. Best of 3 after warmup.

Three trees compared on identical index sets:
**A** = bare windowed cursor (prototype), **B** = the same window wrapped in the
`CapCachedTree` seam that cases 1 and 3 use, **F** = `GR.CellSpaceRTree`.

| window | \|inds\| | build A | build B | build F | join A | join B | join F | **total A** | **total B** | **total F** |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| chunk[2] exact | 2 401 | 0.00000 | 0.00243 | 0.00501 | 0.00688 | 0.00345 | 0.00410 | 0.0069 | **0.0059** | 0.0091 |
| chunk[2] shifted half | 2 401 | 0.00000 | 0.00239 | 0.00479 | 0.01168 | 0.00482 | 0.00389 | 0.0117 | **0.0072** | 0.0087 |
| 10 k unaligned | 10 000 | 0.00000 | 0.01029 | 0.02078 | 0.04492 | 0.01530 | 0.01977 | 0.0449 | **0.0256** | 0.0405 |
| 100 k unaligned | 100 000 | 0.00000 | 0.04601 | 0.20250 | 0.31280 | 0.09286 | 0.08738 | 0.3128 | **0.1389** | 0.2899 |
| 419 k (the cliff) | 419 000 | 0.00000 | 0.11018 | 0.85430 | 1.43553 | 0.39058 | 0.35872 | 1.4355 | **0.5008** | 1.2130 |

Seconds. Build allocations:

| \|inds\| | B (MiB) | F (MiB) | ratio |
|---:|---:|---:|---:|
| 2 401 | 0.09 | 0.39 | 4.3× |
| 10 000 | 0.38 | 1.63 | 4.3× |
| 100 000 | 3.82 | 16.33 | 4.3× |
| 419 000 | 15.99 | 68.41 | 4.3× |

**The nonzero count was identical across A, B and F on every window** — the three
trees produce the same weights, not merely the same leaf set.

### What the numbers say

1. **The fallback's cost is entirely its build, and the build is entirely the
   per-cell cap.** `CellSpaceRTree` calls `_packedcellcap` → `GI.getpoint(getcell(space, i))`
   for every cell (`conservative.jl:39`). Contiguity is irrelevant:

   | \|inds\| | contiguous | stride 7 | random |
   |---:|---:|---:|---:|
   | 2 401 | 0.00479 s | 0.00487 s | 0.00500 s |
   | 10 000 | 0.02085 s | 0.02029 s | 0.02051 s |
   | 100 000 | 0.20012 s | 0.20249 s | 0.20253 s |

   All at identical allocation (0.39 / 1.63 / 16.33 MiB). **A scattered subset costs
   the fallback exactly what a contiguous one does.** This kills option (b) in §5
   before it starts: a "union of chunk ranges" cursor buys nothing that plain
   contiguity does not already buy.

2. **The bare cursor (A) is dominated everywhere.** Its build is free but its join
   is 3–4× worse, because a bare cursor re-derives every node extent from cell
   geometry on each visit. This is exactly the warning at
   `src/cap_cached_tree.jl:100-107`. **Any DGGS `subcursor` must return a
   `CapCachedTree`, not a bare cursor.**

3. **The cached window (B) wins on total time at every size** — 1.54× at chunk
   scale, 2.42× at the cliff — and on memory by a flat 4.3×. The join times for B
   and F are within 9% of each other, so the entire win is the O(1) tree structure
   replacing an O(n) Hilbert pack, with the O(n) cap fill common to both.

4. **DGGS chunk size is a constant 2401 cells at every level.** `_chunklevel`
   targets `DEFAULT_CHUNK_CELLS` and lands on `level - 4` for IGeo7 at L5–L10, so
   `length(ranges[k])` is `7^4` throughout. Two consequences: the `|inds| = 100 k`
   and `419 k` rows describe only the destination-tiling cliff, never a source
   chunk; and `_CHUNK_CAP_CACHE_MAX = 2^16` (`cap_cached_tree.jl:132`) is **dead
   code for IGeo7 at any level** — no chunk ever reaches it.

---

## 3. Design

### The windowing primitive already exists

The natural primitive is not a new cursor type. It is
`Engine.SubtreeIds(grid, first, n)` (`src/engine/partial_grid.jl:10-35`), which
despite its name is a lazy `O(1)` `AbstractVector` view of **any** contiguous
position run of a complete level grid, with `Helpers.strictly_increasing` already
declared true so `PartialGrid`'s O(n) ascent check is elided.

```julia
subcursor(grid, inds) ≈
    CapCachedTree(treeify(PartialGrid(system(grid), level(grid),
                                      SubtreeIds(grid, first(inds), length(inds));
                                      root = enclosing_ancestor(inds))),
                  leafcaps)
```

### Why this is the right shape

`HierarchicalGridCursor` (`src/engine/cursor.jl:28-39`) has two descent arms:

- **Non-`PartialGrid`** (`cursor.jl:132-133`): `_child_window` returns the child's
  full `descendant_range`, **unclamped by the parent's window**.
- **`PartialGrid`** (`cursor.jl:139-149`): clamps to `[first_index, last_index]`
  with two binary searches over the id vector.

This is precisely why `_chunkcursor` (`src/regridding.jl:249-260`) only works on an
*exact* chunk range: a chunk range **is** a subtree's `descendant_range`, so the
children's ranges partition it exactly and the missing clamp never bites. Widen
`inds` by one cell and the unclamped arm emits leaves outside `inds`, violating
`subcursor`'s "must cover exactly the cells at `inds`, no more"
(`src/interface/grid.jl:558-566`).

Routing through a `PartialGrid` gets the clamp for free rather than adding a third
`_child_window` method. It also switches on `_subtree_count`
(`cursor.jl:270-274`), which for a `PartialGrid` compares the stored count against
the true subtree size — enabling the `STORED_UNION_CAP_LIMIT = 64` tightening in
`STI.node_extent` (`cursor.jl:249-262`) that is **dead on a complete level grid**,
where `_subtree_count == _stored_count` always. So a windowed node gets a tighter
cap than the unwindowed hierarchy would give it.

### What it can and cannot serve

| subset shape | cost | how |
|---|---|---|
| whole space | O(1) | already case 1 |
| exact chunk range | O(1) | already case 3 |
| **any contiguous position range** | **O(1) structure + O(\|inds\|) cap fill** | **this design** |
| non-contiguous / `Vector{Int}` inds | — | **must still fall back**; no window exists |
| `!has_sorted_subtrees` system (A5) | — | **must still fall back**; no `descendant_range`, and the selection cursor indexes leaves by selection slot, not grid position (`cap_cached_tree.jl:121`) |

The O(\|inds\|) cap fill is not avoidable and is not the thing being fixed — the
fallback pays it too. What is removed is the Hilbert pack, the boxes, and the
position vector.

### Three things that need care

**(i) The position shift is bidirectional, and this is the one real trap.**
`PartialGrid` positions are `1:length(inds)`, but `subcursor`'s contract is the
parent grid's own positions. Leaves must come *out* shifted by
`first(inds) - 1`, and `getcell` must go back *in* unshifted. In the prototype
this is a `ShiftedCursor` wrapper implementing both directions; getting only one
direction right yields a tree whose leaf set is correct and whose **weights are
silently wrong**. `_ShiftedCaps` (`cap_cached_tree.jl:16-34`) is the existing
one-directional analogue. The wrapper must also carry
`Trees.cell_index_count` (global domain size, for matrix assembly) and
`Trees.split_weight` (restricted population) — the contract asserted at
`lib/GlobalRegridding/test/test_conservative.jl:387-390`.

**(ii) Pentagons are safe, and I verified it rather than reasoning about it.**
The hazard at `src/interface/system.jl:156-172` is that aperture-7 children extend
beyond the parent boundary, so `node_extent` does not nest and child caps can stand
outside a parent's extent. This design never compares a parent's extent against a
child's — it inherits `node_extent(sys, id)` unchanged from the existing cursor and
only *narrows the index window*. Pentagon subtrees are smaller than hexagon ones
(2001 vs 2401 cells at L6; 12 of 492 chunks), but `descendant_range` reports the
true size, and windowing is arithmetic on positions that does not care.

Measured, cached variant, `check3.jl`: **9 of 9 windows exact, 0 mismatches**,
covering the north-pole pentagon subtree, that subtree shifted by one, its
interior, a window straddling two pentagon subtrees, both single-cell poles, the
south-pole tail, a 95 k window spanning 40 chunks, and the whole space. The bare
variant was checked separately on 10 windows against `CellSpaceRTree`'s leaf set:
also 0 mismatches.

**(iii) Rooting.** Descending from the synthetic root fans over all 12 root cells;
the `PartialGrid` clamp prunes empty ones immediately, so it is correct either
way, but rooting at the deepest enclosing ancestor bounds the descent. The
prototype walks `first(levels(sys)):level(grid)` comparing `ancestor(sys, lo, a)`
against `ancestor(sys, hi, a)` — O(depth), and it is what `_tree_root`
(`cursor.jl:64-67`) already consumes.

---

## 4. What this would subsume

A general `subcursor(::HierarchicalLevelGrid, inds)` **replaces `_chunkcursor`
entirely** — an exact chunk range is a contiguous range, and `_cachedchunktree`'s
job becomes the general cached-window path. `src/regridding.jl:249-260` and the
`_CHUNK_CAP_CACHE_MAX` branch (already dead, §2.4) both go. That is a net deletion,
which matters for a plan whose gates are stated in deletions.

The one caveat is measured: at chunk scale the cached window is 1.54× faster than
the fallback but the *existing* case-3 path is faster still (chunk[2] exact: case 3
joins in 0.00323 s vs the window's 0.00345 s, from `check.jl`). Replacing
`_chunkcursor` is a simplification, not a speedup, and must not regress the
already-fast aligned case.

---

## 5. Options

Effort is in the plan's own unit: one card = one focused commit.

| | option | effort | correctness risk | buys |
|---|---|---|---|---|
| **a** | General `subcursor` on the hierarchical cursor, via the lazy-`PartialGrid` window above | **2 cards** — one for the window + shift wrapper, one to retire `_chunkcursor`/`_CHUNK_CAP_CACHE_MAX` | **Medium.** The bidirectional shift is a silent-wrong-weights hazard (§3.i). Everything else is arithmetic on positions and is covered by leaf-set identity tests. Prototype: 19/19 windows exact, nz identical | 1.5–2.4× and 4.3× memory on any contiguous window; removes the cliff; net deletion |
| **b** | Narrower "union of chunk ranges" cursor | 1 card | Low | **Nothing.** §2.1 measures the fallback as contiguity-independent, and the union case is strictly smaller than (a) for the same machinery. Rejected on evidence |
| **c** | Improve `CellSpaceRTree` instead | 1–2 cards | Low | Little. Its cost is the per-cell `getcell` cap (§2.1), which is intrinsic to packing cells it does not know are a hierarchy window. Best case is removing the boxes and position vector — part of the 4.3×, none of the 2.4× |
| **d** | Do nothing | 0 | — | Correct **today**. Leaves the cliff at 4.07× headroom, undocumented and untested |
| **b′** | **Recommended.** Fold (a) into Task B4 as an explicit action, and guard the cliff now | (a)'s 2 cards, absorbed into B4's existing scope; guard is a few lines in an existing card | Medium, same as (a) | (a)'s win where it applies, plus the cliff closed immediately rather than when B4 lands |

### Why (b′)

Task B4 (`plan:337-359`) already owns this exact decision. Its actions read:

> - Add one thin cell-space adapter over `FlexibleRTrees.RTree` […]
> - **Use it only when a native restricted cursor is unavailable.**

B4 assumes native restricted cursors exist where they should. For hierarchical
DGGS they do not, so B4's "only when unavailable" is today vacuous on the DGGS
side. Making the DGGS cursor available is B4 finishing its own sentence, not a new
phase — and B4's verification list already names "native versus fallback
Conservative block identity", which is precisely the proof obligation in §6.

The cliff guard is separable and should not wait for B4: `_defaulttilesizes`
silently degrades from an O(1) hierarchy tree to a 419 k-cell packed R-tree with no
diagnostic. Either clamp destination runs to chunk boundaries when the destination
is chunked, or emit a warning when a destination run is neither whole-space nor
chunk-aligned. Either is a few lines and is worth having regardless of which
option wins.

**Do not create a new phase or card ID for this.** It is B4 scope plus one guard.

---

## 6. Proof obligations

The central assertion, and the one that would have caught the shift trap:

> For every window `inds`, the tree returned by `subcursor` and the tree returned
> by `GR.CellSpaceRTree(space, inds)` yield **the same cell set with the same
> global positions**, and **the same weight block**.

Concretely, three assertions per window:

1. `sort(leafpositions(subcursor(grid, inds))) == collect(inds)` — exact, no more and no fewer.
2. `Trees.cell_index_count(t) == ncells(space)` and `Trees.split_weight(t) == length(inds)`, plus `sum(Trees.split_weight, STI.getchild(t)) == Trees.split_weight(t)`.
3. `build_weights!` over the windowed tree equals `build_weights!` over `CellSpaceRTree` **bit for bit** — nonzero pattern, values, and denominators. This is the one that catches a one-directional shift; assertions 1 and 2 both pass with `getcell` mis-shifted.

Window set — every one of these is exercised in the prototype and passes:

- an exact chunk range (regression: the case-3 path must not change);
- that range shifted by one cell, and its strict interior;
- a window straddling two chunk boundaries;
- **a pentagon subtree, shifted and interior** — at L6, chunks 1, 42, 83, 124, 165, 206, 247, 288, 329, 370, 411, 452 are the 12 non-uniform ones;
- **`1:1` and `n:n`** — the two poles as single cells;
- a window spanning tens of chunks;
- `1:n` (must still take case 1);
- a **scattered** `Vector{Int}` — must still return `CellSpaceRTree`, asserted by type;
- **an A5 space** — must still fall back, asserted by type, since `has_sorted_subtrees` is false.

### Where these belong

| assertion | file |
|---|---|
| leaf-set identity, pentagons, poles, and the fall-through-by-type cases | `test/systems/crosssystem/regrid.jl` — already the home of the DGGS `subtree`-dispatch assertions (`:66-69`, `:341-346`, `:364-371`) |
| the `cell_index_count` / `split_weight` triple | `lib/GlobalRegridding/test/test_conservative.jl`, beside the existing contract test at `:387-390` |
| native-vs-fallback weight-block identity | `lib/GlobalRegridding/test/test_conservative.jl`, mirroring the CopDEM oracle at `test/systems/CopernicusDEM/runtests.jl:1071` |
| the shape of the CopDEM precedent to copy | `test/systems/CopernicusDEM/runtests.jl:995-1039` asserts leaf positions equal `collect(inds)` and that non-rectangular windows fall through to `GR.CellSpaceRTree` — the DGGS version is the same two assertions |

**Add regardless of which option wins:** a test that drives a DGGSpace into case 4
*through the executor* (small `budget`, multi-chunk destination). §1 found no such
test, so the cliff is currently unguarded by anything.

---

## 7. Fit with the plan

**Belongs to Task B4** (`regrid-notes/2026-08-21-regridding-simplification-plan.md:337-359`),
as an added action and two added verification items. Not a new card.

**Blocked by:** B4's prerequisite chain (B3, hence B1/B2 and A1–A4). Nothing in
this design needs the spherical-extent work, but B4 owns the file and the plan
runs one write agent at a time in the shared worktree.

**Blocks:** nothing.

**Interaction with E2** (`plan:504-525`, which removes the `chunktree` bridge —
`lib/GlobalRegridding/src/rastergrid.jl:967-969`): **none.** E2 is raster-side; it
retires `RasterFlatTree`'s last role and the `chunktree` export. The DGGS
`subcursor` touches `src/regridding.jl` and `src/cap_cached_tree.jl` only. The two
can land in either order.

**Standing acceptance laws** (`plan:80-97`): law 1 (no geometric false negatives)
is the §6.1 obligation; law 3 (eager/lazy/differently-chunked equivalence) is §6.3.
Both are already B4 verification items. No law needs amending.

**Phase-5 (W1/W2) allocation work is behind the user's gate and is untouched here.**

---

## 8. Residual uncertainty

1. **One system measured.** All numbers are IGeo7 at L6. H3 and the
   `AbstractQuadFaceGridSystem` family (HEALPix, ISEA4R, S2) also set
   `has_sorted_subtrees = true` and would take the same path, but their
   `descendant_range` implementations differ (H3 has an explicit one-element case
   at `src/systems/H3/system.jl:228-235` specifically because window descent asks
   for it) and I did not benchmark or correctness-check them. H3's icosahedral
   distortion and HEALPix's polar caps are plausible sources of different
   cap-tightening behaviour, not different leaf sets.

2. **The join comparison is single-destination.** All join times use an IGeo7 L4
   destination. A raster destination exercises a different tree pairing, and the
   two trees speak different extent vocabularies at internal nodes — the DGGS
   cursor's are `SphericalCap`s at every level, `CellSpaceRTree`'s internal nodes
   hand back XYZ `Extents.Extent` boxes and only its *leaves* return caps
   (`conservative.jl:74`, `:110-112`). That asymmetry is Task B2's territory
   ("make ConservativeRegridding extent-generic") and could shift the B-vs-F join
   numbers either way once B2 lands.

3. **Threading.** Measured at `--threads=4` on a shared box; `_leafcaps` self-threads
   above ~16 k cells (`cap_cached_tree.jl:158-160`) while the R-tree pack does not,
   so the B/F build ratio at 100 k and 419 k is partly a threading artifact and
   would narrow at `t=1`. The 2401-cell row is single-threaded on both sides and is
   the honest small-window number.

4. **The cliff arithmetic is derived, not observed.** I did not run production at
   `budget = 2^27` to watch it fall into case 4 — that is a 26 h job and the
   Copernicus source is user-gated. The 419 k figure comes from
   `_defaulttilesizes`' arithmetic on the production shape; the 1.21 s / 68 MiB
   cost of a 419 k-cell `CellSpaceRTree` is measured directly.

5. **`_cachedcelltree` has no size guard**, unlike `_cachedchunktree`'s
   `_CHUNK_CAP_CACHE_MAX`. The production path therefore materializes
   `_leafcaps(grid, 1:823_543)` — roughly 26 MB of caps — per work unit, times 40
   workers. This is on the hot path today, is unrelated to the subcursor gap, and
   looks intentional (the `_CACHED_BUCKET_SIZE = 49` seam depends on the cache
   existing), but the asymmetry with the chunk path is unremarked. Flagging, not
   proposing.

---

## Appendix — prototype

`/home/asinghvi17/geo/scratch-subcursor/` (untracked, throwaway, not on any branch):

| file | what |
|---|---|
| `proto2.jl` | the window, the `ShiftedCursor` bidirectional wrapper, timing helper |
| `check.jl` | leaf-set identity vs `CellSpaceRTree`, pentagon subtree sizes, bare-cursor build/join |
| `check2.jl` | the §2 table — bare vs cached vs fallback, build and join, five window sizes |
| `check3.jl` | cached-variant pentagon/pole identity; scattered-vs-contiguous fallback cost |
| `bench.jl` | the first-pass case-3-vs-case-4 comparison that established chunk-size invariance |

Run as `LEVEL=6 julia --project=. --threads=4 <file>` from the repo root.
The prototype is ~60 lines and is deliberately not production-shaped: it takes
`HierarchicalLevelGrid` only, does no `nothing`-returning validation beyond bounds,
and wraps rather than extends the cursor.
