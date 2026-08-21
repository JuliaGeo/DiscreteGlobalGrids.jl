# The perf ladder: seven changes, 2.76x, bit-identical throughout

Companion data: `2026-08-20-perf-ladder.ndjson` (202 measurements).
Spec: `2026-08-20-perf-ladder-spec.md`. Profile it was built from:
`2026-08-20-production-profile.md`. Mid-campaign findings folded in:
`2026-08-20-more-improvements.md`.

Branch: **`claude/perf-ladder`** (pushed, no PR), 7 commits on `47961ea`.
Upstream: **`claude/perf-ladder-predicates`** on GeometryOps.jl (pushed, no PR).

---

## 1. Headline

Over the eight reference columns, single-threaded, bit-identical at every step:

**184.60 -> 67.0 core-seconds, 2.75x** (three independent final runs: 66.66,
67.13, 67.35). Allocation **3.30 -> 1.87 GiB** per column, **-43%**. Full suite
**987 584 pass / 0 fail / 17 broken**.

Projected production ETA at a 24-core budget, against the profile's 20:12Z
snapshot of 5.086e10 cells remaining:

| configuration | cells per core-second | ETA |
|---|---:|---:|
| as running, 22.24 cores | 28 614 | 22.2 h |
| 24 cores, unchanged code | 28 614 | 20.6 h |
| **24 cores, this branch** | **~78 700** | **~7.5 h** |

Scaled by the measured cells-per-core-second ratio, in the `outer` shape that
b2 defaults to. c3's shape is deliberately **not** counted: at a fixed 24-core
budget it costs ~8 % (§2.4).

Three things are worth more than the number.

1. **The spec's ranking was wrong, and measurement is the only reason we know.**
   The two changes the spec led with (c1 at -26%, c2 at -9%) were sized against
   a profile in which `node_extent` was 47.7 % of CPU. That profile was taken
   with the *destination* tree descending to 823 543 single-cell leaves. Once N1
   stops that descent two levels early, most of the calls c1 memoizes never
   happen: **c1 fell from -27.2 % to -4.3 %, and c2 fell to zero.** Run in the
   specced order and stopped there, this campaign would have reported ~35 % and
   left a 2.19x sitting underneath it.
2. **Wall-clock predictions from isolated microbenchmarks mostly did not
   survive; allocation predictions did.** c2 (-9 % predicted, null), N3 (-2.57 %
   predicted, null), c1 (-26 % predicted, -4.3 %) all shrank or vanished, while
   every allocation claim held or beat its estimate. The reason is the same in
   each case: an isolated kernel measurement cannot see that another change
   removed the calls.
3. **Two of the seven changes are correctness improvements, not just speed.**
   The new `SphericalCap` intersection test is *more accurate* than the one it
   replaces (§2.7), and the c2 investigation found that the spec's own
   prescribed design would have silently corrupted self-joins (§2.3).

---

## 2. The ladder, in commit order

Every row measured at `-t 1` on the eight reference columns (three mid-latitude
fed, one single-source 100 %-NaN, four polar), against frozen bit-exact
fingerprints -- NaN count plus a hash of every value's bits. Load 43-70
throughout, recorded per measurement in the ndjson.

| # | commit | change | core-s | Δ | alloc/col |
|---|---|---|---:|---:|---:|
| 0 | `47961ea` | base | 184.60 | — | 3.30 GiB |
| 1 | `9e9e2c6` | **N1** destination bucket size | **84.13** | **-54.4 %** | 2.34 GiB |
| 2 | `c337155` | **c1** source node-extent memo | 80.47 | -4.3 % | 2.34 GiB |
| 3 | `72391eb` | **c2** inline leaf cells | ~80.5 | null | 2.12 GiB |
| 4 | `ef43b0b` | **c3** wave scheduling | no-op at `-t 1` | — | — |
| 5 | `224903b` | **b2** worker sizing | script only | — | — |
| 6 | `f868d7b` | **N3** inline cell boundary | null | null | 1.87 GiB |
| 7 | `7627ee3` | **N2+N4** GeometryOps repin | **66.9** | **-15 %** | 1.87 GiB |

Column 728 alone: **29.11 -> 10.41 core-seconds, 2.80x**; 28.2k -> 79.1k cells
per core-second.

### 2.1 N1 — the destination tree descended one layer too far

`STI.isleaf` for `HierarchicalGridCursor` stops early only when
`bucket_size > 0`, and the default is 0, so a level-12 column rooted at level 5
ended in 823 543 single-cell leaves; the bottom layer alone is 6/7 of the
tree's nodes. Stopping at 49 stored cells (`7^2`, two refinement levels early)
deletes that layer.

Swept on post-c1+c2 code, core-seconds, every point bit-identical:

| column | leaf 1 | leaf 8 | **leaf 49** | leaf 343 |
|---|---:|---:|---:|---:|
| 728 | 20.39 | 13.16 | **12.69** | 18.67 |
| 98241 (all NaN) | 18.27 | 11.66 | **10.97** | 16.60 |
| 115426 (polar) | 20.36 | 13.30 | **12.34** | 18.89 |

**Deviation from the finding as written, and the most important design call of
the campaign.** `2026-08-20-more-improvements.md` proposes the principled home
as a non-zero `bucket_size` default on `PartialGrid`. That is unsafe, and its
own data says so: against a *bare* `HierarchicalGridCursor`, which re-derives a
leaf's cell caps on every visit, the same sweep costs **+25 % at leaf 50 and
+441 % at leaf 350**. The sign of the change depends on whether the caps are
cached. Three uncached paths exist today:

* `GR.celltree(space) = treeify(space.grid)` — bare, no cache;
* `_cachedchunktree` returns the bare cursor above `_CHUNK_CAP_CACHE_MAX` (65 536);
* `_cachedcelltree` falls back to `GR.celltree` for selection cursors.

So the leaf size is attached to **the cap-cached seam, not to the grid**:
`_CACHED_BUCKET_SIZE` is applied by `_bucketed` at exactly the two sites that
return a `CapCachedTree`, and only when the cursor's own `bucket_size` is 0, so
an explicit caller choice still wins and no bare cursor is ever affected. The
production destination (823 543 cells) reaches `_cachedcelltree`, which has no
size limit, so it is cached and gets the benefit.

### 2.2 c1 — memoize the source cursor's derived node extents

`STI.node_extent(::BlockCursor)` re-derived `_node_box` -> `_box_cap` (five
`sincosd` pairs, four spherical distances) for every interior node, every block
build, every column, every worker, over a lattice that never changes.

`MemoBlockCursor` wraps it with a per-task, direct-mapped 1024-slot table in
`task_local_storage`, keyed on `(r0,r1,q0,q1,j0,j1,i0,i1,inpixels)` and cleared
when the task turns to another grid, system or level.

**Deviation:** the spec asked for an LRU keyed by source chunk. The
`MemoRasterTree` shape was used instead — task-local storage is lock-free for
the 21 workers sharing one source `DGGSpace` (the spec's own stated
requirement), memory is O(slots) per task rather than O(tiles), and there is no
eviction policy to get wrong. Verified: `_node_box`/`_leaf_pad` read only `sys`,
`level`, `inpixels` and the eight rectangle fields, so excluding `origin` and
`strategy` from the key is provably correct; and the empty-slot sentinel cannot
collide with a real key because its ninth field is `-1` while `inpixels` is
always 0 or 1.

**-27.2 % standalone, -4.3 % once N1 stands first.**

### 2.3 c2 — hand a leaf's cells back inline

The leaf `(position, cap)` call allocated a fresh `Vector{Tuple{Int,Cap}}` per
leaf: 71.6 % of the base's 3.22 GiB per column.

**Deviation, and a latent-bug finding.** The spec prescribed "a per-task
reusable buffer". Reading the consumer showed that would be wrong:
GeometryOps' `dual_depth_first_search` binds **both** leaves' entry lists and
nests the loops, and `raster_tree_memo.jl` retains the returned vector in a memo
slot. A single shared buffer would therefore silently corrupt any self-join — a
Copernicus grid regridded onto a Copernicus grid, which is supported API.

A lazy view (the overseer's counter-proposal) was implemented and measured at
**-0.9 %**: the join reads the inner leaf's entries once per cell of the outer
leaf, so a view re-derives each cap up to nine times, undoing the hoist the
materialized vector existed to provide.

What shipped is neither: `LeafCells <: AbstractVector`, an `isbits` struct
holding `NTuple{9,Tuple{Int,Cap}}` plus a length. It lives in the caller's
frame, never heap-allocates, is copied by value so aliasing is impossible, and
derives each cap exactly once.

**-46 % allocation standalone; -9.4 % once N1 stands first, with wall inside
noise.** Kept as an allocation change, not a speedup.

### 2.4 c3 — let a narrow wave stand aside

`_fillwave!` wrapped every spawned block build in `@with OUTER_PARALLEL => true`,
suppressing CR's inner threading — measured at **8.7x on 12 threads, 92 %
per-core efficiency** — whenever a column had more than one source chunk, while
the wave itself delivered only 1.05-2.26x.

The fix is in the estimator, not the spawn: `_fillwave!` already has an `i == j`
branch that does *not* set the scoped value, so making `_wavesize` return 1
reaches inner threading with no restructuring. `_wavesize` now estimates
per-chunk cost as `ncells x (cap-overlap(tile, chunk) / chunk cap area)` and
keeps the wave only when its ideal speedup beats `0.73 x nthreads`.

The estimator was validated before being trusted: Spearman 0.70-1.00 against
nnz, where the naive cell-count proxy scores **-0.70, anti-correlated**. A 4x
error in the wide-cap overlap formula was found and fixed against a 4M-sample
Monte Carlo.

**Hazard handled:** when a caller has already set `OUTER_PARALLEL` (the
production worker body does), inner threading is off and the wave is the only
parallelism left, so forcing 1 would serialize the unit. `_wavesize` consults
the scoped value and preserves today's behaviour in that case.

At `-t 12`, `OUTER_PARALLEL` unset:

| col | wall | cores held | cells/core-s | peak RSS |
|---|---|---:|---:|---|
| 728 | 9.56 -> **2.01 s** | 1.48 -> **7.01** | 58.2k -> 58.3k | 2.82 -> **2.21 GiB** |
| 729 | 12.23 -> **2.03 s** | 1.33 -> **7.16** | 50.5k -> 56.6k | 3.09 -> **2.46 GiB** |
| 730 | 10.87 -> **1.88 s** | 1.20 -> **7.16** | 63.1k -> 61.2k | 3.92 -> **2.54 GiB** |
| 115426 | 7.21 -> **1.81 s** | 1.83 -> **7.26** | 62.3k -> 62.7k | 3.71 -> **2.50 GiB** |

Cells per core-second holds, which is what distinguishes real parallelism from
burning cores. Bit-identical at `-t 1` *and* `-t 12`.

**The honest cost.** In the multi-worker regime the new shape spends **~8 % more
CPU per unit of work** (60k vs 65k cells per core-second): inner threading's
92 % efficiency plus pool contention. At a *fixed 24-core budget* the old shape
therefore converts cores into cells slightly faster. c3's value is that it
reaches the same cores with a third of the workers and two thirds of the memory,
and it is the only path past 24 cores. **It is not an ETA improvement at 24
cores and is not counted as one in §1.**

### 2.5 b2 — workers from a core budget

`workers=21`/`threadsper=3` described thread arithmetic the run never had.
`workercount()` now derives `W` from a `cores` budget, and the parallel `shape`
is selectable:

* `outer` — `OUTER_PARALLEL` set in the worker body. 1.06 cores/worker, 65k
  cells per core-second, needs `W ~ cores`, memory scales with `W`.
* `inner` — left unset, so c3's rule applies. Saturates at `W ~ nthreads/4`
  holding ~78 % of the pool, 60k cells per core-second.

`outer` is the default because the standing budget is 24 cores, where the
throughput edge wins; `cores=24 shape=outer` gives **W=23**, independently
reproducing the profile's own `a2` recommendation. Shipping `inner` by default
would have made the production run ~8 % slower at its actual budget.

### 2.6 N3 — inline the Copernicus cell boundary

`cell_boundary` was the only heap polygon left on the hot path; IGeo7's
equivalent already returned an inline `SmallList`. Now `isbits` end to end.

**Allocation -11.6 / -12.0 / -13.0 %** on columns 728 / 98241 / 115426, beating
the -9 % estimate. **Wall: null.** An A/B/A bracket (heap 78.44, 78.17; inline
79.80, 78.72 core-seconds) puts the claimed -2.57 % inside this box's noise.
Kept on allocation alone.

### 2.7 N2+N4 — the GeometryOps repin

Two predicates paid for work they discarded. `spherical_orient` normalized a
cross product and asserted unit-length inputs only to take the sign of a dot
product: **36.5 -> 4.5 ns (8.1x)** once the degeneracy band is tested as
`(n.c)^2 < tol^2 ||n||^2`. The `SphericalCap` intersection test called
`spherical_distance` — an `atan2` — where a squared chord answers:
**16.5 -> 3.4-7.4 ns**.

**-15 % end-to-end, bit-identical, twice**, with allocation unchanged, as
expected of two changes that only remove arithmetic. This is 2.5x the -5.97 %
projected for N2 alone, because N4 attacks the cap-overlap test the profile put
at 24.7 % of CPU.

Both were checked far past "the tests pass":

* `spherical_orient` **keeps** the exact-arithmetic and symbolic-perturbation
  fallback (the report's patch deleted it), now guarded by exactly the criterion
  `robust_cross_product` applies internally. Over **22.2M triples** the two
  implementations show **zero sign flips** — they never disagree about which
  side a point falls on. The 839 708 disagreements are all in deliberately
  constructed boundary suites and all within **±5 % of the tolerance**, which is
  inherent: `tol = 16 eps ~ 3.55e-15` while a three-term dot product's rounding
  noise is ~2e-16.
* The new cap test is **more accurate than the old one**. `spherical_distance`'s
  `cross(p,q)` cancels catastrophically for nearly-equal centers. Adjudicated in
  256-bit arithmetic over the 111 514 knife-edge pairs where the two differ,
  **the old test is wrong on all of them and the new one on none.**
* The sketched N4 patch was benchmarked at **0.71x — slower** — because random
  cap centers sit ~pi/2 apart so its `c <= 1` gate always missed. It was
  restructured to be entirely `sqrt`-free.

GeometryOps' own suite is unchanged: **231 527 pass / 0 fail, per-testset
identical** on both sides.

**N6 was declined** with cause: `find_orthogonal` does not normalize, so the
divisions in `_sh_spherical_intersection` are load-bearing for the `a == b`
case. The principled fix is upstream and independently reviewable, and N6 has no
measurement behind it to justify carrying that risk here.

---

## 3. Deviations from the spec, collected

| # | spec said | what was done | why |
|---|---|---|---|
| 1 | c3/c2 "may live in the pinned CR repo; branch CR and repin" | No CR change at all | `_fillwave!`/`_wavesize` are in `lib/GlobalRegridding`, a DGG workspace member. CR was never touched; its pin is unchanged. |
| 2 | c1: LRU keyed by source chunk | Per-task direct-mapped memo | Lock-free for 21 shared readers, O(slots) not O(tiles), no eviction policy. §2.2 |
| 3 | c2: per-task reusable buffer | `isbits` inline `LeafCells` | The prescribed buffer would corrupt self-joins. §2.3 |
| 4 | N1: non-zero `PartialGrid` default | Attached to the cap-cached seam | A global default regresses uncached cursors by +25 % to +441 %. §2.1 |
| 5 | c4 as a planned step | See §5 | Conditional on the re-profile, per the spec's own wording. |
| 6 | Manifest tree-sha1 beside repo-rev | `Manifest.toml` is gitignored here | The pin travels in `Project.toml`; the tree-sha1 is recorded in the commit message. Method validated by recomputing the *old* pin's tree and matching the existing entry. |

Not attempted, per instruction: c5, N5, and the measured dead ends.

---

## 4. Where the time goes now

`Profile.@profile`, 2 ms sampling, four mid-latitude column builds, 21 896
samples = 43.8 s CPU ~ 10.95 s per column, aggregated topmost-inclusive by
function -- the same methodology as the original profile. Load 48-52.

| frame | orig % | orig s/col | **new %** | new s/col |
|---|---:|---:|---:|---:|
| `intersection_areas` (CR) | 93.0 | 26.97 | **82.95** | 9.08 |
| ├ `assemble_sparse_matrix_coo` / `BlockAreaOperator` | 24.7 | 7.16 | **50.32** | 5.51 |
| │  ├ `_sh_clip_spherical!` (the clipper itself) | — | — | 32.23 | 3.53 |
| │  └ `_memocell`/`cell_polygon`/`cell_boundary` | — | — | 13.77 | 1.51 |
| ├ `get_all_candidate_pairs` / `dual_depth_first_search` | 67.8 | 19.66 | **28.98** | 3.17 |
| │  ├ `child_indices_extents`, both trees | 21.3 | 6.18 | **11.80** | 1.29 |
| │  │   ├ source `BlockCursor` (`LeafCells`) | — | — | 6.17 | 0.68 |
| │  │   └ destination `CapCachedTree` (heap `Vector`) | — | — | 5.63 | 0.62 |
| │  ├ `node_extent`, source `BlockCursor` | 47.7 | 13.83 | **8.28** | 0.91 |
| │  │   └ interior nodes only (what c1 memoizes) | — | — | **1.73** | 0.19 |
| │  │   └ `_box_cap` | 41.9 | 12.15 | **6.26** | 0.69 |
| │  ├ `_intersects` (squared chord, N4) | — | — | 5.64 | 0.62 |
| │  └ `spherical_distance` | 24.7 | 7.16 | **4.72** | 0.52 |
| destination cap-tree build (`_leafcaps`) | ~5 | 1.45 | **13.79** | 1.51 |
| *(cross-cutting)* `snyder_inv_xyz` | — | — | **21.70** | 2.38 |

### The new top-3 sinks

1. **Sparse assembly / `BlockAreaOperator` — 50.3 %** (was 24.7 %). It doubled
   in share while falling 1.7 s absolute. Its internals split cleanly: 32.2 % is
   the Sutherland-Hodgman spherical clipper — genuinely useful work, the thing
   the program exists to do — and 13.8 % is deriving the destination cell
   polygon.
2. **Pair finding — 29.0 %** (was 67.8 %). This is what the campaign aimed at:
   **-16.5 s of 19.7 s, -84 % absolute.** Its internals inverted. The source
   tree's interior-node geometry, formerly the single largest cost in the whole
   program at 47.7 %, is now **1.7 %**.
3. **The destination cap-tree build — 13.8 %** (was ~5 %). Absolutely unchanged
   at 1.51 s/column; it rose in share only because everything around it shrank
   2.7x.

The original headline is now false in its own terms. "Two-thirds of every CPU
second is spent finding candidate pairs, and half of that recomputing spherical
caps for source-tree nodes that never change" reads today as: pair-finding is
**29 %**, and source cap derivation within it is **8 %**.

### What is left, and what changed character

* **`snyder_inv_xyz` — 21.7 % of all CPU** is the largest single remaining
  lever, and it is *not* in either phase this campaign optimized. It is reached
  from two independent places: the destination-polygon memo under
  `BlockAreaOperator` (10.8 %) and `cell_cap` under `_leafcaps` (10.5 %). The
  IGeo7 inverse Snyder projection is now the hot kernel.
* **Runtime dispatch: 0.0137 % self.** Still zero.
* **GC: 10.7 % of samples, but 1.31 -> 1.18 s per column absolute.** Its share
  doubled only because the denominator shrank. (`@timed` reads 4.5-8.4 % on the
  same columns; the profiler figure includes allocation-path frames the `@timed`
  counter excludes. Both are reported rather than picking the flattering one.)
* **`GenericMemory` boot.jl:588 was 7.0-7.5 % self in the original and is now
  below 0.4 %** — c2 and N3 visible directly in the self-cost table.

### Allocation: same dominant type, opposite side of the join

**1.874 GiB per column, 2 443 B per destination cell** (was 3.22 GiB / 4 198 B),
**-42 %**. The dominant type is still
`Memory{Tuple{Int64, SphericalCap{Float64}}}` — but its origin has flipped.
Originally 71.6 % of all bytes were the *source* `BlockCursor`'s per-leaf
vector; that is now **exactly zero**. Every remaining byte of that type comes
from the *destination* `CapCachedTree.child_indices_extents`
(`src/cap_cached_tree.jl:52`), which still builds a heap `Vector` per leaf
visit: **~0.75 GiB/column, ~40 % of all allocation, ~675 000 allocations.**

---

## 5. c4 and `LEAF_CELLS`: both dropped, with numbers

### c4 — drop

Source-side `child_indices_extents` is **6.17 % of CPU** (0.68 s of a 10.95 s
column), of which the memoizable part — the cap derivation, not the tuple
assembly — is **5.07 %**. The original target was -10 to -15 %. A memo that were
perfectly free with a 100 % hit rate could not exceed **-6.2 %**.

The spec's own condition ("only if the re-profile still shows
`child_indices_extents` as the next sink") is **not met**: it fell 21.3 % ->
11.8 %, and its source half now ranks behind the clipper (32.2 %),
`snyder_inv_xyz` (21.7 %), the destination cap-tree build (13.8 %) and the
destination-polygon memo (13.8 %). N1 and c2 already collected most of what c4
was aiming at. **Not implemented; recorded as a measured null.**

**The better-evidenced successor** is the *destination* half: 5.63 % of CPU and
~0.75 GiB/column (40 % of all allocation) from the `Vector` at
`cap_cached_tree.jl:52`. Its caps already sit in a flat vector, so a view or an
inline buffer needs no derivation at all — c2's treatment applied to the other
side of the join. It would not cut much CPU but would remove the single largest
allocation source in the program.

### `LEAF_CELLS` — keep 9

Swept with a `LEAF_CELLS = 9` run first and last to bracket load drift; the two
bracket runs differ by **0.9 %**, so effects above that are real.

| `LEAF_CELLS` | `LeafCells` inline size | Σ min wall | vs 9 | Σ alloc | vs 9 |
|---:|---:|---:|---:|---:|---:|
| 4 | 160 B | 29.32 s | −0.25 % (noise) | 6.87 GiB | +30 % |
| **9 (current)** | **368 B** | **29.26 / 29.52 s** | — | **5.28 GiB** | — |
| 16 | 648 B | 30.77 s | **+4.70 %** | 4.62 GiB | −12.6 % |
| 25 | 1 008 B | 31.27 s | **+6.39 %** | 4.24 GiB | −19.9 % |
| 49 | 1 960 B | 33.78 s | **+14.93 %** | 3.96 GiB | −25.1 % |

Monotone in all three columns independently. **No change recommended.**

**N1's lesson does not transfer to the source side — and the reason confirms the
mechanism rather than contradicting it.** N1 worked because the destination tree
*caches* its leaf caps in a flat vector, so deleting a tree layer removes work
at no per-leaf cost. The source `BlockCursor` *derives* its caps on demand, so a
fatter leaf makes every visit strictly more expensive. This is precisely the
uncached regime where N1's own control measurement showed **+25 % at leaf 50**.
The two curves prove the mechanism: allocation falls monotonically (−25 % at 49,
i.e. materially fewer leaf-pair visits) while wall rises monotonically
(+14.9 %). The loss is per-visit compute, not GC pressure.

**On the `LeafCells` coupling, explicitly:** at 49 the inline struct is 1 960 B
copied by value per visit, and that is part of the +14.9 %. But redesigning
`LeafCells` would unlock nothing, because the allocation curve shows the loss is
not allocation-driven, and because part of it is tree shape — fewer, fatter
leaves prune less and brute-force more, which no container change fixes. A
redesign would have to recover more than 4.7 % just to reach parity with 9.

---

## 6. Reproducing

Harness in `/home/asinghvi17/geo/dggstores/ladder-scratch/`, regenerated from
the production script:

```bash
S=/home/asinghvi17/geo/dggstores/ladder-scratch
head -n -1 /home/asinghvi17/geo/DGG-perf-ladder/scripts/copdem_production.jl > $S/prod.jl
cd /home/asinghvi17/geo/DGG-perf-ladder
RASTERDATASOURCES_PATH=/home/asinghvi17/geo/DiscreteGlobalGrids.jl/bench/data \
  LADDER_TAG=mytag julia --project=bench -t 1 $S/verify.jl
```

| script | what it does |
|---|---|
| `setup.jl` | builds the production world against a scratch store; exports `world()`, `lazycolumn`, `timecolumn`, `withcpu`, `fingerprint`, `emit` |
| `verify.jl` | 8 reference columns vs frozen bit-exact fingerprints; errors on any difference |
| `m0_baseline.jl` | the original baseline run |
| `m1_bucket.jl` | the `bucket_size` sweep (caller-side, no library change) |
| `m2_alloc.jl` | wall + allocation, 3 reps, min wall |
| `baseline-fingerprints.txt` | the bit-identity oracle |

The live production run (PID 2915199), its store, log, done log and column
cache, and the worktree `/home/asinghvi17/geo/DGG-subzone-store` were not
touched at any point.
