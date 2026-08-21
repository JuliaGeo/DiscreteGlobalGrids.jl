# Brief — evidence gathered for the IGEO7 L13 x Copernicus GLO-30 checkpoint

Input document for the optimisation-recommendation pass. Everything here is
measured or quoted code, gathered 2026-08-18 on a 64-core Linux box, Julia
1.12.6, repo at `79e1db0` plus an uncommitted `bench/` workspace member.
Sections marked TODO are filled in when the running campaigns land.

## Reading these numbers: the measurement machine has 751 GB of RAM

**Every figure in this brief was taken with ~280 GB free and no memory pressure**, so the
heap grew whenever it wanted and the GC almost never ran (5.9% of the lazy read, ~0 on the
eager path). Three consequences that must not be lost when reasoning about ordinary hardware:

1. **Allocation churn and peak RSS are different numbers and only one of them was binding
   here.** Baseline eager at the full tile churned **4.89 TB** at a peak RSS of **8.0 GB**.
   On this box the churn was nearly free. On a machine where the heap cannot grow, that same
   churn is collected against a small heap and converts directly into GC time -- so the
   *time* results here are an optimistic bound for smaller machines, by an unmeasured factor.
2. **Peak RSS is the figure that transfers**, and it is the one this brief is thinnest on.
   Where a number below is quoted in GB of *allocation*, it says nothing about whether the
   job fits elsewhere.
3. **Some blockers are absolute and bite everywhere**: the global `DGGSpace` at 29-30 GB
   resident (D4), and the 51 TB global weight matrix. Those do not depend on this machine.

`budget` is the pipeline's nominal memory control (`plans.jl:81-95`, default 2^30, splitting
25% weights / 75% data). **Whether it actually bounds peak RSS at full-tile scale is
measured in the memory-envelope study, not assumed here** -- earlier checkpoints
(`perf-P1.md` L3) found it carried and ignored, which `perf-P2.md` reports as since fixed.

## The workload

| | |
|---|---|
| source | Copernicus GLO-30 tile `N45_00_E010_00`, 45-46N / 10-11E (Lake Garda to the Adamello) |
| | `Matrix{Float32}` 3600 x 3600 = 12 960 000 px, 1 arcsec both axes, no nodata, 7.00-2805.94 m |
| destination | `IGeo7System()` level 13, `MultiOrderCoverage` of the tile footprint |
| | expands to **16 490 405** cells; global level 13 is 968 890 104 072 cells |
| area match | L13 mean cell 526.44 m2 against a 668.69 m2 pixel at 45.5N -- L13 is the
closest level and is 1.27x finer than a pixel. L12 is 5.51x coarser, L14 8.9x finer. |
| harness | `bench/harness.jl`, `bench/moc_probe.jl`, `bench/fetch_tile.jl` in the `bench/` workspace member |

## What is healthy

* **The coverage representation.** `MultiOrderCellSet` holds the whole 16.49 M-cell
  destination in 45 185 mixed-level entries / 728 792 B; `CellVector` adds
  `RangeWindows` at 24 B per window and **0 B per cell** (10 085 windows, 1.05 MB);
  `PartialGrid` shares them by reference (323 984 B). Materialised ids would be
  131 923 280 B. **~0.02 B/cell, a 400x saving**, and `strictly_increasing` is
  short-circuited to `true` (`cell_vector.jl:397`), not an O(n) scan.
* **The MOC query beats the flat query decisively** at depth: L13 is 2.77 s /
  1.77 GiB for the multi-order set against 59.07 s / 20.7 GiB for
  `query(levelgrid, Intersects(box))`. 21x faster, 12x less allocation, 181x
  smaller answer, and the two agree to 41 cells in 16.49 M (the MOC being the
  superset, from IGEO7 non-congruence).
* **Z7 index arithmetic.** `Z7Cell` is a bare `UInt64`, 4 base bits + 20 x 3 digit
  bits = exactly 64. No `Int128`/`BigInt`/`widen` anywhere in `src/`. `NUM_CELLS`
  at level 19 is 1.14e17, **81x inside `Int64`**. `descendant_range` is 20.3 ns,
  zero allocations, `@inferred` clean.
* **Type stability.** Every measured hot entry point (`cell_boundary`,
  `cell_polygon`, `cell_centroid`, `cell_area`, `node_extent`, `descendant_range`,
  `STI.child_indices_extents`, `Trees.getcell`, `treeify`) is concretely inferred;
  `code_warntype` shows no `::Any`, no non-trivial `Union`, no `Core.Box`. The three
  real instabilities are `cellposition -> Union{Nothing,Int}`,
  `STI.getchild(::WindowCursor)` (3-way union, union-split not dispatch), and
  `CopernicusDEM.treeify -> Union{BlockCursor,HierarchicalGridCursor}` (once per
  traversal). Cost here is allocation and redundant work, not dispatch.
* **Clip-window convexity.** Every IGEO7 cell handed to Sutherland-Hodgman is a
  6-vertex CCW spherically-convex hexagon (5 at the twelve pentagons) at **every**
  level 6-14 -- vertex count does not vary with level. The reflex-window hazard
  documented in `examples/copernicus_dem.jl` belongs to HEALPix, not this pipeline.
  Conservation measured on real data: max relative column error **6.09e-10**, total
  closure -2.86e-12, max interior coverage error 1.18e-10.
* **Serial scaling is linear.** Over 64x in problem size the exponent is
  alpha -> 1.034; nnz is exactly linear in source cells at 5.128 nnz/pixel,
  4.0 nnz/destination cell, 1.404 candidate pairs per nonzero (a 71% tree-filter
  hit rate).

## The defects, ranked by measured cost

### D1 -- quadratic merge in the threaded intersection path (UPSTREAM, FIX SUBMITTED)

`ConservativeRegridding/src/utils/MultithreadedDualDepthFirstSearch.jl:66` and
`src/regridder/intersection_areas.jl:218-220` merge per-task results with
`reduce(vcat, x; init = T[])`. Any `init` bypasses
`reduce(::typeof(vcat), ::AbstractVector{<:AbstractVecOrMat})` -- the method that
pre-sizes the output and copies each part once -- and falls back to `mapfoldl`,
which re-copies the accumulator per task: **O(N x ntasks)**.

Measured at the realistic shape (20 736 tasks x 7 pairs): 4 696 ms / 22 966 MiB
with `init`, **0.195 ms / 2.37 MiB** with an emptiness guard. 24 000x.

End-to-end on real IGEO7 x GLO-30 pairs, direct `intersection_areas`:

| src cells | 1 thread | 8 threads | 64 threads |
|---:|---:|---:|---:|
| 5 184 | 3.21 s / 807 MiB | 0.83 s / 2.3 GiB | 1.00 s / 2.4 GiB |
| 20 736 | 10.07 s / 2.0 GiB | 3.78 s / 25.7 GiB | 4.31 s / 26.1 GiB |
| 82 944 | 38.90 s / 7.4 GiB | **52.54 s / 385 GiB** | **62.30 s / 387 GiB** |

Confirmed on the production eager path by `perf-P3.md`: `DirectPlan` at 512 x 512 px is
**17.62 s / 4.64 GiB at one thread and 54.16 s / 191.95 GiB at eight** -- 3.1x slower and
41x more allocated for a bit-identical block. It also means the eager and lazy paths are
**level at one thread** (17.62 s vs 17.47 s), so the eager path's apparent 51 s cost was a
threading regression, not an algorithmic one. The lazy path's blocks are small enough not
to trip it.

Fixed in https://github.com/JuliaGeo/ConservativeRegridding.jl/pull/133, **merged to CR
`main` as `55820f1`**. Measured A/B on this exact workload, patched arm = the pinned tree
plus the PR's two hunks verbatim (so the clipping cache is held constant; `Pkg.resolve`
reported zero version drift). All 30 paired runs bit-identical.

**Eager `DirectPlan`, `-t 8` -- this is where the fix pays:**

| window | base | fixed | speedup |
|---:|---:|---:|---:|
| 512 x 512 | 51.03 s / 191.9 GB | 3.38 s / 4.7 GB | **15.1x / 40.7x** |
| 1024 x 1024 | 247.59 s / 353.6 GB | 11.71 s / 17.3 GB | 21.1x / 20.4x |
| **3600 x 3600 (full tile)** | **3268.75 s / 4892.7 GB** | **140.86 s / 206.5 GB** | **23.2x / 23.7x** |

Fitted exponent vs pixel count goes **1.199 -> 0.939**: super-linear to sub-linear.

**This reverses the standing conclusion that the lazy path is the tractable one.** Patched
eager at the full tile is **140.86 s against lazy's 163.52 s** -- the eager path becomes the
*fastest* option at every size from 512 x 512 up. The cliff was time and allocation churn,
never peak RSS (baseline eager peaked at 8.0 GB even while churning 4.89 TB; nothing OOMed
at any size).

**Lazy path, full tile, `-t 8`: the fix buys ~1% of wall time** (165.49 -> 163.52 s) and 11%
of churn. Its value on the lazy path is at high thread counts, where it removes a genuine
*regression*: baseline speedup peaks at **5.16x at 16 threads then falls back to 4.76x at
64**, while patched is monotone to **5.99x**. Baseline allocation grows 16.7 -> 41.3 GB with
thread count; patched is flat at ~17 GB. At `-t 64`: -13% wall, -49% bytes, peak RSS
**8.44 -> 3.43 GB**.

This refines `perf-P3.md`'s R6, which predicted "no effect on the lazy path": true at
<= 16 threads, false at 32-64.

**Why the two paths differ**, and it is the two patched sites having different drivers:
`intersection_areas.jl:218` partitions at `npartitions = nthreads() * 4`, so its excess is
**thread-driven** and size-independent -- the lazy path's mild high-thread penalty.
`MultithreadedDualDepthFirstSearch.jl:66` spawns per tree-node pair, so its excess is
**problem-size-driven** -- the eager path's blow-up. Isolated, the primitive's cost is exactly
linear in chunk count `k` while the guarded form is flat: at N=1e7, k=8 gives 6x, k=512 gives
344x, k=4096 gives **2238x**.

Second-order contributor, **not** fixed and open for recommendation:
`MultithreadedDualDepthFirstSearch.jl:19-21` spawns a task unconditionally in the
`isleaf` branch, bypassing `should_parallelize` entirely -- measured **one task per
source cell** (20 736 tasks where the tree's own policy asked for ~256).

### D1a -- the dual-DFS spawn test uses `&&`, so one tree is split to its leaves

**Separate from D1's merge defect, in the same file.**
`ConservativeRegridding/src/utils/MultithreadedDualDepthFirstSearch.jl:24`:

```julia
if parallelize1(node1, ext1) && parallelize2(node2, ext2)
```

The conjunction means one side reaching its parallelize granularity does nothing until the
other agrees. When the two trees bottom out at different depths -- exactly what a quadtree
cursor meeting `HierarchicalGridCursor` does, since `_stored_count` barely drops for the first
several levels of a regional window -- the walk keeps splitting the tree that is **already**
fine enough, down to its leaves, where the unconditional `isleaf` spawn at `:19-21` then fires
once per leaf pair. Measured **18 944 tasks for 20 736 source cells**, against a policy target
of `nthreads * 32`.

The `isleaf` branch is the **symptom**; the conjunction is the cause. Fixed by `&&` -> `||` in
https://github.com/JuliaGeo/ConservativeRegridding.jl/pull/134 (open), on the argument that a
node covering <= 1/K of its own tree can only be in <= 1/K of the candidate pairs regardless
of the opposing node.

| threads | tasks before -> after | ms before -> after |
|---:|---:|---:|
| 1 | 18 944 -> 64 | 6 947 -> 4 243 |
| 8 | 18 944 -> **256** | 3 375 -> **686** |
| 64 | 18 944 -> 8 192 | 3 563 -> 1 299 |

Results identical **including order** (20.4 M planar pairs plus the GLO-30 case verified
element-for-element against serial and against `main`). Suite 10 032 / 0 failures.

Known remaining case, documented in the PR rather than papered over: a tree that genuinely
**is** a single leaf still spawns one task for everything (a 2x2 source against 1024x1024 =
one task for ~1 M pairs). Two candidate fixes for it were built and rejected on evidence --
one-sided descent loses pairs (see D1b), and its leaf-restricted variant measured 2.5x worse
(30 722 tasks / 7 466 ms against 18 944 / 2 973 ms).

Behaviour change to weigh: under `||` a policy that fires early caps the task count for
**both** trees, and the `should_parallelize(::Any, ::SphericalCap)` quarter-sphere fallback
fires at the root for any regional tree -- so a foreign tree relying on that fallback paired
with a global grid now spawns a single task at the root.

### D1b -- RESOLVED, NOT A DEFECT: cap non-nesting is a red herring

Verified exhaustively (`bench/nesting-study/`). **The shipping path is sound and loses nothing
with area.**

Cap non-nesting is real -- one-step child-cap escapes are 361/25 452 on the IGEO7 cursor
(max 2.13e-6 rad) and **100%** on `BlockCursor` and `RasterCellTree` by construction, since a
circumscribed box cap of a split rectangle is always wider than half its parent's. But it is
slack against slack: a child's cap over-approximates the child's geometry, so it pokes past a
parent cap tightened to a vertex hull in directions where no vertex lies.

**The invariant the traversal actually needs -- every descendant leaf polygon inside every
ancestor extent -- holds with 0 violations in 928 661 (ancestor, leaf-vertex) checks** across
all four cursor types.

**Why nesting is never asked:** every `predicate` call in GeometryOps' dual DFS
(`dual_depth_first_search.jl`) and CR's threaded variant is **cross-tree**. No extent is ever
compared with an extent from its own tree. The `leaf1`/`leaf2` branches compare a node extent
against an *opposing* child extent -- two independent covering claims, not a nesting claim.

Brute-force ground truth, shipping traversal vs all-pairs:

| window | brute-force candidates | shipping (serial == threaded) | missed | missed with area |
|---:|---:|---:|---:|---:|
| 64x64 | 29 457 | 29 457 | 0 | -- |
| 128x128 | 117 899 | 117 899 | 0 | -- |
| 256x256 | 471 612 | 471 611 | 1 | **0** (true area 0.0) |
| 512x512 | 1 887 405 | 1 887 405 | 0 | -- |

End-to-end through the full lazy stack the weight matrix is **bit-identical** to exhaustive
(`nnz 83987 == 83987`, `max |B-A| = 0.0`), and the exhaustive matrix carries the **identical**
6.09e-10 column residual -- so that residual is entirely Sutherland-Hodgman clipping and area
quadrature, not lost pairs.

Arc soundness confirmed by measurement: the widest non-whole-sphere extent in any real
traversal is **0.783 rad**, far under pi/2, so vertex containment implies polygon containment
everywhere it is used.

**`cap_inflation` is derived, not chosen, and no system is too tight** -- but S2 is the one to
watch:

| system | inflation | descendant/cell ratio | used/extent radius | margin |
|---|---|---:|---:|---|
| IGeo7 | 1.2 | 1.0367 | 0.864 | ~14% |
| H3 | 1.2 | 1.0520 | 0.877 | ~12% |
| A5 | 1.75 | 1.3537 | 0.774 | ~17% |
| HEALPix / ISEA4R | overridden | 1.0000 | 0.928 / 0.935 | ~7% |
| **S2** | overridden | 1.0000 | **1.0000** | **one ULP** |

S2's extent is exactly tight by design (children tile the parent, the max is at a shared
corner, `nextfloat`-guarded, `src/systems/S2/system.jl:242-252`). The argument is rigorous, but
the whole margin is smaller than the relative error of the `acos`-based `spherical_distance`
computing both sides -- the one system where a change to corner or centre arithmetic could flip
a `<=` with no slack to absorb it.

**Conformance recommendations** (`lib/DiscreteGlobalGridsConformanceTesting`): `covering_law_problems`
already tests the right property in the right form (descendant boundary against *every* ancestor
extent on the path, plus a `radius <= pi/2` assertion). Gaps are reach, not shape:
(1) defaults `descent_depth = 3` / `branch_samples = 2` measure ~2% below truth -- the ratio is
still climbing at depth 5->6; raise to >= 6 and >= 3.
(2) sample the **pentagon one-rings**, not the pentagons (see D2b).
(3) **nothing conformance-tests the cursor layer at all** -- `covering_law_problems` tests
`DGG.node_extent(sys, c)` and never `STI.node_extent(::HierarchicalGridCursor)`, `BlockCursor`
or `RasterCellTree`. Promote the walk into a `tree_covering_problems(tree, grid)`.
(4) **do not assert cap nesting** -- it would fail on `BlockCursor`/`RasterCellTree` at 100% of
nodes, correctly, for a benign reason.

### D2 -- IGEO7 has no O(1) `node_extent`, so cell geometry is re-derived ~10.7x per cell

`src/fallbacks/geometry.jl:318-319` falls back to
`inflate_cap(cell_cap(levelgrid(sys, level(c)), c), ...)` and `caps.jl:81`
`cell_cap(grid, c) = points_cap(cell_boundary(grid, c))` -- a full inverse
projection per tree-pruning step. One `cell_boundary` at L13 = ~104 transcendentals
(~14 `atan`, ~20 `acos`, ~26 `sincos`, ~44 `sqrt`, Newton loop averaging 3.36
iterations), **1.215 us, 4 allocations, 416 B**.

**Authoritative count, from `perf-P3.md` on the production lazy path** (1024 x 1024 px
-> 1 331 572 cells, `RasterGrid` source, `chunklevel = 7`, `-t 1`):
**14 206 931 `cell_boundary` calls for 1 331 572 cells = 10.67x**, costing **19.7 s of a
61.17 s read** against **2.204 s** for one pass over every cell (8.96x in time).

| path | share of calls | file:line |
|---|---:|---|
| `cell_cap` | 50.2% | `src/fallbacks/cursor.jl:237-248`, `:265-273` |
| `cell_polygon` via `Trees.getcell` | 45.1% | `src/fallbacks/cursor.jl:291-295` |
| the system's default `node_extent` | 4.3% | `src/fallbacks/geometry.jl:318` |
| `cells_cap` sparse-node tightening | **0.40%** | `src/fallbacks/cursor.jl:243-246` |

**This supersedes an earlier count of 43.2x with "55% via `cells_cap`"**, which was taken
on a different configuration -- a direct `CR.Regridder` build over a much smaller window
(28 800 px -> 49 414 cells) with the source built as `PartialGrid(sys, level, Vector{LevelIndex})`
(a `HierarchicalGridCursor`) rather than the production `RasterGrid`. Both were counted, not
inferred; **P3's is the one that describes the path that actually runs**, and `cells_cap` is
0.14% of wall clock, not a lever. The redundancy factor does still grow with tile size:
`cell_polygon` calls per cell go **1.00 -> 3.76 -> 4.81** across 256/512/1024, tracking tiles
crossing `_TILE_CELL_CACHE_MAX`.

At the leaf, `node_extent` (`cursor.jl:241`) *is* `cell_cap(grid, id)` and
`child_indices_extents` (`cursor.jl:272`) recomputes the same cap for the same cell --
10.42% + 8.58% of the read for one cap. `cursor.jl:188-192` carries the
written-but-unacted-on diagnosis; `cursor.jl:262-263` claims a materialisation "to avoid
recomputing caps during repeated dual-tree passes" that is per-call and so caches nothing
across passes.

The circumradius is closed-form in dev units (`abs(CELL_SCALE[res+1])/SQRT3`, already used
at `z7grid.jl:215`), so an O(1) `node_extent` from centroid + per-level radius is derivable.

### D2a -- `cell_cap` is not an overridable hook, and it carries 50.2% of the calls

**The structural blocker behind D2, and the reason "just give IGEO7 an O(1) `node_extent`"
does not work.** `cell_cap` has **exactly one method in the whole package** (verified via
`methods(Fallbacks.cell_cap)`):

```julia
cell_cap(grid::AbstractGrid, c::AbstractCellIndex) = points_cap(cell_boundary(grid, c))
```
`src/fallbacks/caps.jl:81`.

Per `perf-P3.md`, **50.2% of the 14.2 M `cell_boundary` calls enter through `cell_cap`** and
only **4.3% through `node_extent(sys, c)`**. The 50.2% enters at `src/fallbacks/cursor.jl:241`
(`cursor.level >= cursor.leaf_level && return cell_cap(cursor.grid, cursor.id)`) and `:272`
(`child_indices_extents`). **Neither site can be reached by any system-level `node_extent`
override that exists today** -- not even S2's. Overriding `node_extent` alone therefore fixes
at most 4.3%.

Worse, `cursor.grid` is a `PartialGrid`, and **neither `PartialGrid` nor `AuthalicGrid`
forwards `cell_cap`**. Confirmed by dispatch: all three of `cell_cap(levelgrid, c)`,
`cell_cap(partialgrid, c)` and `cell_cap(authgrid, c)` resolve to `caps.jl:81`. A
specialisation written on `HierarchicalLevelGrid{IGeo7System}` would be **silently bypassed**.

Prerequisites, in order: promote `cell_cap` to a documented grid-level interface method with
the boundary fallback retained; add `cell_cap(::PartialGrid, c) = cell_cap(grid.complete, c)`
(mirroring `partial_grid.jl:142-143`); add `cell_cap(::AuthalicGrid, c)` by the Lipschitz rule.
`cells_cap` is 0.14% of wall clock -- skip it.

### D2b -- the analytic cap constant: K = 1 UNDER-bounds by 13.6%

The obvious O(1) cap is `SphericalCap(cell_centroid(c), K * circumradius(level))` with the
dev-unit circumradius `abs(CELL_SCALE[r+1])/SQRT3` (already used at `z7grid.jl:215`).
**Measured 283.7 ns / 0 B against 1542.0 ns / 416 B for `cell_cap` -- 5.44x and
allocation-free.** 250 of those 284 ns are `cell_centroid` itself; `SphericalCap`
construction is ~33 ns, the radius lookup 2.1 ns.

**But `K = 1` is unsafe.** Ratio `dev_circumradius(level) / max_vertex_arc(centroid)`,
**exhaustive over every cell of levels 0-7** (9.4 M cells at L7) plus sampling to level 19:

| level | ncells | min ratio | max | K needed |
|---:|---:|---:|---:|---:|
| 4 | 24 012 | 0.870781 | 1.029199 | 1.148 |
| **6** | **1 176 492** | **0.863869** | 1.028180 | **1.157583** |
| 7 | 8 235 432 | 0.873278 | 1.115599 | 1.145111 |

Worst case **0.863869 at level 6, `Z7Cell("00463106")` -- a 13.6% under-bound on an ordinary
interior hexagon**. Level-stable from L6 down. **The twelve pentagons are the *safest* cells**
(ratio ~1.1156 at every level >= 3); the defect is Snyder anisotropy on ordinary hexagons, not
the pentagons -- so intuition about the icosahedral defects points the wrong way here.

**Independently confirmed for `K_node`, with a sampling warning:** the covering ratio reaches
~1.179 at depth 6 and is still climbing slowly. Pentagons sampled **alone** give only
0.92-0.94 -- **the worst cells are the hexagons in the pentagons' one-ring**, so any sampling
strategy that reaches for pentagons systematically *under*-measures `K_node`. `K = 1.25` gives
zero child-cap escapes and zero descendant-vertex escapes at every level 0-8; `K = 1.0526`
already fails with 240 escapes at level 1 -- nesting is a threshold property of the extent's
*form*, not something "analytic" confers automatically. Caveat: this fixes nesting only at the
system `node_extent` layer; `STI.node_extent(::HierarchicalGridCursor)` would still tighten
sparse nodes to `cells_cap` and still return `cell_cap` at leaves.

An under-bound **silently drops intersections** with no diagnostic. Safe constants, measured:

* **`K_cell = 1.18`** for the leaf cap (requirement 1.157583, 1.9% margin).
* **`K_node = 1.25`** for the subtree extent -- from the directly measured covering law over
  descendants to depth 6, converged at ~1.185 (increment ~0.001 at d=5->d=6). Here the
  icosahedral defects **do** drive the bound (worst cells `Z7Cell("10000004")`,
  `Z7Cell("04000006")`, all-zero prefixes).

Cost in cap area (~= extra candidate pairs): leaf cap **+23% median / +31% worst**; but
`K_node = 1.25` is **4% tighter than the status quo**, because the status quo is
`1.2 x points_cap` and `points_cap` centres on the normalised vertex mean rather than the
centroid. Net: interior pruning improves, leaves get slightly worse, no boundary is built.

Conformance: add the covering-law check at depth **>= 4** -- depth 1 alone reports 1.151 and
would pass a wrong constant.

### D2c -- authalic: fine at system level, broken at grid level

**System level inherits for free.** `node_extent(sys::AuthalicSystem, c)`
(`src/fallbacks/authalic_grid.jl:284-291`) takes the inner system's cap, re-centres at the
warped centre, multiplies the radius by the Lipschitz constant `authalic_stretch`
(**1.0044899814** for WGS84), and drops to the full sphere past pi/2. It never touches a
boundary. Measured `node_extent(AuthalicSystem(S2System()), c)` = **304.9 ns / 0 B** against
S2's own 243.3 ns -- ~60 ns of wrapper overhead and full inheritance of an O(1) inner extent.

**Grid level falls through.** No `cell_cap(::AuthalicGrid, ...)` exists, so it reaches
`caps.jl:81` -> `cell_boundary(::AuthalicGrid, c) = map(geodetic_point, inner_boundary)`:

| | IGeo7 L13 | Authalic(IGeo7) | S2 L13 | Authalic(S2) |
|---|---:|---:|---:|---:|
| `cell_cap(grid, c)` | 1777 ns / 416 B | **2110 ns / 624 B** | 278 ns / 160 B | **541 ns / 320 B** |
| `node_extent(sys, c)` | 1794 ns / 416 B | 1849 ns / 416 B | 243 ns / 0 B | **305 ns / 0 B** |

The wrapper is a strict improvement at the system level and a strict penalty at the grid
level -- and the grid level carries the 50.2%.

**The generalisable rule for wrapper authors:** a wrapper applying a pointwise map `Phi` keeps
the inner system's extent asymptotics iff it supplies a global Lipschitz constant `L` with
`d(Phi p, Phi q) <= L d(p, q)`; then the extent is `SphericalCap(Phi(centre), L*radius)`,
degrading to the full sphere when `L*radius > pi/2`. Supply `L` as a function proved from the
differential, not a table. `authalic_stretch` is the worked example, and its docstring records
why a multiplicative bound is required: the additive alternative puts a 0.13 degree floor
under every extent, swamping the cell from about level 8 down.

### D2d -- `node_extent_is_expensive` is hard-coded true for every cursor

`src/fallbacks/cursor.jl:192` declares `STI.node_extent_is_expensive(::Type{<:HierarchicalGridCursor}) = true`
unconditionally. GeometryOps consumes it at `dual_depth_first_search.jl:31-34`, allocating a
`Vector` of child extents per internal node -- **15.84% = 9.7 s inclusive** in the P3 profile.
Making it type-parametric on the system (`... where {S} = extent_is_expensive(S)`) is the
trait worth having. **Do not add `has_analytic_extent`** -- it would duplicate an existing
trait that already has a consumer and already folds at compile time.

Per-system extent survey (measured): **S2 is the only genuinely O(1) extent in the package**
(243 ns / 0 B). HEALPix (2254 ns / 848 B) and ISEA4R (6158 ns / 848 B) own `node_extent`
methods that are *more expensive than their own `cell_cap`*, because `CAP_EDGE_SEGMENTS = 8`
means 32 perimeter points and 64 `spherical_distance` calls. IGeo7, H3 and A5 use the
boundary fallback. CopernicusDEM's system method still calls `cell_boundary`
(`systems/CopernicusDEM/system.jl:271`) although its own `BlockCursor` shows the box-only form.

**Marginal payoff, stated honestly.** Analytic extent standalone **-10.3 s of 61.17 s
(1.20x)**, and 2.97 GB of churn to zero. The tile boundary cache (P3 R1) standalone -17.5 s.
They overlap on the same `cell_boundary <- cell_cap` path, so combined is **-18.1 s ->
43.1 s (1.42x)**, and **the analytic extent's marginal value once the cache exists is
-0.6 s, ~1%**. If only one ships, ship the cache. The analytic extent still earns its place
on four grounds wall-clock does not capture: zero allocation vs a 47-528 MiB cache; it works
where the cache is disengaged (already 10 of 22 tiles, and hit rate degrades with problem
size while a closed form does not); it unlocks `node_extent_is_expensive = false`, which a
cached boundary cannot; and wrappers inherit it for free.

### D3 -- `_chunklevel` compares a local cell count against global level counts

`src/regridding.jl:70-77`:
```julia
score = abs(log(n / _levelcells(sys, a)) - log(target))
```
`n = ncells(grid)` is the *subset's* count; `_levelcells(sys, a)` is the *global*
level count. Measured: `chunkcells` of 4096, 65 536 and 524 288 **all** give
`chunklevel 0, 1 chunk` on a 21 659-cell covering, and the full-tile destination
comes out as `DGGSpace(16471968 cells, 1 chunk at level 3)` against a
`DEFAULT_CHUNK_CELLS` of 4096 -- a **4000x miss**. The whole destination becomes one
work unit: no chunk parallelism, no memory bounding, one dual-tree walk over
everything.

Explicit `chunklevel` is not a workaround: `_chunkwindows` (`src/regridding.jl:89`)
loops `for j in 1:ncells(ancestors)` over a **complete global level**. Measured on
the 21 659-cell covering: level 5 (1.7e5 cells) 0.23 s -> 2 chunks, level 6 (1.2e6)
1.6 s -> 2 chunks, level 7 (8.2e6) 11.3 s -> 3 chunks. Level 9 would scan 4.0e8.

### D4 -- a whole-globe `DGGSpace` is not constructible

Same `_chunkwindows` loop. For a complete IGEO7 level-13 destination `_chunklevel`
picks ancestor level 9 -> **403 536 072 chunks**: 72 B/chunk of stored state
(`Z7Cell` 8 + `UnitRange{Int}` 16 + `starts` 8 + `SphericalCap{Float64}` 40) =
**29-30 GB resident before any regridding**, ~210 GB of allocation churn, and
~1 000 s single-threaded just to build the space (extrapolated from 1.72 us/chunk at
level 7 and 3.46 us/chunk at level 8). At level 19 the same code needs 3.4 PB.

Compounding it, `DGGChunkTree` is **flat**: `src/regridding.jl:196-201` declares
`STI.isleaf(::DGGChunkTree) = true` and hands `zip(eachindex(ranges), caps)` to the
search, so chunk-pair discovery is a linear scan over all chunk caps -- O(nsrc_chunks
x ndst_chunks).

### D5 -- `bucket_size = 0` gives one cell per tree leaf

`src/fallbacks/cursor.jl:57` reads `grid.bucket_size`; `PartialGrid`'s default is
**0** (`partial_grid.jl:71`) and `treeify` passes no override, so `STI.isleaf`
(`cursor.jl:194-199`) never buckets and the cursor descends all 13 levels to a
single-cell leaf. `STI.node_extent` is **20.26% = 12.4 s** of the 1024 x 1024 read,
12.36 points of it inside GeometryOps' `_child_extents`. `query` already buckets at
`QUERY_BUCKET_SIZE = 16` (`query.jl:9,563`); `GR.celltree` (`src/regridding.jl:139`)
passes nothing. The source-side `BlockCursor` uses `LEAF_CELLS = 9`.

**Caveat that must not be lost:** bucketing trades interior-node extents for a larger
leaf x leaf candidate cross-product, and clipping is 20% of the read. If the extent
redundancy in D2 is fixed first, extents become nearly free and bucketing may be a net
loss. This needs an A/B at `bucket_size in {0, 7, 16, 49}`, not a blind change.

### D6 -- `NearestCell` with a DGG source does not throw, it hangs

`GR.cellat(::DGGSpace, p)` (`src/regridding.jl:133-137`) reaches
`src/fallbacks/locate.jl:14`, which calls `treeify` and runs a spatial-tree descent
**per destination centroid**. Measured: a 16x16 window (256 px -> 462 destination
cells) spent **196.5 s** in `plan_regrid` -- 0.42 s per cell, 314 M allocations.

### D7 -- `cell_area` re-derives geometry that has a closed form (not on this path)

`src/fallbacks/geometry.jl:194-195` -- IGEO7 defines no `cell_area`, so it is
`GO.area(Spherical, cell_polygon(...))`: **1 560.0 ns, 20 allocations**.
`z7grid.jl:265-269` answers in **4.5 ns** (347x) with zero allocations. **But it takes
zero samples in every P3 profile** -- `Conservative` never calls it. Worth doing for
`Extensive`, for `GR.cellarea`, and for `arealevel`; worth nothing for the measured
workload. Prioritise accordingly.

Not a drop-in substitute: measured `cell_area/equal_area_steradians - 1` converges to a
level-invariant **-11.76% at the twelve pentagons and +2.33% at their neighbours**,
because chart edges near the icosahedral vertex defects are never densified. Total area
still closes to -4.0e-15. `system.jl:322-326` claims this gap narrows with level; it does
not. Related contract gap: `src/interface/grid.jl:64-66` requires `cell_boundary` to
densify chart-curved edges; `system.jl:281-283` acknowledges IGEO7 does not.

### D8 -- traversal primitives duplicate work

* `STI.nchild` (`cursor.jl:201-208`) and `STI.getchild` (`cursor.jl:212-213`) each
  loop `_child_window` over all 7 children -- measured 4.277 us and 4.451 us for the
  same windows. Each `_child_window` probe is two binary searches, each costing a
  `leafposition` (nested `searchsortedfirst`) plus a Z7 decode, ~51 ns per access.
* `STI.getchild(cursor, i)` (`cursor.jl:220-230`) walks the iterator counting --
  O(i), quadratic for an indexed consumer.
* CopernicusDEM `node_extent` (`systems/CopernicusDEM/system.jl:270-273`) decodes the
  cell **3x** per call; `cursor.jl:258-263` has no `sizehint!` (1408 B allocated for
  240 B of payload) and redoes `tilebase`/`ncols` per pixel. Measured full extent walk
  of one tile: **5.21 s, 5.55 M allocations, 1.214 GiB**, extents only.
* `GeometryOps` `_sh_clip_to_edge_spherical!` recomputes
  `robust_cross_product(edge_start, edge_end)` -- with a `normalize`/`sqrt` -- on every
  one of the n x m inner iterations for a loop-invariant value.

### D9 -- other O(global) escape hatches on adjacent paths

* `cursor.jl:52` `selection = has_sorted_subtrees(sys) ? nothing : collect(1:n)` --
  7.75 TB at global L13. IGEO7 escapes; **A5 does not**.
* `cell_vector.jl:331-341` `_grid_positions` allocates `Vector{Int}(undef, ncells)` --
  reached by `PartialGrid(sys, level, ids::Vector)`, the `row_band` shape in
  `examples/copernicus_dem.jl`. 5 TB globally, plus 9.9 TB for the id vector itself.
  `PartialGrid(sys, tilecell, 1)` avoids it via `SubtreeIds`.
* `query.jl:538-550` `query(grid, Disjoint(...))` iterates `1:ncells(grid)`.
* `bands.jl:117-125` `tables()` takes a `ReentrantLock` per call for non-preregistered
  `N` -- 2.37 ns for GLO-30 vs 23.7 ns for an unregistered size, on a per-pixel path.

### D10 -- unenforced contracts

* Sutherland-Hodgman convexity/winding is prose only
  (`GeometryOps/src/methods/clipping/sutherland_hodgman.jl:9`); no check exists
  anywhere. Measured failure: a reflex clip window returns **33.3%** of the true
  area, a clockwise-wound convex window returns **0.0**, and nothing downstream can
  detect either (`should_store_result` only drops exact zeros).
* `Spilled` block storage never garbage-collects its files (`plans.jl:226-229`).
* Nothing in this repo specialises `knownempty`, so the lazy path's empty-chunk drop
  is dead code for every real workload here.

### D11 -- the destination's position -> cell id decode is 8.3% of the read

Unreported before P3. `cellindex(::PartialGrid, i)` (`partial_grid.jl:139`) is documented
O(1) but goes `CellVector` `getindex` (`cell_vector.jl:354-357`) -> `cellindex(sys, l, i)`
(`systems/IGeo7/system.jl:257`) -> `index_to_cell` (`z7grid.jl:625-659`), a 13-step
mixed-radix walk with two integer divisions per step. **51.1 ns a call, 8.34% = 5.1 s** of
the 61.17 s read, ~100 M calls; `div int.jl:301` alone is 2.42% of self time. One pass over
all 1 331 572 ids costs **0.053 s -- 96x redundancy**. Materialising `PartialGrid.ids` as a
`Vector{Z7Cell}` when built from a `CellVector` costs 10.6 MB at 1024 x 1024, 132 MB for the
whole tile.

### D12 -- cap radius calls `spherical_distance` once per vertex instead of once

The cheapest large win found. `spherical_distance` is **20.59% = 12.6 s** of the read:
**15.98 points from `_cornercap`** (`lib/GlobalRegridding/src/rastergrid.jl:610-613`) and
**2.32 from `points_cap`** (`src/fallbacks/caps.jl:67-70`). Both are the same loop,
`r = max(r, spherical_distance(centre, p))`, over 4 raster corners and 6 IGEO7 vertices.
GeometryOps defines `spherical_distance(x, y) = atan(norm(cross(x, y)), x . y)`
(`UnitSpherical/point.jl:144`) -- a cross product, a `sqrt` and an `atan2` each.

The angle is monotone decreasing in the dot product, so
`max_i angle(c, p_i) = angle(c, p_j)` for `j = argmin_i (c . p_i)`: scan the dot products
and call `spherical_distance` **once**. Bit-identical except for ties, which
`_padcap`/`CAP_SLACK` already absorb. **~-12% for about ten lines in two functions.**

### D13 -- the source leaf's cell caps are a generator, so they are re-evaluated per opposing leaf

`STI.child_indices_extents(::RasterCellTree)`
(`lib/GlobalRegridding/src/rastergrid.jl:802-803`) returns a generator, so GeometryOps'
`cie_1` binding (`dual_depth_first_search.jl:48-57`) re-evaluates `_rastercellcap` once per
opposing destination leaf. **14.42% = 8.8 s**, against at most 254 x 16 384 x 189.6 ns =
**0.79 s** if materialised (one pass over all 1 048 576 source cells is 0.184 s).
`RasterFlatTree` (`rastergrid.jl:866`) already shows the materialised form. Cost of the fix:
16 caps x 32 B per leaf, alive for one block build.

### D14 -- what CR PR #131's clipping cache actually buys, and the GeometryOps gate

Measured by isolating the cache with the vcat fix held constant in both arms (three file-copy
CR trees sharing one Manifest; a verbatim `git show a2e5778` revert allocates byte-identically
to a one-method `task_local_operator` stub, so the `_assemble_chunk` split and the second type
parameter cost exactly nothing -- the whole difference is the cache). All arms bit-identical
before any timing compared. GeometryOps confirmed non-stale in every arm (`eRyQ6`,
tree-sha1 `77dcc16f`), and the cache confirmed **engaged**: 1190 instances built and used in one
instrumented 512x512 lazy read.

**~7-10% median wall time, 31-36% less allocation churn, peak RSS unchanged.** Churn is
deterministic and repeats to the megabyte; wall time is not (the box carried load average
55-75 from sibling agents), which is why medians spread 1.01-1.24x.

The budget closes exactly: **1225 B/pair x 1 886 844 clips = 2.152 GiB**, against a measured
whole-workload churn delta of **2.145 GiB**. On time, 0.565 s in the clip loop (3.5%) plus GC
falling 2.000 -> 1.195 s (4.8%) = 8.3%, against 7.3-8.5% measured.

**Attribution** (innermost-ancestor phase split, 1 thread, 512x512 lazy):

| phase | with cache | without |
|---|---:|---:|
| destination geometry (IGEO7 synthesis) | 37.28% | 35.85% |
| unlabelled (GC, cursor bookkeeping) | 30.84% | 29.78% |
| source geometry (RasterGrid) | 11.87% | 11.22% |
| **clip** | **8.99%** | **11.79%** |
| dual-tree descent | 7.78% | 8.44% |
| area | 2.90% | 2.67% |

`Profile.Allocs`: bytes under Sutherland-Hodgman frames **0.501 GiB = 15.79%** with the cache
against **2.540 GiB = 48.71%** without -- i.e. SH scratch was **the single largest allocation
source in the workload**, larger than IGEO7 boundary synthesis. Clip census: 1 886 844 clips,
vertex histogram a **single bucket** -- `(5, 7)` closed rings, 100%, a 4-vertex source quad by
a 6-vertex IGEO7 hexagon.

**Benefit vs vertices per clip pair -- the relative benefit strictly DECREASES:**

| nv | speedup | bytes saved/pair |
|---:|---:|---:|
| 4 | **1.803x** | 1 238 |
| 8 | 1.375x | 1 328 |
| 16 | 1.201x | 4 039 |
| 32 | 1.094x | 4 491 |
| 100 | **1.028x** | 15 609 |

Clipping cost grows ~O(n*m) while scratch overhead grows only linearly. **DGG's 4x6 shape is
the *best* relative case on the whole curve, not the worst** -- the prior microbenchmark's
"the saving grows with vertex count" is true in bytes and **false in time-fraction**. A DGG
HEALPix destination emits 33-point closed rings at level 10 (densified chart edges) against
IGEO7 level 13's 7-point rings, putting it at the 1.09x end -- though clipping would then be a
far larger share of its total run.

**The cache is not defeated.** Worst configuration (lazy, 8 threads, one cache per spawned
chunk per block build): 1190 caches for 1 886 844 clips = **1586 mean / 1259 median / 79 p10**
clips per cache, **zero** unused. Buffers never exceed 9 points, so the win is the four vector
headers plus their initial `Memory` per call, not growth-realloc.

**Verdict on "why hasn't DGG seen much speedup":** the clip loop is only ~12% of wall time.
36-37% is IGEO7 destination geometry through the full inverse projection inside
`Trees.getcell`, which #131 does not touch. 1.195x on 12% is 3.5%; freed GC adds 4.8%. On a
751 GB box where churn is nearly free, a churn fix looks small however well it works -- and on
a memory-constrained machine the same 33% churn reduction would be worth considerably more.

**The GeometryOps release gate is NOT #131-specific and already binds.**
`SutherlandHodgmanCache` came from `b37a02b69` (#414, 2026-08-14), in no tagged release
(newest tag v0.1.43 = `e27b3809`). But `SpatialTreeInterface.node_extent_is_expensive` came
from `9f0d16006` (#462, also untagged) and is **already used** at
`src/fallbacks/position_tree.jl:120` and `lib/GlobalRegridding/src/conservative.jl:87,201`.
Defining a method on it against v0.1.43 is an `UndefVarError` at load, so **DiscreteGlobalGrids
cannot build against the current GeometryOps release today, with or without #131.** Other
post-v0.1.43 commits (`66bfbbe60` spherical `segmentize`, `b003d1207` RelateNG) are not on the
regrid path.

### D15 -- the memory envelope: `budget` works, eager has no control, and `Spilled` on tmpfs is a trap

Measured under **hard cgroup limits** (`systemd-run -p MemoryMax`), not just heap hints -- the
first time this pipeline has been run where OOM is actually reachable. Every prior timing in
this brief was taken with ~280 GB free, so the numbers here are the ones that transfer.

**`budget` bounds the lazy path -- verdict YES, first test at full-tile scale.** Both halves
honoured on every row of a 2^24..2^33 sweep: `storagebytes <= weightbudget` (92-99.9% filled,
never over) and `residency.peakbytes` hits `databudget` **to the byte** (12 582 912 at 2^24).
Peak RSS tracks predictably: **`peak ~= 1.6 GB floor + 1.3 x (0.25 x budget)`**, moving
1.70 -> 4.28 GB over a 512x budget range. `perf-P1.md`'s law-L3 "carried and ignored" finding
is genuinely fixed.

**The eager path has no memory control at all.** `lib/GlobalRegridding/src/api.jl:159`
`_rejectlazykeywords` throws if `budget`, `storage` or `chunks` is passed with `lazy = false`.
There is no way to bound a `DirectPlan`.

**Measured OOM boundary, full tile:**

| path | OOM-killed at | runs at |
|---|---|---|
| eager, 8 threads | 4, 5, 6, 7 GB | 8 GB |
| lazy | 1 GB | **2 GB** |

**Eager needs a 4x larger machine to buy a 7% time saving.** This reverses D-series framing
elsewhere in this brief: post-#133 eager is faster than lazy at the full tile (140.86 s vs
163.52 s, and 138.9 s vs 163.5 s under a 64 GB cap) **only with an unconstrained heap**. Under
any real limit lazy wins decisively.

**Correction -- eager's memory is NOT usefully thread-dependent.** An earlier reading of
"5.60 GB at 1 thread" was a 20 ms sampler missing the peak of a 785 s run; that run's kernel
`maxrss` is **8.03 GB**. Measured peaks are **8.03 / 9.82 / 9.78 / 9.73 / 8.65 GB at
1 / 2 / 4 / 8 / 16 threads** -- eager allocates the same ~206 GB in `plan_regrid` regardless of
thread count, and varies 2.7 GB run-to-run at fixed settings (9.73 vs 7.00 GB at 8 threads).
Confirmed by OOM under a hard 6 GB cap: **killed at 1 thread (837 s), 2 threads (479 s) and
8 threads (182 s)**. **Dropping threads does not rescue eager; it only makes the failure
slower.** Eager is off the table below ~10 GB, and there is no knob to bound it.

**`Spilled` does bound RSS -- but not on tmpfs, and that is a live usability trap.**

| spill target | 2 GB cap | outcome |
|---|---|---|
| `/tmp` (tmpfs = RAM) | 2 GB | **OOM-killed**, `oom_kill 1` |
| `/home` (btrfs SSD) | 2 GB | **ran**: 169.6 s, peak 2.03 GB, `oom_kill 0`, 12 883 reclaims |

Spill files on tmpfs are charged to the cgroup and are not reclaimable; on a real filesystem
they are page cache and the kernel reclaims them instead of killing. So spilling holds RSS down
at full-tile scale for **~2% wall time and 4.1 GB of disk** -- provided the directory is real
storage. **`Spilled(dir)` takes an explicit directory and the obvious idiom is
`Spilled(mktempdir())`, which lands in `/tmp`, which is tmpfs on most Linux systems.** A user
reaching for `Spilled` *because* they are short on memory would silently get no relief and then
be OOM-killed, with nothing in the API hinting why. Cheap fix: detect the filesystem type of
`dir` at construction and warn or refuse on tmpfs/ramfs.

**Measurement correction that matters for any RSS figure quoted elsewhere:** kernel `maxrss`
exceeds a 20 ms-sampled peak by **up to 2.8 GB on the eager path** (5.21 GB sampled vs 8.03 GB
kernel); lazy agrees within 0.25 GB. Sampling misses a transient eager allocation spike the
kernel catches. Prefer `/usr/bin/time -v` / `VmHWM` over in-process sampling for eager.

**`DGGSpace` global construction, measured not extrapolated:** 80.6 B/chunk, flat across
complete levels L6-L12, projecting **~30.3 GB at global L13** -- confirming D4's 29-30 GB
estimate to within 3%.

**Practical guidance implied:** for a memory-constrained user, lazy with `budget` set to taste,
plus `Spilled` on real storage if tight. One Copernicus tile onto IGEO7 L13 fits in **2 GB**
that way, against 8 GB for eager.

## Structural facts about the pipeline

* Plan selection is a **single boolean**: `lazy`, defaulting to `isdiskbacked(data)`
  (`api.jl:64`). No memory heuristic, no automatic escalation from `DirectPlan` to
  `ChunkedPlan` when a whole-domain block would not fit.
* `budget` (default 2^30) does three unrelated jobs: it splits 25% weights / 75% data
  (`plans.jl:81-95`), caps the parallel build wave width (`lazy.jl:464`), and sets the
  derived destination tile size when the space is not tileable (`lazy.jl:98`).
* Threading has exactly two sites: `conservative.jl:293-297` (intra-block, gated only
  on `Threads.nthreads() > 1`, not a keyword) and `lazy.jl:469-484` (inter-block wave).
  The tile loop, the apply loop and the eager slice loop are all serial.
* `BilinearPoint` is structurally invalid with a DGG source --
  `hascellchart(::DGGSpace)` is `false`, and the throw comes from `support_radius` at
  `LazyRegridArray` construction, not at read.
* Measured axis sensitivity at 512 x 512 px: source chunking dominates (1 chunk 118 s,
  16 chunks 6.8 s, 64 chunks 4.5 s); hierarchy-aligned destination tiles beat
  equal-sized runs (8 tiles 4.55 s vs 82 tiles 9.67 s) because `GR.subtree` reuses the
  grid cursor for exact chunk ranges instead of building a `CellCapTree`; storage
  variants are within noise (9.4-10.1 s); `budget` is the only axis that visibly moves
  I/O (2^20 forced 1312 `readblock!` / 82 MiB against 16 reads / 1 MiB at default).

## Extrapolations, and where they are unsafe

Full tile (12.96 M px -> 16.49 M cells): nnz **6.65e7** at the measured 5.128
nnz/pixel; lazy wall time **139-174 s at 8 threads** (measured ladder, exponent 0.86);
eager ~2600 s; peak RSS ~2.3 GiB at 2048 x 2048.

Global GLO-30 -> IGEO7 L13. Source pixel count derived from `bands.jl` banding, not
assumed: `1 296 000 x (100*3600 + 20*2400 + 20*1800 + 20*1200 + 10*720 + 10*360)` =
**620 524 800 000** px against 968 890 104 072 destination cells. nnz ~3.18e12 (a
**lower** bound -- 5.128 nnz/px was measured only at 45.5N, and the ratio rises toward
the equator where a pixel is 953.6 m2 against the same 526.4 m2 cell). The final
`SparseMatrixCSC` alone would be **51 TB**, with ~102 TB of transient `sparse()`
workspace. Index types are not the blocker (`Int64` throughout, nnz 3.18e12 is far
inside range); **memory is**. There is no single-matrix formulation of the global
problem.

Unsafe steps in that chain, explicitly: the 5.128 nnz/px constant is single-latitude;
the measured trees were built via `PartialGrid(sys, level, Vector{LevelIndex})`
(`HierarchicalGridCursor`) rather than the production `PartialGrid(sys, tile, 1)`
(`BlockCursor`), and since ~95% of cost is descent this could move the per-pixel cost
substantially either way; linearity is established over 64x and the global figure is
1.9e6x past the last measured point, so it is an order-of-magnitude statement, not a
prediction.

## TODO -- filled in when the running campaigns land

* Full mode matrix at 1024 x 1024, size ladder to 3600 x 3600, axis sweeps, thread
  sweep, cross-mode correctness, conservation at scale.
* DONE -- `regrid-notes/perf-P3.md` (591 lines). Phase split of the cold lazy 1024 x 1024
  read at one thread: **destination geometry 42.56% (26.0 s), source geometry 23.69%,
  clipping + area 19.97%, tree descent 10.58%, sparse assembly 2.55%, apply 0.45%.**
  Baseline 61.17 s / 16.659 GiB / 5.9% GC, 45.9 us per destination cell; `NearestCell` is
  0.50 us per cell, **93x cheaper**. 8 threads gives **4.32x** (14.17 s) and 57% of that
  profile is idle workers, not work. **No type instability anywhere in either package**
  (13 `@inferred` probes pass, no package frame carries the `runtime_dispatch` bit); the
  only dispatch site is GeometryOps `area.jl:285` at 2.0% self.
  P3's own combined estimate: fixing D2 + D11 + D13 + D12's residual takes the read
  **61.17 s -> ~31 s at one thread (2.0x)** and cuts allocated bytes by well over half.
* DONE -- see D1 above. Results in `bench/results/crfix-*.ndjson` (10 files, 82 rows).

**Pinning caveat that affects every number in this brief.** DiscreteGlobalGrids pins CR at
`rev = "agent/use-geometryops-main-cache"` in three files (root `Project.toml`,
`lib/GlobalRegridding/Project.toml`, `docs/Project.toml`). That branch is **open PR #131**
("Use GeometryOps' SutherlandHodgmanCache in spherical intersection assembly"), currently
**ahead 2 / behind 1** of `main`: it carries the clipping cache (`a2e5778`) which `main` does
not, and lacks the vcat fix (`55820f1`) which `main` has. So the repo cannot simply repin to
`main` without losing the cache. Every measurement in this brief was taken on the pinned,
unfixed tree unless explicitly marked as the patched arm. Resolution is either merging #131,
or merging `main` into that branch.
