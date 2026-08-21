# P3 — performance checkpoint: IGEO7 level-13 destination (report only, no code changed)

Method: `regrid-notes/profile-performance.md` / the `profile-performance` skill —
one long-lived `julia` process per pass (no MCP server on this box), `@timed`
one-shot baselines for a 60-second target, `Profile`/`FlameGraphs` traversed
exhaustively with the `runtime_dispatch`/`gc_event` status bits,
`Profile.Allocs`, `@inferred`/`@code_warntype` behind a function barrier, plus
an instrumented run that **counts** every destination-geometry call. Julia
1.12.6, Linux x86-64 (64 cores, two sibling agents sharing the box), repo `main`
at **79e1db0**, `--project=bench`, ConservativeRegridding
`agent/use-geometryops-main-cache` (v0.2.8, `sLk1Y`), GeometryOps `eRyQ6`.

Every CPU number below is from `-t 1`: single-thread profiles attribute cleanly
and the threaded picture is measured separately in **F**.

Judged against `perf-P1.md` (raster→raster) and `perf-P2.md` (IGEO7 level 7).
Neither's attribution survives unchanged; where I disagree, the numbers are in
**H**.

---

## Headline

1. **The destination is now the expensive side.** 42.6 % of a cold lazy read is
   IGEO7 cell geometry, against 23.7 % for the raster source, 20.0 % for the
   clipping the whole exercise exists to do, and 10.6 % for tree descent. P1's
   raster→raster build was 67 % *source* cap synthesis; P1's edge tables landed
   (`cosd`/`sind` appear nowhere in this profile; `_rastercellcap` is 189.6 ns
   and allocation-free) and the DGGS took its place.
2. **Redundancy measured by counting, not inferred: 14 206 931 `cell_boundary`
   calls for 1 331 572 destination cells — 10.67× — costing 19.7 s of a 61.2 s
   read against 2.20 s for one pass over every cell (8.96× in time).**
   The prior audit's 43× and its "55 % via `cells_cap`" do not reproduce:
   `cells_cap` is **0.40 % of the calls and 0.14 % of the wall clock**. The real
   split is `cell_cap` 50.2 % of calls / `cell_polygon` 45.1 % / the system's
   default `node_extent` 4.3 %.
3. **A second, unreported cost: 8.3 % of the wall clock (5.1 s) is the
   destination's *position → cell id* decode.** `cellindex(::PartialGrid, i)` is
   documented O(1) but goes through a `CellVector` to IGEO7's `index_to_cell`
   mixed-radix walk — 51.1 ns a call, so ≈100 M calls, against 0.053 s for one
   pass over all 1 331 572 ids.
4. **The eager path's 51 s at 512² is not its cost, it is a threading
   regression.** `DirectPlan` at 512² is **17.62 s / 4.64 GiB at one thread** and
   **54.16 s / 191.95 GiB at eight** — 3.1× slower and **41× more allocated**.
   That is CR's `reduce(vcat, ...)` merge, and it means the eager and lazy paths
   are level at one thread (17.62 vs 17.47 s).
5. **The cheapest large win is not a cache.** `spherical_distance` is
   **20.6 % of the read (12.6 s)**, and **18.3 of those points are the cap-radius
   loop** `r = max(r, spherical_distance(centre, p))` over the 4 raster corners
   (`rastergrid.jl:610-613`) and the 6 IGEO7 vertices (`caps.jl:67-70`). The angle is monotone in the dot product, so one
   `spherical_distance` on the `argmin`-dot vertex gives the same cap —
   ≈ −12 % for about ten lines.
6. Nothing in `DiscreteGlobalGrids` or `GlobalRegridding` is type-unstable:
   thirteen `@inferred` probes pass and no package frame carries the
   `runtime_dispatch` bit. The only dispatch site is still GeometryOps'
   `_naive_triangulated_spherical_ring_area` (2.0 % self).

---

## Setup

| | |
|---|---|
| source | GLO-30 tile `N45_00_E010_00`, north-west **1024×1024** window, 1 048 576 `Float32` px, behind `CountingDisk` chunked `(128, 128)` → **64 source chunks**, `RasterGrid` |
| destination | `DGGSpace` over the level-13 IGEO7 covering of the window footprint, **1 331 572 cells**, `chunklevel = 7` → **22 chunks** = 22 output tiles, sizes 1 277 … 117 649, median 60 526 |
| method / policy | `Conservative()`, `Weighted(0.5)` |
| plan | `ChunkedPlan` + `LazyRegridArray`, one `readblock!` over `1:1331572` |
| work | **254 chunk-pair block builds** (`loads = 64`, `hits = 190`); every source chunk read exactly once, `readbytes = 4 194 304` = the window exactly |
| driver | `bench/profiling/p3.jl`, `p3_micro.jl`, `p3_analyze.jl` (new, documented, re-runnable; nothing under `src/`, `lib/` or the rest of `bench/` touched) |

`10 of the 22 tiles exceed _TILE_CELL_CACHE_MAX = 65 536`
(`lib/GlobalRegridding/src/conservative.jl:108`), so for those tiles
`TileCells` (`:165`) keeps on-demand geometry and the destination polygons are
re-synthesized once per connected source chunk.

---

## A — baselines

`@timed` on a fresh plan each time; `@benchmark` is the wrong tool at 60 s.

| configuration | threads | wall | allocated | GC | note |
|---|---:|---:|---:|---:|---|
| **lazy 1024² (the target)** | **1** | **61.17 s** | **16.659 GiB** | **3.61 s = 5.9 %** | repeat 60.87 s; peak RSS 2.00 GiB |
| lazy 1024² | 8 | 14.17 s | 20.085 GiB | 2.91 s = 20.6 % | **4.32×**; matches the harness's 14.07 s |
| lazy 512² | 1 | 15.27 s | 3.907 GiB | 1.16 s = 7.6 % | |
| lazy 512² (`run_case`) | 1 | 17.47 s | 4.13 GiB | 1.01 s | 62 builds |
| lazy 512² (`run_case`) | 8 | 7.51 s | 6.25 GiB | 0.99 s | 2.33×; summed build CPU 17.38 s |
| **eager 512² `DirectPlan`**, chunked source | **1** | **17.62 s** | **4.64 GiB** | ~0 | nnz 1 344 280 |
| eager 512² `DirectPlan`, chunked source | 8 | **54.16 s** | **191.95 GiB** | | **3.1× slower than 1 thread** |
| eager 512² `DirectPlan`, in-memory raster | 1 | 17.59 s | 4.631 GiB | 1.18 s = 6.7 % | apply 0.97 s / 0.146 GiB |
| eager 512² `DirectPlan`, in-memory raster | 8 | 46.51 s | 191.95 GiB | | |
| **`NearestCell()` 1024² eager** | 1 | **0.660 s** | 0.187 GiB | 0.24 % | nnz 1 323 432 |

Derived: **45.9 µs per destination cell** for the lazy conservative read at one
thread; `NearestCell` is **0.50 µs per cell**, i.e. **93× cheaper**.

---

## B — CPU profile, lazy 1024², one thread

`Profile.init(; delay = 0.005, n = 60_000_000)`, one run, **12 151 root
samples** over a 64.2 s profiled read (5 % profiler overhead against the 61.17 s
baseline). Seconds below are the sample share × 61.17 s.

### B1 — phase split

Self samples, attributed to the **innermost** matching ancestor, so the
destination geometry inside the clip loop counts as destination geometry and not
as clipping.

| phase | % | s | what it is |
|---|---:|---:|---|
| **destination geometry (IGEO7)** | **42.56** | **26.0** | `cell_boundary_cartesian` → `dev_to_xyz` → `snyder_inv_xyz`, `points_cap`, `closed_ring`, `index_to_cell` |
| source geometry (`RasterGrid`) | 23.69 | 14.5 | `_rastercellcap`, `_rectcap`/`_sampledcap`, `_cornercap`, `getcell` |
| polygon clipping + area | 19.97 | 12.2 | Sutherland–Hodgman + `GO.area` |
| dual-tree descent | 10.58 | 6.5 | `dual_depth_first_search` bookkeeping, `_child_extents`, the cap predicate |
| sparse assembly | 2.55 | 1.6 | `_fillcoo!`, `sparse`, COO growth |
| apply (read, accumulate, write) | 0.45 | 0.28 | the executor is still a rounding error |
| other | 0.19 | 0.12 | |

### B2 — top self-cost frames (aggregated per `(func, file:line)`, not eyeballed)

| % self | frame |
|---:|---|
| 13.82 | `*` `float.jl:497` |
| 9.16 | `+` `float.jl:495` |
| 6.30 | `<` `float.jl:623` |
| 6.26 | `muladd` `float.jl:500` |
| 4.36 | `Array` `boot.jl:648` *(GC)* |
| 3.88 | `-` `float.jl:496` |
| **3.60** | **`index_to_cell` `int.jl:0`** |
| 3.58 | `GenericMemory` `boot.jl:588` *(GC)* |
| 2.83 | `/` `float.jl:498` |
| 2.42 | `div` `int.jl:301` |
| 2.20 | `<=` `float.jl:624` |
| 2.17 | `_setindex!` `array.jl:991` |
| 1.57 | `sqrt` `math.jl:628` |
| 1.43 | `_naive_triangulated_spherical_ring_area` `array.jl:0` *(dispatch)* |
| 1.41 | `acos` `trig.jl:711` |
| 1.27 | `atan` `trig.jl:0` |

Per function across all lines: `*` 13.59, `+` 9.09, `muladd` 6.27, `<` 6.25,
`-` 4.42, `Array` 4.19, `index_to_cell` 3.52, `GenericMemory` 3.06, `/` 3.01,
`acos` 2.44, `atan` 2.38, `div` 2.32, `sqrt` 1.81,
`dual_depth_first_search` 1.79, `snyder_inv_xyz` 1.26,
`robust_cross_product` 1.10, `_rastercellcap` 0.51, `_rectcap` 0.48.

Self cost is almost entirely generic float arithmetic and two allocation frames;
B3's inclusive totals are what say whose arithmetic it is. The one package
symbol in the top ten is `index_to_cell`.

### B3 — recursion-safe inclusive totals (topmost occurrence, share of the read)

| % | s | frame |
|---:|---:|---|
| 99.00 | 60.6 | `build_weights!` `conservative.jl` |
| 97.48 | 59.6 | `intersection_areas` (CR) |
| 59.25 | 36.2 | `get_all_candidate_pairs` → `dual_depth_first_search` |
| 37.90 | 23.2 | `assemble_sparse_matrix_coo` → `_run_and_store!` (37.73) |
| **32.24** | **19.7** | **`cell_boundary`** |
| 30.58 | 18.7 | `cell_boundary_cartesian` `z7grid.jl:205` |
| 24.54 / 24.29 | 15.0 / 14.9 | `dev_to_xyz` / `snyder_inv_xyz` |
| 20.99 | 12.8 | `_cornercap` `rastergrid.jl:600` |
| 20.62 | 12.6 | `cell_cap` `caps.jl:81` |
| 20.26 | 12.4 | `STI.node_extent` (both trees) |
| 16.04 | 9.8 | `cell_polygon` `geometry.jl:172` |
| 15.84 | 9.7 | `_child_extents` (GO's per-visit destination extent cache) |
| 15.58 | 9.5 | `_intersection_sutherland_hodgman` |
| 14.42 | 8.8 | `_rastercellcap` `rastergrid.jl:593` |
| 8.92 | 5.5 | `child_indices_extents` `cursor.jl:265` |
| **8.34** | **5.1** | **`index_to_cell` `z7grid.jl:625`** |
| 7.88 | 4.8 | `_rectcap` `rastergrid.jl:700` |
| 4.39 | 2.7 | `GO.area` of the clip |
| 2.29 | 1.4 | the cap predicate `_intersects` |
| 0.70 / 0.54 | 0.43 / 0.33 | `_fillcoo!` / `sparse` |
| 0.07 | 0.04 | `STI.isleaf` |

### B4 — who calls the destination geometry

| caller chain | % of read | s |
|---|---:|---:|
| `cell_boundary` ← `cell_cap` (`caps.jl:81`) | **17.74** | 10.9 |
| ⤷ ← `STI.node_extent` (`cursor.jl:237-248`) | 10.42 | 6.4 |
| ⤷ ← `STI.child_indices_extents` (`cursor.jl:265-273`) | 8.58 | 5.2 |
| ⤷ ← `node_extent(sys, c)` (`geometry.jl:318-319`) | 1.62 | 1.0 |
| `cell_boundary` ← `cell_polygon` (`geometry.jl:172`) | **14.36** | 8.8 |
| ⤷ ← `Trees.getcell` (`cursor.jl:291-295`) ← `BlockAreaOperator` (`conservative.jl:251`) | 15.51 (incl. wrappers) | 9.5 |
| `cell_boundary` ← **`cells_cap`** (`cursor.jl:243-246`, the `STORED_UNION_CAP_LIMIT` tightening) | **0.14** | 0.09 |
| `index_to_cell` ← `cellindex` ← `CellVector` `getindex` (`cell_vector.jl:354-357`) | 6.02 | 3.7 |
| `index_to_cell` ← `cellindex` ← `_child_window` (`cursor.jl:137-149`) | 1.88 | 1.1 |
| `_cornercap` ← `_rastercellcap` ← `child_indices_extents` (`rastergrid.jl:802-803`) | 13.77 | 8.4 |
| `_cornercap` ← `_sampledcap` ← `_rectcap` ← `node_extent(::RasterCellTree)` (`rastergrid.jl:786`) | 7.19 | 4.4 |

Two structural facts fall out of this table:

- **The destination leaf's cap is computed twice per visit.** At the leaf level
  `node_extent` (`cursor.jl:241`) *is* `cell_cap(grid, id)`, and
  `child_indices_extents` (`cursor.jl:272`) recomputes the same cap for the same
  cell. 10.42 % + 8.58 % of the wall clock, for one cap.
- **The source leaf's 16 caps are a lazy generator.**
  `STI.child_indices_extents(::RasterCellTree)` (`rastergrid.jl:802-803`) returns
  a generator, so GO's `cie_1` binding re-evaluates `_rastercellcap` on every
  opposing destination leaf (`dual_depth_first_search.jl:48-57`).

### B5 — dispatch and GC

Runtime-dispatch frames, whole profile:

| % self | % total | frame |
|---:|---:|---|
| 1.43 | 1.43 | `_naive_triangulated_spherical_ring_area` GeometryOps `area.jl` (inlined `array.jl:0`) |
| 0.57 | 0.85 | `_naive_triangulated_spherical_ring_area` GeometryOps `src/methods/area.jl:285` |
| 0.04 | 3.88 | `_naive_triangulated_spherical_polygon_area` GeometryOps `area.jl:301` |

**Nothing else — no frame in `DiscreteGlobalGrids` or `GlobalRegridding` carries
the bit.** Same single upstream site P1 filed and P2 re-filed.

GC-flagged self cost: `Array` `boot.jl:648` 4.16 %, `GenericMemory`
`boot.jl:588` 2.93 %, `area.jl:285` 0.57 % — 7.7 % of samples against a `@timed`
`gctime` of 5.9 %.

---

## C — the redundancy, counted rather than inferred

Counting methods added from `Main` (identity wrappers around
`DGG.cell_boundary`, `Fallbacks.cell_cap`, `Fallbacks.cells_cap`,
`DGG.cell_polygon`, `DGG.node_extent`; `bench/profiling/p3_micro.jl counts`).
Instrumented wall time 64.42 s against 61.17 s uninstrumented — 5 % overhead, so
the counts are of the real run.

| window | dst cells | tiles | src chunks | builds | `cell_boundary` | per cell | `cell_polygon` | per cell |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 256² | 84 613 | 4 | 4 | — | 571 129 | **6.75** | 84 613 | **1.00** |
| 512² | 334 521 | 8 | 16 | 62 | 3 205 559 | **9.58** | 1 257 002 | **3.76** |
| **1024²** | **1 331 572** | **22** | **64** | **254** | **14 206 931** | **10.67** | **6 407 973** | **4.81** |

`cell_boundary` calls at 1024², by caller:

| calls | share | caller |
|---:|---:|---|
| 7 133 851 | **50.2 %** | `cell_cap(::PartialGrid, …)` — the cursor's node and leaf extents |
| 6 407 973 | **45.1 %** | `cell_polygon` — the clip loop's destination polygons |
| 608 683 | 4.3 % | `node_extent(sys, c)` on interior nodes (`geometry.jl:318-319`) |
| ≈56 424 over 3 764 calls | **0.40 %** | `cells_cap` (`cursor.jl:243-246`) |

The `cell_polygon` column is the `_TILE_CELL_CACHE_MAX` cliff made visible: at
256² every tile fits the cache and the polygon is built exactly once per cell; at
1024², where 10 of 22 tiles are over the cap, it is built 4.81 times per cell.

### The floor

One pass over every destination cell (`p3_micro.jl ideal 1024`, `-t 1`):

| pass over 1 331 572 cells | time | allocated | per cell |
|---|---:|---:|---:|
| `cellindex` | **0.053 s** | 0.0 MiB | 40 ns |
| `cell_boundary` | **2.204 s** | 528.3 MiB | 1 655 ns |
| `cell_cap` | 2.067 s | 528.3 MiB | 1 552 ns |
| `cell_polygon` | 1.903 s | 894.0 MiB | 1 429 ns |

So `cell_boundary` costs **19.7 s where 2.20 s of work exists — 8.96×** — and
the id decode costs **5.1 s where 0.053 s exists — 96×**.

Source side, same treatment (`p3_micro.jl srcmicro 1024`): `_rastercellcap` is
**189.6 ns, 0 allocations** (P1's edge tables are in place —
`RasterGrid.tables isa LonLatEdgeTables`), `_rectcap` 193.3 ns,
`getcell(::RasterGrid, 1)` 61.9 ns / 240 B. One pass over all 1 048 576 source
cells costs **0.184 s**; the profile spends **8.8 s** in `_rastercellcap`, and
the 254 block builds could not need more than 254 × 16 384 × 190 ns = **0.79 s**.

### Per-call costs (BenchmarkTools minima, `$`-interpolated, level-13 cell)

| call | time | bytes | allocs |
|---|---:|---:|---:|
| `DGG.cell_boundary(grid, c)` | **1 216.1 ns** | 416 | 4 |
| `IGeo7.cell_boundary_cartesian(id)` | 1 200.0 ns | 208 | 2 |
| `Fallbacks.cell_cap(grid, c)` | 1 438.0 ns | 416 | 4 |
| `Fallbacks.points_cap(boundary)` | 216.7 ns | 0 | 0 |
| `cell_polygon(grid, c)` | 1 249.0 ns | 704 | 8 |
| `closed_ring(boundary)` | 18.9 ns | 224 | 2 |
| `cell_area(grid, c)` (fallback via `GO.area`) | **1 560.0 ns** | 1 232 | **20** |
| `IGeo7.cell_area(id)` (closed form, `z7grid.jl:265-269`) | **4.5 ns** | 0 | 0 |
| `node_extent(sys, c₇)` | 1 440.0 ns | 416 | 4 |
| `cells_cap(grid, 64 ids)` | 92 172.0 ns | 44 096 | 264 |
| `cellindex(grid, 1)` | **51.1 ns** | 0 | 0 |
| `STI.nchild(cursor, level 1)` | 8 040.0 ns | 0 | 0 |
| `collect(STI.getchild(cursor, level 1))` | 8 186.7 ns | 512 | 3 |
| `STI.node_extent(cursor, level 1)` | 1 236.0 ns | 352 | 4 |

---

## D — allocations

`Profile.Allocs` at `sample_rate = 0.001` on the lazy **512²** read — 512²
rather than 1024² because the profiler's per-allocation hook makes a
sixty-second run impractical, and the two have the same shape (16.659 GiB at
1024² against 3.907 GiB at 512², a ratio of 4.26 for 3.98× the cells). 52 277
sampled allocations, **3 473 MiB scaled** against a measured 3.890 GiB — the
scaling is honest.

As in P2, **Julia 1.12's allocation profiler resolves no Julia leaf frame**:
every sample's leaf is `maybe_record_alloc_to_profile gc-alloc-profiler.h:43`
and its caller is `jl_gc_alloc_`/`ijl_gc_small_alloc`. The type histogram is what
names the owner.

| MiB | allocs | type | owner |
|---:|---:|---|---|
| **1 583.8** | 11.18 M | `Memory{UnitSphericalPoint{Float64}}` | `cell_boundary` (`system.jl:285-293`), `closed_ring` (`geometry.jl:11-24`), SH clip |
| **491.8** | 3.22 M | `Memory{Tuple{Float64,Float64,Float64}}` | **`cell_boundary_cartesian`** `z7grid.jl:211` — IGEO7 destination geometry, nothing else allocates this type |
| 361.2 | 11.84 M | `Vector{UnitSphericalPoint{Float64}}` | as above |
| 248.7 | 0.005 M | `Profile.Allocs.BufferType` | the profiler itself |
| 153.7 + 134.8 | 9.92 M | `Vector`/`Memory{LinearRing{…}}` | `cell_polygon` `geometry.jl:172-173` |
| 101.0 | 6.62 M | `Float64` (boxed) | |
| 98.8 | 3.24 M | `Vector{Tuple{Float64,Float64,Float64}}` | `cell_boundary_cartesian` |
| 64.9 + 10.2 | | `Memory`/`Vector{SphericalCap{Float64}}` | GO `_child_extents` |
| 59.4 + 58.0 | | `Memory`/`Vector{Vector{UnitSphericalPoint}}` | SH clip |
| 47.1 + 25.6 | | `Memory`/`Vector{Tuple{Int,SphericalCap}}` | `child_indices_extents` `cursor.jl:269` |
| 31.5 | 0.159 M | `Generator{SmallVector{7,Z7Cell}, …}` | `STI.getchild` `cursor.jl:211` |

Cross-check: the counted 3 205 559 `cell_boundary` calls at 512² × the measured
208 B of `cell_boundary_cartesian` = **636 MiB**, against 590.6 MiB of
tuple-typed bytes in the histogram. Destination boundary synthesis is
**≈ 1.2 GiB of the 3.9 GiB** a 512² read allocates, and the same cache that
removes the time removes the bytes.

---

## E — the other two configurations

### E1 — eager `DirectPlan`, 512², one thread (17.59 s, 3 497 samples)

**Same shape, destination share higher**, because a single whole-domain source
chunk amortizes the raster side that the lazy path pays 254 times.

| phase | eager 512² | lazy 1024² |
|---|---:|---:|
| destination geometry | **55.30 %** | 42.56 % |
| clipping + area | 17.56 % | 19.97 % |
| source geometry | 14.53 % | 23.69 % |
| descent | 11.44 % | 10.58 % |
| sparse assembly | 1.09 % | 2.55 % |

`cell_boundary` 38.26 % (via `cell_cap` 22.51, `cell_polygon` 15.36,
**`cells_cap` 0.40**), `index_to_cell` 12.55 %, `node_extent` 21.30 %
(20.65 of it inside `_child_extents`), `_rastercellcap` 12.87 %,
SH clip 13.64 %. Apply is 0.97 s against a 17.59 s build — 18×, and 28
allocations' worth of executor as in P1/P2.

### E2 — `NearestCell()` on a raster source, 1024² (0.660 s, 143 samples)

A different path end to end, and **93× cheaper at the same size** (0.660 s
against 61.17 s; 0.50 µs versus 45.9 µs per destination cell).

| phase | % |
|---|---:|
| destination geometry | 48.25 |
| source geometry | 23.78 |
| sparse assembly | 21.68 |
| other | 6.29 |

There is no descent and no clipping: the cost is
`cellcentroid(::DGGSpace, i)` (`src/regridding.jl:130-131`) →
`cell_centroid` → `_cell_center_xyz` (`z7grid.jl:178`) → `dev_to_xyz` →
`snyder_inv_xyz` (**30.77 % of the run**), plus `cellindex`/`index_to_cell`
(9.79 %), the raster's `searchsortedfirst`/`_edgeindex` locate (12.6 %), and
`sparse` (6.99 %). **The same inverse-Snyder chart dominates, at one point per
cell instead of six** — which is the cleanest possible statement of where IGEO7's
cost lives.

---

## F — threads

| | 1 thread | 8 threads | ratio |
|---|---:|---:|---:|
| **lazy 1024²** | 61.17 s / 16.66 GiB / GC 5.9 % | **14.17 s** / 20.09 GiB / GC 20.6 % | **4.32×** |
| lazy 512² | 17.47 s / 4.13 GiB | 7.51 s / 6.25 GiB | 2.33× |
| **eager 512² `DirectPlan`** | **17.62 s / 4.64 GiB** | **54.16 s / 191.95 GiB** | **0.33× — slower** |
| eager 512², in-memory raster | 17.59 s / 4.63 GiB | 46.51 s / 191.95 GiB | 0.38× |

The 8-thread lazy profile (17 166 samples) is **57.45 % outside any regrid
phase** — `poptask` `task.jl:1216` 19.49 % plus 43.77 % unresolved `:-1`, i.e.
idle workers — while the absolute sample count of real work is unchanged
(destination geometry 5 177 samples at 8 threads against 5 171 at 1). The wave
is not the bottleneck; the threads are idle waiting on it. GC-flagged self cost
rises from 7.7 % to ≈27 % of samples and `@timed gctime` from 5.9 % to 20.6 %.

The eager row is CR's threaded COO merge (`intersection_areas.jl:218-220`,
`reduce(vcat, …)` over `npartitions = 8 × 4` partitions): **41× the allocated
bytes for the same 1 344 280 nonzeros.** This is the item a sibling agent is
fixing; the number to beat is 54.16 s → 17.62 s at 512².

---

## G — inference

Thirteen `@inferred` probes behind a function barrier (concrete arguments, never
a closure over an untyped global — the naive form reports spurious failures):

```
OK cell_boundary(grid,c)   OK cell_cap(grid,c)      OK cell_polygon(grid,c)
OK cell_area(grid,c)       OK node_extent(sys,c7)   OK cells_cap(grid,ids)
OK cellindex(grid,1)       OK STI.node_extent       OK STI.nchild
OK STI.isleaf              OK STI.child_indices_extents(leaf)
OK Trees.getcell(cursor,1) OK GO.area(Spherical, poly)
```

`@code_warntype` on the two hottest interface methods: `STI.node_extent(cursor)`
is `Body::SphericalCap{Float64}`, `GO.area(Spherical, poly)` is `Body::Float64`.

The one instability worth recording is **`STI.getchild(::WindowCursor)`
(`src/fallbacks/cursor.jl:211-212`), whose return type is a 3-way `Union`** of
`Iterators.Filter{…Generator{Vector{Z7Cell}…}}`,
`…{SmallVector{12,Z7Cell}}` and `…{SmallVector{7,Z7Cell}}` — `_child_ids`
(`cursor.jl:122-126`) returns an empty `Vector`, `rootcells` or `children`
depending on the branch. The call sites union-split it (no frame carries the
dispatch bit and the generator is only 31.5 MiB of allocation), so this is a
tidiness item, not a cost.

---

## H — verdict on the six carried-forward claims

| # | claim | verdict |
|---|---|---|
| **1** | IGEO7 destination geometry re-derived ~43× per cell; 55 % via `cells_cap` (`cursor.jl:243-246`), 23 % `cell_polygon`, 22 % `cell_cap`/`node_extent` | **Confirmed in kind, refuted in magnitude and attribution.** Counted: **10.67×** at 1024² (6.75× at 256², 9.58× at 512²), 8.96× measured in time (19.7 s against a 2.20 s floor). Split by calls: `cell_cap` **50.2 %**, `cell_polygon` **45.1 %**, `node_extent(sys,·)` 4.3 %, **`cells_cap` 0.40 %** (0.14 % of wall). The `STORED_UNION_CAP_LIMIT` path is not the problem. Root cause confirmed: IGEO7 defines no `node_extent`, so `src/fallbacks/geometry.jl:318-319` falls back to `inflate_cap(cell_cap(…))`. Per-call cost confirmed to the digit: **1 216.1 ns, 416 B, 4 allocs**. |
| **2** | `PartialGrid` `bucket_size = 0` ⇒ one cell per leaf; `query` uses 16, the regrid path sets nothing | **Confirmed**, structurally and by inspection: default 0 at `src/fallbacks/partial_grid.jl:66,71`, gate at `cursor.jl:194-199`, measured `treeify(grid).bucket_size == 0` with `leaf_level = 13`, `root level = -1`; `QUERY_BUCKET_SIZE = 16` at `src/fallbacks/query.jl:9,563`; `GR.celltree(::DGGSpace)` is a bare `treeify(space.grid)` (`src/regridding.jl:139`). The "19.2 M nodes" figure was not re-derived. |
| **3** | `cell_area` goes through `cell_polygon` + `GO.area` although a closed form exists | **Confirmed per call, refined to worthless on this path.** 1 560.0 ns / 1 232 B / 20 allocs against **4.5 ns / 0 B** — **347×**. But `cell_area` **does not appear anywhere in any of the three profiles**: `Conservative` takes `GO.area` of the *clip*, never of a whole cell. Fix it for `Extensive`/`cellarea` callers, not for this workload. |
| **4** | `STI.nchild` and `STI.getchild` each recompute `_child_window` | **Confirmed per call, refuted as a cost here.** Measured 8 040.0 ns and 8 186.7 ns at a level-1 cursor (about 2× the quoted 4.28/4.45 µs on this machine) and both do walk `_child_window`. But **`nchild` never appears in the profile** — GO's `dual_depth_first_search` (`dual_depth_first_search.jl:45-84`) calls only `isleaf`, `getchild`, `node_extent` and `child_indices_extents`. `isleaf` is 0.07 %; the `_child_window` cost that *is* real reaches the profile as `cellindex ← _child_window` at 1.88 %. |
| **5** | CR's threaded merge is quadratic; expect GC/alloc noise under threads | **Confirmed and quantified.** Eager 512²: **17.62 s / 4.64 GiB at 1 thread vs 54.16 s / 191.95 GiB at 8** for a bit-identical 1 344 280-nonzero block. The lazy path is nearly immune (4.13 → 6.25 GiB at 512²) because its blocks are small. Profiling at `-t 1` was the right call. |
| **6** | P1/P2's raster→raster attributions may not carry over | **They carry over only in part.** P1's headline (`_rastercellcap` 67 % of the build) is **fixed**: the edge tables are in place, `cosd`/`sind` appear nowhere in this profile (P1: `cosd` 10.42 % self), `_rastercellcap` is 189.6 ns and allocation-free, and the source side is down to 23.7 % — what remains of it is `spherical_distance`, not the chart (R2). P2's `cell_boundary` 42.6 % of a cold level-7 read is **32.2 %** here at level 13, and P2's ≈25× redundancy is **10.67×** — the DGGS destination is still the dominant cost, but less dominant and for a different reason (four fifths of it is now `cell_cap` + `cell_polygon`, not the union-cap tightening). P2's Q2 (the lazy path barely threads) has improved from 1.6× to **4.32×**. |

---

## Prioritized proposals — biggest win first, nothing applied

Wins are against the **61.17 s** one-thread lazy 1024² read. R1–R4 overlap where
noted; the combined figure at the end accounts for it.

### R1 — one destination-cell boundary cache for the life of a `subtree`/tile

- **Evidence**: 14 206 931 `cell_boundary` calls for 1 331 572 cells (**10.67×**);
  32.24 % of the read = **19.7 s** against **2.204 s** for one pass over every
  cell (**8.96×**); ≈1.2 GiB of the 3.9 GiB a 512² read allocates is
  `Memory{NTuple{3,Float64}}` + `Memory{UnitSphericalPoint}` from that synthesis;
  50.2 % of the calls come from `cell_cap` and 45.1 % from `cell_polygon`, so
  **one cache serves both**.
- **Where**: `src/fallbacks/cursor.jl:237-248` (`node_extent`), `:265-273`
  (`child_indices_extents`), `:291-295` (`Trees.getcell`) — all three key on a
  grid position and all three end in `cell_boundary(cursor.grid, …)`. A
  `Vector{Vector{USPoint}}` (or a flat 6·n `Float64` buffer) sized to the
  cursor's window, built lazily and dying with the cursor, costs
  **528 MiB for the whole 1 331 572-cell destination** and far less per tile
  (117 649 cells × 416 B = 47 MiB for the largest).
- **Sub-item, free**: at the leaf level `node_extent` (`cursor.jl:241`) *is*
  `cell_cap(grid, id)` and `child_indices_extents` (`cursor.jl:272`) recomputes
  the same cap for the same cell — 10.42 % + 8.58 % of the read, for one cap.
  Even without a cache, letting the leaf's `child_indices_extents` reuse the
  extent the caller already holds removes one of the two.
- **Estimated win**: **−17.5 s, −28.6 % ⇒ 1.40×** at one thread, plus most of
  the 16.66 GiB and the GC (worth 20.6 % at eight threads, not 5.9 %).

### R2 — one `spherical_distance` per cap instead of four (or six)

The smallest change on this list and the second-largest win.

- **Evidence**: `spherical_distance` is **20.59 % = 12.6 s** of the read
  (the remaining 2.25 % is the descent's cap predicate, which this does not
  touch) — **15.98 % from `_cornercap`** (`lib/GlobalRegridding/src/rastergrid.jl:610-613`)
  and **2.32 % from `points_cap`** (`src/fallbacks/caps.jl:67-70`). Both are the
  same loop, `r = max(r, spherical_distance(centre, p))` over 4 raster corners
  and over 6 IGEO7 vertices. GeometryOps defines
  `spherical_distance(x, y) = atan(norm(cross(x, y)), x ⋅ y)`
  (`UnitSpherical/point.jl:144`), so each call is a cross product, a `sqrt` and
  an `atan2` — and the profile agrees: `atan` 2.38 %, `sqrt` 1.81 %, plus the
  float ops feeding them.
- **The fix**: the angle is monotone decreasing in the dot product, so
  `max_i angle(c, p_i) = angle(c, p_j)` where `j = argmin_i (c ⋅ p_i)`. Scan the
  dot products (3 multiplies, 2 adds each) and call `spherical_distance` **once**.
  Bit-for-bit the same cap except for ties, which `_padcap`/`CAP_SLACK` already
  absorb.
- **Estimated win, standalone**: 3/4 of the `_cornercap` cost and 5/6 of the
  `points_cap` cost ⇒ **≈ −7.5 s, −12 %**. After R1 and R3 have removed the
  repeated calls, its residual is the node-extent path `_rectcap` → `_sampledcap`
  → `_cornercap`, still **≈ −2.5 s**.

### R3 — materialize the source leaf's cell caps

- **Evidence**: `STI.child_indices_extents(::RasterCellTree)`
  (`lib/GlobalRegridding/src/rastergrid.jl:802-803`) is a **generator**, so GO's
  `cie_1` binding (`dual_depth_first_search.jl:48-57`) re-evaluates
  `_rastercellcap` for each opposing destination leaf. Cost **14.42 % = 8.8 s**;
  the 254 block builds need at most 254 × 16 384 × 189.6 ns = **0.79 s**
  (one pass over all 1 048 576 source cells is 0.184 s). `RasterFlatTree`
  (`rastergrid.jl:866`) already shows the materialized form.
- **Estimated win**: **−8.0 s, −13.1 %**. Cost: 16 caps × 32 B per leaf, alive
  for one block build.

### R4 — stop decoding the destination's cell ids

- **Evidence**: `index_to_cell` (`src/systems/IGeo7/z7grid.jl:625-659`, a 13-step
  mixed-radix walk with two integer divisions per step) is **8.34 % = 5.1 s**,
  reached through `cellindex(::PartialGrid, i)` (`partial_grid.jl:139`) →
  `CellVector` `getindex` (`cell_vector.jl:354-357`) → `cellindex(sys, l, i)`
  (`src/systems/IGeo7/system.jl:257`). 51.1 ns a call; **one pass over all
  1 331 572 ids costs 0.053 s — 96× redundancy**. `div int.jl:301` alone is
  2.42 % of self time.
- **Where**: materialize `PartialGrid.ids` as a plain `Vector{Z7Cell}` when it is
  constructed from a `CellVector` (10.6 MB at 1024², 132 MB for the whole tile),
  or memoize per cursor window alongside R1's boundary cache.
- **Estimated win**: **−5.0 s, −8.2 %**; overlaps R1 wherever a cached boundary
  makes the decode unnecessary, so count R1+R4 as ≈ −20 s rather than −22.6 s.

### R5 — make the tile geometry cache byte-bounded, not cell-bounded

- **Evidence**: `_TILE_CELL_CACHE_MAX = 65 536`
  (`lib/GlobalRegridding/src/conservative.jl:108`, applied at `:165`) disengages
  for **10 of the 22 tiles** at 1024², and `cell_polygon` calls per cell go
  **1.00 → 3.76 → 4.81** across 256²/512²/1024² exactly as tiles cross the cap.
  The cache stores polygons (704 B/cell); caching **boundaries** (416 B/cell, and
  what R1 needs anyway) halves the cost of the thing the cap exists to bound.
- **Estimated win**: subsumed by R1 if R1 caches at the cursor; standalone it is
  the `cell_polygon` half, ≈ −8.8 s. Either way, a byte budget (the plan already
  carries one) is the right knob, not a cell count.

### R6 — CR's threaded COO merge (upstream, already owned)

- **Evidence**: eager 512² **54.16 s / 191.95 GiB at 8 threads vs 17.62 s /
  4.64 GiB at 1** — 3.1× slower, 41× the bytes, for a bit-identical block.
  `reduce(vcat, getindex.(all_results, k); init = …)` at
  `ConservativeRegridding/src/regridder/intersection_areas.jl:218-220`.
- **Estimated win**: the eager path back to ≈17.6 s at 512² (**3.1×**); no effect
  on the lazy path, whose blocks are too small to trip it.

### R7 — bucket the destination cursor (measure before adopting)

- **Evidence**: `bucket_size = 0` means the cursor descends all 13 levels to one
  cell per leaf; `STI.node_extent` is **20.26 % = 12.4 s**, 12.36 % of it inside
  GO's `_child_extents`. `query` already uses `QUERY_BUCKET_SIZE = 16`
  (`src/fallbacks/query.jl:9,563`) and `GR.celltree` (`src/regridding.jl:139`)
  passes nothing.
- **Caveat, which is why this is R7 and not R1**: bucketing trades interior-node
  extents for a larger leaf×leaf candidate cross-product, and the clip is 20 % of
  the read. With R1 in place the extents are nearly free and this may become a
  net loss. Worth one A/B at `bucket_size ∈ {0, 7, 16, 49}`, not a blind change.

### R8 — upstream `_naive_triangulated_spherical_ring_area`

Still the only runtime-dispatch site anywhere
(`GeometryOps/src/methods/area.jl:285`, 2.0 % self / 3.88 % total through
`_naive_triangulated_spherical_polygon_area` at `:301`). Filed by P1, re-filed by
P2, unchanged. ≈ −2 %.

### R9 — `cell_area`'s closed form (not this workload)

347× per call (1 560.0 ns → 4.5 ns) and 20 allocations → 0, at
`src/fallbacks/geometry.jl:194-195` versus `src/systems/IGeo7/z7grid.jl:265-269`.
**Zero samples in every profile here** — worth doing for `Extensive` and
`cellarea` users, worth nothing for `Conservative`.

**Combined**: R1 + R4 (≈ −20 s) + R3 (−8.0 s) + R2's residual (−2.5 s) takes the
read **61.17 s → ≈ 31 s at one thread (2.0×)**, cuts the allocated bytes by well
over half, and should take the 8-thread GC share down with them — which is where
the remaining 8-thread headroom is, since 57 % of that profile is idle workers,
not work. R2 alone, at perhaps ten lines in two functions, is −12 %.

---

## Reproduction

```
julia -t 1 --project=bench bench/profiling/p3.jl lazy1024      # baseline + CPU profile
julia -t 1 --project=bench bench/profiling/p3.jl eager512
julia -t 1 --project=bench bench/profiling/p3.jl nearest1024
julia -t 1 --project=bench bench/profiling/p3.jl allocs 512 0.001
julia -t 1 --project=bench bench/profiling/p3.jl harness 512   # and -t 8
julia -t 8 --project=bench bench/profiling/p3.jl lazy1024t8
julia -t 1 --project=bench bench/profiling/p3_micro.jl micro|counts 1024|infer|ideal 1024|srcmicro 1024
julia    --project=bench bench/profiling/p3_analyze.jl <dump.jls>
```

`p3.jl` serializes `Profile.retrieve()` to `$P3_DUMP/prof_<pass>.jls` so the
flame graph can be re-traversed without re-running the job; `p3_analyze.jl` does
the phase split, the caller splits and the self-cost tables offline.
`p3_micro.jl counts` adds counting methods to `DGG.cell_boundary`,
`Fallbacks.cell_cap`, `Fallbacks.cells_cap`, `DGG.cell_polygon` and
`DGG.node_extent` **from `Main`** — semantic identity wrappers, 5 % overhead,
nothing under `src/` or `lib/` edited. The only files added are the three under
`bench/profiling/` and this note.
