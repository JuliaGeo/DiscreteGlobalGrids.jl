# Where the eager `plan_regrid` churn goes

Measured 2026-08-18 to attribute the ~206 GiB of allocation churn the eager
(`lazy = false`) path spends building the full-tile IGEO7 L13 x Copernicus GLO-30
plan. Everything below is measured, not modelled; the raw rows are in
`bench/results/alloc-attribution.ndjson`.

## Configuration, stated exactly

| | |
|---|---|
| repo | `79e1db0` + the uncommitted `bench/` workspace member; no tracked file changed |
| Julia | 1.12.6, **1 thread** (`julia --project=bench -t1`) on the 64-core box, load average 25-30 |
| ConservativeRegridding | `0.2.8`, rev `agent/use-geometryops-main-cache`, tree-sha1 **`260f9199ce077283ee76bb0aa2e2e2e8d153f188`**, depot slug `08cPv` |
| GeometryOps | `0.1.43`, rev `main`, tree-sha1 **`77dcc16f35c8271b34a441c923a313c7c2995380`**, depot slug `eRyQ6` (carries `SutherlandHodgmanCache`) |
| workload | `bench/harness.jl`, GLO-30 tile `N45_00_E010_00`, north-west **1024 x 1024** window, source `:chunked` with `srcchunk = (128,128)`, destination IGEO7 **level 13** with `chunklevel = 7`, `Conservative()` / `Weighted(0.5)`, `lazy = false` |
| profiler | `Profile.Allocs.@profile sample_rate = 0.001` around **`plan_regrid` only**, after a 128 x 128 warm-up build so no compilation allocation is sampled |
| samples | 236 614 (cached arm), 318 559 (no-cache arm) |

The source array and the destination `DGGSpace` are built *outside* every timed
and profiled region, so every number below is plan construction alone.

## Step 1 -- the total, and whether 1024^2 really scales to the full tile

`plan_regrid` churn at 1 thread, four windows:

| window | plan churn | s | allocations | nnz |
|---:|---:|---:|---:|---:|
| 128^2 | 0.379 GB (0.353 GiB) | 1.2 | 4.3 M | 83 987 |
| 256^2 | 1.319 GB (1.228 GiB) | 4.3 | 16.0 M | 335 975 |
| 512^2 | 4.973 GB (4.631 GiB) | 17.3 | 62.2 M | 1 344 280 |
| **1024^2** | **18.444 GB (17.18 GiB)** | 62.4 | 236 M | 5 380 753 |

The 1024^2 figure lands on top of the existing 8-thread reference
(`bench/results/crfix-eager-fix.ndjson`: 18.580 GB), confirming the post-#133
eager path allocates the same bytes regardless of thread count -- **0.7% apart**.

Scaling to the full tile: 3600^2 / 1024^2 = **12.359x in pixels**, and the
measured full-tile eager plan churn in that same reference series is
**221.17 GB = 205.98 GiB**, i.e. **11.99x** the measured 1-thread 1024^2 figure
(fitted exponent **0.985** against pixel count). So the ~206 GiB number is not
an extrapolation artefact: linear scaling of the 1024^2 measurement predicts
227.9 GB against a measured 221.2 GB, **3% high**. Everything in this note can
be multiplied by ~12 to reach the full tile.

## Step 2 -- deterministic stage split (no sampling involved)

`bench` script `step3_stages.jl` replays exactly what
`GlobalRegridding.wholeblock` does, stage by stage, with `@timed` around each
stage. The stages sum to **18.434 GB against the 18.444 GB the whole build
churns -- 0.06% unaccounted**, so this split is exact, not sampled.

| stage | bytes | % | s |
|---|---:|---:|---:|
| `subtree` for both spaces (celltree) | 384 | 0.00 | 0.00 |
| `get_all_candidate_pairs` (dual-tree descent) | 5 067 452 872 | **27.49** | 37.2 |
| `assemble_sparse_matrix_coo` (geometry + clip + area) | 12 509 869 112 | **67.86** | 25.0 |
| `SparseArrays.sparse` inside `intersection_areas` | 199 626 024 | 1.08 | 0.12 |
| `GlobalRegridding._fillcoo!` | 446 441 280 | 2.42 | 0.23 |
| `WeightBlock` (`sparse` from the COO) | 210 275 960 | 1.14 | 0.12 |

7 554 895 candidate pairs, 5 380 753 nonzeros (1.404 pairs per nonzero).

The headline here is that **the dual-tree descent is not free**: it is 27.5% of
the plan's churn and 60% of its wall time, and it allocates *before* a single
polygon has been clipped.

## Step 3 -- category decomposition

Sampled bytes were attributed by walking each stack innermost-outward and taking
the first frame that matches a category rule (so a `GeoInterface` wrapper
allocation is charged to whichever geometry synthesised it). The sampled shares
were then re-weighted by the deterministic stage weights above; the sparse
category comes wholly from the deterministic split, because it is a few dozen
very large allocations that a count-based sampler cannot estimate.

| category | bytes | GiB | % of churn | full tile (est.) |
|---|---:|---:|---:|---:|
| **(a) destination cell geometry** | 9.75 GB | 9.08 | **52.9** | ~109 GiB |
| (d) spherical area | 2.75 GB | 2.56 | 14.9 | ~31 GiB |
| (c) Sutherland-Hodgman clip | 2.39 GB | 2.23 | 13.0 | ~27 GiB |
| (b) source pixel geometry | 1.92 GB | 1.79 | 10.4 | ~21 GiB |
| (e) sparse / COO assembly | 0.86 GB | 0.80 | 4.7 | ~10 GiB |
| (f) descent bookkeeping, extents, caps | 0.76 GB | 0.71 | 4.1 | ~9 GiB |

Cross-check: attributing the sampled bytes directly, without the stage
re-weighting, gives 52.5 / 15.2 / 13.2 / 10.6 / -- / 4.0 after allowing for the
4.65% the sampler misses. The two routes agree to within 0.5 points.

Per candidate pair: 1290 B of destination geometry, 364 B of area, 317 B of
clip output, 254 B of source geometry, 113 B of sparse assembly, 101 B of
descent -- **2440 B per pair**, against 24 B of actual weight payload.

### Top allocation sites

Percentages are of total plan churn.

| category | site | % | GiB |
|---|---|---:|---:|
| a | `src/systems/IGeo7/system.jl:287` `cell_boundary` (`Vector{USPoint}` per cell) | 19.1 | 3.28 |
| a | `src/systems/IGeo7/z7grid.jl:211` `cell_boundary_cartesian` (`Vector{NTuple{3,Float64}}`) | 18.8 | 3.22 |
| a | `src/fallbacks/geometry.jl:14` `closed_ring` (third copy of the same 6 corners) | 9.8 | 1.68 |
| a | `src/fallbacks/geometry.jl:172` `cell_polygon` (`[LinearRing]` vector + wrappers) | 2.9 | 0.50 |
| a | `src/fallbacks/cursor.jl:269` `child_indices_extents` (materialised leaf caps) | 1.9 | 0.32 |
| d | `GeometryOps src/methods/area.jl:275` `points = collect(...)` in `_naive_triangulated_spherical_ring_area` | 12.6 | 2.16 |
| d | `GeometryOps src/methods/area.jl:285` (`pop!` path / iterator state) | 2.6 | 0.44 |
| c | `sutherland_hodgman.jl:421` fresh closed result ring | 6.0 | 1.03 |
| c | `sutherland_hodgman.jl:424` `GI.Polygon([result])` | 2.1 | 0.36 |
| c | `sutherland_hodgman.jl:410` **degenerate north-pole polygon for disjoint pairs** | 2.5 | 0.44 |
| b | `lib/GlobalRegridding/src/rastergrid.jl:407` `_cellring` | 7.7 | 1.32 |
| b | `lib/GlobalRegridding/src/rastergrid.jl:376` `getcell` (`[LinearRing]` + wrappers) | 2.9 | 0.50 |
| f | `GeometryOps .../SpatialTreeInterface/dual_depth_first_search.jl:33` `_child_extents` | 2.9 | 0.50 |
| e | `GlobalRegridding src/methods.jl:83-85` via `_fillcoo!` | 2.4 | 0.42 |
| e | `WeightBlock` + CR's `sparse` | 1.8 | 0.36 |

Three observations fall straight out of that table.

**Destination geometry is synthesised twice over, in two different shapes.** Of
its 9.08 GiB, **4.01 GiB is spent inside the dual-tree descent** and 5.07 GiB
inside the clip loop. The descent figure works out to **3232 B per destination
cell** -- `STI.child_indices_extents(::HierarchicalGridCursor)`
(`src/fallbacks/cursor.jl:269`) rebuilds a leaf's cell caps, and each cap costs a
full `cell_boundary` synthesis, once per opposing source leaf that visits it.
The clip-loop figure is **720 B per candidate pair**, which is the whole
`cell_boundary_cartesian` -> `cell_boundary` -> `closed_ring` -> `cell_polygon`
chain: four heap objects carrying the same six corners, re-run from the Snyder
inverse projection for every pair the cell takes part in. At 1024^2 the tile
geometry cache is off by design -- 1 331 572 cells is 20x
`_TILE_CELL_CACHE_MAX` (65 536) -- and the whole-domain eager block is exactly
the case that cache cannot serve.

**The area computation allocates as much as the clipping does.** `area.jl:275`
copies the intersection ring the clipper just built into a *fresh* vector
(`collect(Iterators.map(...))`) before summing triangles: 2.16 GiB, 12.6% of the
whole plan, for a streaming loop over at most nine points.

**2.5% of the plan's churn is spent describing pairs that do not intersect.**
`sutherland_hodgman.jl:410` allocates a fresh `GI.Polygon([[north_pole x 3]])`
for every disjoint pair, and there are 2.17 M of them (7.55 M pairs minus 5.38 M
nonzeros). A shared constant, or a sentinel the operator recognises, removes
0.44 GiB here and the matching `area` call behind it.

## Step 4 -- what the Sutherland-Hodgman cache already saves

The cache was defeated without touching any package source, by overwriting one
method from the profiling script:

```julia
@eval CR task_local_operator(op::DefaultIntersectionOperator{<:GeometryOps.Spherical}) = op
```

which is the single method CR #131 added; with it neutralised, each
`intersection` call constructs its own `SutherlandHodgmanCache` and grows its own
four buffers, exactly as before the PR. Same tree, same manifest, same window,
same sample rate.

| arm | plan churn | s | clip category |
|---|---:|---:|---:|
| cache on (as shipped) | 18.445 GB | 65.8 | 2.23 GiB, **100% result polygons, 0 B of scratch** |
| cache defeated | 27.673 GB | 69.1 | 11.41 GiB (5.92 GiB scratch buffers, 2.25 GiB result, 3.24 GiB clip-edge buffers and per-call cache objects) |

**The cache already removes 9.228 GB -- 33.3% of what this plan would otherwise
churn, a 1.50x reduction.** Per candidate pair that is **1221 B**, which lands on
top of the 1225 B/pair the `perf-P3-brief.md` D14 study derived independently on
the 512^2 lazy path. In the shipped arm the sampler finds **zero** allocations on
any cache-buffer line: the only surviving clip allocations are lines 410, 421 and
424, all of them the returned polygon. There is no scratch left to remove; the
remaining 13.0% of clip churn is the result geometry itself, and only an
area-accumulating operator that never materialises the intersection polygon can
take it.

Scaled to the full tile by the measured 11.99x, the cache is already saving
roughly **111 GB = 103 GiB** of churn on the eager path: without it the 206 GiB
would be about **309 GiB**.

## Where the remaining 206 GiB would have to come from

Ranked by measured size, with the mechanism each one needs:

1. **109 GiB destination geometry** -- memoise or stream `cell_boundary`. The
   descent half (48 GiB) wants leaf caps computed once per leaf rather than once
   per opposing leaf; the clip half (61 GiB) wants either a per-block cell cache
   that scales past `_TILE_CELL_CACHE_MAX` or a non-allocating boundary API
   (fill a caller-owned buffer, skip `closed_ring`, skip the `GI.Polygon`
   wrapper).
2. **31 GiB area** -- one `collect` at `area.jl:275`; a streaming triangulation
   removes essentially all of it.
3. **27 GiB clip output** -- an in-place / area-only intersection operator, plus a
   shared degenerate polygon for the 2.17 M-per-1024^2 disjoint pairs.
4. **21 GiB source geometry** -- same shape as (1) on `rastergrid.jl:376-407`.
5. **10 GiB sparse assembly and 9 GiB descent bookkeeping** -- these are close to
   irreducible; they are proportional to nnz and to the candidate-pair list.

Note that none of this is peak RSS: the eager plan peaks at ~8 GB on the full
tile while churning 206 GiB (`perf-P3-brief.md` D15). On this box the churn is
nearly free; on a machine that cannot grow its heap, it is GC time.
