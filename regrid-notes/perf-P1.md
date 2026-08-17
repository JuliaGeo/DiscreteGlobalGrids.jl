# P1 — performance checkpoint (report only, no code changed)

Method: `regrid-notes/profile-performance.md` (MCP Julia session, BenchmarkTools
baselines with `$`-interpolation, `Profile`/`FlameGraphs` traversed exhaustively,
`Profile.Allocs`, `@inferred`/`@allocated` on the offenders). Julia 1.12.6,
macOS, 8 threads, worktree `../DGG-regrid`, package `lib/GlobalRegridding`,
`ConservativeRegridding` at `agent/use-geometryops-main-cache`, GeometryOps
`eRyQ6`.

## Headline

**One law fails outright and one fails in part.**

- **L2 (locality) FAILS.** Chunk discovery over-reports **2.3–2.6×** at mid size,
  and **any destination chunk that spans a full longitude row connects to 100% of
  the source chunks** — its `SphericalCap` extent degenerates to the whole sphere.
  The "each connected chunk read at most once" half of L2 passes exactly.
- **L3 (bounded residency) FAILS for weights, passes for data.** Data residency is
  one source chunk (verified: 2.38 MiB per warm read against a 49.4 MiB source);
  but `PerChunk()` defaults to unbounded capacity, so every block built stays
  resident (27.0 MiB after one destination chunk, ≈ the whole 158 MiB operator
  after the whole destination). `budget` is carried and ignored — T8's job, but
  today L3 as written does not hold.
- L1, L4, L5 pass cleanly at mid size.

**Plan-build split verdict: the fused accumulating sink is NOT worth building.**
Matrix materialization is **0.24 %** of a conservative plan build (threshold was
5 %). 81 % is tree traversal — and inside it, 67 % of the whole build is
re-synthesizing per-cell spherical caps through the lon/lat chart's `cosd`/`sind`.

**Biggest single win: precompute per-edge `sin`/`cos` on `RasterGrid`.** A scratch
prototype gives a bit-identical cap 5.08× faster (212.0 ns → 41.7 ns), which
removes ≈ 54 % of conservative plan-build wall time.

---

## Setup

| | |
|---|---|
| source array | `DimArray` 3600×1800 `Float64` (X = lon centres, Y = lat centres), 49.4 MiB |
| source space | `RasterGrid`, chunks 512² → 8×4 = **32 chunks**, 6 480 000 cells |
| destination **T** ("tiled") | `RasterGrid` 360×180, chunks 60×60 → 6×3 = **18 chunks**, 64 800 cells |
| destination **B** ("band") | `RasterGrid` 360×180, chunks 360×20 → **9 chunks** × 7200 cells |
| chunked source | `P1Counting <: DiskArrays.AbstractDiskArray` over the same parent, `GridChunks` 512², counting every `readblock!` (a copy of `test/test_lazy.jl`'s `T7Counting`; the test file was not touched) |
| build counter | `P1CountingMethod` wrapping `Conservative()`, counting `build_weights!` calls |

Destination **B** exists because destination **T**'s chunks interleave in
cell-position space and therefore cannot be addressed one at a time at all — see
finding **F5**.

A 1/9-scale pair (1200×600 → 120×60, source chunked 512²) is used where a
measurement needed many repetitions; it is labelled wherever it appears.

---

## A1 — eager whole-domain baselines (mid size)

| measurement | value |
|---|---|
| `@timed plan_regrid(src; to=T, from=src_space, method=Conservative())` | **133.25 s**, 45.5 GiB allocated, 1.28 s GC |
| `Sys.maxrss()` delta over that call | **1335 MiB** (for a 158 MiB result — 8.4×) |
| resulting block | 64 800 × 6 480 000, **nnz = 7 120 800**, 158.1 MiB |
| `@benchmark regrid($data, $plan)` (apply only) | median **23.04 ms**, min 2.00 MiB, **28 allocs**, GC 0 % |
| build ÷ apply | **5783×** |

A one-shot `regrid` is 99.98 % plan build.

## A2 — the conservative plan-build split

### Phase split by wall clock

Measured with a phase-timed replica of `build_weights!(::Conservative, …)` in the
scratch session (no package file edited); mid size, whole domain.

| phase | s | % of build |
|---|---:|---:|
| `subtree` construction (both sides) | 0.003 | 0.00 |
| operator construction | 0.000 | 0.00 |
| **`CR.intersection_areas`** (traversal + clipping + CR's own assembly) | **130.87** | **99.83** |
| `findnz` of CR's matrix (`conservative.jl:271`) | 0.051 | 0.04 |
| COO fill — `addweight!`/`adddenom!` (`conservative.jl:272-277`) | 0.079 | 0.06 |
| `WeightBlock` → `sparse()` (`plans.jl:49`) | 0.094 | 0.07 |
| **total** | **131.10** | |

### Attribution inside `intersection_areas` (flamegraph, thread 1, 9015 samples @ 10 ms)

Percentages are of the whole build; recursion-safe (topmost-occurrence inclusive
totals, so the recursive `dual_depth_first_search` is not double-counted).

| frame | % of build |
|---|---:|
| `get_all_candidate_pairs` — **dual tree descent** | **77.14** |
| ⤷ `_rastercellcap` (`src/rastergrid.jl:579`) — per-cell cap synthesis | **67.24** |
| ⤷⤷ `LonLatToSphere` (`src/rastergrid.jl:23`) | 58.74 |
| ⤷⤷ `spherical_distance` | 4.71 |
| ⤷ `_rectcap` (`src/rastergrid.jl:582`) — node extents | 4.12 |
| `assemble_sparse_matrix_coo` (CR) | 14.69 |
| ⤷ `_run_and_store!` → `BlockAreaOperator` (`src/conservative.jl:207`) | 14.62 |
| ⤷⤷ `getcell` (`src/rastergrid.jl:411`) — cell polygon synthesis | 6.74 |
| ⤷⤷ `_sutherland_hodgman_intersection` (GO) | 4.65 |
| ⤷⤷ `GO.area` of the clip | 3.24 |
| ⤷ CR's own sparse construction (14.69 − 14.62) | 0.07 |
| `SparseArrays.sparse`/`sparse!` anywhere | 0.14 |

`LonLatToSphere` is **66.38 %** of the entire build; within it `cosd` 36.1 %,
`sind` 29.6 %, `rem` 10.5 %.

### The split the plan asked for

| bucket | % of plan build |
|---|---:|
| (a) polygon clipping / intersection (SH clip + `GO.area` + polygon synthesis) | **14.6** |
| (b) tree traversal + extent synthesis | **81.3** |
| (c) `WeightCOO` fill (`addweight!`/`adddenom!`/`push!`) | **0.06** |
| (d) `sparse()` assembly + `WeightBlock` (ours **and** CR's, plus `findnz`) | **0.24** |

**Verdict — sparse assembly is not material (0.24 % ≪ 5 %). Do not build the
fused accumulating sink.** Close that fallback. A perfect fusion would save
0.24 % of a one-shot `regrid`. (The double assembly is still worth removing, but
for peak memory, not time — see **P4**.)

## A3 — the lazy path (mid size, destination **B**, counting source)

| measurement | value |
|---|---|
| `plan_regrid(...; lazy=true)` + `LazyRegridArray` | 10.8 ms, 1.04 MB, **0 source reads, 0 blocks built** |
| **cold** read of band 1 (polar, 7200 cells, 8 connected source chunks) | **59.68 s**, 15.24 GiB allocated, 8 `readblock!` (8 distinct), 8 `build_weights!`, 27.0 MiB of blocks retained |
| **warm** read of the same chunk (blocks cached) | median **5.67 ms**, min 2.38 MiB, **168 allocs**, **0 `build_weights!`**, 8 `readblock!` again |
| cold ÷ warm | ≈ 10 500× |
| values | warm result `==` cold result, and `max|lazy − eager|` = 1.0e-15 at 1/9 scale |

Cold lazy profile has the *same shape* as the eager build:
`build_weights!` 93.2 %, of which descent 81.6 % and `_rastercellcap` **73.3 %**;
clipping 3.1 %, `GO.area` 2.3 %; `_applyslices!`, `blockreference!`,
`WeightBlock`, `cellindices`, `subtree` each **0.0 %**.

---

## Law-by-law

### L1 — construction is free — **PASS**

`plan_regrid(src_chunked; …, lazy=true)` + `LazyRegridArray(...)` at mid size:
`length(counting.reads) == 0`, `counting.bytes == 0`, `nblocks(storage) == 0`,
`method.builds == 0`. 10.8 ms and 1.04 MB, all of it the destination's chunk
extents and spans (`lazy.jl:145-160`, T7 deviation 3).

### L2 — locality — **FAIL** (passes only the "at most once" clause)

Passing half: reading destination band 1 issued **8** `readblock!` calls for **8**
connected source chunks, all distinct — each connected chunk read exactly once,
none outside the connected set.

Failing half — the connected set is far larger than the geometry:

| destination | discovered pairs | geometric truth | over-read |
|---|---:|---:|---:|
| **T** (18 tiles of 60×60) | 205 | 78 | **2.63×** |
| **B** (9 bands of 360×20) | 218 | 96 | **2.27×** |

Per band, discovered vs true source chunks (of 32):

| band | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| discovered | 8 | 16 | **32** | **32** | **32** | **32** | **32** | 18 | 16 |
| true | 8 | 8 | 16 | 8 | 8 | 16 | 8 | 16 | 8 |

Root cause, measured: destination band cap radii are
`[25.1, 56.2, 180.0, 180.0, 180.0, 180.0, 180.0, 56.2, 25.1]°`. `_boxcap`
(`src/rastergrid.jl:575`) and `_mergecaps` (`src/conservative.jl:142`) return
`_WHOLE_SPHERE` once a cap passes π/2, which any 360°-wide box does. Bands 3–7
therefore connect to **every** source chunk. `DilatedIntersects`
(`src/discovery.jl:50`) can only compare caps, so it cannot recover.

This is exactly the shape a global raster is normally chunked in (full rows), and
a DGGS polar chunk will hit the same wall.

### L3 — bounded residency — **PASS for data, FAIL for weights**

Data (T7's claim, **verified**): the warm read of a 7200-cell destination chunk
against 8 source chunks allocates **2.38 MiB / 169 allocations**, of which
2.06 MiB is the single reused source buffer (`_fitbuffer`, `src/lazy.jl:331`,
2 allocations because the last X chunk is 16 wide, not 512) and 0.22 MiB the four
per-destination-chunk accumulators (`src/lazy.jl:249-252`). Peak live-heap delta
over a **cold** read ≈ **13.3 MiB**, against `connected × chunk = 16.0 MiB` and a
whole source of 49.4 MiB. The source array is never materialized, and connected
chunks are loaded one at a time — T7's residency claim holds.

Weights: `PerChunk()` defaults to `capacity = maxbytes = typemax(Int)`
(`src/plans.jl:158-163`), so nothing is ever evicted. After one destination chunk
the plan holds **27.0 MiB** of blocks; the post-read live-heap delta is 33.3 MiB.
Extrapolated over all 218 discovered pairs that is the full 158 MiB operator plus
denominators — i.e. residency grows with the destination, not with the chunk.
`budget` is validated and then ignored (`ChunkedPlan`, `src/plans.jl:280-289`;
T7 "seams for T8").

Eager path for contrast: `Sys.maxrss()` delta **1335 MiB** for a 158 MiB result
(8.4×) — see **P4**.

### L4 — plan reuse — **PASS**

Second read of the same destination chunk: **`build_weights!` calls = 0**,
`nblocks` unchanged at 8, result identical. (Source `readblock!` calls are 8
again — data is deliberately never cached, T7 deviation 4.)

### L5 — hot loops clean — **PASS**

- Eager apply flamegraph (`@bprofile regrid($src, $plan)`, 5757 samples): the only
  runtime-dispatch frame is `wrapoutput` at **`src/executor.jl:394`**, 0.05 %,
  one `DimArray` construction outside every loop. GC self-time 0.14 %.
  The whole profile is `_accumulate!` (`src/executor.jl:166-177`) and the
  `getindex`/`setindex!` it inlines.
- `@inferred applyblock!(...)` and `@inferred finalize!(...)` succeed on the real
  64 800 × 6 480 000 block; `@allocated applyblock!(...) == 0`.
- Warm lazy read `Profile.Allocs` at `sample_rate = 1.0`: **169 allocations,
  2.342 MiB**, all O(chunk) buffers —
  `_fitbuffer` `src/lazy.jl:331` (2 allocs / 2.06 MiB),
  `_fitmatrix` `src/lazy.jl:326` (2),
  `_readdestination!` `src/lazy.jl:251` and `:252` (2),
  `_readsource!` `src/lazy.jl:338` (16 zero-byte `UnitRange` boxes from the
  `sr...` splat). **No per-element allocation anywhere.**
- The `readblock!` frame the profiler flags as dispatch belongs to the
  REPL-defined counting double, not to the package.

---

## C — method baselines on the same pair (eager, mid size)

| method | plan build | allocated | nnz | apply (median) | apply allocs |
|---|---:|---:|---:|---:|---:|
| `Conservative()` | **133.25 s** | 45.5 GiB | 7 120 800 | 23.04 ms | 28 |
| `NearestCell()` | **34.8 ms** | 107.8 MiB | 64 800 | 6.04 ms | 28 |
| `BilinearPoint()` | **38.9 ms** | 124.6 MiB | 259 200 | 6.30 ms | 28 |

Build is 3830× cheaper — orders, as expected, nothing to flag there.

**Flag on apply**: a 110× difference in nonzeros buys only 3.8× in apply time.
Fitting the two points gives **5.88 ms fixed + 2.41 ns per nonzero**; the fixed
part is 0.91 ns × 6 480 000 source columns. `_accumulate!` (`src/executor.jl:169`)
loops `for k in axes(W, 2)` whatever the block contains. See **P7**.

---

## Top self-cost frames — conservative plan build (thread 1)

| % self | frame |
|---:|---|
| 12.61 | `muladd` `float.jl:500` |
| 10.42 | `cosd` |
| 8.95 | `*` `float.jl:497` |
| 7.90 | `-` `float.jl:496` |
| 4.98 | `sind` |
| 3.85 | `isinf` `float.jl:724` |
| 3.11 | `sind` `trig.jl:1209` |
| 2.31 | `rem_internal` |
| 2.23 | `getproperty` |
| 2.22 + 2.16 | `rem` `float.jl:600/601` |
| 1.80 | `_naive_triangulated_spherical_ring_area` (GeometryOps) |

Every one of the top six is inside `LonLatToSphere` / `spherical_distance` under
`_rastercellcap`.

## Runtime-dispatch frames

| % total | frame |
|---:|---|
| 2.06 | `_naive_triangulated_spherical_ring_area`, GeometryOps `src/methods/area.jl:285` |
| 0.05 | `wrapoutput`, `src/executor.jl:394` (eager apply, once per call) |

Nothing else. The package's own kernels are dispatch-free.

## Top allocation sites — plan build

1/9 scale (1200×600 → 120×60), 4.51 GiB total for 789 600 nonzeros = **5.7 KiB
per weight**. `Profile.Allocs` at `sample_rate = 0.001`, scaled ×1000.

| MiB | allocs | site |
|---:|---:|---|
| 1677 | 11.6 M | `_cellring` **`src/rastergrid.jl:425`** (`push!` growth on an empty `USPoint[]`) |
| 644 | 10.5 M | `_naive_triangulated_spherical_ring_area` GO `methods/area.jl:275` |
| 394 | 12.9 M | `getcell` **`src/rastergrid.jl:417`** (`Polygon`/`LinearRing` wrappers; entered through `Trees.getcell` at `:679`) |
| 358 | 7.9 M | `_intersection_sutherland_hodgman` GO `clipping/sutherland_hodgman.jl:410` |
| 336 | 1.5 M | `_cellring` **`src/rastergrid.jl:428`** |
| 197 | 6.5 M | `_cellring` **`src/rastergrid.jl:423`** |
| 176 | 6.6 M | `similar` `broadcast.jl:228` |
| 262 | 5.1 M | `_intersection_sutherland_hodgman` `:421`, `:417`, `:424` |

By type: `Memory{UnitSphericalPoint}` 2.78 GiB, `Vector{UnitSphericalPoint}`
493 MiB, `Vector`+`Memory{LinearRing}` 569 MiB.
`RasterGrid`'s own polygon synthesis (`_cellring` + `getcell`) is **2604 MiB =
56 % of every byte the build allocates**. GC is only ~1 % of wall time, so this is
memory pressure and threading headroom, not latency.

---

## Prioritized proposals (no edits made)

### P1 — precompute per-edge `sin`/`cos` on `RasterGrid`; build caps and corners by table lookup

- **Evidence**: `_rastercellcap` is **67.2 %** of the eager build and **73.3 %** of
  a cold lazy chunk read; `LonLatToSphere` is 66.4 % of the build, of which
  `cosd`+`sind`+`rem` are 76 %. `_boxcap` (`src/rastergrid.jl:553`) evaluates the
  chart 8× per cell cap (4 boundary points, twice — once for the centroid, once
  for the radius) and `_CELL_CAP_SAMPLES = 0` means all four are **corners**, i.e.
  exactly edge-vector entries.
- **Prototype** (scratch only): four `Float64` vectors `cosd.(xedges)`,
  `sind.(xedges)`, `cosd.(yedges)`, `sind.(yedges)`; a corner is
  `(cy*cx, cy*sx, sy)`. Radius bit-identical on a probe cell.
  **212.0 ns → 41.7 ns, 5.08×.**
- **Estimated win**: removes ≈ **54 %** of conservative plan-build wall time
  (131 s → ≈ 60 s at mid size) and ≈ 59 % of a cold lazy chunk read.
- **Cost**: O(nx + ny) memory and construction — exactly what `RasterGrid`'s
  docstring already promises.
- **Where**: `src/rastergrid.jl` (`_boxpoint`, `_boxcap`, `_rastercellcap`,
  `getcell`, `cellcentroid`). Gate on the chart: add a
  `chartedgetables(transform, xe, ye)` hook returning `nothing` for a general
  transform so a projected raster keeps today's path.
- **Absorb into**: an orchestrator fixup before wave 4 (it is `rastergrid.jl`
  only, and T6 already owns seam fixes there).

### P2 — memoize cell caps per chunk in the lazy path

- **Evidence**: the mid-size cold read of destination band 1 costs **75.7 µs per
  weight** against **18.4 µs** for the eager whole-domain build — **4.1×** — with
  an identical profile shape (73 % `_rastercellcap`). A destination chunk's cell
  caps are re-synthesized once per connected source chunk (8× here). At 1/9 scale,
  where a source chunk is 36 % of the source rather than 4 %, lazy and eager are
  level (13.7 vs 16.1 µs/nnz), so the penalty scales with the number of pairs.
- **Estimated win**: up to 4× on cold lazy reads at mid size; with P1 applied the
  lazy-vs-eager penalty should vanish.
- **Where**: `src/rastergrid.jl` (`RasterCellTree.child_indices_extents`,
  `src/rastergrid.jl:673`) or a per-chunk cap vector cached beside the block in
  `src/plans.jl`. **T8.**

### P3 — fix the whole-sphere chunk extent (this is the L2 failure)

- **Evidence**: destination band cap radii `[25.1, 56.2, 180, 180, 180, 180, 180,
  56.2, 25.1]°`; bands 3–7 discover **32/32** source chunks against a truth of
  8–16. Aggregate over-read 2.27× (bands) and 2.63× (tiles).
  `_boxcap` `src/rastergrid.jl:575` and `_mergecaps` `src/conservative.jl:142`
  both bail to the whole sphere past π/2.
- **Estimated win**: 2.3–2.6× fewer block builds *and* source reads on the lazy
  path — which, since block building is 93 % of a cold read, is a 2.3–2.6×
  wall-clock win on lazy materialization; 4× on an equatorial full-row chunk.
- **Options**: (a) let a space override chunk connectivity — `RasterGrid` can
  answer exactly with a lon/lat box-overlap test; (b) represent a chunk extent as
  a small *set* of caps (split a wide box into ≤ 90° longitude runs);
  (c) keep caps and add a cheap box test to `DilatedIntersects`
  (`src/discovery.jl:46-51`).
- **Where**: `src/discovery.jl` + `src/rastergrid.jl`. **T8**, or a fixup — this
  is the largest lazy-path win after P1.

### P4 — remove the double sparse assembly (memory, not time)

- **Evidence**: `build_weights!` (`src/conservative.jl:267-277`) takes CR's
  assembled `SparseMatrixCSC`, `findnz`es it, re-pushes every entry into a
  `WeightCOO`, and `WeightBlock` (`src/plans.jl:49`) `sparse()`s it again. Time
  cost 0.17 %. But at mid size the simultaneously live intermediates are CR's
  matrix (166 MiB) + the `findnz` triple (171 MiB) + the COO (171 MiB, up to
  342 MiB across `push!` growth) + the final matrix (166 MiB) for a 158 MiB
  result. Measured `Sys.maxrss()` delta over the call: **1335 MiB, 8.4× the
  result**.
- **Fix**: `BlockAreaOperator` (`src/conservative.jl:207`) already holds both index
  maps — have it `push!` straight into the `WeightCOO`, or let `WeightBlock` adopt
  CR's matrix and derive `denom` from its nonzeros. Either removes two of four
  copies.
- **Where**: `src/conservative.jl`. Small and safe; good fixup.

### P5 — close the fused-accumulating-sink fallback

Matrix materialization totals **0.24 %** of a conservative plan build (CR's sparse
0.07 + `findnz` 0.04 + COO fill 0.06 + our `sparse()` 0.07). The named fallback in
the P1 task description should be recorded as measured-and-declined. **No work.**

### P6 — `_cellring` allocates by growth

- **Evidence**: `src/rastergrid.jl:423-429` starts from an empty `USPoint[]` and
  `push!`es up to five times; the three lines account for 2.21 GiB of the 4.51 GiB
  a 1/9-scale build allocates (19.6 M allocations). Prototype with exact-size
  allocation: `getcell` **141.1 ns → 118.3 ns, 432 B → 240 B, 5 → 4 allocs**.
- **Estimated win**: ~1 % of build time, ~45 % of the build's allocated bytes.
- **Where**: `src/rastergrid.jl`. Fold into P1 (same function neighbourhood).

### P7 — the apply loop iterates source columns, not nonzeros

- **Evidence**: `_accumulate!` `src/executor.jl:169` (`for k in axes(W, 2)`) walks
  all 6 480 000 columns whatever the block holds. Fitted apply cost =
  **5.88 ms fixed + 2.41 ns/nnz**, so `NearestCell`'s 64 800 nonzeros (worth
  0.16 ms) still cost 6.04 ms. In the lazy path each block is
  `ndst_chunk × nsrc_chunk`, so a 7200 × 262 144 block pays 262 k column
  iterations per source chunk — and P3's over-reported pairs are mostly *empty*
  blocks paying full column cost.
- **Fix**: store the block destination-major (CSC of `Wᵀ`) so accumulation gathers
  over `nnz + ndst`; at minimum hoist the `nzrange` emptiness test above
  `_value(src[k])`.
- **Estimated win**: point-method apply 6 ms → < 1 ms, conservative 23 ms →
  ≈ 18 ms. Matters for N-D reuse (12 monthly slices = 276 ms today) and for the
  lazy path's empty blocks. Low priority against P1–P3.
- **Where**: `src/plans.jl` (`WeightBlock` layout) + `src/executor.jl`.

### P8 — a tiled `RasterGrid` destination has no addressable chunk

- **Evidence**: destination **T** (`RasterGrid` 360×180 chunked 60×60) reports
  `eachchunk = GridChunks(RegularChunks(64800, 0, 64800))` — one chunk for the
  whole cell axis — because its chunks are not contiguous runs of cell positions
  (`_cellchunks`, `src/lazy.jl:183`). Its chunk spans overlap
  (`1:21300`, `61:21360`, …), so `_coveringchunks` (`src/lazy.jl:277`) answers a
  request for tile 1 with **six** destination chunks. Consequently
  `regrid!(dest, data, plan::ChunkedPlan)` (`src/lazy.jl:373`) iterates one chunk
  and materializes the entire destination in a single `readblock!` — there is no
  streaming granularity for the most natural raster destination, and the whole
  205-pair discovery runs in one call.
- **Fix**: address destination chunks by chunk number rather than by
  cell-position span (`_coveringchunks` could intersect `cellindices`; a
  `materialize!` could iterate `1:nchunks(dst_space)` directly), or have `regrid!`
  drive the chunk loop instead of `eachchunk`.
- **Where**: `src/lazy.jl`. Blocks **T8**'s streaming test and **T10**'s
  raster→raster tutorial; DGGS destinations (contiguous `descendant_range`) dodge
  it, so **T9** will not surface it.

---

## Upstream findings (report, do not fix here)

- **GeometryOps `src/methods/area.jl:266-286`** —
  `_naive_triangulated_spherical_ring_area` is commented "streaming iteration (no
  allocation)" but `collect(Iterators.map(...))` at `:273` allocates a
  `Vector{UnitSphericalPoint}` per ring: **644 MiB / 10.5 M allocations** at 1/9
  scale, and `Iterators.peel`'s `Union{Nothing,Tuple}` is the **only real
  runtime-dispatch site in the whole plan build** (2.06 %).
- **GeometryOps `src/methods/clipping/sutherland_hodgman.jl:410/417/421/424`** —
  ≈ 620 MiB allocated per 790 k clips.

## Reproduction

Everything above came from one MCP Julia session on
`lib/GlobalRegridding` with `test/toyspaces.jl` included; the counting disk array
and the counting method are scratch copies of `test/test_lazy.jl`'s `T7Counting`
and `T7CountingMethod`. No file under `lib/GlobalRegridding/src` or
`lib/GlobalRegridding/test` was modified.
