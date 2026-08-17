# P2 — final performance checkpoint (report only, no code changed)

Method: `regrid-notes/profile-performance.md`. MCP Julia sessions on
`<worktree>/test` and `<worktree>/docs` for the interactive work;
`julia --project=... --threads=N` fresh processes for every number that has to
be free of a warm heap or of this session's thread count. Julia 1.12.6, macOS
(Darwin 25.5.0, Apple silicon, 8 CPU threads), worktree `../DGG-regrid` at
`claude/regrid` **6f40cbe**, `lib/GlobalRegridding`, CR
`agent/use-geometryops-main-cache`, GeometryOps git-tree `77dcc16f`.

Judged against `regrid-notes/perf-P1.md` throughout; F1 and F2 verified, not
re-litigated.

---

## Headline — the bad news first

1. **L2's over-report regressed 2.4 % on one shape and improved 31 % on the
   other.** P1's tile-chunked raster destination (360×180 in 60×60 tiles) now
   discovers **210** chunk pairs against a geometric truth of **78** — P1
   measured 205. The cause is measured, not guessed: F2's corner-only
   `_sampledcap` gives a 60°-wide box a 38.50° cap where the sixteen-sample
   construction gives 37.76°, and rebuilding both sides' chunk extents with
   sixteen samples reproduces **exactly 205**. The band-chunked shape, which is
   the one P1 called out, went **218 → 150** against a truth of 96 (2.27× →
   1.56×). Net verdict is still a large win; the regression is small, local, and
   has a one-line fix (see **Q1**).
2. **T10's "12-tile truth" at the south pole is not the truth.** Four of the
   twelve tiles the acceptance file asserts are read produce **zero** weights:
   the destination chunk's cells occupy longitudes `[-180°, -30°) ∪ [90°, 180°)`
   only (2401 of 2401 centroids; none in the 120° gap), so the true fan-in is
   **8** tiles. 1.44 MB of the 4.32 MB that chunk reads is wasted. The
   assertion is still a good under-connection mutant-killer, but it pins the
   current cap rather than the geometry.
3. **A cold lazy chunk read barely threads: 1.01 s at one thread, 0.63 s at
   eight (1.6×)** — against 5.3× for the whole-domain eager build. Blocks are
   built one after another; the only parallelism is inside a single 2401 × 90000
   pair, which is too small to fill eight threads. New since P1 (there was no
   threading at all then).
4. **Destination cell geometry is recomputed ≈25× per cold chunk build.**
   `cell_boundary` is **42.6 %** of a cold south-pole build (0.43 s of 1.01 s);
   computing every one of the chunk's 2401 cell boundaries exactly once costs
   **17.3 ms**. This is the single biggest remaining lever and it is in the main
   package, not in `lib`.

Everything else is green, and the reproductions are clean: **L1–L5 all PASS**,
P1's two failing laws are fixed, and F2's headline numbers reproduce to the
digest.

---

## Law verdict

| law | P1 | P2 | the P2 number |
|---|---|---|---|
| **L1** construction is free | PASS | **PASS** | `plan_regrid` **3.95 ms / 2.62 MB**, `regrid(...; lazy=true)` **4.21 ms / 3.13 MB** on an 8 235 432-cell destination with 3432 chunks. `readblock!` calls **0**, bytes read **0**, `residency(A).loads == 0`, `nblocks == 0`, `build_weights!` **0**. (P1: 10.8 ms / 1.04 MB on an 18-chunk destination.) |
| **L2** locality | **FAIL** (2.27–2.63× over-report; whole-sphere caps) | **PASS**, with a residual superset | *at most once*: exact — over 27 destination chunks read (25 sampled + polar + mid-latitude) `loads == discovered == unique reads` in every case, no read outside the connected set. *over-report*: **1.32×** aggregate at acceptance scale (127 discovered / 96 true over 25 chunks), 1.50× at the pole, 3.0× worst sampled. On P1's own pair: bands **1.56×** (was 2.27×), tiles **2.69×** (was 2.63×). No chunk extent is the whole sphere any more. |
| **L3** bounded residency | **FAIL** for weights, PASS for data | **PASS** both halves | streamed south-pole read at `budget = 2^20`: resident source **360 000 B = exactly one tile** (held: 4 320 000 B = twelve), weight cache **977 320 B in 1 block**. At `budget = 2^24` the byte bound bites exactly: **4 138 648 B / 5 blocks** under a 4 194 304 B share. Values `isequal` at every budget. `Spilled(dir; capacity = 2)`: 12 files, 2 blocks resident, reads 2 and 3 rebuild **nothing** and return `isequal` values. |
| **L4** plan reuse | PASS | **PASS** | second and third read of the south-pole chunk: **0** `build_weights!` (counting method), same values. 3-slice N-D source: **12 builds for 36 chunk-pair applications**; second read 0 builds. |
| **L5** hot loops clean | PASS | **PASS** | warm lazy read **154 allocations / 4.42 MiB, GC 0 %** over 1.08 M source elements; `@allocated applyblock! == 0` on plain, NaN and sentinel paths at 2401 × 90000 / nnz 30 863; `@allocated finalize! == 0` for both policies; **no frame in `GlobalRegridding` or `DiscreteGlobalGrids` carries the `runtime_dispatch` bit** in either the eager apply (2821 thread-1 samples) or the threaded build (16 216 samples). The only real dispatch site is still GeometryOps' `_naive_triangulated_spherical_ring_area`. |

---

## A — law by law at the T10 acceptance scale

Setup is T10's, unchanged: `DEMTiles` 3600 × 1800 `Float32` in **300² tiles**
(12 × 6 = 72, 360 000 B each), synthesized per block and recording every
request; destination `IGeo7System()` area-matched to **level 7**, 8 235 432
cells, **3432 chunks** of 2401. South-pole chunk = **3374**, cells
`8093774:8096174`. `test/systems/crosssystem/regrid_acceptance.jl` runs
**22 pass / 0 fail** in this session (5.1 s).

### L1 — PASS

Both the plan and the lazy array are geometry and bookkeeping. 3.95 ms and
2.62 MB buys the chunk extents and spans of a 3432-chunk destination; the source
is never touched (`length(reads(DEM)) == 0`), no block exists
(`nblocks(plan.storage) == 0`), and `residency(A).loads == 0`.

### L2 — PASS on "at most once", residual superset on "only the connected set"

**The clause that passes exactly.** For every destination chunk measured, the
number of `readblock!` calls equals the number of distinct source chunks equals
the number of discovered pairs. No chunk is read twice within one read and none
outside the connected set is read at all.

**The clause that is a superset.** Truth here is not a cap argument: a
discovered pair whose block has `nnz == 0` contributes nothing, so
`count(nnz > 0)` *is* the geometric truth. Measured by reading the chunk and
inspecting its blocks:

| destination chunk | discovered | true (nnz > 0) | over-report |
|---|---:|---:|---:|
| 3374 (south pole) | 12 | **8** | 1.50× |
| 245 (12°E, 42°N) | 2 | **1** | 2.00× |
| 25 chunks sampled on a 4 × 7 lon/lat lattice | **127** | **96** | **1.32×** |
| worst of the 25 (chunk 2517) | 12 | 4 | 3.00× |

The polar case is worth stating plainly, because T10 asserts the discovered
number. Source tiles 66–69 (longitudes −30° … +90°) return **zero** nonzeros;
tiles 61–65 and 70–72 return 23 730 / 30 863 / 29 654 / 18 274 / 4 152 /
1 589 / 6 372 / 13 677. The chunk's own cells confirm it: centroid latitudes
−89.96 … −85.93, and **0 of 2401** centroids fall in the 120° gap. So the
fan-in is 8 tiles, the file reads 12, and 33 % of the chunk's read volume is
wasted.

**Against P1's own pair** (source 3600×1800 chunked 512² → 32 chunks;
destination 360×180 either as 9 bands of 360×20 or as 18 tiles of 60×60). Truth
was recomputed from the whole-domain 7 120 800-nonzero block by mapping every
nonzero to its chunk pair, and comes out at **78** (tiles) and **96** (bands) —
exactly P1's figures, so the comparison is like for like:

| | P1 discovered | P2 discovered | truth | P1 ratio | P2 ratio |
|---|---:|---:|---:|---:|---:|
| bands (360×20) | 218 | **150** | 96 | 2.27× | **1.56×** |
| tiles (60×60) | 205 | **210** | 78 | 2.63× | **2.69×** |

Per band, discovered: P1 `[8, 16, 32, 32, 32, 32, 32, 18, 16]` →
P2 `[8, 15, 16, 16, 24, 24, 16, 16, 15]`. Band cap radii
`[25.1, 56.2, 180, 180, 180, 180, 180, 56.2, 25.1]°` →
`[20, 40, 60, 80, 100, 80, 60, 40, 20]°`. F1's fix did what it said.

The tile regression is F2's `_sampledcap`, and the attribution is measured:
destination tile caps are **38.50°** by the corner construction against
**37.76°** by `_boxcap(..., _BOX_CAP_SAMPLES = 3)`, source chunk caps **31.36°**
against **30.84°**, and running the pairwise cap test with sixteen-sample caps
on *both* sides gives **205** pairs — P1's number to the pair.

### L3 — PASS, both halves

Budget sweep on the south-pole chunk, values `isequal` at every row:

| `budget` | `weightbudget` | resident weights | `databudget` | peak resident source | loads |
|---|---:|---:|---:|---:|---:|
| 2^20 | 262 144 | **977 320 B / 1 block** | 786 432 | **360 000 B (1.00 tile)** | 12 |
| 2^22 | 1 048 576 | 977 320 B / 1 block | 3 145 728 | 360 000 B (1.00 tile) | 12 |
| 2^24 | 4 194 304 | **4 138 648 B / 5 blocks** | 12 582 912 | 4 320 000 B (12 tiles) | 12 |
| 2^26 | 16 777 216 | 11 154 832 B / 12 blocks | 50 331 648 | 4 320 000 B (12 tiles) | 12 |
| 2^30 | 268 435 456 | 11 154 832 B / 12 blocks | 805 306 368 | 4 320 000 B (12 tiles) | 12 |

Read the first two rows honestly: a single block here is **977 320 B**, larger
than the 262 144 B nominal weight share, and `PerChunk` never evicts the block
it just touched. So the bound below one block's size is not the byte figure but
**one block** — which is exactly what L3 asks for ("one source chunk + one
`WeightBlock` + the output block"). Above one block's size the byte bound is
respected exactly (row 3: 4 138 648 ≤ 4 194 304).

`Sys.maxrss()` across the same read, fresh single-threaded process: **+4.41 MiB**
streamed and **+0.00 MiB** held (the warm-up had already grown the heap). That
is the whole argument for T8's `residency`: `maxrss` cannot resolve a 360 000 B
working set, and `residency` reports it exactly and for free.

**`Spilled` with a bounded deserialized-block LRU.** `Spilled(dir; capacity = 2)`
at `budget = 2^20`, counting method in the plan: first read **0.928 s, 12 builds,
12 files, 2 blocks / 1 837 760 B resident**; second read **0.074 s, 0 builds**;
third **0.030 s, 0 builds**. All three results `isequal` to each other and to
the held (unspilled) answer. Streaming residency stays 360 000 B throughout.

### L4 — PASS

Second and third reads of the south-pole chunk build **zero** weights (counted
through a `P2Method` wrapper around `Conservative`, not through the storage's own
miss counter, so a disk reload cannot be mistaken for a build).

**Slice reuse.** A 3-slice N-D source (3600 × 1800 × 3, chunked 300 × 300 × 1)
read as `A3[POLARCELLS, :]`: **12 `build_weights!` calls for 36 chunk-pair
applications** — one spatial plan, reused across slices. 36 source loads, peak
resident 12 960 000 B (12 tiles × 3 slices) at `budget = 2^30`. Second read: 0
builds, identical values. Slice 1 matches the 2-D reference to Float32 rounding
(max |Δ| 2.4e-4 on ~2000 m elevations) and consecutive slices differ by
9.99976 … 10.00012, which is the synthetic +10 m per slice.

### L5 — PASS

**(i) The warm lazy read loop.** 12 blocks cached, `budget = 2^30`:

- `@benchmark A[POLARCELLS]` — median **27.86 ms**, **4.42 MiB**, **154
  allocations**, GC **0 %**; `build_weights!` still 12 afterwards.
- `Profile.Allocs` at `sample_rate = 1.0`: **158 allocations / 4.41 MiB** — the
  12 held source tiles (4.32 MB) plus the tile accumulators. Over 1.08 M source
  elements and 2401 destination cells that is **no per-element allocation
  anywhere**.
- CPU profile, thread 1, 358 samples: **zero runtime-dispatch frames, zero GC**.
  `readblock!` on the test DEM is **94.4 %** of the read (`elevation` 90.5 %),
  and reading the same 12 tiles with nothing else attached costs **27.33 ms**
  against the read's 27.86 ms. The executor's own work does not clear the 0.3 %
  floor — a warm lazy read at this scale *is* the source.
- On a real block (2401 × 90000, nnz 30 863): `@inferred applyblock!` and
  `@inferred finalize!` succeed; `@allocated` is **0** for `applyblock!` on the
  plain, NaN and sentinel paths and **0** for `finalize!` under both `Weighted`
  and `Extensive`. **2.05 ns/nnz**; `_walknonzeros` correctly routes this block
  (nnz·256 = 7 900 928 > 90 000 columns) to the column walk.

**(ii) The threaded eager path.** Mid-size pair, 64 800 × 6 480 000, nnz
7 120 800.

- *Apply*, thread 1, 2821 samples: 92.1 % `applyblock!` → `_accumulate!`
  (`src/executor.jl:280-286`), 7.3 % `anyinvalid`, 0.2 % `finalize!`; GC 0.21 %.
  **No dispatch-flagged frame in package source at all** — P1's single 0.05 %
  `wrapoutput` (`src/executor.jl:394`) does not even appear. 28 allocations /
  2 098 000 B per apply, identical to P1's 28 / 2.00 MiB.
- *Build*, 8 threads, 16 216 samples across all threads: dispatch-flagged self
  time excluding the unnamed task-root frame is **4.13 %, all of it GeometryOps
  `_naive_triangulated_spherical_ring_area` (`src/methods/area.jl:285`)** — P1
  measured 2.06 % of a serial build; the share of total CPU roughly doubles
  because everything around it got faster. **Zero** dispatch frames in
  `GlobalRegridding` or `DiscreteGlobalGrids`. GC self **6.32 %** of total CPU,
  `poptask` 3.69 %.

**The empty-retry guard (new since P1).**

- *It costs nothing when it does not fire.* On a real pair (nnz 4152)
  `GR._intersectionareas` measures **27.75 ms** against **27.99 ms** for a
  direct threaded `CR.intersection_areas` on the same trees —
  `judge` says **−0.85 % ⇒ invariant**. A `try` with no throw is free.
- *It does not fire on the acceptance case.* All **12** south-pole pairs,
  including the **4** that produce zero nonzeros, return normally from the
  threaded call: an empty *block* is not an empty *reduction*, because the caps
  still meet at node level and tasks still spawn.
- *It does fire on over-reported pairs whose descent spawns nothing.* The
  mid-latitude chunk's spurious source tile raises `empty-reduction` and is
  rebuilt serially: **0.447 ms** through attempt-and-retry against **0.082 ms**
  for a direct serial build (F2 measured 0.465 / 0.020 — reproduced). The total
  cost is `spurious pairs × 0.45 ms`: across the 25 sampled destination chunks
  that is 31 pairs ≈ **14 ms against 4.6 s of reads**. Not a problem; recorded
  because it is a real path, not a theoretical one.

---

## B — reproduced baselines

Every timed comparison below is a **fresh process** at a stated thread count,
except where noted as in-session.

### P1's mid-size eager conservative build (3600×1800 → 360×180)

| | P1 | P2 @1 thread | P2 @8 threads | `judge` vs P1 |
|---|---:|---:|---:|---|
| `plan_regrid` build | **133.25 s** | **37.83 s** | **7.15 s** (in-session 7.851 / 7.688 / 7.760) | −71.6 % / −94.6 % ⇒ improvement |
| allocated | 45.5 GiB | 33.26 GiB | 34.19 GiB | −27 % |
| GC | 1.28 s | 1.38 s | 1.01 s | |
| nnz | 7 120 800 | **7 120 800** | **7 120 800** | identical |
| `Sys.maxrss` delta | 1335 MiB (P1) / 1640.7 MiB (F1, post-fix) | **1548.8 MiB** | **1504.8 MiB** | −5.6 % vs F1 |
| apply median | 23.04 ms, 2.00 MiB, 28 allocs | 20.95 ms | 21.80 ms (in-session 22.2 ms, 2 098 000 B, **28 allocs**) | −9.1 % / −5.4 % ⇒ improvement (i.e. noise) |

`judge(P2 apply in-session, P1 apply)` = **−3.65 % ⇒ invariant**, memory ratio
1.0004, allocations identical. The apply loop has not moved since P1 and was not
meant to.

**New**: the block digest `hash((colptr, rowval, nzval, denom))` is
**17190587724216142268 at one thread and at eight** — F2 asserted threaded ≡
serial bit for bit on a 2832-nonzero pair; it holds at **7 120 800** nonzeros.
(The digest differs from F1's mid-size value because this raster's latitude axis
runs the other way, which permutes cell positions; `nnz` is the invariant that
says it is the same block.)

Threading speedup on this pair: **5.29×** (37.83 → 7.15 s).

### F2's hydrology DEM (Copernicus N46/E010 ×4 → 900×900 → `PartialGrid(IGeo7System, 0061154, 12)`, 823 543 cells)

| | F2 | P2 | `judge` |
|---|---:|---:|---|
| `plan_regrid` @8 threads | 2.5 s (2.49 / 2.53 / 3.22) | **2.327 / 2.360 / 2.349 s** | −6.0 % ⇒ improvement (noise) |
| `plan_regrid` @1 thread, fresh process | 15.5 s | **14.81 / 16.18 / 15.59 s** | **+0.6 % ⇒ invariant** |
| nnz | 2 259 258 | **2 259 258** | identical |
| digest | 904854995744869422 | **904854995744869422** | identical |
| end to end `DGG.regrid(dem; to = grid)` | 2.9 s | **2.42 s** | |
| covered cells / range | 823 056 / 674.2192–3858.329 m | **identical** | |
| `maxrss` delta @1 thread | — | 125.0 MiB | |

F2's numbers reproduce on this machine in a fresh session, at both thread
counts, with a bit-identical operator. Threading speedup 6.6× (F2: 6.2×).

---

## C — the remaining known costs, at acceptance scale

Serial CPU profile (`--threads=1`, fresh process) of a **cold lazy build + read**
of one destination chunk. 9919 samples over 8 south-pole reads and 4673 samples
over 60 mid-latitude reads, 0.5 ms nominal delay. Inclusive, topmost occurrence,
as a share of the whole read.

| frame | south-pole chunk | mid-latitude chunk | owner |
|---|---:|---:|---|
| `build_weights!` | **90.4 %** | 81.1 % | — |
| ⤷ dual descent `get_all_candidate_pairs` | 45.2 % | 43.9 % | GeometryOps |
| ⤷ clip + assembly `assemble_sparse_matrix_coo` → `_run_and_store!` | 43.7 % | 30.8 % | CR |
| **`cell_boundary`** (IGEO7 destination geometry, every path) | **42.6 %** | **41.6 %** | **main pkg** `src/systems/…/z7grid.jl` |
| ⤷ **`snyder_inv_xyz`** | **28.3 %** | **27.7 %** | **main pkg** |
| DGGS-cursor extents: `node_extent` + `child_indices_extents` | 13.7 + 14.6 = **28.3 %** | 19.7 + 15.7 = **35.3 %** | **main pkg** `src/fallbacks/cursor.jl` |
| ⤷ `cell_cap` = `points_cap ∘ cell_boundary` | **20.0 %** | **29.8 %** | **main pkg** `src/fallbacks/caps.jl:81` |
| raster-side extents `_rectcap` + `_rastercellcap` | 7.2 + 5.9 = 13.1 % | 4.8 + 3.7 = 8.5 % | lib |
| `getcell` → `cell_polygon` (destination polygons for the clip) | 30.4 % | 17.5 % | main pkg |
| Sutherland–Hodgman clip | 8.9 % | 9.9 % | GeometryOps |
| `_naive_triangulated_spherical_ring_area` | 4.0 % | 3.0 % | GeometryOps |
| `_fillcoo!` + `sparse()` | 0.15 + 0.51 % | 0.21 + 0.39 % | **lib** |
| `applyblock!` / `_accumulate!` | **0.04 %** | 0.06 % | **lib** |
| `_readsource!` (the IO the whole design is about) | 2.6 % | 6.5 % | lib |
| GC self | **11.4 %** | 12.9 % | — |

`cell_boundary` and `cell_cap` overlap: a cell cap *is* the cap of the cell's
boundary, so the 20.0 % sits inside the 42.6 %. They are listed separately
because they are two different callers of the same synthesis.

**Both of F2's leftovers confirmed, and both are bigger here than on the
hydrology pair.**

1. **The DGGS cursor's node extents — F2's ≈11 %, here `cell_cap` alone is
   20.0 % (pole) / 29.8 % (mid-latitude).** The reason it is worse at level 7 is
   structural: the south-pole chunk's `HierarchicalGridCursor` subtree has
   **2401 leaves for 2401 cells — one cell per leaf**. `child_indices_extents`
   (`src/fallbacks/cursor.jl:265`) therefore allocates a one-element
   `Vector{Tuple{Int,Cap}}` and computes one `cell_cap` on *every visit to that
   leaf*, and GeometryOps' `_child_extents` cache — which covers only the second
   tree — has nothing to amortize. `cell_cap` itself
   (`src/fallbacks/caps.jl:81`) is `points_cap(cell_boundary(grid, c))`, i.e. a
   full inverse-Snyder boundary synthesis per call.
2. **`snyder_inv_xyz` — F2's ≈23 %, here 28.3 %.** IGEO7 cell boundaries, and
   still paid identically by `CR.Regridder`, so it is not a gap against the
   baseline — it is the ceiling on this case.

**Where they live.** Neither is in `lib/GlobalRegridding`. `cell_cap` and
`cells_cap` are `src/fallbacks/caps.jl`; `node_extent`/`child_indices_extents`
for the cursor are `src/fallbacks/cursor.jl`; `snyder_inv_xyz` and
`cell_boundary_cartesian` are the IGEO7 system under `src/systems/`. The lib's
own share of a cold build is `_rectcap` + `_rastercellcap` (13.1 %) plus 0.7 %
of COO/sparse assembly; its executor is 0.04 %.

**What a fix looks like (proposal only, nothing applied).**
*Memoize `cell_boundary` per cell for the lifetime of one `subtree` / block
build.* Evidence, measured on the south-pole chunk's own 2401 cells:

| one pass over all 2401 destination cells | time | bytes |
|---|---:|---:|
| `cell_boundary` | **17.32 ms** | 3.02 MiB |
| `cell_cap` | 14.37 ms | 3.20 MiB |
| `cell_polygon` | 14.55 ms | 3.72 MiB |

A cold south-pole read takes **1.01 s** serially, of which `cell_boundary` is
42.6 % ≈ **0.43 s** — a redundancy factor of **≈25×**. One boundary cache keyed
by cell position, living on the cursor built by `subtree` and dying with it,
serves `cell_cap`, `cell_polygon` and the clip from a single pass. It costs
**3 MiB per 2401-cell chunk**, which is the granularity `budget` already
governs, and it needs no change in `lib`. A vectorized or cached inverse
projection would attack the same 28.3 % from below and is the only other lever
on `snyder_inv_xyz`.

---

## D — new since P1 and F2

### D1 — the lazy path barely threads (1.6× of the 5.3× the same machine gives the eager build)

Fresh processes, same script, cold read of one destination chunk:

| | 1 thread | 8 threads | speedup |
|---|---:|---:|---:|
| south-pole chunk (12 pairs) | 1.024 / 0.925 / 1.008 / 1.012 s | 0.632 / 0.625 / 0.646 / 0.625 s | **1.6×** |
| mid-latitude chunk (2 pairs) | 0.055 – 0.067 s | 0.039 – 0.050 s | **1.4×** |
| whole-domain eager build (for contrast) | 37.83 s | 7.15 s | **5.3×** |

F2 put threading *inside* one `build_weights!`. That is the right place for a
`DirectPlan`, which is one block, but a `ChunkedPlan` destination chunk is a
sequence of 2–12 small blocks and each of them is too small to fill eight
threads. `build_weights!`'s docstring already promises "no state outside the
call, so concurrent builds of different chunk pairs are independent", and
`PerChunk` already documents that "concurrent readers of different keys are
safe" — so the pairs of one destination tile are spawnable as they stand.
The interactions to think about are `SourceHold` (which is per `readblock!`
call and not currently task-safe) and the fact that two levels of `@spawn` would
oversubscribe, so the inner `True()` would want to become conditional on the
outer fan-out. **Owner: `lib/GlobalRegridding` (`src/lazy.jl`,
`src/conservative.jl`).** Estimated win 2–3× on cold lazy materialization, which
is the operation the whole design exists to make cheap.

### D2 — the tile-shape discovery regression, and its exact cause

Covered under L2 above. Restating as an actionable item: `_rectcap` is called
in two very different regimes — **once per chunk** when `chunkextents` builds
the discovery tree (72 + 3432 calls per plan, cold), and **once per visited node
pair** inside a block build (hundreds of thousands, hot). F2's corner-only
`_sampledcap` is exactly right for the second and measurably loose for the
first. Passing the sample count down, or giving `chunkextent` the
sixteen-sample construction unconditionally, recovers P1's 205 pairs for an
O(nchunks) cost that is invisible against a 2.3-second build. **Owner:
`lib/GlobalRegridding/src/rastergrid.jl`.**

### D3 — GC is 11–13 % of a cold lazy chunk build, and the bytes are destination geometry

P1 recorded GC at ~1 % of the eager build's wall time and did not report a GC
share for a cold lazy read. At acceptance scale a cold south-pole read spends
**11.35 %** of its samples in GC (mid-latitude 12.9 %), and
`Profile.Allocs` at `sample_rate = 0.01` puts **≈486 MiB** through the allocator
for that single 2401-cell destination chunk:

| MiB | type |
|---:|---|
| 241.3 | `Memory{UnitSphericalPoint{Float64}}` |
| 75.5 | `Memory{Tuple{Float64,Float64,Float64}}` |
| 56.8 | `Vector{UnitSphericalPoint{Float64}}` |
| 23.9 + 22.5 | `Vector`/`Memory{LinearRing{…}}` |
| 15.3 | `Vector{Tuple{Float64,Float64,Float64}}` |

That is the same object population C's boundary cache would stop re-creating —
the ≈25× redundancy and the 11 % GC are one finding seen twice. (`Array
boot.jl:648` and `jl_alloc_genericmemory_unchecked` are the only sites the 1.12
allocation profiler resolves; the type histogram is what names the owner.)

### D4 — the upstream dispatch site is unchanged and now a larger share

`GeometryOps` `src/methods/area.jl:285`,
`_naive_triangulated_spherical_ring_area`: P1 reported it as the only real
runtime-dispatch site in the whole build at 2.06 %. It is now **4.13 %** of the
threaded build's total CPU and 4.0 % of a cold lazy build. Reported, not fixed —
it is upstream, and P1 already filed it.

---

## Prioritized proposals (nothing applied)

| # | proposal | evidence | est. win | owner |
|---|---|---|---|---|
| **Q1** | Cache `cell_boundary` per cell for the life of one `subtree`/block build; serve `cell_cap`, `cell_polygon` and the clip from it | `cell_boundary` 42.6 % of a 1.01 s cold build ≈ 0.43 s, against 17.3 ms for one pass over the chunk's 2401 cells (**≈25×**); 486 MiB and 11.4 % GC on the same read | up to ~40 % of a cold lazy chunk build, and most of the GC | **main pkg** `src/fallbacks/{cursor,caps}.jl` |
| **Q2** | Spawn the chunk pairs of one destination tile; make the inner `True()` conditional on the outer fan-out | cold chunk read 1.01 s @1t vs 0.63 s @8t (1.6×) against 5.3× for the eager build | 2–3× on cold lazy materialization | **lib** `src/lazy.jl`, `src/conservative.jl` |
| **Q3** | Use the sixteen-sample cap for **chunk** extents only; keep the corner cap for tree nodes | 210 pairs vs 205 with sixteen samples on both sides; truth 78. Chunk extents are O(nchunks) per plan, node extents are O(node pairs) per build | 2.4 % fewer pairs on tile-chunked rasters, no build-time cost | **lib** `src/rastergrid.jl` |
| **Q4** | Loosen T10's polar assertion from "12 tiles read" to "≥ the 8 that carry weight, none outside the southern row" | 4 of 12 tiles yield `nnz == 0`; 0 of 2401 chunk centroids lie in their longitude gap | correctness of the acceptance claim, 33 % of that chunk's IO | test file (`test/systems/crosssystem/regrid_acceptance.jl`) |
| **Q5** | Upstream: `_naive_triangulated_spherical_ring_area`'s `Iterators.peel` union | 4.13 % of the threaded build's CPU, the only dispatch site anywhere | ~4 % | GeometryOps (already filed by P1) |

`_fillcoo!` + `sparse()` are **0.66 %** of a cold build and 0.24 % of the eager
one; P1's declined fused-accumulating-sink stays declined.

---

## State of performance

The implementation phase closes with all five laws holding at the acceptance
scale, and with the two that P1 recorded as failures fixed at their root rather
than papered over: chunk extents no longer degenerate to the whole sphere
(band-chunked over-report 2.27× → 1.56×, aggregate 1.32× at acceptance scale),
and residency is now bounded on both halves — one source tile and one weight
block at a tiny budget, with the byte bound respected exactly wherever it
exceeds a single block, and identical values at every budget and through the
spill path. The hot loops are what P1 found them to be and slightly better: no
package frame anywhere carries a runtime-dispatch flag, `applyblock!` and
`finalize!` allocate zero on every path including the sentinel one, and a warm
lazy read of a 2401-cell chunk against 1.08 M source elements costs 154
allocations, 98 % of which are the source's own buffers. The eager path is
**18.6×** faster than P1's baseline at eight threads and **3.5×** at one, on a
bit-identical operator at 7.1 M nonzeros, and F2's hydrology numbers reproduce
to the digest in a fresh session at both thread counts. What is left is not in
this package: a third of a cold chunk build is the IGEO7 destination's own cell
geometry, recomputed about twenty-five times per chunk because nothing caches a
cell boundary between the cap that prunes it and the polygon that clips it —
one bounded, chunk-sized cache in the main package is the last large lever, and
spawning a destination chunk's independent block builds is the second. Neither
is a correctness risk and both are additive; the regridding path as landed is
fast enough to be used and clean enough to be optimized further without moving
a number.

---

## Reproduction

```
cd <worktree>
julia --project=test -e 'using Test; using DiscreteGlobalGrids;
    include("test/systems/crosssystem/regrid_acceptance.jl")'
T = Main.RegridAcceptanceTests
```

Scratch scripts used for the fresh-process numbers (scratchpad only, nothing in
the repo): `eager_fresh.jl` (P1's mid-size build, 1 and 8 threads),
`hydro_fresh.jl` (F2's hydrology build; needs `Rasters.checkmem!(false)` when
`Sys.free_memory()` is small), `profile_acceptance2.jl` (the serial cold-build
profile and allocation profile), `coldread_time.jl` (cold-read wall time by
thread count), `maxrss_polar.jl` (the `maxrss` versus `residency` comparison).

No file under `lib/GlobalRegridding/src`, `lib/GlobalRegridding/test`, `src/` or
`test/` was modified. The counting method (`P2Method`) and the 3-D counting
source (`DEM3`) are session-local copies of the patterns in
`lib/GlobalRegridding/test/test_lazy.jl` and
`test/systems/crosssystem/regrid_acceptance.jl`.
