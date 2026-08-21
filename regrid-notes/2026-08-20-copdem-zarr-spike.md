# Copernicus DEM → IGEO7 level 12 → chunked Zarr: spike record

2026-08-20. Worktree `/home/asinghvi17/geo/DGG-copdem-spike`, branch
`claude/copdem-zarr-spike`, based on `origin/main` `92ef101` (merge of PR #60).
Script: `scripts/copdem_to_igeo7_zarr.jl`. Julia 1.12.6, 64 cores, 751 GB RAM,
`-t 8 --gcthreads=8,1`.

Priority per the user mid-spike: **conservative first, bilinear best-effort**.
Both landed; see §7.

---

## 1. Architecture as built

Chunking is tree-aligned on both sides. No latitude-band mosaics, no ad-hoc
regional rasters.

| | unit | how |
|---|---|---|
| **Source chunk** | one Copernicus DEM 1°×1° tile | `TiledDEM <: DiskArrays.AbstractDiskArray{Float32,1}` over the **complete** level-1 (pixel) grid — 68 947 200 000 pixels — with `IrregularChunks` equal to each level-0 tile's `descendant_range`. Passed as `from = DGGSpace(levelgrid(sys, 1); chunklevel = 0)` → exactly 64 800 source chunks. |
| **Destination chunk** | one IGEO7 level-`La` ancestor subtree | `DGGSpace(dstgrid; chunklevel = La)`. |
| **Write chunk** | the same subtree run | `dggwrite(...; chunks = :auto, chunk_target = 7^(12-La))`. |

The one decision that makes all three coincide exactly:

> **Query the coverage AT the ancestor level, not at level 12.**
> `query(IGeo7System(), MultiOrderCoverage(region); level = La)` emits cells at
> level `La` *or coarser*. Expanding that set with `CellVector(set; level = 12)`
> therefore yields only **complete** subtrees, so every level-`La` run is exactly
> `7^(12-La)` cells. One uniform Zarr chunk length then lands on *every* run
> boundary and `dggwrite`'s manifest reports `ancestor_aligned = true`.
>
> Querying at level 12 instead — the obvious reading of "the region's cells" —
> would leave the boundary ancestors partial and the runs unequal. By the rule in
> `write.jl`, `:auto` would then take `chunklength = best[searchsortedlast(best,
> target)]`, the cumulative length of the first few *partial* runs, and report
> `aligned = false`. No uniform length lands on unequal boundaries.
>
> Price: over-coverage. The destination is the smallest union of whole subtrees
> containing the region, not the region.

`chunks` is never passed to `regrid` (it defeats empty-pair pruning); the tiling
comes from the destination `DGGSpace` alone.

**Level 12 rationale.** IGEO7 level 12 cell ≈ 3685 m² against a GLO-90
equatorial pixel of 8548 m² — the closest level under the closest-area rule.
GLO-30 would want level 13.

**Synthetic fill.** Only tiles named by `real=` are read; every other tile is
generated at its posts from
`z(λ,φ) = 1000 sin 3λ cos 2φ + 500 cos 7λ sin 5φ + 100` (radians in, metres out).
Four real GLO-90 tiles, two of them newly downloaded (2 downloads total):
`S85_00_E000`, `S89_00_E000`, `S90_00_E000`, `S90_00_E001`.

---

## 2. Environment: two broken pins (fixed on the branch, commit `0b7b56c`)

`origin/main`'s `[sources]` revisions are **branch names**, and both branches had
moved under them. A fresh instantiate of the worktree produced code that would
not precompile, and then code that failed at the first weight block:

1. `ConservativeRegridding = {rev = "claude/budget-frontier"}` resolved to tree
   `eb71b25a` — commit `89ec5ca5`, *"Bump version from 0.2.8 to 0.2.9 (#135)"* —
   which predates the frontier work. `GlobalRegridding` needs
   `Trees.split_weight`, added later on that branch.
   → `UndefVarError: split_weight not defined in ConservativeRegridding.Trees`;
   `GlobalRegridding`, `DiscreteGlobalGrids` and `DiscreteGlobalGridsZarrExt` all
   fail to precompile.
2. `GeometryOps = {rev = "as/intersection_area"}` resolved to tree `c30ec910` —
   commit `02750768`, *"Bump patch version to v0.1.44"* — which predates
   `intersection_area` itself.
   → `UndefVarError: intersection_area not defined in GeometryOps`, thrown from
   `lib/GlobalRegridding/src/intersection_area.jl:31`, only once a weight block
   is actually built.

Fixed by pinning both to SHAs across `Project.toml`,
`lib/GlobalRegridding/Project.toml`, `docs/Project.toml`, `test/Project.toml`,
`lib/DiscreteGlobalGridsConformanceTesting/Project.toml`:

* CR → `6a4b997ab2e66dea45b5b62b5b2f32f0a3b279b0` ("Conserve the frontier's work
  estimate across splits") — the tip the earlier campaigns measured.
* GO → `2825c476647cfe6791dbb89fffb8592697b46995` (`as/intersection_area` tip).

**`Pkg.resolve()` does not re-read `[sources]` for a git dep; `Pkg.update(name)`
does.** The first attempt at the fix looked like it had failed for this reason.

This is a live hazard for `main`, not a spike artefact: any fresh checkout of
`origin/main` today can resolve to a non-compiling manifest. The pins want to be
SHAs, or the branches want to stop moving. (`regrid-notes/cr-pin-state` and
`go-pin-state` still apply — these are still temporary pins, just immovable ones.)

`bench` was also made a workspace member with a `Zarr` dependency, in the same
commit.

---

## 3. `_rawids`: verdict, and what it means for one global store

**Survey A is right.** `ext/DiscreteGlobalGridsZarrExt/write.jl`:

```julia
function _cellaxis(src)
    for d in DD.dims(src)
        lk = DD.val(d)
        lk isa Union{CellLookup,ChunkedCellLookup} || continue
        grid = levelgrid(system(lk), level(lk))
        return d, grid, _rawids(grid, lk)      # <- always
    end
    ...
end

_rawids(grid, cells) = (I = idtype(grid); I[convert(I, rawid(c)) for c in cells])
```

There is **one** method of `_cellaxis` and **one** of `_rawids`. The raw-id
vector is materialised before anything else runs, for **every** encoding —
`:ranges`, `:dense` and `:implicit` alike. Everything downstream (`_encoding`
eligibility, `_chunkplan`'s `_runends`, `_coordinate`) reads that array.

Survey B did not notice because a whole-level-6 write is
`ncells(IGeo7, 6) = 1 176 492` ids = 9.4 MB. The cost only bites at this scale.

> **This is orthogonal to the on-disk encoding.** `:ranges` (§5) makes the axis
> *on disk* essentially free — 0.2-0.6 % of a dense axis, and 0.0035 % with
> `merge = :rank`.
> `_rawids` is about **RAM during the write**, and the ranges coordinate is
> computed *from* the dense vector, not instead of it. Choosing `:ranges` or
> `:implicit` does not avoid the allocation.

**Consequence.** A single flat global IGEO7 level-12 store is **not writable** by
today's `dggwrite`:

* `ncells(IGeo7, 12) = 138 412 872 012`
* × 8 B (`Z7Cell` wraps a `UInt64`) = **1.107 TB** for the id vector alone, on a
  751 GB machine. `:implicit` does not help.
* It is also O(n) in *time* — 1.4·10¹¹ `rawid`/`convert` calls before a byte is
  written.

Two ways out, neither attempted here:

* **Shard the store on an ancestor boundary.** One store per level-3 IGEO7
  ancestor is 3432 stores of `7⁹ = 40 353 607` cells (323 MB f32 layer, 323 MB id
  vector); per level-2 ancestor is 492 stores of `7¹⁰ = 282 475 249` cells
  (1.13 GB layer, 2.26 GB id vector). Comfortable either way, and the shard
  boundary is the same tree alignment one level up.
* **Make the axis streamable.** `_rawids` could return a lazy id accessor:
  `_runends` walks it once sequentially, `RangesEncoding`'s `_coordinate`
  compresses it to ranges, `_fill!` already streams. Only `DenseEncoding`'s
  coordinate genuinely needs the whole vector, and that is the encoding nobody
  would pick at this scale. This is the fix worth doing.

---

## 4. Runs

All: GLO-90, IGEO7 level 12, ancestor level 7 (`7⁵ = 16 807` cells/chunk,
65.7 KB f32), `budget = 2^30`, `-t 8`, `merge = :step` unless noted.

`ancestor = 7` rather than the suggested 5 or 6 because the spike needed small
regions and the smallest possible destination is one whole chunk: at level 5
that is 823 543 cells.

| run | region | cells | chunks | src tiles per chunk | write wall | cells/s | peak RSS | store |
|---|---|---|---|---|---|---|---|---|
| **A0** | `cap:-89.9`, all synthetic | 201 684 | 12 | 234 / 245 / 216 / 360 | 131.1 s | 1539 | 2.83 GiB | 0.6 MiB |
| **A** | `cap:-89.9`, 4 real tiles | 201 684 | 12 | 234 / 245 | 128.1 s | 1575 | 2.93 GiB | 0.6 MiB |
| **B** | `box:10,11,45,46` (45 °N) | 2 890 804 | 172 | 2 / 0 | **truncated** at 96/172 chunks by my own 45-min timeout | ~613 | — | — |
| **B2** | `box:10,10.4,45,45.4`, both methods, **uncontended** | 605 052 | 36 | 2 / 0 / 0 | 902.8 s cons, **5.1 s bilinear** | **670 / 118 637** | 2.63 GiB | 1.5 + 1.4 MiB |
| **C** | `box:0,2,-86,-84` (band edge) | 1 815 156 | 108 | 9 / 3 / 2 | 586.0 s | 3098 | 2.46 GiB | 4.1 MiB |
| **D** | `box:10,10.2,45,45.2` both methods | 218 491 | 13 | 2 / 0 / 1 / 0 | 234.2 s cons, 4.2 s bilinear | 933 / 52 000 | 2.48 GiB | 0.53 + 0.51 MiB |
| **E** | `box:10,10.05,45,45.05` both methods | 50 421 | 3 | 3 | 47.0 s cons, 4.0 s bilinear | 1073 / 12 600 | 2.16 GiB | 0.12 MiB |

A0/A/B2/C/D/E: **ALL CHECKS PASSED** (16 checks in the both-methods runs). D and E
ran concurrently with B and C, so their per-cell rates are depressed; **B2 is the
one clean, uncontended mid-latitude measurement** and is the number to quote.
B was killed by a `timeout 2700` I set, not by any failure; its 96 completed
chunks give the ~613 cells/s that B2 then confirmed at 670.

### Store layout, verified on disk (run A0)

```
elevation/.zarray  {"chunks":[16807], "shape":[201684], "dtype":"<f4",
                    "compressor":{"id":"blosc","cname":"lz4","clevel":5,"shuffle":1}}
elevation/0 .. elevation/11              # 12 files; 201684 / 16807 = 12 exactly
cell_chunk_manifest/.zattrs
  {"dggs_chunk_manifest":{"ancestor_aligned":true,"ancestor_level":7,
                          "chunk_length":16807,"grid":"igeo7","level":12,...}}
```

`dggread` round trip: cell count and ids identical in every run; one destination
chunk recomputed from the lazy cube is **bit-identical** to the store in every
run.

### 4.1 The pairing, and the parallelism it does or does not give

Logged per destination chunk (tiles newly decoded):

* **Pole** (`cap:-89.9`): **234, 245, 216, 360** tiles per level-7 subtree. All
  360 tiles of the `S90` row meet at the pole, so one small subtree there
  straddles most of them.
* **Band edge** (−86…−84): 9, 3, 2.
* **Mid-latitude** (45 °N): **2, 0, 1, 0** — the second and fourth chunks needed
  no new tile at all, i.e. the LRU already held everything.

That is exactly the locality the tile chunking is supposed to buy, and it is also
the thing that limits `-t 8`:

> **Parallelism in the lazy plan is over SOURCE CHUNKS within one destination
> tile.** At mid-latitude a destination subtree touches 2–3 tiles, so 8 threads
> can use at most 2–3 of them. Measured: run B sits at **193 % CPU** for its
> whole write, and run B2 — alone on a 64-core box, nothing else of mine running
> — at **185 %** (`/usr/bin/time -v`, 18 m 38 s wall). That is 1.85 of 8 threads
> requested and 1.85 of 64 available. The polar runs, with 200–360 source tiles
> per destination tile, reached ~470 % (15 m 32 user / 3 m 18 real).
>
> So the ideal source chunking for *throughput* is the opposite of the ideal one
> for *locality*, and the mid-latitude case — the bulk of the globe — is the one
> that loses. Parallelising over destination tiles instead (or as well) is the
> obvious lever, and is not something this script can do from outside.

### 4.2 Mid-latitude (runs B, B2)

B2 is the clean measurement: 605 052 cells in 36 chunks, **902.8 s** conservative
(**670 cells/s**), 185 % CPU, peak RSS 2.63 GiB, all 16 checks passed. B — the
same latitude at 4.8× the area — was killed by a `timeout 2700` after 96 of 172
chunks, giving ~613 cells/s, consistent.

Mid-latitude is **the slow case**, at a fifth of the band-edge rate. The reason
is the DEM's own band table: a 45 °N tile is 1200 × 1200 = 1.44 M pixels against
240 × 1200 = 288 k at −85 and 120 × 1200 = 144 k at −90. The source pixels a
destination cell must be intersected against go *up* as the tiles stop being
longitudinally reduced — and mid-latitude is most of the globe.

**Scaling summary, per destination cell**, holding the ancestor level fixed:
work is roughly (source pixels the subtree overlaps) × (cells in the subtree),
and the DEM's band table makes the first factor 10× bigger at the equator than
at the pole while longitude convergence makes the *tile count* 100× bigger at
the pole. The two do not cancel; the mid-latitude case is the expensive one.

### 4.3 Band edge (run C) — no seam

`box:0,2,-86,-84` straddles the 85° band edge, so tiles with 240 columns
(band `[80,85)`, the `S85` row) and 120 columns (band `[85,90)`, the `S86` row)
feed the same ancestor chunks — chunk 1 decoded nine `S86` tiles, chunk 2 three
`S85` tiles.

* **0 non-finite cells** in 1 815 156.
* **0 cells** outside the analytic field's own range over their cell.
* max |value − field(centre)| = **4.6 × 10⁻² m**, RMS 1.3 × 10⁻² m — and the
  worst cell's error is 0.11× the field's own spread across that cell, i.e. the
  residual is the cell-mean-vs-centre-value gap, not error.
* The real `S85_00_E000` tile sits inside this region; the real/synthetic seam
  produced no holes (all cells finite, and the 82.1 % of cells classified
  synthetic-only all satisfy the bracket).

---

## 5. The id axis: `:auto` picks ranges, and what that costs

Verified on every store. The group attribute is
`"dggs": {"coordinate":"cell_id_ranges", "compression":"ranges", ...}` — so
`encoding = :auto` did resolve to `RangesEncoding`, and the store carries
`cell_id_ranges` rather than a dense `cell_ids`.

Measured on run C's store (1 815 156 cells; a dense `UInt64` axis is
**14 521 248 B** raw):

| variant | coordinate rows | coordinate on disk | vs dense raw | total store |
|---|---|---|---|---|
| `:auto` → ranges, `merge = :step` (**default, shipped**) | 259 308 | **26 494 B** | 0.18 % | 4 336 094 B |
| ranges, `merge = :rank` | **44** | **512 B** | 0.0035 % | 4 310 104 B |
| `encoding = :dense` | — | 138 433 B | 0.95 % | 4 447 943 B |
| no compressor at all | 259 308 | 4 149 179 B | 28.6 % | 11 415 774 B |

Layer bytes are 4 304 570 B in all four (7 260 624 B raw f32, 59 %).

**Why 259 308 rows and not 44 under `:step`.** `merge = :step` merges only ids
adjacent **as integers**. An IGEO7 `Z7Cell` packs each base-7 digit into four
bits, so only seven of every sixteen integer values name a cell: **seven siblings
is the longest integer-contiguous run that exists**, and the coordinate is
`ncells / 7` rows however tidy the cell set is. 259 308 × 7 = 1 815 156 exactly.
Run E confirms it at another size: 7203 rows for 50 421 cells.

`merge = :rank` merges runs of consecutive *cells* and collapses the same
destination to **44 rows** — the number of maximal contiguous groups of level-7
subtrees in it, which is what the ancestor-snapped design earns. It is read back
correctly only by a rank-aware reader (this package), so the script leaves the
default at `:step` and exposes `merge=rank` as a flag.

Per-run, as logged by the script (`:step`, the shipped default):

| run | cells | rows | coordinate on disk | dense `UInt64` axis | ratio |
|---|---|---|---|---|---|
| E | 50 421 | 7 203 | 2 372 B | 403 368 B | 0.59 % |
| B2 | 605 052 | 86 436 | 10 157 B | 4 840 416 B | 0.21 % |
| C | 1 815 156 | 259 308 | 26 494 B | 14 521 248 B | 0.18 % |

`rows = cells / 7` exactly in all three.

Even at `:step` the coordinate is a rounding error on disk (0.2–0.6 % of the layer),
because the deltas are perfectly regular and Blosc's shuffle+lz4 eats them. The
practical reading: **ranges vs dense matters far less than it looks; `:step` vs
`:rank` is the real lever, and only for a rank-aware reader.**

`:implicit` writes no coordinate at all, but "needs a whole level" — it is the
right answer for a *global* level-12 store and is unavailable to any of these
partial ones. It does **not** relieve the `_rawids` allocation (§3).

**Naming collision, flagged:** the `dggs` attribute called `"compression"` holds
the **id-axis encoding** (`"ranges"` / `"none"` / `"implicit"`), not any byte
codec. The byte codec lives in each array's `.zarray` `compressor` field and is
Zarr.jl's default Blosc lz4 clevel 5 with byte shuffle. Two different things,
one word. (FYI only, from a measurement made before the item was cancelled:
`BloscCompressor(cname="zstd", clevel=9, shuffle=1)` gives 3 743 176 B for run
C's layer against 4 304 570 B for the lz4 default, −13 %; a *bare*
`ZstdCompressor(level=19)` is **worse** at 5 252 866 B, because Blosc's byte
shuffle is doing most of the work on `Float32` and the standalone codec has no
filter in front of it. `dggwrite` has no knob for this and was deliberately left
without one.)

---

## 6. `ancestor_level` in the manifest is the COARSEST level that fits

Run E (three level-7 subtrees, each the only one under its own level-6 parent)
came back with `ancestor_level = 6, aligned = true, chunk_length = 16807`, not
`ancestor_level = 7`. That is `_chunkplan` doing what it documents: it walks
`A = L-1` upwards while `_maxrun(ends) <= target` and keeps the **coarsest**
level that still fits. Where each level-`La` run happens to also be a complete
level-`La-1` run, it reports the parent — a stronger statement, not a wrong one.

A checker that asserts `ancestor_level == La` will fail intermittently on small
or sparse regions. The right assertion is
`ancestor_level <= La && aligned && chunk_length == 7^(L-La)`.

---

## 7. Bilinear: it works, and it is 177× faster than conservative

`BilinearPoint` requires `hascellchart(source)`, which only `RasterGrid` answers
`true` to, so the DGGSpace source cannot be reused. The script keeps the chunk
*shape* instead of the source:

* walk the same destination ancestor subtrees;
* partition each subtree's cells by the DEM tile their **centre** falls in
  (`cellat(levelgrid(sys, 0), cell_centroid(...))`) — so no cell is ever sampled
  from a raster it lies outside, which is what stops `BilinearPoint` from
  clamp-extrapolating;
* build that one tile as a small `RasterGrid` with a **one-pixel halo** taken
  from its neighbours, and regrid it onto a `PartialGrid` of just that subset;
* expose the whole thing as `BilinearCells <: AbstractDiskArray{Float32,1}` whose
  chunks are the ancestor runs, so `dggwrite` streams it exactly as it streams
  the conservative `LazyRegridArray` — same write unit, no full materialisation.

Halo rule: east/west neighbours always share the tile's band and column count, so
that halo is exact. North/south neighbours only sometimes do; where the band
changes the posts do not line up and the edge row is replicated, degrading to
clamping across that one border. `BAND_EDGE_CLAMPS` counts it — **0** in every
run here, because the mid-latitude test boxes sit inside one band.

Results (run B2, 605 052 cells, all synthetic, uncontended):

* store chunks `(16807,)`, ancestor-aligned, `dggread` round trip clean;
* **0 non-finite cells**; 0 cells outside their cell's analytic field range;
* max |value − field(centre)| = **2.9 × 10⁻³ m**;
* **bilinear vs conservative: RMS 2.6 × 10⁻³ m, max 4.8 × 10⁻³ m,
  correlation 0.99999994** over all 605 052 cells;
* **5.1 s against 902.8 s** for the same cells — **177×**. (Runs D and E, at
  smaller sizes and under contention, gave 56× and 12×; B2 is the clean one.)

The speed gap is the honest headline: conservative is paying for exact spherical
polygon intersection against ~1.4 M pixels per tile, bilinear for four samples
per cell. Whether that trade is acceptable is a data question (bilinear is not
conservative of mass, and a DEM at 61 m cells from 90 m posts is close to a
point-sample regime anyway), not a pipeline question.

Two caveats on the 177×. The bilinear path here is `lazy = false` per
(subtree, tile) pair, so it uses **one thread** — the 118 637 cells/s is a
*serial* number against conservative's 1.85 cores, which makes the per-core gap
larger still. And the tiles it samples are generated in process; a real global
run decodes 26 475 COGs, but at ~0.1 s each that is ~45 min total against days of
compute either way.

**Not covered by these runs:** bilinear across a band edge (the `S85`/`S86`
boundary), and bilinear at the pole — the centre-in-tile partition is defined
there but the per-tile raster's clamp behaviour at ±90 was not exercised. The
band-edge clamp path is written and counted but never fired.

---

## 8. Problems found, in order of how much they matter

1. **`_rawids` caps a single store at ~90 G cells on this machine** (§3). Global
   level 12 needs 1.1 TB for the axis alone. Blocking for the stated goal;
   fixable either by sharding or by streaming the axis.
2. **Both `[sources]` branch pins were broken on `main`** (§2). Anyone
   instantiating `origin/main` today gets a non-compiling manifest.
3. **Thread utilisation collapses exactly where the data is** (§4.1). 185 % of
   8 threads at mid-latitude — 1.85 of the box's 64 — because parallelism is over
   source chunks within a destination tile and a mid-latitude subtree touches
   2–3 tiles. The agreed fix is outer parallelism over destination tiles,
   designed in §10 and deliberately not built here.
4. **`DGGSpace(...; chunklevel = La)` scans the whole global ancestor level.**
   Measured: level 5 → 0.2 s, level 6 → 1.6 s, **level 7 → 12.1 s**, every time a
   space is built, whatever the region's size. For a 50 421-cell destination that
   is 12 s of the 60 s run. `_chunkwindows` iterates `1:ncells(ancestors)`; for a
   `PartialGrid` it could instead walk the grid's own windows.
5. **`ancestor_level` in the manifest is not the level you asked for** (§6).
6. **The word "compression" means the id encoding, not the byte codec** (§5).
7. **A naive analytic oracle is wrong near a pole** (§9) — worth knowing before
   anyone writes a polar regression test against a smooth field.

---

## 9. The oracle: why the first version of it failed, and what replaced it

The first check compared each destination cell against `SYNTHETIC` at the cell
**centre**. On the all-synthetic polar run it failed: max 617 m, RMS 2.97 m,
p99 = 1.0 m, p99.9 = 10 m — with a perfectly correct mean of 100.00 m and the
field's exact range.

That is the oracle, not the regrid. A level-12 cell is ~61 m across, but a degree
of longitude at 89.9° is 194 m and at 89.99° is 19 m, so a cell there spans a
third of a *degree* of longitude, over which `SYNTHETIC` genuinely swings tens of
metres. The cell's mean is not its centre's value. Confirming the diagnosis:
every one of the 18 cells with error > 100 m sits within 10 cell radii of the
pole, and the two worst (617 m, at lon −168.75 and lon 11.25) are 0.64 cell radii
from it.

Replaced with a check that is valid everywhere: **a conservative cell mean is an
average of the field over the cell, so it must lie between the field's smallest
and largest values there.** Sampled at the cell's boundary vertices plus its
centre, with tolerance `1e-3 + 0.05 × spread`. Zero violations in every run,
polar included, and the residual is reported as a fraction of the cell's own
field spread (worst: 0.11× at the band edge, 0.10× mid-latitude, 0.06×
bilinear).

The same geometry broke the *classifier* for "which cells are fed only by
synthetic tiles". Sampling the centroid and the boundary vertices misses a real
tile that slices through a cell without containing a vertex — routine near a
pole, where a tile is narrower than a cell. Fixed by densifying each edge to 16
samples and treating any cell whose cap reaches a pole as touching all 360 tiles
of that row.

---

## 10. Outer parallelism: the agreed follow-up, designed but not built

Not implemented in this spike. This section is the design and the arithmetic;
the only code that landed is `worker_groups`, the partition function, which is
computed and logged but never executed against.

### 10.1 Why it is the first lever

Measured (§4.1): a mid-latitude write sits at **193 % CPU with `-t 8`**, because
the lazy plan's parallelism is over **source chunks within one destination
tile** and a mid-latitude destination subtree touches only 2–3 DEM tiles. Six of
eight threads have nothing to do. The polar case reaches ~470 % because a polar
subtree touches 200–360 tiles — the one place the existing parallelism has
material to work with, and it is 0.4 % of the globe.

Destination tiles, by contrast, are **embarrassingly parallel**: each is an
independent weight build over an independent slice of the cell axis, and this
spike's chunking already makes them disjoint on the write side.

### 10.2 Shape

`W` workers × `T` threads. Workers are **in-process** — a `Threads.@spawn` per
group, each calling `regrid`/`dggwrite` on its own destination subset — not
`Distributed`. In-process matters for the memory arithmetic below: there is one
Julia base image, one set of loaded packages, one set of band tables, and the
per-worker cost is only the working set.

Each worker takes a **group of destination ancestor chunks in the same
geographic area**. `worker_groups(dstspace, sys7, La - k)` gives exactly that:
group the level-`La` chunks by their level-`(La-k)` ancestor. IGEO7 ids sort into
contiguous subtrees, so a coarse ancestor's descendants are neighbours on the
sphere *and* a contiguous run of the cell axis — one partition buys both
properties at once:

* **LRU locality.** A worker walking one area re-reads the same 2–3 DEM tiles for
  every chunk in its group. Measured in run B2: chunk 1 decoded 2 tiles, chunks
  2 and 3 decoded **0**.
* **Write disjointness.** A group is itself an ancestor run, so it is exactly one
  store shard and no two workers ever touch the same Zarr chunk.

`k` is the sizing knob. At `La = 6` and `k = 2`, a group is `7² = 49` chunks of
117 649 cells = 5.8 M cells; the globe has `10·7⁴+2 = 24 012` such groups, which
is a healthy queue depth for any `W`.

### 10.3 Sizing rule

Two independent constraints:

```
CPU:  W · T  ≤  available CPU threads
RAM:  base + W · (budget + LRU + working set)  ≤  0.7 · RAM
```

Calibration from these runs, all at `budget = 2^30` (1 GiB) and `cache = 128`
tiles: **peak RSS 2.16–2.95 GiB for one whole process**, of which ~0.6 GiB is the
Julia base and loaded packages. So

> **per worker ≈ 2.5 GiB** at `budget = 2^30`, `cache = 128`.

The LRU term is worth spelling out because it is counter-intuitive: a
mid-latitude tile is 1200 × 1200 × 4 B = **5.76 MB** and a polar one is
120 × 1200 × 4 B = **0.576 MB**, so 128 cached tiles is **0.72 GiB at 45 °N and
0.07 GiB at the pole**. The LRU is *ten times smaller* at the pole, not larger —
what grows there is the number of *distinct* tiles a chunk wants (360), so a
polar worker wants a bigger `cache` in entries but not in bytes: 360 polar tiles
is 0.2 GiB.

There is also a third, softer constraint from §4.1: **`T` above ~3 is wasted at
mid-latitude.** Raising `T` past the number of source tiles a destination subtree
touches buys nothing, so a split that maximises `W` at small `T` beats an even
one — except at the pole, where `T = 8` is genuinely used.

### 10.4 Worked example 1 — this box, 64 threads / 751 GB

**CPU binds; RAM is free by two orders of magnitude.**

| split | threads used | expected cores busy | RAM |
|---|---|---|---|
| `W=8, T=8` (the safe, even split) | 64 | 8 × 1.85 = **14.8** | 0.6 + 8·2.5 = **20.6 GiB** (2.7 %) |
| `W=21, T=3` (matched to the measured ceiling) | 63 | 21 × min(1.85, 3) ≈ **39** | 0.6 + 21·2.5 = **53 GiB** (7 %) |

`W=21, T=3` is the better pick *on the measurement*: `T=8` hands each worker six
threads it demonstrably cannot use at mid-latitude. `W=8, T=8` is the right pick
for a **polar-heavy** workload, where a single worker really does use ~4.7 cores.
A latitude-aware split — few fat workers near the poles, many thin ones in the
band — is the obvious refinement and needs no new machinery, only a different
`(W, T)` per group batch.

Since RAM is nowhere near binding here, the free upgrade is to **spend it**:
`budget = 4 GiB` and `cache = 512` puts a worker at ~8.5 GiB, and `W = 21` still
only reaches ~180 GiB (24 %).

### 10.5 Worked example 2 — AMD 5900X, 12 cores / 24 threads / 64 GB

**CPU binds here too; RAM has headroom at the default budget and would bind if
the budget were raised.**

| split | threads | expected cores busy | RAM | verdict |
|---|---|---|---|---|
| `W=3, T=8` | 24 | 3 × 1.85 = 5.6 | 8.1 GiB | wastes 6 threads per worker |
| `W=4, T=6` | 24 | 4 × 1.85 = 7.4 | 10.6 GiB | still wasteful |
| **`W=8, T=3`** | 24 | 8 × 1.85 = 14.8 nominal, **~12 real** | 0.6 + 8·2.5 = **20.6 GiB** (32 %) | **recommended** — demand just exceeds the 12 physical cores, which is where you want it |
| `W=12, T=2` | 24 | 12 × 1.85 ≈ 22 nominal, ~12 real | 30.6 GiB (48 %) | no faster than `W=8` once capped, and 50 % more RAM |

The 5900X has 12 physical cores behind 24 SMT threads, and the spherical
intersection kernel is FP-heavy, so **12 is the real ceiling** however the
threads are spelled. `W=8, T=3` reaches it with the least memory; anything
beyond `W=8` is buying queue depth, not throughput.

Recommended: **`W = 8, T = 3`, `budget = 2^30`, `cache = 128`** → 20.6 GiB of
64 GB, a third of RAM, with the other two thirds left for page cache over the
COGs. Do **not** raise `budget` on this box: at `budget = 4 GiB` and
`cache = 512`, a worker is ~8.5 GiB and `W = 8` needs ~68 GiB — over the 64 GB
before headroom. That is the point at which RAM, not CPU, becomes the binding
constraint on the 5900X, and it never does on the 64-core box.

Script side, the memory knobs already exist as config args: `budget=` and
`cache=`. `workers=`/`threads=` would be the two to add when the outer loop is
built; the partition they would consume is `worker_groups`, which is already
there and logged.

### 10.6 How the writes compose

Two options, as posed:

**(a) One store shard per worker, at a coarser ancestor.** Each worker calls
`dggwrite` on its own group into its own directory. No writer change, no locking,
no shared mutable state; shard boundaries are ancestor boundaries so the shards
concatenate in canonical order and a reader can treat them as one axis.

**(b) Concurrent disjoint-chunk writes into one pre-created store.** Would need
`dggwrite` split into a "create and stamp" phase and a "fill this chunk range"
phase, the second callable concurrently. Zarr v2 chunks are separate files, so
the storage layer already permits it; what does not exist is the API.

**Recommendation: (a) now, (b) as the eventual design — and (b) is where the
`_rawids` fix belongs.**

The reason is §3. A single global store is not writable at all today because its
axis is 1.107 TB in RAM, so option (b) *as it stands* does not even work for the
case it exists to serve: pre-creating one global store still runs `_cellaxis`
over 1.4·10¹¹ cells. Sharding is therefore not merely the convenient parallel
answer, it is the same answer §3 already had to give for a completely different
reason, and the two align exactly — the shard level and the worker-group level
can be the same level. Option (b) only becomes worth building **after** the axis
is streamable, at which point a "fill this chunk range" call needs only that
range's ids and the memory problem dissolves with it. Building (b) first would
buy nothing and would still hit the wall.

### 10.7 Expected speedup

Linear in `W`, because the workers are independent, CPU-bound, and neither RAM
nor disk is close to binding (§10.4/10.5; the whole 4.2 MiB/2.9 M-cell store rate
is ~1.5 KB/s of I/O).

From the measured **670 cells/s at 185 % CPU** (mid-latitude, run B2,
uncontended):

| configuration | cores busy | cells/s | global level 12 (1.384·10¹¹ cells) |
|---|---|---|---|
| today, `-t 8` | 1.85 | 670 | 6.5 years |
| `W=8, T=8` | 14.8 | ~5 400 | 297 days |
| `W=21, T=3` | ~39 | ~14 100 | **114 days** |
| BilinearPoint alone, 1 thread (§7) | 1 | 118 637 | 13.5 days |
| `W=21` × BilinearPoint (1 thread each) | 21 | ~2.5·10⁶ | **~15 hours** |

So outer parallelism is worth ~20× on the conservative path and turns "not
viable" into "a long batch job"; the bilinear path is worth another ~180× and
turns it into overnight. Neither is a change to the chunking, the store, or the
axis — all three of which already work at this scale.

Note the last row's shape: bilinear needs no `T` at all in this implementation,
so on a 64-thread box it wants `W = 64, T = 1`, not `W = 21, T = 3`. The
`(W, T)` split is a property of the *method*, not of the machine alone.

---

## 11. Full-scale plan sketch

**Disk.** Global IGEO7 level 12 = 138 412 872 012 cells.

| | |
|---|---|
| dense `Float32` layer | 553.7 GB |
| measured Blosc lz4 ratio on this data | 59–75 % → **330–415 GB** |
| id coordinate, `:implicit` (whole level) | 0 B |
| id coordinate, `:ranges`/`:step` if sharded | ~0.6 % of layer |
| free on `/home` | 2.5 TB — comfortable |

Chunking: at `La = 5`, 168 072 chunks of 823 543 cells (3.29 MB f32);
at `La = 6`, 1 176 492 chunks of 117 649 cells (470 KB). `La = 6` is the better
Zarr chunk size; `La = 5` halves the number of destination-space builds.

**Shape.** One store per level-2 or level-3 ancestor (§3): 492 or 3432 stores,
each written by one `dggwrite` call whose axis fits in 2.26 GB or 323 MB. Shard
boundaries are ancestor boundaries, so the sharding is the same tree alignment
one level up and a reader can concatenate without resorting anything. `:implicit`
becomes available per shard only if a shard is a whole level, which it is not —
so `:ranges` with `merge = :rank` (44-ish rows per shard) is the right pick.

**Time.** Honest extrapolation is uncomfortable. Measured `-t 8` conservative
write rates: 1539 cells/s (pole), 3098 cells/s (−85), **670 cells/s (45 °N,
uncontended, run B2)**. The globe is mostly mid-latitude and mid-latitude is the
slow case, so 670 cells/s is the number to plan against, not the average:

* at 670 cells/s and `-t 8`: 1.384 × 10¹¹ / 670 = **6.5 years**.
* with the §10 outer parallelism at `W=21, T=3`: **114 days**.
* with `BilinearPoint` and `W=64, T=1`: **~15 hours**.

So the pipeline is *correct* at full scale and *not viable* at full scale as it
stands. The gap is not the store, the chunking or the axis — those all work and
are cheap. It is the conservative weight build, at 10⁴–10⁵ source pixels per
destination cell over 1.4 × 10¹¹ cells. Levers, in the order they would pay:

1. **Parallelise over destination tiles**, not only source chunks within a tile
   — designed in §10; the measured 1.85 cores → ~39 is a ~21× that costs nothing
   in accuracy.
2. **Use `BilinearPoint`** where mass conservation is not required: measured
   **177×** faster on identical cells with 4.8 × 10⁻³ m maximum disagreement and
   correlation 0.99999994 (§7).
3. **Reuse weights across shards.** Nothing here persists a weight block, per the
   brief. At full scale the same (tile → subtree) geometry recurs identically for
   every longitude in a band; a band-relative weight cache is the obvious
   structural win and was out of scope for a one-shot spike.

---

## Reproducing

```
export RASTERDATASOURCES_PATH=/home/asinghvi17/geo/DiscreteGlobalGrids.jl/bench/data
cd /home/asinghvi17/geo/DGG-copdem-spike
julia --project=bench -t 8 --gcthreads=8,1 scripts/copdem_to_igeo7_zarr.jl \
    region=cap:-89.9 ancestor=7 method=both out=/tmp/spike
```

`region` is `cap:<max latitude>` or `box:<w>,<e>,<s>,<n>`; `real` is `auto`
(every tile in the cache), `none` (all synthetic), or a list of stems; `merge` is
`step` or `rank`; `pairs` is how many destination chunks to log the source-tile
pairing for. Logs and stores from these runs are under the session scratchpad
(`runA0.log`, `runA.log`, `runB.log`, `runB2.log`, `runC.log`, `runD.log`,
`runE.log`, `spike*/`). `runB2.log` is the clean uncontended mid-latitude run and
the one to quote; `runB.log` was cut short by my own `timeout`.

Branch `claude/copdem-zarr-spike` is pushed to origin, three commits, no PR
opened:

* `0b7b56c` the two SHA pins (§2),
* `075f8ac` the spike script,
* `f1aab8a` `worker_groups`, the §10 partition, computed but not run.

---

## 12. Conservative source-path diagnosis (2026-08-20)

The 670 cells/s of run B2 is not the cost of conservative regridding. It is the
cost of feeding the regridder a **`DGGSpace` source**. Holding everything else
fixed — same destination, same source pixels, same chunk granularity, same
thread count, bit-identical output — a `RasterGrid` source over the same four
DEM tiles runs **44.5× faster**. The mechanism is named to the line in §12.2-3;
§12.4 prototypes the one-method library fix and measures it at **60.8×**, at
parity with the raster path.

Worktree `/home/asinghvi17/geo/DGG-copdem-spike`, branch
`claude/copdem-zarr-spike`, `-t 8 --gcthreads=8,1`, Julia 1.12.6, uncontended.

### 12.1 The A/B

Region `box:10,10.05,45,45.05` (the spike's run E), IGEO7 level 12, ancestor 7:
3 destination chunks, **50 421 cells**. The destination touches four GLO-90
tiles — (45,9), (45,10), (44,9), (44,10) — 5 760 000 pixels. All three arms use
the *same* `dstspace = DGGSpace(dstgrid; chunklevel = 7)`, the same
`Conservative()`/`Weighted(0.5)`, `lazy = true`, `budget = 2^30`, and are timed
warm (second repetition, fresh `regrid` call each time).

| arm | source | wall | CPU | cells/s | vs arm 1 |
|---|---|---|---|---|---|
| **1** | `DGGSpace(levelgrid(CopernicusDEMSystem(90), 1); chunklevel = 0)` + `TiledDEM` — **the spike's path** | **61.42 s** | 309 % | 821 | — |
| **2b** | `RasterGrid` over the 2400×2400 mosaic of those four tiles, chunks **1200×1200** = one tile, i.e. *identical* chunk granularity | **1.38 s** | 308 % | 36 567 | **44.5×** |
| **2** | the same `RasterGrid`, chunks 128×128 | **0.96 s** | 143 % | 52 286 | **63.7×** |

* **The outputs are bit-identical.** `RMS difference 0.000e+00 m, max 0.000e+00 m`
  over all 50 421 cells, 0 non-finite in either. This is purely a performance
  difference, not a semantic one.
* **The work is the same.** A counting wrapper at the `GR.subtree` seam (the
  `CountingSpace` pattern from `lib/GlobalRegridding/test/toyspaces.jl`, since
  PR #51 removed `cellcaptree_builds()`) counts **10 source subtree builds over
  14 400 000 source cells** for arm 1 and **11 builds over 15 840 000 cells**
  for arm 2b. Same number of `(destination tile × source chunk)` block builds,
  same 1.44 M-cell chunks.
* **Thread utilisation is the same** — 309 % vs 308 %. §4.1's "parallelism is
  over source chunks" observation is real but it is *not* what separates these
  two arms. The 44.5× is pure per-core work.
* Arm 2 shows that finer source chunks are worth only a further **1.43×**. The
  granularity of the tile chunking is not the problem.

Two corrections to §4 while we are here:

* Run E's 47.0 s is **not** a regrid measurement. `pairs = 4` reads
  `cube[dstspace.ranges[k]]` for the first four chunks *before* the write clock
  starts, and run E has only three chunks — so every block was already in the
  `PerChunk` storage (10 blocks ≈ 10 MB against a 2^30 budget: nothing was
  evicted) and `writestore` recomputed nothing. Runs A0/D/E are all affected;
  B2, at 36 chunks, is understated by ~11 % and remains the honest number.
* Correspondingly, arm 1's 821 cells/s and run B2's 670 cells/s are the same
  measurement to within the region-size effect.

### 12.2 Where arm 1 spends the time

`Profile` at 2 ms over one arm-1 pass, and at 1 ms over 20 arm-2b passes,
normalised to thread-seconds inside `buildblock` (`plans.jl:420`):

| | arm 1 (DGGSpace) | arm 2b (RasterGrid) | ratio |
|---|---|---|---|
| whole block build | **116.0 s** | **1.40 s** | 83× |
| `_intersectionareas` (`conservative.jl:393`) | 106.4 | 0.656 | 162× |
| — `get_all_candidate_pairs` / `dual_depth_first_search` | **106.1** | **0.466** | **228×** |
| — clipping + areas (`intersection_areas.jl:273,277`) | **0.32** | **0.19** | **1.7×** |
| `GR.subtree` → source tree construction | **9.4** | ≈ 0 | — |
| destination `indexmap` (`shared.jl:25`, via `conservative.jl:331`) | ≈ 0.4 | 0.68 | ≈ 1 |

The last two rows are the point. **The actual spherical clipping — the work the
method exists to do — costs the same in both arms (0.32 s against 0.19 s, on
identical weights).** Everything else is overhead in *finding* the candidate
pairs. Two readings follow immediately:

1. **Building the source tree is only 8.1 % of arm 1.** `CellCapTree`'s eager
   per-pixel polygon + cap pass (`conservative.jl:37-41`, `_cellcap` at
   `conservative.jl:82-84` → `getcell` → `cell_polygon` → `cell_boundary` →
   `closed_ring`, three heap allocations per source pixel, 4.3 M small
   allocations per 1.44 M-cell chunk) is real but **eliminating it entirely
   would buy 1.09×, not 44×.**
2. **The 44× lives in the dual-tree candidate search**, which is 228× more
   expensive for the DGG source. Inside it, the single largest line is the
   *destination* side: `STI.child_indices_extents(cursor::HierarchicalGridCursor)`
   at `src/engine/cursor.jl:272` → `Fallbacks.cell_cap` (`src/fallbacks/caps.jl:81`)
   → IGEO7 `cell_boundary` (`src/systems/IGeo7/system.jl:287` →
   `src/systems/IGeo7/z7grid.jl:205`), **25 074 of 57 978 samples = 43.2 %** of
   arm 1's block time. That cost is *induced*: the destination tree is
   re-descended once per unprunable source node, and arm 1 has far more of them.
   In arm 2b, with the same destination, the same lines cost ~5 %.

### 12.3 Why the source tree cannot prune: two independent defects

`GR.subtree(space::DGGSpace, inds)` (`src/regridding.jl:155-160`) tries three
things in order: the whole-space cached tree, `_chunkcursor`, then
`GR.CellCapTree`. `_chunkcursor` (`src/regridding.jl:163-173`) demands

```julia
root = treeify(space.grid)
root isa HierarchicalGridCursor && root.selection === nothing || return nothing
```

and `DGG.treeify` for a CopernicusDEM level grid returns a **`BlockCursor`**
(`src/systems/CopernicusDEM/cursor.jl:309`). So the fast path is dead code for
this source and **every** block falls through to `CellCapTree`. What that tree
gives us, measured directly on one 1200×1200 GLO-90 tile at 45 °N against the
`RasterCellTree` the raster arm gets for the same 1.44 M cells:

| depth | node cells | `CellCapTree` cap radius, median (max) | `RasterCellTree`, median |
|---|---|---|---|
| 0 | 1 440 000 | **274.9 km** | **68.1 km** |
| 3 | 180 000 | 224.2 (226.2) | 29.5 |
| 6 | 22 500 | 204.8 (220.0) | 8.5 |
| 8 | 5 625 | 194.8 (218.4) | 4.24 |
| 11 | 703 | **121.0 (203.9)** | **1.82** |
| 13 | 176 | 5.76 (176.5) | 0.93 |
| 17 | 11 | 0.39 (71.2) | 0.24 |

A level-7 IGEO7 destination subtree has a cap radius of **6.556 km**. The
`RasterCellTree` median drops below that at depth 8; the `CellCapTree` median
not until depth 13, and its *maximum* stays above 70 km all the way down. A
single-cap descent against that one destination cap visits **11 729 nodes /
20 978 leaf cells** in the `CellCapTree` against **2 743 nodes / 15 325 cells**
in the `RasterCellTree`. In the real dual-tree join the destination tree is
re-descended at every one of those extra nodes, so the two multiply — hence 228×
rather than 4×.

Two separable causes, both in `lib/GlobalRegridding/src/conservative.jl`:

* **(D1) `_cellcapnode` (`conservative.jl:43-56`) bisects the linear index
  range** — `mid = (lo + hi) >> 1`, no spatial sort. CopernicusDEM level-1 cell
  order is raster row-major, so every internal node above ~depth 11 is a set of
  complete pixel *rows*: a strip spanning the tile's full 1° of longitude
  (78.7 km at 45 °N), whose bounding cap is ~39 km wide however thin the strip
  is. Longitude pruning is therefore impossible over the top ~11 of 17 levels.
  `RasterCellTree` bisects the longer side of the 2-D index rectangle
  (`rastergrid.jl:825-836`) and does not have this problem.
* **(D2) `_mergecaps` (`conservative.jl:62-79`) is not a bounding cap.** It
  centres on the *mean* of the child cap centres and takes
  `max(distance + radius)`, so slop compounds bottom-up. Measured inflation on
  the same point set: root **274.9 km against a true 68.1 km — 4.0×**; at depth
  11 (one pixel row, true half-width 39.4 km) **121 km — 3.1×**. Every factor of
  2 in radius costs about two levels of pruning depth, so D2 alone is worth ~2
  levels ≈ 4× in visited nodes.

For completeness, the destination-side item is a third, smaller defect:
`STI.child_indices_extents(cursor::HierarchicalGridCursor)`
(`src/engine/cursor.jl:265-275`) materialises a fresh `Vector{Tuple{Int,Cap}}`
and recomputes every leaf cap from `cell_boundary` on **each** visit. Its
docstring says the result "is materialized to avoid recomputing caps during
repeated dual-tree passes", but the materialisation lives only inside one call —
there is no memo across visits. `MemoRasterTree`
(`lib/GlobalRegridding/src/raster_tree_memo.jl:62-86`) is exactly that memo for
rasters and has no DGG counterpart. `TileCells` (`conservative.jl:120-189`)
caches destination *polygons*, not caps, and only below
`_TILE_CELL_CACHE_MAX = 1 << 16` (`conservative.jl:111`) — La = 7 (16 807 cells)
keeps it, La = 6 (117 649) loses it.

### 12.4 The fix, prototyped and measured

The proposed change (F1 below) was prototyped by overwriting one method — no
library edit, no commit:

```julia
function GR.subtree(space::DGG.DGGSpace, inds::AbstractUnitRange{<:Integer})
    # exact chunk range of a CopernicusDEM level-1 grid chunked at level 0?
    k = searchsortedfirst(space.starts, Int(first(inds)))
    if k <= length(space.ranges) && space.ranges[k] == inds
        r, q, _, _ = CD.decode(sys, space.chunkids[k])
        return CD.BlockCursor(grid, sys, CD.Bisected(), 1, Int64(-1), r, r, q, q,
                              0, N - 1, 0, nc - 1, true)   # the tile's whole raster
    end
    ...
end
```

Same region, same destination, same 10 block builds — **all 10 took the
`BlockCursor`** — and the same output (`sum = -1.0228710316057205e6`, identical
to arms 1 and 2b to every digit):

| arm | warm wall | cells/s | vs arm 1 |
|---|---|---|---|
| 1 — `CellCapTree` (today) | 61.42 s | 821 | — |
| **1F — windowed `BlockCursor`** | **1.01 s** | **49 774** | **60.8×** |
| 2b — `RasterGrid`, same process | 1.23 s | 41 033 | 50× |

So the fix does not merely narrow the gap, it **closes it**: a `DGGSpace` source
over CopernicusDEM ends up at parity with — here marginally ahead of — the
`RasterGrid` source over the same pixels, because a `BlockCursor` node is an
index rectangle with an O(1) exact box cap, exactly as `RasterCellTree` is.
(1F and 2b ran in one process against each other while another benchmark held
8 of the box's 64 cores; the clean arm-2b figure is the 1.38 s of §12.1.)

That is the measurement that turns the work plan below from an estimate into a
plan: **F1 is worth ~60× on its own and is the smallest change in the table.**

### 12.5 Calibration: the reference GLO-30 case on this machine

The A/B above uses synthetic 90 m tiles and a 50 421-cell destination, so it
wants anchoring against the campaign configuration everyone else quotes: the
**cached real GLO-30 tile** `Copernicus_DSM_COG_10_N45_00_E010_00_DEM.tif`
(3600×3600) onto its covering IGEO7 **level-13** MOC, 16 490 405 cells,
`Weighted(0.5)`, `-t 8 --gcthreads=8,1`, via `bench/harness.jl`'s `run_case`.
Box load average 20-22 throughout.

| configuration | wall | cells/s | reference |
|---|---|---|---|
| **eager, in-memory raster**, `chunkcells = 4096` — *the reference config* | **51.5 s** cold, **43.1 / 44.0 s** warm (plan + run) | 320 k – 383 k | `regrid-notes/2026-08-19-post-stack-profile.md`: **50.8–51.8 s** warm at t8 |
| lazy, `CountingDisk` 128×128 source chunks, `chunklevel = 7` | 183.6 s | 89 843 | `bench/results/crfix-size-*.ndjson` `size_3600`: **163.5–165.5 s** at t8 |
| lazy, 128×128 source chunks, `chunklevel = 6` | 126.8 s | 130 090 | — |
| lazy, 512×512 source chunks, `chunklevel = 6` | 243.7 s¹ | 67 663 | — |

¹ that row is the first regrid in its process and so carries first-call
compilation; treat it as an upper bound.

`nnz = 66 728 744` in the eager rows, against post-stack's 66 129 312 for a
coverage of 16 181 892 cells — same operator, marginally larger MOC here.

**The reference case reproduces**: 43–44 s warm against the recorded 50.8–51.8 s
and the ~56 s expected, i.e. this tree is if anything slightly faster than the
last recorded baseline. So the raster arm of §12.1 is the known-good
configuration, not a strawman, and the 44.5×/60.8× is measured against a path
that is behaving exactly as the campaigns say it should.

Two incidental readings: **source-chunk size is a real lever on the raster path
too** (512² → 128² measures 1.9× at La = 6, an upper bound per note 1), and **coarser destination chunking wins**
(La = 6 beats La = 7 by 1.45×) even though La = 6's 117 649 cells/tile exceed
`_TILE_CELL_CACHE_MAX = 1 << 16` (`conservative.jl:111`) and so lose the
`TileCells` polygon cache. The spike's La = 7 keeps that cache; if a future
design moves to La = 6 for Zarr-chunk reasons, raising or removing that constant
should be measured at the same time.

### 12.6 Work plan: making a DGGSpace source fast in the library

Ranked by measured share of the gap. Every item has a raster-side analogue that
already exists and works; none of them is a change to the spike, the chunking,
the store or the axis.

| # | change | closes | difficulty | lands in |
|---|---|---|---|---|
| **F1** | **Let `_chunkcursor` window a `BlockCursor`.** `_chunkcursor` (`src/regridding.jl:163-173`) requires a `HierarchicalGridCursor`; CopernicusDEM's `treeify` gives a `BlockCursor` (`src/systems/CopernicusDEM/cursor.jl:309`), which is already an O(1) lazily-bisecting tree over **2-D tile/pixel rectangles** with an exact box cap (`STI.node_extent` → `_box_cap`, `cursor.jl:245-248`) and a correct `Trees.split_weight` (`cursor.jl:289`) — structurally the same object `RasterCellTree` is for rasters. `_block_cursor` (`cursor.jl:340-377`) already constructs the exact node we want ("one tile, whole raster rows", `cursor.jl:364-369`); it just takes a `PartialGrid` rather than a chunk of a `LevelGrid`. Needs a small hook (e.g. `DGG.subcursor(root::BlockCursor, chunkid, inds)`) plus the `isa` relaxation. | **measured 60.8× (§12.4)**: kills both the 8.1 % construction and the 228× candidate-search blowup at once, because `Bisected` (`cursor.jl:26-30,164-165`) splits the longer axis exactly as `RasterCellTree` does, so the node caps become §12.3's `RasterCellTree` column. The residual against arm 2b is F4: `child_indices_extents(::BlockCursor)` (`cursor.jl:255-286`) still allocates and recomputes leaf caps per visit, where `MemoRasterTree` memoises them. | low–moderate | DGG `src/regridding.jl` + `src/systems/CopernicusDEM/cursor.jl` |
| **F2** | **Fix `_mergecaps` to be an actual bounding cap** (`conservative.jl:62-79`) — a minimal (or Ritter-style) enclosing cap instead of a mean-centre one. | ~4× fewer visited nodes for *every* system that falls back to `CellCapTree` (H3, ISEA, IGEO7 at non-chunk ranges). Measured: removes a 3–4× radius inflation. Does nothing once F1 lands for CopernicusDEM, everything for the others. | low | `lib/GlobalRegridding/src/conservative.jl` |
| **F3** | **Split `CellCapTree` spatially, not by index** (`_cellcapnode`, `conservative.jl:43-56`) — median-split the node's cap centres on their widest axis, permuting `ix`/`caps` (leaves keep global positions, so nothing downstream changes). | the rest of D1 for the generic fallback: turns the 11 729-node descent into something near the 2 743-node one | moderate | `lib/GlobalRegridding/src/conservative.jl` |
| **F4** | **A `MemoRasterTree` analogue for DGG cursors** — memoise `node_extent` and leaf `child_indices_extents` in task-local direct-mapped slots keyed on `(level, id)` instead of the index rectangle. Copy `raster_tree_memo.jl:26-88` almost verbatim. | 43.2 % of arm-1 block time *today*; ~5 % once F1 removes the surplus descents. Helps every DGG **destination**, which is every run. | low | DGG `src/engine/cursor.jl` (or a generic wrapper in `lib/GlobalRegridding`) |
| **F5** | **Cheap per-cell caps in the fallback.** `_cellcap` (`conservative.jl:82-84`) synthesises the full polygon and merges per-vertex caps; `Fallbacks.cell_cap` (`src/fallbacks/caps.jl:81`) is the tight-cap hook `cap_cached_tree.jl:93` already uses. Add a `GR.cellcap(space, i)` hook with a `DGGSpace` method. | at most the 8.1 % construction share (1.09×), and 0 once F1 lands | low | `lib/GlobalRegridding/src/conservative.jl` + DGG `src/regridding.jl` |
| **F6** | **Memoise the source subtree across block builds** — the `TileCells` analogue on the source side. Today `TileCells` wraps only the destination (`plans.jl:398`, `plans.jl:416`, `lazy.jl:415`) and the source tree is rebuilt per `(dst tile × src chunk)` pair (measured: 10 builds for 3 destination chunks over 4 tiles, so each tile's tree is built 2–3×). | ⅔ of the 8.1 %; ~0 once F1 makes the build O(1) | low | `lib/GlobalRegridding/src/plans.jl` |
| **F7** | **`CachedCellTree` has no `Trees.split_weight`** (`conservative.jl:197-224`) and its `ncells` forwards to the wrapped tree's, which for a `CellCapTree` is `ncells(tree.space)` — the whole space (`conservative.jl:99`). The frontier's default `split_weight` then reads a global cell count. Planning-only, no correctness impact, but it will over-split. | scheduling quality only | low | `lib/GlobalRegridding/src/conservative.jl` |

**Order: F1 → F4 → F2 → F3 → F6/F5/F7.** F1 is measured at 60.8× (§12.4) and is
the smallest change in the table; F4 is the general destination-side win
and is independent; F2/F3 are what make *every other* DGG system's source path
fast, and they are the ones to do if the fallback is ever meant to be the answer
rather than the fallback. Nothing here needs an upstream ConservativeRegridding
change: `intersection_areas` is reached identically by both arms
(`conservative.jl:391-395` → `ConservativeRegridding/.../intersection_areas.jl:249`),
and it is the trees handed to it that differ.

**Regression guard worth adding with F1:** the A/B in §12.1 is a two-line test —
the same destination fed a `DGGSpace` source and an equivalent `RasterGrid`
source must agree bit for bit (it already does) *and* be within a small constant
factor in `GR.subtree` node-visit count.

### 12.7 What this does to §11's arithmetic

`670 cells/s` was never the conservative kernel's speed; it was the cost of
`GR.subtree`'s choice of source tree. Applying the measured multiplier —
**44.5×** from the A/B, **60.8×** from the F1 prototype; the table uses the
conservative end of that range — to run B2's uncontended mid-latitude rate:

| configuration | cells/s | global IGEO7 level 12 (1.384 × 10¹¹ cells) |
|---|---|---|
| as measured today, `-t 8` (§11) | 670 | 6.5 years |
| **with F1, `-t 8`, one process** | **3.0 – 4.1 × 10⁴** | **39 – 53 days** |
| with F1 + §10 outer parallelism (`W = 21, T = 3`) | 6.3 – 8.6 × 10⁵ | **1.9 – 2.5 days** |
| `BilinearPoint`, `W = 64, T = 1` (§7, unchanged) | ≈ 2.5 × 10⁶ | ≈ 15 hours |

The middle row is corroborated independently by §12.5: a GLO-30 tile onto
level 13 runs at 89 843 cells/s (La = 7) / 130 090 cells/s (La = 6) through the
*raster* path on this machine, so tens of thousands of cells per second is the
right order for the fixed DGG path at GLO-90 / level 12 with tile-sized source
chunks.

So §11's conclusion — "correct at full scale and not viable at full scale" — was
right about the arithmetic and wrong about the cause. The cause is not "the
conservative weight build at 10⁴–10⁵ source pixels per destination cell"; it is
one dead `isa` check in `_chunkcursor` and a bounding cap that is 4× too big.
With F1 the conservative path lands within ~4× of the bilinear path *with the
same outer parallelism*, which is a different design conversation from the 177×
in §7.

**The tree-aligned chunking survives unchanged.** Both arms use the identical
`DGGSpace(dstgrid; chunklevel = La)` destination and the identical one-tile
source chunk; nothing about §1's three-way alignment, the ancestor-snapped
coverage query, or the Zarr chunk plan is implicated. The fix is entirely inside
`GR.subtree`'s choice of source tree.

---

## 13. Source-path fixes landed (2026-08-20)

§12's work plan is implemented on branch `claude/dgg-source-subtree` (worktree
`/home/asinghvi17/geo/DGG-source-fast`), three commits on top of `origin/main`
`92ef101` plus the cherry-picked SHA pins. F1, F4 and F2 all landed; F3, F5, F6
and F7 did not.

### 13.1 What landed

**F1 — `subcursor`, the windowed source tree** (`4d55dd4`). A new grid seam,
`DGG.subcursor(grid, inds) -> tree or nothing` (`src/interface/grid.jl`), whose
default is `nothing` and which `GR.subtree(::DGGSpace, inds)`
(`src/regridding.jl`) tries before `_chunkcursor` — dispatch only, so it costs
grids without a method nothing. `CopernicusDEM` implements it
(`src/systems/CopernicusDEM/cursor.jl`) by lifting `_block_cursor`'s
run-to-rectangle test off the whole grid and onto the window: `_window_cursor`
checks that the positions `inds` are one contiguous id run (`hi.index -
lo.index + 1 == length(inds)`) and `_run_cursor` then applies `treeify`'s own
three rectangle rules. `_block_cursor` is now `_window_cursor(grid, strategy,
1:ncells(grid))`, so there is one rule set, not two.

Two consequences worth naming, both of which fall out rather than being
designed in:

* **Any window of whole raster rows qualifies**, not just an exact chunk range.
  A range that is not one lattice rectangle — a mid-row window, a run spanning
  two tiles — still takes `GR.CellCapTree`, which is the documented behaviour.
* **A sparse holding works out of the box.** Copernicus ships land tiles only
  (~26 k of 64 800), so production will build the source over a `PartialGrid`.
  Its *whole-grid* tree is still the generic `HierarchicalGridCursor` — the id
  run is not contiguous across missing tiles — but every per-tile chunk of it
  is one contiguous run and windows normally. `DGGSpace(pg; chunklevel = 0)`
  produces exactly one chunk per held tile. Nothing else is needed.

**F4 — a chunk's cursor carries its own leaf caps** (`ecd4f34`).
`CapCachedTree` (`src/cap_cached_tree.jl`) existed only for the whole space;
`_chunkcursor` handed back a bare cursor whose `node_extent` at a leaf and
whose `child_indices_extents` both re-derive every cap through
`Fallbacks.cell_cap` → `cell_boundary`, once per opposing node — §12.2's 43 %
line. Give `CapCachedTree` an index offset (position `p` reads `caps[p -
offset]`, `0` for the whole space and `first(inds) - 1` for a chunk) and wrap
chunk cursors up to `_CHUNK_CAP_CACHE_MAX = 1 << 16` cells. On the lazy path
the destination reaches this through `GR.TileCells`, which holds one tree per
tile, so the caps are filled once and read by every source chunk paired with
that tile.

This is the eager-per-chunk shape rather than §12.6's suggested per-task
direct-mapped memo. It reuses a wrapper the package already has, it is exact
rather than best-effort, and `TileCells` already gives it the reuse a memo
would have to earn; the memo's advantage — bounded memory — is bought here by
the size cap instead. The memo remains the better answer if a destination tile
ever has to exceed 65 536 cells.

The cap lands in the same place §12.5's does: La = 7 (16 807 cells) keeps the
cache, La = 6 (117 649) loses it, exactly as it loses `TileCells`' polygon cache
at `_TILE_CELL_CACHE_MAX = 1 << 16`. If a future design moves to La = 6 for
Zarr-chunk reasons, both constants want raising together and both want
measuring at that size.

**F2 — `_mergecaps` is now a bounding cap** (`a30ca26`,
`lib/GlobalRegridding/src/conservative.jl`). Fold the smallest cap containing
two caps — the construction `GO.UnitSpherical._merge` and
`DGG.Fallbacks.merge_caps` already use — instead of centring on the mean of the
inputs' centres. §12.3's table, remeasured on the same tile (GLO-90, 45 °N,
1 440 000 cells), radius in km:

| depth | §12.3 median (max) | with F2, median (max) | `RasterCellTree` median |
|---|---|---|---|
| 0 | 274.9 | **94.6** | 68.1 |
| 3 | 224.2 (226.2) | **45.9 (46.2)** | 29.5 |
| 6 | 204.8 (220.0) | **39.9 (40.2)** | 8.5 |
| 8 | 194.8 (218.4) | **39.2 (39.6)** | 4.24 |
| 11 | 121.0 (203.9) | **38.8 (39.4)** | 1.82 |
| 13 | 5.76 (176.5) | 5.76 **(39.4)** | 0.93 |
| 17 | 0.39 (71.2) | 0.40 **(39.4)** | 0.24 |

The ~39 km floor that remains at every depth is exactly D1 and nothing else:
39.4 km is a tile's full degree of longitude at 45 °N, so a node that is a band
of complete pixel rows cannot be narrower. D2 is gone.

### 13.2 Measured

`-t 8 --gcthreads=8,1`, Julia 1.12.6, box load average 19–20 throughout, so
comparable to §12.5's conditions and not to §12.4's. Both arms run in one
process off one script, warm (second repetition), with the fallback arm forced
at the `GR.subtree` seam by a test-local wrapper — which, with F4 and F2 backed
out, reproduces §12.1's arm 1 to within 1 % (60.9 and 61.9 s on two runs against
61.42 s), so it is a faithful stand-in for `main`.

Region `box:10,10.05,45,45.05` (§12.1's A/B), IGEO7 level 12, ancestor 7,
50 421 cells, 3 destination chunks, **10 source subtree builds** in both arms —
the same count §12.1 measured:

| tree | wall | cells/s | vs arm 1 |
|---|---|---|---|
| `CellCapTree`, before any of this (= §12.1 arm 1) | 61.9 s | 815 | — |
| + F4 | 32.8 s | 1 536 | 1.9× |
| + F2 | **13.0 s** | 3 890 | **4.8×** |
| **windowed `BlockCursor` (F1), all three landed** | **1.39–1.50 s** | **33 600–36 200** | **~42×** |

Region `box:10,10.2,45,45.2` — the spike's **run D**, 218 491 cells, 13 chunks,
33 source subtree builds:

| | wall | cells/s |
|---|---|---|
| spike run D as recorded in §4 (contended, includes the Zarr write) | 234.2 s | 933 |
| forced `CellCapTree` on this branch (F4 + F2 already in it) | 87.0 s | 2 512 |
| **this branch** | **7.63–7.81 s** | **~28 600** |

The rate holds within 15 % across a 4.3× change in destination size (33 600
against 28 600 cells/s), so §12.7's arithmetic can use the branch's rate
directly rather than a ratio.

**Every arm's output is bit-identical to every other arm's** — `isequal` over
all 50 421 / 218 491 cells, zero non-finite in any of them. That is the bar
§12.4 set and it holds for F2 and F4 as well as F1: caps prune, they do not
weigh.

F4 and F2 also touch the *raster*-source path — F4 through the destination
chunk tree, F2 through any non-rectangular range — so §12.5's reference
configuration is the regression guard. `bench/harness.jl` `run_case(size =
3600, source = :chunked, srcchunk = (128, 128), chunklevel = 7, lazy = true)`,
the cached real GLO-30 tile onto its 16 490 405-cell level-13 MOC:

| | wall | cells/s |
|---|---|---|
| §12.5 as recorded | 183.6 s | 89 843 |
| this branch, two warm repetitions | 179.2 / 187.6 s | 92 021 / 87 879 |

Unchanged, with the same `nnz`.

### 13.3 Tests

* `test/systems/CopernicusDEM/runtests.jl`, new section (k2), 49 assertions:
  every source chunk of the complete level grid and of a scattered
  `PartialGrid` gets a `BlockCursor`; its leaves are the chunk's positions
  exactly; whole-row windows qualify and mid-row and cross-tile ones fall back;
  and for each chunk of the sparse holding the intersection-area matrix against
  a shared IGEO7 level-7 destination equals the `GR.CellCapTree` arm's, entry
  for entry. Whole file: 16 162 pass, 3 broken, 0 fail.
* `test/systems/crosssystem/regrid.jl`: the chunk tree is now a
  `CapCachedTree`; its caps, leaf entries and cells are `_chunkcursor`'s bit
  for bit (`checktree`), its cap vector is chunk-length and its offset is
  `first(inds) - 1`. 34 pass in that testset.
* `lib/GlobalRegridding/test/test_conservative.jl`: a new testset asserts
  `_mergecaps` both **covers** (every node's cap contains every vertex of every
  cell below it) and is **tight** (root radius within 1.3× of half the widest
  vertex separation, which Jung's bound puts the true minimum below 1.16× of).
  Whole suite: 2 202 pass, 1 broken, 0 fail.
* Full `test/runtests.jl` at `-t 4`, Julia 1.12.6: **987 238 pass, 17 broken,
  0 fail**, 12 m 54 s. The broken count is `main`'s; the pass count is
  `main`'s 987 162 plus the 76 assertions added here.

### 13.4 Not done

Referenced to §12.6, in the order they are worth doing:

* **F3 — split `CellCapTree` spatially.** The whole of what is left of D1: the
  ~39 km floor in §13.1's table. A median split of the node's cap centres on
  their widest axis; nothing downstream reads `ix` in order, so permuting `ix`
  and `caps` together suffices. A comment at `_cellcapnode` now says so. Worth
  the remaining ~4× for every system that has no windowed cursor — which is
  every system except Copernicus.
* **F6 — memoise the source subtree across block builds.** Still 33 builds for
  13 destination chunks in run D. Now that a build is O(1) this is scheduling
  polish, not time.
* **F5 — a `GR.cellcap(space, i)` hook.** Zero once F1 lands for Copernicus;
  it is `CellCapTree`'s construction cost for everyone else, and F3 is the
  bigger half of that.
* **F7 — `CachedCellTree` has no `Trees.split_weight`.** Unchanged, and still
  planning-only.

One thing §12 did not list: `Trees.ncells` of a chunk tree answers for the
chunk while its leaf indices are global, so a chunk tree cannot be fed to
`CR.intersection_areas` directly — only through a block build's index maps.
That is why §13.3's chunk-level comparison is `checktree` rather than a weight
matrix.

## 14. Ancestor-subzone 2-D store layout (2026-08-20)

Built on branch `claude/subzone-store` (worktree
`/home/asinghvi17/geo/DGG-subzone-store`), three commits on top of
`origin/claude/dgg-source-subtree` `660640f`:

| commit | |
|---|---|
| `39b1094` | the layout's arithmetic and vocabulary (`src/io/subzones.jl`, 624 lines) |
| `05f8af0` | the Zarr half: store, incremental writer, lazy view (`ext/DiscreteGlobalGridsZarrExt/subzones.jl`, 626 lines) |
| `c784f23` | the store-IO reference section |
| `86f3174` | one more test: the stack form of the one-shot write |

### 14.1 Why a second layout at all

§3 and §5 are about a one-dimensional cell axis cut into equal chunks. Zarr v2
chunks are uniform by format, so a chunk grid that follows the TREE is not
expressible that way: a level-La pentagon's subtree holds
`p(d) = (5·7^d + 1)/6` cells and a hexagon's holds `7^d`, and no single chunk
length lands on both. Every consequence in §3 — `_rawids` materializing a dense
id vector, the ranges coordinate, the chunk manifest — is a consequence of
keeping the axis in one dimension.

The ancestor-subzone layout spends a dimension on the tree instead:

```
dim 1 (fastest, Julia dim 1)   subzone position 1:capacity within one subtree
dim 2                          the ancestor: the i-th level-La cell, ascending
chunks                         (capacity, 1) — one chunk IS one subtree
fill_value                     NaN
```

Three things follow, and they are the whole reason for it:

* **Sparsity is free.** A land-only global store writes only the columns it
  covers; an unwritten ancestor is a chunk that was never stored and reads back
  as NaN at zero bytes. No mask, no index, no compaction.
* **A column is one file.** Production writes columns from as many tasks as it
  likes with no coordination, because nothing shared is rewritten per column —
  not the attrs, not `.zmetadata`, not a manifest. The suite asserts exactly
  that (one new file per column write, every pre-existing file untouched to the
  mtime).
* **The reader gets the tree's chunking back.** The cell-axis view publishes
  `DiskArrays.IrregularChunks` of the real column lengths — `7^d` for a hexagon
  ancestor, `p(d)` for the twelve pentagons — which is the chunking Zarr's own
  grid could not hold.

No level-L id vector is ever materialized: the mapping both ways is
`ancestor`/`cellposition`/`descendant_range`, O(level) digit arithmetic, so the
`_rawids` finding of §3 simply does not apply to this path. It is also
system-generic rather than Z7-specific — it works on any system with
`has_sorted_subtrees` (HEALPix included; A5 is refused by name).

### 14.2 Public API

```julia
# one shot
dggwrite(dest, cube; layout = :subzones, ancestor_level = 6,
         capacity = nothing, fill_value = NaN, ancestor_coordinate = true,
         compressor = Zarr.BloscCompressor())            -> dest

# incremental / parallel: the production shape
store = subzonestore(dest, system, level; ancestor_level, layers,
                     capacity = nothing, fill_value = NaN,
                     ancestor_coordinate = true, attrs = Dict{String,Any}(),
                     compressor = Zarr.BloscCompressor())   # create
store = subzonestore(dest)                                  # reopen writeable
dggwrite!(store, ancestor_cell_or_column, values; var = the only layer)
dggwrite!(store, ancestor_cell_or_column, (elevation = v1, slope = v2))
dggwrite!(store, cube)          # every complete column of a cube
keys(store)                     # layer names

# read
dggread(path)                          -> DimStack over the WHOLE level
dggread(path; ancestors = cells_or_columns)   -> only those columns
dggread(path, :elevation)              -> DimArray
```

`layers` is `name => eltype` pairs (`("elevation" => Float32,)`, a `Dict`, a
`NamedTuple`). `values` is exactly as long as that ancestor's subtree really
is — `columnlength(layout, column)` — and a pentagon column takes `p(d)`, not
`7^d`; passing `7^d` to a pentagon is an `ArgumentError` rather than a silent
overwrite of the padding.

The layout descriptor and its arithmetic are exported/`public` from the package
proper and need no Zarr:

```julia
SubzoneLayout(system, level, ancestor_level; gridname, capacity)
DGG.subzone_capacity(system, ancestor_level, level)   # measured longest subtree
DGG.columncell(l, i) / DGG.columnindex(l, ancestor)
DGG.columnpositions(l, i) / DGG.columnlength(l, i)
DGG.subzoneindex(l, cell) -> (column, row)  |  DGG.positionindex(l, p)
DGG.subzone_runs(l, cells) -> Vector{SubzoneRun}    # cube axis -> columns
DGG.subzone_cellvector(l, columns)                  # columns -> cell axis
DGG.subzone_attrs(l; ...) / DGG.subzone_layout(attrs) / DGG.issubzonestore(attrs)
```

### 14.3 On-disk vocabulary

One group attribute, `dggs`, with the spec's own keys honest and everything the
layout adds nested under `subzone_layout`:

```json
{"dggs": {"name": "igeo7", "refinement_level": 12, "indexing_scheme": "z7int",
  "subzone_layout": {
    "layout": "ancestor_subzone", "version": 1, "writer": "DiscreteGlobalGrids.jl",
    "ancestor_level": 6, "ancestor_dimension": "ancestor",
    "ancestor_count": 1176492, "ancestor_order": "ascending_id",
    "ancestor_coordinate": "ancestor_cell_ids",
    "subzone_dimension": "subzone", "subzone_count": 117649,
    "subzone_order": "ascending_id",
    "padding": "trailing_fill", "padding_fill_value": "NaN",
    "chunk_shape": [1, 117649], "variables": ["elevation"]}}}
```

Deliberately **no `zarr_conventions` declaration**: this is not the
one-dimensional layout `zarr-conventions/dggs` describes, and declaring it would
send a convention-aware reader down a path that cannot open the store. A foreign
DGGS reader gets "no convention detected", which is true. Which twelve columns
are pentagons is derivable from the grid and is not recorded. `subzone_order`
is OGC API-DGGS's term for position within a parent zone.

`ancestor_cell_ids` (optional, on by default) is the level-La ids in order,
carrying `grid_name`/`level` — so xdggs sees the ancestor axis for the level-La
cell axis it really is, and a store read by anything but this package still
resolves its columns. This package never reads it: the column axis is implicit,
position IS the ancestor.

The declared `ancestor_count`/`subzone_count` are CHECKED against the grid's own
arithmetic on read, not believed — a store written against another grid
definition would otherwise read every column at the wrong offset, silently.

### 14.4 Pentagons, verified

`p(d) = (5·7^d + 1)/6`, confirmed empirically through `descendant_range` at
La = 2, L = 4: exactly 12 of the 492 columns are short, all of them 41 cells
against a capacity of 49, and they are exactly the digit-0 chains from the
twelve root cells (`first(children(sys, c))` iterated La times). The tests
assert both the count and the identity, and the store test asserts the padding:
41 real values then 8 NaN, with the axis that comes back holding 41 cells.

### 14.5 The read side: a 2-D store faked into a 1-D cell axis

`dggread` recognizes the layout from its own attribute BEFORE the conventions
are asked (its data arrays share no cell dimension, so there is nothing for a
`StoreDescription` to describe) and hands back an ordinary `DimStack`: one
`Cells` dimension carrying a real `CellLookup`, one lazy `SubzoneCellArray` per
layer. The wrapper is a `DiskArrays.AbstractDiskArray{T,1}`:

* `size` is the cell axis length; `readblock!` walks columns, reading
  `z[rows, col]` per column, so a read inside one column is one chunk read and
  one spanning columns is one per column, in order, padding dropped.
* `haschunks` is `Chunked()` and `eachchunk` is
  `GridChunks(IrregularChunks(chunksizes = column lengths))`, built on demand
  and memoized — a whole-level view at La = 6 would be 1 176 492 integers, and a
  read that never asks how the store is chunked never builds them.
* The default view is the WHOLE level (position in the view == position in the
  level grid, so no tables at all: constant memory whatever the level).
  `ancestors = [...]` restricts it, and then the offsets are one vector over the
  selected columns.
* Writing through the view is refused by name, not by `MethodError`: a column is
  written whole, padding rule included.

Stack metadata carries `"source"`, `"attrs"` and `"layout"` (a `SubzoneLayout`);
no `"description"`/`"conventions"`/`"encoding"`, none of which have anything to
say about a two-dimensional store.

### 14.6 Deliberately not built

* **No multi-dimensional layers.** A time or band axis would want a third
  dimension; the layout spends both of its own on the cell axis. A 2-D cube is
  refused with a message saying to use the one-dimensional layout or one store
  per step.
* **No written-column discovery.** There is no index of which columns exist and
  no listing of chunk keys; the whole-level view simply reads unwritten columns
  as fill, which is Zarr's own semantics and costs nothing. A production script
  that wants the covered set should keep it (it knows it — it wrote it).
* **No chunk cache in the wrapper.** Repeated reads of one column re-read its
  chunk; `DiskArrays.cache` is the answer where one is wanted.
* **No `CellEncoding` registration.** The layout is not an encoding: encodings
  choose the shape of the cell COORDINATE, and this chooses the shape of the
  store. It is `layout = :subzones` on `dggwrite`, orthogonal to `encoding`.
* **No remote writes** (same rule as `dggwrite`), and no `description` escape
  hatch on read: the shape is in the store's own attributes.

### 14.7 Tests

`test/io/subzones.jl` (110 assertions, no Zarr) and `test/io/subzone_store.jl`
(86 assertions, self-skipping on Zarr), both wired into `test/io/runtests.jl`.
The whole io suite passes with them in: **937 assertions, 0 failures**, of which
196 are new, and the package's whole suite is green on the branch:
**987 456 pass / 17 broken / 0 fail** (the same 17 broken as before). Covered: a
mid-latitude hexagon batch, a pentagon column, unwritten columns reading NaN,
incremental == one-shot (array bytes AND group attrs), the incomplete-subtree
refusal at three shapes, the irregular chunk grid, a column write touching one
file, several layers, the error surface, a dense store re-written as subzones
through a `ChunkedCellLookup`, and an integration test that regrids a synthetic
15° global raster onto three complete level-2 subtrees, writes them, and reads
back values identical to the direct regrid output.

### 14.8 For the production script

* **`capacity` is measured, once, at store creation**: one `descendant_range`
  per level-La cell, ~1.2 M of them at La = 6 (a fraction of a second). Pass
  `capacity = 7^d` explicitly to skip it if a run creates many stores.
* **Chunk size is `capacity × sizeof(T)`**, and that is the only thing the
  layout constrains La by: GLO-90 at L12/La6 is `d = 6`, capacity 117 649,
  470 KB Float32 chunks over 1 176 492 columns. GLO-30 at L13 is `d = 7` at
  La = 6 — 3.3 MB chunks over the same 1.18 M columns — or `d = 6` at La = 7,
  470 KB chunks over 8 235 432 columns and therefore that many files. See
  `2026-08-20-la-choice.md` for the regridding side of the same choice; the
  store side has no opinion beyond those two numbers.
* **Writing columns in parallel is safe on a `DirectoryStore`** as long as they
  are disjoint. Create the store once from the driver, hand the same
  `SubzoneStore` to every task, `dggwrite!(store, ancestor, values)`.
* **`ancestor_coordinate = true` writes 1.2 M UInt64 (9.4 MB)** at La = 6 —
  once, at creation. Pass `false` where that matters; nothing in this package
  reads it.
* **Coverage must be ancestor-snapped.** A cube covering part of a subtree is
  refused (`DGGSFormatError(check = :incomplete_subtree)`), which is what
  ensures no data cell is indistinguishable from padding. Regridding onto
  `subzone_cellvector(layout, columns)` — or any covering expanded from level-La
  cells — satisfies it by construction.
* **Environment note.** The committed `Manifest.toml` on
  `claude/dgg-source-subtree` records a `git-tree-sha1` for
  ConservativeRegridding that does not belong to the pinned rev `6a4b997`
  (it is the tree of `89ec5ca`, which has no `Trees.split_weight`, so
  GlobalRegridding fails to precompile). `Manifest.toml` is gitignored, so the
  fix is local: set `git-tree-sha1 = "873cc732e77638eb44d7d92e06880585a31200a9"`
  beside that `repo-rev` and re-instantiate.

## 15. Production dress rehearsal, synthetic (2026-08-20)

The whole land surface, end to end, from one script into one store. Branch
`claude/copdem-production` (worktree `/home/asinghvi17/geo/DGG-subzone-store`),
two commits on top of `claude/subzone-store` `86f3174`:

| commit | |
|---|---|
| `648aca6` | `_chunkwindows` narrows a ROOTED `PartialGrid` to its own root |
| `47961ea` | `scripts/copdem_production.jl` |

No PR.

### 15.1 What is real and what is synthetic

The DEM is never downloaded. What IS downloaded, and is the whole point, is the
**tile list** — `https://copernicus-dem-90m.s3.amazonaws.com/tileList.txt`,
26 475 names, cached at `bench/data/CopernicusDEM/tileList-glo90.txt`
(gitignored, ~900 KB). Copernicus ships land tiles only: 26 475 of the 64 800 a
1x1-degree lattice has, **40.9 %**. A tile off that list does not exist, so the
source is a `PartialGrid` over the listed tiles' pixels and an over-covered
open-ocean destination pairs with no source chunk at all.

Two further real inputs, both small:

* the four GLO-90 GeoTIFFs already in `bench/data/CopernicusDEM/90m/` (S85/E000,
  S89/E000, S90/E000, S90/E001), which decode for real;
* Natural Earth 10m land (`ne_10m_land.shp`, 3.3 MB, cached under
  `bench/data/naturalearth/`), rasterised into the nodata mask below.

Everything else is `SYNTHETIC(lon, lat) = 1000 sin(3L) cos(2P) + 500 cos(7L)
sin(5P) + 100`, the spike's field, which is smooth on the scale of a level-12
cell by four orders of magnitude and is therefore its own oracle.

### 15.2 Ocean pixels inside listed tiles are NODATA

A listed tile is a *land* tile, not an all-land tile: it is a degree square with
land somewhere in it. Making every pixel of it valid would have skipped the
missing-data machinery entirely, which at this scale is most of what is worth
rehearsing. So the run rasterises the Natural Earth coastlines once into a global
Bool lattice and sets every ocean post of a synthesised tile to `NaN32` — the
regridder's own invalid sentinel (`GR._isvalid(x::AbstractFloat) = !isnan(x)`),
so a synthetic tile's nodata is indistinguishable in kind from a real one's.

The rasteriser is a scanline even-odd fill with an active edge list
(`rasterize_land`), edges bucketed by the first row they reach: **6838 rings /
446 175 vertices onto 86 400 x 43 200 at 15 arcsec in 0.5 s**, 445 MiB as a
`BitMatrix`, one copy shared by all workers. Even-odd handles interior rings as
the holes they are, and Natural Earth is clipped to [-180, 180] so no edge wraps
the antimeridian. It is cheap enough that the run rebuilds it every start rather
than caching it. Land is **33.09 %** of the lattice by cell — plate-carree, so
Antarctica is over-weighted against the 29 % of area.

Measured consequence, from the smoke tests: **57 % of the pixels inside listed
tiles are valid** in a coastal region, and destination columns come out anywhere
from 0 % to 100 % NaN.

### 15.3 Configuration

Everything below is `bench`'s environment on the pinned CR `6a4b997` / GO
`2825c47`, Julia 1.12.6.

| | |
|---|---|
| source | GLO-90 (N = 1200), 26 475 listed tiles, 2.513e10 pixels |
| source space | `DGGSpace(PartialGrid(TileIds); chunklevel = 0)`, 26 475 chunks, **built in 0.1 s** |
| destination | IGEO7 level 12, ancestor level **La = 5** |
| work unit | one level-5 column = `7^7` = 823 543 cells = one Zarr chunk = one file |
| columns | **66 178** cover the tiles (of 168 072 global), computed in 10.8 s |
| total | **5.4500e10** level-12 cells, 203.0 GiB dense f32 |
| method | `Conservative()`, `Weighted(0.5)`, `lazy = true`, `budget = 2^30`; `chunks` NEVER passed |
| store | `/home/asinghvi17/geo/dggstores/copdem90-igeo7-l12-synthetic.zarr`, subzone layout, `elevation` Float32, fill NaN, `capacity = 7^7` |
| parallelism | W = 21 workers, `julia -t 63 --gcthreads=8,1`, batch = 8 columns |
| box | 64 cores / 751 GB, SHARED; 1-min load 19.9 at launch (the brief's drop-to-W=14 threshold is 30) |

`La = 5` is the la-choice verdict and nothing here revisits it.

### 15.4 The lazy source id vector, and why the sparse holding is free

`PartialGrid` stores its ids **by reference** and never materialises them, and
`subtree`'s `SubtreeIds` already shows the shape a lazy one takes. `TileIds` is
the multi-range version: an offsets table of 26 475 integers, `getindex` a
binary search, and `Helpers.strictly_increasing` answered `true` by construction
rather than by the O(n) scan the constructor would otherwise run over 2.5e10
elements.

`_chunkwindows` at `chunklevel = 0` then visits the 64 800 level-0 ancestors,
binary-searching that vector twice each, and the tiles that are not held produce
empty windows and are skipped. **0.1 s for the whole global source space.**
Section 13.1 predicted this would work out of the box; it does, at full scale.

### 15.5 The destination-side fix: a rooted subset should not scan its level

A work unit's destination is `subtree(sys7, column, 12)` — a `PartialGrid`
ROOTED at the column's own level-5 cell. `DGGSpace(...; chunklevel = 5)` then ran
`_chunkwindows` over all **168 072** level-5 ancestors to discover that exactly
one of them is non-empty: the la-choice note measured that at **0.30 s per
space**, region-independent, and this run builds 66 178 of them — ~5.5 CPU-hours
spent rediscovering the answer it was handed.

`648aca6` narrows the visit: a rooted `PartialGrid` holds nothing outside its
root's subtree, so only that root's level-`a` descendants can qualify, and a root
DEEPER than the chunk level puts the whole grid under one ancestor found by
arithmetic. Exact, not a heuristic — `PartialGrid`'s constructor checks the
endpoint ancestry. The test asserts the narrow and the wide answers agree
ancestor for ancestor and range for range at every chunk level between the root
and the grid; `test/systems/crosssystem/regrid.jl` is +17 assertions and green.

### 15.6 Shape of the run

Per section 10, with one deviation named below.

* **W x T in process.** 21 worker tasks in one Julia process, `-t 63`. There is
  no explicit `T`: `_wavesize` reads `Threads.nthreads()`, so a worker's inner
  wave is `min(63, source chunks meeting this destination tile, budget fit)`,
  which at mid-latitude is the 2-4 the measurements have always shown and at the
  pole is the 200-360 that is the one place the inner parallelism has material.
  `threadsper` is therefore config that is REPORTED, not enforced; `-t W*T` is
  what sets the cap.
* **`GR.OUTER_PARALLEL` is set inside each worker body.** It disables
  `_innerthreaded` — CR's per-weight-build threading — and does NOT disable the
  lazy plan's own source-chunk wave, which spawns without consulting the scoped
  value. That is exactly the split section 10 wanted.
* **Deviation: batches are PULLED, not dealt.** The brief says workers take
  contiguous column groups; a static W-way split of the sorted column list is
  contiguous but badly balanced, because la-choice measured a polar column at
  ~5x a mid-latitude one and Antarctica is a contiguous run. So the columns are
  cut into contiguous batches of 8 and workers pull batches off one atomic
  cursor. A batch is still geographically connected — Z7 ids sort into subtrees,
  so neighbouring column indices are neighbours on the sphere and a worker's tile
  cache stays hot within a batch — but the tail cannot set the wall clock.
* **The tile cache is striped.** One `TiledDEM` is shared by all 21 workers, so a
  single lock over one LRU would serialise them for the ~15 ms a 1200x1200
  synthesis takes. 64 stripes, 64 locks, and the tile is built OUTSIDE its lock;
  a race builds it twice and drops the loser.
* **Skip pruning is structural.** The column list is built FROM the tiles
  (`covering_columns`), so every enqueued column meets at least one listed tile
  and a column with no source is never enqueued at all.

### 15.7 Resume

`--resume` is the default and takes the **union** of two sources, because
neither is complete on its own:

* the append-only done log `<store>.done.ndjson`, one flushed line per column;
* the store's own chunk listing, `<store>/elevation/<column-1>.0`.

The log can be lost or truncated. And **Zarr stores no chunk for a column whose
every value is the fill value** — which on the ocean side of the covering is a
normal outcome, measured at 1 of 10 columns in the first smoke test — so an
all-NaN column leaves no file behind and is indistinguishable on disk from one
nobody computed. Taking the union recomputes nothing and leaves no hole; the
disagreement is reported rather than silently trusted.

### 15.8 Smoke tests

Two, both with `checks=true`, both green after one fix.

**A. Coastal Adriatic + Antarctic, 10 columns, W=8 on `-t 24`.** 8.235e6 cells in
50 s, 160 589 cells/s aggregate, peak RSS 5.4 GiB, 64.3 % of written cells NaN.
Per-column walls 3.4-31.7 s and per-column NaN fractions 0.0 % to 100.0 % — the
spread the mask exists to produce. Verified on four columns read back from disk:
every value finite or NaN; unfed cells NaN; fully-fed cells within
**1.7e-2 m** of `SYNTHETIC` at the cell centre (RMS 2.2e-3 m); coastal cells
inside the field's own range across their own cell; a column nobody wrote reads
back 686 286 of 686 286 NaN.

**B. Antarctic only, 18 columns, W=6 on `-t 16`,** to reach the real tiles. All
four decode and all four are finite at their tile centres, **both S90 pole tiles
included** — 2761.96 m at S90/E000 and 2766.74 m at S90/E001, against the S85
tile's real range of 2424-2580 m. 100 % of pixels valid, as an all-Antarctica
region should be.

**The fix the smoke found** was in the ORACLE, not the pipeline. The first
version asserted "the mask says land, therefore the cell is finite", which fails
wherever Natural Earth calls a point land that Copernicus never listed a tile for
— a coastline the two datasets disagree about, or, in a regional smoke, a tile
the region filter cut. 229 of 401 sampled cells in one column, 370 of 401 in
another. The predicate the oracle actually needs is `SourceMask`: a point has a
valid source pixel when its tile is on the list AND the mask calls the post land.
With that, both smoke tests are clean.

A second, smaller correction: "some sampled cell straddles the coast" is a claim
about a GLOBAL run and is asserted only there; an all-Antarctica box has no coast
and nothing to renormalise, and saying so is not the same as failing.

### 15.9 Launch

```
2026-08-20T18:32:30Z   launched, PID 2915199 (in /home/asinghvi17/geo/dggstores/copdem90-run.pid)
2026-08-20T18:32:51Z   21 workers up, 66 178 columns to do
```

```bash
cd /home/asinghvi17/geo/DGG-subzone-store
nohup env RASTERDATASOURCES_PATH=/home/asinghvi17/geo/DiscreteGlobalGrids.jl/bench/data \
  julia --project=bench -t 63 --gcthreads=8,1 scripts/copdem_production.jl \
  config=/home/asinghvi17/geo/dggstores/copdem90-run.conf \
  >> /home/asinghvi17/geo/dggstores/copdem90-run.log 2>&1 &
```

The config file is `/home/asinghvi17/geo/dggstores/copdem90-run.conf` and is the
whole of section 15.3; the run log is
`/home/asinghvi17/geo/dggstores/copdem90-run.log`.

**First five minutes, watched.** Startup to the first column was 14 s (package
load, then the first regrid carrying compile). Nothing else on the startup path
cost anything measurable: the source space is 0.1 s, the column list came off its
cache, and the store reopened in 3 s.

| elapsed (worker time) | columns | aggregate cells/s | RSS | process CPU |
|---:|---:|---:|---:|---:|
| 62 s | 32 | 4.25e5 | 10.3 GiB | 1683 % |
| 92 s | 56 | 5.01e5 | 15.4 GiB | 1889 % |
| 180 s | 122 | 5.58e5 | 11.9 GiB | 2080 % |
| 254 s | 173 | 5.61e5 | 19.8 GiB | 2155 % |
| 290 s | 199 | 5.65e5 | 20.2 GiB | 2187 % |

and the run's own heartbeats:

```
18:37:52  HEARTBEAT  210/66178 columns (0 skipped) | 1.728e+08 cells | 575412 cells/s
                     | elapsed 0.08 h | ETA 26.21 h | RSS 19.7 GiB | NaN 38.01%
18:42:52  HEARTBEAT  427/66178 columns (0 skipped) | 3.515e+08 cells | 585498 cells/s
                     | elapsed 0.17 h | ETA 25.68 h | RSS 22.1 GiB | NaN 26.25%
18:47:54  HEARTBEAT  644/66178 columns (0 skipped) | 5.302e+08 cells | 587423 cells/s
                     | elapsed 0.25 h | ETA 25.51 h | RSS 24.1 GiB | NaN 23.91%
18:52:55  HEARTBEAT  854/66178 columns (0 skipped) | 7.032e+08 cells | 584236 cells/s
                     | elapsed 0.33 h | ETA 25.57 h | RSS 25.8 GiB | NaN 18.95%
```

Across the first twenty minutes the rate rises and then settles
(575 412 -> 585 498 -> 587 423 -> 584 236 cells/s) and the ETA converges from
below and holds (26.21 -> 25.68 -> 25.51 -> 25.57 h) as the compile of the first
columns washes out of the average. Zero errors, zero skips, and the
columns-per-heartbeat is steady at 217 / 217 / 210.

**Settled figure: ~5.85e5 cells/s aggregate, ~25.5 h, completing about
2026-08-21T20:00Z.**

RSS climbs slowly (19.7 -> 22.1 -> 24.1 -> 25.8 GiB) as the 3072-entry tile cache
fills, and should plateau near 17.7 GiB of cache plus the workers' budgets —
comfortably inside the 214 GiB the box had free at launch. If it ever approaches
that, `cache=1024` on a resume is the knob.

The NaN fraction falling monotonically (38.0 -> 26.3 -> 23.9 -> 19.0 %) is the
covering leaving the first base cell: the earliest columns are the ocean-heavy
edge of the holding, and the figure should settle toward the whole-run number as
the workers spread. It is worth re-reading at the end, because the whole-run
value is the one that says how much of the land-tile covering is actually sea —
the number the mask and the tile list exist to produce, and the one a run without
either would have reported as zero.

Zero errors. Store 206 MB after 199 columns. **38.0 % of written cells are NaN**
so far, which is the mask and the tile list doing their job — a figure that would
be 0 % without either.

**Projection: 5.75e5 cells/s aggregate over 5.4500e10 cells = ~26 hours**, i.e.
completion around **2026-08-21T20:45Z**. That is 7x the brief's "~78k cells/s
order of magnitude" and 7x the la-choice single-process figure, which is what W
was for. It is also ~3x SLOWER than la-choice's `W=21 x T=3 -> 9.0 h`
extrapolation, and the reason is visible in the CPU column.

**Disk: ~1.03 MB per column measured, so ~68 GB for the store** against 203 GiB
dense f32 — a 3x Blosc ratio, better than the spike's 59-75 % because the ocean
side of the covering is long runs of NaN. 2.5 TB free, so it is not a
constraint.

### 15.9.1 The utilisation gap, named rather than fixed

The process holds **~21.9 cores of 64**, i.e. **~1.04 cores per worker**, where
section 10.4 projected `21 x min(1.85, 3) ~ 39`. The run is producing and
correct, so this is recorded rather than chased, but three candidates are worth
naming for whoever picks it up:

1. **The early sample is one neighbourhood.** All 21 workers pull from the front
   of the same sorted column list, so at five minutes every worker is inside base
   cell 0 within a few hundred column indices of the others. Their inner waves
   are therefore competing for the same handful of source tiles — good for the
   cache, and possibly bad for the wave width, since a wave is bounded by the
   number of source chunks meeting the destination. The first heartbeats past the
   base-cell boundary will say whether the figure moves.
2. **Source discovery now has 26 475 chunks, not 24.** Every column tests itself
   against the whole listed holding. The smoke tests had a source space three
   orders of magnitude smaller and cannot have seen this. A 100 %-NaN column
   costing 9-35 s of wall is the shape to be suspicious of: it computes nothing
   and still pays.
3. **The box is shared.** 1-min load was 19.9 at launch and 51 at five minutes,
   of which ~25 is not this run.

Note that `_wavesize` reads `Threads.nthreads()`, which is 63 for the whole
process rather than a per-worker `T`. That is the section-10 design as the code
actually stands, and it means `T` is emergent, not enforced: a worker whose
destination meets 3 source tiles gets a wave of 3 whatever `threadsper` says.

### 15.9.2 The live store, read while it is being written

The concurrent-write claim of section 14.8 is checked on the PRODUCTION store
rather than only on a test one: with 21 workers writing, an independent process
opened the store read-only and pulled one column back.

```
column 33 (ancestor Z7Cell("0000055")): 823543 values, 787155 finite, 36388 NaN
  |value - SYNTHETIC(centre)| over 501 sampled finite cells:
      max 7.83e-3 m, RMS 2.89e-3 m
```

So the store is readable, correct, and partial all at once — a column that has
been written reads its values, a column that has not reads `NaN`, and neither
the writers nor the reader needed to know about the other. `dggread` is the
right entry point for this; `subzonestore(path)` is NOT, because it opens the
group writeable.


### 15.10 Monitoring, resuming, reading

```bash
# progress
tail -f  /home/asinghvi17/geo/dggstores/copdem90-run.log
grep HEARTBEAT /home/asinghvi17/geo/dggstores/copdem90-run.log | tail -5
wc -l    /home/asinghvi17/geo/dggstores/copdem90-igeo7-l12-synthetic.zarr.done.ndjson
du -sh   /home/asinghvi17/geo/dggstores/copdem90-igeo7-l12-synthetic.zarr

# health
ps -o pid,etime,rss,%cpu -p $(cat /home/asinghvi17/geo/dggstores/copdem90-run.pid)
grep -c ERROR /home/asinghvi17/geo/dggstores/copdem90-run.log

# resume after a kill: the SAME command. `resume=true` is the default and the
# union of the done log and the chunk listing is skipped.
cd /home/asinghvi17/geo/DGG-subzone-store
nohup env RASTERDATASOURCES_PATH=/home/asinghvi17/geo/DiscreteGlobalGrids.jl/bench/data \
  julia --project=bench -t 63 --gcthreads=8,1 scripts/copdem_production.jl \
  config=/home/asinghvi17/geo/dggstores/copdem90-run.conf \
  >> /home/asinghvi17/geo/dggstores/copdem90-run.log 2>&1 &

# verify when it finishes (reads columns back off disk, computes nothing)
julia --project=bench -t 8 scripts/copdem_production.jl \
  config=/home/asinghvi17/geo/dggstores/copdem90-run.conf checks=true checkcolumns=24

# read
julia --project=bench -e 'import DiscreteGlobalGrids as DGG, Zarr
  st = DGG.dggread("/home/asinghvi17/geo/dggstores/copdem90-igeo7-l12-synthetic.zarr")'
```

If the machine's other tenants grow, kill the PID, re-launch with `workers=14`
appended; resume makes the restart cost one partial column.
