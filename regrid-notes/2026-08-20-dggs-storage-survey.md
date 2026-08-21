# Full-globe hexagonal DGGS datasets at scale — storage survey (2026-08-20)

> Opus research-agent deliverable for the CopDEM→IGeo7 Zarr design. Companion to `2026-08-20-copdem-zarr-spike.md` (esp. §on pentagon alignment and the global implicit-axis store).

Research complete. Here is the deliverable.

---

# Storage of full-globe DGGS cell-indexed data at scale — survey

*Scope: what has actually been published, verified against primary sources (specs, store metadata, S3 listings) wherever possible. Numbers I computed myself are marked **[calc]**.*

## Direct answers first

**(a) Has anyone published a ≥1e10-cell dense hexagonal-DGGS array?**
**No.** The largest published *hexagonal* DGGS dataset is ~1e10 cells and it is **tabular Parquet with an explicit id column**, not an array. The largest published *dense array* on any spherical DGGS is HEALPix — and it does exceed 1e10 (see NICAM below), so HEALPix is not merely "the only family with a track record", it is the only family where dense ≥1e10 has shipped at all.

**(b) Closest existing solution to the pentagon-alignment problem?**
**H3 pads the digit space** — the canonical prior art. It reserves all 7 digit values everywhere and simply declares pentagon digit-1 subtrees invalid ("deleted subsequence"). Nobody chunks per-base-cell into 12 arrays. The formal alternative is **OGC API-DGGS "sub-zone order"**: data for one parent zone at relative depth *d* delivered as a 1-D array whose length is allowed to differ for pentagons vs hexagons — i.e. ragged per-parent shards, standardized.

**(c) Does xdggs / DGGS-Zarr assume a materialized id coordinate, and has the scale problem been flagged?**
**xdggs: yes historically, and it was flagged loudly with numbers.** The new **`zarr-conventions/dggs`** convention: **no** — an omitted id coordinate is a first-class mode.

---

## 1. The closest existing thing to your project (verified)

**[walkthru.earth DEM-Terrain](https://source.coop/walkthru-earth/dem-terrain)** — a global DEM (GEDTM30, which itself fuses Copernicus DEM + ALOS + ICESat-2/GEDI) regridded onto **H3** (aperture-7 hexagonal icosahedral — same family as IGeo7).

| | |
|---|---|
| Grid / levels | H3 res 1–10; **8.96e9 cells at res 10**, ~10.45e9 total |
| Format | Parquet, ZSTD-3, 1,000,000-row row groups |
| Id axis | **explicit `h3_index BIGINT` column**, table **sorted by it** (for delta encoding + range queries) |
| Partitioning | Hive by `h3_res=` only — **one file per resolution**, no spatial sharding |
| Size (I listed the bucket) | res 10 = **167.94 GB in a single Parquet file**; total v2 = **195.68 GB** |
| Density | **18.7 bytes/cell** for int64 id + 5 float32 values **[calc]** |
| Notable | v1→v2 **dropped `geometry`/`lat`/`lon` because they are derivable from the id** — 262.7 GB → 167.9 GB at res 10. Dropping the id column itself is the same move, one step further. |

This is the single best calibration point you have. It is one order of magnitude below your L12 and two below L13, and it stops at res 10 (~65 m) precisely where a table stops being tractable.

---

## 2. The id axis: what published datasets actually do

### The four encodings, and who uses each

| Encoding | Used by | Evidence |
|---|---|---|
| **Implicit / positional (no id stored)** | HEALPix FITS `INDXSCHM='IMPLICIT'`; nextGEMS / DestinE / easy.gems Zarr; NICAM 220 m Zarr | [HEALPix FITS spec v0.6.0](https://healpix.sourceforge.io/data/examples/healpix_fits_specs.pdf); easy.gems repr literally reads **"Dimensions without coordinates: cell"** ([starter](https://easy.gems.dkrz.de/Processing/healpix/healpix_starter.html)) |
| **Id in the *path*, not the array** | HiPS (`NorderK/DirD/NpixN.ext`, `D = ⌊N/10000⌋·10000`); HATS (`Norder=k/Dir=d/Npix=p.parquet`) | [IVOA HiPS 1.0 REC](https://www.ivoa.net/documents/HiPS/20170519/REC-HIPS-1.0-20170519.pdf); [IVOA HATS Note 2025](https://www.ivoa.net/documents/Notes/HATS/20250822/NOTE-hats-ivoa-1.0-20250822.pdf) |
| **Dense explicit id column** | every H3 dataset ever published; HEALPix `INDXSCHM='EXPLICIT'` (`PIXEL` must be column 1); xdggs default | dem-terrain, Kontur, AusPIX |
| **Ranges / compaction** | IVOA MOC (NUNIQ + RANGE); H3 `compactCells`; S2 cell unions; OGC API-DGGS compact zone lists | [MOC 2.0 REC](https://www.ivoa.net/documents/MOC/20220727/REC-MOC-2.0-20220727.html) |

### The astronomy precedent is exactly your problem, solved in 2020

The HEALPix-in-FITS standard makes the choice a metadata keyword:

- `INDXSCHM = 'IMPLICIT'` — *"the pixel index is not given, but that p-th data value read correspond to the pixel p"*. Required for full-sky. `FIRSTPIX = 0`, `LASTPIX = 12·Nside²−1`.
- `INDXSCHM = 'EXPLICIT'` — with the rationale spelled out: *"it would be wasteful to store them in full sky maps with all unobserved pixels set to the sentinel value"*. Adds `OBS_NPIX` and a mandatory first column `TTYPE1 = 'PIXEL'`.

The spec's own worked example flips to EXPLICIT at ~10% sky coverage. The break-even is just "index width × coverage fraction vs 1".

### Real datasets, cell counts, storage

| Dataset | Grid / level | Cells | Storage | Id axis |
|---|---|---:|---|---|
| **NICAM 220 m** ([zarr](https://nowake.nicam.jp/files/220m/data_healpix_15.zarr)) | HEALPix zoom 15 | **12,884,901,888** | Zarr v2, 10 vars, `(1, 12884901888)` f32, **chunks `(1, 1048576)` = 4¹⁰**, blosc/lz4. 12,288 chunks/var; **51.5 GB/var**, ~515 GB total | **none** — no coords at all, empty root attrs (I read `.zmetadata`) |
| PanSTARRS DR1 HiPS | HEALPix order 20, 76% sky | **~1.0e13** | 47 M FITS tiles, **>20 TB/band**; tile = 512×512 = 4⁹ cells | in the path |
| FLAMINGO lightcones | HEALPix Nside 16384 | 3.22e9/map | HDF5, ~12.9 GB/map, 600 TB of maps | implicit |
| Planck | Nside 2048 | 5.03e7 | FITS, 0.6–1.9 GB | implicit |
| nextGEMS/DestinE ICON | HEALPix zoom ≤10 | ≤1.26e7 | Zarr, chunks `{"time":32,"cell":4⁷}` / `4⁸`, **~8 MB target chunk** | implicit (global) / `np.arange` (limited-area) |
| **dem-terrain** | H3 res 10 | 8.96e9 | Parquet 195.7 GB | explicit int64, sorted |
| Kontur Population | H3 res 8, populated only | subset of 6.9e8 | GeoPackage/SQLite, 6.6 GB (2.44 GB gz); H3 as **15-char hex string** + redundant WKB geometry per row | explicit string |
| AusPIX Major Geographies | rHEALPix L10, **Australia only** | **4.3e8 rows** | PostgreSQL on AWS + CSV | explicit |
| DES Y6 Gold | HEALPix Nside 16384 | 3.2e9 addressable | **healsparse** (see §4) | coverage-map offsets |
| Eco-ISEA3H | ISEA3H | ~1e6–1e7 | tabular | explicit |
| DGGS flow routing, Amazon/Yukon ([ESSD 2025](https://essd.copernicus.org/articles/17/2035/2025/)) | ISEA3H L10–13 | regional | GeoJSON/GeoPackage/GeoParquet, **9 MB** | explicit |

**Negative results, searched thoroughly:** Overture Maps does **not** key on H3 — it uses bbox structs + a **geohash-15 spatial sort** ([discussion #91](https://github.com/OvertureMaps/data/discussions/91)). Foursquare OS Places has **no H3 column** (lat/lon + WKB). Uber publishes the H3 *library*, no H3 *dataset* (Movement is decommissioned). WorldPop and GHSL are regular-grid GeoTIFF only. There is **no published DGGRID ISEA3H/ISEA7H global dataset at any scale**. OpenEAGGR is dormant and its ISEA3H mode has no working hierarchical index. PYXIS publishes nothing. **No DEM-on-DGGS other than dem-terrain.**

---

## 3. Chunking: how others handle hierarchy alignment

### The universal answer where the grid permits it: chunk = subtree

HEALPix NESTED makes this free — children of *p* are `4p…4p+3`, so **any chunk of 4ⁿ aligned to 4ⁿ is exactly one subtree**, and `12·4^k` is always divisible. Everyone independently landed there:

- NICAM: chunk = **4¹⁰** cells
- easy.gems: chunk = **4⁷ = 16,384** (3-D) / **4⁸ = 65,536** (2-D), *"We aim at a chunk size of about 8 MB"*
- HiPS: tile = **4⁹ = 262,144** cells
- healsparse: block = `2^bit_shift`, always a power of 4
- easy.gems limited-area: cuts the region **on whole chunks, not the AOI boundary** — *"maintains the proximity by keeping full chunks that follow the 4ⁿ relation"*

Your aperture-7 translation is `7ⁿ`: **7⁵ = 16,807** and **7⁶ = 117,649** bracket exactly the sizes the HEALPix community chose. At float32 those are 67 kB / 470 kB — you'll want a second axis or sharding to reach the 1–10 MB target.

### The formal ragged-shard model: OGC API-DGGS sub-zone order

[OGC API - DGGS Part 1: Core](https://docs.ogc.org/is/21-038r1/21-038r1.html) (approved standard) is the only spec that confronts pentagons head-on. A DGGRS defines a **deterministic sub-zone order**; a data packet is `dggrs` + one `zoneId` (the parent) + `depths` + **bare value arrays with zero per-value ids**:

```json
{ "dggrs": ".../ISEA3H", "zoneId": "C0-2B-A", "depths": [0,1],
  "values": {"t":[{"depth":1,"shape":{"count":7,"subZones":7},"data":[25.1,...]}]}}
```

And it explicitly permits ragged lengths: *"For most DGGRSs, that count of sub-zones is a pre-determined value consistent for all root zones **of a given type (e.g., hexagons or pentagons for hexagonal DGGRS)**"* ([CF discussion #447](https://github.com/orgs/cf-convention/discussions/447), Jérôme St-Louis). Its netCDF/Zarr profiles come in matched pairs: `zarr2-dggs` (*"one axis corresponds to local sub-zone indices"*) and `zarr2-dggs-zoneids` (global ids). **ISEA7H is a registered DGGRS**, and Annex B names IGEO7/Z7 as an alternate indexing of the same hierarchy.

### Zarr v3 now supports genuinely ragged chunk grids

The **[rectilinear chunk grid extension](https://github.com/zarr-developers/zarr-extensions/tree/main/chunk-grids/rectilinear)** (maintainer @d-v-b) lets each axis carry a list of per-chunk edge lengths, **run-length encoded as `[V, N]`**:

```json
"chunk_grid": {"name":"rectilinear","configuration":{"kind":"inline",
  "chunk_shapes": [[[16807, 686285], 14006]]}}
```

Merged into zarr-python in [#3802](https://github.com/zarr-developers/zarr-python/pull/3802), shipping **3.2 behind `zarr.config.set({"array.rectilinear_chunks": True})`**; also in `zarrs` and Icechunk. It supports **rectilinear shard boundaries with regular inner chunks**. [Earthmover's announcement](https://www.earthmover.io/blog/zarr-variable-length-chunks/) names your case: *"Hierarchical spatial schemes like HEALPix, S2, and H3 group cells by parent tile, with counts that vary across the surface… the natural read-and-write unit is the group, not a fixed-size slab."*

**[calc]** If you lay out each base-cell subtree as `[5 hex blocks] ++ [pentagon subtree]` recursively, the chunk sequence per base cell is exactly `[[7^k, M], P(k)]` — **one ragged chunk per base cell, 24 JSON entries for the whole L12 axis**:

| chunk | regular chunks / base cell | ragged tail | total chunks (L12) |
|---|---:|---:|---:|
| 7⁵ = 16,807 | 686,285 | 14,006 | 8,235,432 |
| 7⁶ = 117,649 | 98,040 | 98,041 | 1,176,492 |
| 7⁷ = 823,543 | 14,005 | 686,286 | 168,072 |

### The tabular-world equivalent

Warehouses cluster on the raw id (`CLUSTER BY h3`, `ZORDER BY h3`), but **cannot reason through it**: Google's [BigQuery spatial clustering guide](https://cloud.google.com/blog/products/data-analytics/best-practices-for-spatial-clustering-in-bigquery) states that `H3_ToParent(idx, 7) = @parent` **scans the whole table** because the bit ops defeat the query planner; you must use `BETWEEN H3_CellRangeStart(...) AND H3_CellRangeEnd(...)`, or **materialize a coarse parent column as the partition key**. That is a per-parent-cell shard key by another name. [raster2dggs](https://github.com/manaakiwhenua/raster2dggs) does exactly that on disk: `example.pq/h3_03=83bb09fffffffff/part.0.parquet`.

---

## 4. Sparsity (land-only)

**Zarr's native mechanism is simply not writing the chunk** — an absent chunk decodes as `fill_value`; `fill_value` is a *required* v3 metadata field. Zarr has **no native sparse-array support** ([zarr-specs #48](https://github.com/zarr-developers/zarr-specs/issues/48), [#245](https://github.com/zarr-developers/zarr-specs/issues/245) both open); everyone rolls a two-level scheme by hand.

**[calc]** For CopDEM on IGeo7 with global land ≈ 29%:

| chunk | cell area | ~chunk width | total chunks | retained |
|---|---:|---:|---:|---:|
| L12, 7⁵ | 62 km² | 7.9 km | 8.24 M | **~32%** |
| L12, 7⁶ | 434 km² | 21 km | 1.18 M | ~36% |
| L12, 7⁷ | 3034 km² | 55 km | 0.17 M | ~46% |

So plain missing-chunk elision gets you ~3× at 7⁵–7⁶ with **zero id storage and zero extra machinery**.

**The gold-standard sparse design is [healsparse](https://github.com/LSSTDESC/healsparse)** (Rykoff & Sanchez, LSST DESC — repo/Zenodo only, no paper). Two levels: a dense `int64` **coverage map** at a coarse nside (typically 32 = 12,288 entries) storing *pre-adjusted offsets, not booleans*, plus packed fine blocks of `nfine_per_cov = 2^bit_shift` cells. Lookup is one add:

```
index = pix_nest + cov_map[pix_nest >> bit_shift]
```

Empty init sets `cov_map[:] = -arange(npix_cov) * nfine_per_cov` so every absent block resolves into **block 0, permanently sentinel-filled** — branchless missing path, no search, no hash, no id column anywhere. Serializations: FITS (with **tile compression `blocksize = nfine_per_cov`**), **Parquet (`iopix=###/###.parquet` + `_coverage.parquet` mapping coverage pixel → row-group index, one row group per coverage pixel)**, and HDF5 (data reshaped 2-D as `(ncov, nfine_per_cov)`, one chunk per coverage pixel). All three make **the storage unit one coarse cell's subtree**. The Parquet layout is a drop-in template for a Zarr layout with a sidecar coverage array.

**How well does coverage compress?** [h3o-zip](https://github.com/HydroniumLabs/h3o-zip) gives real numbers: **mainland France at H3 res 11 = 267,532,208 cells → compact + CHT = 100.93 KiB.** Paris res 11 (54,812 cells): raw 438 kB → 933 B (0.14 bits/index). But: *"H3 compaction completely falls apart"* on sparse/linear sets (cycle lanes: 63.95 bits/index — no saving at all). Coastlines are more like a boundary than a blob, so expect the compacted form to be dominated by the coastal fringe, not the interior.

---

## 5. Conventions status — the decisive finding for you

### `zarr-conventions/dggs` — omitting the id array is a first-class mode

**[github.com/zarr-conventions/dggs](https://github.com/zarr-conventions/dggs)**, UUID `7b255807-…`, owner **@keewis**, maturity **"Pilot"**, listed by [GeoZarr](https://geozarr.org/conventions.html) under conventions under consideration. A `dggs` JSON object on a group or array:

| field | required | |
|---|---|---|
| `name`, `refinement_level`, `spatial_dimension` | ✓ | |
| `ellipsoid` | ✗ | projjson-shaped |
| **`coordinate`** | **✗** | *"If not provided, **the entire domain must be covered** and the `refinement_level` MUST NOT be `null`."* |
| **`compression`** | conditional | `"none"` \| `"compacted"` \| **`"ranges"` — shape `(n_ranges, 2)`** |

Those three `compression` values are *precisely* your options list for design problem 1, and its second worked example is HEALPix `refinement_level: 16` (5.15e10 cells) with **no coordinate at all**.

Caveats: the `v1` tag its own schema URLs point at does not exist yet; and xdggs's reader ([`conventions/zarr.py`](https://github.com/xarray-contrib/xdggs/blob/main/xdggs/conventions/zarr.py)) currently `NotImplementedError`s on *both* the missing-coordinate and the compressed-coordinate branches. **The convention permits what you need; no reader implements it.** Only HEALPix has a standardized parameter extension so far — an `igeo7` entry would be a natural PR.

### The scale problem was flagged, with numbers

**[xdggs#143 "Design/Datamodel Decision: Very large datasets"](https://github.com/xarray-contrib/xdggs/issues/143)** (jbusecke, May 2025, still open) is the thread:

> *"The current implementation basically faces a hard scaling stop … since we rely on a labelled xarray coordinate ('cell_ids'), which will by default be loaded to memory."*

- **keewis**: the `pandas.Index` is what forces materialization; fixing it needs *"a lazy index (e.g. using dask) or a new index implementation based on compacted cell ids / MOCs. Both options appear to be quite a bit of work."*
- **d70-t (Tobias Kölling, MPI-M / easy.gems)** — the "we don't want to store the ids" statement: *"One simple solution for global data we often use, is to not have a coordinate at all… where `cell_ids` would be the equivalent of `np.arange(n_cells)`, there's no real point in having them listed individually."*
- Same author, in [CF#433](https://github.com/cf-convention/cf-conventions/issues/433#issuecomment-2899274927): *"when experimenting with the 220m NICAM output… at `refinement_level = 15`, we would have about **100GB for the cell index coordinate alone**… That size rendered this approach pretty unusable: **xdggs currently basically relies on this approach**."*
- The related dask bug is literally `da.arange(12 * 4**15)`: [dask/dask#11997](https://github.com/dask/dask/issues/11997) — "96 GB, chunks of 256 MB", workers dying.

**Resolution:** [xdggs#151](https://github.com/xarray-contrib/xdggs/pull/151) merged a `HealpixMocIndex` (`index_kind="moc"`, backed by `RangeMOCIndex` from `healpix-geo`/`cds-healpix-rust`). It short-circuits when `array.size == 12·4**level` to `RangeMOCIndex.full_domain(level)` **without reading the ids**. But it is **HEALPix-nested only**, `index_kind="pandas"` is still the default, and you still must hand it an array object.

**For your grid specifically:** aperture-7 ISEA support in xdggs is a third-party plugin, [`xdggs-dggrid4py`](https://github.com/LandscapeGeoinformatics/xdggs-dggrid4py) (`@register_dggs("igeo7")`). Its `IGeo7Index` is a plain `PandasIndex` — **no lazy/MOC/range path at all** — and `from_variables` does `isinstance(var.values[0], str)`, forcing a compute. Its own cell-count table hard-codes L12 = 138,412,872,012.

### CF-1.13 went the other way, and only covers HEALPix

[PR #605](https://github.com/cf-convention/cf-conventions/pull/605) merged Dec 2025: `grid_mapping_name = healpix` with `indexing_scheme` and `refinement_level`, and standard name `healpix_index`. It **recommends materializing**: *"it is recommended to store strictly monotonic HEALPix indices in a coordinate variable, which allows an application to assess their monotonicity without expensive checking."* d70-t's implicit-coordinate proposal (`grid_reference = "dimension"`) did not survive. Generic DGGS in CF was explicitly deferred — davidhassell: *"solving the general use case would require many more degrees of freedom."* So **CF-1.13 does not bind an IGeo7 product**; the escape hatches it does sanction are compression-by-gathering (dense kept-index list) and compression-by-coordinate-subsampling / tie points (ranges).

---

## 6. Aperture-7 specifics: Z7, pentagons, and the padding trade

From [Kmoch, Sahr, Chan & Uuemaa, *IGEO7*, AGILE GIScience Series 6:32 (2025)](https://agile-giss.copernicus.org/articles/6/32/2025/) (I extracted the PDF text):

> *"A Z7 index is a 64-bit unsigned integer. The first four bits … indicate the base cell number (0 to 11), and the remaining 60 bits encode the digits for each resolution, using 3 bits per resolution. Each digit has a value from 0 to 6, **with a value of 7 used for digits greater than the resolution of the cell being indexed**."*

Base cells are the **12 pentagons** of the dual dodecahedron (unlike H3, whose 122 res-0 cells come from non-hierarchical aperture-3/4 refinements — this gives Z7 five extra indexed levels). Its Table 1 confirms `N(L) = 10·7^L + 2` with exactly 12 pentagons at every level: **L12 = 138,412,872,012**, **L13 = 968,890,104,072**. Its only published application is a 1 km-rasterized suitability model at **resolution 9**, stored **"in tabular form, indexing it by their Z7 cell identifiers … in the database."**

DGGRID also exposes **SEQNUM**, documented as *"linear address (1 to size-of-DGG)"* — a genuine dense implicit axis — with `seqnum_to_z7` / `z7_to_seqnum` converters ([duck_dggs](https://duckdb.org/community_extensions/extensions/duck_dggs.html), [webDggrid](https://github.com/am2222/webDggrid)). But SEQNUM derives from Q2DI (quad + i,j scanline), **not** the Z7 tree — the existence of the converters is itself the proof they are different orders. So SEQNUM buys density at the cost of tree alignment.

### The padding trade, quantified **[calc]**

H3's approach ([Cell mode docs](https://h3geo.org/docs/library/index/cell/)): *"in the case of the 12 pentagonal cells the indexing hierarchy produced by sub-digit 1 is removed at all resolutions"* — digit 1 after a run of 0s from a pentagon base cell returns `E_DELETED_DIGIT`. The digit space is **padded**, and the invalid slots are simply never used.

The cost of padding differs sharply between the two grids, because **all 12 of IGeo7's base cells are pentagons whereas only 12 of H3's 122 are**:

| grid | actual cells | padded (12 or 122 · 7^L) | waste |
|---|---:|---:|---:|
| H3 res 11 | 237,279,209,162 | 241,233,862,646 | **1.67 %** |
| **IGeo7 L12** | 138,412,872,012 | **166,095,446,412** | **20.0 %** |
| IGeo7 L13 | 968,890,104,072 | 1,162,668,124,884 | 20.0 % |

**But 20 % is the wrong number to fear, because the holes are contiguous and tree-aligned.** With the mixed-radix index `base·7^L + Σ dᵢ·7^(L−i)`, the invalid set per base cell is exactly *L* contiguous runs of sizes `7^(L−1), 7^(L−2), …, 7⁰`, each aligned to its own power of 7. So with a regular 7ᵏ chunk grid:

| L12 chunk | chunk slots | **never written (all fill)** | wasted cells *inside written chunks* |
|---|---:|---:|---:|
| 7⁵ | 9,882,516 | 1,647,084 (16.7 %) | 33,612 (0.000024 %) |
| 7⁶ | 1,411,788 | 235,296 (16.7 %) | 235,296 (0.00017 %) |
| 7⁷ | 201,684 | 33,612 (16.7 %) | 1,647,084 (0.0012 %) |

**A padded uniform 12·7^L axis with regular 7ᵏ chunks costs essentially nothing on disk** — 16.7 % of chunks are simply never written (Zarr's missing-chunk-is-fill_value), and real waste inside written chunks is under one part in 10⁵. You get: chunk index = the Z7 digit prefix, no rectilinear extension needed, no ragged tail, no id array. That is, as far as I can tell, the cleanest available answer to your design problem 2 — and it is the same trade H3 already made, just with a worse constant that Zarr absorbs for free.

**One documentation gap:** neither the IGEO7 paper nor igeo7.org states *which* digit value is deleted for Z7 pentagon descendants (H3 documents digit 1). The Table 1 counts imply exactly one child deleted per pentagon per level; you'd need DGGRID source or your own implementation to confirm the digit.

**A nuance worth knowing:** DGGAL distinguishes the *index* hierarchy (Z7: exact partition, 7 children per hex / 6 per pentagon) from the *geometric* sub-zone set (its `RI7H` `getSubZonesCount` is a nontrivial scanline computation, because aperture-7 hex children only approximately tile their parent). For array storage the index hierarchy is what matters — but don't take geometric sub-zone counts as chunk sizes.

---

## 7. Prior art worth reading before you design

1. **[moczarr](https://github.com/espg/moczarr)** (v0.5.0, PyPI, by @espg) — I verified this exists. *"Sparse-DGGS xarray reader for morton-hive zarr stores: MOC-declared domains, arithmetic shard paths, lazy dense views."* A digit tree of self-describing Zarr v3 leaves keyed by Morton ids, a static manifest plus hierarchical `coverage.moc` files, so *"a reader intersects an area of interest arithmetically instead of listing objects"*; **"NESTED is fabricated, never stored"**; registers with xdggs as a `"morton"` grid. Someone has already built roughly the thing you are about to design — on aperture 4.
2. **[healsparse file spec](https://github.com/LSSTDESC/healsparse/blob/main/docs/filespec.rst)** — the coverage-offset trick and the three block-aligned serializations.
3. **[IVOA MOC 2.0](https://www.ivoa.net/documents/MOC/20220727/REC-MOC-2.0-20220727.html)** — NUNIQ (`uniq = 4·4^order + npix`, self-delimiting via top set bit) and RANGE (sorted non-overlapping `[lo,hi)` at max depth, so set algebra is a linear merge). **Note the one place aperture 7 genuinely costs you**: 7ᵏ gives no "position of top set bit" trick, so a Z7 NUNIQ analogue needs an explicit order field or a digit-count sentinel — which Z7's trailing-7 padding already effectively provides.
4. **[IVOA HATS](https://www.ivoa.net/documents/Notes/HATS/20250822/NOTE-hats-ivoa-1.0-20250822.pdf)** — HiPS for tables, with **adaptive tile depth by data density** (split when a tile's bytes exceed a threshold → uniform *bytes*, non-uniform tree depth) and a **margin cache** holding points just outside each tile boundary so cross-matching never touches neighbours. Directly relevant to your halo/boundary-ring problem.
5. **[HiPS 1.0](https://www.ivoa.net/documents/HiPS/20170519/REC-HIPS-1.0-20170519.pdf)** — the 10,000-per-directory fanout rule that has run 47 M files in production since 2019.

---

## 8. Honest negative results

- **No dense hexagonal-DGGS array of any size has been published.** Every hex-keyed dataset found is tabular.
- **No hexagonal DGGS dataset above ~1e10 cells exists**, in any format. Above that, only HEALPix (NICAM 1.29e10 dense; PanSTARRS HiPS ~1e13 tiled).
- **No DGGS DEM product other than dem-terrain**, and it stops at res 10.
- **No published ISEA7H/ISEA3H global dataset** at any scale from DGGRID.
- **Nobody chunks per-base-cell into 12 separate arrays.** The only formal ragged-shard model is OGC's sub-zone order, and no bulk product implements it.
- **No published measurement of MOC/compaction size for a global land mask** at high order — my §4 reasoning is extrapolation from h3o-zip's France/Paris figures.
- **Zarr has no native sparse arrays**; missing-chunk-is-fill is the whole mechanism.
- **A correction to one line I'd otherwise have relayed:** it is sometimes said the largest dense HEALPix array is FLAMINGO at Nside 16384 (3.2e9). NICAM 220 m at zoom 15 is **4× larger and dense in Zarr** — I read its `.zmetadata` directly to confirm.

---

## 9. What I'd steal, in priority order

1. **Don't write an id array on the dense path.** `zarr-conventions/dggs` sanctions omitting it; HEALPix FITS has required it be omitted since 2020; the largest deployed HEALPix stores have no coordinate at all. At L12 an int64 axis is **1.107 TB** against a **554 GB** float32 payload — the index would cost 2× the data.
2. **Use the padded uniform `12·7^L` axis with regular `7ᵏ` chunks.** Holes are tree-aligned and contiguous; 16.7 % of chunks are never written; real waste < 1e-5. Chunk index = Z7 digit prefix. No rectilinear extension needed. Fall back to rectilinear + RLE (24 JSON entries) only if you must have an exactly-dense axis.
3. **Size chunks at 7⁶ = 117,649 and bundle a second axis (or shard) to ~8 MB** — the figure easy.gems, healsparse, HiPS and Dask all converged on independently.
4. **Get ocean elision for free** by not writing all-fill chunks (~3× at 7⁵–7⁶), before reaching for a coverage map. If you later need one, copy healsparse's offset-not-boolean coverage array with a permanently-sentinel block 0.
5. **Consider putting the level-*k* Z7 prefix in the chunk key / path** rather than a flat linear index — that is the HiPS/HATS pattern, it makes coarse↔fine joins arithmetic, and it survives at 47 M objects.
6. **Two upstream places your numbers would land well**: [xdggs#143](https://github.com/xarray-contrib/xdggs/issues/143) (an IGeo7 L12 = 1.1 TB id-axis datapoint would be the largest in that thread) and [zarr-conventions/dggs](https://github.com/zarr-conventions/dggs/issues) (proposing an `igeo7` parameter extension and confirming `"ranges"` semantics for Z7). You would be the first non-HEALPix implementation of the omit-coordinate mode.
---

## Appendix: sub-agent details not fully surfaced in the main report

Two sub-researchers fed the report above; incremental specifics worth keeping:

- **DGGRID practice ceiling**: the manual warns global binning "may fail at higher DGG resolutions due to memory restrictions" and that even reported successes at very high res "should be checked… to make sure they are not degenerate"; workshop examples top out at ISEA7H res 8 / **IGEO7 res 11**. dggridR/FRK ship precomputed grids only to res 6. So the generator ecosystem itself has never operated at our L12/L13.
- **No DEM on any DGGS exists** beyond walkthru.earth's dem-terrain (H3, res ≤ 10): the ISEA3H side has methods papers and code only.
- **dem-terrain measured density**: 10,450,894,746 cells, 195,675,804,510 B ≈ **18.7 B/cell** (int64 id + 5 float32 layers); v1→v2 saved 36 % purely by deleting id-derivable columns (geometry/lat/lon). Query pattern is `WHERE h3_index BETWEEN ? AND ?` with Parquet row-group pruning.
- **Kontur Population is GeoPackage (SQLite), not Parquet**, with a 15-char hex-string h3 column and ~165 B of redundant WKB polygon per row — an anti-pattern catalogue for our purposes. Res-8 row count unpublished.
- **Warehouse practice**: h3 always an explicit column (often STRING), clustered/Z-ordered; BigQuery has no native H3 and needs manual range predicates (25 MB vs 5.77 GB scan difference); no vendor documents a parent-cell partition key — raster2dggs' `h3_03=…` hive layout is the only on-disk parent-shard example found.
- **easy.gems/DKRZ dense HEALPix Zarr** is the cleanest production statement of the implicit-axis mode: `Dimensions without coordinates: cell`, NESTED ids exactly `0..12·4^z−1` — the contiguity H3's packed ids lack and Z7's padded digit space restores.
- **Uber Movement is decommissioned** (verified); Uber has never published an H3-keyed dataset. Overture geohash-sorts and rejected H3/S2 for its spatial key.

### From the H3-focused sub-report (new facts only)

- **H3 v4.1 added `cellToChildPos` / `childPosToCell` / `cellToChildrenSize`** — an official per-parent ordinal bijection (pentagon deletions handled). PR uber/h3#719's stated motivation includes verbatim: *"dense data association: mapping data arrays directly to H3 indexes by position, eliminating the need to store H3 indexes themselves."* This is the strongest single precedent for the positional/padded layout: the flagship hex-DGGS library built rank/inverse-rank for exactly this purpose. (DGG's `z7grid.jl` rank machinery is the same primitive.)
- **H3's deleted pentagon digit is digit 1** ("the indexing hierarchy produced by sub-digit 1 is removed at all resolutions", `E_DELETED_DIGIT`). The corresponding Z7 digit remains to be confirmed from DGGRID source — GBT comments in DGG's own `gbt.jl` suggest deleted digit `2` for northern / another for southern pentagons; verify before hardcoding the padded-axis hole pattern.
- **OGC API-DGGS explicitly blesses positional encoding**: zone-data Zarr responses may use "local indices into the deterministic zone order" instead of global ids; DGGAL's ISEA7H_Z7 adopts Z7 ordering for DGGRID/IGEO7 interop.
- **healsparse addresses Nside 131,072 ≈ 2×10^11 pixels** (DES 5000 deg² mask < 4 GB) — an address space the same order as IGeo7 L12, handled by coverage-offsets + dense blocks.
- H3's own id space is deliberately sparse/non-bijective (padding digits set to 7, 122/128 base slots, deleted pentagon subtrees) — "not a bijection onto 0..N−1" is normal for this family; nobody found it a reason to store ids.
