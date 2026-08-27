# Storing Partial-Coverage DGGS Datasets

*A research synthesis on representing regional and global datasets on Discrete
Global Grid Systems (DGGS), with a focus on HEALPix and quadtree-on-polyhedron
grids.*

---

## The problem

We want to store datasets on a DGGS. The datasets may be fully global or may
cover only part of the globe — and the covered region (e.g. a country) often
does not fit neatly inside a single face's `(i, j)` matrix and may cross face
boundaries.

Two starting points, each with drawbacks:

- **Raw cell IDs (the xarray / `xdggs` approach).** Tree acceleration is weak:
  an index must be rebuilt on load, costs memory, and you lose the data
  contiguity of a face-matrix layout.
- **Dense face matrices (e.g. `(i, j, face)`).** Natural contiguity and
  ready-made 2D tiles, but a regional dataset wastes space on uncovered cells
  and the region may straddle faces.

The question: how do we actually store this — a ragged representation, or
something else? Has this been solved before?

---

## The key reframe

For quadtree-on-polyhedron DGGS — HEALPix, rHEALPix, S2, ISEA aperture-4 on
quads — **the `(face, i, j)` cube and a 1-D nested-order cell-ID array are the
same object.** The nested ID is just *face bits* concatenated with a
Morton/Hilbert interleave of `(i, j)`.

S2 makes this explicit: a cell ID is 3 face bits followed by the Hilbert-curve
position, the subdivision level is readable from the lowest set bit, and
containment reduces to a prefix comparison.

So the real design decisions are **not** "IDs vs. matrices." They are:

1. **Which linearization is canonical?** (Use *nested* / Morton / Hilbert order,
   never ring order or raw lat/lon slicing, which scramble spatial proximity.)
2. **What do you physically materialize for uncovered cells?** (Nothing, a
   sentinel fill value, or an explicit list.)
3. **Where does the index live?** (Arithmetic, implicit in sort order, or a
   small persisted manifest.)

The "rebuild the tree every time" complaint is entirely about (3), and it
dissolves once the index is either arithmetic or a tiny static table.

A consequence used throughout: an **aligned dyadic `2^k x 2^k` square in
`(i, j)` is exactly the descendant set of one parent cell.** Morton interleaving
means *subtrees and aligned squares are the same thing.* This is what lets both
storage models decompose a query into the same parent tiles.

> **Caveat — hexagonal systems.** The clean nesting / range / face-matrix
> properties assume **aperture-4 quad hierarchies**. For aperture-3 or -4
> *hexagonal* grids (H3, ISEA3H, ISEA4H) the parent–child relationship is
> ambiguous and hexes do not tile a face matrix exactly. There, the sorted-ID +
> prefix-range approach is the realistic path; the dense-cube option does not
> apply cleanly.

---

## Prior art

This exact problem has been solved at least three times, in three communities,
and the solutions converge.

### Astronomy — partial-sky HEALPix, standardized for decades

- **FITS partial-sky convention.** Full-sky files use *implicit* indexing (row
  number = pixel number). Partial-sky files are marked `OBJECT=PARTIAL` /
  `INDXSCHM=EXPLICIT`, and the first column holds the pixel indices. The
  gamma-ray-astronomy formats extend this with `LOCAL` and `SPARSE` variants for
  partial maps and cubes. This is precisely the "store covered cells + their
  IDs" model.
- **MOC (Multi-Order Coverage maps), IVOA standard.** A first-class object for
  *coverage itself*. Two packagings: **NUNIQ**, which packs an `(order, index)`
  pair into a single integer; and **RANGE**, which expresses everything as
  continuous intervals of indices at the deepest order. The same machinery was
  generalized to attach a value per pixel — this is how LIGO/Virgo distributes
  variable-resolution gravitational-wave sky localizations.
- **HiPS.** Exploits nested ordering to store tiles as fully-filled square
  arrays in a per-order directory tree — effectively a sparse tile pyramid.
- **Range sets in relational DBs (Healpix-Alchemy / SkyPortal).** Each tile is
  represented as the range of order-29 cells it contains; PostgreSQL range-set
  types make union/intersection and point-in-region trivial. For
  relational-database stacks this range representation is described as the only
  viable implementation for rapid spatial queries.

### ML / climate — dense faces

- **DLWP-HPX** arranges data as the twelve HEALPix faces as 2-D arrays and
  handles cross-face stencils with a *padding* operation that aligns and rotates
  neighbouring faces, with a special averaging rule for the missing corner
  pixel. The consistent east–west cell orientation of HEALPix lets them use
  *location-invariant convolution kernels*, unlike the cubed sphere.
- **km-scale climate (DKRZ "easy.gems", ICON, nextGEMS).** Hierarchical
  output = multiple resolution copies on the HEALPix grid, where zoom level `z`
  gives `12 * 4^z` cells. Establishes nested-order + multi-zoom as a working
  convention.

### Spatial databases — sorted IDs are the index

Because the Hilbert curve gives numerically close IDs to spatially close cells,
a plain **sorted index on cell IDs** supports spatial queries: a region query
becomes a *covering* set of cells, each of which is one contiguous ID range — a
range scan. (Google S2 on Spanner is the canonical write-up.)

### Earth observation — the active standardization effort

**GRID4EARTH (the EOPF-DGGS activity)** is building an end-to-end operational
pipeline for a DGGS-based EO/climate format: **HEALPix on Zarr v3**, validated
against Sentinel and DestinE products, feeding tools into services like
DestinE's DeltaTwin. HEALPix was already adopted in DestinE's Climate Digital
Twin, making it the natural foundation.

Two notable threads:

- **Ellipsoidal HEALPix.** HEALPix is defined on the sphere, but EO data is
  referenced to the WGS84 ellipsoid. The mismatch causes non-negligible area
  distortion at Copernicus resolutions and can bias zonal/regional analyses.
  GRID4EARTH extends HEALPix to WGS84 via the associated *authalic sphere*.
- **`healpix-geo`** (Python, built on the `cdshealpix` Rust crate) provides
  HEALPix indexing on both the sphere and the WGS84/GRS80 ellipsoid, with cover
  requests and MOC support, and no astronomy dependencies.

> **Maturity note.** This is prototype-grade plumbing. The Zarr v3 *rectilinear
> chunk grid* used below landed via a recent `zarr-python` PR, and the public
> demo still requires a fork branch of xarray. The **patterns** are stable; the
> **APIs** are moving.

---

## The design options, concretely

### A. Logically global, physically sparse dense cube *(good default for raster/EO/ML)*

Store `(face, i, j, ...)`, chunked, and **don't write empty chunks.** Zarr
supports this natively: a chunk equal to the fill value need not be written and
is assumed empty at read time (`write_empty_chunks=False`).

- A country on a global grid costs only its intersecting chunks.
- Waste is partially-filled boundary chunks, proportional to *perimeter / area*;
  fill-value runs compress to almost nothing.
- **No index at all:** cell ID → `(face, i, j)` → chunk key is pure bit
  arithmetic. Absent chunk = fill.
- Global and regional datasets share one schema and align trivially.
- Add a **coverage sidecar** (a MOC or compressed bitmap) so consumers can
  intersect footprints without touching data.
- Cross-face stencils need halo gathers — precompute these once per
  grid/resolution as static metadata (the DLWP-HPX pattern).

### B. Ragged per-face crops

Per touched face, store the minimal `(i, j)` rectangle plus its offset (e.g. a
DataTree of per-face arrays). Workable and xarray-friendly, but it is really
option A with irregular hand-rolled tiles, and it degrades for elongated or
fragmented regions (Chile; France + overseas territories). Pick it only if a
downstream consumer needs exact dense rectangles.

### C. Sorted nested-ID table (COO done right) *(good for tabular/join-heavy work)*

The `xdggs` layout — cell IDs as an index — but with nested order **mandated and
sorted**, stored in Zarr/Parquet/Arrow with per-chunk min/max ID statistics.

- The **sort order is the index.** Point lookup = binary search. Region query =
  covering cells → each a contiguous ID range → merged against chunk stats.
- Nothing is rebuilt; the "index" is a tiny static manifest.
- This is the astronomy `EXPLICIT` map and the S2-in-a-database pattern.
- Refinement: drop the ID column, store coverage as a compressed
  bitmap/interval list (≈ MOC RANGE) plus values packed in grid order, and
  address by rank — O(1) lookup, kilobytes of persisted index.

### Recommendation

Make **nested order canonical.** Keep the **dense Zarr cube (A)** as source of
truth for array/stencil/ML workloads with a coverage sidecar; **derive the
sorted table (C)** when you need joins. Converting between them is a
*permutation, not a resampling.* Adopt parent-tile-aligned chunking so the two
remain mutually derivable for the price of a reshuffle. Track the EOPF-DGGS work,
since they are settling these conventions in public now.

---

## The GRID4EARTH sample layout

GRID4EARTH published a sample store
(`data-taos.ifremer.fr/GRID4EARTH/no_chunk_healpix.zarr`) that is now the demo
dataset in the zarr-python docs. It is option C ("sorted COO done right"):

- A Zarr group with two **aligned 1-D arrays**: `cell_ids` (int64) and the data
  variable (float32), both shape `(222_442,)`.
- `cell_ids` carries a `level` attribute (10) and uses **nested ordering** via
  `healpix-geo`.
- Only covered cells are stored — 222k of the ~12.6M global zoom-10 cells
  (~1.8%) — sorted in nested order. No `(face, i, j)` axes, no fill values.

**The chunk grid is the spatial index.** Cells are grouped by their parent tile
at a coarser level (zoom 4 = level − 6), so each chunk holds all stored cells
under one parent — at most `4^6 = 4096`, fewer where coverage is partial. This
yields variable-sized chunks (e.g. `25, 645, 1510, ..., 4096`), serialized with
the Zarr v3 **rectilinear chunk-grid extension**, which stores run-length-encoded
chunk edge lengths inline in `zarr.json`:

```json
{"name": "rectilinear",
 "configuration": {"kind": "inline",
   "chunk_shapes": [[25, 645, 1510, 2363, 3203, 74, 769, 3963, 4096,
                     233, 1603, 2450, [4096, 2], 3327, 4047, [4096, 2],
                     "..."]]}}
```

Chunk boundaries coincide with parent-cell boundaries on the nested curve, so
the whole "tree" collapses to a tiny static table — read once, never rebuilt.

---

## Worked example: zonal mean over Germany

Suppose a **level-15** HEALPix store (~200 m cells), parent-aligned chunks at
zoom 9, and we want the mean NDVI over Germany.

A pleasant accident: Germany (47–55°N, 6–15°E) sits entirely above the 41.8°N
polar/equatorial face transition and does not cross the 0° meridian, so it lives
on a **single base face** — no face-boundary logic at all. (Mainland France
crosses 0° and splits across two polar faces: still just two key sets, no
special handling, because zonal reductions are permutation-invariant. The
face-alignment/rotation padding machinery is needed only for *stencils* across
face edges, not reductions.)

### On the sorted-ID layout (C)

1. **Polygon → ranges.** Cover Germany's boundary at level 15
   (`healpix-geo`/`cdshealpix`; on the ellipsoid to avoid the area bias).
   Internally a MOC; flattened to level 15 it is a sorted list of disjoint
   `[lo, hi)` nested-ID intervals — hundreds of ranges, kilobytes.
2. **Ranges → chunks.** Parent of cell `c` at Δ = 6 is `c >> 12`; a parent `p`
   owns `[p * 4^6, (p+1) * 4^6)`. Bit-shift each query range to a parent-tile set,
   intersect with the store's chunk table (binary search on per-chunk start IDs).
   Germany ≈ 9M level-15 cells ≈ ~2,200 full parent tiles + ~300 partial boundary
   tiles. **No I/O yet — pure arithmetic on metadata.**
3. **Fetch + slice.** Read only those chunks. Within each, `searchsorted` the
   range bounds against `cell_ids` to mask non-Germany cells in boundary tiles
   (~10–15% overscan, discarded after decode).
4. **Aggregate.** Equal-area cells ⇒ the zonal mean is an **unweighted** mean —
   no cos(lat) weights, no reprojection. Compute fractional-area weights only for
   the O(perimeter) border-straddling cells if needed.
5. **Many zones at once** (all NUTS-3 regions): rasterize zone polygons to a
   `zone_id` array on the same cell dimension once, then it is a plain groupby.
   Once data is on the DGGS, spatial joins become attribute joins on cell IDs,
   avoiding geometric intersection entirely. Hierarchical zoom copies let you
   preview at zoom 8 in milliseconds before touching zoom 15.

### On the dense cube layout (A)

Layout — choose `2^k x 2^k` chunks aligned to the grid so an aligned dyadic
square *is* one parent cell. A `(1, 1024, 1024)` chunk at zoom 15 is a zoom-5
HEALPix cell.

```
/ndvi      (face: 12, i: 32768, j: 32768)   float32, chunks (1, 1024, 1024)
           fill_value: NaN
/zone_id   same shape/chunks, uint16        # optional, written once
attrs: {grid: healpix, zoom: 15, order: nested -> (f, Morton(i,j))}
```

1. **Polygon → chunk keys.** Cover Germany at zoom 5 (= 15 − log2 1024) — the
   *same* cover call, coarser depth. The resulting cell IDs, bit-shifted, *are*
   the chunk keys: ~9 interior tiles + ~19 boundary tiles. No store listing, no
   metadata read — keys are computed; tiles outside coverage return fill.
2. **Read.** ~28 GETs of ~4 MB ≈ 110 MB raw vs. a 36 MB payload → overscan ≈ 3×.
   With 256² chunks (zoom-7 tiles): ~215 GETs ≈ 55 MB, overscan ≈ 1.5×. Zarr v3
   **sharding** keeps small logical chunks inside larger objects with byte-range
   fetches, mostly dissolving the tradeoff.
3. **Mask.** Interior tiles need nothing. For boundary tiles, rasterize the
   polygon onto cell centers, or precompute a boolean mask layer once (the dense
   analogue of the MOC; compresses to near-nothing). Semantic difference from C:
   there, "not Germany's data" is *structurally absent*; here, absence is a NaN
   sentinel, so cells inside Germany but outside acquisition coverage are
   excluded by value, not existence.
4. **Aggregate.** Equal-area ⇒ unweighted `nanmean`. Repeated zonal stats →
   write `zone_id` once, then `groupby(zone_id).mean()` over computed keys.

```python
tiles = cover(germany_poly, depth=zoom-10)          # zoom-5 parents
for t in tiles:                                     # f, i, j via bit-deinterleave
    block = ndvi[t.face, t.i0:t.i0+1024, t.j0:t.j0+1024]
    m = interior(t) or rasterize(germany_poly, t)
    acc += np.nansum(block[m]); n += np.isfinite(block[m]).sum()
mean = acc / n
```

---

## Where each model pays

| Dimension      | Dense cube (A)                                   | Sorted-ID table (C)                                  |
|----------------|--------------------------------------------------|------------------------------------------------------|
| Index          | None — pure arithmetic, absent chunk = fill      | Tiny chunk-start table in `zarr.json`                |
| Both           | Static, **nothing rebuilt** (original complaint dissolves either way)                                  |
| Bytes read     | 55–110 MB across 30–215 objects (bandwidth-bound, tunable) | ~40 MB across ~2,500 tiny objects (request-bound) |
| Storage        | Boundary-tile slack (fill compresses away), zero ID overhead; **wins for contiguous coverage** | Packs exactly; ~8 B/cell ID overhead (delta-compresses); **wins for scattered coverage** (ship tracks, point retrievals) |
| Compute shape  | 2-D tiles ready for convolution, pyramids, GPU   | Packed 1-D runs ideal for joins, groupbys, hierarchical scans |
| Alignment      | Two datasets at same zoom align by construction  | Need a (cheap, sorted) ID intersection first         |

---

## The deep takeaway

For the zonal-stats access pattern, **the two layouts execute the same
algorithm** — polygon → parent-tile cover → fetch → mask → reduce — through the
same cover machinery. Converting between them is a *permutation, not a
resampling.*

So the choice is workload-driven:

- **Raster-like EO/ML processing** → dense cube as canonical + coverage sidecar.
- **Tabular analytics / joins** → sorted nested-ID table.
- **Either way** → mandate nested order and parent-tile-aligned chunking, which
  makes the index static (killing the rebuild problem) and keeps the two
  representations mutually derivable.

---

## References / pointers

- **xdggs** — `github.com/xarray-contrib/xdggs` (design doc, cell IDs as an
  xarray index).
- **GRID4EARTH / EOPF-DGGS** — `eopf-dggs.github.io`; sample store at
  `data-taos.ifremer.fr/GRID4EARTH/no_chunk_healpix.zarr`.
- **healpix-geo** — `healpix-geo.readthedocs.io` (sphere + WGS84/GRS80 indexing,
  cover requests, MOC).
- **Zarr rectilinear chunk grids** — zarr-python docs "Rectilinear chunks"
  example; `zarr-developers/zarr-python` PR #3802; extension spec under
  `zarr-developers/zarr-extensions`.
- **MOC** — IVOA MOC 2.0 recommendation (`ivoa.net/documents/MOC/`); NUNIQ and
  RANGE packagings.
- **FITS partial-sky** — healpy `read_map`/`write_map`; gamma-ray-astronomy data
  formats (IMPLICIT / EXPLICIT / LOCAL / SPARSE).
- **DLWP-HPX** — Karlbauer et al. 2024, *JAMES*; `CognitiveModeling/dlwp-hpx`
  (face arrangement + padding).
- **DKRZ easy.gems** — `easy.gems.dkrz.de/Processing/healpix/` (hierarchical
  HEALPix output, limited-area HEALPix).
- **S2** — `s2geometry.io` (cell hierarchy, Hilbert encoding); S2-on-Spanner
  spatial indexing write-ups.
- **Healpix-Alchemy / SkyPortal** — range-set HEALPix in PostgreSQL.

*Note: this synthesis draws on web sources gathered during research and reflects
prototype-stage tooling as of mid-2026; verify current API status before
building.*
