# Lazy Copernicus GLO-90 source and graph net-download count

Branch: `claude/copdem-lazy-source`, based at `claude/perf-ladder` `0d0f069`
and plain-merged with `claude/chunk-dag` (`0419988`). Worktree:
`/home/asinghvi17/geo/DGG-copdem-lazy`.

## 1. Design

`scripts/copdem_production.jl` now selects the elevation source with
`CONFIG.source = :real | :synthetic`. The committed production configuration is
`:real` and points at a distinct `copdem90-igeo7-l12-real.zarr` store path.

The existing source seam was retained:

- `TileIds` and `TiledDEM` still expose one disk-array chunk per listed 1-degree
  tile, in the Copernicus pixel position order.
- ArchGDAL remains the decoder used by the six-local-GeoTIFF path. No dependency
  was added: `ArchGDAL` and stdlib `Downloads` were already in `bench/Project.toml`.
- Local files selected by the existing `real = :auto | :none | [stems...]`
  setting remain overrides. All other listed tiles use `LazyCopernicusTiles` in
  real mode; all other listed tiles use `synthetic_tile` in synthetic mode.
- The Natural Earth mask remains exactly where it was: the synthetic tile
  generator applies it, while a real COG supplies its own values/nodata.

`LazyCopernicusTiles` is inert at construction. First access to a listed tile:

1. checks the persistent final cache name;
2. acquires that listed tile's own `ReentrantLock` (the production executor is
   multiple worker tasks in one process), then checks the cache again;
3. downloads with `Downloads.download` to `<stem>.tif.part`;
4. atomically renames the completed file to `<stem>.tif` on the same filesystem.

Only the final name is trusted, so a killed run can leave a `.part` but cannot
leave a truncated file that resume accepts. The default cache is
`bench/data/CopernicusDEM/tiles/` and can be changed with `COPDEM_TILE_CACHE` or
`CONFIG.tilecache`; it is gitignored. Transient network/HTTP 408/429/5xx errors
get four total attempts with exponential backoff. A 403 or 404 for a tile that
is in the official list is surfaced immediately with the stem, status, and URL.

Unlisted tiles return a correctly sized all-`NaN32` tile from the provider
without constructing a URL or making a request. In the production space they
are even cheaper: the `PartialGrid` contains only the official 26,475 listed
tiles, so open-ocean source chunks do not exist and downstream cells stay nodata.

Decoded COG dimensions are checked against the system lattice before a tile is
accepted: 1,200 rows and the latitude-band-dependent 1,200/800/600/400/240/120
columns.

The production store metadata now records `REAL` versus `SYNTHETIC` consistently.
Two bounded supporting scripts were added:

- `scripts/copdem_lazy_smoke.jl`: three-tile single-flight/decode/ocean checks
  and one materialized level-5-column regrid.
- `scripts/copdem_download_count.jl`: both production dependency graphs plus an
  18-HEAD transfer-size estimate; no tile bodies are downloaded.

## 2. AWS layout verification

The conventional URL was verified before being encoded, with an HTTP HEAD on
an exact entry from `tileList-glo90.txt`:

```text
https://copernicus-dem-90m.s3.amazonaws.com/
  Copernicus_DSM_COG_30_N00_00_E006_00_DEM/
  Copernicus_DSM_COG_30_N00_00_E006_00_DEM.tif
```

On 2026-08-21 S3 returned `HTTP/1.1 200 OK`, `Content-Type: image/tiff`,
`Accept-Ranges: bytes`, and `Content-Length: 490119`. The encoded layout is
therefore:

```text
https://copernicus-dem-90m.s3.amazonaws.com/<TILENAME>/<TILENAME>.tif
```

## 3. Bounded real-tile and regrid evidence

Exactly three tile bodies were downloaded through the new provider (3.653 MB
total on disk), followed by a cache-resume run that made zero GETs:

| tile | latitude | decoded dimensions | bytes | finite values | elevation range (m) |
|---|---:|---:|---:|---:|---:|
| `Copernicus_DSM_COG_30_N00_00_E006_00_DEM` | 0..1 N | 1200 x 1200 | 490,119 | 1,440,000 / 1,440,000 | -1.976 .. 1969.091 |
| `Copernicus_DSM_COG_30_N60_00_E010_00_DEM` | 60..61 N | 600 x 1200 | 2,700,509 | 720,000 / 720,000 | 62.000 .. 1006.599 |
| `Copernicus_DSM_COG_30_S90_00_E000_00_DEM` | 90..89 S | 120 x 1200 | 462,262 | 144,000 / 144,000 | 2711.177 .. 2836.327 |

Four tasks concurrently requesting the initially uncached `N00 E006` tile
returned the same cache path and incremented the successful-download counter
once. A direct access to the first unlisted lattice tile returned all NaN and
did not change that counter. No `.part` remained after any successful download.

The normal conservative lazy path was run with `N00 E006` as one source chunk
and one complete IGeo7 level-5 column (column position 54,998) as the rooted
level-12 destination:

| output cells | finite | NaN | range (m) | cached-run materialization |
|---:|---:|---:|---:|---:|
| 823,543 (`7^7`) | 823,543 | 0 | 0.000 .. 1590.069 | 6.39 s |

## 4. Graph executor: distinct source tiles

The spaces match `scripts/copdem_production.jl` and the chunk-DAG proof:
26,475 listed GLO-90 source chunks and 66,178 IGeo7 level-12 destination chunks,
each rooted at level 5. Conservative support radius is zero. Four Julia threads.

The explicit intersection with the official list is an identity for this
production `PartialGrid`, but it was retained and asserted in the counting
script because that is the executor's net-I/O boundary.

| graph | build time | edges | source chunks with >=1 consumer | intersect official list: net downloads |
|---|---:|---:|---:|---:|
| cap-conservative (`refine = nothing`) | 0.122 s | 326,392 | 26,475 | **26,475** |
| lon/lat-box narrow phase | 0.172 s | 250,821 | 26,475 | **26,475** |

The box phase removes 75,571 false-positive dependency edges, but no source tile
entirely loses its destination consumers. It therefore improves residency and
work scheduling, not the number of distinct COGs transferred for the closed
global production workload.

### Net tile distribution

Tile bands use absolute tile-center latitude.

| band | cap tiles | cap share | refined tiles | refined share |
|---|---:|---:|---:|---:|
| \|lat\| 0-30 | 7,318 | 27.6% | 7,318 | 27.6% |
| \|lat\| 30-60 | 7,027 | 26.5% | 7,027 | 26.5% |
| \|lat\| 60-90 | **12,130** | **45.8%** | **12,130** | **45.8%** |
| **total** | **26,475** | **100%** | **26,475** | **100%** |

Nearly half the object count is polar, although the longitude-reduced rasters
have fewer columns there.

### HTTP HEAD sample

Eighteen listed tiles were selected across latitude, longitude, and hemisphere:
six per requested broad band. Sizes are compressed COG object sizes, hence the
large terrain/content variation even at a fixed raster width.

| band | tile | center lat | pixel columns | Content-Length (bytes) | MB |
|---|---|---:|---:|---:|---:|
| \|lat\| 0-30 | `Copernicus_DSM_COG_30_S26_00_W131_00_DEM` | -25.5 | 1200 | 53,611 | 0.054 |
| \|lat\| 0-30 | `Copernicus_DSM_COG_30_S15_00_W077_00_DEM` | -14.5 | 1200 | 521,430 | 0.521 |
| \|lat\| 0-30 | `Copernicus_DSM_COG_30_S06_00_W036_00_DEM` | -5.5 | 1200 | 3,717,377 | 3.717 |
| \|lat\| 0-30 | `Copernicus_DSM_COG_30_N04_00_E029_00_DEM` | 4.5 | 1200 | 4,771,196 | 4.771 |
| \|lat\| 0-30 | `Copernicus_DSM_COG_30_N14_00_E093_00_DEM` | 14.5 | 1200 | 75,253 | 0.075 |
| \|lat\| 0-30 | `Copernicus_DSM_COG_30_N24_00_E153_00_DEM` | 24.5 | 1200 | 51,754 | 0.052 |
| \|lat\| 30-60 | `Copernicus_DSM_COG_30_S56_00_E158_00_DEM` | -55.5 | 800 | 50,267 | 0.050 |
| \|lat\| 30-60 | `Copernicus_DSM_COG_30_S45_00_W076_00_DEM` | -44.5 | 1200 | 169,855 | 0.170 |
| \|lat\| 30-60 | `Copernicus_DSM_COG_30_S35_00_W054_00_DEM` | -34.5 | 1200 | 752,884 | 0.753 |
| \|lat\| 30-60 | `Copernicus_DSM_COG_30_N35_00_E032_00_DEM` | 35.5 | 1200 | 736,316 | 0.736 |
| \|lat\| 30-60 | `Copernicus_DSM_COG_30_N44_00_E090_00_DEM` | 44.5 | 1200 | 4,356,594 | 4.357 |
| \|lat\| 30-60 | `Copernicus_DSM_COG_30_N55_00_E155_00_DEM` | 55.5 | 800 | 1,628,597 | 1.629 |
| \|lat\| 60-90 | `Copernicus_DSM_COG_30_S85_00_W151_00_DEM` | -84.5 | 240 | 988,613 | 0.989 |
| \|lat\| 60-90 | `Copernicus_DSM_COG_30_S76_00_W091_00_DEM` | -75.5 | 400 | 1,456,288 | 1.456 |
| \|lat\| 60-90 | `Copernicus_DSM_COG_30_S65_00_W057_00_DEM` | -64.5 | 600 | 131,203 | 0.131 |
| \|lat\| 60-90 | `Copernicus_DSM_COG_30_N64_00_E030_00_DEM` | 64.5 | 600 | 2,770,622 | 2.771 |
| \|lat\| 60-90 | `Copernicus_DSM_COG_30_N74_00_E090_00_DEM` | 74.5 | 400 | 1,872,042 | 1.872 |
| \|lat\| 60-90 | `Copernicus_DSM_COG_30_N81_00_E099_00_DEM` | 81.5 | 240 | 49,894 | 0.050 |

### Byte model and extrapolation

The model is the arithmetic mean `Content-Length` per broad latitude band,
multiplied by the exact net-tile count in that band. GB is decimal (`10^9` B).

| band | sample n | mean bytes/tile | cap tiles | cap GB | refined tiles | refined GB |
|---|---:|---:|---:|---:|---:|---:|
| \|lat\| 0-30 | 6 | 1,531,770 | 7,318 | 11.209 | 7,318 | 11.209 |
| \|lat\| 30-60 | 6 | 1,282,419 | 7,027 | 9.012 | 7,027 | 9.012 |
| \|lat\| 60-90 | 6 | 1,211,444 | 12,130 | 14.695 | 12,130 | 14.695 |
| **total** | **18** | — | **26,475** | **34.916 GB** | **26,475** | **34.916 GB** |

The point estimate is **34.916 GB for either graph** (about 32.52 GiB). With
only six very high-variance compressed objects per broad band it should be used
as a planning estimate, not a billing-grade total.

### Land-mask-restricted count

Skipped. The production source/destination spaces are built before the
15-arcsecond Natural Earth mask is applied per decoded source pixel, and neither
space carries a destination-column `nontrivial` flag. Deriving that subset is a
separate polygon/raster-to-level-5 spatial pass rather than a cheap property of
the constructed spaces.

## 5. Regression tests

All commands used at most eight Julia threads and no concurrent GC sweep thread.

| check | result |
|---|---|
| `test/systems/CopernicusDEM/runtests.jl` | **16,258 pass**, 3 expected broken, 57.0 s |
| `lib/GlobalRegridding/test/runtests.jl` from the bench workspace environment | **2,295 pass**, 1 expected broken, 81.7 s |
| first lazy smoke | all checks pass; 3 successful GETs |
| cache-resume lazy smoke | all checks pass; 0 successful GETs; regrid materialized in 6.39 s |

The full root suite was not requested and remains pending. A deliberately bad
listed-tile GET was not issued merely to exercise 403/404 handling, preserving
the network-request cap; the response classification is covered directly by the
provider logic.

## 6. Open questions

1. The 34.916 GB estimate is sensitive to COG compression/content. A later
   production planner could HEAD several hundred net tiles cheaply and fit the
   native Copernicus width bands (0-50, 50-60, 60-70, 70-80, 80-85, 85-90)
   instead of the requested three broad bands.
2. The single-flight lock is per tile and per Julia process, matching the
   production executor's worker-task model. If several independent Julia
   processes are ever pointed at one cache concurrently, the final-name rename
   remains atomic/trustworthy but cross-process duplicate GET suppression would
   need a stale-safe filesystem lock.
3. A land-mask-restricted destination count needs an explicit definition for
   “non-trivial” and a separate spatial construction; it is not present in the
   current production spaces.
