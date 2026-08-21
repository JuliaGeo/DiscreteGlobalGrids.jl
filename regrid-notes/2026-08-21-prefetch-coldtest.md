# CopDEM synthetic-source guard and cold-network prefetch acceptance

**Date:** 2026-08-21  
**Pipeline:** Copernicus DEM GLO-90 -> IGeo7 level 12, level-5 destination columns  
**Work branch:** `claude/perf-ladder`  
**Task-1 commit:** `083b2175173e6c274b866a85810d652e1089966d`  
**Cold-test harness commit:** `893fe09` (commands updated after the concurrent
`bench` -> `benchmark` project rename in `fe16005`)  
**Final verdict:** source-purity and network-path acceptance **PASS**. There is
no correctness blocker for the real GLO-90 run. The one warning is performance:
the Himalayan S3 objects took 40.7--65.3 seconds each in this environment. The
prefetch pool initiated every GET before demand and overlapped them strongly,
but a three-column all-cold start still exposed 50.3 seconds versus the warm
control.

Raw rows are in
`regrid-notes/2026-08-21-prefetch-coldtest.ndjson`. Scratch outputs were kept at:

- `/home/asinghvi17/geo/scratch-stores/prefetch-coldtest-cache/`
- `/home/asinghvi17/geo/scratch-stores/prefetch-coldtest-cold.zarr`
- `/home/asinghvi17/geo/scratch-stores/prefetch-coldtest-warm.zarr`
- `/home/asinghvi17/geo/scratch-stores/prefetch-coldtest-results.ndjson`

## 1. Task 1: `source=:synthetic` is absolute

The hybrid incident had two independent source selectors. `CONFIG.source`
controlled whether the lazy AWS provider existed, but `CONFIG.real` was still
passed to `realtiles`; its old `:auto` default allowed local GeoTIFF overrides
to replace analytic tiles even when `source=:synthetic`.

Commit `083b217` makes the contract explicit and early:

- the checked synthetic default is now `real=:none`;
- `effective_realspec(:synthetic, :none)` returns `:none`;
- every other synthetic real-tile setting, including `:auto` and an explicit
  stem list, throws `ArgumentError` before grids, stores, or caches are opened;
- real mode returns its existing `:auto`, `:none`, or explicit list unchanged.

The dependency-free seam is `scripts/copdem_source_mode.jl`, exercised by
`test/scripts/copdem_source_mode.jl` and included by the production driver.
Focused result: **7/7 tests passed**. A direct production-entry check produced:

```text
ArgumentError: source=:synthetic requires real=:none; got real=:auto.
Synthetic runs must never read cached real GeoTIFF overrides.
default: source=synthetic real=none
```

This verification was done while cached Copernicus files existed under the
worktree's `bench/data/CopernicusDEM/tiles/`; the absolute synthetic selection
does not construct a provider and admits no local override dictionary.

## 2. Pre-download plan

Planning ran before the cache path existed. It built the same unrefined
cap/broad-phase dependency graph used by the executor and issued HEAD requests
only for the predicted `Content-Length`; no GeoTIFF body was read and the fresh
cache path remained absent.

Three Himalayan real-land level-5 columns were selected:

| column | centroid (approximately) | graph tiles |
|---:|---|---:|
| 70493 | 88.72 E, 28.75 N | 4 |
| 73037 | 84.39 E, 28.20 N | 5 |
| 73054 | 86.45 E, 28.03 N | 5 |

There are 14 graph edges and 13 unique listed tiles because N28 E085 is shared
by columns 73037 and 73054. Predicted transfer was **64,371,900 bytes**
(64.372 decimal MB), safely below the 2 GB cap.

| predicted tile | HEAD bytes | observed latency (s) |
|---|---:|---:|
| `Copernicus_DSM_COG_30_N27_00_E083_00_DEM` | 5,290,704 | 43.138 |
| `Copernicus_DSM_COG_30_N27_00_E084_00_DEM` | 5,408,240 | 44.213 |
| `Copernicus_DSM_COG_30_N27_00_E085_00_DEM` | 5,471,895 | 42.081 |
| `Copernicus_DSM_COG_30_N27_00_E086_00_DEM` | 5,315,373 | 44.534 |
| `Copernicus_DSM_COG_30_N28_00_E083_00_DEM` | 5,218,401 | 49.406 |
| `Copernicus_DSM_COG_30_N28_00_E084_00_DEM` | 5,126,732 | 65.343 |
| `Copernicus_DSM_COG_30_N28_00_E085_00_DEM` | 4,731,595 | 41.595 |
| `Copernicus_DSM_COG_30_N28_00_E086_00_DEM` | 4,662,520 | 46.165 |
| `Copernicus_DSM_COG_30_N28_00_E087_00_DEM` | 4,636,116 | 55.861 |
| `Copernicus_DSM_COG_30_N28_00_E088_00_DEM` | 4,505,954 | 40.741 |
| `Copernicus_DSM_COG_30_N28_00_E089_00_DEM` | 4,594,604 | 44.232 |
| `Copernicus_DSM_COG_30_N29_00_E088_00_DEM` | 4,687,172 | 46.272 |
| `Copernicus_DSM_COG_30_N29_00_E089_00_DEM` | 4,722,594 | 46.432 |
| **total / distribution** | **64,371,900** | **sum 610.013; min 40.741; median 44.534; p90 54.570; max 65.343** |

## 3. Controlled cold run

Invocation used Julia 1.12.6, `nice -n 10`, eight Julia threads, the safe
single-field `--gcthreads=4` (`gcsweep=0`), three outer workers, refcount cache,
affinity order, tapered batches, prefetch depth 32, and fetch concurrency 16.
Source was real with `real=:none`, so all selected tiles went through the lazy
provider and none through local override selection.

Predicted versus actual:

- predicted: 13 unique stems / 64,371,900 bytes;
- actual cache finals: the same 13 stems / 64,371,900 bytes;
- missing predicted tiles: 0;
- unexpected tiles: 0;
- HEAD length versus final file size mismatches: 0.

Provider/cache counters:

| counter | result |
|---|---:|
| successful downloads | 13 |
| demand-cold downloads | **0** |
| prefetch requests issued | 13 |
| source decodes / cache loads | 13 |
| cache hits | 7 |
| joined cache loads | 0 |
| double-loaded source chunks | 0 |
| uncredited demands | 0 |
| live cache tiles at exit | 0 |

The successful-download count equals the number of unique final paths and each
source chunk's cache attempt count stayed at one, so the concurrent path made no
duplicate GET. There were no retry warnings. All 13 `.part` paths were observed
while downloads were in flight; **zero `.part` files remained at exit**.

### Latency hiding

All 13 `.part` files appeared at about 6.94 seconds after the driver timer
started, showing that the depth-32/concurrency-16 pool filled the whole bounded
frontier concurrently. The tile latencies summed to **610.013 seconds**, while
the complete cold driver wall was **72.709 seconds**: concurrency compressed
about 8.39 aggregate tile-seconds into each wall second.

It did not hide the all-cold critical path completely. The separate-process
warm control took **22.433 seconds**, so the cold run exposed **50.276 seconds**.
The per-column driver times provide a worker-stall proxy:

| column | cold (s) | warm (s) | cold excess (s) |
|---:|---:|---:|---:|
| 70493 | 45.5 | 11.1 | 34.4 |
| 73037 | 63.1 | 10.9 | 52.2 |
| 73054 | 53.6 | 9.4 | 44.2 |

There is no direct timer around a worker waiting on the provider's per-tile
lock. `cache_joined_loads=0` therefore does not mean zero wait: `prepare` owns
the provider lock before the cache flight is installed. The strongest facts are
that demand initiated zero downloads, every network request belonged to the
prefetch path, and the cold-minus-warm control quantifies the exposed start-up
cost.

The roughly 0.99 MB/s aggregate transfer rate over this small, high-relief
64.4 MB set would correspond to about 9.8 hours for the earlier 34.916 GB global
HEAD estimate if sustained. That estimate is illustrative, not an extrapolation
guarantee. Since the real regrid was previously sized as a much longer compute
run and its cold start occurs once, this is a performance warning rather than a
correctness blocker. A longer steady-state pilot would be appropriate if the
full run's go/no-go depends on depth 32 hiding all network time.

## 4. Output sanity

Every selected column wrote all `7^7 = 823,543` cells as finite real elevation:

| column | finite | NaN | minimum (m) | maximum (m) |
|---:|---:|---:|---:|---:|
| 70493 | 823,543 | 0 | 3,949.385 | 6,100.409 |
| 73037 | 823,543 | 0 | 269.160 | 5,644.099 |
| 73054 | 823,543 | 0 | 1,355.110 | 8,146.419 |
| **total** | **2,470,629** | **0** | **269.160** | **8,146.419** |

These ranges are plausible for the Nepal/Tibet/Himalaya footprint, including
the high-elevation column around 86.45 E.

## 5. Warm rerun and byte identity

The warm control ran in a fresh Julia process with the same columns, graph,
workers, prefetch depth, and concurrency. It reused the retained scratch cache
and wrote a second store rather than deleting the cold result.

- successful downloads: **0**;
- demand-cold downloads: **0**;
- decoded cache loads: 13, all from final cached TIFFs;
- stray `.part` files: 0;
- wall: 22.433 seconds;
- output finite counts and ranges: identical to cold;
- encoded `elevation` chunk files: **3/3 byte-identical**.

The paired chunk SHA-256 values were:

```text
70492.0  e0c8ace22418b85a2ef6972b004b6a5a8a63dedbcfdbdfddc4b3fd5ff73bcd27
73036.0  c84c19c385686cbbeccf94ce4753d6ffcc0490a3065922a24870379fece62850
73053.0  76bda3ee3eec08b355a0bf402d9eef8fcdeb233a2b1744de4d6b01add82293d5
```

## 6. Go/no-go

**GO on correctness.** The real network seam exercised the intended
ensure-downloaded path, prefetch owned all cold GETs, predicted and actual sets
matched exactly, single-flight and atomic rename held, the warm cache made zero
downloads, and output was physically sane and byte-reproducible.

**Performance caution, not a blocker:** the 40.7--65.3 second object latencies
are much higher than the 1--3 second latency used to size the original lookahead
model. On this deliberately tiny all-cold frontier the pool overlapped requests
but could not get ahead of workers before they started. For the full run, watch
the demand-cold counter and early throughput; pause if it rises above zero or if
the aggregate download rate fails to keep up with compute.
