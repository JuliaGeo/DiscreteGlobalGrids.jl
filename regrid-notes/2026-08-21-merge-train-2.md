# Merge train 2: CR PR 138, DAG driver, and lazy Copernicus source

Date: 2026-08-21

Target: `claude/perf-ladder`, starting at `0d0f0690381b12b270bc333c8e45441734ceec72`

Final local tip before push: `22cf2cf8a3212044b7424296f02199371c84d88e`

## Merge order and commits

1. `claude/cr-pr138` at `5c3960b4b19a13c8dc0c1dcaf2421c94952647fa`
   fast-forwarded the target. This contributed:
   - `ea0df210ff8eb6b20f6735aad68f1a99940fe45b` — Pin ConservativeRegridding to PR 138 head
   - `5c3960b4b19a13c8dc0c1dcaf2421c94952647fa` — Advance the ConservativeRegridding pin to `claude/cached-dual-dfs`
   The resulting CR source SHA is
   `66ed54cbe8621018fc3c1df936c32f3f420cba57`. There is no merge commit for
   this step because the requested first merge was a fast-forward.
2. `claude/dag-driver` at `92737066cd1ff0a605e336c9e2b7b698ee0a2256`
   merged cleanly as `b5aa39bc6c57cda2819985ed8a4e7b9b3ae35e7b`.
3. `claude/copdem-lazy-source` at
   `9acbeda91d70a07ac22d638721ba8f5351a9ac96` merged as
   `6d342f3da93612e172d1a26b6f3eb4e690e5c084` after resolving the production
   script overlap described below.
4. Glue commit `22cf2cf8a3212044b7424296f02199371c84d88e` — Wire lazy tile downloads into DAG prefetch.

All three supplied branch tips are ancestors of the final tip. No branch was
rebased or rewritten.

## Conflicts and resolution

Git reported one content conflict: `scripts/copdem_production.jl` while merging
`claude/copdem-lazy-source`. `.gitignore` and `scripts/copdem_store.jl`
auto-merged without conflicts. The store edit retains source-dependent REAL vs
SYNTHETIC metadata, and `.gitignore` retains the persistent tile-cache rule.

The production-script conflict was resolved as follows:

- Kept the DAG driver's generic `SubtreeIds`, `TileBuilder`, `TiledDEM`, and
  policy-cache split. `TiledDEM` remains the disk-array wrapper over a cache;
  `TileBuilder` remains the sole loading seam.
- Added the lazy provider to `TileBuilder`. A build uses a named local GeoTIFF
  override first, the synthetic generator when no provider is active, or
  `loadtile(provider, ordinal)` for the lazy real source. The lazy path retains
  its validated raster dimensions and GDAL serialization.
- Kept the complete lazy-source API: `LazyCopernicusTiles`, `tilepath!`,
  `tileurl`, `tilecachepath`, `loadtile`, `readtile`, `listedtiles`, and
  `realtiles`, plus URL/retry/backoff/timeout/cache configuration.
- Kept `source = :synthetic` as the default and restored the synthetic
  production store path. `source = :real` constructs the lazy provider.
- Kept `refinegraph = false`. The unsafe lon/lat narrow phase was not enabled.
- Kept the DAG driver's affinity order, guided/tapered schedule, refcount and
  striped-LRU policies, prefetcher, and `reportcache` checks.
- Kept driver-reporting's session heartbeat and separate session/store closing
  banners. Lazy-provider decode/download counts were added alongside the cache
  report rather than replacing those banners.
- Adapted the newly added `scripts/copdem_lazy_smoke.jl` from the old monolithic
  `TiledDEM` constructor to `TileBuilder` + `StripedLRUCache` + the new
  `TiledDEM` wrapper. This was not a textual conflict, but it was required by
  the selected DAG-driver structure.

No feature from either branch was intentionally dropped.

## Glue design

The separate glue commit gives `Prefetcher` an optional `prepare(source)` seam.
Real mode translates the graph's source-chunk number through `tiles[source]`
and calls `tilepath!(provider, ordinal; demand = false)` before the existing
speculative `gettile!` call. Thus the prefetch pool ensures the GeoTIFF is on
disk before decode/cache publication. Synthetic mode supplies no preparation
callback and follows the pre-existing speculative cache-load path unchanged.

`LazyCopernicusTiles.ncold` counts a demand call that reaches an absent final
file, acquires that tile's single-flight lock, rechecks the file, and must start
the network fetch itself. A demand that waits behind an already-running
prefetch sees the final file after acquiring the lock and does not increment the
counter. The counter appears in real-mode heartbeats and the final source/cache
report. A focused policy test verifies that preparation happens before
speculative decode; the lazy smoke checks the cold-counter behavior for
concurrent demand.

## Validation

Commands used no more than eight Julia threads. Long commands ran under
`nice -n 10`; GC options, where supplied, used `--gcthreads=4` with no second
field.

1. Root resolve/instantiate:
   - The first `Pkg.instantiate()` detected the expected stale manifest and
     precompile could not load `Graphs` from `GlobalRegridding`.
   - The required fresh `Pkg.resolve(); Pkg.instantiate()` succeeded, adding
     `Graphs v1.14.0` and its resolver-selected transitive dependencies to the
     ignored manifest. `GlobalRegridding` and `DiscreteGlobalGrids`
     precompiled successfully. The CR test environments reported source
     `66ed54c`.
2. Synthetic byte identity, `scripts/copdem_dag_validate.jl`:
   - The first launch stopped before computation because this worktree lacked
     its untracked tile-list fixture. The tile list and Natural Earth shapefile
     components were copied read-only from the main checkout, and later removed
     from the target worktree after validation.
   - Rerun passed: 11 selected chunks, 8 stored files, 11.98 MiB,
     byte-identical. Both legacy and DAG sides had zero failures. The DAG side
     loaded 389 source tiles exactly once, with zero uncredited demands and
     complete retirement.
3. Lazy smoke, `scripts/copdem_lazy_smoke.jl`:
   - Copied the three named cached GeoTIFFs from the lazy-source worktree into
     this worktree's ignored `bench/data/CopernicusDEM/tiles/` cache.
   - Passed every smoke check with 0 fresh GETs. Raster dimensions were
     1200x1200, 600x1200, and 120x1200. The real-tile regrid wrote 823,543
     finite cells and 0 NaNs.
4. `GlobalRegridding` suite:
   - 2,295 pass, 0 fail, 1 broken (2,296 total), 1m35.5s.
5. Full root suite:
   - 1,043,846 pass, 0 fail, 17 broken (1,043,863 total), 17m21.0s.
   - The public-store network suite was skipped because
     `DGG_IO_NETWORK_TESTS` was unset, as in the normal local test run.

## Push

`git push origin claude/perf-ladder` succeeded, advancing the remote from
`0d0f069` to `22cf2cf` without force.

## Open questions

- `source = :synthetic` still honors any local files selected by
  `real = :auto` before falling back to the synthetic generator. Both input
  branches had local override precedence, so the merge preserved it as the
  smallest compatible resolution. If `source` is intended to be an absolute
  mode switch, a follow-up should make synthetic mode ignore local real-tile
  overrides (or default `real = :none`). The byte-identity validation used zero
  local overrides.
- The required lazy smoke used the existing three-file disk cache, so it did
  not exercise a fresh network transfer through the real-mode prefetch seam.
  Preparation ordering and provider single-flight/cold accounting are covered
  separately, but a controlled uncached integration run remains a useful
  follow-up. This is not a merge blocker.
