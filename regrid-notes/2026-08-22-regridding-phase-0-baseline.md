# Regridding simplification Phase 0 baseline

Date: 2026-08-22

This is the first baseline taken against the real dependency trees required by
the simplification plan. It supersedes results produced from the stale ignored
Manifest.

## Resolved sources

- GeometryOps rev `32c60581afa09f19aeaefee26446d95693ec52c4`, tree
  `632748249d0ba30f1faa997670f5bda2ec447551`.
- ConservativeRegridding rev `66ed54cbe8621018fc3c1df936c32f3f420cba57`,
  tree `9469fdc5e6eb9a1c304f61b906f2be778ae45437`.
- `GeometryOps.intersection_area`, `ConservativeRegridding.Trees.split_weight`,
  and `ConservativeRegridding.cached_dual_depth_first_search` are all present.

The ignored root Manifest was deleted and regenerated with `Pkg.instantiate()`.

## Correctness gates

- `GlobalRegridding`: 2,295 pass, 1 expected broken before adding the cache-bound
  test below.
- `test/systems/crosssystem/regrid.jl`: 148 pass.
- `test/systems/crosssystem/regridding_conservation.jl`: 76 pass, 12 expected
  broken.
- `test/systems/crosssystem/regrid_acceptance.jl`: 22 pass.
- `test/systems/CopernicusDEM/runtests.jl`: 16,258 pass, 3 expected broken.

## Planning baseline

Command:

```sh
julia -t 4 --project=benchmark benchmark/regridding_plan_baseline.jl
```

The live official GLO-90 list produces 26,475 source chunks. With the current
geometry tree it covers 66,175 IGeo7 level-5 destination chunks and produces
326,386 graph edges. Destination-space construction took 0.912 s. The graph
occupied 3,352,520 bytes and built in a 0.0589 s median over five samples
(0.0584 s minimum), allocating 10,654,960 bytes.

## Conservative production-tree baseline

`benchmark/conservative_block_baseline.jl` builds IGeo7 level 5 against a
360x180 global RasterGrid. Both direct and chunked modes produce 487,174
nonzeros, a 9,657,976-byte block, and sum to `4pi`.

- direct: 3.36 s, 799,699,328 allocated bytes, 754.5 MiB peak RSS;
- chunked: 2.31 s, 744,407,184 allocated bytes, 755.4 MiB peak RSS.

## Sparse round-trip baseline

`benchmark/conservative_roundtrip_baseline.jl` exactly reproduces the earlier
P1/P4 3600x1800 to 360x180 RasterGrid workload. Both paths produce 7,120,800
nonzeros, a 166,291,416-byte block, and sum to `4pi`.

- direct CSC adoption: 3.14 s, 12,333,000,896 allocated bytes, 2,640.2 MiB peak
  RSS;
- chunked CSC-to-COO-to-CSC: 3.16 s, 13,143,046,416 allocated bytes, 2,784.4 MiB
  peak RSS.

The live round-trip penalty is therefore 810,045,520 allocated bytes and
144.2 MiB peak RSS. The old 1,335 MiB figure included an earlier `findnz`-based
copy path and is not the current baseline.

## Cache-bound artifact

`test_conservative.jl` now pins `_TILE_CELL_CACHE_MAX` on both sides of the
boundary: exactly 65,536 destination cells are synthesized once and reused;
65,537 cells allocate no tile-wide geometry vector and remain on demand.
