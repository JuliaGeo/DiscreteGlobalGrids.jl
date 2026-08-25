# P3 — fused tile-weight execution, measured

- Date: 2026-08-25
- Card: `regrid-notes/2026-08-23-barycentric-regridding-plan.md`, "P3 — fused
  tile-weight execution"
- Commits: `8c05a2f` (Build point weights per destination tile), `2800ba3`
  (Read exact source chunks for point tiles)
- Benchmark: `benchmark/point_tile_baseline.jl`

## 1. What the gate asks

> Work is O(destination cells + nonzeros), not O(destination cells × candidate
> chunks), and budget/spill tests pass.

Two things, then: the shape of the work, and that the tile survives the cache
policy it is now a unit of.

## 2. Which bullet each commit answers

| P3 card bullet | where it holds |
|---|---|
| the tile-level build/cache unit and exact dependency manifest | `8c05a2f`: `TileWeights(sourcechunks, blocks)` in `lib/GlobalRegridding/src/plans.jl`, built by `tileweights`, cached and locked by tile number through `gettile!` |
| route `BarycentricPoint` through it in lazy execution, retaining the old path for conservative and unopted custom methods | `8c05a2f`: `tilesampler(plan)` answers a sampler where `outputsampling` is `Points()` and `sampler(method, src_space)` returns one; `blockfor` builds the tile and answers a pair from it, `_readdestination!` asks for the tile before any source read. No method type takes part, so `NearestCell`, `BilinearPoint` and any method with no sampler keep the chunk-pair build |
| partition once by source local index, apply exact source blocks in ascending order | `8c05a2f`: every nonzero entry is filed under `chunkat(src_space, i)` at that chunk's chunk-local index in the one destination pass, blocks are finalized ascending, and the read applies them in that order |
| preserve eager/lazy and source-rechunking equivalence | `8c05a2f`, testset `point tile values match eager, lazy and another chunking`; `2800ba3`, testset `a point tile's manifest tracks the source's chunking` — the stencils do not move with the chunking and the manifest does |
| call-count assertion: one stencil query per requested destination | `8c05a2f`, testset `a point tile is one build with an exact manifest`: driving every pair the relation names leaves `method.placed[] == length(tile)`, and `< length(tile) * k`, which is what the same counting fixture records on the pair route in `df1d8b8` |
| read-count assertion: no source chunk outside the stencil union is read | `2800ba3`, testset `a point tile reads exactly the chunks its stencils name`: the recorded `readblock!` ranges are exactly the manifest's chunk ranges, a strict subset of the relation's row |

The gate's second half is `8c05a2f`'s testset
`a point tile builds once, evicts and spills whole`: two tasks asking for one
tile cause one build; a `PerChunk` budget too small for two tiles evicts one and
rebuilds it once, the same tile it was; and a `Spilled` tile comes back off disk
with its manifest and no destination pass. Its neighbour
`the conservative build unit is still one chunk pair` holds the other side.

## 3. The measurement

One run of `benchmark/point_tile_baseline.jl`, both routes timed in that one
process, so every ratio below is same-session by construction.

```
arm                                  time        allocated
one pass over the tile          67.801 ms      101,422,496 B
pair loop (4 candidates)       324.564 ms      221,472,320 B
plan construction                0.066 ms            4,320 B
cold tile read                 331.528 ms      254,579,792 B
warm tile read                   4.478 ms       19,275,648 B
tile weights, one pass         105.406 ms      212,468,720 B
plan construction, tile route    0.062 ms           78,304 B
cold read, tile route          105.509 ms      245,848,432 B
warm read, tile route            3.021 ms       19,349,536 B
```

The first five arms are `BilinearPoint` on the chunk-pair route; the last four
are `BarycentricPoint`, which supplies a `sampler` and is therefore built one
destination tile at a time. `one pass` is one `buildweights!` over the whole
source index range — the point-location floor, one unpartitioned COO.

| quantity | chunk-pair route | fused tile route | ratio |
|---|---|---|---|
| build of one tile's weights | 324.564 ms (4 pair builds) | 105.406 ms (1 `TileWeights`) | 3.08x |
| cold read of the tile | 331.528 ms | 105.509 ms | 3.14x |
| warm read of the tile | 4.478 ms | 3.021 ms | 1.48x |
| point locations | 1,500,000 | 375,000 | 4.00x |
| stencil entries kept | 1,500,000 | 1,500,000 | 1.00x |
| source chunks read | 4, each once | 4, each once | 1.00x |

Accounting for the same tile, from one further instrumented read of a clean
plan on each route:

```
                                  chunk-pair route     fused tile route
destination cells                          375,000              375,000
candidate source chunks (k)                      4                    4
manifest length                                  -                    4
point locations                          1,500,000              375,000
weight nonzeros held                     1,500,000            1,500,000
entries per destination cell                  4.00                 4.00
cache entries held                     4 pair blocks    1 tile, 4 blocks
weight bytes, plan accounting           40,147,488           40,147,584
weight bytes, summarysize               28,147,912           28,148,112
source chunks read                               4                    4
readblock! calls                                 4                    4
peak source bytes (residency)            4,147,200            4,147,200
```

Point locations on the fused route are counted, not asserted: the accounting
build runs through a sampler that records every `weightsat!` query and answers
with `BarycentricPoint`'s own stencil, and its manifest is checked against the
manifest of the tile the plan built.

### What the numbers say

- The work changed shape. 375,000 destinations are located once rather than
  once per candidate chunk, and the same 1,500,000 entries come out: the build
  is O(destination cells + nonzeros) where the pair loop was O(destination
  cells × candidates). The 4.00x fall in locations is exactly `k`. Time falls
  3.08x rather than 4.00x because only the location repeats on the pair
  route — each pair emits its own share, so the 1,500,000 entries are assembled
  once in total either way — and the fused pass additionally asks `chunkat` of
  every entry it files.
- Weight construction still dominates a cold read on both routes — 98.6% of
  the pair route's, 97.1% of the fused route's — so the 3.14x on the cold read
  is the build's 3.08x, not source I/O.
- Reads did not change, and were not expected to. The four candidates all carry
  weights on this fixture, so the exact manifest and the relation's row are the
  same four chunks, each read once. What removes a read is a candidate carrying
  no stencil entry; what P0 said was left to remove here was the repeated point
  location, and that is what fell.
- The tile costs 96 bytes more than the four pair blocks: its manifest.
- One cache entry replaces four. A budget bounds a tile as it bounds a pair, and
  the tile is what the budget evicts and the spill file holds.

### Metadata

- Julia 1.12.6 on macOS, Apple M4 Pro, `-t auto` = 8 threads and 8 GC threads of
  12 CPUs, `pmset -g` reports `powermode 2`.
- Source: a global 0.125-degree raster, 2880 x 1440 = 4,147,200 cells in 8 x 4 =
  32 chunks of 360 x 360, behind a counting `DiskArrays` source.
- Destination: one 750 x 500 = 375,000-cell tile spanning 60 degrees of
  longitude and 40 of latitude. `k = 4` candidate source chunks, read from the
  plan's relation rather than assumed.
- Statistic: every arm warmed once, then the minimum of three samples.
- Allocations: `@allocated` on one further call of the same arm.
- Weight memory: the plan's own `PerChunk` accounting (`storagebytes` over
  `_blockbytes`, manifest included), with `Base.summarysize` of the same blocks
  beside it. Storage is unbounded, so the bytes are what the weights occupy.
- Source memory: `residency(::LazyRegridArray)`; per-chunk read counts from the
  `readblock!` ranges the counting source records.

These absolutes belong to that machine state. Only the ratios above are
portable, and only because both routes were measured in one process.

## 4. P0's record, superseded

The file's header previously recorded the five chunk-pair arms alone, from an
earlier session on the same machine and power mode:

```
one pass over the tile          72.950 ms       96,081,312 B
pair loop (4 candidates)       320.903 ms      244,901,440 B
plan construction                0.058 ms            4,096 B
cold tile read                 308.418 ms      252,154,736 B
warm tile read                   2.834 ms       19,275,648 B
```

That session reported the pair loop at 4.40x one pass and weight construction at
99.1% of a cold read; this one reports 4.79x and 98.6% for the same unchanged
code. The spread between the two is why no ratio in section 3 is taken across
sessions.

## 5. What is not shown here

- Nothing measures a fixture where a candidate carries no weight, which is where
  an exact manifest removes a read rather than confirming one. `2800ba3`'s
  testsets assert that case; the benchmark does not price it.
- The 3.08x is one tile, one destination shape and `k = 4`. The ratio grows with
  `k` by construction and is not evidence about any other `k`.
- `BarycentricPoint` and `BilinearPoint` are different kernels, so the two
  routes' absolute times are not a controlled comparison of the routing alone;
  the counted point locations, the entry count and the manifest are.
- P4's conforming-DGGS source is untouched: every number here comes from a chart
  source, whose stencils are the chart-axis brackets.
