# N1 — NearestCell on the point seam

- Date: 2026-08-25
- Card: this note. N1 is outside the barycentric and simplification plans; it
  applies their point seam to the method that was left on the chunk-pair route.
- Commit: `Sample nearest cells per destination tile`
- Benchmark: `benchmark/point_tile_baseline.jl`

## 1. What changed

`NearestCell` supplies a [`sampler`](../lib/GlobalRegridding/src/interpolation.jl),
so a chunked plan builds one `TileWeights` per destination tile and locates each
destination sample site once, where it previously located each site once per
candidate source chunk.

- `NearestSampler(space)` holds the source space and nothing else. It is not
  `BarycentricPoint`'s `Sampler`: that one carries prepared chart state, and a
  nearest query prepares nothing.
- `weightsat!(row, s::NearestSampler, p)` clears the row, asks
  `cellat(s.space, p)`, and either leaves one entry `(i, 1.0)` in the source's
  local index with `WeightsMapped`, or leaves the row empty with
  `WeightsOutside` — the status the seam already uses for a point the source
  covers nowhere. Zero bytes allocated once warm.
- `supportradius(::NearestCell, ::RegridSpace) = 0.0`, now stated rather than
  inherited from the default.
- `buildweights!(::NearestCell, …)` is unchanged. It is still the eager
  whole-domain build and the un-cached pair build (`weightblock` → `pairblock`),
  and a plan that has no sampler for its method still takes it.

`plans.jl` needed no change: `tilesampler`, `tileweights` and `_blockfor` name
no method type, and `tileweights` needs only `weightsat!`, `samplesites` and
`chunkat`.

`WeightRow`, `WeightStatus`, `ismapped` and `_addentry!` moved from
`barycentric.jl` to `methods.jl`, unchanged but for one docstring line. They are
the seam's vocabulary rather than barycentric's, `methods.jl` already holds the
other weight container (`WeightCOO`), and it is included before
`interpolation.jl`, which is what lets `NearestCell`'s sampler sit beside its
pair-route build instead of in a third file. None of the four names is public.

## 2. Eager and chunked are bit-identical, not merely close

A nearest stencil is one entry of exactly `1.0`. Nothing accumulates across
blocks, so no chunking can reassociate a sum, and `Weighted` divides by a valid
weight of exactly one. The two routes therefore agree bit for bit, and the tests
and the comparator assert `isequal`/`==` rather than `≈`. This is what separates
`NearestCell` from every other point method: `BarycentricPoint` and
`BilinearPoint` sum two to four entries, and a source rechunking splits that sum
across a different number of blocks, which moves the last ulp.

The one assumption underneath is the documented one:
`samplesites(space)[i] == cellcentroid(space, i)`. The eager build reads
`cellcentroid`, the tile build reads `samplesites`; a space where those differ
would break the identity here before it broke anything else.

## 3. Why the radius is zero

`supportradius` is a bound on the stencil, and chunk discovery is cap overlap
plus that bound. A destination sample site lies inside its own cell, so inside
the covering cap of the destination chunk that owns it. The source cell `cellat`
names contains that same point, so the source chunk owning it covers the point
too. Two covers sharing a point overlap, so a radius of zero already relates
every source chunk a destination tile can read. That is the containment every
method's discovery already rests on — the area route included — and not a second
assumption.

This is why the check in `_selectsource!` never fires for this method: the
manifest a tile emits names only chunks the relation already holds, including
for tiles straddling a source chunk seam. A chart stencil is the case that
needs a radius, because it reaches out to sample sites in cells the destination
never touches; a nearest stencil never leaves the point.

`samplesites` placing a site outside its own cell would break this. No shipped
space does: the DGG side is `cell_centroid`, pinned by
`cellat(grid, cell_centroid(grid, c)) == c`; `RasterGrid` and the toy space use
the cell centre; and `RasterGrid`'s longitude wrap maps a point to a 360-degree
representative of itself, so a cell it names still contains the point. A space
whose sample site can leave its cell would have to declare a radius bounding
that distance.

## 4. What it verifies

`lib/GlobalRegridding/test/test_interpolation.jl` and `test_lazy.jl`, 43 new
assertions; `test/systems/crosssystem/regrid.jl`, 7.

- **The route, counted.** `T6LocateCount` wraps a source space and counts
  `cellat`, delegating every other verb. On one 200-cell destination tile
  against a 72-chunk source whose relation names `k = 10` candidates, driving
  every pair the relation names leaves `located == length(tile)` on the tile
  route and `length(tile) * k` through `buildblock`, which is still the pair
  build. The tile is the one cache entry and no chunk pair is one. Kills a
  method that reports `Points()` but keeps locating per candidate chunk, and a
  `buildblock` that stopped being the pair route.
- **The manifest.** It is exactly `chunkat(src, cellat(src, site))` over the
  tile — four of the ten candidates — computed from `cellat` directly rather
  than from any build, with one entry of weight exactly `1.0` per destination.
  Kills a manifest that follows the relation rather than the stencils, and a
  weight that is not one.
- **The reads.** On the S3 fixture the recorded `readblock!` ranges are exactly
  the manifest's chunk ranges, on each of two tiles that straddle a source chunk
  seam, and the manifest is a strict subset of the relation's row. Kills a wrong
  bound: with a radius too small the manifest check would raise instead, and
  with the row deciding the read the count would be larger.
- **Outside coverage.** A destination reaching past the source on both axes
  takes no weight there and the missing policy decides it, identically on the
  eager, chunked, lazy and rechunked routes (`isequal`, so one missing equals
  another). Kills a sampler that emits a bogus entry for a point `cellat`
  refuses.
- **Parity.** The existing eager/chunked parity testset now also compares a lazy
  read and a differently chunked source, for both shipped point methods and two
  destinations, exactly for `NearestCell` and approximately for
  `BilinearPoint`. The DGG bridge testset does the same across the bridge:
  `NearestCell` off a rechunked raster into an S2 level-3 grid is bit-identical
  to the eager answer, where `BarycentricPoint` differs in 20 of 768 values by
  an ulp.
- **The radius.** `supportradius(NearestCell(), space) == 0.0` on a chart source
  whose `BilinearPoint` radius is large, so the method cannot be forwarding it.

### The values comparator

`n1_values.jl` (a toy pair, a coarse-from-fine toy pair, and a `RasterGrid`
pair) and `n1_dgg_values.jl` (a global raster into an S2 level-3 grid, and that
grid back into a regional raster), both local capture scripts rather than
tracked files, extend `s1_values.jl`. They capture, for `Conservative`,
`NearestCell` and `BilinearPoint`: the eager whole-domain block, eager values under both missing
policies, every `(destination chunk, source chunk)` pair `buildblock` produces,
every block `blockfor` answers, lazy values under both policies with the source
read ranges, the same under a second source chunking, and a lazy read under a
one-byte budget.

196 entries captured on the untouched worktree and again after the change (126
and 70): **188 identical, 8 changed, and all eight are read ranges, not
values.** Every weight, every block and every value is unchanged, on both
routes, in both directions across the bridge.

| entry | before | after |
|---|---|---|
| `raster/NearestCell/reads/*` | 6 chunk reads | 2 |
| `raster/NearestCell/rechunked-reads/*` | 6 | 2 |
| `dgg-to-raster/NearestCell/reads/*` | 12 | 8 |
| `dgg-to-raster/NearestCell/rechunked-reads/*` | 36 | 24 |

That is the exact manifest replacing the relation's row: a candidate chunk
holding no cell of the tile is no longer read. No other method's reads moved,
and neither did the two toy cases, whose destinations cover their whole source.

## 5. The measurement

One run of `benchmark/point_tile_baseline.jl`; both nearest routes are timed in
that one process, so the ratios are same-session by construction. `PairNearest`
forwards `NearestCell`'s `buildweights!` and supplies no sampler, so it is the
same weights on the chunk-pair route — the routing is the only difference
between the two arms.

```
arm                                  time        allocated
nearest one pass over the tile  44.263 ms       31,404,960 B
nearest pair loop (4)          226.055 ms       88,399,360 B
nearest tile weights            72.156 ms       88,399,600 B
nearest cold read, pair route  226.120 ms      107,678,432 B
nearest cold read, tile route   72.811 ms      107,680,240 B
nearest warm read, pair route    1.727 ms       19,275,648 B
nearest warm read, tile route    1.889 ms       19,275,552 B
```

| quantity | chunk-pair route | tile route | ratio |
|---|---|---|---|
| point locations | 1,500,000 | 375,000 | 4.00x |
| build of one tile's weights | 226.055 ms | 72.156 ms | 3.13x |
| cold read of the tile | 226.120 ms | 72.811 ms | 3.11x |
| warm read of the tile | 1.727 ms | 1.889 ms | 0.91x |
| stencil entries kept | 375,000 | 375,000 | 1.00x |
| weight bytes, plan accounting | 22,147,488 | 22,147,584 | +96 B |
| source chunks read | 4, each once | 4, each once | 1.00x |

The fall in locations is exactly `k`. Time falls 3.13x rather than 4.00x because
only the location repeats on the pair route — the 375,000 entries are assembled
once in total either way — and the fused pass additionally asks `chunkat` of
every entry it files. Weight construction is 99.2% of the pair route's cold read
and 97.4% of the tile route's, so the cold-read ratio is the build's. The warm
read is unchanged: one block application per destination either way. The 96
extra bytes are the manifest.

Point locations are counted, not asserted: the pair route through a method that
records what each build was asked to place, the tile route through a sampler
that records every `weightsat!` query and answers with `NearestCell`'s own
stencil. The counted build's manifest is checked against the plan's.

### Metadata

- Julia 1.12.6 on macOS, M-series, `-t auto` = 8 threads and 8 GC threads of 12
  CPUs, `pmset -g` reports `powermode 2`.
- Source: a global 0.125-degree raster, 2880 x 1440 = 4,147,200 cells in 8 x 4 =
  32 chunks of 360 x 360, behind a counting `DiskArrays` source.
- Destination: one 750 x 500 = 375,000-cell tile spanning 60 degrees of
  longitude and 40 of latitude. `k = 4` candidate source chunks, read from the
  plan's relation rather than assumed.
- Statistic: every arm warmed once, then the minimum of three samples;
  allocations from `@allocated` on one further call.
- Weight memory: the plan's own `PerChunk` accounting, manifest included. Source
  memory: `residency(::LazyRegridArray)`, with per-chunk read counts from the
  `readblock!` ranges the counting source records.

Absolutes belong to that machine state. Only the ratios are portable, and only
because both routes were measured in one process.

## 6. What this does not show

- **No read was removed on the benchmark fixture.** All four candidates hold
  cells of the tile there, so the manifest is the row. Where an exact manifest
  removes a read is the comparator's raster and DGG-source cases (6 → 2, 12 →
  8), which are correctness fixtures and are not timed.
- **`k = 4`, one tile, one destination shape.** The location ratio is `k` by
  construction; the time ratio is not evidence about any other `k`.
- **A raster source only.** Every timed number comes from `cellat` on a
  `RasterGrid`, which is two binary searches on prepared edge vectors.
- A production nearest run — Copernicus DEM tiles into IGeo7 — will exercise
  what the toy benchmark cannot: `cellat` on a DGG system, which is a hierarchy
  descent rather than a bracketed lookup and by far the dominant cost per point,
  so the ratio there is a ratio of DGG locations, not of raster ones; a
  destination tiling of thousands of tiles rather than one, where the tile cache
  and its budget do the evicting and a tile may be rebuilt; source chunks that
  are real files behind real I/O, where a chunk dropped from a manifest is a
  file not opened; a destination whose cells are much finer than a source chunk,
  so a tile's manifest is one or two chunks out of a much larger candidate row;
  and the pole, where the caps a relation is built from are widest and a
  candidate row is longest. None of those is measured here.
