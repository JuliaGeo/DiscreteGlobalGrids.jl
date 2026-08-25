# Concurrent destination tile builds, measured

2026-08-25. `lib/GlobalRegridding` plus two methods in the DGGS package;
benchmarks `benchmark/copdem_nearest.jl` and `benchmark/point_tile_baseline.jl`.

Two independent levers on the point-method tile route: `TilePrefetch` builds the
weights of upcoming destination tiles in background tasks while the tile in hand
has its sources read, and `localindex(grid, p)` removes a duplicated membership
search from the per-cell constant those builds are made of. A `TileWeights` build
reads geometry and no source data, so tiles have no cross-tile dependency, and it
is the larger half of a cold read: 356.4 against 644.6 ns per destination cell.

## Where it lives, and why not elsewhere

- **On the array.** `LazyRegridArray` holds a `TilePrefetch` and
  `_readdestination!` takes its tile through it: one field, one call site. The
  queue starts after two consecutive tiles, so an array read once is unchanged.
- **Not in one read.** `eachchunk` makes a destination tile a chunk, so a sweep
  asks for one tile per `readblock!` and a read-scoped queue would see nothing.
- **Not in the storage.** `gettile!` is handed one tile number under a lock and
  holds neither the tiling nor the sampler a build needs.

## The bound

    depth = min(nthreads() - 1, weightlimit(storage) ÷ largest - 1)

- `weightlimit` is new: `maxbytes` for `PerChunk`, its memory tier for
  `Spilled`, zero elsewhere, which then reads tile by tile.
- Queued tiles plus the one in hand are at most `weightlimit` bytes, so the
  prefetch puts no more weight bytes in flight than the cache could hold. Depth
  is 0 before the first tile is built, at one thread, under `OUTER_PARALLEL`, and
  where the limit holds fewer than two tiles.
- The queue holds a `StableTask{CachedTile}` per tile: an evicted queued tile is
  served rather than rebuilt, and a fetched tile is typed without an assertion.
  The wave of source-block builds holds `StableTask{CachedBlock}`s the same way.
- A `ReentrantLock` guards each of a take's two regions, both of which leave
  `queue[k]` naming tile `served + k`; the wait for a tile is between them,
  holding nothing. Readers whose tiles interleave empty the queue on each other.

## Per cell

One 823,543-cell tile, each phase compiled on its own, medians of three
interleaved runs of each state.

| phase, ns/cell | before | after |
|---|---:|---:|
| `cellat(space, p)` on the DGGS space | 159.6 | **102.5** |
| its floor: one locate + one index lookup | 101.9 | 103.4 |
| the sample site `sites[i]` | 121.4 | 123.5 |
| the same sites through an in-order cursor | 126.9 | 127.1 |
| whole `tileweights` call | 356.4 | **310.7** |
| cold read of the tile | 644.6 | **593.4** |

- `cellat` on a subset searched its sorted ids twice: for membership while
  locating, then to read the located cell back as an index. `localindex(grid, p)`
  keeps the index the membership test produced. Semantics are unchanged.
- The sample site was measured and left alone: an in-order cursor over a
  `CellVector`'s windows is slower, the decode being 8.5 ns of a 121.4 ns site.
- Flat profile after, of ~177 working samples: `searchsortedfirst` 32 %,
  `getindex` 9 %, `gc_sweep_pool_page` 8 %, nothing else above 5 %.

## The gate

Two checkouts on one manifest at base `409c3e3`, arms interleaved rep by rep,
three reps each, `-t auto`, each arm waiting for an idle machine; medians. Every
arm verified 1,000 destination cells, none wrong, same mapped/outside split.

| arm | cold before | cold after | ratio | warm ratio | RSS ratio |
|---|---:|---:|---:|---:|---:|
| GLO-90 1 tile (2.31 M cells, 8 tiles) | 0.987 s | 0.461 s | **0.467** | 1.000 | 1.047 |
| GLO-90 strip (4.62 M, 13 tiles) | 2.407 s | 0.760 s | **0.316** | 0.970 | 1.047 |
| GLO-90 2×2 (9.16 M, 21 tiles) | 6.258 s | 1.407 s | **0.225** | 0.186 | 1.032 |
| GLO-30 1 tile (16.18 M, 33 tiles, L13) | 7.817 s | 3.292 s | **0.421** | 0.396 | 1.217 |

- CPU over wall on the cold sweep: 1.39→2.49×, 1.24→4.26×, 1.47→3.12× (1-tile,
  2×2, GLO-30). `point_tile_baseline.jl`'s nearest arms 0.971–1.011×: one tile.
- A sweep overlaps every tile but the first two: 6 of 8, 11 of 13, 19 of 21,
  31 of 33. The 2×2 and GLO-30 warm sweeps hold more tiles than the budget.
- **RSS cost.** GLO-30 rises 1.217×. Weights stay inside the budget (475.7 MB of
  536.9): it is the transient allocation of up to eight concurrent builds, and
  15.7 % of that sweep is GC against 8.1 % serially.
- **Eviction order.** On the 2×2 arm the cache left after the cold sweep is 14
  tiles / 501.7 MB against 15 / 530.2, a queued tile being inserted before it is
  served. Both build 21 tiles for 21, load the same 64 chunks, and answer the
  same 9,142,229 finite cells.
- Sanity arms on the rebased tree, GLO-90 one tile, two reps each: 0.469/0.462 s
  with the queue lock against 0.487/0.510 s without it, and 0.466/0.455 s once
  the task handles are typed.

## Suites

Rebased tree, all `--depwarn=yes`.

| suite | before | after | added |
|---|---|---|---:|
| `lib/GlobalRegridding`, `-t auto` | 4,757 / 0 / 0 / 1 | **5,002 / 0 / 0 / 1** | 245 |
| `lib/GlobalRegridding`, one thread | — | **4,998 / 0 / 0 / 1** | — |
| bridge cross-check | 383 / 0 / 0 / 12 | **390 / 0 / 0 / 12** | 7 |
| root | 1,047,337 / 0 / 0 / 17 | **1,047,344 / 0 / 0 / 17** | 7 |

Zero deprecation warnings in the library suite; the root suite has the same two
pre-existing `globalindices` ones. 54 serialized eager and chunked weight and
value entries are all `==` and `isequal` to the reference taken before.

What the added tests pin: concurrent and serial reads of a 16-tile destination
give `isequal` values and identical `LazyStats`, `builds`, resident tiles and
manifests; the depth arithmetic, directly and against a plan; sweeps at
a one-byte, a partial and an unbounded budget each build every tile once;
`OUTER_PARALLEL` is seen inside queued builds, not the first two; the queue's and
the wave's handle types, with `@inferred` fetches; a queued build raises its own
exception; `cellat` matches the two-step form at interior, boundary and outside
points on every system; and two tasks reading disjoint halves of one array, 24
rounds, agree with a serial read — without the lock that one fails in four runs.

## The floor

The build is 52.4 % of a cold read, against 55.3 % before. Prefetching moves it
off the critical path rather than removing it, with a ceiling of `(n - 2) / n`
for `n` tiles. What is left: 123.5 ns of sample site, 102.5 of location, 82.5 of
filing entries — none of them duplicated work any more.
