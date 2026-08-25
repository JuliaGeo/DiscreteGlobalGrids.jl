# S2 — a point method's build unit is one destination tile

- Date: 2026-08-25
- Card: `regrid-notes/2026-08-21-regridding-simplification-plan.md`,
  "Phase 9 — point-method admission", Task S2
- Commit: `Build point weights per destination tile`

## 1. The unit

`TileWeights(sourcechunks, blocks)` in `lib/GlobalRegridding/src/plans.jl` holds
one destination tile's weights: the exact ascending source chunks the tile's
stencils name, and one `WeightBlock` per chunk in the same order. Rows are
chunk-local destination indices within the tile and columns chunk-local source
indices within that chunk, the convention `WeightCOO` documents. A candidate
chunk no stencil names has no block and is not in the manifest, so the manifest
bounds the tile's reads exactly. A point tile carries no denominator, as a point
sample has no coverage.

`tileweights(method, dst_space, dinds, src_space, sampler)` builds one in a
single pass over the tile's sample sites: one `weightsat!` per destination cell,
then every nonzero entry filed under `chunkat(src_space, i)` at that chunk's
chunk-local index — `i - first(ownedindices) + 1` where the chunk is contiguous
and a lookup table where it is not, both through `indexmap`/`localindex`.
Entries naming one source cell twice are summed into one. Blocks are finalized in
ascending source-chunk order.

Temporary storage is one `WeightRow`, one slot per source chunk number, and one
accumulator per contributing chunk holding the entries that reached it. Nothing
is sized by destination cells times candidate chunks, and no source chunking
changes a stencil: the whole stencil is found before it is partitioned.

## 2. Which methods take it, and where

`tilesampler(plan)` answers a `Sampler` when the plan's method reports `Points()`
and `sampler(method, src_space)` returns one, and `nothing` otherwise. No trait,
no method-type union and no concrete method name takes part. The default
`sampler(::AbstractRegriddingMethod, ::RegridSpace)` is `nothing`, so
`NearestCell`, `BilinearPoint` and a third-party point method keep the chunk-pair
build byte for byte.

- `blockfor(plan, (d, s), dinds)` finds or builds tile `d`'s `TileWeights` once
  and answers `s` from it, or with an empty block of the right shape when the
  manifest does not hold `s`.
- `buildblock` is always the single pair, through `weightblock`, which is also
  what the eager whole domain takes: a tile is worth building only where it is
  kept.
- `_readdestination!` asks for the tile once, before any source read, and applies
  its blocks in ascending source-chunk order.

## 3. Cache, budget and spill

The tile is the cache and locking unit. `PerChunk` keeps tiles in a dictionary of
its own beside the chunk-pair one, keyed by tile number; the two share the
recency clock, the entry count and the byte budget, and evict each other by
recency alone. `gettile!` builds a tile once: a request that meets a build in
flight waits on the storage's condition instead of starting a second, and the
build runs outside the storage lock.

A tile's bytes are its blocks' `_blockbytes`, reference vectors included, plus
its manifest, so `weightbudget` bounds a tile as it bounds a pair, and an evicted
tile is rebuilt in one pass. `Spilled` writes one `gr-<tag>-t<tile>.tile` file
holding the manifest and then the blocks, and reads both back without running a
destination pass.

`Conservative`'s `(destination tile, source chunk)` key, its `PerChunk` eviction
order and its `.blk` files are untouched. The two file kinds share the block
body's serialization; the bytes a `.blk` file holds are unchanged.

## 4. What it verifies

`lib/GlobalRegridding/test/test_interpolation.jl`, 44 assertions over a counting
point method with a sampler whose stencil is the two source sample sites
bracketing each destination in longitude:

- on the counting fixture the pair route walks — one tile, an 18-chunk source —
  the tile route performs `length(tile)` point locations rather than
  `length(tile) * k`, produces two entries per destination, and leaves one cache
  entry rather than `k`;
- the manifest is ascending, strictly increasing, exactly the chunks owning a
  nonzero entry, a strict subset of the tile's relation row, and every block in
  it is nonempty; every entry sits at its owner's chunk-local column, and some
  stencil straddles a source-chunk seam;
- eager, `regrid!`, a lazy read and a differently chunked source give one answer,
  which is the bracket's own interpolation; the stencils do not move with the
  chunking and the manifest does;
- a stencil naming one source cell twice leaves one entry of the summed weight;
- two tasks requesting one tile cause one build; a budget too small for two tiles
  evicts one and rebuilds it once, the same tile it was; a spilled tile comes
  back off disk with its manifest and no destination pass;
- an area method has no sampler, keeps the pair key, and holds no tile.

## 5. What S3 needs

- The manifest is on the cached tile, `plan.storage.tiles[d].sourcechunks`, and
  on the `TileWeights` a build returns. It is ascending, strictly increasing and
  a subset of `sourcesof(dependencies(plan), d)`.
- A tile-route read already loads only chunks the manifest holds:
  `_readdestination!` intersects the tile's manifest with the row
  `_connectedsource!` produced, after that function's own `knownempty` filtering,
  and applies a block for every chunk that survives. What is left is the
  selection itself — `_connectedsource!` still computes the relation's row first,
  nothing asserts a per-tile read count, and the relation is not yet documented
  as a superset that never decides a read.
- A derived tile spanning several destination chunks is keyed by its tile number,
  the number the pair cache already used, and its manifest is the tile's own
  rather than a union of chunk rows.
