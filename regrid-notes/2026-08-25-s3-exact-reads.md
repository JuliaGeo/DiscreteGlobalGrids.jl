# S3 — a point tile reads exactly the chunks its stencils name

- Date: 2026-08-25
- Card: `regrid-notes/2026-08-21-regridding-simplification-plan.md`,
  "Phase 9 — point-method admission", Task S3
- Commit: `Read exact source chunks for point tiles`

## 1. What decides a read

`_readdestination!` (`lib/GlobalRegridding/src/lazy.jl`) obtains a tile's
weights before it selects source chunks, and `_connectedsource!` takes the
selection from whichever object the tile's build unit is.

- **Chunk pair.** Every area method, and a `Points` method that supplies no
  `sampler`, has no `TileWeights`. Its selection is the plan's relation: a tile
  that is destination chunk `d` takes row `d`, and a derived tile spanning
  several destination chunks takes the ascending union of their rows.
- **Destination tile.** A tile with `TileWeights` takes `sourcechunks`, exactly
  — the ascending union of the chunks owning its nonzero stencil entries. The
  two are not intersected, and the row decides nothing; it is read only to
  check that the manifest sits inside it.

`knownempty` filtering runs after the selection on both routes, unchanged: it is
data-dependent, so it may drop a chunk the selection holds and may never add one
the selection does not.

Nothing buffers. No card dilates a destination cap, joins neighbouring chunks or
widens a manifest to approximate a stencil's reach.

## 2. What the relation is for

A chunked plan owns one relation whatever its build unit is, and a
`LazyRegridArray` still requires one. On the tile route it orders tiles, costs
waves from the caps it carries, and carries refcounts and prefetch. Its rows are
a **superset** of every tile's manifest and decide no read of a tile that has
one. That is the documented contract on `dependencies(plan)`,
`dependencies(A::LazyRegridArray)`, `_connectedsource!` and
`chunk_dependency_graph`.

## 3. The declared reach, and the check that enforces it

Cap overlap alone is **not** a superset of a point stencil's reach. A
destination cell finer than a source cell is bracketed by sample sites whose own
cells it never touches: a 0.1° destination cell centred at 0.35° brackets sites
at 0° and 1°, the 1° source cell's polygon begins at 0.5°, and the two extents
never overlap. The chunk owning that site is absent from a radius-free row, so
an intersection of manifest and row would silently drop its weight.

`supportradius(method, src_space)` is therefore a **bound on the stencil, not an
approximation of it**. A method that supplies a `sampler` declares how far its
stencils reach so that the relation stays a superset; `BarycentricPoint` answers
the larger chart-axis spacing on a source it queries through a cell chart.
`methods.jl` says so on the hook.

The second half of standing law 7 is checked rather than assumed. When a tile's
manifest names a source chunk no row of that tile holds, `_connectedsource!`
throws an `ArgumentError` naming the tile, the chunk, `supportradius`, the
method and the radius it answered. A silent wrong answer is the failure mode
this replaces, and the residency and refcount code is never asked to hold a
chunk the relation never named.

## 4. What it verifies

`lib/GlobalRegridding/test/test_lazy.jl`, three new testsets, 60 assertions.

- **Per-tile read count.** On the counting point fixture — a 20×10 destination
  in two chunks against a 36×18 source in 18 chunks, with stencils bracketing
  each destination across chunk seams — each tile's manifest is the chunks its
  brackets own, the recorded `readblock!` ranges are exactly that manifest's
  chunk ranges, and the manifest is a strict subset of the relation's row (2
  chunks against 4). Reading the whole destination in one call, which shares one
  source hold across the tiles, touches nothing outside the two manifests, and
  their chunks are fewer than the rows'.
- **A stencil reaching past cap overlap.** A 2×2 destination of 0.5° cells,
  wholly inside one 10° source cell and just past its sample site, against a
  source where every cell is its own chunk. Its stencils name four chunks; the
  destination's cap overlaps exactly one of them. With `BarycentricPoint`'s
  declared radius the relation holds all four, the read is all four and the
  values match eager. With a test method whose `supportradius` is the default
  `0.0`, the relation holds one, and the read throws an `ArgumentError` naming
  `supportradius` and the method rather than answering with three weights
  dropped.
- **A rechunked raster source.** A 72×36 raster under 4×2 and 6×4 chunkings,
  against a finer 17×11 destination in three tiles straddling both chunk seams.
  Each tile's manifest is the chunks owning its share of the whole-domain
  operator the chunk-pair builder produced, the reads are those chunks, both
  chunkings agree with the eager operator to one ulp, and the manifests differ.
  The ulp is the pre-existing reassociation of a stencil's partial sums across
  a different number of blocks, not a weight that moved.

The mutants these kill: reading a row chunk that carries no block; skipping a
manifest chunk; the intersection or the row deciding a read; a declared radius
of zero accepted on a chart source; a manifest that does not follow the source's
chunking.

Conservative read counts, residency statistics and values are unchanged, pinned
by the existing lazy testsets: a conservative tile's `_connectedsource!` still
answers its relation row exactly, and `L2 — locality` still reads each connected
chunk once.

## 5. What follows from it

- A space or method that answers dual cells of its own, rather than through a
  cell chart, must bound them in `supportradius` before a chunked read of it
  will run. Until it does, its stencils are found by cap overlap alone and the
  check refuses the first tile that reaches further.
- The exact manifest is what a point tile reads; the relation's rows remain the
  only thing a whole-run scheduler can order and prefetch against before any
  tile is built. Nothing here makes a scheduler exact, and nothing here lets one
  buffer adjacent chunks to pretend it is.
