# W1 — one reference vector, and it lives in the block

- Date: 2026-08-25
- Card: `regrid-notes/2026-08-21-regridding-simplification-plan.md`,
  "Phase 5 — one final weight block representation", Task W1
- Commit: `Store reference weights in WeightBlock`

## 1. The representation

A final `WeightBlock` (`lib/GlobalRegridding/src/plans.jl`) holds three fields:
`weights`, the optional `denom`, and `reference`, the per-destination weight the
block's values are normalized against.

- On a **denominated** block `reference === denom`. It is the same object, not a
  copy of one, so there is no second opinion of a denominator to keep in step.
- On a block reporting **no denominator** — every point sample — `denom` is
  `nothing` and `reference` is the row sums of `weights`, computed once when the
  block is built, sparse and dense alike.
- `denom === nothing` therefore still records *which* of the two a reference is,
  which is all the deleted `hasdenom` predicate ever answered.

Both constructor forms are unchanged for a caller: `WeightBlock(weights, denom)`
and `WeightBlock(coo, ndst, nsrc)`. The reference is filled at construction, so
a third-party builder that emits through `WeightCOO` gets one without knowing it
exists.

`WeightCOO.denom` is `Union{Nothing,Vector{Float64}}` and is `nothing` until a
builder declares denominators, through `markdenominated!` or the first
`adddenom!`, either of which allocates the zero-filled vector. A point-method
COO never allocates one. `markdenominated!` before any `adddenom!` is what keeps
a builder that reports coverage for no destination — a conservative pair whose
cells are disjoint — producing a denominated block of zeros rather than a block
of row sums.

`CachedBlock` is `(block, bytes, used)`. It carries no reference of its own, and
neither `PerChunk` nor `CachedTile` nor an evicted-and-rebuilt entry copies
numerical state: a cache hit hands back the entry, and therefore the identical
`reference` object the block was built with. `blockreference!` is gone; nothing
needs a scratch reference anywhere, because every place that wanted one now has
the block's. `addreference!` remains, and is the accumulation of one block's
reference into a destination's running total.

The eager path applies a block through its stored reference: `applyplan!` takes
`block.reference` rather than filling a fresh vector, so a second application of
one plan allocates its two accumulators and nothing else. The lazy path's
`_applygroup!` reads the same field, so every slice of a chunk group shares it.

`Spilled` stores weights and the optional denominator, exactly as before — the
`.blk` and `.tile` bytes are unchanged and `SPILL_VERSION` stays `0x01`. The
reference is reconstructed on read, aliased to the denominator where the file
holds one and recomputed as row sums where it does not.

`_blockbytes(block)` counts the aliased pair once: a denominated block is its
weights plus one vector, and so is a block holding row sums instead. On an 18×18
diagonal fixture — 440 bytes of weights, a 144-byte vector — a denominated block
is 648 bytes where the separately cached reference made it 792, and a block with
no denominator is 648 as it always was. The saving is one `8 · ndst` vector per
resident denominated block, whatever the block's size. `CachedTile.bytes` and
`PerChunk.bytes` are sums of it and follow.

Nothing about weights, the meaning of denominators, missing policies or values
changes.

## 2. What it verifies

`lib/GlobalRegridding/test/test_executor.jl`, 24 assertions:

- a denominated block's `reference` **is** its `denom`; a block from a builder
  reporting none has `denom === nothing` and row sums as its reference, over the
  same weights;
- that builder's `WeightCOO` holds `nothing` before and after the build, so the
  empty-denominator path allocates no denominator at any point;
- dense and sparse weights of one operator give equal references, from separate
  vectors;
- a zero-coverage row has reference `0`, and the policy decides what that
  destination becomes — `Weighted` blanks it, `Extensive` returns the sum;
- empty sides build and apply: no destinations gives an empty reference, no
  sources a zero reference of full length, and a builder declaring denominators
  for an empty destination still references its own empty vector;
- `_blockbytes` is weights plus one vector for both kinds of block;
- a second eager application of one plan allocates fewer bytes than three
  reference vectors would occupy — two accumulators and no third array.

`lib/GlobalRegridding/test/test_lazy.jl`, 19 assertions:

- a `PerChunk` hit returns the same entry and therefore the identical
  `reference` object, with no rebuild; a `CachedTile` entry references its own
  block's;
- a denominated block and a block with none over identical weights differ on
  disk by exactly `8 · ndst` bytes, the denominator alone, so no reference is
  serialized;
- both round-trip with the reference reconstructed — aliased to the denominator
  where there is one, recomputed where there is not — with identical weights,
  denominators and byte accounting, and a spilled tile's blocks the same way.

The mutants these kill: a constructor that copies the denominator into a
separate reference; a COO that allocates a denominator no builder asked for; a
row-sum path that handles only one matrix layout; a reference left zero where a
method reports no denominator, which would blank every `Weighted` destination; a
cache that keeps a reference of its own and answers an equal copy; a reference
written to the spill file; byte accounting that counts an aliased vector twice.

## 3. What W2 needs

- The seam a unified builder must produce is `WeightBlock(weights, denom)`:
  supply the assembled matrix and either a denominator or `nothing`, and the
  reference follows. There is no third argument to compute and no reference to
  pass along a build path.
- A Conservative build that adopts its assembled CSC directly still needs a
  denominator vector of length `ndst`; handing that vector to the constructor is
  what makes it the block's reference, and copying it first would undo this
  card.
- The empty-side asymmetry the generic path carries is unchanged and still
  observable: an empty `dst_inds` returns before `markdenominated!`, so its
  block reports no denominator, while an empty `src_inds` returns after it and
  reports a denominator of zeros.
- `_blockbytes` takes the block alone. Any measurement of chunked residency
  compares against blocks that hold one vector, not two.
