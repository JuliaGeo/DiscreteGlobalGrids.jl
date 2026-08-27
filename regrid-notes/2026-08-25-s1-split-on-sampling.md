# S1 — the shared builders choose a build path by output sampling

- Date: 2026-08-25
- Card: `regrid-notes/2026-08-21-regridding-simplification-plan.md`,
  "Phase 9 — point-method admission", Task S1
- Commit: `Split weight construction on output sampling`

## 1. The seam

Every weight block, eager or chunked, is now assembled by one function in
`lib/GlobalRegridding/src/plans.jl`:

```julia
weightblock(method, dst_space, dst_inds, src_space, src_inds) -> WeightBlock
```

It dispatches once, on `outputsampling(method)`, into two branches:

- `weightblock(::DD.Lookups.Sampling, method, …)` — the area path. `Conservative`
  and every method reporting `Intervals` land here.
- `weightblock(::DD.Lookups.Points, method, …)` — the point path. `NearestCell`,
  `BilinearPoint` and any third-party point method land here.

Both call `pairblock`, which fills one `WeightCOO` through `buildweights!` and
returns a `WeightBlock` sized `(length(dst_inds), length(src_inds))` — the body
the two builders held before. No concrete method type takes part in the choice,
and no new trait or method-type union exists.

The two builders that reach it:

- `buildblock(plan::ChunkedPlan, dinds, sinds, dst_space)`
  (`lib/GlobalRegridding/src/plans.jl`) — one `(destination tile, source chunk)`
  pair. `blockfor`, its `PerChunk`/`Spilled` keys and its eviction order are
  untouched.
- `wholeblock(method, dst_space, src_space)`
  (`lib/GlobalRegridding/src/api.jl`) — the eager whole domain, as
  `1:ncells(dst_space)` against `1:ncells(src_space)`.

`wholeblock(::Conservative, …)`, which adopts the assembled sparse matrix
without a `WeightCOO` round trip, is a separate method and is unchanged; its
degenerate-side `invoke` reaches the generic `wholeblock` and therefore the area
branch.

The point branch's body is today the area branch's body: a point method that
supplies no weights of its own is served entirely by `buildweights!`, so
`NearestCell`, `BilinearPoint` and third-party point methods build
byte-identical weights.

## 2. What it verifies

`lib/GlobalRegridding/test/test_interpolation.jl`, testset
`output sampling selects the build path`, using two methods whose sampling is a
field so one method type can be put on either path:

- `T4SplitMethod` supplies weights only through `buildweights!`. On the point
  path it reaches that hook once per block and keeps every entry, gives the same
  block as the same method on the area path, and gives the same values eagerly
  and through `regrid!`.
- `T4TileMethod` supplies a constant block by defining
  `weightblock(::Points, ::T4TileMethod, …)` and no `buildweights!`. Reporting
  `Points()` it answers from both builders; reporting `Intervals(Center())` the
  identical type falls to `buildweights!` and throws. This is what pins the
  choice to the trait rather than to the method type.
- `NearestCell` and `BilinearPoint` give one answer eagerly and by chunk on a
  coarse-from-fine pair whose stencils cross source-chunk seams.

## 3. What S2 needs from it

S2 replaces the point branch's body. `TileWeights` enters at
`weightblock(::DD.Lookups.Points, …)` and, for a `ChunkedPlan`, at the point
branch of `blockfor`/`buildblock` that is lifted from it: a point method builds
one destination tile's stencils in a single pass, splits them by
`chunkat(src_space, i)`, and returns the block for the source chunk `src_inds`
names. The area branch and `pairblock` stay where they are, so a point method
that supplies no tile weights keeps the `buildweights!` fallback.

Deleting the `Points` branch outright is behaviour-neutral today, because it
forwards to the same assembly the area branch uses. It exists so S2 has one line
to rewrite and one place the point contract is documented.
