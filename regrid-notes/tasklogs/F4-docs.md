# F4 — documentation fixups

Text only. Nothing under `lib/` or `src/` touched, no Julia run.

## A. `README.md:172-175` — `from` takes the source directly

The prose taught `DGGSpace` as required; a bare grid, `CellVector` or
`CellLookup` is a source as it stands. Comment reworded to mirror how `to` is
described above it, and the example is now `from = grid`.

## B. `README.md` — `lib/GlobalRegridding` acknowledged

- `README.md:222` — Layout table row: the generic regridding engine (spaces,
  weights, plans, lazy executor), consumed by the main package for
  `regrid`/`plan_regrid`. Placed above the conformance-testing row.
- `README.md:311-312` — one sentence after the Tests paragraph: the subpackage
  carries its own suite, not run by the root `Pkg.test()`, runnable with
  `julia --project=lib/GlobalRegridding -e 'using Pkg; Pkg.test()'`.

**The `945,225 assertions, ~2m55s warm` figures at `README.md:307` are stale and
were deliberately left alone** — the orchestrator re-measures them after the
parallel code pass.

## C. `docs/src/tutorials/hydrology.jl` — scale comment, decorative timing

- `:47` — `leaf = 12` was annotated `≈ 1.1 km cells`, wrong by three levels.
  IGEO7 level 12 is `≈ 65 m` across, consistent with the page's own 2.31 M cells
  over the ~8,500 km² tile. Comment corrected; the level itself is unchanged.
- `:130-131` — `@time` removed from `GM.topographic_position_index` and
  `GM.flowaccumulation`. Neither makes a performance point, and a rendered page
  would show compilation-dominated numbers.
- No prose on the page calls the cells kilometre-scale, so nothing else changed.
  The `@time` on `DGG.regrid` at `:77` stays.

## D. `docs/src/tutorials/multiorder.jl:156`

`@time` removed from `A = DGG.regrid(box; to = lk)`; the page makes no timing
point.

## E. `examples/dimensionaldata.jl:190-199`

Section 5's point is that the regridder consumes the cube's own axis, so
`from = DGG.DGGSpace(DGG.PartialGrid(lk))` became `from = lk`, and the sentence
that named `DGGSpace` now says the lookup is a source as it stands.

No `DGGSpace` reference and no decorative `@time` remains in the four files.
