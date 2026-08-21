# F5 — test honesty

Files touched: `lib/GlobalRegridding/test/{test_conservative,test_lazy,
test_rastergrid}.jl`, `test/systems/crosssystem/regrid_acceptance.jl`,
`.github/workflows/CI.yml`.

## Item 1 — the threaded bit-for-bit test compares two things

`test_conservative.jl`, "threaded and serial builds agree bit for bit". The
testset built only a serial reference and compared it to `conservative_block`,
which routes through `_intersectionareas` and short-circuits to serial at
`Threads.nthreads() == 1` — so single-threaded CI compared serial to serial.

The explicit `CR.intersection_areas` call is now a local `assemble(threaded)`
that takes the `GOCore.False()`/`GOCore.True()` flag, with a fresh
`BlockAreaOperator` per call because the inner operator carries a mutable
clipping cache. `assemble(GOCore.True())` is the threaded arm; both it and the
production `block` are compared against `assemble(GOCore.False())` in the
testset's existing style (`colptr`/`rowval` `==`, `nzval` and `denom` `.===`).
The production block stays tied to the pair at any thread count: at one thread
it is the serial side and the explicit `True()` arm supplies the threaded one;
above one thread it is a threaded side and the serial reference is the contrast.

`.github/workflows/CI.yml` gains `JULIA_NUM_THREADS: '2'` as a job-level `env`
on the `test` job, so `julia-actions/julia-runtest` and the conformance-package
step both run two-threaded and the threaded arm sees real concurrency.

## Item 2 — "disjoint chunks and mismatched manifolds"

Renamed to **"disjoint chunks keep zero denominators"**; no mismatch arm added.

Branch taken: rename only. Grepped the lib suite — the only occurrence of
"manifold" in a test title was this one, and nothing anywhere asserts the
`plan_regrid` ArgumentError from `src/api.jl:144`. But every space the lib
tests can construct is hard-wired to one manifold: `manifold(::ToyLonLatSpace)`,
`manifold(::RasterGrid)` and the file's own `manifold(::DensifiedCellSpace)` all
return `GOCore.Spherical(; radius = 1.0)` with no radius kwarg or field to vary.
Constructing a mismatched pair needs a new space type or a new field on
`ToyLonLatSpace`, i.e. new test machinery, so per the plan the testset is
renamed to what it tests and the mismatch arm is left unwritten.

## Item 3 — law numbering legend

`test_lazy.jl` gains a five-line legend under the file header, each line derived
from what the corresponding testset asserts (construction is free; locality;
residency; plan reuse; chunking invariance). Titles are now one number per law:
`L1 — construction is free`, `L2 — locality`, `L3 — the budget bounds
residency`, `L4 — plan reuse`, and `law 5 — chunking invariance` →
`L5 — chunking invariance`. The second `L4` became the unnumbered descriptive
`non-spatial slices reuse blocks`, since it is a corollary of L4.

## Item 4 — acceptance level pinned

`regrid_acceptance.jl:126` re-derived the level that built `DST`. Now
`@test DGG.level(DST.grid) == 7`. Confirmed by running `arealevel(IGeo7System(),
RasterGrid(DEM))` on the file's own 3600×1800 axes (prints `7`) and by the
acceptance run passing with the literal in place.

## Item 5 — structural spill count

`length(GR.spilledfiles(storage)) == heldstats.loads` compared against a load
counter from a different run. Now `== NX ÷ CHUNK`: the destination chunk has one
block per source tile it connects to, and those tiles are the southernmost tile
row, `NX ÷ CHUNK` of them — the same expression the read-count assertion above
already uses. The surrounding comment now says the count is the chunk pairs.

## Item 6 — `reference_boxpoint`'s dead `m`

Branch taken: **dropped the parameter**. The fast paths under test produce
corners only — `_rastercellcap` either calls `_cornercap` on the four tabulated
corners or `_boxcap(..., _CELL_CAP_SAMPLES)` with `_CELL_CAP_SAMPLES == 0`,
which is `m = 1`, four samples, all corners. Nothing under test densifies a cell
edge, so an `m = 8` arm would compare against a fast path that does not exist.
(The densified `_BOX_CAP_SAMPLES` path belongs to chunk caps, which the "wide
chunk extents" testset covers by containment, not bitwise.)

`reference_boxpoint` is replaced by `reference_corners(space, ix, iy)`, which
returns the four corners in ascending native order straight off the chart —
exactly `_cellcorners`' untabulated form. `reference_cellcap` loops over the
tuple, and `reference_cell` now calls the same helper instead of open-coding the
identical tuple. Sum order is unchanged, so the bitwise comparisons are
unchanged; dropping the `xlo + 0.0 * (xhi - xlo)` form also removes a latent
`-0.0`-vs-`0.0` discrepancy at a box edge that lands on negative zero.

## Tests

| command | result |
|---|---|
| `julia --project=lib/GlobalRegridding -t 2 -e 'using Pkg; Pkg.test()'` | 427 pass, 1 broken |
| `julia --project=lib/GlobalRegridding -t 8 -e 'using Pkg; Pkg.test()'` | 427 pass, 1 broken |
| `julia --project=test -t 8 … regrid_acceptance.jl` | 22 pass |

No new failures. The one broken arm is the pre-existing non-convex destination
clip in `test_conservative.jl` ("non-convex cells"), unchanged. The lib count
moves 423 → 427: item 1's threaded arm adds exactly four comparisons; item 6
neither adds nor removes assertions, and items 2 and 3 are titles and comments.
The acceptance count is unchanged at 22 — items 4 and 5 replace right-hand
sides, one for one.
