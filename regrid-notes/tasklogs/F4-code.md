# F4 — code fixes from the review

Files touched: `lib/GlobalRegridding/src/{GlobalRegridding,shared,api,methods,
conservative,interpolation,rastergrid,executor,lazy,plans}.jl`,
`lib/GlobalRegridding/Project.toml`, `lib/GlobalRegridding/test/test_executor.jl`,
`src/regridding.jl`, `test/systems/crosssystem/regrid.jl`. `shared.jl` is new.

## Item 1 — source-aware target resolution replaces the type piracy

`plan_regrid` now resolves the source space first and passes it to a new
three-argument `_asspace(to, "to", src_space)`, whose fallback delegates to the
two-argument form, so existing extensions are unaffected. The main package's
pirate `GR.plan_regrid(::DD.AbstractDimArray; …)` and `_resolvetarget` are gone,
replaced by `GR._asspace(::AbstractHierarchicalGridSystem, name, src_space)`.
Both former bugs are fixed: a plain vector with a measurable `from` plans, and
`from = levelgrid(sys, l)` resolves through `_asspace` before reaching
`arealevel`. The two-argument bare-system method stays as the `from` error and
now says a bare system names no cells until a level is chosen — matched to the
source's areas as a destination, spelled `levelgrid(sys, l)` as a source.

## Item 2 — the dead `hasdenom` argument

`finalize!` lost its trailing `::Bool` from both implementations and its
docstring, and the six-argument `WeightBlock` convenience (no callers) is
deleted. `lazy.jl` no longer threads `denominated` through `_readdestination!`
and `_writechunk!`; `executor.jl` no longer computes `hasdenom` in
`applyplan!`. `hasdenom(::WeightBlock)` itself is kept as an accessor on an
exported type.

## Item 5 — merged helpers

New `src/shared.jl`, included first, holds one chunk-local index map pair
(`OffsetIndexMap`/`LookupIndexMap`, built by `indexmap`, queried by
`localindex`, `0` on a miss), `_WHOLE_SPHERE`, `_CELL_TREE_LEAF`, and `_padcap`.
`conservative.jl`'s `_FULL_SPHERE`, `_SUBTREE_LEAFSIZE`, `_cachedposition` and
`interpolation.jl`'s `_RangeLocalIndex`/`_DictLocalIndex`/`_localindex` are
gone. `BlockAreaOperator` now checks both local indices for `0` before clipping
and skips the pair — a block emits weights only for pairs both chunks contain,
which is what `build_weights!` already promises. The five copies of
`nextfloat(r * 1.0001 + 1e-12)` call `_padcap`.

`_rectcap` and `_chunkcap` are **not** identical — `_rectcap` goes through
`_sampledcap` (four tabulated corners for lon/lat boxes under 180°) where
`_chunkcap` always takes the sixteen-sample `_boxcap`, which is F3's deliberate
split. Both kept.

`_outputgrid`'s open-coded chunk-shape validation is replaced by one call to
api.jl's `_checkchunks`, so the messages exist once. The call is kept rather
than deleted outright because `ChunkedPlan` is exported and its constructors
bypass `plan_regrid`.

## Item 7 — `Spilled` docstring

States that the caller owns `dir` and its lifetime, and that the per-instance
tag makes the files unreadable garbage once the plan is dropped.

## Item 8 — interpolation guards

`_require_pointlocation` and `_require_centroids` (reflection probes that always
passed, and would have named the wrapper type) are deleted; `cellat` and
`cellcentroid` now raise their own errors. `_chart_required` lost its dead
`method` parameter and its "must implement …" lecture. What remains is the trait
check (`_require_chart`) plus the throwing chart fallbacks — the two idioms that
actually execute.

## Item 9 — lazy-only keywords behave alike

`budget` defaults to `nothing` in `plan_regrid`, `regrid` and `regrid!`; the
lazy path resolves it to `2^30`. An eager plan now rejects any explicitly set
`chunks`, `budget` or `storage`, naming the offending ones in one sentence.
Docstrings list `budget` as bytes for lazy reads and weights, default `2^30`,
and say the three keywords apply only to `lazy = true`.

## Item 12

`markdenominated!(coo)` sits next to `adddenom!` and replaces
`adddenom!(coo, 1, 0.0)`. `TileCells.cells::Any` becomes `initialized::Bool`
plus `cells::Union{Nothing,Vector}`, with `_cachedtree` still the function
barrier. `_indexrect` gained an `AbstractUnitRange` method that returns in O(1);
the scan survives only for scattered index sets. The `public` list drops
`isvalidvalue` (deleted), `weightbudget` and `databudget`, and gains `cellarea`.
The `_sourcespace` fallback and spill-version errors are one sentence each. The
module docstring says the source and destination spaces implement `RegridSpace`.

`LinearAlgebra` moved from `[deps]` to `[extras]` and the test target; `[compat]`
kept; `import LinearAlgebra` removed from the module.

## Tests

| suite | result |
|---|---|
| `julia --project=lib/GlobalRegridding -t 8 -e 'using Pkg; Pkg.test()'` | 423 pass, 1 broken |
| `test/systems/crosssystem/regrid.jl` | 68 pass |
| `test/systems/crosssystem/regrid_acceptance.jl` | 22 pass |

Two tests changed, both required by the items. `regrid.jl`'s
"every spelling of `to`" now asserts that `to = SYS, from = SRC` on a plain
vector plans at `arealevel(SYS, SRC) == LEVEL`, and that `from = GRID` plans at
the same level. `test_executor.jl`'s API-surface testset asserts that an eager
plan names the lazy-only keywords it was handed, and moved the
`budget = 0` positivity assertion onto `lazy = true`, the path that reads it.
