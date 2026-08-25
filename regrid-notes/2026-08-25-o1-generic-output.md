# O1 — a DGG destination's axis is generic output labelling

- Date: 2026-08-25
- Card: `regrid-notes/2026-08-21-regridding-simplification-plan.md`,
  "Phase 6 — generic DGG output and API forwarding", Task O1
- Commit: `Use generic output dimensions for DGG regridding`

## 1. The seam

A regrid labels its result in one place: `wrapoutput` on the eager route
(`lib/GlobalRegridding/src/executor.jl`) and `wraplazy` on the lazy one
(`lazy.jl`), both over `destinationdims(plan)`, which is the destination space's
`destinationdims(space, sampling)`. A `DGGSpace` answers that hook and supplies
nothing else:

```julia
GR.destinationdims(space::DGGSpace, ::DD.Lookups.Sampling) =
    (Cells(CellLookup(space.grid)),)
```

That one line is the whole of the DGG side of output. The `_DirectToDGG` and
`_ChunkedToDGG` plan aliases, the two `GlobalRegridding.regrid` methods
specialized on them, and `_ascube` — which took the parent out of the generic
result and relabelled it — are gone. This package defines no `regrid` method
at all now: a plan whose destination is a `DGGSpace` reaches
`GlobalRegridding`'s own application on both routes, and the eager route no
longer builds a `Cell`-labelled cube only to discard it.

The lookup is the one `_ascube` built, `CellLookup(dst_space.grid)`, so every
spelling of `to` — a level grid, a `CellVector`, a `CellLookup`, a
`MultiOrderCellSet`, a `PartialGrid`, a bare system — keeps the axis it had,
and `regrid!`'s accepted destination shapes are unchanged, the space's own
axis being the flat cell axis on a one-axis destination.

The hook ignores `sampling`. A cell holds one value however that value was
measured, so `BarycentricPoint`'s `Points` and `Conservative`'s `Intervals` name
the same cells the same way, and a `CellLookup` carries no sampling to vary.

## 2. What the generic path gained

The cells are the whole of a one-axis destination's shape, so neither route
needs to reshape anything to put a result on that axis. Both wrappers now say
so, generically — any space whose `destinationdims` names a single axis gets
them, and two or more axes still reshape, so a `RasterGrid` result is untouched:

- `wraplazy` labels the `LazyRegridArray` directly instead of building a
  `ShapedRegridArray` view whose shape is the array's own. Without it a DGG lazy
  result would have gained a reshaping view between it and its reads that it
  never had.
- `wrapoutput` labels the written array itself rather than `reshape`-ing it to
  the size it already has. `reshape` returns a fresh array header even when the
  dimensions match exactly, so this keeps `parent(result)` the very array
  `applyplan!` filled, which is what `_ascube` handed back.

`DimArray` construction still checks that the axis lengths match the array, so
a destination declaring an axis of the wrong length is refused where
`ShapedRegridArray` used to refuse it.

## 3. What it verifies

`test/systems/crosssystem/regrid.jl`, 9 assertions net:

- the axis a result carries *is* `destinationdims(plan)`, for the eager and the
  chunked plan alike, and every `regrid` method applicable to either plan type
  is defined in `GlobalRegridding`;
- a one-axis lazy result's parent is the `LazyRegridArray` itself, on the
  conservative and on the point route;
- the lazy result's dimensions — the cells and the pass-through month, with
  their lookups — equal the eager result's;
- `BarycentricPoint`, which asks the output for `Points` sampling where
  `Conservative` asks for `Intervals`, lands on the same axis, and its eager and
  lazy answers are equal cell for cell.

The mutants these kill: a DGG-side wrapper that labels a result but leaves
`destinationdims` answering `nothing`, which would leave `regrid!`'s size check
and every other reader of the hook disagreeing with the result; any `regrid`
method specialized on a plan with a DGG destination; a generic lazy path that
wraps a one-axis destination in a `ShapedRegridArray` anyway; a generic path
that drops or reorders the source's non-spatial dimensions or loses their
lookups; and a hook that lets `sampling` change which cells the axis names.

No test was deleted. `_ascube` had no test of its own — what it did is what
"the destination axis is the cells" asserts, and that testset is kept and
extended.

## 4. What is unchanged

A comparator regridded to DGG destinations on every route the bridge offers —
`RasterGrid` and DGG sources; `Conservative`, `BarycentricPoint` and
`BilinearPoint`; eager and lazy; a plan applied twice; `regrid!` into a
preallocated array; a complete level, a `CellVector` region and a `PartialGrid`
destination; a source with two non-spatial dimensions, a source with none, and a
plain matrix carrying no dimensions at all — and serialised each result's
dimension names and types, lookup types and values, parent type, element type,
size and values, before the change and after.

Twenty-five of the twenty-seven cases are identical field for field. Every
parent type is what it was, the lazy ones included, and every value is
bit-identical. The two cases that differ are the change itself:

- `destinationdims(plan)` answers the `Cells` dimension where it answered
  `nothing`, for the direct and the chunked plan;
- the `regrid` method applicable to a plan with a DGG destination is
  `GlobalRegridding`'s, where it was this package's.

## 5. What the phase's second card needs

- The bridge's remaining API surface in `src/regridding.jl` is target and source
  resolution alone: `regridgrid`, the `RegridTarget` union, the `_asspace`
  methods and `dimsource`. Nothing there labels output.
- This package defines no `regrid`, `regrid!` or `plan_regrid` method, so moving
  keyword defaults and forwarding `kwargs...` is entirely inside
  `GlobalRegridding` and has no DGG-side application to keep in step.
- Output labelling reads no keyword, `sampling` included, so a default that
  moves into `plan_regrid` cannot move a DGG result's axis.
