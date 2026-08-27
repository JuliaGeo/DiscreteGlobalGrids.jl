# O2 — one method per target spelling, one owner for the keywords

- Date: 2026-08-25
- Card: `regrid-notes/2026-08-21-regridding-simplification-plan.md`,
  "Phase 6 — generic DGG output and API forwarding", Task O2
- Commit: `Simplify regridding target and keyword resolution`

## 1. A spelling is a method

`GlobalRegridding` asks one question of a package that supplies spaces:
`_asspace(target, name[, src_space])`, "what space is this `to` or `from`?".
This package answers it once per spelling it accepts:

| spelling | space |
|---|---|
| `AbstractGrid` | `DGGSpace(grid)` |
| `AbstractCellLookup` | `DGGSpace(PartialGrid(lk))` |
| `AbstractCellVector` | `DGGSpace(PartialGrid(cv))` |
| `MultiOrderCellSet` | `DGGSpace(PartialGrid(CellVector(set)))` |
| `AbstractHierarchicalGridSystem` as `to` | `DGGSpace(levelgrid(sys, levelfor(sys, src_space)))` |
| `AbstractHierarchicalGridSystem` as `from` | an `ArgumentError` |

`regridgrid` — an open generic with four one-line conversions — and the
`RegridTarget` union that named its domain are gone. They were one indirection
with no second caller: `_asspace` was the only thing that asked, and the union
existed only to give that single `_asspace` method a signature. A grid already
stands for itself, so three of the four conversions were the constructor call
that now stands in the method body, and the fourth was the identity.

`PartialGrid` is an `AbstractGrid`, so a partial grid destination reaches the
first row rather than a conversion, exactly as before.

The bare-system rows are the two that carry behaviour. As a destination the
level is chosen from the source — `levelfor(sys, src_space)`, the level whose
cell size is closest in ratio to the source's median cell — which is the one
spelling that reads `src_space` and the reason the three-argument form of the
hook exists. As a source there is nothing to match a level against, and the
error says so, names the keyword it came in as, and names the spelling that
fixes it:

```
`from = S2System()` names no cells until a level is chosen. As a destination
the level is matched to the source's cell areas, but as a source you must name
it with `levelgrid(sys, l)`.
```

## 2. One owner for the keywords

`plan_regrid` states every keyword, every default and every check. `regrid` and
`regrid!` are two lines each: forward `kwargs...` to `plan_regrid`, apply the
plan it returns. Neither declares a keyword of its own, so neither can restate
a default or skip a check.

Their docstrings keep the full keyword list with the defaults spelled out —
what a caller reads is explicit even though what the code does is forward.

Three keywords do not forward. `dependencies`, `refine` and `narrow` describe a
plan somebody keeps: a relation to adopt, the narrow phase to build one with,
and the name that phase goes by. A one-shot regrid builds its plan and drops
it, so `regrid` and `regrid!` refuse them, saying which one was passed and
where it belongs.

The one default that is a value rather than a sentinel, the lazy budget, is now
written once as `GlobalRegridding.DEFAULT_BUDGET`. `plan_regrid` resolves the
API keyword against that name and the `ChunkedPlan` keyword constructor defaults
to it, where both spelled `2^30` before.

## 3. What it verifies

`test/systems/crosssystem/regrid.jl`, in "every spelling of `to` names the same
cells":

- every spelling resolves to a `DGGSpace` over exactly the cells it names — the
  grid's own cell list, not a matching cell count — through the `to` form and
  the `from` form alike, with a complete level and a subset for the collection
  spellings;
- the bare system's level follows the source for a raster source and for a DGG
  source, and its `from` error keeps its type and both halves of its message;
- `regridgrid` and `RegridTarget` are not defined.

`lib/GlobalRegridding/test/test_integration.jl`, in "one keyword surface, owned
by `plan_regrid`":

- `regrid` and `regrid!` declare no keyword but the splat, and `plan_regrid`
  declares all thirteen;
- the lazy budget default is `DEFAULT_BUDGET` through `plan_regrid` and through
  the `ChunkedPlan` constructor;
- the one-shot form produces what the two-step form produces;
- six bad keyword values raise the *same message* through `regrid`, `regrid!`
  and `plan_regrid`, and that message is an `ArgumentError` rather than
  whatever a call that quietly succeeded would return.

`lib/GlobalRegridding/test/test_chunkgraph.jl` keeps its "a narrow phase cannot
be supplied after the plan exists" testset and now checks all three relation
keywords through both convenience methods.

`lib/GlobalRegridding/test/test_integration.jl`, in "a one-axis destination is
labelled, not reshaped", pins O1's short circuits on this package's own toy
space rather than only through the bridge: `ToyCellAxisSpace` wraps a toy space
and answers `destinationdims` with a single axis, and for a dimensional source
with one non-spatial axis the eager wrapper hands back the very array it was
given (`===`) and the lazy wrapper the `LazyRegridArray` itself.

The mutants these kill: a spelling routed to the wrong grid — a lookup or a
multi-order set resolving to its complete level rather than its own cells; a
bare-system error downgraded to a `MethodError` or stripped of the spelling
that fixes it; a default restated on `regrid` and drifting from
`plan_regrid`'s; a check `plan_regrid` performs that the one-shot form walks
around; a relation keyword silently accepted by a form that drops the plan; and
a one-axis destination that gets a `reshape` or a `ShapedRegridArray` between
the result and the array that was written.

## 4. What a caller sees differently

Nothing about a result, a default, a keyword name or a target spelling. Two
error shapes change, both on misuse:

- `regrid(...; refine = f)` and its `dependencies`/`narrow` siblings raise an
  `ArgumentError` naming the keyword and pointing at `plan_regrid`, where they
  raised a `MethodError` listing every keyword the method did have.
- A misspelt keyword still raises a `MethodError`, but it now names
  `plan_regrid` — the function that owns the keyword — rather than `regrid`.

## 5. What Phase 7 needs

- The bridge's regridding surface is now spaces and hooks alone. There is no
  `regrid`, `regrid!` or `plan_regrid` method, no target-conversion generic and
  no output labelling outside `destinationdims`, so a destination-geometry
  change has no package-side API to keep in step.
- Destination cells reach a block builder as `dst_inds` on a space, through
  `subtree` and the destination polygons; nothing on the API surface reads
  them. `DGGSpace.subtree` is the one hook a prepared destination tile would
  consume from this side.
- The `_asspace` methods construct a `DGGSpace` and nothing else, so a prepared
  destination cannot be built at target resolution time and does not need to
  be: planning reads no geometry.
