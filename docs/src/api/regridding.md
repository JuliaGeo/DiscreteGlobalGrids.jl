# Regridding calls and plans

```@meta
CurrentModule = DiscreteGlobalGrids
```

Use `regrid` for one transfer, `plan_regrid` to reuse the same spatial mapping,
and `regrid!` to write into an existing destination. These functions belong to
GlobalRegridding and are re-exported by DiscreteGlobalGrids.
Choose the meaning of the transfer in [Choosing a regridding method](regridding-methods.md).

## One transfer

This example transfers a constant field between two HEALPix levels. A plain
vector needs an explicit source grid because it carries no spatial metadata.

```@example regridding-api
import DiscreteGlobalGrids as DGG

source = DGG.levelgrid(DGG.HEALPixSystem(), 1)
target = DGG.levelgrid(DGG.HEALPixSystem(), 2)
values = fill(3.0, DGG.ncells(source))
result = DGG.regrid(values; from=source, to=target, method=DGG.NearestCell())
@assert length(result) == DGG.ncells(target)
@assert all(==(3.0), result)
size(result)
```

A dimensional source can describe its own spatial axes. DGGS targets accept
complete or partial grids, cell vectors, cell lookups, and mixed-level sets.
A bare system as `to` chooses a level from the source resolution. Use an
explicit grid when the level must be fixed.

Source spatial dimensions must come first and follow the source space's cell
order. Other dimensions follow them in their original order. Dimensional
results use the destination's axes; DGGS destinations provide a `Cells` axis.
Plain array inputs return arrays. Results are floating point, including
containing-cell transfers of integer categories.

## Reuse a plan or destination

A plan fixes the spaces, cell order, method, and weight policy. Reuse it for
fields with that same spatial layout; changing geometry or cell order requires
a new plan.

```@example regridding-api
plan = DGG.plan_regrid(values; from=source, to=target,
    method=DGG.NearestCell())
dest = zeros(DGG.ncells(target))
DGG.regrid!(dest, 2 .* values, plan)
@assert dest == 2 .* result
sum(dest)
```

`regrid(values, plan)` allocates a result. `regrid!(dest, values, plan)` uses
existing storage. The destination must have the expected shape and an element
type that can hold its missing-value sentinel. A dimensional source can contain
time or band slices; one plan serves those slices without rebuilding weights.

## Lazy execution and memory

| Control | Meaning |
| --- | --- |
| `lazy=true` | Compute destination data on demand; defaults to true for chunked sources |
| `chunks` | Lazy output tiling; separate from the spatial chunks in `DGGSpace` |
| `budget` | Target bytes for lazy reads and weights, not a strict total-process memory limit |
| `storage=PerChunk()` | Cache lazy weight blocks; this explicit form has no cache size limit |
| `storage=Spilled(directory)` | Store lazy weight blocks on disk |
| `DGGSpace(grid; chunklevel, chunkcells)` | Control the DGGS space's ancestor chunks |

Without an explicit storage policy, lazy plans use a budget-limited cache.
`chunks`, `budget`, and `storage` apply to lazy execution. `sampling` controls
lookup sampling for eager execution. See [Assigning chunks to workers](partitioning.md)
for partitioning an existing lazy plan's dependencies.

Missing-value normalization and output sentinels are separate choices; read
[Missing data](regridding-methods.md#missing-data) before reusing a plan across
fields with different coverage.

## Callable reference

```@docs
regrid
regrid!
plan_regrid
DGGSpace
```

## Point methods and weight storage

`Conservative`, `Weighted`, and `Extensive` are documented on the
[method-choice page](regridding-methods.md).

```@docs
NearestCell
DirectNearest
BarycentricPoint
PerChunk
Spilled
```

## Index

```@index
Pages = ["api/regridding.md"]
```
