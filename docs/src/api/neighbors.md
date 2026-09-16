# Neighbours and stencils

```@meta
CurrentModule = DiscreteGlobalGrids
```

Use this API to find cells around a cell, build an adjacency table, or compute
new values from a neighbourhood. For the edge of an entire region, use
[Region boundaries](boundaries.md). The [stencil tutorial](../tutorials/stencils.md)
works through smoothing, edge detection and graph traversal.

## Find neighbours around a cell

`neighbors(grid, cell, k)` returns cells within `k` adjacency steps, excluding
the centre. `ring(grid, cell, k)` returns only the cells at step `k`. Each ring
runs counter-clockwise as seen from outside the sphere.

```@example neighbors
import DiscreteGlobalGrids as DGG

grid = DGG.levelgrid(DGG.HEALPixSystem(), 3)
cell = DGG.cellat(grid, 8.5, 47.4)

(; first_ring = length(DGG.neighbors(grid, cell)),
   second_ring = length(DGG.ring(grid, cell, 2)),
   within_two = length(DGG.neighbors(grid, cell, 2)))
```

Pass a cell id to receive cell ids, or a local integer index to receive local
indices. On a subset, results include only members of that subset.

```@docs
neighbors
ring
```

## Choose connectivity

The default `Vertex()` includes cells touching at a corner. `Edge()` requires
a shared edge. These choices often coincide on hexagonal grids and differ on
quadrilateral grids.

```@docs
Connectivity
Vertex
Edge
neighborcount
```

## Build an adjacency table

`adjacency(region)` caches a neighbour list for each cell. Row `i` contains
local indices for the neighbours of cell `i`, in ring order. Reuse the table
when an algorithm repeatedly traverses the same cells.

The `halo` keyword controls neighbours outside the region:

| Value | Row entries |
|---|---|
| `0` (default) | Only neighbours belonging to the region |
| `1` | Indices into a combined region-and-halo buffer |
| `:mark` | Complete neighbour slots, with `0` for an outside neighbour |

For `halo = 1`, `haloindices(table)` and `halocells(table)` identify the cells
in the buffer's halo portion. See [Region boundaries](boundaries.md) for the
relationship between that halo and the region's border.

```@docs
adjacency
AdjacencyTable
halocells
haloindices
```

## Compute with neighbourhoods

`mapneighbors` applies a one-ring kernel to each cell and collects the results.
`mapneighbors!` writes into an existing destination, and `foreachneighbors`
runs a callback without collecting its return values.

All three sweeps currently supply one-ring neighbors. They do not accept a
radius-k or exact-ring selector. A wider chunk halo does not change that limit.
Repeated one-ring passes perform diffusion, not one radius-k convolution.

For a dimensional array, choose the callback shape with `pass`:

| Mode | Callback | Use it when |
| --- | --- | --- |
| `Neighbors()` | `f(cell, neighbors)` with indexed cell handles | The callback needs cell identity or controls data access |
| `Values()` | `f(cell, value, neighbors)` with values at one non-spatial position | Each time or band slice needs the same scalar kernel |
| `NeighborSlices()` | `f(cell, center_slice, neighbor_slices)` | One call needs each cell's whole time series or other non-spatial slice |

`Values()` preserves the input dimensions and visits each non-spatial position.
`NeighborSlices()` makes one call per cell and collects one return value per cell.
See its reference below for the exact slice and output contracts.

```@example neighbor-modes
import DiscreteGlobalGrids as DGG
import DimensionalData as DD

grid = DGG.levelgrid(DGG.HEALPixSystem(), 1)
lookup = DGG.CellLookup(grid)
values = [Float64(t + c) for t in 1:3, c in 1:length(lookup)]
cube = DD.DimArray(values, (DD.Dim{:time}(1:3), DGG.Cells(lookup)))

per_time = DGG.mapneighbors((c, x, ns) -> x + sum(ns), cube;
    pass=DGG.Values(), threaded=false)
per_cell = DGG.mapneighbors((c, xs, ns) -> sum(xs) + sum(sum, ns), cube;
    pass=DGG.NeighborSlices(), threaded=false)
@assert size(per_time) == size(cube)
@assert size(per_cell) == (length(lookup),)
@assert parent(per_cell) ≈ vec(sum(parent(per_time); dims=1))
(; per_time=size(per_time), per_cell=size(per_cell))
```

Use `pass = Values()` to receive the centre value and neighbouring values.
For kernels that also need geometry, indices or multiple variables, see
[Requesting neighbour fields](neighbor-fields.md). For stored data, see
[chunked sweeps](chunk-sweep.md).

```@docs
mapneighbors
mapneighbors!
foreachneighbors
Values
Neighbors
NeighborSlices
DiscreteGlobalGrids.StorageOrder
DiscreteGlobalGrids.NeighborCallbackError
```

## Neighbours across levels

`member_neighbors` finds adjacent members of a mixed-level `MultiOrderCellSet`.
Use it when the cells themselves have different resolutions; ordinary
`neighbors` queries a collection at one level.

```@docs
member_neighbors
```

## Neighbour bounds and ordering declarations

Grid implementations use these declarations to describe neighbourhood size and
ring orientation. The query functions above handle the resulting traversal.

```@docs
maxneighbors
maxring
winding
Winding
CounterClockwise
Clockwise
CustomOrder
Unordered
```

## Index

```@index
Pages = ["api/neighbors.md"]
```
