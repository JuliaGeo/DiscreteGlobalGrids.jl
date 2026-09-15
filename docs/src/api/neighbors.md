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

`mapneighbors` applies a kernel to each cell and collects the results.
`mapneighbors!` writes into an existing destination, and `foreachneighbors`
runs a callback without collecting its return values.

Use `pass = Values()` to receive the centre value and neighbouring values.
For kernels that also need geometry, indices or multiple variables, see
[Requesting neighbour fields](neighbor-fields.md). For stored data, see
[chunked sweeps](chunk-sweep.md).

The `neighborhood` keyword chooses which cells each visit hands the kernel.
`Disc(k)` hands it `neighbors(grid, cell, k)`, every cell within `k` steps;
`Ring(k)` hands it `ring(grid, cell, k)`, the cells at exactly `k` steps. The
default is `Disc(1)`, the one-ring. Every callback form keeps its arity: only
the ring argument widens. The same keyword drives the one-argument
[`neighbors`](@ref) iterator and the chunked [`mapneighbors!`](@ref), whose
halo follows the radius.

```@example neighbors
cells = DGG.CellVector(grid)
i = DGG.localindex(cells, cell)

within3 = DGG.mapneighbors((c, nbrs) -> [DGG.localindex(h) for h in nbrs], cells;
    neighborhood = DGG.Disc(3))
at3 = DGG.mapneighbors((c, nbrs) -> [DGG.localindex(h) for h in nbrs], cells;
    neighborhood = DGG.Ring(3))

(; disc_is_neighbors = within3[i] == DGG.neighbors(cells, i, 3),
   ring_is_ring = at3[i] == DGG.ring(cells, i, 3),
   one_ring_leads = within3[i][1:length(DGG.neighbors(cells, i))] == DGG.neighbors(cells, i))
```

The order is the one [`neighbors`](@ref) and [`ring`](@ref) fix: a `Disc(k)`
ring holds ring 1, then ring 2, out to ring `k`, each counter-clockwise, and a
subset drops its non-members in place. The ring boundaries carry no marker, so
a kernel that weights by distance sweeps `Ring(j)` once per `j`.

One pass at radius `k` is a different computation from `k` one-ring passes.
The disc pass weights every cell within `k` steps once; repeated one-ring
passes compound, reaching a cell at distance two through every path of length
two, so their weights fall off with distance.

```@example neighbors
mean3(c, x, nbrs) = (x + sum(nbrs)) / (1 + length(nbrs))
values = sin.(eachindex(cells) ./ 7)

radius3 = DGG.mapneighbors(mean3, cells, values; neighborhood = DGG.Disc(3))
threepasses = foldl((v, _) -> DGG.mapneighbors(mean3, cells, v), 1:3; init = values)
maximum(abs, radius3 .- threepasses)
```

`Disc(0)` is legal and hands every cell an empty ring; `Ring(0)` throws,
because the sweep already hands the centre to the kernel as its first
argument. On systems whose [`winding`](@ref) is `CounterClockwise` or
`Clockwise` a `k ≥ 2` visit is one shell walk, and H3 answers from libh3's
own automaton, both on the stack, so a sequential sweep through `Disc(3)` and
`Ring(3)` allocates nothing there; on `CustomOrder` or undeclared-winding
systems it sorts each
shell by centroid azimuth, so a wide sweep there costs more than `k` one-ring
sweeps. [`adjacency`](@ref) stays a one-ring table; widen the
region with [`grow`](@ref) when a wider table is wanted.

```@docs
mapneighbors
mapneighbors!
foreachneighbors
Values
Neighbors
Neighborhood
Disc
Ring
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
