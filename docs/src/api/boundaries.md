# Region boundaries

```@meta
CurrentModule = DiscreteGlobalGrids
```

Use these operations to find the cells along a region's edge, fetch values just
outside it, or enlarge it by a few layers of cells. A region is a collection at
one level, such as a [`PartialGrid`](@ref), [`CellVector`](@ref) or
[`CellLookup`](@ref).

| Operation | Cells returned |
|---|---|
| [`border`](@ref) | Cells in the region with a neighbour outside it |
| [`interior`](@ref) | Cells in the region whose neighbours all belong to it |
| [`halo`](@ref) | Cells outside the region that touch it, at the same level |
| [`grow`](@ref) | The region plus a requested number of neighbouring layers |

`border` and `interior` partition the region. Its halo supplies the outside
values needed to compute a stencil at the border. Holes also have boundaries:
cells just inside a hole can belong to the halo. A complete global grid has an
empty border and halo, and all of its cells belong to the interior.

These operations return cells or their indices. To obtain polygon coordinates,
use [`cell_boundary`](@ref) for a single cell or [`cell_polygons`](@ref) for a
collection. For neighbours around individual cells and adjacency tables, see
[Neighbours and stencils](neighbors.md).

## Find the inside and outside of an edge

This example uses the descendants of one IGeo7 cell as a region. The border
and interior divide its cells; the halo contains the touching cells outside it.

```@example boundaries
import DiscreteGlobalGrids as DGG

sys = DGG.IGeo7System()
root = first(DGG.CellVector(DGG.levelgrid(sys, 1)))
region = DGG.subtree(sys, root, 3)

edge = collect(DGG.border(region; cells = true))
inside = collect(DGG.interior(region; cells = true))
outside = collect(DGG.halo(region; cells = true))

(; border = length(edge), interior = length(inside), halo = length(outside),
   region = DGG.ncells(region))
```

`connectivity = Edge()` requires a shared edge; the default `Vertex()` also
counts corner contact. Use the same connectivity for the border, halo and
stencil that will consume them.

## Use the indices with data

`border`, `interior` and `halo` return lazy iterators. By default they yield
indices in ascending order; `cells = true` requests cell ids.

| Operation | Default index space |
|---|---|
| `border(region)` and `interior(region)` | Local positions in the region's data array |
| `halo(region)` | Positions in the complete level grid, where the outside cells live |

Use local border indices to read or update values stored on the region. Use
halo indices to fetch context from the complete grid. For data on disk,
[chunked sweeps](chunk-sweep.md) manage these reads and index translations.

`collect(iterator)` materialises the indices or ids; `Set(iterator)` makes a
membership set. Some iterators have no cheap exact `length`, so collect them
when you need to count the results. Contributor details are in
[Boundary traversal engines](../internals/boundary-engines.md).

```@docs
halo
border
interior
```

## Enlarge a region

`grow(region, n)` includes the original region and `n` layers of neighbours.
Use it when a calculation needs a wider margin around an area of interest.

```@example boundaries
grown = DGG.grow(region, 1)
length(grown) == DGG.ncells(region) + length(outside)
```

Growth uses a single level. Expand a [`MultiOrderCellSet`](@ref) to the desired
level with [`expand`](@ref) before applying it.

```@docs
grow
```

## Construct a region from an ancestor

[`subtree`](@ref) selects all descendants of one cell at a target level. It is
useful for regions aligned with the grid hierarchy, as in the example above.

```@docs
subtree
```

## Index

```@index
Pages = ["api/boundaries.md"]
```
