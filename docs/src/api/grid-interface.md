# The grid interface

```@meta
CurrentModule = DiscreteGlobalGrids
```

Use this interface to locate cells, read their geometry and navigate between
levels. An [`AbstractGrid`](@ref) is a finite collection of cells; an
[`AbstractHierarchicalGridSystem`](@ref) defines grids and their hierarchy.
[`levelgrid`](@ref) returns one complete level of a system.

A bare `Int` addresses a local position in `1:ncells(grid)`. An
[`AbstractCellIndex`](@ref) identifies a cell and carries its level. Keep
that distinction when working with subsets, where local positions change.
[Grids and cell indices](../abstractions.md) introduces these concepts.
The [architecture guide](../architecture.md) explains the design.

## What a grid answers

Locate a cell with longitude and latitude in degrees. Geometry methods return
unit-sphere points, and `cell_area` returns steradians. Convert a point to
geographic coordinates when the next consumer expects longitude and latitude:

```@example grid-basics
import DiscreteGlobalGrids as DGG
import GeometryOps as GO
sys = DGG.HEALPixSystem()
grid = DGG.levelgrid(sys, 2)
cell = DGG.cellat(grid, 8.5, 47.4)
togeographic = GO.GeographicFromUnitSphere()
centroid_lonlat = togeographic(DGG.cell_centroid(grid, cell))
boundary_lonlat = togeographic.(DGG.cell_boundary(grid, cell))
(; centroid_lonlat, boundary_points=length(boundary_lonlat))
```

The converted pairs are `(longitude, latitude)` in degrees, with longitude
in `[-180, 180]`. Coordinate conversion does not change the latitude frame;
see [Choosing a grid](../tutorials/choosing_a_grid.md) for authalic and geodetic data.

```@docs
ncells
levelgrid
cell_boundary
cell_corners
cell_centroid
cell_polygon
cell_area
cell_extent
cellat
system
level
```

## Cell size, and choosing a level

`cellsize` measures a typical cell width in metres; `levelfor` finds the
closest level for a requested width or another dataset's resolution. Both take
an `over` area of interest, within which `levelfor` measures both sides.

```@docs
cellsize
levelfor
```

## The hierarchy

Larger level numbers mean finer cells. Choose the output needed by your task:

| Operation | Result |
| --- | --- |
| `parent(sys, cell)` / `children(sys, cell)` | Immediate coarser / finer relatives |
| `ancestor(sys, cell, l)` | One ancestor at `l <= level(cell)` |
| `descendants(sys, cell, l)` | All descendants at `l >= level(cell)` |
| `descendant_range(sys, cell, l)` | Their complete-level positions, only with sorted subtrees |
| `subtree(sys, cell, l)` | A regional grid holding those descendants |

```@example grid-basics
coarse = parent(sys, cell)
leaves = DGG.descendants(sys, coarse, 3)
positions = DGG.descendant_range(sys, coarse, 3)
regional = DGG.subtree(sys, coarse, 3)
@assert length(leaves) == length(positions) == DGG.ncells(regional)
(; parent_level=DGG.level(coarse), fine_cells=length(leaves))
```

`descendant_range` requires `has_sorted_subtrees(sys)`. A5 does not provide
that property. Hierarchical descendants need not geometrically tile their
parent; that stronger property is `has_congruent_refinement(sys)`.
See [Region boundaries](boundaries.md) to compute on a subtree region.


```@docs
levels
maxlevel
rootcells
children
Base.parent(::AbstractHierarchicalGridSystem, ::AbstractCellIndex)
ancestor
descendants
descendant_range
```

## Trees over a grid

[`treeify`](@ref ConservativeRegridding.Trees.treeify) exposes a grid as a spatial tree for queries and regridding.
Hierarchical grids can use their existing parent/child structure and compute
node geometry as the traversal needs it.

```@docs
treeify
getcell
```

## Identifiers

```@docs
cellindex
localindex
globalindex
LevelIndex
cellid
rawid
reindex
cellindextype
cellindextypes
```

## Grids that stand for a system

[`levelgrid`](@ref) returns a [`HierarchicalLevelGrid`](@ref).
[`AuthalicSystem`](@ref) adapts supported systems to geodetic latitude while
preserving cell ids and hierarchy. See [Choosing a
grid](../tutorials/choosing_a_grid.md) for coordinate guidance and the
[gallery](../all_dggs.md) for the available systems.

```@docs
HierarchicalLevelGrid
AuthalicSystem
AuthalicGrid
```

## The types the contract is stated in

```@docs
AbstractGrid
PartialGrid
AbstractHierarchicalGridSystem
AbstractQuadFaceGridSystem
AbstractCellIndex
```

Implementation traits and traversal helpers are in the
[Grid extension reference](../internals/grid-contracts.md).

## Index

```@index
Pages = ["api/grid-interface.md"]
```
