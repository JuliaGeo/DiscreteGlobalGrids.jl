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
The [architecture guide](../architecture.md) explains the design.

## What a grid answers

Four of these are the whole implementor's contract — [`ncells`](@ref),
[`cellindex`](@ref), [`cell_boundary`](@ref) and [`cell_centroid`](@ref) — and
the rest have working defaults built on them.

```@docs
ncells
levelgrid
cell_boundary
cell_centroid
cell_polygon
cell_area
cell_extent
cellat
system
level
```

## Cell size, and choosing a level

`cellsize` measures a typical cell width in metres. `levelfor` finds the
closest available level for a requested width or another dataset's resolution.

```@docs
cellsize
levelfor
```

## The hierarchy

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

[`treeify`](@ref) exposes a grid as a spatial tree for queries and regridding.
Hierarchical grids can use their existing parent/child structure and compute
node geometry as the traversal needs it.

```@docs
treeify
getcell
DiscreteGlobalGrids.Engine.node_cell
DiscreteGlobalGrids.Engine.node_indices
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

## System capabilities

These declarations let generic algorithms use a system's hierarchy, cell
ordering and point-location methods efficiently.

```@docs
has_sorted_subtrees
has_congruent_refinement
has_direct_location
node_extent
```

## Geometry and traversal helpers

These helpers support grid implementations: coordinate transforms, boundary
rings, spatial bounds and traversal engines.

```@docs
DiscreteGlobalGrids.one_ring
DiscreteGlobalGrids.cap_inflation
DiscreteGlobalGrids.authalic_sphere
DiscreteGlobalGrids.Fallbacks.authalic_stretch
DiscreteGlobalGrids.Fallbacks.authalic_shift
DiscreteGlobalGrids.Helpers.AuthalicTransform
DiscreteGlobalGrids.Helpers.authalic_radius
DiscreteGlobalGrids.Helpers.EllipsoidShapeError
DiscreteGlobalGrids.Fallbacks.closed_ring
DiscreteGlobalGrids.border_engine
DiscreteGlobalGrids.lattice_decode
Trees.AbstractCurvilinearGrid
Trees.cell_range_extent
```

## Index

```@index
Pages = ["api/grid-interface.md"]
```
