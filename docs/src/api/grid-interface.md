# The grid interface

```@meta
CurrentModule = DiscreteGlobalGrids
```

Every verb on this page is asked of one of two things, and
[Architecture](../architecture.md) says why they are disjoint: an
[`AbstractGrid`](@ref) is a finite collection of cells with a length, and an
[`AbstractHierarchicalGridSystem`](@ref) is a rule for tessellating the sphere
with no length at all. [`levelgrid`](@ref) is the bridge — a complete level of a
system, answered as a grid.

A cell is named twice over. A bare `Int` is always a *local index* into
`1:ncells(grid)`; an [`AbstractCellIndex`](@ref) is a typed identity that carries
its own level. The two are separated by method signature alone, so a verb that
takes either takes both.

## What a grid answers

Four of these are the whole implementor's contract — [`ncells`](@ref),
[`cellindex`](@ref), [`cell_boundary`](@ref) and [`cell_centroid`](@ref) — and
the rest have working defaults built on them.

```@docs
ncells
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

The two directions of one question: how coarse is this level, and which level is
that coarse.

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
```

## Trees over a grid

[`treeify`](@ref) is how the regridding and point-location machinery walks a
grid: it hands back a `ConservativeRegridding.Trees` cursor whose nodes are
subtrees of the system's own hierarchy, so no tree is built and no geometry is
touched until a query descends.

```@docs
treeify
getcell
DiscreteGlobalGrids.Engine.node_cell
DiscreteGlobalGrids.Engine.node_indices
```

## Identifiers

```@docs
cellid
rawid
reindex
cellindextype
cellindextypes
```

## Neighbour counts and ring order

The bounds and the ordering declaration the neighbourhood machinery reads.
[`neighbors`](@ref) and [`ring`](@ref) themselves are documented with the
[region boundaries](boundaries.md).

```@docs
neighborcount
maxring
winding
Winding
CounterClockwise
Clockwise
CustomOrder
Unordered
```

## Grids that stand for a system

No shipped system defines a grid type of its own: [`levelgrid`](@ref) returns a
[`HierarchicalLevelGrid`](@ref), and [`AuthalicSystem`](@ref) wraps any system so
that geometry reads at geodetic latitude, leaving ids, indices, hierarchy and
ordering exactly as they were. The concrete systems themselves are drawn and
compared in the [gallery](../all_dggs.md).

```@docs
HierarchicalLevelGrid
AuthalicSystem
AuthalicGrid
```

## The types the contract is stated in

```@docs
AbstractGrid
AbstractHierarchicalGridSystem
AbstractQuadFaceGridSystem
AbstractCellIndex
```

## The internals these link to

Listed for the same reason as the corresponding section of
[Region boundaries](boundaries.md): the docstrings above cross-reference them,
and these entries are what make the links resolve.

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
