# Grid extension reference

```@meta
CurrentModule = DiscreteGlobalGrids
```

Use these contracts when implementing a grid or optimizing a system.
[Writing a grid system](../extending.md) provides a complete example.
User operations remain in [The grid interface](../api/grid-interface.md).

## Tree node access

```@docs
DiscreteGlobalGrids.Engine.node_cell
DiscreteGlobalGrids.Engine.node_indices
```

## Hierarchy traits

These declarations let generic algorithms use a system's hierarchy, cell
ordering and point-location methods efficiently.

```@docs
has_sorted_subtrees
has_congruent_refinement
has_direct_location
node_extent
```

## Geometry and traversal hooks

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
