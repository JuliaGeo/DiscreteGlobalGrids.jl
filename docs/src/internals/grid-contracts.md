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

A tree node also answers the hooks GridInterface's generic verbs read;
[Trees](../architecture.md#Trees) lists what each cursor returns.

| hook | answers |
|---|---|
| `treeify(grid, n::Integer)` | the grid's own tree with up to `n` cells per leaf; the one builder every intent and `subtree(grid, cells)` reach |
| `GridInterface.position_span(node)` | a bound on the grid positions beneath the node, which the generic `subtree` prunes by |
| `SpatialTreeInterface.node_extent_is_expensive(::Type{<:Node})` | whether a dual-tree sweep should carry child extents across visits |

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
```
