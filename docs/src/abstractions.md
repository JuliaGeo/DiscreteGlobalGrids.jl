# Grids and cell indices

DiscreteGlobalGrids separates a grid's geometry from the values stored on it.
These concepts also apply to regional grids and grids without a hierarchy.

## Grids and systems

A **grid** is a finite collection of cells. An H3 grid at level 12 is a grid;
so is a selected region of that level. Geometry, point queries and
neighbourhood operations act on grids.

A **system** defines related grids at several levels, including their parent
and child relationships. `levelgrid(sys, level)` selects one complete level.
A grid can also supply its geometry directly, as an ocean-model grid might.

The [grid interface](api/grid-interface.md) documents both abstractions.
[Writing a grid system](extending.md) implements a complete example, and the
[architecture guide](architecture.md) describes the generic algorithms.

## Cell ids and array positions

A **cell id** identifies a cell independently of the array that stores its
value. A **local index** identifies a position in one collection. Selecting
a region changes local positions while preserving cell ids.

A `Cells` dimension connects those identities to array values. Use it to
[select data by point or region](api/selecting-cells.md), while ordinary
integer indexing selects positions in the current array.
