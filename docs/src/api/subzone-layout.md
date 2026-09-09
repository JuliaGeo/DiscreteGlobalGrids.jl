# The ancestor-subzone layout

Use this layout to group cells by ancestor, so each storage chunk covers a
complete subtree. It stores data in a two-dimensional array: subzone
index within an ancestor's subtree, then the ancestor. One column is one
subtree, and its row order is the OGC API-DGGS **sub-zone order**, here ascending
cell id. This gives three useful properties:

  - an ancestor nobody wrote is a chunk that was never stored, so a land-only
    global store costs nothing for the ocean and reads back as `fill_value`;
  - a column is one file, so production tasks can write columns independently;
  - a reader gets the tree's own irregular chunking back.

IGEO7 pentagon subtrees contain fewer cells than hexagon subtrees. The layout
pads the shorter columns to fit Zarr's uniform chunks, and the reader removes
that padding from the cell axis.

Write it with `dggwrite(dest, cube; layout = :subzones, ancestor_level = k)` or
incrementally through [`subzonestore`](@ref) and [`dggwrite!`](@ref). Read it
with [`dggread`](@ref), which returns a `Cells` dimension over a lazy
`DiskArrays` view, drops pentagon padding, and publishes subtree chunk
boundaries.

Each column is written whole. A cube whose coverage stops inside a subtree
raises `DGGSFormatError(check = :incomplete_subtree)`.

```@docs
SubzoneLayout
subzonestore
dggwrite!
DiscreteGlobalGrids.SubzoneRun
DiscreteGlobalGrids.subzone_capacity
DiscreteGlobalGrids.subzone_depth
DiscreteGlobalGrids.subzone_runs
DiscreteGlobalGrids.subzone_cellvector
DiscreteGlobalGrids.subzone_columns
DiscreteGlobalGrids.columncell
DiscreteGlobalGrids.columnindex
DiscreteGlobalGrids.columnindices
DiscreteGlobalGrids.columnlength
DiscreteGlobalGrids.subzoneindex
DiscreteGlobalGrids.columnrow
```

## Attributes

The `dggs` attribute object carries the grid name and level, with layout fields
nested under `subzone_layout`. It omits `zarr_conventions` because that field
describes the one-dimensional layout.

```@docs
DiscreteGlobalGrids.subzone_attrs
DiscreteGlobalGrids.subzone_layout
DiscreteGlobalGrids.issubzonestore
DiscreteGlobalGrids.subzone_coordinate
DiscreteGlobalGrids.SUBZONE_LAYOUT
DiscreteGlobalGrids.SUBZONE_ORDER
DiscreteGlobalGrids.SUBZONE_PADDING
DiscreteGlobalGrids.gridnamefor
```
