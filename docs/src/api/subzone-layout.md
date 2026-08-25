# The ancestor-subzone layout

!!! warning "Interim"
    Zarr chunks are uniform by format, so a chunk grid that follows the tree
    exactly is not expressible on a one-dimensional cell axis: an IGEO7
    pentagon's subtree holds `p(d) = (5·7^d + 1)/6` cells where a hexagon's
    holds `7^d`. This layout works around that and goes away once Zarr supports
    variable chunk sizes.

The layout buys tree-aligned chunking by spending a dimension on it. Data arrays
are two-dimensional — subzone index within one ancestor's subtree, then the
ancestor — chunked one column per ancestor, so a chunk *is* a subtree. Index
within a column is what OGC API-DGGS calls the **sub-zone order**, here ascending
cell id. So:

  - an ancestor nobody wrote is a chunk that was never stored, so a land-only
    global store costs nothing for the ocean and reads back as `fill_value`;
  - a column is one file, so a production run writes columns from as many tasks
    as it likes, with no coordination and no shared file to rewrite;
  - a reader gets the tree's own irregular chunking back.

Written with `dggwrite(dest, cube; layout = :subzones, ancestor_level = k)` or
incrementally through [`subzonestore`](@ref) and [`dggwrite!`](@ref); read by
[`dggread`](@ref) like any other store, which hands back a `Cells` dimension over
a lazy `DiskArrays` view of the two-dimensional arrays, with the pentagon padding
dropped and the subtree chunk boundaries published.

A column is written whole: a cube whose coverage stops inside a subtree is
refused with `DGGSFormatError(check = :incomplete_subtree)`.

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

A `dggs` attribute object carrying the grid name and the level, with everything
the layout adds nested under `subzone_layout`. It declares no `zarr_conventions`:
this is not the one-dimensional layout that convention describes.

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
