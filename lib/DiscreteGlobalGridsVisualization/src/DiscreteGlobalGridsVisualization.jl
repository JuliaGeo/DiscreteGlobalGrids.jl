"""
    DiscreteGlobalGridsVisualization

Fast Makie plotting for discrete global grids.

The package draws a set of DGGS cells as **one mesh** rather than as one polygon
per cell.  A cell boundary is already a small convex ring, so it needs no
general-purpose triangulator: fanning it from its first vertex is exact, costs
`n - 2` triangles, and can be done for every cell in parallel.  The result is a
single `mesh!` — one GPU buffer, one draw call — instead of `ncells` polygons
that Makie has to tessellate, concatenate and re-upload whenever anything
changes.

The entry point is [`dggpoly`](@ref), which reads like Makie's `poly`:

```julia
dggpoly(cells; color = values)
```

`cells` is anything that names a set of DGGS cells — an `AbstractGrid`, a
`CellVector`, a `CellLookup`, a `MultiOrderCellSet`, or a `(system, ids)` pair.

## Drawing the field instead of the cells

[`dggsurface`](@ref) takes the same arguments and draws the same cells as the
continuous field they sample: a vertex at each cell's **centroid**, carrying
that cell's value, joined by the triangles of the grid's dual.  The value varies
smoothly between cell centres rather than jumping at cell edges — the picture
for elevation or temperature, not for a categorical field.

It is the smaller mesh, one vertex and about two triangles per cell against six
and four.  Its triangles come from adjacency alone, under a rule that gives each
grid corner one owner, so none is emitted twice and none has to be de-duplicated
afterwards; see `surface.jl`.

It also takes heights, the way Makie's own `surface` does — `dggsurface(cells,
zs)` lifts the vertex over each cell to that cell's height, as a third
coordinate on a flat map and as a height above the ellipsoid on a globe.  A
one-dimensional `DimArray` over a cell dimension carries both, so
`dggsurface(A)` draws `A` as relief.

## Drawing less than you were given

[`dggresample`](@ref) takes the same arguments and answers the other half of the
question.  Past a certain level every extra cell lands under a pixel another
cell already owns, so instead of drawing the set it descends the system's own
hierarchy — keeping only branches that are on screen and hold data — and stops
at the level whose cells come out a few pixels across, colouring each of them by
the leaf cell under its centre.  It follows the camera, so zooming in refines
and zooming out coarsens, and neither costs anything proportional to the number
of cells handed in.

Given heights as well, `dggresample(cells, zs)` draws the frame as a surface
rather than as patches, so relief has the same answer to "more cells than
pixels" that a flat field does.

## Taking the mesh away

`GeometryBasics.mesh(cellregion(cells)[, elevation]; attributes...)` builds the
same mesh outside any plot, as a `GeometryBasics.Mesh` — to write out as a mesh
file, to hand to something else, or to plot with Makie's own `mesh`.  Values
arrive per cell and leave per vertex.  `cellset` in place of `cellregion` gives
the flat patches instead.

## Where the mesh lives

A DGGS cell is a spherical polygon, and where its corners land on screen depends
on what the axis does with them.  The recipe therefore reads the axis's own
transformation off the plot (`:transform_func`) and turns it into a
[`PlotTarget`](@ref):

  * [`GlobeTarget`](@ref) — a `GeoMakie.GlobeAxis`.  The mesh is built in the
    axis's 3D space directly from the unit-sphere corners, and no cell has to be
    split, because the sphere has no cut.
  * [`PlanarTarget`](@ref) — a plain `Axis` (identity transform) or a
    `GeoAxis` (a `Proj.Transformation`).  The mesh is built in longitude and
    latitude, cells that straddle the map's cut meridian are split against it,
    and the whole vertex buffer is then projected in one bulk call.

Both are ordinary types with ordinary methods, so a new axis kind is a new
[`plot_target`](@ref) method, and a new grid system needs nothing at all: cell
geometry is read through `DiscreteGlobalGrids.cell_boundary`, which every system
implements.

!!! warning "Proof of concept"
    This package is an experiment living in `lib/`.  Names and behaviour are not
    stable, and the pieces that earn their keep are meant to move into
    `DiscreteGlobalGrids`' own Makie extension.
"""
module DiscreteGlobalGridsVisualization

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
import GeometryBasics
using GeometryBasics: Point2d, Point3d, GLTriangleFace
import Makie
using Makie: @recipe, on

export dggpoly, dggpoly!
export dggsurface, dggsurface!
export dggresample, dggresample!

include("chunks.jl")
include("targets.jl")
include("cellsets.jl")
include("tessellate.jl")
include("recipe.jl")
include("surface.jl")
include("surface_recipe.jl")
include("export.jl")
include("pyramid.jl")
include("resample.jl")
include("resample_recipe.jl")
include("dimensionaldata.jl")

end # module
