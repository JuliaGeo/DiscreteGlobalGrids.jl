"""
    DiscreteGlobalGridVisualizations

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
module DiscreteGlobalGridVisualizations

import DiscreteGlobalGrids as DGG
import GeometryBasics
using GeometryBasics: Point2d, Point3d, GLTriangleFace
import Makie
using Makie: @recipe, on

export dggpoly, dggpoly!
export dggresample, dggresample!

include("targets.jl")
include("cellsets.jl")
include("tessellate.jl")
include("recipe.jl")
include("pyramid.jl")
include("resample.jl")
include("resample_recipe.jl")

end # module
