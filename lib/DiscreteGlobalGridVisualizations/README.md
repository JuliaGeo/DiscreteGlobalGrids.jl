# DiscreteGlobalGridVisualizations

Fast Makie plotting for discrete global grids — a proof of concept for
[JuliaGeo/DiscreteGlobalGrids.jl#12](https://github.com/JuliaGeo/DiscreteGlobalGrids.jl/issues/12).

```julia
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGridVisualizations, GLMakie, GeoMakie

sys = DGG.IGeo7System()
region = DGG.query(sys, DGG.MultiOrderCoverage(extent); level = 12)
cells = DGG.CellVector(region)

dggpoly(cells; color = elevation)                       # a plain Axis
dggpoly!(GeoAxis(f[1, 1]; dest = "+proj=moll"), cells)  # a projected map
dggpoly!(GlobeAxis(f[1, 1]), cells)                     # a globe
```

`dggpoly` reads like `poly`: cells in, one filled patch per cell out, coloured
by a vector as long as the cell set.  It accepts anything that names a set of
DGGS cells — an `AbstractGrid`, a `CellVector`, a `CellLookup`, a
`MultiOrderCellSet`, or a system with a vector of ids — so a partial grid is the
normal case rather than a special one.

## Why it is faster

`poly` on `n` cells builds `n` polygons, triangulates each with a
general-purpose triangulator, and concatenates the results; it also builds the
outline of every cell whether or not a stroke is drawn.  All of that is redone
whenever anything about the plot changes.

A DGGS cell is a convex ring of five to ten corners, so it needs no triangulator
at all — fanning it from its first corner is exact and costs `n - 2` triangles —
and no cell depends on any other, so the whole job is a parallel loop into
disjoint buffers.  The result is one mesh: one buffer, one draw call.

Colour is kept separate from geometry.  Each vertex records the cell it came
from, so a new colour vector is one gather and the geometry never moves.

## Where the mesh lives

Where a cell's corners land on screen is the axis's business, and the axis says
so in one place: its transform function.  The recipe reads that off the plot and
turns it into a *plot target*:

| axis | transform function | target |
| --- | --- | --- |
| `Axis` | `(identity, identity)` | `PlanarTarget`, cut at the antimeridian |
| `GeoAxis` | `Proj.Transformation` | `PlanarTarget`, cut read from the projection |
| `GlobeAxis` | `GeoMakie.GlobeTransform` | `GlobeTarget` |

* **On a globe** there is no seam and no projection to apply: a unit-sphere
  corner becomes an earth-centred vertex with three multiplies and a square
  root.  The ellipsoid is *measured* out of the axis's own transform by asking
  where two known points go, so a globe on WGS84, on a sphere, or in kilometres
  all come out where the axis would have put them.
* **On a map** the cells are built in longitude and latitude, cut, and then the
  whole vertex buffer is projected in one strided `proj_trans_generic` call
  rather than one `ccall` per vertex.  The seam is the meridian opposite the
  projection's `lon_0`, which PROJ publishes; failing that it is recovered by
  asking where the origin of the projected plane came from.

Supporting a new axis is a new `plot_target` method.  Supporting a new grid
system needs nothing: cell geometry is read through
`DiscreteGlobalGrids.cell_boundary`, which every system implements.

## Longitude near a pole

Three things can happen to a cell on a map, and the mesh builder tells them
apart by how far longitude turns in going once around the cell's ring:

* **`0`, and no long step** — an ordinary cell.  If it straddles the seam it is
  clipped there and the outside piece is moved a full turn, so it reappears on
  the far edge of the map.
* **`±360`** — the cell *contains* a pole.  Longitude cannot cut such a cell,
  but it describes it perfectly as the band between its ring and the pole, so it
  is drawn spanning the map and closed along the top or bottom edge.
* **`0`, but with a step of more than a quarter turn** — an edge running over or
  beside a pole, where the short way round is not obviously the right way.  The
  edge is walked, halving until each step is unambiguous; an edge crossing a pole
  exactly gets the two pole corners inserted.  IGEO7 puts each pole on an edge
  shared by two cells, so this is not a corner case there — it is what makes the
  top and bottom of a global IGEO7 map close.

## Backends

`primitive = automatic` gives CairoMakie one filled path per cell and every
other backend the mesh.  A vector renderer draws a path faster than a fan of
triangles and without the hairline seams that show between triangles of the same
cell; a GPU backend wants the opposite.  Both paths share the same tessellation,
so the cut and the pole handling do not depend on the backend.  A globe plot is
three-dimensional and gets the mesh whatever the backend is.

## Measurements

IGEO7 cells covering a 1° tile over the Alps, Julia 1.12, 8 threads, GLMakie
under `xvfb` (software rasterisation, so the *plot* column's second half is
pessimistic).  Reproduce with `bench/bench.jl`.

```
IGEO7 level 7   n=     169   convert    0.003 ->   0.030 s (  0.1x)   plot    0.385 ->   0.073 s (  5.3x)
IGEO7 level 8   n=    1035   convert    0.010 ->   0.003 s (  3.7x)   plot    0.040 ->   0.042 s (  0.9x)
IGEO7 level 9   n=    6945   convert    0.059 ->   0.016 s (  3.8x)   plot    0.090 ->   0.053 s (  1.7x)
IGEO7 level 10  n=   47659   convert    0.419 ->   0.073 s (  5.7x)   plot    0.441 ->   0.151 s (  2.9x)
IGEO7 level 11  n=  331491   convert    3.709 ->   0.541 s (  6.9x)   plot    3.629 ->   0.832 s (  4.4x)
IGEO7 level 12  n= 2313802   convert   26.146 ->   3.889 s (  6.7x)   plot   27.236 ->   6.588 s (  4.1x)
```

`convert` is the part this package replaces and is backend-independent.

The docs' hydrology tutorial at its level 12, run stage by stage with
`bench/hydrology_level.jl`:

| stage | `poly` | `dggpoly` |
| --- | --- | --- |
| build the elevation figure | 37.9 s | 19.4 s |
| save the elevation figure | 13.2 s | 10.3 s |
| build and save the terrain figures | 60.2 s | 15.8 s |
| peak RSS for the page | 15.4 GiB | 6.6 GiB |

Plotting stops being what sizes that page: regridding, at 5.2 GiB, becomes its
peak.  At level 13 (16.2 M cells) the page reaches 33.6 GiB and the recipe is no
longer what stops it — regridding, adjacency and the terrain analysis together
reach 18.8 GiB before the last two figures are built.

## Status

This is a proof of concept living in `lib/`; names and behaviour are not stable.
The pieces that earn their keep are meant to move into `DiscreteGlobalGrids`'
own Makie extension, at which point plain `poly` on a cell set can take this
path and existing code gets faster without changing.

Known gaps:

* A projection is applied to cell corners, not along cell edges, so an edge is
  drawn straight in the projected plane.  This is what `poly` does today as well,
  and it shows only for cells large enough for the projection to bend an edge
  visibly.
* Nothing here reduces the *number* of things drawn.  Past a few million cells
  the cost is the vertex buffer itself, and the answer is to resample cells onto
  the screen rather than to draw them all — a second recipe, not a faster mesh.
