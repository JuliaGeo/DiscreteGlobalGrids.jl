# DiscreteGlobalGrids.jl

Six discrete global grid systems — IGEO7, H3, HEALPix, A5, S2, ISEA4R — behind
one small interface, with every algorithm written against the interface exactly
once.

```@example index
import DiscreteGlobalGrids as DGG
using GeoMakie
using WGLMakie
using DiscreteGlobalGridsVisualization

WGLMakie.activate!()

figure = Figure(size = (640, 640), figure_padding = 2)
axis = GlobeAxis(
    figure[1, 1];
    source = "+proj=longlat +R=1",
    dest = GeoMakie.Geodesy.Ellipsoid(; a = "1", b = "1"),
    camera_longlat = (10, 25), camera_altitude = 1.9,
)
meshimage!(axis, -180 .. 180, -90 .. 90, fill("#f0faea", 1, 1);
    zlevel = -0.05, npoints = 300)
dggpoly!(axis, DGG.levelgrid(DGG.IGeo7System(), 3);
    color = "#dcf5d7", strokecolor = "#2c7a1e", strokewidth = 0.8)
lines!(axis, GeoMakie.coastlines(); color = ("#212529", 0.7), linewidth = 1.0,
    zlevel = 0.002)
figure
```

That is IGEO7 at level 3: hexagons, with twelve pentagons where an
icosahedron's vertices fall. A grid is a system and a level and nothing else:

```@example index
grid = DGG.levelgrid(DGG.IGeo7System(), 4)
```

```@example index
DGG.cellat(grid, 8.5, 47.4)     # the cell under Zürich, as a typed id
```

```@example index
DGG.cellsize(grid)              # a typical cell's width, in metres
```

```@example index
import Extents
DGG.query(grid, DGG.Intersects(Extents.Extent(X = (5, 12), Y = (45, 50))))
```

Swap `IGeo7System()` for any of the six and nothing else changes.
[Choosing a grid](tutorials/choosing_a_grid.md) is how to decide which one you
want; the [README](https://github.com/JuliaGeo/DiscreteGlobalGrids.jl) walks
the whole surface.

## The mental model

Two tiers. A **grid** is one finite collection of cells on the sphere — a
complete level, or a regional subset of one — and geometry, stencils and
queries are all answered there. A **system** adds the parent/child hierarchy
across levels, always as a fast path: hierarchy is an optimisation, never a
semantic. A bare `Int` is always an **index** in `1:ncells(grid)`, a local
index into that collection's own storage; a typed cell id knows its own level.

## Installation

Neither package is in the General registry yet, so both install from the
repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaGeo/DiscreteGlobalGrids.jl")
Pkg.add(url = "https://github.com/JuliaGeo/DiscreteGlobalGrids.jl",
        subdir = "lib/DiscreteGlobalGridsVisualization")
```

The second one is the plotting companion, and it is what the tutorials draw
with: `dggpoly`, `dggpoly!` and `dggsurface` live there. `using Makie` on its
own teaches Makie how to read a cell set, but does not give you those verbs.

`Pkg.develop(url = ...)` instead for a checkout you intend to edit. Julia 1.11
or newer is required. Two more capabilities arrive with the package that
provides them: `using Makie` (or a backend) makes cells drawable, and
`using Zarr` gives `dggread` and `dggwrite` real methods — without it they
exist only to report that Zarr is missing.

## Where to go next

The [DGGS gallery](all_dggs.md) draws every system. The tutorials run in order:

  - [Choosing a grid](tutorials/choosing_a_grid.md) — hexagons or
    quadrilaterals, equal-area or not, how fine, and on which sphere.
  - [Regridding: getting data onto a grid](tutorials/regridding.md) — a real
    monthly climatology moved conservatively onto a DGGS, and back to a raster.
  - [Stencil operations](tutorials/stencils.md) — smoothing, Laplacians and
    diffusion from each cell's neighbourhood.
  - [Zonal statistics](tutorials/zonal.md) — reduce a field over regions, with
    spatial queries.
  - [Multi-order coverage](tutorials/multiorder.md) — one region at every
    resolution at once.
  - [Moving between DGGS](tutorials/between_grids.md) — across systems and
    across resolutions, and whether your values are areas or points.
  - [Hydrology: a DEM on an IGEO7 grid](tutorials/hydrology.md) — elevation
    data on a regional subset, and flow routing across it.
  - [A round trip through a DGGS store](tutorials/store_io.md) — `dggwrite` and
    `dggread` over a Zarr store.
  - [Out of core](tutorials/out_of_core.md) — a stencil over a stored cube that
    never holds more than one chunk and its halo.
  - [The sky in HEALPix](tutorials/healpix_astronomy.md) — nested order, cone
    searches, and a galactic-plane cut.

Six API pages sit under the tutorials, one for each verb that needs more room
than a worked example gives it:

  - [Choosing a regridding method](api/regridding-methods.md)
  - [Region boundaries](api/boundaries.md) — `halo`, `border`, `interior`,
    `adjacency`, and the engines behind them
  - [Reading and writing DGGS stores](api/store-io.md)
  - [Sweeping a cube along its chunk lines](api/chunk-sweep.md)
  - [Requesting neighbour fields](api/neighbor-fields.md)
  - [The ancestor-subzone layout](api/subzone-layout.md)

[Writing a grid system](extending.md) is the other direction: what the
interface asks of a new system, built up on a plain lon/lat grid. There is a
second worked example in the repository — `DGG.CopernicusDEMSystem` implements
the same interface for the Copernicus DEM raster lattice, which is a raster
rather than a DGGS and so is not one of the six, and
`examples/copernicus_dem/copernicus_dem.jl` puts it to work.
