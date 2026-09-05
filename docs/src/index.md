# DiscreteGlobalGrids.jl

DiscreteGlobalGrids.jl lets you analyse data on global grids: move rasters
onto cells, select regions, compute with neighbours, and read or write Zarr
stores. The same operations work across IGEO7, H3, HEALPix, A5, S2 and ISEA4R.

A discrete global grid divides the Earth's surface into cells. Each cell has
an identifier, a boundary and neighbours, so you can work with it much as you
would with a raster pixel.

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

The globe shows IGEO7 at level 3. To work with a grid, choose a system and a
resolution level:

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

These calls also work with the other systems. [Choosing a
grid](tutorials/choosing_a_grid.md) compares cell shapes, sizes and coordinate
conventions.

## Grids, cells and data

| Term | Meaning |
|---|---|
| System | A family of grids at different resolutions, such as HEALPix |
| Grid | A collection of cells at one level, covering the globe or a region |
| Cell id | A typed identifier for a cell, including its level |
| Local index | A cell's position in a particular collection or array |
| `Cells` dimension | The link between an array's values and its grid cells |

Regridding a monthly raster produces an array with `Cells` and time dimensions.
Spatial selectors act on `Cells`; ordinary Julia indexing and reductions work
on the result. See the [grid interface](api/grid-interface.md) for the full
reference.

## Installation

Install the package and its plotting companion from the repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaGeo/DiscreteGlobalGrids.jl")
Pkg.add(url = "https://github.com/JuliaGeo/DiscreteGlobalGrids.jl",
        subdir = "lib/DiscreteGlobalGridsVisualization")
```

Julia 1.11 or newer is required. The tutorials import the packages needed
for each example. Plotting uses `DiscreteGlobalGridsVisualization` with a
Makie backend; store I/O requires `using Zarr` to load `dggread` and `dggwrite`.
Use `Pkg.develop(url = ...)` for a checkout you intend to edit.

## Choose a tutorial

Start with **Choosing a grid** and **Regridding** if you are new to DGGS data.
Each tutorial includes its own setup, so you can then follow the capability
you need:

| Task | Tutorial |
|---|---|
| Choose cell geometry, resolution and latitude convention | [Choosing a grid](tutorials/choosing_a_grid.md) |
| Bring a monthly raster onto a DGGS and export it back | [Regridding](tutorials/regridding.md) |
| Change grid system or resolution; compare interpolation methods | [Moving between DGGS](tutorials/between_grids.md) |
| Select cells by region and calculate regional means | [Zonal statistics](tutorials/zonal.md) |
| Represent a region compactly and regrid onto it | [Multi-order coverage](tutorials/multiorder.md) |
| Smooth a field, detect edges and traverse the cell graph | [Stencil operations](tutorials/stencils.md) |
| Use Geomorphometry and write a custom terrain kernel | [Hydrology](tutorials/hydrology.md) |
| Save data, reopen it lazily and read a region | [DGGS stores](tutorials/store_io.md) |
| Run a neighbourhood kernel over stored chunks | [Out of core](tutorials/out_of_core.md) |
| Work with HEALPix vectors, sky masks and cone searches | [The sky in HEALPix](tutorials/healpix_astronomy.md) |

The Earth-data examples use simple spherical setups to demonstrate the
operations. For work that requires alignment with geodetic data, follow the
[coordinate guidance](tutorials/choosing_a_grid.md#match-the-ellipsoid-of-the-source)
when choosing your grid.

## Reference and extension

The API pages cover [grids](api/grid-interface.md),
[spatial selection](api/selecting-cells.md),
[regridding methods](api/regridding-methods.md),
[region boundaries](api/boundaries.md),
[neighbours and stencils](api/neighbors.md),
[neighbour fields](api/neighbor-fields.md),
[store I/O](api/store-io.md),
[chunked computation](api/chunk-sweep.md) and
[subzone storage](api/subzone-layout.md).

To add a grid, follow [Writing a grid system](extending.md). The
[architecture guide](architecture.md) explains how grids, cell collections
and algorithms fit together.
