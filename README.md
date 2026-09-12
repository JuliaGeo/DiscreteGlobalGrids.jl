# DiscreteGlobalGrids.jl

Work with data on global grids through a common interface.

A [discrete global grid system (DGGS)](https://en.wikipedia.org/wiki/Discrete_global_grid) divides the Earth's surface into cells
at different resolutions. Each cell has an identifier, a boundary, and
neighbors, much like a raster pixel.

DiscreteGlobalGrids brings six grid systems into the Julia geo ecosystem.
The same interface lets you:

- Locate cells, select regions, and calculate cell geometry.
- Compute with neighboring cells and follow parent–child relationships.
- Regrid data between rasters and global grids with
  [GlobalRegridding.jl](lib/GlobalRegridding/README.md).
- Work with dimensional arrays, process data in chunks, and read or write
  Zarr stores.

[DimensionalData.jl](https://github.com/rafaqz/DimensionalData.jl) connects array
values to grid cells through the `Cells` dimension. The companion
[DiscreteGlobalGridsVisualization](lib/DiscreteGlobalGridsVisualization/README.md)
package provides plotting with Makie.

## Quick start

Use Julia 1.11 or later. In the REPL, install the package from this repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaGeo/DiscreteGlobalGrids.jl")
```

Choose a system for its cell geometry and compatibility with your data:

| System | Cell shape |
| --- | --- |
| `IGeo7System()` | Hexagons, with twelve pentagons |
| `H3System()` | Hexagons, with twelve pentagons; uses H3 identifiers |
| `HEALPixSystem()` | Equal-area curvilinear diamonds |
| `A5System()` | Equal-area pentagons |
| `S2System()` | Geodesic quadrilaterals |
| `ISEA4RSystem()` | Equal-area rhombi |

A system describes a family of grids. `levelgrid` selects one resolution level;
higher levels have smaller cells. This example locates a cell near Zürich and
finds its neighbors:

```julia
import DiscreteGlobalGrids as DGG

sys = DGG.HEALPixSystem()
grid = DGG.levelgrid(sys, 4)
DGG.ncells(grid)                         # 3072

cell = DGG.cellat(grid, 8.5, 47.4)        # longitude, latitude in degrees
DGG.cell_area(grid, cell)                # area in steradians
DGG.neighbors(grid, cell)                # adjacent cell ids
DGG.ring(grid, cell, 2)                  # cells exactly two neighbor steps away

coarser = parent(sys, cell)              # the parent at level 3
finer = DGG.children(sys, cell)          # children at level 5
```

Replace `HEALPixSystem()` to use another system with the same operations.
Level numbers represent different cell sizes across systems. Use
`DGG.levelfor(sys, 100_000)` to choose a level with cells roughly 100 km across.
The [grid selection tutorial](docs/src/tutorials/choosing_a_grid.jl) explains
cell sizes, connectivity, and latitude conventions.

## Regridding data

Install `DimensionalData` with `Pkg.add("DimensionalData")` for this example.
Using the grid above, map a longitude/latitude raster onto HEALPix:

```julia
import DimensionalData as DD

lon = collect(-175.0:10.0:175.0)
lat = collect(-85.0:10.0:85.0)
raster = DD.DimArray([cosd(y) * cosd(x) for x in lon, y in lat],
                    (DD.X(lon), DD.Y(lat)))

values = DGG.regrid(raster; to = grid, method = DGG.Conservative())
size(values)  # (3072,)

# Reuse the weights for other fields on the same grids.
plan = DGG.plan_regrid(raster; to = grid, method = DGG.Conservative())
values = DGG.regrid(raster, plan)
```

`Conservative()` weights values by cell overlap area. `NearestCell()` samples
the source cell containing each destination centroid. `BarycentricPoint()`
interpolates between source samples. Point methods do not preserve integrals.

The result has a `Cells` dimension. Inputs with time or other non-spatial
dimensions retain those dimensions. Regional grids also work as destinations.
See [GlobalRegridding](lib/GlobalRegridding/README.md) for missing-data policies,
lazy execution, and weight caching.

## API

| API | Purpose |
| --- | --- |
| `levelgrid(sys, level)` | Select a complete grid at one resolution. |
| `cellat(grid, lon, lat)` | Find the cell containing a location in degrees. |
| `cell_polygon(grid, cell)`, `cell_area(grid, cell)` | Read a cell's geometry or area. |
| `neighbors(grid, cell)`, `ring(grid, cell, k)` | Find adjacent cells or cells at distance `k`. |
| `query(grid, Intersects(geometry))` | Select cells that intersect a region. |
| `regrid(data; to, method, ...)`, `plan_regrid(data; to, ...)` | Regrid data or prepare reusable weights. |

These names are available through `DGG` in the examples. See the
[grid interface](docs/src/api/grid-interface.md) for the full reference.

## How it works

Grids implement `AbstractGrid`, which describes a finite collection of cells.
The core interface provides cell counts, identifiers, boundaries, and centroids.
Generic algorithms use this information for spatial queries, geometry, and
neighborhood operations. Hierarchical systems also provide parent and child
relationships to accelerate searches and represent regions compactly.

A typed cell id identifies a cell and records its level. An integer index is
its position in a particular collection. `CellVector` represents cell collections,
including regions stored as compressed index ranges. `CellLookup` connects
these collections to dimensional arrays through the `Cells` dimension.

Geometry uses the unit sphere internally, and cell areas are in steradians.
Longitude/latitude entry points use degrees. For geodetic data, follow the
[coordinate guidance](docs/src/tutorials/choosing_a_grid.jl) when choosing your grid.

## Going further

The [tutorials](https://juliageo.org/DiscreteGlobalGrids.jl/dev/tutorials/) cover regional statistics, neighborhood
operations, regridding, and Zarr storage. Store I/O requires `using Zarr`.
To add a grid, follow [Writing a grid system](https://juliageo.org/DiscreteGlobalGrids.jl/dev/extending/) and use the
[conformance tests](lib/DiscreteGlobalGridsConformanceTesting/README.md).

## License and attribution

This package uses the [MIT license](LICENSE.md). The IGEO7 adjacency kernel
includes code from Alexander Kmoch's IGEO7.jl, with permission to relicense it
under MIT; the [source header](src/systems/IGeo7/gbt.jl) records the attribution.
H3 uses libh3 through `H3_jll`, and A5 includes arithmetic ported from upstream A5.

## AI disclosure

This package was created with the help of AI agents, including Claude and Codex,
and will continue to be developed with these agents.
