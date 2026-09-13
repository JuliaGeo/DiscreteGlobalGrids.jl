# DiscreteGlobalGridsVisualization.jl

Plot global grid cells and their data with Makie.

This companion to [DiscreteGlobalGrids.jl](../../README.md) draws complete grids,
regional cell collections, and data with a `Cells` dimension. It supports plain
longitude/latitude maps, projected maps, and globes through GeoMakie.

The package provides three plot types:

- `dggpoly` draws filled cells, with optional outlines. Use it to show individual
  cell values or categories.
- `dggsurface` interpolates between values at cell centers. Use it for continuous
  fields such as temperature, or supply heights to draw terrain.
- `dggresample` adjusts the displayed resolution as you zoom. Use it when the
  data contains more cells than the screen can show.

The plotting API is experimental and may change.

## Quick start

Use Julia 1.11 or later. In the REPL, install the packages and a Makie backend:

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaGeo/DiscreteGlobalGrids.jl")
Pkg.add(url = "https://github.com/JuliaGeo/DiscreteGlobalGrids.jl",
        subdir = "lib/DiscreteGlobalGridsVisualization")
Pkg.add("CairoMakie")
```

This example colors HEALPix cells by the sine of their latitude and saves a map:

```julia
import DiscreteGlobalGrids as DGG
using DiscreteGlobalGridsVisualization
using CairoMakie

CairoMakie.activate!()
grid = DGG.levelgrid(DGG.HEALPixSystem(), 3)

# The third coordinate of a unit-sphere centroid is sin(latitude).
values = [DGG.cell_centroid(grid, c)[3] for c in DGG.CellVector(grid)]

figure, axis, plot = dggpoly(grid;
    color = values, colormap = :viridis, strokewidth = 0.3,
    axis = (; xlabel = "Longitude (°)", ylabel = "Latitude (°)",
             aspect = DataAspect()))
Colorbar(figure[1, 2], plot; label = "sin(latitude)")
save("global-grid.png", figure)
figure
```

`color` accepts one color or a vector with one value per cell, in the collection's
order. `dggpoly` omits cells whose color is `missing` or `NaN`.
For a dimensional array `A`, use `dggpoly(A; color = A)` to color cells by its values.

Use GLMakie or WGLMakie for interactive plots. With GeoMakie loaded, the `!`
forms accept a `GeoAxis` for a projected map or a `GlobeAxis` for a globe.
The plot reads the projection from the axis.

## API

| API | Purpose |
| --- | --- |
| `dggpoly(cells; color, ...)` | Draw filled cells and optional outlines. |
| `dggsurface(cells; color, ...)` | Draw a continuous field between cell centers. |
| `dggsurface(cells, heights; color, ...)` | Draw a surface with a height per cell. |
| `dggresample(cells; color, ...)` | Draw cells at a resolution suited to the current view. |
| `dggpoly!(axis, cells; ...)` | Add cells to an existing axis; the other plot types also have `!` forms. |

The plot types accept grids, `CellVector`s, `CellLookup`s, and multi-order cell
sets. They also accept a system followed by a vector of cell ids.
For `dggsurface(A)`, a one-dimensional array's values supply heights;
`color = A` also uses them for color.

## How it works

The recipes read geometry through the DiscreteGlobalGrids interface, so the
same plotting code works across grid systems. On maps, they split geometry
at the projection's cut meridian and handle cells around the poles. On globes,
they place vertices directly in the axis's three-dimensional space.

`dggpoly` builds a combined triangle mesh for GPU backends. For flat maps,
CairoMakie uses filled polygon paths. The recipes keep color data separate from
geometry.
`dggsurface` connects neighboring cell centers into triangles, which interpolate
the values between centers. At a region's edge, the surface ends at the
outermost centers.

`dggresample` follows the grid hierarchy to select visible cells at a suitable
resolution. By default, each displayed cell samples one underlying value.
Pass `aggregate = mean`, with `using Statistics`, to average the underlying
values instead. Aggregation controls the displayed values.

See the [plot recipe](src/recipe.jl), [surface recipe](src/surface_recipe.jl),
and [resampling recipe](src/resample_recipe.jl) docstrings for further options.

## License

This package uses the repository's [MIT license](LICENSE.md).

## AI disclosure

This package was created with the help of AI agents, including Claude and Codex,
and will continue to be developed with these agents.
