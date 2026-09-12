# GlobalRegridding.jl

Remap data from one grid to another, using a variety of methods.

GlobalRegridding currently focuses on spherical regridding, with
plans to generalize to planar regridding. It supports conservative area
weighting, point sampling, and interpolation. It computes results eagerly or
reads them on demand in chunks.

This is an extension of the concept from [ConservativeRegridding.jl](https://github.com/JuliaGeo/ConservativeRegridding.jl),
but meant to also support:
- derivative-aware and interpolating regridding, which requires neighbourhood information,
- chunk-aware regridding for big data,
- lazy and cached regridding for reuse across multidimensional datasets.


Grids are described by the `RegridSpace` interface.  In here, we define the `RasterGrid`, which assumes longitude and latitude in degrees.
Other packages like [DiscreteGlobalGrids.jl](../../README.md) define other spaces for their grids.

## Quick start

Use Julia 1.11 or later.  In the REPL, type:
```julia
using Pkg
Pkg.add("GlobalRegridding")
```
to add the package.  Then, load via `using GlobalRegridding`.



Choose a method with the `method` keyword.  Currently, three methods are implemented,
though we'll add more in the future:

- `Conservative()` constructs weights by the overlap of the source and destination cells.  This is the same as ConservativeRegridding's default method.
- `NearestCell()` samples the source cell containing the destination cell's centroid.  This is a nearest neighbor interpolation.
- `BarycentricPoint()` interpolates the value at the destination cell's centroid by a weighted average of the source cell values, following [Barycentric interpolation](https://en.wikipedia.org/wiki/Barycentric_coordinate_system).

Point methods do not preserve integrals. This example uses `Conservative()` to
coarsen a global raster from 10° to 20° cells:

```julia
using GlobalRegridding
using DimensionalData

lon = collect(-175.0:10.0:175.0)
lat = collect(-85.0:10.0:85.0)
data = DimArray([cosd(y) * cosd(x) for x in lon, y in lat],
                (X(lon), Y(lat)))
target = RasterGrid(DimArray(zeros(18, 9),
    (X(collect(-170.0:20.0:170.0)), Y(collect(-80.0:20.0:80.0)))))

result = regrid(data; to = target, method = Conservative())
size(result)  # (18, 9)

# Keep the weights for repeated use on the same grids.
plan = plan_regrid(data; to = target, method = Conservative())
result = regrid(data, plan)
```

`RasterGrid` assumes longitude and latitude in degrees; projected coordinates
require an explicit transform. For plain arrays, supply a source space with
`from`. Put spatial dimensions first, in the source space's cell order.
Dimensional results use the destination's axes and retain non-spatial dimensions,
such as time.

## API

| API | Purpose |
| --- | --- |
| `regrid(data; to, from, method, ...)` | Build a plan and return the regridded data. Only `to` is required. |
| `plan_regrid(data; to, ...)` | Build a reusable plan without reading source values. |
| `regrid(data, plan)` | Apply an existing plan to compatible data. |
| `regrid!(dest, data, plan)` | Write results into a preallocated array. |
| `RasterGrid(data_or_dims; ...)` | Describe raster cells from a dimensional array or a tuple of dimensions. |

See the [API docstrings](src/api.jl) for all keywords.

## How it works

1. Source and destination spaces describe cell geometry, cell order, and spatial
   indexes on a common sphere.
2. The method builds weights from geometry. `Conservative()` uses spherical
   polygon intersection areas; point methods use source samples.
3. A plan applies these weights to each field or time slice. Eager plans store
   one sparse weight block. Lazy plans build and cache blocks as needed.

The default `Weighted(0.5)` policy returns means normalized by valid source
weight. It marks results as missing when valid weight falls below 50% of the
source-covered weight, or no valid weight remains. `missing`, `NaN`, and declared
nodata values do not contribute. `Extensive()` returns unnormalized weighted
sums; with `Conservative()`, these represent integrals over valid source coverage.

Chunked sources default to lazy execution; set `lazy = true` to request it
explicitly. A chunk dependency graph identifies which source chunks each
destination chunk needs. For lazy execution, `chunks` sets destination tiling
and `budget` sets a memory target in bytes. Use `storage = Spilled(dir)` to store
weights on disk.
