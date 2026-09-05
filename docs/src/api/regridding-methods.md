# Choosing a regridding method

```@meta
CurrentModule = DiscreteGlobalGrids
```

Choose a regridding method according to what each source value represents:
an average over a cell, a sample at a point, or a value to copy from the
containing cell. Pass the method with `regrid(...; method = ...)`.

| method | what a destination cell gets | conserves the integral |
|---|---|---|
| `Conservative()` | the mean of the source values over the ground the cell covers | yes |
| `BarycentricPoint()` | a sample at the cell's centre, interpolated between the source sample sites around it | no |
| `NearestCell()` | the value of the source cell the centre falls in | no |
| `DirectNearest()` | the same value, without building an operator for it | no |

`Conservative()` is the default. Its weights describe area overlaps;
`missingpolicy` controls how those weights apply when coverage is incomplete.

## Area means

Use `Conservative()` for values that represent cell averages, such as model
output or observations supplied as averages over pixels. It weights each source
value by its overlap with the destination cell. With complete coverage and
consistent overlap geometry, this preserves the area integral. Missing-data
thresholds and normalization affect that property at partially covered cells.

## Point samples

`BarycentricPoint()` treats each source value as sitting at one point, its
cell's sample site, and interpolates between the sites surrounding the
destination point: tensor Q1 on a raster's four surrounding samples, barycentric
coordinates on a triangle of them, mean-value coordinates on a larger convex
polygon. The weights are non-negative and sum to one, so the output stays inside
the range of the sources it came from, and a destination point landing exactly
on a source sample site reproduces that source exactly.

`NearestCell()` copies the value from the source cell containing the destination
centre. It preserves the selected value, with discontinuities at source-cell
boundaries. `DirectNearest()` produces the same result by direct lookup.

- Use `DirectNearest()` for a transfer that needs little plan reuse.
- Use `NearestCell()` when several fields share a plan or you need its operator.

These point methods do not preserve area integrals.

## Which one your data wants

Check the source dataset's definition of a pixel. Elevations published at
coordinates, often called **posts**, suit `BarycentricPoint()`. Elevations
published as averages over pixel footprints suit `Conservative()`.

The [hydrology tutorial](../tutorials/hydrology.md) uses the default area
method to demonstrate the workflow. Choose the point method when your analysis
needs to retain the interpretation of DEM values as samples. When aligning
Earth data, also check the [latitude convention](../tutorials/choosing_a_grid.md#match-the-ellipsoid-of-the-source).

A finer destination makes the difference especially visible: area averaging
repeats a source value across destination cells wholly inside it, while point
interpolation varies between sample sites. [Moving between
DGGS](../tutorials/between_grids.md) compares both methods at two resolutions.

## Missing data

`missingpolicy` decides what a destination cell holds when some of the source
under it is missing. `Weighted(t)` divides by the valid weight and blanks a cell
whose valid weight falls below `t` of the total.

For example, `Weighted(0.5)` marks a coastal cell missing when only 30% of
its total weight comes from valid source values. For point interpolation,
`Weighted(1)` requires the full interpolation weight to be valid. Use it when
you want missing samples to leave gaps; lower thresholds allow normalization
over the valid samples.

The policy and conservative-method reference:

```@docs
Weighted
Extensive
Conservative
```

## The sentinel a blanked cell holds

`missingpolicy` decides *which* cells are blanked. `missingval` decides what
they hold, and it is the source's own nodata convention unless you say
otherwise. A raster in is a raster out, declaring what it was handed:

| the source | the result | a blanked cell holds |
|---|---|---|
| a `Raster` with `missingval = -9999.0` | a `Raster` declaring `-9999.0`, over a concrete float array | `-9999.0` |
| a `Raster` with `missingval = missing` | a `Raster` declaring `missing`, over a `Union{Missing, Float64}` array | `missing` |
| a `Raster` declaring no `missingval` | a `Raster` declaring `NaN` | `NaN` |
| a `DimArray`, or a bare array | the same wrapper it came in as, with no nodata to declare | `missing` when the element type holds it, `NaN` otherwise |

The `missingval` keyword names the nodata of the whole regrid: invalid on the
way in, written on the way out. Give it a value the destination element type
holds and the result stays concrete, which is what makes

```julia
DGG.regrid(ras; to = grid, missingval = NaN)
```

the fast reading of a `Union{Missing, Float64}` raster — a `Float64` array,
smaller in memory and quicker in every operation that follows.

`regrid!` uses the missing-value convention declared by `dest`. Its element
type must be able to hold that value.

## Rims, degeneracies, and poles

Point interpolation requires a polygon of source sample sites containing the
destination point. The following cases remain unmapped:

  - the strip between the outermost sample sites and the source's own boundary,
    where no such polygon exists;
  - a destination whose stencil names a cell the source collection does not hold
    — a partial grid's rim, where nothing is substituted from elsewhere;
  - a degenerate polygon: repeated sites, no area, or a fold.

Copernicus DEM adds one more. Its posts stop a fraction of a pixel short of the
poles, and the natural polygon covering that gap has the entire polemost post
row as its corners — tens of thousands of them. A point there takes the policy
the method's `poles` keyword names. `BarycentricPoint(poles = NearestCell())`,
the default, gives it the nearest post of that row: one entry, weight one, and a
discontinuity where the fallback begins. `BarycentricPoint(poles = nothing)`
leaves those points unmapped instead. A source whose own sample sites reach the
poles — every conforming grid — has no such region and ignores the keyword.

## An example

A coarse source and a finer destination, which is where the two readings differ
most.

```@example methods
import DiscreteGlobalGrids as DGG

src = DGG.levelgrid(DGG.IGeo7System(), 3)
dst = DGG.levelgrid(DGG.HEALPixSystem(), 5)

peak(p) = 2000 * p[3] + 800 * p[1] * p[2]
elevation = [peak(DGG.cell_centroid(src, c)) for c in DGG.CellVector(src)]

area = DGG.regrid(elevation; to = dst, from = src)
point = DGG.regrid(elevation; to = dst, from = src,
    method = DGG.BarycentricPoint(), missingpolicy = DGG.Weighted(1))
count(isfinite, area), count(isfinite, point)
```

The counts show how many destination cells each method can fill. A point
remains unmapped when no valid source polygon contains it.

The source in this example samples a known function. Evaluate that function
at the destination centres to compare the interpolation errors:

```@example methods
both = isfinite.(area) .& isfinite.(point)
truth = [peak(DGG.cell_centroid(dst, c)) for c in DGG.CellVector(dst)]

(maximum(abs, area[both] .- point[both]),
 maximum(abs, area[both] .- truth[both]),
 maximum(abs, point[both] .- truth[both]))
```

This comparison tests agreement with point samples of `peak`. Testing area
averages would require averaging `peak` over each destination footprint.
Use the test that matches the meaning of your data.

## Implementing point interpolation for a source space

To add point interpolation to a `GlobalRegridding.RegridSpace`, describe the
polygon of source sample sites surrounding a destination point. Implement
these three hooks:

  - `GlobalRegridding.hasdualcells(space)` returns `true`. A space declaring
    neither this nor `hascellchart` is refused rather than silently answering
    every point with nothing.
  - `GlobalRegridding.samplerstate(space)` returns whatever the lookup needs
    prepared once — axis locators, topology tables. It is read concurrently
    during a sweep and must not be written.
  - `GlobalRegridding.dualcellat(sampler, p)` returns a
    `GlobalRegridding.DualCell(indices, nodes, kind)`: the source *local*
    indices of the sites in cyclic order, their coordinates in one plane, and a
    `GlobalRegridding.BasisKind` — `Bilinear` for isoparametric Q1 on a
    quadrilateral, `MeanValue` for a convex polygon, which on a triangle is that
    triangle's barycentric coordinates. Return a cell with no nodes where none
    contains `p`.

`sampler` is `GlobalRegridding.Sampler`, which those two dispatch on; it carries
the space, its sample sites and that state. A space with a cell chart writes its
nodes in that chart and needs nothing more; a space without one also defines
`GlobalRegridding.chartat(sampler, p)` to say which plane the nodes are in.
Chunking is never an input: a stencil is built whole and partitioned afterwards.
