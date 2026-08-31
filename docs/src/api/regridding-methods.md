# Choosing a regridding method

```@meta
CurrentModule = DiscreteGlobalGrids
```

`regrid` takes a `method`, and the choice is a statement about what the source
*values mean*, not about how accurate you would like the answer to be. There are
three readings of the data, and four names for them.

| method | what a destination cell gets | conserves the integral |
|---|---|---|
| `Conservative()` | the mean of the source values over the ground the cell covers | yes |
| `BarycentricPoint()` | a sample at the cell's centre, interpolated between the source sample sites around it | no |
| `NearestCell()` | the value of the source cell the centre falls in | no |
| `DirectNearest()` | the same value, without building an operator for it | no |

`Conservative()` is the default, and the only one whose output is an area
quantity.

## Area means

`Conservative()` clips every source cell against every destination cell it
meets and weights by the overlap area. If the source values are averages over
their own cells — most model output, most gridded observations, a DEM
distributed as area means — this is the reading that keeps the field's integral,
and no other method does.

## Point samples

`BarycentricPoint()` treats each source value as sitting at one point, its
cell's sample site, and interpolates between the sites surrounding the
destination point: tensor Q1 on a raster's four surrounding samples, barycentric
coordinates on a triangle of them, mean-value coordinates on a larger convex
polygon. The weights are non-negative and sum to one, so the output stays inside
the range of the sources it came from, and a destination point landing exactly
on a source sample site reproduces that source exactly.

`NearestCell()` is the degenerate case of the same idea — one source, weight
one. It costs less and it is not continuous. `DirectNearest()` answers exactly
what `NearestCell()` answers, element for element, by looking the source cell up
and copying the value instead of assembling a matrix of ones and multiplying by
it: prefer it where a plan is applied about as often as it is built, and
`NearestCell()` where one plan serves many different sources or the operator
itself is wanted.

None of them is conservative. Sampling a field at points and calling the result an
area mean is the mistake this page exists to prevent.

## Which one your data wants

Copernicus DEM is the clear case for points: its pixels are **posts**, elevations
published *at* a coordinate rather than averaged over a footprint, so
`BarycentricPoint()` is the faithful reading of them and `Conservative()`
quietly turns a post into a cell average it never was. A DEM distributed as area
means, and anything whose documentation says "cell average", wants
`Conservative()`.

When source and destination are close in size the two agree closely. They part
company when the destination is *finer* than the source: an area mean has
nothing finer to say, so it repeats one source value across the cells inside it,
while a point sample interpolates between the sites around each one.

## Missing data

`missingpolicy` decides what a destination cell holds when some of the source
under it is missing. `Weighted(t)` divides by the valid weight and blanks a cell
whose valid weight falls below `t` of the total.

An area weight row is a spectrum: a coastal cell may be 30 % covered and
`Weighted(0.5)` blanks it. A point row is different — it is complete or it is
empty, four posts or none — so the threshold is a yes/no switch on whether an
absent source may be interpolated across. `Weighted(1)` is the strict choice and
the one to prefer for point output: a stencil naming a post with no value blanks
the cell instead of renormalising over the posts that do.

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

`regrid!` reads the buffer you hand it rather than the source, so blanked cells
take the sentinel `dest` declares. A buffer that cannot hold the sentinel says
so instead of guessing.

## Rims, degeneracies, and poles

Point interpolation is defined inside the polygons whose corners are source
sample sites, and nowhere else. Three things therefore stay unmapped rather than
being invented:

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

The point field leaves a couple of cells unmapped: no polygon of source sample
sites contains them, and the alternative to a blank would be extrapolation.

Where both placed a value, they differ by tens of metres on a field whose range
is thousands, and the smooth field they were both built from says which is
closer to it:

```@example methods
both = isfinite.(area) .& isfinite.(point)
truth = [peak(DGG.cell_centroid(dst, c)) for c in DGG.CellVector(dst)]

(maximum(abs, area[both] .- point[both]),
 maximum(abs, area[both] .- truth[both]),
 maximum(abs, point[both] .- truth[both]))
```

The area column is not wrong: over a destination cell coarser than the source it
is the better answer, and it is the only one that keeps the total. On a
destination finer than the source it has nothing finer to say.

## Making a source space point-samplable

A `GlobalRegridding.RegridSpace` answers point queries by describing, for one
destination point, the polygon of its own sample sites containing it. Three
hooks and a small closed set of names:

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
