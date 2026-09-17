```@meta
CurrentModule = DiscreteGlobalGrids
```

# Rasterization, extraction, and zonal statistics

These ten verbs are owned by DiscreteGlobalGrids and share the names Rasters.jl
uses, so call them qualified as `DGG.rasterize`, `DGG.extract`, `DGG.zonal` and
so on when both packages are loaded. They work in memory on a DGGS cell axis.

Every verb accepts the same geometry inputs: GeoInterface geometries, features,
feature collections, nested iterables of those, Tables tables, longitude/latitude
`Extents.Extent` regions, and `GO.UnitSpherical.SphericalCap` caps. Table inputs
read `:geometry` by default; `geometrycolumn=:geom` selects another column and
`geometrycolumn=(:longitude, :latitude)` reads point coordinates.

Input coordinates follow the package geometry contract: longitude/latitude for
ordinary GeoInterface points, or `UnitSphericalPoint` coordinates. Edges are
great-circle arcs, and these functions do not transform projected coordinates.
Cell intersection uses each system's published [`cell_boundary`](@ref),
including its approximation of curved edges.

## Boundary rules

| `boundary` | Polygon membership |
|:--|:--|
| `:center` (default) | The [`CentroidCovered`](@ref) rule: the polygon covers the cell's canonical interior representative, [`cell_centroid`](@ref). |
| `:intersects` | Any intersection, including shared edges or vertices. |
| `:touches` | Alias for `:intersects`; different from the DE9IM [`Touches`](@ref) predicate. |
| `:inside` | The entire cell lies within the polygon. Coincident polygon boundaries are permitted. |

Lines select intersected cells and points use deterministic containing-cell
lookup. `shape=:point`, `:line`, or `:polygon` reinterprets a geometry as that
kind — the rings of a polygon become lines, its vertices become points — and
plural aliases are accepted. The canonical cell representative is not necessarily
the mathematical area centroid.

## Destinations

`to` accepts a grid, partial grid, cell vector or lookup, a cell dimension,
dimension tuple, or an existing dimensional array or stack. A system requires a
level: `to=DGG.H3System(), level=5`. A `MultiOrderCellSet` mixes levels and is
not a destination; pass `CellVector(set)` or a level grid instead. Core
allocation returns `DimArray`/`DimStack`, and Rasters templates keep their
wrappers through the optional Rasters extension. The cell dimension can occur
anywhere in a cube; spatial values broadcast over the other dimensions.

## Rasterize values

```julia
import DiscreteGlobalGrids as DGG
import GeoInterface as GI
using Statistics

grid = DGG.levelgrid(DGG.HEALPixSystem(), 5)
points = [GI.Point((10.0, 25.0)), GI.Point((10.0, 25.0))]
counts = DGG.rasterize(count, points; to=grid)
sums = DGG.rasterize(sum, points; to=grid, fill=[2, 4])
DGG.rasterize!(sums, points; op=+, fill=1)
```

`mean` is `sum ./ count`, as in Rasters: `init` joins the sum, not the count.

Here is the vector-valued pattern from Rasters' crazy rasterization tutorial:

```julia
regions = [GI.Polygon([[(-5.0, 5.0), (15.0, 5.0), (15.0, 30.0),
                       (-5.0, 30.0), (-5.0, 5.0)]]),
           GI.Polygon([[(5.0, 15.0), (25.0, 15.0), (25.0, 40.0),
                       (5.0, 40.0), (5.0, 15.0)]])]
ids = DGG.rasterize(regions; to=grid, op=vcat,
    fill=[[i] for i in eachindex(regions)], eltype=Vector{Int},
    init=Int[], missingval=Int[], boundary=:intersects)
```

Every cell owns its mutable state, and overlapping cells hold region IDs in input
order. Mutating binary operations such as `append!` also work.

## Extract labelled slices

```julia
rows = DGG.extract(sums, points; id=true, index=true)
```

`index=true` reports the local cell-axis position, not a global cell ID. Points
retain their input coordinates; polygon and line rows report cell representative
longitude/latitude. Unsampled time and band dimensions stay labelled slices, and
`skipmissing=true` drops a row when any element of its slice is missing.

## Reduce zones without repeating selection per slice

```julia
stats = DGG.zonal(mean, sums; of=regions, emptyval=NaN)
```

The default `spatialslices=true` reduces the cell dimension independently for
every nonspatial slice, so `Ti × Cells × Band` becomes `Ti × Band × Zone`. This
extends Rasters' whole-zone behavior. `spatialslices=false` reduces the whole
selected cube and answers an array with a plain `Vector`, one entry per zone, and
a stack with a `NamedTuple` of those vectors. A single geometry omits the Zone
dimension, and stack layers retain their own nonspatial dimensions.

## Reference

```@docs
rasterize
rasterize!
extract
zonal
mask
mask!
boolmask
boolmask!
missingmask
missingmask!
```

## Compatibility and execution

The initial compatibility reference is Rasters v0.15. Geometry selection uses a
spherical grid/edge dual-tree traversal with prepared point location,
conservative subtree acceptance, and exact leaf predicates. Common reducers
stream accumulators; arbitrary iterable reducers gather ordered values. Threading
uses disjoint output ownership: built-in reducers and operations run threaded,
and a custom `op`, reducer, or fill function runs serially unless
`threadsafe=true`.

`progress` and `verbose` are accepted compatibility controls; this version
displays no progress bars. File output (`filename`, `suffix`, `force`), raster
`res`/`size`, `crs`/`mappedcrs`, reprojection, fractional coverage, and chunk
execution are outside the in-memory API, and unsupported options raise errors.
Future chunk execution can reuse the partitioning API once selection plans are
stable.
