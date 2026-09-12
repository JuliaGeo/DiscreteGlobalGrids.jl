# Rasterization, extraction, and zonal statistics

These operations are owned by DiscreteGlobalGrids. Call them qualified as
`DGG.rasterize`, `DGG.extract`, and `DGG.zonal` when also using Rasters.
They work in memory on a DGGS cell axis and accept GeoInterface geometries,
features, collections, nested inputs, and Tables tables. Zonal operations also
accept longitude/latitude `Extents.Extent` regions and spherical caps. Table inputs use
`:geometry` by default; `geometrycolumn=:geom` selects another column and
`geometrycolumn=(:longitude, :latitude)` reads point coordinates.

Input coordinates follow the package geometry contract: longitude/latitude
for ordinary GeoInterface points, or `UnitSphericalPoint` coordinates. Edges
are great-circle arcs. These functions do not transform projected coordinates.
Cell intersection uses each system's published `cell_boundary`, including
its approximation of curved edges.

## Boundary rules

| `boundary` | Polygon membership |
|:--|:--|
| `:center` (default) | The polygon covers the cell's canonical interior representative, `cell_centroid`. |
| `:intersects` | Any intersection, including shared edges or vertices. |
| `:touches` | Alias for `:intersects`; different from the DE9IM `Touches` predicate. |
| `:inside` | The entire cell lies within the polygon. Coincident polygon boundaries are permitted. |

Lines select intersected cells and points use deterministic containing-cell
lookup. `shape=:point`, `:line`, or `:polygon` overrides geometry interpretation;
plural aliases are also accepted. The canonical cell representative is not
necessarily the mathematical area centroid.

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

`to` accepts a grid, partial grid, cell vector or lookup, a cell dimension,
dimension tuple, or existing dimensional array/stack. A system requires
`to=DGG.H3System(), level=5`. A `MultiOrderCellSet` keeps its stored cells;
use `CellVector(set)` explicitly to expand it. Core allocation returns
`DimArray`/`DimStack`; Rasters templates retain their wrappers through the
optional Rasters extension. The cell dimension can occur anywhere in a cube;
spatial values are broadcast over the other dimensions.

`fill` can be scalar, iterable, a feature-property/column symbol, a tuple of
symbols for several layers, a NamedTuple of layer fills, or a function updating
the current cell. Except for `count`, `fill` is required. Multiple features
need a reducer, binary `op`, or function fill. Input order is preserved.
`init`, `eltype`, and `missingval` can be supplied for custom cell types.

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

Every cell owns its mutable state. Overlapping cells hold region IDs in input
order. Mutating binary operations such as `append!` are also supported.

`DGG.boolmask`, `DGG.missingmask`, and their mutating variants use the same
membership rules. `DGG.mask(A; with=regions)` and `mask!` replace unselected
values; `invert=true` reverses the selection.

## Extract labelled slices

```julia
rows = DGG.extract(sums, points; id=true, index=true)
```

Rows are NamedTuples containing selected layer values and, by default,
`:geometry`. Points retain their input coordinates; polygon and line rows
report cell representative longitude/latitude. `index=true` reports the local
cell-axis position, not a global cell ID. `name` selects stack layers.
`skipmissing=true` drops missing rows; `flatten=false` groups polygon/line rows
by feature. Unsampled time/band dimensions remain labelled slices.

## Reduce zones without repeating selection per slice

```julia
stats = DGG.zonal(mean, sums; of=regions, emptyval=NaN)
```

The default `spatialslices=true` reduces the cell dimension independently for
every nonspatial slice. `Ti × Cells × Band` becomes `Ti × Band × Zone`.
This default extends Rasters' whole-zone behavior. `spatialslices=false`
reduces the whole selected cube, and a tuple chooses dimensions to reduce
(including the cell dimension). A single geometry omits the Zone dimension.
Stack layers retain their own nonspatial dimensions.

Zones outside a regional holding return `missing`, retaining nonspatial shape.
For zones within the holding, missing filtering and `emptyval` apply per slice. Without `emptyval`, the
reducer receives an empty iterator and retains its ordinary empty-input
behavior, including errors. Means are unweighted cell means; area weighting
must be requested through an explicit reducer.

## Compatibility and execution

The initial compatibility reference is Rasters v0.15. Geometry selection uses
a spherical grid/edge dual-tree traversal with prepared point location,
conservative subtree acceptance, and exact leaf predicates. Common reducers
stream accumulators; arbitrary iterable reducers gather ordered values.
Threading uses disjoint output ownership and serial fallback for unsafe custom
operations. `threadsafe=true` opts a custom in-place operation into threading.

`progress` and `verbose` are accepted compatibility controls; this version does
not display progress bars. File output (`filename`, `suffix`, `force`), raster
`res`/`size`, `crs`/`mappedcrs`, reprojection, fractional coverage, and chunk execution are outside
the in-memory API and unsupported options raise errors. Future chunk execution
can reuse the partitioning API after selection plans are stabilized.
