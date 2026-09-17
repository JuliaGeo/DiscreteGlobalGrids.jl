# Selecting cells

```@meta
CurrentModule = DiscreteGlobalGrids
```

Select cells by their spatial relationship to a point, polygon or spherical
cap. [`query`](@ref) returns cells matching a predicate, such as intersection
or containment. [`covering`](@ref) finds cells covering a region, and
[`covering_indices`](@ref) returns positions in a cell collection.

The relation is a DE9IM predicate — [`Intersects`](@ref), [`Within`](@ref),
[`Disjoint`](@ref) and the rest — read with spherical semantics rather than
planar ones. The result depends on the query form:

| Query | Result |
| --- | --- |
| `query(grid, predicate)` or `query(sys, predicate; level=l)` | Sorted `Vector` of typed cell IDs at one level |
| `query(sys, MultiOrderCoverage(target); level=l)` or `maxcells=n` | `MultiOrderCellSet` with cells at several levels |

A query result contains identities, not field values. Use the collection
conversions below to attach values or apply region operations.

## Cell collections

A `CellVector` holds cells at one level, and a `CellLookup` connects them to a
`Cells` dimension. `region` obtains the cell collection used by regional
operations, including for a stored axis.

| You have | You need | Conversion |
| --- | --- | --- |
| A complete or partial grid `g` | A cell-ID vector | `CellVector(g)` |
| Sorted unique IDs `ids` at level `l` | A compressed one-level collection | `CellVector(sys, l, ids)` |
| A `CellVector` `cv` | Grid geometry and queries | `PartialGrid(cv)` |
| A grid or `CellVector` | A dimensional array axis | `CellLookup(g)` or `CellLookup(cv)`, wrapped in `Cells` |
| A mixed-level set `set` | A one-level collection | `CellVector(set; level=l)` |
| A supported region or stored lookup | Its regional collection | `region(x)` |
| A `CellVector` | Its recorded source set | `cellset(cv)` |

`CellVector` preserves canonical cell order. Keep values in that same order;
these conversions do not sort an unrelated data array for you. Expansion from
a mixed-level set uses hierarchy membership, which need not equal the union of
its drawn polygons.

`PartialGrid(sys, l, ids)` requires canonical, strictly ascending, unique IDs
at the stated level. It retains `ids` by reference: later mutation can invalidate
the region. It does not provide an owning constructor that normalizes ordinary
unsorted input. `CellVector(sys, l, ids)` also requires strictly ascending IDs.

```@docs
AbstractCellVector
CellVector
AbstractCellLookup
CellLookup
region
```

## Query functions

```@docs
query
covering
covering_indices
cellset
```

## The predicates

The DE9IM types are DE9IM.jl's and `CentroidCovered` is this package's; the
semantics throughout are this package's, evaluated on the sphere. Every one of
them takes the target geometry as its argument, and `Base.parent` gives it back.

Read `Predicate(target)` as **cell RELATION target**. The DE9IM predicates
test the cell footprint; `CentroidCovered` tests the cell's centroid.

| Predicate | Keep a cell when… | Geometry or lon/lat extent | Spherical cap |
| --- | --- | --- | --- |
| `Intersects` | Cell and target share a point | Yes | Yes |
| `Disjoint` | Cell and target share no point | Yes | Yes |
| `Within` | The cell lies within the target | Yes | Yes; cap boundary contact is accepted |
| `Contains` | The target lies within the cell | Yes | No |
| `CoveredBy` | The target covers the cell, including its boundary | Yes | No |
| `Covers` | The cell covers the target, including its boundary | Yes | No |
| `Touches` | Boundaries meet without interior intersection | Yes | No |
| `Overlaps` | Same-dimensional interiors partly overlap | Yes | No |
| `Equals` | Cell and target are topologically equal | Yes | No |
| `Crosses` | Not implemented by this query engine | No | No |
| `CentroidCovered` | The cell's centroid lies on or inside the target | Yes | Yes |

For polygon targets, `Within` and `Contains` require the relevant interiors to
intersect. `CoveredBy` and `Covers` include boundary-only containment.
Thus `Within(zone)` asks for cells inside a zone; `Contains(zone)` asks for
cells large enough to contain the zone. `CoveredBy(zone)` and `Covers(zone)`
reverse direction in the same way.

`CentroidCovered(zone)` selects cells by centroid containment, the raster zonal
rule `boundary = :center` of [Rasterize, extract, and zonal](raster-work.md);
`Within(zone)` selects a subset of it and `Intersects(zone)` a superset.

Unsupported query predicates raise `ArgumentError`.

```@docs
DE9IMPredicate
Intersects
Disjoint
Contains
Within
Covers
CoveredBy
Touches
Crosses
Overlaps
Equals
CentroidCovered
```

## Multi-order answers

`MultiOrderCoverage` queries a region using cells at several levels and returns
a `MultiOrderCellSet`. See [Multi-order coverage](../tutorials/multiorder.md)
for a worked example.

```@docs
MultiOrderCoverage
MultiOrderCellSet
level_ranges
cellindices
iscontained
coarsest_contained
cell_polygons
```

## Region algebra

Use `expand` to obtain cells at a chosen level and `compact` to merge complete
sibling groups. [`grow`](@ref), documented with [region boundaries](boundaries.md),
adds neighbouring cells to a region.

```@docs
expand
compact
```

## The DimensionalData layer

A cube over cells carries a [`Cells`](@ref DiscreteGlobalGrids.CellLookups.Cells) dimension, whose lookup is a
[`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup). [`Covering`](@ref) is the selector that turns a geometry
into a selection on that dimension, so `A[Cells(Covering(geom))]` is the cube
spelling of [`query`](@ref).

```@docs
Cells
Covering
predicate_indices
```

## Index

```@index
Pages = ["api/selecting-cells.md"]
```
