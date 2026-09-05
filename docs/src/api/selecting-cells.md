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
planar ones. The answer comes back as a [`CellVector`](@ref) at one level, or as
a [`MultiOrderCellSet`](@ref) spanning several when the query is allowed to
coarsen.

## Cell collections

A `CellVector` holds cells at one level, and a `CellLookup` connects them to a
`Cells` dimension. `region` obtains the cell collection used by regional
operations, including for a stored axis.

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

The types are DE9IM.jl's; the semantics are this package's, evaluated on the
sphere. Every one of them takes the target geometry as its argument, and
`Base.parent` gives it back.

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

A cube over cells carries a [`Cells`](@ref) dimension, whose lookup is a
[`CellLookup`](@ref). [`Covering`](@ref) is the selector that turns a geometry
into a selection on that dimension, so `A[Cells(Covering(geom))]` is the cube
spelling of [`query`](@ref).

```@docs
Cells
Covering
predicate_indices
NeighborSlices
```

## Index

```@index
Pages = ["api/selecting-cells.md"]
```
