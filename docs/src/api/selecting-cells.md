# Selecting cells

```@meta
CurrentModule = DiscreteGlobalGrids
```

Selection has one verb and three answers. [`query`](@ref) asks a grid or a
system which of its cells stand in a named relation to a geometry;
[`covering`](@ref) is the same question phrased as "give me the cells", and
[`covering_indices`](@ref) is it phrased as "give me their indices".

The relation is a DE9IM predicate — [`Intersects`](@ref), [`Within`](@ref),
[`Disjoint`](@ref) and the rest — read with spherical semantics rather than
planar ones. The answer comes back as a [`CellVector`](@ref) at one level, or as
a [`MultiOrderCellSet`](@ref) spanning several when the query is allowed to
coarsen.

## Asking

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

A coverage that may span levels, and the verbs that read one. The set types
themselves — [`MultiOrderCoverage`](@ref) and [`MultiOrderCellSet`](@ref) — are
documented with the [region boundaries](boundaries.md).

```@docs
level_ranges
cellindices
iscontained
coarsest_contained
cell_polygons
```

## Region algebra

[`grow`](@ref) is documented with the boundary verbs it is built on; these two
complete the trio.

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
