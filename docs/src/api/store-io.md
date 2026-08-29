# Reading and writing DGGS stores

```@meta
CurrentModule = DiscreteGlobalGrids
```

A DGGS store is a Zarr group holding data variables over one cell axis, plus
attributes saying which grid, which level and how the cell ids are laid out.
[`dggread`](@ref) turns that into plain DimensionalData — one `Cells` dimension
shared by every layer, lazy arrays behind it — and [`dggwrite`](@ref) turns it
back. Neither introduces a container type: what a producer hands over is a
`DimStack`, and what a consumer gets is a `DimStack`.

Zarr.jl loads both methods through an extension. Before `using Zarr`, the two
names are stubs that explain how to load their implementations. The reference
below renders the stub docstrings; `?dggread` displays the longer extension
docstring after Zarr loads.

[`StoreDescription`](@ref) connects store attributes to the cube. It records the
grid name, level, encoding, array names and grid parameters. Reading follows
attributes → description → axis; writing follows axis → description →
attributes.

Two independent abstractions support those pipelines:

  - [`DGGSConvention`](@ref) interprets metadata in a [`StoreSnapshot`](@ref),
    which contains attribute dictionaries and an array listing.
  - [`CellEncoding`](@ref) implements cell-id storage layout while the grid
    implements id arithmetic.

This separation lets conventions run in tests without store chunks or values,
and lets every grid and encoding combination share the same implementation.

Both are registries rather than closed sets. A downstream package that speaks a
fourth dialect calls [`register_convention!`](@ref); one that stores its ids
some fifth way calls [`register_encoding!`](@ref). Neither requires a change
here.

[`coarsen`](@ref) builds a mixed-level axis. [`dggwrite`](@ref) stores that
axis in the `compacted` layout as two aligned columns, `cell_ids` and
`cell_levels`, under `refinement_level: null`; the `refinement_levels`
attribute names the second column. [`dggread`](@ref) restores it as a
[`MultiOrderLookup`](@ref).

This package defines `compacted` as an extension to v1 of
`zarr-conventions/dggs`. Version 1 specifies `compression: "none"` when
`refinement_level` is null, so this package supplies the reader for the extended
layout. [`expand`](@ref) presents mixed-level data at the single level required
by the other encodings.

A [`ChunkedCellLookup`](@ref) resolves `At`, `Contains` and `Covering` through a
[`ChunkManifest`](@ref), the chunk grid expressed in cells. Stores written by
[`dggwrite`](@ref) and stores with arithmetic axes open without reading a
coordinate chunk. Selections then fetch only the chunks they touch. Opening a
foreign dense store scans its ids once to validate them and build the manifest.

A stored axis also supports `halo`, `border`, `interior` and `adjacency` as a
region. The first [`region`](@ref) call builds and caches a compressed
[`CellVector`](@ref) representation:

  - Ranges and implicit encodings build it with rank arithmetic.
  - Dense encoding reads the ids once in chunk order.

Both paths preserve index order, so region results align directly with the
stored axis. Sweeping a store along its own chunk lines is
[its own page](@ref "Sweeping a cube along its chunk lines").

[`DGGSFormatError`](@ref) names the failed check for an unknown grid, conflicting
conventions, an invalid cell id or a length mismatch. To read an attribute-free
store, pass `dggread(store; description = StoreDescription(...))`; the supplied
description asserts the grid and skips detection.

## Reading and writing

```@docs
dggread
dggwrite
```

## The stored axis

```@docs
ChunkedCellLookup
DiscreteGlobalGrids.ChunkedCellVector
DiscreteGlobalGrids.axisindex
ChunkManifest
DiscreteGlobalGrids.chunkmanifest
nchunks(::ChunkManifest)
chunkof
chunkbounds
```

## Describing a store

```@docs
StoreDescription
StoreSnapshot
ArrayEntry
describe_store
Detection
DiscreteGlobalGrids.Ellipsoid
DiscreteGlobalGrids.GridOrientation
DiscreteGlobalGrids.ellipsoid
DiscreteGlobalGrids.DEFAULT_ELLIPSOID
DiscreteGlobalGrids.ellipsoid_attrs
```

## Conventions

A convention is a dialect of store attributes. Several may fire on one store —
the published stores are stamped twice on purpose — and where they do, their
descriptions must agree field for field.

```@docs
DGGSConvention
ZarrDGGSConvention
XdggsConvention
LegacyHealpixConvention
DKRZConvention
CONVENTION_REGISTRY
DEFAULT_WRITE_CONVENTIONS
register_convention!
```

A new dialect is a subtype of [`DGGSConvention`](@ref) with a `detect`, a
`decode` and — if it is to be written and not only read — an `encode!`. These
names stay qualified: they are generic enough that exporting them would collide
with half the ecosystem.

```@docs
DiscreteGlobalGrids.detect
DiscreteGlobalGrids.decode
DiscreteGlobalGrids.encode!
DiscreteGlobalGrids.conventionname
DiscreteGlobalGrids.gridname
```

## Encodings and grid names

An encoding determines how cell ids reach disk. The grid reference table maps a
stored grid name to a system this package can compute on. Both lookup tables
require registered names. Their exported types and registration functions let a
downstream package register an encoding or grid; implementation verbs remain
qualified.

```@docs
CellEncoding
DenseEncoding
RangesEncoding
ImplicitEncoding
CompactedEncoding
ENCODING_REGISTRY
register_encoding!
GridReference
GRID_REFERENCE
register_grid!
DiscreteGlobalGrids.gridreference
```

An encoding implements three things here: build the axis, declare write
eligibility, and name itself for the store's vocabulary. It asks the grid for
everything about the ids themselves, which is the layering rule — a grid that
answers the five id functions below works under every encoding, and an encoding
written against them works on every grid.

The Zarr extension adds the store-specific methods. `storedaxis` chooses the
array and read size before calling [`cellaxis`](@ref). Writable encodings also
implement the four private write-pipeline verbs. Missing methods raise
`DGGSFormatError(check = :unsupported_encoding)` with the encoding name.

```@docs
DiscreteGlobalGrids.cellaxis
DiscreteGlobalGrids.write_eligible
DiscreteGlobalGrids.encodingname
DiscreteGlobalGrids.idrank
DiscreteGlobalGrids.idselect
DiscreteGlobalGrids.idcount_between
DiscreteGlobalGrids.idvalid
DiscreteGlobalGrids.idcell
```

## Errors

One exception type, defined layer-neutrally so that the encoding and lookup
layers — which never learn what a store is — can throw it too. The store URL
and the conventions that fired are optional context, added at the boundary by
the layer that does know them.

```@docs
DGGSFormatError
DiscreteGlobalGrids.with_store_context
DiscreteGlobalGrids.store_context
```

## Index

```@index
Pages = ["api/store-io.md"]
```
