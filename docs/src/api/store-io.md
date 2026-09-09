# Reading and writing DGGS stores

```@meta
CurrentModule = DiscreteGlobalGrids
```

A DGGS store is a Zarr group containing variables over one cell axis and the
metadata needed to interpret that axis. [`dggwrite`](@ref) writes a
DimensionalData `DimStack`; [`dggread`](@ref) reopens it as a `DimStack` with a
shared `Cells` dimension and lazy arrays.

The methods live in the Zarr.jl extension, which `using Zarr` loads. The
extension supplies the store-aware methods documented here.

[`StoreDescription`](@ref) carries plain data between store attributes and the
cube: grid name, level, encoding, array names, and grid parameters. Reading
follows attrs → description → axis; writing follows axis → description → attrs.
[`StoreSnapshot`](@ref) supplies metadata-only input to
[`DGGSConvention`](@ref), while [`CellEncoding`](@ref) supplies the identifier
layout. The grid retains id arithmetic, so encodings work across systems.

Both are registries. A downstream package can add a metadata dialect with
[`register_convention!`](@ref) or an id layout with
[`register_encoding!`](@ref).

The reader returns a [`ChunkedCellLookup`](@ref). It answers `At`, `Contains`
and `Covering` through the [`ChunkManifest`](@ref), which describes the chunk
grid in cells. Arithmetic and range encodings can open without coordinate
reads; a foreign dense store reads its ids once at open to validate them.

A stored axis is also a **region** and answers `halo`,
`border`, `interior` and `adjacency` with the same code as an in-memory axis. It
does so through [`region`](@ref), which is the axis's compressed
[`CellVector`](@ref) twin, built on the first call and kept. What that
conversion costs is the encoding's and not the axis's length: a ranges or
implicit store converts by arithmetic alone, because a stored interval is a run
of consecutive ranks and a rank plus one is an index; a dense store reads its
ids once, in the order that touches each chunk once. Index order is
preserved either way, which is what lets a result computed through the twin be
written back against the store's own axis with no permutation. Sweeping a store
along its own chunk lines is
[its own page](@ref "Sweeping a cube along its chunk lines").

Validation is strict. An unknown grid, conflicting metadata, invalid id, or
inconsistent length raises [`DGGSFormatError`](@ref) with the failed check.
Supply `description = StoreDescription(...)` to assert the metadata explicitly
when a store has no attributes.

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

## Encodings and grid references

An encoding maps cell ids to disk storage; a grid reference maps the store's grid
name to a system. Both are lookup tables with registration functions and strict
recognition rules. The exported types and tables let downstream packages add
entries.

```@docs
CellEncoding
DenseEncoding
RangesEncoding
ImplicitEncoding
ENCODING_REGISTRY
register_encoding!
GridReference
GRID_REFERENCE
register_grid!
DiscreteGlobalGrids.gridreference
```

An encoding builds the axis, declares write eligibility, and supplies its store
name. It asks the grid for id arithmetic, so one encoding works across systems
and one system works with every encoding.

The Zarr extension uses `storedaxis` to open an axis and dispatches four write
operations for encodings that support output. An incomplete registration raises
`DGGSFormatError(check = :unsupported_encoding)`.

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
