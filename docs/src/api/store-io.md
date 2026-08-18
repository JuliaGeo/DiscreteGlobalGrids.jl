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

Both methods live in an extension on Zarr.jl: until `using Zarr` loads it the
two names are stubs, and calling one gets an error saying exactly that. What is
rendered below is the stub's docstring; the extension's methods carry a longer
one of their own, which `?dggread` prints once Zarr is loaded.

Between the store's attributes and the cube stands one plain-data value, the
[`StoreDescription`](@ref): grid name, level, encoding, array names, grid
parameters. Reading is attrs → description → axis and writing is
axis → description → attrs, and the two halves are deliberately separate
because they fail differently. Deciding **what a store says** is metadata
logic over a [`StoreSnapshot`](@ref) — attribute dictionaries and an array
listing, no chunks and no values — which is why a
[`DGGSConvention`](@ref) is written and tested without a store existing at all.
Deciding **how the ids are stored** is the other axis entirely, and it is a
[`CellEncoding`](@ref): encodings own layout mechanics and the grid owns id
arithmetic, so a new grid makes every encoding work and a new encoding works on
every grid.

Both are registries rather than closed sets. A downstream package that speaks a
fourth dialect calls [`register_convention!`](@ref); one that stores its ids
some fifth way calls [`register_encoding!`](@ref). Neither requires a change
here.

What a reader gets back is a [`ChunkedCellLookup`](@ref): the axis a store
wrote, which answers `At`, `Contains` and `Covering` the way an ordinary
`CellLookup` does but resolves them through the [`ChunkManifest`](@ref) — the
chunk grid described in cells — rather than by scanning ids. That is the whole
point of the layer: a store written by [`dggwrite`](@ref), or one whose ids are
arithmetic rather than stored, opens without reading a single coordinate chunk
however many cells it holds, and a selection over it fetches the chunks it lands
in and no others. A foreign dense store pays one pass over its ids at open, and
that pass is what proves them.

Refusing to guess is policy. A grid name in no registry, two conventions that
disagree about the level, an id that names no cell, a length that does not
check out: all of them raise [`DGGSFormatError`](@ref) naming the check that
failed, rather than being decoded into something plausible. The escape is
explicit — `dggread(store; description = StoreDescription(...))` asserts the
grid yourself and skips detection, which is also how an attribute-less store is
read.

## Reading and writing

```@docs
dggread
dggwrite
```

## The stored axis

```@docs
ChunkedCellLookup
DiscreteGlobalGrids.ChunkedCellVector
DiscreteGlobalGrids.axisposition
ChunkManifest
DiscreteGlobalGrids.chunkmanifest
nchunks
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

An encoding is how the cell ids reach the disk; the grid reference table is how
a store's spelling of a grid becomes a system this package can compute on. Both
are lookup tables with a registration function, and both are deliberately
strict about what they do not recognise. The types and the two tables are
exported, so a downstream package registers an encoding or a grid without
qualifying a name; the verbs it implements stay qualified, as a convention's do.

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

An encoding implements three things here: build the axis, declare write
eligibility, and name itself for the store's vocabulary. It asks the grid for
everything about the ids themselves, which is the layering rule — a grid that
answers the four id functions below works under every encoding, and an encoding
written against them works on every grid.

Reaching a *store* takes two more, both in the Zarr extension and both still
private: one `storedaxis` method, which says which array to open and how much of
it to read before handing it to [`cellaxis`](@ref), and — to be written as well
as read — the four verbs the write pipeline dispatches on. An encoding that
registers without them is refused by name, with
`DGGSFormatError(check = :unsupported_encoding)`, rather than by a `MethodError`
about a private function.

```@docs
DiscreteGlobalGrids.cellaxis
DiscreteGlobalGrids.write_eligible
DiscreteGlobalGrids.encodingname
DiscreteGlobalGrids.idrank
DiscreteGlobalGrids.idselect
DiscreteGlobalGrids.idcount_between
DiscreteGlobalGrids.idvalid
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
