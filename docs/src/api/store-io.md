# Reading and writing DGGS stores

```@meta
CurrentModule = DiscreteGlobalGrids
```

A DGGS store is a Zarr group containing variables over one cell axis and the
metadata needed to interpret that axis. [`dggwrite`](@ref DiscreteGlobalGrids.dggwrite) writes a
DimensionalData `DimStack`; [`dggread`](@ref DiscreteGlobalGrids.dggread) reopens it as a `DimStack` with a
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

The reader returns a [`ChunkedCellLookup`](@ref DiscreteGlobalGrids.ChunkedLookups.ChunkedCellLookup). It answers `At`, `Contains`
and `Covering` through the [`ChunkManifest`](@ref), which describes the chunk
grid in cells. Arithmetic and range encodings can open without coordinate
reads; a foreign dense store reads its ids once at open to validate them.

A mixed-level axis — the one [`coarsen`](@ref) builds — is stored in the `compacted`
layout instead: `cell_ids` and `cell_levels` as aligned columns under
`refinement_level: null`, with the `refinement_levels` attribute naming the
level column. [`dggread`](@ref) restores that axis as a
[`MultiOrderLookup`](@ref). The layout extends v1 of `zarr-conventions/dggs`,
which specifies `compression: "none"` when `refinement_level` is null, so this
package supplies its reader; [`expand`](@ref) presents the same data at the one
level every other encoding needs. The [multi-order storage
tutorial](../tutorials/moc_storage.md) writes and reads such a store.

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

## Writing a store for xdggs

[xdggs](https://xdggs.readthedocs.io)'s default convention reads a
one-dimensional `cell_ids` coordinate whose attributes are the fields of its
grid-info dataclass, for a grid its registry knows. `target = :xdggs` writes
that layout and checks the store's description against it before the group is
created:

```julia
using DiscreteGlobalGrids, Zarr
dggwrite("tas.zarr", cube; target = :xdggs)
```

The store then opens in Python with no arguments beyond the path:

```python
import xarray as xr, xdggs
ds = xr.open_dataset("tas.zarr", engine="zarr").pipe(xdggs.decode)
ds.dggs.cell_centers()
```

What the target chooses and checks, and why:

| Choice | Reason |
|---|---|
| `encoding = :dense` | `:auto` prefers the ranges encoding, which stores `(n, 2)` intervals and no `cell_ids` array; xdggs finds nothing to decode there. |
| the coordinate is `cell_ids` and carries a level | `xdggs.decode` looks the coordinate up by that name, and its grid info has no default level. |
| grid in [`XDGGS_GRIDS`](@ref) | `healpix` and `h3` ship with xdggs; `igeo7` needs the `xdggs-dggrid4py` plugin from its main branch, whose grid info takes the attributes [`XdggsConvention`](@ref) writes. |

Two properties every dense store from this writer already has matter to xdggs:
the coordinate's attributes are grid keys only, because xdggs forwards every
attribute but `grid_name` to its dataclass constructor and a stray `units` is
a `TypeError` there; and `cell_ids` has no fill value, because xarray reads a
Zarr v2 fill value as a mask and promotes the ids to `Float64`.

Both write conventions are stamped, so `ds.dggs.decode(convention="zarr")`
opens the same store through the `zarr-conventions/dggs` group attributes.
`xdggs.decode` strips the grid attributes off the coordinate as it builds its
index, so a dataset edited in Python is written back with
`ds.dggs.encode("xdggs").to_zarr(...)`.

The chunk manifest sidecar rides along as a data variable on dimensions of its
own, which xdggs ignores.

```@docs
DiscreteGlobalGrids.require_xdggs_readable
DiscreteGlobalGrids.XDGGS_GRIDS
DiscreteGlobalGrids.xdggs_ellipsoid_attrs
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

A compacted store restores a mixed-level axis instead. `coarsen` chooses the
levels and `aggregate` reduces values onto them; the remaining verbs query the
mixed-level container the store holds. [Multi-order
coverage](../tutorials/multiorder.md) is the tutorial for them.

```@docs
MultiOrderLookup
MultiOrderVector
coarsen
aggregate
covering_index
complement
reference_level
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
CompactedEncoding
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
