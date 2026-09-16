# Workflow execution details

```@meta
CurrentModule = DiscreteGlobalGrids
```

This reference holds detailed storage, scheduling, and metadata contracts.
Start with [Neighbours and stencils](../api/neighbors.md),
[Requesting neighbour fields](../api/neighbor-fields.md), or
[Reading and writing DGGS stores](../api/store-io.md) for normal use.

All neighborhood sweeps described here use one-ring callbacks. Input halo
width does not change that reach. Cache behavior is an implementation detail;
cell identity and callback/output alignment remain the public contract.

## Cell-vector sweep execution

mapneighbors(f, cv; order = StorageOrder(), threaded = true,
                 connectivity = Vertex())
    mapneighbors(f, cv, data::AbstractVector; ...)
    mapneighbors(f, cv; needs = (Value(data), Centroid()), ...)

Apply `f` to each cell and its clipped one-ring. `cv` may be a
[`CellVector`](@ref), [`PartialGrid`](@ref), or [`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup).

Without `data`, `f(cell, nbrs)` receives the same indexed handles yielded
by the one-argument [`neighbors`](@ref) iterator. With a vector laid out
against the subset, `f(cell, value, values)` receives the cell value and its
neighbour values in the counter-clockwise order [`neighbors`](@ref) states, so
slot `j` of the callback's ring names a direction.

`needs` names the per-neighbour fields the kernel reads — a tuple of `Cell`,
`Index`, `Value` and `Centroid` requests — and the callback becomes
`f(center, rings)`: one entry per need for the visited cell, and one ring per
need for its clipped neighbours. The rings are field-major, `rings[j]` being
need `j`'s value for every neighbour, with slot `i` of every ring naming the
same neighbour; a caller who wants one record per neighbour writes
`zip(rings...)`. `Index(Local())` is the index in the collection passed here,
and it stays that index however the sweep is run: [`mapneighbors!`](@ref DiscreteGlobalGrids.mapneighbors!)
answers the same request chunk by chunk and translates each chunk's own
numbering back to this collection's, so its result is this one's cell for
cell. `Centroid()` is answered from a bounded working set kept per task and
keyed by the local index, so a centroid several neighbourhoods name is
computed once wherever the visit order keeps them close in that index — the
default storage order does,
and a random permutation `order` does not. A field request and a positional
`data` vector are exclusive.

Results are stored in subset index order. A concrete tuple result produces
a tuple of vectors, one per component. `order` accepts [`StorageOrder`](@ref)
or a permutation of `1:length(cv)`; invalid permutations throw
`ArgumentError`.

When `threaded` is true, contiguous ranges run in separate tasks and write to
disjoint output indices — legal exactly when `f` is order-independent, and
the results are then identical to the sequential ones. A callback that throws
there raises one [`NeighborCallbackError`](@ref) naming the cell and index
it failed at, not one exception per task.
[`foreachneighbors`](@ref) provides the side-effecting form and defaults to
sequential execution.

## Dimensional-array sweep execution

mapneighbors(f, A::AbstractDimArray; spatialdim = nothing, pass = Neighbors(),
                 order = StorageOrder(), threaded = true, connectivity = Vertex())
    mapneighbors(f, A::AbstractDimArray; needs = (Value(a), Centroid()), ...)

Apply `f` to each cell and its neighbors. The result uses `A`'s wrapper and
lookups. If `f` returns a concrete tuple, each component becomes an array.

`spatialdim` accepts any selector supported by `DimensionalData.dims`, and
`mapneighbors(f, A, dims; kw...)` is the same request spelled positionally.
By default, the first dimension with a [`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup) is used. An
array without one, or a selector that misses or names a non-cell dimension,
is an `ArgumentError`.

`pass` controls the callback arguments and output shape:

- [`Neighbors`](@ref): indexed handles; one result per cell, on the cell
  dimension.
- [`Values`](@ref): scalar values; the same dimensions as `A`.
- [`NeighborSlices`](@ref): views across the other dimensions; one result per
  cell, on the cell dimension.

`needs` names the per-neighbour fields the kernel reads instead of `pass`, and
the callback becomes `f(center, rings)` — the contract is the [`CellVector`](@ref)
method's, and the requests are `Cell`, `Index`, `Value` and `Centroid`. Values
reach the callback through the request's `Value` entries rather than from `A`
itself, so an array of any dimensionality gives one result per cell, on the
cell dimension. Any `pass` other than the default alongside `needs` is an
`ArgumentError`.

With `pass = Values()` or a `needs` request, and `order = StorageOrder()`, a
cube whose data is chunked on disk is swept along those chunks rather than cell
by cell — see [`chunkplan`](@ref DiscreteGlobalGrids.chunkplan). The result is identical either way; what
changes is that each stored chunk is decoded once instead of once per scalar
read, and that a request's stored `Value`s are read the same way.
`Index(Local())` still answers this cube's cell-axis index, never a chunk's
own. A permutation `order` names a visit order over the whole axis and a
chunked sweep visits by chunk, so the two cannot both be honoured and the
permutation wins.

## Stored one-ring execution

mapneighbors!(dest, f, A::AbstractDimArray; halo = 1, chunks = :auto,
                  spatialdim = nothing, connectivity = Vertex(), threaded = true)
    mapneighbors!(dest, f, A, plan::MapChunkPlan;
                  needs = (Value(a), Centroid()), threaded = true)

Apply `f` to each cell of `A` and its neighbours, chunk by chunk, writing one
result per cell into `dest`.

This is [`mapneighbors`](@ref)'s out-of-core form, and the difference is where
the results go: `mapneighbors` collects them, which needs one array of them in
memory, and this writes them into `dest` a chunk at a time, which does not.
`dest` is anything indexable along the cell dimension by a range — an `Array`,
another cube, or a lazy array being written back to a store.

`f` is called as `f(cell, value, neighbor_values)`, [`Values`](@ref)' form: a
chunked sweep is exactly the case where the values must flow through the
traversal rather than be fetched by the callback.

`needs` names the per-neighbour fields the kernel reads instead, and the
callback becomes `f(center, rings)` — [`mapneighbors`](@ref)'s field-request
contract, here on the chunk route. The request is stated once, about the cube
the caller passed, and translated onto each chunk: `Index(Local())` answers
that cube's cell-axis index on every chunk, never a chunk-local one, and a
[`Value`](@ref) over a stored array is read along its own storage chunks the
way the swept data is. `dest` then holds one result per cell.

The result is the whole-axis sweep's, cell for cell, in either form. Each
chunk's halo carries every axis neighbour of every cell the chunk owns, so a
ring computed on a chunk is the ring computed on the axis — clipped
identically, and in the same order.

`threaded` threads WITHIN each chunk. To run chunks themselves in parallel,
`split` a [`chunkplan`](@ref DiscreteGlobalGrids.chunkplan) and call this on the pieces. Build the plan before
splitting: doing so is what fills the axis's [`region`](@ref) memo, so the
pieces share one conversion rather than each repeating it.

## Store read mechanics

dggread(store; vars = All(), lazy = true, validate = :strict,
            conventions = CONVENTION_REGISTRY, description = nothing) -> DimStack
    dggread(store, var::Symbol; kwargs...) -> DimArray

Requires `using Zarr`. The following details supplement the concise public function contract.

Read a DGGS store into plain DimensionalData: one `Cells` dimension shared by
every layer, carrying a [`ChunkedCellLookup`](@ref DiscreteGlobalGrids.ChunkedLookups.ChunkedCellLookup) — the lookup over an axis a
store wrote, which resolves a cell without scanning it. The grid SYSTEM is in
that lookup's type, the level is a field of the grid it holds, and what is
neither — orientation, ellipsoid, the layout the store keeps its axis in — rides
in the [`StoreDescription`](@ref) under the stack's `metadata["description"]`.

`store` is a `Zarr.ZGroup`, a `Zarr.AbstractStore`, a local path, or a URL
(`gs://`, `s3://`, `https://`); `s3://` also requires `using AWSS3`. `vars = All()` reads every data variable, or
name the `Symbol`s to read. Data arrays are lazy by default; `lazy = false`
materializes them. The detected convention, the verbatim original attributes and
the source encoding ride in the stack's `metadata`, which is enough to
regenerate a value-identical store.

`validate = :strict` checks that every stored id names a cell of the declared
level and `:lazy` samples instead. Neither reaches a store carrying a chunk
manifest this package wrote: that axis is built from the manifest and no id is
scanned, which is what opens a store of tens of millions of cells at all.
`validate = :scan` declines the manifest and runs the full scan on any store.

`description` bypasses detection: pass a [`StoreDescription`](@ref) and the
caller asserts grid, level, encoding and array names, leaving only the
mechanical checks. That is how an attribute-less store is read.

## Store write mechanics

dggwrite(dest, stack_or_array; encoding = :auto,
             conventions = DEFAULT_WRITE_CONVENTIONS, chunks = :auto,
             merge = :step, chunk_target = 1_000_000) -> dest

Requires `using Zarr`. The following details supplement the concise public function contract.

Write a `DimStack` or `DimArray` over a `Cells` dimension to a Zarr v2 directory
store. The cell dimension has to carry a `CellLookup` or a
[`ChunkedCellLookup`](@ref DiscreteGlobalGrids.ChunkedLookups.ChunkedCellLookup) — this package's way of saying the axis is still
sorted, unique and at one level. `dest` is a local directory path or a writeable
`Zarr.ZGroup`; a remote URL is refused rather than half-written — write locally
and upload.

`encoding = :auto` picks ranges where the axis is eligible — sorted, unique,
single-level — and dense otherwise; `:dense` is the interop escape for readers
that cannot expand ranges, `:ranges` forces the compact form, and `:implicit`
writes no cell coordinate at all, which needs a whole level. `conventions`
stamps the store, dual by default so that both a convention-aware reader and
xdggs can open it.

`merge` picks the ranges run rule: `:step` (default) merges unit-increment ids,
which a structural reader also counts correctly; `:rank` merges rank-adjacent
cells for the fewest rows, and is read back correctly only by a rank-aware
reader such as this package. `chunks = :auto` aims each chunk at
`chunk_target` as a whole number of complete coarse-ancestor subtree runs; an
`Integer` fixes the chunk length in cells instead. `chunk_target` counts the
ELEMENTS of a chunk — cells times the extents of the non-cell dimensions, which
are one chunk each — so a layer with a 40-step time axis gets a fortieth of the
cells per chunk.

Each layer's `metadata` is written as its array attributes and the stack's
`metadata["attrs"]` as the group's, the two places `dggread` puts them, so a
store read and rewritten keeps its `units`, `long_name` and group vocabulary;
convention-generated keys are stamped over the producer's. A round trip
normalizes two things: layers are written in alphabetical order, and each
layer's attributes carry the `_ARRAY_DIMENSIONS` this writer stamps.

`layout` chooses the SHAPE of the store rather than the shape of its cell
coordinate: `:cells` (the default) is everything above, a one-dimensional cell
axis; `:subzones` is the two-dimensional [`SubzoneLayout`](@ref), which takes an
`ancestor_level` and none of the keywords above it.

## Zarr reader metadata and validation

dggread(store; vars = All(), lazy = true, validate = :strict,
            conventions = CONVENTION_REGISTRY, description = nothing) -> DimStack
    dggread(store, var::Symbol; kwargs...) -> DimArray

Read a DGGS store into plain DimensionalData: one `Cells` dimension shared by
every layer, carrying a [`ChunkedCellLookup`](@ref DiscreteGlobalGrids.ChunkedLookups.ChunkedCellLookup) that resolves a cell to a
index without scanning the axis. Any other dimension of the store — time,
bands — becomes an ordinary `Dim`, with the values of the like-named coordinate
array where the store has one.

`store` is a `Zarr.ZGroup`, a `Zarr.AbstractStore`, a local path, or a URL. A
`gs://BUCKET/PATH` URL is read as `https://storage.googleapis.com/BUCKET/PATH`,
which works for public-read buckets; an `s3://` URL needs Zarr's own AWSS3
extension (`using AWSS3`) and says so otherwise.

  - `vars`: `All()` for every data variable, or the `Symbol`s to read. An
    unknown name raises, listing what the store holds.
  - `lazy`: leave the data as the store's own chunked arrays (the default), or
    materialize them.
  - `validate`: `:strict` (the default) checks that every stored id names a cell
    of the declared level; `:lazy` checks `LAZY_SAMPLES` per chunk
    instead. Sortedness, uniqueness and the length checks are not optional and
    run either way. None of the three reaches a store carrying a chunk manifest
    this package wrote: its axis is built from the manifest without a scan, so
    no id is checked at all, and what bounds that trust is a per-chunk
    comparison of first id, last id and length as the ids are read
    (`persistedmanifest`). `:scan` is how that trust is declined — the
    sidecar is ignored, the ids are scanned, and every check runs on every
    store.
  - `conventions`: the conventions to try, in order.
  - `ancestors`: an ancestor-subzone store only — the level-`ancestor_level`
    cells (or column indices) to restrict the cell axis to. The default reads
    the WHOLE level, since a column nobody wrote is not absent from such a store,
    it reads back as fill.
  - `description`: a [`StoreDescription`](@ref) that bypasses detection. The
    caller then asserts grid, level, encoding and array names, and only the
    mechanical checks — the id scan, the closed-form counts — still run. This
    is how an attribute-less store is read.

The stack's `metadata` carries the provenance a value-identical rewrite needs:

| key | |
|---|---|
| `"source"` | the store URL or path |
| `"conventions"` | the conventions that fired, in order |
| `"encoding"` | the cell-axis layout, as the store spells it |
| `"attrs"` | the group attributes verbatim |
| `"description"` | the [`StoreDescription`](@ref) everything was read through |

Each layer keeps its own array attributes as its metadata.

**The ancestor-subzone layout** is recognized from its own attributes and read
by its own path: the layers come back over a `Cells` dimension as ever, backed
by a lazy `SubzoneCellArray` whose chunks are the store's subtrees. Such
a stack's metadata carries `"layout"` — a [`SubzoneLayout`](@ref) — in place of
`"description"`, `"conventions"` and `"encoding"`, none of which have anything
to say about a two-dimensional store.


## Zarr writer metadata and validation

dggwrite(dest, stack_or_array; encoding = :auto,
             conventions = DEFAULT_WRITE_CONVENTIONS, chunks = :auto,
             merge = :step, chunk_target = DEFAULT_CHUNK_TARGET) -> dest

Write a `DimStack` or `DimArray` over a cell axis to a **Zarr v2 directory
store**, consolidated metadata included. `dest` is a local path or an open
writeable `Zarr.ZGroup`; a `gs://`/`s3://`/`https://` URL is refused rather than
half-written.

The cell dimension must carry an `AbstractCellLookup`, which is
this package's way of saying the axis is still sorted, unique and at one level;
`reverse` and friends degrade it to a `Categorical`, and that is refused.

  - `encoding = :auto` writes ranges where the axis is eligible and dense
    otherwise. `:dense` is the interop escape for readers that cannot expand
    ranges, `:ranges` forces the compact form, and `:implicit` writes no cell
    coordinate at all — the index is the cell — which needs a whole level.
  - `merge = :step` merges only ids adjacent as integers, so no interval can
    enclose an id that names no cell — what a structural-count reader needs, and
    what the published IGEO7 range stores hold. `merge = :rank` merges runs of
    consecutive CELLS instead, giving the fewest rows, and is read back correctly
    only by a rank-aware reader such as this package; see `idranges`.
  - `chunks = :auto` groups whole coarse-ancestor subtree runs into chunks of
    about `chunk_target` elements; an integer is a fixed chunk length in CELLS.
    See `WriteChunkPlan` for what that guarantees and what it only aims at.
    `chunk_target` counts the elements of a chunk — cells times the extents of
    the non-cell dimensions, which are one chunk each — so a layer with a
    40-step time axis gets a fortieth of the cells per chunk.
  - `conventions` stamps the store, `zarr-conventions/dggs` plus xdggs by
    default, so both a convention-aware reader and xdggs can open it.

The chunk grid is persisted as a `(n_chunks, 2)` sidecar array of per-chunk
first and last id, so a reader need not scan the axis to rebuild it.

**Attributes.** Each layer's `metadata` is written as its array attributes and
the stack's `metadata["attrs"]` as the group's — the two places `dggread` puts
them, so a store read and rewritten keeps its `units`, its `long_name` and its
group vocabulary. Convention-generated keys are stamped OVER the producer's: a
`_ARRAY_DIMENSIONS` or a `dggs` object carried in from another layout would
describe this store wrongly. Other stack metadata is not written; the cell
coordinate is regenerated by the encoding and carries no producer attributes.

**Two documented normalizations of a round trip.** Layers are written in
alphabetical order, so that is the order they come back in whatever order went
in; and a layer's attributes include the `_ARRAY_DIMENSIONS` this writer stamps,
so a stack read back carries it in each layer's `metadata`.

Layers are never overwritten: a `ZGroup` destination that already holds an array
this write would create raises before anything is stamped.
