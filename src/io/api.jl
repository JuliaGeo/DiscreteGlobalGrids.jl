# `dggread`/`dggwrite` are defined and exported here as stubs. Their methods
# live in `DiscreteGlobalGridsZarrExt`, so the store types, keyword defaults
# and IO all stay behind the Zarr weak dependency.

# Everything that matches no method of the extension lands back on the stub, so
# what it says has to depend on whether the extension is there: telling a caller
# to run `using Zarr` when Zarr is already loaded sends them to fix the one thing
# that is not wrong.
_needs_zarr(f) = error(_no_zarr_method(f,
    Base.get_extension(DiscreteGlobalGrids, :DiscreteGlobalGridsZarrExt) !== nothing))

function _no_zarr_method(f, loaded::Bool)
    loaded && return """
        no `$f` method matches these arguments; the Zarr extension is loaded, so \
        this is an argument-type problem and not a missing package. A store is a \
        `Zarr.ZGroup`, a `Zarr.AbstractStore`, a path or a URL, and the cube \
        `dggwrite` takes is a `DimArray` or `DimStack` over a cell dimension."""
    return """
        `$f` requires Zarr.jl. Run

            using Zarr

        to load `DiscreteGlobalGridsZarrExt`, which provides the store methods."""
end

"""
    dggread(store; vars = All(), lazy = true, validate = :strict,
            conventions = CONVENTION_REGISTRY, description = nothing) -> DimStack
    dggread(store, var::Symbol; kwargs...) -> DimArray

**Requires `using Zarr`.** The methods live in `DiscreteGlobalGridsZarrExt`;
until it loads, this name is a stub whose only behaviour is to say so, and the
extension's docstring — not this one — is the full keyword reference.

Read a DGGS store into plain DimensionalData: one `Cells` dimension shared by
every layer, carrying a [`ChunkedCellLookup`](@ref) — the lookup over an axis a
store wrote, which resolves a cell without scanning it. The grid SYSTEM is in
that lookup's type, the level is a field of the grid it holds, and what is
neither — orientation, ellipsoid, the layout the store keeps its axis in — rides
in the [`StoreDescription`](@ref) under the stack's `metadata["description"]`.

`store` is a `Zarr.ZGroup`, a `Zarr.AbstractStore`, a local path, or a URL
(`gs://`, `s3://`, `https://`). `vars = All()` reads every data variable, or
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
"""
dggread(args...; kwargs...) = _needs_zarr("dggread")

"""
    dggwrite(dest, stack_or_array; encoding = :auto,
             conventions = DEFAULT_WRITE_CONVENTIONS, chunks = :auto,
             merge = :step, chunk_target = 1_000_000) -> dest

**Requires `using Zarr`.** The methods live in `DiscreteGlobalGridsZarrExt`;
until it loads, this name is a stub whose only behaviour is to say so, and the
extension's docstring — not this one — is the full keyword reference.

Write a `DimStack` or `DimArray` over a `Cells` dimension to a Zarr v2 directory
store. The cell dimension has to carry a `CellLookup` or a
[`ChunkedCellLookup`](@ref) — this package's way of saying the axis is still
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
"""
dggwrite(args...; kwargs...) = _needs_zarr("dggwrite")

"""
    subzonestore(dest, system, level; ancestor_level, layers, kwargs...) -> SubzoneStore
    subzonestore(dest) -> SubzoneStore

**Requires `using Zarr`.** The methods live in `DiscreteGlobalGridsZarrExt`,
whose docstring is the full keyword reference.

Create — or reopen — an ancestor-subzone store for INCREMENTAL writing: the
group, its arrays and its attributes are stamped once, and the columns are
filled afterwards, one [`dggwrite!`](@ref) at a time. That is the production
shape of a global regrid, where one ancestor subtree is one task's worth of work
and the tasks finish in no particular order.

A column is one Zarr chunk and therefore one file, so writes of DISJOINT columns
do not touch a byte in common and need no coordination. Nothing shared is
rewritten by a column write — not the attributes, not the consolidated metadata,
not a manifest — which is what makes that true.

See [`SubzoneLayout`](@ref) for the layout itself and [`dggwrite`](@ref)'s
`layout = :subzones` for the one-shot form.
"""
subzonestore(args...; kwargs...) = _needs_zarr("subzonestore")

"""
    dggwrite!(store::SubzoneStore, ancestor, values; var = the only layer) -> store
    dggwrite!(store::SubzoneStore, cube) -> store

**Requires `using Zarr`.** The methods live in `DiscreteGlobalGridsZarrExt`.

Fill columns of a store [`subzonestore`](@ref) has already created: one ancestor
cell's subtree from a vector in ascending cell id, or every complete column of a
cube over a cell axis.

`values` is as long as that ancestor's subtree really is — `7^d` for a hexagon
and `(5*7^d + 1)/6` for a pentagon — and the rest of the column stays fill.
"""
dggwrite!(args...; kwargs...) = _needs_zarr("dggwrite!")
