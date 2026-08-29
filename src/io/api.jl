# Stubs keep Zarr types and IO behind the optional Zarr extension.

# The loaded-extension check distinguishes unsupported arguments from setup errors.
_needs_zarr(f) = error(_no_zarr_method(f,
    Base.get_extension(DiscreteGlobalGrids, :DiscreteGlobalGridsZarrExt) !== nothing))

function _no_zarr_method(f, loaded::Bool)
    loaded && return """
        no `$f` method matches these arguments. The Zarr extension is loaded; \
        check the argument types. A store is a `Zarr.ZGroup`, a \
        `Zarr.AbstractStore`, a path or a URL. `dggwrite` accepts a `DimArray` \
        or `DimStack` over a cell dimension."""
    return """
        `$f` requires Zarr.jl. Run

            using Zarr

        to load `DiscreteGlobalGridsZarrExt`, which provides the store methods."""
end

"""
    dggread(store; vars = All(), lazy = true, validate = :strict,
            conventions = CONVENTION_REGISTRY, description = nothing) -> DimStack
    dggread(store, var::Symbol; kwargs...) -> DimArray

**Requires `using Zarr`.** Loading Zarr activates the methods and their full
keyword reference from `DiscreteGlobalGridsZarrExt`.

Read a DGGS store into plain DimensionalData with one `Cells` dimension shared
by every layer. The dimension carries one of two lookups:

  - [`ChunkedCellLookup`](@ref) resolves cells on a stored single-level axis
    without scanning it.
  - [`MultiOrderLookup`](@ref) represents the mixed-level cells of a
    `compacted` store.

The lookup type identifies the grid system, and a single-level lookup's grid
stores its level. The [`StoreDescription`](@ref) in the stack's
`metadata["description"]` stores the orientation, ellipsoid and cell-axis
layout.

`store` accepts a `Zarr.ZGroup`, `Zarr.AbstractStore`, local path or URL
(`gs://`, `s3://`, `https://`). The main keywords are:

  - `vars = All()` reads every data variable; a collection of `Symbol`s selects
    specific variables.
  - `lazy = true` preserves store-backed arrays; `false` materializes them.
  - `validate = :strict` checks every id in scanned coordinates and trusts a
    package-written chunk manifest.
  - `validate = :lazy` samples each scanned coordinate chunk.
  - `validate = :scan` ignores a trusted manifest and checks every id.
  - `description = StoreDescription(...)` supplies grid, level, encoding and
    array names directly and skips convention detection.

The stack metadata records the detected conventions, source encoding and
original attributes needed for a value-identical rewrite. A trusted manifest
lets the default validation open large package-written stores without scanning
their coordinate ids. A supplied [`StoreDescription`](@ref) enables reading an
attribute-free store while retaining the mechanical checks.
"""
dggread(args...; kwargs...) = _needs_zarr("dggread")

"""
    dggwrite(dest, stack_or_array; encoding = :auto,
             conventions = DEFAULT_WRITE_CONVENTIONS, chunks = :auto,
             merge = :step, chunk_target = 1_000_000) -> dest

**Requires `using Zarr`.** Loading Zarr activates the methods and their full
keyword reference from `DiscreteGlobalGridsZarrExt`.

Write a `DimStack` or `DimArray` over a `Cells` dimension to a Zarr v2 directory
store. A `CellLookup` or [`ChunkedCellLookup`](@ref) identifies a sorted,
unique, single-level axis; a [`MultiOrderLookup`](@ref) identifies a mixed-level
axis. `dest` is a local directory path or a writable `Zarr.ZGroup`. URL
destinations are rejected; write locally and upload, or pass an already-open
writable remote group.

`encoding` selects the cell-axis layout:

  - `:auto` selects compacted for a mixed-level axis, ranges for an eligible
    single-level axis, and dense otherwise.
  - `:dense` writes every id for broad reader compatibility.
  - `:ranges` writes the compact single-level range representation.
  - `:implicit` writes a complete level with no cell coordinate.
  - `:compacted` writes the aligned id and level columns of a mixed-level axis.

Single-level encodings require [`expand`](@ref) to present mixed-level data at
one level. `conventions` stamps the store with both default conventions for
single-level layouts. Compacted stores carry the compatible DGGS convention
metadata; xdggs attributes describe only a single-level coordinate.

The remaining layout controls are:

  - `merge = :step` merges unit-increment ids for compatibility with structural
    readers; `:rank` merges rank-adjacent cells for fewer rows and requires a
    rank-aware reader.
  - `chunks = :auto` groups complete coarse-ancestor subtree runs near the
    `chunk_target`; an integer fixes the chunk length in cells.
  - `chunk_target` counts all elements in a chunk, including the extents of
    non-cell dimensions.

The writer restores each layer's metadata as array attributes and
`metadata["attrs"]` as group attributes. Convention-generated keys take
precedence. A round trip sorts layers alphabetically and adds
`_ARRAY_DIMENSIONS` to each layer's metadata.

`layout` selects the store shape. `:cells` uses the one-dimensional cell axis
described above. `:subzones` uses the two-dimensional [`SubzoneLayout`](@ref)
and takes an `ancestor_level`.
"""
dggwrite(args...; kwargs...) = _needs_zarr("dggwrite")

"""
    subzonestore(dest, system, level; ancestor_level, layers, kwargs...) -> SubzoneStore
    subzonestore(dest) -> SubzoneStore

**Requires `using Zarr`.** The methods live in `DiscreteGlobalGridsZarrExt`,
whose docstring is the full keyword reference.

Create or reopen an ancestor-subzone store for incremental writing. Creation
stamps the group, arrays and attributes once. Each later [`dggwrite!`](@ref)
fills one column, represented by one chunk and one file. Tasks can therefore
write disjoint columns independently.

See [`SubzoneLayout`](@ref) for the layout itself and [`dggwrite`](@ref)'s
`layout = :subzones` for the one-shot form.
"""
subzonestore(args...; kwargs...) = _needs_zarr("subzonestore")

"""
    dggwrite!(store::SubzoneStore, ancestor, values; var = the only layer) -> store
    dggwrite!(store::SubzoneStore, cube) -> store

**Requires `using Zarr`.** The methods live in `DiscreteGlobalGridsZarrExt`.

Fill a store created by [`subzonestore`](@ref). A vector supplies one ancestor
cell's subtree in ascending cell-id order; a cube supplies every complete column
over its cell axis.

`values` is as long as that ancestor's subtree really is — `7^d` for a hexagon
and `(5*7^d + 1)/6` for a pentagon — and the rest of the column stays fill.
"""
dggwrite!(args...; kwargs...) = _needs_zarr("dggwrite!")
