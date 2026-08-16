# `dggread`/`dggwrite` are defined and exported here as stubs. Their methods
# live in `DiscreteGlobalGridsZarrExt`, so the store types, keyword defaults
# and IO all stay behind the Zarr weak dependency.

_needs_zarr(f) = error("""
    `$f` requires Zarr.jl. Run

        using Zarr

    to load `DiscreteGlobalGridsZarrExt`, which provides the store methods.""")

"""
    dggread(store; vars, lazy = true, validate = :strict,
            conventions = CONVENTION_REGISTRY, description = nothing) -> DimStack
    dggread(store, var::Symbol; kwargs...) -> DimArray

Read a DGGS store into plain DimensionalData: one `Cells` dimension shared by
every layer, carrying a `CellLookup` that holds grid, level and orientation in
its types.

`store` is a `Zarr.ZGroup`, a local path, or a URL (`gs://`, `s3://`,
`https://`). Data arrays are lazy by default; `lazy = false` materializes them.
The detected convention, the verbatim original attributes and the source
encoding ride in the stack's `metadata`, which is enough to regenerate a
value-identical store.

`description` bypasses detection: pass a [`StoreDescription`](@ref) and the
caller asserts grid, level, encoding and array names, leaving only the
mechanical checks. That is how an attribute-less store is read.

Requires `using Zarr`.
"""
dggread(args...; kwargs...) = _needs_zarr("dggread")

"""
    dggwrite(dest, stack_or_array; encoding = :auto,
             conventions = DEFAULT_WRITE_CONVENTIONS, chunks = :auto,
             merge = :rank, chunk_target = 1_000_000) -> dest

Write a `DimStack` or `DimArray` over a `Cells` dimension to a DGGS store.
`dest` is a local directory path or a writeable `Zarr.ZGroup`; remote stores
are not written in v1 — write locally and upload.

`encoding = :auto` picks ranges where the axis is eligible — sorted, unique,
single-level — and dense otherwise; `encoding = :dense` is the interop escape
for readers that cannot expand ranges. `conventions` stamps the store, dual by
default so that both a convention-aware reader and xdggs can open it.

`merge` picks the ranges run rule: `:rank` (default) merges rank-adjacent cells
for the fewest rows; `:step` merges unit-increment ids, which a structural
reader also counts correctly. `chunks = :auto` aims each chunk at
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

Requires `using Zarr`.
"""
dggwrite(args...; kwargs...) = _needs_zarr("dggwrite")
