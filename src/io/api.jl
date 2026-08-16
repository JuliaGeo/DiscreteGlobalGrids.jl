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
             conventions = DEFAULT_WRITE_CONVENTIONS, chunks = :auto) -> dest

Write a `DimStack` or `DimArray` over a `Cells` dimension to a DGGS store.

`encoding = :auto` picks ranges where the axis is eligible — sorted, unique,
single-level — and dense otherwise; `encoding = :dense` is the interop escape
for readers that cannot expand ranges. `conventions` stamps the store, dual by
default so that both a convention-aware reader and xdggs can open it.

Requires `using Zarr`.
"""
dggwrite(args...; kwargs...) = _needs_zarr("dggwrite")
