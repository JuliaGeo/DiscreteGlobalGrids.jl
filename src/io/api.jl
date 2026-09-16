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

Read a DGGS store as a dimensional cube. Requires `using Zarr`.
The stack shares a [`Cells`](@ref DiscreteGlobalGrids.CellLookups.Cells) axis backed by [`ChunkedCellLookup`](@ref DiscreteGlobalGrids.ChunkedLookups.ChunkedCellLookup).
Data stay lazy unless `lazy=false`; the single-variable form returns a `DimArray`.

`store` accepts a local path, URL, `Zarr.ZGroup`, or `Zarr.AbstractStore`.
Public `gs://` URLs use HTTPS; `s3://` additionally requires `using AWSS3`.
`vars` selects data variables. Metadata retain source attributes and the detected
grid description for a later rewrite.

`validate=:strict` validates IDs when scanning an axis. A stored chunk manifest
can avoid that scan; `:scan` forces it. `:lazy` samples IDs instead.
A supplied [`StoreDescription`](@ref) bypasses metadata detection, not mechanical
validation. Invalid formats raise [`DGGSFormatError`](@ref).

See [Reading and writing DGGS stores](@ref) for formats and examples, and
[Workflow execution details](@ref) for validation and metadata rules.
"""
dggread(args...; kwargs...) = _needs_zarr("dggread")

"""
    dggwrite(dest, stack_or_array; encoding = :auto,
             conventions = DEFAULT_WRITE_CONVENTIONS, chunks = :auto,
             merge = :step, chunk_target = 1_000_000) -> dest

Write a `DimArray` or `DimStack` with a cell lookup to a Zarr v2 store.
Requires `using Zarr`. `dest` is a local directory or writable `Zarr.ZGroup`;
remote URL writing is not supported.

`encoding=:auto` chooses an eligible ranges encoding or dense IDs. `:dense`
stores each ID; `:ranges` stores intervals; `:implicit` requires a complete level.
`merge=:step` joins integer-adjacent IDs. `merge=:rank` joins consecutive valid
cells and requires a rank-aware reader.

`chunks` is a cell chunk length or `:auto`. `chunk_target` counts all elements
per chunk, including non-cell dimensions. Layer metadata become array attributes;
group attributes come from `metadata["attrs"]`. Generated convention keys take
precedence, and layer order is normalized alphabetically.

`layout=:subzones` selects the separate ancestor-subzone writer and requires
`ancestor_level`. See [Reading and writing DGGS stores](@ref),
[The ancestor-subzone layout](@ref), and [Workflow execution details](@ref)
for layout-specific options and metadata rules.
"""
dggwrite(args...; kwargs...) = _needs_zarr("dggwrite")

"""
    subzonestore(dest, system, level; ancestor_level, layers, kwargs...) -> SubzoneStore
    subzonestore(dest) -> SubzoneStore

**Requires `using Zarr`.** The methods live in `DiscreteGlobalGridsZarrExt`,
whose docstring is the full keyword reference.

Create — or reopen — an ancestor-subzone store for incremental writing: the
group, its arrays and its attributes are stamped once, and the columns are
filled afterwards, one [`dggwrite!`](@ref DiscreteGlobalGrids.dggwrite!) at a time. A column is one chunk and
therefore one file, and a column write rewrites nothing shared, so tasks writing
disjoint columns need no coordination.

See [`SubzoneLayout`](@ref) for the layout itself and
[`dggwrite`](@ref DiscreteGlobalGrids.dggwrite)'s
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
