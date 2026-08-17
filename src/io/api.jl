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
    dggread(store; vars, lazy = true, validate = :strict,
            conventions = CONVENTION_REGISTRY, description = nothing) -> DimStack
    dggread(store, var::Symbol; kwargs...) -> DimArray

Read a DGGS store into plain DimensionalData: one `Cells` dimension shared by
every layer, carrying a [`ChunkedCellLookup`](@ref) — the lookup over an axis a
store wrote, which resolves a cell without scanning it — or, for a `compacted`
store, a [`MultiOrderLookup`](@ref) over its mixed-level cells. The grid SYSTEM
is in that lookup's type, the level is a field of the grid it holds, and what is
neither — orientation, ellipsoid, the layout the store keeps its axis in — rides
in the [`StoreDescription`](@ref) under the stack's `metadata["description"]`.

`store` is a `Zarr.ZGroup`, a local path, or a URL (`gs://`, `s3://`,
`https://`). Data arrays are lazy by default; `lazy = false` materializes them.
The detected convention, the verbatim original attributes and the source
encoding ride in the stack's `metadata`, which is enough to regenerate a
value-identical store.

`validate = :strict` checks that every stored id names a cell of the declared
level and `:lazy` samples instead. Neither reaches a store carrying a chunk
manifest this package wrote: that axis is built from the manifest and no id is
scanned, which is what opens a store of tens of millions of cells at all.
`validate = :scan` declines the manifest and runs the full scan on any store.

`description` bypasses detection: pass a [`StoreDescription`](@ref) and the
caller asserts grid, level, encoding and array names, leaving only the
mechanical checks. That is how an attribute-less store is read.

Requires `using Zarr`.
"""
dggread(args...; kwargs...) = _needs_zarr("dggread")

"""
    dggwrite(dest, stack_or_array; encoding = :auto,
             conventions = DEFAULT_WRITE_CONVENTIONS, chunks = :auto,
             merge = :step, chunk_target = 1_000_000) -> dest

Write a `DimStack` or `DimArray` over a `Cells` dimension to a DGGS store.
`dest` is a local directory path or a writeable `Zarr.ZGroup`; remote stores
are not written in v1 — write locally and upload.

`encoding = :auto` picks compacted for a mixed-level axis (a
[`MultiOrderLookup`](@ref)), ranges where a single-level axis is eligible —
sorted, unique — and dense otherwise; `encoding = :dense` is the interop
escape for readers that cannot expand ranges. A single-level encoding
requested for a mixed-level axis is refused; present the cube at one level
with [`expand`](@ref) to write it that way. `conventions` stamps the store,
dual by default so that both a convention-aware reader and xdggs can open it;
a compacted store carries no xdggs stamp, which cannot say mixed-level.

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

Requires `using Zarr`.
"""
dggwrite(args...; kwargs...) = _needs_zarr("dggwrite")
