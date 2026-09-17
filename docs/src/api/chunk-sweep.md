# Sweeping a cube along its chunk lines

```@meta
CurrentModule = DiscreteGlobalGrids
```

This API runs neighbourhood kernels over chunked cubes. It batches the cells
owned by each storage chunk with the halo needed by their rings, reducing the
repeated chunk decoding caused by scalar-at-a-time access.

The traversal follows the cube's existing chunk grid. Each callback receives an
ordinary in-memory cube containing the owned cells and the halo cells reached by
their rings. A foreign chunk can serve several halos and may be read again.

The API separates a **plan** from its **runner**. The plan records the read
order, chunk partitioning, halo, and estimated work before data access begins.

```julia
A    = dggread("dem.zarr")[:elevation]
plan = chunkplan(A; halo = 1)          # no data is read here

slope(c, v, vs) = maximum(abs(v - u) for u in vs; init = 0.0)

out = zeros(Float64, size(A))
mapneighbors!(out, slope, A, plan)
```

## The plan

[`chunkplan`](@ref DiscreteGlobalGrids.chunkplan) reads boundaries from the data's chunk grid, including
irregular layouts such as one chunk per ancestor subtree. It finds each halo by
walking its boundary through [`halo`](@ref), so planning reads metadata and
performs CPU work without loading data chunks.

[`split`](@ref Base.split(::MapChunkPlan, ::Integer)) divides a plan into
contiguous pieces with similar chunk counts. For weighted work, unequal worker
capacities, or grouping chunks by shared inputs, use the
[partitioning API](partitioning.md). Build the plan before assigning work so
its [`region`](@ref) conversion is shared by in-process tasks.

```@docs
chunkplan
MapChunkPlan
MapChunk
ownedindices
chunkhalo
nchunks(::MapChunkPlan)
halowidth
Base.split(::MapChunkPlan, ::Integer)
```

## Running it

[`foreachchunk`](@ref) hands the callback a [`ChunkCube`](@ref) containing the
chunk's cells and halo in memory, represented as an ordinary cube over a
[`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup). Package operations work on it as they do on any one-level
cell cube; ownership metadata remains alongside the cube.

A chunk is a partial grid, so its indices are **chunk-local**. The owned cells
form a contiguous range in the block, making [`localindices`](@ref) a range;
[`ownedindices`](@ref) maps that range to the caller's cell axis.
[`axisindices`](@ref) names **every** cell of the block, halo included, in that
same axis, so `axisindices(cc)[localindices(cc)] == ownedindices(cc)`. It is
what lets a sweep over a chunk report numbers the caller can use, and it is how
a [field request](@ref "Requesting neighbour fields") is translated onto a
chunk.

```@docs
foreachchunk
ChunkCube
chunkcube
localindices
axisindices
globalindices
```

## The sweeps built on it

Halo width describes available input data. The `neighborhood` selector
chooses which cells the callback receives. Repeating a one-ring sweep
is a different operation from one sweep over a wider neighborhood.

[`mapneighbors!`](@ref DiscreteGlobalGrids.mapneighbors!) keeps that invariant itself. Its `halo` defaults to
the radius of the `neighborhood` selector, `k` for both [`Disc`](@ref)`(k)`
and [`Ring`](@ref)`(k)` and one ring for `Disc(0)`, and a supplied plan must carry at least that many
rings: [`halowidth`](@ref)`(plan)` narrower than the radius throws an
`ArgumentError` naming both widths. A wider halo than the radius is allowed.
With `A` a stored cube, `out` a destination of the same shape and `kernel` a
`(cell, value, values)` function:

```julia
plan = chunkplan(A; halo = 3)
mapneighbors!(out, kernel, A, plan; neighborhood = Ring(3))   # halowidth(plan) ≥ 3
mapneighbors!(out, kernel, A; neighborhood = Disc(3))         # plans halo = 3 itself
```

[`mapneighbors!`](@ref DiscreteGlobalGrids.mapneighbors!) is the streaming form: results are written into `dest` a
chunk at a time, so neither the input nor the output has to fit in memory.
[`mapneighbors`](@ref) with `pass = Values()` takes the same route by itself
whenever the cube's data is chunked, and collects the results.

A [field request](@ref "Requesting neighbour fields") takes it too.
`mapneighbors!(dest, f, A, plan; needs = (Value(dem), Centroid()))` states the
request once, about the cube that was passed, and the route translates it onto
each chunk: `Index(Local())` keeps answering the caller's cell-axis index,
never a chunk-local one, and each stored `Value` is read along its own chunk
grid the way the swept data is. `mapneighbors` and `foreachneighbors` reach for
it by themselves under the rule `Values()` uses — a chunked parent and
`order = StorageOrder()`.

```julia
plan = chunkplan(A; halo = 1)
out  = zeros(Float64, size(A))
mapneighbors!(out, steepest, A, plan; needs = (Value(A), Centroid()))
```

[`Neighbors`](@ref) uses a different callback contract: it supplies cell handles
and the callback reads values from the original array. Use `Values()` when a
chunked sweep should stream the fields through the traversal.

The [neighbourhood API](neighbors.md#compute-with-neighbourhoods) documents
these kernels and callback forms.

## Index

```@index
Pages = ["api/chunk-sweep.md"]
```
