# Sweeping a cube along its chunk lines

```@meta
CurrentModule = DiscreteGlobalGrids
```

A cell-at-a-time pass over a lazy cube decodes one storage chunk **per scalar
read**. A cell near a chunk boundary has ring members in the next chunk over, so
that chunk is decoded, dropped, and decoded again for the cell after it — and
again for the one after that. The cost is not the arithmetic; it is the same
compressed bytes being inflated hundreds of times.

The traversal here follows the chunk grid the store already has. For each chunk
it reads the cells the chunk owns once, reads the cells their rings reach
outside it once per *distinct* foreign chunk, and hands the pair to the caller
as an ordinary in-memory cube.

It is spelled as a **plan** and a **runner**, because three decisions are made
on a plan and none of them can be made inside a single opaque call: what order
to read in, how to cut the work up, and what it will cost before anything is
read.

```julia
A    = dggread("dem.zarr")[:elevation]
plan = chunkplan(A; halo = 1)          # no data is read here

slope(c, v, vs) = maximum(abs(v - u) for u in vs; init = 0.0)

out = zeros(Float64, size(A))
mapneighbors!(out, slope, A, plan)
```

## The plan

[`chunkplan`](@ref) reads the chunk boundaries from the data's own chunk grid,
so a store chunked irregularly — one chunk per ancestor subtree, say — is
planned on its real boundaries rather than on a nominal chunk length. Each
chunk's halo is found by walking its boundary through [`halo`](@ref), which is
CPU and no IO: a plan over a store of tens of millions of cells is built without
touching a single data chunk.

[`split`](@ref Base.split(::MapChunkPlan, ::Integer)) is how a sweep is
parallelised, which is why the runner has no `threaded` keyword for the chunks
themselves. Pieces of one plan own disjoint ranges of the axis, so tasks running
at once write disjoint ranges of the destination and need no coordination —
including when the destination is a store. Build the plan *before* splitting:
that call fills the axis's [`region`](@ref) memo, so the pieces share one
conversion instead of each repeating it.

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

[`foreachchunk`](@ref) hands the callback a [`ChunkCube`](@ref): the chunk's
cells **and its halo**, in memory, as an ordinary cube over a
[`CellLookup`](@ref). Every verb in this package then works on it unmodified,
because there is nothing special about it — it is a cube over cells at one
level, and the fact that some of those cells are context rather than results is
carried beside it, not inside it.

A chunk is a partial grid, so an index into that cube is **chunk-local** and
means nothing outside it. Two accessors carry the translation. The owned cells
are contiguous within the block — every halo cell is by definition outside the
chunk's own run, so sorting the two into one axis leaves the owned run unbroken
— which makes [`localindices`](@ref) a range, and [`ownedindices`](@ref) says
where those results belong in the caller's cell axis.
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

Because a chunk's halo carries every axis neighbour of every cell the chunk
owns, a stencil reaching no further than the plan's halo width computes on a
chunk exactly what it computes on the whole axis — clipped identically, in the
same order. That is what lets the two sweep forms be built on the plan without
qualifying their results.

[`mapneighbors!`](@ref) is the streaming form: results are written into `dest` a
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

[`Neighbors`](@ref) deliberately does **not** take this route. Its callback is
handed cell handles and reaches back into the original array for values, so a
sweep over blocks would hand it block indices to index the whole cube with.
A chunked sweep is exactly the case where the values must flow through the
traversal, which is what `Values()` means.

```@docs
mapneighbors!
mapneighbors
foreachneighbors
Values
Neighbors
```

## Index

```@index
Pages = ["api/chunk-sweep.md"]
```
