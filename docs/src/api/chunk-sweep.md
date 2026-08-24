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
globalindices
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

The owned cells are contiguous within the block: every halo cell is by
definition outside the chunk's own run, so sorting the two into one axis leaves
the owned run unbroken. [`localindices`](@ref) is therefore a range, and
[`globalindices`](@ref) says where those results belong in the full axis.

```@docs
foreachchunk
ChunkCube
chunkcube
localindices
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
whenever the cube's data is chunked, and collects the results as it always has.

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
