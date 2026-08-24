# Sweeping a cube chunk by chunk.
#
# A cell-at-a-time pass over a lazy cube decodes one storage chunk PER SCALAR
# READ: a cell, then each of its ring, then the next cell. The ring of a cell
# near a chunk boundary lands in a neighbouring chunk, so the same chunk is
# decoded, dropped and decoded again, cell after cell.
#
# The fix is to make the traversal follow the chunk lines the store already has.
# For each chunk: read its cells once, read the cells its rings reach outside it
# once per foreign chunk, and hand the pair to the caller as an ordinary
# in-memory cube. Everything in this package then works on it unmodified,
# because it IS an ordinary cube — a `CellLookup` over cells at one level.
#
# The plan is separate from the run for the reason the user of one asks for:
# what to read, in what order, and how much of it is a value you can look at,
# reorder, and cut into pieces that run on different tasks.

import DimensionalData as DD
import DiskArrays

# `nchunks` already counts the chunks of a `ChunkManifest`; a plan is the other
# thing with chunks to count, so it answers the same verb.
import .ChunkedLookups: nchunks

# ===========================================================================
# The plan
# ===========================================================================

"""
    MapChunk

One chunk of a [`MapChunkPlan`](@ref): the axis positions it owns, and the axis
positions outside it that its cells' rings reach.

`chunkrange(mc)` is the first, `chunkhalo(mc)` the second; both are positions in
the cube's own cell axis, and both ascend. The halo holds only positions the
axis really has — a ring that leaves the axis altogether is clipped here exactly
as [`neighbors`](@ref) clips it.
"""
struct MapChunk
    index::Int
    range::UnitRange{Int}
    halo::Vector{Int}
end

"""
    chunkrange(mc::MapChunk) -> UnitRange{Int}
    chunkrange(cc::ChunkCube) -> UnitRange{Int}

The positions in the CUBE'S CELL AXIS that this chunk owns. Contiguous, because
a chunk is a contiguous block of the axis.
"""
chunkrange(mc::MapChunk) = mc.range

"""
    chunkhalo(mc::MapChunk) -> Vector{Int}

The axis positions outside [`chunkrange`](@ref) that this chunk's cells reach,
ascending. These are the cells a sweep over the chunk has to read but does not
produce a result for.
"""
chunkhalo(mc::MapChunk) = mc.halo

Base.length(mc::MapChunk) = length(mc.range)

Base.show(io::IO, mc::MapChunk) = print(io, "MapChunk(", mc.index, ", ",
    length(mc.range), " cells, ", length(mc.halo), " halo)")

"""
    MapChunkPlan

What a chunk-following sweep will read, in what order — a value, not a running
traversal.

Build one with [`chunkplan`](@ref) and run it with [`foreachchunk`](@ref). It is
worth being a value because three separate decisions are made on it:

  - **order**: `plan[i]` is the `i`th chunk to visit. Reorder by rebuilding the
    plan from a permutation of its chunks.
  - **parallelism**: `split(plan, n)` cuts it into `n` plans over disjoint
    chunks. Running those on separate tasks is what parallelises a sweep;
    there is no `threaded` keyword here because the split IS the decision.
  - **cost**: `nchunks(plan)` and the per-chunk `chunkhalo` lengths say how much
    will be read before anything is.

The halo is computed when the plan is built, which costs one region walk per
chunk and no IO. Chunk boundaries come from the data's own chunk grid, so a
store whose chunks are irregular — one per ancestor subtree, say — is planned on
its real boundaries and not on an assumed uniform length.
"""
struct MapChunkPlan{L<:AbstractCellLookup,K<:Connectivity}
    lookup::L
    chunks::Vector{MapChunk}
    width::Int
    connectivity::K
end

"""
    nchunks(plan::MapChunkPlan) -> Int

How many chunks the plan visits.
"""
ChunkedLookups.nchunks(plan::MapChunkPlan) = length(plan.chunks)

Base.length(plan::MapChunkPlan) = length(plan.chunks)
Base.getindex(plan::MapChunkPlan, i::Integer) = plan.chunks[i]
Base.iterate(plan::MapChunkPlan, s=1) =
    s > length(plan.chunks) ? nothing : (plan.chunks[s], s + 1)
Base.eltype(::Type{<:MapChunkPlan}) = MapChunk

"""
    halowidth(plan::MapChunkPlan) -> Int

How many rings of context each chunk's halo carries. `1` is what a one-ring
stencil needs and is the default.
"""
halowidth(plan::MapChunkPlan) = plan.width

"""
    split(plan::MapChunkPlan, n::Integer) -> Vector{MapChunkPlan}

Cut the plan into at most `n` plans over disjoint chunks, in order.

This is how a chunked sweep is parallelised: run the pieces on separate tasks.
Each chunk writes only the positions it owns, so pieces running at once cannot
collide however the cut falls — and the result does not depend on how many
pieces there were.
"""
function Base.split(plan::MapChunkPlan, n::Integer)
    nc = length(plan.chunks)
    k = clamp(Int(n), 1, max(nc, 1))
    nc == 0 && return [plan]
    return [MapChunkPlan(plan.lookup,
        plan.chunks[(1+div((t - 1) * nc, k)):div(t * nc, k)],
        plan.width, plan.connectivity) for t in 1:k]
end

function Base.show(io::IO, plan::MapChunkPlan)
    h = isempty(plan.chunks) ? 0 : sum(length(chunkhalo(c)) for c in plan.chunks)
    print(io, "MapChunkPlan(", nchunks(plan), " chunks, ",
        length(plan.lookup), " cells, halo=", plan.width, ", ", h, " halo reads)")
end

Base.show(io::IO, ::MIME"text/plain", plan::MapChunkPlan) = show(io, plan)

"""
    chunkplan(A::AbstractDimArray; halo = 1, chunks = :auto, spatialdim = nothing,
              connectivity = Vertex()) -> MapChunkPlan
    chunkplan(lk::AbstractCellLookup, bounds; halo = 1, connectivity = Vertex())

Plan a chunk-following sweep of `A`: which cells each chunk owns, and which
cells outside it its rings reach.

`chunks = :auto` takes the boundaries from the data's own chunk grid
(`DiskArrays.eachchunk`), which is the whole point — reading along them is what
makes each chunk decode once. An `Integer` fixes a uniform chunk length in cells
instead, and an in-memory array, having no chunk grid, is one chunk.

`halo = n` carries `n` rings of context, which is what a stencil reaching `n`
rings needs; `1` is the one-ring default. `spatialdim` and `connectivity` are
[`mapneighbors`](@ref)'s.

The second form plans against a cell axis alone, with `bounds` either a chunk
length or the chunk ranges themselves — the form to reach for when the cube is
not the thing being read.

Building the plan walks each chunk's boundary once. That is CPU, not IO: no
chunk of the data is touched until [`foreachchunk`](@ref) runs the plan.

```julia
plan = chunkplan(A; halo = 1)
nchunks(plan)                             # what it will read
foreachchunk(A, plan) do cc               # ... and reading it
    r = mapneighbors(f, chunkcube(cc))
    out[chunkrange(cc)] = parent(r)[chunkpositions(cc)]
end
```
"""
function chunkplan(A::DD.AbstractDimArray; halo::Integer=1, chunks=:auto,
        spatialdim=nothing, connectivity::Connectivity=Vertex())
    dnum = CellLookups._cells_dimnum(A, spatialdim)
    lk = DD.lookup(A, dnum)
    return chunkplan(lk, _chunkbounds(parent(A), dnum, chunks, length(lk));
        halo=halo, connectivity=connectivity)
end

function chunkplan(lk::AbstractCellLookup, bounds; halo::Integer=1,
        connectivity::Connectivity=Vertex())
    w = Int(halo)
    w >= 1 || throw(ArgumentError(
        "a chunk's halo carries at least one ring, not $w; a sweep with no " *
        "context is `mapneighbors` on each chunk separately"))
    cv = region(lk)
    ranges = _asranges(bounds, length(cv))
    chunks = Vector{MapChunk}(undef, length(ranges))
    for (i, r) in enumerate(ranges)
        chunks[i] = MapChunk(i, r, _chunkhalo(cv, r, w, connectivity))
    end
    return MapChunkPlan(lk, chunks, w, connectivity)
end

# The chunk grid, from whichever of the three sources was named.
_chunkbounds(data, dnum::Int, chunks::Integer, n::Int) = Int(chunks)

function _chunkbounds(data, dnum::Int, chunks::Symbol, n::Int)
    chunks === :auto || throw(ArgumentError(
        "chunks is :auto or a positive chunk length in cells, not $(repr(chunks))"))
    DiskArrays.haschunks(data) isa DiskArrays.Chunked || return n
    grid = DiskArrays.eachchunk(data)
    return [UnitRange{Int}(r) for r in grid.chunks[dnum]]
end

_asranges(bounds::AbstractVector{UnitRange{Int}}, n::Int) = bounds

function _asranges(chunklength::Integer, n::Int)
    cl = Int(chunklength)
    cl >= 1 || throw(ArgumentError(
        "a chunk length is at least one cell, not $cl"))
    n == 0 && return UnitRange{Int}[]
    return [lo:min(lo + cl - 1, n) for lo in 1:cl:n]
end

# The axis positions a chunk's cells reach but does not own.
#
# Asked of the REGION rather than cell by cell: `halo` walks a boundary in
# `O(border)`, where enumerating every cell's ring and discarding the ones that
# land back inside is `O(cells x degree)`. The answer ascends without sorting,
# because a complete-level position and an axis position both ascend in cell id,
# so the map between them is monotone.
function _chunkhalo(cv::CellVector, r::UnitRange{Int}, width::Int,
        conn::Connectivity)
    out = Int[]
    isempty(r) && return out
    sub = cv[r]
    if width == 1
        for p in halo(sub; connectivity=conn)
            q = Engine.windowposition(cv.windows, p)
            q === nothing || push!(out, q)
        end
    else
        # Wider than one ring has no single boundary walk, so the region is
        # grown and what it gained outside the chunk is the halo.
        grown = grow(sub, width; connectivity=conn)
        for p in Engine.leafpositions(grown.windows)
            q = Engine.windowposition(cv.windows, p)
            (q === nothing || first(r) <= q <= last(r)) && continue
            push!(out, q)
        end
    end
    return out
end

# ===========================================================================
# One loaded chunk
# ===========================================================================

"""
    ChunkCube

One chunk of a cube, read: its own cells and its halo, in memory, as an ordinary
`DimArray` over a [`CellLookup`](@ref).

This is what [`foreachchunk`](@ref) hands its callback, and the reason it is
worth handing over rather than hiding: every verb in this package already works
on a cube, so a chunk-following pass is written in the same vocabulary as the
whole-cube pass it replaces.

| accessor | what it gives |
|---|---|
| [`chunkcube`](@ref) | the in-memory cube: the chunk's cells AND its halo |
| [`chunkpositions`](@ref) | the positions IN THAT CUBE the chunk owns |
| [`chunkrange`](@ref) | the same cells' positions in the FULL axis |

A result computed on the cube is meaningful for the owned positions only:
a halo cell is there to complete its neighbours' rings, and its own ring is
missing whatever fell outside the block. `r[chunkpositions(cc)]` selects the
part to keep and `chunkrange(cc)` says where it belongs.

Because the halo holds every axis neighbour of every owned cell, a stencil that
reaches no further than the plan's [`halowidth`](@ref) computes exactly what the
whole-axis sweep computes for those cells.
"""
struct ChunkCube{C<:DD.AbstractDimArray}
    cube::C
    positions::UnitRange{Int}
    range::UnitRange{Int}
    index::Int
end

"""
    chunkcube(cc::ChunkCube) -> AbstractDimArray

The chunk's cells and its halo, in memory, over a [`CellLookup`](@ref).
"""
chunkcube(cc::ChunkCube) = cc.cube

"""
    chunkpositions(cc::ChunkCube) -> UnitRange{Int}

The positions in [`chunkcube`](@ref) that the chunk OWNS — the ones a result is
kept for. Contiguous: a chunk is a run of the axis, and every halo cell is
outside that run, so the owned cells stay together when the two are sorted into
one axis.
"""
chunkpositions(cc::ChunkCube) = cc.positions

chunkrange(cc::ChunkCube) = cc.range

"""
    chunkhalo(cc::ChunkCube) -> Vector{Int}

The positions in [`chunkcube`](@ref) that are HALO — read to complete the owned
cells' rings, and not results to keep. The complement of
[`chunkpositions`](@ref) over the cube, so the two together are all of it.

These are positions in the block, not in the axis; `chunkhalo` of the plan's
[`MapChunk`](@ref) gives the axis positions they came from.
"""
chunkhalo(cc::ChunkCube) = vcat(1:(first(cc.positions)-1),
    (last(cc.positions)+1):_cellcount(cc))

# The block cube's cell dimension, and how long it is. Read off the block rather
# than carried, so the accessors work on a `ChunkCube` alone; on an N-D cube the
# cell count is that dimension's length and not the cube's element count.
_celldim(cc::ChunkCube) = CellLookups._cells_dimnum(cc.cube, nothing)
_cellcount(cc::ChunkCube) = size(cc.cube, _celldim(cc))

Base.parent(cc::ChunkCube) = cc.cube

Base.show(io::IO, cc::ChunkCube) = print(io, "ChunkCube(", cc.index, ", ",
    length(cc.range), " cells + ", _cellcount(cc) - length(cc.range), " halo)")

# ===========================================================================
# Running a plan
# ===========================================================================

"""
    foreachchunk(f, A::AbstractDimArray, plan = chunkplan(A; kw...); kw...)

Run `plan` over `A`, calling `f(cc::ChunkCube)` once per chunk.

Each call gets the chunk's cells and its halo as an in-memory cube — one
contiguous read for the chunk itself, and one read per foreign chunk the halo
touches, never more. `f`'s return value is discarded; this is the primitive the
result-producing forms are built on, and the one to reach for when what a pass
produces is not one value per cell.

Chunks are visited in the plan's order, on the calling task. To use more than
one, `split` the plan and run the pieces — see [`MapChunkPlan`](@ref).

```julia
plan = chunkplan(A; halo = 1)
@sync for p in Base.split(plan, Threads.nthreads())
    Threads.@spawn foreachchunk(A, p) do cc
        write!(out, chunkrange(cc), mine(mapneighbors(f, chunkcube(cc)), cc))
    end
end
```
"""
function foreachchunk(f::F, A::DD.AbstractDimArray, plan::MapChunkPlan;
        spatialdim=nothing) where {F}
    dnum = CellLookups._cells_dimnum(A, spatialdim)
    cv = region(DD.lookup(A, dnum))
    # The storage chunk grid is the array's and does not change between chunks,
    # so it is read once here rather than per gather.
    bounds = _storagebounds(parent(A), dnum)
    for mc in plan.chunks
        f(_loadchunk(A, dnum, cv, mc, bounds))
    end
    return nothing
end

foreachchunk(f::F, A::DD.AbstractDimArray; spatialdim=nothing, kw...) where {F} =
    foreachchunk(f, A, chunkplan(A; spatialdim=spatialdim, kw...);
        spatialdim=spatialdim)

# Read one chunk and its halo into an ordinary cube.
#
# The three pieces are read separately and joined rather than gathered in one
# scattered index, so the chunk's own cells stay a single contiguous read and
# each foreign chunk the halo reaches is read exactly once.
function _loadchunk(A::DD.AbstractDimArray, dnum::Int, cv::CellVector,
        mc::MapChunk, bounds::Vector{UnitRange{Int}})
    r = mc.range
    below = _lessthan(mc.halo, first(r))
    above = @view mc.halo[(length(below)+1):end]
    data = parent(A)
    parts = Any[]
    isempty(below) || push!(parts, _readpositions(data, dnum, below, bounds))
    push!(parts, _readrange(data, dnum, r))
    isempty(above) || push!(parts, _readpositions(data, dnum, above, bounds))
    block = length(parts) == 1 ? only(parts) : cat(parts...; dims=dnum)
    positions = vcat(collect(below), collect(r), collect(above))
    lk = CellLookup(cv[positions])
    dims = ntuple(i -> i == dnum ? Cells(lk) : DD.dims(A)[i], Val(ndims(A)))
    own = (length(below)+1):(length(below)+length(r))
    return ChunkCube(DD.rebuild(A; data=block, dims=dims), own, r, mc.index)
end

_lessthan(halo::Vector{Int}, lo::Int) =
    @view halo[1:(searchsortedfirst(halo, lo)-1)]

function _readrange(data, dnum::Int, r::AbstractUnitRange)
    inds = ntuple(i -> i == dnum ? r : Colon(), Val(ndims(data)))
    return Array(data[inds...])
end

# Scattered positions, read one covering range per storage chunk they land in.
# Within a chunk the covering range is decoded as a unit anyway, so reading it
# whole costs nothing over reading part of it — and reading the runs separately
# is what stops a halo that straddles a distant chunk from pulling in everything
# between.
function _readpositions(data, dnum::Int, idx::AbstractVector{Int},
        bounds::Vector{UnitRange{Int}})
    parts = Any[]
    for run in _positionruns(idx, bounds)
        block = _readrange(data, dnum, run)
        take = [p - first(run) + 1 for p in idx if first(run) <= p <= last(run)]
        inds = ntuple(i -> i == dnum ? take : Colon(), Val(ndims(block)))
        push!(parts, block[inds...])
    end
    return length(parts) == 1 ? only(parts) : cat(parts...; dims=dnum)
end

# Group ascending positions by the storage chunk they fall in, so each chunk a
# gather reaches is read once and no two reads overlap.
function _positionruns(idx::AbstractVector{Int}, bounds::Vector{UnitRange{Int}})
    runs = UnitRange{Int}[]
    i = firstindex(idx)
    while i <= lastindex(idx)
        stop = last(_boundsof(bounds, idx[i]))
        j = i
        while j < lastindex(idx) && idx[j+1] <= stop
            j += 1
        end
        push!(runs, idx[i]:idx[j])
        i = j + 1
    end
    return runs
end

_storagebounds(data, dnum::Int) =
    DiskArrays.haschunks(data) isa DiskArrays.Chunked ?
    [UnitRange{Int}(r) for r in DiskArrays.eachchunk(data).chunks[dnum]] :
    [1:size(data, dnum)]

# Binary search, not a scan: the chunk grid of a store of tens of millions of
# cells has thousands of rows and this is asked once per gathered run.
function _boundsof(bounds::Vector{UnitRange{Int}}, p::Int)
    lo, hi = firstindex(bounds), lastindex(bounds)
    while lo < hi
        mid = (lo + hi) >>> 1
        p > last(bounds[mid]) ? (lo = mid + 1) : (hi = mid)
    end
    return bounds[lo]
end

# ===========================================================================
# The result-producing forms
# ===========================================================================

# --- the automatic route ---------------------------------------------------

# `mapneighbors(f, A; pass = Values())` on a cube whose data is chunked runs the
# plan. Only in storage order: a permutation names a visit order over the whole
# axis, and a chunked sweep visits by chunk, so the two cannot both be honoured.
CellLookups._chunkedvalues(A::DD.AbstractDimArray, dnum::Int, ::StorageOrder,
    conn::Connectivity) =
    DiskArrays.haschunks(parent(A)) isa DiskArrays.Chunked ?
    chunkplan(DD.lookup(A, dnum),
        _chunkbounds(parent(A), dnum, :auto, size(A, dnum));
        halo=1, connectivity=conn) : nothing

# The results are collected in memory, as `mapneighbors` promises: what the plan
# saves here is the READ side. `mapneighbors!` is the form that streams the
# write side too.
function CellLookups._map_values_chunked(f::F, A, dnum::Int,
        plan::MapChunkPlan, threaded, conn::Connectivity) where {F}
    cv = region(DD.lookup(A, dnum))
    cap = Engine._capacity(system(cv), conn)
    H = Engine.SubsetPositionedCell{eltype(cv)}
    T = Base.promote_op(f, H, eltype(A), Engine._ringtype(cap, eltype(A)))
    sz = size(A)
    dest = T <: Tuple && isconcretetype(T) ?
           ntuple(j -> Array{fieldtype(T, j)}(undef, sz), fieldcount(T)) :
           Array{T}(undef, sz)
    mapneighbors!(dest, f, A, plan; threaded=threaded)
    return dest
end

"""
    mapneighbors!(dest, f, A::AbstractDimArray; halo = 1, chunks = :auto,
                  spatialdim = nothing, connectivity = Vertex(), threaded = true)

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

The result is the whole-axis sweep's, cell for cell. Each chunk's halo carries
every axis neighbour of every cell the chunk owns, so a ring computed on a chunk
is the ring computed on the axis — clipped identically, and in the same order.

`threaded` threads WITHIN each chunk. To run chunks themselves in parallel,
`split` a [`chunkplan`](@ref) and call this on the pieces. Build the plan before
splitting: doing so is what fills the axis's [`region`](@ref) memo, so the
pieces share one conversion rather than each repeating it.
"""
function mapneighbors!(dest, f::F, A::DD.AbstractDimArray; halo::Integer=1,
        chunks=:auto, spatialdim=nothing, connectivity::Connectivity=Vertex(),
        threaded=true) where {F}
    plan = chunkplan(A; halo=halo, chunks=chunks, spatialdim=spatialdim,
        connectivity=connectivity)
    return mapneighbors!(dest, f, A, plan; spatialdim=spatialdim,
        threaded=threaded)
end

function mapneighbors!(dest, f::F, A::DD.AbstractDimArray, plan::MapChunkPlan;
        spatialdim=nothing, threaded=true) where {F}
    foreachchunk(A, plan; spatialdim=spatialdim) do cc
        cube = chunkcube(cc)
        out = mapneighbors(f, cube; pass=Values(), threaded=threaded,
            connectivity=plan.connectivity)
        _storechunk!(dest, out, cc)
    end
    return dest
end

# Keep the owned rows and put them where the chunk came from. A tuple result is
# one output per component, so `dest` is then a tuple of the same length.
function _storechunk!(dest, out, cc::ChunkCube)
    dnum = CellLookups._cells_dimnum(out, nothing)
    src = _slice(parent(out), dnum, chunkpositions(cc))
    _destview(dest, dnum, chunkrange(cc)) .= src
    return nothing
end

function _storechunk!(dest::Tuple, out::Tuple, cc::ChunkCube)
    map((d, o) -> _storechunk!(d, o, cc), dest, out)
    return nothing
end

_slice(data, dnum::Int, idx) =
    data[ntuple(i -> i == dnum ? idx : Colon(), Val(ndims(data)))...]

_destview(dest::DD.AbstractDimArray, dnum::Int, r) = _destview(parent(dest), dnum, r)
_destview(dest, dnum::Int, r) =
    view(dest, ntuple(i -> i == dnum ? r : Colon(), Val(ndims(dest)))...)
