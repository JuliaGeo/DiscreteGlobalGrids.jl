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

One chunk of a [`MapChunkPlan`](@ref): the axis indices it owns, and the axis
indices outside it that its cells' rings reach.

`ownedindices(mc)` is the first, `chunkhalo(mc)` the second; both are indices in
the cube's own cell axis, and both ascend. The halo holds only indices the
axis really has — a ring that leaves the axis altogether is clipped here exactly
as [`neighbors`](@ref) clips it.
"""
struct MapChunk
    index::Int
    range::UnitRange{Int}
    halo::Vector{Int}
end

"""
    ownedindices(mc::MapChunk) -> UnitRange{Int}
    ownedindices(cc::ChunkCube) -> UnitRange{Int}

The indices in the CELL AXIS OF THE CUBE THE CALLER PASSED that this chunk
owns — the cells it produces results for. Contiguous, because a chunk is a
contiguous block of that axis.

This is the axis index, not the complete level's numbering: on a cube over a
subset of a level the two differ, and it is the caller's axis a result is
written back to. [`axisindices`](@ref) names every cell of a loaded chunk in
the same space, halo included.
"""
ownedindices(mc::MapChunk) = mc.range

"""
    chunkhalo(mc::MapChunk) -> Vector{Int}

The axis indices outside [`ownedindices`](@ref) that this chunk's cells reach,
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
Each chunk writes only the indices it owns, so pieces running at once cannot
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
    out[ownedindices(cc)] = parent(r)[localindices(cc)]
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

# The axis indices a chunk's cells reach but does not own.
#
# Asked of the REGION rather than cell by cell: `halo` walks a boundary in
# `O(border)`, where enumerating every cell's ring and discarding the ones that
# land back inside is `O(cells x degree)`. The answer ascends without sorting,
# because a complete-level index and an axis index both ascend in cell id,
# so the map between them is monotone.
function _chunkhalo(cv::CellVector, r::UnitRange{Int}, width::Int,
        conn::Connectivity)
    out = Int[]
    isempty(r) && return out
    sub = cv[r]
    if width == 1
        for p in halo(sub; connectivity=conn)
            q = Engine.windowindex(cv.windows, p)
            q === nothing || push!(out, q)
        end
    else
        # Wider than one ring has no single boundary walk, so the region is
        # grown and what it gained outside the chunk is the halo.
        grown = grow(sub, width; connectivity=conn)
        for p in Engine.leafindices(grown.windows)
            q = Engine.windowindex(cv.windows, p)
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

A chunk is a partial grid, so an index into its cube is CHUNK-LOCAL and means
nothing to the caller. Its accessors name cells in the two index spaces:

| accessor | what it gives |
|---|---|
| [`chunkcube`](@ref) | the in-memory cube: the chunk's cells AND its halo |
| [`localindices`](@ref) | the indices IN THAT CUBE the chunk owns |
| [`ownedindices`](@ref) | those same owned cells' indices in the CALLER'S CELL AXIS |
| [`axisindices`](@ref) | EVERY cube cell's index in the caller's cell axis, halo included |

The owned pair is one set of cells in the two numberings, so
`axisindices(cc)[localindices(cc)] == ownedindices(cc)` always.

A result computed on the cube is meaningful for the owned indices only:
a halo cell is there to complete its neighbours' rings, and its own ring is
missing whatever fell outside the block. `r[localindices(cc)]` selects the
part to keep and `ownedindices(cc)` says where it belongs.

Because the halo holds every axis neighbour of every owned cell, a stencil that
reaches no further than the plan's [`halowidth`](@ref) computes exactly what the
whole-axis sweep computes for those cells.
"""
struct ChunkCube{C<:DD.AbstractDimArray}
    cube::C
    localindices::UnitRange{Int}
    range::UnitRange{Int}
    axisindices::Vector{Int}
    index::Int
end

"""
    chunkcube(cc::ChunkCube) -> AbstractDimArray

The chunk's cells and its halo, in memory, over a [`CellLookup`](@ref).
"""
chunkcube(cc::ChunkCube) = cc.cube

"""
    localindices(cc::ChunkCube) -> UnitRange{Int}

The indices in [`chunkcube`](@ref) that the chunk OWNS — the ones a result is
kept for. Contiguous: a chunk is a run of the axis, and every halo cell is
outside that run, so the owned cells stay together when the two are sorted into
one axis.
"""
localindices(cc::ChunkCube) = cc.localindices

ownedindices(cc::ChunkCube) = cc.range

"""
    axisindices(cc::ChunkCube) -> Vector{Int}

Where every cell of [`chunkcube`](@ref) came from in the CALLER'S CELL AXIS,
in the cube's own cell order: entry `k` is the axis index of the cube's cell
`k`, halo cells included. Ascending, since the halo below the chunk, the chunk
itself and the halo above it are stored in that order.

This is the translation between the two numberings, and it is what lets a
sweep over a chunk report indices the caller can use: `Index(Local())` on a
chunk is `Value(axisindices(cc))`, so a request answered chunk by chunk names
the same cells the whole-axis sweep would have.
"""
axisindices(cc::ChunkCube) = cc.axisindices

# `globalindices` is the old name of `ownedindices` and forwards to it, so a
# call of the old name answers the same with a deprecation warning. It was
# never the complete level's numbering that `globalindex` answers — the word
# meant "the whole axis", which is now said as the axis.

"""
    globalindices(mc::MapChunk) -> UnitRange{Int}
    globalindices(cc::ChunkCube) -> UnitRange{Int}

Deprecated. Use [`ownedindices`](@ref), which this forwards to, so existing
calls keep their old behaviour exactly.
"""
function globalindices end

@deprecate globalindices(mc::MapChunk) ownedindices(mc) false
@deprecate globalindices(cc::ChunkCube) ownedindices(cc) false

"""
    chunkhalo(cc::ChunkCube) -> Vector{Int}

The indices in [`chunkcube`](@ref) that are HALO — read to complete the owned
cells' rings, and not results to keep. The complement of
[`localindices`](@ref) over the cube, so the two together are all of it.

These are indices in the block, not in the axis; `chunkhalo` of the plan's
[`MapChunk`](@ref) gives the axis indices they came from.
"""
chunkhalo(cc::ChunkCube) = vcat(1:(first(cc.localindices)-1),
    (last(cc.localindices)+1):_cellcount(cc))

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
        write!(out, ownedindices(cc), mine(mapneighbors(f, chunkcube(cc)), cc))
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

# The three runs of axis indices a chunk's cube holds, in the cube's own cell
# order: the halo below the chunk, the chunk's own contiguous block, the halo
# above it. Named from the plan, before anything is read.
function _chunkparts(mc::MapChunk)
    below = _lessthan(mc.halo, first(mc.range))
    return below, mc.range, (@view mc.halo[(length(below)+1):end])
end

# The same three runs, named from a chunk already loaded. `localindices` is the
# owned run inside the cube, so what lies before and after it is the halo.
function _chunkparts(cc::ChunkCube)
    idx, own = cc.axisindices, cc.localindices
    below = @view idx[1:(first(own)-1)]
    above = @view idx[(last(own)+1):end]
    return below, cc.range, above
end

# Read one chunk's three runs out of a cell-axis array and join them.
#
# The pieces are read separately and joined rather than gathered in one
# scattered index, so the chunk's own cells stay a single contiguous read and
# each foreign chunk the halo reaches is read exactly once. Both the cube being
# swept and any array a field request names go through here, so a `Value` over
# a stored array costs one read per storage chunk it touches, never one per
# scalar.
function _readcells(data, dnum::Int, below, r::UnitRange{Int}, above,
        bounds::Vector{UnitRange{Int}})
    parts = Any[]
    isempty(below) || push!(parts, _readindices(data, dnum, below, bounds))
    push!(parts, _readrange(data, dnum, r))
    isempty(above) || push!(parts, _readindices(data, dnum, above, bounds))
    return length(parts) == 1 ? only(parts) : cat(parts...; dims=dnum)
end

# Read one chunk and its halo into an ordinary cube.
function _loadchunk(A::DD.AbstractDimArray, dnum::Int, cv::CellVector,
        mc::MapChunk, bounds::Vector{UnitRange{Int}})
    below, r, above = _chunkparts(mc)
    block = _readcells(parent(A), dnum, below, r, above, bounds)
    indices = vcat(collect(below), collect(r), collect(above))
    lk = CellLookup(cv[indices])
    dims = ntuple(i -> i == dnum ? Cells(lk) : DD.dims(A)[i], Val(ndims(A)))
    own = (length(below)+1):(length(below)+length(r))
    return ChunkCube(DD.rebuild(A; data=block, dims=dims), own, r, indices,
        mc.index)
end

_lessthan(halo::Vector{Int}, lo::Int) =
    @view halo[1:(searchsortedfirst(halo, lo)-1)]

function _readrange(data, dnum::Int, r::AbstractUnitRange)
    inds = ntuple(i -> i == dnum ? r : Colon(), Val(ndims(data)))
    return Array(data[inds...])
end

# Scattered indices, read one covering range per storage chunk they land in.
# Within a chunk the covering range is decoded as a unit anyway, so reading it
# whole costs nothing over reading part of it — and reading the runs separately
# is what stops a halo that straddles a distant chunk from pulling in everything
# between.
function _readindices(data, dnum::Int, idx::AbstractVector{Int},
        bounds::Vector{UnitRange{Int}})
    parts = Any[]
    for run in _indexruns(idx, bounds)
        block = _readrange(data, dnum, run)
        take = [p - first(run) + 1 for p in idx if first(run) <= p <= last(run)]
        inds = ntuple(i -> i == dnum ? take : Colon(), Val(ndims(block)))
        push!(parts, block[inds...])
    end
    return length(parts) == 1 ? only(parts) : cat(parts...; dims=dnum)
end

# Group ascending indices by the storage chunk they fall in, so each chunk a
# gather reaches is read once and no two reads overlap.
function _indexruns(idx::AbstractVector{Int}, bounds::Vector{UnitRange{Int}})
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
# A field request, translated onto a chunk
# ===========================================================================
#
# A request names quantities in the collection the CALLER passed. A chunk is a
# different collection — a partial grid whose `1:ncells` is its own — so a
# sweep over one has to be asked for the same quantities in ITS terms and give
# back answers in the caller's. That translation is this section, and it is one
# `map` over the request: each need either does not depend on the collection
# at all, or is a `Value` whose vector is restricted to the chunk's cells.
#
# `Index(Local())` is the whole point. The caller's axis index of every cell of
# the chunk's cube is `axisindices(cc)`, so the request becomes a `Value` of
# that vector and the sweep's existing machinery answers it — for the visited
# cell and for every ring slot — with the number the whole-axis sweep would
# have reported. No per-need code, and the pin holds by construction.

_chunkneeds(needs::Tuple, cc::ChunkCube) = map(n -> _chunkneed(n, cc), needs)

# A cell id, the complete level's numbering and a centroid are properties of
# the cell, not of the collection it was found in: a chunk answers them
# unchanged.
_chunkneed(n::Engine.Cell, ::ChunkCube) = n
_chunkneed(n::Engine.Index{Engine.Global}, ::ChunkCube) = n
_chunkneed(n::Engine.Index{Type{T}}, ::ChunkCube) where {T} = n
_chunkneed(n::Engine.Centroid, ::ChunkCube) = n

_chunkneed(::Engine.Index{Engine.Local}, cc::ChunkCube) =
    Engine.Value(axisindices(cc))

_chunkneed(n::Engine.Value, cc::ChunkCube) = Engine.Value(_restrict(n.data, cc))

# The chunk cube's own collection: what a sweep over the chunk runs on, and
# what a field restricted onto the chunk is built over. The lookup holds it, so
# every call answers the same object and a field built here is accepted by the
# sweep's `===` check without falling back to comparing windows.
_chunkcv(cc::ChunkCube) = region(DD.lookup(cc.cube, _celldim(cc)))

_rawdata(a) = a
_rawdata(a::DD.AbstractDimArray) = parent(a)

_ischunked(a) = DiskArrays.haschunks(_rawdata(a)) isa DiskArrays.Chunked

# A vector laid out on the caller's axis, restricted to the chunk's cells and
# laid out on the chunk cube's. In memory that is a gather; on disk it is the
# same three-part read the chunk's own data took, so the request costs one read
# per storage chunk it touches rather than one per scalar.
#
# The gather COPIES rather than viewing. On time the two are the same: over a
# 117,649-cell in-memory sweep a view ran at 0.994 of the copy's median, inside
# the run-to-run spread either way. What decides it is the type — a copy and
# the disk-backed read below both answer a `Vector{eltype(a)}`, so the
# translated request's types do not depend on which branch ran, which is what
# lets the automatic route size its result from the whole axis before a chunk
# exists — and lifetime: a view would keep the caller's whole array reachable
# through every chunk's request and read each entry through one more
# indirection on each of its `degree + 1` touches, against 8 bytes per block
# cell copied once.
function _restrict(a::AbstractVector, cc::ChunkCube)
    v = _ischunked(a) ? _readalong(a, cc) : a[axisindices(cc)]
    return convert(Vector{eltype(a)}, v)
end

function _readalong(a, cc::ChunkCube)
    data = _rawdata(a)
    below, r, above = _chunkparts(cc)
    # A `Value` names one entry per cell, so its cell axis is its only one.
    return vec(_readcells(data, 1, below, r, above, _storagebounds(data, 1)))
end

# A cell field is a function of the cell plus whatever the caller already knew,
# so a chunk rebuilds it over the chunk cube's collection. `known` follows the
# rule its own constructor states: a dense vector is laid out on the collection
# and is restricted with it, and a subset cube says which cells it holds by
# cell id, which no chunk changes.
_restrict(a::Engine.CellField, cc::ChunkCube) =
    Engine.cellfield(Engine._fieldfunction(a), _chunkcv(cc);
        known=_restrictknown(Engine._fieldknown(a), cc))

_restrictknown(::Nothing, ::ChunkCube) = nothing
_restrictknown(k::AbstractVector, cc::ChunkCube) = _restrict(k, cc)
_restrictknown(k::Engine._SubsetKnown, ::ChunkCube) = k

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
    H = Engine.SubsetIndexedCell{eltype(cv)}
    T = Base.promote_op(f, H, eltype(A), Engine._ringtype(cap, eltype(A)))
    sz = size(A)
    dest = T <: Tuple && isconcretetype(T) ?
           ntuple(j -> Array{fieldtype(T, j)}(undef, sz), fieldcount(T)) :
           Array{T}(undef, sz)
    mapneighbors!(dest, f, A, plan; threaded=threaded)
    return dest
end

function CellLookups._map_needs_chunked(f::F, A, dnum::Int,
        plan::MapChunkPlan, needs, threaded, conn::Connectivity) where {F}
    cv = region(DD.lookup(A, dnum))
    # What the caller wrote is checked against the collection they wrote it
    # about, before a chunk is read, so an error names the request as written.
    Engine._checkneeds(needs, cv)
    # The translated request's element types are the request's own: every rule
    # in `_chunkneed` preserves the need's element type — `Index(Local())`
    # becomes a `Value` of `Int`s, a `Value` keeps its array's eltype, a field
    # keeps the one it was built with, and the rest are unchanged — so the
    # callback's argument types are the same on every chunk AND the same as
    # the whole-axis sweep's. They can therefore be derived here, from the
    # whole axis, before any chunk exists.
    rn = Engine._resolveneeds(needs, cv)
    cap = Engine._capacity(system(cv), conn)
    T = Base.promote_op(f, Engine._centertype(rn, cv),
        Engine._ringstype(rn, cv, cap))
    n = length(cv)
    dest = T <: Tuple && isconcretetype(T) ?
           ntuple(j -> Array{fieldtype(T, j)}(undef, n), fieldcount(T)) :
           Array{T}(undef, n)
    _mapchunks!(needs, dest, f, A, plan, nothing, threaded)
    return dest
end

# The side-effecting form. A sweep over a chunk visits the halo too — those
# cells are there to complete the owned cells' rings — and a callback that is
# only called for its effects must not see them, so the owned range is handed
# down as the visit filter. `mapneighbors!` needs no such filter: it drops the
# halo rows when it stores.
function CellLookups._foreach_needs_chunked(f::F, A, dnum::Int,
        plan::MapChunkPlan, needs, threaded, conn::Connectivity) where {F}
    Engine._checkneeds(needs, region(DD.lookup(A, dnum)))
    foreachchunk(A, plan) do cc
        Engine._foreachneighbors(f, _chunkcv(cc), _chunkneeds(needs, cc),
            StorageOrder(), threaded, plan.connectivity,
            Base.Fix2(in, localindices(cc)))
    end
    return nothing
end

"""
    mapneighbors!(dest, f, A::AbstractDimArray; halo = 1, chunks = :auto,
                  spatialdim = nothing, connectivity = Vertex(), threaded = true)
    mapneighbors!(dest, f, A, plan::MapChunkPlan;
                  needs = (Value(a), Centroid()), threaded = true)

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

`needs` names the per-neighbour fields the kernel reads instead, and the
callback becomes `f(center, rings)` — [`mapneighbors`](@ref)'s field-request
contract, here on the chunk route. The request is stated once, about the cube
the caller passed, and translated onto each chunk: `Index(Local())` answers
that cube's cell-axis index on every chunk, never a chunk-local one, and a
[`Value`](@ref) over a stored array is read along its own storage chunks the
way the swept data is. `dest` then holds one result per cell.

The result is the whole-axis sweep's, cell for cell, in either form. Each
chunk's halo carries every axis neighbour of every cell the chunk owns, so a
ring computed on a chunk is the ring computed on the axis — clipped
identically, and in the same order.

`threaded` threads WITHIN each chunk. To run chunks themselves in parallel,
`split` a [`chunkplan`](@ref) and call this on the pieces. Build the plan before
splitting: doing so is what fills the axis's [`region`](@ref) memo, so the
pieces share one conversion rather than each repeating it.
"""
function mapneighbors!(dest, f::F, A::DD.AbstractDimArray; halo::Integer=1,
        chunks=:auto, spatialdim=nothing, needs=nothing,
        connectivity::Connectivity=Vertex(), threaded=true) where {F}
    plan = chunkplan(A; halo=halo, chunks=chunks, spatialdim=spatialdim,
        connectivity=connectivity)
    return mapneighbors!(dest, f, A, plan; spatialdim=spatialdim, needs=needs,
        threaded=threaded)
end

function mapneighbors!(dest, f::F, A::DD.AbstractDimArray, plan::MapChunkPlan;
        spatialdim=nothing, needs=nothing, threaded=true) where {F}
    _mapchunks!(needs, dest, f, A, plan, spatialdim, threaded)
    return dest
end

# No field request: the `Values()` form, untouched — `needs = nothing` reaches
# it by dispatch, not by a branch.
function _mapchunks!(::Nothing, dest, f::F, A, plan::MapChunkPlan, spatialdim,
        threaded) where {F}
    foreachchunk(A, plan; spatialdim=spatialdim) do cc
        cube = chunkcube(cc)
        out = mapneighbors(f, cube; pass=Values(), threaded=threaded,
            connectivity=plan.connectivity)
        _storechunk!(dest, out, cc)
    end
    return nothing
end

# A field request. The sweep runs on the chunk cube's own collection rather
# than the cube — a request never reads the array it is called on, so there is
# nothing for the cube's other dimensions to do — and the answer is one result
# per cell, of which the owned rows are kept.
function _mapchunks!(needs, dest, f::F, A, plan::MapChunkPlan, spatialdim,
        threaded) where {F}
    dnum = CellLookups._cells_dimnum(A, spatialdim)
    Engine._checkneeds(needs, region(DD.lookup(A, dnum)))
    foreachchunk(A, plan; spatialdim=spatialdim) do cc
        out = mapneighbors(f, _chunkcv(cc); needs=_chunkneeds(needs, cc),
            threaded=threaded, connectivity=plan.connectivity)
        _storeneeds!(dest, out, cc)
    end
    return nothing
end

# Keep the owned rows and put them where the chunk came from. A tuple result is
# one output per component, so `dest` is then a tuple of the same length.
function _storechunk!(dest, out, cc::ChunkCube)
    dnum = CellLookups._cells_dimnum(out, nothing)
    src = _slice(parent(out), dnum, localindices(cc))
    _destview(dest, dnum, ownedindices(cc)) .= src
    return nothing
end

function _storechunk!(dest::Tuple, out::Tuple, cc::ChunkCube)
    map((d, o) -> _storechunk!(d, o, cc), dest, out)
    return nothing
end

# A field request gives one result per cell, so its result is a plain vector
# over the chunk cube's cells and `dest` is indexed along the cell axis alone.
function _storeneeds!(dest, out::AbstractVector, cc::ChunkCube)
    _destview(dest, 1, ownedindices(cc)) .= @view out[localindices(cc)]
    return nothing
end

function _storeneeds!(dest::Tuple, out::Tuple, cc::ChunkCube)
    map((d, o) -> _storeneeds!(d, o, cc), dest, out)
    return nothing
end

_slice(data, dnum::Int, idx) =
    data[ntuple(i -> i == dnum ? idx : Colon(), Val(ndims(data)))...]

_destview(dest::DD.AbstractDimArray, dnum::Int, r) = _destview(parent(dest), dnum, r)
_destview(dest, dnum::Int, r) =
    view(dest, ntuple(i -> i == dnum ? r : Colon(), Val(ndims(dest)))...)
