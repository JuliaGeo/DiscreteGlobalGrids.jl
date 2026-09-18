# The DEM as one lazy vector: ids, an in-memory tile cache, the disk array, and
# the destination chunks that cover a set of tiles.

"""
    SubtreeIds(sys, parents, level)

The level-`level` descendants of `parents` as one lazy ascending vector of cell
ids, subtree `k` being `parents[k]`. `parents` must ascend.

Storage is two integers per parent, so the 3.8e10 pixels of GLO-90 cost a table
of 26 475 entries. The vector declares itself sorted, which spares a
`PartialGrid` its O(n) scan.
"""
struct SubtreeIds{G,ID} <: AbstractVector{ID}
    complete::G
    starts::Vector{Int}      # first level-`level` index of each subtree
    offsets::Vector{Int}     # cells before each subtree; `offsets[1] == 0`
    n::Int
end

function SubtreeIds(sys, parents::AbstractVector, level::Int)
    complete = DGG.levelgrid(sys, level)
    starts = Vector{Int}(undef, length(parents))
    offsets = Vector{Int}(undef, length(parents))
    acc = 0
    for (k, p) in enumerate(parents)
        r = DGG.descendant_range(sys, p, level)
        starts[k] = Int(first(r))
        offsets[k] = acc
        acc += length(r)
    end
    # The cell type differs by system (`LevelIndex`, `Z7Cell`, ...).
    ID = typeof(DGG.cellindex(complete, first(starts)))
    return SubtreeIds{typeof(complete),ID}(complete, starts, offsets, acc)
end

"The pixel (level-1) cell ids of the level-0 `tiles`."
TileIds(sys, tiles::AbstractVector{<:Integer}) =
    SubtreeIds(sys, [DGG.LevelIndex(0, t) for t in tiles], 1)

Base.size(v::SubtreeIds) = (v.n,)
Base.IndexStyle(::Type{<:SubtreeIds}) = IndexLinear()

"Which subtree holds index `p`, and the offset within it."
@inline function tileat(v::SubtreeIds, p::Int)
    k = searchsortedlast(v.offsets, p - 1)
    return k, p - v.offsets[k]
end

@inline function Base.getindex(v::SubtreeIds, i::Int)
    @boundscheck (1 <= i <= v.n) || throw(BoundsError(v, i))
    k, off = tileat(v, i)
    return DGG.cellindex(v.complete, v.starts[k] + off - 1)
end

DGG.Helpers.strictly_increasing(::SubtreeIds) = true

"Cells in subtree `k`."
subtreelength(v::SubtreeIds, k::Int) =
    (k < length(v.offsets) ? v.offsets[k + 1] : v.n) - v.offsets[k]

"""
    StripedLRUCache{V}(load; slots = 256, stripes = 16)

A bounded in-memory cache, callable as `cache(k)`: `stripes` independent LRUs of
`slots ÷ stripes` entries, one lock apiece, keyed by integer.

`load(k)` runs outside the stripe lock, so two racing callers may both load the
same key; the first to publish wins and the other copy is dropped. Peak
residency is `slots` values whatever the access pattern.
"""
struct StripedLRUCache{V,F}
    load::F
    caches::Vector{Dict{Int,V}}
    orders::Vector{Vector{Int}}
    locks::Vector{ReentrantLock}
    per::Int
    loads::Threads.Atomic{Int}
    hits::Threads.Atomic{Int}
end

function StripedLRUCache{V}(load; slots::Integer = 256, stripes::Integer = 16) where {V}
    ns = max(1, Int(stripes))
    return StripedLRUCache{V,typeof(load)}(load,
        [Dict{Int,V}() for _ in 1:ns], [Int[] for _ in 1:ns],
        [ReentrantLock() for _ in 1:ns], max(1, Int(slots) ÷ ns),
        Threads.Atomic{Int}(0), Threads.Atomic{Int}(0))
end

function (c::StripedLRUCache)(k::Integer)
    ki = Int(k)
    st = mod(ki, length(c.locks)) + 1
    lk, cache, order = c.locks[st], c.caches[st], c.orders[st]
    hit = lock(lk) do
        v = get(cache, ki, nothing)
        v === nothing && return nothing
        push!(order, splice!(order, findfirst(==(ki), order)))
        return v
    end
    if hit !== nothing
        Threads.atomic_add!(c.hits, 1)
        return hit
    end
    Threads.atomic_add!(c.loads, 1)
    v = c.load(ki)
    return lock(lk) do
        existing = get(cache, ki, nothing)
        existing === nothing || return existing
        length(cache) >= c.per && delete!(cache, popfirst!(order))
        cache[ki] = v
        push!(order, ki)
        return v
    end
end

"`(; loads, hits, live, bytes)` for a [`StripedLRUCache`](@ref)."
cachestats(c::StripedLRUCache) = (loads = c.loads[], hits = c.hits[],
    live = sum(length, c.caches),
    bytes = sum(cache -> sum(sizeof, values(cache); init = 0), c.caches))

"""
    TiledDEM(ids::SubtreeIds, gettile)
    TiledDEM(source, tiles; slots = 256, stripes = 16)

Every listed tile's posts as one lazy `Float32` vector in the source grid's own
index order. Chunk `k` is tile `k`, so reads are tile aligned and a regridder's
source chunks are the DEM's own tiles.

`gettile(k)` returns chunk `k`'s posts and owns all caching policy. The second
form caches `loadtile(source, tiles[k])` in a [`StripedLRUCache`](@ref).

```julia
sys   = DGG.CopernicusDEMSystem(90)
tiles = listedtiles(sys, tilelist(datadir))
dem   = TiledDEM(CopernicusTiles(sys, tiles; cachedir), tiles)
space = DGG.DGGSpace(DGG.PartialGrid(sys, 1, dem.ids); chunklevel = 0)
```
"""
struct TiledDEM{I<:SubtreeIds,F,C} <: DiskArrays.AbstractDiskArray{Float32,1}
    ids::I
    gettile::F
    chunks::C
end

function TiledDEM(ids::SubtreeIds, gettile)
    widths = [subtreelength(ids, k) for k in eachindex(ids.offsets)]
    chunks = DiskArrays.GridChunks(DiskArrays.IrregularChunks(; chunksizes = widths))
    return TiledDEM(ids, gettile, chunks)
end

function TiledDEM(source, tiles::AbstractVector{<:Integer}; slots = 256, stripes = 16)
    cache = StripedLRUCache{Vector{Float32}}(k -> loadtile(source, Int(tiles[k]));
        slots, stripes)
    return TiledDEM(TileIds(source.sys, tiles), cache)
end

Base.size(A::TiledDEM) = (A.ids.n,)
DiskArrays.eachchunk(A::TiledDEM) = A.chunks
DiskArrays.haschunks(::TiledDEM) = DiskArrays.Chunked()

function DiskArrays.readblock!(A::TiledDEM, out, r::AbstractUnitRange)
    p = first(r)
    while p <= last(r)
        k, off = tileat(A.ids, p)
        stop = min(last(r), A.ids.offsets[k] + subtreelength(A.ids, k))
        seg = (p - first(r) + 1):(stop - first(r) + 1)
        v = A.gettile(k)
        out[seg] .= view(v, off:(off + length(seg) - 1))
        p = stop + 1
    end
    return out
end

"""
    covering_chunks(dstsys, sys, tiles, level; nthreads = 1) -> Vector{Int}

The level-`level` cells of `dstsys` that meet any of `tiles`, as ascending local
indices of `DGG.levelgrid(dstsys, level)`.

Each tile's 1x1-degree extent is queried at `level`, so the answer is whole
cells and each one is a complete subtree below it: one work unit, one store
chunk. Every cell returned meets a listed tile, so ocean-only work is never
enqueued.
"""
function covering_chunks(dstsys, sys, tiles::AbstractVector{<:Integer}, level::Int;
        nthreads::Integer = 1)
    g = DGG.levelgrid(dstsys, level)
    parts = [Set{Int}() for _ in 1:nthreads]
    Threads.@sync for w in 1:nthreads
        Threads.@spawn for k in w:nthreads:length(tiles)
            lat, lon = CD.tilecorner(sys, DGG.LevelIndex(0, tiles[k]))
            ex = Extents.Extent(X = (Float64(lon), Float64(lon) + 1),
                Y = (Float64(lat), Float64(lat) + 1))
            set = DGG.query(dstsys, DGG.MultiOrderCoverage(ex); level)
            for c in DGG.CellVector(set; level)
                push!(parts[w], DGG.localindex(g, c))
            end
        end
    end
    return sort!(collect(union(parts...)))
end
