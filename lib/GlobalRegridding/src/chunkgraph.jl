# The chunk dependency graph: a materialized, conservative bipartite relation
# between destination chunks and the source chunks that may contribute to them.
#
# `discovery.jl` answers "which source chunks does this destination chunk need?"
# one query at a time. That is the right shape for a lazy read, but a scheduler
# needs the whole relation at once: refcounts need every consumer of a source
# chunk, prefetch needs lookahead, and an affinity order needs both directions.
# This file materializes the relation once, from chunk extents alone.

"""
    ChunkDependencyGraph{T} <: Graphs.AbstractGraph{T}

A conservative bipartite dependency relation between the chunks of a destination
space and the chunks of a source space, materialized as bidirectional CSR.

An edge `(s, d)` means source chunk `s` *may* contribute to destination chunk
`d`. The relation is a **superset**: it never omits a contributing pair, but it
may include pairs whose geometries do not actually overlap, because it is
computed from covering spherical caps rather than cell geometry. Treat an edge
as "must be loaded", never as "has nonzero weight".

# Vertex numbering

The graph implements the Graphs.jl `AbstractGraph` interface over a single
vertex range, so that ecosystem algorithms work directly on it:

| vertices | meaning |
|---|---|
| `1 : nsrc` | source chunk `v` |
| `nsrc+1 : nsrc+ndst` | destination chunk `v - nsrc` |

Source vertices are numbered first so that a destination's neighbours are
already source *chunk* numbers with no offset arithmetic. Use
[`sourcesof`](@ref) and [`consumersof`](@ref) for the chunk-numbered views, and
[`srcvertex`](@ref)/[`dstvertex`](@ref)/[`srcchunk`](@ref)/[`dstchunk`](@ref) to
convert. The graph is undirected: every edge has one source and one destination
endpoint, so `has_edge` and `edges` are role-agnostic.

The spaces themselves are deliberately *not* stored. The graph is a plain
combinatorial object identified by `(nsrc, ndst, radius)`, which keeps it small,
serializable, and free of the space type parameters.

# Identity

The relation depends on the regridding method only through its support radius,
so [`dependency_radius`](@ref) is part of the graph's identity: two methods with
the same support radius share a graph, and a graph built at one radius is not
valid at a larger one.

See [`chunk_dependency_graph`](@ref).
"""
struct ChunkDependencyGraph{T<:Integer} <: Graphs.AbstractGraph{T}
    nsrc::Int
    ndst::Int
    radius::Float64
    # Destination-major CSR: sources of destination chunk `d` are
    # `srcof[dstoff[d]:dstoff[d+1]-1]`, ascending source chunk numbers.
    dstoff::Vector{Int}
    srcof::Vector{T}
    # Source-major CSR: consumers of source chunk `s` are
    # `dstof[srcoff[s]:srcoff[s+1]-1]`, ascending destination chunk numbers.
    srcoff::Vector{Int}
    dstof::Vector{T}
end

Base.show(io::IO, g::ChunkDependencyGraph) =
    print(io, "ChunkDependencyGraph(", g.nsrc, " source × ", g.ndst,
        " destination chunks, ", length(g.srcof), " edges, radius ", g.radius, ")")

"""
    nsourcechunks(g::ChunkDependencyGraph) -> Int
    ndestinationchunks(g::ChunkDependencyGraph) -> Int

Return the size of each side of the bipartition.
"""
nsourcechunks(g::ChunkDependencyGraph) = g.nsrc
ndestinationchunks(g::ChunkDependencyGraph) = g.ndst

"""
    dependency_radius(g::ChunkDependencyGraph) -> Float64

Return the support radius, in radians, the relation was built at. A graph is
valid for any method whose [`support_radius`](@ref) is at most this value.
"""
dependency_radius(g::ChunkDependencyGraph) = g.radius

# Vertex/chunk conversion

"""
    srcvertex(g, srcchunk) -> Int
    dstvertex(g, dstchunk) -> Int

Return the graph vertex for a source or destination chunk number.
"""
function srcvertex(g::ChunkDependencyGraph, s::Integer)
    i = Int(s)
    1 <= i <= g.nsrc || throw(BoundsError(g, s))
    return i
end

function dstvertex(g::ChunkDependencyGraph, d::Integer)
    i = Int(d)
    1 <= i <= g.ndst || throw(BoundsError(g, d))
    return g.nsrc + i
end

"""
    issrcvertex(g, v) -> Bool
    isdstvertex(g, v) -> Bool

Return whether vertex `v` is on the source or the destination side.
"""
issrcvertex(g::ChunkDependencyGraph, v::Integer) = 1 <= Int(v) <= g.nsrc
isdstvertex(g::ChunkDependencyGraph, v::Integer) = g.nsrc < Int(v) <= g.nsrc + g.ndst

"""
    srcchunk(g, v) -> Int
    dstchunk(g, v) -> Int

Return the chunk number of vertex `v`, throwing if `v` is on the other side of
the bipartition. The role check is the point: silently reinterpreting a
destination vertex as a source chunk is the easiest way to corrupt a schedule.
"""
function srcchunk(g::ChunkDependencyGraph, v::Integer)
    issrcvertex(g, v) || throw(ArgumentError(
        "vertex $v is not a source vertex of $g; source vertices are 1:$(g.nsrc)"))
    return Int(v)
end

function dstchunk(g::ChunkDependencyGraph, v::Integer)
    isdstvertex(g, v) || throw(ArgumentError(
        "vertex $v is not a destination vertex of $g; destination vertices are " *
        "$(g.nsrc + 1):$(g.nsrc + g.ndst)"))
    return Int(v) - g.nsrc
end

"""
    srcvertices(g) -> UnitRange{Int}
    dstvertices(g) -> UnitRange{Int}

Return each side's vertex range.
"""
srcvertices(g::ChunkDependencyGraph) = 1:g.nsrc
dstvertices(g::ChunkDependencyGraph) = (g.nsrc+1):(g.nsrc+g.ndst)

# Adjacency, in chunk numbers

"""
    sourcesof(g, dstchunk) -> AbstractVector{<:Integer}

Return the ascending source chunk numbers that may contribute to destination
chunk `dstchunk`. The result is a read-only view into the graph; do not mutate
it.

This is the executor's query: the set of source chunks that must be resident to
compute this destination chunk.
"""
function sourcesof(g::ChunkDependencyGraph, d::Integer)
    i = Int(d)
    1 <= i <= g.ndst || throw(BoundsError(g, d))
    return view(g.srcof, g.dstoff[i]:(g.dstoff[i+1]-1))
end

"""
    consumersof(g, srcchunk) -> AbstractVector{<:Integer}

Return the ascending destination chunk numbers that may need source chunk
`srcchunk`. The result is a read-only view into the graph; do not mutate it.

`length(consumersof(g, s))` is the initial refcount for source chunk `s` when
the whole destination space is being computed. See the eviction notes in
[`chunk_dependency_graph`](@ref) before using it as one.
"""
function consumersof(g::ChunkDependencyGraph, s::Integer)
    i = Int(s)
    1 <= i <= g.nsrc || throw(BoundsError(g, s))
    return view(g.dstof, g.srcoff[i]:(g.srcoff[i+1]-1))
end

"""
    sourcedegree(g, dstchunk) -> Int
    consumerdegree(g, srcchunk) -> Int

Return adjacency counts without materializing the rows. `sourcedegree` is a
cheap per-destination cost and I/O estimate; `consumerdegree` is a source
chunk's consumer count.
"""
sourcedegree(g::ChunkDependencyGraph, d::Integer) =
    (i = Int(d); 1 <= i <= g.ndst || throw(BoundsError(g, d));
     g.dstoff[i+1] - g.dstoff[i])
consumerdegree(g::ChunkDependencyGraph, s::Integer) =
    (i = Int(s); 1 <= i <= g.nsrc || throw(BoundsError(g, s));
     g.srcoff[i+1] - g.srcoff[i])

# Graphs.jl `AbstractGraph` interface
#
# Read-only by design. Mutating methods (`add_edge!`, `rem_vertex!`, ...) are
# deliberately not defined: they could not preserve the bipartition, and the
# relation is a derived geometric fact, not something to edit.

Graphs.is_directed(::Type{<:ChunkDependencyGraph}) = false
Graphs.is_directed(::ChunkDependencyGraph) = false
Graphs.edgetype(::ChunkDependencyGraph{T}) where {T} = Graphs.SimpleEdge{T}
Graphs.nv(g::ChunkDependencyGraph{T}) where {T} = T(g.nsrc + g.ndst)
Graphs.ne(g::ChunkDependencyGraph) = length(g.srcof)
Graphs.vertices(g::ChunkDependencyGraph{T}) where {T} = Base.OneTo(Graphs.nv(g))
Graphs.has_vertex(g::ChunkDependencyGraph, v::Integer) = 1 <= Int(v) <= g.nsrc + g.ndst

Base.zero(::Type{ChunkDependencyGraph{T}}) where {T} =
    ChunkDependencyGraph{T}(0, 0, 0.0, [1], T[], [1], T[])
Base.zero(g::ChunkDependencyGraph{T}) where {T} = zero(ChunkDependencyGraph{T})

# A destination's neighbours are source chunk numbers, which *are* source vertex
# numbers, so that side needs no translation. A source's neighbours are
# destination chunk numbers, which need the `nsrc` offset; hand back a lazy
# shifted view rather than allocating a shifted copy on every query.
struct ShiftedChunks{T,V<:AbstractVector{T}} <: AbstractVector{T}
    parent::V
    shift::T
end

Base.size(v::ShiftedChunks) = size(v.parent)
Base.IndexStyle(::Type{<:ShiftedChunks}) = IndexLinear()
Base.@propagate_inbounds Base.getindex(v::ShiftedChunks, i::Int) = v.parent[i] + v.shift

function Graphs.outneighbors(g::ChunkDependencyGraph{T}, v::Integer) where {T}
    Graphs.has_vertex(g, v) || throw(BoundsError(g, v))
    issrcvertex(g, v) && return ShiftedChunks(consumersof(g, v), T(g.nsrc))
    return sourcesof(g, Int(v) - g.nsrc)
end

# Undirected: in- and out-neighbours coincide.
Graphs.inneighbors(g::ChunkDependencyGraph, v::Integer) = Graphs.outneighbors(g, v)

function Graphs.has_edge(g::ChunkDependencyGraph, u::Integer, v::Integer)
    (Graphs.has_vertex(g, u) && Graphs.has_vertex(g, v)) || return false
    # Exactly one endpoint must be a source vertex; same-side pairs are never
    # adjacent in a bipartite relation.
    if issrcvertex(g, u) && isdstvertex(g, v)
        s, d = Int(u), Int(v) - g.nsrc
    elseif issrcvertex(g, v) && isdstvertex(g, u)
        s, d = Int(v), Int(u) - g.nsrc
    else
        return false
    end
    # Search the shorter row.
    if sourcedegree(g, d) <= consumerdegree(g, s)
        return insorted(s, sourcesof(g, d))
    else
        return insorted(d, consumersof(g, s))
    end
end

# Every edge runs from a source vertex to a destination vertex, and source
# vertices are numbered first, so emitting from the source-major CSR yields each
# undirected edge exactly once with `src(e) < dst(e)`.
struct ChunkEdgeIter{T}
    graph::ChunkDependencyGraph{T}
end

Base.eltype(::Type{ChunkEdgeIter{T}}) where {T} = Graphs.SimpleEdge{T}
Base.length(it::ChunkEdgeIter) = Graphs.ne(it.graph)
Graphs.edges(g::ChunkDependencyGraph) = ChunkEdgeIter(g)

function Base.iterate(it::ChunkEdgeIter{T}, state::Tuple{Int,Int} = (1, 1)) where {T}
    g = it.graph
    s, k = state
    # Skip source chunks with no consumers.
    while s <= g.nsrc && k >= g.srcoff[s+1]
        s += 1
        k = s <= g.nsrc ? g.srcoff[s] : k
    end
    s > g.nsrc && return nothing
    e = Graphs.SimpleEdge{T}(T(s), T(g.nsrc) + g.dstof[k])
    return e, (s, k + 1)
end

# Construction

"""
    chunk_dependency_graph(dst_space, src_space; radius = 0.0, refine = nothing,
                           prefilter = true) -> ChunkDependencyGraph
    chunk_dependency_graph(plan::ChunkedPlan; refine = nothing) -> ChunkDependencyGraph

Build the conservative bipartite dependency relation between the chunks of
`dst_space` and those of `src_space`, from chunk extents alone. No cell
geometry is built and no data is read.

The `plan` form takes the radius from the plan's method, which is the safe entry
point; the space form is for callers that know the radius they want. Argument
order matches the rest of `discovery.jl`: destination first.

# Keywords

- `radius`: the method's angular support radius in radians. Two caps are
  connected when their centres are within the sum of their radii plus this. Must
  be finite and non-negative.
- `refine`: an optional conservative narrow phase, `refine(dstchunk, srcchunk)
  -> Bool`. Returning `false` asserts that the pair *cannot* contribute at this
  radius, and drops the edge. Returning `true` keeps it. The default keeps every
  candidate the cap test accepts. A wrong `refine` silently corrupts results, so
  it must only ever reject pairs it can prove disconnected.
- `prefilter`: use the latitude-band prefilter (default `true`). Setting it to
  `false` tests every destination against every source chunk; this is the
  brute-force reference the prefilter is checked against, and is only useful for
  small spaces and tests.

# Method

Both sides are reduced to their covering spherical caps via
[`chunkextents`](@ref) — the one seam every `RegridSpace` provides — so the
builder is generic over any source/destination pair the regridder accepts.

Pairs are then found by a **cap join** rather than by the dual-tree descent of
[`connectedchunkpairs`](@ref). Source chunks are sorted by cap-centre latitude
once; each destination cap then binary-searches the latitude band it could
possibly reach and applies the exact [`DilatedIntersects`](@ref) test to the
candidates. The band is a valid bound because the latitude gap between two
points is a lower bound on their spherical distance, so a pair separated by more
than `maxsrcradius + dstradius + radius` in latitude cannot pass the cap test.
Every candidate still gets the exact predicate, so the prefilter changes speed,
not results.

This is deliberate: a space is free to return a *flat* chunk tree — a single
node holding every chunk cap — and several do, which collapses a dual-tree
descent into a full quadratic scan. The cap join is uniformly
`O(nsrc log nsrc + ndst log nsrc + candidates)` regardless of tree shape.

Destination rows are built in parallel over blocks of destination chunks and
concatenated in chunk order, so the result is bit-identical regardless of thread
count.

# Using it for eviction and ordering

Two hazards, both measured:

- **Refcount eviction is not a bounded cache policy.** Holding a source chunk
  until its last consumer finishes guarantees each source is loaded exactly once
  — under refcount eviction the load count equals the source chunk count, by
  construction — but it puts *no* bound on residency. A locality-destroying
  order can hold most of the source space at once.
- **Do not order destinations by degree.** Degree is a good cost estimate and a
  terrible schedule: descending-degree order plus refcount eviction was measured
  at 85 GiB peak residency on the Copernicus DEM workload, against 1.43 GiB for
  a Morton-order spatial-affinity schedule. Order by spatial affinity and use
  degree only to estimate cost.

Any ordering derived from this graph is a *priority order* for a dynamic
pull-cursor, not a static partition of work; workers must still pull, or the
cost skew across the destination space stops being absorbed.

# Example

```julia
graph = chunk_dependency_graph(dst_space, src_space; radius = 0.0)
sourcesof(graph, 1)              # source chunks needed by destination chunk 1
length(consumersof(graph, 7))    # initial refcount for source chunk 7
Graphs.connected_components(graph)   # independent groups of work
```
"""
function chunk_dependency_graph(dst_space::RegridSpace, src_space::RegridSpace;
        radius::Real = 0.0, refine = nothing, prefilter::Bool = true)
    r = Float64(radius)
    (isfinite(r) && r >= 0) || throw(ArgumentError(
        "radius must be finite and non-negative, got $radius"))
    return _chunkgraph(chunkextents(dst_space), chunkextents(src_space), r,
        refine, prefilter)
end

chunk_dependency_graph(plan::ChunkedPlan; refine = nothing, prefilter::Bool = true) =
    chunk_dependency_graph(plan.dst_space, plan.src_space;
        radius = support_radius(plan.method, plan.src_space), refine, prefilter)

# Latitude of a unit-sphere point. `atan(z, hypot(x, y))` is the numerically
# robust spelling: `asin(z)` loses precision near the poles.
_caplatitude(c::SphericalCap) = atan(c.point[3], hypot(c.point[1], c.point[2]))

function _chunkgraph(dstcaps::Vector{<:SphericalCap}, srccaps::Vector{<:SphericalCap},
        radius::Float64, refine, prefilter::Bool)
    T = Int32
    ndst, nsrc = length(dstcaps), length(srccaps)
    (nsrc <= typemax(T) && ndst <= typemax(T)) || throw(ArgumentError(
        "chunk counts $nsrc and $ndst exceed the $(T) vertex numbering"))
    # Sort source chunks by cap-centre latitude once. `order[j]` is the source
    # chunk number of the j-th band entry.
    lats = map(_caplatitude, srccaps)
    order = sortperm(lats)
    sortedlats = lats[order]
    maxsrcradius = isempty(srccaps) ? 0.0 : maximum(c -> c.radius, srccaps)

    rows = Vector{Vector{T}}(undef, ndst)
    # Block the destination range so each task owns a contiguous span; rows are
    # written by index, so the concatenation below is order-independent.
    nblocks = max(1, min(ndst, 8 * Threads.nthreads()))
    blocks = _blockranges(ndst, nblocks)
    Threads.@sync for b in blocks
        let b = b
            Threads.@spawn begin
                # One scratch buffer per task, not per thread, so task migration
                # cannot alias it.
                buf = T[]
                for d in b
                    _fillrow!(buf, rows, d, dstcaps[d], srccaps, order,
                        sortedlats, maxsrcradius, radius, refine, prefilter)
                end
            end
        end
    end

    dstoff = Vector{Int}(undef, ndst + 1)
    dstoff[1] = 1
    for d in 1:ndst
        dstoff[d+1] = dstoff[d] + length(rows[d])
    end
    nedges = dstoff[end] - 1
    srcof = Vector{T}(undef, nedges)
    for d in 1:ndst
        copyto!(srcof, dstoff[d], rows[d], 1, length(rows[d]))
    end

    srcoff, dstof = _transpose(dstoff, srcof, nsrc, ndst, T)
    return ChunkDependencyGraph{T}(nsrc, ndst, radius, dstoff, srcof, srcoff, dstof)
end

# Split `1:n` into `k` contiguous ranges of near-equal length.
function _blockranges(n::Int, k::Int)
    out = UnitRange{Int}[]
    n == 0 && return out
    lo = 1
    for i in 1:k
        hi = lo + cld(n - lo + 1, k - i + 1) - 1
        push!(out, lo:hi)
        lo = hi + 1
        lo > n && break
    end
    return out
end

# Fill destination row `d` with the source chunks that may contribute to it.
function _fillrow!(buf::Vector{T}, rows::Vector{Vector{T}}, d::Int, dcap::SphericalCap,
        srccaps::Vector, order::Vector{Int}, sortedlats::Vector{Float64},
        maxsrcradius::Float64, radius::Float64, refine, prefilter::Bool) where {T}
    empty!(buf)
    intersects = DilatedIntersects(radius)
    band = maxsrcradius + dcap.radius + radius
    n = length(order)
    # A threshold at or above pi reaches every cap, and the latitude bound
    # degenerates; skip the prefilter rather than trusting the search bounds.
    lo, hi = if !prefilter || !(band < Float64(pi))
        1, n
    else
        dlat = _caplatitude(dcap)
        # Widen by one ulp on each side so an exactly-equal boundary case is
        # never dropped by rounding in the comparison.
        (searchsortedfirst(sortedlats, prevfloat(dlat - band)),
         searchsortedlast(sortedlats, nextfloat(dlat + band)))
    end
    @inbounds for j in lo:hi
        s = order[j]
        intersects(dcap, srccaps[s]) || continue
        (refine === nothing || refine(d, s)::Bool) || continue
        push!(buf, T(s))
    end
    sort!(buf)
    rows[d] = copy(buf)
    return nothing
end

# Counting-sort transpose of the destination-major CSR into source-major CSR.
function _transpose(dstoff::Vector{Int}, srcof::Vector{T}, nsrc::Int, ndst::Int,
        ::Type{T}) where {T}
    # Bin each edge's count one slot past its source, so the prefix sum turns
    # `srcoff[s]` into "1 + edges belonging to sources before s" directly.
    srcoff = zeros(Int, nsrc + 1)
    @inbounds for k in eachindex(srcof)
        srcoff[srcof[k]+1] += 1
    end
    srcoff[1] = 1
    @inbounds for s in 2:(nsrc+1)
        srcoff[s] += srcoff[s-1]
    end
    dstof = Vector{T}(undef, length(srcof))
    cursor = copy(srcoff)
    # Walking destinations in ascending order leaves every source row ascending.
    @inbounds for d in 1:ndst
        for k in dstoff[d]:(dstoff[d+1]-1)
            s = srcof[k]
            dstof[cursor[s]] = T(d)
            cursor[s] += 1
        end
    end
    return srcoff, dstof
end
