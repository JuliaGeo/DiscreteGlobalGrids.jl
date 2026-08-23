# The chunk dependency graph: a materialized, conservative bipartite relation
# between destination chunks and the source chunks that may contribute to them.
#
# `discovery.jl` answers "which source chunks does this destination chunk need?"
# one query at a time. That is the right shape for a lazy read, but a scheduler
# needs the whole relation at once: refcounts need every consumer of a source
# chunk, prefetch needs lookahead, and an affinity order needs both directions.
# This file materializes the relation once, by asking `discovery.jl` the same
# question for every destination chunk at once. Building it from the chunk caps
# directly would be faster and wrong: a space's own index need not test the caps
# `chunkextents` reports, and a graph that crosses the executor's relation
# instead of dominating it retires sources that are still going to be demanded.

# --------------------------------------------------------------------------
# Identity
# --------------------------------------------------------------------------

"""
    SpaceStamp

A cheap, serializable snapshot of the parts of a [`RegridSpace`](@ref) that a
dependency relation is derived from: the name of the space's type, its cell and
chunk counts, and a digest over its covering chunk caps.

Two spaces with different stamps cannot produce the same relation, so a stamp
mismatch is a *proof* that a graph does not belong to a pair of spaces. The
converse does not hold. See [`spacestamp`](@ref) for exactly what equal stamps
do and do not establish.

The fields are a `Symbol`, two `Int64`s and a `UInt64`, so a stamp survives
`Serialization` unchanged and carries no space type parameters.
"""
struct SpaceStamp
    tag::Symbol
    ncells::Int64
    nchunks::Int64
    digest::UInt64
end

Base.show(io::IO, s::SpaceStamp) =
    print(io, "SpaceStamp(", s.tag, ", ", s.ncells, " cells, ", s.nchunks,
        " chunks, 0x", string(s.digest; base = 16, pad = 16), ")")

const _EMPTY_STAMP = SpaceStamp(Symbol(""), 0, 0, UInt64(0))

# Domain-separating seed, so a cap digest cannot collide with an unrelated hash
# that happens to run over the same Float64s.
const _STAMP_SEED = 0x63_68_75_6e_6b_67_72_66 % UInt

"""
    spacestamp(space::RegridSpace) -> SpaceStamp

Fingerprint `space` for dependency-graph identity.

The digest covers the *name* of the space's type, [`ncells`](@ref),
[`nchunks`](@ref), and every covering cap [`chunkextents`](@ref) reports, in
chunk order. Cost is one `chunkextents` call plus `O(nchunks)` hashing; on a
`DGGSpace` the extents are a stored field, on a `RasterGrid` they are `nchunks`
rectangle-to-cap conversions.

# What a stamp match does and does not establish

**It proves nothing on its own.** Equal stamps are *evidence* that two spaces
are interchangeable for graph purposes, and a mismatch is proof that they are
not. Specifically, equal stamps still permit:

- A 64-bit digest collision. Unlikely, and silent when it happens.
- Two spaces of the same type *name* with different type **parameters** — a
  different grid system, element type or index type behind the same
  `DGGSpace`/`RasterGrid` name. Only the name is stamped, because rendering a
  full parametric type is [measured] 0.79 ms on a `RasterGrid`, which is more
  than building the relation. In practice a parameter that changes the relation
  changes a cap or a count as well, but the type parameters themselves are not
  checked.
- Two spaces with identical type name, counts and chunk caps but different
  **cell geometry**. The destination side of the relation is a function of the caps
  alone, so this is harmless there; on the source side it is not, because the
  relation comes from [`chunkindex`](@ref) — a native hierarchy that need not
  test the caps `chunkextents` reports (that divergence is the whole reason PR
  #69 exists). Two source spaces with equal caps and different hierarchies can
  therefore produce different relations and identical stamps.
- In-place mutation of a space after the graph was built. A stamp is a snapshot,
  not a live binding.

What it *does* catch, which is the failure this exists to prevent, is a graph
built against one pair of spaces being handed to a caller working on a different
pair: a different destination grid, a different source resolution, a different
chunking of either side, a swapped source and destination, or a space of an
entirely different type. Every one of those moves a count, the type tag, or a
cap.

Nothing cheaper is sound, and nothing sound is cheap: the only exact check is to
rebuild the relation and compare it, which is the work the identity exists to
avoid.
"""
spacestamp(space::RegridSpace) = _spacestamp(space, chunkextents(space))

function _spacestamp(space, caps::AbstractVector{<:SphericalCap})
    # `nameof`, not the full parametric type: rendering a `RasterGrid`'s type to
    # a string costs [measured] 0.79 ms — more than building the whole relation
    # on a small case — because `show(::Type)` searches every module for an
    # alias. The parameters that could change the relation change a cap or a
    # count anyway.
    tag = nameof(typeof(space))
    h = hash(tag, _STAMP_SEED)
    nc = Int64(ncells(space))
    h = hash(nc, h)
    h = hash(length(caps), h)
    @inbounds for c in caps
        p = c.point
        h = hash(c.radius, hash(p[3], hash(p[2], hash(p[1], h))))
    end
    return SpaceStamp(tag, nc, Int64(length(caps)), h)
end

"""
    DependencyIdentity

Everything a [`ChunkDependencyGraph`](@ref) must agree with a caller about
before the caller may reuse it: the two space stamps, the support radius the
relation was built at, and a tag naming the narrow phase that was applied.

All four fields are serializable and none of them holds a space, a closure or a
type parameter, so an identity can be written to disk beside a graph and
compared after a restart.

See [`dependency_identity`](@ref), [`validate_dependencies`](@ref).
"""
struct DependencyIdentity
    dst::SpaceStamp
    src::SpaceStamp
    radius::Float64
    narrow::Symbol
end

DependencyIdentity() = DependencyIdentity(_EMPTY_STAMP, _EMPTY_STAMP, 0.0, :none)

Base.show(io::IO, id::DependencyIdentity) =
    print(io, "DependencyIdentity(dst ", id.dst.nchunks, " chunks, src ",
        id.src.nchunks, " chunks, radius ", id.radius, ", narrow ",
        repr(id.narrow), ")")

"""
    dependency_identity(dst_space, src_space; radius = 0.0, narrow = :none)
        -> DependencyIdentity
    dependency_identity(graph::ChunkDependencyGraph) -> DependencyIdentity

The identity a graph over these spaces, at this radius, with this narrow phase
*would* carry — or the identity a graph does carry.

The two-space form is the one a caller uses to decide whether a graph it already
has is the graph it wants; it costs two [`spacestamp`](@ref) calls and builds no
relation. Callers that check repeatedly should compute the stamp once and keep
it, rather than restamping per column.
"""
dependency_identity(dst_space::RegridSpace, src_space::RegridSpace;
        radius::Real = 0.0, narrow::Symbol = :none) =
    DependencyIdentity(spacestamp(dst_space), spacestamp(src_space),
        _checkedradius(radius), narrow)

_checkedradius(radius::Real) = let r = Float64(radius)
    (isfinite(r) && r >= 0) || throw(ArgumentError(
        "radius must be finite and non-negative, got $radius"))
    r
end

"""
    GlobalRegridding.UNNAMED_NARROW

The narrow-phase tag (`:unnamed`) stamped on a graph built with a `refine` its
caller did not name. [`validate_dependencies`](@ref) rejects it unconditionally,
including against itself: an anonymous closure has no identity, so a relation
narrowed by one cannot be certified as the relation any later caller wants. Pass
`narrow` alongside `refine` to make such a graph reusable.
"""
const UNNAMED_NARROW = :unnamed

# A narrow-phase tag is a *claim about the relation*, so the tag and the closure
# must agree: `:none` means "every candidate the cap test accepted is here", and
# a name means "these pairs, and only these, survived the phase that name refers
# to". Either half lying makes reuse unsafe, so both mismatches are errors. The
# unsupplied tag is `nothing`, not `:none`, precisely so that claiming `:none`
# over a `refine` is distinguishable from not claiming anything.
_narrowtag(refine, ::Nothing) = refine === nothing ? :none : UNNAMED_NARROW

function _narrowtag(refine, narrow::Symbol)
    if refine === nothing
        narrow === :none || throw(ArgumentError(
            "narrow = $(repr(narrow)) names a narrow phase, but no `refine` was " *
            "given; a graph that applied no narrow phase must be tagged :none"))
        return :none
    end
    narrow === :none && throw(ArgumentError(
        "narrow = :none claims no narrow phase was applied, but a `refine` was " *
        "given; tag it with the name of the phase it implements, or leave " *
        "`narrow` unset to record it as $(repr(UNNAMED_NARROW))"))
    return narrow
end

"""
    ChunkDependencyGraph{T} <: Graphs.AbstractGraph{T}

A conservative bipartite dependency relation between the chunks of a destination
space and the chunks of a source space, materialized as bidirectional CSR.

An edge `(s, d)` means source chunk `s` *may* contribute to destination chunk
`d`. The relation is a **superset**: it never omits a contributing pair, but it
may include pairs whose geometries do not actually overlap, because it is
computed from covering spherical caps rather than cell geometry. Treat an edge
as "must be loaded", never as "has nonzero weight".

What makes it usable as a *refcount* is a second, sharper property: the relation
is the one [`candidatechunks!`](@ref) answers on the source space's own chunk
index, so a lazy read can never demand a source chunk the graph did not predict.
See [`chunk_dependency_graph`](@ref).

# Vertex numbering

The graph implements the Graphs.jl `AbstractGraph` interface over a single
vertex range, so that ecosystem algorithms work directly on it:

| vertices | meaning |
|---|---|
| `1 : nsrc` | source chunk `v` |
| `nsrc+1 : nsrc+ndst` | destination **row** `v - nsrc` |

Source vertices are numbered first so that a destination's neighbours are
already source *chunk* numbers with no offset arithmetic. Use
[`sourcesof`](@ref) and [`consumersof`](@ref) for the chunk-numbered views, and
[`srcvertex`](@ref)/[`dstvertex`](@ref)/[`srcchunk`](@ref)/[`dstchunk`](@ref) to
convert. The graph is undirected: every edge has one source and one destination
endpoint, so `has_edge` and `edges` are role-agnostic.

Destination *rows* and destination *chunks* coincide on a graph built over a
whole destination space, and diverge on a row view; see [`restrict`](@ref) and
[`globaldestination`](@ref).

The spaces themselves are deliberately *not* stored. The graph carries a
[`DependencyIdentity`](@ref) instead, which keeps it small, serializable, and
free of the space type parameters.

# Chunk extents

A graph built from spaces also retains the two covering-cap vectors it was
derived from — [`chunkextents`](@ref) of each side — because they are its own
inputs and it would otherwise throw them away for every consumer to rebuild.
[`destinationextents`](@ref), [`sourceextents`](@ref) and their singular forms
read them; [`hasextents`](@ref) says whether a graph carries them at all, which
a graph assembled from bare CSR arrays does not.

This is where per-chunk cap metadata *lives*. A consumer that needs a chunk's
extent — the lazy executor weighs a wave's blocks by how much of each source
chunk a destination tile can reach — reads it off the relation both sides
already agree on, rather than keeping a private copy per array, per worker or
per column. The vectors are shared by reference: on a `DGGSpace` they are the
space's own stored field, and a [`restrict`](@ref) row view shares its parent's.

# Identity

A graph is only reusable against the inputs it was built from, and it carries
enough to say so: both space stamps, the support radius, and the narrow-phase
tag. See [`dependency_identity`](@ref) for the record and
[`validate_dependencies`](@ref) for the check. The radius is part of that
identity because the relation depends on the regridding method only through its
support radius: two methods with the same support radius share a graph, and a
graph built at one radius is not valid at a larger one.

See [`chunk_dependency_graph`](@ref).
"""
struct ChunkDependencyGraph{T<:Integer} <: Graphs.AbstractGraph{T}
    id::DependencyIdentity
    nsrc::Int
    ndst::Int
    # Destination-major CSR, indexed by GLOBAL destination chunk: sources of
    # destination chunk `c` are `srcof[dstoff[c]:dstoff[c+1]-1]`, ascending
    # source chunk numbers. A row view shares both arrays with its parent
    # unchanged; that sharing is the whole point of `restrict`.
    dstoff::Vector{Int}
    srcof::Vector{T}
    # Local row `d` -> global destination chunk, ascending. EMPTY means the
    # identity map, which is the case for every graph built over a whole
    # destination space; `_row` is the one place that distinction is read.
    dstrows::Vector{Int}
    # Source-major CSR over THIS graph's rows: consumers of source chunk `s` are
    # `dstof[srcoff[s]:srcoff[s+1]-1]`, ascending LOCAL destination rows. Always
    # private, because a row view's refcounts are not its parent's.
    srcoff::Vector{Int}
    dstof::Vector{T}
    # Whether the two cap vectors below are populated. An explicit flag, not an
    # `isempty` test: a zero-chunk side has no caps and still carries them, and
    # a subspace view's `dstcaps` is its parent's, whose length matches neither
    # of this graph's own counts.
    extents::Bool
    # The covering caps of each side, in chunk order, as `chunkextents` reported
    # them while the relation was built. `dstcaps` is indexed by GLOBAL
    # destination chunk, exactly like `dstoff`, so a row view shares both by
    # reference and `_row` indexes both.
    dstcaps::Vector{Cap}
    srccaps::Vector{Cap}
end

"""
    ChunkDependencyGraph(id::DependencyIdentity, dstoff, srcof, srcoff, dstof;
                         dstcaps = Cap[], srccaps = Cap[])

Assemble a graph over a whole destination space from prebuilt bidirectional CSR
arrays. `dstoff` and `srcoff` are offset vectors of length `ndst+1` and
`nsrc+1`. For callers that build the relation themselves; production goes
through [`chunk_dependency_graph`](@ref).

`dstcaps`/`srccaps` are the two sides' covering chunk extents. Supply both or
neither: a graph assembled without them answers `hasextents(g) === false` and
cannot serve a consumer that needs per-chunk geometry, such as the lazy
executor's wave costing.
"""
function ChunkDependencyGraph(id::DependencyIdentity, dstoff::Vector{Int},
        srcof::Vector{T}, srcoff::Vector{Int}, dstof::Vector{T};
        dstcaps::AbstractVector{<:SphericalCap} = Cap[],
        srccaps::AbstractVector{<:SphericalCap} = Cap[]) where {T<:Integer}
    length(srcof) == length(dstof) || throw(ArgumentError(
        "the two CSR directions hold $(length(srcof)) and $(length(dstof)) edges"))
    nsrc, ndst = length(srcoff) - 1, length(dstoff) - 1
    have = !(isempty(dstcaps) && isempty(srccaps))
    have || return ChunkDependencyGraph{T}(id, nsrc, ndst, dstoff, srcof, Int[],
        srcoff, dstof, false, Cap[], Cap[])
    return ChunkDependencyGraph{T}(id, nsrc, ndst, dstoff, srcof, Int[],
        srcoff, dstof, true, _graphcaps(dstcaps, ndst, "destination"),
        _graphcaps(srccaps, nsrc, "source"))
end

# Extents are stored as a plain `Vector{Cap}` so the graph stays free of space
# type parameters. `convert` aliases when the vector already has that type,
# which is the case on every shipped space, so nothing is copied.
function _graphcaps(caps::AbstractVector{<:SphericalCap}, n::Int, role::AbstractString)
    length(caps) == n || throw(ArgumentError(
        "$(length(caps)) $role chunk extents for $n $role chunks; supply both " *
        "sides' extents or neither"))
    return convert(Vector{Cap}, caps)
end

function Base.show(io::IO, g::ChunkDependencyGraph)
    print(io, "ChunkDependencyGraph(", g.nsrc, " source × ", g.ndst,
        " destination ", isrestricted(g) ? "rows" : "chunks", ", ",
        length(g.dstof), " edges, radius ", g.id.radius)
    g.id.narrow === :none || print(io, ", narrow ", repr(g.id.narrow))
    isrestricted(g) && print(io, ", row view of ", g.id.dst.nchunks)
    print(io, ")")
end

"""
    nsourcechunks(g::ChunkDependencyGraph) -> Int
    ndestinationchunks(g::ChunkDependencyGraph) -> Int

Return the size of each side of the bipartition. On a row view
`ndestinationchunks` is the number of rows the view holds, not the destination
space's chunk count; that one is `dependency_identity(g).dst.nchunks`.
"""
nsourcechunks(g::ChunkDependencyGraph) = g.nsrc
ndestinationchunks(g::ChunkDependencyGraph) = g.ndst

"""
    dependency_identity(g::ChunkDependencyGraph) -> DependencyIdentity

Return the graph's identity record: both space stamps, the radius, and the
narrow-phase tag.
"""
dependency_identity(g::ChunkDependencyGraph) = g.id

# --------------------------------------------------------------------------
# Chunk extents: the relation's own inputs, kept
# --------------------------------------------------------------------------

"""
    hasextents(g::ChunkDependencyGraph) -> Bool

Return whether `g` carries the covering chunk extents of both sides.

True for every graph built from spaces — [`chunk_dependency_graph`](@ref) and
plan construction — because those are the relation's own inputs. False only for
a graph assembled from bare CSR arrays, which never had them.
"""
hasextents(g::ChunkDependencyGraph) = g.extents

"""
    destinationextents(g::ChunkDependencyGraph) -> Vector{Cap}
    sourceextents(g::ChunkDependencyGraph) -> Vector{Cap}

Return the covering spherical cap of every chunk of each side, in chunk order,
as [`chunkextents`](@ref) reported it when the relation was built.

The result is the graph's own vector, shared by reference and **not** to be
mutated. `destinationextents` is indexed by the *destination space's* chunk
number, so on a row view it is the parent space's whole vector; index it with
[`globaldestination`](@ref), or use [`destinationextent`](@ref), which does.
"""
function destinationextents(g::ChunkDependencyGraph)
    _requireextents(g)
    return g.dstcaps
end

function sourceextents(g::ChunkDependencyGraph)
    _requireextents(g)
    return g.srccaps
end

"""
    destinationextent(g::ChunkDependencyGraph, d) -> Cap
    sourceextent(g::ChunkDependencyGraph, s) -> Cap

Return one chunk's covering cap. `d` is a destination **row** of `g`, as
everywhere else on this type; `s` is a source chunk number.
"""
function destinationextent(g::ChunkDependencyGraph, d::Integer)
    _requireextents(g)
    i = Int(d)
    1 <= i <= g.ndst || throw(BoundsError(g, d))
    return @inbounds g.dstcaps[_row(g, i)]
end

function sourceextent(g::ChunkDependencyGraph, s::Integer)
    _requireextents(g)
    i = Int(s)
    1 <= i <= g.nsrc || throw(BoundsError(g, s))
    return @inbounds g.srccaps[i]
end

@inline _requireextents(g::ChunkDependencyGraph) =
    g.extents || throw(ArgumentError(
        "$g carries no chunk extents. It was assembled from bare CSR arrays " *
        "rather than built from spaces, so the per-chunk geometry this asks " *
        "for was never part of it; build the relation with " *
        "`chunk_dependency_graph` or a plan's `dependencies` keyword."))

# The one place the "empty means identity" encoding of `dstrows` is read.
@inline _row(g::ChunkDependencyGraph, d::Int) =
    isempty(g.dstrows) ? d : @inbounds(g.dstrows[d])

"""
    dependency_radius(g::ChunkDependencyGraph) -> Float64

Return the support radius, in radians, the relation was built at. A graph is
valid for any method whose [`support_radius`](@ref) is at most this value.
"""
dependency_radius(g::ChunkDependencyGraph) = g.id.radius

"""
    narrowphase(g::ChunkDependencyGraph) -> Symbol

Return the serializable tag naming the narrow phase the relation was built with:
`:none` when no `refine` was applied, [`UNNAMED_NARROW`](@ref) when one was
applied without a name, and the caller's own tag otherwise.
"""
narrowphase(g::ChunkDependencyGraph) = g.id.narrow

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
    sourcesof(g, d) -> AbstractVector{<:Integer}

Return the ascending source chunk numbers that may contribute to destination row
`d`. The result is a read-only view into the graph; do not mutate it.

This is the executor's query: the set of source chunks that must be resident to
compute this destination chunk. On a row view `d` is the view's own row index
and [`globaldestination`](@ref) maps it back to the destination space's chunk
number; on a whole-space graph the two coincide.
"""
function sourcesof(g::ChunkDependencyGraph, d::Integer)
    i = Int(d)
    1 <= i <= g.ndst || throw(BoundsError(g, d))
    r = _row(g, i)
    return view(g.srcof, g.dstoff[r]:(g.dstoff[r+1]-1))
end

"""
    consumersof(g, srcchunk) -> AbstractVector{<:Integer}

Return the ascending destination rows that may need source chunk `srcchunk`. The
result is a read-only view into the graph; do not mutate it.

`length(consumersof(g, s))` is the initial refcount for source chunk `s` over
exactly the rows this graph holds — so a row view's refcounts are the row view's
own, not its parent's. See the eviction notes in
[`chunk_dependency_graph`](@ref) before using them as one.
"""
function consumersof(g::ChunkDependencyGraph, s::Integer)
    i = Int(s)
    1 <= i <= g.nsrc || throw(BoundsError(g, s))
    return view(g.dstof, g.srcoff[i]:(g.srcoff[i+1]-1))
end

"""
    sourcedegree(g, d) -> Int
    consumerdegree(g, srcchunk) -> Int

Return adjacency counts without materializing the rows. `sourcedegree` is a
cheap per-destination cost and I/O estimate; `consumerdegree` is a source
chunk's consumer count over the rows this graph holds.
"""
sourcedegree(g::ChunkDependencyGraph, d::Integer) =
    (i = Int(d); 1 <= i <= g.ndst || throw(BoundsError(g, d));
     r = _row(g, i); g.dstoff[r+1] - g.dstoff[r])
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
# `srcof` is the parent's whole edge array on a row view, so the edge count is
# read from the private source-major side, which holds exactly this graph's rows.
Graphs.ne(g::ChunkDependencyGraph) = length(g.dstof)
Graphs.vertices(g::ChunkDependencyGraph{T}) where {T} = Base.OneTo(Graphs.nv(g))
Graphs.has_vertex(g::ChunkDependencyGraph, v::Integer) = 1 <= Int(v) <= g.nsrc + g.ndst

Base.zero(::Type{ChunkDependencyGraph{T}}) where {T} =
    ChunkDependencyGraph{T}(DependencyIdentity(), 0, 0, [1], T[], Int[], [1], T[],
        true, Cap[], Cap[])
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
    chunk_dependency_graph(dst_space, src_space; radius = 0.0)
        -> ChunkDependencyGraph

Build the conservative bipartite dependency relation between the chunks of
`dst_space` and those of `src_space`, from destination chunk extents and the
source space's chunk index alone. No cell geometry is built and no data is read.

Argument order matches the rest of `discovery.jl`: destination first.

This is the low-level builder, for a caller that knows the radius it wants and
owns the result itself. A caller that has a plan should not use it: a
[`ChunkedPlan`](@ref) owns exactly one relation, built or validated when the
plan was constructed, and [`dependencies`](@ref) is how to read it. There is
deliberately no `plan` method here and no narrow-phase keyword: a narrow phase
is an argument to plan construction ([`plan_regrid`](@ref)) and to nothing else,
so that a relation a plan hands out can never have been thinned behind the
plan's back. See [`dependencies`](@ref).

# Keywords

- `radius`: the method's angular support radius in radians. Two caps are
  connected when their centres are within the sum of their radii plus this. Must
  be finite and non-negative.

# Method

The destination side is reduced to its covering spherical caps via
[`chunkextents`](@ref). Each destination cap is then handed to
[`candidatechunks!`](@ref) on [`chunkindex`](@ref)`(src_space)` — *the same
query, on the same index, that a lazy read issues for that destination*.

That identity is the point, and it is a correctness property rather than a
convenience. A space's chunk index need not test the caps `chunkextents`
reports: a hierarchy derives node extents its own way, and a covering cap
derived from a node rectangle is neither contained in nor containing the cap
derived from the chunk's own boundary. Building the graph from `chunkextents`
directly therefore produced a relation that *crosses* the executor's rather than
dominating it — measured on the Copernicus DEM × IGeo7 pair as 71 demanded
pairs the graph did not hold, against 437 it held and no read ever asked for.
Refcounts derived from such a graph retire a source that is still going to be
demanded, and the demand then reloads it. Querying the index closes that gap by
construction, for every space, with no per-space invariant to maintain.

Destination rows are built in parallel over blocks of destination chunks and
written by index, so the result is identical regardless of thread count. The
source index is built once and shared; `candidatechunks!` is already called
concurrently from the lazy executor's workers and must be safe for that.

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
chunk_dependency_graph(dst_space::RegridSpace, src_space::RegridSpace;
    radius::Real = 0.0) =
    _builddependencies(dst_space, src_space, radius, nothing, nothing)

# The narrow phase enters here and nowhere a caller can reach with a plan in
# hand: `plans.jl` calls this once, during construction, and stores the result.
# `refine(dstchunk, srcchunk) -> Bool` returning `false` asserts the pair cannot
# contribute at this radius and drops the edge; `narrow` is the Symbol that
# names the phase in the graph's identity. See `plan_regrid`.
function _builddependencies(dst_space::RegridSpace, src_space::RegridSpace,
        radius::Real, refine, narrow)
    r = _checkedradius(radius)
    tag = _narrowtag(refine, narrow)
    # Both cap vectors are the graph's own input — the destination's directly,
    # the source's through its stamp — so take each once and keep it, rather
    # than paying `chunkextents` twice here and again in every consumer.
    dstcaps = chunkextents(dst_space)
    srccaps = chunkextents(src_space)
    id = DependencyIdentity(_spacestamp(dst_space, dstcaps),
        _spacestamp(src_space, srccaps), r, tag)
    return _chunkgraph(id, dstcaps, chunkindex(src_space), srccaps, r, refine)
end

function _chunkgraph(id::DependencyIdentity,
        dstcaps::AbstractVector{<:SphericalCap}, srcindex,
        srccaps::AbstractVector{<:SphericalCap}, radius::Float64, refine)
    nsrc = length(srccaps)
    T = Int32
    ndst = length(dstcaps)
    (nsrc <= typemax(T) && ndst <= typemax(T)) || throw(ArgumentError(
        "chunk counts $nsrc and $ndst exceed the $(T) vertex numbering"))

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
                buf = Int[]
                for d in b
                    _fillrow!(buf, rows, d, dstcaps[d], srcindex, radius, refine)
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
    return ChunkDependencyGraph{T}(id, nsrc, ndst, dstoff, srcof, Int[],
        srcoff, dstof, true, convert(Vector{Cap}, dstcaps),
        convert(Vector{Cap}, srccaps))
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
# `candidatechunks!` returns them ascending and duplicate-free, which is what the
# CSR rows and every `insorted` lookup on them rely on.
function _fillrow!(buf::Vector{Int}, rows::Vector{Vector{T}}, d::Int,
        dcap::SphericalCap, srcindex, radius::Float64, refine) where {T}
    candidatechunks!(buf, srcindex, convert(Cap, dcap); radius)
    refine === nothing || filter!(s -> refine(d, s)::Bool, buf)
    rows[d] = T.(buf)
    return nothing
end

_transpose(dstoff::Vector{Int}, srcof::Vector{T}, nsrc::Int, ndst::Int,
    ::Type{T}) where {T} = _transpose(dstoff, srcof, Base.OneTo(ndst), nsrc, T)

# Counting-sort transpose of a *selection* of destination-major rows into a
# source-major CSR whose destination numbers are the selection's own positions.
# `rows` holds global destination chunk numbers in ascending order; passing
# `OneTo(ndst)` transposes the whole relation, which is the whole-space case.
function _transpose(dstoff::Vector{Int}, srcof::Vector{T},
        rows::AbstractVector{<:Integer}, nsrc::Int, ::Type{T}) where {T}
    # Bin each edge's count one slot past its source, so the prefix sum turns
    # `srcoff[s]` into "1 + edges belonging to sources before s" directly.
    srcoff = zeros(Int, nsrc + 1)
    @inbounds for r in rows
        for k in dstoff[r]:(dstoff[r+1]-1)
            srcoff[srcof[k]+1] += 1
        end
    end
    srcoff[1] = 1
    @inbounds for s in 2:(nsrc+1)
        srcoff[s] += srcoff[s-1]
    end
    dstof = Vector{T}(undef, srcoff[end] - 1)
    cursor = copy(srcoff)
    # Walking destinations in ascending order leaves every source row ascending.
    @inbounds for (d, r) in enumerate(rows)
        for k in dstoff[r]:(dstoff[r+1]-1)
            s = srcof[k]
            dstof[cursor[s]] = T(d)
            cursor[s] += 1
        end
    end
    return srcoff, dstof
end

# --------------------------------------------------------------------------
# Row views
# --------------------------------------------------------------------------

"""
    isrestricted(g::ChunkDependencyGraph) -> Bool

Return whether `g` is a row view over part of a destination space rather than a
graph over the whole of one. See [`restrict`](@ref).
"""
isrestricted(g::ChunkDependencyGraph) = !isempty(g.dstrows)

"""
    globaldestinations(g::ChunkDependencyGraph) -> AbstractVector{Int}

Return the destination space's chunk number for each of `g`'s rows, ascending.
On a whole-space graph this is `1:ndestinationchunks(g)`.

This is what makes a row view usable for refinement: the view renumbers its rows
`1:k` for Graphs.jl's sake, but never loses which chunk of the destination space
each row is.

On a [`subspace_dependencies`](@ref) view the numbers are the *parent* relation's
destination chunks, not the sub-space's own — the sub-space's chunks are its
rows, `1:k`. That is the provenance a caller wants from a re-stamped view: which
chunks of the whole destination these rows came from.
"""
globaldestinations(g::ChunkDependencyGraph) =
    isrestricted(g) ? g.dstrows : Base.OneTo(g.ndst)

"""
    globaldestination(g::ChunkDependencyGraph, d) -> Int

Return the destination space's chunk number of `g`'s row `d`.
"""
function globaldestination(g::ChunkDependencyGraph, d::Integer)
    i = Int(d)
    1 <= i <= g.ndst || throw(BoundsError(g, d))
    return _row(g, i)
end

"""
    localdestination(g::ChunkDependencyGraph, chunk) -> Union{Int,Nothing}

Return the row of `g` that holds destination chunk `chunk`, or `nothing` if the
view does not hold it. `O(log k)`: rows are ascending in the destination chunk
number.
"""
function localdestination(g::ChunkDependencyGraph, chunk::Integer)
    c = Int(chunk)
    if !isrestricted(g)
        return 1 <= c <= g.ndst ? c : nothing
    end
    i = searchsortedfirst(g.dstrows, c)
    return (i <= length(g.dstrows) && g.dstrows[i] == c) ? i : nothing
end

"""
    restrict(g::ChunkDependencyGraph, destinations) -> ChunkDependencyGraph

Return the row view of `g` over `destinations`, an ascending, duplicate-free
collection of rows of `g`.

The result is a `ChunkDependencyGraph` like any other — same accessors, same
bidirectional CSR, same Graphs.jl interface — over `length(destinations)`
destination rows and the same source side. Its rows are renumbered `1:k`;
[`globaldestination`](@ref) maps them back, and its
[`dependency_identity`](@ref) still stamps the **whole** destination space, so a
view knows what it is a view *of*.

That last point decides who may adopt one, and [`subspace_dependencies`](@ref)
is the way to change it. A view is the relation of *these rows
of that space*, so [`validate_dependencies`](@ref) certifies it only against
that space, with `destinations` naming the rows. It is **not** the relation of a
smaller space built over the same cells: a plan whose destination is its own
one-chunk grid has a different [`spacestamp`](@ref) and is refused. To hand a
view to a plan over such a sub-space, re-stamp it with
[`subspace_dependencies`](@ref), which checks cap for cap that the sub-space's
chunks are the selected chunks and stamps the destination half against it.

# What is shared and what is not

The destination-major direction — the offsets and the whole edge array — is the
parent's, by reference, not a copy. Only the source-major direction is rebuilt,
because a view's refcounts are genuinely different numbers: `consumersof` on a
view must count only the rows the view holds, or a refcount derived from it
retires nothing. That rebuild is a counting-sort over the selected rows' edges,
`O(selected edges + nsourcechunks)`, with no spatial index, no cap test and no
query of either space.

That is the saving, and it is [measured] 19–104× for a single destination and
33–370× for a column of a sixteenth of the destination space, across the whole
`benchmark/chunk_graph_gates.jl` matrix. On the production Copernicus × IGeo7
pair: [measured] 26× for one destination, 80× and 16× fewer bytes for a
4136-chunk column. What a rebuild pays that a view does not is the destination
caps, the source [`chunkindex`](@ref), and one `candidatechunks!` query per row —
plus, on the generic path, a fresh packed R-tree over every source chunk extent,
though neither shipped native space takes that path.

The `O(nsourcechunks)` term is not free. On a source with many chunks and a view
with very few rows, the transpose's offset array can cost more than re-querying
one row would: [measured] on the production pair, a rebuild of a *single*
destination's rows without its identity and CSR is 24 µs against `restrict`'s
39 µs. Restrict a column, not a row, and the term disappears into the edges.

# Ordering

`destinations` must be strictly ascending. A schedule is a separate permutation
applied by whoever walks the rows, exactly as it is for a whole-space graph; the
graph is never built in walk order. Requiring ascending rows is what keeps
`globaldestinations` searchable and both CSR directions sorted.

# Example

```julia
graph = chunk_dependency_graph(dst_space, src_space)
column = restrict(graph, 17:24)          # eight destination chunks
sourcesof(column, 1) == sourcesof(graph, 17)
globaldestination(column, 1) == 17
```
"""
function restrict(g::ChunkDependencyGraph{T}, destinations) where {T}
    rows = _restrictionrows(g, destinations)
    srcoff, dstof = _transpose(g.dstoff, g.srcof, rows, g.nsrc, T)
    return ChunkDependencyGraph{T}(g.id, g.nsrc, length(rows), g.dstoff, g.srcof,
        rows, srcoff, dstof, g.extents, g.dstcaps, g.srccaps)
end

"""
    subspace_dependencies(g::ChunkDependencyGraph, subspace, destinations)
        -> ChunkDependencyGraph

Return [`restrict`](@ref)`(g, destinations)` re-stamped as the relation of
`subspace` — a destination space whose chunks *are* `g`'s destination chunks
`destinations`, in that order.

This is what lets a plan over a piece of a destination space adopt the relation
built over the whole of it, instead of rebuilding its own. The result is a whole
destination-space relation *for `subspace`*: `validate_dependencies(view,
subspace, src_space; radius)` certifies it with no `destinations` argument,
because for `subspace` there are no other rows.

# Why this is sound, and what is checked

The destination half of the relation is a function of the destination chunk caps
alone: a row is one [`candidatechunks!`](@ref) query of the source index against
one cap. So if `subspace`'s chunk `k` has *the same cap* as `g`'s destination
chunk `destinations[k]`, the row `subspace` would produce for `k` is the row `g`
already holds — the same query, the same index, the same cap, the same answer.
The source half is untouched: it is literally the parent's, and its stamp is
carried over unchanged.

That equality is therefore checked, cap for cap, and a mismatch is an
`ArgumentError` rather than a re-stamp. Two caps are compared exactly, because
they are two computations of the same covering cap and only a difference in how
the two spaces derive it can separate them; a re-stamp on nearly-equal caps
would be a relation for caps that are not these. The count must match too:
`nchunks(subspace) == length(destinations)`.

What is **not** checked, and cannot be, is that `subspace`'s chunk `k` holds the
same *cells* as the parent's chunk `destinations[k]`. Equal caps do not imply
equal cell sets. That is the same hole [`spacestamp`](@ref) documents, in the
one place where it matters most, and it is why the destination stamp of the
result is computed from `subspace` itself rather than asserted: a caller that
hands over an unrelated space with coincidentally equal caps gets a certified
relation for a destination it did not mean. Passing a space that is genuinely a
sub-space of the graph's destination is the caller's obligation.

`g` must carry chunk extents ([`hasextents`](@ref)); a graph assembled from bare
CSR arrays has nothing to compare against.

# Cost

One `chunkextents(subspace)` call, `length(destinations)` cap comparisons, one
`SpaceStamp` over the sub-space, and `restrict`'s `O(selected edges +
nsourcechunks)` transpose. No query of either space's chunk index, and no cap
test against the source.

# Example

```julia
graph  = chunk_dependency_graph(dst_space, src_space)
piece  = subspace_of(dst_space, 17:24)               # caller's own sub-space
view   = subspace_dependencies(graph, piece, 17:24)
sourcesof(view, 1) == sourcesof(graph, 17)
plan   = plan_regrid(data; to = piece, from = src_space, lazy = true,
                     dependencies = view)
```
"""
function subspace_dependencies(g::ChunkDependencyGraph{T}, subspace::RegridSpace,
        destinations) where {T}
    _requireextents(g)
    rows = _restrictionrows(g, destinations)
    subcaps = chunkextents(subspace)
    length(subcaps) == length(rows) || throw(ArgumentError(
        "$subspace has $(length(subcaps)) chunks, but $(length(rows)) " *
        "destination chunks of $g were selected; a sub-space relation must " *
        "cover the sub-space exactly"))
    @inbounds for k in eachindex(rows)
        subcaps[k] == g.dstcaps[rows[k]] || throw(ArgumentError(
            "chunk $k of $subspace has extent $(subcaps[k]), but destination " *
            "chunk $(rows[k]) of $g has $(g.dstcaps[rows[k]]); a sub-space may " *
            "only adopt rows whose destination cap it reproduces exactly, " *
            "because the cap is the whole of what those rows were derived from"))
    end
    id = DependencyIdentity(_spacestamp(subspace, subcaps), g.id.src,
        g.id.radius, g.id.narrow)
    srcoff, dstof = _transpose(g.dstoff, g.srcof, rows, g.nsrc, T)
    return ChunkDependencyGraph{T}(id, g.nsrc, length(rows), g.dstoff, g.srcof,
        rows, srcoff, dstof, true, g.dstcaps, g.srccaps)
end

# Translate a selection of this graph's rows into global destination chunk
# numbers, which are also indices into the shared destination-major offsets.
# Validation happens here, before anything is built, so an ill-formed selection
# never produces a graph.
function _restrictionrows(g::ChunkDependencyGraph, destinations)
    rows = Vector{Int}(undef, 0)
    sizehint!(rows, length(destinations))
    last = 0
    for d in destinations
        i = Int(d)
        1 <= i <= g.ndst || throw(ArgumentError(
            "destination $i is outside the 1:$(g.ndst) rows of $g"))
        i > last || throw(ArgumentError(
            "restrict needs strictly ascending destinations, got $i after $last; " *
            "a walk order is a separate permutation, not a row order"))
        last = i
        push!(rows, _row(g, i))
    end
    return rows
end

# --------------------------------------------------------------------------
# Validated reuse
# --------------------------------------------------------------------------

"""
    validate_dependencies(g::ChunkDependencyGraph, dst_space, src_space;
                          radius = 0.0, narrow = :none, destinations = nothing)
        -> g

Return `g` if it is a relation a caller regridding `src_space` onto `dst_space`
at `radius` may reuse, and throw an `ArgumentError` naming the mismatch
otherwise. This is the check a plan runs *before* adopting a graph it did not
build, so that an invalid reuse fails at construction rather than as a wrong
answer later.

Five things must agree.

1. **The destination space.** [`spacestamp`](@ref)`(dst_space)` must equal the
   stamp the graph carries.
2. **The source space.** Likewise. A swapped pair fails here, because the two
   stamps are compared in their own roles.
3. **The radius.** `radius` must be at most [`dependency_radius`](@ref)`(g)`. The
   relation grows monotonically with the radius, so a graph built wider is still
   a conservative superset; one built narrower is not, and is rejected.
4. **The narrow phase.** `narrow` must equal the graph's tag exactly, and
   [`UNNAMED_NARROW`](@ref) never matches. A caller expecting the full candidate
   relation must not silently receive one somebody else thinned, and a caller
   expecting a specific thinning must not silently receive the full one.
5. **The rows.** With `destinations === nothing`, `g` must hold one row per chunk
   of the destination space it stamps. That is every whole-space graph, and also
   a [`subspace_dependencies`](@ref) view, which stamps the sub-space it covers.
   Pass `destinations` — the ascending chunk numbers of the graph's own
   destination the caller intends to compute — to validate a [`restrict`](@ref)
   row view instead; they must be exactly the view's
   [`globaldestinations`](@ref).

What this does **not** establish is that `g` is *the* relation these spaces
produce: the stamps are fingerprints, not proofs. [`spacestamp`](@ref) documents
precisely what equal stamps leave open. The only exact check is to rebuild the
relation, which is the work this exists to avoid.

Cost is two `spacestamp` calls. A caller validating repeatedly should stamp once
and compare [`dependency_identity`](@ref) records directly.
"""
function validate_dependencies(g::ChunkDependencyGraph, dst_space::RegridSpace,
        src_space::RegridSpace; radius::Real = 0.0, narrow::Symbol = :none,
        destinations = nothing)
    r = _checkedradius(radius)
    _checkstamp(g, :destination, g.id.dst, spacestamp(dst_space))
    _checkstamp(g, :source, g.id.src, spacestamp(src_space))
    r <= g.id.radius || throw(ArgumentError(
        "$g was built at support radius $(g.id.radius) and cannot be reused at " *
        "the wider radius $r; the relation is only monotone the other way"))
    narrow === UNNAMED_NARROW && throw(ArgumentError(
        "$(repr(UNNAMED_NARROW)) is not a narrow phase a caller can ask for; it " *
        "is the tag for a `refine` nobody named"))
    g.id.narrow === UNNAMED_NARROW && throw(ArgumentError(
        "$g was built with an unnamed `refine`, so which pairs it dropped is not " *
        "recorded and it cannot be reused; rebuild it passing `narrow`"))
    g.id.narrow === narrow || throw(ArgumentError(
        "$g was built with narrow phase $(repr(g.id.narrow)), but the caller " *
        "expects $(repr(narrow))"))
    if destinations === nothing
        # A graph is whole for the space it stamps when it has a row per chunk
        # of it. That is every unrestricted graph; it is also a
        # `subspace_dependencies` view, whose stamp is the sub-space's own. A
        # plain row view has fewer rows than the space it stamps and lands in
        # the error, which is what refuses a fraction of a destination.
        g.ndst == Int(g.id.dst.nchunks) || throw(ArgumentError(
            "$g holds $(g.ndst) rows for a $(g.id.dst.nchunks)-chunk " *
            "destination; pass `destinations` to validate it as a row view, " *
            "validate its parent, or re-stamp it onto the space it covers " *
            "with `subspace_dependencies`"))
    else
        wanted = collect(Int, destinations)
        globaldestinations(g) == wanted || throw(ArgumentError(
            "$g holds destination chunks $(_summarize(globaldestinations(g))), " *
            "but the caller asked for $(_summarize(wanted))"))
    end
    return g
end

function _checkstamp(g, role::Symbol, have::SpaceStamp, want::SpaceStamp)
    have == want && return nothing
    reason = have.tag != want.tag ?
             "space type $(want.tag) against $(have.tag)" :
             have.nchunks != want.nchunks ?
             "$(want.nchunks) chunks against $(have.nchunks)" :
             have.ncells != want.ncells ?
             "$(want.ncells) cells against $(have.ncells)" :
             "the same type, cell and chunk counts but different chunk extents"
    throw(ArgumentError(
        "$g was not built against this $role space: $reason. A dependency " *
        "relation is only valid for the spaces it was derived from."))
end

_summarize(v::AbstractVector{Int}) = length(v) <= 6 ? string(v) :
    string(length(v), " chunks from ", first(v), " to ", last(v))
_summarize(v::AbstractUnitRange{Int}) = string(v)
