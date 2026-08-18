# ---------------------------------------------------------------------------
# Stencils on subsets.
#
# `neighbors` / `ring` on a SUBSET of a complete level are the complete level's
# answer clipped to membership — `ring(sub, c, k) == filter(in(sub), ring(complete,
# c, k))`, stated in `src/interface/grid.jl`. Two consequences shape this file:
#
#   * the per-system automata are reached through the complete level grid, so a
#     subset gets the fast path instead of the geometric tree walk; and
#   * distance is the SYSTEM's, so nothing here does a breadth-first walk of its
#     own and a hole in the subset cannot lengthen a path around itself.
#
# The position forms are the same answers read as indices into a data vector,
# and `halo_table` is those for a whole subset at once. On a rooted subtree it
# treats interior and rim cells separately: an interior cell has no neighbour outside the subtree by
# definition, so its row needs no membership test at all.
# ---------------------------------------------------------------------------

# ===========================================================================
# Membership, and the error for a cell that has none
# ===========================================================================

@noinline function _not_a_member(sub, c)
    throw(ArgumentError(
        "$c is not a cell of $sub. On a subset `neighbors` and `ring` are the " *
        "complete level's answer clipped to membership, so a cell the subset " *
        "does not hold has no neighbourhood here — asking the complete level " *
        "would answer about cells this collection does not contain"))
end

@inline function _member_or_throw(sub, c::AbstractCellIndex)
    p = cellposition(sub, c)
    p === nothing && _not_a_member(sub, c)
    return p
end

# ===========================================================================
# The clipped ids
# ===========================================================================

# One shape for both verbs and both containers: run the complete level's answer
# and keep what the subset holds, in the order it came in. That order IS the
# rotational contract, so clipping preserves it for free.
#
# The membership check is the CALLER's, not this function's: an argument is
# evaluated before the call, so checking here would compute a whole `k == 3`
# disc for a cell the subset does not hold and throw afterwards.
# The clip keeps the container it was handed. A clipped ring is a subset of the
# unclipped one, so `k == 1` stays in the `SmallVector` the system answered in
# rather than moving to the heap for the membership test; `k >= 2` has no static
# bound, arrives in a `Vector`, and leaves in one.
function _clip(sub, cells::SmallCollections.SmallVector{N,T}) where {N,T}
    out = SmallCollections.SmallVector{N,T}()
    for nb in cells
        cellposition(sub, nb) === nothing ||
            (out = SmallCollections.push(out, nb))
    end
    return out
end

function _clip(sub, cells)
    out = Vector{eltype(cells)}()
    for nb in cells
        cellposition(sub, nb) === nothing || push!(out, nb)
    end
    return out
end

"""
    neighbors(pg::PartialGrid, c, k = 1; connectivity = Vertex())
    ring(pg::PartialGrid, c, k; connectivity = Vertex())

The complete level's neighbourhood of `c`, clipped to the subset — which is
what [`neighbors`](@ref) means on a subset. The complete grid is
`levelgrid(system(pg), level(pg))`, so this reaches the system's own automaton
and never the geometric tree walk; the clip is one `cellposition` per candidate.

`c` outside `pg` throws an `ArgumentError` rather than answering about cells
`pg` does not hold.
"""
function neighbors(pg::PartialGrid, c::AbstractCellIndex, k::Integer = 1;
        connectivity::Connectivity = Vertex())
    _member_or_throw(pg, c)
    return _clip(pg, neighbors(pg.complete, c, Int(k); connectivity))
end

function ring(pg::PartialGrid, c::AbstractCellIndex, k::Integer;
        connectivity::Connectivity = Vertex())
    _member_or_throw(pg, c)
    return _clip(pg, ring(pg.complete, c, Int(k); connectivity))
end

"""
    neighbors(cv::CellVector, c, k = 1; connectivity = Vertex())
    ring(cv::CellVector, c, k; connectivity = Vertex())

[`neighbors`](@ref) and [`ring`](@ref) on the compressed collection, with the
same clipped-to-membership meaning a [`PartialGrid`](@ref) has. Membership is
the window search — `O(log #windows)` — so the clip costs nothing the vector
was not already able to answer.
"""
function neighbors(cv::CellVector, c::AbstractCellIndex, k::Integer = 1;
        connectivity::Connectivity = Vertex())
    _member_or_throw(cv, c)
    return _clip(cv, neighbors(cv.grid, c, Int(k); connectivity))
end

function ring(cv::CellVector, c::AbstractCellIndex, k::Integer;
        connectivity::Connectivity = Vertex())
    _member_or_throw(cv, c)
    return _clip(cv, ring(cv.grid, c, Int(k); connectivity))
end

# ===========================================================================
# Neighbour counts default to the length of the clipped one-ring. Systems with
# structural degree information may specialize this method.
# ===========================================================================

neighborcount(grid::AbstractGrid, c::AbstractCellIndex;
        connectivity::Connectivity = Vertex()) =
    length(neighbors(grid, c, 1; connectivity))

neighborcount(cv::CellVector, c::AbstractCellIndex;
        connectivity::Connectivity = Vertex()) =
    length(neighbors(cv, c, 1; connectivity))

# ===========================================================================
# The position forms
#
# A bare `Int` is a position, so these are the same two verbs read as indices
# into a data vector laid out against the collection — the id answer mapped
# through `cellposition`, element for element. The ROTATIONAL order is carried
# across: slot `i` of a row names the same direction as slot `i` of the ring,
# which is the whole content of the order contract and the thing an oriented
# stencil reads.
# ===========================================================================

# The generic route, for a grid with no faster membership test than its own
# `cellposition`. A complete level grid lands here too, where the clip is a
# no-op.
function _positions(sub, cells)
    out = Int[]
    for nb in cells
        p = cellposition(sub, nb)
        p === nothing || push!(out, p)
    end
    return out
end

"""
    neighbors(grid::AbstractGrid, p::Int, k = 1; connectivity = Vertex()) -> Vector{Int}
    ring(grid::AbstractGrid, p::Int, k; connectivity = Vertex()) -> Vector{Int}

The neighbourhood of the cell at **position** `p`, as in-set positions in the
same counter-clockwise order the id form answers in. `p` outside
`1:ncells(grid)` is a `BoundsError`, from [`cellindex`](@ref).

See [`neighbors`](@ref) for the contract these keep, and [`halo_table`](@ref)
for the whole grid's worth at once.
"""
neighbors(grid::AbstractGrid, p::Int, k::Integer = 1;
        connectivity::Connectivity = Vertex()) =
    _positions(grid, neighbors(grid, cellindex(grid, p), Int(k); connectivity))

ring(grid::AbstractGrid, p::Int, k::Integer;
        connectivity::Connectivity = Vertex()) =
    _positions(grid, ring(grid, cellindex(grid, p), Int(k); connectivity))

"""
    neighbors(sub::PartialGrid, p::Int, k = 1; connectivity = Vertex()) -> Vector{Int}
    ring(sub::PartialGrid, p::Int, k; connectivity = Vertex()) -> Vector{Int}
    neighbors(cv::CellVector, p::Int, k = 1; connectivity = Vertex()) -> Vector{Int}
    ring(cv::CellVector, p::Int, k; connectivity = Vertex()) -> Vector{Int}

The position forms on the two subset collections: positions in the subset, in
the ring order [`neighbors`](@ref) states.

The complete level's candidates are clipped ONCE here rather than by going
through the id form and testing membership again on the way back — a position
already names a cell of the subset, so there is nothing for the id form's
out-of-set check to decide.
"""
neighbors(pg::PartialGrid, p::Int, k::Integer = 1;
        connectivity::Connectivity = Vertex()) =
    _positions(pg, neighbors(pg.complete, cellindex(pg, p), Int(k); connectivity))

ring(pg::PartialGrid, p::Int, k::Integer;
        connectivity::Connectivity = Vertex()) =
    _positions(pg, ring(pg.complete, cellindex(pg, p), Int(k); connectivity))

neighbors(cv::CellVector, p::Int, k::Integer = 1;
        connectivity::Connectivity = Vertex()) =
    _positions(cv, neighbors(cv.grid, cv[p], Int(k); connectivity))

ring(cv::CellVector, p::Int, k::Integer;
        connectivity::Connectivity = Vertex()) =
    _positions(cv, ring(cv.grid, cv[p], Int(k); connectivity))

# ===========================================================================
# The halo table
# ===========================================================================

halo_table(grid::AbstractGrid, k::Integer = 1;
        connectivity::Connectivity = Vertex()) =
    [neighbors(grid, p, Int(k); connectivity) for p in 1:ncells(grid)]

# Use the cursor sweep for one-ring rows unless the vector retains a rooted
# grid. Other radii use the per-cell implementation.
function halo_table(cv::CellVector, k::Integer = 1;
        connectivity::Connectivity = Vertex(), threaded = true)
    steps = Int(k)
    steps == 1 ||
        return [neighbors(cv, p, steps; connectivity) for p in 1:length(cv)]
    b = cv.backing
    if b isa PartialGrid
        r = _whole_subtree_range(b)
        r === nothing || return _rooted_halo(b, r, connectivity)
    end
    return _swept_rows(cv, connectivity, GOCore.booltype(threaded))
end

# The rooted-subtree fast path. A subtree's INTERIOR is exactly the descendants
# with no neighbour outside it, so an interior cell's whole row is in-set by
# definition and its positions are the complete grid's minus the block offset —
# no membership test, no search. Only the rim, which is O(sqrt(subtree)) of the
# cells, has anything to decide.
#
# Three conditions, all load-bearing:
#
#   * the grid must hold the WHOLE subtree, not a subset of one, or "interior"
#     stops meaning "every neighbour is present";
#   * the iterators must be built with the SAME connectivity the halo is asked
#     for, since `Edge()`'s interior is strictly larger than `Vertex()`'s and
#     using the wrong one would admit a neighbour that is not there; and
#   * `k` must be 1. The interior guarantee applies only to the one-ring:
#     a cell two steps inside the rim can still reach outside at `k == 2`, and
#     there is no "k-interior" iterator to ask.
function halo_table(pg::PartialGrid, k::Integer = 1;
        connectivity::Connectivity = Vertex(), threaded = true)
    steps = Int(k)
    steps == 1 ||
        return [neighbors(pg, p, steps; connectivity) for p in 1:ncells(pg)]
    r = _whole_subtree_range(pg)
    # Unrooted subsets use the position-identical vector sweep.
    r === nothing &&
        return _swept_rows(CellVector(pg), connectivity, GOCore.booltype(threaded))
    return _rooted_halo(pg, r, connectivity)
end

# The subtree's position block, or `nothing` when this grid is not one. The
# count decides it: the ids are validated ascending and inside the root's range,
# so holding as many of them as the range is long means holding all of it.
function _whole_subtree_range(pg::PartialGrid)
    _is_rooted(pg) || return nothing
    has_sorted_subtrees(pg.system) || return nothing
    r = descendant_range(pg.system, pg.root_id, pg.level)
    return length(r) == ncells(pg) ? r : nothing
end

# Separate loops keep each iterator state concrete. Interleaving the two
# iterator types would introduce union-typed dispatch on every step; the
# separate passes instead pay one `cellposition` lookup per cell.
function _rooted_halo(pg::PartialGrid, r::UnitRange{Int}, connectivity::Connectivity)
    sys, root, complete = pg.system, pg.root_id, pg.complete
    l = pg.level
    out = Vector{Vector{Int}}(undef, ncells(pg))
    _inner_rows!(out, complete, InnerCellIterator(sys, root, l; connectivity),
        first(r), last(r), connectivity)
    _rim_rows!(out, complete, EdgeCellIterator(sys, root, l; connectivity),
        first(r), last(r), connectivity)
    return out
end

# `cellposition(complete, ·)` cannot be `nothing` for either the walked cell or
# its neighbours: both are cells OF the complete level, which is the grid being
# asked. The `::Int` says so where a `nothing` would otherwise be silent.
#
# The interior pass takes "every neighbour is inside the block" on the
# iterator's word, and that word is the ONLY thing keeping it from writing a
# wrong row: a position outside the block still lands somewhere in `out` after
# the offset. So it is checked rather than assumed — an interior cell with an
# outside neighbour means the split is broken, and this says so at the cell.
# Both passes check the ROW index the same way, for the same reason: what they
# are told to write is only as trustworthy as the walk that named it.
@noinline _escaped(d, nb) = error(
    "interior cell $d has neighbour $nb outside its own subtree: the " *
    "interior/rim split is wrong, and the halo table would be silently bad")

@noinline _outside(d) = error(
    "walked cell $d is outside the subtree it was walked from: the " *
    "interior/rim split is wrong, and the halo table would be silently bad")

function _inner_rows!(out, complete, inner, lo::Int, hi::Int,
        connectivity::Connectivity)
    for d in inner
        row = Int[]
        for nb in neighbors(complete, d, 1; connectivity)
            q = cellposition(complete, nb)::Int
            @boundscheck (lo <= q <= hi) || _escaped(d, nb)
            push!(row, q - lo + 1)
        end
        p = cellposition(complete, d)::Int
        @boundscheck (lo <= p <= hi) || _outside(d)
        out[p-lo+1] = row
    end
    return out
end

function _rim_rows!(out, complete, rim, lo::Int, hi::Int, connectivity::Connectivity)
    for d in rim
        row = Int[]
        for nb in neighbors(complete, d, 1; connectivity)
            q = cellposition(complete, nb)::Int
            (q < lo || q > hi) && continue
            push!(row, q - lo + 1)
        end
        p = cellposition(complete, d)::Int
        @boundscheck (lo <= p <= hi) || _outside(d)
        out[p-lo+1] = row
    end
    return out
end

# ===========================================================================
# The subset halo
#
# The OUTSIDE face of a subset, next to `halo_table`, which is its inside face
# read as positions. The two are deliberately adjacent and deliberately not the
# same verb: `halo_table` says which of a subset's OWN cells each of its cells
# touches, and `halo` says which cells the subset does not hold touch it. A
# stencil needs both — the table to gather from local storage, the halo to name
# the extra fetch list — and neither can be derived from the other.
# ===========================================================================

"""
    halo(pg::PartialGrid; connectivity = Vertex())

Return cells immediately outside a subset grid, lazily. A rooted grid containing
a complete subtree delegates to [`SubtreeHaloIterator`](@ref); other subsets use
the outside-first subset walk. `_whole_subtree_range` determines the path.

The subset walk prunes by the subset's own [`subset_span`](@ref)s — the root, if
there is one, is not consulted at all. Construction is `O(1)`: nothing about the
answer, and nothing about the input either, is computed before the first
`iterate`.
"""
function halo(pg::PartialGrid; connectivity::Connectivity = Vertex())
    _whole_subtree_range(pg) === nothing ||
        return SubtreeHaloIterator(pg.system, pg.root_id, pg.level; connectivity)
    return SubsetHaloIterator(pg, connectivity, subset_halo_engine(pg.system, pg,
        pg.complete, pg.level, connectivity))
end

"""
    halo(cv::CellVector; connectivity = Vertex())

The cells immediately outside the compressed collection, lazily. See
[`halo`](@ref).

A [`CellVector`](@ref) always uses the subset walk because it stores windows but
no root ancestor. Build a grid from the root cell to enable the subtree engine.

Membership is the window search, `O(log #windows)`, so the walk's per-candidate
cost is lower here than on a `PartialGrid` over a bare id vector — and so is its
per-NODE cost, since [`subset_span`](@ref) reads a stored run as a block the
vector holds entire.
"""
halo(cv::CellVector; connectivity::Connectivity = Vertex()) =
    SubsetHaloIterator(cv, connectivity, subset_halo_engine(system(cv), cv,
        cv.grid, level(cv), connectivity))

# ===========================================================================
# The chunk-plus-halo stencil: the two faces above, addressed together
#
# A chunked stencil addresses subset and halo cells in one `[chunk; halo]`
# buffer. `halo_table` clips rows to the subset, while `halo` names exterior
# cells without mapping them back to rows.
#
# `stencil_table` therefore has a separate contract:
#
#   * rows are COMPLETE, never clipped, and that is checked rather than hoped;
#   * entries index the CONCATENATED buffer, not the subset;
#   * the result is CSR, not a vector of row vectors.
#
# Rows carry the same counter-clockwise order every other neighbour verb here
# does; what is special is that they are never short, so the cycle is not just
# ordered but closed.
# ===========================================================================

"""
    StencilTable <: AbstractVector

The stencil of a chunk read together with its halo: `t[i]` is the FULL
neighbourhood of the chunk's `i`-th cell, as indices into the concatenated
`[chunk; halo]` buffer — `1:nchunk` for a cell of the chunk, `nchunk+1 :
nchunk+nhalo` for one of the halo. Built by [`stencil_table`](@ref), which is
where the contract is stated.

CSR, and the two arrays are public: `t.offsets` is `nchunk + 1` long with row `i`
occupying `t.offsets[i] : t.offsets[i+1] - 1`, and `t.indices` is the flat
neighbour array those slices cut. `t.nchunk` and `t.nhalo` are the two halves'
lengths, so `t.nchunk + t.nhalo` is the buffer length a row can address and the
one to check a fetch against.

Indexing gives the row as a `view`, so `t[i]` allocates nothing and a stencil
pass can be written either way — over rows, or over `t.indices` directly with
`t.offsets` as the row bounds, which is the loop a kernel wants.
"""
struct StencilTable <: AbstractVector{SubArray{Int,1,Vector{Int},Tuple{UnitRange{Int}},true}}
    offsets::Vector{Int}
    indices::Vector{Int}
    nchunk::Int
    nhalo::Int
end

Base.size(t::StencilTable) = (t.nchunk,)
Base.IndexStyle(::Type{StencilTable}) = Base.IndexLinear()

# `offsets` is `nchunk + 1` long by construction, so `i + 1` is in range as soon
# as `i` is, and reading the two bounds unchecked is safe. The VIEW is left
# CHECKED, and deliberately: both arrays are public fields, so a hand-assembled
# table would otherwise be able to read past the flat array's end, and one range
# comparison per row is nothing beside reading the row.
Base.@propagate_inbounds function Base.getindex(t::StencilTable, i::Int)
    @boundscheck checkbounds(t, i)
    lo = @inbounds t.offsets[i]
    hi = @inbounds t.offsets[i+1] - 1
    return view(t.indices, lo:hi)
end

Base.show(io::IO, t::StencilTable) = print(io, "StencilTable(nchunk=", t.nchunk,
    ", nhalo=", t.nhalo, ", entries=", length(t.indices), ")")
Base.show(io::IO, ::MIME"text/plain", t::StencilTable) = show(io, t)

# Which slot of the buffer's FIRST half a neighbour occupies, or `0` — free as a
# sentinel, since a position is one-based — when the chunk does not hold it.
#
# Two shapes, chosen once in `stencil_table` and passed by value, so the inner
# loop is monomorphic on whichever it got. This is the same split
# `_whole_subtree_range` makes for `halo_table` and `halo`, deliberately read
# from the same function: the halo a caller was handed came out of `halo(pg)`,
# which branched on that predicate, so a second notion of "is a subtree" here
# could hand back a table addressed against a different halo than the one that
# was fetched.
struct BlockChunk
    lo::Int
    hi::Int
end

# The unnamed arguments are the ones the shape does not read, as elsewhere in
# this package: the block decides by POSITION alone, the membership form by the
# CELL alone, and naming what each ignores would suggest it consults it.
@inline _chunk_slot(b::BlockChunk, ::AbstractCellIndex, q::Int) =
    (b.lo <= q <= b.hi) ? q - b.lo + 1 : 0

# The general subset: a hole, a scattered id set, or a system with no
# `descendant_range` to give the block. `O(log nchunk)` per neighbour against
# the block form's `O(1)`, which is the whole reason the block form exists.
struct MemberChunk{S}
    subset::S
end

@inline function _chunk_slot(m::MemberChunk, nb::AbstractCellIndex, ::Int)
    p = cellposition(m.subset, nb)
    return p === nothing ? 0 : p
end

@noinline _stencil_k(k) = throw(ArgumentError(
    "stencil_table: k must be 1, got $k. This is the ONE-RING table: `halo` " *
    "and `subtree_halo` produce a width-1 margin, so a one-ring is what a " *
    "fetched halo can complete, and answering k = 2 from it would give SHORT " *
    "rows — the defect this verb exists to remove. `halo_table` is the in-set " *
    "table at any k"))

@noinline _stencil_unsorted() = throw(ArgumentError(
    "stencil_table: the halo positions must be strictly ascending. That is " *
    "the order `halo` and `subtree_halo` emit in, and it is what the binary " *
    "search that addresses the halo half of the buffer rests on — an " *
    "unsorted list would silently misaddress rows rather than fail"))

@noinline _stencil_incomplete(d, nb, q) = throw(ArgumentError(
    "stencil_table: cell $d has neighbour $nb at level position $q, which is " *
    "in neither the chunk nor the halo, so its row cannot be completed. The " *
    "halo passed is not this chunk's one-ring halo under this connectivity — " *
    "the usual causes are a halo collected under Vertex() and a table asked " *
    "for under Edge() or the reverse, a halo of a different chunk, or a " *
    "truncated fetch list"))

"""
    stencil_table(pg::PartialGrid, halo_positions::AbstractVector{<:Integer}, k = 1;
                  connectivity = Vertex()) -> StencilTable

The positional stencil of a chunk read together with its halo. Row `i` is the
**complete** one-ring of `cellindex(pg, i)`, in the system's rotational order,
as indices into the concatenated `[chunk; halo]` buffer: `1:ncells(pg)` names a
cell of the chunk and `ncells(pg)+j` names `halo_positions[j]`.

This is the addressing [`halo_table`](@ref) does not give. That verb is IN-SET,
so at the rim its rows are short; here the missing neighbours are exactly the
halo, and pointing at them is the whole content of this function.

`halo_positions` is the halo's positions on `levelgrid(system(pg), level(pg))`,
strictly ascending — the fetch list itself:

```julia
grid  = levelgrid(sys, l)
pg    = PartialGrid(sys, root, l)                       # the chunk
hpos  = [cellposition(grid, x) for x in halo(pg)]       # the margin to fetch
buf   = vcat(read(store, descendant_range(sys, root, l)), read(store, hpos))
table = stencil_table(pg, hpos)
mean_of_neighbours = [sum(@view buf[row]) / length(row) for row in table]
```

# The three decisions, and why

**POSITIONS, MATERIALISED, NOT THE ITERATOR.** [`halo`](@ref) is lazy by design
and this function refuses it (loudly — there is a method whose only job is to
say so). A table of buffer indices presupposes a buffer, and there is no buffer
until the halo has been walked and read; a caller who has one already owns this
array, so taking the iterator would mean either walking the halo twice or
hiding an `nhalo`-sized `collect` inside a function whose caller already paid
for it. Positions rather than ids for the same reason: the buffer is indexed by
integer, `cellposition` is what turns the walk into a fetch list, and searching
positions asks nothing of the id order beyond what the fetch already asked.

Ascending is not sorted for and not assumed: it is CHECKED, in one allocation-free
pass, because the binary search that addresses the halo half rests on it and a
violation would misaddress rows rather than fail.

The list need not be MINIMAL, only ascending and covering. A caller who fetched
a two-ring margin gets a one-ring table addressed into it, with the outer ring's
slots simply unreferenced; what is refused is a margin that is too *small*.

**CSR, NOT A VECTOR OF ROWS.** `nchunk` separate row vectors is `nchunk`
allocations to build a table whose rows are read once each. The flat form is one
`offsets` array and one `indices` array, both output-sized, and
[`StencilTable`](@ref) exposes them; `t[i]` is a `view`, so reading by row costs
nothing either.

**`k == 1` ONLY.** [`halo`](@ref) and [`subtree_halo`](@ref) produce a width-1
margin — nothing in this package produces a wider one — so a one-ring is what a
fetched halo can complete, and anything else is an `ArgumentError` rather than a
short row. The refusal is deliberately in front of the completeness check rather
than behind it: the check would catch a `k == 2` request against a one-ring
margin, but only after the caller had built a workflow around a shape no verb
here supplies. Accepting a wider margin (above) and answering a wider table are
different things, and only the first is offered.

# Completeness, which is checked

Every neighbour of every chunk cell is looked up, and one that is in neither
half is an `ArgumentError` naming the cell, the neighbour and its position. So a
`Vertex()` halo handed to an `Edge()` call, a halo of the wrong chunk, or a
truncated fetch list fails at the first row that would have been short. Nothing
here returns a partial row.

# The two paths

A [`PartialGrid`](@ref) holding a whole rooted subtree has its chunk positions
as one contiguous block, and membership is then two integer comparisons and a
subtraction. Everything else — a hole, a forgotten root, an arbitrary id list —
resolves membership through `cellposition`, `O(log nchunk)`. The split is
`_whole_subtree_range`, the same predicate [`halo`](@ref) and
[`halo_table`](@ref) branch on, read from the same function so that "is a
subtree" cannot come to mean two things: the halo a caller passes here came out
of `halo(pg)`, which asked that question first.

**A HOLE IS ADDRESSED, NOT AN ERROR.** `halo(pg)` counts a punched-out interior
cell as a halo cell, so a holed chunk's rows complete through the halo half like
any other. That is the hole law of [`halo`](@ref) paying for itself.

**A5 IS SUPPORTED**, on the membership path, at every root and depth — it is the
one system with no [`descendant_range`](@ref), so it never takes the block path.
Its subtrees do in fact occupy contiguous position blocks, but no primitive in
this package promises that, and nothing here depends on it.

ONE CONTAINER, DELIBERATELY. [`halo`](@ref) answers on three collections; this
verb takes only a [`PartialGrid`](@ref), because a chunk is a thing you READ —
one block of storage plus one margin — and that is the shape a `PartialGrid`
carries. A [`CellVector`](@ref) reaches it as `PartialGrid(system(cv), level(cv),
cv)`, which copies nothing, though membership then goes through the vector's
`getindex` rather than its own `O(log #windows)` window search; a
[`CellLookup`](@ref) unwraps to its `CellVector` first.

`O(nchunk · degree · log nhalo)` time. The only allocation of this function's own
is the two output arrays: no intermediate is sized by the halo, there is no hash
map anywhere, and the halo is neither copied nor sorted — which is what the
ascending requirement buys. Whatever the system's own `neighbors` allocates per
call it allocates here too, unchanged.

See [`halo_table`](@ref) for the in-set table, [`halo`](@ref) for the fetch list
this addresses, and [`StencilTable`](@ref) for the layout.
"""
function stencil_table(pg::PartialGrid, halo_positions::AbstractVector{<:Integer},
        k::Integer = 1; connectivity::Connectivity = Vertex())
    Int(k) == 1 || _stencil_k(k)
    Helpers.strictly_increasing(halo_positions) || _stencil_unsorted()
    r = _whole_subtree_range(pg)
    chunk = r === nothing ? MemberChunk(pg) : BlockChunk(first(r), last(r))
    return _stencil_rows(pg, chunk, halo_positions, connectivity)
end

# The refusal `halo` invites. Written as a method rather than left to
# `MethodError` because passing the iterator is the one mistake the lazy design
# makes natural, and the answer to it is a sentence, not a signature.
@noinline stencil_table(::PartialGrid,
        ::Union{SubtreeHaloIterator,SubsetHaloIterator}, ::Vararg{Any}; kw...) =
    throw(ArgumentError(
        "stencil_table needs the halo MATERIALISED as ascending positions, not " *
        "as the iterator: its rows are indices into a `[chunk; halo]` buffer, " *
        "and that buffer does not exist until the halo has been walked and " *
        "read. Pass `[cellposition(levelgrid(sys, l), x) for x in halo(pg)]` — " *
        "the same fetch list the read used"))

# One pass, in POSITION order, which is what lets the offsets be filled as the
# rows are produced. `halo_table`'s rooted path splits the walk into
# `InnerCellIterator` and `EdgeCellIterator` instead; that split is not reused
# here, for two reasons. It would cost the single pass — the two walks emit
# interleaved, and CSR rows have to be appended in row order — and it would buy
# nothing, because `_chunk_slot` already retires every interior cell's
# neighbours before the halo search is reached, so an interior row does no work
# a "this cell is interior" flag could save. The split that IS reused is
# `_whole_subtree_range`, one function above, which is the one that has to agree
# with `halo`.
function _stencil_rows(pg::PartialGrid, chunk, halo_positions::AbstractVector,
        connectivity::Connectivity)
    complete = pg.complete
    n = ncells(pg)
    m = length(halo_positions)
    offsets = Vector{Int}(undef, n + 1)
    @inbounds offsets[1] = 1
    indices = Int[]
    # The degree ceiling, which is the row length almost everywhere: this is the
    # output's own size to within the pentagons and seam cells, not a guess, and
    # asking for it once is what keeps the append from reallocating on the way
    # up. Deliberately not a second pass to count exactly — that pass would have
    # to compute every neighbourhood twice.
    sizehint!(indices, n * _hint_degree(_capacity(pg.system, connectivity)))
    for i in 1:n
        d = cellindex(pg, i)
        for nb in neighbors(complete, d, 1; connectivity)
            # Cannot be `nothing`: `nb` is a cell OF the complete level, which is
            # the grid it came from. The `::Int` says so where a `nothing` would
            # otherwise flow silently into the search below.
            q = cellposition(complete, nb)::Int
            slot = _chunk_slot(chunk, nb, q)
            if slot == 0
                j = Helpers.sorted_index(halo_positions, q)
                j == 0 && _stencil_incomplete(d, nb, q)
                slot = n + j
            end
            push!(indices, slot)
        end
        @inbounds offsets[i+1] = length(indices) + 1
    end
    return StencilTable(offsets, indices, n, m)
end

# ===========================================================================
# Cross-level adjacency on a multi-order set
# ===========================================================================

"""
    member_neighbors(set::MultiOrderCellSet, c; connectivity = Vertex()) -> Vector

The members of `set` that share a boundary with the member `c` — its neighbours
in a set whose cells sit at different levels, so a neighbour may be COARSER than
`c` (an ancestor of one of its system-neighbours) or FINER (a cell deep inside
one), and this is the one call that finds both.

`c` must be a member; anything else is an `ArgumentError`. `c` itself is never
in the result. The order is the package's one order — counter-clockwise about
`cell_centroid(system(set), c)` seen from outside, as [`neighbors`](@ref)
states — measured on each neighbour's own centroid, whatever level it sits at.
The start is the smallest-`(level, position)` neighbour, and exact azimuth ties
break the same way, so the answer is deterministic and independent of how the
walk found them.

# The algorithm

A member's subtree is a region, and two regions of the hierarchy touch exactly
when two of their cells touch at a common depth. The set already names one:
its REFERENCE LEVEL `L`, the depth its covering guarantee is stated at and no
shallower than any member. So the question becomes a question about level `L`,
and it is answered without expanding anything:

 1. walk `c`'s subtree RIM at `L` — [`EdgeCellIterator`](@ref), `O(rim)` time
    and `O(depth)` memory, and the rim is where every contact with the outside
    lives by its own definition;
 2. take each rim cell's one-ring at `L`, under the requested connectivity;
 3. map each of those cells to the member that contains it, if any. Members are
    disjoint subtrees, so there is at most one, and on a sorted-subtree system
    it is found by binary search over the set's own curve keys: the keys are the
    members' reference-level range starts, ascending and disjoint, so
    `searchsortedlast` plus one range test decides it. A neighbour inside `c`'s
    own subtree maps back to `c` and is dropped.

`O(|rim(c, L)| · degree · log|set|)` time and `O(|answer|)` memory beyond the
rim walk's own `O(depth)`. The rim is the square root of the subtree, so a
member many levels above `L` costs the perimeter of its block and never its
area.

# Exactness, and where it is the hierarchy's answer rather than geometry's

Two members are neighbours here when their level-`L` cells are, which is
*geometric* boundary sharing exactly where the system's refinement is congruent
— HEALPix, S2 and ISEA4R, whose four children tile their parent, so a member's
footprint is the union of its level-`L` descendants'. Under `Edge()` the same
equivalence holds for shared edges, so vertex-only contact is excluded rather
than approximated away.

Where children do not tile their parent — IGEO7 and H3 (aperture 7) and A5 —
a member's footprint is NOT its descendants' union, and no level-`L` statement
can be a statement about the drawn polygons; see [`MultiOrderCoverage`](@ref)
for the size of that gap. The answer there is the hierarchy's, which is the
same relation [`subtree_border`](@ref) is defined by, and it is consistent with
every other subtree verb in this package. That is a carve-out about the
SYSTEMS, not about this walk.

!!! note "A5 pays for its missing primitives here too"
    Without [`has_sorted_subtrees`](@ref) there are no curve keys to binary
    search, so the member lookup is a set built per call, `O(|set|)`; and the
    rim iterator materialises the subtree rather than walking it. The answer
    is unchanged.
"""
function member_neighbors(set::MultiOrderCellSet, c::AbstractCellIndex;
        connectivity::Connectivity = Vertex())
    sys = set.system
    L = set.reference_level
    grid = levelgrid(sys, L)
    lookup = _member_lookup(set)
    home = _home_index(lookup, set, L, c)
    home === nothing && throw(ArgumentError(
        "$c is not a member of this MultiOrderCellSet, so it has no neighbours " *
        "in it; the verb is about members, not about cells of a level"))
    hits = Int[]
    seen = Set{Int}((home,))
    for d in EdgeCellIterator(sys, c, L; connectivity)
        for nb in neighbors(grid, d, 1; connectivity)
            i = _member_of(lookup, set, grid, nb)
            (i === nothing || i in seen) && continue
            push!(seen, i)
            push!(hits, i)
        end
    end
    length(hits) <= 1 && return [set.cells[i] for i in hits]
    # `(level, key)` IS ascending `(level, position)`: within one level the keys
    # are that level's own order, whichever branch built them. It is the tie
    # break and the start anchor, not the order — the order is the ring's.
    rank(i) = (level(set.cells[i]), set.keys[i])
    sort!(hits; by = rank)
    # Each member is read on its own level's grid, because a coarse neighbour
    # has no centroid at `L`.
    centroid(x) = cell_centroid(levelgrid(sys, level(x)), x)
    centre = centroid(c)
    anchor = centroid(set.cells[first(hits)])
    e1, e2 = _tangent_basis(centre, anchor)
    spoke = _azimuth(centre, e1, e2, anchor)
    sort!(hits; by = i -> (_phase(_azimuth(centre, e1, e2,
            centroid(set.cells[i])) - spoke), rank(i)))
    return [set.cells[i] for i in hits]
end

# The curve keys are the members' reference-level range starts — ascending and
# disjoint — so no index has to be built at all.
struct CurveKeyLookup end

# A5 has no ranges to search, so the members go into a dictionary and the walk
# up from a level-`L` cell asks it once per generation.
struct AncestorLookup{ID}
    index::Dict{ID,Int}
end

_member_lookup(set::MultiOrderCellSet) =
    has_sorted_subtrees(set.system) ? CurveKeyLookup() :
    AncestorLookup(Dict(c => i for (i, c) in enumerate(set.cells)))

# Is `c` a member, and which one — the same two structures answering about a
# cell at its OWN level rather than at the reference level. Binary search again
# rather than a scan, so the membership check is not the expensive step.
function _home_index(::CurveKeyLookup, set::MultiOrderCellSet, L::Int,
        c::AbstractCellIndex)
    level(c) <= L || return nothing
    i = searchsortedfirst(set.keys, first(descendant_range(set.system, c, L)))
    (i <= length(set.cells) && set.cells[i] == c) || return nothing
    return i
end

_home_index(lk::AncestorLookup, ::MultiOrderCellSet, ::Int, c::AbstractCellIndex) =
    get(lk.index, c, nothing)

function _member_of(::CurveKeyLookup, set::MultiOrderCellSet, grid::AbstractGrid,
        nb::AbstractCellIndex)
    p = cellposition(grid, nb)
    p === nothing && return nothing
    i = searchsortedlast(set.keys, p)
    i >= 1 || return nothing
    return p <= last(descendant_range(set.system, set.cells[i], level(grid))) ?
           i : nothing
end

function _member_of(lk::AncestorLookup, set::MultiOrderCellSet, grid::AbstractGrid,
        nb::AbstractCellIndex)
    top = first(levels(set.system))
    a = nb
    while true
        i = get(lk.index, a, nothing)
        i === nothing || return i
        level(a) <= top && return nothing
        a = Base.parent(set.system, a)
    end
end
