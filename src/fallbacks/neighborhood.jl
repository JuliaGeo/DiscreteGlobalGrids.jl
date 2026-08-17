# ---------------------------------------------------------------------------
# The neighborhood iterator: every cell's clipped one-ring, in storage order,
# with positions resolved as a side effect of the walk rather than by one
# search per neighbour.
#
# Stencil kernels over a subset spend their time translating ids back to
# positions — ~20 `cellposition` calls per cell for a one-ring metric, each an
# `O(log #windows)` search on a `CellVector`. Walking the subset in storage
# order makes almost all of those searches redundant: the walked cell's window
# is known, and ring members overwhelmingly land in that window or the one the
# previous member landed in. The iterator keeps that state — a window cursor —
# and yields cells already paired with their positions, as handles.
# ---------------------------------------------------------------------------

"""
    SubsetPositionedCell{C}

A cell paired with its position in the subset collection that minted it — the
element type of the one-arg [`neighbors`](@ref) iterator, never
user-constructed. The position is relative to the minting collection's axis
and is **trusted unconditionally** wherever a handle is accepted: indexing a
`Cells`-lookup raster with one reads storage at that position directly, with
no membership check. Against any other axis the position is meaningless; the
bare cell is the portable identity — [`cellid`](@ref) unwraps it, and
`cellposition(other, cellid(h))` re-resolves it, so correct cross-axis use is
one accessor away.

`==`, `hash`, `show` and `convert(::Type{C})` delegate to the cell, so a
handle compares, hashes and prints as the cell it names.
[`cellposition`](@ref)`(h)` reads the minted position back.
"""
struct SubsetPositionedCell{C<:AbstractCellIndex}
    cell::C
    position::Int
end

"""
    cellid(h::SubsetPositionedCell)
    cellid(c::AbstractCellIndex)

The bare cell: a positioned handle unwrapped, an already-bare cell unchanged.
The escape from a handle's single-axis contract — `cellposition(cv,
cellid(h))` resolves against any axis, where the handle's own position is
valid only against the one that minted it.
"""
cellid(h::SubsetPositionedCell) = h.cell
cellid(c::AbstractCellIndex) = c

"""
    cellposition(h::SubsetPositionedCell) -> Int

The position the handle was minted with — valid only against the collection
that minted it, and trusted there without a check.
"""
cellposition(h::SubsetPositionedCell) = h.position

level(h::SubsetPositionedCell) = level(h.cell)

Base.:(==)(a::SubsetPositionedCell, b::SubsetPositionedCell) = a.cell == b.cell
Base.:(==)(a::SubsetPositionedCell, b::AbstractCellIndex) = a.cell == b
Base.:(==)(a::AbstractCellIndex, b::SubsetPositionedCell) = a == b.cell
Base.hash(h::SubsetPositionedCell, u::UInt) = hash(h.cell, u)
Base.show(io::IO, h::SubsetPositionedCell) = show(io, h.cell)
Base.convert(::Type{C}, h::SubsetPositionedCell{C}) where {C<:AbstractCellIndex} =
    h.cell

# ===========================================================================
# The window cursor
#
# Three integers of iteration state: the storage position `k`, the index `wj`
# of the window holding cell `k`, and the index `hj` of the last window a ring
# member resolved into. `k` ascends one at a time, so `wj` advances by
# comparison, never by search; a ring member checks `wj`, then `hj`, and only
# then pays the binary search — updating `hj` so the next member of a
# cross-window run gets the hit for free. `0` is the "not held" sentinel, as
# in `_chunk_slot`: a position is one-based, so it is free.
# ===========================================================================

@inline _wbase(w::RangeWindows, j::Int) = j == 1 ? 0 : @inbounds w.offsets[j-1]

@inline function _advance(w::RangeWindows, wj::Int, k::Int)
    @inbounds while k > w.offsets[wj]
        wj += 1
    end
    return wj
end
@inline _advance(::PositionWindows, wj::Int, k::Int) = k

@inline _leaf_at(w::RangeWindows, wj::Int, k::Int) =
    @inbounds w.starts[wj] + (k - _wbase(w, wj)) - 1
@inline _leaf_at(w::PositionWindows, ::Int, k::Int) = @inbounds w.positions[k]

# Leaf position -> (concatenation position or 0, updated hint).
@inline function _cursor_find(w::RangeWindows, wj::Int, hj::Int, p::Int)
    @inbounds if w.starts[wj] <= p <= w.stops[wj]
        return _wbase(w, wj) + (p - w.starts[wj]) + 1, hj
    end
    @inbounds if hj != wj && w.starts[hj] <= p <= w.stops[hj]
        return _wbase(w, hj) + (p - w.starts[hj]) + 1, hj
    end
    j = searchsortedfirst(w.stops, p)
    (j <= length(w.stops) && p >= @inbounds w.starts[j]) || return 0, hj
    return _wbase(w, j) + (p - @inbounds w.starts[j]) + 1, j
end

# A position list has no runs to be inside of; the hint is the one entry the
# previous member resolved to, worth one comparison before the search.
@inline function _cursor_find(w::PositionWindows, ::Int, hj::Int, p::Int)
    ps = w.positions
    @inbounds if hj <= length(ps) && ps[hj] == p
        return hj, hj
    end
    j = searchsortedfirst(ps, p)
    (j <= length(ps) && @inbounds(ps[j]) == p) || return 0, hj
    return j, j
end

# The clip with positions: each ring member that survives membership, paired
# with the position the cursor resolved. The output capacity is the SYSTEM's
# `max_neighbors` bound rather than the incoming container's parameter,
# because the incoming type is a union over the ring verb's `k` — the systems
# answer `k >= 2` in a `Vector` — and a concrete output is what lets the
# yielded tuple stay off the heap. The bound is the same static capacity the
# systems' own one-ring containers carry. The complete level holds every
# same-level cell, so its `cellposition` cannot be `nothing`.
@inline function _positioned(cv::CellVector, wj::Int, hj::Int, cells,
        ::Val{M}) where {M}
    w = cv.windows
    out = SmallCollections.SmallVector{M,SubsetPositionedCell{eltype(cells)}}()
    for nb in cells
        q, hj = _cursor_find(w, wj, hj, cellposition(cv.grid, nb)::Int)
        q == 0 || (out = SmallCollections.push(out, SubsetPositionedCell(nb, q)))
    end
    return out, hj
end

# ===========================================================================
# The iterator
# ===========================================================================

"""
    neighbors(cv::CellVector; connectivity = Vertex())
    neighbors(pg::PartialGrid; connectivity = Vertex())

Iterate the whole subset in storage order, yielding `(cell, nbrs)` with both
sides as [`SubsetPositionedCell`](@ref) handles: `cell` is the `k`-th cell
paired with `k`, and `nbrs` is its one-ring clipped to membership — the same
cells the two-arg `neighbors(cv, cell)` answers, in the same order and
container discipline — each paired with its own position.

The point of the form is amortized position resolution. Iteration state is a
window cursor: cells are walked without any search, and each ring member is
tried against the walked cell's window and the last window that answered
before paying the `O(log #windows)` search the two-arg form pays every time.
The whole sweep allocates nothing beyond what the system's own ring does.

A `PartialGrid` iterates as `CellVector(pg)` — `O(1)` for a rooted subtree or
a vector-backed grid, one position-resolution pass for a bare id vector — and
the minted positions are the grid's own, since the vector preserves them.
"""
neighbors(cv::CellVector; connectivity::Connectivity = Vertex()) =
    NeighborhoodIterator(cv, connectivity)

neighbors(pg::PartialGrid; connectivity::Connectivity = Vertex()) =
    NeighborhoodIterator(CellVector(pg), connectivity)

# `M` is `max_neighbors(system(cv), connectivity)` — a constant of the two
# singleton types, so construction infers — and is the static capacity of
# every yielded ring.
struct NeighborhoodIterator{M,CV<:CellVector,CN<:Connectivity}
    cv::CV
    connectivity::CN
end

NeighborhoodIterator(cv::CellVector, conn::Connectivity) =
    NeighborhoodIterator{max_neighbors(system(cv), conn),typeof(cv),typeof(conn)}(
        cv, conn)

Base.length(it::NeighborhoodIterator) = length(it.cv)
Base.IteratorSize(::Type{<:NeighborhoodIterator}) = Base.HasLength()
Base.IteratorEltype(::Type{<:NeighborhoodIterator}) = Base.EltypeUnknown()

Base.show(io::IO, it::NeighborhoodIterator) =
    print(io, "neighbors(", it.cv, "; connectivity=", it.connectivity, ")")
Base.show(io::IO, ::MIME"text/plain", it::NeighborhoodIterator) = show(io, it)

Base.iterate(it::NeighborhoodIterator) = _neighborhood_next(it, 1, 1, 1)
Base.iterate(it::NeighborhoodIterator, st::NTuple{3,Int}) =
    _neighborhood_next(it, st[1], st[2], st[3])

function _neighborhood_next(it::NeighborhoodIterator{M}, k::Int, wj::Int,
        hj::Int) where {M}
    cv = it.cv
    k > length(cv) && return nothing
    w = cv.windows
    wj = _advance(w, wj, k)
    c = cellindex(cv.grid, _leaf_at(w, wj, k))
    nbrs, hj = _positioned(cv, wj, hj,
        neighbors(cv.grid, c, 1; connectivity = it.connectivity), Val(M))
    return (SubsetPositionedCell(c, k), nbrs), (k + 1, wj, hj)
end

# ===========================================================================
# The closure form: the grid owns the loop
#
# `mapneighbors` and `foreachneighbors` run the same cursor sweep the
# iterator runs, but the FRAME owns traversal, chunking and the writes: `f`
# returns a cell's value(s) and the frame stores them at the position it
# already knows. Threaded chunks are disjoint position ranges, so a pure `f`
# is race-free by construction and never touches an output index.
# ===========================================================================

# The window holding storage position `k`, by search — the once-per-range (or
# once-per-visit, under a permutation) initialization the cursor's
# comparison-only `_advance` starts from.
@inline _window_at(w::RangeWindows, k::Int) = searchsortedfirst(w.offsets, k)
@inline _window_at(::PositionWindows, k::Int) = k

"""
    StorageOrder()

The default traversal of [`mapneighbors`](@ref) and
[`foreachneighbors`](@ref): cells in storage order, the walk the window
cursor is built for. The alternative is an explicit permutation of
`1:length(cv)` — see `order` in [`mapneighbors`](@ref).
"""
struct StorageOrder end

# The engine: call `g(k, cell, nbrs)` for every storage position in `r`, in
# order, with one cursor for the whole range. `g` is the frame's store — the
# user's `f` never sees a position it could misuse.
function _sweep!(g::G, cv::CellVector, conn::Connectivity, r::UnitRange{Int},
        ::Val{M}) where {G,M}
    isempty(r) && return nothing
    w = cv.windows
    wj = _window_at(w, first(r))
    hj = wj
    for k in r
        wj = _advance(w, wj, k)
        c = cellindex(cv.grid, _leaf_at(w, wj, k))
        nbrs, hj = _positioned(cv, wj, hj,
            neighbors(cv.grid, c, 1; connectivity = conn), Val(M))
        g(k, SubsetPositionedCell(c, k), nbrs)
    end
    return nothing
end

# The permutation engine: visit `perm[j]` for `j` in `r`. Successive visits
# share no locality, so no hint is carried between cells — each visit pays
# one window search for its OWN window, which its ring members still resolve
# against before falling back to the full search.
function _sweep_perm!(g::G, cv::CellVector, conn::Connectivity,
        perm::AbstractVector{<:Integer}, r::UnitRange{Int}, ::Val{M}) where {G,M}
    w = cv.windows
    for j in r
        k = Int(@inbounds perm[j])
        wj = _window_at(w, k)
        c = cellindex(cv.grid, _leaf_at(w, wj, k))
        nbrs, _ = _positioned(cv, wj, wj,
            neighbors(cv.grid, c, 1; connectivity = conn), Val(M))
        g(k, SubsetPositionedCell(c, k), nbrs)
    end
    return nothing
end

# Length + in-range + no duplicate IS a permutation. Checked up front rather
# than trusted: the frame stores at `perm[j]`, so a duplicate would silently
# overwrite one cell's output and leave another's slot undefined.
function _check_permutation(perm::AbstractVector{<:Integer}, n::Int)
    length(perm) == n || throw(ArgumentError(
        "order must be a permutation of 1:$n, got length $(length(perm))"))
    seen = falses(n)
    for x in perm
        k = Int(x)
        1 <= k <= n || throw(ArgumentError(
            "order must be a permutation of 1:$n, got entry $k"))
        @inbounds seen[k] && throw(ArgumentError(
            "order must be a permutation of 1:$n, but it visits $k twice"))
        @inbounds seen[k] = true
    end
    return nothing
end

# Contiguous ranges, one task and one cursor each. The split point does not
# touch the answer — every chunk writes only its own positions — so no
# window alignment is needed; a chunk starting mid-window pays one
# `_window_at` search and walks on. The ranges partition `1:n` in order, so
# anything assembled chunk-by-chunk in range order is the sequential result.
function _chunk_ranges(n::Int)
    step = cld(n, min(n, 2 * Threads.nthreads()))
    return [lo:min(lo + step - 1, n) for lo in 1:step:n]
end

function _foreach_chunk(body!::F, n::Int, ::GOCore.True) where {F}
    n == 0 && return nothing
    @sync for r in _chunk_ranges(n)
        Threads.@spawn body!(r)
    end
    return nothing
end

_foreach_chunk(body!::F, n::Int, ::GOCore.False) where {F} = body!(1:n)

_run!(g::G, cv::CellVector, conn::Connectivity, ::StorageOrder, thr,
    ::Val{M}) where {G,M} =
    _foreach_chunk(r -> _sweep!(g, cv, conn, r, Val(M)), length(cv), thr)

function _run!(g::G, cv::CellVector, conn::Connectivity,
        perm::AbstractVector{<:Integer}, thr, ::Val{M}) where {G,M}
    _check_permutation(perm, length(cv))
    return _foreach_chunk(r -> _sweep_perm!(g, cv, conn, perm, r, Val(M)),
        length(perm), thr)
end

@noinline _run!(g, cv::CellVector, conn, order, thr, v) = throw(ArgumentError(
    "order must be StorageOrder() or a permutation of 1:length(cv), " *
    "got $(typeof(order))"))

# --- the output frame -------------------------------------------------------

# One vector for a plain return, one vector PER COMPONENT for a concrete
# tuple return. The split is decided by the inferred return type, once: a
# kernel that sometimes returns a tuple and sometimes not infers abstract and
# lands in the single-vector shape, which can hold anything.
_outputs(::Type{T}, n::Int) where {T} =
    T <: Tuple && isconcretetype(T) ?
    ntuple(j -> Vector{fieldtype(T, j)}(undef, n), fieldcount(T)) :
    Vector{T}(undef, n)

@inline _store!(out::Vector, k::Int, v) = (@inbounds out[k] = v; nothing)
@inline function _store!(outs::Tuple, k::Int, v::Tuple)
    map((o, x) -> (@inbounds o[k] = x), outs, v)
    return nothing
end

# The ring's values, gathered by the frame from the positions it minted.
@inline function _gather(data::AbstractVector,
        nbrs::SmallCollections.SmallVector{M,H}) where {M,H}
    out = SmallCollections.SmallVector{M,eltype(data)}()
    for h in nbrs
        out = SmallCollections.push(out, @inbounds data[h.position])
    end
    return out
end

@noinline _data_mismatch(nd, n) = throw(ArgumentError(
    "data must be laid out against the collection: expected a vector with " *
    "axis 1:$n, got axis $(nd)"))

_check_data(data::AbstractVector, n::Int) =
    eachindex(data) == Base.OneTo(n) || _data_mismatch(eachindex(data), n)

"""
    mapneighbors(f, cv; order = StorageOrder(), threaded = true,
                 connectivity = Vertex())
    mapneighbors(f, cv, data::AbstractVector; ...)

Apply `f` to every cell of the subset and its clipped one-ring, with the
frame owning the loop: traversal order, chunking, threading, and the writes.
`f` returns the cell's value — or a tuple of values — and the frame stores it
at the position it already knows, so a pure kernel never resolves a position
and never indexes an output. `cv` is a [`CellVector`](@ref), a
[`PartialGrid`](@ref) (iterated as its vector), or a [`CellLookup`](@ref).

Two calling forms:

  * `f(cell, nbrs)` — both arguments as [`SubsetPositionedCell`](@ref)
    handles, exactly the pairs the one-arg [`neighbors`](@ref) iterator
    yields.
  * with `data`, a vector laid out against the subset's axis:
    `f(cell, value, values)` — the cell's own entry and its neighbours'
    entries, gathered by the frame, aligned with the clipped ring. Pure
    local metrics read their samples without touching a handle at all.

Returns a `Vector` of `f`'s results in position order; a kernel whose
inferred return is a concrete tuple gets one vector per component,
returned as a tuple — the multi-output frame.

`order` is [`StorageOrder`](@ref) (the window walk, cursor at full
efficiency) or a permutation of `1:length(cv)` (e.g. `sortperm(elevation)`);
outputs land by position either way, so for a pure `f` the order changes
only the visit sequence. A non-permutation is an `ArgumentError`, since a
duplicate visit would silently overwrite one cell's output and leave
another's undefined.

`threaded` (`Bool` or GeometryOps' `True()`/`False()`) chunks the traversal
into contiguous ranges, one cursor per task, disjoint writes — legal
exactly when `f` is order-independent. Results are identical to the
sequential ones. [`foreachneighbors`](@ref) is the side-effecting form and
defaults to sequential.
"""
function mapneighbors(f::F, cv::CellVector; order = StorageOrder(),
        threaded = true, connectivity::Connectivity = Vertex()) where {F}
    M = max_neighbors(system(cv), connectivity)
    H = SubsetPositionedCell{eltype(cv)}
    T = Base.promote_op(f, H, SmallCollections.SmallVector{M,H})
    outs = _outputs(T, length(cv))
    return _mapstore!(f, outs, cv, connectivity, order, GOCore.booltype(threaded),
        Val(M))
end

function mapneighbors(f::F, cv::CellVector, data::AbstractVector;
        order = StorageOrder(), threaded = true,
        connectivity::Connectivity = Vertex()) where {F}
    _check_data(data, length(cv))
    M = max_neighbors(system(cv), connectivity)
    H = SubsetPositionedCell{eltype(cv)}
    T = Base.promote_op(f, H, eltype(data),
        SmallCollections.SmallVector{M,eltype(data)})
    outs = _outputs(T, length(cv))
    return _mapstore!(f, outs, cv, data, connectivity, order,
        GOCore.booltype(threaded), Val(M))
end

# Function barriers: `_outputs`' shape is runtime, so the store closure is
# built where the output type is concrete again.
function _mapstore!(f::F, outs::O, cv::CellVector, conn::Connectivity, order,
        thr, ::Val{M}) where {F,O,M}
    _run!((k, c, nbrs) -> _store!(outs, k, f(c, nbrs)), cv, conn, order, thr,
        Val(M))
    return outs
end

function _mapstore!(f::F, outs::O, cv::CellVector, data::AbstractVector,
        conn::Connectivity, order, thr, ::Val{M}) where {F,O,M}
    _run!((k, c, nbrs) -> _store!(outs, k,
            f(c, (@inbounds data[k]), _gather(data, nbrs))),
        cv, conn, order, thr, Val(M))
    return outs
end

"""
    foreachneighbors(f, cv; order = StorageOrder(), threaded = false,
                     connectivity = Vertex())
    foreachneighbors(f, cv, data::AbstractVector; ...)

[`mapneighbors`](@ref) without the output frame: `f` is called for its side
effects and its return is discarded. Same calling forms, same `order`
contract. `threaded` defaults OFF — a side-effecting `f` is order-dependent
until its author says otherwise — and turning it on is the caller asserting
order-independence, not the frame checking it.
"""
function foreachneighbors(f::F, cv::CellVector; order = StorageOrder(),
        threaded = false, connectivity::Connectivity = Vertex()) where {F}
    M = max_neighbors(system(cv), connectivity)
    _run!((k, c, nbrs) -> (f(c, nbrs); nothing), cv, connectivity, order,
        GOCore.booltype(threaded), Val(M))
    return nothing
end

function foreachneighbors(f::F, cv::CellVector, data::AbstractVector;
        order = StorageOrder(), threaded = false,
        connectivity::Connectivity = Vertex()) where {F}
    _check_data(data, length(cv))
    M = max_neighbors(system(cv), connectivity)
    _run!((k, c, nbrs) -> (f(c, (@inbounds data[k]), _gather(data, nbrs)); nothing),
        cv, connectivity, order, GOCore.booltype(threaded), Val(M))
    return nothing
end

mapneighbors(f::F, pg::PartialGrid; kw...) where {F} =
    mapneighbors(f, CellVector(pg); kw...)
mapneighbors(f::F, pg::PartialGrid, data::AbstractVector; kw...) where {F} =
    mapneighbors(f, CellVector(pg), data; kw...)
foreachneighbors(f::F, pg::PartialGrid; kw...) where {F} =
    foreachneighbors(f, CellVector(pg); kw...)
foreachneighbors(f::F, pg::PartialGrid, data::AbstractVector; kw...) where {F} =
    foreachneighbors(f, CellVector(pg), data; kw...)

# ===========================================================================
# The materialized sweep: the one-ring table in CSR
# ===========================================================================

"""
    HaloTable(cv; connectivity = Vertex(), threaded = true)
    HaloTable(pg::PartialGrid; connectivity = Vertex(), threaded = true)
    HaloTable(lk::CellLookup; connectivity = Vertex(), threaded = true)

The one-arg [`neighbors`](@ref) sweep materialized as positions, in CSR:
`t[p]` is the clipped one-ring of the subset's `p`-th cell as in-set
positions, **in the system's rotational ring order** — exactly the positions
the iterator's handles carry, row for row. Two output-sized arrays and
nothing else: `t.offsets` is `length(cv) + 1` long with row `p` occupying
`t.offsets[p] : t.offsets[p+1] - 1`, and `t.nbrs` is the flat neighbour
array those slices cut. `t[p]` is a `view`, so reading by row allocates
nothing; a kernel can equally loop `t.nbrs` directly with `t.offsets` as the
row bounds.

This is [`halo_table`](@ref)'s content in the layout a consumer that reads
every row wants: one flat array instead of a vector per cell, and ring order
instead of ascending — position `t[p][i]` is the cell's `i`-th surviving
ring member, so direction survives the materialization. `halo_table` remains
the ascending row-vector form, and the only one answering `k != 1`.

`threaded` (`Bool` or GeometryOps' `True()`/`False()`) builds the table in
[`mapneighbors`](@ref)' contiguous chunks, one cursor per task: each chunk
sweeps its range once into a local CSR, then one offset shift and one copy
per chunk stitch the pieces. Rows never interleave, so `offsets` and `nbrs`
come out identical to the sequential build's, whatever the chunking.
"""
struct HaloTable <: AbstractVector{SubArray{Int,1,Vector{Int},Tuple{UnitRange{Int}},true}}
    offsets::Vector{Int}
    nbrs::Vector{Int}
end

Base.size(t::HaloTable) = (length(t.offsets) - 1,)
Base.IndexStyle(::Type{HaloTable}) = Base.IndexLinear()

# Same shape and same reasoning as `StencilTable`'s indexing: the offset
# reads are safe by construction, the view stays checked because the fields
# are public.
Base.@propagate_inbounds function Base.getindex(t::HaloTable, p::Int)
    @boundscheck checkbounds(t, p)
    lo = @inbounds t.offsets[p]
    hi = @inbounds t.offsets[p+1] - 1
    return view(t.nbrs, lo:hi)
end

# The two arrays determine the rows and the rows determine the two arrays,
# so field equality is row equality without cutting a view per row.
Base.:(==)(a::HaloTable, b::HaloTable) =
    a.offsets == b.offsets && a.nbrs == b.nbrs

Base.show(io::IO, t::HaloTable) =
    print(io, "HaloTable(ncells=", length(t), ", entries=", length(t.nbrs), ")")
Base.show(io::IO, ::MIME"text/plain", t::HaloTable) = show(io, t)

HaloTable(cv::CellVector; connectivity::Connectivity = Vertex(),
    threaded = true) =
    _halo_table(cv, connectivity, GOCore.booltype(threaded))

function _halo_table(cv::CellVector, conn::Connectivity, ::GOCore.False)
    n = length(cv)
    offsets = Vector{Int}(undef, n + 1)
    @inbounds offsets[1] = 1
    nbrs = Int[]
    M = max_neighbors(system(cv), conn)
    # The degree ceiling as the append hint, as in `_stencil_rows`: the
    # output's own size to within the clip, asked for once.
    sizehint!(nbrs, n * M)
    _sweep!(cv, conn, 1:n, Val(M)) do k, c, ring
        for h in ring
            push!(nbrs, h.position)
        end
        @inbounds offsets[k+1] = length(nbrs) + 1
    end
    return HaloTable(offsets, nbrs)
end

# The chunked build. Each task sweeps its contiguous range with its own
# cursor into a LOCAL flat array, writing each of its rows' end positions
# into the global `offsets` as if its chunk began the table. Rings are
# computed once, here; what remains is arithmetic and copying.
function _halo_table(cv::CellVector, conn::Connectivity, ::GOCore.True)
    n = length(cv)
    n == 0 && return _halo_table(cv, conn, GOCore.False())
    M = max_neighbors(system(cv), conn)
    ranges = _chunk_ranges(n)
    parts = Vector{Vector{Int}}(undef, length(ranges))
    offsets = Vector{Int}(undef, n + 1)
    @inbounds offsets[1] = 1
    @sync for (i, r) in enumerate(ranges)
        Threads.@spawn @inbounds parts[i] = _halo_chunk(cv, conn, r, offsets,
            Val(M))
    end
    # One shift per chunk closes the seams: adding the entries that precede a
    # chunk turns its local end positions into global ones, and because the
    # ranges partition `1:n` in order, the result is the sequential build's
    # `offsets` exactly.
    total = 0
    @inbounds for (i, r) in enumerate(ranges)
        if total != 0
            for k in r
                offsets[k+1] += total
            end
        end
        total += length(parts[i])
    end
    # The stitch: the local arrays back to back, in range order.
    nbrs = Vector{Int}(undef, total)
    pos = 1
    for part in parts
        copyto!(nbrs, pos, part, 1, length(part))
        pos += length(part)
    end
    return HaloTable(offsets, nbrs)
end

function _halo_chunk(cv::CellVector, conn::Connectivity, r::UnitRange{Int},
        offsets::Vector{Int}, ::Val{M}) where {M}
    loc = Int[]
    sizehint!(loc, length(r) * M)
    _sweep!(cv, conn, r, Val(M)) do k, c, ring
        for h in ring
            push!(loc, h.position)
        end
        @inbounds offsets[k+1] = length(loc) + 1
    end
    return loc
end

HaloTable(pg::PartialGrid; kw...) = HaloTable(CellVector(pg); kw...)

# `halo_table`'s generic one-ring rows, produced by the sweep instead of one
# resolved call per cell: the same ascending `Vector{Int}` rows, with the
# per-neighbour window search amortized by the cursor. Every row is its own
# allocation landing in its own slot, so the chunked sweep needs no stitch:
# `thr` only decides how many cursors walk.
function _swept_rows(cv::CellVector, conn::Connectivity, thr = GOCore.False())
    n = length(cv)
    out = Vector{Vector{Int}}(undef, n)
    M = max_neighbors(system(cv), conn)
    _foreach_chunk(n, thr) do r
        _sweep!(cv, conn, r, Val(M)) do k, c, ring
            row = Vector{Int}(undef, length(ring))
            for (i, h) in enumerate(ring)
                @inbounds row[i] = h.position
            end
            @inbounds out[k] = sort!(row)
        end
    end
    return out
end
