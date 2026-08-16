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
