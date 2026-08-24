# ---------------------------------------------------------------------------
# Iterate clipped one-rings in storage order while reusing index lookups.
# ---------------------------------------------------------------------------

"""
    SubsetIndexedCell{C}

A cell and its local index in the subset that created the handle. The one-argument
[`neighbors`](@ref) iterator yields these handles.

The index is valid only for that subset and is trusted unconditionally
wherever a handle is accepted: indexing a `Cells`-lookup raster with one reads
storage at that index directly, with no membership check. Use
[`cellid`](@ref) and resolve the bare cell when working with another
collection.

`==`, `hash` and `convert(::Type{C})` delegate to the cell, so a handle
compares and hashes as its cell; `show` prints the wrapper and the index
as well. [`localindex`](@ref)`(h)` returns the stored index.

The read-only cell verbs — [`cell_boundary`](@ref), [`cell_centroid`](@ref),
[`cell_polygon`](@ref), [`cell_area`](@ref), [`cell_extent`](@ref),
[`node_extent`](@ref), [`localindex`](@ref), [`globalindex`](@ref),
[`neighbors`](@ref), [`ring`](@ref), [`neighborcount`](@ref),
[`reindex`](@ref), [`level`](@ref) and [`rawid`](@ref) — accept a handle
wherever they accept a cell, and answer for the cell.
"""
struct SubsetIndexedCell{C<:AbstractCellIndex}
    cell::C
    index::Int
end

"""
    cellid(h::SubsetIndexedCell)
    cellid(c::AbstractCellIndex)

Return the bare cell from an indexed handle, or return a bare cell unchanged.
The escape from a handle's single-axis contract: `localindex(other,
cellid(h))` resolves against any axis, where the handle's own index is valid
only against the collection that minted it.
"""
cellid(h::SubsetIndexedCell) = h.cell
cellid(c::AbstractCellIndex) = c

"""
    localindex(h::SubsetIndexedCell) -> Int

Return the handle's local index in the collection that created it.
"""
localindex(h::SubsetIndexedCell) = h.index

level(h::SubsetIndexedCell) = level(h.cell)
rawid(h::SubsetIndexedCell) = rawid(h.cell)

# The read-only cell verbs answer for the underlying cell, so a handle works
# wherever the cell it prints does. The collection stays untyped because grids,
# systems, `CellVector`s and `CellLookup`s all take a cell in this slot.
localindex(x, h::SubsetIndexedCell) = localindex(x, h.cell)
globalindex(x, h::SubsetIndexedCell) = globalindex(x, h.cell)
cell_boundary(x, h::SubsetIndexedCell) = cell_boundary(x, h.cell)
cell_centroid(x, h::SubsetIndexedCell) = cell_centroid(x, h.cell)
cell_polygon(x, h::SubsetIndexedCell) = cell_polygon(x, h.cell)
cell_area(x, h::SubsetIndexedCell) = cell_area(x, h.cell)
cell_extent(x, h::SubsetIndexedCell) = cell_extent(x, h.cell)
node_extent(x, h::SubsetIndexedCell) = node_extent(x, h.cell)
neighbors(x, h::SubsetIndexedCell, k::Integer = 1;
    connectivity::Connectivity = Vertex()) = neighbors(x, h.cell, k; connectivity)
ring(x, h::SubsetIndexedCell, k::Integer;
    connectivity::Connectivity = Vertex()) = ring(x, h.cell, k; connectivity)
neighborcount(x, h::SubsetIndexedCell;
    connectivity::Connectivity = Vertex()) = neighborcount(x, h.cell; connectivity)
reindex(T::Type{<:AbstractCellIndex}, sys, h::SubsetIndexedCell) =
    reindex(T, sys, h.cell)

Base.:(==)(a::SubsetIndexedCell, b::SubsetIndexedCell) = a.cell == b.cell
Base.:(==)(a::SubsetIndexedCell, b::AbstractCellIndex) = a.cell == b
Base.:(==)(a::AbstractCellIndex, b::SubsetIndexedCell) = a == b.cell
Base.hash(h::SubsetIndexedCell, u::UInt) = hash(h.cell, u)
# Print what it is: a bare cell here would hide both the wrapper and the
# index, which is the pair every handle question turns on.
Base.show(io::IO, h::SubsetIndexedCell) =
    print(io, "SubsetIndexedCell(", h.cell, ", index ", h.index, ")")
Base.show(io::IO, ::MIME"text/plain", h::SubsetIndexedCell) = show(io, h)
Base.convert(::Type{C}, h::SubsetIndexedCell{C}) where {C<:AbstractCellIndex} =
    h.cell

# ===========================================================================
# The window cursor
#
# The cursor tracks the storage index, its window, and the last window that
# contained a neighbour. A returned index of zero means the cell is absent.
# ===========================================================================

@inline _wbase(w::RangeWindows, j::Int) = j == 1 ? 0 : @inbounds w.offsets[j-1]

@inline function _advance(w::RangeWindows, wj::Int, k::Int)
    @inbounds while k > w.offsets[wj]
        wj += 1
    end
    return wj
end
@inline _advance(::IndexWindows, wj::Int, k::Int) = k

@inline _leaf_at(w::RangeWindows, wj::Int, k::Int) =
    @inbounds w.starts[wj] + (k - _wbase(w, wj)) - 1
@inline _leaf_at(w::IndexWindows, ::Int, k::Int) = @inbounds w.indices[k]

# Leaf index -> (concatenation index or 0, updated hint).
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

# Check the previous match before searching an index list.
@inline function _cursor_find(w::IndexWindows, ::Int, hj::Int, p::Int)
    ps = w.indices
    @inbounds if hj <= length(ps) && ps[hj] == p
        return hj, hj
    end
    j = searchsortedfirst(ps, p)
    (j <= length(ps) && @inbounds(ps[j]) == p) || return 0, hj
    return j, j
end

# ===========================================================================
# Ring buffer capacity
#
# `maxneighbors` is the system's static degree bound. A system that declares
# one buys stack-allocated `SmallVector` rings; a system that declares none
# (the trait defaults to `nothing`) gets heap `Vector` rings with the same
# contents and the same order. The capacity travels as a singleton VALUE —
# `Val(M)` or `nothing` — so every function below specializes on which of the
# two containers it is building, and the declared path stays exactly what it
# was.
# ===========================================================================

_capacity(sys, conn::Connectivity) =
    static_capacity(maxneighbors(sys, conn), cellindextype(sys))

_capacity(M, ::Type{T}) where {T} = static_capacity(M, T)

_ringtype(::Val{M}, ::Type{T}) where {M,T} = SmallCollections.SmallVector{M,T}
_ringtype(::Nothing, ::Type{T}) where {T} = Vector{T}

@inline _newring(::Val{M}, ::Type{T}) where {M,T} =
    SmallCollections.SmallVector{M,T}()

@inline _pushring(v::SmallCollections.SmallVector, x) =
    SmallCollections.push(v, x)

# Row-length guess for `sizehint!`. With no declared bound this is only a
# reservation hint, never a capacity: the buffers grow if it is wrong.
_hint_degree(::Val{M}) where {M} = M
_hint_degree(::Nothing) = 8

# Clip a ring to the subset and attach each surviving cell's local index. A
# declared degree bound keeps the result statically sized and off the heap.
@inline function _indexed(cv::CellVector, wj::Int, hj::Int, cells, cap)
    w = cv.windows
    out = _newring(cap, SubsetIndexedCell{eltype(cells)})
    for nb in cells
        q, hj = _cursor_find(w, wj, hj, globalindex(cv.grid, nb)::Int)
        q == 0 || (out = _pushring(out, SubsetIndexedCell(nb, q)))
    end
    return out, hj
end

# Undeclared: the unclipped ring is an exact upper bound on the clipped one, so
# the row is allocated once at that size and trimmed. `push!`-from-empty would
# reallocate on the way to a degree-12 neighbourhood.
@inline function _indexed(cv::CellVector, wj::Int, hj::Int, cells, ::Nothing)
    w = cv.windows
    out = Vector{SubsetIndexedCell{eltype(cells)}}(undef, length(cells))
    m = 0
    for nb in cells
        q, hj = _cursor_find(w, wj, hj, globalindex(cv.grid, nb)::Int)
        if q != 0
            m += 1
            @inbounds out[m] = SubsetIndexedCell(nb, q)
        end
    end
    m == length(out) || resize!(out, m)
    return out, hj
end

# ===========================================================================
# The iterator
# ===========================================================================

"""
    neighbors(cv::CellVector; connectivity = Vertex())
    neighbors(pg::PartialGrid; connectivity = Vertex())

Iterate a subset in storage order as `(cell, nbrs)`. Both the cell and its
one-ring, clipped to membership, are [`SubsetIndexedCell`](@ref) handles:
`nbrs` holds the same cells `neighbors(cv, cell)` answers, in the same order.

The iterator reuses window lookups across the sweep and allocates nothing
beyond the system's native one-ring operation. A `PartialGrid` uses its
`CellVector` representation without changing indices.
"""
neighbors(cv::CellVector; connectivity::Connectivity = Vertex()) =
    NeighborhoodIterator(cv, connectivity)

neighbors(pg::PartialGrid; connectivity::Connectivity = Vertex()) =
    NeighborhoodIterator(CellVector(pg), connectivity)

# `CAP` is the ring capacity witness: `Val{M}` for a declared bound, `Nothing`
# for a system that declares none.
struct NeighborhoodIterator{CAP,CV<:CellVector,CN<:Connectivity}
    cv::CV
    connectivity::CN
end

function NeighborhoodIterator(cv::CellVector, conn::Connectivity)
    cap = _capacity(system(cv), conn)
    return NeighborhoodIterator{typeof(cap),typeof(cv),typeof(conn)}(cv, conn)
end

Base.length(it::NeighborhoodIterator) = length(it.cv)
Base.IteratorSize(::Type{<:NeighborhoodIterator}) = Base.HasLength()
Base.IteratorEltype(::Type{<:NeighborhoodIterator}) = Base.EltypeUnknown()

Base.show(io::IO, it::NeighborhoodIterator) =
    print(io, "neighbors(", it.cv, "; connectivity=", it.connectivity, ")")
Base.show(io::IO, ::MIME"text/plain", it::NeighborhoodIterator) = show(io, it)

Base.iterate(it::NeighborhoodIterator) = _neighborhood_next(it, 1, 1, 1)
Base.iterate(it::NeighborhoodIterator, st::NTuple{3,Int}) =
    _neighborhood_next(it, st[1], st[2], st[3])

function _neighborhood_next(it::NeighborhoodIterator{CAP}, k::Int, wj::Int,
        hj::Int) where {CAP}
    cv = it.cv
    k > length(cv) && return nothing
    w = cv.windows
    wj = _advance(w, wj, k)
    c = cellindex(cv.grid, _leaf_at(w, wj, k))
    nbrs, hj = _indexed(cv, wj, hj,
        neighbors(cv.grid, c, 1; connectivity = it.connectivity), CAP())
    return (SubsetIndexedCell(c, k), nbrs), (k + 1, wj, hj)
end

# ===========================================================================
# Closure-based neighbourhood sweeps.
# ===========================================================================

# Find the window containing storage index `k`.
@inline _window_at(w::RangeWindows, k::Int) = searchsortedfirst(w.offsets, k)
@inline _window_at(::IndexWindows, k::Int) = k

"""
    StorageOrder()

Traverse cells in storage order. This is the default for
[`mapneighbors`](@ref) and [`foreachneighbors`](@ref).
"""
struct StorageOrder end

# Call `g(k, cell, nbrs)` in storage order with one cursor for the range.
function _sweep!(g::G, cv::CellVector, conn::Connectivity, r::UnitRange{Int},
        cap::CAP) where {G,CAP}
    isempty(r) && return nothing
    w = cv.windows
    wj = _window_at(w, first(r))
    hj = wj
    for k in r
        wj = _advance(w, wj, k)
        c = cellindex(cv.grid, _leaf_at(w, wj, k))
        nbrs, hj = _indexed(cv, wj, hj,
            neighbors(cv.grid, c, 1; connectivity = conn), cap)
        g(k, SubsetIndexedCell(c, k), nbrs)
    end
    return nothing
end

# Visit the indices selected by `perm`; each visit starts with a fresh window
# lookup because permutations provide no locality guarantee.
function _sweep_perm!(g::G, cv::CellVector, conn::Connectivity,
        perm::AbstractVector{<:Integer}, r::UnitRange{Int},
        cap::CAP) where {G,CAP}
    w = cv.windows
    for j in r
        k = Int(@inbounds perm[j])
        wj = _window_at(w, k)
        c = cellindex(cv.grid, _leaf_at(w, wj, k))
        nbrs, _ = _indexed(cv, wj, wj,
            neighbors(cv.grid, c, 1; connectivity = conn), cap)
        g(k, SubsetIndexedCell(c, k), nbrs)
    end
    return nothing
end

# Reject orders that do not visit every output index exactly once.
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

# Split `1:n` into contiguous ranges, one cursor per task.
function _chunk_ranges(n::Int)
    step = cld(n, min(n, 2 * Threads.nthreads()))
    return [lo:min(lo + step - 1, n) for lo in 1:step:n]
end

"""
    NeighborCallbackError

A [`mapneighbors`](@ref) or [`foreachneighbors`](@ref) callback that threw
during a threaded sweep, naming the `cell` and subset `index` it was called
with. The callback's own exception is `err` and is shown as the cause.

One callback failure raises one of these: the sweep waits for every chunk, then
reports the failure at the lowest index and drops the rest. The sequential
path lets the callback's exception through untouched.
"""
struct NeighborCallbackError <: Exception
    cell::AbstractCellIndex
    index::Int
    err::Any
    backtrace::Any
end

function Base.showerror(io::IO, e::NeighborCallbackError)
    print(io, "NeighborCallbackError: the callback failed at cell ", e.cell,
        ", subset index ", e.index, ".",
        "\nRerun with `threaded = false` to raise it on its own.",
        "\ncaused by: ")
    showerror(io, e.err, e.backtrace)
    return nothing
end

# Name the cell a failing callback was called with. Only the threaded path is
# wrapped, so the sequential sweep keeps its bare stack trace.
struct _ReportingCallback{G}
    g::G
end

@inline function (cb::_ReportingCallback)(k::Int, c, nbrs)
    try
        return cb.g(k, c, nbrs)
    catch err
        err isa NeighborCallbackError && rethrow()
        throw(NeighborCallbackError(cellid(c), k, err, catch_backtrace()))
    end
end

_reporting(g::G, ::GOCore.False) where {G} = g
_reporting(g::G, ::GOCore.True) where {G} = _ReportingCallback(g)

# Every chunk is waited on before anything is reported, so no task outlives the
# call and the reported failure is the earliest in index order rather than
# the first to be scheduled.
function _foreach_chunk(body!::F, n::Int, ::GOCore.True) where {F}
    n == 0 && return nothing
    ranges = _chunk_ranges(n)
    tasks = Vector{Task}(undef, length(ranges))
    for (i, r) in enumerate(ranges)
        @inbounds tasks[i] = Threads.@spawn body!(r)
    end
    failed = 0
    for (i, t) in enumerate(tasks)
        try
            wait(t)
        catch
            failed == 0 && (failed = i)
        end
    end
    failed == 0 && return nothing
    return _report_chunk_failure(@inbounds tasks[failed])
end

@noinline function _report_chunk_failure(t::Task)
    stack = Base.current_exceptions(t)
    if !isempty(stack)
        err = last(stack).exception
        err isa NeighborCallbackError && throw(err)
    end
    # Machinery failures are not the callback's to explain.
    return wait(t)
end

_foreach_chunk(body!::F, n::Int, ::GOCore.False) where {F} = body!(1:n)

function _run!(g::G, cv::CellVector, conn::Connectivity, ::StorageOrder, thr,
        cap::CAP) where {G,CAP}
    h = _reporting(g, thr)
    return _foreach_chunk(r -> _sweep!(h, cv, conn, r, cap), length(cv), thr)
end

function _run!(g::G, cv::CellVector, conn::Connectivity,
        perm::AbstractVector{<:Integer}, thr, cap::CAP) where {G,CAP}
    _check_permutation(perm, length(cv))
    h = _reporting(g, thr)
    return _foreach_chunk(r -> _sweep_perm!(h, cv, conn, perm, r, cap),
        length(perm), thr)
end

@noinline _run!(g, cv::CellVector, conn, order, thr, v) = throw(ArgumentError(
    "order must be StorageOrder() or a permutation of 1:length(cv), " *
    "got $(typeof(order))"))

# --- output storage ---------------------------------------------------------

# Concrete tuple results use one output vector per component.
_outputs(::Type{T}, n::Int) where {T} =
    T <: Tuple && isconcretetype(T) ?
    ntuple(j -> Vector{fieldtype(T, j)}(undef, n), fieldcount(T)) :
    Vector{T}(undef, n)

@inline _store!(out::Vector, k::Int, v) = (@inbounds out[k] = v; nothing)
@inline function _store!(outs::Tuple, k::Int, v::Tuple)
    map((o, x) -> (@inbounds o[k] = x), outs, v)
    return nothing
end

# Gather neighbour values in ring order, into the same container shape the
# handles arrived in.
@inline function _gather(data::AbstractVector,
        nbrs::SmallCollections.SmallVector{M}) where {M}
    out = SmallCollections.SmallVector{M,eltype(data)}()
    for h in nbrs
        out = SmallCollections.push(out, @inbounds data[h.index])
    end
    return out
end

@inline function _gather(data::AbstractVector, nbrs::Vector)
    out = Vector{eltype(data)}(undef, length(nbrs))
    for (i, h) in enumerate(nbrs)
        @inbounds out[i] = data[h.index]
    end
    return out
end

@noinline _data_mismatch(nd, n) = throw(ArgumentError(
    "data must be laid out against the collection: expected a vector with " *
    "axis 1:$n, got axis $(nd)"))

# `axes`, not `eachindex`: what the contract asks is that index `k` of the
# data is index `k` of the collection, and a lazy array satisfies that while
# reporting a chunked `eachindex` that is not a `OneTo`. Passing one is legal
# and slow — see `foreachchunk` for the traversal that makes it fast.
_check_data(data::AbstractVector, n::Int) =
    axes(data) == (Base.OneTo(n),) || _data_mismatch(axes(data, 1), n)

"""
    mapneighbors(f, cv; order = StorageOrder(), threaded = true,
                 connectivity = Vertex())
    mapneighbors(f, cv, data::AbstractVector; ...)

Apply `f` to each cell and its clipped one-ring. `cv` may be a
[`CellVector`](@ref), [`PartialGrid`](@ref), or [`CellLookup`](@ref).

Without `data`, `f(cell, nbrs)` receives the same indexed handles yielded
by the one-argument [`neighbors`](@ref) iterator. With a vector laid out
against the subset, `f(cell, value, values)` receives the cell value and its
neighbour values in the counter-clockwise order [`neighbors`](@ref) states, so
slot `j` of the callback's ring names a direction.

Results are stored in subset index order. A concrete tuple result produces
a tuple of vectors, one per component. `order` accepts [`StorageOrder`](@ref)
or a permutation of `1:length(cv)`; invalid permutations throw
`ArgumentError`.

When `threaded` is true, contiguous ranges run in separate tasks and write to
disjoint output indices — legal exactly when `f` is order-independent, and
the results are then identical to the sequential ones. A callback that throws
there raises one [`NeighborCallbackError`](@ref) naming the cell and index
it failed at, not one exception per task.
[`foreachneighbors`](@ref) provides the side-effecting form and defaults to
sequential execution.
"""
function mapneighbors(f::F, cv::CellVector; order = StorageOrder(),
        threaded = true, connectivity::Connectivity = Vertex()) where {F}
    cap = _capacity(system(cv), connectivity)
    H = SubsetIndexedCell{eltype(cv)}
    T = Base.promote_op(f, H, _ringtype(cap, H))
    outs = _outputs(T, length(cv))
    return _mapstore!(f, outs, cv, connectivity, order, GOCore.booltype(threaded),
        cap)
end

function mapneighbors(f::F, cv::CellVector, data::AbstractVector;
        order = StorageOrder(), threaded = true,
        connectivity::Connectivity = Vertex()) where {F}
    _check_data(data, length(cv))
    cap = _capacity(system(cv), connectivity)
    H = SubsetIndexedCell{eltype(cv)}
    T = Base.promote_op(f, H, eltype(data), _ringtype(cap, eltype(data)))
    outs = _outputs(T, length(cv))
    return _mapstore!(f, outs, cv, data, connectivity, order,
        GOCore.booltype(threaded), cap)
end

# Build the store closure after the output container type is known.
function _mapstore!(f::F, outs::O, cv::CellVector, conn::Connectivity, order,
        thr, cap::CAP) where {F,O,CAP}
    _run!((k, c, nbrs) -> _store!(outs, k, f(c, nbrs)), cv, conn, order, thr,
        cap)
    return outs
end

function _mapstore!(f::F, outs::O, cv::CellVector, data::AbstractVector,
        conn::Connectivity, order, thr, cap::CAP) where {F,O,CAP}
    _run!((k, c, nbrs) -> _store!(outs, k,
            f(c, (@inbounds data[k]), _gather(data, nbrs))),
        cv, conn, order, thr, cap)
    return outs
end

"""
    foreachneighbors(f, cv; order = StorageOrder(), threaded = false,
                     connectivity = Vertex())
    foreachneighbors(f, cv, data::AbstractVector; ...)

Call `f` for each cell and clipped one-ring, discarding its return value. The
calling forms and `order` contract match [`mapneighbors`](@ref). Threading is
disabled by default; enabling it requires `f` to be order-independent.
"""
function foreachneighbors(f::F, cv::CellVector; order = StorageOrder(),
        threaded = false, connectivity::Connectivity = Vertex()) where {F}
    _run!((k, c, nbrs) -> (f(c, nbrs); nothing), cv, connectivity, order,
        GOCore.booltype(threaded), _capacity(system(cv), connectivity))
    return nothing
end

function foreachneighbors(f::F, cv::CellVector, data::AbstractVector;
        order = StorageOrder(), threaded = false,
        connectivity::Connectivity = Vertex()) where {F}
    _check_data(data, length(cv))
    _run!((k, c, nbrs) -> (f(c, (@inbounds data[k]), _gather(data, nbrs)); nothing),
        cv, connectivity, order, GOCore.booltype(threaded),
        _capacity(system(cv), connectivity))
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
