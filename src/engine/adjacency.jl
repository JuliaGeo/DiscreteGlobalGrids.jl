# ---------------------------------------------------------------------------
# The cached one-ring table.
#
# `border`/`interior`/`halo` answer about a region's boundary one cell at a
# time; this answers about every cell at once and keeps the answer. It is also
# the only verb here that threads: the lazy walks are sequential by contract,
# and a tens-of-millions-cell sweep is what this is for.
#
# Three row shapes, one builder. The shape decides what happens to a ring member
# the region does not hold — dropped, marked with 0, or addressed into a halo
# buffer — and each is a separate singleton so the inner loop stays monomorphic.
# ---------------------------------------------------------------------------

"""
    AdjacencyTable <: AbstractVector

The product of [`adjacency`](@ref): `t[p]` is the one-ring of in-region position
`p` as a non-allocating `view`, counter-clockwise, and `length(t)` is the region
size — the quantity an entry is compared against to tell a region slot
(`1:length(t)`) from a halo slot (`length(t) + j`).

CSR, and the two arrays are public: `t.offsets` is `length(t) + 1` long with row
`p` occupying `t.offsets[p] : t.offsets[p+1] - 1`, and `t.indices` is the flat
array those slices cut, so a kernel can loop over it directly instead of over
rows.

[`halopositions`](@ref) and [`halocells`](@ref) give the buffer's second half.
"""
struct AdjacencyTable{G} <: AbstractVector{SubArray{Int,1,Vector{Int},Tuple{UnitRange{Int}},true}}
    offsets::Vector{Int}
    indices::Vector{Int}
    nregion::Int
    halo::Vector{Int}
    grid::G
    mode::Symbol
end

Base.size(t::AdjacencyTable) = (t.nregion,)
Base.IndexStyle(::Type{<:AdjacencyTable}) = Base.IndexLinear()

# `offsets` is `nregion + 1` long by construction, so reading both bounds
# unchecked is safe once `p` is. The VIEW stays checked: both arrays are public
# fields, and one range comparison per row is nothing beside reading the row.
Base.@propagate_inbounds function Base.getindex(t::AdjacencyTable, p::Int)
    @boundscheck checkbounds(t, p)
    lo = @inbounds t.offsets[p]
    hi = @inbounds t.offsets[p+1] - 1
    return view(t.indices, lo:hi)
end

Base.:(==)(a::AdjacencyTable, b::AdjacencyTable) =
    a.offsets == b.offsets && a.indices == b.indices && a.halo == b.halo

Base.show(io::IO, t::AdjacencyTable) = print(io, "AdjacencyTable(halo = ",
    t.mode === :buffer ? 1 : t.mode === :mark ? ":mark" : 0, ", ncells=",
    t.nregion, ", halocells=", length(t.halo), ", entries=", length(t.indices), ")")
Base.show(io::IO, ::MIME"text/plain", t::AdjacencyTable) = show(io, t)

"""
    halopositions(t::AdjacencyTable) -> Vector{Int}
    halocells(t::AdjacencyTable) -> Vector{<:AbstractCellIndex}

The second half of the buffer `t`'s entries address, as complete-level positions
or as cell ids: slot `length(t) + j` of a row is element `j` of these, which is
element `j` of [`halo`](@ref)`(region)` under the table's connectivity.

Both are empty unless the table was built with `halo = 1` — the clipped and
marked forms name no cell outside the region, so they walk no halo. `halocells`
resolves ids from the positions on each call.
"""
halopositions(t::AdjacencyTable) = t.halo

@doc (@doc halopositions)
halocells(t::AdjacencyTable) = [cellindex(t.grid, p) for p in t.halo]

# ===========================================================================
# The three row shapes
# ===========================================================================

struct ClippedRows end
struct MarkedRows end

struct BufferedRows{V<:AbstractVector{<:Integer}}
    positions::V
    nregion::Int
end

@inline function _adj_push!(loc::Vector{Int}, ::ClippedRows, q::Int, ::Int,
        ::Any, ::Any)
    q == 0 || push!(loc, q)
    return nothing
end

# A position is one-based, so `0` is free as the out-of-region sentinel and the
# marked row is the clipped one with nothing dropped.
@inline function _adj_push!(loc::Vector{Int}, ::MarkedRows, q::Int, ::Int,
        ::Any, ::Any)
    push!(loc, q)
    return nothing
end

@inline function _adj_push!(loc::Vector{Int}, b::BufferedRows, q::Int, p::Int,
        c, nb)
    if q == 0
        j = Helpers.sorted_index(b.positions, p)
        j == 0 && _adjacency_incomplete(c, nb, p)
        push!(loc, b.nregion + j)
    else
        push!(loc, q)
    end
    return nothing
end

@noinline _adjacency_deep(n) = throw(ArgumentError(
    "adjacency: halo = $n is not available. A row exists only for an in-region " *
    "position, so a receptive field wider than one ring is a wider REGION: " *
    "build the table over `grow(region, $(n - 1))` and read the middle rows"))

@noinline _adjacency_mode(s) = throw(ArgumentError(
    "adjacency: halo must be 0, 1 or :mark, got :$s"))

@noinline _adjacency_unsorted() = throw(ArgumentError(
    "adjacency: the halo positions must be strictly ascending. That is the " *
    "order `halo(region)` emits in, and it is what the search addressing the " *
    "buffer's second half rests on — an unsorted list would silently " *
    "misaddress rows rather than fail"))

@noinline _adjacency_lazy() = throw(ArgumentError(
    "adjacency: the second argument is the halo as a materialised vector of " *
    "ascending positions, not the lazy walk. `collect(halo(region; " *
    "connectivity))` is that vector, and `adjacency(region; halo = 1)` walks " *
    "it for you"))

@noinline _adjacency_incomplete(c, nb, p) = throw(ArgumentError(
    "adjacency: $c has neighbour $nb at level position $p, which is in neither " *
    "the region nor the halo, so its row cannot be completed. The positions " *
    "passed are not this region's one-ring halo under this connectivity — the " *
    "usual causes are a halo walked under Vertex() and a table asked for under " *
    "Edge() or the reverse, a halo of a different region, or a truncated list"))

# ===========================================================================
# The builder
# ===========================================================================

_unwrap_task(e) = e
_unwrap_task(e::CompositeException) = _unwrap_task(first(e.exceptions))
_unwrap_task(e::TaskFailedException) = _unwrap_task(e.task.result)

# One pass in storage order, so the offsets fill as the rows are produced.
# Threaded chunks write local flat arrays and local row ends, both shifted by
# the number of preceding entries once every chunk is in.
function _adjacency_rows(cv::CellVector, conn::Connectivity, shape, thr)
    n = length(cv)
    offsets = Vector{Int}(undef, n + 1)
    @inbounds offsets[1] = 1
    n == 0 && return (offsets, Int[])
    cap = _capacity(system(cv), conn)
    ranges = thr isa GOCore.True ? _chunk_ranges(n) : [1:n]
    parts = Vector{Vector{Int}}(undef, length(ranges))
    if thr isa GOCore.True
        # A row that cannot be completed is the caller's error either way, so
        # the task wrapper is unwrapped rather than surfaced.
        try
            @sync for (i, r) in enumerate(ranges)
                Threads.@spawn @inbounds parts[i] =
                    _adjacency_chunk(cv, conn, r, shape, offsets, cap)
            end
        catch e
            throw(_unwrap_task(e))
        end
    else
        for (i, r) in enumerate(ranges)
            @inbounds parts[i] = _adjacency_chunk(cv, conn, r, shape, offsets, cap)
        end
    end
    total = 0
    @inbounds for (i, r) in enumerate(ranges)
        if total != 0
            for k in r
                offsets[k+1] += total
            end
        end
        total += length(parts[i])
    end
    indices = Vector{Int}(undef, total)
    pos = 1
    for part in parts
        copyto!(indices, pos, part, 1, length(part))
        pos += length(part)
    end
    return (offsets, indices)
end

function _adjacency_chunk(cv::CellVector, conn::Connectivity, r::UnitRange{Int},
        shape, offsets::Vector{Int}, cap::CAP) where {CAP}
    w = cv.windows
    wj = _window_at(w, first(r))
    hj = wj
    loc = Int[]
    sizehint!(loc, length(r) * _hint_degree(cap))
    for k in r
        wj = _advance(w, wj, k)
        c = cellindex(cv.grid, _leaf_at(w, wj, k))
        for nb in neighbors(cv.grid, c, 1; connectivity = conn)
            # Cannot be `nothing`: `nb` is a cell OF the complete level, which
            # is the grid it came from.
            p = cellposition(cv.grid, nb)::Int
            q, hj = _cursor_find(w, wj, hj, p)
            _adj_push!(loc, shape, q, p, c, nb)
        end
        @inbounds offsets[k+1] = length(loc) + 1
    end
    return loc
end

# ===========================================================================
# The verb
# ===========================================================================

function adjacency(region::Region; halo::Union{Integer,Symbol} = 0,
        connectivity::Connectivity = Vertex(), threaded = true)
    cv = CellVector(region)
    thr = GOCore.booltype(threaded)
    if halo isa Symbol
        halo === :mark || _adjacency_mode(halo)
        offsets, indices = _adjacency_rows(cv, connectivity, MarkedRows(), thr)
        return AdjacencyTable(offsets, indices, length(cv), Int[], cv.grid, :mark)
    end
    width = Int(halo)
    if width == 0
        offsets, indices = _adjacency_rows(cv, connectivity, ClippedRows(), thr)
        return AdjacencyTable(offsets, indices, length(cv), Int[], cv.grid, :clip)
    end
    width == 1 || _adjacency_deep(width)
    return adjacency(region, _halo_positions(region, connectivity);
        connectivity, threaded)
end

function adjacency(region::Region, hpos::AbstractVector{<:Integer};
        connectivity::Connectivity = Vertex(), threaded = true)
    Helpers.strictly_increasing(hpos) || _adjacency_unsorted()
    cv = CellVector(region)
    shape = BufferedRows(hpos, length(cv))
    offsets, indices = _adjacency_rows(cv, connectivity, shape,
        GOCore.booltype(threaded))
    return AdjacencyTable(offsets, indices, length(cv), collect(Int, hpos),
        cv.grid, :buffer)
end

# Passing the walk itself is the mistake the lazy design invites, so it is a
# message rather than a `MethodError`.
adjacency(::Region, ::Union{SubtreeHaloIterator,SubsetHaloIterator,
    HaloPositionIterator}; kw...) = _adjacency_lazy()
