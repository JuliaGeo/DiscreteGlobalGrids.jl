# ---------------------------------------------------------------------------
# Lazy subtree id vectors, and neighbor stepping over them
#
# `DGGSGlobeIds` is the whole sphere at one level, computed instead of stored.
# A tile is the other half of the same idea: the complete set of `level`
# descendants of one cell, which is what a regional workflow actually holds —
# a DEM over one catchment, a regridding destination, a stencil domain.
#
# Two structural facts make it a *vector* rather than a set:
#
#   * `descendant_range` is a two-sided interval on the encoded id space, so
#     membership of a leaf-level id is two comparisons; and
#   * `cell_to_ordinal` is strictly monotone in the id (`kernel.jl`, "Dense
#     ordinals"), so a subtree's descendants occupy a *contiguous* ordinal
#     interval. Position within the tile is therefore arithmetic, not a search:
#     `cell_to_ordinal(id) - offset`.
#
# The second fact is the one that matters for stencils. A partial lookup
# resolves a neighbor id to an array position with `searchsortedfirst` over the
# stored ids — O(log n) with a random memory probe per neighbor per cell. Over
# a tile it is a fixed-cost integer computation with no memory traffic at all.
#
# On top of that sits the neighbor stepper: the per-tile object a stencil
# sweep asks for a cell's neighbors. It exists because the three ways to answer
# have wildly different cost profiles and the right choice is not universal —
# see `AbstractNeighborStepper` below.
# ---------------------------------------------------------------------------

"""
    DGGSSubtreeIds(system, root_level, root_id, level) -> AbstractVector{cell_id_type(system)}

Every `level` descendant of cell `(root_level, root_id)`, in ascending
canonical-id order, computed on demand instead of stored — the regional
counterpart of [`DGGSGlobeIds`](@ref).

The whole vector is a handful of words however many cells it names, so it is
what makes a tile-shaped `<X>Lookup` (and hence a tile-shaped `DimArray`
dimension) free to build:

```julia
tile = DGGSSubtreeIds(IGEO7DGGS(), 5, 0x0c4d9fffffffffff, 13)   # 5.76e6 cells
lookup = IGeo7Lookup(tile, 13, Dict{String,Any}())              # O(1)
```

Requires [`has_descendant_ranges`](@ref); systems without it (A5) have no
two-sided interval to test membership against and throw here.

Beyond indexing it answers two things in O(1), both used by the stencil path:

  * `id in tile` — membership, two comparisons against the
    [`descendant_range`](@ref) endpoints;
  * [`subtree_position`](@ref)`(tile, id)` — the array position of a cell,
    from the ordinal contract's monotonicity, with no search.
"""
struct DGGSSubtreeIds{S<:AbstractDGGS,ID<:Integer} <: AbstractVector{ID}
    system::S
    root_level::Int
    root_id::ID
    level::Int
    lo::ID              # first leaf descendant, `descendant_range` low end
    hi::ID              # last leaf descendant, high end
    offset::Int         # ordinal of `lo`, minus one
    n::Int

    function DGGSSubtreeIds(system::AbstractDGGS, root_level::Integer, root_id,
            level::Integer)
        rl, lvl = Int(root_level), Int(level)
        0 <= rl <= lvl || throw(ArgumentError("expected 0 <= root_level <= level"))
        limit = max_level(system)
        limit === nothing || lvl <= limit ||
            throw(ArgumentError("$(system_name(system)) level must be in 0:$limit"))
        has_descendant_ranges(system) || throw(ArgumentError(
            "$(system_name(system)) has no descendant ranges, so a subtree is not " *
            "an id interval; build the ids explicitly with cell_descendants"))
        ID = cell_id_type(system)
        lo, hi = descendant_range(system, rl, root_id, lvl)
        o1 = Int(cell_to_ordinal(system, lvl, lo))
        o2 = Int(cell_to_ordinal(system, lvl, hi))
        return new{typeof(system),ID}(system, rl, ID(root_id), lvl, ID(lo), ID(hi),
            o1 - 1, o2 - o1 + 1)
    end
end

Base.size(t::DGGSSubtreeIds) = (t.n,)
Base.IndexStyle(::Type{<:DGGSSubtreeIds}) = Base.IndexLinear()
Base.getindex(t::DGGSSubtreeIds, i::Int) = ordinal_to_cell(t.system, t.level, t.offset + i)

# Same reasoning as `Helpers.strictly_increasing(::DGGSGlobeIds)`: ascending
# ordinals enumerate ascending ids by the ordinal contract, so the lookups'
# O(n) verification pass has nothing to discover.
Helpers.strictly_increasing(::DGGSSubtreeIds) = true

Base.show(io::IO, t::DGGSSubtreeIds) = print(io, "DGGSSubtreeIds(", t.system, ", ",
    t.root_level, ", ", repr(t.root_id), ", ", t.level, ")")
Base.show(io::IO, ::MIME"text/plain", t::DGGSSubtreeIds) = show(io, t)

"""
    Base.in(id, tile::DGGSSubtreeIds) -> Bool

`true` when the **`tile.level`** id `id` is one of the tile's cells. Two
comparisons against the [`descendant_range`](@ref) endpoints — restricted to one
level that interval *is* the tile, which is why this is not a search.

The level restriction is load-bearing, and it is not symmetric. A coarser id is
correctly rejected: an ancestor's encoding pads with the sentinel digit, which
sorts above every active digit, so it sorts above `hi`. A **finer** id is not —
a descendant of one of the tile's own cells shares the prefix and lands inside
the interval, so it tests `true` for a cell the tile does not contain. This is
[`descendant_range`](@ref)'s documented shape (a prefix interval, not a cell
set) surfacing; the callers here — [`subtree_position`](@ref) and the neighbor
steppers — only ever present ids that [`cell_neighbors`](@ref) produced at
`tile.level`, so the restriction costs them nothing.
"""
@inline Base.in(id, t::DGGSSubtreeIds) = t.lo <= id <= t.hi

"""
    subtree_position(tile, id) -> Int

The position of leaf-level `id` within `tile`, or `0` when it is outside.
Arithmetic, not a search: descendants occupy a contiguous ordinal interval, so
this is [`cell_to_ordinal`](@ref) minus the tile's offset.
"""
@inline function subtree_position(t::DGGSSubtreeIds, id)
    id in t || return 0
    return Int(cell_to_ordinal(t.system, t.level, id)) - t.offset
end

# --------------------------------------------------------------------------
# Rim and interior
# --------------------------------------------------------------------------

"""
    subtree_border_positions(tile) -> Vector{Int}

Positions of the tile cells with at least one edge neighbor *outside* the
tile — [`subtree_border`](@ref) mapped through [`subtree_position`](@ref),
ascending.

This is the set a stencil cannot evaluate from the tile alone: its neighborhood
is truncated, so a reduction over it sees fewer values than an interior cell
does. Naming it is the point — an unnoticed truncated neighborhood is how a
flow-routing sweep grows spurious pits along a tile edge.

The materialized form of [`edge_cells`](@ref); the rim is small enough that
the difference rarely matters, unlike the interior's.
"""
subtree_border_positions(t::DGGSSubtreeIds) = collect(edge_cells(t))

"""
    subtree_interior_positions(tile, border=subtree_border_positions(tile)) -> Vector{Int}

Positions of the tile cells whose whole edge neighborhood is inside the tile —
the complement of [`subtree_border_positions`](@ref). Pass `border` back in if
it is already computed.

This materializes the interior, which is almost the whole tile; prefer
[`interior_cells`](@ref), which is the same sequence without the array.
"""
subtree_interior_positions(t::DGGSSubtreeIds,
    border::AbstractVector{Int}=subtree_border_positions(t)) =
    collect(interior_cells(t, border))

# --------------------------------------------------------------------------
# Rim and interior, without the arrays
#
# The two sets are wildly different sizes, and the useful laziness follows the
# asymmetry rather than treating them alike.
#
# The rim is small — `3^(d+1) - 3` cells against the `7^d` a hexagon subtree
# holds, 2.4% of it at depth 8 — and `subtree_border` already enumerates it in
# `O(result)`. So there is nothing to save by not holding it; what
# `edge_cells` drops is only the *second* array, mapping ids to positions on
# access instead of into a fresh `Vector{Int}`.
#
# The interior is the other 97.6%, and materializing it is close to writing
# down `1:n` — 6.4 MB of `Int` at depth 8 to say "nearly everything". That one
# is worth not building: the complement of an ascending set is a merge walk
# with no state beyond a cursor, and its `k`-th element is a binary search.
#
# Both are `AbstractVector{Int}`, not just iterators, so they drop into
# `subtree_stencil`'s `border` / `interior` keywords (which is where they are
# now the defaults) and into anything else expecting an indexable position
# list.
# --------------------------------------------------------------------------

"""
    edge_cells(tile) -> AbstractVector{Int}

Positions of the tile cells with at least one edge neighbor *outside* the
tile — the rim — ascending, computed lazily from
[`subtree_border`](@ref)'s ids.

This is the set a stencil cannot evaluate from the tile alone: its
neighborhood is truncated, so a reduction over it sees fewer values than an
interior cell does. Naming it is the point — an unnoticed truncated
neighborhood is how a flow-routing sweep grows spurious pits along a tile edge.

`tile[i]` recovers the cell id at position `i`.

    for i in edge_cells(tile)
        halo_exchange!(data, i)
    end

See [`interior_cells`](@ref) for the complement, and
[`subtree_border_positions`](@ref) for the materialized form.
"""
struct SubtreeEdgeCells{T<:DGGSSubtreeIds,V<:AbstractVector} <: AbstractVector{Int}
    tile::T
    ids::V
end

edge_cells(t::DGGSSubtreeIds) =
    SubtreeEdgeCells(t, subtree_border(t.system, t.root_level, t.root_id, t.level))

Base.size(e::SubtreeEdgeCells) = (length(e.ids),)
Base.IndexStyle(::Type{<:SubtreeEdgeCells}) = IndexLinear()

# Ascending without a sort: ids come out ascending and ordinals are monotone
# in the id, so the positions inherit the order.
Base.@propagate_inbounds function Base.getindex(e::SubtreeEdgeCells, k::Int)
    @boundscheck checkbounds(e, k)
    return subtree_position(e.tile, @inbounds e.ids[k])
end

"""
    interior_cells(tile, edge=edge_cells(tile)) -> AbstractVector{Int}

Positions of the tile cells whose whole edge neighborhood is inside the tile —
the complement of [`edge_cells`](@ref) — ascending, allocating nothing beyond
`edge` itself. These are the cells a stencil can evaluate from tile-local data
alone, so their loop needs no membership branch and their neighbor container
is always full.

Pass `edge` back in when it is already computed; the two are usually wanted
together.

    for i in interior_cells(tile)
        data[i] = f(data[i])
    end

Iteration is a merge against `edge` (`O(1)` per cell); `getindex` is a binary
search over it (`O(log |edge|)`).
"""
struct SubtreeInteriorCells{T<:DGGSSubtreeIds,E<:AbstractVector{Int}} <: AbstractVector{Int}
    tile::T
    edge::E
    n::Int          # length(tile)
    m::Int          # length(edge)
end

function interior_cells(t::DGGSSubtreeIds, edge::AbstractVector{Int}=edge_cells(t))
    return SubtreeInteriorCells(t, edge, length(t), length(edge))
end

Base.size(v::SubtreeInteriorCells) = (v.n - v.m,)
Base.IndexStyle(::Type{<:SubtreeInteriorCells}) = IndexLinear()

# The k-th position that `edge` does not contain. With `edge` ascending and
# distinct, `g(j) = edge[j] - j` is non-decreasing and counts the gap before
# `edge[j]`; the answer sits after exactly `j*` edge entries, where `j*` is the
# last `j` with `g(j) < k`. Then the answer is `k + j*`.
Base.@propagate_inbounds function Base.getindex(v::SubtreeInteriorCells, k::Int)
    @boundscheck checkbounds(v, k)
    lo, hi = 0, v.m                  # invariant: g(lo) < k, and lo is feasible
    @inbounds while lo < hi
        mid = (lo + hi + 1) >>> 1
        if v.edge[mid] - mid < k
            lo = mid
        else
            hi = mid - 1
        end
    end
    return k + lo
end

# Iteration walks the complement directly rather than re-searching per element:
# `(next position, how many edge entries are already behind it)`.
@inline function Base.iterate(v::SubtreeInteriorCells, state::Tuple{Int,Int}=(1, 1))
    pos, j = state
    @inbounds while j <= v.m && v.edge[j] == pos
        pos += 1
        j += 1
    end
    pos > v.n && return nothing
    return pos, (pos + 1, j)
end

# --------------------------------------------------------------------------
# Neighbor steppers
#
# Three ways to answer "which positions are cell `i`'s neighbors", with
# genuinely different cost profiles — hence a small type hierarchy rather than
# one function:
#
#   generic   `cell_neighbors` + `subtree_position`. Correct for every wired
#             system. Pays the system's own neighbor construction per cell,
#             which on a system whose `cell_neighbors` is geometric means
#             geometry — but IGEO7's is integer digit arithmetic
#             (`src/IGeo7/gbt_neighbors.jl`), so there it is already cheap.
#   twiddle   integer digit arithmetic on the *tile*, no geometry and no
#             memory: a stepper that reads positions straight off the digit
#             suffix, skipping even `subtree_position`. Per-system; no system
#             wires one at present. IGEO7 did until its `cell_neighbors`
#             became integer too and made the generic path the faster of the
#             two (see `src/IGeo7/IGeo7Kernel.jl`, "Tile stencils").
#   table     any of the above, materialized once into `6 x n` `Int32`.
#             O(1) per lookup afterwards, at 4 bytes per cell per direction.
#
# `neighbor_stepper(tile)` picks the best *computed* one; `neighbor_table`
# wraps whatever that is into the materialized form. The trade is setup and
# memory against per-access cost, and it does not resolve the same way twice:
# a repeated sweep over a tile that fits in RAM wants the table, a single pass
# (or a tile that does not fit, or a sparse walk touching a few cells) wants
# the computed stepper.
# --------------------------------------------------------------------------

"""
    AbstractNeighborStepper

Per-tile object answering [`step_neighbors`](@ref)`(stepper, i)`: the positions
of tile cell `i`'s edge neighbors, ascending, with off-tile neighbors **absent**
rather than zero-padded.

That container convention is [`stencil`](@ref)'s: neighbors outside the stored
coverage are simply not there, so a reduction skips them. It also makes the
rim self-describing — `length(step_neighbors(s, i)) < max_neighbors(system)`
is exactly "cell `i` is on the tile border or a pentagon".

Build one with [`neighbor_stepper`](@ref); materialize it with
[`neighbor_table`](@ref).
"""
abstract type AbstractNeighborStepper end

"""
    step_neighbors(stepper, i) -> SmallVector{N,Int}

Positions of tile cell `i`'s edge neighbors, ascending, off-tile neighbors
absent. See [`AbstractNeighborStepper`](@ref).
"""
function step_neighbors end

"""
    neighbor_stepper(tile) -> AbstractNeighborStepper

The best *computed* neighbor stepper for `tile` — no materialization, O(1)
memory. Systems that can decide adjacency from the id arithmetic alone
specialize this; the fallback is [`cell_neighbors`](@ref) resolved through
[`subtree_position`](@ref).
"""
neighbor_stepper(t::DGGSSubtreeIds) = GenericNeighborStepper(t)

"""
    GenericNeighborStepper(tile)

[`cell_neighbors`](@ref) resolved through [`subtree_position`](@ref) — correct
for every kernel-wired system, and the reference the per-system twiddles are
tested against.
"""
struct GenericNeighborStepper{T<:DGGSSubtreeIds} <: AbstractNeighborStepper
    tile::T
end

@inline function step_neighbors(s::GenericNeighborStepper, i::Int)
    t = s.tile
    id = @inbounds t[i]
    nbrs = cell_neighbors(t.system, t.level, id)
    out = _empty_like(nbrs, Int)
    for nb in nbrs                      # ascending by id => ascending by position
        pos = subtree_position(t, nb)
        pos == 0 || (out = SmallCollections.push(out, pos))
    end
    return out
end

# The neighbor container's capacity is a type parameter of whatever
# `cell_neighbors` returned, so the position container inherits it rather than
# re-deriving `max_neighbors` at runtime — same trick `_neighbor_positions` uses.
@inline _empty_like(::SmallVector{N}, ::Type{T}) where {N,T} = SmallVector{N,T}()

"""
    TableNeighborStepper(tile, source=neighbor_stepper(tile))

A materialized halo table: `source` evaluated once for every cell and stored as
a `max_neighbors x n` `Matrix{Int32}`, zero-padded, so
[`step_neighbors`](@ref) becomes a contiguous read.

Costs 4 bytes per cell per direction — 138 MB for an IGEO7 level-13 tile of
5.76e6 cells, and 6.7x that at level 15. Worth it when the same tile is swept
many times (flow accumulation, iterative fill) and it fits; the computed
stepper is the one that scales.

Built threaded over `Threads.nthreads()`.
"""
struct TableNeighborStepper{N} <: AbstractNeighborStepper
    table::Matrix{Int32}        # N x n, zero-padded, ascending within a column
    counts::Vector{UInt8}       # how many of each column are real neighbors
end

function TableNeighborStepper(t::DGGSSubtreeIds, source::AbstractNeighborStepper=neighbor_stepper(t))
    n = length(t)
    N = Int(max_neighbors(t.system))
    table = zeros(Int32, N, n)
    counts = zeros(UInt8, n)
    Threads.@threads for i in 1:n
        v = step_neighbors(source, i)
        @inbounds counts[i] = length(v)
        @inbounds for k in eachindex(v)
            table[k, i] = Int32(v[k])
        end
    end
    return TableNeighborStepper{N}(table, counts)
end

@inline function step_neighbors(s::TableNeighborStepper{N}, i::Int) where {N}
    out = SmallVector{N,Int}()
    @inbounds for k in 1:s.counts[i]
        out = SmallCollections.push(out, Int(s.table[k, i]))
    end
    return out
end

"""
    neighbor_table(tile, source=neighbor_stepper(tile)) -> TableNeighborStepper

Materialize `source` over the whole tile. See
[`TableNeighborStepper`](@ref) for when that is the right trade.
"""
neighbor_table(t::DGGSSubtreeIds, source::AbstractNeighborStepper=neighbor_stepper(t)) =
    TableNeighborStepper(t, source)

# --------------------------------------------------------------------------
# The sweep
# --------------------------------------------------------------------------

"""
    subtree_stencil(f, data, tile; stepper=neighbor_stepper(tile),
                    border=edge_cells(tile), interior=interior_cells(tile, border))

Apply `f(center_value, neighbor_values::SmallVector)` over every cell of
`tile`, sweeping the interior and the border as separate loops.

The split is not cosmetic. Interior cells have a full neighborhood, so their
loop carries no membership branch and the container is always full; border
cells are the ones whose result is computed from *fewer* values than the
interior sees, which is a fact about the answer and not just about the loop.
Keeping them apart lets a caller weight, mask, or halo-exchange the border
without re-deriving which cells it is.

Returns a `Vector` in tile order. `stepper`, `border` and `interior` are
keyword arguments so a repeated sweep pays for each of them once.
"""
function subtree_stencil(f, data::AbstractVector, t::DGGSSubtreeIds;
        stepper::AbstractNeighborStepper=neighbor_stepper(t),
        border::AbstractVector{Int}=edge_cells(t),
        interior::AbstractVector{Int}=interior_cells(t, border))
    length(data) == length(t) || throw(DimensionMismatch(
        "data has $(length(data)) elements, tile has $(length(t)) cells"))
    return _subtree_sweep(f, data, stepper, border, interior)
end

# Function barrier, exactly as `_stencil_sweep`: the stepper's container width
# enters as a type parameter here, so neither loop touches the heap.
function _subtree_sweep(f, data::AbstractVector, stepper, border, interior)
    T = eltype(data)
    probe = step_neighbors(stepper, 1)
    out = similar(data, Base.promote_op(f, T, typeof(_values_like(probe, T))))
    @inbounds for i in interior
        out[i] = f(data[i], _gather(data, step_neighbors(stepper, i)))
    end
    @inbounds for i in border
        out[i] = f(data[i], _gather(data, step_neighbors(stepper, i)))
    end
    return out
end

@inline _values_like(::SmallVector{N}, ::Type{T}) where {N,T} = SmallVector{N,T}()

@inline function _gather(data::AbstractVector{T}, positions::SmallVector{N,Int}) where {N,T}
    values = SmallVector{N,T}()
    @inbounds for j in positions
        values = SmallCollections.push(values, data[j])
    end
    return values
end
