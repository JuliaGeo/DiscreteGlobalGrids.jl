# ---------------------------------------------------------------------------
# The region verbs.
#
# A region is a subset of one complete level — `PartialGrid`, `CellVector`,
# `CellLookup` — or a complete level itself. `halo` walks what is outside it,
# `border` and `interior` split what is inside. All three are lazy and serial;
# `adjacency` is the cached, threaded product.
#
# The currency is INDICES. A halo cell has none in the region, so it is named
# by the complete level's global index; a border or interior cell has one, so
# it is named by the region's own local index.
# ---------------------------------------------------------------------------

const Region = Union{AbstractGrid,CellVector}

# ===========================================================================
# The outside
# ===========================================================================

Base.@constprop :aggressive function halo(region::Region;
        connectivity::Connectivity = Vertex(), cells::Bool = false)
    walk = _halo_walk(region, connectivity)
    return cells ? walk : HaloIndexIterator(walk, _halo_grid(walk))
end

# A rooted grid holding a complete subtree keeps its system's specialization;
# everything else takes the outside-first walk against membership.
function _halo_walk(pg::PartialGrid, connectivity::Connectivity)
    _whole_subtree_range(pg) === nothing ||
        return SubtreeHaloIterator(pg.system, pg.root_id, pg.level; connectivity)
    return SubsetHaloIterator(pg, connectivity, subset_halo_engine(pg.system, pg,
        pg.complete, pg.level, connectivity))
end

# A `CellVector` stores windows but no root ancestor, so it always walks; its
# membership is the window search, `O(log #windows)`.
_halo_walk(cv::CellVector, connectivity::Connectivity) =
    SubsetHaloIterator(cv, connectivity, subset_halo_engine(system(cv), cv,
        cv.grid, level(cv), connectivity))

_halo_walk(grid::AbstractGrid, connectivity::Connectivity) =
    _halo_walk(CellVector(grid), connectivity)

# Used by `adjacency`, where the keyword `halo` shadows this verb.
_halo_indices(region::Region, connectivity::Connectivity) =
    collect(halo(region; connectivity))

# ===========================================================================
# The two insides
# ===========================================================================

"""
    RegionSide

The lazy walk [`border`](@ref) and [`interior`](@ref) return. Its engine yields
`(index, cell)` pairs and the wrapper projects them to whichever of the two
the caller asked for, so `cells = true` costs no second pass.
"""
struct RegionSide{E,C}
    engine::E
end

Base.IteratorSize(::Type{<:RegionSide}) = Base.SizeUnknown()
Base.eltype(::Type{<:RegionSide{E,C}}) where {E,C} = C

@inline _side_value(::Type{Int}, pc) = pc[1]
@inline _side_value(::Type{C}, pc) where {C} = pc[2]

Base.iterate(it::RegionSide) = _side_step(it, iterate(it.engine))
Base.iterate(it::RegionSide, state) = _side_step(it, iterate(it.engine, state))

@inline function _side_step(::RegionSide{E,C}, r) where {E,C}
    r === nothing && return nothing
    return (_side_value(C, r[1]), r[2])
end

Base.show(io::IO, it::RegionSide{E,C}) where {E,C} =
    print(io, "RegionSide(", it.engine, "; cells = ", !(C === Int), ")")

# --- the rooted-subtree engine ---------------------------------------------

# The system's own `O(border)` automaton plus the block offset that turns a
# walked cell into an in-region index. The block is contiguous and ascending,
# so the walk's canonical order is in-region index order.
struct SubtreeSideEngine{W,G}
    walk::W
    complete::G
    lo::Int
end

Base.iterate(e::SubtreeSideEngine) = _subtree_side(e, iterate(e.walk))
Base.iterate(e::SubtreeSideEngine, state) = _subtree_side(e, iterate(e.walk, state))

@inline function _subtree_side(e::SubtreeSideEngine, r)
    r === nothing && return nothing
    c, s = r
    return ((globalindex(e.complete, c)::Int - e.lo + 1, c), s)
end

Base.show(io::IO, e::SubtreeSideEngine) = print(io, e.walk)

# --- the general engine -----------------------------------------------------

# One storage-order pass, keeping the cells whose clipped one-ring is short
# (`border`) or full (`interior`). `SHORT` picks which, so each walk is
# monomorphic. Membership rides the window cursor the sweeps use, so the pass
# costs one native one-ring and `O(degree)` cursor hits per cell.
struct ScanSideEngine{CV,K,SHORT}
    cv::CV
    connectivity::K
end

_scan_side(cv::CellVector, conn::Connectivity, ::Val{S}) where {S} =
    ScanSideEngine{typeof(cv),typeof(conn),S}(cv, conn)

struct ScanSideState
    k::Int
    wj::Int
    hj::Int
end

Base.iterate(e::ScanSideEngine) = iterate(e, ScanSideState(1, 1, 1))

function Base.iterate(e::ScanSideEngine{CV,K,SHORT},
        s::ScanSideState) where {CV,K,SHORT}
    cv = e.cv
    w = cv.windows
    n = length(cv)
    k, wj, hj = s.k, s.wj, s.hj
    while k <= n
        wj = _advance(w, wj, k)
        c = cellindex(cv.grid, _leaf_at(w, wj, k))
        short = false
        for nb in neighbors(cv.grid, c, 1; connectivity = e.connectivity)
            q, hj = _cursor_find(w, wj, hj, globalindex(cv.grid, nb)::Int)
            if q == 0
                short = true
                break
            end
        end
        k += 1
        short == SHORT && return ((k - 1, c), ScanSideState(k, wj, hj))
    end
    return nothing
end

Base.show(io::IO, e::ScanSideEngine{CV,K,SHORT}) where {CV,K,SHORT} =
    print(io, SHORT ? "border(" : "interior(", e.cv, ")")

# --- the verbs --------------------------------------------------------------

Base.@constprop :aggressive border(region::Region;
        connectivity::Connectivity = Vertex(), cells::Bool = false) =
    _region_side(region, connectivity, cells, Val(true))

Base.@constprop :aggressive interior(region::Region;
        connectivity::Connectivity = Vertex(), cells::Bool = false) =
    _region_side(region, connectivity, cells, Val(false))

function _region_side(region::Region, connectivity::Connectivity, cells::Bool,
        short::Val)
    engine = _side_engine(region, connectivity, short)
    return cells ? RegionSide{typeof(engine),_celltype(region)}(engine) :
           RegionSide{typeof(engine),Int}(engine)
end

_celltype(cv::CellVector) = eltype(cv)
_celltype(grid::AbstractGrid) = cellindextype(system(grid))

function _side_engine(pg::PartialGrid, connectivity::Connectivity, ::Val{true})
    r = _whole_subtree_range(pg)
    r === nothing && return _scan_side(CellVector(pg), connectivity, Val(true))
    return SubtreeSideEngine(EdgeCellIterator(pg.system, pg.root_id, pg.level;
        connectivity), pg.complete, first(r))
end

function _side_engine(pg::PartialGrid, connectivity::Connectivity, ::Val{false})
    r = _whole_subtree_range(pg)
    r === nothing && return _scan_side(CellVector(pg), connectivity, Val(false))
    return SubtreeSideEngine(InnerCellIterator(pg.system, pg.root_id, pg.level;
        connectivity), pg.complete, first(r))
end

_side_engine(cv::CellVector, connectivity::Connectivity, short::Val) =
    _scan_side(cv, connectivity, short)

_side_engine(grid::AbstractGrid, connectivity::Connectivity, short::Val) =
    _scan_side(CellVector(grid), connectivity, short)

# ===========================================================================
# The cheap estimate
# ===========================================================================

"""
    sizehint(walk) -> Union{Int,Nothing}

An approximate count of what `walk` will yield, or `nothing` when no cheap bound
exists. Suitable only for `sizehint!`.

    h = sizehint(w)
    out = eltype(w)[]
    h === nothing || sizehint!(out, h)
    for x in w; push!(out, x); end

Unlike `length`, an estimate may over- or undershoot: `Base.IteratorSize` has no
inexact slot, and this is the inexact answer. Engines with an exact count return
it; seam bands use `4·side + 8`, hexagonal walks `3^(d+1) + 3`, and the scanning
and outside-first walks return `nothing` because no general perimeter bound is
available.
"""
function sizehint end

sizehint(it::SubtreeHaloIterator) = _halo_sizehint(it.engine)
sizehint(it::SubsetHaloIterator) = _halo_sizehint(it.engine)
sizehint(it::HaloIndexIterator) = sizehint(it.halo)
sizehint(it::RegionSide) = _side_sizehint(it.engine)

_side_sizehint(::ScanSideEngine) = nothing
_side_sizehint(e::SubtreeSideEngine) =
    Base.IteratorSize(typeof(e.walk)) isa Base.HasLength ? length(e.walk) : nothing

# `collect` on the lazy walks reserves against the estimate; the exact-count
# engines are still validated by `collect_subtree`.
Base.collect(it::SubtreeHaloIterator) = collect_subtree(it, sizehint(it))
Base.collect(it::SubsetHaloIterator) = collect_subtree(it, sizehint(it))
Base.collect(it::HaloIndexIterator) = collect_subtree(it, sizehint(it))
Base.collect(it::RegionSide) = collect_subtree(it, sizehint(it))
