# ===========================================================================
# GridInterface: the answers this package's grids give in GridInterface's
# spelling where that differs from their own.
#
# The verbs themselves are GridInterface's bindings (see the import block in
# `DiscreteGlobalGrids.jl`), so a method such as `ncells(::HierarchicalLevelGrid)`
# is already a GridInterface method. What is left here is
#
#   - the traits a GridInterface consumer reads (`cellshape`, `isglobal`, …),
#     where the per-system declarations live in `systems/*/system.jl`;
#   - the position (`::Integer`) forms of the geometry verbs, which this
#     package spells with a typed id;
#   - the two verbs whose GridInterface meaning differs from this package's
#     meaning of the same name: `cell_extent` (a `SphericalCap` there, a lon/lat
#     `Extents.Extent` here) and `cellat` (a position there, a typed id here).
#     Those two stay this package's own functions; GridInterface's get their
#     methods here;
#   - `MultiOrderCellSet` as a GridInterface multi-level grid.
# ===========================================================================

import GridInterface: cellshape, edgekind, tessellation, location_is_analytic,
    neighbors_are_analytic, extent_is_expensive, cell_geometry_is_expensive,
    isglobal, gridrelation, cell_geometry, cell_contains,
    vertices_are_analytic, vertex_numbering_is_analytic, nvertices, vertex_key,
    vertex_index, vertex_cells, vertex_point, VertexKey,
    ConvexCells, SimpleCells, GeneralCells, GeodesicEdges, CurvedEdges,
    Tessellating, Overlapping

# --- identity ---------------------------------------------------------------

# Every grid here names its cells with typed ids, never with its positions, so
# an `Integer` is not an id of any declared scheme.
localindex(::AbstractGrid, ::Integer) = nothing

# A non-`Int` position is the `Int` one (a grid that answers no `Int` form at
# all gets the `MethodError` it would have had, not a loop).
cellindex(grid::AbstractGrid, i::Integer) =
    i isa Int ? throw(MethodError(cellindex, (grid, i))) : cellindex(grid, Int(i))

# Scheme conversion is the system's.
function reindex(::Type{T}, grid::AbstractGrid, c::AbstractCellIndex) where {T<:AbstractCellIndex}
    c isa T && return c
    sys = system(grid)
    sys === nothing && throw(ArgumentError(
        "$(typeof(grid)) has no system, so no conversion to $T"))
    return reindex(T, sys, c)
end

# --- traits -----------------------------------------------------------------

isglobal(::HierarchicalLevelGrid) = true
isglobal(grid::AuthalicGrid) = isglobal(grid.grid)

location_is_analytic(grid::AbstractGrid) =
    (sys = system(grid); sys !== nothing && has_direct_location(sys))
neighbors_are_analytic(grid::AbstractGrid) = system(grid) !== nothing
extent_is_expensive(grid::AbstractGrid) = !(Fallbacks.cell_cap_is_cheap(grid) isa Val{true})
# Every boundary here is computed (a projection or a library call) except the
# lon/lat lattice's.
cell_geometry_is_expensive(grid::AbstractGrid) = !(system(grid) isa CopernicusDEMSystem)

# The per-system shape declarations (`systems/*/system.jl`) name the level
# grid; the subset grids and the authalic wrapper carry their parent's.
cellshape(grid::PartialGrid) = cellshape(grid.complete)
cellshape(grid::AuthalicGrid) = cellshape(grid.grid)
edgekind(grid::PartialGrid) = edgekind(grid.complete)
edgekind(grid::AuthalicGrid) = edgekind(grid.grid)
tessellation(::HierarchicalLevelGrid) = Tessellating()
tessellation(grid::PartialGrid) = tessellation(grid.complete)
tessellation(grid::AuthalicGrid) = tessellation(grid.grid)

# --- geometry, by position ---------------------------------------------------

for f in (:cell_boundary, :cell_corners, :cell_centroid, :cell_area)
    @eval $f(grid::AbstractGrid, i::Integer) = $f(grid, cellindex(grid, Int(i)))
end
cell_geometry(grid::AbstractGrid, i::Integer) = cell_polygon(grid, cellindex(grid, Int(i)))
cell_geometry(grid::AbstractGrid, c::AbstractCellIndex) = cell_polygon(grid, c)

# GridInterface's extent is the pruning cap, which this package calls `cell_cap`.
GridInterface.cell_extent(grid::AbstractGrid, i::Integer) =
    Fallbacks.cell_cap(grid, cellindex(grid, Int(i)))
GridInterface.cell_extent(grid::AbstractGrid, c::AbstractCellIndex) = Fallbacks.cell_cap(grid, c)

# --- location, by position ---------------------------------------------------

# GridInterface's `cellat` answers a position: through the system's own
# location where it has one, otherwise through GridInterface's tree walk.
function GridInterface.cellat(grid::AbstractGrid, p)
    m = GOCore.best_manifold(grid)
    q = GridInterface.check_point(m, p)
    location_is_analytic(grid) || return GridInterface.fallback_cellat(m, grid, q)
    c = cellat(grid, q)
    c === nothing && return nothing
    return localindex(grid, c)
end

# --- regions of a complete level ---------------------------------------------

# A complete level has no outside and no border, and its interior is every
# position: GridInterface's answers for a grid with no parent, given here so the
# whole-grid forms are the same values a wrapper of the grid gives (the lazy
# region iterators are for subsets, where there is something to walk).
halo(grid::HierarchicalLevelGrid; connectivity::Connectivity=Vertex(), cells::Bool=false) =
    cells ? cellindextype(system(grid))[] : Int[]
border(grid::HierarchicalLevelGrid; connectivity::Connectivity=Vertex(), cells::Bool=false) =
    cells ? cellindextype(system(grid))[] : Int[]
interior(grid::HierarchicalLevelGrid; connectivity::Connectivity=Vertex(), cells::Bool=false) =
    cells ? (cellindex(grid, i) for i in 1:ncells(grid)) : Base.OneTo(ncells(grid))

# The authalic wrapper is the same complete level under a latitude warp.
halo(grid::AuthalicGrid; kw...) = halo(grid.grid; kw...)
border(grid::AuthalicGrid; kw...) = border(grid.grid; kw...)
interior(grid::AuthalicGrid; kw...) = interior(grid.grid; kw...)

# --- rooted chunks: this package's O(border) subtree engines --------------------

# GridInterface's primitive `halo(parent, positions, k)` / `border(parent, positions)`
# on a level grid: when `positions` is the descendant range of one cell (a
# rooted chunk, which is what `chunks(grid)` hands out), the subtree engine
# answers in O(border) instead of the generic O(|set|) neighbour merge. Any other
# set, radius or connectivity falls through to the generic method. (Not
# CopernicusDEM: its 1° tiles hold about a million cells and the subtree halo
# engine does not finish there; the generic merge does.)
function _rooted_subtree(lg::HierarchicalLevelGrid, r::AbstractUnitRange{<:Integer})
    sys = system(lg)
    (isempty(r) || !has_sorted_subtrees(sys) || sys isa CopernicusDEMSystem) && return nothing
    c = cellindex(lg, Int(first(r)))
    for l in levels(sys)
        l > lg.level && break
        a = ancestor(sys, c, l)
        descendant_range(sys, a, lg.level) == r && return subtree(sys, a, lg.level)
    end
    return nothing
end
const _GENERIC_HALO = Tuple{Any,AbstractVector{<:Integer},Integer}
const _GENERIC_BORDER = Tuple{Any,AbstractVector{<:Integer}}
function halo(lg::HierarchicalLevelGrid, positions::AbstractUnitRange{<:Integer}, k::Integer=1;
        connectivity::Connectivity=Vertex())
    pg = k == 1 && connectivity isa Vertex ? _rooted_subtree(lg, positions) : nothing
    pg === nothing && return invoke(halo, _GENERIC_HALO, lg, positions, k; connectivity)
    return sort!(collect(Int, halo(pg)))
end
function border(lg::HierarchicalLevelGrid, positions::AbstractUnitRange{<:Integer};
        connectivity::Connectivity=Vertex())
    pg = connectivity isa Vertex ? _rooted_subtree(lg, positions) : nothing
    pg === nothing && return invoke(border, _GENERIC_BORDER, lg, positions; connectivity)
    return sort!(collect(Int, border(pg)))
end

# --- vertices through the authalic wrapper -----------------------------------

# A vertex is combinatorial: the wrapper has its parent's keys and numbering,
# and a vertex point is the parent's corner through the same warp
# `cell_corners` applies, so the two agree bit for bit. Only parents that
# declare analytic vertices are forwarded; the rest keep GridInterface's
# geometric table over the warped corners.
const _AnalyticVertexParent = HierarchicalLevelGrid{<:HEALPixSystem}
const _AuthalicAnalyticVertices = AuthalicGrid{<:_AnalyticVertexParent}
vertices_are_analytic(grid::_AuthalicAnalyticVertices) = vertices_are_analytic(grid.grid)
vertex_numbering_is_analytic(grid::_AuthalicAnalyticVertices) = vertex_numbering_is_analytic(grid.grid)
nvertices(grid::_AuthalicAnalyticVertices) = nvertices(grid.grid)
vertex_key(grid::_AuthalicAnalyticVertices, i::Integer, k::Integer) = vertex_key(grid.grid, i, k)
vertex_key(grid::_AuthalicAnalyticVertices, v::Integer) = vertex_key(grid.grid, v)
vertex_index(grid::_AuthalicAnalyticVertices, key::VertexKey) = vertex_index(grid.grid, key)
vertex_cells(grid::_AuthalicAnalyticVertices, key::VertexKey) = vertex_cells(grid.grid, key)
vertex_point(grid::_AuthalicAnalyticVertices, key::VertexKey) =
    Fallbacks.geodetic_point(grid.transform, vertex_point(grid.grid, key))

# --- chunks -------------------------------------------------------------------

# The level whose cells hold about 4096 of the grid's; a partition by that
# level's ancestor is GridInterface's default chunking of a level grid.
function _chunk_level(sys::AbstractHierarchicalGridSystem, l::Int)
    n = Float64(ncells(levelgrid(sys, l)))
    best, bestd = l, Inf
    for lc in levels(sys)
        lc > l && break
        d = abs(log(n / ncells(levelgrid(sys, lc))) - log(4096))
        d < bestd && ((best, bestd) = (lc, d))
    end
    return best
end

function chunks(grid::AbstractGrid)
    sys = system(grid)
    (sys === nothing || !has_sorted_subtrees(sys) || level(grid) === nothing) && return nothing
    return GridInterface.hierarchy_chunks(grid, _chunk_level(sys, Int(level(grid))))
end

# --- queries ------------------------------------------------------------------

# GridInterface's four predicates are this package's DE-9IM ones, answered as
# ascending positions.
_dgg_predicate(p::GridInterface.Intersects) = DE9IM.Intersects(p.region)
_dgg_predicate(p::GridInterface.Within) = DE9IM.Within(p.region)
_dgg_predicate(p::GridInterface.Contains) = DE9IM.Contains(p.region)
_dgg_predicate(p::GridInterface.CentroidIn) = CentroidCovered(p.region)

function query(grid::AbstractGrid, pred::GridInterface.QueryPredicate)
    ncells(grid) == 0 && return Int[]
    return Engine._query_indices(grid, _dgg_predicate(pred), Engine._query_target(pred.region))
end

# --- MultiOrderCellSet as a multi-level grid ---------------------------------
#
# A coverage result is a grid of typed ids drawn from several levels of one
# system: GridInterface's `MultiLevelSet` shape (parent system, no
# `globalindex`). On a sorted system the set already holds the ascending
# reference-level keys that shape stores, so GridInterface's own set is an O(1)
# view sharing the vectors, and location and the trees are its. On A5 (no
# sorted subtrees) the set answers identity itself and the trees are the
# generic ones.

GridInterface.isgrid(::MultiOrderCellSet) = true
ncells(set::MultiOrderCellSet) = length(set.cells)
cellindex(set::MultiOrderCellSet, i::Integer) = set.cells[i]
cellindextypes(::MultiOrderCellSet{S,ID}) where {S,ID} = (ID,)
gridrelation(::MultiOrderCellSet) = GridInterface.ParentSystem()
localindex(::MultiOrderCellSet, ::Integer) = nothing

_sorted_layout(set::MultiOrderCellSet) = has_sorted_subtrees(set.system)
_refgrid(set::MultiOrderCellSet) = levelgrid(set.system, set.reference_level)

"""
    _multilevel_view(set::MultiOrderCellSet) -> GridInterface.MultiLevelSet

GridInterface's multi-level set over the same cells, in the same order, sharing
the vectors; sorted systems only. The set's keys are already the first
reference-level descendant of each member, ascending, which is that layout.
"""
function _multilevel_view(set::MultiOrderCellSet{S,ID}) where {S,ID}
    _sorted_layout(set) || throw(ArgumentError(
        "$(typeof(set.system)) has no sorted subtrees, so no multi-level view"))
    return GridInterface.MultiLevelSet{S,ID}(set.system, set.cells, set.keys, Int[],
        set.reference_level, true)
end

# Membership: on a sorted system one descendant range and a binary search in
# the keys; otherwise a scan.
function localindex(set::MultiOrderCellSet, c::AbstractCellIndex)
    c isa eltype(set.cells) || return nothing
    if _sorted_layout(set)
        level(c) <= set.reference_level || return nothing
        k = first(descendant_range(set.system, c, set.reference_level))
        i = searchsortedfirst(set.keys, k)
        (i <= length(set.cells) && @inbounds(set.cells[i]) == c) || return nothing
        return i
    end
    i = findfirst(==(c), set.cells)
    return i === nothing ? nothing : Int(i)
end

# The member whose subtree holds reference-level position `p`, or nothing.
@inline function _member_at_reference(set::MultiOrderCellSet, p::Int)
    i = searchsortedlast(set.keys, p)
    i >= 1 || return nothing
    r = descendant_range(set.system, @inbounds(set.cells[i]), set.reference_level)
    return p <= last(r) ? i : nothing
end

GOCore.best_manifold(set::MultiOrderCellSet) = GOCore.best_manifold(_refgrid(set))
cellshape(set::MultiOrderCellSet) = cellshape(_refgrid(set))
edgekind(set::MultiOrderCellSet) = edgekind(_refgrid(set))
tessellation(set::MultiOrderCellSet) =
    has_congruent_refinement(set.system) ? Tessellating() : Overlapping()
cell_geometry_is_expensive(set::MultiOrderCellSet) = cell_geometry_is_expensive(_refgrid(set))
extent_is_expensive(set::MultiOrderCellSet) = extent_is_expensive(_refgrid(set))
location_is_analytic(set::MultiOrderCellSet) =
    _sorted_layout(set) && location_is_analytic(_refgrid(set))
neighbors_are_analytic(::MultiOrderCellSet) = false
# Members of different levels meet at hanging vertices: no dual, no alignment.
GridInterface.conformity(::MultiOrderCellSet) = GridInterface.NonConforming()

# Geometry reads from the member's own level grid (the typed-id forms the set
# already answers for `cell_boundary`/`cell_centroid`/`cell_polygon`).
cell_corners(set::MultiOrderCellSet, c::AbstractCellIndex) =
    cell_corners(levelgrid(set.system, level(c)), c)
cell_area(set::MultiOrderCellSet, c::AbstractCellIndex) =
    cell_area(levelgrid(set.system, level(c)), c)
cell_geometry(set::MultiOrderCellSet, c::AbstractCellIndex) = cell_polygon(set, c)
GridInterface.cell_extent(set::MultiOrderCellSet, c::AbstractCellIndex) =
    Fallbacks.cell_cap(levelgrid(set.system, level(c)), c)
for f in (:cell_boundary, :cell_corners, :cell_centroid, :cell_area, :cell_geometry)
    @eval $f(set::MultiOrderCellSet, i::Integer) = $f(set, cellindex(set, Int(i)))
end
GridInterface.cell_extent(set::MultiOrderCellSet, i::Integer) =
    GridInterface.cell_extent(set, cellindex(set, Int(i)))
function cell_contains(set::MultiOrderCellSet, i::Integer, p)
    c = cellindex(set, Int(i))
    g = levelgrid(set.system, level(c))
    return cell_contains(g, localindex(g, c), p)
end

# Location through the reference level, then up to the member (the hierarchy's
# answer, exact under congruent refinement); the tree walk otherwise.
function GridInterface.cellat(set::MultiOrderCellSet, p)
    m = GOCore.best_manifold(set)
    q = GridInterface.check_point(m, p)
    if location_is_analytic(set)
        pos = GridInterface.cellat(_refgrid(set), q)
        pos === nothing && return nothing
        return _member_at_reference(set, pos)
    end
    return GridInterface.fallback_cellat(m, set, q)
end

# Trees and chunks: GridInterface's own over the view, which emits set positions.
treeify(set::MultiOrderCellSet, n::Integer) = _sorted_layout(set) ?
    treeify(_multilevel_view(set), n) : GridInterface._default_treeify(set, Int(n))
treeify(set::MultiOrderCellSet, hint::GridInterface.LeafHint) =
    treeify(set, GridInterface.leaf_capacity(set, hint))
treeify(set::MultiOrderCellSet; leaf_capacity=GridInterface.Balanced()) =
    treeify(set, leaf_capacity)
subtree(set::MultiOrderCellSet, cells::AbstractVector{<:Integer}) = _sorted_layout(set) ?
    subtree(_multilevel_view(set), cells) : GridInterface._default_subtree(set, cells)
GridInterface.hierarchy_chunks(set::MultiOrderCellSet, l::Integer) =
    GridInterface.hierarchy_chunks(_multilevel_view(set), l)
GridInterface.expand_to_reference(set::MultiOrderCellSet) =
    GridInterface.expand_to_reference(_multilevel_view(set))
GridInterface.aggregation_pairs(set::MultiOrderCellSet) =
    GridInterface.aggregation_pairs(_multilevel_view(set))
GridInterface.check_antichain(set::MultiOrderCellSet) = _sorted_layout(set) ?
    GridInterface.check_antichain(_multilevel_view(set)) :
    GridInterface.check_antichain(GridInterface.MultiLevelSet(set.system, set.cells; check=false))
