"""
    MultiOrderGrid(mov::MultiOrderVector)

Create an O(1) grid view over a [`MultiOrderVector`](@ref)'s stored cells.
Geometry delegates to each cell's level grid, and point location uses the
container's covering-cell lookup. `level(g)` returns `nothing` because the
cells span several levels. Topology operations are unavailable because
coarse/fine boundaries form a nonconforming tiling.
"""
struct MultiOrderGrid{ID,S<:AbstractHierarchicalGridSystem} <: AbstractGrid
    cells::MultiOrderVector{ID,S}
end

"""
    cellset(g::MultiOrderGrid)

Return the backing [`MultiOrderVector`](@ref).
"""
cellset(g::MultiOrderGrid) = g.cells

Base.show(io::IO, g::MultiOrderGrid) =
    print(io, "MultiOrderGrid(", typeof(g.cells.system).name.name, ", ",
        length(g.cells), " cells", isempty(g.cells) ? "" :
        ", levels $(minimum(level, g.cells)):$(maximum(level, g.cells))", ")")

Base.show(io::IO, ::MIME"text/plain", g::MultiOrderGrid) = show(io, g)

# --- the base grid interface ------------------------------------------------

ncells(g::MultiOrderGrid) = length(g.cells)
cellindex(g::MultiOrderGrid, i::Int) = g.cells[i]
system(g::MultiOrderGrid) = g.cells.system

level(::MultiOrderGrid) = nothing

_ownlevel(g::MultiOrderGrid, c::AbstractCellIndex) = levelgrid(g.cells.system, level(c))

cell_boundary(g::MultiOrderGrid, c::AbstractCellIndex) = cell_boundary(_ownlevel(g, c), c)
cell_centroid(g::MultiOrderGrid, c::AbstractCellIndex) = cell_centroid(_ownlevel(g, c), c)
cell_polygon(g::MultiOrderGrid, c::AbstractCellIndex) = cell_polygon(_ownlevel(g, c), c)
cell_extent(g::MultiOrderGrid, c::AbstractCellIndex) = cell_extent(_ownlevel(g, c), c)

# Level-grid area methods preserve exact answers for curvilinear cells.
cell_area(g::MultiOrderGrid, c::AbstractCellIndex) = cell_area(_ownlevel(g, c), c)

# Forward analytic caps and their cost trait at the cell's own level.
Fallbacks.cell_cap(g::MultiOrderGrid, c::AbstractCellIndex) =
    Fallbacks.cell_cap(_ownlevel(g, c), c)

# The cost trait is system-wide, so the reference level represents all cells.
Fallbacks.cell_cap_is_cheap(g::MultiOrderGrid) =
    Fallbacks.cell_cap_is_cheap(levelgrid(g.cells.system, g.cells.reference_level))

# --- membership and location ------------------------------------------------

localindex(g::MultiOrderGrid, c::AbstractCellIndex) = localindex(g.cells, c)

Base.in(c::AbstractCellIndex, g::MultiOrderGrid) = localindex(g, c) !== nothing

"""
    localindex(g::MultiOrderGrid, p::GO.UnitSphericalPoint) -> Union{Int,Nothing}
    cellat(g::MultiOrderGrid, p::GO.UnitSphericalPoint)

Return the covering stored cell of a point, or `nothing` for a container hole.
The lookup combines reference-level point location with one binary search.
"""
localindex(g::MultiOrderGrid, p::GO.UnitSphericalPoint) = localindex(g.cells, p)

cellat(g::MultiOrderGrid, p::GO.UnitSphericalPoint) = cellat(g.cells, p)

# --- the spatial tree -------------------------------------------------------

"""
    treeify(::Manifold, g::MultiOrderGrid) -> IndexTreeNode

Build an O(stored cells) [`IndexTree`](@ref) from cell caps. This cap-based tree
supports the grid's mixed levels.
"""
treeify(::GOCore.Manifold, g::MultiOrderGrid) = IndexTreeNode(IndexTree(g), 1)

# --- unavailable single-level operations -----------------------------------

@noinline function _nolevelindex(g::MultiOrderGrid, c::AbstractCellIndex)
    throw(ArgumentError(
        "$c has no global index in a MultiOrderGrid, whose cells span " *
        "$(isempty(g.cells) ? "no levels" :
           "levels $(minimum(level, g.cells)):$(maximum(level, g.cells))"). " *
        "Use `localindex(grid, c)` for the stored index or " *
        "`globalindex(levelgrid(system(grid), level(c)), c)` for the cell's " *
        "level-specific index."))
end

globalindex(g::MultiOrderGrid, c::AbstractCellIndex) = _nolevelindex(g, c)

@noinline function _nonconforming(verb::AbstractString)
    throw(ArgumentError(
        "`$verb` requires a conforming tiling, but MultiOrderGrid boundaries " *
        "contain coarse/fine T-junctions. Use `member_neighbors(set, c)` on a " *
        "MultiOrderCellSet or `$verb(levelgrid(system(grid), level(c)), c)` " *
        "for level-specific topology."))
end

neighbors(::MultiOrderGrid, ::AbstractCellIndex, ::Integer = 1;
    connectivity::Connectivity = Vertex()) = _nonconforming("neighbors")

neighbors(::MultiOrderGrid, ::AbstractCellIndex, ::Val;
    connectivity::Connectivity = Vertex()) = _nonconforming("neighbors")

ring(::MultiOrderGrid, ::AbstractCellIndex, ::Integer;
    connectivity::Connectivity = Vertex()) = _nonconforming("ring")

ring(::MultiOrderGrid, ::AbstractCellIndex, ::Val;
    connectivity::Connectivity = Vertex()) = _nonconforming("ring")

one_ring(::MultiOrderGrid, ::AbstractCellIndex, ::Connectivity = Vertex()) =
    _nonconforming("one_ring")

halo(::MultiOrderGrid; connectivity::Connectivity = Vertex(), cells::Bool = false) =
    _nonconforming("halo")

border(::MultiOrderGrid; connectivity::Connectivity = Vertex(), cells::Bool = false) =
    _nonconforming("border")

interior(::MultiOrderGrid; connectivity::Connectivity = Vertex(), cells::Bool = false) =
    _nonconforming("interior")

adjacency(::MultiOrderGrid; halo::Union{Integer,Symbol} = 0,
    connectivity::Connectivity = Vertex(), threaded = true) =
    _nonconforming("adjacency")

adjacency(::MultiOrderGrid, ::AbstractVector{<:Integer};
    connectivity::Connectivity = Vertex(), threaded = true) =
    _nonconforming("adjacency")
