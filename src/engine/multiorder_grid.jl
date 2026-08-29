# ---------------------------------------------------------------------------
# `MultiOrderGrid`: the grid of a `MultiOrderVector`'s STORED cells, each at its
# own level. One cell per stored cell, never one per leaf, which is the whole
# point — a container that stores 25 cells over 3,073 reference-level leaves
# presents 25 here.
#
# Geometry is read per cell from `levelgrid(sys, level(c))`, exactly as
# `cell_polygons(mov)` already does. Point location is the container's own:
# locate at the reference level, then take the covering ancestor, both binary
# searches over the interval index. Nothing here queries a coverage.
#
# The tiling is NON-CONFORMING: a coarse cell's edge carries T-junctions where
# finer neighbours subdivide it, and on a hexagonal hierarchy parent and child
# vertices do not coincide at all. So every topology verb refuses rather than
# returning a geometrically-derived ring that would be plausible and wrong.
# ---------------------------------------------------------------------------

"""
    MultiOrderGrid(mov::MultiOrderVector)

The grid of a mixed-level container's **stored** cells, each at its own level.

  - `ncells` is the stored count and `cellindex(g, i)` is `mov[i]`, so a
    regridding source built on this scales with what was stored, not with the
    leaf count its [`reference_level`](@ref) would name.
  - Geometry — boundary, centroid, area, cap — is forwarded to each cell's own
    `levelgrid`, never derived, so a system that overrides area on a level grid
    keeps that answer here.
  - `level(g)` is `nothing`: the cells are at several levels, so there is no one
    level and no [`globalindex`](@ref) space either.
  - `localindex(g, p)` is the container's covering-ancestor lookup, so a point
    inside a coarse cell resolves to that cell. `localindex(g, c)` stays exact
    membership.
  - The tiling does not conform, so `neighbors`, `ring`, `halo`, `border`,
    `interior` and `adjacency` throw rather than guess across a T-junction.

Construction is `O(1)` and holds the container itself.
"""
struct MultiOrderGrid{ID,S<:AbstractHierarchicalGridSystem} <: AbstractGrid
    cells::MultiOrderVector{ID,S}
end

"""
    cellset(g::MultiOrderGrid)

The container the grid is a face of.
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

# Several levels, so there is no one level to report. `DGGSpace` reads this and
# takes its single-chunk fallback, which is what keeps chunking out of here.
level(::MultiOrderGrid) = nothing

# The grid of the level `c` itself sits at — the one authority on its geometry.
_ownlevel(g::MultiOrderGrid, c::AbstractCellIndex) = levelgrid(g.cells.system, level(c))

cell_boundary(g::MultiOrderGrid, c::AbstractCellIndex) = cell_boundary(_ownlevel(g, c), c)
cell_centroid(g::MultiOrderGrid, c::AbstractCellIndex) = cell_centroid(_ownlevel(g, c), c)
cell_polygon(g::MultiOrderGrid, c::AbstractCellIndex) = cell_polygon(_ownlevel(g, c), c)
cell_extent(g::MultiOrderGrid, c::AbstractCellIndex) = cell_extent(_ownlevel(g, c), c)

# Forwarded, not derived, for `PartialGrid`'s reason: HEALPix and ISEA4R have
# curvilinear edges and override area on their level grid with the exact
# `4pi/ncells`, so the polygon ring is not the authority.
cell_area(g::MultiOrderGrid, c::AbstractCellIndex) = cell_area(_ownlevel(g, c), c)

# What `IndexTree` reads to place a cell, and what a system may answer
# analytically. Forward both the cap and the trait, per level.
Fallbacks.cell_cap(g::MultiOrderGrid, c::AbstractCellIndex) =
    Fallbacks.cell_cap(_ownlevel(g, c), c)

# Every level of one system answers the same, so the reference level speaks for
# all of them; an empty container has no cell to ask about at all.
Fallbacks.cell_cap_is_cheap(g::MultiOrderGrid) =
    Fallbacks.cell_cap_is_cheap(levelgrid(g.cells.system, g.cells.reference_level))

# --- membership and location ------------------------------------------------

# Exact membership: a stored ancestor or descendant of `c` is not `c`.
localindex(g::MultiOrderGrid, c::AbstractCellIndex) = localindex(g.cells, c)

Base.in(c::AbstractCellIndex, g::MultiOrderGrid) = localindex(g, c) !== nothing

"""
    localindex(g::MultiOrderGrid, p::GO.UnitSphericalPoint) -> Union{Int,Nothing}
    cellat(g::MultiOrderGrid, p::GO.UnitSphericalPoint)

The **covering** stored cell of a point: one location at the container's
reference level, then one binary search for the stored ancestor. `nothing`
where the container covers the point nowhere — a hole is a real answer, and no
nearest stored cell stands in for it.

This is what makes a nearest-cell regrid off a `MultiOrderGrid` bit-identical
to one off the reference-level expansion: both decide the point at the same
level with the same call, then read the same stored value.
"""
localindex(g::MultiOrderGrid, p::GO.UnitSphericalPoint) = localindex(g.cells, p)

cellat(g::MultiOrderGrid, p::GO.UnitSphericalPoint) = cellat(g.cells, p)

# --- the spatial tree -------------------------------------------------------

"""
    treeify(::Manifold, g::MultiOrderGrid) -> IndexTreeNode

The `O(stored)` [`IndexTree`](@ref) over the stored cells.

Stated explicitly rather than inherited: the generic sends any grid with a
system to [`HierarchicalGridCursor`](@ref), whose index windows assume one leaf
level. `IndexTree` sorts cell caps and has no level in it at all.
"""
treeify(::GOCore.Manifold, g::MultiOrderGrid) = IndexTreeNode(IndexTree(g), 1)

# --- what a mixed-level tiling cannot answer --------------------------------

@noinline function _nolevelindex(g::MultiOrderGrid, c::AbstractCellIndex)
    throw(ArgumentError(
        "$c has no global index in a MultiOrderGrid: its cells sit at " *
        "$(isempty(g.cells) ? "no level at all" :
           "levels $(minimum(level, g.cells)):$(maximum(level, g.cells))"), so " *
        "there is no one complete level for an index to be carved out of. Use " *
        "`localindex(grid, c)` for the grid's own index, or " *
        "`globalindex(levelgrid(system(grid), level(c)), c)` for the index at " *
        "the cell's own level."))
end

globalindex(g::MultiOrderGrid, c::AbstractCellIndex) = _nolevelindex(g, c)

@noinline function _nonconforming(verb::AbstractString)
    throw(ArgumentError(
        "`$verb` has no answer on a MultiOrderGrid: cells at mixed levels do " *
        "not tile conformingly — a coarse cell's edge carries T-junctions " *
        "where finer neighbours subdivide it, and on a hexagonal hierarchy " *
        "parent and child vertices do not coincide at all — so any ring " *
        "derived from shared vertices would be plausible and wrong. Ask the " *
        "container instead: `member_neighbors(set, c)` on a MultiOrderCellSet, " *
        "or `$verb(levelgrid(system(grid), level(c)), c)` at the cell's own " *
        "level."))
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
