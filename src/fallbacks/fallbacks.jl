# ---------------------------------------------------------------------------
# T2 — Fallback substrate.
#
# Every generic implementation of the interface declared in `src/interface/`
# lives under this directory, split into the files `include`d below. The layers,
# bottom up:
#
#   caps.jl          spherical-cap arithmetic; the one extent vocabulary
#   geometry.jl      cell_polygon / cell_area / cell_extent / node_extent
#   identity.jl      cellposition, reindex, ancestor, descendants
#   subtree.jl       subtree_border / subtree_interior
#   level_grid.jl    `HierarchicalLevelGrid`, the one complete-level grid type
#   partial_grid.jl  `PartialGrid`, the one subset-grid type
#   authalic_grid.jl `AuthalicGrid`/`AuthalicSystem`, the ellipsoid wrapper
#   cursor.jl        `HierarchicalGridCursor` — the hierarchy IS the tree
#   position_tree.jl the fallback tree for grids with no system
#   locate.jl        `cellat`, geometric `neighbors` / `ring`
#   query.jl         the DE9IM query engine (tree prune + rim sandwich)
#   multiorder.jl    `MultiOrderCoverage` -> `MultiOrderCellSet`
#   cell_vector.jl   `CellVector` — a coverage read as a lazy id vector
#
# The parent's generics are extended from in here by importing the bindings and
# defining methods on them, so a method added here is the package's method.
#
# ---------------------------------------------------------------------------
# What a system owes this substrate
#
# Everything below is written against the interface contracts, and each of the
# following is one of those contracts that the generic code *depends on* rather
# than merely documents. A system that gets one wrong produces wrong answers
# here rather than an error, so they are collected in one place:
#
#   * `levelgrid(sys, l)` is the COMPLETE grid at `l`, and its positions ascend
#     in canonical id order. The cursor treats any non-`PartialGrid` grid with a
#     system as that grid, and reads a child's position window straight off
#     `descendant_range`. Its default is a `HierarchicalLevelGrid`, which turns
#     the four grid contracts below into the system-level methods
#     `ncells(sys, l)`, `cellindex(sys, l, i)`, `cellposition(sys, c)`,
#     `cell_boundary(sys, c)` and `cell_centroid(sys, c)`.
#   * `descendant_range(sys, c, l)` is in that same position space, and answers
#     the `l == level(c)` case with the cell's own one-element range — window
#     descent asks for it at the level above the leaves.
#   * `children(sys, c)` and `rootcells(sys)` are in ascending canonical order:
#     selection-mode descent binary-searches the child list.
#   * `cell_centroid(grid, c)` is STRICTLY interior to the cell. The query
#     engine's fast accept is "centroid inside the target implies the cell
#     meets the target", which is only true if it is.
#   * `cell_boundary(grid, c)` is implicitly closed and counter-clockwise seen
#     from outside; area, containment and the rim sandwich all read its winding.
#   * `node_extent(sys, c)` obeys the covering law. The default here (the cell's
#     cap inflated by `cap_inflation`) is available and is O(1) at every depth.
#
# Worth overriding for speed, in this order: `cellposition` on the level grid
# (the generic is a linear scan), `neighbors` (the generic is geometric), and
# `cellat` (the generic is a tree descent plus point-in-polygon).
# ---------------------------------------------------------------------------

module Fallbacks

import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, LevelIndex, Connectivity, Vertex, Edge,
    ncells, cellindex, cell_boundary, cell_centroid,
    cellposition, rawid, reindex, cellindextypes,
    cell_polygon, cell_area, cell_extent, getcell,
    cellat, neighbors, ring, treeify, query, system, level,
    cellindextype, levels, max_level, levelgrid, rootcells, children,
    node_extent, cap_inflation, max_neighbors, has_sorted_subtrees,
    ancestor, descendants, descendant_range,
    subtree_border, subtree_interior
import ..DiscreteGlobalGrids: Helpers

import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeoInterface as GI
import Extents
import DE9IM
import ConservativeRegridding: Trees
import GeometryOps: SpatialTreeInterface as STI
# Only for the `getindex` tie-break in `cell_vector.jl`: a neighbour list is one
# of SmallCollections' vectors, and indexing a cell vector by one is otherwise
# ambiguous against that package's own method.
import SmallCollections

# The unit-sphere vocabulary, spelled once. Every extent in this package is a
# `SphericalCap` and every coordinate a `UnitSphericalPoint`.
#
# `SphericalCap` is aliased as the UnionAll, never as `SphericalCap{Float64}`:
# the two-argument `(point, radius)` constructor is a method of the UnionAll,
# and a parametrised alias would silently reach the three-field default
# constructor instead.
const US = GO.UnitSpherical
const USPoint = GO.UnitSphericalPoint{Float64}
const SphericalCap = GO.UnitSpherical.SphericalCap
const Cap = GO.UnitSpherical.SphericalCap{Float64}

include("caps.jl")
include("geometry.jl")
include("identity.jl")
include("subtree.jl")
include("level_grid.jl")
include("partial_grid.jl")
include("authalic_grid.jl")
include("cursor.jl")
include("position_tree.jl")
include("locate.jl")
include("query.jl")
include("multiorder.jl")
include("cell_vector.jl")

end # module Fallbacks
