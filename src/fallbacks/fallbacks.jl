# Generic implementations of the interfaces declared in `src/interface/`.
#
# Required system contracts:
#
#   * `levelgrid(sys, l)` is complete and ordered by ascending canonical id.
#     Its default is a `HierarchicalLevelGrid`, which forwards the grid
#     contracts to `ncells(sys, l)`, `cellindex(sys, l, i)`,
#     `cellposition(sys, c)`, `cell_boundary(sys, c)`, `cell_centroid(sys, c)`.
#   * `descendant_range(sys, c, l)` is in that same position space, and answers
#     `l == level(c)` with the cell's one-element position range.
#   * `children(sys, c)` and `rootcells(sys)` are in ascending canonical order:
#     selection-mode descent binary-searches the child list.
#   * `cell_centroid(grid, c)` is strictly interior to the cell.
#   * `cell_boundary(grid, c)` is implicitly closed and counter-clockwise seen
#     from outside; area, containment and the rim sandwich all read its winding.
#   * `node_extent(sys, c)` obeys the covering law.
#
# Performance-sensitive overrides are `cellposition`, `neighbors`, and `cellat`.

module Fallbacks

import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, LevelIndex, Connectivity, Vertex, Edge,
    ncells, cellindex, cell_boundary, cell_centroid,
    cellposition, rawid, reindex, cellindextypes,
    cell_polygon, cell_area, cell_extent, getcell,
    cellat, neighbors, ring, halo_table, neighborcount, treeify, query, system, level,
    cellindextype, levels, max_level, levelgrid, rootcells, children,
    node_extent, cap_inflation, max_neighbors, has_sorted_subtrees,
    ancestor, descendants, descendant_range,
    subtree_border, subtree_interior,
    rim_engine, interior_engine, halo_engine,
    lattice_decode, lattice_cell, face_orientation,
    hex_child_direction, seeded_rim_engine
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

# Unit-sphere aliases. Keep `SphericalCap` as the UnionAll so its two-argument
# constructor remains available.
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
# After `authalic_grid.jl`: the wrapper forwards both engines, so it must be a
# type by the time those methods are defined.
include("subtree_iterators.jl")
# Straight after: the halo is the outside face of the boundary the file above
# walks the inside of, and it reuses that file's stack vocabulary. Its geometry
# provider calls `_match_tolerance` / `_shared_vertices` from `locate.jl`, which
# is included below — a forward reference between function bodies in one module,
# which Julia resolves at call time.
include("halo.jl")
include("cursor.jl")
include("position_tree.jl")
include("locate.jl")
include("query.jl")
include("multiorder.jl")
include("cell_vector.jl")
# Last: the stencil layer reads every collection above it — the subset grid, the
# compressed vector, the multi-order set — and the lazy rim walkers besides.
include("stencil.jl")
# The positioned iterator depends on the stencil and window helpers.
include("neighborhood.jl")

end # module Fallbacks
