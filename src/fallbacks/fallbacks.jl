# The overridable generic defaults of the interfaces declared in
# `src/interface/`: identity, location, geometry, the subtree walkers, and the
# level-grid and authalic wrappers. Each is a method a system may replace with
# a faster one of its own. The machinery no system overrides — the region
# containers and the walks over them — is `Engine`, in `src/engine/`.
#
# Required system contracts:
#
#   * `levelgrid(sys, l)` is complete and ordered by ascending canonical id.
#     Its default is a `HierarchicalLevelGrid`, which forwards the grid
#     contracts to `ncells(sys, l)`, `cellindex(sys, l, i)`,
#     `globalindex(sys, c)`, `cell_boundary(sys, c)`, `cell_centroid(sys, c)`.
#   * `descendant_range(sys, c, l)` is in that same index space, and answers
#     `l == level(c)` with the cell's one-element index range.
#   * `children(sys, c)` and `rootcells(sys)` are in ascending canonical order:
#     selection-mode descent binary-searches the child list.
#   * `cell_centroid(grid, c)` is strictly interior to the cell.
#   * `cell_boundary(grid, c)` is implicitly closed and counter-clockwise seen
#     from outside; area, containment and the border sandwich all read its winding.
#   * `node_extent(sys, c)` obeys the covering law.
#
# Performance-sensitive overrides are `localindex`, `globalindex`, `neighbors`,
# and `cellat`.

module Fallbacks

import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractQuadFaceGridSystem,
    AbstractCellIndex, LevelIndex, Connectivity, Vertex, Edge,
    Winding, CounterClockwise, Clockwise, CustomOrder, Unordered,
    ncells, cellindex, cell_boundary, cell_centroid,
    localindex, globalindex, rawid, reindex, cellindextypes,
    cell_polygon, cell_area, cell_extent, getcell,
    cellat, cellindices, neighbors, ring, one_ring, neighborcount,
    treeify, query,
    system, level,
    cellindextype, levels, maxlevel, levelgrid, rootcells, children,
    node_extent, cap_inflation, maxneighbors, maxring, winding,
    static_capacity, STATIC_RING_CAP, STATIC_RING_BYTES,
    has_sorted_subtrees, has_congruent_refinement, has_direct_location,
    ancestor, descendants, descendant_range,
    border_engine, interior_engine, halo_engine,
    lattice_decode, lattice_cell, face_orientation,
    hex_child_direction, seeded_border_engine
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
include("authalic_grid.jl")
include("subtree_iterators.jl")
# After the border and interior engines: the radix-4 quad-face family wires them
# to one supertype.
include("quad_face.jl")
include("locate.jl")

end # module Fallbacks
