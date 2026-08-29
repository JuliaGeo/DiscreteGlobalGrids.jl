module Engine

import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractQuadFaceGridSystem,
    AbstractCellIndex, LevelIndex, AbstractCellVector, Connectivity, Vertex, Edge,
    Winding, CounterClockwise, Clockwise, CustomOrder, Unordered,
    ncells, cellindex, cell_boundary, cell_centroid,
    localindex, globalindex, rawid, reindex, cellindextypes,
    cell_polygon, cell_area, cell_extent, getcell,
    cellat, cellindices, neighbors, ring, one_ring, neighborcount,
    halo, border, interior, adjacency,
    treeify, query, subcursor,
    raster_tiles, raster_shape, raster_localindex, raster_cap,
    system, level,
    cellindextype, levels, maxlevel, levelgrid, rootcells, children,
    node_extent, cap_inflation, maxneighbors, maxring, winding,
    static_capacity, STATIC_RING_CAP, STATIC_RING_BYTES,
    has_sorted_subtrees, has_direct_location,
    ancestor, descendants, descendant_range,
    subtree,
    border_engine, interior_engine, halo_engine,
    lattice_decode, lattice_cell, face_orientation,
    hex_child_direction, seeded_border_engine
import ..DiscreteGlobalGrids: Helpers
import ..Fallbacks
import ..Fallbacks: AuthalicGrid, AuthalicSystem, HierarchicalLevelGrid,
    EdgeCellIterator, InnerCellIterator, collect_subtree,
    MortonCurve, quadrant_step, _SQUARE_CAP,
    checked_id, subtree_curve, nside,
    cap_contains, cell_cap, cell_cap_is_cheap, cells_cap, full_sphere_cap,
    points_cap, lonlat, unit_point, query_point,
    point_in_cell, open_ring, closed_ring,
    _canonical, _match_tolerance, _shared_vertices,
    _tangent_basis, _azimuth, _phase

import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeoInterface as GI
import Extents
import DE9IM
import ConservativeRegridding: Trees
import GeometryOps: SpatialTreeInterface as STI
# Resolve `getindex` dispatch for SmallCollections neighbor vectors.
import SmallCollections

# The UnionAll alias preserves the two-argument `SphericalCap` constructor.
const US = GO.UnitSpherical
const USPoint = GO.UnitSphericalPoint{Float64}
const SphericalCap = GO.UnitSpherical.SphericalCap
const Cap = GO.UnitSpherical.SphericalCap{Float64}

include("partial_grid.jl")
# Halo walkers reuse the subtree-boundary stack types.
include("halo.jl")
include("cursor.jl")
include("index_tree.jl")
# Per-task tables share derived extents across cursors.
include("extent_memo.jl")
# The tiled raster tree packs tiles the way the index tree packs cells.
include("tiled_raster.jl")
include("query.jl")
include("multiorder.jl")
include("cell_vector.jl")
include("multiorder_vector.jl")
include("multiorder_grid.jl")
include("aggregate.jl")
include("stencil.jl")
# The indexed iterator depends on the stencil and window helpers.
include("neighborhood.jl")
include("cellfield.jl")
# The field requests resolve against the clip the sweep above already made.
include("needs.jl")
# The region verbs read every container above and the cursor the sweeps use.
include("region.jl")
# The cached table reads the region verbs and the same cursor.
include("adjacency.jl")
include("region_algebra.jl")

end # module Engine
