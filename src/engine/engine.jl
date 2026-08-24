# The machinery no system overrides: the region containers — `PartialGrid`,
# `CellVector`, `MultiOrderCellSet` — the tree cursor and position tree, the
# query planner, and the halo, adjacency and neighbourhood walks over all of
# them. One implementation each, reached through the exported verbs.
#
# `Fallbacks` is the other half: the overridable defaults a system replaces.
# The dependency runs one way, `Engine` on `Fallbacks`, with two methods
# deliberately living here rather than there because their signatures name a
# type defined here — `Fallbacks._check_wrappable(::PartialGrid)` and
# `halo_engine(::AbstractQuadFaceGridSystem, ...)`.

module Engine

import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractQuadFaceGridSystem,
    AbstractCellIndex, LevelIndex, AbstractCellVector, Connectivity, Vertex, Edge,
    Winding, CounterClockwise, Clockwise, CustomOrder, Unordered,
    ncells, cellindex, cell_boundary, cell_centroid,
    cellposition, rawid, reindex, cellindextypes,
    cell_polygon, cell_area, cell_extent, getcell,
    cellat, cellindices, neighbors, ring, one_ring, neighborcount,
    halo, border, interior, adjacency,
    treeify, query,
    system, level,
    cellindextype, levels, maxlevel, levelgrid, rootcells, children,
    node_extent, cap_inflation, maxneighbors, maxring, winding,
    static_capacity, STATIC_RING_CAP, STATIC_RING_BYTES,
    has_sorted_subtrees,
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
    cap_contains, cell_cap, cells_cap, full_sphere_cap,
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

include("partial_grid.jl")
# The halo is the outside face of the boundary `Fallbacks`' subtree iterators
# walk the inside of, and it reuses that file's stack vocabulary.
include("halo.jl")
include("cursor.jl")
include("position_tree.jl")
include("query.jl")
include("multiorder.jl")
include("cell_vector.jl")
# The stencil layer reads every collection above it — the subset grid, the
# compressed vector, the multi-order set — and the lazy border walkers besides.
include("stencil.jl")
# The positioned iterator depends on the stencil and window helpers.
include("neighborhood.jl")
# The region verbs read every container above and the cursor the sweeps use.
include("region.jl")
# The cached table reads the region verbs and the same cursor.
include("adjacency.jl")
# The algebra composes everything above it: the halo walk grows a region, the
# window helpers merge one, and the multi-order set is what compaction answers.
include("region_algebra.jl")

end # module Engine
