"""
    DiscreteGlobalGrids

Discrete global grid systems for the Julia geo ecosystem.

The package centers on five abstractions:

  - [`AbstractGrid`](@ref) represents a finite global or regional collection.
  - [`AbstractHierarchicalGridSystem`](@ref) provides parent/child structure
    and level grids through [`levelgrid`](@ref).
  - [`PartialGrid`](@ref) and [`CellVector`](@ref) represent regions at one
    level.
  - [`MultiOrderCellSet`](@ref) and [`MultiOrderVector`](@ref) represent
    mixed-level regions.
  - [`CellLookup`](@ref), [`MultiOrderLookup`](@ref), and [`Cells`](@ref)
    connect these collections to DimensionalData.

An `Int` is a local index in `1:ncells(grid)`. An
[`AbstractCellIndex`](@ref) is a typed cell identity that records its level.
Geometry uses the exported `UnitSphericalPoint` type internally and converts
longitude/latitude values at API boundaries.

Region topology is available through [`neighbors`](@ref), [`ring`](@ref),
[`halo`](@ref), [`border`](@ref), [`interior`](@ref), and
[`adjacency`](@ref). [`member_neighbors`](@ref) provides cross-level adjacency
for mixed-level cell sets.

[`query`](@ref) uses DE9IM predicates with spherical semantics. [`regrid`](@ref),
[`regrid!`](@ref), and [`plan_regrid`](@ref) implement GlobalRegridding's space
contract for package grids and cell collections.
"""
module DiscreteGlobalGrids

import GeometryOps as GO
# Export the point type used by every geometry signature.
using GeometryOps: UnitSphericalPoint
import GeometryOpsCore as GOCore
import GeoInterface as GI
import Extents
import ConservativeRegridding
import ConservativeRegridding: Trees
import GeometryOps: SpatialTreeInterface as STI
# `maxneighbors` bounds the fixed-capacity topology buffers.
import SmallCollections
using SmallCollections: SmallVector

# Qualification disambiguates DE9IM.jl from GeometryOps.DE9IM.
import DE9IM
using DE9IM: DE9IMPredicate,
    Intersects, Disjoint, Contains, Within, Covers, CoveredBy,
    Touches, Crosses, Overlaps, Equals

# Share the original tree bindings with ConservativeRegridding.
import ConservativeRegridding.Trees: treeify, ncells, getcell

# Extend GlobalRegridding's bindings so both packages share one function per name.
import GlobalRegridding: cellat, regrid, regrid!, plan_regrid
# Re-export regridding policies used directly in keyword arguments.
using GlobalRegridding: Conservative, NearestCell, DirectNearest,
    BarycentricPoint, Weighted, Extensive, PerChunk, Spilled

include("Helpers/Helpers.jl")

# GeometryOps adapters stay outside the dependency-free `Helpers` module.
include("core/manifolds.jl")

# Interface declarations and trait defaults.
include("interface/types.jl")
include("interface/grid.jl")
include("interface/system.jl")

# Both fallback and engine modules extend this shared binding.
function cellindices end

# The overridable generic defaults, then the machinery that reads them.
include("fallbacks/fallbacks.jl")
include("engine/engine.jl")

# Bind fallback types before including systems that extend them.
using .Fallbacks: HierarchicalLevelGrid, AuthalicGrid, AuthalicSystem,
    EdgeCellIterator, InnerCellIterator

using .Engine: PartialGrid,
    HierarchicalGridCursor, TiledRasterCursor,
    MultiOrderCoverage, MultiOrderCellSet, level_ranges,
    iscontained, coarsest_contained, cell_polygons,
    CellVector, cellset, covering, covering_indices, covering_index,
    reference_level,
    grow, expand, compact, member_neighbors,
    SubtreeHaloIterator, SubsetHaloIterator, HaloIndexIterator, RegionSide,
    halo_indices, sizehint,
    AdjacencyTable, halocells, haloindices,
    SubsetIndexedCell, cellid,
    mapneighbors, foreachneighbors, StorageOrder,
    NeighborCallbackError,
    AbstractNeed, Cell, Index, Local, Global, Value, Centroid,
    cellfield,
    MultiOrderVector, MultiOrderGrid, aggregate, coarsen, complement

using .Fallbacks: collect_subtree,
    MortonCurve, quadrant_step, SquareBorderEngine, SquareInteriorEngine,
    adjacency_shells, shell_ring, shell_disc, checked_steps, _ring_frame, _wind!
using .Engine: SquareBandEngine, square_halo_engine, generic_halo_engine,
    check_halo_level, HexChildHaloEngine, HexArcHaloEngine, hex_halo_engine

using .Fallbacks: nbasefaces, systemname, idname,
    subtree_curve, subtree_orientation,
    nside, checked_id, chart_perimeter, sampled_cap,
    morton_encode, morton_decode

# The Snyder/icosahedron basis IGeo7 and ISEA4R share, before either of them.
include("systems/ISEA/ISEA.jl")
include("systems/IGeo7/IGeo7.jl")
include("systems/H3/H3.jl")
include("systems/HEALPix/HEALPix.jl")
include("systems/A5/A5.jl")
include("systems/S2/S2.jl")
include("systems/ISEA4R/ISEA4R.jl")

# CopernicusDEM has nonuniform resolution and stays outside `systems()`.
include("systems/CopernicusDEM/CopernicusDEM.jl")

using .IGeo7: IGeo7System, Z7Cell, RelativeZ7Cell, directioncode
using .H3: H3System, H3Cell
using .HEALPix: HEALPixSystem, HEALPixRingIndex
using .A5: A5System, A5Cell
using .S2: S2System
using .ISEA4R: ISEA4RSystem
using .CopernicusDEM: CopernicusDEMSystem

# DimensionalData wrappers over the dependency-free `Fallbacks.CellVector`.
include("dimensionaldata.jl")

using .CellLookups: AbstractCellLookup, CellLookup, Cells, Covering, Neighbors,
    Values, NeighborSlices, MultiOrderLookup

# Definitions precede the I/O modules that import their types.
include("io/errors.jl")
include("io/encodings.jl")
include("io/chunked_lookup.jl")

using .Encodings: CellEncoding, DenseEncoding, RangesEncoding, ImplicitEncoding,
    CompactedEncoding,
    ENCODING_REGISTRY, encodingname, register_encoding!, cellaxis,
    idrank, idselect, idcount_between, idvalid, idcell, idtype,
    idranges, write_eligible, validate_ranges
using .ChunkedLookups: ChunkManifest, nchunks, chunkof, chunkbounds,
    ChunkedCellVector, axisindex, chunkmanifest, ChunkedCellLookup

include("io/description.jl")
include("io/conventions.jl")
include("io/subzones.jl")
include("io/api.jl")

include("chunks.jl")

include("regridding.jl")
include("cap_cached_tree.jl")

# Point methods load after the shared regridding space they extend.
Base.include(CopernicusDEM, joinpath(@__DIR__, "systems", "CopernicusDEM", "point.jl"))

# Sizing dispatch depends on both raster and regrid-space target types.
include("sizing.jl")
include("deprecated.jl")

# Cross-system sweeps require uniform cell sizes, which CopernicusDEM lacks.
"""
    systems() -> Tuple{Vararg{AbstractHierarchicalGridSystem}}

Return the package's grid systems as a stable-order tuple of singletons.

    julia> using DiscreteGlobalGrids

    julia> systems()
    (IGeo7System(), H3System(), HEALPixSystem(), A5System(), S2System(), ISEA4RSystem())

The registry supports enumeration; interface dispatch remains open to external
systems. [`CopernicusDEMSystem`](@ref) stays separate because its cells have
nonuniform resolution.

# Shipped systems

The size column reports [`cellsize`](@ref) in kilometers at levels 0, 4, and 8.
Use [`levelfor`](@ref) to select a level for a target resolution.

| system | cells at level `l` | ≈ size (km) at `l` = 0 / 4 / 8 | cell shape | equal-area |
|---|---|---|---|---|
| [`IGeo7System`](@ref) | `10·7^l + 2` | 6520 / 146 / 3.0 | hexagons + 12 pentagons | by construction; see `IGeo7.equal_area_steradians` |
| [`H3System`](@ref) | `120·7^l + 2` | 2026 / 42.6 / 0.87 | hexagons + 12 pentagons | no (libh3's gnomonic faces) |
| [`HEALPixSystem`](@ref) | `12·4^l` | 6520 / 408 / 25.5 | curvilinear diamonds | yes, exactly `4π/(12·4^l)` |
| [`A5System`](@ref) | `12`, `60`, then `60·4^(l-1)` | 6524 / 365 / 22.8 | pentagons (Cairo-style) | yes |
| [`S2System`](@ref) | `6·4^l` | 9220 / 584 / 36.1 | geodesic quadrilaterals | no; ~2.08× within-level spread |
| [`ISEA4RSystem`](@ref) | `10·4^l` | 7142 / 446 / 27.9 | rhombi on ten diamonds | yes, exactly `4π/(10·4^l)` |

Key traits differ by system:

| Trait | Contract |
|---|---|
| [`maxneighbors`](@ref) | Bounds the variable neighbor degree for [`Vertex`](@ref) and [`Edge`](@ref) connectivity. |
| [`node_extent`](@ref) | S2, ISEA4R, and HEALPix provide exact subtree caps; IGeo7, H3, and A5 provide conservative inflated caps. |
| [`has_sorted_subtrees`](@ref) | All listed systems except A5 provide contiguous [`descendant_range`](@ref)s. |
| [`has_direct_location`](@ref) | Every listed system locates a point analytically from its coordinates. |
| [`border`](@ref), [`interior`](@ref) | Five systems provide O(border) subtree walkers; A5 uses the O(subtree) fallback. |
| [`halo`](@ref) | Subtree halos use resumable O(depth)-memory iterators. Square and hexagonal hierarchies provide specialized boundary walks; A5 scans the target level. |
| [`adjacency`](@ref) | `halo = 0` clips rows, `halo = 1` addresses a `[region; halo]` buffer, and `halo = :mark` preserves slots with zero sentinels. |
| [`member_neighbors`](@ref) | HEALPix, S2, and ISEA4R use geometric boundary sharing; IGeo7, H3, and A5 use the hierarchy relation. |

# Interoperability caveats

  - ISEA4R interchange requires a fixture-derived permutation because this
    package defines its own face pairing, diamond numbering, and orientation.
  - S2 uses the scaffold ordinal as its canonical identifier; native
    `s2_cellid` values are outside the available [`reindex`](@ref) schemes.

Use [`levels`](@ref) and [`levelgrid`](@ref) to construct queryable grids. Each
system module documents its identifier codec and optimized operations.
"""
systems() = (IGeo7System(), H3System(), HEALPixSystem(),
             A5System(), S2System(), ISEA4RSystem())

# --- Type vocabulary -------------------------------------------------------
export AbstractGrid, AbstractHierarchicalGridSystem, AbstractCellIndex
export AbstractQuadFaceGridSystem
export AbstractCellVector, AbstractCellLookup
export LevelIndex
export Connectivity, Vertex, Edge
export Winding, CounterClockwise, Clockwise, CustomOrder, Unordered
export UnitSphericalPoint

# --- Base grid interface ---------------------------------------------------
export ncells, cellindex, cell_boundary, cell_centroid
export localindex, globalindex
# `cellposition` stays exported for the deprecation shim in `deprecated.jl`.
export cellposition, rawid, reindex, cellindextypes
export cell_polygon, cell_area, cell_extent, getcell
export cellat, neighbors, ring, neighborcount
export treeify, query
export system, level

# --- Cell size and level choice --------------------------------------------
export cellsize, levelfor

# --- Hierarchical system interface -----------------------------------------
export cellindextype, levels, maxlevel, levelgrid, rootcells, children
export node_extent, maxneighbors, maxring, winding, has_sorted_subtrees
export has_direct_location
export ancestor, descendants, descendant_range
export subtree
export cellid
export mapneighbors, foreachneighbors

# --- The region verbs ------------------------------------------------------
export halo, border, interior
export adjacency, AdjacencyTable, halocells, haloindices
# The container those four are answered as, and the conversion into it.
export region

# --- Following a stored cube's chunk lines ---------------------------------
export chunkplan, foreachchunk, mapneighbors!
export MapChunkPlan, MapChunk, ChunkCube
export chunkcube, localindices, ownedindices, axisindices, chunkhalo, halowidth
# `globalindices` remains a compatibility alias for `ownedindices`.
export globalindices

# --- Reachable by name, not exported ---------------------------------------
public EdgeCellIterator
public InnerCellIterator
public SubtreeHaloIterator
public SubsetHaloIterator
public HaloIndexIterator
public RegionSide
# The inexact size estimate `Base.IteratorSize` has no slot for.
public sizehint
# The indices view of an id halo walk; `halo` already answers in indices.
public halo_indices
public SubsetIndexedCell
public HierarchicalGridCursor
public TiledRasterCursor
public StorageOrder
public AbstractNeed, Cell, Index, Local, Global, Value, Centroid
public cellfield
public cap_inflation
public NeighborCallbackError

# --- Query predicates (DE9IM.jl types, our semantics) ----------------------
export DE9IMPredicate
export Intersects, Disjoint, Contains, Within, Covers, CoveredBy
export Touches, Crosses, Overlaps, Equals

# --- Fallback substrate ----------------------------------------------------
export HierarchicalLevelGrid, PartialGrid
export AuthalicGrid, AuthalicSystem
export MultiOrderCoverage, MultiOrderCellSet, level_ranges, cellindices
export iscontained, coarsest_contained, cell_polygons, member_neighbors

# --- The compressed cell collection ----------------------------------------
export CellVector, covering, covering_indices, cellset

# --- Region algebra --------------------------------------------------------
export grow, expand, compact

# --- Multi-order storage ----------------------------------------------------
export MultiOrderVector, MultiOrderGrid, coarsen, covering_index, complement

public aggregate
public reference_level

# --- The DimensionalData layer ---------------------------------------------
export CellLookup, Cells, Covering
export Neighbors, Values, NeighborSlices
export MultiOrderLookup

# --- Grid systems ----------------------------------------------------------
export systems
export IGeo7System, Z7Cell, RelativeZ7Cell
# The Z7 development-frame direction of a child, which only IGeo7 code reads.
public directioncode
export H3System, H3Cell
export HEALPixSystem, HEALPixRingIndex
export A5System, A5Cell
export S2System
export ISEA4RSystem
export CopernicusDEMSystem

# --- Manifolds -------------------------------------------------------------
public authalic_sphere

# --- Regridding ------------------------------------------------------------
export regrid, regrid!, plan_regrid, DGGSpace
export Conservative, NearestCell, DirectNearest, BarycentricPoint
export Weighted, Extensive
export PerChunk, Spilled

# --- Store IO --------------------------------------------------------------
export dggread, dggwrite, dggwrite!
export Detection, DGGSFormatError
export DGGSConvention, ZarrDGGSConvention, XdggsConvention,
    LegacyHealpixConvention, DKRZConvention
export register_convention!
export CellEncoding, DenseEncoding, RangesEncoding, ImplicitEncoding,
    CompactedEncoding
export register_encoding!
export register_grid!
export describe_store
export ChunkedCellLookup, nchunks, chunkof, chunkbounds

# --- The ancestor-subzone layout -------------------------------------------
export SubzoneLayout, subzonestore
public SUBZONE_LAYOUT, SubzoneRun
public subzone_attrs, subzone_capacity, subzone_cellvector, subzone_columns,
    subzone_coordinate, subzone_depth, subzone_layout, subzone_runs,
    issubzonestore
public SUBZONE_ORDER, SUBZONE_PADDING
public columncell, columnindex, columnlength, columnindices, columnrow,
    subzoneindex, gridnamefor

public StoreSnapshot
public StoreDescription
public ArrayEntry
public ChunkManifest
public GridReference
# The three mutable registries, reached through their `register_*!` verbs.
public CONVENTION_REGISTRY
public DEFAULT_WRITE_CONVENTIONS
public ENCODING_REGISTRY
public GRID_REFERENCE

end # module DiscreteGlobalGrids
