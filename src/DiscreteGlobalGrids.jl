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
  - [`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup),
    [`MultiOrderLookup`](@ref), and [`Cells`](@ref DiscreteGlobalGrids.CellLookups.Cells)
    connect these collections to DimensionalData.

An `Int` is a local index in `1:ncells(grid)`. An
[`AbstractCellIndex`](@ref) is a typed cell identity that records its level.

Region topology is available through [`neighbors`](@ref), [`ring`](@ref),
[`halo`](@ref), [`border`](@ref), [`interior`](@ref), and
[`adjacency`](@ref). [`member_neighbors`](@ref) provides cross-level adjacency
for mixed-level cell sets.

Internal geometry uses `UnitSphericalPoint` — `GeometryOps`', re-exported here
so the type the contract asks an implementor to return needs no module path;
explicitly named wrappers convert longitude and latitude at API boundaries.

# Layout

  - `src/interface/`: abstract types and generic contracts.
  - `src/fallbacks/`: the overridable generic defaults — identity, location,
    geometry, the subtree walkers, and the level-grid and authalic wrappers.
  - `src/engine/`: the machinery no system overrides — the region containers,
    the cursor and index tree, the query planner, and the walks over them.
  - `src/dimensionaldata.jl`: the cube face of [`CellVector`](@ref) —
    [`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup), [`Cells`](@ref DiscreteGlobalGrids.CellLookups.Cells), [`Covering`](@ref).
  - `src/systems/`: grid-system implementations.
  - `src/core/`: the authalic manifold pair.
  - [`Helpers`](@ref): shared allocation-free primitives.
  - `lib/DiscreteGlobalGridsConformanceTesting/`: the test-only package whose
    `test_grid_interface` / `test_hierarchical_system` suites make these
    contracts executable; a new grid or system is expected to pass them.

[`query`](@ref) uses DE9IM.jl predicate types with spherical semantics defined
here. `treeify`, `ncells`, and `getcell` extend and re-export
`ConservativeRegridding.Trees` bindings.

[`regrid`](@ref), [`regrid!`](@ref) and [`plan_regrid`](@ref) are
`GlobalRegridding`'s own, keywords and all; `src/regridding.jl` adds the target
resolution that lets a grid, a [`CellVector`](@ref), a [`CellLookup`](@ref DiscreteGlobalGrids.CellLookups.CellLookup), a
[`MultiOrderCellSet`](@ref), or a bare system spell a destination — a
mixed-level container expands to its reference level. `cellat` is
that package's binding for the same reason the `Trees` ones are. The rest of
the space contract
`src/regridding.jl` fills in — `nchunks`, `ownedindices`, `chunkat`,
`cellcentroid`, `celltree`, `chunkextents`, `chunkindex`, `candidatechunks!`,
`chunkranges`, `subtree` and `destinationdims` — is extended under
`GlobalRegridding`'s own name rather than imported; `samplesites` is left to
that package's own centroid vector.
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
    BarycentricPoint, Weighted, Extensive, PerChunk, Spilled,
    TreeLocator, AnalyticLocator

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
    HierarchicalGridCursor, TiledRasterCursor, node_cell,
    MultiOrderCoverage, MultiOrderCellSet, level_ranges,
    iscontained, coarsest_contained, cell_polygons,
    CellVector, cellset, covering, covering_indices, covering_index,
    reference_level, predicate_indices,
    CentroidCovered, QueryPredicate,
    grow, expand, compact, member_neighbors,
    SubtreeHaloIterator, SubsetHaloIterator, HaloIndexIterator, RegionSide,
    halo_indices, sizehint,
    AdjacencyTable, halocells, haloindices,
    SubsetIndexedCell, cellid,
    mapneighbors, foreachneighbors, StorageOrder,
    Neighborhood, Disc, Ring,
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
    nside, checked_id, chart_perimeter, sampled_cap, corner_cap,
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
include("partitioning.jl")
include("partitioning_backends.jl")
include("cap_cached_tree.jl")

# Point methods load after the shared regridding space they extend.
Base.include(CopernicusDEM, joinpath(@__DIR__, "systems", "CopernicusDEM", "point.jl"))

# Sizing dispatch depends on both raster and regrid-space target types.
include("sizing.jl")
include("deprecated.jl")

# Package-owned raster verbs: keyword destinations cannot dispatch Rasters verbs.
include("raster/common.jl")
include("raster/selection.jl")
include("raster/rasterize.jl")
include("raster/extract.jl")
include("raster/zonal.jl")

# CopernicusDEM is deliberately absent: registering a system enrols it in every
# cross-system sweep, whose hardcoded cases and level choices assume a globally
# uniform cell size. Reach for it by name: `DGG.CopernicusDEMSystem(90)`.
# (Kept ABOVE the docstring — a comment between docstring and function detaches it.)
"""
    systems() -> Tuple{Vararg{AbstractHierarchicalGridSystem}}

Return the six registered global systems in a stable-order tuple:
`IGeo7System()`, `H3System()`, `HEALPixSystem()`, `A5System()`, `S2System()`,
and `ISEA4RSystem()`.

This registry is not an exhaustive list of shipped or externally defined
systems, and it does not control dispatch. [`CopernicusDEMSystem`](@ref) is
also exported. Construct it directly with `CopernicusDEMSystem(30)` or `(90)`;
its latitude-dependent raster lattice does not fit the uniform-size assumptions
of the registry's cross-system sweeps.

Use [`levelfor`](@ref) to compare physical resolutions and [`levelgrid`](@ref)
to construct a grid. Level numbers are not comparable between systems.
See [Choosing a grid](@ref) for the user comparison and
[System capabilities and traversal costs](@ref) for traits and algorithm costs.

```jldoctest
julia> systems()
(IGeo7System(), H3System(), HEALPixSystem(), A5System(), S2System(), ISEA4RSystem())
```
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
export ncells, cellindex, cell_boundary, cell_corners, cell_centroid, cell_polygon
export localindex, globalindex
# `cellposition` stays exported for the deprecation shim in `deprecated.jl`.
export cellposition, rawid, reindex, cellindextypes
export cell_area, cell_extent, getcell
export cellat, neighbors, ring, neighborcount
export treeify, query
export system, level

# --- Cell size and level choice --------------------------------------------
export cellsize, levelfor

# --- Hierarchical system interface -----------------------------------------
export cellindextype, levels, maxlevel, levelgrid, rootcells, children
export node_extent, maxneighbors, maxring, winding, has_sorted_subtrees
export has_congruent_refinement, has_direct_location
export ancestor, descendants, descendant_range
export subtree
export cellid
export mapneighbors, foreachneighbors
export Disc, Ring
public Neighborhood

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
# The cell a cursor node stands for; `nothing` at the synthetic root.
public node_cell
public TiledRasterCursor
public StorageOrder
public AbstractNeed, Cell, Index, Local, Global, Value, Centroid
public cellfield
public cap_inflation
public NeighborCallbackError

# The raster verbs. Rasters.jl owns the same ten names, so callers spell them
# qualified and neither package shadows the other on `using`.
public rasterize
public rasterize!
public extract
public zonal
public mask
public mask!
public boolmask
public boolmask!
public missingmask
public missingmask!

# --- Query predicates (DE9IM.jl types, our semantics) ----------------------
export DE9IMPredicate
export Intersects, Disjoint, Contains, Within, Covers, CoveredBy
export Touches, Crosses, Overlaps, Equals
export CentroidCovered

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
export TreeLocator, AnalyticLocator

# --- Chunk partitioning ----------------------------------------------------
export AbstractPartitioningAlgorithm, WeightedContiguous, MetisPartition, KaHyParPartition, ScotchPartition
export PartitionProblem, ChunkPartition, PartitionBackendUnavailable
export partitionproblem, partitionlabels, partition
export npartitions, partindices, partchunks, partsources, partweights

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
public XDGGS_GRIDS, require_xdggs_readable, xdggs_ellipsoid_attrs

function __init__()
    Base.Experimental.register_error_hint(PartitionBackendUnavailable) do io, err
        err.backend in (:Metis, :KaHyPar_jll, :Scotch) || return
        backend = string(err.backend)
        algorithm = nameof(typeof(err.algorithm))
        print(io, "\nLoad `$backend` with `using $backend` to enable $algorithm. " *
            "If it is not installed, run `import Pkg; Pkg.add(\"$backend\")`.")
    end
end

end # module DiscreteGlobalGrids
