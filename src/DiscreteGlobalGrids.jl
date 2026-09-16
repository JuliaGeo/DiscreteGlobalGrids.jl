"""
    DiscreteGlobalGrids

Discrete global grid systems for the Julia geo ecosystem.

[`AbstractGrid`](@ref) represents a finite global or regional cell collection.
Implementations provide [`ncells`](@ref), [`cellindex`](@ref),
[`cell_boundary`](@ref), and [`cell_centroid`](@ref); generic methods provide
geometry, lookup, topology, trees, and queries.

[`AbstractHierarchicalGridSystem`](@ref) adds analytic parent/child structure,
subtree operations, and hierarchy-based query acceleration. Overrides must
preserve the generic methods' semantics.

A system needs no grid type of its own. [`levelgrid`](@ref) defaults to
[`HierarchicalLevelGrid`](@ref), which holds `(system, level)` and forwards the
base grid interface to the five level-grid primitives a system writes instead.

A bare `Int` is a local index in `1:ncells(grid)`. An
[`AbstractCellIndex`](@ref) is a typed cell identity that records its level.

A **region** is a subset of one complete level — [`PartialGrid`](@ref), any
[`AbstractCellVector`](@ref) or [`AbstractCellLookup`](@ref), so a cube axis read
from a store as readily as one built in memory — or a complete level itself, and
[`subtree`](@ref) is the reifier that makes a subtree one. On a region,
[`neighbors`](@ref) and [`ring`](@ref) return complete-level adjacency clipped
to membership, and four verbs answer about the region as a whole:
[`halo`](@ref) walks what is outside it, including the cells punched out of its
middle; [`border`](@ref) and [`interior`](@ref) split what is inside;
[`adjacency`](@ref) tables every one-ring at once, clipped, completed over a
`[region; halo]` buffer, or marked in place. The three walks are lazy, serial
and `O(depth)` in memory; `adjacency` is the cached, threaded product, and the
one that keeps its halo ([`haloindices`](@ref)).
[`member_neighbors`](@ref) asks the adjacency question across the levels of a
[`MultiOrderCellSet`](@ref), which has no `halo` because it has no single level
to answer at.

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
[`MultiOrderCellSet`](@ref), or a bare system spell a destination. `cellat` is
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
# The unit-sphere point type is in every geometry signature the interface asks
# an implementor to write, so it is re-exported rather than left behind `GO`.
using GeometryOps: UnitSphericalPoint
import GeometryOpsCore as GOCore
import GeoInterface as GI
import Extents
import ConservativeRegridding
import ConservativeRegridding: Trees
import GeometryOps: SpatialTreeInterface as STI
# `neighbors` and `ring` use variable-length, fixed-capacity containers bounded
# by `maxneighbors`.
import SmallCollections
using SmallCollections: SmallVector

# Keep the module qualified because GeometryOps also defines a `DE9IM` name.
# This package imports only DE9IM.jl's predicate types and defines their methods.
import DE9IM
using DE9IM: DE9IMPredicate,
    Intersects, Disjoint, Contains, Within, Covers, CoveredBy,
    Touches, Crosses, Overlaps, Equals

# Extend the original `Trees` bindings so re-exporting them does not introduce
# wrapper functions or binding ambiguity.
import ConservativeRegridding.Trees: treeify, ncells, getcell

# `GlobalRegridding` owns the regridding verbs and the space contract, and has no
# dependency on this package; `src/regridding.jl` implements the contract for the
# grids here. `cellat` is its binding for the same reason the `Trees` ones are:
# extending it rather than shadowing it keeps one function per name for a
# session that holds both surfaces. Every other hook `src/regridding.jl` fills
# in is extended qualified, so those names stay out of this namespace.
import GlobalRegridding: cellat, regrid, regrid!, plan_regrid
# The method and policy names ride along re-exported: they appear in the verbs'
# keyword arguments, so a session that can call `regrid` can also spell
# `method = Conservative()` without a second import.
using GlobalRegridding: Conservative, NearestCell, DirectNearest,
    BarycentricPoint, Weighted, Extensive, PerChunk, Spilled

include("Helpers/Helpers.jl")

# GeometryOps adapters stay outside the dependency-free `Helpers` module.
include("core/manifolds.jl")

# Interface declarations and trait defaults.
include("interface/types.jl")
include("interface/grid.jl")
include("interface/system.jl")

# Declared here rather than beside its method because both modules included
# below import it; the method expands a `MultiOrderCellSet` to a level
# (`src/engine/multiorder.jl`).
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
    CellVector, cellset, covering, covering_indices, predicate_indices,
    grow, expand, compact, member_neighbors,
    SubtreeHaloIterator, SubsetHaloIterator, HaloIndexIterator, RegionSide,
    halo_indices, sizehint,
    AdjacencyTable, halocells, haloindices,
    SubsetIndexedCell, cellid,
    mapneighbors, foreachneighbors, StorageOrder,
    Neighborhood, Disc, Ring,
    NeighborCallbackError,
    AbstractNeed, Cell, Index, Local, Global, Value, Centroid,
    cellfield

# Internal extension points for system-specific subtree walkers and shell
# winding.
using .Fallbacks: collect_subtree,
    MortonCurve, quadrant_step, SquareBorderEngine, SquareInteriorEngine,
    adjacency_shells, shell_ring, shell_disc, checked_steps, _ring_frame, _wind!
using .Engine: SquareBandEngine, square_halo_engine, generic_halo_engine,
    check_halo_level, HexChildHaloEngine, HexArcHaloEngine, hex_halo_engine

# The radix-4 quad-face family: the declarations its members write, and the
# shared arithmetic and geometry their own files call.
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

# The seventh system, included but not registered — see the comment above
# `systems()` for why it stays out of the tuple.
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
    Values, NeighborSlices

# The store-IO layer. Encodings and the chunked lookup own layout mechanics;
# conventions are plain-data metadata logic with no Zarr and no arrays.
# Order matters: errors.jl defines DGGSFormatError for the submodules'
# `import ..DGGSFormatError`, and description.jl types its encoding field
# with encodings.jl's CellEncoding.
include("io/errors.jl")
include("io/encodings.jl")
include("io/chunked_lookup.jl")

using .Encodings: CellEncoding, DenseEncoding, RangesEncoding, ImplicitEncoding,
    ENCODING_REGISTRY, encodingname, register_encoding!, cellaxis,
    idrank, idselect, idcount_between, idvalid, idcell, idtype,
    idranges, write_eligible, validate_ranges
using .ChunkedLookups: ChunkManifest, nchunks, chunkof, chunkbounds,
    ChunkedCellVector, axisindex, chunkmanifest, ChunkedCellLookup

include("io/description.jl")
include("io/conventions.jl")
# The two-dimensional ancestor-subzone layout: arithmetic and vocabulary only,
# read after the conventions whose grid reference table it spells names out of.
include("io/subzones.jl")
include("io/api.jl")

# Following the chunk lines of a stored cube: the plan, the runner, and the
# out-of-core neighbourhood sweep built on them. Reads the cube layer above and
# the region verbs below it.
include("chunks.jl")

# Last: the regridding face reads the grids, the compressed collection, and the
# cube axis alike.
include("regridding.jl")
include("partitioning.jl")
include("partitioning_backends.jl")
include("cap_cached_tree.jl")

# Copernicus DEM answers point queries from its own row arithmetic. The methods
# dispatch on the space above, so the file is read into that system's module
# here rather than where the module is defined.
Base.include(CopernicusDEM, joinpath(@__DIR__, "systems", "CopernicusDEM", "point.jl"))

# After it: a target resolution may be spelled as a raster or a regrid space.
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
# The region contract and its cube face: what code meaning "a set of cells at
# one level" dispatches on, whichever backing produced it.
export AbstractCellVector, AbstractCellLookup
export LevelIndex
export Connectivity, Vertex, Edge
export Winding, CounterClockwise, Clockwise, CustomOrder, Unordered
# `GeometryOps.UnitSphericalPoint`, re-exported: every boundary and centroid
# method in the contract is written in it.
export UnitSphericalPoint

# --- Base grid interface ---------------------------------------------------
export ncells, cellindex, cell_boundary, cell_centroid
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
# `Base.parent(sys, c)` belongs to this list and is absent from it deliberately:
# the hierarchy's parent is a method on Base's function, not a name to re-export.
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
# A region is a subset of one complete level, or a complete level itself: the
# outside, the two insides, and the whole adjacency at once.
export halo, border, interior
export adjacency, AdjacencyTable, halocells, haloindices
# The container those four are answered as, and the conversion into it.
export region

# --- Following a stored cube's chunk lines ---------------------------------
export chunkplan, foreachchunk, mapneighbors!
export MapChunkPlan, MapChunk, ChunkCube
export chunkcube, localindices, ownedindices, axisindices, chunkhalo, halowidth
# Deprecated: the old name of `ownedindices`.
export globalindices

# --- Reachable by name, not exported ---------------------------------------
# The lazy walk types: an argument of the verbs above, never a name a caller
# spells to get an answer.
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
# A traversal order, not a traversal.
public StorageOrder
# The fields a neighbourhood sweep can be asked to stream, and the two index
# spaces `Index` names.
public AbstractNeed, Cell, Index, Local, Global, Value, Centroid
# The vector a `Value` reads when the quantity is computed rather than stored.
public cellfield
# A tuning knob for the default `node_extent`, read by no caller that does not
# implement a system.
public cap_inflation
# Caught, not called.
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

# --- Fallback substrate ----------------------------------------------------
# These fallback types are bound before system modules extend them.
export HierarchicalLevelGrid, PartialGrid
export AuthalicGrid, AuthalicSystem
export MultiOrderCoverage, MultiOrderCellSet, level_ranges, cellindices
export iscontained, coarsest_contained, cell_polygons, member_neighbors

# --- The compressed cell collection ----------------------------------------
# `CellVector` is the DimensionalData-independent compressed collection.
export CellVector, covering, covering_indices, cellset

# --- Region algebra --------------------------------------------------------
# Growth, bulk level movement, and compaction over the region types; `union`,
# `vcat`, `intersect` and `issubset` are Base's and carry no name of their own.
export grow, expand, compact

# --- The DimensionalData layer ---------------------------------------------
# Do not re-export DimensionalData's `Contains`; it conflicts with the DE9IM
# geometry predicate exported above.
export CellLookup, Cells, Covering
export Neighbors, Values, NeighborSlices

# --- Grid systems ----------------------------------------------------------
# Export system types rather than modules whose names collide with packages.
# All seven systems return `HierarchicalLevelGrid` from `levelgrid`.
export systems
export IGeo7System, Z7Cell, RelativeZ7Cell
# The Z7 development-frame direction of a child, which only IGeo7 code reads.
public directioncode
export H3System, H3Cell
export HEALPixSystem, HEALPixRingIndex
export A5System, A5Cell
export S2System
export ISEA4RSystem
# Exported and reachable by name, but not in `systems()` — see the comment there.
export CopernicusDEMSystem

# --- Manifolds -------------------------------------------------------------
# The manifold pair behind `AuthalicSystem`; reached through the wrapper, not
# by name.
public authalic_sphere

# --- Regridding ------------------------------------------------------------
# The verbs are `GlobalRegridding`'s, extended for this package's targets;
# `DGGSpace` is the space they resolve to, and the place chunking is tuned.
# Methods, policies, and storage flavors are re-exported so the verbs' keyword
# arguments are spellable without importing `GlobalRegridding`.
export regrid, regrid!, plan_regrid, DGGSpace
export Conservative, NearestCell, DirectNearest, BarycentricPoint
export Weighted, Extensive
export PerChunk, Spilled

# --- Chunk partitioning ----------------------------------------------------
export AbstractPartitioningAlgorithm, WeightedContiguous, MetisPartition, KaHyParPartition, ScotchPartition
export PartitionProblem, ChunkPartition, PartitionBackendUnavailable
export partitionproblem, partitionlabels, partition
export npartitions, partindices, partchunks, partsources, partweights

# --- Store IO --------------------------------------------------------------
# `detect`, `decode`, `encode!` and `gridname` stay qualified: they are
# extension points, and the names are too generic to export.
export dggread, dggwrite, dggwrite!
export Detection, DGGSFormatError
export DGGSConvention, ZarrDGGSConvention, XdggsConvention,
    LegacyHealpixConvention, DKRZConvention
export register_convention!
export CellEncoding, DenseEncoding, RangesEncoding, ImplicitEncoding
export register_encoding!
export register_grid!
export describe_store
export ChunkedCellLookup, nchunks, chunkof, chunkbounds

# --- The ancestor-subzone layout -------------------------------------------
# The layout descriptor and the store handle its incremental writer hands back.
# The arithmetic around them stays qualified: `columnindices` and friends are
# names a production script spells once, not vocabulary for every user.
export SubzoneLayout, subzonestore
public SUBZONE_LAYOUT, SubzoneRun
public subzone_attrs, subzone_capacity, subzone_cellvector, subzone_columns,
    subzone_coordinate, subzone_depth, subzone_layout, subzone_runs,
    issubzonestore
public SUBZONE_ORDER, SUBZONE_PADDING
public columncell, columnindex, columnlength, columnindices, columnrow,
    subzoneindex, gridnamefor

# The store description vocabulary: what `describe_store` hands back and what a
# convention writer reads, never a name a reader or writer has to spell.
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
