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
    [`CellLookup`](@ref), [`Cells`](@ref), [`Covering`](@ref).
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
`GlobalRegridding`'s, extended in `src/regridding.jl` so that a grid, a
[`CellVector`](@ref), a [`CellLookup`](@ref), a [`MultiOrderCellSet`](@ref), or a
bare system spells a destination. `cellat` is that package's binding for the
same reason the `Trees` ones are. The rest of the space contract
`src/regridding.jl` fills in — `nchunks`, `ownedindices`, `chunkat`,
`cellcentroid`, `samplesites`, `celltree`, `chunkextents`, `chunkindex`,
`candidatechunks!`, `chunkranges`, `subtree` and `destinationdims` — is
extended under `GlobalRegridding`'s own name rather than imported.
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
using GlobalRegridding: Conservative, NearestCell, BilinearPoint,
    Weighted, Extensive, PerChunk, Spilled

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
    HierarchicalGridCursor, MultiOrderCoverage, MultiOrderCellSet, level_ranges,
    iscontained, coarsest_contained, cell_polygons,
    CellVector, cellset, covering, covering_indices,
    grow, expand, compact, member_neighbors,
    SubtreeHaloIterator, SubsetHaloIterator, HaloIndexIterator, RegionSide,
    halo_indices, sizehint,
    AdjacencyTable, halocells, haloindices,
    SubsetIndexedCell, cellid,
    mapneighbors, foreachneighbors, StorageOrder,
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
include("cap_cached_tree.jl")

# After it: a target resolution may be spelled as a raster or a regrid space.
include("sizing.jl")
include("deprecated.jl")

# CopernicusDEM is deliberately absent: registering a system enrols it in every
# cross-system sweep, whose hardcoded cases and level choices assume a globally
# uniform cell size. Reach for it by name: `DGG.CopernicusDEMSystem(90)`.
# (Kept ABOVE the docstring — a comment between docstring and function detaches it.)
"""
    systems() -> Tuple{Vararg{AbstractHierarchicalGridSystem}}

The grid systems shipped by this package, as a stable-order tuple of singletons.

    julia> using DiscreteGlobalGrids

    julia> systems()
    (IGeo7System(), H3System(), HEALPixSystem(), A5System(), S2System(), ISEA4RSystem())

This registry does not include externally defined systems and is not used for
interface dispatch. Tuple order has no semantic meaning.

# The six, and how they differ

The size column is [`cellsize`](@ref) in km at levels 0, 4 and 8 — the side of a
square with the median cell's area, which is the only size a system whose cells
differ in shape and area has. Ask [`levelfor`](@ref) for the level a given
resolution wants rather than reading a level off this table.

| system | cells at level `l` | ≈ size (km) at `l` = 0 / 4 / 8 | cell shape | equal-area |
|---|---|---|---|---|
| [`IGeo7System`](@ref) | `10·7^l + 2` | 6520 / 146 / 3.0 | hexagons + 12 pentagons | by construction; see `IGeo7.equal_area_steradians` |
| [`H3System`](@ref) | `120·7^l + 2` | 2026 / 42.6 / 0.87 | hexagons + 12 pentagons | no (libh3's gnomonic faces) |
| [`HEALPixSystem`](@ref) | `12·4^l` | 6520 / 408 / 25.5 | curvilinear diamonds | yes, exactly `4π/(12·4^l)` |
| [`A5System`](@ref) | `12`, `60`, then `60·4^(l-1)` | 6524 / 365 / 22.8 | pentagons (Cairo-style) | yes |
| [`S2System`](@ref) | `6·4^l` | 9220 / 584 / 36.1 | geodesic quadrilaterals | no; ~2.08× within-level spread |
| [`ISEA4RSystem`](@ref) | `10·4^l` | 7142 / 446 / 27.9 | rhombi on ten diamonds | yes, exactly `4π/(10·4^l)` |

Important cross-system traits:

  - **Neighbour degree varies, and [`Vertex`](@ref)/[`Edge`](@ref) need not
    coincide.** [`maxneighbors`](@ref) is an upper bound. IGeo7/H3 have degree
    6 except for 12 pentagons of degree 5, with identical connectivities. A5's
    4-valent corners give
    `maxneighbors(A5System(), Vertex()) == 11` against
    `maxneighbors(A5System(), Edge()) == 5`; at resolution 1, some cells have
    11 vertex-neighbours and 3 edge-neighbours. Above level 1, ISEA4R has ten
    degree-9 cells at icosahedral vertices 0 and 11, thirty degree-7 cells, and
    degree 8 elsewhere; level 0 is 6-regular.
  - **[`node_extent`](@ref).** S2, ISEA4R, and HEALPix provide exact uninflated
    subtree caps. IGeo7, H3, and A5 use inflated caps; A5 sets
    [`cap_inflation`](@ref) to `1.75`.
  - **[`has_sorted_subtrees`](@ref).** True except for A5, whose canonical order
    has not established the two-sided [`descendant_range`](@ref) contract.
  - **[`border`](@ref) on a subtree region.** IGeo7, H3, HEALPix, ISEA4R, and
    S2 provide `O(border)` walkers; A5 uses the `O(subtree)` fallback. Each is a
    resumable [`EdgeCellIterator`](@ref) / [`InnerCellIterator`](@ref) in
    `O(depth)` memory, which is what `border` and [`interior`](@ref) read.
  - **[`halo`](@ref) on a subtree region.** The same boundary from outside, as a resumable
    [`SubtreeHaloIterator`](@ref) in `O(depth)` memory. `l == level(c)` is the
    cell's own one-ring on every system, emitted ascending without a sort.
    Deeper, five of the six take a specialization and A5 does not:

      + **HEALPix, S2 and ISEA4R** walk the width-one band around their square
        block, one pruned quadtree descent per face the halo touches, taken in
        face order — which is canonical order, so the merge is concatenation.
        `O(halo + depth)` time. A block nowhere flush with its face edge has a
        closed-form count (`4·side + 4` under [`Vertex`](@ref), `4·side` under
        [`Edge`](@ref)) and therefore a `length`; a block that crosses a seam
        has a conservative rectangle band that every candidate is filtered
        through the native one-ring before yielding, and declares
        `SizeUnknown()`.
      + **IGeo7 and H3** approach the halo from the root's own same-level
        neighbours. One level down the calibration is already the answer; deeper,
        each neighbour's subtree-border automaton is seeded with an exposed-direction
        arc calibrated by observation and walked to the target. Neighbouring
        subtrees occupy disjoint [`descendant_range`](@ref)s, so walking them in
        range order is canonical without a heap or a seen-set, and every
        candidate goes through the native one-ring. Neither hexagonal engine
        declares a `length`: the observed counts (`3^(d+1) + 3` around a
        hexagon, `5(3^d + 1)/2` around a pentagon) are validated by enumeration,
        not derived from the transition recurrence.
      + **A5** has no [`descendant_range`](@ref) to prune a descent by and no
        validated boundary automaton, so it scans the target level in `O(1)`
        memory and `O(ncells)` time. Its aperture and Hilbert-like indexing are
        not evidence of a square fast path, and none is inferred from them.

    The generic outside-first hierarchy walk — canonical-order candidates with
    [`node_extent`](@ref) cap pruning — is still what every specialization's
    guards return to and what a newly registered system inherits.
  - **[`halo`](@ref) on any other region.** The same question about an
    arbitrary subset, and always an iterator. A rooted grid holding a complete
    subtree delegates to
    [`SubtreeHaloIterator`](@ref) and keeps its system's specialization;
    everything else — a hole, a forgotten root, an arbitrary id list — takes an
    outside-first walk against membership, pruned by the subset's own index
    spans rather than by geometry: a block the subset holds entire is retired by
    one lookup, and a block no NEIGHBOUR of which it touches is retired by the
    coarse-containment law. The walk follows the subset's boundary, so its cost
    is the halo's and not the subset's. A cell punched
    out of the middle of a subset is outside it and touches it, so it joins the
    halo. A5 is again the exception to the delegation: without
    [`has_sorted_subtrees`](@ref) there is no way to recognise a held subtree, so
    even its rooted complete grid takes the subset walk — to the same answer.
  - **[`adjacency`](@ref).** Every one-ring of a region at once, CSR and
    counter-clockwise, with the out-of-region members dropped (`halo = 0`),
    addressed into a `[region; halo]` buffer (`halo = 1`), or marked with `0`
    in place (`halo = :mark`). The two complete-width shapes preserve slot
    indices against the canonical `one_ring`, so a direction code is a property
    of the cell; the clipped shape preserves order only. Rows exist for
    in-region indices alone, so `halo` above 1 throws and points at
    [`grow`](@ref).
  - **Cross-level adjacency ([`member_neighbors`](@ref)).** Boundary sharing in
    the geometric sense on HEALPix, S2 and ISEA4R, whose four children tile
    their parent exactly; the hierarchy's own relation on IGEO7, H3 and A5,
    where they do not and a member's footprint is not its descendants' union.

# Interoperability caveats

  - ISEA4R's face pairing, diamond numbering, and `(x, y)` orientations are
    package-defined and anchored on vertices `(0, 11)`. Compatibility with
    external ISEA4R identifiers, including DGGAL, is not claimed; establish a
    fixture-derived permutation before interchange.
  - S2's native 64-bit `s2_cellid` is not an available [`reindex`](@ref) scheme.
    The scaffold ordinal is the canonical identifier.

Use [`levels`](@ref) and [`levelgrid`](@ref) to construct queryable grids. Each
system module documents its identifier codec and optimized operations.
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
export cell_polygon, cell_area, cell_extent, getcell
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
export ancestor, descendants, descendant_range
export subtree
export cellid
export mapneighbors, foreachneighbors

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
export Conservative, NearestCell, BilinearPoint, Weighted, Extensive
export PerChunk, Spilled

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

end # module DiscreteGlobalGrids
