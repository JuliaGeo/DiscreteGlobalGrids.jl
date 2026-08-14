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
four required grid methods to system-level counterparts.

A bare `Int` is a position in `1:ncells(grid)`. An
[`AbstractCellIndex`](@ref) is a typed cell identity that records its level.

Adjacency is a verb on subsets, not only on complete levels. On a
[`PartialGrid`](@ref), a [`CellVector`](@ref) or a [`CellLookup`](@ref),
[`neighbors`](@ref) and [`ring`](@ref) are the system's own answer clipped to
membership, in ids or in positions, and [`halo_table`](@ref) is a whole
subset's stencil in one call. [`member_neighbors`](@ref) asks the same question
across the levels of a [`MultiOrderCellSet`](@ref).

Internal geometry uses `GeometryOps.UnitSphericalPoint`; explicitly named
wrappers convert longitude and latitude at API boundaries.

# Layout

  - `src/interface/`: abstract types and generic contracts.
  - `src/fallbacks/`: generic implementations and wrapper types.
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
"""
module DiscreteGlobalGrids

import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeoInterface as GI
import Extents
import ConservativeRegridding
import ConservativeRegridding: Trees
import GeometryOps: SpatialTreeInterface as STI
# The neighbour containers of `neighbors` / `ring`: fixed capacity (the
# `max_neighbors` trait), variable length, no allocation. Only the type is
# brought in; the non-mutating verbs stay qualified (`SmallCollections.push`).
import SmallCollections
using SmallCollections: SmallVector

# Imported qualified, then re-exported by name below: GeometryOps has an
# unrelated internal `DE9IM` matrix struct, so the module name must never be
# `using`-ed into this namespace. DE9IM.jl supplies the predicate *types* only
# — every semantic behind them is implemented in this package.
import DE9IM
using DE9IM: DE9IMPredicate,
    Intersects, Disjoint, Contains, Within, Covers, CoveredBy,
    Touches, Crosses, Overlaps, Equals

# `import`, not `using`: these are extended for every `AbstractGrid` in
# `src/fallbacks/`. They stay `Trees`' own bindings rather than wrappers, so
# re-exporting them here cannot make a `using ConservativeRegridding` alongside
# this package ambiguous.
import ConservativeRegridding.Trees: treeify, ncells, getcell

include("Helpers/Helpers.jl")

# GeometryOps manifolds -> the authalic transform and the authalic-sphere
# compute manifold. Lives here rather than in `Helpers` because `Helpers` is a
# deliberately dependency-free leaf module.
include("core/manifolds.jl")

# The interface: types, then the base grid contract, then the hierarchical
# system contract. Declarations and trait defaults only — no algorithms.
include("interface/types.jl")
include("interface/grid.jl")
include("interface/system.jl")

# Generic implementations of everything the interface declares.
include("fallbacks/fallbacks.jl")

# The concrete types the generic layer ships: the one complete-level grid, the
# one subset grid, the ellipsoid wrapper pair, the one cursor, and the
# multi-order coverage pair. Systems define none of these. Bound here rather
# than beside their exports below because the system modules build on the first
# of them.
using .Fallbacks: HierarchicalLevelGrid, PartialGrid, AuthalicGrid, AuthalicSystem,
    HierarchicalGridCursor, MultiOrderCoverage, MultiOrderCellSet, level_ranges,
    cellindices, is_contained, coarsest_contained, cell_polygons,
    CellVector, cellset, covering, covering_positions,
    EdgeCellIterator, InnerCellIterator, member_neighbors

# The lazy subtree walkers' extension point and the parts a system builds one
# from. Not exported — a caller reaches the iterators, a system reaches these.
using .Fallbacks: collect_subtree,
    MortonCurve, quadrant_step, SquareRimEngine, SquareInteriorEngine

# Grid systems, all six ported. Include order never matters: the two ISEA-family
# systems (IGeo7, ISEA4R) share `src/systems/ISEA/`, and whichever is included
# first defines it behind an `isdefined` guard.
include("systems/IGeo7/IGeo7.jl")
include("systems/H3/H3.jl")
include("systems/HEALPix/HEALPix.jl")
include("systems/A5/A5.jl")
include("systems/S2/S2.jl")
include("systems/ISEA4R/ISEA4R.jl")

using .IGeo7: IGeo7System, Z7Cell, RelativeZ7Cell,
    directioncode, trytranslate
using .H3: H3System, H3Cell
using .HEALPix: HEALPixSystem, HEALPixRingIndex
using .A5: A5System, A5Cell
using .S2: S2System
using .ISEA4R: ISEA4RSystem

# The DimensionalData layer. In-package rather than a package extension because
# DimensionalData is a hard dependency, and because a cube axis is the shape
# most consumers meet this package in — not an optional garnish. Included after
# the systems so its cross-system tests can be written against `systems()`; it
# depends on nothing any of them define.
#
# It is a thin face over `Fallbacks.CellVector`, above: the compression itself
# is DimensionalData-free, and everything that is not a cube — regridding,
# chunking, plain arrays — reaches it without going through this file.
include("dimensionaldata.jl")

using .CellLookups: CellLookup, Cells, Covering

"""
    systems() -> Tuple{Vararg{AbstractHierarchicalGridSystem}}

The grid systems shipped by this package, as a stable-order tuple of singletons.

    julia> using DiscreteGlobalGrids

    julia> systems()
    (IGeo7System(), H3System(), HEALPixSystem(), A5System(), S2System(), ISEA4RSystem())

This registry does not include externally defined systems and is not used for
interface dispatch. Tuple order has no semantic meaning.

# The six, and how they differ

| system | cells at level `l` | cell shape | equal-area |
|---|---|---|---|
| [`IGeo7System`](@ref) | `10·7^l + 2` | hexagons + 12 pentagons | by construction; see `IGeo7.equal_area_steradians` |
| [`H3System`](@ref) | `120·7^l + 2` | hexagons + 12 pentagons | no (libh3's gnomonic faces) |
| [`HEALPixSystem`](@ref) | `12·4^l` | curvilinear diamonds | yes, exactly `4π/(12·4^l)` |
| [`A5System`](@ref) | `12`, `60`, then `60·4^(l-1)` | pentagons (Cairo-style) | yes |
| [`S2System`](@ref) | `6·4^l` | geodesic quadrilaterals | no; ~2.08× within-level spread |
| [`ISEA4RSystem`](@ref) | `10·4^l` | rhombi on ten diamonds | yes, exactly `4π/(10·4^l)` |

Important cross-system traits:

  - **Neighbour degree varies, and [`Vertex`](@ref)/[`Edge`](@ref) need not
    coincide.** [`max_neighbors`](@ref) is an upper bound. IGeo7/H3 have degree
    6 except for 12 pentagons of degree 5, with identical connectivities. A5's
    4-valent corners give
    `max_neighbors(A5System(), Vertex()) == 11` against
    `max_neighbors(A5System(), Edge()) == 5`; at resolution 1, some cells have
    11 vertex-neighbours and 3 edge-neighbours. Above level 1, ISEA4R has ten
    degree-9 cells at icosahedral vertices 0 and 11, thirty degree-7 cells, and
    degree 8 elsewhere; level 0 is 6-regular.
  - **[`node_extent`](@ref).** S2, ISEA4R, and HEALPix provide exact uninflated
    subtree caps. IGeo7, H3, and A5 use inflated caps; A5 sets
    [`cap_inflation`](@ref) to `1.75`.
  - **[`has_sorted_subtrees`](@ref).** True except for A5, whose canonical order
    has not established the two-sided [`descendant_range`](@ref) contract.
  - **[`subtree_border`](@ref).** IGeo7, H3, HEALPix, ISEA4R, and S2 provide
    `O(rim)` walkers; A5 uses the `O(subtree)` fallback. Each is a resumable
    [`EdgeCellIterator`](@ref) / [`InnerCellIterator`](@ref) in `O(depth)`
    memory, of which [`subtree_border`](@ref) and [`subtree_interior`](@ref) are
    the `collect` forms.
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
export LevelIndex
export Connectivity, Vertex, Edge

# --- Base grid interface ---------------------------------------------------
export ncells, cellindex, cell_boundary, cell_centroid
export cellposition, rawid, reindex, cellindextypes
export cell_polygon, cell_area, cell_extent, getcell
export cellat, neighbors, ring, halo_table
export treeify, query
export system, level

# --- Hierarchical system interface -----------------------------------------
export cellindextype, levels, max_level, levelgrid, rootcells, children
export node_extent, cap_inflation, max_neighbors, has_sorted_subtrees
export ancestor, descendants, descendant_range
export subtree_border, subtree_interior
export EdgeCellIterator, InnerCellIterator

# --- Query predicates (DE9IM.jl types, our semantics) ----------------------
export DE9IMPredicate
export Intersects, Disjoint, Contains, Within, Covers, CoveredBy
export Touches, Crosses, Overlaps, Equals

# --- Fallback substrate ----------------------------------------------------
# `using`-ed above the system includes, because `HierarchicalLevelGrid` is what
# all six of them return from `levelgrid` and attach their fast paths to.
export HierarchicalLevelGrid, PartialGrid, HierarchicalGridCursor
export AuthalicGrid, AuthalicSystem
export MultiOrderCoverage, MultiOrderCellSet, level_ranges, cellindices
export is_contained, coarsest_contained, cell_polygons, member_neighbors

# --- The compressed cell collection ----------------------------------------
# A coverage read as a lazy id vector at one level, the region selector over
# one, and the accessor for what either was built from. None of the three
# involves DimensionalData: this is the compression, and `CellLookup` below is
# its cube-shaped face. `covering_positions` is the position-space sibling of
# `covering` and is reached as `DiscreteGlobalGrids.covering_positions`.
export CellVector, covering, cellset

# --- The DimensionalData layer ---------------------------------------------
# A lookup and the dimension it goes in, plus the one selector DimensionalData
# does not already have a spelling for. `At` and `Contains` are
# DimensionalData's own and are not re-exported here: this package already
# exports DE9IM's `Contains`, a predicate about geometries rather than a
# selector about positions, and the two must never end up as the same name in a
# caller's namespace.
export CellLookup, Cells, Covering

# --- Grid systems ----------------------------------------------------------
# System modules are not exported because their names collide with registered
# packages. No system exports a grid type: all six return
# `HierarchicalLevelGrid` from `levelgrid`. S2 and ISEA4R use `LevelIndex` over
# their scaffold ordinals.
export systems
export IGeo7System, Z7Cell, RelativeZ7Cell
export directioncode, trytranslate
export H3System, H3Cell
export HEALPixSystem, HEALPixRingIndex
export A5System, A5Cell
export S2System
export ISEA4RSystem

# --- Manifolds -------------------------------------------------------------
export authalic_sphere

end # module DiscreteGlobalGrids
