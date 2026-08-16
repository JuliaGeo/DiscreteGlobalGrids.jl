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

On [`PartialGrid`](@ref), [`CellVector`](@ref), and [`CellLookup`](@ref),
[`neighbors`](@ref) and [`ring`](@ref) return complete-level adjacency clipped
to membership. [`halo`](@ref) lazily returns absent cells adjacent to the
subset, including holes. [`stencil_table`](@ref) addresses complete adjacency
rows into a concatenated `[subset; halo]` buffer.
[`member_neighbors`](@ref) asks the adjacency question across the levels
of a [`MultiOrderCellSet`](@ref), which has no `halo` because it has no single
level to answer at.

[`subtree_border`](@ref) and [`subtree_interior`](@ref) are its inside;
[`subtree_halo`](@ref) is its outside — the level-`l` cells that are not
descendants but touch one. All three are `collect` of a resumable
`O(depth)`-memory iterator ([`EdgeCellIterator`](@ref),
[`InnerCellIterator`](@ref), [`SubtreeHaloIterator`](@ref)). Halos are emitted
in ascending target-grid position. [`halo_positions`](@ref) streams those
positions without materializing ids, and [`halo_sizehint`](@ref) provides an
optional allocation hint where exact length is unavailable.

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
# `neighbors` and `ring` use variable-length, fixed-capacity containers bounded
# by `max_neighbors`.
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

include("Helpers/Helpers.jl")

# GeometryOps adapters stay outside the dependency-free `Helpers` module.
include("core/manifolds.jl")

# Interface declarations and trait defaults.
include("interface/types.jl")
include("interface/grid.jl")
include("interface/system.jl")

# Generic interface implementations.
include("fallbacks/fallbacks.jl")

# Bind fallback types before including systems that extend them.
using .Fallbacks: HierarchicalLevelGrid, PartialGrid, AuthalicGrid, AuthalicSystem,
    HierarchicalGridCursor, MultiOrderCoverage, MultiOrderCellSet, level_ranges,
    cellindices, is_contained, coarsest_contained, cell_polygons,
    CellVector, cellset, covering, covering_positions,
    EdgeCellIterator, InnerCellIterator, member_neighbors,
    SubtreeHaloIterator, SubsetHaloIterator, HaloPositionIterator,
    subtree_halo, halo, halo_positions, halo_sizehint,
    StencilTable, stencil_table

# Internal extension points for system-specific subtree walkers.
using .Fallbacks: collect_subtree,
    MortonCurve, quadrant_step, SquareRimEngine, SquareInteriorEngine,
    SquareBandEngine, square_halo_engine, generic_halo_engine, check_halo_level,
    HexChildHaloEngine, HexArcHaloEngine, hex_halo_engine

# IGeo7, ISEA4R and ISEAGrids share guarded definitions in `src/systems/ISEA/`,
# so include order does not matter.
include("systems/IGeo7/IGeo7.jl")
include("systems/H3/H3.jl")
include("systems/HEALPix/HEALPix.jl")
include("systems/A5/A5.jl")
include("systems/S2/S2.jl")
include("systems/ISEA4R/ISEA4R.jl")
include("systems/ISEAGrids/ISEAGrids.jl")
include("systems/RHEALPix/RHEALPix.jl")
include("systems/IVEARTEA/IVEARTEA.jl")

using .IGeo7: IGeo7System, Z7Cell
using .H3: H3System, H3Cell
using .HEALPix: HEALPixSystem, HEALPixRingIndex
using .A5: A5System, A5Cell
using .S2: S2System
using .ISEA4R: ISEA4RSystem
using .ISEAGrids: ISEA3HSystem, Z3Cell, ISEA4HSystem, ISEA4TSystem
using .RHEALPix: RHEALPixSystem, AusPIXSystem, RHEALPixCell
using .IVEARTEA: IVEA4RSystem, IVEA9RSystem, RTEA4RSystem, RTEA9RSystem

# DimensionalData wrappers over the dependency-free `Fallbacks.CellVector`.
include("dimensionaldata.jl")

using .CellLookups: CellLookup, Cells, Covering

"""
    systems() -> Tuple{Vararg{AbstractHierarchicalGridSystem}}

The grid systems shipped by this package, as a stable-order tuple of singletons.

    julia> using DiscreteGlobalGrids

    julia> systems()
    (IGeo7System(), H3System(), HEALPixSystem(), A5System(), S2System(),
     ISEA4RSystem(), ISEA3HSystem(), ISEA4HSystem(), ISEA4TSystem(),
     RHEALPixSystem(), AusPIXSystem(), IVEA4RSystem(), IVEA9RSystem(),
     RTEA4RSystem(), RTEA9RSystem())

This registry does not include externally defined systems and is not used for
interface dispatch. Tuple order has no semantic meaning.

# The registered systems, and how they differ

| system | cells at level `l` | cell shape | equal-area |
|---|---|---|---|
| [`IGeo7System`](@ref) | `10·7^l + 2` | hexagons + 12 pentagons | by construction; see `IGeo7.equal_area_steradians` |
| [`H3System`](@ref) | `120·7^l + 2` | hexagons + 12 pentagons | no (libh3's gnomonic faces) |
| [`HEALPixSystem`](@ref) | `12·4^l` | curvilinear diamonds | yes, exactly `4π/(12·4^l)` |
| [`A5System`](@ref) | `12`, `60`, then `60·4^(l-1)` | pentagons (Cairo-style) | yes |
| [`S2System`](@ref) | `6·4^l` | geodesic quadrilaterals | no; ~2.08× within-level spread |
| [`ISEA4RSystem`](@ref) | `10·4^l` | rhombi on ten diamonds | yes, exactly `4π/(10·4^l)` |
| [`ISEA3HSystem`](@ref) | `10·3^l + 2` | hexagons + 12 pentagons | yes |
| [`ISEA4HSystem`](@ref) | `10·4^l + 2` | hexagons + 12 pentagons | yes |
| [`ISEA4TSystem`](@ref) | `20·4^l` | triangles | yes, exactly `4π/(20·4^l)` |
| [`RHEALPixSystem`](@ref) | `6·9^l` | curvilinear quads/darts/caps | yes, exactly `2π/(3·9^l)` |
| [`AusPIXSystem`](@ref) | `6·9^l` | WGS84 rHEALPix profile | yes on the authalic sphere |
| `IVEA4RSystem` / `RTEA4RSystem` | `10·4^l` | rhombi | yes |
| `IVEA9RSystem` / `RTEA9RSystem` | `10·9^l` | rhombi | yes |

Important cross-system traits:

  - **Neighbour degree varies, and [`Vertex`](@ref)/[`Edge`](@ref) need not
    coincide.** [`max_neighbors`](@ref) is an upper bound. IGeo7/H3 have degree
    6 except for 12 pentagons of degree 5, with identical connectivities. A5's
    4-valent corners give
    `max_neighbors(A5System(), Vertex()) == 11` against
    `max_neighbors(A5System(), Edge()) == 5`; at resolution 1, some cells have
    11 vertex-neighbours and 3 edge-neighbours. Above level 1, ISEA4R has ten
    degree-9 cells at icosahedral vertices 0 and 11, thirty degree-7 cells, and
    degree 8 elsewhere; level 0 is 6-regular. Among the literature families,
    ISEA3H/4H are degree 6 with twelve degree-5 pentagons and the two
    connectivities coincide; ISEA4T is 3-regular under `Edge()` and reaches 12
    under `Vertex()`; rHEALPix/AusPIX and the IVEA/RTEA rhombi are 4-regular
    under `Edge()` and reach 8 and 9 respectively under `Vertex()`.
  - **[`node_extent`](@ref).** Systems with spatially nested chart trees provide
    tight caps. The ISEA3H/4H prefix trees are deliberately non-spatial and may
    require caps wider than a hemisphere; their documentation states that
    query-pruning tradeoff.
  - **[`has_sorted_subtrees`](@ref).** False for A5 and the IVEA/RTEA rhombic
    systems; their canonical orders do not establish the two-sided
    [`descendant_range`](@ref) contract.
  - **[`subtree_border`](@ref).** IGeo7, H3, HEALPix, ISEA4R, and S2 provide
    `O(rim)` walkers; the remaining systems use the `O(subtree)` fallback. Each is a resumable
    [`EdgeCellIterator`](@ref) / [`InnerCellIterator`](@ref) in `O(depth)`
    memory, of which [`subtree_border`](@ref) and [`subtree_interior`](@ref) are
    the `collect` forms.
  - **[`subtree_halo`](@ref).** The same boundary from outside, as a resumable
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
        each neighbour's subtree-rim automaton is seeded with an exposed-direction
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
  - **[`halo`](@ref).** The same question about a SUBSET rather than a subtree,
    on [`PartialGrid`](@ref), [`CellVector`](@ref) and [`CellLookup`](@ref), and
    always an iterator. A rooted grid holding a complete subtree delegates to
    [`SubtreeHaloIterator`](@ref) and keeps its system's specialization;
    everything else — a hole, a forgotten root, an arbitrary id list — takes an
    outside-first walk against membership, pruned by the subset's own position
    spans rather than by geometry: a block the subset holds entire is retired by
    one lookup, and a block no NEIGHBOUR of which it touches is retired by the
    coarse-containment law. The walk follows the subset's boundary, so its cost
    is the halo's and not the subset's. A cell punched
    out of the middle of a subset is outside it and touches it, so it joins the
    halo. A5 is again the exception to the delegation: without
    [`has_sorted_subtrees`](@ref) there is no way to recognise a held subtree, so
    even its rooted complete grid takes the subset walk — to the same answer.
  - **[`stencil_table`](@ref).** The chunk-plus-halo addressing, and the only
    verb here whose rows are guaranteed COMPLETE: given a subset and its halo
    materialised as ascending positions, CSR rows of full one-rings indexed into
    the concatenated `[chunk; halo]` buffer, in rotational order. `k == 1` only,
    since a width-one halo completes nothing wider, and a neighbour in neither
    half is an error rather than a short row. A rooted complete subtree resolves
    membership by integer range and everything else — a hole, a scattered id
    set, all of A5 — by `cellposition`.
  - **Cross-level adjacency ([`member_neighbors`](@ref)).** Boundary sharing in
    the geometric sense on HEALPix, S2, ISEA4R, ISEA4T, rHEALPix/AusPIX, and
    the rhombic IVEA/RTEA systems, whose children tile their parent; the
    hierarchy's own relation on IGEO7, H3, A5, ISEA3H and ISEA4H, where a
    member's footprint is not its descendants' union.

# Interoperability caveats

  - ISEA4R's face pairing, diamond numbering, and `(x, y)` orientations are
    package-defined and anchored on vertices `(0, 11)`. Compatibility with
    external ISEA4R identifiers, including DGGAL, is not claimed; establish a
    fixture-derived permutation before interchange.
  - S2's native 64-bit `s2_cellid` is not an available [`reindex`](@ref) scheme.
    The scaffold ordinal is the canonical identifier.
  - ISEA4T uses a package-defined face/path index; the sealed DGGRID corpus does
    not establish a SEQNUM crosswalk. ISEA3H uses the documented Z3 prefix tree.
    Because ISEA3H/4H prefix parents are non-spatial, fixed-level
    `MultiOrderCoverage` is supported but budget mode throws rather than claim
    a covering that the hierarchy cannot guarantee.
  - IVEA/RTEA currently expose the aperture-4 and aperture-9 rhombic systems.
    Exact 3H/7H atlas parity, pentagon mapping and corrected odd-level RI7
    lookup require deeper post-fix oracle vectors and are not claimed here.
  - ISEA3H/4H `cell_boundary` chooses a common Snyder face for each paired edge
    when one exists, with eight samples and a great-circle fallback across the
    five-face development cut. It is a finite polygon approximation; analytic
    `cell_area` is the equal-area value and is independent of it.
  - IVEA/RTEA `cell_boundary` densifies each chart edge to thirty-two segments
    at every level — 128 points per cell, four times the density of the other
    congruent systems — which holds midpoint-to-chord deviation under 2e-5 rad
    at a level-0 root and far under it below. The count is level-INDEPENDENT on
    purpose: chording a parent and its children at different steps opens a lens
    between the two polylines that a mixed-level `MultiOrderCoverage` reads as a
    sliver. `cell_area` is the analytic equal-area value and does not read the
    polygon.

Use [`levels`](@ref) and [`levelgrid`](@ref) to construct queryable grids. Each
system module documents its identifier codec and optimized operations.
"""
systems() = (IGeo7System(), H3System(), HEALPixSystem(),
             A5System(), S2System(), ISEA4RSystem(),
             ISEA3HSystem(), ISEA4HSystem(), ISEA4TSystem(),
             RHEALPixSystem(), AusPIXSystem(),
             IVEA4RSystem(), IVEA9RSystem(), RTEA4RSystem(), RTEA9RSystem())

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
export node_extent, cap_inflation, coarse_probe_rings, max_neighbors,
    has_sorted_subtrees
export ancestor, descendants, descendant_range
export subtree_border, subtree_interior
export EdgeCellIterator, InnerCellIterator
export SubtreeHaloIterator, SubsetHaloIterator, HaloPositionIterator
export subtree_halo, halo, halo_positions, halo_sizehint
export StencilTable, stencil_table

# --- Query predicates (DE9IM.jl types, our semantics) ----------------------
export DE9IMPredicate
export Intersects, Disjoint, Contains, Within, Covers, CoveredBy
export Touches, Crosses, Overlaps, Equals

# --- Fallback substrate ----------------------------------------------------
# These fallback types are bound before system modules extend them.
export HierarchicalLevelGrid, PartialGrid, HierarchicalGridCursor
export AuthalicGrid, AuthalicSystem
export MultiOrderCoverage, MultiOrderCellSet, level_ranges, cellindices
export is_contained, coarsest_contained, cell_polygons, member_neighbors

# --- The compressed cell collection ----------------------------------------
# `CellVector` is the DimensionalData-independent compressed collection.
export CellVector, covering, cellset

# --- The DimensionalData layer ---------------------------------------------
# Do not re-export DimensionalData's `Contains`; it conflicts with the DE9IM
# geometry predicate exported above.
export CellLookup, Cells, Covering

# --- Grid systems ----------------------------------------------------------
# Export system types rather than modules whose names collide with packages.
export systems
export IGeo7System, Z7Cell
export H3System, H3Cell
export HEALPixSystem, HEALPixRingIndex
export A5System, A5Cell
export S2System
export ISEA4RSystem
export ISEA3HSystem, Z3Cell, ISEA4HSystem, ISEA4TSystem
export RHEALPixSystem, AusPIXSystem, RHEALPixCell
export IVEA4RSystem, IVEA9RSystem, RTEA4RSystem, RTEA9RSystem

# --- Manifolds -------------------------------------------------------------
export authalic_sphere

end # module DiscreteGlobalGrids
