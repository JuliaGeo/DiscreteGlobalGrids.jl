"""
    DiscreteGlobalGrids

Discrete global grid systems for the Julia geo ecosystem, built around one
small interface that every grid — hierarchical or not, global or regional —
implements, and against which every algorithm is written exactly once.

# The two tiers

  - [`AbstractGrid`](@ref) is **one finite collection of cells on the sphere**:
    a complete DGGS level, a regional subset of one, or a standalone structured
    grid with no hierarchy at all. Four required methods
    ([`ncells`](@ref), [`cellindex`](@ref), [`cell_boundary`](@ref),
    [`cell_centroid`](@ref)) buy the whole generic surface — geometry,
    [`cellat`](@ref), [`neighbors`](@ref), [`treeify`](@ref), [`query`](@ref).
  - [`AbstractHierarchicalGridSystem`](@ref) adds **analytic parent/child
    structure** and powers the fast paths: tree pruning under the covering law
    of [`node_extent`](@ref), subtree ranges, sublinear queries. Hierarchy is
    always an optimisation, never a semantic.

A system does not normally write a grid type at all. [`levelgrid`](@ref)
defaults to [`HierarchicalLevelGrid`](@ref), which stores `(system, level)` and
forwards the four required grid methods to system-level counterparts, so the
four are answered once per system rather than once per grid type.

Implementors write primitives; consumers get contracts. A system that overrides
a generic for speed and changes an answer is wrong, and the separate
`DiscreteGlobalGridsConformanceTesting` package is how that is caught.

# Position vs identity

A bare `Int` is always a **position** in a grid's canonical dense order
`1:ncells(grid)` — the storage coordinate data arrays and regridding matrices
are laid out against. A typed [`AbstractCellIndex`](@ref) is always an
**identity** — a name relative to a system, meaningful with no grid in hand,
and self-describing about its level. Ids are never bare integers, so the two
never collide.

All internal geometry is on the **unit sphere**
(`GeometryOps.UnitSphericalPoint`); longitude and latitude appear only at the
edges, in wrappers that say so.

# Layout

  - `src/interface/` — the type vocabulary and every generic's contract.
  - `src/fallbacks/` — the generic implementations: the cursor, the trees,
    [`HierarchicalLevelGrid`](@ref), `PartialGrid`,
    `AuthalicGrid`/`AuthalicSystem`, `MultiOrderCellSet`, the query engine.
  - `lib/DiscreteGlobalGridsConformanceTesting/` — the separate test-only
    package whose property suites make the contracts executable.
  - `src/systems/{IGeo7,H3,HEALPix,A5,S2,ISEA4R}/` — one directory per grid
    system; see [`systems()`](@ref) for what distinguishes them. `src/systems/ISEA/`
    is the Snyder/icosahedron basis IGeo7 and ISEA4R share.
  - [`Helpers`](@ref) — shared allocation-free primitives (`SmallList`,
    `sorted_index`, the `AuthalicTransform`).

Predicate types for [`query`](@ref) come from DE9IM.jl and are re-exported
here (`Intersects`, `Covers`, `Touches`, ...); this package implements their
semantics on the sphere. `treeify`/`ncells`/`getcell` are
`ConservativeRegridding.Trees`' own bindings, extended and re-exported, so a
grid is a regridding source with no imports and no wrapper.
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
    cellindices, is_contained, coarsest_contained, cell_polygons

# Grid systems, all six ported. Include order never matters: the two ISEA-family
# systems (IGeo7, ISEA4R) share `src/systems/ISEA/`, and whichever is included
# first defines it behind an `isdefined` guard.
include("systems/IGeo7/IGeo7.jl")
include("systems/H3/H3.jl")
include("systems/HEALPix/HEALPix.jl")
include("systems/A5/A5.jl")
include("systems/S2/S2.jl")
include("systems/ISEA4R/ISEA4R.jl")

using .IGeo7: IGeo7System, Z7Cell
using .H3: H3System, H3Cell
using .HEALPix: HEALPixSystem, HEALPixRingIndex
using .A5: A5System, A5Cell
using .S2: S2System
using .ISEA4R: ISEA4RSystem

"""
    systems() -> Tuple{Vararg{AbstractHierarchicalGridSystem}}

Every grid system this package ships, as a tuple of singletons.

    julia> using DiscreteGlobalGrids

    julia> systems()
    (IGeo7System(), H3System(), HEALPixSystem(), A5System(), S2System(), ISEA4RSystem())

Written for the two things a caller actually does with such a list: run one
piece of code across all of them (a conformance sweep, a benchmark, a
comparison table), and discover what is available without reading the source.
Order is stable but carries no meaning.

This is a **registry**, not an interface generic: nothing in the package
dispatches on it, and a system defined outside this package is a first-class
system that simply is not in this tuple. It replaces the old `all_systems()`,
which returned metadata-only singletons for systems that had no working
implementation; every entry here is fully ported and passes both
`DiscreteGlobalGridsConformanceTesting` suites.

# The six, and how they differ

| system | cells at level `l` | cell shape | equal-area |
|---|---|---|---|
| [`IGeo7System`](@ref) | `10·7^l + 2` | hexagons + 12 pentagons | by construction; see `IGeo7.equal_area_steradians` |
| [`H3System`](@ref) | `120·7^l + 2` | hexagons + 12 pentagons | no (libh3's gnomonic faces) |
| [`HEALPixSystem`](@ref) | `12·4^l` | curvilinear diamonds | yes, exactly `4π/(12·4^l)` |
| [`A5System`](@ref) | `12`, `60`, then `60·4^(l-1)` | pentagons (Cairo-style) | yes |
| [`S2System`](@ref) | `6·4^l` | geodesic quadrilaterals | no; ~2.08× within-level spread |
| [`ISEA4RSystem`](@ref) | `10·4^l` | rhombi on ten diamonds | yes, exactly `4π/(10·4^l)` |

Traits worth knowing before writing generic code across them:

  - **Neighbour degree is not uniform, and [`Vertex`](@ref)/[`Edge`](@ref) do
    not always coincide.** [`max_neighbors`](@ref) is the container bound, not
    the typical degree. IGeo7/H3 are 6 with 12 pentagons at 5, and their two
    connectivities *do* coincide. A5 is the exception the interface docs used
    to deny: its Cairo-style tiling has **4-valent** corners, so
    `max_neighbors(A5System(), Vertex()) == 11` against
    `max_neighbors(A5System(), Edge()) == 5`, and at resolution 1 a cell really
    has 11 vertex-neighbours and 3 edge-neighbours. ISEA4R is 8 in the lattice
    interior but **9** at the two icosahedral vertices 0 and 11, where five
    diamond corners meet and the diagonal offset resolves to two cells — do not
    assume 8.
  - **[`node_extent`](@ref).** S2 and ISEA4R ship the cell's own *exact,
    uninflated* cap: children tile the parent (a geodesic quad, a chart
    rectangle) so the tight cap is already a legal covering. HEALPix likewise.
    IGeo7, H3 and A5 take the generic inflated-cap default; A5 raises
    [`cap_inflation`](@ref) to `1.75`.
  - **[`has_sorted_subtrees`](@ref).** True for every system but A5, whose
    resolution-0 → resolution-1 quintant numbering is a *rotation* of a5's own
    segment walk, leaving the two-sided [`descendant_range`](@ref) contract
    unverified. A5 therefore treeifies to a *selection-mode* cursor — correct,
    just not range-pruned.
  - **[`subtree_border`](@ref) automata.** IGeo7, H3, HEALPix and ISEA4R ship
    `O(rim)` rim walkers. A5 and S2 keep the `O(subtree)` generic fallback: A5
    because its Hilbert children cover the parent's *area* but not its
    footprint, so there is no digit predicate to read a rim off; S2 because
    nobody has written it yet (its quad lattice plus sorted subtrees is exactly
    the shape that would benefit — an open future-work item, not a defect).

# Interoperability caveats

  - **ISEA4R's ten-diamond numbering has no external oracle.** The pairing of
    the twenty icosahedral faces, the numbering of the ten diamonds, and the
    orientation of the `(x, y)` square inside each are *this package's own
    canonical choice*, anchored on the vertex pair `(0, 11)`. Identifier
    compatibility with any external ISEA4R product — DGGAL included — is
    deliberately not claimed and must not be inferred; interop needs a
    permutation fitted against fixtures first.
  - **S2's native 64-bit `s2_cellid` is not shipped** as an alternate
    [`reindex`](@ref) scheme. The codec is one step from the canonical scaffold
    ordinal, but this repository carries no s2geometry fixtures, so shipping it
    would publish an interoperability claim nothing checks. Purely additive
    later; the scaffold ordinal is canonical either way.

See [`levels`](@ref) and [`levelgrid`](@ref) for turning one of these into a
grid you can query, and each system's own module docstring
(`?DiscreteGlobalGrids.A5`) for its id codec and fast paths.
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
export cellat, neighbors, ring
export treeify, query
export system, level

# --- Hierarchical system interface -----------------------------------------
export cellindextype, levels, max_level, levelgrid, rootcells, children
export node_extent, cap_inflation, max_neighbors, has_sorted_subtrees
export ancestor, descendants, descendant_range
export subtree_border, subtree_interior

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
export is_contained, coarsest_contained, cell_polygons

# --- Grid systems ----------------------------------------------------------
# One singleton and one canonical id type per system, plus the registry that
# lists them. No system exports a grid type: all six return the package's
# `HierarchicalLevelGrid` from `levelgrid` and hang their fast paths off
# `HierarchicalLevelGrid{TheSystem}`. The submodules themselves
# (`DiscreteGlobalGrids.H3` and friends) are deliberately NOT exported: `H3`,
# `HEALPix`, `A5` and `S2` are also the names of registered packages, and a bare
# `using DiscreteGlobalGrids` must not shadow them. Anything past this list —
# `IGeo7.equal_area_steradians`, `IGeo7.z7_string`, `H3.H3Native`,
# `A5.A5Native` — is reached through the qualified module.
#
# S2 and ISEA4R contribute no id type: both are canonically indexed by the
# interface's own `LevelIndex` (exported above), over the scaffold ordinal
# `face * 4^level + hilbert` and `diamond * 4^level + morton` respectively.
export systems
export IGeo7System, Z7Cell
export H3System, H3Cell
export HEALPixSystem, HEALPixRingIndex
export A5System, A5Cell
export S2System
export ISEA4RSystem

# --- Manifolds -------------------------------------------------------------
export authalic_sphere

end # module DiscreteGlobalGrids
