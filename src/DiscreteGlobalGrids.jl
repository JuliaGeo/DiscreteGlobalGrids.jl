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

Implementors write primitives; consumers get contracts. A system that overrides
a generic for speed and changes an answer is wrong, and
`test_grid_interface` / `test_hierarchical_system` are how that is caught.

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
    `PartialGrid`, `MultiOrderCellSet`, the query engine.
  - `src/conformance.jl` — the property suites that make the contracts
    executable.
  - `src/systems/{IGeo7,H3,HEALPix}/` — one directory per grid system.
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

# The executable form of the contracts above.
include("conformance.jl")

# Grid systems.
include("systems/IGeo7/IGeo7.jl")
include("systems/H3/H3.jl")
include("systems/HEALPix/HEALPix.jl")

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

# --- Query predicates (DE9IM.jl types, our semantics) ----------------------
export DE9IMPredicate
export Intersects, Disjoint, Contains, Within, Covers, CoveredBy
export Touches, Crosses, Overlaps, Equals

# --- Conformance harness ---------------------------------------------------
using .Conformance: test_grid_interface, test_hierarchical_system
export test_grid_interface, test_hierarchical_system

# --- Manifolds -------------------------------------------------------------
export authalic_sphere

end # module DiscreteGlobalGrids
