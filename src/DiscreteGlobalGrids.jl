"""
    DiscreteGlobalGrids

Discrete global grid systems (DGGS) for the Julia geo ecosystem: native cell
indexing/geometry per system, `DimensionalData.Lookups` integration, and
`ConservativeRegridding.Trees` spatial-tree adapters.

Everything common lives in the package namespace: the abstract DGGS interface
(`AbstractDGGS`, the trait functions `aperture`, `root_count`, ..., and
`cell_polygon`), the per-system registry (`all_systems`), the per-system
operations kernel, and the generic grid/tree family built on it (`DGGSGrid`,
`DGGSPartialGrid`, `subtree_grid`, `DGGSCursor`). The shared chart + ordering
face-grid layer lives there too (`FaceGridSystem` and the `FaceGridSpace` /
`FaceGrid` / `FaceChartGrid` / `FaceGridRoot` family, `src/core/face_grid.jl`);
the HEALPix, S2 and ISEA4R submodules each supply one system singleton, their
orderings, and nothing else. The `DimensionalData` lookups the systems define
share one supertype here too (`AbstractDGGSLookup`), and one lazy id vector
(`DGGSGlobeIds`) that makes a globe-complete dimension cost two words rather
than one id per cell.

A grid becomes a spatial tree in one call — `treeify(grid)`, no manifold, no
`Trees` import: `treeify`, `ncells` and `getcell` are re-exported from
`ConservativeRegridding.Trees`, and `intersects_cap` supplies the query
predicate. Individual grid systems are submodules:

- [`Helpers`](@ref) — shared allocation-free primitives (`SmallList`,
  `sorted_index`, `to_uint64_id`, ...).
- [`ISEA`](@ref) — shared icosahedral machinery: spherical icosahedron
  (vertices, faces, `Orientation`) and the Snyder equal-area projection
  charts. Common to the ISEA family (ISEA7H/IGEO7 today; ISEA3H/4H/9H later).
- [`A5`](@ref) — A5 pentagonal DGGS (pure-Julia native layer).
- [`H3`](@ref) — H3 hexagonal DGGS (backed by `H3_jll`).
- [`HEALPix`](@ref) — EOPF-convention nested HEALPix. Note the capitalization:
  the submodule is `HEALPix`, distinct from the registered Healpix.jl package
  it wraps.
- [`IGeo7`](@ref) — clean-room IGEO7 (ISEA7H + Z7). See `docs/IGeo7/` for the
  clean-room provenance record.
- [`ISEA4R`](@ref) — ISEA4R diamond chart grids: the ten rhombus charts over
  `ISEA`'s Snyder machinery, dense face grids under swappable orderings. The
  ten-diamond layout is a package convention with no external oracle — see
  `docs/design/isea4r_diamond_layout.md`.
- [`S2`](@ref) — S2 cube-face chart grids: the six closed-form charts with the
  quadratic ST↔UV transform, dense face grids under swappable orderings.
  Closed forms only — no s2geometry dependency.

Each system module contains its native layer plus, where ported, an
`<X>Lookups` integration module and an `<X>Kernel.jl` wiring file — the latter
is what puts the system on the generic kernel (S2 and ISEA4R wire geometry
only; no `<X>Lookups` module yet). Systems share generic vocabulary (`cell_center`,
`lonlat_to_cell`, ...), so system-level names are never re-exported here.
"""
module DiscreteGlobalGrids

import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeoInterface as GI
import ConservativeRegridding
import ConservativeRegridding: Trees
import GeometryOps: SpatialTreeInterface as STI
# Imported, never `using`-ed: `DimensionalData` exports a large generic
# vocabulary (`dims`, `rebuild`, `format`, ...) that would collide with this
# namespace's own. Core needs it only for the lookup supertype below.
import DimensionalData as DD
# The neighbor containers of `cell_neighbors` / `neighbor_indices` / `stencil`:
# fixed capacity (the `max_neighbors` trait), variable length, no allocation.
# Only the type is brought in; the non-mutating verbs stay qualified
# (`SmallCollections.push`, `SmallCollections.insert`).
import SmallCollections
using SmallCollections: SmallVector

# Re-exported below so tree consumers never need `Trees` in scope. These are
# `Trees`' own bindings, not wrappers: `treeify` stays the single function whose
# `(::Spherical, ::DGGSGrid)` methods live in `core/generic_cursor.jl`, and
# because the binding is shared, a `using ConservativeRegridding` in the same
# namespace (which exports `ncells`/`getcell` too) cannot make either ambiguous.
using ConservativeRegridding.Trees: treeify, ncells, getcell

include("Helpers/Helpers.jl")

# GeometryOps manifolds -> the authalic transform and the authalic-sphere
# compute manifold. Lives here rather than in `Helpers` because `Helpers` is a
# deliberately dependency-free leaf module.
include("core/manifolds.jl")

# Common DGGS interface and metadata registry (formerly the DGGSScaffold
# staging module).
include("core/interface.jl")
# Per-system operations kernel and the generic tree family built on it.
include("core/kernel.jl")
# The lazy globe-complete id vector and the lookup supertype the per-system
# `<X>Lookups` modules subtype. Both sit above the kernel: `DGGSGlobeIds` is
# `ordinal_to_cell` seen as a vector.
include("core/globe_ids.jl")
# The regional counterpart of `DGGSGlobeIds` — one subtree as a lazy id vector —
# and the neighbor steppers a stencil sweep over such a tile resolves through.
include("core/subtree_ids.jl")
include("core/lookups.jl")
# Lookup-level operations built on the kernel: the neighbor halo table,
# `stencil`, and `zonal`, generic over any lookup whose system wires
# `cell_neighbors` (stencil) / `descendant_range` (zonal).
include("core/lookup_ops.jl")
include("core/grid_types.jl")
include("core/generic_cursor.jl")
# Shared chart + ordering face-grid layer (the FaceGridSystem contract);
# per-system instances live in the HEALPix / S2 / ISEA4R submodules.
include("core/face_grid.jl")
include("core/systems/h3.jl")
include("core/systems/s2.jl")
include("core/systems/a5.jl")
include("core/systems/igeo7.jl")
include("core/systems/isea3h.jl")
include("core/systems/isea4h.jl")
include("core/systems/isea4t.jl")
include("core/systems/isea4r_isea9r.jl")
include("core/systems/ivea.jl")
include("core/systems/rtea.jl")
include("core/systems/rhealpix.jl")
include("core/systems/auspix.jl")
include("core/systems/healpix.jl")

export AbstractDGGS, NotPortedError, OrdinalRangeError
export authalic_sphere
export all_systems, system_name, grid_family, base_solid, cell_shape
export is_equal_area, aperture, canonical_index_name, max_level
export root_count, radix, supports_prefix_ranges
export child_ids, leaf_interval, leaf_count
export cell_polygon, cell_extent

# Operations kernel + generic tree family.
#
# Three kernel generics stay unexported because a system submodule still
# exports the same name for its own native, system-specific function — a
# `using DiscreteGlobalGrids` alongside `using ...H3.H3Native` would make the
# name ambiguous. Reach for them qualified (`DiscreteGlobalGrids.num_cells`):
#
#   num_cells      A5Native, H3Native, IGeo7, IGeo7Lookups
#   cell_boundary  A5Native, A5Lookups, H3Native, H3Lookups, IGeo7, IGeo7Lookups
#   cell_center    A5Lookups, H3Native, H3Lookups, IGeo7, IGeo7Lookups
#
# `cell_polygon_unitsphere` was in that group only while `A5Trees` exported it;
# that module is gone and no submodule claims the name, so it is exported here.
export cell_id_type, has_ordinal_ids, has_descendant_ranges, has_exact_subtree_cap
export has_congruent_geometry
export root_ids, cell_children, cell_parent, cell_descendants
export subtree_leaf_count, cell_to_ordinal, ordinal_to_cell, descendant_range
export max_neighbors, cell_neighbors, subtree_border
export cell_polygon_unitsphere, cell_cap, cells_cap, subtree_cap, cell_cap_inflation
export subtree_polygon_unitsphere
export intersects_cap
export DGGSGrid, DGGSPartialGrid, subtree_grid, DGGSCursor, node_level, node_id
export node_indices
export DGGSGlobeIds, AbstractDGGSLookup
export DGGSSubtreeIds, subtree_position
export subtree_border_positions, subtree_interior_positions
export AbstractNeighborStepper, GenericNeighborStepper, TableNeighborStepper
export neighbor_stepper, step_neighbors, neighbor_table, subtree_stencil
export dggs_system, dggs_level
# `zonal` and `stencil` are also exported from `HealpixLookups`, which owned
# them before they went generic — but as *these* bindings, imported back and
# re-exported, so a `using` of both namespaces cannot make either ambiguous
# (the `treeify` pattern above).
export neighbor_indices, stencil, zonal
export grid_manifold
export treeify, ncells, getcell

export H3DGGS, S2DGGS, A5DGGS, IGEO7DGGS, ISEA3HDGGS, ISEA4HDGGS
export ISEA4TDGGS, ISEA4RDGGS, ISEA9RDGGS, IVEADGGS, RTEADGGS
export RHEALPixDGGS, AusPIXDGGS, HEALPixDGGS

# Shared icosahedral / Snyder equal-area projection machinery
include("ISEA/ISEA.jl")

# Grid systems
include("A5/A5.jl")
include("H3/H3.jl")
include("HEALPix/HEALPix.jl")
include("IGeo7/IGeo7.jl")
include("ISEA4R/ISEA4R.jl")
include("S2/S2.jl")

const HexIndex = IGeo7.HexIndex
const IGEO7Index = IGeo7.IGEO7Index
const RelativeIGEO7Index = IGeo7.RelativeIGEO7Index
const cellarea = IGeo7.cellarea
const cellbearing = IGeo7.cellbearing
const celldistance = IGeo7.celldistance
const cell_to_position = IGeo7.cell_to_position
const directioncode = IGeo7.directioncode
const edges = IGeo7.edges
const get_resolution = IGeo7.get_resolution
const neighbors = IGeo7.neighbors
const position_to_cell = IGeo7.position_to_cell
const trytranslate = IGeo7.trytranslate
export HexIndex, IGEO7Index, RelativeIGEO7Index
export cellarea, cellbearing, celldistance, cell_to_position, directioncode, edges, get_resolution
export neighbors
export position_to_cell
export trytranslate

end # module DiscreteGlobalGrids
