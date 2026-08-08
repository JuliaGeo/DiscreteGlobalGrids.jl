"""
    S2DGGS()

S2 — cube-based Hilbert quadtree of spherical quadrilaterals, indexed by
`s2_cellid`, levels 0–30. Not equal-area.

Report section 1.2. Storage model `:sorted_ids_or_coverings`; native tree
strategy `:hilbert_quadtree` (6 roots, radix 4, prefix ranges).

# What is built today

**Geometry, over the scaffold ordinal.** `DiscreteGlobalGrids.S2` owns the six
closed-form cube-face charts and dense whole-sphere face grids at any
`nside >= 1` (`S2.S2FaceGrid(nside; ordering)` → `treeify` → `Regridder`), and
`S2Kernel.jl` answers the package's geometry generics from those same charts:
`cell_boundary`, `cell_center`, `cell_cap`, `cell_polygon_unitsphere` and
`cell_polygon(S2DGGS(), level, id)` all work. The `id` they take is the
**scaffold ordinal** `face * 4^level + hilbert_position` below — equivalently
`S2.HilbertOrder`'s data position minus one at `nside = 2^level` — and the two
paths are bitwise identical, sharing one evaluation of the chart.

**Dense geometry enumeration is wired, but the hierarchy is not.** `num_cells`
and `ordinal_to_cell` enumerate the scaffold ids accepted by the geometry
methods. `cell_children`, `cell_parent`, `cell_descendants`, `cell_to_ordinal`,
`descendant_range` and `root_ids` still throw `NotPortedError`, so there is no
`DGGSGrid(S2DGGS(), level)` either. The blocker is an id-model *decision*, not
code: `canonical_index_name` here is `:s2_cellid`, the native 64-bit
face/Hilbert/sentinel encoding, which is not ported (see the first note below).
Over scaffold ordinals the whole radix-4 group would be exact and free —
children `4p:4p+3`, ancestor `p ÷ 4^Δ`, descendants
`[p * 4^Δ, (p + 1) * 4^Δ)`, all derivable from `has_ordinal_ids = true` plus the
already-wired `root_count = 6` / `radix = 4` — but declaring it would silently
make `S2DGGS` ids mean something other than what `canonical_index_name` says.
The two must be settled together.

Tree notes:
- The scaffold ordinal is face * 4^level + hilbert_position.
- Native S2CellId uses a 64-bit face/Hilbert/sentinel encoding.

Sources:
- https://s2geometry.io/devguide/s2cell_hierarchy.html
- https://github.com/google/s2geometry
- https://pkg.go.dev/github.com/golang/geo/s2

Local references:
- global_grid_systems_report.md#12-s2-geometry
- docs/dggs-partial-coverage-storage.md

Notes:
- Port S2CellId orientation/range_min/range_max before using native ids for ranges. That port is also what unblocks the hierarchy group above; the Hilbert tables in `S2/chart.jl` transcribe s2geometry's conventions and are *intended* to agree with native ids, but there are no s2geometry fixtures here to verify it against.
- Cell polygons are geodesic quadrilaterals on the unit sphere — exactly, not to 4-gon accuracy, since a chart line `u = const` is a central plane section.
"""
struct S2DGGS <: AbstractDGGS end

system_name(::S2DGGS) = :S2
grid_family(::S2DGGS) = :cube_quadtree
base_solid(::S2DGGS) = :cube
cell_shape(::S2DGGS) = :spherical_quadrilateral
is_equal_area(::S2DGGS) = false
aperture(::S2DGGS) = 4
canonical_index_name(::S2DGGS) = :s2_cellid
max_level(::S2DGGS) = 30
supports_prefix_ranges(::S2DGGS) = true
root_count(::S2DGGS) = 6
radix(::S2DGGS) = 4
