"""
    H3DGGS()

H3 — icosahedral aperture-7 hexagonal grid (hexagons plus 12 pentagons),
indexed by the 64-bit `h3_index`, levels 0–15. Not equal-area.

Report section 1.1. Storage model `:sorted_ids`; native tree strategy
`:logical_aperture7` (122 roots, radix 7, prefix ranges).

Tree notes:
- Prefix ranges are logical H3 hierarchy ranges, not exact geographic containment.
- Pentagon deleted subsequences and Class II/Class III rotations must be handled by H3 math.

Sources:
- https://h3geo.org/docs/core-library/overview/
- https://h3geo.org/docs/highlights/indexing/
- https://h3geo.org/docs/core-library/coordsystems/
- https://github.com/uber/h3

Local references:
- global_grid_systems_report.md#11-h3

Notes:
- Implemented under `src/H3/` with H3_jll-backed native calls, a DimensionalData lookup, and the generic kernel wiring that feeds the ConservativeRegridding trees.
- Final polygon predicates must use exact H3 cell boundaries; logical parent prefixes are only hierarchy filters.
"""
struct H3DGGS <: AbstractDGGS end

system_name(::H3DGGS) = :H3
grid_family(::H3DGGS) = :icosahedral_hex
base_solid(::H3DGGS) = :icosahedron
cell_shape(::H3DGGS) = :hexagon_with_12_pentagons
is_equal_area(::H3DGGS) = false
aperture(::H3DGGS) = 7
canonical_index_name(::H3DGGS) = :h3_index
max_level(::H3DGGS) = 15
supports_prefix_ranges(::H3DGGS) = true
root_count(::H3DGGS) = 122
radix(::H3DGGS) = 7
