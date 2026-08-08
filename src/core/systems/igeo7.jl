"""
    IGEO7DGGS()

IGEO7 (ISEA7H + Z7) — icosahedral aperture-7 equal-area hexagonal grid, indexed
by `ri7h_u64`, levels 0–19.

Report section 1.4. Storage model `:sorted_ri7h_ids`; native tree strategy
`:ri7h_aperture7` (12 roots, radix 7, prefix ranges).

Tree notes:
- Native implementation lives in `src/IGeo7/` (`z7.jl`, `engine.jl`, `grid.jl`).
- Uses DGGAL-compatible RI7H fields internally; public Z7 digit-string compatibility against recorded oracle output still needs cross-edge fixtures.
- Snyder inverse projection is adapted from SphericalSpatialTrees.NativeISEA because the package currently has a DimensionalData compat conflict here.

Sources:
- https://agile-giss.copernicus.org/articles/6/32/2025/
- https://igeo7.org/
- an independent external implementation, used only as a black-box validation oracle (its recorded output lives in the test fixtures)

Local references:
- global_grid_systems_report.md#14-igeo7--isea7h_z7

Notes:
- Add DGGAL and recorded-oracle fixture comparison for public Z7 text ids and edge/polar ordering.
- Consider replacing the vendored Snyder slice with a SphericalSpatialTrees dependency once its DimensionalData compat permits it.
"""
struct IGEO7DGGS <: AbstractDGGS end

system_name(::IGEO7DGGS) = :IGEO7_ISEA7H_Z7
grid_family(::IGEO7DGGS) = :isea_hex
base_solid(::IGEO7DGGS) = :icosahedron
cell_shape(::IGEO7DGGS) = :hexagon
is_equal_area(::IGEO7DGGS) = true
aperture(::IGEO7DGGS) = 7
canonical_index_name(::IGEO7DGGS) = :ri7h_u64
# The Z7 encoding carries 19 digit slots after the base cell (`IGeo7.z7.jl`'s
# `MAX_RESOLUTION`); slot 20 is the all-7 padding sentinel, which has no
# geometry, so 19 is the deepest addressable level.
max_level(::IGEO7DGGS) = 19
supports_prefix_ranges(::IGEO7DGGS) = true
# Verified: `IGeo7.res0_cells()` enumerates exactly twelve ascending base cells,
# one per icosahedron vertex, pinned against `test/IGeo7/vectors/res0_cells.csv`
# and `num_cells.csv` (`10 * 7^r + 2`, which is 12 at r = 0).
root_count(::IGEO7DGGS) = 12
radix(::IGEO7DGGS) = 7
