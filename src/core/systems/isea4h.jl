"""
    ISEA4HDGGS()

ISEA4H — icosahedral aperture-4 equal-area hexagonal grid; the canonical index
is implementation-defined and the level is unbounded.

Report section 1.6. Storage model `:sorted_ids`; native tree strategy
`:hex_aperture4` (root count not verified, radix 4, no prefix ranges).

Tree notes:
- Aperture-4 area scaling is clean, but hex cells are not a dense quadtree matrix.

Sources:
- an independent external implementation, used only as a black-box validation oracle (its recorded output lives in the test fixtures); R and Python wrappers around that same oracle exist too

Local references:
- global_grid_systems_report.md#16-isea4h

Notes:
- Use sorted ids plus exact implementation compact/uncompact.
"""
struct ISEA4HDGGS <: AbstractDGGS end

system_name(::ISEA4HDGGS) = :ISEA4H
grid_family(::ISEA4HDGGS) = :isea_hex
base_solid(::ISEA4HDGGS) = :icosahedron
cell_shape(::ISEA4HDGGS) = :hexagon
is_equal_area(::ISEA4HDGGS) = true
aperture(::ISEA4HDGGS) = 4
canonical_index_name(::ISEA4HDGGS) = :implementation_defined
max_level(::ISEA4HDGGS) = nothing
supports_prefix_ranges(::ISEA4HDGGS) = false
# No `root_count` method: the root zone count is not verified yet, so it falls
# back to the `AbstractDGGS` method and throws `NotPortedError`.
radix(::ISEA4HDGGS) = 4
