"""
    ISEA3HDGGS()

ISEA3H — icosahedral aperture-3 equal-area hexagonal grid (hexagons with
pentagons); the canonical index is implementation-defined and the level is
unbounded.

Report section 1.5. Storage model `:sorted_ids`; native tree strategy
`:hex_aperture3` (root count not verified, radix 3, no prefix ranges).

Tree notes:
- Aperture-3 hex orientation alternates; use implementation parent/child functions.

Sources:
- https://github.com/mocnik-science/geogrid
- https://github.com/ecere/dggal
- an independent external implementation, used only as a black-box validation oracle (its recorded output lives in the test fixtures)

Local references:
- global_grid_systems_report.md#15-isea3h

Notes:
- Prefer DGGAL/geogrid and recorded oracle fixtures before native Snyder Equal Area port.
"""
struct ISEA3HDGGS <: AbstractDGGS end

system_name(::ISEA3HDGGS) = :ISEA3H
grid_family(::ISEA3HDGGS) = :isea_hex
base_solid(::ISEA3HDGGS) = :icosahedron
cell_shape(::ISEA3HDGGS) = :hexagon_with_pentagons
is_equal_area(::ISEA3HDGGS) = true
aperture(::ISEA3HDGGS) = 3
canonical_index_name(::ISEA3HDGGS) = :implementation_defined
max_level(::ISEA3HDGGS) = nothing
supports_prefix_ranges(::ISEA3HDGGS) = false
# No `root_count` method: the root zone count is not verified yet, so it falls
# back to the `AbstractDGGS` method and throws `NotPortedError`.
radix(::ISEA3HDGGS) = 3
