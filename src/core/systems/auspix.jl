"""
    AusPIXDGGS()

AusPIX — the Geoscience Australia rHEALPix profile: equal-area aperture-9
squares on a cube-like base, indexed by `auspix_id_with_rhealpix_ordinal`,
unbounded level.

Report section 1.12. Storage model `:dense_faces_or_sorted_ids`; native tree
strategy `:auspix_rhealpix_base9` (6 roots, radix 9, prefix ranges).

Tree notes:
- Geometry follows rHEALPix; AusPIX adds naming, metadata, and dataset conventions.

Sources:
- https://github.com/GeoscienceAustralia/AusPIX_DGGS
- https://github.com/GeoscienceAustralia/AusPIX-DGGS-dataset

Local references:
- global_grid_systems_report.md#112-auspix

Notes:
- Implement as an rHEALPix wrapper with an AusPIX id codec.
"""
struct AusPIXDGGS <: AbstractDGGS end

system_name(::AusPIXDGGS) = :AusPIX
grid_family(::AusPIXDGGS) = :rhealpix_auspix
base_solid(::AusPIXDGGS) = :cube_like
cell_shape(::AusPIXDGGS) = :square
is_equal_area(::AusPIXDGGS) = true
aperture(::AusPIXDGGS) = 9
canonical_index_name(::AusPIXDGGS) = :auspix_id_with_rhealpix_ordinal
max_level(::AusPIXDGGS) = nothing
supports_prefix_ranges(::AusPIXDGGS) = true
root_count(::AusPIXDGGS) = 6
radix(::AusPIXDGGS) = 9
