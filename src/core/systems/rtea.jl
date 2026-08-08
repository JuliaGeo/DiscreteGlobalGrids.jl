"""
    RTEADGGS(variant::Symbol = :RTEA_family)

RTEA — the DGGAL family of rhombic-triacontahedron equal-area grids, indexed by
`dggal_zone`, unbounded level. `variant` selects a member of the family and is
the system's [`system_name`](@ref); it drives [`aperture`](@ref)
(`:RTEA4R` → 4, `:RTEA9R` → 9, `:RTEA3H` → 3, `:RTEA7H`/`:RTEA7H_Z7` → 7,
anything else → `:family`), [`radix`](@ref), and [`cell_shape`](@ref)
(`:rhomb` for the rhombic variants, `:hexagon_or_rhomb` otherwise). The default
`:RTEA_family` describes the family as a whole and pins no aperture.

Report section 1.10. Storage model `:sorted_ids`; native tree strategy
`:dggal_rtea` (30 roots, radix = aperture when the variant pins one, no prefix
ranges).

Tree notes:
- The geometric solid has 30 rhombic faces; DGGAL data order still needs to be authoritative.

Sources:
- https://dggal.org/docs/html/dggal.html
- https://github.com/ecere/dggal

Local references:
- global_grid_systems_report.md#110-rtea-family

Notes:
- Rhombic variants should use prefix partial trees after DGGAL zone ordering is pinned.
"""
struct RTEADGGS <: AbstractDGGS
    variant::Symbol
end

RTEADGGS() = RTEADGGS(:RTEA_family)

system_name(system::RTEADGGS) = system.variant
grid_family(::RTEADGGS) = :rtea
base_solid(::RTEADGGS) = :rhombic_triacontahedron
cell_shape(system::RTEADGGS) =
    system.variant in (:RTEA4R, :RTEA9R) ? :rhomb : :hexagon_or_rhomb
is_equal_area(::RTEADGGS) = true

function aperture(system::RTEADGGS)
    variant = system.variant
    return variant === :RTEA4R ? 4 :
           variant === :RTEA9R ? 9 :
           variant === :RTEA3H ? 3 :
           variant in (:RTEA7H, :RTEA7H_Z7) ? 7 : :family
end

canonical_index_name(::RTEADGGS) = :dggal_zone
max_level(::RTEADGGS) = nothing
supports_prefix_ranges(::RTEADGGS) = false
root_count(::RTEADGGS) = 30

function radix(system::RTEADGGS)
    a = aperture(system)
    a isa Int || throw(NotPortedError(system_name(system), :radix,
        "The child radix is not a fixed verified integer for this system."))
    return a
end
