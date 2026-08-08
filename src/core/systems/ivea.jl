"""
    IVEADGGS(variant::Symbol = :IVEA_family)

IVEA — the DGGAL family of icosahedral equal-area grids, indexed by
`dggal_zone`, unbounded level. `variant` selects a member of the family and is
the system's [`system_name`](@ref); it drives [`aperture`](@ref)
(`:IVEA4R` → 4, `:IVEA9R` → 9, `:IVEA3H` → 3, `:IVEA7H`/`:IVEA7H_Z7` → 7,
anything else → `:family`), [`radix`](@ref), and [`cell_shape`](@ref)
(`:rhomb` for the rhombic variants, `:hexagon_or_rhomb` otherwise). The default
`:IVEA_family` describes the family as a whole and pins no aperture.

Report section 1.9. Storage model `:sorted_ids`; native tree strategy
`:dggal_ivea` (root count not verified, radix = aperture when the variant pins
one, no prefix ranges).

Tree notes:
- IVEA uses Icosahedral Vertex-oriented Great-circle Equal Area projection.
- Use DGGAL zone ids and geometry as the reference contract.

Sources:
- https://dggal.org/docs/html/dggal.html
- https://github.com/ecere/dggal

Local references:
- global_grid_systems_report.md#19-ivea-family

Notes:
- Share aperture-specific tree code with ISEA once projection and root layout are pinned.
"""
struct IVEADGGS <: AbstractDGGS
    variant::Symbol
end

IVEADGGS() = IVEADGGS(:IVEA_family)

system_name(system::IVEADGGS) = system.variant
grid_family(::IVEADGGS) = :ivea
base_solid(::IVEADGGS) = :icosahedron
cell_shape(system::IVEADGGS) =
    system.variant in (:IVEA4R, :IVEA9R) ? :rhomb : :hexagon_or_rhomb
is_equal_area(::IVEADGGS) = true

function aperture(system::IVEADGGS)
    variant = system.variant
    return variant === :IVEA4R ? 4 :
           variant === :IVEA9R ? 9 :
           variant === :IVEA3H ? 3 :
           variant in (:IVEA7H, :IVEA7H_Z7) ? 7 : :family
end

canonical_index_name(::IVEADGGS) = :dggal_zone
max_level(::IVEADGGS) = nothing
supports_prefix_ranges(::IVEADGGS) = false

# No `root_count` method: the DGGAL root layout is not verified yet, so it falls
# back to the `AbstractDGGS` method and throws `NotPortedError`.

function radix(system::IVEADGGS)
    a = aperture(system)
    a isa Int || throw(NotPortedError(system_name(system), :radix,
        "The child radix is not a fixed verified integer for this system."))
    return a
end
