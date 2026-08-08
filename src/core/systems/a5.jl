"""
    A5DGGS()

A5 — dodecahedral equal-area pentagonal grid, indexed by `a5_u64`, levels 0–30.
The refinement aperture is implementation-defined.

Report section 1.3. Storage model `:sorted_ids`; native tree strategy
`:external_a5_hierarchy` (12 roots, no fixed radix, no prefix ranges).

Tree notes:
- Hierarchy from upstream serialization: world cell, 12 res0 pentagons, 5 res1 quintants per pentagon, then 4-way Hilbert children from res2 onward.
- Generic fixed-radix prefix ranges do not apply until numeric descendant ranges are fixture-verified.
- Current Julia integration lives in src/A5/A5Native.jl, src/A5/A5Lookups.jl, and src/A5/A5Kernel.jl.
- Authoritative source lives in felixpalmer/a5 modules/core, modules/lattice, modules/geometry, modules/projections.

Sources:
- https://a5geo.org/
- https://github.com/felixpalmer/a5
- https://github.com/felixpalmer/a5/tree/main/modules/core
- https://github.com/felixpalmer/a5/tree/main/modules/projections

Local references:
- global_grid_systems_report.md#13-a5

Notes:
- Native projection, boundary, point lookup, serialization, parent, and children math are ported.
- Cell ids are bigint/u64 in the upstream implementation.
"""
struct A5DGGS <: AbstractDGGS end

system_name(::A5DGGS) = :A5
grid_family(::A5DGGS) = :dodecahedral_pentagonal
base_solid(::A5DGGS) = :dodecahedron
cell_shape(::A5DGGS) = :pentagon
is_equal_area(::A5DGGS) = true
aperture(::A5DGGS) = :implementation_defined
canonical_index_name(::A5DGGS) = :a5_u64
max_level(::A5DGGS) = 30
supports_prefix_ranges(::A5DGGS) = false
root_count(::A5DGGS) = 12
# No `radix` method: A5 has no fixed verified radix, so it falls back to the
# `AbstractDGGS` method and throws `NotPortedError`.
