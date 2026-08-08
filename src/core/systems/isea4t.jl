"""
    ISEA4TDGGS()

ISEA4T — icosahedral aperture-4 equal-area triangular grid, indexed by
`face_quadtree_ordinal`, unbounded level.

Report section 1.7. Storage model `:sorted_ids`; native tree strategy
`:triangular_quadtree` (20 roots, radix 4, prefix ranges).

Tree notes:
- 20 icosahedron triangle roots; each triangle has four child triangles.
- Snyder Equal Area inverse projection still needs to be ported.

Sources:
- https://github.com/riskaware-ltd/open-eaggr
- an independent external implementation, used only as a black-box validation oracle (its recorded output lives in the test fixtures)

Local references:
- global_grid_systems_report.md#17-isea4t

Notes:
- Good candidate for the generic prefix partial tree once polygon math is ported.
"""
struct ISEA4TDGGS <: AbstractDGGS end

system_name(::ISEA4TDGGS) = :ISEA4T
grid_family(::ISEA4TDGGS) = :isea_triangle
base_solid(::ISEA4TDGGS) = :icosahedron
cell_shape(::ISEA4TDGGS) = :triangle
is_equal_area(::ISEA4TDGGS) = true
aperture(::ISEA4TDGGS) = 4
canonical_index_name(::ISEA4TDGGS) = :face_quadtree_ordinal
max_level(::ISEA4TDGGS) = nothing
supports_prefix_ranges(::ISEA4TDGGS) = true
root_count(::ISEA4TDGGS) = 20
radix(::ISEA4TDGGS) = 4
