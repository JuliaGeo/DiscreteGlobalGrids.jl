"""
    HEALPixDGGS()

HEALPix — equal-area curvilinear quadrilaterals over the 12 HEALPix base faces,
`nested` indexing, levels 0–29.

Report section 1.13. Storage model `:dense_faces_or_sorted_ids`; native tree
strategy `:twelve_root_quadtree` (12 roots, radix 4, prefix ranges).

Tree notes:
- nside = 2^level; global leaf count is 12 * 4^level.
- Nested pixel p has children 4p, 4p+1, 4p+2, 4p+3.

Sources:
- https://healpix.sourceforge.io/html/intro_Geometric_Algebraic_Propert.htm
- https://healpix.sourceforge.io/
- https://github.com/healpy/healpy

Local references:
- src/HEALPix/HealpixLookups.jl
- src/HEALPix/HealpixKernel.jl
- src/HEALPix/chart.jl
- src/HEALPix/face_grid.jl
- ConservativeRegridding.jl/ext/ConservativeRegriddingHealpixExt.jl
- ConservativeRegridding.jl/ext/ConservativeRegriddingRingGridsExt/healpix.jl

Notes:
- Most complete local implementation today.
- Exact cell polygons come from the package's own chart kernel
  (`src/HEALPix/chart.jl`): `pixel_corners` / `pixel_center` are closed forms,
  free of Healpix.jl, and are the single evaluation shared by the id-hierarchy
  cells here and the dense face grid in `src/HEALPix/face_grid.jl`.
"""
struct HEALPixDGGS <: AbstractDGGS end

system_name(::HEALPixDGGS) = :HEALPix
grid_family(::HEALPixDGGS) = :healpix
base_solid(::HEALPixDGGS) = :healpix_12_faces
cell_shape(::HEALPixDGGS) = :curvilinear_quadrilateral
is_equal_area(::HEALPixDGGS) = true
aperture(::HEALPixDGGS) = 4
canonical_index_name(::HEALPixDGGS) = :nested
# The nested id space is unbounded in principle, but the package addresses it
# with `Int64`: `leaf_count = 12 * 4^level` already overflows at level 30
# (12 * 4^30 > typemax(Int64)), and the kernel's `radix^delta` range/ordinal
# arithmetic wraps silently rather than throwing, so 29 is the bound.
max_level(::HEALPixDGGS) = 29
supports_prefix_ranges(::HEALPixDGGS) = true
root_count(::HEALPixDGGS) = 12
radix(::HEALPixDGGS) = 4
