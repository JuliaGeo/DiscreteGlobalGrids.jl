"""
    AbstractISEARDGGS

Supertype of the two rhombic ISEA systems of report section 1.8,
[`ISEA4RDGGS`](@ref) and [`ISEA9RDGGS`](@ref).
"""
abstract type AbstractISEARDGGS <: AbstractDGGS end

"""
    ISEA4RDGGS()

ISEA4R — icosahedral aperture-4 equal-area rhombic grid, indexed by
`isea4r_ordinal`, unbounded level.

Report section 1.8. Storage model `:sorted_ids_or_dense_faces`; native tree
strategy `:rhombic_aperture4` (10 roots, radix 4, prefix ranges).

# What is built today

**The dense diamond grid, and geometry over the canonical ordinal.**
`DiscreteGlobalGrids.ISEA4R` owns the ten rhombus charts over `ISEA`'s Snyder
machinery plus whole-sphere grids at any `nside >= 1`
(`ISEA4R.Isea4rFaceGrid(nside; ordering)` → `treeify` → `Regridder`), and
`Isea4rKernel.jl` answers the package's geometry generics from those same
charts: `cell_boundary`, `cell_center`, `cell_cap`,
`cell_polygon_unitsphere` and `cell_polygon(ISEA4RDGGS(), level, id)` all work.
The `id` they take is the `isea4r_ordinal` — equivalently
`ISEA4R.MortonOrder`'s data position minus one at `nside = 2^level` — and the
two paths are bitwise identical, sharing one evaluation of the chart.

`ISEA4R.MortonOrder` realizes the `isea4r_ordinal` shape
`diamond * 4^level + position` at `nside = 2^level`, and **`position` is hereby
pinned to the Morton (Z-order) code with `ix` in the even bit positions and
`iy` in the odd ones**, in this package's own diamond layout. That layout has no
external oracle behind it, so *no identifier compatibility with any external
ISEA4R/ISEA9R product, DGGAL included, is claimed* — see
`docs/design/isea4r_diamond_layout.md`, which records the numbering, the
derivation rule, and exactly what a fixture-based permutation would have to
pin.

**Dense geometry enumeration is wired; the id hierarchy is deferred.**
`num_cells` and `ordinal_to_cell` enumerate the canonical Morton ids accepted
by geometry. `cell_children`, `cell_parent`, `cell_descendants`,
`cell_to_ordinal`, `descendant_range` and `root_ids` still throw
`NotPortedError`, so there is no `DGGSGrid(ISEA4RDGGS(), level)` either. Unlike
S2 there is no open question of
which id space is canonical — this one is, and the radix-4 arithmetic over it is
exact (the Morton code drops two bits per level up and the lattice nesting is
bit-exact, `fl(ix/n) === fl(2ix/2n)`), so children `4p:4p+3`, ancestor
`p ÷ 4^Δ` and descendant interval `[p * 4^Δ, (p + 1) * 4^Δ)` would all follow
from `has_ordinal_ids = true` plus the already-wired `root_count = 10` /
`radix = 4`. What that one line additionally *owes* is the kernel-test battery
the resulting grids deserve and a `has_exact_subtree_cap` decision, which is why
it lands as its own unit rather than as a side effect of the geometry wiring.

Tree notes:
- SphericalSpatialTrees.jl current GitHub main implements an ISEA10 diamond aperture-4 layout with 10 roots and 2^r x 2^r cells per root. `ISEA4R.RowMajorOrder` is layout-*isomorphic* to SST's `LinearIndices((2^r, 2^r, 10))` order, but the diamond numbering and per-diamond axis orientation are not cross-pinned; SST interop would be a permutation of that ordering fitted against SST fixtures, addable as a third ordering with zero tree changes.
- DGGAL describes ISEA9R as a 5x6 Cartesian equal-area square-zone layout; do not apply the SST 10-root layout to ISEA9R without fixtures.

Sources:
- https://dggal.org/docs/html/dggal.html
- https://github.com/ecere/dggal
- https://github.com/meggart/SphericalSpatialTrees.jl/blob/main/src/iseatree.jl
- https://github.com/meggart/SphericalSpatialTrees.jl/blob/main/src/nativeisea.jl
- an independent external implementation, used only as a black-box validation oracle (its recorded output lives in the test fixtures)

Local references:
- global_grid_systems_report.md#18-isea4r--isea9r
- docs/design/isea4r_diamond_layout.md
- ConservativeRegridding.jl/examples/isea20_sst.jl

Notes:
- The concrete backend is the package's own `ISEA4R` submodule (Snyder charts from `ISEA`, ten-diamond layout derived here), not SST.ISEACircleTree; verify naming against DGGAL/SST or recorded oracle fixtures before declaring identifier compatibility.
"""
struct ISEA4RDGGS <: AbstractISEARDGGS end

"""
    ISEA9RDGGS()

ISEA9R — icosahedral aperture-9 equal-area rhombic grid, indexed by
`isea9r_ordinal`, unbounded level.

Report section 1.8. Storage model `:sorted_ids_or_dense_faces`; native tree
strategy `:rhombic_aperture9` (root count not verified, radix 9, no prefix
ranges).

# What an ISEA9R face grid would need, and what blocks it

*If* the ten-diamond layout applies to ISEA9R, the geometry is already written:
the identical `ISEA4R.DIAMONDS` table and rhombus chart, which are
aperture-agnostic (the chart admits any `nside`). The delta over the shipped
ISEA4R face grid would be one ordering type of about thirty lines — `nside =
3^level`, `data_index = diamond * 9^level + base3_interleave(ix, iy) + 1`, or
row-major for arbitrary `nside`, restricted to `3^k` by `validate_ordering`
exactly as `ISEA4R.MortonOrder` restricts to `2^k` — plus its tests.

What blocks it is the layout question, not the code: DGGAL describes ISEA9R as
a 5x6 Cartesian equal-area square-zone layout, which may not be the ten-diamond
layout at all. Per the registry note below, the SST/ISEA4R 10-root layout must
not be applied to ISEA9R without DGGAL or other external fixtures. Until then
`root_count` stays unwired (`NotPortedError`), `supports_prefix_ranges` stays
`false`, and no ISEA9R face grid is built. If fixtures do pin the ten-diamond
reading, the package's own ISEA4R layout
(`docs/design/isea4r_diamond_layout.md`) is the geometry it would reuse
verbatim — with that note's compatibility caveat still applying to identifiers.

Tree notes:
- Rhombic variants should become prefix-range trees after the root layout is pinned against DGGAL or other external fixtures.
- DGGAL describes ISEA9R as a 5x6 Cartesian equal-area square-zone layout; do not apply the SST 10-root layout to ISEA9R without fixtures.

Sources:
- https://dggal.org/docs/html/dggal.html
- https://github.com/ecere/dggal
- https://github.com/meggart/SphericalSpatialTrees.jl/blob/main/src/iseatree.jl
- https://github.com/meggart/SphericalSpatialTrees.jl/blob/main/src/nativeisea.jl
- an independent external implementation, used only as a black-box validation oracle (its recorded output lives in the test fixtures)

Local references:
- global_grid_systems_report.md#18-isea4r--isea9r
- ConservativeRegridding.jl/examples/isea20_sst.jl

Notes:
- Verify root count/order before enabling prefix ranges.
"""
struct ISEA9RDGGS <: AbstractISEARDGGS end

system_name(::ISEA4RDGGS) = :ISEA4R
grid_family(::ISEA4RDGGS) = :isea_rhombic
base_solid(::ISEA4RDGGS) = :icosahedron
cell_shape(::ISEA4RDGGS) = :rhomb
is_equal_area(::ISEA4RDGGS) = true
aperture(::ISEA4RDGGS) = 4
canonical_index_name(::ISEA4RDGGS) = :isea4r_ordinal
max_level(::ISEA4RDGGS) = nothing
supports_prefix_ranges(::ISEA4RDGGS) = true
root_count(::ISEA4RDGGS) = 10
radix(::ISEA4RDGGS) = 4

system_name(::ISEA9RDGGS) = :ISEA9R
grid_family(::ISEA9RDGGS) = :isea_rhombic
base_solid(::ISEA9RDGGS) = :icosahedron
cell_shape(::ISEA9RDGGS) = :rhomb
is_equal_area(::ISEA9RDGGS) = true
aperture(::ISEA9RDGGS) = 9
canonical_index_name(::ISEA9RDGGS) = :isea9r_ordinal
max_level(::ISEA9RDGGS) = nothing
supports_prefix_ranges(::ISEA9RDGGS) = false
# No `root_count` method for ISEA9R: the root layout is not verified yet, so it
# falls back to the `AbstractDGGS` method and throws `NotPortedError`.
radix(::ISEA9RDGGS) = 9
