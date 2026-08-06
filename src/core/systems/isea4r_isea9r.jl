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
- The ten-diamond *layout* is shared with ISEA9R: DGGAL's "5x6 Cartesian equal-area square-zone" phrasing describes a container CRS, and OGC 21-038r1 Annex B.2 names ten root rhombuses. `ISEA9R` therefore imports this system's charts unchanged. Identifier compatibility with DGGAL remains unclaimed for both — see [`ISEA9RDGGS`](@ref) and `docs/design/isea9r_layout.md`.

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
strategy `:rhombic_aperture9` (10 roots, radix 9, prefix ranges).

# What is built today

**The dense diamond grid, and geometry over the canonical ordinal** — the same
two things `ISEA4RDGGS` has, over the same ten charts.
`DiscreteGlobalGrids.ISEA9R` owns whole-sphere grids at any `nside >= 1`
(`ISEA9R.Isea9rFaceGrid(nside; ordering)` → `treeify` → `Regridder`) and
`Isea9rKernel.jl` answers the package's geometry generics from the charts:
`cell_boundary`, `cell_center`, `cell_cap`, `cell_polygon_unitsphere` and
`cell_polygon(ISEA9RDGGS(), level, id)` all work. The `id` they take is the
`isea9r_ordinal` — equivalently `ISEA9R.MortonOrder`'s data position minus one
at `nside = 3^level` — and the two paths are bitwise identical, sharing one
evaluation of the chart.

The charts are not a second implementation: `ISEA9R` imports `ISEA4R`'s
`xyd_to_point` / `cell_corners` by name, the rhombus chart being
aperture-agnostic (it takes continuous `(x, y)` and quantises nothing). So an
ISEA9R cell and an ISEA4R cell at the same lattice position and the same
`nside = 3^k` are the same four `Float64`-triples, bitwise.

`ISEA9R.MortonOrder` realizes the `isea9r_ordinal` shape
`diamond * 9^level + position` at `nside = 3^level`, and **`position` is hereby
pinned to the base-9 Morton (Z-order) code**: base-9 digit `k` of the code is
`ix_k + 3 * iy_k`, the base-3 digits of the two coordinates interleaved. That is
this package's canonical choice in this package's own diamond numbering, and is
*not* DGGAL's index — see the compatibility section below.

**The id hierarchy is deferred, not blocked**, on exactly the ISEA4R sibling's
line and as one decision with it. `cell_children`, `cell_parent`,
`cell_descendants`, `cell_to_ordinal`, `ordinal_to_cell`, `descendant_range`,
`num_cells` and `root_ids` still throw `NotPortedError`, so there is no
`DGGSGrid(ISEA9RDGGS(), level)` either. The radix-9 arithmetic over these
ordinals is exact (the base-9 Morton code drops one digit per level up and the
lattice nesting is bit-exact, `fl(ix/n) === fl(3ix/3n)`), so children `9p:9p+8`,
ancestor `p ÷ 9^Δ` and descendant interval `[p * 9^Δ, (p + 1) * 9^Δ)` would all
follow from `has_ordinal_ids = true` plus the already-wired `root_count = 10` /
`radix = 9`. What that one line additionally owes is the kernel-test battery the
resulting grids deserve and a `has_exact_subtree_cap` decision.

# The layout question is settled; the DGGAL-compatibility question is not

The registry note this entry used to carry — *"DGGAL describes ISEA9R as a 5x6
Cartesian equal-area square-zone layout; do not apply the SST 10-root layout to
ISEA9R without fixtures"* — conflated two questions. The first is **resolved**
from primary sources:

> "The ten root rhombuses are formed by combining two icosahedron triangles at
> their base." — OGC 21-038r1 (*OGC API — DGGS Part 1: Core*), Annex B.2,
> Listing B.2, <https://docs.ogc.org/is/21-038r1/21-038r1.html#isea9r-dggrs>

and DGGAL's `RhombicIcosahedral9R::countZones(level)` returns `10 * 9^level`
(`src/dggrs/RI9R.ec`), i.e. ten zones at level 0. The 5×6 space is a *container*
CRS chosen so the grid is also an OGC 2D Tile Matrix Set: the sphere occupies
ten of its thirty unit cells in a diagonal staircase, and each of those unit
squares is one icosahedral rhombus — two triangular faces glued along an edge,
which is precisely the ten-diamond layout. So the ten-diamond chart *is* ISEA9R
geometry, and it is wired.

The second question stands unchanged: **no DGGAL/DGGRID identifier or geometry
compatibility is claimed, and fixtures are still required before any could be.**
Three separate deltas would each have to close:

1. *Identifiers.* DGGAL's zone id is `{LevelChar}{RootRhombus}-{HexIndex}` with
   the index row-major within the rhombus; its root numbering runs
   north/south/north/… around the staircase where this package numbers the five
   northern diamonds first, and its in-square axes are the transpose of this
   package's. Both the permutation and the transpose are *derivations*, not
   oracle output.
2. *Orientation.* DGGAL places the icosahedron at 11.20°E; this package's
   DGGRID-standard placement is 11.25°E — about 5.6 km at the equator.
3. *Ellipsoid.* DGGAL converts geodetic↔authalic latitude at the WGS84
   boundary; `ISEA` works on the authalic sphere and does not.

`docs/design/isea9r_layout.md` carries the citations, the derived permutation,
and the seven-item fixture dump that would settle all three.

Tree notes:
- 10 roots and radix 9, both wired; `supports_prefix_ranges` is true of the package-canonical base-9 Morton ordinal, not of DGGAL's zone id.
- DGGAL's "5x6 Cartesian equal-area square-zone" phrasing describes a container CRS, not a 30-root decomposition (OGC 21-038r1 Annex B.2; DGGAL RI9R.ec).

Sources:
- https://docs.ogc.org/is/21-038r1/21-038r1.html#isea9r-dggrs
- https://dggal.org/docs/html/dggal.html
- https://github.com/ecere/dggal
- https://github.com/meggart/SphericalSpatialTrees.jl/blob/main/src/iseatree.jl
- https://github.com/meggart/SphericalSpatialTrees.jl/blob/main/src/nativeisea.jl
- an independent external implementation, used only as a black-box validation oracle (its recorded output lives in the test fixtures)

Local references:
- global_grid_systems_report.md#18-isea4r--isea9r
- docs/design/isea9r_layout.md
- docs/design/isea4r_diamond_layout.md
- ConservativeRegridding.jl/examples/isea20_sst.jl

Notes:
- The concrete backend is the package's own `ISEA9R` submodule over `ISEA4R`'s charts, not DGGAL; verify naming against DGGAL fixtures before declaring identifier compatibility.
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

"""
    supports_prefix_ranges(::ISEA9RDGGS) -> true

True **of the package-canonical `isea9r_ordinal`**, and of nothing else.

That ordinal is `diamond * 9^level + morton9_position` with `position` the
base-9 Morton code (`ISEA9R.MortonOrder`). Over it, a cell's `leaf_level`
descendants are the contiguous interval `[p * 9^Δ, (p + 1) * 9^Δ)`: the Morton
code drops its low base-9 digit per level up, and the within-diamond lattice
nesting is bit-exact (`fl(ix/n) === fl(3ix/3n)`, since all four integers are
exactly representable as `Float64` and the real quotients are identical), so
the interval arithmetic `leaf_interval` performs is exact rather than
approximate.

It is emphatically **not** true of DGGAL's ISEA9R identifiers. DGGAL packs a
zone as `level << 59 | row << 30 | col` over the global 5×6 lattice (its textual
form `{LevelChar}{RootRhombus}-{HexIndex}` reorders the same information), and
the descendants of a zone are a *rectangular block* in `(row, col)` — never a
contiguous integer interval in either encoding. DGGAL itself compacts by
explicit nine-child set membership rather than by interval arithmetic
(`compactI9RZones`, `src/dggrs/RI9R.ec`). Anything reading this trait as a
statement about DGGAL ids is reading it wrong; see
`docs/design/isea9r_layout.md`.
"""
supports_prefix_ranges(::ISEA9RDGGS) = true

"""
    root_count(::ISEA9RDGGS) -> 10

Ten root rhombuses, each two icosahedron faces glued along an edge. Normative:

> "The ten root rhombuses are formed by combining two icosahedron triangles at
> their base." — OGC 21-038r1 (*OGC API — Discrete Global Grid Systems — Part 1:
> Core*), Annex B.2 "ISEA9R DGGRS definition", Listing B.2,
> <https://docs.ogc.org/is/21-038r1/21-038r1.html#isea9r-dggrs>

and implemented: DGGAL's `RhombicIcosahedral9R::countZones(level)` returns
`10 * 9^level` (`src/dggrs/RI9R.ec`), so `countZones(0) == 10`, and its zone-id
codec admits exactly `root ∈ 0:9`.

The "5×6 Cartesian equal-area square-zone" phrasing in DGGAL's own
documentation — which an earlier version of this registry entry read as a
possible thirty-root decomposition — describes a *container* CRS for tiling, not
a root decomposition: the sphere occupies ten of the thirty unit cells in a
diagonal staircase and the other twenty are empty. See
`docs/design/isea9r_layout.md`.
"""
root_count(::ISEA9RDGGS) = 10

radix(::ISEA9RDGGS) = 9
