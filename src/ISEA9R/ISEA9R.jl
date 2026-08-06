"""
    ISEA9R

ISEA9R diamond chart grids: the icosahedron's ten rhombi at aperture 9, and
dense whole-sphere grids built on them, under a swappable data ordering.

Neither the projection nor the chart is new here. The projection is
[`ISEA`](@ref)'s Snyder equal-area machinery and the chart is
[`ISEA4R`](@ref)'s ten-diamond rhombus chart, **imported unchanged** — the
rhombus chart takes continuous `(x, y)` and quantises nothing, so it carries no
aperture, and the same ten charts serve `nside = 3^level` exactly as they serve
`nside = 2^level` (`chart.jl` records the delegation and its basis;
`test/ISEA9R/test_delegation.jl` checks it bitwise). What this module adds is
the aperture-9 half: the base-9 Morton index maps, the two orderings, the system
singleton, and geometry for `ISEA9RDGGS` over the canonical `isea9r_ordinal`.

# Why the ten-diamond layout, and on whose authority

The question `docs/design/isea4r_diamond_layout.md` §7 left open — whether
DGGAL's "5×6 Cartesian equal-area square-zone" ISEA9R is the ten-diamond layout
at all — is resolved, from primary sources, in favour of ten roots:

> "The ten root rhombuses are formed by combining two icosahedron triangles at
> their base." — OGC 21-038r1 (*OGC API — DGGS Part 1: Core*), Annex B.2,
> Listing B.2, <https://docs.ogc.org/is/21-038r1/21-038r1.html#isea9r-dggrs>

and DGGAL's own implementation computes `10 * 9^level` zones at every level
(`src/dggrs/RI9R.ec`, `RhombicIcosahedral9R::countZones`). The 5×6 space is a
*container* CRS for tiling, not a root decomposition: the sphere occupies ten of
its thirty unit cells in a diagonal staircase, and each of those unit squares is
one icosahedral rhombus. The full record, with the citations and with everything
this package deliberately does *not* claim, is
`docs/design/isea9r_layout.md`.

# Scope, honestly

**The dense diamond grid, plus chart geometry for `ISEA9RDGGS`.**
`face_grid.jl` — [`Isea9rFaceSpace`](@ref) and [`Isea9rFaceGrid`](@ref) — owns
dense whole-sphere grids at one resolution, where the data vector *is* the grid:
any `nside >= 1`, laid out in row-major or base-9 Morton order by a swappable
[`AbstractIsea9rOrdering`](@ref) component that no grid or tree code has to know
about. `Isea9rKernel.jl` then answers the package's geometry generics for
`ISEA9RDGGS` over the canonical `isea9r_ordinal`
`diamond * 9^level + morton9_position`, so `cell_polygon(ISEA9RDGGS(), level,
id)` returns a real polygon — bitwise the one the diamond grid emits for the
same ordinal.

The `isea9r_ordinal` id *hierarchy* is still deferred: `cell_children`,
`cell_parent`, `descendant_range`, `num_cells` and the rest of that group throw
`NotPortedError`, exactly as for the ISEA4R sibling and for the same reason
(scope, not obstruction — the radix-9 arithmetic over these ordinals is exact,
which is what `supports_prefix_ranges(ISEA9RDGGS())` asserts at the interface
level). What [`MortonOrder`](@ref) does assert today is that it is a bijection
and that it realizes the registry's ordinal `diamond * 9^level + position` with
`position := morton9(ix, iy)`.

# Provenance, and what is NOT claimed

The *layout* — ten roots, each a rhombus of two icosahedron faces — is the
standard's, cited above. The *numbering* of those ten, the in-diamond axis
orientation and the in-diamond index are this package's own conventions,
inherited from `ISEA4R` and pinned by nothing external. **DGGAL identifier and
geometry compatibility is deliberately not claimed and must not be inferred**,
and three separate deltas would each have to be closed before it could be:

1. DGGAL numbers its roots in a staircase that alternates north/south, and
   indexes within a root in row-major order over the *transpose* of this
   package's square. Both the permutation and the transpose are derivations
   with no fixture behind them.
2. DGGAL's icosahedron sits at 11.20°E; this package's DGGRID-standard
   placement sits at 11.25°E — about 5.6 km at the equator, six orders of
   magnitude above any tolerance in `test/ISEA9R/`.
3. DGGAL converts geodetic↔authalic latitude at the WGS84 boundary; `ISEA`
   works on the authalic sphere and does not.

`docs/design/isea9r_layout.md` §7 lists the seven-item fixture dump that would
settle all three.

Corner rings are counterclockwise as seen from outside the sphere. That is a
contract, not a convention: the convex-clip kernel behind spherical
intersections clips a clockwise ring to EMPTY, so a reversed ring yields silent
zero areas instead of an error. Include order follows the dependency — the
delegation record and index maps, the face grid, then the kernel wiring built on
them.
"""
module ISEA9R

# No `using ..ISEA`: this module reaches the projection only through `ISEA4R`'s
# chart, which `chart.jl` imports name by name.

# Where the chart comes from (`ISEA4R`, unchanged) + the base-9 index maps.
include("chart.jl")
# Dense per-resolution diamond grids built on those charts.
include("face_grid.jl")
# Geometry-only operations-kernel wiring for `ISEA9RDGGS` (see
# `src/core/kernel.jl`); last, because it evaluates the chart and documents
# itself against the grid.
include("Isea9rKernel.jl")

# The dense face-grid layer's contract surface. System vocabulary stays in the
# submodule — nothing here is re-exported from `DiscreteGlobalGrids` — so reach
# for it as `using DiscreteGlobalGrids.ISEA9R`. The chart functions themselves
# stay unexported (`ISEA9R.xyd_to_morton`, ...): they are the kernel these types
# are built out of, not the API. `cell_polygon`, `cell_center` and `num_cells`
# stay unexported for the additional reason that the package namespace already
# has those names.
export Isea9rFaceSpace, Isea9rFaceGrid
export AbstractIsea9rOrdering, RowMajorOrder, MortonOrder

end # module ISEA9R
