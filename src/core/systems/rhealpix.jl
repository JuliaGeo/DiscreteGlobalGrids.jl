"""
    RHEALPixDGGS()

rHEALPix — rectified HEALPix: an equal-area aperture-9 square grid on a
cube-like base, indexed by `root_base9_digits`, unbounded level.

Report section 1.11. Storage model `:dense_faces_or_sorted_ids`; native tree
strategy `:six_root_base9` (6 roots, radix 9, prefix ranges).

# What an rHEALPix face grid would need, and why it is deferred

Registry-only today: no geometry and no hierarchy are wired, so every kernel
operation throws `NotPortedError`. That is a disposition, not an oversight — an
rHEALPix face grid *is* expressible on the `FaceGridSystem` contract
(`src/core/face_grid.jl`), and this note records what it would take.

The shape is already there: 6 faces, and the four equatorial charts are the
HEALPix equatorial charts verbatim (`src/HEALPix/chart.jl`), since rHEALPix
rearranges the polar caps and leaves the equatorial belt alone. Each cap face is
a *piecewise* chart that assembles the four HEALPix polar-cap triangles into one
square by rigid planar motions — the same piecewise-assembly pattern
`ISEA4R.xyd_to_point` already carries for its two half-triangles, with the same
obligation to make the seam ownership deterministic so shared lattice points
come out bit-identical.

Three things block it, in order of weight:

1. **No pinned placement parameters.** The assembly is parameterized by
   `(north_square, south_square)` — which quadrant each cap's triangles are
   rotated into — and those must be pinned against `rhealpixdggs-py` fixtures
   before any claim about `root_base9_digits` identifiers can be made. Without
   that, a grid would be a plausible-looking layout with no interoperability,
   which is exactly the trap `docs/design/isea4r_diamond_layout.md` records for
   ISEA4R.
2. **Cap-face block caps are unsound across the diagonal seams.** The four
   corners of an index block that straddles a cap's diagonal do not bound the
   block, so `FourCornerCap` is not available there — while the current contract
   makes `cap_policy` a per-*system* choice, not a per-face one. An rHEALPix
   instance should therefore either take `PerimeterWalkCap()` globally (sound,
   available today, zero new contract surface, merely slower on the four
   equatorial faces that would not have needed it) or the contract needs a
   per-face policy extension. The first is the recommended first cut.
3. **rHEALPix is officially ellipsoidal** (WGS84 authalic sphere). The
   unit-sphere face grid is the right first target — it is what
   `ConservativeRegridding` consumes and what every other grid here emits — with
   the ellipsoidal boundary added later as a pointwise latitude
   reparameterization of the same charts, not as a second geometry.

Tree notes:
- Common cell names are N, O, P, Q, R, S with 3x3 child digits.
- Use rHEALPix projection for exact ellipsoidal cell boundaries.

Sources:
- https://cdnsciencepub.com/doi/10.1139/geomat-2018-0008
- https://proj.org/en/stable/operations/projections/rhealpix.html
- https://github.com/manaakiwhenua/rhealpixdggs-py

Local references:
- global_grid_systems_report.md#111-rhealpix

Notes:
- Prefix partial tree works for ordinal ids: root * 9^level + base9_digits.
"""
struct RHEALPixDGGS <: AbstractDGGS end

system_name(::RHEALPixDGGS) = :rHEALPix
grid_family(::RHEALPixDGGS) = :healpix_rectified
base_solid(::RHEALPixDGGS) = :cube_like
cell_shape(::RHEALPixDGGS) = :square
is_equal_area(::RHEALPixDGGS) = true
aperture(::RHEALPixDGGS) = 9
canonical_index_name(::RHEALPixDGGS) = :root_base9_digits
max_level(::RHEALPixDGGS) = nothing
supports_prefix_ranges(::RHEALPixDGGS) = true
root_count(::RHEALPixDGGS) = 6
radix(::RHEALPixDGGS) = 9
