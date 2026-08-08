"""
    ISEA4R

ISEA4R diamond chart grids: the ten continuous charts `[0, 1]² → S²` over the
icosahedron's ten rhombi, and dense whole-sphere grids built on them, under a
swappable data ordering.

The projection is not new here — it is [`ISEA`](@ref)'s Snyder equal-area
machinery, used unchanged. What this module adds is the *layout*: which two
icosahedron faces form each diamond, how the ten are numbered, and how the
`(x, y)` square is oriented inside each. `diamonds.jl` derives that table at
load time from `ISEA.FACE_TRIPLES` / `ISEA.NBRS_CCW` / `ISEA.NEIGHBORS` and
asserts it against pinned literals; `chart.jl` assembles the two face triangles
into one piecewise-affine, exactly equal-area square chart; `face_grid.jl` turns
the charts into a `SpatialTreeInterface` tree a
`ConservativeRegridding.Regridder` consumes directly.

# Why not inside `ISEA`

`ISEA` is charter-bound to aperture-agnostic, oracle-locked machinery: every
constant and convention in it is pinned against recorded oracle output or a
published spec. The ten-diamond layout is deliberately *not* oracle-locked — it is a
package convention with no external fixture behind it (see below) — so it lives
in a consumer module that does `using ..ISEA`, exactly as [`IGeo7`](@ref) does.

# Scope, honestly

**The dense diamond grid, plus chart geometry for `ISEA4RDGGS`.**
`face_grid.jl` — [`Isea4rFaceSpace`](@ref) and [`Isea4rFaceGrid`](@ref) — owns
dense whole-sphere grids at one resolution, where the data vector *is* the
grid: any `nside >= 1`, laid out in row-major or Morton order by a swappable
[`AbstractIsea4rOrdering`](@ref) component that no grid or tree code has to know
about. `Isea4rKernel.jl` then answers the package's geometry generics for
`ISEA4RDGGS` over the canonical `isea4r_ordinal`
`diamond * 4^level + morton_position`, so `cell_polygon(ISEA4RDGGS(), level,
id)` returns a real polygon — bitwise the one the diamond grid emits for the
same ordinal.

`num_cells` and `ordinal_to_cell` enumerate those geometry ids directly. The
`isea4r_ordinal` id *hierarchy* is still deferred: `cell_children`,
`cell_parent`, `descendant_range` and the rest of that group throw
`NotPortedError`. Nothing blocks it — the radix-4 arithmetic over these
ordinals is exact — it is simply out of this milestone's geometry scope. What
[`MortonOrder`](@ref) does assert today is that it is a bijection and that it
realizes the registry's ordinal `diamond * 4^level + position` with
`position := morton(ix, iy)`.

# Provenance, and what is NOT claimed

This numbering, pairing, and in-diamond orientation is this package's canonical
choice, derived at build time from `ISEA`'s tables in the standard ISEA
placement (identity `Orientation`), anchored on the vertex pair `(0, 11)`.
**There is no external oracle pinning it: identifier compatibility with any
external ISEA4R/ISEA9R product, DGGAL included, is deliberately not claimed
and must not be inferred.** The layout coincides in shape with
SphericalSpatialTrees.jl's `ISEACircleTree` (10 × 2^r × 2^r), but the diamond
numbering and per-diamond axis orientation have not been cross-pinned against
SST either. Anyone needing DGGAL/SST or any other external identifier interop
must first pin a permutation against fixtures. The full record is
`docs/design/isea4r_diamond_layout.md`.

Corner rings are counterclockwise as seen from outside the sphere. That is a
contract, not a convention: the convex-clip kernel behind spherical
intersections clips a clockwise ring to EMPTY, so a reversed ring yields silent
zero areas instead of an error. Include order follows the dependency — layout
table, chart, face grid, then the kernel wiring built on them.
"""
module ISEA4R

using ..ISEA

# The ten-diamond layout table: derived, snapped, asserted against pinned
# literals at load time.
include("diamonds.jl")
# The piecewise-affine rhombus chart + the row-major / Morton index maps.
include("chart.jl")
# Dense per-resolution diamond grids built on the charts.
include("face_grid.jl")
# Geometry-only operations-kernel wiring for `ISEA4RDGGS` (see
# `src/core/kernel.jl`); last, because it evaluates the chart and documents
# itself against the grid.
include("Isea4rKernel.jl")

# The dense face-grid layer's contract surface. System vocabulary stays in the
# submodule — nothing here is re-exported from `DiscreteGlobalGrids` — so reach
# for it as `using DiscreteGlobalGrids.ISEA4R`. The chart functions themselves
# stay unexported (`ISEA4R.xyd_to_point`, ...): they are the kernel these types
# are built out of, not the API. `cell_polygon`, `cell_center` and `num_cells`
# stay unexported for the additional reason that the package namespace already
# has those names.
export Isea4rFaceSpace, Isea4rFaceGrid
export AbstractIsea4rOrdering, RowMajorOrder, MortonOrder

end # module ISEA4R
