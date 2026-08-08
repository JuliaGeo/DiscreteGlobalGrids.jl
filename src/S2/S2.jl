"""
    S2

S2 cube-face chart grids: the six closed-form charts `[0, 1]² → S²` and dense
whole-sphere grids built on them, under a swappable data ordering.

`chart.jl` is the kernel, and it is s2geometry-*convention* code with zero
external dependency — the face frames and the quadratic `ST → UV` transform of
`s2coords.h`, and the Hilbert tables of `s2cell_id.h`, transcribed from the
published conventions rather than vendored or linked. Everything else is one
evaluation of `stf_to_point` on the `(s, t)` lattice, so neighbouring cells
share bit-identical corners and the tessellation is exact rather than consistent
to rounding. Because there is no S2 oracle in this repository, the tests here
are correspondingly oracle-free: analytic invariants and internal consistency,
not a reference implementation (contrast `test/HEALPix/test_chart.jl`, which has
Healpix.jl to check against).

# Scope, honestly

**The dense face grid, plus chart geometry for `S2DGGS`.** `face_grid.jl` —
[`S2FaceSpace`](@ref) and [`S2FaceGrid`](@ref) — owns dense whole-sphere grids
at one resolution, where the data vector *is* the grid: any `nside >= 1`, laid
out in row-major or Hilbert order by a swappable
[`AbstractS2Ordering`](@ref) component that no grid or tree code has to know
about. `S2Kernel.jl` then answers the package's geometry generics for
`S2DGGS` over the *scaffold ordinal* `face * 4^level + hilbert_position`, so
`cell_polygon(S2DGGS(), level, id)` returns a real polygon — bitwise the one
the face grid emits for the same ordinal.

`num_cells` and `ordinal_to_cell` enumerate those geometry ids directly. The
`s2_cellid` id hierarchy is still a later milestone: `cell_children`,
`cell_parent`, `descendant_range` and the rest of that group throw
`NotPortedError`, because the canonical id is the native 64-bit `s2_cellid` and
declaring the hierarchy over scaffold ordinals would answer in the wrong
coordinate system. Native-id alignment claims wait for that port too. What the
Hilbert ordering does assert today is that it is a bijection, local, and nests
across resolutions — which is what makes it the S2 analogue of HEALPix's
`NestedOrder`.

# Recorded decisions

- **The quadratic transform.** `st_to_uv` is s2geometry's
  `S2_QUADRATIC_PROJECTION`, not the simpler linear `u = 2s - 1`. Only the
  quadratic map puts the `(s, t)` lattice at `nside = 2^level` on exact
  canonical S2 cell boundaries, which is the property the eventual native-id
  wiring will need; it also narrows the within-level cell-area spread from ~5.2× to
  ~2.08×. It does not make the grid equal-area — S2 is
  `is_equal_area == false`, and each face carries exactly `4π/6` for *any*
  monotone transform anyway, by cube symmetry.
- **No `MortonOrder`.** S2 has no Morton convention to be compatible with,
  `RowMajorOrder` already covers every `nside`, and the ordering contract makes
  a Morton ordering about fifteen lines of user code.

Corner rings are counter-clockwise as seen from outside the sphere. That is a
contract, not a convention: the convex-clip kernel behind spherical
intersections clips a clockwise ring to EMPTY, so a reversed ring yields silent
zero areas instead of an error. Include order follows the dependency — chart,
face grid, then the kernel wiring built on both.
"""
module S2

# Pure closed-form face charts + index maps; no s2geometry dependency.
include("chart.jl")
# Dense per-resolution grids built on the charts.
include("face_grid.jl")
# Geometry-only operations-kernel wiring for `S2DGGS` (see `src/core/kernel.jl`);
# last, because it evaluates the chart and documents itself against the grid.
include("S2Kernel.jl")

# The dense face-grid layer's contract surface. System vocabulary stays in the
# submodule — nothing here is re-exported from `DiscreteGlobalGrids` — so reach
# for it as `using DiscreteGlobalGrids.S2`. The chart functions themselves stay
# unexported (`S2.stf_to_point`, ...): they are the kernel these types are built
# out of, not the API. `cell_polygon`, `cell_center` and `num_cells` stay
# unexported for the additional reason that the package namespace already has
# those names.
export S2FaceSpace, S2FaceGrid
export AbstractS2Ordering, RowMajorOrder, HilbertOrder

end # module S2
