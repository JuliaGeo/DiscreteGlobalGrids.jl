"""
    HEALPix

HEALPix as a DGGS, in two representations over one geometry kernel. Named
`HEALPix` — not `Healpix` — so the registered Healpix.jl package these modules
import is never shadowed.

`chart.jl` is that kernel, and the only Healpix.jl-free layer: the twelve
closed-form face charts `[0, 1]² → S²` of Górski et al. 2005
(DOI:10.1086/427976), plus the two index maps over the `nside × nside` lattice
— RING (1-based data position, any `nside >= 1`) and NESTED (0-based EOPF id, a
Morton code, hence `nside = 2^k` only). All pixel geometry is one chart
evaluation at lattice points, so neighbouring pixels share bit-identical
corners and the tessellation is exact rather than consistent to rounding.

The two representations are the same grid seen two ways: nested ids are one
*ordering* of HEALPix face space, and the DGGS id hierarchy is that ordering's
quadtree. The **id hierarchy** — `HealpixLookups.jl` and `HealpixKernel.jl`,
which put `HEALPixDGGS` on the package's generic kernel and tree family — owns
what only ids can express: sparse coverage of an arbitrary stored subset
(`HealpixLookup`), `subtree_grid` chunking, cross-level parent/child descent.
It needs `nside = 2^k`. The **face grid** — `face_grid.jl`,
[`HealpixFaceSpace`](@ref) and [`HealpixFaceGrid`](@ref) — owns dense
whole-sphere grids at one resolution, where the data vector *is* the grid: any
`nside >= 1`, laid out in RING or NESTED order (or, later, a permuted or
reduced layout) by a swappable [`AbstractHealpixOrdering`](@ref) component that
no grid or tree code has to know about.

Both routes evaluate the same `pixel_corners` / `pixel_center`, so a pixel's
polygon is bitwise identical whichever one produced it — and on both, corner
rings are counter-clockwise as seen from outside the sphere. That is a
contract, not a convention: the convex-clip kernel behind spherical
intersections clips a clockwise ring to EMPTY, so a reversed ring yields silent
zero areas instead of an error. Include order follows the dependency — chart,
face grid, then the Healpix.jl-using layers.
"""
module HEALPix

import ..Helpers

# Pure closed-form face charts + index maps; no Healpix.jl dependency.
include("chart.jl")
# Dense per-resolution grids built on the charts; also Healpix.jl-free.
include("face_grid.jl")
include("HealpixLookups.jl")
# Operations-kernel wiring for `HEALPixDGGS` (see `src/core/kernel.jl`).
include("HealpixKernel.jl")

# The dense face-grid layer's contract surface. System vocabulary stays in the
# submodule — nothing here is re-exported from `DiscreteGlobalGrids` — so reach
# for it as `using DiscreteGlobalGrids.HEALPix`. The chart functions themselves
# stay unexported (`HEALPix.xyf_to_point`, ...): they are the kernel these
# types are built out of, not the API.
export HealpixFaceSpace, HealpixFaceGrid
export AbstractHealpixOrdering, RingOrder, NestedOrder

end # module HEALPix
