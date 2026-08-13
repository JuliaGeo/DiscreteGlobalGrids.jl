# ---------------------------------------------------------------------------
# T6 — HEALPix.
#
# Nested HEALPix on the grid interface. Everything here is a closed form on the
# twelve equal-area face charts; nothing depends on the `Healpix.jl` package.
# ---------------------------------------------------------------------------

"""
    DiscreteGlobalGrids.HEALPix

The [HEALPix](https://healpix.sourceforge.io) discrete global grid system:
twelve equal-area base pixels on the sphere, each refined by aperture-4
quadrant subdivision, levels `0:29` (`nside = 2^level`, `npix = 12 * 4^level`).

  - [`HEALPixSystem`](@ref) — the system singleton.
  - [`HEALPixGrid`](@ref) — one complete level, from
    `levelgrid(HEALPixSystem(), l)`.
  - `LevelIndex` — the canonical id: `LevelIndex(level, index)` with `index` the
    **NESTED** pixel number, **0-based**, as the standard defines it.
  - [`HEALPixRingIndex`](@ref) — the alternate scheme, the **RING** pixel number
    and **1-based** (the convention Healpix.jl and SpeedyWeather's `RingGrids`
    both use). Convert with [`reindex`](@ref) in either direction.

# The layers, bottom up

| file | what it is |
|---|---|
| `chart.jl` | the twelve equal-area face charts and the nested/ring/xyf codecs — pure closed forms, plus the chart's analytic inverse (`point_to_xyf`). Copied wholesale from the pre-redesign `src/HEALPix/chart.jl`. |
| `neighbors.jl` | the 3x3 lattice neighbourhood on the nested id, lifted off the old lookup layer. |
| `system.jl` | `HEALPixSystem`, `HEALPixGrid`, `HEALPixRingIndex`, and every interface method. |
| `border.jl` | the Morton rim walk — a subtree's border with no neighbour queries. |

# Fast paths over the generic fallbacks

| operation | how |
|---|---|
| [`cellat`](@ref) | `point_to_nested`, the chart's analytic inverse — no tree descent, no point-in-polygon |
| [`cellindex`](@ref) / [`cellposition`](@ref) | the identity, up to the 0-based-id / 1-based-position `± 1` |
| [`descendant_range`](@ref) | `4^Δ`-wide shift of the nested id: subtrees are contiguous, hence `has_sorted_subtrees` |
| [`node_extent`](@ref) | the pixel's own bounding cap, uninflated — nested parents *are* the union of their children, so `cap_inflation` is never consulted |
| [`cell_area`](@ref) | `4π / npix` in closed form, **exactly** equal-area, rather than the densified polygon's area |
| [`neighbors`](@ref) / [`ring`](@ref) | the lattice one-ring on the Morton code, under both `Vertex()` (8) and `Edge()` (4) connectivity |
| `subtree_border` | the Morton rim walk in `border.jl` |

# Note the capitalisation

This submodule is `HEALPix`, distinct from the registered `Healpix.jl` package.
Nothing here depends on that package; it is used only by this system's TESTS,
as an independent oracle for the codecs and for point location.
"""
module HEALPix

import ..DiscreteGlobalGrids as DGG

import GeometryOps as GO
const US = GO.UnitSpherical
# The UnionAll, never `SphericalCap{Float64}`: the two-argument
# `(point, radius)` constructor is a method of the UnionAll, and a parametrised
# alias would silently reach the three-field default constructor instead.
const SphericalCap = US.SphericalCap

import SmallCollections
using SmallCollections: SmallVector

include("chart.jl")
include("neighbors.jl")
include("system.jl")
include("border.jl")

export HEALPixSystem, HEALPixGrid, HEALPixRingIndex

end # module HEALPix
