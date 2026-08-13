# ---------------------------------------------------------------------------
# T11 — S2.
#
# S2 on the grid interface. Everything here is a closed form on the six
# cube-face charts; nothing depends on the s2geometry library, and no
# s2geometry code is vendored.
# ---------------------------------------------------------------------------

"""
    DiscreteGlobalGrids.S2

The [S2](https://s2geometry.io) discrete global grid system: six cube-face
charts, each refined by aperture-4 quadrant subdivision and ordered along a
Hilbert curve, levels `0:30` (`nside = 2^level`, `6 * 4^level` cells).

  - [`S2System`](@ref) — the system singleton.
  - [`S2Grid`](@ref) — one complete level, from `levelgrid(S2System(), l)`.
  - `LevelIndex` — the canonical id: `LevelIndex(level, index)` with `index` the
    **scaffold ordinal** `face * 4^level + hilbert_position`, **0-based**.

# The layers, bottom up

| file | what it is |
|---|---|
| `chart.jl` | the six cube-face charts, the quadratic `ST ↔ UV` transform, the Hilbert and row-major codecs, and the chart's analytic inverse (`point_to_xyf`). Ported from the pre-redesign `src/S2/chart.jl`, plus the inverse the old port never wired. |
| `neighbors.jl` | the 3×3 lattice neighbourhood, and the cube-edge seam table it crosses faces through — derived from the face frames, not transcribed. |
| `system.jl` | `S2System`, `S2Grid`, and every interface method. |

# What S2 trades, and what it buys

S2 is **not equal-area** — the quadratic `ST → UV` projection narrows the
within-level cell-area spread to about 2.08×, and that is all it does about it.
What it buys is that every cell edge is an exact great-circle arc, so an S2
cell is a spherical quadrilateral *exactly*:

  - [`cell_boundary`](@ref) is four vertices, not a densification (HEALPix
    needs 32);
  - `cell_area` needs no override — the generic spherical polygon area of that
    ring is the true cell area;
  - [`node_extent`](@ref) is the cell's own four-corner cap, exactly, because
    the cell is the geodesic convex hull of those corners and children tile
    their parent.

# Fast paths over the generic fallbacks

| operation | how |
|---|---|
| [`cellat`](@ref) | `point_to_xyf`, the chart's analytic inverse — no tree descent, no point-in-polygon |
| [`cellindex`](@ref) / [`cellposition`](@ref) | the identity, up to the 0-based-id / 1-based-position `± 1` |
| [`descendant_range`](@ref) | `4^Δ`-wide shift of the scaffold ordinal: subtrees are contiguous, hence `has_sorted_subtrees` |
| [`node_extent`](@ref) | the cell's own four-corner cap, uninflated |
| [`neighbors`](@ref) / [`ring`](@ref) | the lattice one-ring plus the cube-edge seam table, under both `Vertex()` (8, or 7 in a face corner) and `Edge()` (4) |

# Deferred

The native 64-bit `s2_cellid` (face bits, Hilbert bits, lsb sentinel) as an
alternate scheme via [`reindex`](@ref). The codec is one step from
`xyf_to_hilbert`, but this repository carries no s2geometry fixtures, so
shipping it would publish an interoperability claim nothing checks. The
scaffold ordinal is canonical either way, so adding the scheme later is purely
additive.
"""
module S2

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

export S2System, S2Grid

end # module S2
