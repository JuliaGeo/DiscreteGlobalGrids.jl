# ---------------------------------------------------------------------------
# S2 grid implementation using closed-form cube-face charts, without a runtime
# dependency on s2geometry.
# ---------------------------------------------------------------------------

"""
    DiscreteGlobalGrids.S2

The [S2](https://s2geometry.io) discrete global grid system: six cube-face
charts, each refined by aperture-4 quadrant subdivision and ordered along a
Hilbert curve, levels `0:30` (`nside = 2^level`, `6 * 4^level` cells).

  - [`S2System`](@ref) — the system singleton.
  - `levelgrid(S2System(), l)` — one complete level, as the package's
    [`HierarchicalLevelGrid`](@ref). S2 ships no grid type of its own.
  - `LevelIndex` — the canonical id: `LevelIndex(level, index)` with `index` the
    **scaffold ordinal** `face * 4^level + hilbert_position`, **0-based**.

S2 is not equal-area; its quadratic `ST → UV` projection gives an approximately
2.08× within-level area spread. Every edge is an exact great-circle arc, so an
S2 cell is exactly its four-corner spherical quadrilateral:

  - [`cell_boundary`](@ref) contains four vertices without densification;
  - generic `cell_area` returns the true cell area;
  - [`node_extent`](@ref) is the cell's four-corner cap because children tile
    their parent.

# Fast paths over the generic fallbacks

| operation | how |
|---|---|
| [`cellat`](@ref) | `point_to_xyf`, the chart's analytic inverse — no tree descent, no point-in-polygon |
| [`cellindex`](@ref) / [`cellposition`](@ref) | the identity, up to the 0-based-id / 1-based-position `± 1` |
| [`descendant_range`](@ref) | `4^Δ`-wide shift of the scaffold ordinal: subtrees are contiguous, hence `has_sorted_subtrees` |
| [`node_extent`](@ref) | the cell's own four-corner cap, uninflated |
| [`neighbors`](@ref) / [`ring`](@ref) | the lattice one-ring plus the cube-edge seam table, under both `Vertex()` (8, or 7 in a face corner) and `Edge()` (4) |

Native 64-bit `s2_cellid` is not an available [`reindex`](@ref) scheme. Its
compatibility is not verified against s2geometry fixtures; the scaffold ordinal
is canonical.
"""
module S2

import ..DiscreteGlobalGrids as DGG

import GeometryOps as GO
const US = GO.UnitSpherical
# Preserve the UnionAll so the two-argument `(point, radius)` constructor applies;
# `SphericalCap{Float64}` would silently reach the three-field default instead.
const SphericalCap = US.SphericalCap

import SmallCollections
using SmallCollections: SmallVector

include("chart.jl")
include("neighbors.jl")
include("system.jl")

export S2System

end # module S2
