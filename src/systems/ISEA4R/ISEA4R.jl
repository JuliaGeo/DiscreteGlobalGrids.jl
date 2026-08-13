# ---------------------------------------------------------------------------
# T12 — ISEA4R.
#
# The shared ISEA basis (icosahedron tables + Snyder charts) is included here
# behind the same guard `IGeo7.jl` uses: it belongs to the ISEA *family*, not to
# any one system, so whichever system is included first defines it and the
# second include is a no-op. The include order of `src/systems/` never matters.
# ---------------------------------------------------------------------------

isdefined(@__MODULE__, :ISEA) || include("../ISEA/ISEA.jl")

"""
    DiscreteGlobalGrids.ISEA4R

The ISEA4R discrete global grid system: the icosahedron's twenty faces paired
into **ten rhombi** ("diamonds"), each carrying a square chart `[0, 1]² → S²`
that is refined by aperture-4 quadrant subdivision. Levels `0:29`
(`nside = 2^level`, `10 * 4^level` cells).

  - [`ISEA4RSystem`](@ref) — the system singleton.
  - `levelgrid(ISEA4RSystem(), l)` — one complete level, as the package's
    [`HierarchicalLevelGrid`](@ref). ISEA4R ships no grid type of its own.
  - `LevelIndex` — the canonical id: `LevelIndex(level, index)` with
    `index = diamond * 4^level + morton(ix, iy)`, **0-based** and dense, so a
    cell's position in the level grid is `index + 1`.

# The layers, bottom up

| file | what it is |
|---|---|
| `diamonds.jl` | the ten-diamond layout table — which two faces pair, how the ten are numbered, how the `(x, y)` square sits in each. Derived at load time from `ISEA`'s tables and asserted against pinned literals. Ported from the pre-redesign `src/ISEA4R/diamonds.jl`. |
| `chart.jl` | the piecewise-affine, exactly equal-area rhombus chart and the row-major / Morton codecs, plus the chart's analytic inverse `point_to_xyd`. Ported from `src/ISEA4R/chart.jl`; the inverse is new. |
| `topology.jl` | the rhombus lattice one-ring and the ten-diamond seam topology — the edge-adjacency and vertex-fan tables, derived from the layout. New. |
| `system.jl` | `ISEA4RSystem` and every interface method, system-level and grid-level both. New. |
| `border.jl` | the lattice-block rim walk — a subtree's border with no neighbour queries. New. |

# Fast paths over the generic fallbacks

| operation | how |
|---|---|
| [`cellat`](@ref) | `point_to_morton`, the chart's analytic inverse — no tree descent, no point-in-polygon |
| [`cellindex`](@ref) / [`cellposition`](@ref) | the identity, up to the 0-based-id / 1-based-position `± 1` |
| [`descendant_range`](@ref) | `4^Δ`-wide shift of the Morton id: subtrees are contiguous, hence `has_sorted_subtrees` |
| [`node_extent`](@ref) | the cell's own bounding cap, uninflated — children tile the parent rhombus in chart space, so the tight cap is already legal |
| [`cell_area`](@ref) | `4π / (10 * 4^level)` in closed form, **exactly** equal-area |
| [`neighbors`](@ref) / [`ring`](@ref) | the lattice one-ring plus the seam tables, under both `Vertex()` (9) and `Edge()` (4) connectivity |
| [`subtree_border`](@ref) | the boundary ring of the subtree's lattice block, `O(rim)` rather than `O(subtree)` |

# Provenance — read this before claiming interoperability

The projection is entirely `ISEA`'s Snyder machinery, unchanged. What this
module adds is the **layout**: the diamond pairing, the numbering of the ten,
and the orientation of the `(x, y)` square inside each. That layout is this
package's own canonical choice, derived from `ISEA.FACE_TRIPLES` /
`ISEA.NBRS_CCW` / `ISEA.NEIGHBORS` in the standard ISEA placement and anchored
on the vertex pair `(0, 11)`. **There is no external oracle pinning it:
identifier compatibility with any external ISEA4R/ISEA9R product, DGGAL
included, is deliberately not claimed and must not be inferred.** The layout
coincides in shape with SphericalSpatialTrees.jl's `ISEACircleTree`
(`10 × 2^r × 2^r`), but neither the diamond numbering nor the per-diamond axis
orientation has been cross-pinned against it. Anyone needing external
identifier interop must first fit a permutation against fixtures.

ISEA9R (the aperture-9 sibling over the same charts) is not part of this
module: it would want a radix-9 id space over the same layout, and nothing here
is aperture-specific except `system.jl`.
"""
module ISEA4R

import ..DiscreteGlobalGrids as DGG
using ..ISEA

import GeometryOps as GO
const US = GO.UnitSpherical
# The UnionAll, never `SphericalCap{Float64}`: the two-argument
# `(point, radius)` constructor is a method of the UnionAll, and a parametrised
# alias would silently reach the three-field default constructor instead.
const SphericalCap = US.SphericalCap

import SmallCollections
using SmallCollections: SmallVector

# Dependency order: the layout table, the chart built on it, the topology
# derived from the layout, then the interface wiring over all three.
include("diamonds.jl")
include("chart.jl")
include("topology.jl")
include("system.jl")
include("border.jl")

# The system's contract surface. The package-level exports are T13's; these
# make the names reachable as `DiscreteGlobalGrids.ISEA4R.<name>`.
export ISEA4RSystem

end # module ISEA4R
