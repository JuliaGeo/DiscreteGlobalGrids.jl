# ---------------------------------------------------------------------------
# T5 — H3.
#
# Uber's H3 grid on the redesigned interface. The geometry, the location
# inverse, the hierarchy arithmetic and the adjacency all come from libh3
# itself through `H3_jll`, so what this module contributes is wiring, an order,
# and two things libh3 does not have: a dense position numbering and a subtree
# border walk.
# ---------------------------------------------------------------------------

"""
    DiscreteGlobalGrids.H3

The [H3](https://h3geo.org) discrete global grid system: aperture-7 hexagons on
an icosahedron, twelve pentagons, resolutions `0:15`.

  - [`H3System`](@ref) — the system singleton.
  - [`H3Cell`](@ref) — the canonical id, a `UInt64` with the resolution in
    bits 52-55.
  - [`H3Grid`](@ref) — one complete resolution, from `levelgrid(H3System(), l)`.
  - [`H3Native`](@ref) — the raw libh3 ccall layer, if you want it directly.

# The canonical order

Cells are numbered **base cell by base cell** (the `0:121` order `getRes0Cells`
returns), and **within a base cell by H3's own child position** — the digit path
read as a number, with the digit paths pentagons delete simply absent. That
order is identical to unsigned comparison of the raw index at a fixed
resolution, which is what [`H3Cell`](@ref)'s `isless` is, and it makes every
subtree a contiguous run of positions — hence `has_sorted_subtrees` and the
closed-form [`descendant_range`](@ref).

# Fast paths over the generic fallbacks

| operation | how |
|---|---|
| [`cellat`](@ref) | `latLngToCell`, a closed-form inverse projection |
| [`cellindex`](@ref) / [`cellposition`](@ref) | base-cell prefix sums + `childPosToCell` |
| [`descendant_range`](@ref) | `childPosToCell` + `cellToChildrenSize`, two calls |
| [`neighbors`](@ref) / [`ring`](@ref) | `gridRingUnsafe` shell walks, azimuth-ordered `gridDisk` at pentagons |
| [`ancestor`](@ref) / [`descendants`](@ref) | `cellToParent` / `cellToChildren` across any level gap |
| `subtree_border` | the digit-arc automaton in `border.jl` |

[`node_extent`](@ref) is the generic default — the cell's cap inflated by
[`cap_inflation`](@ref) — which is sound here by measurement; see that
docstring.
"""
module H3

# Exactly the generics this module defines methods on, plus the types those
# methods dispatch on. `node_extent` is deliberately absent: H3 takes the
# generic default (an inflated cell cap), and declaring `cap_inflation` is the
# whole of its say in the matter.
import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, Connectivity, Vertex, Edge,
    ncells, cellindex, cell_boundary, cell_centroid,
    cellposition, rawid,
    cellat, neighbors, ring, system, level,
    cellindextype, levels, levelgrid, rootcells, children,
    cap_inflation, max_neighbors, has_sorted_subtrees,
    ancestor, descendants, descendant_range

import GeometryOps as GO
import SmallCollections
using SmallCollections: SmallVector

# The unit-sphere vocabulary, spelled the same way the fallback substrate
# spells it.
const USPoint = GO.UnitSphericalPoint{Float64}

# The libh3 ccall layer first: everything below is a wiring of it.
include("native.jl")

include("cell.jl")
include("system.jl")
include("geometry.jl")
include("neighbors.jl")
include("border.jl")

end # module H3
