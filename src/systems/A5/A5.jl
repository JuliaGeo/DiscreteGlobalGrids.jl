# A5 interface over the upstream-compatible arithmetic in `native.jl`.

"""
    DiscreteGlobalGrids.A5

The [A5](https://a5geo.org) equal-area pentagonal dodecahedral grid at
resolutions `0:29`. [`A5Cell`](@ref) uses upstream a5's self-describing
`UInt64` encoding; [`A5Native`](@ref) exposes the low-level arithmetic.

A5 is not fixed-radix:

| level | cells | refinement |
|---|---|---|
| 0 | 12 | the dodecahedron's faces |
| 1 | 60 | each face cut into 5 triangular quintants |
| ≥ 2 | `60·4^(l-1)` | 4 Hilbert children per cell |

Cells are ordered by quintant, then Hilbert state `S`:

    index = quintant · 4^(level-1) + S + 1,     quintant = 5·origin + segment

This equals ascending raw-`UInt64` order within a level. Level 0 contains the
twelve faces in origin order.

[`has_sorted_subtrees`](@ref) is `false` because contiguity is not established
across the level-0 fan-out. Location and topology use A5's native arithmetic;
neighbours are wound counter-clockwise. [`node_extent`](@ref) uses the generic
cell cap with inflation `1.75`; the subtree border also uses its fallback.
"""
module A5

# `node_extent` and `descendant_range` use the generic implementations.
import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, Connectivity, Vertex, Edge, HierarchicalLevelGrid,
    ncells, cellindex, cell_boundary, cell_centroid,
    globalindex, rawid,
    cellat, neighbors, ring, one_ring, system, level,
    cellindextype, levels, levelgrid, rootcells, children,
    cap_inflation, maxneighbors, has_sorted_subtrees,
    ancestor, descendants

import GeometryOps as GO
import SmallCollections
using SmallCollections: SmallVector

# Unit-sphere point type used by the public geometry methods.
const USPoint = GO.UnitSphericalPoint{Float64}

# Low-level A5 arithmetic must be loaded before the interface definitions.
include("native.jl")

include("cell.jl")
include("system.jl")
include("geometry.jl")
include("neighbors.jl")

end # module A5
