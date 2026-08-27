# H3 interface over libh3, with dense grid indices and subtree-border walks.

"""
    DiscreteGlobalGrids.H3

The [H3](https://h3geo.org) aperture-7 icosahedral grid: hexagons, twelve
pentagons, and resolutions `0:15`. [`H3Cell`](@ref) wraps libh3's `UInt64` id;
[`H3Native`](@ref) exposes the low-level calls.

Canonical order is base-cell-major, then H3 child position with deleted
pentagon paths omitted. It matches raw-id order within a resolution and makes
subtrees contiguous, enabling exact [`descendant_range`](@ref) values. Geometry,
location, hierarchy, and adjacency use libh3; `border.jl` implements the
digit-arc subtree border. [`node_extent`](@ref) uses the generic inflated cap.
"""
module H3

# `node_extent` uses the generic inflated cell cap.
import ..DiscreteGlobalGrids as DGG
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, Connectivity, Vertex, Edge, HierarchicalLevelGrid,
    ncells, cellindex, cell_boundary, cell_centroid,
    globalindex, rawid,
    cellat, neighbors, ring, one_ring, system, level,
    cellindextype, levels, levelgrid, rootcells, children,
    cap_inflation, maxneighbors, has_sorted_subtrees, has_direct_location,
    ancestor, descendants, descendant_range

import GeometryOps as GO
import SmallCollections
using SmallCollections: SmallVector

# Unit-sphere point type used by the public geometry methods.
const USPoint = GO.UnitSphericalPoint{Float64}

# Load the libh3 call layer before the interface definitions.
include("native.jl")

include("cell.jl")
include("system.jl")
include("geometry.jl")
include("neighbors.jl")
include("border.jl")

end # module H3
