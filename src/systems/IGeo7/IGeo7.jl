"""
    IGeo7

An aperture-7 hexagonal hierarchy on the icosahedron (ISEA7H) with Z7 indexing.
Use [`IGeo7System`](@ref) for the system and [`Z7Cell`](@ref) for canonical
cell identifiers. Integer lattice arithmetic and Z7 codecs are combined with
the shared [`ISEA`](@ref) geometry.
Subtree borders are generated in `O(rim)`; `subtree_border_count` returns their
size in closed form.
"""
module IGeo7

import ..DiscreteGlobalGrids as DGG
# Interface methods qualify package-level generics with `DGG.`.
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, Connectivity, Vertex, Edge, level, rawid,
    subtree_border
import ..Helpers
using ..ISEA

import GeometryOps as GO
import SmallCollections
using SmallCollections: SmallVector

const USPoint = GO.UnitSphericalPoint{Float64}

# Dependency order: identifiers, lattice arithmetic, geometry, adjacency, then
# the public system interface.
include("z7.jl")
include("engine.jl")
include("z7grid.jl")
include("gbt.jl")
include("system.jl")
include("indexing.jl")

export IGeo7System,
    Z7Cell,
    RelativeZ7Cell,
    InvalidZ7Error,
    RelativeZ7Error,
    MAX_RESOLUTION,
    equal_area_steradians,
    directioncode,
    is_pentagon,
    is_valid_cell,
    subtree_border,
    subtree_border_count,
    z7_hex,
    z7_string

end # module IGeo7
