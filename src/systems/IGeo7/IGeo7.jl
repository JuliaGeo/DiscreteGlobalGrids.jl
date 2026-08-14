# Load the shared ISEA geometry once, independent of system include order.

isdefined(@__MODULE__, :ISEA) || include("../ISEA/ISEA.jl")

"""
    IGeo7

An aperture-7 hexagonal hierarchy on the icosahedron (ISEA7H) with Z7 indexing.
Use [`IGeo7System`](@ref) for the system and [`Z7Cell`](@ref) for canonical
cell identifiers. Integer lattice arithmetic and Z7 codecs are combined with
the shared [`ISEA`](@ref) geometry; the fitted digit and chirality conventions
are pinned by the sealed DGGRID oracle vectors in `test/systems/IGeo7/vectors/`.
Subtree borders are generated in `O(rim)`; `subtree_border_count` returns their
size in closed form.
"""
module IGeo7

import ..DiscreteGlobalGrids as DGG
# Native geometry names remain unimported; interface methods use `DGG.`.
import ..DiscreteGlobalGrids: AbstractGrid, AbstractHierarchicalGridSystem,
    AbstractCellIndex, Connectivity, Vertex, Edge, level, rawid,
    subtree_border
import ..Helpers
using ..ISEA

import GeometryOps as GO
import SmallCollections
using SmallCollections: SmallVector

const USPoint = GO.UnitSphericalPoint{Float64}

# Dependency order: identifiers, lattice arithmetic, geometry, interface.
include("z7.jl")
include("engine.jl")
include("z7grid.jl")
include("system.jl")

# Public IGeo7 API.
export IGeo7System,
    Z7Cell,
    InvalidZ7Error,
    MAX_RESOLUTION,
    equal_area_steradians,
    is_pentagon,
    is_valid_cell,
    subtree_border,
    subtree_border_count,
    z7_hex,
    z7_string

end # module IGeo7
