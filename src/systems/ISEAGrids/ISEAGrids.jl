# The three central-place/triangular ISEA grids share the package's Snyder
# projection.  This module deliberately owns only the lattice, hierarchy and
# index layers; the projection remains in ../ISEA.
isdefined(@__MODULE__, :ISEA) || include("../ISEA/ISEA.jl")

module ISEAGrids

import ..DiscreteGlobalGrids as DGG
using ..ISEA

import GeometryOps as GO
const US = GO.UnitSpherical
const USPoint = GO.UnitSphericalPoint
const LevelGrid{S} = DGG.HierarchicalLevelGrid{S}

include("z3.jl")
include("hex.jl")
include("triangles.jl")

export ISEA3HSystem, ISEA4HSystem, ISEA4TSystem, Z3Cell,
    z3_string, is_pentagon, equal_area_steradians

end # module ISEAGrids
