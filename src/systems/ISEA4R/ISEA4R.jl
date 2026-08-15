# Load the shared ISEA geometry once, independent of system include order.

isdefined(@__MODULE__, :ISEA) || include("../ISEA/ISEA.jl")

"""
    DiscreteGlobalGrids.ISEA4R

An aperture-4 hierarchy over ten equal-area ISEA rhombus charts. Level `l`
uses `nside = 2^l` and contains `10*4^l` cells. Canonical `LevelIndex`
identifiers use `diamond*4^l + morton(ix, iy)`.

The diamond pairing, numbering, and axis orientations are package-specific.
They are not asserted compatible with external ISEA4R or ISEA9R identifiers;
external interoperability requires a fixture-derived permutation.
"""
module ISEA4R

import ..DiscreteGlobalGrids as DGG
using ..ISEA

import GeometryOps as GO
const US = GO.UnitSpherical
# Preserve the UnionAll so the two-argument `(point, radius)` constructor is
# selected.
const SphericalCap = US.SphericalCap

import SmallCollections
using SmallCollections: SmallVector

# Dependency order: layout, chart, topology, interface, border.
include("diamonds.jl")
include("chart.jl")
include("topology.jl")
include("system.jl")
include("border.jl")

export ISEA4RSystem

end # module ISEA4R
