"""
    DiscreteGlobalGrids.IVEARTEA

Slice-and-dice equal-area icosahedral grids in a 5-by-6 atlas derived from the
published DGGAL tables.
IVEA uses an icosahedron vertex as the radial vertex of each of the 120
fundamental triangles; RTEA uses the adjacent edge midpoint.  The projection
kernel follows van Leeuwen and Strebe (2006).  Atlas and interoperability
tables are derived from BSD-3-Clause DGGAL commit `e16cea7`; only level-zero
geometry has currently been checked against the sealed DGGAL corpus.

The concrete systems currently implemented are the rhombic aperture-4 and
aperture-9 grids.  Their canonical `LevelIndex` is a dense, root-major,
row-major lattice index. Deep DGGAL/ZIRS identifier and exact seam-tie
compatibility remain unclaimed until deeper post-fix oracle vectors exist.
"""
module IVEARTEA

import ..DiscreteGlobalGrids as DGG
import GeometryOps as GO
const US = GO.UnitSpherical
const SphericalCap = US.SphericalCap

abstract type AbstractSliceDiceRhombicSystem <: DGG.AbstractHierarchicalGridSystem end
"""`IVEA4RSystem`: IVEA projection with aperture-4 rhombic refinement (10 roots)."""
struct IVEA4RSystem <: AbstractSliceDiceRhombicSystem end
"""`IVEA9RSystem`: IVEA projection with aperture-9 rhombic refinement (10 roots)."""
struct IVEA9RSystem <: AbstractSliceDiceRhombicSystem end
"""`RTEA4RSystem`: RT(S)EA projection with aperture-4 rhombic refinement (10 roots)."""
struct RTEA4RSystem <: AbstractSliceDiceRhombicSystem end
"""`RTEA9RSystem`: RT(S)EA projection with aperture-9 rhombic refinement (10 roots)."""
struct RTEA9RSystem <: AbstractSliceDiceRhombicSystem end

include("projection.jl")
include("atlas.jl")
include("topology.jl")
include("system.jl")

export IVEA4RSystem, IVEA9RSystem, RTEA4RSystem, RTEA9RSystem,
    IVEAProfile, RTEAProfile, forward_5x6, inverse_5x6

end # module IVEARTEA
