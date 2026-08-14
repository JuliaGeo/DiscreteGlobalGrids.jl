"""
    DiscreteGlobalGrids.RHEALPix

The rHEALPix equal-area discrete global grid: the HEALPix projection with its
four northern and four southern polar triangles assembled into two squares,
followed by a row-major `3 × 3` hierarchy.  The six roots are `N O P Q R S`
and canonical [`RHEALPixCell`](@ref) identifiers use the published prefix SUID
grammar (`N`, `R7`, `S08`, ...).

[`AusPIXSystem`](@ref) is the WGS84 authalic profile of the default `(0, 0)`
rHEALPix layout.  It deliberately reuses the same cell type and identifiers;
AusPIX does not define a second ordinal codec.
"""
module RHEALPix

import ..DiscreteGlobalGrids as DGG
import GeometryOps as GO
import SmallCollections
using SmallCollections: SmallVector

const US = GO.UnitSpherical
const SphericalCap = US.SphericalCap

include("projection.jl")
include("system.jl")
include("neighbors.jl")

export RHEALPixSystem, AusPIXSystem, RHEALPixCell
export healpix_forward, healpix_inverse, rhealpix_forward, rhealpix_inverse
export auspix_forward, auspix_inverse, in_rhealpix_image
export suid, parse_suid

end # module RHEALPix
