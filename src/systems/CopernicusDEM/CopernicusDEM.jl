# Copernicus DEM's two-level AWS Open Data COG lattice.

"""
    DiscreteGlobalGrids.CopernicusDEM

The [Copernicus DEM](https://spacedata.copernicus.eu/collections/copernicus-digital-elevation-model)
GLO-30 and GLO-90 lattices in the **AWS Open Data COG** convention: pixel-is-point,
with longitude spacing reduced at latitude bands 50/60/70/80/85.

`CopernicusDEMSystem(30)` selects GLO-30 and `CopernicusDEMSystem(90)` GLO-90.
Level 0 is a global lattice of 1°x1° tiles; level 1 is pixels in north-to-south,
west-to-east raster order. `levelgrid` returns a [`HierarchicalLevelGrid`](@ref).

!!! note "Cannot be wrapped in `AuthalicSystem`"
    Coordinates are geodetic WGS84-G1150 (EPSG:4326), so [`AuthalicSystem`](@ref)
    would reinterpret and warp them. It refuses this system for that reason.

!!! note "DGED, not DTED"
    This implements the AWS DGED profile. It does not implement the DTED
    (`.dt1`/`.dt2`) five-band profile.

!!! note "A holding of any tiles is still this lattice"
    A holding is a `PartialGrid` over the tiles' pixels, in any arrangement, and
    `treeify` gives it a tiled raster tree. System counts, geometry, and nesting
    stay global.

!!! note "`refine` and `coarsen` here are module-local"
    Call `DiscreteGlobalGrids.CopernicusDEM.refine` or `.coarsen`. They relate
    resolutions at one level; `children` and `parent` relate hierarchy levels.
"""
module CopernicusDEM

import ..DiscreteGlobalGrids as DGG
# The tiled raster tree, the leaf container and the rectangle split are the
# engine's; this module answers the lattice hooks they read.
import ..DiscreteGlobalGrids: Engine
# The inline boundary storage IGeo7 already publishes its rings in.
import ..Helpers

import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeometryOps: SpatialTreeInterface as STI
import ConservativeRegridding: Trees
const US = GO.UnitSpherical
# The two-argument constructor belongs to the UnionAll.
const SphericalCap = US.SphericalCap
const Cap = US.SphericalCap{Float64}

const TO_SPHERE = US.UnitSphereFromGeographic()
const FROM_SPHERE = US.GeographicFromUnitSphere()

include("bands.jl")
include("system.jl")
include("nesting.jl")
include("cursor.jl")
include("cursor_memo.jl")

export CopernicusDEMSystem

end # module CopernicusDEM
