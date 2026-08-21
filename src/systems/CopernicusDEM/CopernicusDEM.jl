# Copernicus DEM's two-level AWS Open Data COG lattice.

"""
    DiscreteGlobalGrids.CopernicusDEM

The [Copernicus DEM](https://spacedata.copernicus.eu/collections/copernicus-digital-elevation-model)
GLO-30 and GLO-90 lattices in the **AWS Open Data COG** convention: pixel-is-point,
with longitude spacing reduced at latitude bands 50/60/70/80/85.

`CopernicusDEMSystem(30)` selects GLO-30 and `CopernicusDEMSystem(90)` GLO-90.
Level 0 is a global lattice of 1°x1° tiles; level 1 is pixels in north-to-south,
west-to-east raster order. `levelgrid` returns a [`HierarchicalLevelGrid`](@ref).

!!! note "Do not wrap in `AuthalicSystem`"
    Coordinates are geodetic WGS84-G1150 (EPSG:4326). [`AuthalicSystem`](@ref)
    would reinterpret and warp them.

!!! note "DGED, not DTED"
    This implements the AWS DGED profile. It does not implement the DTED
    (`.dt1`/`.dt2`) five-band profile.

!!! note "A holding without the far-south tiles is still this lattice"
    Rows `S90`–`S85` are the contiguous tail of both level orders, so a holding
    without them is a `PartialGrid` over one id run, and `treeify` still gives it
    the block cursor. System counts, geometry, and nesting stay global.

!!! note "`refine` and `coarsen` here are module-local"
    Call `DiscreteGlobalGrids.CopernicusDEM.refine` or `.coarsen`. They relate
    resolutions at one level; `children` and `parent` relate hierarchy levels.
"""
module CopernicusDEM

import ..DiscreteGlobalGrids as DGG
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
