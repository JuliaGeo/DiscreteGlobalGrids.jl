# The Copernicus DEM lattice as a two-level grid system: 1-degree tiles over
# variable-longitude-spacing pixel rasters, in the AWS Open Data COG convention.

"""
    DiscreteGlobalGrids.CopernicusDEM

The [Copernicus DEM](https://spacedata.copernicus.eu/collections/copernicus-digital-elevation-model)
GLO-30 and GLO-90 lattices, in the **AWS Open Data COG** convention: pixel-is-point,
3600 (GLO-30) or 1200 (GLO-90) latitude intervals per degree, and a longitude
spacing that steps down by a reduction factor at latitude bands 50/60/70/80/85.

[`CopernicusDEMSystem`](@ref) is the system singleton; `CopernicusDEMSystem(30)` is
GLO-30 and `CopernicusDEMSystem(90)` is GLO-90. Level 0 is one 1°x1° tile (64 800 of
them, covering the globe whether or not AWS ships that tile); level 1 is one pixel.
Within a tile, pixels are in raster order — north row first, west to east — which
makes a tile's pixels one contiguous position window, so `has_sorted_subtrees` is
`true` and `descendant_range` is closed form.

This module ships no grid type: `levelgrid(sys, l)` returns the package's
[`HierarchicalLevelGrid`](@ref).

!!! note "Do not wrap in `AuthalicSystem`"
    Copernicus DEM coordinates are geodetic (WGS84-G1150, EPSG:4326) already.
    [`AuthalicSystem`](@ref) re-reads a grid's geometry as if its latitudes were
    authalic, so wrapping this system would warp coordinates that are not warped.

!!! note "DGED, not DTED"
    The band table below is the DGED profile, which is what both AWS buckets ship.
    The DTED profile (`.dt1`/`.dt2`) has a different table — five bands with factors
    1, 2, 3, 4, 6 — and this system does not implement it.

!!! note "A holding without the far-south tiles is still this lattice"
    Tiles are numbered north to south, so the far-south rows `S90`–`S85` are
    the contiguous tail of the id order (at both levels). A holding that omits
    them is a `PartialGrid` over one contiguous run, and `treeify` still gives
    it the block cursor. Such a subset changes which cells the grid names, not
    what a cell is: `ncells(sys, l)`, [`cell_box`](@ref) and the nesting pair
    keep describing the whole globe.

!!! note "`refine` and `coarsen` here are module-local"
    Reachable only as `DiscreteGlobalGrids.CopernicusDEM.refine` / `.coarsen`,
    so every mention must carry the module. They relate two systems at a fixed
    level (GLO-30 inside GLO-90), not levels within a hierarchy — that relation
    is `children`/`parent`.
"""
module CopernicusDEM

import ..DiscreteGlobalGrids as DGG

import GeometryOps as GO
import GeometryOpsCore as GOCore
import GeometryOps: SpatialTreeInterface as STI
import ConservativeRegridding: Trees
const US = GO.UnitSpherical
# The UnionAll, never `SphericalCap{Float64}`: the two-argument `(point, radius)`
# constructor is a method of the UnionAll.
const SphericalCap = US.SphericalCap
const Cap = US.SphericalCap{Float64}

const TO_SPHERE = US.UnitSphereFromGeographic()
const FROM_SPHERE = US.GeographicFromUnitSphere()

include("bands.jl")
include("system.jl")
include("nesting.jl")
include("cursor.jl")

export CopernicusDEMSystem

end # module CopernicusDEM
