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

!!! note "Dropping the far-south tiles is a lookup, not a lattice"
    The 64 800 tiles are numbered north to south, row-major:
    `tileordinal(r, q) = r*360 + q` with `r = 89 - lat_s`. The six southernmost
    tile rows — `S90` through `S85`, `lat_s in -90:-85` — are therefore the last
    `6*360 = 2 160` ids, the contiguous TAIL of the order, and level 1 inherits
    that because its ids are tile-major in the same order.

    So a holding that does not carry those tiles needs no second lattice and no
    renumbering. It is this lattice minus ONE contiguous window — which a
    `PartialGrid` names directly, and which `treeify` still recognises as a
    rectangle and gives the block cursor:

    ```julia
    import DiscreteGlobalGrids as DGG

    sys = DGG.CopernicusDEMSystem(90)
    south = [DGG.CopernicusDEM.tilecell(sys, lat_s, lon_w)
             for lat_s in -90:-85 for lon_w in -180:179]
    cut = minimum(c.index for c in south)      # the tail runs cut : ncells(sys, 0) - 1
    kept = DGG.PartialGrid(sys, 0, [DGG.LevelIndex(0, i) for i in 0:(cut - 1)])
    DGG.treeify(kept) isa DGG.CopernicusDEM.BlockCursor     # true
    ```

    Keep the two apart: this changes which cells a grid NAMES, not what a cell
    IS. `ncells(sys, l)`, [`cell_box`](@ref), the nesting pair and every id in
    this module go on describing the whole globe.

!!! note "`refine` and `coarsen` here are module-local"
    They are defined in this module and reachable only as
    `DiscreteGlobalGrids.CopernicusDEM.refine` / `.coarsen`: not exported, and
    not methods of any `DiscreteGlobalGrids` generic. An unqualified `refine(…)`
    will not resolve for a reader who has only `using DiscreteGlobalGrids`, so
    every mention of them — here, in docs, at a call site — must carry the
    module.

    Flag for whoever adds a package-level `refine`/`coarsen` later: the names
    collide, and this pair would have to become methods of it or be renamed.
    They relate two SYSTEMS at a fixed level (GLO-30 inside GLO-90), which is
    not the relation a hierarchy `refine` would name — that one is `children`
    and `parent`, which this system already has.
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
