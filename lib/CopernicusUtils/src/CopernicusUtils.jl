"""
Copernicus DEM as a data source for DiscreteGlobalGrids: the tile list, real
tiles fetched lazily from the public AWS bucket, a synthetic stand-in with an
analytic oracle, and the whole DEM as one lazy chunked vector.
"""
module CopernicusUtils

import DiscreteGlobalGrids as DGG
import GeometryOps as GO
import GeoInterface as GI
import ArchGDAL
import DiskArrays
import Downloads
import Extents

const CD = DGG.CopernicusDEM
const US = GO.UnitSpherical

export stemtile, tilestem, tilelist, listedtiles
export CopernicusTiles, SyntheticTiles, loadtile, tilepath!, tilecachepath, tileurl
export LandMask, NOMASK, landmask, island, synthetic_elevation, synthetic_tile
export SourceMask, hassource, cellsource
export SubtreeIds, TileIds, tileat, TiledDEM, StripedLRUCache, covering_chunks

include("tiles.jl")
include("sources.jl")
include("synthetic.jl")
include("tileddem.jl")

end # module CopernicusUtils
