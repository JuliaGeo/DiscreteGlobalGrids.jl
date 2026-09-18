# CopernicusUtils

Copernicus DEM (GLO-90 / GLO-30) as a data source for DiscreteGlobalGrids.jl.

| Piece | What it is |
|---|---|
| `tilelist`, `landtiles` | the bucket's tile list, and the level-0 ordinals of the tiles it names; every other tile is ocean |
| `CopernicusTiles` | real COGs, cached on disk and fetched from AWS on first touch |
| `TileDirectory` | the default tile locator: AWS file names under a prefix, flat or as a bucket mirror |
| `SyntheticTiles` | an analytic elevation field with an optional land mask, for oracles |
| `TiledDEM` | every land tile as one lazy `DiskArrays` vector, one chunk per tile |
| `covering_chunks` | the destination cells that meet a set of tiles |

```julia
import DiscreteGlobalGrids as DGG
using CopernicusUtils

sys   = DGG.CopernicusDEMSystem(90)
tiles = landtiles(sys, tilelist(datadir), [(6.0, 8.0, 45.0, 47.0)])
dem   = TiledDEM(CopernicusTiles(sys, tiles; cachedir = joinpath(datadir, "tiles")), tiles)
space = DGG.DGGSpace(DGG.PartialGrid(sys, 1, dem.ids); chunklevel = 0)
# `dem` and `space` are the `data` and `from` of a lazy `GlobalRegridding.regrid`.
```

## Tiles you already have

`CopernicusTiles` finds files through a locator, `locate(lat, lon)`, which takes
a tile's whole-degree south-west corner and returns its path, or `nothing` for a
tile that does not exist. `cachedir = dir` is the locator `TileDirectory(dir, 90)`.
Any callable serves, and `landtiles(sys, locate)` returns the tiles it finds on
disk:

```julia
locate(lat, lon) = joinpath("/mnt/dem", lat < 0 ? "south" : "north", tilestem(90, lat, lon) * ".tif")
tiles = landtiles(sys, locate)
dem   = TiledDEM(CopernicusTiles(sys, tiles; locate, download = false), tiles)
```

`scripts/copdem_production.jl` and `scripts/dagger_regrid/` are the full
GLO-90 -> IGeo7 drivers built on it.
