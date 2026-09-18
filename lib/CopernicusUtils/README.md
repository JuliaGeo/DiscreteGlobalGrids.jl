# CopernicusUtils

Copernicus DEM (GLO-90 / GLO-30) as a data source for DiscreteGlobalGrids.jl.

| Piece | What it is |
|---|---|
| `tilelist`, `listedtiles` | the bucket's tile list, and the level-0 tile ordinals it names |
| `CopernicusTiles` | real COGs, cached on disk and fetched from AWS on first touch |
| `SyntheticTiles` | an analytic elevation field with an optional land mask, for oracles |
| `TiledDEM` | every listed tile as one lazy `DiskArrays` vector, one chunk per tile |
| `covering_chunks` | the destination cells that meet a set of tiles |

```julia
import DiscreteGlobalGrids as DGG
using CopernicusUtils

sys   = DGG.CopernicusDEMSystem(90)
tiles = listedtiles(sys, tilelist(datadir), [(6.0, 8.0, 45.0, 47.0)])
dem   = TiledDEM(CopernicusTiles(sys, tiles; cachedir = joinpath(datadir, "tiles")), tiles)
space = DGG.DGGSpace(DGG.PartialGrid(sys, 1, dem.ids); chunklevel = 0)
# `dem` and `space` are the `data` and `from` of a lazy `GlobalRegridding.regrid`.
```

`scripts/copdem_production.jl` and `scripts/dagger_regrid/` are the full
GLO-90 -> IGeo7 drivers built on it.
