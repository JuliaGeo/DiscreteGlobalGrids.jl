# Copernicus DEM onto IGeo7

`regrid.jl` regrids Copernicus DEM GLO-30 (the default) or GLO-90 onto IGeo7
and writes one Zarr store. It runs the same way on a laptop and on a Slurm
cluster: Dagger spreads the work over whatever Julia workers exist.

```sh
julia --project=scripts/copdem_regrid -e 'using Pkg; Pkg.instantiate()'   # once

julia --project=scripts/copdem_regrid -p 4 -t 2 scripts/copdem_regrid/regrid.jl   # one machine
sbatch scripts/copdem_regrid/regrid.sbatch                                        # Slurm
```

## How it works

1. The land tiles come from the bucket's tile list, or from a scan of the tile
   directory when downloading is off.
2. The IGeo7 cells at `chunklevel` that meet a land tile are the chunks. A chunk
   is one unit of work and one Zarr chunk, so workers write disjoint files.
3. Each worker opens the DEM as one lazy array, `CopernicusUtils.TiledDEM`,
   which reads and caches a tile the first time a regrid touches it.
4. Dagger runs batches of neighbouring chunks on the workers, one per thread.
5. The coordinator appends finished chunks to `<store>.done.txt`. A rerun
   skips them, so a job that hits its time limit is resubmitted unchanged.

## Settings

The settings block at the top of `regrid.jl` holds everything. These
environment variables override the ones a job usually changes:

| Variable | Default | Meaning |
|---|---|---|
| `COPDEM_RES` | `30` | `30` for GLO-30 onto IGeo7 level 13, `90` for GLO-90 onto level 12 |
| `COPDEM_DATA` | `~/copdem` | parent of the defaults below; the tile list is cached here |
| `COPDEM_TILES` | `$COPDEM_DATA/glo<res>` | tile GeoTIFFs, flat or laid out like the AWS bucket |
| `COPDEM_DOWNLOAD` | `1` | `1` fetches missing tiles from AWS; `0` reads `COPDEM_TILES` as complete |
| `COPDEM_STORE` | `$COPDEM_DATA/glo<res>-igeo7-l<level>.zarr` | output store |
| `COPDEM_REGION` | the globe | `west,east,south,north` in whole degrees |
| `COPDEM_SYNTHETIC` | `0` | `1` fabricates every tile from an analytic field and reads no DEM |
| `COPDEM_BATCH` | `8` | chunks per Dagger task |

With `COPDEM_DOWNLOAD=0` a tile with no file on disk reads as ocean, so check
the mirror is complete first: GLO-30 has 26 450 tiles and GLO-90 has 26 475.

## Tiles in your own layout

`COPDEM_TILES` covers tiles named as AWS names them, flat
(`<dir>/<stem>.tif`) or laid out like the bucket (`<dir>/<stem>/<stem>.tif`).
For any other layout, replace `locator` in the "Where the tiles are" section of
`regrid.jl`. It returns a function from a tile's south-west corner, in whole
degrees, to the tile's path, or to `nothing` for a tile that does not exist:

```julia
# /mnt/dem/north/45/Copernicus_DSM_COG_10_N45_00_E006_00_DEM.tif
@everywhere locator(s) = (lat, lon) -> joinpath("/mnt/dem", lat < 0 ? "south" : "north",
    string(abs(lat)), tilestem(s.res, lat, lon) * ".tif")
```

`tilestem(res, lat, lon)` is the AWS name; a layout with its own names builds
the file name from `lat` and `lon` directly. With `COPDEM_DOWNLOAD=1` a missing
tile is fetched to the path the locator gives.

## On a cluster

- Put `COPDEM_TILES` and `COPDEM_STORE` on a filesystem every node sees.
- Mirror the tiles before the job, since compute nodes are often offline:
  `aws s3 sync --no-sign-request s3://copernicus-dem-30m/ $COPDEM_TILES`.
- Start from 4 GB of memory per CPU. A regrid holds up to 1 GiB of weights, and
  each worker caches 32 decoded tiles of 52 MB at 30 m, or 256 of 6 MB at 90 m.
- Try a small synthetic box first. It needs no tiles and takes about a minute
  on two workers: `COPDEM_SYNTHETIC=1 COPDEM_REGION=6,6,45,45`.

`scripts/copdem_production.jl` is the older single-node threaded driver, with a
dependency-graph tile cache and a synthetic source for verification.
