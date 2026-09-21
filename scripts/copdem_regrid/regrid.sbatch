#!/bin/bash
# Submit from the repository root:  sbatch scripts/copdem_regrid/regrid.sbatch
#
# Each Slurm task becomes one Julia worker with --cpus-per-task threads, and
# every thread regrids one chunk at a time. Several small workers per node beat
# one large one: GeoTIFF decoding is serial within a process.

#SBATCH --job-name=copdem-igeo7
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=4
#SBATCH --mem=1000G
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --partition=big
#SBATCH --output=copdem-%j.log

export COPDEM_RES=30
export COPDEM_TILES=/Net/Groups/BGI/data/DataStructureMDI/DATA/Incoming/Copernicus/cop30dem/tiles             # shared filesystem; for another layout edit `locator` in regrid.jl
export COPDEM_STORE=/Net/Groups/BGI/work_5/scratch/fgans/glo30-igeo7.zarr # shared filesystem
export COPDEM_DOWNLOAD=0                              # compute nodes read the tiles, never fetch them
export JULIA_WORKER_TIMEOUT=120
# The coordinator only schedules, so it runs here in the batch step; the
# script starts the workers itself with `srun`.
julia --project=scripts/copdem_regrid scripts/copdem_regrid/regrid.jl
