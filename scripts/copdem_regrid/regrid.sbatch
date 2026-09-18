#!/bin/bash
# Submit from the repository root:  sbatch scripts/copdem_regrid/regrid.sbatch
#
# Each Slurm task becomes one Julia worker with --cpus-per-task threads, and
# every thread regrids one chunk at a time. Several small workers per node beat
# one large one: GeoTIFF decoding is serial within a process.

#SBATCH --job-name=copdem-igeo7
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=4G
#SBATCH --time=24:00:00
#SBATCH --output=copdem-%j.log

export COPDEM_RES=90
export COPDEM_TILES=/path/to/copdem/glo90             # shared filesystem
export COPDEM_STORE=/path/to/output/glo90-igeo7.zarr  # shared filesystem
export COPDEM_DOWNLOAD=0                              # compute nodes read the tiles, never fetch them

# The coordinator only schedules, so it runs here in the batch step; the
# script starts the workers itself with `srun`.
julia --project=scripts/copdem_regrid scripts/copdem_regrid/regrid.jl
