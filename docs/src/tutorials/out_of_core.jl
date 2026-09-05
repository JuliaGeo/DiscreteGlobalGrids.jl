# # Out of core: a stencil over a stored cube
#
# This tutorial applies a neighbourhood kernel to a DGGS cube stored on disk.
# Each sequential step loads a chunk and its neighbouring cells (the halo),
# so the computation can process a large field with bounded working memory.

ENV["RASTERDATASOURCES_PATH"] = mkpath(get(ENV, "RASTERDATASOURCES_PATH", joinpath(tempdir(), "rasterdatasources")))

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
using Rasters, RasterDataSources
import NCDatasets
using Zarr
using GLMakie, GeoMakie
using DiscreteGlobalGridsVisualization: dggpoly, dggpoly!
GLMakie.activate!(inline = true)

# ## Get the data onto a store
#
# July soil moisture on IGEO7 level 5, written to a Zarr store with 4096 cells
# per chunk. `dggread` opens the store as a lazy, cell-indexed cube.

soil = Raster(RasterDataSources.getraster(CPCSoil; period = "1981-2010"); name = :soilw)
grid = DGG.levelgrid(DGG.IGeo7System(), 5)
july = DGG.regrid(view(soil, Ti = 7); to = grid, missingval = NaN32)
path = DGG.dggwrite(joinpath(mktempdir(), "soil.zarr"), july; chunks = 4096)
A = DGG.dggread(path)[:soilw]
A = DD.rebuild(A; metadata = DD.NoMetadata())

# `A` is a lazy `DimArray` over a `Cells` dimension. Opening it reads the axis;
# data values arrive as the sweep requests their chunks.
#
# ## Write the kernel
#
# The kernel receives a cell, its value, and the values of its ring. It returns
# the largest difference between the centre and its neighbours; a missing soil
# value remains `NaN`.

function roughness(_, value, values)
    isnan(value) && return NaN32
    biggest = 0.0f0
    for v in values
        isnan(v) || (biggest = max(biggest, abs(value - v)))
    end
    return biggest
end

# ## Apply it chunk by chunk
#
# `mapneighbors!` writes one result per cell into `dest`, here a `DimArray` over
# the same cell axis. This first pass provides an in-memory result for checking
# the computation.

dest = similar(A)
DGG.mapneighbors!(dest, roughness, A)

# Write to a Zarr array to keep the output on disk too. This array stores only
# values; a reusable DGGS store also needs the cell axis and metadata written
# by `dggwrite` (see [DGGS stores](store_io.md)). The equality check below loads
# this small result into memory for comparison:

out = Zarr.zcreate(Float32, length(A); path = joinpath(mktempdir(), "rough.zarr"),
    chunks = (4096,), fill_value = NaN32)
DGG.mapneighbors!(out, roughness, A)
isequal(out[:], parent(dest))

# Roughness is high where the field changes fast: coasts, desert margins, the
# edge of the boreal forest. `dggpoly!` draws the cells that hold a value, so the
# `NaN` ocean stays empty.

fig = Figure(size = (860, 440))
ax = GeoAxis(fig[1, 1]; dest = "+proj=moll", title = "one-ring roughness, IGEO7 level 5")
plt = dggpoly!(ax, dest; color = dest, colormap = :magma, colorrange = (0, 150),
    highclip = :white)
Colorbar(fig[1, 2], plt; label = "mm")
fig

# ## Run the chunks in parallel
#
# `chunkplan` records each owned chunk and the halo width it needs. Set
# `halo = n` when a kernel reaches `n` rings.

plan = DGG.chunkplan(A; halo = 1)

# `split` cuts the plan into pieces over disjoint chunks. Each piece writes a
# disjoint range of `dest`, so the pieces can run as separate tasks. Each
# active task needs memory for its own chunk and halo.

pieces = Base.split(plan, 4)
threaded = similar(A)
@sync for piece in pieces
    Threads.@spawn DGG.mapneighbors!(threaded, roughness, A, piece)
end
isequal(threaded, dest)

# A chunk's halo supplies the neighbours of every owned cell, so the chunked
# result matches the whole-cube result cell for cell. The sweep may read a
# foreign chunk again when it belongs to several halos; the plan controls those
# reads while keeping each callback's working set bounded. [Sweeping a cube
# along its chunk lines](../api/chunk-sweep.md) covers plans, `foreachchunk`,
# and field requests (`needs`).
