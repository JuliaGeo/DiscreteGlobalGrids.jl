# # Out of core: a stencil over a stored cube
#
# `mapneighbors!` runs a neighbourhood kernel over a cube stored on disk, one
# chunk at a time. Each chunk is read once together with the ring of cells
# around it, the kernel runs on that block, and the block's results go to the
# destination. Memory holds one chunk and its halo, whatever the store's size.

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
# to a chunk. `dggread` opens any DGGS store the same way.

soil = Raster(RasterDataSources.getraster(CPCSoil; period = "1981-2010"); name = :soilw)
grid = DGG.levelgrid(DGG.IGeo7System(), 5)
july = DGG.regrid(view(soil, Ti = 7); to = grid, missingval = NaN32)
path = DGG.dggwrite(joinpath(mktempdir(), "soil.zarr"), july; chunks = 4096)
A = DGG.dggread(path)[:soilw]
A = DD.rebuild(A; metadata = DD.NoMetadata())   # the store's header is noise here

# `A` is a lazy `DimArray` over a `Cells` dimension: opening it read only the
# cell axis.
#
# ## Write the kernel
#
# The kernel receives a cell, its value, and the values of its ring. Roughness
# is the largest difference between a cell and any member of its ring. Soil
# moisture is a land field, so an ocean cell stays `NaN`.

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
# `mapneighbors!` writes one result per cell into `dest`. `dest` is anything
# indexable along the cell axis — here a `DimArray` over the same cells.

dest = similar(A)
DGG.mapneighbors!(dest, roughness, A)

# A Zarr array as `dest` streams the results to disk, so the whole pass runs
# with one chunk and its halo in memory:

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
# `chunkplan` lists the chunks a sweep visits and the halo each one reads.
# `halo = n` carries `n` rings of context, for a kernel that reaches `n` rings.

plan = DGG.chunkplan(A; halo = 1)

# `split` cuts the plan into pieces over disjoint chunks. The pieces write
# disjoint ranges of `dest`, so they run on separate tasks with no coordination.

pieces = Base.split(plan, 4)
threaded = similar(A)
@sync for piece in pieces
    Threads.@spawn DGG.mapneighbors!(threaded, roughness, A, piece)
end
isequal(threaded, dest)

# A chunk's halo holds every neighbour of every cell it owns, so the chunked
# result matches the whole-cube result cell for cell. [Sweeping a cube along its
# chunk lines](../api/chunk-sweep.md) covers the plan, `foreachchunk` and field
# requests (`needs`).
