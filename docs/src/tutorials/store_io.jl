# # A round trip through a DGGS store
#
# This tutorial shows how to persist a cell-indexed cube, reopen it lazily, and
# select only the stored cells needed for a region. It also compares the two
# cell-id encodings available to a DGGS store.
#
# `dggwrite` and `dggread` are provided by the Zarr.jl extension, loaded by
# `using Zarr`.

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
using Zarr
using DiscreteGlobalGridsVisualization: dggpoly, dggpoly!
using GLMakie, GeoMakie
GLMakie.activate!(inline = true)

# ## A cube to write
#
# The example uses two level-1 IGEO7 cells and their level-4 descendants as a
# small regional store. `CellVector` names the cells, `CellLookup` turns them
# into a one-level axis, and `Cells` makes that axis a cube dimension.

sys = DGG.IGeo7System()
roots = DGG.CellVector(DGG.levelgrid(sys, 1))[1:2]
cells = sort!(reduce(vcat, [collect(DGG.descendants(sys, c, 4)) for c in roots]))
lookup = DGG.CellLookup(DGG.CellVector(sys, 4, cells))

# The two layers share the cell axis. Their values encode axis positions, which
# makes the later selection checks easy to read.

n = length(cells)
elevation = Float32.(1:n)
slope = Float32.(0.5 .* (1:n))
cube = DD.DimStack((; elevation, slope), (DGG.Cells(lookup),))

# ## Write the cube and read it back
#
# `dggwrite` returns the path that `dggread` opens. `chunks = 128` gives this
# small cube five chunks; `:auto` would place the whole example in one chunk.

path = DGG.dggwrite(joinpath(mktempdir(), "demo.zarr"), cube; chunks = 128)
store = DGG.dggread(path)
#
DD.metadata(store)["description"]

# `dggread` reconstructs the grid description from store attributes. The
# returned arrays stay lazy until a value is requested; `collect` makes an
# explicit in-memory comparison with the original cube:

axis = DD.lookup(store[:elevation], DGG.Cells)
#
collect(axis) == cells, collect(parent(store[:elevation])) == elevation

# ## Selecting a region out of a store
#
# The axis is a `ChunkedCellLookup`. It supports the same selectors as a
# `CellLookup` and uses the chunk manifest to locate the required id data. Three
# selectors name a single cell:
#
# - `At(cell)` — the cell itself;
# - `Contains(cell)` — the same cell;
# - `Contains((lon, lat))` — the cell holding a point.

c = cells[300]
at = store[:elevation][DGG.Cells(DD.At(c))]
contains_cell = store[:elevation][DGG.Cells(DD.Contains(c))]
contains_point = store[:elevation][DGG.Cells(DD.Contains((67.5, 66.7)))]
at, contains_cell, contains_point

# `Covering(target)` selects every stored cell reached by the coverage of
# `target`. Here `target` is the extent of a level-2 ancestor, so the selection
# demonstrates the small spill beyond a region's exact boundary.

target = DGG.cell_extent(DGG.levelgrid(sys, 2), DGG.ancestor(sys, cells[300], 2))
region = store[:elevation][DGG.Cells(DGG.Covering(target))]

# The selected result has a compressed in-memory `CellLookup`; the source
# remains a `ChunkedCellLookup`.
#
# ## Which chunks a selection touches
#
# A `ChunkManifest` describes the chunk grid in cell-axis positions. It records
# each chunk's bounds and maps an axis position to its chunk.

manifest = DGG.chunkmanifest(axis, 128)
#
DGG.nchunks(manifest), length(manifest)
#
DGG.chunkbounds(manifest, 5)

# The elevation values equal their positions, so selected values reveal which
# positions were fetched. `chunkof` maps each position to its chunk:

selected = Int.(collect(region))
sort(unique(DGG.chunkof.(Ref(manifest), selected)))

# The figure shows stored cells coloured by chunk, the target as a dashed box,
# and the selected cells outlined. A compact spherical region can span several
# file chunks.

chunk = DGG.chunkof.(Ref(manifest), 1:n)
corners = [(target.X[1], target.Y[1]), (target.X[2], target.Y[1]),
           (target.X[2], target.Y[2]), (target.X[1], target.Y[2])]
box = [c1 .+ t .* (c2 .- c1) for (c1, c2) in zip(corners, circshift(corners, -1))
       for t in range(0, 1; length = 30)]

fig = Figure(size = (780, 470))
ax = GeoAxis(fig[1, 1]; dest = "+proj=laea +lon_0=44 +lat_0=64",
    limits = ((-16.0, 104.0), (45.0, 83.0)),
    xticks = 0:20:100, yticks = 50:10:80,
    title = "Stored cells by chunk; the cells Covering selects, outlined")
plt = dggpoly!(ax, cube[:elevation]; color = chunk,
    colormap = cgrad(:Set2, 5; categorical = true), colorrange = (0.5, 5.5))
lines!(ax, GeoMakie.coastlines(); color = ("#212529", 0.55), linewidth = 0.6)
dggpoly!(ax, region; color = :transparent, strokecolor = :black, strokewidth = 0.8)
lines!(ax, box; color = :black, linewidth = 2, linestyle = :dash)
Colorbar(fig[1, 2], plt; label = "chunk", ticks = 1:5)
fig

# How to make that outline land in fewer chunks is the subject of [Subzone
# layout](../api/subzone-layout.md).
#
# ## Choosing how the cell ids are stored
#
# An *encoding* describes how the store lays out cell ids. `encoding = :auto`
# chooses from the axis shape:
#
# | `encoding` | Stores | `:auto` picks it when |
# |---|---|---|
# | `:ranges` | `(n, 2)` inclusive `[start, stop]` id intervals | the axis is sorted, unique and one level |
# | `:dense` | one id per cell | otherwise; also the interop choice for readers without interval support |
#
# A ranges axis opens without reading coordinate data: length, chunk boundaries
# and selectors use rank/select arithmetic over its intervals. This store uses
# `RangesEncoding`, with one row per interval:

size(Zarr.zopen(path)["cell_id_ranges"], 2)

# `merge` chooses what one interval may span:
#
# | `merge` | A run is | Rows | Read back correctly by |
# |---|---|---|---|
# | `:step` (default) | ids adjacent as integers | more (91 here) | any reader that counts ids, grid-aware or not |
# | `:rank` | consecutive cells | fewest (1 here) | a rank-aware reader such as this package |
#
# This axis is a single run of consecutive cells, so under `:rank` it is one row:

ranks = DGG.dggwrite(joinpath(mktempdir(), "ranks.zarr"), cube; chunks = 128,
                     merge = :rank)
size(Zarr.zopen(ranks)["cell_id_ranges"], 2)

# ## Reading a store by URL
#
# `dggread` also opens a public `gs://`, `s3://` or `https://` store in place.
# A selection then fetches the chunks it needs. For example:
#
# ```julia
# pori = DGG.dggread("https://storage.googleapis.com/geo-assets/igeo7-zarr/pori_z7_r10.zarr")
# ```
#
# `dggwrite` writes to a local path or an open `Zarr.ZGroup`; publishing a
# remote store means uploading the directory it produced.
#
# [Out of core](out_of_core.md) sweeps a kernel over a store chunk by chunk,
# starting from a store like the one written here.
