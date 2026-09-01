# # A round trip through a DGGS store
#
# A cube over a cell axis is an ordinary `DimArray` or `DimStack` whose one
# spatial dimension is `Cells`. `dggwrite` puts such a cube in a Zarr store and
# `dggread` opens one back.
#
# Both methods live in the Zarr.jl extension: `using Zarr` loads them.

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
using Zarr
using DiscreteGlobalGridsVisualization: dggpoly, dggpoly!
using GLMakie, GeoMakie
GLMakie.activate!(inline = true)

# ## A cube to write
#
# The first two level-1 IGEO7 cells, expanded to their level-4 descendants: a
# regional store in miniature. `CellVector` names the cells, `CellLookup` reads
# them as a one-level cell axis, and `Cells` makes that axis a cube dimension.

sys = DGG.IGeo7System()
roots = DGG.CellVector(DGG.levelgrid(sys, 1))[1:2]
cells = sort!(reduce(vcat, [collect(DGG.descendants(sys, c, 4)) for c in roots]))
lookup = DGG.CellLookup(DGG.CellVector(sys, 4, cells))

# Two layers over that one axis, as a `DimStack`. Elevation doubles as a
# position marker: value `k` sits at axis position `k`. The chunk section reads
# positions back out of a selection that way.

n = length(cells)
elevation = Float32.(1:n)
slope = Float32.(0.5 .* (1:n))
cube = DD.DimStack((; elevation, slope), (DGG.Cells(lookup),))

# ## Write the cube and read it back
#
# `dggwrite` returns its destination, ready to hand to `dggread`. `chunks = 128`
# fixes the chunk length in cells and gives this small cube five chunks to look
# at; the default `:auto` targets a million elements per chunk and would put the
# whole cube in one.

path = DGG.dggwrite(joinpath(mktempdir(), "demo.zarr"), cube; chunks = 128)
store = DGG.dggread(path)
#
DD.metadata(store)["description"]

# The description is read from the store's attributes: grid system, level, and
# the layout of the ids. The arrays `dggread` returns are lazy; `collect` pulls
# the axis and the values into memory to compare them with what went in:

axis = DD.lookup(store[:elevation], DGG.Cells)
#
collect(axis) == cells, collect(parent(store[:elevation])) == elevation

# ## Selecting a region out of a store
#
# The axis is a `ChunkedCellLookup`. It answers the same selectors a
# `CellLookup` does and resolves each through the chunk manifest, touching at
# most one chunk of ids. Three selectors name a single cell:
#
# - `At(cell)` — the cell itself;
# - `Contains(cell)` — the same cell;
# - `Contains((lon, lat))` — the cell holding a point.

c = cells[300]
at = store[:elevation][DGG.Cells(DD.At(c))]
contains_cell = store[:elevation][DGG.Cells(DD.Contains(c))]
contains_point = store[:elevation][DGG.Cells(DD.Contains((67.5, 66.7)))]
at, contains_cell, contains_point

# `Covering(target)` selects a region: every stored cell the coverage of
# `target` lands on. The target is the extent of a level-2 cell, an ancestor of
# one stored cell, which makes it a small piece of what was written. A coverage
# is a superset of its region, and the selection spills a little past the box.

target = DGG.cell_extent(DGG.levelgrid(sys, 2), DGG.ancestor(sys, cells[300], 2))
region = store[:elevation][DGG.Cells(DGG.Covering(target))]

# The `dims` line of the result shows a `CellLookup`: a selection materialises
# the cells it names and carries them as the package's compressed in-memory
# axis. The `ChunkedCellLookup` stays with the store.
#
# ## Which chunks a selection touches
#
# A `ChunkManifest` describes the store's chunk grid in cells: for each chunk,
# its first id, its last id and its length. `nchunks` counts chunks, `chunkof`
# maps an axis position to its chunk, and `chunkbounds` maps a chunk to the
# positions it holds.

manifest = DGG.chunkmanifest(axis, 128)
#
DGG.nchunks(manifest), length(manifest)
#
DGG.chunkbounds(manifest, 5)

# The selection's values are stored positions, by construction of the elevation
# layer. `chunkof` maps each position to the chunk a reader fetches for it:

selected = Int.(collect(region))
sort(unique(DGG.chunkof.(Ref(manifest), selected)))

# The figure shows the read: stored cells coloured by chunk, the level-2 target
# as a dashed box, the cells `Covering` returned outlined. The outline crosses
# three colours — a compact region on the sphere spans three chunks of the file.

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
# An *encoding* is how the store lays out its cell ids. `encoding = :auto` picks
# one by the shape of the axis:
#
# | `encoding` | Stores | `:auto` picks it when |
# |---|---|---|
# | `:ranges` | `(n, 2)` inclusive `[start, stop]` id intervals | the axis is sorted, unique and one level |
# | `:dense` | one id per cell | otherwise; also the interop choice for readers without interval support |
#
# A ranges axis opens with zero data IO at any size: its length, its chunk
# boundaries and every selector are closed-form rank/select arithmetic over the
# intervals. The description above named this store's encoding,
# `RangesEncoding`; its cost is one row per interval:

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
# `dggread` opens a public store in place from a `gs://`, `s3://` or `https://`
# URL as readily as from a local path or a `Zarr.ZGroup`; a selection fetches
# only the chunks it touches. A published IGEO7 store over Pori, Finland, reads
# like this (displayed only: the docs build offline):
#
# ```julia
# pori = DGG.dggread("https://storage.googleapis.com/geo-assets/igeo7-zarr/pori_z7_r10.zarr")
# ```
#
# `dggwrite` writes to a local path or an open `Zarr.ZGroup`; publish a remote
# store by uploading the directory it produced.
#
# [Out of core](out_of_core.md) sweeps a kernel over a store chunk by chunk,
# starting from a store like the one written here.
