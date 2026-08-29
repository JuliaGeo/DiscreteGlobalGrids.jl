# # A round trip through a DGGS store
#
# A cube over a cell axis is an ordinary `DimArray` whose one spatial dimension
# is `Cells`. `dggwrite` puts such a cube in a Zarr store and `dggread` opens
# one back, and neither invents a container type: what goes in is
# DimensionalData and what comes out is DimensionalData. The grid SYSTEM is in
# the type of the lookup that comes back, the level is a field of the grid it
# holds, and what is neither — orientation, ellipsoid, the layout the ids were
# stored in — rides in the stack's metadata.
#
# The methods live in an extension on Zarr.jl, so `using Zarr` is what turns the
# two stubs into functions.

import DiscreteGlobalGrids as DGG
import DimensionalData as DD
using Zarr

# ## A cube to write
#
# The first two level-1 IGEO7 cells, expanded to their level-4 descendants: 629
# cells, which is a regional store in miniature. It is not 686 because one of
# the two roots is a pentagon, and a pentagon's subtree is short. `CellVector`
# names the cells, `CellLookup` reads them as a one-level cell axis, and `Cells`
# makes that axis a cube dimension.

sys = DGG.IGeo7System()
roots = DGG.CellVector(DGG.levelgrid(sys, 1))[1:2]
cells = sort!(reduce(vcat, [collect(DGG.descendants(sys, c, 4)) for c in roots]))
lookup = DGG.CellLookup(DGG.CellVector(sys, 4, cells))

# Deterministic values make the two-layer `DimStack` round trip exact and
# reproducible.

n = length(cells)
elevation = Float32.(1:n)
slope = Float32.(0.5 .* (1:n))
cube = DD.DimStack((; elevation, slope), (DGG.Cells(lookup),))

# ## Out and back
#
# `dggwrite` returns its destination, so it composes. `chunks = 128` fixes the
# chunk length in cells; the default `:auto` aims each chunk at a million
# elements instead, which for 629 cells is one chunk and nothing to look at.

path = DGG.dggwrite(joinpath(mktempdir(), "demo.zarr"), cube; chunks = 128)
store = DGG.dggread(path)

# The axis came back as the cells that went in, and the values with it. What
# `dggread` hands over is a lazy store-backed array, so `collect` is what
# forces the comparison.

axis = DD.lookup(store[:elevation], DGG.Cells)
(; ncells = length(axis), level = DGG.level(axis),
   axis_ok = collect(axis) == cells,
   values_ok = collect(parent(store[:elevation])) == elevation)

# A `ChunkedCellLookup` resolves these selectors through the store's chunk grid:
#
#   - `At` and `Contains` select cells.
#   - `Contains` also selects a longitude/latitude point.
#   - `Covering` selects a region.
#
# The `ChunkManifest` records every chunk's first id, last id and cell count.
# `nchunks`, `chunkof` and `chunkbounds` query that manifest.

manifest = DGG.chunkmanifest(axis, 128)
(; chunks = DGG.nchunks(manifest), cells = length(manifest))

# ## Selecting a region out of a store
#
# `Covering(target)` runs a coverage of `target` against the axis and keeps the
# indices it lands on. The cheapest interesting target here is the extent of a
# coarse ancestor of one of the stored cells: a level-2 cell, a small piece of
# what was written.

target = DGG.cell_extent(DGG.levelgrid(sys, 2), DGG.ancestor(sys, cells[300], 2))
region = store[:elevation][DGG.Cells(DGG.Covering(target))]

# A subset is no longer a stored axis: the cells it names are materialised, and
# the result is the package's own compressed `CellLookup` over them.

(; ncells = length(region), lookup = nameof(typeof(DD.lookup(region, DGG.Cells))))

# ## Encodings
#
# The *encoding* controls the cell-id layout. `encoding = :auto` selects:
#
#   - `RangesEncoding` for a sorted, unique, single-level axis;
#   - `CompactedEncoding` for a mixed-level axis; and
#   - `DenseEncoding` for every other single-level axis.
#
# Ranges store `[start, stop]` intervals and reconstruct the axis with
# rank/select arithmetic, so opening them requires no cell-id reads.
# `encoding = :dense` writes one id per cell for readers without range support.
#
# `merge` controls range formation. The default `:step` merges integer-adjacent
# ids, giving structural and cell-aware readers the same counts. `:rank` merges
# consecutive cells for the fewest rows and requires a rank-aware reader. The
# eligible axis above uses the default:

(; encoding = DD.metadata(store)["encoding"],
   rows = size(Zarr.zopen(path)["cell_id_ranges"], 2))

# ## Reading a store by URL
#
# `dggread` accepts a `Zarr.ZGroup`, local path or URL. A public remote store
# stays lazy, and selections fetch only the chunks they touch. The Pori stores
# that established this format are available over HTTPS; this example remains
# unevaluated so the documentation builds offline:
#
# ```julia
# pori = DGG.dggread("https://storage.googleapis.com/geo-assets/igeo7-zarr/pori_z7_r10.zarr")
# ```
#
# Remote publication starts with a local `dggwrite` and then uploads the result.
#
# `dggwrite` stores a mixed-level array from the
# [multi-order storage tutorial](moc_storage.md) in the `compacted`
# layout: two aligned columns, `cell_ids` and `cell_levels`, under
# `refinement_level: null`. `dggread` restores a
# [`MultiOrderLookup`](@ref) axis. Single-level encodings require
# [`expand`](@ref) to present the data at one level.
