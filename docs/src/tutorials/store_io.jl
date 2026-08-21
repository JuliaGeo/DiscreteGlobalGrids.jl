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

# Two layers over that one axis, as a `DimStack`. The values are deterministic
# so the round trip below is checkable rather than plausible.

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

# The lookup is a `ChunkedCellLookup`: it answers the same selectors as a
# `CellLookup` — `At` and `Contains` on a cell, `Contains` on a lon/lat point,
# `Covering` on a region — but resolves them against the store's chunk grid
# rather than by scanning it. That chunk grid is the `ChunkManifest`, described
# in cells: for every chunk, the first and last id it holds and how many. It is
# what says which chunks a selection touches, and `nchunks`, `chunkof` and
# `chunkbounds` are how it is asked.

manifest = DGG.chunkmanifest(axis, 128)
(; chunks = DGG.nchunks(manifest), cells = length(manifest))

# ## Selecting a region out of a store
#
# `Covering(target)` runs a coverage of `target` against the axis and keeps the
# positions it lands on. The cheapest interesting target here is the extent of a
# coarse ancestor of one of the stored cells: a level-2 cell, a small piece of
# what was written.

target = DGG.cell_extent(DGG.levelgrid(sys, 2), DGG.ancestor(sys, cells[300], 2))
region = store[:elevation][DGG.Cells(DGG.Covering(target))]

# A subset is no longer a stored axis: the cells it names are materialised, and
# the result is the package's own compressed `CellLookup` over them.

(; ncells = length(region), lookup = nameof(typeof(DD.lookup(region, DGG.Cells))))

# ## Encodings
#
# How the cell ids are laid out in the store is the *encoding*, and
# `dggwrite`'s default `encoding = :auto` chooses one: `RangesEncoding` where
# the axis is eligible — sorted, unique, and all one level — and
# `DenseEncoding` otherwise, which is always eligible and so makes `:auto`
# total. Ranges store `[start, stop]` intervals instead of ids, and the axis is
# recovered from them by rank/select arithmetic, so opening one costs no data
# IO at all however many cells it names. `encoding = :dense` is the interop
# escape for readers that cannot expand intervals; it writes one id per cell.
#
# `merge` chooses what an interval is allowed to hold. The default `:step`
# merges ids that are adjacent AS INTEGERS, so no interval encloses an id naming
# no cell and a reader that counts well-formed ids rather than cells recovers
# the same axis. `merge = :rank` merges consecutive CELLS instead: the fewest
# rows — this axis is rank-contiguous, so a single `[start, stop]` row — read
# back correctly only by a rank-aware reader. The axis above is eligible, and
# went to disk under the default:

(; encoding = DD.metadata(store)["encoding"],
   rows = size(Zarr.zopen(path)["cell_id_ranges"], 2))

# ## Reading a store by URL
#
# `dggread` takes a `Zarr.ZGroup`, a local path, or a URL, so a public store is
# read where it lives and only the chunks a selection touches are fetched. The
# Pori stores this format was reverse-engineered from are readable over HTTPS —
# not run here, because these pages build without a network:
#
# ```julia
# pori = DGG.dggread("https://storage.googleapis.com/geo-assets/igeo7-zarr/pori_z7_r10.zarr")
# ```
#
# Writing goes the other way round: `dggwrite` writes locally, and a remote
# store is an upload of what it produced.
#
# A store axis is single-level. The mixed-level arrays of the multi-order
# storage page reach a store through `expand`, which presents them at one
# level; `dggwrite` refuses a `MultiOrderLookup` axis directly.
