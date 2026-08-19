# DiscreteGlobalGrids.jl

Six discrete global grid systems — IGEO7, H3, HEALPix, A5, S2, ISEA4R — behind
one small interface, with every algorithm written against the interface exactly
once.

`DGG.CopernicusDEMSystem` implements that same interface for the Copernicus DEM
raster lattice, but is a raster rather than a DGGS and so stays outside
`DGG.systems()`; `examples/copernicus_dem.jl` puts it to work.

The mental model is two tiers. A **grid** is one finite collection of cells on
the sphere — a complete level, or a regional subset of one — and geometry,
stencils and queries are all answered there. A **system** adds the parent/child
hierarchy across levels, always as a fast path: hierarchy is an optimisation,
never a semantic. A bare `Int` is a position in `1:ncells(grid)`; a typed cell
id knows its own level.

```@example index
import DiscreteGlobalGrids as DGG

grid = DGG.levelgrid(DGG.IGeo7System(), 4)  # the complete level: system + level, nothing else
c = DGG.cellat(grid, 8.5, 47.4)             # the cell under Zürich, as a typed id
DGG.cell_area(grid, c)                      # steradians, on the unit sphere
```

```@example index
import Extents
DGG.query(grid, DGG.Intersects(Extents.Extent(X = (5, 12), Y = (45, 50))))
```

Swap `IGeo7System()` for any of the six and nothing else changes. The
[README](https://github.com/JuliaGeo/DiscreteGlobalGrids.jl) walks the whole
surface; the docstring of `DGG.systems()` is the comparison table.

## Installation

The package is not in the General registry yet, so it installs from the
repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/JuliaGeo/DiscreteGlobalGrids.jl")
```

`Pkg.develop(url = ...)` instead for a checkout you intend to edit. Julia 1.11
or newer is required. Two capabilities ride in extensions and are loaded by
their package: `using Makie` (or a backend) draws cells, and `using Zarr` turns
`dggread`/`dggwrite` from stubs that only report their absence into methods.

## Where to go next

The [DGGS gallery](all_dggs.md) draws every system. Each tutorial is the
shortest honest path to one result:

  - [Stencil operations](tutorials/stencils.md) — smoothing, Laplacians and
    diffusion from each cell's neighbourhood, on a whole level and on a subset.
  - [Zonal statistics](tutorials/zonal.md) — reduce a field over regions, with
    spatial queries.
  - [Regridding a time series](tutorials/regridding.md) — conservative
    regridding between a lon/lat raster and a DGGS.
  - [Multi-order coverage](tutorials/multiorder.md) — one region at every
    resolution at once: coarse cells inside, leaf cells along the boundary.
  - [Hydrology: a DEM on an IGEO7 grid](tutorials/hydrology.md) — elevation
    data on a regional subset, and flow routing across it.
  - [The sky in HEALPix](tutorials/healpix_astronomy.md) — nested order, cone
    searches, and a galactic-plane cut.
  - [A round trip through a DGGS store](tutorials/store_io.md) — `dggwrite` and
    `dggread` over a Zarr store, with the grid carried in the lookup type.

Two API pages are rendered so far.
[Subtree and subset boundaries](api/boundaries.md) covers the boundary family —
`subtree_border`, `subtree_halo`, `halo`, `halo_positions`, `halo_sizehint` and
the engines behind them — and
[Reading and writing DGGS stores](api/store-io.md) covers `dggread`, `dggwrite`
and the conventions, encodings and chunked axis they are built on.
