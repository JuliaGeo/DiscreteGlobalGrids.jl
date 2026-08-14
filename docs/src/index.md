# DiscreteGlobalGrids.jl

Six discrete global grid systems — IGEO7, H3, HEALPix, A5, S2, ISEA4R — behind
one small interface, with every algorithm written against the interface exactly
once.

```@example index
import DiscreteGlobalGrids as DGG
import Extents

# A complete level is the entry point. It stores the system and the level and
# nothing else, so building one is free.
grid = DGG.levelgrid(DGG.HEALPixSystem(), 4)
DGG.ncells(grid)
```

A bare `Int` is always a **position** in `1:ncells(grid)`. A typed
`AbstractCellIndex` is always an **identity**, self-describing about its level,
so no call passes a level and an id side by side. The two are one bijection:

```@example index
c = DGG.cellindex(grid, 1000)          # position -> identity
DGG.cellposition(grid, c), DGG.level(c)
```

Geometry is on the unit sphere; the `(lon, lat)` wrappers take degrees.

```@example index
DGG.cell_area(grid, c),                # steradians
DGG.cellat(grid, 8.5, 47.4),           # the cell at a point
length(DGG.neighbors(grid, c))         # ring 1, CCW seen from outside
```

Spatial queries take DE9IM predicate types and answer with sorted typed ids:

```@example index
DGG.query(grid, DGG.Intersects(Extents.Extent(X = (5, 12), Y = (45, 50))))
```

The hierarchy hangs off the system rather than the grid, because an id already
knows its level:

```@example index
sys = DGG.HEALPixSystem()
DGG.children(sys, c), DGG.descendant_range(sys, c, 6), length(DGG.subtree_border(sys, c, 6))
```

Swapping `HEALPixSystem()` for `IGeo7System()`, `H3System()`, `A5System()`,
`S2System()` or `ISEA4RSystem()` changes nothing else, and `AuthalicSystem`
wraps any of them to read geometry at geodetic latitude. `DGG.systems()` lists
all six, and its docstring is the comparison table.

## Where to go next

The [DGGS gallery](all_dggs.md) draws every system. The tutorials are each the
shortest honest path to one result:

  - [Stencil operations](tutorials/stencils.md) — `halo_table` on a whole level
    and on a subset of one, then smoothing, Laplacians and diffusion.
  - [Zonal statistics](tutorials/zonal.md) — `query` with `Within` and
    `Intersects`, and the multi-order coverage that compresses the answer.
  - [Regridding a time series](tutorials/regridding.md) — conservative
    regridding onto a DGGS, with no adapter between the two packages.
  - [Multi-order coverage](tutorials/multiorder.md) — one region at every
    resolution at once: giant cells inside, leaf cells along the coastline.
  - [Hydrology: a DEM on an IGEO7 grid](tutorials/hydrology.md) —
    `MultiOrderCoverage`, `PartialGrid` over one subtree, and flow routing off
    its `halo_table`.
  - [The sky in HEALPix](tutorials/healpix_astronomy.md) — nested order,
    cone searches against a `SphericalCap`, and a galactic-plane cut.
