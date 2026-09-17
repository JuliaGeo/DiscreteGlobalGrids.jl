# The pyramid layout

Use this layout to store every level of a hierarchy at once, with no cell
coordinate at all: each level is one array, laid out on the level's **slot
space** — its id space with the addresses naming no cell left in — and chunked
at a power of the system's aperture.

```
elevation/level0     12 slots,     chunked 7^0
elevation/level1     84 slots,     chunked 7^1
…
elevation/level13    12·7^13 slots, chunked 7^6
```

Two properties follow from that, and the layout exists for both.

  - **A chunk is a subtree.** Slots `(c-1)·7^k+1 : c·7^k` of level `J` are
    exactly the descendants of slot `c` of level `J-k`. A region written into
    the store touches only the chunks its subtrees fall in, and Zarr never
    stores the rest, so half an earth costs half an earth.
  - **The coarse level is the chunk index.** `level{J-k}` has one value per
    chunk of `level{J}`, at the same position, and a parent is fill exactly when
    every cell under it is. A reader learns which chunks of a deep level exist
    by reading a shallow one — no `LIST`, no per-chunk probe. That is why the
    overviews are not optional here.

The price is the addresses that name no cell: a sixth of IGEO7's slots, the
twelve pentagons' deleted branches. They fall in aligned blocks of `7^(l-1) …
7^0` per base, so all but the smallest are whole chunks that are never written,
and a full store really carries `2·(7^k − 1)` fill values per level — 684 at the
default `k = 6`.

Write it with `dggwrite(dest, cube; layout = :pyramid)` and read it with
[`dggread`](@ref DiscreteGlobalGrids.dggread), which returns a
[`StorePyramid`](@ref DiscreteGlobalGrids.StorePyramid) — every level, lazily.
Reading one level is `pyr[7]`; reading one variable is `pyr[:elevation]`.

```julia
import DiscreteGlobalGrids as DGG
using Zarr

DGG.dggwrite("dem.zarr", elevation; layout = :pyramid)

pyr = DGG.dggread("dem.zarr")[:elevation]
DGG.levels(pyr)                     # 0:13
DGG.holdsdata(pyr, cell)            # is anything under this cell?
DGG.cellvalues(pyr, 7, cells)       # a frame's worth, one read per chunk
```

## Fill is the presence marker

An absent cell and a deleted address read back as the same value, and that is
deliberate: a reader asking "is there data here" gets one answer, and the coarse
level can then be the chunk index. Two consequences:

  - a cube whose element type has no `NaN` must name a `fill_value`, because a
    store whose absent cells read back as zero cannot tell a written zero from
    an unwritten cell;
  - the writer decides fill by **presence**, not by what the aggregate returned.
    A parent with at least one present child is written non-fill, and an
    `aggregate` that returns the fill value raises rather than marking a
    populated subtree empty.

## Aggregating up

Every level carries the cube's own element type, so `aggregate` has to land back
in it. It defaults to `mean` for a floating-point cube; a cube of another
element type has to name one — `mode` for a class, `vs -> round(Int16, mean(vs))`
for a rounded mean.

A parent is reduced from the children that carry a value, not from all seven, so
a region's coarse levels stay usable right up to its boundary. A cell on that
boundary is a mean over three of seven children with nothing recording it; where
that matters, store the count as a second variable.

## Compared with the other layouts

|  | `:cells` | `:subzones` | `:pyramid` |
| --- | --- | --- | --- |
| cell coordinate | stored | implicit column axis | none at all |
| levels per store | one | one | every level to the root |
| chunk | an element target | one ancestor subtree | one subtree, `a^k` slots |
| partial chunks | yes | refused | yes |
| xdggs can read it | yes | partly | no |

The pyramid layout writes no `zarr_conventions` declaration and no array of ids,
so nothing generic can fingerprint it; it is read by a reader that knows the
`dggs.pyramid_layout` block below.

```@docs
DiscreteGlobalGrids.StorePyramid
DiscreteGlobalGrids.PyramidLayout
DiscreteGlobalGrids.AbstractPyramid
DiscreteGlobalGrids.holdsdata
DiscreteGlobalGrids.cellvalues
```

## The slot space

```@docs
DiscreteGlobalGrids.slotcount
DiscreteGlobalGrids.slotindex
DiscreteGlobalGrids.slotcell
```

## Chunks and subtrees

```@docs
DiscreteGlobalGrids.levelname
DiscreteGlobalGrids.levelfromname
DiscreteGlobalGrids.levelslots
DiscreteGlobalGrids.chunkexponent
DiscreteGlobalGrids.chunklength
DiscreteGlobalGrids.chunkrootlevel
DiscreteGlobalGrids.chunkcount
DiscreteGlobalGrids.slotchunk
DiscreteGlobalGrids.chunkroot
DiscreteGlobalGrids.chunkindices
DiscreteGlobalGrids.chunkslots
DiscreteGlobalGrids.cellslot
DiscreteGlobalGrids.PyramidRun
DiscreteGlobalGrids.pyramid_runs
DiscreteGlobalGrids.pyramid_index_runs
DiscreteGlobalGrids.pyramid_reduce
```

## Attributes

The `dggs` attribute object carries the grid name and the finest level, with
layout fields nested under `pyramid_layout`.

```@docs
DiscreteGlobalGrids.pyramid_attrs
DiscreteGlobalGrids.pyramid_level_attrs
DiscreteGlobalGrids.pyramid_dimension
DiscreteGlobalGrids.pyramid_layout
DiscreteGlobalGrids.ispyramidstore
DiscreteGlobalGrids.PYRAMID_LAYOUT
DiscreteGlobalGrids.PYRAMID_LEVEL_PREFIX
DiscreteGlobalGrids.PYRAMID_ORDER
DiscreteGlobalGrids.PYRAMID_PADDING
DiscreteGlobalGrids.DEFAULT_PYRAMID_CHUNK_EXPONENT
```
