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

## Reading a region back

A pyramid level is the COMPLETE level — at level 13 that is 968 890 104 072
cells — so nothing reads one whole. A region comes back by asking for the parts
of it you want, and how you ask matters more here than in a flat layout, because
a region is many runs over few chunks. The hydrology tutorial's tile is 9887
runs living in 76 chunks:

| reading its 16.17 M cells back | |
| --- | --- |
| run by run, no cache | 2.09 s |
| run by run, `cache = true` | 0.01 s |
| walking the runs, holding the current chunk | 0.04 s |

Each run visits a chunk, and without a cache each visit decompresses it again —
about 130 times per chunk here. Either fix works: pass `cache` to
[`dggread`](@ref DiscreteGlobalGrids.dggread), or walk the runs in order and keep
the chunk you are in. `cache` is off by default so that a level array shared
between tasks stays stateless.

[`cellvalues`](@ref DiscreteGlobalGrids.cellvalues) needs neither: it groups its
cells by chunk itself, which is why a frame costs one read per chunk touched.

## Compared with the other layouts

|  | `:cells` | `:subzones` | `:pyramid` |
| --- | --- | --- | --- |
| cell coordinate | stored | implicit column axis | none at all |
| levels per store | one | one | every level to the root |
| chunk | an element target | one ancestor subtree | one subtree, `a^k` slots |
| partial chunks | yes | refused | yes |
| xdggs can read it | yes | partly | no |

On the hydrology tutorial's tile — 16 172 725 cells at level 13, 64.7 MB of
`Float32` before compression — the three come out like this:

| | write | on disk | of which index | files |
| --- | --- | --- | --- | --- |
| `:cells`, `encoding = :dense` | 1.02 s | 48.0 MB | 0.87 MB | 35 |
| `:cells`, `encoding = :ranges` | 1.05 s | 47.4 MB | 0.27 MB | 19 |
| `:pyramid` | 0.74 s | 55.7 MB | none | 227 |
| `:subzones` | — | refused: the coverage is not snapped to subtrees | | |

Two things worth reading off that table. **The explicit coordinate is cheap**:
16.17 million `UInt64` ids are 129 MB raw and 0.87 MB stored, because a cell axis
ascends and the compressor is very good at that — dropping it saves under two
per cent. **The overviews are not free**: the thirteen of them add 8.4 MB to the
leaf level's 47.3 MB, eighteen per cent, which is the `1/6` the aperture implies.

So the pyramid is not the small option. What it buys for that eighteen per cent
is every level of the hierarchy and an index that needs no probing — and the
ability to store this cube at all where `:subzones` cannot, because a
`MultiOrderCoverage` of a lon/lat box refines at its boundary and leaves ancestor
subtrees partly covered.

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
