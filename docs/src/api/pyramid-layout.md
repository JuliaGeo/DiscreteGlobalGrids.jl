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

## Reading a subset, which is what plotting asks for

The table above reads everything. A view reads a part, and that is where the
layouts separate. Same tile, same machine, against `:cells` with the ranges
encoding:

| | `:cells` | `:pyramid` | |
| --- | --- | --- | --- |
| open the store | 0.192 s | 0.0016 s | 120× |
| a 0.1° window at level 13 (163 098 cells) | 0.082 s | 0.0078 s | 10× |
| the same window at level 8 (1 038 cells) | 9.84 s | 0.0017 s | 5 700× |

Three different reasons, worth separating.

**Opening** is `O(n)` for a stored coordinate and `O(1)` for an implicit one:
`:cells` has to read its cell axis and check it before it can answer anything,
while a pyramid's axis is arithmetic on the id. That gap widens with the store.

**Locating** a window is arithmetic either way, but `:cells` resolves it against
the stored axis and a pyramid against the slot number, so the pyramid reads only
the chunks the window touches and nothing about where they are. The cost of "which
subtree covers this cell range" is a handful of base-7 digit operations — it never
touches the store.

**The coarse view** is not a speedup at all, it is a different algorithm. A
single-level store has no level 8; getting one means reading all 16 million
leaves and aggregating them, every time. The pyramid reads the overview that was
computed once at write time. This is the whole point of the layout, and it is
the number that grows with the data: at ten times the cells the first column is
ten times slower and the second is unchanged.

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

## Where it degenerates

Every one of these layouts assumes some spatial coherence, and each loses it
differently. With `n` cells stored, `r` maximal runs of consecutive indices, `C`
chunks touched, `W = a^k` slots per chunk and a value type of `v` bytes, they lay
down, before the compressor:

| | logical bytes | fails when |
| --- | --- | --- |
| `:cells`, dense | `v·n + 8n` | never — but the index is `Θ(n)` forever, twice the data at `Float32` |
| `:cells`, ranges | `v·n + 16r` | `r > n/2`, i.e. the mean run drops below two; then it is worse than dense |
| `:subzones` | `v·C'·capacity` | any subtree is covered in part — it refuses outright |
| `:pyramid` | `v·C·W·7/6` | the fill `n/(C·W)` collapses; note the cost does not mention `n` at all |

That last row is the one to internalise. **A pyramid pays for the chunks it
touches, full or not**, so its cost is flat in the number of cells and its enemy
is scatter rather than size. `benchmark/pyramid_sparsity.jl` measures the whole
range at IGEO7 level 13, on disk, with incompressible values:

| case | cells | runs | chunks | fill | dense | ranges | pyramid |
| --- | --- | --- | --- | --- | --- | --- | --- |
| tile, 1°×1° | 16 181 892 | 9 909 | 170 | 80.9% | 58.1 MB / 35 | 57.5 MB / 19 | 67.1 MB / 227 |
| straddling a base cell | 16 177 029 | 10 502 | 167 | 82.3% | 58.1 MB / 41 | 57.5 MB / 22 | 67.1 MB / 233 |
| 1°×1° at the north pole | 210 558 | 5 239 | 20 | 8.9% | 0.8 MB / 3 | 0.8 MB / 3 | 1.1 MB / 42 |
| thin strip, 1° × 0.01° | 165 737 | 5 250 | 14 | 10.1% | 0.6 MB / 3 | 0.6 MB / 3 | 0.9 MB / 38 |
| every 7th cell | 2 311 699 | 2 311 699 | 170 | 11.6% | 8.4 MB / 7 | 8.5 MB / 5 | 35.2 MB / 227 |
| every 49th cell | 330 243 | 330 243 | 169 | 1.7% | 1.2 MB / 3 | 1.2 MB / 3 | 10.2 MB / 226 |
| 1 000 scattered in the tile | 1 000 | 1 000 | 157 | ~0 | 0.0 MB / 3 | 0.0 MB / 3 | 0.5 MB / 212 |
| 1 000 scattered on the earth | 1 000 | 1 000 | 1 000 | ~0 | 0.0 MB / 3 | 0.0 MB / 3 | 10.6 MB / 5 413 |

Five things the arithmetic above does not tell you.

  - **A rectangle is fine wherever you put it.** Straddling a Z7 base-cell
    boundary costs nothing measurable — 82.3% fill against 80.9% mid-face. At
    the chunk root level a one-degree box already spans about 170 cells, so it
    crosses subtree boundaries wherever it sits, and the level-0 ones are not
    special. What hurts is a large perimeter for the area: the polar box and the
    thin strip fall to 9–10% fill.
  - **Compression rescues most of the amplification.** The earth-scatter case is
    549 MB of logical slots and 10.6 MB on disk, because solid fill is a constant
    and the compressor eats it.
  - **It does not rescue interleaved sparsity.** Every 7th cell costs 4.2× dense
    and every 49th 8.5×, because a chunk of alternating fill and value has no
    runs to squeeze — unlike a chunk that is all fill.
  - **The real wall is file count.** A thousand cells scattered over the earth
    cost 5 413 files: five per stored value, one object-store request each. The
    coarse levels contribute most of them, because scattered cells stay scattered
    all the way up.
  - **`ranges` really does cross over.** Every 7th cell makes every run length
    one, and it lands at 8.5 MB against dense's 8.4 — the predicted flip at mean
    run length two, muted on disk because both halves compress.

The chunk exponent trades bytes against files, and which way depends on whether
the sparsity is finer or coarser than a chunk:

| | `7^2` | `7^3` | `7^4` | `7^6` |
| --- | --- | --- | --- | --- |
| the whole tile | 82.0 MB / 387 441 | 73.1 MB / 55 914 | 69.8 MB / 8 213 | 67.1 MB / 227 |
| every 49th cell | 32.3 MB / 386 000 | 19.9 MB / 55 758 | 14.3 MB / 8 187 | 10.2 MB / 226 |
| 1 000 over the earth | 0.7 MB / 9 363 | 1.0 MB / 8 375 | 1.1 MB / 7 387 | 10.6 MB / 5 411 |

Thinning touches every chunk whatever the chunk size, so a bigger chunk is
strictly better there. Isolated points waste one chunk each, so a smaller chunk
is better on bytes and worse on files. **If your sparsity is scatter rather than
shape, use `:cells`** — that is what its explicit index is for, and it costs
under two per cent.

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
