# Locating a point in a partial grid

- Date: 2026-08-25
- Benchmark: `benchmark/copdem_nearest.jl`, unchanged
- Companion: `2026-08-25-copdem-nearest-performance.md`, which measured the
  defect this fixes. Its sections 3-7, 9 and 11 are the "before" column here.

`PartialGrid` defined no `cellat`, so a subset of a level over a system that
names the cell containing a point itself still went to the generic locator: a
spatial-tree query over the subset's cells, a candidate sort, and an exact
point-in-polygon test per candidate. This note measures what defining that
method costs and what it removes.

The headline is that a Copernicus DEM source cell lookup falls from 6,085 ns and
347 bytes to 166 ns and no bytes, that the four-orders-of-magnitude cliff for a
tile set crossing a latitude row is gone from the nearest path entirely, and
that the 2x2-degree box the earlier note called unreachable now sweeps
9,160,874 destination cells in 6.04 seconds.

## 1. The defect, as measured

The companion note's section 4 priced two routes to the same answer on the
Copernicus source over one 1x1-degree tile, 100,000 sample sites:

| route | ns per call | bytes per call |
|---|---|---|
| `cellat(levelgrid(sys, 1), p)` then `localindex(grid, c)` | 102.2 | 0 |
| `cellat(space, p)` on the `DGGSpace` over that grid, as shipped | 6,366.6 | 346.7 |

`CopernicusDEMSystem` has a closed-form `cellat` on its complete level grid — a
latitude-band lookup and two floors — and so do IGeo7, H3, A5, S2, HEALPix and
ISEA4R on theirs. None of them was reached, because the grid the source is
built on is a `PartialGrid` and `cellat(grid::AbstractGrid, p)` in
`src/fallbacks/locate.jl` was the only method that matched it.

Two consequences followed. Location was 84.7% of a cold destination-tile read
(companion section 9), which made `NearestCell` 2.1x slower than `Conservative`
on the same data (section 5) although a nearest stencil is one entry of one. And
because the tree the locator queried is the one `treeify` gives, the locator
inherited `treeify(::PartialGrid{<:CopernicusDEMSystem})`'s split: a tile set
that is one contiguous id run gets the block cursor, anything else gets the
generic hierarchical cursor, and one `cellat` there cost 662 ms instead of
7 us. Every tile set crossing a latitude row is on the wrong side, which is
nearly every production source, and a 2x2-degree box could not finish a single
destination-tile read in 30 seconds of profiling.

## 2. The change

`cellat(grid::PartialGrid, p)` asks the complete level for the cell and answers
it when it is a member of `ids`, `nothing` otherwise. That is the whole method.
The generic locator in `locate.jl`, `treeify`, and the Copernicus cursor are
untouched; `DGGSpace`'s `cellat` and `CellVector`'s already reduced to the
grid's, so both inherit the change without a line.

Composing is only right where the complete level itself locates directly. Where
the complete level would fall to the same generic search, a tree over the whole
level is strictly worse than a tree over the subset, so such a partial grid
keeps the fallback. That is a fact about the system, and the code base had no
name for it, so the change adds one:

    has_direct_location(sys::AbstractHierarchicalGridSystem) -> Bool

`false` by default, declared per system beside `has_sorted_subtrees`, which is
the trait it is modelled on. It is a per-type constant, so the branch folds at
compile time and the composed path costs no dispatch and no reflection. A
declaration can drift from the methods it stands for, so the cross-system suite
asserts, for every registered system, that the trait says `true` exactly when
`cellat` on that system's complete level grid selects a method of the system's
own rather than the interface-wide one.

### Which side each system lands on

| system | `has_direct_location` | why |
|---|---|---|
| `IGeo7System` | true | `cellat` on its level grid: the gnomonic inverse and the Z7 digit walk name the cell from the point |
| `H3System` | true | `cellat` on its level grid: libh3's own `latLngToCell` |
| `A5System` | true | `cellat` on its level grid: A5's `lonlat_to_cell` inverse projection |
| `S2System` | true | inherited from `AbstractQuadFaceGridSystem`; its own `cellat` is `point_to_xyf`, the chart's analytic inverse |
| `HEALPixSystem` | true | inherited from `AbstractQuadFaceGridSystem`; its own `cellat` is the chart inverse |
| `ISEA4RSystem` | true | inherited from `AbstractQuadFaceGridSystem`; its own `cellat` is the chart inverse |
| `CopernicusDEMSystem{N}` | true | `cellat` on its level grid: a latitude-band lookup and two floors |
| `AuthalicSystem(sys)` | `has_direct_location(sys)` | its grid's `cellat` warps the query point and forwards, so it locates directly exactly when the system underneath does |
| any other system | false (the default) | no `cellat` on its complete level grid, so the generic search answers there too |

The family declaration for `AbstractQuadFaceGridSystem` is not a shortcut: the
family contract already obliges every subtype to supply `cellat`, alongside
`levels`, `cell_boundary`, `cell_centroid`, `node_extent` and `one_ring`.

Every system this package ships therefore locates directly, and every partial
grid over one of them composes. The false side is reachable and is exercised:
the substrate suite's mock systems implement no `cellat` at all, and the
conformance library's `GenericFallbackSystem` hides a real system's fast paths
by type.

### The semantics that changed

The composed locator answers what the complete level answers, restricted to the
subset. A point the complete level assigns to a cell outside the subset is
outside this grid's coverage **even when it lies on the boundary of a member
cell** — the two cells share that boundary and the complete level's tie rule
decides which of them owns it. The tree fallback tests member cells only, and
`point_in_cell` reports a point on a boundary as inside, so it would hand back
the member. The set of points where the two differ has measure zero.

It is constructible deterministically, and the suite pins it on every registered
system: take a member cell's boundary vertex that the complete level awards to a
non-member neighbour, and the composed locator answers `nothing` where the
generic search answers the member cell.

## 3. Where the time goes now

One destination tile of 823,543 cells, GLO-90 N46 E010 into IGeo7 level 12,
each phase in its own compiled function, one shot after a warm-up on a
1,024-cell slice. The "before" column is the companion note's section 3.

| phase | ns per destination cell, before | after | bytes, before | after |
|---|---|---|---|---|
| `CellVector` decode alone | 60.3 | 59.9 | 0 | 0 |
| `sites[i]` — decode + IGeo7 centroid | 128.3 | 128.6 | 0 | 0 |
| `weightsat!` — site read + `cellat` | 6,373.5 | **335.7** | 285,488,416 | **0** |
| `cellat` alone, sites precomputed | 6,085.2 | **166.1** | 285,488,416 | **0** |
| `chunkat` over the located indices | 1.1 | 1.0 | 0 | 0 |
| whole `tileweights` | 6,290.1 | **356.2** | 401,288,656 | 126,204,080 |
| residual | 75.6 | 60.5 | | |
| cold read of the same tile | 6,621.0 | **669.0** | 511,855,280 | 219,258,032 |
| warm read of the same tile | 16.8 | 8.4 | 40,408,096 | 34,032,464 |

`cellat` is **36.6x faster and allocates nothing**. Weights fall from 95.0% of a
cold read to 53.2% of it, and the read as a whole from 6,621 to 669 ns per
destination cell.

The companion note's section 11 predicted 640 ns per destination cell for
exactly this change, by adding up its own measured parts. The measurement is
669 ns.

## 4. The two routes, now within 58 ns of each other

The same 100,000 sample sites of the same tile:

| route | ns per call, before | after | bytes, after |
|---|---|---|---|
| `cellat(levelgrid(sys, 1), p)` then `localindex(grid, c)` | 102.2 | 100.2 | 0 |
| `cellat(space, p)` on the `DGGSpace` over that grid | 6,366.6 | **158.3** | **0** |
| `treeify(grid)` on its own | 47.1 | 55.5 | 352 |

**62x becomes 2x.** The remaining 58 ns is `DGGSpace`'s own translation: its
`cellat` calls the grid's, which already resolves membership with a binary
search over `ids`, and then calls `localindex` on the same grid again to turn
the cell into a local index. That is one redundant search per call and is the
next thing to remove on this path; it is 8.7% of a cold read.

## 5. `NearestCell` is no longer slower than `Conservative`

The same destination tile against the same source chunk, one session, 20,000
destination cells each:

| method, `buildweights!` on one (tile, chunk) pair | ns per cell, before | after | bytes, after |
|---|---|---|---|
| `Conservative()` | 3,186.5 | 3,119.8 | 87,754,448 |
| `NearestCell()` | 6,758.9 | **364.7** | 1,372,608 |

`Conservative` is unchanged, as it must be: it never called `cellat`. Nearest
was 2.1x the clip and is now **8.6x faster than it**, which is the ordering a
one-entry stencil against a polygon intersection should have had all along.

## 6. The cliff is gone from the nearest path

Locating a point in the source, over up to 200,000 destination sites per shape.
The "before" column is the companion note's section 6, which measured the same
call on the same six shapes.

| source tiles | chunks | tree `treeify` gives | us per `cellat`, before | after |
|---|---|---|---|---|
| N46 E010 | 1 | `MemoBlockCursor` | 7.0 | **0.160** |
| N46 E010-E011 | 2 | `MemoBlockCursor` | 8.0 | **0.240** |
| N46 E010-E012 | 3 | `MemoBlockCursor` | 7.0 | **0.306** |
| N46 E010-E013 | 4 | `MemoBlockCursor` | 7.0 | **0.375** |
| N46-N47 E010 | 2 | `HierarchicalGridCursor` | **662,300** | **0.237** |
| N46-N47 E010-E011 | 4 | `HierarchicalGridCursor` | **731,700** | **0.375** |

Every shape allocates nothing per call. Which cursor `treeify` gives no longer
appears in the numbers at all, because no location consults a tree: the two
`HierarchicalGridCursor` shapes cost what the `MemoBlockCursor` shapes of the
same chunk count cost, to within 2%. The residual growth with chunk count is the
membership search over a longer id vector and the wider destination box, not the
tree.

### What the tree still costs on that shape

The cursor split is untouched, so the shape still gets the generic cursor and
the paths that do consult a tree still pay for it. Measured on the 2x2 source
(N46-N47 x E010-E011) against its own destination at level 12:

| what | cost |
|---|---|
| `treeify(source grid)`, first call | 160 B, below timer resolution |
| the same call again, memoized | 0 B, below timer resolution |
| `candidatechunks!` per destination chunk, 609 chunks | 4.5 us each, 2.7 ms total, 96 B total |
| `Conservative` `buildweights!`, one destination chunk x one source chunk, 2,000 cells | **449,867 ns per destination cell**, 83,974,944 B |
| `NearestCell` `buildweights!`, the same pair and cells | 477 ns per destination cell, 144 B |

Building the tree is free and the chunk-candidate search is 2.7 ms for the whole
destination, so neither is a reason to change the cursor. What is left is the
per-cell tree descent that `Conservative` does: 449,867 ns per destination cell
on this shape against 3,120 ns on a one-tile source, a factor of 144. The cliff
did not disappear; it stopped being on the nearest path. Anyone regridding this
shape conservatively still meets it.

## 7. The arms

`benchmark/copdem_nearest.jl` unchanged, one fresh process per arm, machine
otherwise idle, `-t auto`. "Before" for the first arm is a run of the same
command on the unchanged library in this same session, so the two are comparable
without carrying a machine state across days; the companion note's own figure is
given beside it.

### GLO-90 into IGeo7 level 12

| arm | cold sweep | cold cells/s | warm sweep | warm cells/s | destination arm | peak RSS |
|---|---|---|---|---|---|---|
| 1 tile, before (this session) | 15.172 s | 152,507 | 0.018 s | 130,056,462 | 1.152 s | 1.37 GiB |
| 1 tile, before (companion note) | 15.809 s | 146,364 | 0.039 s | 59,300,001 | 1.163 s | 1.27 GiB |
| **1 tile, after** | **1.040 s** | **2,224,097** | 0.018 s | 128,300,137 | 1.135 s | 1.45 GiB |
| 2 tiles, before (companion note) | 31.216 s | 148,153 | 0.048 s | 96,922,022 | 1.543 s | 1.48 GiB |
| **2 tiles, after** | **2.473 s** | **1,869,897** | 0.085 s | 54,540,585 | 1.509 s | 1.57 GiB |
| 2x2 block, before | **did not finish** | | | | | |
| **2x2 block, after** | **6.040 s** | **1,516,684** | 6.103 s | 1,501,012 | 2.026 s | 1.98 GiB |

The one-tile arm is **14.6x faster end to end** against the same-session
baseline, and the two-tile arm 12.6x against the companion note's.

### GLO-30 into IGeo7 level 13

| arm | cold sweep | cold cells/s | warm sweep | warm cells/s | destination arm | peak RSS |
|---|---|---|---|---|---|---|
| 1 tile, before | 119.084 s | 135,886 | 118.608 s | 136,432 | 5.678 s | 2.05 GiB |
| **1 tile, after** | **8.182 s** | **1,977,639** | 7.661 s | 2,112,362 | 5.675 s | 2.38 GiB |

**14.6x**, the same factor as GLO-90, which is what a change that removes a
fixed per-cell cost should give: the companion note showed a nine-fold finer
source cost only 13%, because the locator's cost was the tree and the polygon,
not the pixel count.

### The 2x2 block

The box `10,12,46,48` is four tiles across two latitude rows — the shape the
companion note's section 14 listed as unreachable, at about 75 days for a
9.3-million-cell destination. It now reads 9,160,874 destination cells in
21 tiles in 6.040 seconds, and its 1,000-cell check is 998 mapped, 2 outside
coverage, 0 wrong.

### Correctness

`check=1000` per arm: the destination cell's centroid is located in the source
space, the pixel it names is read directly, and the regridded value must equal
it bit for bit.

| arm | cells checked | mapped | outside coverage | wrong |
|---|---|---|---|---|
| GLO-90, 1 tile | 1,000 | 996 | 4 | **0** |
| GLO-90, 2 tiles | 1,000 | 997 | 3 | **0** |
| GLO-90, 2x2 block | 1,000 | 998 | 2 | **0** |
| GLO-30, level 13, 1 tile | 1,000 | 999 | 1 | **0** |

The one-tile GLO-90 arm reports the same 996 mapped and 4 outside coverage as
the before run, over the same cells, and the same value range
[227.63, 3867.07] m and the same 2,309,219 finite of 2,313,802.

## 8. The profile

Two `Profile` listings over the same four destination tiles (526,893 cells) of
the two-tile GLO-30 strip, the same shape and delay the companion note's
section 9 used.

**Cold read, weights built here.** 135 samples reached this package, of 2,480
over all threads — against 1,273 of 22,032 before. At an unchanged sampling
delay that is a **9.4x drop in the time this package is on the stack** for the
same work.

The inclusive listing no longer contains a `locate.jl` frame at all. What it
contains, below the 100% read frames:

| inclusive | frame |
|---|---|
| 9.63% | `cell_centroid` through `partial_grid.jl` and `level_grid.jl` |
| 9.63% | `getindex` on the source's tile ids |
| 6.67% | `_readdestination!` |
| 5.93% | `cellindex` at `level_grid.jl:70` |

Top self time is 3.70% (`getindex`, and GC marking), and the Copernicus closed
form itself is 2.22% at `system.jl:307`. Before, `cellat` at `locate.jl:14` was
84.68% inclusive and the top self frames were `sind` at 12.25% and `cosd` at
8.32%. Nothing dominates the read now.

**Warm read, weights cached**: 23 samples of 368, against 11 of 208 before. Too
cheap to profile at this delay either way.

## 9. What a full globe would take

The production GLO-90 destination is 66,175 IGeo7 level-5 chunks of
823,543 cells, 5.45e10 cells; at level 13 the same land is 3.81e11. The
extrapolation assumes the measured rate holds, which sections 6 and 7 establish
over four-fold growth in box size and across the tile shape that used to decide
everything, and assumes no scheduling loss.

| run | destination cells | one core, before | one core, after | at 40 workers, after |
|---|---|---|---|---|
| GLO-90 to level 12 | 5.45e10 | 145,000 cells/s -> 104 hours | 1.52-2.22e6 cells/s -> **6.8-10.0 hours** | 0.17-0.25 hours |
| GLO-30 to level 13 | 3.81e11 | 130,000 cells/s -> 815 hours | 1.98e6 cells/s -> **53.5 hours** | 1.3 hours |

The lower end of the GLO-90 range is the 2x2 arm's rate and the upper end the
one-tile arm's; the difference is not the locator but the weight budget, below.

The companion note said neither was reachable, because a production destination
chunk pairs with a multi-row tile set and every such chunk paid 0.7 s per
destination cell. Section 6 is why that sentence no longer holds: the multi-row
shapes locate in 0.24 us, indistinguishable from the single-row ones.

Three walls remain, and none of them is location.

- **The weight budget.** The 2x2 arm holds 530,198,752 B of weights against a
  536,870,912 B budget and keeps 15 of its 21 destination tiles, so its warm
  sweep costs what its cold sweep did — 6.103 s against 6.040 s. The GLO-30 arm
  does the same, 7.661 s warm against 8.182 s cold, keeping 4 tiles of 33. The
  one-tile and two-tile GLO-90 arms fit entirely and sweep warm 60x and 29x
  faster than cold. This is the companion note's section 7 boundary, unchanged,
  and it is now the largest remaining lever on a repeated sweep.
- **Source re-copying.** The GLO-30 arm read its one source chunk 33 times and
  copied 427,680,000 pixels for 16,181,892 destination cells. Nothing caches a
  decoded source across destination tiles, so a chunk is re-copied for every
  tile that names it.
- **Data.** The GLO-90 globe is roughly 132 GB of GeoTIFF and GLO-30 roughly
  1.1 TB. Neither fits this machine, which is why every arm here is a box.

## 10. Metadata

- Julia 1.12.6, macOS 26.5, M-series, `-t auto` = 8 threads and 8 GC threads of
  12 CPUs, `pmset -g` reports `powermode 2`, 24 GiB of RAM.
- Every timed arm ran with the machine otherwise idle, checked before each arm.
- One shot per arm after a compilation warm-up, as in the companion note. The
  before/after pairs for the first GLO-90 arm are two runs of the same command
  in one session on the same machine state, differing only in the library; the
  other before columns are the companion note's own figures, taken on the same
  machine on the same day.
- Source tiles are the same local GeoTIFFs the companion note used, passed
  through the benchmark's `data=` argument. Nothing was downloaded for this
  note and nothing was written into the package.

## 11. What this does not show

- **Nothing about the `Conservative` path.** Section 6 measures that the tree
  descent it does per destination cell still costs 449,867 ns on a source that
  crosses a latitude row. Fixing that is the cursor's problem, not the
  locator's, and nothing here changes the cursor.
- **No full-globe or production-shaped run.** Section 9 is arithmetic on
  measured rates over boxes of up to four tiles.
- **Nothing about accuracy beyond the bit-identity check.** The regridded values
  are identical to the before run on the arm where both were measured; no arm
  compares against an independent implementation.
- **Nothing about the boundary-tie change at scale.** The set of points where
  the composed locator and the tree fallback differ has measure zero, and the
  suite pins one such point per system by construction. No arm samples enough
  boundary points to say anything statistical about them, and none needs to.
