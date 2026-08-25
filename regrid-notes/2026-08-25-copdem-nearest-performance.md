# Nearest regridding of Copernicus DEM tiles to IGeo7

> Superseded for the locator cost by `2026-08-25-partialgrid-cellat.md`.

- Date: 2026-08-25
- Benchmark: `benchmark/copdem_nearest.jl`
- Nothing in the library changed. This is what `NearestCell` costs on real
  Copernicus DEM data, and where that cost is.

`benchmark/point_tile_baseline.jl` priced the two nearest routes on a toy
raster, where `cellat` is two binary searches on prepared edge vectors. This
note prices the same method on the workload that toy stands in for: real
1x1-degree Copernicus DEM GeoTIFFs into IGeo7, where the source is a
`PartialGrid` over a `CopernicusDEMSystem` and the destination is a `PartialGrid`
over IGeo7 cells.

The headline is not a throughput number. It is that `cellat` on the source
`PartialGrid` costs 6,367 ns and 347 bytes per call where the same answer is
available in 102 ns and no bytes, and that the same call costs about 0.7 s
whenever the source tiles do not form one contiguous run. Sections 3 and 4
establish both, section 5 shows what that does to a regrid, and section 14 says
which measurements the second effect makes unreachable.

## 1. The workload

The source is built exactly as `scripts/copdem_production.jl` builds it — the
benchmark includes that file rather than reimplementing it. One `TiledDEM` over
a `PartialGrid` of the listed tiles' pixels, chunk `k` being tile `k`, decoded
through `TileBuilder` from local GeoTIFFs and held in a `StripedLRUCache` with
one slot per tile. The source space is `DGGSpace(grid; chunklevel = 0)`, so a
source chunk is exactly a Copernicus tile.

The destination is the IGeo7 cells a box covers at the destination level —
`query(sys, MultiOrderCoverage(box); level)` resolved to a `CellVector` and a
`PartialGrid` — tiled at `level - 7`, so a full destination tile is
7^7 = 823,543 cells, the chunk the production run writes. Taking whole ancestor
subtrees the way a subzone store must would add cells outside the box that no
source can reach, and those would dilute every per-cell number here.

Method `NearestCell()`, missing policy `Weighted(0.5)`, `lazy = true`, the
default 2 GiB budget, `-t auto`.

## 2. The level the bare-system rule picks

`to = IGeo7System()` names no cells until a level is chosen, and as a
destination it takes `levelfor(sys, src_space)` — the level whose median cell
area is closest to the source's.

| product | rows per degree | median pixel over the box | IGeo7 level | IGeo7 cell | cells globally |
|---|---|---|---|---|---|
| GLO-90 | 1200 | 76.88 m | **12** | 60.71 m | 138,412,872,012 |
| GLO-30 | 3600 | 25.63 m | **13** | 22.94 m | 968,890,104,072 |

`scripts/copdem_production.jl` fixes `level = 12` for GLO-90, which is the level
the rule picks anyway, so for GLO-90 the two agree and there is one arm. For
GLO-30 they differ, and both are measured: level 13 is the rule's, level 12 is
the production setting applied to 30 m data.

A bare system as destination names the WHOLE globe at that level, which is
1.4e11 cells at level 12 and 9.7e11 at level 13. No arm here sweeps that; every
arm names the IGeo7 cells of its box at the rule's level, and section 11
extrapolates.

## 3. Where the time goes: `cellat` on the source

One destination tile of 823,543 cells, GLO-90 N46 E010 into IGeo7 level 12,
every phase in its own compiled function so a top-level loop over globals cannot
price boxing instead of work. Each phase is one shot after a warm-up on a
1,024-cell slice.

| phase | s | ns per destination cell | bytes |
|---|---|---|---|
| `CellVector` decode alone | 0.050 | 60.3 | 0 |
| `sites[i]` — decode + IGeo7 centroid | 0.106 | 128.3 | 0 |
| `weightsat!` — site read + `cellat` | 5.249 | 6,373.5 | 285,488,416 |
| `cellat` alone, sites precomputed | 5.011 | **6,085.2** | 285,488,416 |
| `chunkat` over the located indices | 0.001 | 1.1 | 0 |
| **whole `tileweights`** | **5.180** | **6,290.1** | 401,288,656 |
| residual (`tileweights` less the three above) | 0.062 | 75.6 | |
| cold read of the same tile | 5.453 | 6,621.0 | 511,855,280 |
| warm read of the same tile | 0.014 | 16.8 | 40,408,096 |

The phases sum: 128.3 + 6,085.2 + 1.1 + 75.6 = 6,290.2 ns of `tileweights`, and
`tileweights` is 95.0% of the cold read; the remaining 331 ns per cell is the
GeoTIFF decode, the source load and the weight application together. A warm read
of the same tile is 16.8 ns per cell, so applying the weights is 0.25% of a cold
read and every question about this workload is a question about `cellat`.

Allocation follows the same line: 622 bytes per destination cell over the cold
read, 487 of them in `tileweights` and 347 in `cellat` alone. A nearest stencil
is one entry of one; none of those bytes is the stencil.

The process ran 46.55 s of wall against 52.55 s of user time, a ratio of 1.13,
and `tileweights` is a serial loop over the tile's cells. The sweep is
single-threaded plus GC threads, so per-core throughput and aggregate throughput
are the same number here.

## 4. The same answer costs 102 ns

`CopernicusDEMSystem` defines a closed-form `cellat(::LevelGrid, p)` —
`src/systems/CopernicusDEM/system.jl` — that is a latitude-band lookup and two
floors. `DGGSpace`'s `cellat` reduces whatever the grid answers through
`localindex` anyway. Composing those two directly, over 100,000 sample sites of
the same tile:

| route | ns per call | bytes per call |
|---|---|---|
| `cellat(levelgrid(sys, 1), p)` then `localindex(grid, c)` | **102.2** | **0** |
| `cellat(space, p)` on the `DGGSpace` over that grid, as shipped | **6,366.6** | **346.7** |
| `treeify(grid)` on its own | 47.1 | 352 |

**62x, for the same answer.** `treeify` is not the cost — it is memoized and
answers in 47 ns. The cost is what `cellat(grid::AbstractGrid, p)` in
`src/fallbacks/locate.jl` does with the tree it gets: an `STI.query` over the
tile's 1,440,000 pixel cells, a `sort!` of the surviving candidates into a fresh
vector, and then an exact `point_in_cell` against a `cell_boundary` polygon
built for each candidate. That is the generic locator for a grid with no
closed-form method, and `PartialGrid` defines none, so a partial grid over a
system that HAS one still takes it.

## 5. `NearestCell` is slower than `Conservative` on this data

On the same destination tile against the same source chunk, in one session,
20,000 destination cells each:

| method, `buildweights!` on one (tile, chunk) pair | ns per destination cell | bytes |
|---|---|---|
| `Conservative()` | 3,186.5 | 84,690,640 |
| `NearestCell()` | 6,758.9 | 8,307,376 |

`Conservative` clips polygons and `NearestCell` looks up one cell, and the
lookup is **2.1x the clip**. `Conservative` never calls `cellat`; it goes
through the source's spatial tree once per destination cell and pays the tree
descent once, where the nearest locator pays a tree query, a sort and a polygon
test for a point it could have computed arithmetically. This inversion is the
clearest statement of the defect in section 4: nothing about nearest is
expensive, and the method is nonetheless the slow one.

## 6. The cliff: a source that is not one contiguous run

`treeify(::PartialGrid{<:CopernicusDEMSystem})` gives a `BlockCursor` only when
the grid's cells are one rectangle of the Copernicus lattice held as one
contiguous id run — documented at `src/systems/CopernicusDEM/cursor.jl` — and a
`HierarchicalGridCursor` otherwise. Tiles in ONE latitude row are one such run;
tiles spanning two rows are not. Because `cellat` is built on the tree, the
locator inherits that split.

200 `cellat` calls per shape on the fast shapes, 3 on the slow ones:

| source tiles | chunks | tree `treeify` gives | ms per `cellat` |
|---|---|---|---|
| N46 E010 | 1 | `MemoBlockCursor` | 0.007 |
| N46 E010–E011 | 2 | `MemoBlockCursor` | 0.008 |
| N46 E010–E012 | 3 | `MemoBlockCursor` | 0.007 |
| N46 E010–E013 | 4 | `MemoBlockCursor` | 0.007 |
| N46–N47 E010 | 2 | `HierarchicalGridCursor` | **662.3** |
| N46–N47 E010–E011 | 4 | `HierarchicalGridCursor` | **731.7** |

Four orders of magnitude, decided by whether the tile set crosses a latitude
row. A 996-cell destination against the 2x2 block did not finish a single read
in 30 seconds of profiling, and 95.6% of the samples that reached the package
were under `tileweights` → `weightsat!`.

This is not a corner case. Every production destination is a chunk of a global
run, and the source it pairs with is whatever tiles its cap reaches — a 55 km
IGeo7 level-5 cell straddles a latitude row wherever it crosses a whole degree
of latitude. A 2x2-degree box is the smallest interesting one and it is already
on the wrong side.

## 7. Throughput, GLO-90 and GLO-30

Only east-west strips within one latitude row are reachable, for the reason in
section 6, so the scaling arms grow the box along longitude. Boxes are 1x1
degree tiles in latitude row 46, growing east. Each arm is a fresh process;
`cold sweep` reads every destination tile once in order on a plan that has never
been read, `warm sweep` repeats it against the plan that sweep left behind.

### GLO-90 into IGeo7 level 12

| source tiles | destination cells | destination tiles | candidates per tile (med/max) | destination arm | cold sweep | cold cells/s | warm sweep | warm cells/s | peak RSS |
|---|---|---|---|---|---|---|---|---|---|
| 1 (E010) | 2,313,802 | 8 | 1.0 / 1 | 1.163 s | 15.809 s | 146,364 | 0.039 s | 59,300,001 | 1.27 GiB |
| 2 (E010-E011) | 4,624,703 | 13 | 2.0 / 2 | 1.543 s | 31.216 s | 148,153 | 0.048 s | 96,922,022 | 1.48 GiB |
| 3 (E010-E012) | 6,935,841 | 18 | 2.0 / 3 | 1.717 s | 48.295 s | 143,613 | 0.079 s | 87,629,025 | 1.65 GiB |
| 4 (E010-E013) | 9,246,816 | 21 | 2.0 / 3 | 2.070 s | 64.063 s | 144,338 | **63.974 s** | **144,540** | 1.67 GiB |

Reads and weights for the same four arms:

| source tiles | source chunk reads | pixels copied | GeoTIFF decodes | resident weight bytes | budget | destination tiles cached |
|---|---|---|---|---|---|---|
| 1 | 8 | 11,520,000 | 1 | 147,619,072 | 536,870,912 | 8 of 8 |
| 2 | 16 | 23,040,000 | 2 | 309,841,656 | 536,870,912 | 13 of 13 |
| 3 | 25 | 36,000,000 | 3 | 486,996,576 | 536,870,912 | 18 of 18 |
| 4 | 33 | 47,520,000 | 4 | 525,094,312 | 536,870,912 | **16 of 21** |

`Conservative()` over the identical one-tile destination: 6.284 s and **368,199
cells/s**, against nearest's 15.809 s and 146,364 — `Conservative` is **2.52x
faster end to end**, the whole-pipeline form of section 5.

### GLO-30 into IGeo7

| destination | source tiles | destination cells | destination tiles | candidates (med/max) | destination arm | cold sweep | cold cells/s | warm sweep | warm cells/s | peak RSS |
|---|---|---|---|---|---|---|---|---|---|---|
| level 13 (the rule's) | 1 | 16,181,892 | 33 | 1.0 / 1 | 5.678 s | 119.084 s | 135,886 | 118.608 s | 136,432 | 2.05 GiB |
| level 13 | 2 | 32,351,955 | 59 | 1.0 / 2 | 7.058 s | 251.166 s | 128,807 | 246.933 s | 131,015 | 2.15 GiB |
| level 13 | 3 | 48,523,116 | 81 | 2.0 / 3 | 7.652 s | 371.226 s | 130,710 | 369.284 s | 131,398 | 2.37 GiB |
| level 12 (production's) | 1 | 2,313,802 | 8 | 1.0 / 1 | 1.147 s | 18.218 s | 127,009 | 17.876 s | 129,438 | 1.97 GiB |

The level-12 arm is the same destination as every GLO-90 arm above, fed by a
source nine times finer. It runs at 127,009 cells/s where GLO-90 runs at
146,364, so a nine-fold finer source costs 13% — `cellat`'s cost is the tree
query and the polygon test, not the number of pixels it discriminates between.

**Every GLO-30 arm sweeps warm at cold speed**, and the reason is in the weight
bytes: 451,737,104 B for the level-12 destination against GLO-90's 147,619,072
for the identical destination cells. A weight block carries a reference vector
sized by the SOURCE chunk it names, so a 12,960,000-pixel source tile makes a
tile's weights three times the bytes a 1,440,000-pixel one does. Only 4 of the 8
destination tiles then fit the 512 MiB weight budget, and an in-order sweep
against an LRU cache holding less than the working set reuses nothing. The
pair-route GLO-30 arm says this outright: it placed 16,181,892 destination cells
on the cold sweep and **16,181,892 again on the warm one**, so not one weight
survived to be reused.

The same boundary produces the GLO-90 four-tile warm sweep: 21 destination tiles
needing 525,094,312 B against a 536,870,912 B budget keeps 16, and the warm
sweep costs what the cold one did (63.974 s against 64.063 s). Three tiles fit
and sweep warm in 0.079 s, eight hundred times faster. The cliff is at the
budget, not at the box.

Source reads amplify with the destination tiling rather than the box: the
level-13 three-tile arm read its 3 source chunks 98 times and copied
1,270,080,000 pixels — 5.1 GB of memcpy for 128 MB of GeoTIFF decoded three
times. Data is not cached across reads, so a source chunk is re-copied for every
destination tile that names it.

Throughput is flat across box size — 146,364, 148,153, 143,613 and 144,338
cells/s over a four-fold growth on GLO-90, and 135,886, 128,807 and 130,710 over
a three-fold growth on GLO-30 — so the cost is per destination cell and nothing
about the box is superlinear. That flatness is what makes section 3's
attribution the whole story: at 145,000 cells/s a destination cell costs 6.9 us,
and section 3 accounts for 6.6 us of it in one `cellat`.

What does grow is destination-space construction: 1.163, 1.543, 1.717 and
2.070 s for 2.3 to 9.2 million cells at level 12, allocating 760 MB to 1.70 GB,
and 5.678 to 7.652 s for 16.2 to 48.5 million at level 13. That is the
`MultiOrderCoverage` query, the `CellVector` and the `PartialGrid` over it. It
is sublinear in cells, because the coverage set compacts whole subtrees and only
the boundary refines. It is not part of any regrid, but a caller pays it.

`plan_regrid` never exceeded 8.5 MB or a measurable fraction of a second on any
arm. The relation is not a cost on this workload.

## 8. The two routes on real data

`PairNearest` is `NearestCell`'s own `buildweights!` with no `sampler`, so a plan
around it builds one `(destination tile, source chunk)` pair at a time. Nothing
but the routing differs, and each pair of arms ran on the same data.

| product | source tiles | candidates (med/max) | pair locations | tile locations | location ratio | pair cold sweep | tile cold sweep | time ratio |
|---|---|---|---|---|---|---|---|---|
| GLO-90 | 1 | 1.0 / 1 | 2,313,802 | 2,313,802 | 1.00x | 15.111 s | 15.809 s | 0.96x |
| GLO-90 | 2 | 2.0 / 2 | 8,561,833 | 4,624,703 | **1.85x** | 56.876 s | 31.216 s | **1.82x** |
| GLO-30 | 1 | 1.0 / 1 | 16,181,892 | 16,181,892 | 1.00x | 117.745 s | 119.084 s | 0.99x |
| GLO-30 | 2 | 1.0 / 2 | 50,372,164 | 32,351,955 | **1.56x** | 377.992 s | 251.166 s | **1.51x** |

At one source tile every destination tile has exactly one candidate, the two
routes do identical work, and the pair route is 1-4% ahead — the tile route's
`chunkat` per entry and its manifest, and nothing else. At two source tiles the
pair route performs 1.85 and 1.56 locations per destination cell where the tile
route performs one, and it costs 1.82x and 1.51x the time. **The location ratio
and the time ratio agree to within 3% in both cases.** That is the statement
`benchmark/point_tile_baseline.jl` makes on the toy raster (4.00x locations, 3.13x time) made again on real data:
what the tile route removes is the repeated location and only that, and here,
where the location is 92% of the read, the two ratios coincide instead of
diverging.

The multiplicity is small because a destination tile is SMALLER than a source
tile — an IGeo7 level-5 cell is about 55 km against a 1x1-degree tile's 111 km
by 77 km — so most destination tiles sit inside one source tile and only the
straddlers have two. A destination tiling coarser than the source, or a source
chunked finer, would raise it.

## 9. The profile

Two `Profile` listings over the same four destination tiles (526,893 cells) of
the two-tile GLO-30 strip, taken on a plan of its own after every source tile has
been decoded and made resident, so neither is a profile of GDAL. Self time is the
frame a sample was taken in; the inclusive table below it credits a frame once
per sample it appears in anywhere, so those nest rather than sum. Only samples
whose backtrace reaches this package are counted — under `-t auto` the idle
threads park in `__psynch_cvwait` and are 94% of all samples.

The east-west strip stands in for the 2x2-degree box: section 6's cliff makes
that box unprofilable, since a 996-cell destination against it did not finish a
single read in 30 seconds of sampling and every sample would be the same frame.

**Cold read, weights built here** (1,273 samples reaching the package, of 22,032):

| inclusive | frame |
|---|---|
| 98.90% | `readblock!` at `lazy.jl:422` |
| 97.96% | `_readdestination!` → `tilefor` → `gettile!` |
| 94.27% | `tileweights` at `plans.jl:1000` |
| 92.54% | `weightsat!` at `interpolation.jl:76` |
| 89.95% | `cellat` at `regridding.jl:158` (the `DGGSpace` method) |
| **84.68%** | **`cellat` at `locate.jl:14` (the generic fallback)** |
| 26.16% | `node_extent` at `cursor.jl:247` |
| 22.94% | `child_indices_extents` at `cursor_memo.jl:101` |
| 21.45% | `_box_cap` at `cursor.jl:200` |
| 20.42% | `_leafcell` at `cursor.jl:289` |
| 6.91% | `cap_contains` at `caps.jl:25` |

Top self time: `sind` 12.25%, `cosd` 8.32% across two entries, `_taskmemo`
5.89% across three, `node_extent` 3.46%, `getchild` 1.65%. Trigonometry and
node-extent bookkeeping — the spherical caps the tree query tests against —
with no frame belonging to a nearest stencil, because a nearest stencil is one
entry of one.

The chain reads plainly: 92.5% of the cold read is `weightsat!`, 90% is the
`cellat` it calls, and **84.7% is inside the generic locator alone** — the tree
query, the candidate sort and the exact polygon test that section 4 shows a
closed-form lookup would replace with 102 ns.

**Warm read, weights cached** (11 samples reaching the package, of 208): the
whole read of 526,893 cells produced eleven samples. `_applygroup!` and
`anyinvalid` at 54.55%, `_readsource!` and the `TiledDEM` `readblock!` at
27.27%, `applyblock!` at 9.09%. Applying a nearest tile's weights and copying
its source is too cheap to profile at this delay, which is the same statement as
the 16.8 ns per cell in section 3.

## 10. Correctness on real data

`check=1000` verifies random destination cells against the source read directly:
the sample site is the destination cell's centroid, `cellat` on the source space
names the pixel, and the regridded value must equal that pixel BIT for bit
(`isequal`, so a NaN pixel must give a NaN cell). Where `cellat` finds nothing
the cell must be missing.

A nearest stencil is one entry of exactly 1.0 and `Weighted` divides by exactly
one, so bit identity is the right assertion and approximate equality would be a
weaker test.

| arm | cells checked | mapped | outside coverage | wrong |
|---|---|---|---|---|
| GLO-90, 1 tile | 1,000 | 996 | 4 | **0** |
| GLO-90, 2 tiles | 1,000 | 997 | 3 | **0** |
| GLO-90, 3 tiles | 1,000 | 998 | 2 | **0** |
| GLO-90, 4 tiles | 1,000 | 996 | 4 | **0** |
| GLO-90, 1 tile, pair route | 1,000 | 996 | 4 | **0** |
| GLO-90, 2 tiles, pair route | 1,000 | 997 | 3 | **0** |
| GLO-30, level 13, 1 tile | 1,000 | 999 | 1 | **0** |
| GLO-30, level 13, 2 tiles | 1,000 | 997 | 3 | **0** |
| GLO-30, level 13, 3 tiles | 1,000 | 997 | 3 | **0** |
| GLO-30, level 13, 1 tile, pair route | 1,000 | 999 | 1 | **0** |
| GLO-30, level 13, 2 tiles, pair route | 1,000 | 997 | 3 | **0** |
| GLO-30, level 12 | 1,000 | 996 | 4 | **0** |

Twelve thousand cells, none wrong. The cells outside coverage are the box's own
boundary: the coverage query keeps every IGeo7 cell the box touches, so a cell
whose centroid falls just outside a tile's 1-degree extent has no source, and
every one of them came back missing under `Weighted(0.5)` as the policy
requires. Finite fractions over the whole destination run 99.7% to 99.9%, and
the elevations span 22.5 m to 3,874.1 m, which is the Alps.

## 11. What a full globe would take

The production GLO-90 destination is 66,175 IGeo7 level-5 chunks of 823,543
cells, 5.45e10 cells; at level 13 the same land is seven times that, 3.81e11.
Extrapolating the measured per-cell cost assumes the rate stays flat, which
section 7 establishes over a four-fold growth in box size and no further, and
that a full run pays the same per-cell cost per core with no scheduling loss.

| run | destination cells | at the measured rate, one core | at 40 workers, perfect scaling |
|---|---|---|---|
| GLO-90 to level 12 | 5.45e10 | 145,000 cells/s -> **104 hours** | 2.6 hours |
| GLO-30 to level 13 | 3.81e11 | 130,000 cells/s -> **815 hours** | 20 hours |

**Neither is reachable today**, because both rates belong to a source that is one
contiguous run of tiles. A production destination chunk pairs with whatever
tiles its cap reaches, and a 55 km level-5 cell crosses a whole degree of
latitude nearly everywhere, so nearly every chunk gets the `HierarchicalGridCursor`
path and 0.7 s per destination cell. At that rate GLO-90 alone is 3.8e10 seconds
of core time, which is about 1,200 years.

With `cellat` composed from the closed form section 4 measures, a destination
cell would cost 128 ns for its site, 102 ns for the lookup, 1 ns for `chunkat`
and 76 ns for the row and accumulator work — 307 ns of weights — plus the 331 ns
of source read and application the cold read shows, so about 640 ns, or
1.57 million cells/s. That puts GLO-90 at roughly 9.6 hours on one core and
GLO-30 at 67 hours, and it removes the cliff along with the cost, because the
closed form does not consult a tree at all. Those two numbers are arithmetic on
the measured parts, not a measurement of a changed library.

Data is the other wall. The GLO-90 globe is 26,475 tiles of about 5 MB, roughly
132 GB; GLO-30 is the same 26,475 tiles at about 43 MB, roughly 1.1 TB. Neither
fits the 23 GiB free on this machine, which is why every arm here is a box.


## 12. Data footprint

Downloaded into the directory the `data=` argument named, outside the
repository, and safe to delete:

    copdem-data/    280 MB total, 288,835,243 bytes of GeoTIFF

- `tiles-glo90/` — 77 MB, 15 tiles downloaded (N46–N49 x E010–E013) plus one
  symlink to a tile already in `$RASTERDATASOURCES_PATH`.
- `tiles-glo30/` — 203 MB, 5 tiles downloaded (E012 for N46–N47, and N48 x
  E010–E012) plus four symlinks to tiles already in `$RASTERDATASOURCES_PATH`.

The four GLO-30 and one GLO-90 tiles the user already had are symlinked rather
than copied, so they cost nothing and are untouched. Nothing was written into
`$RASTERDATASOURCES_PATH`. The 6 GiB download cap was never approached: time,
not disk, bounded the arms, and section 6 bounded them harder than either.

`/Users/anshul/Downloads/Copernicus_DSM_10_N46_00_E010_00_DEM.tif` was read but
not used. It is 3601x3601, the ESA distribution shape, where the AWS COG for the
same tile is 3600x3600 and is what `CopernicusDEMSystem(30)` and
`scripts/copdem_production.jl`'s `readtile` require — that reader raises on a
raster size it does not expect. The COG was fetched instead.


## 13. Metadata

- Julia 1.12.6, macOS 26.5, M-series, `-t auto` = 8 threads and 8 GC threads of
  12 CPUs, `pmset -g` reports `powermode 2`, 24 GiB of RAM.
- Every timed arm ran with the machine otherwise idle: the driver waits until no
  other Julia process matching a benchmark or scratch path is running, and each
  arm's log records the count of such processes before and after it.
- Statistic: one shot per arm, after a compilation warm-up. Ratios inside one
  process are the portable part; the absolutes belong to this machine state and
  this power mode.
- Source reads are counted by a `DiskArrays` wrapper that records every
  `readblock!` range against the tile it lands in; GeoTIFF decodes and decoded
  pixels come from `TileBuilder`'s own atomics; source residency from
  `residency(::LazyRegridArray)`; resident weight bytes from the plan's own
  `PerChunk` accounting with the cached tiles' manifest lengths beside them;
  peak RSS from `Sys.maxrss()`, which covers the whole process, warm-up
  included.
- Profile listings are self time — the frame a sample was taken in — over the
  samples whose backtrace reaches this package. Under `-t auto` the idle threads
  park in `__psynch_cvwait` and are more than nine samples in ten; both totals
  are printed so the discarded fraction is visible.

## 14. What this does not show

- **No box that crosses a latitude row.** Section 6 is why. The 2x2-degree box
  the plan called for, and every larger block, is unreachable: at 0.7 s per
  destination cell a 9.3-million-cell destination is about 75 days.
- **No full-globe or production-shaped run.** A production destination chunk is
  a level-5 subtree of the global grid, not a box, and its source is whatever
  tiles the chunk's cap reaches — which is a multi-row tile set for most chunks.
- **One machine state, one process per arm.** Every arm is a one-shot
  measurement of a phase that runs for seconds; there is no min-of-n. What is
  removed instead is compilation, by one complete miniature regrid at the same
  level with the same method before any arm runs.
- **Nothing about accuracy beyond the bit-identity check in section 10.** No arm
  compares against an independent implementation.
- **Nothing about the destination being an IGeo7 grid specifically.** The
  destination side costs 128 ns per cell here (section 3) and is not what
  dominates; a different destination system would move that number and leave the
  conclusion.

