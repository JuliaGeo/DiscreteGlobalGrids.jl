# P0 — the point-method baseline

- Date: 2026-08-25
- Card: `regrid-notes/2026-08-23-barycentric-regridding-plan.md`, "P0 — baseline
  and instrumentation"
- Commit: `Record the point-method baseline`

Nothing in production changed. What follows is what is now observable, and
where.

## 1. What `BilinearPoint` does today, pinned

`lib/GlobalRegridding/test/test_interpolation.jl`, testset
`BilinearPoint sampling policy`. Weights, stencil sizes and row sums are
asserted, not merely the absence of an error.

- **Interior.** A source of 20-degree cells reproduces the whole Q1 space
  `a + bx + cy + dxy` on the chart, not merely affine fields, and every interior
  destination has exactly four entries summing to one. Reproducing the cross
  term is what separates a tensor-product stencil from P1 on a triangulated
  dual quadrilateral.
- **Periodic longitude seam.** A destination 35 degrees east of the last source
  centre, on a 90-degree global lattice, interpolates across the seam to the
  first centre: longitude weights 11/18 and 7/18, four entries, row sum one.
- **Non-periodic rim.** Past the last centre of one axis the stencil drops to
  that axis' rim line and the other axis keeps interpolating: two entries
  carrying the other axis' weights (0.75/0.25 across a latitude rim, 0.35/0.65
  across a longitude rim), row sum one. Both axes past their span give one
  entry of weight one.
- **Outside coverage.** The chart coordinate is located but never tested for
  coverage. A destination 40 degrees north of the source's northmost cell, whose
  `cellat` is `nothing`, is served by the rim line with a row summing to one;
  `NearestCell`, which asks `cellat`, emits nothing for the same point. The same
  holds along a non-periodic longitude axis, on whatever branch the chart
  reports.

Decision 8 of the plan makes rims and points outside coverage unmapped instead.
These are the numbers P2 must document each intentional difference against.

## 2. The repeated work, counted

`lib/GlobalRegridding/test/test_interpolation.jl`, testset
`point lookup repeats once per candidate source chunk`. `T4PlaceCount` reports
`Points()`, emits `BilinearPoint` stencils, and records how many destination
cells each `buildweights!` call is asked to place.

One 200-cell destination tile against an 18-chunk source. The plan's relation
names `k` candidate source chunks for the tile; `k` is read from
`sourcesof(dependencies(plan), 1)`, never written down. Driving every pair the
relation names:

- one block per candidate, `calls == k`;
- **`placed == length(tile) * k`** — every pair walks every destination cell of
  the tile, so a point location is performed once per candidate chunk;
- the stencils those passes produce are `4 * length(tile)` entries in total;
- they land on fewer chunks than were searched, so some candidates are walked
  to produce nothing.

The middle assertion is the one P3 inverts to `== length(tile)`. It is a single
number in a single line.

## 3. What one large tile costs

`benchmark/point_tile_baseline.jl`. A global 0.125-degree raster,
2880 x 1440 = 4,147,200 cells in 32 chunks of 360 x 360, against one
750 x 500 = 375,000-cell destination tile spanning 60 by 40 degrees. Four arms
over the same tile: one `buildweights!` pass over the whole source index range
(the point-location floor), the `k` pair builds the plan performs, a cold tile
read on a fresh plan, and a warm read against cached weights.

Recorded numbers, the metadata they belong to, and how to read them are in the
file's own header. In summary: with `k = 4` the pair loop costs 4.40x one pass,
weight construction is 99.1% of a cold tile read, and the four source chunks are
each read exactly once. The gap to close is the repeated point location, not
source I/O and not the weight application.

Source reads are instrumented by a counting `DiskArrays` source that records
every `readblock!` range; source residency comes from
`residency(::LazyRegridArray)`; weight residency is the plan's own `PerChunk`
block accounting with `Base.summarysize` of the same blocks beside it.

The timings are absolutes belonging to one machine state, power mode included.
P3's comparison must be a same-session ratio — re-run this file on the unchanged
tree and again on the changed one, on the same machine in the same state — and
never a comparison against the numbers recorded here.
