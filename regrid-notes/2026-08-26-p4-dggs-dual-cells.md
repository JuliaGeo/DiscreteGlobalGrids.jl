# P4 — dual cells behind a `DGGSpace`

- Date: 2026-08-26
- Plan: `regrid-notes/2026-08-23-barycentric-regridding-plan.md`, task P4
- Scope: `src/dual_cells.jl`, the chart kernel and the refusal trait in
  `lib/GlobalRegridding/src/barycentric.jl`

## 1. What was built

`BarycentricPoint()` on a `DGGSpace` source used to answer `WeightsOutside`
everywhere and produce an all-missing result without a word. It now interpolates
on the dual complex of the source grid. `samplerstate` holds the complete level
rings and sites are read on, which for a subset is its `complete`. Per query:

1. `localindex(grid, p)` locates the host cell, the point location every point
   method already pays for;
2. the host's `Vertex()` one-ring gives every cell touching it, and its `Edge()`
   one-ring marks, as bit positions in that same ring, where one primal vertex's
   fan ends and the next begins;
3. each run between consecutive marks, with the host's site in front, is a
   candidate — one per primal vertex of the host, and their union covers the
   host cell, so the one holding the point is among them.

No cell polygon is matched, no primal vertex is numbered, no table is built.
Buffers are fixed-capacity `SmallVector`s sized by `MAX_DUAL_NODES = 12` and the
ring marks are a `UInt16`, so a query touches no heap. `dualcellat` and
`weightsat!` are one search; the second also wants the reason there is none.

## 2. The chart

`TangentChart(centre)` is the azimuthal-equidistant (log map) chart, centred on
the **queried point**, not on the host cell. Chosen because:

- a node's chart radius and bearing are its true geodesic ones, so the
  coordinates read the sphere's own angles and distances where evaluated;
- the point is the origin, so a destination sitting on a source sample site is a
  node hit and reproduces that source exactly — a coordinate within `1e-12`
  radians of the centre snaps to the origin, which is what makes it exact;
- it is injective everywhere but the antipode, where gnomonic diverges at a
  quarter turn. A node past the chart's `reach`, a quarter turn by default, has
  no coordinates and the cell holding it is `WeightsDegenerate`.

Adjacent dual cells are read in the same chart at the same point, so they
partition it consistently. Raw longitude and latitude appear nowhere.

## 3. The basis rule

Every dual cell says `MeanValue`, at any node count. On the three nodes a
hexagonal source gives, mean-value coordinates *are* that triangle's barycentric
coordinates, so P1 is exact there. A four-node cell is **not** promoted to
`Bilinear` — a DGGS dual quadrilateral is no tensor-product cell — and a
five-node one, ISEA4R's poles, is a pentagon, never silently triangulated.

## 4. Cost model, and chunk invariance

`benchmark/dggs_dual_cells.jl`, `julia -t 1`, one session, IGeo7 level 5 source
(168,072 cells), 1,536 destination points off the source sites:

| arm | ns/query | bytes/query |
|---|---:|---:|
| `NearestCell`: the point location every point method pays | 313 | 0 |
| that location plus one `cell_centroid` | 384 | 0 |
| that location plus the host's two one-rings | 530 | 0 |
| `BarycentricPoint`: the dual cell and its coordinates | 1,658 | 0 |

That is 5.3× the location, and the 1,345 ns surcharge splits about evenly
between the seven sites a query reads, the two ring walks, and the chart with
the containment tests on top. `supportradius` is twice the widest chunk cover:
a cell sits inside its own chunk's cover, so twice that radius bounds its width.

A stencil is a function of the space and the point, and nothing on the path is
told about chunks. Verified rather than asserted: two chunkings of one IGeo7
level give identical rows entry for entry and an operator exactly equal to the
eager whole-domain one; eager, chunked and lazy reads agree to `1e-13`; and every
chunk pair the stencils name is inside the plan's dependency relation, so the
lazy path's exact-manifest check has nothing to refuse.

## 5. Gate

> topology lookup is bounded local work, uses no per-query polygon matching in
> production, and is invariant to chunking.

| clause | result |
|---|---|
| bounded local work | **pass** — one point location, two one-rings, at most one containment test per host vertex; 1,658 ns and 0 bytes per query |
| no per-query polygon matching | **pass** — `cell_boundary` is never called on the query path; coordinate matching appears only in the test oracle |
| invariant to chunking | **pass** — identical rows and an exactly equal operator across two chunkings, eager/chunked/lazy agreeing to `1e-13` |

## 6. What it verifies

`test/systems/crosssystem/regrid_dual.jl`, 91 assertions: every dual cell of a
whole level tallied by node count against `E - F + 2` primal vertices — IGeo7 980
triangles including its 12 pentagons, S2 90 quadrilaterals and 8 cube-corner
triangles, HEALPix 186 and 8, ISEA4R 150, 10 and 2 pentagons, all `MeanValue`; a
coordinate oracle at genuine primal vertices including every defect on four
systems; the point laws on 2,000 random points per system (all mapped,
nonnegative, rows to `1e-14`, the chart field to `1e-15`, every source site
reproducing itself, zero bytes); continuity across primal edges and face seams;
chunk invariance and the discovery superset; and a rooted `PartialGrid` whose
border cells are `WeightsRim` while no stencil names one outside it.

`lib/GlobalRegridding/test/test_barycentric.jl` gains 26: the chart's origin,
geodesic radii, subtended angle, right-handedness, reach, inference and zero
bytes; `containspoint`; and the refusal — `sampler(BarycentricPoint(), space)` on
a source reporting neither `hascellchart` nor the new `hasdualcells` throws an
`ArgumentError` naming the type rather than mapping nothing.

## 7. Open

- **The reach bound is coarse.** Twice the widest chunk cover is sound but tens
  of cell widths wide, so a chunked plan's relation is a loose superset. Reads
  stay exact — the tile manifest decides those — but a tighter bound wants either
  a cheap cell-width bound at a level or the deferred `chunksupport`.
- **The sample sites cost more than the topology.** A query reads seven
  `cell_centroid`s. A materialised centroid field behind `samplesites`, which
  `cellfield`'s `known=` already accepts, would take most of that back.
- **A5 allocates** about 200 bytes a query where every other system here is zero,
  and at level 1 declares an `Edge()` ring that is not its cells' shared-edge
  adjacency, so its dual complex there is that ring's; it is not tested here.
- **`BilinearPoint` compatibility** is untouched and undecided, as the plan says.
