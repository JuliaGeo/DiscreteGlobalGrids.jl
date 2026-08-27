# Copernicus DEM point stencils

`src/systems/CopernicusDEM/point.jl`, section (n) of its tests, `benchmark/copdem_point.jl`.

`BarycentricPoint` over a Copernicus DEM source used to produce an all-missing
result: the space offered neither a chart nor dual cells. It now answers from
the lattice's rows — pixels are posts, so a query interpolates between them.

## The stencil, per case

`N` is the latitude intervals per degree, `J` a post row over `0:180N-1` from
north to south, `K` a column of a row with `nc` columns to the degree. A post
stands at `lon = -180 + K/nc`, `lat = 90 - J/N`, the two pole rows at their
clamped box's midpoint `90 - 1/(4N)` and `-90 + 3/(4N)`; longitude is carried as
`x = lon + 180`, so a cell edge stands at `(2K ± 1)/(2nc)`. `JA` is the northern
of the rows bracketing the query, `a` and `b` their columns, `v` its latitude
fraction `(latA-lat)/(latA-latB)`.

| case | test | posts | coordinates |
|---|---|---|---|
| ordinary rows | `a == b` | `(JA,K)`, `(JA,K+1)`, `(JA+1,K+1)`, `(JA+1,K)`, `K = ⌊x·a⌋` | tensor Q1 in `u = x·a - K` and `v` |
| band edge, shared cell edge | `(2K-1)b == (2L-1)a` | `(JA,K-1)`, `(JA,K)`, `(JA+1,L)`, `(JA+1,L-1)` | Q1 on the trapezoid: `v`, and `u` along the parallel at `v` |
| band edge, north edge alone | `(2K-1)b > (2L-1)a` | `(JA,K-1)`, `(JA,K)`, `(JA+1,L)` | barycentric, `meanvalueweights!` on the triangle |
| band edge, south edge alone | `(2L-1)a > (2K-1)b` | `(JA,K)`, `(JA+1,L-1)`, `(JA+1,L)` | as above |
| antimeridian | none — columns wrap `mod 360nc`, chart longitudes stay unwrapped | ordinary | ordinary |
| 1-degree tile seam | none — a seam changes a row's base, not its length | ordinary | ordinary |
| holding without a post | a post carrying weight has no index in the collection | none | `WeightsRim`, empty row |
| poleward of the outermost row | `lat` outside the extreme rows' latitudes | the nearest post of that row | weight 1.0, see below |

Cell edges are compared in units of `1/(2ab)` degrees, so "shared" is integer
equality and no floating ring is intersected. The dual cell at a band edge is
found by walking to the first segment joining a north post to the post facing it
at or east of the query; the cell is what that segment and its western neighbour
bound. The walk is capped at 8 steps — the rows differ by at most one column in
two — and a query it cannot settle is degenerate, which none in the suites is.

A shared cell edge exists only where, with `g = gcd(a,b)`, both `a/g` and `b/g`
are odd. As the poleward row's columns against the equatorward row's, the five
edges are 2:3 at 50 degrees, 3:4 at 60, 2:3 at 70, 3:5 at 80 and 1:2 at 85 — so
on every Copernicus lattice the trapezoid is the 80-degree edge alone. Both
kinds are checked in both hemispheres against a brute-force oracle.

Nodes are built on the complete level and each exchanged for its index in the
collection sampled. A node carrying no weight is never looked up, so a
destination on a post maps even where its neighbours are absent. Nothing is
substituted and nothing renormalised.

## Cost and reads

One process per thread count, `res=90 box=10,11,46,47`, IGeo7 level 12 over
2,313,802 destination cells against one 1,440,000-pixel GLO-90 tile, Julia 1.12.6
on 12 CPUs. Cold reads every destination tile once on a fresh plan; warm repeats it.

| threads | method | cold ns/cell | cold cells/s | warm ns/cell | warm cells/s | cold vs nearest |
|---|---|---:|---:|---:|---:|---:|
| 8 | `NearestCell` | 192.2 | 5,203,961 | 7.6 | 132,291,280 | 1.00x |
| 8 | `BarycentricPoint` | 367.9 | 2,717,801 | 30.8 | 32,513,297 | **1.91x** |
| 8 | `Conservative` | 3,234.3 | 309,186 | 9.3 | 107,828,502 | 16.83x |
| 1 | `NearestCell` | 424.3 | 2,357,002 | 9.5 | 104,735,628 | 1.00x |
| 1 | `BarycentricPoint` | 793.1 | 1,260,829 | 35.2 | 28,372,266 | **1.87x** |
| 1 | `Conservative` | 14,449.8 | 69,205 | 19.5 | 51,173,233 | 34.06x |

Peak RSS 1.91 GiB at 8 threads and 1.71 GiB at 1, warm-up included. **A point
regrid costs about twice a nearest one, and 0.11x a conservative one at 8
threads or 0.05x at 1.** Warm it costs about four times nearest, its four
entries against nearest's one. The stencil alone, serial, one `weightsat!` per
destination site into a reused row with no read and no weight stored:

| the stencil is written in | ns per query | bytes per query | entries per placed cell |
|---|---:|---:|---:|
| the one-tile holding | 519.7 | 0.0 | 4.000 |
| the complete GLO-90 level | 219.8 | 0.0 | 4.000 |

Allocation is zero in steady state: the suite asserts `@allocated == 0` on a
warm query and the sweeps round to 0.0 B over 2.3 million. The 300 ns between
the two lines is `localindex` on the holding — four binary searches of 1,440,000
ascending ids a query, the one lever left; the other 220 ns is the site and the
row arithmetic.

GLO-30 over the same box and level, a nine times finer source, 8 threads: 500.8,
811.9 and 6,218.7 ns per cell for the three methods, at 2.61 GiB resident — a
weight block's reference vector is sized by the source chunk it names.

Per destination tile, the source chunks read equal the sorted union of the
chunks owning that tile's nonzero entries — asserted for two holdings, a
destination tiled as its source is and one destination tile inside three held
source tiles. Planning and building the lazy array read nothing, and the direct
and chunked lazy routes agree. `supportradius` is the largest post-to-point
diagonal over the six bands, doubled, checked against every post named.

## The poles

A point poleward of the outermost post row lies in the dual cell of the pole,
whose nodes are that row entire: 43,200 posts on GLO-90, 129,600 on GLO-30 —
linear in the row against a constant four everywhere else. So the query takes
the **nearest post of that row**: one entry, weight 1.0, `WeightsMapped`. The
region is `lat > 90 - 1/(4N)` and `lat < -90 + 3/(4N)`, a quarter and three
quarters of an arcsecond on GLO-30. The policy is named on the method:
`BarycentricPoint()` carries `poles = NearestCell()`; `poles = nothing` leaves
the region `WeightsDegenerate` with an empty row — a construction not settled,
not `WeightsOutside`, which would claim the source does not cover the point.

Every post of a row stands at one latitude, so the nearest in longitude is the
nearest on the sphere: the cosine of the distance grows with the cosine of the
longitude difference and with nothing else. The column is one rounding, wrapping
at the antimeridian, and the post is exchanged for its index in the collection
sampled — absent, it is `WeightsRim`, as anywhere else. The interpolant steps
where the fallback begins, which is what the constant stencil costs.

## The gate, per clause

- **every query is O(1)** — two floors and four products, at most 8 walk steps
  at a band edge, one rounding in a pole region; 0.0 B and `@allocated == 0`.
- **emits at most four entries, one in a pole region** — asserted over every
  sweep on all three lattices, and 4.000 per placed cell on the benchmark.
- **reads exactly its owning chunks** — asserted per destination tile against
  the union recomputed from the stencils.
- **agrees across seams and band transitions** — a 1-degree seam and the
  antimeridian reproduce the weights of the same query a degree west, and every
  band edge agrees with the oracle in both hemispheres.
